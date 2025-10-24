import time

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.widgets as widgets
import os
import rasterio
import tkinter as tk
import tkinter.messagebox as messagebox
from tkinter import filedialog
from tkinterdnd2 import DND_FILES, TkinterDnD
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

# Helper to load hyperspectral data

def nan_median(a, axis=None):
    """
    Compute the median of an array along a given axis, ignoring NaNs.

    Parameters
    ----------
    a : array_like
        Input array or object that can be converted to an array.
    axis : int or None, optional
        Axis along which the medians are computed. The default (None)
        is to compute the median of the flattened array.

    Returns
    -------
    median : float or ndarray
        Median of the array elements along the specified axis, ignoring NaNs.
    """
    a = np.asanyarray(a)

    # Mask NaNs
    mask = np.isnan(a)
    if not np.any(mask):
        # no NaNs -> normal median
        return np.median(a, axis=axis)

    # For each axis slice, select non-nan values
    if axis is None:
        # flatten and drop NaNs
        return np.median(a[~mask])
    else:
        # Move target axis to first position for easier looping
        a = np.moveaxis(a, axis, 0)
        mask = np.moveaxis(mask, axis, 0)
        result = np.full(a.shape[1:], np.nan)

        for idx in np.ndindex(result.shape):
            valid = a[(~mask[:, idx]), idx]
            if valid.size > 0:
                result[idx] = np.median(valid)
        return result



def load_hyperspectral_image(path):
    ext = os.path.splitext(path)[1].lower()
    wavelengths = None
    if ext == '.npy':
        data = np.load(path)
    elif ext == '.img':
        with rasterio.open(path) as src:
            data = src.read()  # (bands, height, width)
            data[data >= 1] = np.nan  # Handle no-data values
            wavelengths = []
            for wl in src.descriptions:
                wavelengths.append(float(wl.split(" ")[0]))
    else:
        raise ValueError(f"Unsupported file format: {ext}")
    return data, wavelengths


# Main function
class HyperspectralViewer(TkinterDnD.Tk):
    def __init__(self):
        super().__init__()
        self.title('Hyperspectral Movie Viewer')
        self.geometry('900x700')
        self.drop_label = tk.Label(self, text='Drop a hyperspectral image file here', font=('Arial', 14))
        self.drop_label.pack(pady=10)
        # Register drop only if tkinterdnd2 is available
        if DND_FILES is not None:
            self.drop_label.drop_target_register(DND_FILES)
            self.drop_label.dnd_bind('<<Drop>>', self.on_drop)
        else:
            # Fallback: bind a left-click to open file dialog
            def _open_dialog(event=None):
                path = filedialog.askopenfilename(title='Select hyperspectral image', filetypes=[('ENVI img', '*.img'), ('NumPy', '*.npy'), ('All files', '*.*')])
                if path:
                    self.load_image(path)
            self.drop_label.bind('<Button-1>', _open_dialog)
        self.fig, self.ax = plt.subplots()
        self.fig.subplots_adjust(bottom=0.2)
        self.canvas = FigureCanvasTkAgg(self.fig, master=self)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        # Ensure the canvas has keyboard focus so Left/Right work
        self.canvas.get_tk_widget().focus_set()
        # Bind left/right arrow keys on the canvas widget so arrows work when it has focus
        widget = self.canvas.get_tk_widget()
        widget.bind('<Left>', self.on_left)
        widget.bind('<Right>', self.on_right)
        self.slider = None
        self.text = None
        self.data = None
        self.wavelengths = None
        self.bands = None
        self.im = None
        self.current_band = 0
        self.slider_is_touched = False

    def run_animation(self, sleep_time=0.5):
        for band in range(self.bands):
            self.step_band(1)
            time.sleep(sleep_time)

    def on_drop(self, event):
        path = event.data.strip().replace('{', '').replace('}', '')
        self.load_image(path)
        # start non-blocking animation (0.5s per frame by default)
        self.start_animation(sleep_time=0.1)

    def start_animation(self, sleep_time=0.1):
        """Start a non-blocking animation using Tk.after."""
        if self.data is None or self.bands is None:
            return
        self.animation_delay_ms = int(sleep_time * 1000)
        self.animating = True
        # Start from current band (or set self.current_band = -1 to start at 0)
        self._animation_step()

    def stop_animation(self):
        """Stop the running animation."""
        self.animating = False

    def _animation_step(self):
        """One animation step; schedules the next via after()."""
        if not getattr(self, 'animating', False):
            return
        # advance one band; step_band will update display/slider
        self.step_band(1)
        # if we've reached the last band, stop; remove this check to loop instead
        if self.current_band >= (self.bands - 1):
            self.animating = False
            return
        # schedule next step without blocking the GUI
        self.after(self.animation_delay_ms, self._animation_step)


    def load_image(self, path):
        self.ax.clear()
        self.data, self.wavelengths = load_hyperspectral_image(path)
        if self.data.ndim == 3:
            self.bands = self.data.shape[0]
        elif self.data.ndim == 2:
            _, _, self.bands = self.data.shape
            self.data = np.transpose(self.data, (2, 0, 1))
        else:
            messagebox.showerror('Error', 'Unsupported data shape. Expected 3D array.')
            return
        if not self.wavelengths or len(self.wavelengths) != self.bands:
            self.wavelengths = [f"Band {i+1}" for i in range(self.bands)]
        self.im = self.ax.imshow(self.data[0], cmap='gray', vmin=0, vmax=1)
        self.text = self.ax.text(0.5, 0.95, '', transform=self.ax.transAxes, ha='center', va='top', fontsize=10, color='red',
                                 bbox=dict(facecolor='black', alpha=0.6, boxstyle='round,pad=0.1'))
        self.ax.axis('off')
        self.setup_slider()
        self.update_band(0)
        self.canvas.draw()

        self.save_frames(outdir=None, overwrite=False)

    def setup_slider(self):
        if self.slider:
            self.slider.ax.clear()
        ax_slider = self.fig.add_axes([0.2, 0.08, 0.6, 0.04])
        self.slider = widgets.Slider(ax_slider, 'Band', 1, self.bands, valinit=1, valstep=1)
        self.slider.on_changed(self.on_slider_change)

    def on_left(self, event):
        # step one band lower
        self.step_band(-1)
        self.stop_animation()

    def on_right(self, event):
        # step one band higher
        self.step_band(1)
        self.stop_animation()

    def step_band(self, delta):
        """Move current band by delta (±1). Updates slider if present."""
        if self.data is None or self.bands is None:
            return
        new_band = int(np.clip(self.current_band + delta, 0, self.bands - 1))
        if new_band == self.current_band:
            return
        self.current_band = new_band
        if self.slider:
            # Slider in this UI is 1-based
            try:
                self.slider.set_val(new_band + 1)
            except Exception:
                # fallback: update directly
                self.update_band(new_band)
                self.canvas.draw()
        else:
            self.update_band(new_band)
            self.canvas.draw()

    def on_slider_change(self, val):
        band = int(self.slider.val) - 1
        self.update_band(band)
        self.canvas.draw()

    def save_frames(self, outdir=None, overwrite=False):
        """
        Save every band as a PNG into a `bands_single` folder.
        Each PNG receives a small red wavelength label in the image.
        """
        if self.data is None or self.bands is None:
            messagebox.showwarning("No data", "No image loaded to export.")
            return

        base = os.path.dirname(getattr(self, 'current_path', '')) or os.getcwd()
        if outdir is None:
            outdir = os.path.join(base, 'bands_single')
        os.makedirs(outdir, exist_ok=True)

        for i in range(self.bands):
            band = self.data[i]
            band_min = np.nanmin(band)
            band_center = np.nanmedian(band)
            normed = band - band_min
            if band_center != 0 and not np.isnan(band_center):
                normed = normed / band_center
            arr = np.nan_to_num(normed, nan=0.0)

            p1, p99 = np.nanpercentile(arr, [1, 99])
            if p99 > p1:
                arr = (arr - p1) / (p99 - p1)
            else:
                arr = np.clip(arr, 0.0, 1.0)
            arr = np.clip(arr, 0.0, 1.0)

            fname = os.path.join(outdir, f'band_{i + 1:03d}.png')
            if os.path.exists(fname) and not overwrite:
                continue

            # Prepare wavelength label
            wl = None
            if isinstance(self.wavelengths, (list, tuple, np.ndarray)) and i < len(self.wavelengths):
                wl = self.wavelengths[i]
            wl_text = f"{wl} nm" if isinstance(wl, (int, float)) else (str(wl) if wl is not None else f"Band {i + 1}")

            # Create figure sized to image (keep reasonable DPI)
            h, w = arr.shape
            dpi = 100
            fig = plt.figure(figsize=(w / dpi, h / dpi), dpi=dpi)
            ax = fig.add_axes([0, 0, 1, 1])
            ax.imshow(arr, cmap='gray', vmin=0, vmax=1, aspect='auto')
            ax.axis('off')

            # Small red wavelength label in the lower-right corner
            ax.text(0.98, 0.02, wl_text, transform=ax.transAxes,
                    ha='right', va='bottom', color='red', fontsize=8)

            # Save and close
            fig.savefig(fname, dpi=dpi, bbox_inches='tight', pad_inches=0)
            plt.close(fig)

    def update_band(self, frame):
        # Normalize and display the requested band. Keep current_band in sync.
        frame = int(frame)
        if self.bands is not None:
            frame = int(np.clip(frame, 0, self.bands - 1))
        self.current_band = frame
        band_data = self.data[frame]
        normed = (band_data - np.nanmin(band_data)) / np.nanmedian(band_data)
        self.im.set_data(normed)
        self.im.set_clim(np.nanpercentile(normed, 1), np.nanpercentile(normed, 99))
        wl = self.wavelengths[frame] if isinstance(self.wavelengths, (list, np.ndarray)) else f"Band {frame+1}"
        if self.text is None:
            self.text = self.ax.text(0.5, 0.95, '', transform=self.ax.transAxes, ha='center', va='top', fontsize=10, color='red', bbox=dict(facecolor='black', alpha=0, boxstyle='round,pad=0.3'))
        self.text.set_text(f"Band {frame+1} | Wavelength: {wl}nm")

if __name__ == "__main__":
    app = HyperspectralViewer()
    app.mainloop()
