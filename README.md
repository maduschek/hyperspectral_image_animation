# Hyperspectral Image Animation

Simple Tkinter GUI to view hyperspectral image stacks, step bands, run a non-blocking animation, and export single-band PNGs with a small red wavelength label.

## Files
- `hs_animation.py` — main viewer script.
- Exports saved to `bands_single` (created next to the opened file by default).

## Requirements
- Python 3
- Python packages: `numpy`, `matplotlib`, `rasterio`, `tkinter` (system), `tkinterdnd2` (optional)

Install common Python deps:
`pip install numpy matplotlib rasterio`

Install optional drag\-and\-drop support:
`pip install tkinterdnd2`

## Usage
Run the viewer:
`python3 hs_animation.py`

Load a file:
- If `tkinterdnd2` is available: drag & drop a raster (e.g., ENVI `.img`) onto the label.
- Otherwise: click the label to open a file dialog.

Supported input formats must be readable by `rasterio` (e.g., ENVI `.img`, GeoTIFF). The script also accepts NumPy `.npy` arrays.

## Controls
1. Left / Right arrow keys: step one wavelength bin lower / higher (also stops any running animation).
2. Slider: jump to a specific band.
3. Press `s` (if bound in the UI) or use the export routine to save frames.

## Animation
- Animation runs non\-blocking via Tkinter's `after()` so the GUI stays responsive.
- On file load the script may optionally start an animation stepping frames at a configurable delay.

## Export
- Each band can be saved as a PNG into `bands_single` (next to the opened file or current working directory).
- Exported images are normalized using the display normalization (1st/99th percentile) and include a small red wavelength label.
- Files are named `band_001.png`, `band_002.png`, etc. Existing files are skipped unless overwritten.

## Notes & Troubleshooting
- If `tkinterdnd2` is not installed you still can open files via the click dialog.
- Ensure `rasterio` can open your file if you expect wavelengths/metadata to be read.
- If arrow keys do not work, click the image/canvas so it has keyboard focus.
