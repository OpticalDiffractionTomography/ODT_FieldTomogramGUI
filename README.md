# ODT_FieldTomogramGUI

## Overview

**ODT_FieldTomogramGUI** is a MATLAB-based graphical user interface (GUI) designed for **field retrieval** and **3D tomogram reconstruction** in **optical diffraction tomography (ODT)**. It enables users to visualize and reconstruct refractive index distributions from off-axis holographic measurements with minimal programming effort.

---

## System Requirements

### Hardware
- A computer equipped with an **NVIDIA GPU** is recommended to enable GPU-accelerated tomogram reconstruction.
- Minimum requirement: GPU with **4 GB of memory**.

### Software
- MATLAB **R2018a** or later
- Required toolboxes:
  - **Parallel Computing Toolbox**
  - **Image Processing Toolbox**

## Installation Guide

1. Download and place all files (`Field_Tomogram_Reconstruction.fig`, `Field_Tomogram_Reconstruction.m`, and auxiliary scripts) into the same directory.
2. In MATLAB, set the path to this directory using:
   ```matlab
   addpath('path_to_your_folder');
   ```
3. Launch the GUI by running the main `.m` script (e.g., `Field_Tomogram_Reconstruction.m`).

## Input Data

The software requires at least two `.mat` files located in the same folder:

- `sample*.mat`: Off-axis holograms of samples (use consistent `sample` naming).
- `bg*.mat`: Background (empty field-of-view) holograms.

Each `.mat` file should include the following variables required for field retrieval and tomogram reconstruction:
- `res`: Pixel resolution (in µm). This is the image-plane resolution, calculated as camera pixel size divided by total magnification.
- `NA`: Numerical aperture of the detection objective lens.
- `lambda`: Illumination wavelength (in µm).
- `tomogMap`: A 3D array of raw holograms (x, y, θ) from 150 illumination angles. The hologram corresponding to normal illumination should be at index 49.
 
---

## How to Use

The workflow consists of three steps:

### 1. Sample Selection (Gray Box)

1. Click **"Browse"** to select the folder containing the sample and background `.mat` files.
2. The surrounding medium’s refractive index is preset to **1.337** by default. Update this value if using a different medium.
3. If multiple background files are available, you can select one using the **"Background #"** field.

### 2. Field Retrieval (Red Box)

1. Choose the **Analysis Option** to retrieve the field from either a single hologram or all sample holograms.
2. Click **"Retrieval"** to begin. The progress is displayed in real time.
3. Once completed:
   - The phase delay map at normal incidence is shown.
   - Retrieved fields are saved as `Field_*.mat` in the `field_retrieval/` folder.
   - Six representative phase maps and a central cross-section from 150 angles are also saved as `Field_*.png`.

### 3. Tomogram Reconstruction (Blue Box)

1. After field retrieval, proceed to tomogram reconstruction.
2. Use **"Field Inspection"** to examine and optionally exclude corrupted frames (e.g., those affected by dust or artifacts).
   - Set thresholds for:
     - **Phase value** (red line)
     - **Phase gradient** (green line)
   - Frames exceeding these thresholds will be excluded for tomogram reconstruction.
3. Click **"Reconstruct"** to begin reconstruction. Progress is displayed live.
4. Once completed:
   - A cross-sectional slice (x-y, y-z, or x-z) of the latest sample is displayed.
   - Navigate through different planes using the GUI.
   - Tomograms are saved as `Tomogram_*.mat` in the `field_retrieval/` folder.
   - Cross-sectional slices and maximum projection images are also saved.
5. Export to TIFF:
   - Click **"Export TIFF"** to generate 16-bit stackable TIFF files readable by ImageJ/FIJI.
   - Refractive index (RI) values are scaled by 10,000 (e.g., RI 1.344 → value 13440).
   - Axial and lateral resolutions differ and can be checked from the variables `res3` (lateral) and `res4` (axial) in the `.mat` file.
   - Because the voxel size is anisotropic, be sure to adjust scaling properly in ImageJ/FIJI.

---

## Citation

If you use **ODT_FieldTomogramGUI** in your research or publication, please cite it as:

> Kyoohyun Kim (2025). *ODT_FieldTomogramGUI version X.X.X: MATLAB GUI for Field Retrieval and Tomogram Reconstruction in Optical Diffraction Tomography* [Software]. Available at: https://github.com/OpticalDiffractionTomography/ODT_FieldTomogramGUI
