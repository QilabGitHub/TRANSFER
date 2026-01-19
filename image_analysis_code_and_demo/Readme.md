# Colocalization Analysis Script

This MATLAB script calculates a colocalization score between fluorescent
channels from multi-channel TIFF microscopy images.

The MATLAB code has been tested on MATLAB 9.13.0.2105380 (R2022b) Update 2.

---

## Required Input Files

### 1. Image files
- All `.tif` images must be placed in the **same folder**
- All images must share the **same channel acquisition order**

### 2. Channel sequence file
- An Excel file named: Channel_sequence.xlsx

---

## ROI Selection

- The user will be prompted to manually select a **region of interest (ROI)**
representing the **recipient cell** based on nuclear fluorescent protein signals.
- Donor cells labeled by a dye of different color should be excluded from the ROI.
- All downstream analysis is restricted to the selected ROI.

---

## Thresholding and Mask Generation

- Intensity thresholds are calculated based on pixel-intensity statistics
within the selected ROI.
- Binary masks are generated as follows:

pixel intensity < threshold → 0
pixel intensity ≥ threshold → 1


---

## Colocalization Definition

- Pixels overlapping with the **binary mask of endosome markers** are defined
  as **colocalized pixels**.

---

## Output Files

### 1. ROI mask file

Cell_roi_mask.mat

- Stores all manually selected ROIs
- Automatically loaded if the script has been run previously

---

### 2. Endosome mask file

OMM_mask.mat

- endosome binary masks
- corresponding image filenames

---

### 3. Colocalization score file

colocalization.csv (or intensity.csv if analyzing intensity on the mask as well)

- The first column (Filenames) specifies the image names.
- For colocalization score analysis, the second column (**RNA_colocal**) specifies the percentage of mCherry or GFP signal co-localized with the endosome marker.
- For analysis of integrated intensity of fluorescent signal on the mCherry mask, the column (**Intensity_RNA**) specifies the total intensity of Sephluorin on the mask,
- and the column (**Intensity_OMM**) specifies the total intensity of mCherry on the mask.

---

## Re-running the Analysis

If you wish to test different threshold parameters:

- Simply rerun the script
- Do not modify the saved ROI files

---

## Demo

Three folders:
### 1. Quantifying mCherry on endosome marker
### 2. Quantifying GFP on endosome marker
### 3. Quantifying fluorescence intensity in recipient cells

Each containing Demo data (.tif files) under the Demo folder, MATLAB code for analysis, and Channel intensity file.
Example output files can be found in the Demo folder as well.

