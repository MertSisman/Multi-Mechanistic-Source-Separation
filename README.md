# Multi-Mechanistic-Source-Separation
that maps sub-voxel brain microstructure. By combining biophysical gradient-echo signal simulations with a PyTorch deep learning model, MMSS translates standard mGRE and QSM data into specific biological biomarkers, including myelin susceptibility, iron susceptibility, myelin volume fraction, and free water fraction.

4. **Run the Script:** Execute `run_MMSS_pipeline.m`. 
   
The script will automatically handle path routing, normalize your in vivo data, call the Python deep learning environment for inference, and reconstruct the output maps back in MATLAB. Intermediate files are automatically cleaned up upon completion.

---

## Outputs

The pipeline produces a single output file, `MMSS_results.mat`, saved directly in your `data/` directory. It contains the following coregistered 3D volumetric maps:

*   **`chi_neg`**: Negative (diamagnetic) susceptibility map (myelin dominated).
*   **`chi_pos`**: Positive (paramagnetic) susceptibility map (iron dominated).
*   **`MVF`**: Myelin Volume Fraction.
*   **`fwf`**: Free Water Fraction.
*   **`theta`**: Estimated fiber orientation angle.

---

## Acknowledgments and Theoretical Background

The biophysical dictionary generation within this framework utilizes equations derived in:
> Wharton, S. and Bowtell, R., 2012. *Fiber orientation-dependent white matter contrast in gradient echo MRI.* Proceedings of the National Academy of Sciences, 109(45), pp.18559-18564.

## Citation

If you use the MMSS framework in your research, please cite our corresponding manuscript:

*(Citation details to be updated upon publication)* 


***

# Multi-Mechanistic Source Separation (MMSS)

**Multi-Mechanistic Source Separation (MMSS)** is a computational MRI framework that maps sub-voxel brain microstructure. By combining biophysical gradient-echo signal simulations with a PyTorch deep learning model, MMSS translates standard multi-echo gradient echo (mGRE) and Quantitative Susceptibility Mapping (QSM) data into specific biological biomarkers, including myelin susceptibility, iron susceptibility, myelin volume fraction, and free water fraction.

## Table of Contents
- [Prerequisites](#prerequisites)
- [Installation](#installation)
- [Repository Structure](#repository-structure)
- [Usage & Path Configuration](#usage--path-configuration)
- [Outputs](#outputs)
- [Citation](#citation)

---

## Prerequisites

To run the full MMSS pipeline, you will need:
1. **MATLAB**: For dictionary simulation, data preparation, and output reconstruction.
2. **Python 3.x**: For the deep learning Multi-Layer Perceptron (MLP) inference.
3. **MEDI Toolbox**: **Important:** All mGRE data preprocessing and QSM dipole inversion *must* be performed using the [Cornell MEDI Toolbox](http://weill.cornell.edu/mri/pages/qsm.html) prior to running the MMSS pipeline. Ensure your `iField.mat`, `QSM.mat`, and `RDF.mat` files are properly formatted.

---

## Installation

1. **Clone the repository:**
   ```bash
   git clone [https://github.com/YourUsername/MMSS.git](https://github.com/YourUsername/MMSS.git)
   cd MMSS
