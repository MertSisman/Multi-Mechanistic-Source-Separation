# Multi Mechanistic Source Separation (MMSS)

**Multi Mechanistic Source Separation (MMSS)** is a computational MRI framework that maps sub voxel brain microstructure. By combining biophysical gradient echo signal simulations with a PyTorch deep learning model, MMSS translates standard multi echo gradient echo (mGRE) and Quantitative Susceptibility Mapping (QSM) data into specific biological biomarkers, including myelin susceptibility, iron susceptibility, myelin volume fraction, and free water fraction.

## Table of Contents
* Prerequisites
* Installation
* Repository Structure
* Usage & Path Configuration
* Outputs
* Citation

***

## Prerequisites

To run the full MMSS pipeline, you will need:
1. **MATLAB**: For dictionary simulation, data preparation, and output reconstruction.
2. **Python 3.x**: For the deep learning Multi Layer Perceptron (MLP) inference.
3. **MEDI Toolbox**: **Important:** All mGRE data preprocessing and QSM dipole inversion *must* be performed using the [Cornell MEDI Toolbox](http://weill.cornell.edu/mri/pages/qsm.html) prior to running the MMSS pipeline. Ensure your `iField.mat`, `QSM.mat`, and `RDF.mat` files are properly formatted.

***

## Installation

1. **Clone the repository:**
   ```bash
   git clone https://github.com/MertSisman/Multi-Mechanistic-Source-Separation.git
   cd MMSS
   ```
Install PyTorch:
Because PyTorch installation depends heavily on your specific hardware (CPU vs. CUDA/GPU), we highly recommend installing it first via the official instructions at pytorch.org.

Install remaining Python dependencies:
Once PyTorch is installed, install the rest of the required Python libraries using the provided requirements file:

   ```bash
   pip install numpy scipy mat73 matplotlib
   ```
Repository Structure
```Plaintext
MMSS/
* data/                   # Place your preprocessed MEDI outputs here (iField, QSM, RDF)
* dictionary/             # Pre simulated dictionary (.mat files)
* scripts/                
  * run_MMSS_pipeline.m   # Main executable wrapper for the entire pipeline
* src/                    
  * simulation/           # Biophysical field modeling (compute_field, estimate_chi_neg, etc.)
  * utils/                # Helper functions (dipole_kernel, polar_mesh)
  * ml/                 
    * train_and_process.py  # PyTorch deep learning architecture and training loop
* requirements.txt        # Python dependencies
* README.md
```
Usage & Path Configuration
The MMSS pipeline is designed to be highly automated. You only need to configure the root paths in the main execution script.

Place your MEDI preprocessed data (iField.mat, QSM.mat, RDF.mat) into the data/ directory.

2. Open scripts/run_MMSS_pipeline.m in MATLAB.

3. Configure your paths: At the very top of the script, update the REPO_PATH and PYTHON_EXEC variables to match your local machine. This is the only code you need to change.

```Matlab
%% ========================================================================
%  USER CONFIGURATION
%  ========================================================================
% 1. Define the absolute path to the main MMSS repository folder.
REPO_PATH = '/path/to/your/MMSS'; 
% 2. Define the path to your Python executable. 
% (e.g., '/home/username/miniconda3/envs/MMSS_env/bin/python3' or simply 'python')
PYTHON_EXEC = 'python3';

```


4. **Run the Script:** Execute `run_MMSS_pipeline.m`. 
   
The script will automatically handle path routing, normalize your in vivo data, call the Python deep learning environment for inference, and reconstruct the output maps back in MATLAB. Intermediate files are automatically cleaned up upon completion.

***

## Outputs

The pipeline produces a single output file, `MMSS_results.mat`, saved directly in your `data/` directory. It contains the following coregistered 3D volumetric maps:

* **`chi_neg`**: Negative (diamagnetic) susceptibility map (myelin dominated).
* **`chi_pos`**: Positive (paramagnetic) susceptibility map (iron dominated).
* **`MVF`**: Myelin Volume Fraction.
* **`fwf`**: Free Water Fraction.
* **`theta`**: Estimated fiber orientation angle.

***

## Acknowledgments and Theoretical Background

The biophysical dictionary generation within this framework utilizes equations derived in:
> Wharton, S. and Bowtell, R., 2012. *Fiber orientation dependent white matter contrast in gradient echo MRI.* Proceedings of the National Academy of Sciences, 109(45), pp.18559-18564.

## Citation

If you use the MMSS framework in your research, please cite our corresponding manuscript:

*(Citation details to be updated upon publication)* 

**Author:** Mert Sisman  
**Institution:** Weill Cornell Medicine
