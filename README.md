# Energy Finance

## Abstract
This study presents a framework for modeling electricity prices using a Geometric model incorporating an Ornstein-Uhlenbeck (OU) process driven by a Variance Gamma jump process. 

The spot price dynamics are defined by mean-reverting stochastic processes to capture market behaviors. Forward price expressions under the risk-neutral measure are derived using the Esscher transform, and swap prices are estimated based on discrete forward prices over the delivery period. 

Model parameters are calibrated using market data, employing Mean Squared Error minimization via MATLAB’s optimization tools. 
The impact of standardization on calibration accuracy is examined, revealing that normalization improves parameter estimation and reduces error metrics. 

The study also explores an alternative model using only a Gaussian OU process, demonstrating its limitations in capturing electricity price volatility. 

Finally, implied volatility curves are computed, showing characteristics consistent with commodity markets. 
The results suggest that incorporating jump processes enhances model accuracy, making it more suitable for energy finance applications.

## Repository Structure
* [Report](Report.pdf)
  contains a deeper description of the study and a detailed explanation of the methodologies.

* **`GeometricModel_GAUSS_MAIN.m`** contains the main script that executes the core logic of the project involving the Gaussian OU process.

* **`GeometricModel_VG_MAIN.m`** contains the main script that executes the core logic of the project involving the OU model driven by a Variance Gamma jump process.
  
* **`Data/`** contains the original dataset used in the project.

* **`Functions/`** includes auxiliary functions for processing, analysis, and other operations.


