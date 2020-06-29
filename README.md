# Microdialysis (qµD)
![qµD simulation showing chamber compound concentration as a function of K<sub>D</sub>](2020-06-26_Figure1-KDvsConcentrations_SS1.svg)

Functions for simulation of quantitative microdialysis (qµD) experiments, as described and used in the recently submitted MDPI Methods and Protocols manuscript 
:

"*Quantitative Microdialysis for Rapid Affinity Determination of Small
Molecules to Target Proteins and for Exclusion of Compounds with Poor
Physicochemical Properties.*"

## Installation

The main simulation code is contained within microdialysis_equations.py. It is not set up to be a Python module, but rather, useable within different projects without installation.

## Programs
Whilst microdialysis_equations.py contains code to integrate simulations into custom processes, the following demonstration programs are available, and also produce the plots used in the submitted publication.
#### 01_simulate_qud_concentrations.py
This program allows single point simulation of qµD experiments.  Looking through the code should allow the reader to familiarise themselves with the available qµD functions.

#### 01_simulate_qud_concentrations.py

#### 02_plot_KDvsConcentrations.py
![qµD simulation showing chamber compound concentration as a function of K<sub>D</sub>](2020-06-26_Figure1-KDvsConcentrations_SS1.svg)
Program to simulate compound concentration in the red and white chambers as a function of K<sub>D</sub>.
#### 03_deriveKD_from_lred.py
#### 03_deriveKD_from_lwhite.py
#### 03_deriveKD_from_pt.py
#### 04_deriveKD_from_multiple_lred.py
#### 04_deriveKD_from_multiple_lwhite.py
#### 04_deriveKD_from_multiple_pt.py
#### 05_plot_KDvsPt.py
![qµD simulation showing p<sub>t</sub> as a function of K<sub>D</sub>](2020-06-26_Figure2-KDvsPt_SS1.svg)
Program to simulate *p<sub>t</sub>* as a function of K<sub>D</sub>.
#### 06_calc_derivatives.py
#### 07_plot_conc_to_kd_accuracy.py

#### microdialysis_equations.py
![qµD simulation showing 5% readout error over different K<sub>D</sub> ranges](2020-06-26_Figure3-KDvsPtWith5pctError_SS1.svg)

## Contributing
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

## License
[MIT](https://choosealicense.com/licenses/mit/)