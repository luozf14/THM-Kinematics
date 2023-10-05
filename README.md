# Kinematics calculation for THM
This repo provides three scripts to calculate the knematics for Trojan Horse Method. 

## Prerequisites
- ROOT6
- Bound state wave function from FRESCO/RADCAD/...

## How to use
### GetPsi_p.C
This script converts the bound state wave function obtained from FRESCO from coordinate representation to momentum representation. It will generate ``Psi_p2.txt``.

### GetAngCorrel.C
This script reads the ``Psi_p2.txt``. It will draw the energy vs angle and angle vs angle figures. 

### GetAngCorrelMT.C
Almost the same as ``GetAngCorrel.C`` except it is multithreads. The number of threads used is defined in ``nWorkers``.

### Analysis
The output data file is ``myFile.root``.
