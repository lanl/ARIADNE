# ARIADNE
O# (O5151) 

The code ARIADNE performs uncertainty quantification of experimental data for nuclear data. It also performs nuclear data evaluations with the algorithm generalized squares.

Experimental uncertainty quantification is performed for experimental data that is usually pulled from the EXFOR database. The EXFOR database is an open-access database with
nuclear reaction experiments that have been measured. EXFOR data are used as examples for demonstrating ARIADNE. These data are published under a CC BY 4.0 license (https://nds.iaea.org/nrdc/exfor-master/). ARIADNE ingests these data together with user-defined information on uncertainties that are either retrieved from the open-access database EXFOR or the open literature. It produces as output total covariances for these experimental data based on uncertainty quantification algorithms proposed in the open literature. 

Nuclear data evaluations are performed based on user-defined input; that is experimental data and covariances provided by ARIADNE and user-defined prior input. Prior input is
either produced via the LANL open-source codes CoH3 or CGMF. It produces evaluated nuclear data that are then released to U.S. open-access nuclear data libraries such as ENDF/B-VII.0 and later ENDF/B libraries. These Evaluated Nuclear Data File (ENDF/B) nuclear data libraries are publicly released (open-access and with unlimited distribution) by the NNDC, BNL (https://www.nndc.bnl.gov/endf-releases/).

© 2026. Triad National Security, LLC. All rights reserved.

## Authors
D. Neudecker, N.A.W. Walton, A. Khatiwada

## Requirements
- python >=3.9,
- numpy>=1.26,
- scipy>=1.13.1,
- jupyter,
- lxml>=5.2,
- matplotlib>=3.8.

## Installation
- Download or cloned ARIADNE,
- cd ARIADNE,
- Open the jupyter example notebook: /Examples/LoadReferenceData/LoadNewData.ipynb and run it. That installs the reference library on your computer necessary for uncertainty quantification.
- Try running the example python notebooks in the "Examples/" folder to see how the notebooks work.

## License 
This program is Open-Source under the BSD-3 License.

Copyright (c) 2026, Los Alamos National Laboratory
Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
- Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
- Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
- Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

 

(End of Notice)
