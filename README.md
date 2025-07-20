# Boostlets_SparsityAnalysis

**1. INTRODUCTION**

This code can be used to reproduce the sparsity analysis in the paper:
  E. Zea, M. Laudato, J. Andén, "A boostlet transform for wave-based acoustic signal processing in space--time", arXiv:2403.11362, DOI: [10.48550/arXiv.2403.11362](https://doi.org/10.48550/arXiv.2403.11362), 2025. 

**2. INSTALLATION**

Download the .zip file and extract it in your preferred folder. You do not need to install files, as the dependencies are added (and removed) automatically when running the main script. 

OBS: The user must have the Wavelet Toolbox installed in Matlab. 

**3. EXAMPLE OF USAGE**

  To reproduce the sparsity analysis in Figure 6 of the manuscript, run the command: 

      > Boostlets_SparsityAnalysis;

  The output produces a figure (Figure 1) corresponding to Figure 6 in the paper. 
  A table including the l1-norm of the 20.0000 largest coefficients, corresponding 
  to Table I in the manuscript, is output in the Command Window. See below. 

This is the output of Figure 6:

![Figure 6](https://github.com/eliaszea/Boostlets_SparsityAnalysis/blob/main/Fig6.jpg)

This is what the Table looks like in the Command Window: 

![Table 1](https://github.com/eliaszea/Boostlets_SparsityAnalysis/blob/main/Table1.jpg)

**4. DEPENDENCIES**

Folders:    

- 'Datasets':          measured acoustic fields in three rooms[^1]
- 'Toolboxes':         curvelets (CurveLab-2.1.3)[^2], wave atoms (WaveAtom-1.1.1)[^3], and shearlets (Fast Finite Shearlet Transform, FFST)[^4]

[^1]: [E. Zea, “Compressed sensing of impulse responses in rooms of unknown properties and contents,” J. Sound Vib. 459, 114871 (2019)](https://doi.org/10.1016/j.jsv.2019.114871).
[^2]: [E. Candès, L. Demanet, D. Donoho, L. Ying, “Fast Discrete Curvelet Transforms,” Multiscale Modeling & Simulation 5(3), 861–899 (2006)](https://doi.org/10.1137/05064182X).
[^3]: [L. Demanet, L. Ying, 'Wave atoms and sparsity of oscillatory patterns', Appl. Comput. Harmon. Anal. 23(3), 368–387 (2007)](https://doi.org/10.1016/j.acha.2007.03.003).
[^4]: [S. Häuser, G. Steidl, “Fast finite shearlet transform: a tutorial,” ArXiv 1202.1773, 1-41 (2012)](https://arxiv.org/abs/1202.1773).

**5. RELEASE HISTORY**

	Release #1	 Version 001 	E. Zea	2024-03-13
 	Release #2	 Version 002    E. Zea  2025-07-20 	Renamed repo

**6. FEEDBACK & CONTACT INFORMATION**

Your questions, suggestions, and feedback can help improve the quality of this software. Feel free to contact me at

	Elias Zea (zea@kth.se)
	Marcus Wallenberg Laboratory for Sound and Vibration Research
	Department of Engineering Mechanics
	KTH Royal Institute of Technology
	Teknikringen 8
	10044 Stockholm, SWEDEN

**7. LEGAL INFORMATION**

Copyright 2024 Elias Zea, Marco Laudato, Joakim Andén

This software was written by Elias Zea. 

The provided code is free software. You can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version. If not stated otherwise, this applies to all files contained in this package and its sub-directories. 

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
