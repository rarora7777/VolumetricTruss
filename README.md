This repository accompanies the SCF 2019 paper "**Volumetric Michell Trusses for Parametric Design & Fabrication**".

The code takes in as input a) a 3D tetrahedral mesh as input, and b) Dirichlet (fixed points) and Neumann (static loads) boundary conditions, and produces a 3D truss structure whose elements follow the stress field induced by the given boundary conditions.

The solution involves solving for a global **R**<sup>3</sup>-valued parametrization on the input mesh, whose isolines are traced to generate the resulting truss. In effect, it means that the truss is composed of end-to-end curves from three mutually-orthogonal familes of curves. Please look at the paper for details. [https://doi.org/10.1145/3328939.3328999](https://doi.org/10.1145/3328939.3328999)

# Installation and dependencies

Simply clone the repository to install.

`git clone https://github.com/rarora7777/VolumetricTruss`

The provided MATLAB code depends on the following packages.
1. [GAUSS](https://github.com/dilevin/GAUSS)
2. [gptoolbox](https://github.com/alecjacobson/gptoolbox/)
3. MATLAB. We have primarily tested with the version R2018b.

# License

This code is available under the MIT license.

# Usage

Please see `batchResults.m` for example usage.

# Citation

If you utilize this code or dataset for a publication, please cite

```
@inproceedings{Arora:michell:scf:2019,
 author = {Arora, Rahul and Jacobson, Alec and Langlois, Timothy R. and Huang, Yijiang 
 and Mueller, Caitlin and Matusik, Wojciech and Shamir, Ariel and 
 Singh, Karan and Levin, David I.W.},
 title = {Volumetric Michell Trusses for Parametric Design \& Fabrication},
 booktitle = {Proceedings of the 3rd ACM Symposium on Computation Fabrication},
 series = {SCF '19},
 year = {2019},
 location = {Pittsburgh, PA, USA},
 numpages = {13},
 publisher = {ACM},
 address = {New York, NY, USA}
}
```
