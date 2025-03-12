# Power series approach to solving continuum backstepping kernel equations

These codes contain a MATLAB implementation of the power series-based algorithm to solving continuum kernel equations as described in the manusript "On Computation of Approximate Solutions to Large-Scale Backstepping Kernel Equations via Continuum Approximation" by Jukka-Pekka Humaloja and Nikolaos Bekiaris-Liberis. The journal version of the paper can be found [here](https://www.sciencedirect.com/science/article/pii/S0167691124002706) and a preprint is available in [arXiv](https://arxiv.org/abs/2406.13612).

The mathematical background and the documentation for the codes are desribed in the file "codesdoc.pdf.".

## Requirements

The codes require MATLAB with the Symbolic Math Toolbox.

## Usage

The algorithm for solving continuum kernel equations is implemented in the function file "kernelsolver.m". The file "solverdemo.m" contains an example of using the function. The syntax of the "kernelsover.m" function is described in the file "codesdoc.pdf".

## License

Copyright Jukka-Pekka Humaloja 2024. See LICENSE.txt for licensing information.

## Acknowledgements

Funded by the European Union (ERC, C-NORA, 101088147). Views and opinions expressed are however those of the authors only and do not necessarily reflect those of the European Union or the European Research Council Executive Agency. Neither the European Union nor the granting authority can be held responsible for them.

## Cite this work
Journal version:
```
@article{HumBek25,
	author = {J-P. Humaloja and N. Bekiaris-Liberis},
	journal = {Syst. Control Lett.},
	pages = {105982},
	title = {On computation of approximate solutions to large-scale backstepping kernel equations via continuum approximation},
	volume = {196},
	year = {2025}}
```
ArXiv preprint:
```
@Unpublished{HumBek24arxivb,
  author = {J.-P. Humaloja and N. Bekiaris-Liberis},
  note   = {{arXiv}, 2406.13612, 2024},
  title  = {On computation of approximate solutions to large-scale backstepping kernel equations via continuum approximation},
}
```
