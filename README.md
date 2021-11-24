# Communication coordination assumptions in network controllability

[![License](https://img.shields.io/badge/License-Apache_2.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)

<!-- [![arXiv](https://img.shields.io/badge/arXiv-1234.56789-b31b1b.svg)](https://arxiv.org/abs/1234.56789) -->

Code for the paper on "Communication coordination assumptions in network controllability", van den Heuvel and Nys (2021)

## Content

* controllability_analyses.ipynb: reproduces all the analyses in the paper. Note that if fully run, it will remake all the files that are present in the 'simulation_output' folder in the main folder. Since these simulations take a long time to complete, we have included our simulation-output in the folder 'simulation_output'. So if you want to reproduce the exact figures from the paper, move the simulation files into the same folder as the notebook.
* folder 'data' includes the network files for the: 
  * consultancy network (source: http://opsahl.co.uk/tnet/datasets/Cross_Parker-Consulting_info.txt)
  * Ownership-USCorp network (source: http://vlado.fmf.uni-lj.si/pub/networks/data/econ/Eva/Eva.htm)


## How to use

To run the example in controllability_analyses.ipynb
* install anaconda, networkx, and graph-tool (https://git.skewed.de/count0/graph-tool/-/wikis/installation-instructions)
* folder 'modules' contains the code to find maximal matching, which is used to calculate the driver node fraction n<sub>D</sub>
* run cells in controllability_analyses.ipynb

## Cite

```python
@article{heuvel2021communication,
      title={Communication coordination in network controllability}, 
      author={Milan van den Heuvel and Jannes Nys},
      year={2021},
      eprint={2105.04164},
      archivePrefix={arXiv}
}
```
