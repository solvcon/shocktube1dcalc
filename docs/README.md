[![PRs Welcome](https://img.shields.io/badge/PRs-welcome-brightgreen.svg?style=flat-square)](http://makeapullrequest.com)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

# 1D shocktube caculator
This tool provides 1D Shock Tube solutions. You may output analytic solutions, or numeric solutions backed by CESE
method, *Sin-Chung Chang, “The Method of Space-Time Conservation Element and Solution Element – A New Approach for
Solving the Navier-Stokes and Euler Equations”, Journal of Computational Physics, Volume 119, Issue 2, July 1995, Pages 295-324. doi: 10.1006/jcph.1995.1137*.

## Getting Started

### Prerequisites
* [Python](https://www.python.org/downloads/)
* numpy
* scipy

## Usage
### Analytic Solution
```python
from shocktube1dcalc import solver_analytic

# by default it will create a the shock tube based on Sod's classic condition.
shocktube = solver_analytic.ShockTube()

import numpy as np
mesh = np.linspace(-0.5, 0.5, 50)

analytic_solution = shocktube.get_analytic_solution(
    mesh, t=0.4
)
```


You may customize the physical status of the shocktube via:
```python
shocktube = solver_analytic.ShockTube(rho_left=1.0, u_left=0.0, p_left=1.0, rho_right=0.125, u_right=0.0, p_right=0.1)
```


### Numeric Solution
```python
from shocktube1dcalc import cese

elapsed_time = 0.4
cese_grid_size_t = 0.004
# multiply 2 for half grids, so total iteration number should be double
# the iteration number is always less than 1 by the grid number
iteration_number = round(elapsed_time / cese_grid_size_t * 2) - 1

shocktube = cese.ShockTube(iteration=iteration_number, grid_size_t=cese_grid_size_t)
shocktube.run_cese_iteration()

numeric_solution = shocktube.data.solution
```

## Contributing
See [Contributing](contributing.md)

## Authors
Taihsiang Ho (tai271828) <tai271828@solvcon.net>


Created from [Lee-W/cookiecutter-python-template](https://github.com/Lee-W/cookiecutter-python-template/tree/0.7.1) version 0.7.1
