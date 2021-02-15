[![PRs Welcome](https://img.shields.io/badge/PRs-welcome-brightgreen.svg?style=flat-square)](http://makeapullrequest.com)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

# 1D shocktube caculator
This tool provdes 1D Shock Tube analytic solutions.

## Getting Started

### Prerequisites
* [Python](https://www.python.org/downloads/)
* numpy
* scipy

## Usage
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

## Contributing
See [Contributing](contributing.md)

## Authors
Taihsiang Ho (tai271828) <tai271828@gmail.com>


Created from [Lee-W/cookiecutter-python-template](https://github.com/Lee-W/cookiecutter-python-template/tree/0.7.1) version 0.7.1
