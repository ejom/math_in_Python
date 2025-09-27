# Math in Python: ODE Analysis

This project demonstrates numerical analysis of ordinary differential equations (ODEs) using Python's scientific computing libraries. It solves three different ODEs and compares numerical solutions with analytical solutions where available.

## Overview

The project analyzes three different ODEs:

1. **Damped Oscillator**: `y'' + y' + y = 0`
2. **Airy Equation**: `y'' - t*y = 0`
3. **Third-order Nonlinear ODE**: `(t²+1)y''' + ty'' + t(y')² + y = 0`

## Features

- Numerical solution of ODEs using `scipy.integrate.solve_ivp`
- Comparison with analytical solutions where available
- Error analysis and visualization
- High-precision calculations with custom tolerances
- Jacobian matrix implementation for improved numerical stability

## Requirements

The project requires the following Python packages (see `requirements.txt`):

- `numpy` (2.3.3) - Numerical computing
- `scipy` (1.16.2) - Scientific computing and ODE solving
- `matplotlib` (3.10.6) - Plotting and visualization
- `sympy` (1.14.0) - Symbolic mathematics

## Installation

1. Clone or download this repository
2. Install the required dependencies:

```bash
pip install -r requirements.txt
```

## Usage

Run the main analysis script:

```bash
python ODE_analyze.py
```

## License

This project is for educational purposes. Feel free to use and modify as needed.
