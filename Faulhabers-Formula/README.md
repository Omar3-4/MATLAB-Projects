# Faulhaber's Formula: Sum of Powers

This folder contains a custom MATLAB function, `p2N.m`, developed to analytically compute the polynomial expression for the sum of the first N integers raised to a specific power p.

## Background and Motivation
This script was developed as a mathematical and computational exercise. Instead of relying on built-in polynomial fitting tools that may yield floating-point approximations, the goal was to derive the general form of Faulhaber's formula from first principles and translate it into efficient MATLAB code. 

## Methodology
The function constructs a Vandermonde-like matrix to build a system of linear equations representing the power series. It then solves this system utilizing MATLAB's Symbolic Math Toolbox. This approach ensures the extraction of the exact fractional coefficients for the resulting polynomial, maintaining mathematical rigor.

## Usage
The function `p2N` takes two arguments:
- `p`: The power to which the integers are raised.
- `N_value`: The specific integer N at which you want to evaluate the sum.

**Example:**
To find the polynomial equivalent for the sum of squares (p=2) up to N=10, run the following in your Command Window:
```matlab
[polynomial_formula, total_sum] = p2N(2, 10)
```
The function will output the exact symbolic polynomial expression as a function of N, along with the evaluated numeric sum.

## Requirements 
- MATLAB
- Symbolic Math Toolbox
