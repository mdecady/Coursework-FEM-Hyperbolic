# Finite Element Method for Hyperbolic Equation

Course project on numerical methods.

## Problem

Numerical solution of a second-order hyperbolic equation using

- Finite Element Method in space
- Three-layer implicit time scheme
- Uniform and non-uniform time grids
- MSG iterative solver with LU preconditioning

## Equation

χ u_tt + σ u_t - div(λ grad u) + γu = f

## Features

- Quadratic triangular finite elements
- 7-point Gaussian quadrature
- Three-layer implicit time discretization
- Uniform and non-uniform time grids
- Convergence study
- Approximation order study

## Project structure

src/
mesh/
report/
results/

## Build

Visual Studio 2022

or

g++ *.cpp -O2
