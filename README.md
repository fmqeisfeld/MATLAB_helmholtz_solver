# MATLAB Helmholtz-Solver

## Description
This Matlab script demonstrates an application for a simple 1D-Helmholtz-solver for inhomogeneous media by means of the Transfer-Matrix method.
Note: The code is not utilizing Matlab's vectorization capabilities to increase performance, since it served as prototype for a C-port.
Its main purpose is to provide you an idea of how laser light of a specific wavelength, polarization und angle of incidence 
interacts with some material or several layers of differenct materials.
![example](<example.png>)

## Features
- user defined, arbitrary density profile
- you can adjust the angle of incidence and the polarization of the incoming light (s- or p-polarization)
- plots the profile for the calculated absorbed power density
- outputs the integral absorption, reflection and transmission
- the material's permittivity can be modeled as a function of wavelenght, density and temperature

## Usage

In order to compute the absolute power density in units of W/m^3, simply multiply the relative power density by the incident laser intensity.
The main input-parameters are the following:
- *m_polar*: either 1 (for s-polarization) or 2 (for p-polarization)
- *lambda*: wavelength of the incident radiation in units of meters
- *theta*: angle of incidence
- *nelements*: nr of piecewise constant material elements
- *delta*: the widht of each element in units of nano-meters
- *dprof*-vector: an adjustable vector containing the density profile of the material. The same can be done for the temperature.
- *getEpsilon(lambda,Te,Ti,rho)*: a function, giving the relative permittivity as a function of electron temperature, ion-temperature, density and wavelength. 
