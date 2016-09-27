#!/bin/bash
gfortran -c parameters.f90 procedures.f90 
gfortran -g parameters.f90 procedures.f90 labchoice.f90 minpack.f90 main.f95 -o seniors
