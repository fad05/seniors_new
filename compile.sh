#!/bin/bash
gfortran -c parameters.f90 nrtype.f90 nrutil.f90  && gfortran -g import_data.f90 import_matrices.f90 import_vector.f90 labchoice.f90 linspace.f90 minpack.f90 main.f95 -o seniors
