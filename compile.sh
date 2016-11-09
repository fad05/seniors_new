#!/bin/bash
gfortran -c parameters.f90 procedures.f90 random.f90 csv_file.f90
gfortran -g parameters.f90 procedures.f90 csv_file.f90 truncated_normal.f90 random.f90 labchoice.f90 main.f95 -o seniors
