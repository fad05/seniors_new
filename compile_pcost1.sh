#!/bin/bash
gfortran -c  parameters.f90 procedures.f90 random.f90 csv_file.f90 nrtype.f90 nrutil.f90
gfortran -g  parameters.f90 procedures.f90 nrtype.f90 nrutil.f90 csv_file.f90 amoeba.f90 truncated_normal.f90 random.f90 labchoice_partcost1.f90 main_partcost1.f95 -o seniors_partcost1

