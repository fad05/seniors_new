#!/bin/bash
gfortran -c  parameters.f90 procedures.f90 random.f90 csv_file.f90 nrtype.f90 nrutil.f90
gfortran -g parameters.f90 procedures.f90 nrtype.f90 nrutil.f90 csv_file.f90 amoeba.f90 truncated_normal.f90 random.f90 labchoice_cd.f90 main_cd.f95 -o seniors_cd
