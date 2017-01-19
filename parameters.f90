module parameters
!	1. Parameters for calibration
	integer :: ltot = 3800, kappa = 900, d = 100000, hmin = 0	

    integer, parameter :: grid_asset = 11, grid_ss = 11 ,nwagebins = 5, nmedbins = 5, ntypes = 8, nhealth = 2
    integer, parameter :: lifespan = 41
    integer, parameter :: statesize = nwagebins*nmedbins*nhealth
    integer, parameter :: nsim = 5000 !number of simulations
    integer, parameter :: tmom(2) = (/34, 18/) !number of periods available for method of simulated moments: 34 periods (57 to 90) for c-0, 18 periods (50 to 67) for c-1
    integer, parameter :: init_workyears = 25	
    integer testparam	
    									
    real (kind = 8), parameter ::	beta = 0.95d0,	delta = 0.2d0, eta = 1.8d0, r = 2.0d-2, sigma = 3.75d0, ugamma = 0.5d0
    real (kind = 8), parameter ::	xi = 1.4d0, nu = 2.0d0, age0 = 60
    real (kind = 8), parameter :: 	tol = 1.0d-6, asset_min = 1.0d3, asset_max = 5.0d5, cmin = 2600
    
    !inflation adjustment : http://data.bls.gov/cgi-bin/cpicalc.pl
    !I calculate cap in the following way: calculate PIA at NRA taking base income from the table and increase PIA actuarially until age 70
    !Year 2003: max income 87000 (interpret it as AIME); bendpoints are  7272$ and 43836$ => PIA = 7272*0.9+0.32*(43836-7272)+0.15*(87000-43836) = 24720$ (2060$ monthly); NRA is 66.
    !at age 70 it's going to be B70 = 24720*(1+4*0.08) = 32630 ~ 32700 (this is our ss_max)
    real (kind = 8), parameter :: ss_min = 1.0d3, ss_max = 3.467d4 !minimal social security: 1000USD/YEAR (arbitrary), maximal (real) ~ 2060USD/MONTH~32700USD/YEAR(*CPI2005/2003 = 1.06) = 34670$  to adjust for inflation
    real (kind = 8), parameter :: aime_max = 9.0d4	!Contribution And Benefit Base for 2005 is 90000USD (namely, a base for benefit calc cant be larger) 
    real (kind = 8), parameter :: scale_factor = 1.0d3
    real (kind = 8), parameter :: nra(2) = (/65.0d0,66.0d0/)
    real (kind = 8), parameter :: credit_delret(2) = (/4.25d-2,8.0d-2/) !averages for the first cohort (3% to 5%, and 8% to second cohort)
    real (kind = 8), parameter :: penalty_er3 = 6.6666d-2, penalty_long = 5.0d-2
    
    
    !AIME bendpoints for two cohorts; column denote cohort;
    !first entry in column (1,x): first bendpoint, second entry in column (2,x) - second bendpoint
    
    real (kind = 8), dimension(2,2), parameter :: aime_bend = reshape((/6.395d3,3.8496d4,7.719d3,4.6528d4/),(/2,2/))  !bendpoints are for year 1987 and 2003; adjusted by CPI to year 2005 http://data.bls.gov/cgi-bin/cpicalc.pl
    !benefit caps also exist, according to SSA Contribution And Benefit Base: https://www.ssa.gov/oact/cola/cbb.html
    !I take again years 1987 ans 2003
    !I calculate cap in the following way: calculate PIA at NRA taking base income from the table
    !Example 2003: max income 87000 (interpret it as AIME); bendpoints are  7272$ and 43836$ => PIA = 7272*0.9+0.32*(43836-7272)+0.15*(87000-43836) = 24720$ (2060$ monthly); NRA is 66. Using cpi, in 2005$ is  
    !Example 1987: max income 43800; bendpoints are 3720 and 22392 => PIA at NRA (65) = 3720*0.9+0.32*(22392-3720)+0.15*(43800-22392) = 12534$; using cpi, in 2005$ it is 21548$
    real (kind = 8), dimension(2), parameter :: sscap_nra = (/2.1548d4,2.6238d4/) !Adjusted to 2005USD using cpi1987 and cpi2003 http://data.bls.gov/cgi-bin/cpicalc.pl
    
    !Taxation and ss taxation
    !tax formula is of the form: tax(y) = b*(1-(s*y**p+1)**(-1/p)); numbers are taken from Guner etal, 2014(RED)
    real (kind = 8), parameter :: btax = 2.64d-1, stax = 1.36d-2, ptax = 9.64d-1 !stax is adjusted: stax  = s(from paper)/f^(-p); f = CPI2005/CPI2000 = 1.13416; s = 0.012
    !taxation of benefits
    real (kind = 8), parameter :: first_treshold(2) = (/2.5d4,3.2d4/), second_treshold(2) = (/3.4d4,4.4d4/)
    
    
    !Earning test exempt amounts
    real (kind = 8), parameter :: et_2005 = 12000 !exempt amount of young cohort: 62-66
    real (kind = 8), parameter :: et_1994(2) = (/10593,14704/) !two exempt amounts for old cohort: 62-65 and 66-70 years of age
     save
end module parameters
