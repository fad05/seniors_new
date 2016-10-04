module parameters
    integer, parameter :: grid_asset = 10, grid_ss = 10 ,nwagebins = 3, nmedbins = 5, ntypes = 8, nhealth = 2
    integer, parameter :: lifespan = 41, ltot = 4800, kappa = 900, d = 10000
    integer, parameter :: statesize = nwagebins*nmedbins*nhealth
    real (kind = 8), parameter ::  beta = 0.95d0, delta = 0.2d0, eta = 1.8d0, r = 2.0d-2, sigma = 3.75d0, ugamma = 0.5d0
    real (kind = 8), parameter :: tol = 1.0d-6, asset_min = 0.0d0, asset_max = 5.0d5, cmin = 2600
    real (kind = 8), parameter :: ss_min = 2.4d3, ss_max = 7.2d4 !minimal social security: 200USD/MONTH = 2400USD/YEAR, maximal = 6000USD/MONTH=72000USD/YEAR
    real (kind = 8), parameter :: scale_factor = 1.0d3
    real (kind = 8), parameter :: nra(2) = (/65.0d0,67.0d0/)
    real (kind = 8), parameter :: credit_delret(2) = (/4.25d-2,8.0d-2/)
    real (kind = 8), parameter :: penalty_er3 = 6.6666d-2, penalty_long = 5.0d-2
    !AIME bendpoints for two cohorts; column denote cohort;
    !first entry in column (1,x): first bendpoint, second entry in column (2,x) - second bendpoint
    real (kind = 8), dimension(2,2), parameter :: aime_bend = reshape((/3.72d3,2.2392d4,7.272d3,4.3836d4/),(/2,2/))  
     save
end module parameters
