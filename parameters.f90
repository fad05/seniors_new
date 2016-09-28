module parameters
    integer, parameter :: grid_asset = 10, grid_ss = 10 ,nwagebins = 3, nmedbins = 5, ntypes = 8, nhealth = 2
    integer, parameter :: lifespan = 41, ltot = 4800, kappa = 900, d = 10000
    integer, parameter :: statesize = nwagebins*nmedbins*nhealth
    real (kind = 8), parameter ::  beta = 0.95d0, delta = 0.2d0, eta = 1.8d0, r = 2.0d-2, sigma = 3.75d0, ugamma = 0.5d0
    real (kind = 8), parameter :: tol = 1.0d-6, asset_min = 0.0d0, asset_max = 5.0d5, cmin = 2600
    real (kind = 8), parameter :: ss_min = 2.4d3, ss_max = 7.2d4 !minimal social security: 200USD/MONTH = 2400USD/YEAR, maximal = 6000USD/MONTH=72000USD/YEAR
    real (kind = 8), parameter :: scale_factor = 1.0d3
     save
end module parameters
