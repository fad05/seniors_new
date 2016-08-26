module parameters
    integer, parameter :: grid_asset = 20, nwagebins = 3, nmedbins = 5, ntypes = 8
    integer, parameter :: lifespan = 41, ltot = 4800, kappa = 900, d = 100000
    integer, parameter :: statesize = nwagebins*nmedbins*2
    real (kind = 8), parameter ::  beta = 0.95d0, delta = 0.2d0, eta = 1.8d0, r = 2.0d-2, sigma = 3.75d0, ugamma = 0.5d0
    real (kind = 8), parameter :: tol = 1.0d-5, asset_min = 0.0d0, asset_max = 1.0d5
    save
end module parameters
