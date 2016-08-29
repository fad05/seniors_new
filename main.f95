! Updated Fortran code for "seniors" project
! By Alexey Filatov
program main

use parameters
implicit none


interface
subroutine import_data(datafile,outdata)
    implicit none
    character*200, intent(in) :: datafile ! address of a data file to read
    real (kind = 8), dimension(:,:), intent(out) :: outdata
end subroutine
subroutine import_matrices(datafile,outdata)
    implicit none
    character*200, intent(in) :: datafile ! address of a data fie to read
    real (kind = 8), dimension(:,:,:), intent(out) :: outdata
end subroutine import_matrices
subroutine import_vector(datafile,outdata)
    implicit none
    character*200, intent(in) :: datafile ! address of a data file to read
    real (kind = 8), dimension(:), intent(out) :: outdata
end subroutine import_vector
subroutine  linspace(a,b,n,outvec)
    real (kind = 8), intent(in) :: a, b !starting point, endpoint
    integer, intent(in) :: n !number of points
    real (kind = 8), dimension(n), intent(out) :: outvec !output vector
end subroutine linspace
subroutine valfun(wage,medexp,transmat,surv,ms,vf,pf,lchoice)
    use parameters
    implicit none
    real (kind = 8), dimension(:,:), intent(in) :: wage,surv
    real (kind = 8), dimension(:,:,:), intent(in) :: medexp, transmat
    integer, dimension(statesize,3), intent(in) :: ms
    real (kind = 8), dimension(:,:,:), intent(out) :: vf
    integer, dimension(:,:,:), intent(out) :: pf
    real (kind = 8), dimension(:,:,:), intent(out) :: lchoice
end subroutine valfun
subroutine labchoice (lchoice, acur,anext, h, wage, mexp)
    use parameters
    implicit none
    real (kind = 8), intent(in) :: acur, anext, wage,mexp !current assets, next period assets, current wage, mexp and current wage and mexp shocks
    integer, intent(in) :: h !health status
    real (kind = 8), intent(out) :: lchoice !optimal labor choice
end subroutine labchoice
end interface

!PROGRAM VARIABLES
character*200 datafile !address of the data file on the disk
character*30, names(5)
!wageprof is a matrix that holds coefficients of the wage profile for each type;
!dimensons are:
!1: 8 types (2 cohorts,2 sexes,2 college),
!2: cohort sex college, profile const, profile coefficient
!3: health (bad(1) or good(2)) (for med coefs)
!Common factors for wage and med profiles are cohort and college.
!Wage is further conditioned by sex, and med is further conditioned by health
real (kind = 8) wagecoefs(8,5), wages(lifespan,nwagebins,ntypes), wage_z(nwagebins,1),wagetrans(nwagebins,nwagebins,lifespan-1)
real (kind = 8) medcoefs(8,5,2), medexp(lifespan,2,nmedbins,ntypes), mexp_z(nmedbins,1), medtrans(nmedbins,nmedbins,lifespan-1)
real (kind = 8), dimension(8,3) :: typemat
real (kind = 8) healthtrans(lifespan-1,2,ntypes), surv(lifespan,2,ntypes)
real (kind = 8) htrans_full(2,2,lifespan-1,ntypes)
real (kind = 8) transmat(statesize,statesize,lifespan-1,ntypes)
real (kind = 8) vf(grid_asset,lifespan+1,statesize) !lifespan + additional period of bequests (leave bequest at age 91)
integer pf(grid_asset,lifespan,statesize) !lifespan (no policy for age 91)
real (kind = 8) lchoice(grid_asset,lifespan,statesize) !optimal labor choice given current asset, current state and current age
integer i,j,k,t,w,wn,m,mn,h,hn,ind_type
integer :: cohort, sex, college, health, age(lifespan) = (/(i,i=50,90)/)
integer ms(statesize,3)



!PROGRAM BODY

!-----------------------------------------------------------------------
!1. Import files
!1.1. Import coefficients of wage profiles
datafile = 'data_inputs/wage_profiles_coefs.csv'
open(unit = 1, file = datafile, action = "read")
read (1,*) names(:)
!print *, names
do i = 1,ntypes
    read (1,*) wagecoefs(i,:)
    print *, wagecoefs(i,:)
enddo
close(1)

!1.2. Import coefficients of medical expenditure profiles
datafile = 'data_inputs/medical_expenses_profiles_coefs.csv'
open(unit = 1, file = datafile, action = "read")
read (1,*) names(:)
!print *, names
do i = 1,ntypes/2
    read (1,*) medcoefs(2*i-1,:,1) !bad health, sex 0
    medcoefs(2*i,:,1) = medcoefs(2*i-1,:,1) !coefficients at bad health, sex 1 are equal to bad health, sex 0
    medcoefs(2*i,3,1) = 1 !setting sex to 1
    read (1,*) medcoefs(2*i-1,:,2) !good health, sex 0
    medcoefs(2*i,:,2) = medcoefs(2*i-1,:,2) !good health, sex 1 is equal to good health, sex 0
    medcoefs(2*i,3,2) = 1 !setting sex to 1
    print *, medcoefs(2*i-1,:,1)
    print *, medcoefs(2*i,:,1)
    print *, medcoefs(2*i-1,:,2)
    print *, medcoefs(2*i,:,2)
enddo
close(1)

!1.3 Import bins for wage and mexp

datafile = 'data_inputs/wage_means_by_bins.csv'
call import_data(datafile,wage_z)
datafile = 'data_inputs/mexp_means_by_bins.csv'
call import_data(datafile,mexp_z)

!1.4 Import transition matrices
!1.4.1 Wage transitions
datafile = 'data_inputs/transmat_wage.csv'
call import_matrices(datafile,wagetrans)
!1.4.2 Med expenses transitions
datafile = 'data_inputs/transmat_mexp.csv'
call import_matrices(datafile,medtrans)
!1.4.3 Health transitions, survival and type matrix
do cohort = 0,1
    do sex = 0,1
        do college = 0,1
            do health = 0,1
                !type matrix
                !!!! VERY IMPORTANT!!!
                !this has to be consistent with stata formulas!!!
                ind_type = 1+sex+2*college+4*cohort
                typemat(ind_type,:) = (/cohort,college,sex/)

                !import health transitions: for each type and health status, has length of lifespan-1 (at age 90 no transition to any health status)
                WRITE(datafile,'(a,i1,a,i1,a,i1,a,i1,a)') &
                    "data_inputs/healthtrans_coh",cohort,"_sex",sex,"_college",college,"_health",health,".csv"
                call import_vector(datafile,healthtrans(:,health+1,ind_type)) !"short" health transition matrix: only probabilities of transition to a good state
                !print *, size(healthtrans(:,health+1,ind_type))
                !print *, healthtrans(1,:,ind_type)
                !import survival rates: for each type w/o health and health status, has length of lifespan, but survival rate at age 90 equals to zero
                WRITE(datafile,'(a,i1,a,i1,a,i1,a,i1,a)') &
                    "data_inputs/survival_coh",cohort,"_sex",sex,"_college",college,"_health",health,".csv"
                call import_vector(datafile,surv(1:lifespan-1,health+1,ind_type))
                surv(lifespan,health+1,ind_type) = 0
                !print *, surv(1,:,ind_type)
            enddo
        enddo
    enddo
enddo

!-----------------------------------------------------------------------
!2.!construct wages and medical expenses life cycle profiles for each type
do i = 1,ntypes
    do k = 1,nwagebins
        !print *, wagecoefs(i,4)
        !print *, wagecoefs(i,5)
        !print *, wage_z(k,1)
        !print *, wagecoefs(i,4) + wagecoefs(i,5)*50 + wage_z(k,1)
        !print *, exp(wagecoefs(i,4) + wagecoefs(i,5)*50 + wage_z(k,1))
        wages(:,k,i) = exp(wagecoefs(i,4) + wagecoefs(i,5)*age + wage_z(k,1))
        wages(26:lifespan,k,i) = 0 ! once agent reaches 76, no wage available
        !print *, wages(:,k,i)
    enddo

    do k = 1,nmedbins
        do j = 1,2 !over health statuses
            medexp(:,j,k,i) = exp(medcoefs(i,4,j) + medcoefs(i,5,j)*age + mexp_z(k,1))
            !print *, medexp(:,j,k,i)
        enddo
    enddo
    !print *, medexp(:,3,i)
enddo

!-----------------------------------------------------------------------
!3. Construct matrix of the correspondence of exostate to shocks
do w = 1,nwagebins
    do m = 1,nmedbins
        do h = 1,2
            k = h+2*(m-1)+nmedbins*2*(w-1)
            ms(k,:) = (/w,m,h/) !matrix that translates combination of wage shocks, med exp shocks and health shocks into a single number
            print *, ms(k,:)
        enddo
    enddo
enddo

!-----------------------------------------------------------------------
!4. construct Kronecker product of transition matrices for each individual type and each age;
!result is stored in transmat
do k = 1,ntypes
    do t = 1,lifespan-1
        !construct full health trans matrix: the difference from reduced form is
        !in probabilities of transition from both good and bad states,
        !it is just expansion of reduced matrix
        htrans_full(1,1,t,k) = 1-healthtrans(t,1,k) !bad-to-bad health
        htrans_full(1,2,t,k) = healthtrans(t,1,k) !bad-to-good health
        htrans_full(2,1,t,k) = 1-healthtrans(t,2,k) !good-to-bad health
        htrans_full(2,2,t,k) = healthtrans(t,2,k) !good-to-good health
        !print *, wagetrans(:,:,t)
        !print *, medtrans(:,:,t)
        !print *, htrans_full(:,:,t,z)
        do j =  1,statesize !rows of transition matrix
            do i = 1,statesize !columns of transition matrix
                w = ms(i,1) !current wage  shock: 1 to grid_wage
                wn = ms(j,1) !next period wage
                m = ms(i,2) !current medical shock: 1 to grid_med
                mn = ms(j,2) !next period med
                h = ms(i,3) !current health status: 1(unhealthy) or 2(healthy)
                hn = ms(j,3) !next period med shock
                transmat(i,j,t,k) = wagetrans(w,wn,t)*medtrans(m,mn,t)*htrans_full(h,hn,t,k)
            enddo
        enddo
    enddo
enddo
!-----------------------------------------------------------------------
!5. Pick a type between 1 and 8, and calculate optimal labor choice matrix and value function for a chosen type
ind_type = 1
!labor choice
lchoice = -10
!value function
call valfun(wages(:,:,ind_type),medexp(:,:,:,ind_type),transmat(:,:,:,ind_type),surv(:,:,ind_type),ms,vf,pf,lchoice)
end
!end main
!========================================================================
subroutine valfun(wage,medexp,transmat,surv,ms,vf,pf,lchoice)
    use parameters
    implicit none
!INPUTS AND OUTPUTS
    !dimensions, wage: lifespan x nwagebins
    !dimensions, surv: lifespan x 2 (health)
    real (kind = 8), dimension(:,:), intent(in) :: wage, surv
    !dimensions: lifespan x 2 (health) x nmedbins
    real (kind = 8), dimension(:,:,:), intent(in) :: medexp, transmat
    !dimensions: statesize x 3
    integer, dimension(statesize,3), intent(in) :: ms
    !dimensions: asset x lifespan+1 x state
    real (kind = 8), dimension(:,:,:), intent(out) :: vf
    !dimensions: asset x lifespan x state
    integer, dimension(:,:,:), intent(out) :: pf
    !dimensions: asset x lifespan x state
    real (kind = 8), dimension(:,:,:), intent(out) :: lchoice
!INTERNAL VARIABLES
    real (kind = 8) assets(grid_asset), beq(grid_asset), cons, util, opt_labor, val_nomax(grid_asset,grid_asset)
    real (kind = 8) lab_cur_next(grid_asset,grid_asset) !labor choice given current and next period asset
    integer i,j,k,t,s,w,m,h, age
!BODY OF VALFUN
!1. Check input and output matrices' dimensions
    if (size(wage,1) .NE. lifespan) then
        print *, 'Wage matrix "lifespan" size is incorrect!'
        vf = -1
        return
    elseif (size(wage,2) .NE. nwagebins) then
        print *, 'Wage matrix "wage bins" size is incorrect!'
        vf = -2
        return
    endif

    if (size(medexp,1) .NE. lifespan) then
        print *, 'Med matrix "lifespan" size is incorrect!'
        vf = -3
        return
    elseif (size(medexp,2) .NE. 2) then
        print *, 'Med matrix "health state" size is incorrect!'
        vf = -4
        return
    elseif (size(medexp,3) .NE. nmedbins) then
        print *, 'Med matrix "med bins" size is incorrect!'
        vf = -5
        return
    endif

    if (size(transmat,1) .NE. statesize) then
        print *, 'Transition matrix "current state" size is incorrect!'
        vf = -6
        return
    elseif (size(transmat,2) .NE. statesize) then
        print *, 'Transition matrix "next-period state" size is incorrect!'
        vf = -7
        return
    elseif (size(transmat,3) .NE. lifespan-1) then
        print *, 'Transition matrix "lifespan-1" size is incorrect!'
        vf = -8
        return
    endif

    if (size(vf,1) .NE. grid_asset) then
        print *, 'Value function "asset" size is incorrect!'
        vf = -9
        return
    elseif (size(vf,2) .NE. lifespan+1) then
        print *, 'Value function "period of life" size is incorrect!'
        vf = -10
        return
    elseif (size(vf,3) .NE. statesize) then
        print *, 'Value function "state" size is incorrect!'
        vf = -11
        return
    endif

    if (size(surv,1) .NE. lifespan) then
        print *, 'Survival matrix "lifespan" size is incorrect!'
        vf = -12
        return
    elseif (size(surv,2) .NE. 2) then
        print *, 'Survival matrix "health state" size is incorrect!'
        vf = -13
        return
    endif

    if (size(pf,1) .NE. grid_asset) then
        print *, 'Policy function "asset" size is incorrect!'
        vf = -14
        return
    elseif (size(pf,2) .NE. lifespan) then
        print *, 'Policy function "period of life" size is incorrect!'
        vf = -15
        return
    elseif (size(pf,3) .NE. statesize) then
        print *, 'Policy function "state" size is incorrect!'
        vf = -16
        return
    endif

    if (size(lchoice,1) .NE. grid_asset) then
        print *, 'Labchoice "asset" size is incorrect!'
        vf = -17
        return
    elseif (size(lchoice,2) .NE. lifespan) then
        print *, 'Labchoice "period of life" size is incorrect!'
        vf = -18
        return
    elseif (size(lchoice,3) .NE. statesize) then
        print *, 'Labchoice "state" size is incorrect!'
        vf = -19
        return
    endif

!2. IF INPUTS ARE CORRECT, CALCULATE VALUE FUNCTION BY BACKWARDS INDUCTION
    !create assets grid
    call linspace(asset_min,asset_max,grid_asset,assets)
    !bequests depending on assets
    beq = eta*(assets+d)**(1-sigma)/(1-sigma)
    !Value function next period has two dimensions: possible CURRENT assets and CURRENT state of the world:
    !a combination of wage shock, med expenditure shock and health.
    !Every particular combination of current asset and current shock in equilibrium define an optimal behavior.
    !Namely, optimal choice of next period assets, labor supply and consumption. These choices define current utility, and,
    !together with transition probabilities from current state, next period expected utility. Combination of these two define value function.

    !First, calculate expected value function upon death (at t = lifespan(= 90)).
    !Value function upon death is just utility of bequest, which uniquely defined by assets that agent chooses to leave as bequest.
    !No uncertainty here, since everybody dies with certainty after age of 90.
    !Similarly, value function upon death is the same as expected value function upon death (since no uncertainty).
    !value function upon death
    do s = 1,statesize
        vf(:,lifespan+1,s) = beq !at age 91, when dead
    enddo
    !MAIN ITERATION CYCLE
    do t = lifespan-16,1,-1 !for each age 90 to 50
		age = t+49
        do s = 1,statesize !for every possible state
            !here, determine wage, med and health shock corresponding to the current state
            w = int(ms(s,1))
            m = int(ms(s,2))
            h = int(ms(s,3)) !h is either 1 or 2
            !over all possible combinations of assets:
            do j = 1,grid_asset     !next-period assets
                do i = 1,grid_asset !current assets 
                    !First, calculate optimal labor choice
					call labchoice(opt_labor, assets(i),assets(j),h,wage(t,w),medexp(t,h,m))
                    !print out if optimal labor <0
                    if (opt_labor<0) then
                        print *, opt_labor
                        print *, 'Labor choice is out of bounds:', opt_labor
                        call sleep(60)
                    elseif (opt_labor>ltot-kappa) then
						print *, opt_labot
						print *, 'Labor choice is out of bounds:', opt_labor
                        call sleep(60)
                    endif
                    lab_cur_next(i,j) = opt_labor
                    !given optimal choice of labor: calculate consumption and utility, different cases
                    if (assets(j) == 0) then !if next-period asset is zero, then govt transfer to cover minimal consumption is possible
                    cons = (1+r)*assets(i) + wage(t,w)*opt_labor-assets(j)-medexp(t,h,m) !uncompensated consumption
						if (cons < cmin) then
							util = ((1+delta*(h-1))*(cons**ugamma*ltot**(1-ugamma))**(1-sigma))/(1-sigma)
						endif
                    endif
                    !unmaxed "value function"
                    
                    val_nomax(i,j) = util + beta*(surv(t,h)*dot_product(vf(j,t+1,:),transmat(s,:,t)) &
                                + (1-surv(t,h))*beq(j))
                enddo
            enddo
            vf(:,t,s) = maxval(val_nomax,2)
            pf(:,t,s) = maxloc(val_nomax,2)
            lchoice(:,t,s) = (/(lab_cur_next(k,pf(k,t,s)),k=1,grid_asset)/)
        enddo
    enddo
    return
end subroutine valfun

