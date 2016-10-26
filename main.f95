! Updated Fortran code for "seniors" project
! By Alexey Filatov
program main

use parameters
use procedures
implicit none


interface
subroutine valfun(cohort,wage,medexp,transmat,surv,ms,vf,vf_na,pf,pf_na,lchoice,lchoice_na,app_policy)
    use parameters
    use procedures
    implicit none
	integer, intent(in) :: cohort
    real (kind = 8), dimension(lifespan,nwagebins), intent(in) :: wage
    real (kind = 8), dimension(lifespan,nhealth), intent(in) :: surv
    real (kind = 8), dimension(lifespan,nhealth,nmedbins), intent(in) :: medexp
    real (kind = 8), dimension(statesize,statesize,lifespan-1), intent(in) :: transmat
    integer, dimension(statesize,3), intent(in) :: ms
    real (kind = 8), dimension(grid_asset,statesize,grid_ss,lifespan+1), intent(out) :: vf
    integer, dimension(grid_asset,statesize,grid_ss,lifespan), intent(out) :: pf
    real (kind = 8), dimension(grid_asset,statesize,grid_ss,lifespan), intent(out) :: lchoice
	real (kind = 8), dimension(grid_asset,statesize,grid_ss,lifespan+1,2), intent(out) :: vf_na
	integer, dimension(grid_asset,statesize,grid_ss,lifespan,2), intent(out) :: pf_na
	real (kind = 8), dimension(grid_asset,statesize,grid_ss,lifespan,2), intent(out) :: app_policy
	real (kind = 8), dimension(grid_asset,statesize,grid_ss,lifespan,2), intent(out) :: lchoice_na
end subroutine valfun
subroutine labchoice (lchoice, cons, util, acur,anext, h, wage, mexp)
    use parameters
    use procedures
    implicit none
    real (kind = 8), intent(in) :: acur, anext, wage,mexp !current assets, next period assets, current wage, mexp
    integer, intent(in) :: h !health status, 1 or 2
    real (kind = 8), intent(out) :: lchoice !optimal labor choice
    real (kind = 8), intent(out) :: cons !consumption under optimal choice
    real (kind = 8), intent(out) :: util !value of utility under optimal choice
end subroutine labchoice
subroutine lc_simulation(cohort,in_asset,in_aime,years_aime,wage,mexp,ms,transmat,vf,vf_na,pf,pf_na,lchoice,lchoice_na,app_policy, &
						life_assets,life_labor,app_age)
	use parameters
	use procedures
	integer, intent(in) :: cohort
	real (kind = 8), intent(in) :: in_asset, in_aime, years_aime
	real (kind = 8), dimension(lifespan,nwagebins), intent(in) :: wage
	real (kind = 8), dimension(lifespan,nhealth,nmedbins), intent(in) :: mexp
	integer, dimension(statesize,3), intent(in) :: ms
	real (kind = 8), dimension(statesize,statesize,lifespan-1), intent(in) :: transmat
	real (kind = 8), intent(in) :: vf(grid_asset,statesize,grid_ss,lifespan+1) 
	integer, intent(in) :: pf(grid_asset,statesize,grid_ss,lifespan)
	real (kind = 8), intent(in) :: lchoice(grid_asset,statesize,grid_ss,lifespan) 
	real (kind = 8), dimension(grid_asset,statesize,grid_ss,lifespan+1,2), intent(in) :: vf_na
	integer, dimension(grid_asset,statesize,grid_ss,lifespan,2), intent(in) :: pf_na
	real (kind = 8), dimension(grid_asset,statesize,grid_ss,lifespan,2), intent(in) :: app_policy
	real (kind = 8), dimension(grid_asset,statesize,grid_ss,lifespan,2), intent(in) :: lchoice_na
	real (kind = 8), intent(out) :: life_assets(lifespan+1), life_labor(lifespan)
	integer, intent(out) :: app_age
end subroutine lc_simulation
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
real (kind = 8) wagecoefs(8,5), wages(lifespan,nwagebins,ntypes), wage_z(nwagebins,1),wagetrans(nwagebins,nwagebins,lifespan-1), &
				mlwages(lifespan,ntypes) !mean log wages for each type

real (kind = 8) medcoefs(8,5,2), medexp(lifespan,2,nmedbins,ntypes), mexp_z(nmedbins,1), medtrans(nmedbins,nmedbins,lifespan-1)
real (kind = 8), dimension(8,3) :: typemat
real (kind = 8) healthtrans(lifespan-1,2,ntypes), surv(lifespan,2,ntypes)
real (kind = 8) htrans_full(2,2,lifespan-1,ntypes)
real (kind = 8) transmat(statesize,statesize,lifespan-1,ntypes)
!for social security applicants
real (kind = 8) vf(grid_asset,statesize,grid_ss,lifespan+1) !lifespan + additional period of bequests (leave bequest at age 91)
integer pf(grid_asset,statesize,grid_ss,lifespan) !lifespan (no policy for age 91)
real (kind = 8) lchoice(grid_asset,statesize,grid_ss,lifespan) !optimal labor choice given current asset, current state and current age
!for non-applicants
real (kind = 8), dimension(grid_asset,statesize,grid_ss,lifespan+1,2) :: vf_na
integer, dimension(grid_asset,statesize,grid_ss,lifespan,2) :: pf_na
real (kind = 8), dimension(grid_asset,statesize,grid_ss,lifespan,2) :: app_policy
real (kind = 8), dimension(grid_asset,statesize,grid_ss,lifespan,2) :: lchoice_na
!results of simulation
real (kind = 8) life_assets(lifespan+1), life_labor(lifespan)
integer app_age
	
	
integer i,j,k,l,t,w,wn,m,mn,h,hn,ind_type
integer :: cohort, sex, college, health, age(lifespan) = (/(i,i=50,90)/)
integer ms(statesize,3)



!PROGRAM BODY

!user inputs the type

print *, "Enter an integer 1:8, which defines an individual type"
print *, "cohort college sex"
print *, "1: 0	0	0"
print *, "2: 0 	0	1"
print *, "3. 0	1	0"
print *, "4. 0	1	1"
print *, "5. 1	0	0"
print *, "6. 1	0	1"
print *, "7. 1	1	0"
print *, "8. 1	1	1"

read (*,*) ind_type

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
				print *, typemat(ind_type,:) 
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
	mlwages(:,i) = wagecoefs(i,4) + wagecoefs(i,5)*age !mean log wage profiles for each type; we're going to need them for simulation
    do k = 1,nwagebins
        !print *, wagecoefs(i,4)
        !print *, wagecoefs(i,5)
        !print *, wage_z(k,1)
        !print *, wagecoefs(i,4) + wagecoefs(i,5)*50 + wage_z(k,1)
        !print *, exp(wagecoefs(i,4) + wagecoefs(i,5)*50 + wage_z(k,1))
        wages(:,k,i) = exp(wagecoefs(i,4) + wagecoefs(i,5)*age + wage_z(k,1))

		!QUESTIONABLE: TRUNCATED WAGES --------------------------------------
        !wages(26:lifespan,k,i) = 0 ! once agent reaches 76, no wage available
		!QUESTIONABLE: TRUNCATED WAGES --------------------------------------

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
!		TESTING GROUNDS

!5. Pick a type between 1 and 8, and calculate optimal labor choice matrix and value function for a chosen type
ind_type = 1
!labor choice
lchoice = -10
!value function
call valfun(1,wages(:,:,ind_type),medexp(:,:,:,ind_type),transmat(:,:,:,ind_type),surv(:,:,ind_type),ms, &
vf,vf_na,pf,pf_na,lchoice,lchoice_na,app_policy)

cohort = 1

call lc_simulation(cohort,1.5d4,0.8d3,2.5d1, wages(:,:,ind_type),medexp(:,:,:,ind_type),ms, &
transmat(:,:,:,ind_type),vf,vf_na,pf,pf_na,lchoice,lchoice_na,app_policy, life_assets,life_labor,app_age)

!		END TESTING
!-----------------------------------------------------------------------
end
!End of main program
!========================================================================
!Local subroutines
subroutine valfun(cohort,wage,medexp,transmat,surv,ms,vf,vf_na,pf,pf_na,lchoice,lchoice_na,app_policy)
    use parameters
    use procedures
    implicit none
!INPUTS AND OUTPUTS
	integer, intent(in) :: cohort
    real (kind = 8), dimension(lifespan,nwagebins), intent(in) :: wage
    real (kind = 8), dimension(lifespan,nhealth), intent(in) :: surv
    real (kind = 8), dimension(lifespan,nhealth,nmedbins), intent(in) :: medexp
    real (kind = 8), dimension(statesize,statesize,lifespan-1), intent(in) :: transmat
    integer, dimension(statesize,3), intent(in) :: ms
    !vf, pf and labchoice for benefit applicants ("first vf")
    real (kind = 8), dimension(grid_asset,statesize,grid_ss,lifespan+1), intent(out) :: vf
    integer, dimension(grid_asset,statesize,grid_ss,lifespan), intent(out) :: pf
    real (kind = 8), dimension(grid_asset,statesize,grid_ss,lifespan), intent(out) :: lchoice
    !vf, pf, labchoice and application policy for non-applicants ("second vf")
    !note the last dimension of the size 2, which signifies the way we update next period benefits (more or less than 35 working years)
	real (kind = 8), dimension(grid_asset,statesize,grid_ss,lifespan+1,2), intent(out) :: vf_na !value function for those not applied to social security benefits
	integer, dimension(grid_asset,statesize,grid_ss,lifespan,2), intent(out) :: pf_na !policy function for those not applied to social security benefits
	real (kind = 8), dimension(grid_asset,statesize,grid_ss,lifespan,2), intent(out) :: app_policy !each entry is 0-1; application policy!
	real (kind = 8), dimension(grid_asset,statesize,grid_ss,lifespan,2), intent(out) :: lchoice_na !labor choice if not applied
!INTERNAL VARIABLES
    real (kind = 8) assets(grid_asset), beq(grid_asset), cons, util, opt_labor, val_nomax(grid_asset,grid_asset)
    real (kind = 8) val_nomax_na(grid_asset,grid_asset)
    real (kind = 8) ss(grid_ss) !ss: spans on ALL levls of ss INCLUDING 0
    real (kind = 8) lab_cur_next(grid_asset,grid_asset), lab_cur_next_na(grid_asset,grid_asset)  !labor choice given current and next period asset
    integer i,j,k,l,t,s,w,m,h, age, st, bcalc
    real (kind = 8) b, bnext, aime, aime_next !variables to contain current value of ss (on the grid) and next period (calculated by formula), as well as AIMEs
    integer ibprev, ibnext !indices of the interval points on ss grid within which lies calculated next-period benefit 
    real (kind = 8) vf_na_interp(statesize) !container for interpolated value function
    real (kind = 8) val_max_na(grid_asset) !technical variable, maxed vf_na BEFORE comparison vith vf
    

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
    elseif (size(vf,2) .NE. statesize) then
        print *, 'Value function "state" size is incorrect!'
        vf = -10
        return
    elseif (size(vf,3) .NE. grid_ss) then
        print *, 'Value function "social security" size is incorrect!'
        vf = -11
        return
    elseif (size(vf,4) .NE. lifespan+1) then
        print *, 'Value function "lifespan+1" size is incorrect!'
        vf = -12
        return
    endif

    if (size(surv,1) .NE. lifespan) then
        print *, 'Survival matrix "lifespan" size is incorrect!'
        vf = -13
        return
    elseif (size(surv,2) .NE. 2) then
        print *, 'Survival matrix "health state" size is incorrect!'
        vf = -14
        return
    endif

    if (size(pf,1) .NE. grid_asset) then
        print *, 'Policy function "asset" size is incorrect!'
        vf = -15
        return
    elseif (size(pf,2) .NE. statesize) then
        print *, 'Policy function "state" size is incorrect!'
        vf = -16
        return
    elseif (size(pf,3) .NE. grid_ss) then
        print *, 'Policy function "social security" size is incorrect!'
        vf = -17
    elseif (size(pf,4) .NE. lifespan) then
        print *, 'Policy function "lifespan" size is incorrect!'
        vf = -18
        return
    endif

    if (size(lchoice,1) .NE. grid_asset) then
        print *, 'Labchoice "asset" size is incorrect!'
        vf = -19
        return
    elseif (size(lchoice,2) .NE. statesize) then
        print *, 'Labchoice "state" size is incorrect!'
        vf = -20
        return
    elseif (size(lchoice,3) .NE. grid_ss) then
        print *, 'Labchoice "social security" size is incorrect!'
        vf = -21
        return
    elseif (size(lchoice,4) .NE. lifespan) then
        print *, 'Labchoice "lifespan" size is incorrect!'
        vf = -22
        return
    endif

!2.    !create assets grid and social security grid
    assets(1)	= 0
    ss(1)		= 0
    call logspace(asset_min,asset_max,grid_asset-1,assets(2:grid_asset))
    call logspace(ss_min,ss_max,grid_ss-1,ss(2:grid_ss))


!4. IF INPUTS ARE CORRECT, CALCULATE TWO VALUE FUNCTIONS BY BACKWARDS INDUCTION
!	1)	THE FIRST VALUE FUNCTION (SPANS FROM AGE 62 TO AGE 90) IS CALCULATED FOR THOSE WHO HAVE APPLIED FOR SOCIAL
!		SECURITY AND RECEIVE FIXED AMOUNT EVERY YEAR. WE HAVE TO CALCULATE IT FIRST, AS IT IS NEEDED TO CALCULATE THE SECOND ONE.
!	2)	THE SECOND VALUE FUNCTION IS FOR THSE WHO ARE STILL NOT APPLIED FOR SS BENEFITS. IT HAS TO CONTAIN ADDITIONAL DIMENSION:
!		FOR EVERY COMBINATION OF CURRENT SHOCK (WAGE, HEALTH, MEDICINE), CURRENT ASSET LEVEL, AND CURRENT "POTENTIAL" BENEFIT LEVEL
!		(NAMELY, THE AMOUNT OF BENEFITS ONE WILL RECEIVE IF APPLIES TODAY) IT CONTAINS A 0-1 DECISION WHETHER TO APPLY FOR SS OR NOT.
!		THIS ACTS AS AN INDICATOR OF MOVING TO ANOTHER REGIME ("REGIME WITH BENEFITS"). tHIS VALUE FUNCTION SPANS FROM AGE 50 TO AGE 90,
!		BUT ONE CAN APPLY NOT EARLIER THAN AGE 62.

!	bequests depending on assets, scaled properly
	call bequest_value(assets,beq)
    
!	1) FIRST VALUE FUNCTION
    !Value function next period has two dimensions: possible CURRENT assets and CURRENT state of the world:
    !a combination of wage shock, med expenditure shock and health, and current social security level.
    !For the second value function, I'm going to interpolate between levels to estimate social security at points off the grid.
    
    !Every particular combination of current asset, current shock and ss level shock in equilibrium define an optimal behavior.
    !Namely, optimal choice of next period assets, labor supply and consumption. These choices define current utility, and,
    !together with transition probabilities from current state, next period expected utility. Combination of these two define value function.

    !First, calculate expected value function upon death (at t = lifespan(= 90)).
    !Value function upon death is just utility of bequest, which uniquely defined by assets that agent chooses to leave as bequest.
    !No uncertainty here, since everybody dies with certainty after age of 90.
    !Similarly, value function upon death is the same as expected value function upon death (since no uncertainty).
    !value function upon death
    
    !initialize vf and pf; if we will observe these numbers after the execution, means something wrong
    
    vf = -101
    pf = -102
    
    do k = 1,grid_ss
		do s = 1,statesize
			vf(:,s,k,lifespan+1) = beq !at age 91, when dead; for every state and every level of social security is the same
		enddo
    enddo
    
    !MAIN ITERATION CYCLE
    vf(:,:,:,1:12) = -1.0d5 !no values until age 62
    do t = lifespan,13,-1 !for each age 90 to 62
		age = t+49
		do k = 1,grid_ss !for each grid level of social security
	        do s = 1,statesize !for every possible state
	            !here, determine wage, med and health shock corresponding to the current state
	            w = int(ms(s,1))
	            m = int(ms(s,2))
	            h = int(ms(s,3)) !h is either 1 or 2
	            !over all possible combinations of assets:
	            do j = 1,grid_asset     !next-period assets
	                do i = 1,grid_asset !current assets
	                    !First, calculate optimal labor choice
!	                    if (age == 65 .AND. k == 4 .AND. i == 5 .AND. w == 2) then
!							print *, 'stop'
!						endif
						call labchoice(opt_labor, cons, util, assets(i),assets(j),h,wage(t,w),medexp(t,h,m),ss(k))
	                    lab_cur_next(i,j) = opt_labor
	                    !print *, opt_labor

	                    !unmaxed "value function"
	                    val_nomax(i,j) = util + beta*(surv(t,h)*dot_product(vf(j,:,k,t+1),transmat(s,:,t)) &
	                                + (1-surv(t,h))*beq(j))
	                enddo
	            enddo
	            vf(:,s,k,t) = maxval(val_nomax,2)
	            pf(:,s,k,t) = maxloc(val_nomax,2)
!	            if (age <= 60) then
!		            do l = 1,grid_asset
!						print *, val_nomax(l,1), val_nomax(l,2), val_nomax(l,3), val_nomax(l,4), val_nomax(l,5), val_nomax(l,6), &
!								val_nomax(l,7), val_nomax(l,8),  val_nomax(l,9),val_nomax(l,10)
!		            enddo
!		            do l = 1,grid_asset
!						print *, vf(l,s,k,t)
!					enddo
!		            do l = 1,grid_asset
!						print *, pf(l,s,k,t)
!					enddo
!				endif
	            lchoice(:,s,k,t) = (/(lab_cur_next(l,pf(l,s,k,t)),l=1,grid_asset)/)
	        enddo
        enddo
    enddo
!	FIRST VALUE FUNCTION FINISHED

!	2)NOW GIVEN THE FIRST VALUE FUNCTION, WE CAN CALCULATE THE SECOND ONE ALONG WITH THE DECISION RULE:
!	WHETHER TO APPLY TO SOCIAL SECURITY, A CRUCIAL OBJECT WE'RE INTERESTED IN.
!	THE SECOND VALUE FUNCTION HAS MORE ELABORATE STRUCTURE. WE NEED A GRID ON SOCIAL SECURITY HERE TO CALCULATE 
!	CURRENT DECISION BY COMPARISON. THE DIFFERENCE FROM THE FIRST VALUE FUNCTON IS: IN THE 1ST, WE HAD LEVEL OF SOCIAL SECURITY THAT 
!	WERE FIXED (IF YOU HAVE 1000USD THIS YEAR, YOU'LL HAVE 1000USD NEXT YEAR). IN THE 2D ONE, WE HAVE TO HAVE ALL LEVELS AS WELL.
!	INTERPRETATION IS THIS: EACH LEVEL SHOWS A "POSSIBLE" BENEFIT THAT INDIVIDUAL WILL RECEIVE IF HE APPLIES RIGHT NOW. IT IS POSSIBLE THAT 
!	UNDER CERTAIN CONDITIONS INDIVIDUALS HAVE 0 CURRENT POSSIBLE BENEFITS (SAY, THEY DIDN'T WORK ENOUGH WORKING YEARS TO QUALIFY FOR THEM).
!	SO WHAT WE DO IS FOR EVERY SS GRID POINT, EVERY CURRENT LEVEL OF ASSET, EVERY STATE AND EVERY YEAR WE ESTIMATE WEIGTHED AVERAGE OF
!	NEXT-YEAR POSSIBILITIES.
!	EXAMPLE:   assume current shock is S, current asset is A, current "possible" level of benefits B, current age is AGE. Then, the utility in the case
!	i choose to apply for the benefit is just vf(A,S,B,AGE). But which is the utility if i choose not to apply? Well, current utility is sum of contemporaneous 
!	utility and continuation value. Continuation value is just discounted value function next period. 
!	Note that in the last period T everybody applies to ss no mater what. Then, the two value functions are identical in the last period for every possible B (social security). 
!	In the period T-1, however, there are two possibilities: apply now with current B and switch to the first value function with this B, or don't apply
!	and next period receive B' = F(B) (B'  is a function of B which is a function of AIME; AIME' =  AIME + max{0,(wl-AIME)/35} if a person worked more than 35 years, and
!	AIME' = (AIME+wl)/35 if individual worked less; for now i assume that at age 62 everybody worked at least 35 years). We calculate both these scenarios (the "application" 
!	scenario is already calculated in the first vf), and choose the one that grants higher utility (setting ss policy to 0 or 1 accordingly), and thus second VF ths period is calculated.  
	
	vf_na(:,:,:,lifespan:lifespan+1,1) = vf(:,:,:,lifespan:lifespan+1) !upon death and in the last period, the two value functions are identical (everyone applies at age 90 with certainty)
	vf_na(:,:,:,lifespan:lifespan+1,2) = vf(:,:,:,lifespan:lifespan+1)
	app_policy = 0 !initialize application policy
	app_policy(:,:,:,lifespan,:) = 1 ! We claim that in the last period of life everybody applies for ss with certainty
	
	do bcalc = 1,2 
		!bcalc = 1: we calculate "next period" benefits as if person has less than 35 working years
		!bcalc = 2: calculate bnext as if an individual has more than 35 working years
		do t = lifespan-1,1,-1 !age 89 to 50
			age = t+49 !real age
			do k = 1,grid_ss
				do s = 1,statesize !for every possible state
					!here, determine wage, med and health shock corresponding to the current state
		            w = int(ms(s,1))
		            m = int(ms(s,2))
		            h = int(ms(s,3)) !h is either 1 or 2
					!over all possible combinations of assets:
		            do j = 1,grid_asset     !next-period assets
						do i = 1,grid_asset !current assets
							!First, calculate optimal labor choice WITH 0 social security: "what if I don't apply?"
							call labchoice(opt_labor, cons, util, assets(i),assets(j),h,wage(t,w),medexp(t,h,m),0.0d0)
		                    lab_cur_next_na(i,j) = opt_labor
		                    
		                    !NOW, HERE ARE THE BIG DIFFERENCES BETWEEN THIS VF AND THE FIRST ONE
		                    !First, we calculate next-period B (it's going to be off the grid), and use linear interpolation to calculate
		                    !value functions between points of B. Then we have to take expectation of these interpolations.
		                    !On top of that, we need to control for application age: an individual can't apply if she's younger than 62
		                    
		                    !current benefit, in monetary terms
		                    b = ss(k) 
		                    !current AIME, in monetary terms		                    
		                    call aime_calc(cohort,age,b,aime)
		                    
		                    !update AIME for next period
		                    if (bcalc == 1) then !individual worked less than 35 years
								aime_next = aime + wage(t,w)*opt_labor/35.0 		      !next period aime grows              
		                    else !individual aready has more than 35 working years
								aime_next = aime + max(0.0,(wage(t,w)*opt_labor-aime)/35.0) !next-period aime doesn't decline
		                    endif
		                    
		                    !calculate next period benefit using proper age and aime NEXT PERIOD
		                    call benefit_calc(cohort,age+1,aime_next,bnext)
		                    
		                    !bnext lies between ss(ibprev) and ss(ibnext), we need to find these
		                    ibprev = locate_greater(ss,bnext)-1 !index of a previous gridpoint
		                    ibnext = ibprev+1 !index of a next gridpoint
		                    
		                    !Calculate set of interpolations for next-period vf_na: one interpolation for every state s
		                    do st = 1,statesize
								vf_na_interp(st) = linproj(bnext,ss(ibprev),ss(ibnext),vf_na(j,st,ibprev,t+1,bcalc),vf_na(j,st,ibnext,t+1,bcalc)) !vf_na(j - next period asset, st - next-period state,ib.. - indices of ss, t+1 - priod)
		                    enddo
	
		                    val_nomax_na(i,j) = util + beta*(surv(t,h)*dot_product(vf_na_interp,transmat(s,:,t)) &
		                                + (1-surv(t,h))*beq(j))
						enddo
		            enddo
		            !Now, calculate vf_na; notice the difference!
					val_max_na = maxval(val_nomax_na,2)
					pf_na(:,s,k,t,bcalc) = maxloc(val_nomax_na,2)
					lchoice_na(:,s,k,t,bcalc) = (/(lab_cur_next_na(l,pf_na(l,s,k,t,bcalc)),l=1,grid_asset)/)
					
!					if (age == 90 .AND. k == 3 .AND. s == 7) then
!						read (*,*)
!					endif
					!If age is larger or equal than 62, an individual CAN apply for socal security.

					if (age>=62) then
!						print *, vf(:,s,k,t)
!						print *, val_max_na
!						print *, app_policy(:,s,k,t,bcalc)
!						print *, lchoice_na(:,s,k,t,bcalc)
!						print *, lchoice(:,s,k,t)
!						print *, vf_na(:,s,k,t,bcalc)
						where (vf(:,s,k,t) >= val_max_na) !ss application grants higher utility
							vf_na(:,s,k,t,bcalc) = vf(:,s,k,t) 
							app_policy(:,s,k,t,bcalc) = 1	!individual applies for ss, indicator of regime swithing
							lchoice_na(:,s,k,t,bcalc) = lchoice(:,s,k,t)
						elsewhere
							vf_na(:,s,k,t,bcalc) = val_max_na 
							app_policy(:,s,k,t,bcalc) = 0 !individual doesn't apply
						endwhere
		            else !individual doesn't apply
						vf_na(:,s,k,t,bcalc) = val_max_na
						app_policy(:,s,k,t,bcalc) = 0 
		            endif
!					print *, vf(:,s,k,t)
!					print *, val_max_na
!					print *, app_policy(:,s,k,t,bcalc)
!					print *, lchoice_na(:,s,k,t,bcalc)
!					print *, lchoice(:,s,k,t)
!					print *, vf_na(:,s,k,t,bcalc)
				enddo !do s = 1,statesize
			enddo !do k = 1,grid_ss
		enddo !do t = lifespan,1,-1 n
	enddo ! do bcalc = 1,2
	
!	SECOND VALUE FUNCTION FINISHED	
    return
end subroutine valfun

subroutine lc_simulation(cohort,in_asset,in_aime,years_aime,wage,mexp,ms,transmat,vf,vf_na,pf,pf_na,lchoice,lchoice_na,app_policy, &
						life_assets,life_labor,app_age)
!this subroutine simulates the life cycle of an individual given:
!-cohort
!-initial assets
!-initial AIME
!-number of working years up to now (used to calculate AIME)
!-wage matrix
!-value functions
!-labor choice functions
!-policy functions
use parameters
use procedures
integer, intent(in) :: cohort
real (kind = 8), intent(in) :: in_asset, in_aime, years_aime
real (kind = 8), dimension(lifespan,nwagebins), intent(in) :: wage
real (kind = 8), dimension(lifespan,nhealth,nmedbins), intent(in) :: mexp
integer, dimension(statesize,3), intent(in) :: ms
real (kind = 8), dimension(statesize,statesize,lifespan-1), intent(in) :: transmat
real (kind = 8), intent(in) :: vf(grid_asset,statesize,grid_ss,lifespan+1) 
integer, intent(in) :: pf(grid_asset,statesize,grid_ss,lifespan)
real (kind = 8), intent(in) :: lchoice(grid_asset,statesize,grid_ss,lifespan) 
real (kind = 8), dimension(grid_asset,statesize,grid_ss,lifespan+1,2), intent(in) :: vf_na
integer, dimension(grid_asset,statesize,grid_ss,lifespan,2), intent(in) :: pf_na
real (kind = 8), dimension(grid_asset,statesize,grid_ss,lifespan,2), intent(in) :: app_policy
real (kind = 8), dimension(grid_asset,statesize,grid_ss,lifespan,2), intent(in) :: lchoice_na
real (kind = 8), intent(out) :: life_assets(lifespan+1), life_labor(lifespan)
integer, intent(out) :: app_age !the variable wil contain age of social security application



!LOCAL VARIABLES
real (kind = 8) life_earnings(lifespan), life_medexp(lifespan)
real (kind = 8) acur, anext, workyears, app, apn, apx, aime, lpp,lpn,lpx, bcur, vfp, vfn, vfx1,vfx2
real (kind = 8) wage_tn, med_tn !wage and medical expenses drawn from truncated normal
integer iaprev, ianext, ibprev, ibnext, iwprev,iwnext,imprev,imnext !indices of previous and next gridpoints
real (kind = 8) rndnum
real (kind = 8) assets(grid_asset), ss(grid_ss)
integer t, curstate, bcalc, age, seed_size
integer w,m,h !wage, med, health index
integer, allocatable :: seed(:)
integer :: tn_seed = 10000
real (kind = 8) :: tn_sigma = 1.5d0, tn_mu = 0.0d0, tn_samp

!determine initial state
call random_seed(seed_size)
allocate(seed(seed_size))
seed = 100 !fix seed to have reproducible results
call random_seed(put = seed) !start new random sequence
call random_number(rndnum) !get a random number
deallocate(seed)

!initialize vectors
life_earnings 	= -1.0d0
life_assets 	= -1.0d0
life_labor		= -1.0d0
life_medexp 	= -1.0d0

!   create assets grid and social security grid
assets(1)	= 0
ss(1)		= 0
call logspace(asset_min,asset_max,grid_asset-1,assets(2:grid_asset))
call logspace(ss_min,ss_max,grid_ss-1,ss(2:grid_ss))


!test truncated normal

call truncated_normal_ab_sample(tn_mu,tn_sigma,-2.0d0,1.d0,tn_seed,tn_samp)
!------------------------------------------
!For now, initial state is just state 1; 
!I will update it soon
nextstate = 1
!------------------------------------------

acur = in_asset !input asset
aime = in_aime !as well as starting AIME
workyears = years_aime !and years worked
app_age = 1000 !we first set it to very high number, just for initialization


do t = 1,lifespan
	age = t+49
	curstate = nextstate !update state
	!determine wage shock and med shock corresponding to the current state; we're going to need it for AIME update
	!that is, which bin we're at now for wage and for med exp
	w = int(ms(curstate,1))
	m = int(ms(curstate,2))
	h = int(ms(curstate,3))
	!given the bin, make two draws from truncated normal for both wage and med exp
	if (w /= 1 .AND. w /= nwagebins) then
		call truncated_normal_ab_sample(tn_mu,tn_sigma,-2.0d0,1.d0,tn_seed,tn_samp)
	elseif (w == 1) then
	else
	endif
	
	!check the next period benefit calculation regime
	if (workyears<35) then
		bcalc = 1
	else
		bcalc = 2
	endif
	!put current asset into personal history
	life_assets(t) = acur	
	!given current aime, calculate current "potential" benefit bcur (the amount one will receive IF 1) eligible 2)chooses to apply)
	call benefit_calc(cohort,age,aime,bcur)
	!given current benefit, calculate two closest locations on benefit grid
	ibprev = locate_greater(ss,bcur)-1 !index of closest asset on the grid from below
	ibnext = ibprev+1 !index of closest asset on the grid from above
	!given current asset, calculate two closest locations on asset grid
	iaprev = locate_greater(assets,acur)-1 !index of closest asset on the grid from below
	ianext = iaprev+1 !index of closest asset on the grid from above	
	
	!check whether individual applies to social security this period.
	!Three possible cases:
	!1)	individual is eligible but not yet applied (age<app_age) and 
	!	at the same time it is unsure whether to apply or not: app_policy(iaprev,curstate,ibprev,t,bcalc) /= app_policy(ianext,curstate,ibnext,t,bcalc).
	!2) individual definitely applies this period: this means app_policy(iaprev,curstate,ibprev,t,bcalc) = 1 and app_policy(ianext,curstate,ibnext,t,bcalc) = 1.
	!3) individual doesn't apply this period: app_policy(iaprev,curstate,ibprev,t,bcalc) = 0 and app_policy(ianext,curstate,ibnext,t,bcalc) = 0
	!Let's check it case by case.
	if (app_age>age) then !individual not yet applied	
		if (app_policy(iaprev,curstate,ibprev,t,bcalc) /= app_policy(ianext,curstate,ibnext,t,bcalc)) then !1st case
			!we have to check both vf and vf_na at exact point (acur,bcur)
			!first, check vf_na, result stored in vfx1
			vfp = vf_na(iaprev,curstate,ibprev,t,bcalc)
			vfn = vf_na(ianext,curstate,ibnext,t,bcalc)
			vfx1 = vfp + (vfn-vfp)*distance2d(acur,bcur,assets(iaprev),ss(ibprev))/ &
				distance2d(assets(ianext),ss(ibnext),assets(iaprev),ss(iaprev))
			!same for vf, store in vfx2
			vfp = vf(ianext,curstate,ibnext,t)
			vfn = vf(ianext,curstate,ibnext,t)
			vfx2 = vfp + (vfn-vfp)*distance2d(acur,bcur,assets(iaprev),ss(ibprev))/ &
				distance2d(assets(ianext),ss(ibnext),assets(iaprev),ss(iaprev))
			!now check which one is higher; and take the application decision correspondingly
			if (vfx2>=vfx1) then
				app_age = age !age of application is current age
			endif
		elseif (app_policy(iaprev,curstate,ibprev,t,bcalc) == 1 .AND. app_policy(ianext,curstate,ibnext,t,bcalc) == 1) then !2d case
		!2) second case, individual applies
			app_age = age
		endif !if not these cases, individual doesn't apply (3d case)
	endif
	
	!given all this, calculate optimal asset choice next period and labor choice today (using linear interpolation).
	!update relavant variables
	if (age<62 .OR. (age>=62 .AND. age<app_age)) then !before individual becomes eligible to social security OR is eligible but doesn't apply for it
		!For the simulation for age 50-61, we only use value function of "not applied" individual vf_na.
		app = assets(pf_na(iaprev,curstate,ibprev,t,bcalc)) !optimal choice next period on the gridpoint (iaprev,ibprev)
		apn = assets(pf_na(ianext,curstate,ibnext,t,bcalc)) !optimal choice next period on the gridpoint (ianext,ibnext)		
		lpp = lchoice_na(iaprev,curstate,ibprev,t,bcalc) !optimal current labor choice on the gridpoint (iaprev,ibprev)
		lpn = lchoice_na(ianext,curstate,ibnext,t,bcalc) !optimal current labor choice on the gridpoint (ianext,ibnext)
		!--------------------------------------------------------------------------------------------------------------------
		!now use pythagorean theorem and similar triangles to calculate an approximate points in (asset,benefit,asset')-space and (asset,benefit,lchoice)-space (linearly projected)
		!assets next period (linear interpolation)
		apx = app + (apn-app)*distance2d(acur,bcur,assets(iaprev),ss(ibprev))/ &
				distance2d(assets(ianext),ss(ibnext),assets(iaprev),ss(iaprev))	!linear intrapolation between points for the optimal choice next period
		
		!given this, update acur for the next iteration
		acur = apx !notice, that we will put this value into the vector of lifetime assets in the beginning of next cycle
		
		!current optimal labor choice: linear interpolation
		lpx = lpp + (lpn-lpp)*distance2d(acur,bcur,assets(iaprev),ss(ibprev))/ &
				distance2d(assets(ianext),ss(ibnext),assets(iaprev),ss(iaprev)) !linear interpolation for optimal choice of labor given current asset (off the grid)	
		!put this into personal working history
		life_labor(t) = lpx
		!update number of working years (if individual worked)
		if (lpx>0) then
			workyears = workyears+1
		endif	
		
		!--------------------------------------------------------------------------------------------------------------------		
		!We need to update aime
		
		if (bcalc == 1) then !individual worked less than 35 years
			aime = aime + wage(t,w)*lpx/35.0 		      !next period aime grows              
		else !individual aready has more than 35 working years
			aime = aime + max(0.0,(wage(t,w)*lpx-aime)/35.0) !next-period aime doesn't decline
		endif		
		!--------------------------------------------------------------------------------------------------------------------		
	else	
		!individual applied for social security before OR this period: (age>=app_age)
		!thus, for calculations we use vf and pf for those who applied
		app = assets(pf(iaprev,curstate,ibprev,t))
		apn = assets(pf(ianext,curstate,ibnext,t))		
		lpp = lchoice(iaprev,curstate,ibprev,t)
		lpn = lchoice(ianext,curstate,ibnext,t) 
		apx = app + (apn-app)*distance2d(acur,bcur,assets(iaprev),ss(ibprev))/ &
				distance2d(assets(ianext),ss(ibnext),assets(iaprev),ss(iaprev))	!linear intrapolation between points for the optimal choice next period		
		!given this, update acur for the next iteration
		acur = apx	!notice, that we will put this value into the vector of lifetime assets in the beginning of next cycle	
		!current optimal labor choice: linear interpolation
		lpx = lpp + (lpn-lpp)*distance2d(acur,bcur,assets(iaprev),ss(ibprev))/ &
				distance2d(assets(ianext),ss(ibnext),assets(iaprev),ss(iaprev)) !linear interpolation for optimal choice of labor given current asset (off the grid)	
		!put this into personal working history
		life_labor(t) = lpx
		!here, we don't care about total working years and benefit updating anymore
		!--------------------------------------------------------------------------------------------------------------------
	endif
	
	life_earnings(t) = 	wage(t,w)*life_labor(t)
	life_medexp(t) = mexp(t,h,m)
	
	!now, we need to update state, using transition matrix
	call random_number(rndnum)
	nextstate = locate_greater(transmat(curstate,:,t),rndnum) - 1
enddo


end subroutine lc_simulation




