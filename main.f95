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
subroutine lc_simulation(cohort,in_asset,in_aime,years_aime,wage,logwage,wzmean,wzsd,mexp,logmexp,mzmean,mzsd, &
						ms,transmat,vf,vf_na,pf,pf_na,lchoice,lchoice_na,app_policy,life_assets,life_labor,app_age)
	use parameters
	use procedures
	integer, intent(in) :: cohort
	real (kind = 8), intent(in) :: in_asset, in_aime, years_aime
	real (kind = 8), dimension(lifespan,nwagebins), intent(in) :: wage
	real (kind = 8), dimension(lifespan), intent(in) :: wzmean, wzsd,logwage
	real (kind = 8), dimension(lifespan,nhealth,nmedbins), intent(in) :: mexp
	real (kind = 8), dimension(lifespan,nhealth), intent(in) :: mzmean, mzsd,logmexp
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
real (kind = 8) statvec(statesize), cumstvec(statesize) !stationary distribution of states, and cumulative stationary distribution
!-------------------------------
!new variables related to wage
!for each age, each type the entry of logwagedata will contain three values: (1. predicted mean log wage, 2. mean and 3. sd of parenting normal distribution of residuals at this age)
!Thus, to compute a wage at particular age, I only need to know:
!	1.  predicted mean log wage
!	2.	bin we're at;
!	4.	rule of how bins are organised : -inf<mu-1.5sd<mu-0.5sd<m+0.5sd<mu+1.5sd<+inf; 
!	3.	mean and as of parenting normal. Then, given the rule the bins are organised as above, I make a necessary draw from truncated normal!
real (kind = 8) logwagedata(lifespan,3,ntypes), logwage(lifespan,ntypes), wzmean(lifespan,ntypes), wzsd(lifespan,ntypes) 
real (kind = 8) wagegrid(lifespan,nwagebins,ntypes) !wage grid points
real (kind = 8) wagetrans(nwagebins,nwagebins,lifespan-1) !transition matrix between wage bins
!-------------------------------
!new variables related to medical expenditure
!for each age, each type, and health the entry of logmexpdata will contain three values: (1. predicted mean log wage, 2. mean and 3. sd of parenting normal distribution of residuals at this age)
!Thus, to compute a wage at particular age, I only need to know:
!	1.  predicted mean log wage
!	2.	bin we're at;
!	4.	rule of how bins are organised : -inf<mu-1.5sd<mu-0.5sd<m+0.5sd<mu+1.5sd<+inf; 
!	3.	mean and as of parenting normal. Then, given the rule the bins are organised as above, I make a necessary draw from truncated normal!
real (kind = 8) logmexpdata(lifespan,3,nhealth,ntypes), logmexp(lifespan,nhealth,ntypes), &
				mzmean(lifespan,nhealth,ntypes), mzsd(lifespan,nhealth,ntypes) 
real (kind = 8) mexpgrid(lifespan,nhealth,nmedbins,ntypes)
real (kind = 8) medtrans(nmedbins,nmedbins,lifespan-1)
!-------------------------------
real (kind = 8) bin_mean !intermediate variable
!-------------------------------

real (kind = 8), dimension(8,3) :: typemat
real (kind = 8) healthtrans(lifespan-1,nhealth,ntypes), surv(lifespan,nhealth,ntypes)
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

!-----------------------------------------------------------------------
!1. Import files

!1.1 Import transition matrices
!1.1.1 Wage transitions
datafile = 'data_inputs/transmat_wage.csv'
call import_matrices(datafile,wagetrans)
!1.1.2 Med expenses transitions
datafile = 'data_inputs/transmat_mexp.csv'
call import_matrices(datafile,medtrans)
!1.1.3 Health transitions, survival and type matrix
do cohort = 0,1
	do college = 0,1
		do sex = 0,1        
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

!2.	Construct wage and medical expenses grids
!Import a wage profile, profile of sd and mean of parenting normal distribution (which we are going to use to draw a wage)
!		And med exp profile, profile of sd and mean of parenting normal distribution (which we are going to use to draw a med exp), additionally conditional on health
do i = 1,ntypes
	!import from a proper file
	WRITE(datafile,'(a,i1,a)') "data_inputs/lwage_smooth_",i,".csv"
	call import_data(datafile,logwagedata(:,:,i))
	logwage(:,i) = logwagedata(:,1,i)
	wzmean(:,i) = logwagedata(:,2,i)
	wzsd(:,i) = logwagedata(:,3,i)
	!construct a wage grid: 
	!for every type, for every age, construct a mean of corresponding wage bin from truncated normal distribition
	!remember, residuals are distributed normally, but they are devidations from the log mean
	!we do it in a few steps:
	do k = 1,nwagebins !for every bin
		do t = 1,lifespan !for every age
			!1. calculate a mean of residual in that bin
			if (k == 1) then !first bin
				call truncated_normal_b_mean (wzmean(t,i), wzsd(t,i), wzmean(t,i)-1.5*wzsd(t,i),bin_mean) !calculate mean of the bin
			elseif (k == nwagebins) then !last bin
				call truncated_normal_a_mean (wzmean(t,i), wzsd(t,i), wzmean(t,i)+1.5*wzsd(t,i),bin_mean)
			else !interior bins
				call truncated_normal_ab_mean (wzmean(t,i), wzsd(t,i), wzmean(t,i)+wzsd(t,i)*(k-3.5), &
												wzmean(t,i)+wzsd(t,i)*(k-2.5), bin_mean)
			endif
			wagegrid(t,k,i) = exp(logwage(t,i)+bin_mean)  !calculate real wage (exponent of mean log wage plus bin average)
		enddo 
	enddo
	
!Same for medical expenses, additionally condition on health	
	do j = 1,nhealth
		WRITE(datafile,'(a,i1,i1,a)') "data_inputs/lmexp_smooth_",i,j-1,".csv"
		call import_data(datafile,logmexpdata(:,:,j,i))
		logmexp(:,j,i) = logmexpdata(:,1,j,i)
		mzmean(:,j,i) = logmexpdata(:,2,j,i)
		mzsd(:,j,i) = logmexpdata(:,3,j,i)
		do k = 1,nmedbins !for every bin
			do t = 1,lifespan !for every age
				!1. calculate a mean of residual in that bin
				if (k == 1) then !first bin
					call truncated_normal_b_mean (mzmean(t,j,i), mzsd(t,j,i), mzmean(t,j,i)-1.5*mzsd(t,j,i),bin_mean) !calculate mean of the bin
				elseif (k == nmedbins) then !last bin
					call truncated_normal_a_mean (mzmean(t,j,i), mzsd(t,j,i), mzmean(t,j,i)+1.5*mzsd(t,j,i),bin_mean)
				else !interior bins
					call truncated_normal_ab_mean (mzmean(t,j,i), mzsd(t,j,i), mzmean(t,j,i)+mzsd(t,j,i)*(k-3.5), &
													mzmean(t,j,i)+mzsd(t,j,i)*(k-2.5), bin_mean)
				endif
				mexpgrid(t,j,k,i) = exp(logmexp(t,j,i)+bin_mean)  !calculate real med exp (exponent of mean log med exp plus bin average)
			enddo
		enddo
	enddo
enddo

!-----------------------------------------------------------------------
!3. Construct matrix of the correspondence of exostate to shocks
do w = 1,nwagebins
    do m = 1,nmedbins
        do h = 1,nhealth
            k = h+nhealth*(m-1)+nmedbins*nhealth*(w-1)
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

!5. Pick a type between 1 and 8, and calculate optimal labor choice matrix and value function for a chosen type (choice have been made in the beginning)
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
print *, "Chosen type:", ind_type
cohort = (ind_type-1)/4 + 1 !cohort is 1/2 and not 0/1 because of indexing style; it is passed to valfun, where it is passed to aime_calc and benefit_calc
!sex = mod(ind_type,2)
!Calculate stationary distribution of states at age 50 for that type
call stationary_dist(transmat(:,:,1,ind_type),statvec)
call cumsum(statvec,cumstvec)
!labor choice
lchoice = -10
!value function
call valfun(cohort,wagegrid(:,:,ind_type),mexpgrid(:,:,:,ind_type),transmat(:,:,:,ind_type),surv(:,:,ind_type),ms, &
vf,vf_na,pf,pf_na,lchoice,lchoice_na,app_policy)

call lc_simulation(cohort,1.5d4,0.8d3,2.5d1, wagegrid(:,:,ind_type),logwage(:,ind_type),wzmean(:,ind_type),wzsd(:,ind_type), &
mexpgrid(:,:,:,ind_type),logmexp(:,:,ind_type),mzmean(:,:,ind_type),mzsd(:,:,ind_type), ms, &
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

!	print *, size(wage,1)
!	print *, size(wage,2)
!	print *, size(medexp,1)
!	print *, size(medexp,2)
!	print *, size(medexp,3)
!	print *, size(transmat,1)
!	print *, size(transmat,2)
!	print *, size(transmat,3)
!	print *, size(vf,1)
!	print *, size(vf,2)
!	print *, size(vf,3)
!	print *, size(vf,4)
!	print *, size(pf,1)
!	print *, size(pf,2)
!	print *, size(pf,3)
!	print *, size(pf,4)
	
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
    elseif (size(medexp,2) .NE. nhealth) then
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
    
    vf = -101.0d0
    pf = -102.0d0
    
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

subroutine lc_simulation(cohort,in_asset,in_aime,years_aime,wage,logwage,wzmean,wzsd,mexp,logmexp,mzmean,mzsd, &
						ms,transmat,vf,vf_na,pf,pf_na,lchoice,lchoice_na,app_policy,life_assets,life_labor,app_age)
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
implicit none
integer, intent(in) :: cohort
real (kind = 8), intent(in) :: in_asset, in_aime, years_aime
real (kind = 8), dimension(lifespan,nwagebins), intent(in) :: wage
real (kind = 8), dimension(lifespan), intent(in) :: wzmean, wzsd,logwage
real (kind = 8), dimension(lifespan,nhealth,nmedbins), intent(in) :: mexp
real (kind = 8), dimension(lifespan,nhealth), intent(in) :: mzmean, mzsd,logmexp
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
real (kind = 8) life_earnings(lifespan), life_wage(lifespan), life_medexp(lifespan)
real (kind = 8) anext, workyears, aime, apx, lpx, acur, wcur, mcur, bcur, vx1,vx2
real (kind = 8) wage_tn, med_tn !wage and medical expenses drawn from truncated normal
integer iaprev, ianext, ibprev, ibnext
integer iwprev,iwnext,imprev,imnext !indices of previous and next gridpoints in (w,m) - space


real (kind = 8) rndnum, statvec(statesize), cumstvec(statesize)
real (kind = 8) assets(grid_asset), ss(grid_ss)
integer t, curstate, nextstate , bcalc, age, seed_size
integer w,m,h !wage, med, health index
integer, allocatable :: seed(:)
integer :: tn_seed = 10000
real (kind = 8) :: tn_sigma = 1.5d0, tn_mu = 0.0d0, tn_samp

!Value function polygone vertices for linear interpolation for vf, assets, labor
real (kind = 8) v0000,v0001,v0010,v0011,v0100,v0101,v0110,v0111, &
				v1000,v1001,v1010,v1011,v1100,v1101,v1110,v1111
real (kind = 8) a0000,a0001,a0010,a0011,a0100,a0101,a0110,a0111, &
				a1000,a1001,a1010,a1011,a1100,a1101,a1110,a1111		
real (kind = 8) l0000,l0001,l0010,l0011,l0100,l0101,l0110,l0111, &
				l1000,l1001,l1010,l1011,l1100,l1101,l1110,l1111				
!Value function polygone vertices for linear interpolation for vf_na, assets_na, labor_na		
real (kind = 8) vn0000,vn0001,vn0010,vn0011,vn0100,vn0101,vn0110,vn0111, &
				vn1000,vn1001,vn1010,vn1011,vn1100,vn1101,vn1110,vn1111
real (kind = 8) an0000,an0001,an0010,an0011,an0100,an0101,an0110,an0111, &
				an1000,an1001,an1010,an1011,an1100,an1101,an1110,an1111		
real (kind = 8) ln0000,ln0001,ln0010,ln0011,ln0100,ln0101,ln0110,ln0111, &
				ln1000,ln1001,ln1010,ln1011,ln1100,ln1101,ln1110,ln1111	
!Corresponding state values in s-space
integer state00, state01,state10,state11

!determine initial state
call random_seed(seed_size)
allocate(seed(seed_size))
seed = 100 !fix seed to have reproducible results
call random_seed(put = seed) !start new random sequence
deallocate(seed)

!initialize vectors
life_earnings 	= -1.0d0
life_assets 	= -1.0d0
life_labor		= -1.0d0
life_medexp 	= -1.0d0
life_wage		= -1.0d0

!   create assets grid and social security grid
assets(1)	= 0
ss(1)		= 0
call logspace(asset_min,asset_max,grid_asset-1,assets(2:grid_asset))
call logspace(ss_min,ss_max,grid_ss-1,ss(2:grid_ss))


!test truncated normal

!call truncated_normal_ab_sample(tn_mu,tn_sigma,-2.0d0,1.d0,tn_seed,tn_samp)
!------------------------------------------
!calculate initial stationry distribution of states at age 50 and determine initial state by random draw
call stationary_dist(transmat(:,:,1),statvec)
call cumsum(statvec,cumstvec)
call random_number(rndnum) !get a random number from uniform (0,1) distribution
nextstate = locate_greater(cumstvec,rndnum) !see where this number falls in cumulative current distribution
!------------------------------------------

acur = in_asset !input asset
aime = in_aime !as well as starting AIME
workyears = years_aime !and years worked
app_age = 1000 !we first set it to very high number, just for initialization


do t = 1,lifespan
	age = t+49
	curstate = nextstate !update state
	
	!	A) ASSETS HISTORY UPDATE
	!put current asset into personal history
	life_assets(t) = acur
	
	!determine wage shock and med shock corresponding to the current state; we're going to need it for AIME update
	!that is, which bin we're at now for wage and for med exp
	w = int(ms(curstate,1))
	m = int(ms(curstate,2))
	h = int(ms(curstate,3))
	!given the bin, make two draws from truncated normal for both wage and med exp
	!First, draw wage residual
	if (w /= 1 .AND. w /= nwagebins) then
		call truncated_normal_ab_sample(wzmean(t),wzsd(t),wzmean(t)+wzsd(t)*(w-3.5), &
										wzmean(t)+wzsd(t)*(w-2.5),tn_seed,wage_tn)
	elseif (w == 1) then
		call truncated_normal_ab_sample(wzmean(t),wzsd(t),wage(t,1), &
										wzmean(t)-1.5*wzsd(t),tn_seed,wage_tn)
	elseif (w == nwagebins) then
		call truncated_normal_ab_sample(wzmean(t),wzsd(t),wzmean(t)+1.5*wzsd(t), &
										wage(t,nwagebins),tn_seed,wage_tn)
	else
		print *, "Something's wrong with current wage state"
		read (*,*)
		stop
	endif
	! Calculate wage
	wcur = exp(logwage(t)+wage_tn)
	
	!	B) WAGE HISTORY UPDATE
	life_wage(t) = wcur
	
	! Find nearest wage grid points
	if (life_wage(t)>=wage(t,w)) then
		iwprev = w
		iwnext = w + 1
	else
		iwprev = w - 1
		iwnext = w
	endif
	! Second draw med exp residuals
	if (m /= 1 .AND. m /= nmedbins) then
		call truncated_normal_ab_sample(mzmean(t,h),mzsd(t,h),mzmean(t,h)+mzsd(t,h)*(m-3.5), &
										mzmean(t,h)+mzsd(t,h)*(m-2.5),tn_seed,med_tn)
	elseif (m == 1) then
		call truncated_normal_ab_sample(mzmean(t,h),mzsd(t,h),mexp(t,h,1), &
										mzmean(t,h)-1.5*mzsd(t,h),tn_seed,med_tn)
	elseif (m == nmedbins) then
		call truncated_normal_ab_sample(mzmean(t,h),mzsd(t,h),mzmean(t,h)+1.5*mzsd(t,h), &
										mexp(t,h,nmedbins),tn_seed,med_tn)
	else
		print *, "Something's wrong with current med exp state"
		read (*,*)
		stop
	endif
	! Calculate med exp
	mcur = exp(logmexp(t,h)+med_tn)
	
	!	C) MEDICAL EXPENSES HISTORY UPDATE
	life_medexp(t) = mcur
	
	! Nearest med grid points
	if (life_medexp(t)>=mexp(t,h,m)) then
		imprev = m
		imnext = m + 1
	else
		imprev = m - 1
		imnext = m
	endif
	!check the next period benefit calculation regime
	if (workyears<35) then
		bcalc = 1
	else
		bcalc = 2
	endif
		
	!given current aime, calculate current "potential" benefit bcur (the amount one will receive IF 1) eligible 2)chooses to apply)
	call benefit_calc(cohort,age,aime,bcur)
	!given current benefit, calculate two closest locations on benefit grid
	ibprev = locate_greater(ss,bcur)-1 !index of closest asset on the grid from below
	ibnext = ibprev+1 !index of closest asset on the grid from above
	!given current asset, calculate two closest locations on asset grid
	iaprev = locate_greater(assets,acur)-1 !index of closest asset on the grid from below
	ianext = iaprev+1 !index of closest asset on the grid from above	
	
	
	!	Now we have two ends of 4-tuple: (iaprev,ibprev,iwprev,imprev); (ianext,ibnext,iwnext,imnext)
	!	However, value function at a given period t is 3-dimensional: 
	!	(imprev, iwprev) -> (imprev, iwprev,h) -> sprev 
	!	(imnext, iwnext) -> (imnext, iwnext,h) -> snext 
	!	Therefre, I have to be VERY careful when calculating 4-tuple: i have to go from v(a,s,b) to v(a,w,m,b)
	!	Therefore we can calculate value functons v(iaprev,sprev,ibprev,t) == v0000 == v(iaprev,(iwprev,imprev),ibprev,t) 
	!	to v(ianext,snext,ibnext,t) == v1111 == v(ianext,(iwnext,imnext),ibnext,t)
	!	ON TOP OF THAT, we have TWO value functions, and we have to calculate it for BOTH, that's why I have vn0000 to vn1111
	!	LET's DO IT EXTRA CAREFULLY
	!	16 value functions v0000 to v1111
	!	and 16 vf_na: vn0000 to vn1111
	
	!states (0,0) -> (iwprev,imprev); (1,1) -> (iwnext,imnext)
	state00 = h+nhealth*(imprev-1)+nmedbins*nhealth*(iwprev-1)
	state01 = h+nhealth*(imnext-1)+nmedbins*nhealth*(iwprev-1)
	state10 = h+nhealth*(imprev-1)+nmedbins*nhealth*(iwnext-1)
	state11 = h+nhealth*(imnext-1)+nmedbins*nhealth*(iwnext-1)
		!0000
	v0000 	= vf			(iaprev,state00,ibprev,t)
	vn0000 	= vf_na			(iaprev,state00,ibprev,t,bcalc)
	a0000 	= assets(pf		(iaprev,state00,ibprev,t))
	an0000 	= assets(pf_na	(iaprev,state00,ibprev,t,bcalc))
	l0000 	= lchoice		(iaprev,state00,ibprev,t)
	ln0000 	= lchoice_na	(iaprev,state00,ibprev,t,bcalc)
		!0001
	v0001 	= vf			(iaprev,state00,ibnext,t)
	vn0001 	= vf_na			(iaprev,state00,ibnext,t,bcalc)
	a0001 	= assets(pf		(iaprev,state00,ibnext,t))
	an0001 	= assets(pf_na	(iaprev,state00,ibnext,t,bcalc))
	l0001 	= lchoice		(iaprev,state00,ibnext,t)
	ln0001 	= lchoice_na	(iaprev,state00,ibnext,t,bcalc)
		!0010
	v0010	= vf			(iaprev,state01,ibprev,t)
	vn0010	= vf_na			(iaprev,state01,ibprev,t,bcalc)
	a0010 	= assets(pf		(iaprev,state01,ibprev,t))
	an0010 	= assets(pf_na	(iaprev,state01,ibprev,t,bcalc))
	l0010 	= lchoice		(iaprev,state01,ibprev,t)
	ln0010 	= lchoice_na	(iaprev,state01,ibprev,t,bcalc)
		!0011
	v0011	= vf			(iaprev,state01,ibnext,t)
	vn0011	= vf_na			(iaprev,state01,ibnext,t,bcalc)
	a0011 	= assets(pf		(iaprev,state01,ibnext,t))
	an0011 	= assets(pf_na	(iaprev,state01,ibnext,t,bcalc))
	l0011 	= lchoice		(iaprev,state01,ibnext,t)
	ln0011 	= lchoice_na	(iaprev,state01,ibnext,t,bcalc)
		!0100
	v0100	= vf			(iaprev,state10,ibprev,t)
	vn0100	= vf_na			(iaprev,state10,ibprev,t,bcalc)
	a0100 	= assets(pf		(iaprev,state10,ibprev,t))
	an0100 	= assets(pf_na	(iaprev,state10,ibprev,t,bcalc))
	l0100 	= lchoice		(iaprev,state10,ibprev,t)
	ln0100 	= lchoice_na	(iaprev,state10,ibprev,t,bcalc)
		!0101
	v0101	= vf			(iaprev,state10,ibnext,t)
	vn0101	= vf_na			(iaprev,state10,ibnext,t,bcalc)
	a0101 	= assets(pf		(iaprev,state10,ibnext,t))
	an0101 	= assets(pf_na	(iaprev,state10,ibnext,t,bcalc))
	l0101 	= lchoice		(iaprev,state10,ibnext,t)
	ln0101 	= lchoice_na	(iaprev,state10,ibnext,t,bcalc)
		!0110
	v0110	= vf			(iaprev,state11,ibprev,t)
	vn0110	= vf_na			(iaprev,state11,ibprev,t,bcalc)
	a0110 	= assets(pf		(iaprev,state11,ibprev,t))
	an0110 	= assets(pf_na	(iaprev,state11,ibprev,t,bcalc))
	l0110 	= lchoice		(iaprev,state11,ibprev,t)
	ln0110 	= lchoice_na	(iaprev,state11,ibprev,t,bcalc)
		!0111
	v0111	= vf			(iaprev,state11,ibnext,t)
	vn0111	= vf_na			(iaprev,state11,ibnext,t,bcalc)
	a0111 	= assets(pf		(iaprev,state11,ibnext,t))
	an0111 	= assets(pf_na	(iaprev,state11,ibnext,t,bcalc))
	l0111 	= lchoice		(iaprev,state11,ibnext,t)
	ln0111 	= lchoice_na	(iaprev,state11,ibnext,t,bcalc)
		!1000
	v1000 	= vf			(ianext,state00,ibprev,t)
	vn1000 	= vf_na			(ianext,state00,ibprev,t,bcalc)
	a1000 	= assets(pf		(ianext,state00,ibprev,t))
	an1000 	= assets(pf_na	(ianext,state00,ibprev,t,bcalc))
	l1000 	= lchoice		(ianext,state00,ibprev,t)
	ln1000 	= lchoice_na	(ianext,state00,ibprev,t,bcalc)
		!1001
	v1001 	= vf			(ianext,state00,ibnext,t)
	vn1001 	= vf_na			(ianext,state00,ibnext,t,bcalc)
	a1001 	= assets(pf		(ianext,state00,ibnext,t))
	an1001 	= assets(pf_na	(ianext,state00,ibnext,t,bcalc))
	l1001 	= lchoice		(ianext,state00,ibnext,t)
	ln1001 	= lchoice_na	(ianext,state00,ibnext,t,bcalc)
		!1010
	v1010	= vf			(ianext,state01,ibprev,t)
	vn1010	= vf_na			(ianext,state01,ibprev,t,bcalc)
	a1010 	= assets(pf		(ianext,state01,ibprev,t))
	an1010 	= assets(pf_na	(ianext,state01,ibprev,t,bcalc))
	l1010 	= lchoice		(ianext,state01,ibprev,t)
	ln1010 	= lchoice_na	(ianext,state01,ibprev,t,bcalc)
		!1011
	v1011	= vf			(ianext,state01,ibnext,t)
	vn1011	= vf_na			(ianext,state01,ibnext,t,bcalc)
	a1011 	= assets(pf		(ianext,state01,ibnext,t))
	an1011 	= assets(pf_na	(ianext,state01,ibnext,t,bcalc))
	l1011 	= lchoice		(ianext,state01,ibnext,t)
	ln1011 	= lchoice_na	(ianext,state01,ibnext,t,bcalc)
		!1100
	v1100	= vf			(ianext,state10,ibprev,t)
	vn1100	= vf_na			(ianext,state10,ibprev,t,bcalc)
	a1100 	= assets(pf		(ianext,state10,ibprev,t))
	an1100 	= assets(pf_na	(ianext,state10,ibprev,t,bcalc))
	l1100 	= lchoice		(ianext,state10,ibprev,t)
	ln1100 	= lchoice_na	(ianext,state10,ibprev,t,bcalc)
		!1101
	v1101	= vf			(ianext,state10,ibnext,t)
	vn1101	= vf_na			(ianext,state10,ibnext,t,bcalc)
	a1101 	= assets(pf		(ianext,state10,ibnext,t))
	an1101 	= assets(pf_na	(ianext,state10,ibnext,t,bcalc))
	l1101 	= lchoice		(ianext,state10,ibnext,t)
	ln1101 	= lchoice_na	(ianext,state10,ibnext,t,bcalc)
		!1110
	v1110	= vf			(ianext,state11,ibprev,t)
	vn1110	= vf_na			(ianext,state11,ibprev,t,bcalc)
	a1110 	= assets(pf		(ianext,state11,ibprev,t))
	an1110 	= assets(pf_na	(ianext,state11,ibprev,t,bcalc))
	l1110 	= lchoice		(ianext,state11,ibprev,t)
	ln1110 	= lchoice_na	(ianext,state11,ibprev,t,bcalc)
		!1111
	v1111	= vf			(ianext,state11,ibnext,t)
	vn1111	= vf_na			(ianext,state11,ibnext,t,bcalc)
	a1111 	= assets(pf		(ianext,state11,ibnext,t))
	an1111 	= assets(pf_na	(ianext,state11,ibnext,t,bcalc))
	l1111 	= lchoice		(ianext,state11,ibnext,t)
	ln1111 	= lchoice_na	(ianext,state11,ibnext,t,bcalc)
	
	!check whether individual is: 
	!1) age<62: can't apply; use only vf_na (not vf)
	!2) age>=62 AND age<app_age: not yet applied, but we must check whether applies: compare values of two quadrilinear interpolations and use higher one (change app_age if it is the case)
	!3) age>=62 AND age>age_app: already applied, use only vf (not vf_na)
	
	!Before calculating interpolated assets and labor choice, we check if individual is eligible and didn't apply yet. If it's the case, we
	!switch app_age to a current age
	
	if (age>=62 .AND. age<app_age) then 
		!Here we have to check both cases: applicant and non_applicant; 
		!if applicant's value is higher, she/he applies, and we swithc app_age to current age and use applicants vf and lchoice;
		!otherwise, we use vf_na
		!(x,y,z,q) == (a,w,m,b)
		vx1 = quadlin(	acur,wcur,mcur,bcur, &
						assets(iaprev),wage(t,iwprev),mexp(t,h,imprev),ss(ibprev), &
						assets(ianext),wage(t,iwnext),mexp(t,h,imnext),ss(ibnext), &
						v0000,v0001,v0010,v0011,v0100,v0101,v0110,v0111, &
						v1000,v1001,v1010,v1011,v1100,v1101,v1110,v1111)
		vx2 = quadlin(	acur,wcur,mcur,bcur, &
						assets(iaprev),wage(t,iwprev),mexp(t,h,imprev),ss(ibprev), &
						assets(ianext),wage(t,iwnext),mexp(t,h,imnext),ss(ibnext), &
						vn0000,vn0001,vn0010,vn0011,vn0100,vn0101,vn0110,vn0111, &
						vn1000,vn1001,vn1010,vn1011,vn1100,vn1101,vn1110,vn1111)
		if (vx1 >= vx2) then !it is time to apply for ss
			!SWITCH APPLICATION AGE, MAIN THING HERE!!!
			app_age = age 
		endif !otherwise do nothing: individual didn't apply
	endif
	
	!NOW, we calculate interpolated assets and labor
	if (age<62 .OR. (age>=62 .AND. age<app_age)) then !individual 1) can't apply or 2) can but doesn't (we checked it before)
		!Here, we only use vf_na to construct current assets and labor choice; quadlin
		!(x,y,z,q) == (a,w,m,b)
		apx = quadlin(	acur,wcur,mcur,bcur, &
						assets(iaprev),wage(t,iwprev),mexp(t,h,imprev),ss(ibprev), &
						assets(ianext),wage(t,iwnext),mexp(t,h,imnext),ss(ibnext), & 
						an0000,an0001,an0010,an0011,an0100,an0101,an0110,an0111, &
						an1000,an1001,an1010,an1011,an1100,an1101,an1110,an1111)
		lpx = quadlin(	acur,wcur,mcur,bcur, &
						assets(iaprev),wage(t,iwprev),mexp(t,h,imprev),ss(ibprev), &
						assets(ianext),wage(t,iwnext),mexp(t,h,imnext),ss(ibnext), & 
						ln0000,ln0001,ln0010,ln0011,ln0100,ln0101,ln0110,ln0111, &
						ln1000,ln1001,ln1010,ln1011,ln1100,ln1101,ln1110,ln1111)
						
		!we have to update AIME
		if (bcalc == 1) then !individual worked less than 35 years
			aime = aime + wage(t,w)*lpx/35.0 		      !next period aime grows              
		else !individual aready has more than 35 working years
			aime = aime + max(0.0,(wage(t,w)*lpx-aime)/35.0) !next-period aime doesn't decline
		endif
		
		!update number of working years (if individual worked)
		if (lpx>0) then
			workyears = workyears+1
		endif
	elseif (age>=62 .AND. age>=app_age) then !individual applied		
		!calculate interpolated assets and labor choice using functions of applicants
		apx = quadlin(	acur,wcur,mcur,bcur, &
						assets(iaprev),wage(t,iwprev),mexp(t,h,imprev),ss(ibprev), &
						assets(ianext),wage(t,iwnext),mexp(t,h,imnext),ss(ibnext), & 
						a0000,a0001,a0010,a0011,a0100,a0101,a0110,a0111, &
						a1000,a1001,a1010,a1011,a1100,a1101,a1110,a1111)
		lpx = quadlin(	acur,wcur,mcur,bcur, &
						assets(iaprev),wage(t,iwprev),mexp(t,h,imprev),ss(ibprev), &
						assets(ianext),wage(t,iwnext),mexp(t,h,imnext),ss(ibnext), & 
						l0000,l0001,l0010,l0011,l0100,l0101,l0110,l0111, &
						l1000,l1001,l1010,l1011,l1100,l1101,l1110,l1111)
		!In this case, we don't care about the AIME and working years anymore, it doesn't play any role
	else !Excluded case, just fishing for errors
		print *, "Wrong scenario, check SS applicaton code in lc_simulation"
		stop 
	endif
	
	!given this, update acur for the next iteration
	acur = apx	!notice, that we will put this value into the vector of lifetime assets in the beginning of next cycle	
	!put labor into personal working history
	
	!	D) LABOR HISTORY UPDATE
	life_labor(t) = lpx
	!Fill in individual history	
	
	!	E) EARNINGS HISTORY UPDATE
	life_earnings(t) = 	wcur*life_labor(t)
	
	
	!now, we need to update state, using transition matrix
	call random_number(rndnum)
	!nextstate = locate_greater(transmat(curstate,:,t),rndnum) - 1 !wrong, I have to use cumulative sum
enddo


end subroutine lc_simulation




