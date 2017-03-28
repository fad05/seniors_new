! Updated Fortran code for "seniors" project
! By Alexey Filatov
program main

use parameters
use procedures
use csv_file
use random
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
subroutine labchoice (lchoice, cons, util, acur,anext, h, wage, mexp,ss,age,cohort)
	use procedures
    use parameters
    implicit none
    real (kind = 8), intent(in) :: acur, anext, wage,mexp,ss !current assets, next period assets, current wage, mexp, social scurity
    integer, intent(in) :: h, age, cohort !health status, age and cohort
    real (kind = 8), intent(out) :: lchoice !optimal labor choice
    real (kind = 8), intent(out) :: cons !consumption under optimal choice
    real (kind = 8), intent(out) :: util !value of utility under optimal choice
end subroutine labchoice
subroutine lc_simulation(cohort, in_asset,in_aime,in_wage,in_mexp,in_health,years_aime, &
						wage,logwage,wzmean,wzsd,wbinborders,mexp,logmexp,mzmean,mzsd, &
						mbinborders,ms,transmat,vf,vf_na,pf,pf_na,lchoice,lchoice_na,app_policy, surv, &
						life_assets,life_labor,life_earnings, life_wage, life_medexp,life_health,app_age,tn_seed)
	use parameters
	use procedures
	integer, intent(in)	:: cohort, in_health
	real (kind = 8)		:: in_asset, in_aime,in_wage,in_mexp
	integer, intent(in) :: years_aime
	real (kind = 8), dimension(lifespan,nwagebins), intent(in) :: wage
	real (kind = 8), dimension(lifespan), intent(in) :: wzmean, wzsd,logwage
	real (kind = 8), intent(in) :: wbinborders(lifespan,nwagebins+1), mbinborders(lifespan,nmedbins+1)
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
	real (kind = 8), dimension(lifespan,nhealth), intent(in) :: surv
	real (kind = 8), dimension(grid_asset,statesize,grid_ss,lifespan,2), intent(in) :: lchoice_na
	real (kind = 8), intent(out) :: life_assets(lifespan)
	real (kind = 8), dimension(lifespan), intent(out) :: life_labor, life_earnings, life_wage, life_medexp, life_health
	integer, intent(out) :: app_age
	integer, intent(inout) :: tn_seed

end subroutine lc_simulation

subroutine amoeba(p,y,ftol,func,iter,info)
	use nrtype; use nrutil, only: assert_eq, imaxloc, iminloc, nrerror, swap
	implicit none
	integer(I4B), intent(out) :: iter,info
	real(DP), intent(in) :: ftol
	real(DP), dimension(:), intent(inout) :: y
	real(DP), dimension(:,:), intent(inout) :: p
	interface
		function func(x)
		use nrtype
		implicit none
		real(DP), dimension(:), intent(inout) :: x
		real(DP) :: func
		end function func
	end interface
end subroutine amoeba
end interface

!PROGRAM VARIABLES

!	1.
real (kind = 8) structparams(nparams_cd), spsimplex(nparams_cd+1,nparams_cd) !STRUCTURAL PARAMETERS OF THE MODEL, THE ONES TO BE CALIBRATED; and corresponding simplex for Nelder-Mead algorithm
!	The parameters are:	ltot kappa d cmin beta delta eta sigma cgamma lgamma xi nu
real (kind = 8) y(nparams_cd), funevals(nparams_cd+1)
real (kind = 8), parameter :: usual_delta =  0.1, zero_term_delta =  0.025
integer iter, info

!	2.
character*200 datafile !address of the data file on the disk
character*200 command !bash command
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
real (kind = 8) lwagetemp(16,ntypes)
real (kind = 8) wagegrid(lifespan,nwagebins,ntypes), wgupd(lifespan, nwagebins) !wage grid points, on bin means
real (kind = 8) wagetrans(nwagebins,nwagebins,lifespan-1) !transition matrix between wage bins
real (kind = 8) wbinborders(lifespan,nwagebins+1) !n bins => n+1 borders
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
real (kind = 8) mexpgrid(lifespan,nhealth,nmedbins,ntypes) !med grid, on bin means
real (kind = 8) medtrans(nmedbins,nmedbins,lifespan-1)
real (kind = 8) mbinborders(lifespan,nmedbins+1) !n bins => n+1 borders

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
!results of simulation: single simulation and averaged over "nsim" simulations
real (kind = 8)  life_assets(lifespan), mean_lifeassets(lifespan)
real (kind = 8), dimension(lifespan) :: life_labor, life_earnings, life_wage, life_medexp, life_lfp, worker_wage, life_health
real (kind = 8), dimension(lifespan) :: mean_lifelabor,mean_lifeearnings,mean_lifewage,mean_lifemedexp,mean_lfp,mean_workerwage, &
										mean_applied
real (kind = 8), dimension(lifespan) :: delta_wage, initial_wage, updated_wage, upd_logwage
integer, 		 dimension(lifespan) :: nalive, nwork !for every age, number of individuals alive at that period (to average properly) and no of indivs who work
integer app_age


integer i,j,k,l,n,t,w,wn,m,mn,h,hn,ind_type, ind_subtype, counter, itype
integer :: icoh, cohort, sex, college, health, age(lifespan) = (/(i,i=50,90)/)
integer ms(statesize,3)

!variables necessary to use "random" module;
!to draw from bivariate normal
real fchol3(9), fchol4(16) !LOWER TRIANGULAR DECOMPOSITION OF VARIANCE MATRIX, 3-dim and 4-dim
real draw3(3), draw4(4) !a draw from mv normal
integer ier !    ier = 1 if the input covariance matrix is not +ve definite
!    FIRST = .TRUE. IF THIS IS THE FIRST CALL OF THE ROUTINE
!    OR IF THE DISTRIBUTION HAS CHANGED SINCE THE LAST CALL OF THE ROUTINE.
!    OTHERWISE SET TO .FALSE.
!            (INPUT,LOGICAL)
logical first

!Variables necessary to make random draws
integer :: tn_seed
integer seed_size
integer, allocatable :: seed(:)
real (kind = 8) rndnum

real (kind = 8) :: 	in_aime, in_asset, in_wage, in_mexp ! initial values of aime, assets, wage and mexp drawn from distribution
integer				in_health

!variables that contain parameters of joint distributions of AIME, assets, wage and medical expenses (conditional on asses being zero)
real (kind = 8) prob_ah(4,4), cdf_ah(4,4) 		!4 = male/femele * college/noncollege; 4 = a0h0/a0h1/a1h0/a1h1
real (kind = 8) joint_awma0(4,9) 				!4 = male/femele * college/noncollege; 9 = mu_a, mu_w, mu_m + varcov mat (6 values)
real (kind = 8) joint_awma(4,14)				!4 = male/femele * college/noncollege; 14 = mu_a, mu_w, mu_m, mu_assets + varcov mat (10 values)

!Variables containing data moments
real (kind = 8), dimension(:), allocatable ::	datalfp, datahours, datassets, dataret, sdlfp, sdhours,sdassets,sdret
real (kind = 8), dimension(:), allocatable ::	modellfp, modelhours, modelassets, modelret
real (kind = 8), dimension(:,:), allocatable ::	moments

!Miscellaneous varables
real (kind = 8) ssr, conv_crit



!PROGRAM BODY

!	Initialize structural parameters
!	ltot kappa d cmin beta delta eta sigma cgamma lgamma xi nu
structparams = (/5.0d3, 9.0d2, 1.0d5, 2.6d3, 0.95d0, 0.2d0,  1.8d0, 1.75d0, 0.5d0, 0.5d0, 2.0d0, 2.3d0/)
!structparams = (/3.8d0,870.07d0, 98691.4d0, 2584.5d0, 0.98d0, 0.2d0, 1.8d0, 1.78d0,  0.31d0,  0.39d0,1.9d0, 1.5d0/)
!read (*,*)
!parameter values:   4329.5102084019554        870.07494196261041        98691.438175558287        2584.4577746812288       0.98089214831123517       0.19956086722606781        1.7998490549485420        1.7784006858736456       0.30618202318229315       0.38606688813185541        1.3594207172878821        1.9883658329442573

ltot 	= structparams(1)
kappa 	= structparams(2)
d 		= structparams(3)
cmin 	= structparams(4)
beta 	= structparams(5)
delta	= structparams(6)
eta 	= structparams(7)
sigma 	= structparams(8)
cgamma 	= structparams(9)
lgamma 	= structparams(10)
xi 		= structparams(11)
nu 		= structparams(12)

! Pick a type between 1 and 8, and calculate optimal labor choice matrix and value function for a chosen type (choice have been made in the beginning)
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

read (*,*) itype
print *, "Chosen type:", itype
if (itype < 5) then
	ind_subtype = 	itype
else
	ind_subtype	=	itype-4
endif

icoh = (itype-1)/4 + 1 !cohort is 1/2 and not 0/1 because of indexing style; it is passed to valfun, where it is passed to aime_calc and benefit_calc
!sex = mod(ind_type,2)


!Load inputs

call inputs(1) ! (1) means initial inputs
!-----------------------------------------------------------------------
!		TESTING GROUNDS

!----------------------
!once the type is chosen, load data profiles for comparison: lfp, hours, assets and SS application timing
!	1.allocate arrays
!	data
allocate(datalfp(tmom(icoh)))
allocate(datahours(tmom(icoh)))
allocate(datassets(tmom(icoh)))
allocate(dataret(tmom(icoh)))
allocate(sdlfp(tmom(icoh)))
allocate(sdhours(tmom(icoh)))
allocate(sdassets(tmom(icoh)))
allocate(sdret(tmom(icoh)))
allocate(moments(4*tmom(icoh),1))
!	model
allocate(modellfp(tmom(icoh)))
allocate(modelhours(tmom(icoh)))
allocate(modelassets(tmom(icoh)))
allocate(modelret(tmom(icoh)))
!	2.load data
write(datafile,'(a,i1.1,a)') "data_inputs/lfp_data_",itype,".csv"
call import_vector(datafile,datalfp)
write(datafile,'(a,i1.1,a)') "data_inputs/sd_lfp_data_",itype,".csv"
call import_vector(datafile,sdlfp)
write(datafile,'(a,i1.1,a)') "data_inputs/hours_data_",itype,".csv"
call import_vector(datafile,datahours)
write(datafile,'(a,i1.1,a)') "data_inputs/sd_hours_data_",itype,".csv"
call import_vector(datafile,sdhours)
write(datafile,'(a,i1.1,a)') "data_inputs/assets_data_",itype,".csv"
call import_vector(datafile,datassets)
write(datafile,'(a,i1.1,a)') "data_inputs/sd_assets_data_",itype,".csv"
call import_vector(datafile,sdassets)
write(datafile,'(a,i1.1,a)') "data_inputs/applied_",itype,".csv"
call import_vector(datafile,dataret)

sdret = sd(dataret)
!----------------------




!!Calculate stationary distribution of states at age 50 for that type
!call stationary_dist(transmat(:,:,1,itype),statvec)
!call cumsum(statvec,cumstvec)

!MAIN SIMULATION CYCLE
!first, initialize a randomizer
call random_seed(seed_size)
allocate(seed(seed_size))
seed = 100 !fix seed to have reproducible results
call random_seed(put = seed) !start new random sequence
deallocate(seed)
!seed truncated normal
tn_seed = 100

!---!---BIG CALIBRATION CYCLE---!---!
!	Sequence of actions:
!	1. Correct for wage selection bias correction (see full procedure in French 2005): the difference between mean model wages of workers vs all. Use guess on structural parameters (from literature).
!	We need to recover true (hidden) wage here.

!	!	Citation from French:
!	If, for example, the fixed-effects wage profiles overstate average wages at age 60 by 10% in
!	the simulated sample, then it is likely that wages have been overestimated at age 60 by 10%
!	in the PSID data. Therefore, the candidate for the unobserved average wage at age 60 is the
!	fixed-effects estimate from the PSID data, less 10%. This new candidate wage profile is fed
!	into the model and the procedure is repeated. If, for example, the fixed-effects profile using
!	simulated data still indicates a 1% upward bias, the candidate true wage profile is reduced
!	by an additional 1%. This iterative process is continued until a fixed point is found.

!	2.	Feed the updated profiles to the model and calculate set of parameters that minimizes moment conditons.
!	3. 	With the updated new set of structural parameters I calculate wage correction for the profile I calculated previously in (1.), observe the difference.
!	4.	Calclate new set of structural parameters.
!	5.  Repeat 1-4 until wage profiles AND structural parameters change very little with every consecutive iteration.



!initialize convergence criterion
conv_crit = 10
!initialize vector of percentage differences in wages
delta_wage = 10
!initialize initial wage profile (which is currently corresponds to the biased wage profile in the data
initial_wage = exp(logwage(:,icoh))
!initialize iteration counter
l = 0
do while (conv_crit>10*tol)
	l = l+1
	if (l == 50) then
		read (*,*)
	endif
	print *, 'iteration no ', l
	!	1. Correct for wage selection bias correction
    !	1.1. Solve for wage update
    call moment_wage_calc(structparams, ssr, delta_wage)
    
    !Use simulated data to update grid on wage and transition probabilities
    write(command,'(a,i1.1,a)') "cd data_output/stata && stata simulwage_update.do ",itype," && cd ../.."
    
    call execute_command_line (command)  
    call inputs(2)
	
!	initial_wage = updated_wage
!	conv_crit = sum(sqrt(delta_wage**2))
!enddo
!	2.	Feed the updated profiles to the model and calculate set of parameters that minimizes moment conditons.
!I use updated wage grid it to calibrate parameters: so far, to fit lfp, labor supply and assets. Amoeba procedure.
!I put no weight on assets moment for now, since the differences are huge.


	!Construct initial simplex to feed into 'amoeba' procedure; algorithm is borrowed from Matlab "fminsearch"
	print *, 'Total number of simplex vectors:', nparams_cd+1
	spsimplex(1,:) = structparams
	funevals(1) = mom_calc(structparams)
!	!$OMP PARALLEL DO
	do n = 1,nparams_cd
		print *, 'Simplex vectot number', n+1
		y = structparams
		if (y(n) /= 0.) then
			y(n) = (1+usual_delta)*y(n)
		else
			y(n) = zero_term_delta
		endif
		spsimplex(n+1,:) = y
		funevals(n+1) = mom_calc(y)
	enddo
!	!$OMP END PARALLEL DO
	!Solving for minimal ssr
	call amoeba(spsimplex,funevals,tol,mom_calc,iter,info)
	!update structural parameters: choose one of the simplex values
	structparams = spsimplex(1,:)

	initial_wage = updated_wage
	
	conv_crit = sum(sqrt((delta_wage)**2))
	print *, 'convergence criterion', conv_crit
enddo

!		END TESTING
!-----------------------------------------------------------------------

!dellocation
deallocate(datalfp)
deallocate(datahours)
deallocate(datassets)
deallocate(dataret)
deallocate(sdlfp)
deallocate(sdhours)
deallocate(sdassets)
deallocate(sdret)
deallocate(moments)
deallocate(modellfp)
deallocate(modelhours)
deallocate(modelassets)
deallocate(modelret)


contains
	subroutine inputs(mode)
		integer, intent(in) :: mode
		character(len = 10) expr
		if (mode == 1) then
			expr = ''
		else
			expr = '_simul'
		endif
			!-----------------------------------------------------------------------
		!1. Import files

		!1.1 Import transition matrices
		!1.1.1 Wage transitions
		datafile = 'data_inputs/transmat_wage'//trim(expr)//'.csv'
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
						!print *, typemat(ind_type,:)
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

		!1.1.4. Initial joint distribution of young cohort at age 50 of AIME, assets, wage and medical expenses; joint probabilities of bad health and zero assets
		!		Every imported matrix has four rows, which denote:
		!					1 row: male 0, college 0
		!					2 row: male 0, college 1
		!					3 row: male 1, college 0
		!					4 row: male 1, college 1

		!1.1.4.1. 	Joint probabilities of having good/bad health and zero/nonzero assets
		!			4 columns, denote following:
		!			1: Pr(assets = 0, health = 0)
		!			2: Pr(assets = 0, health = 1)
		!			3: Pr(assets = 1, health = 0)
		!			4: Pr(assets = 1, health = 1)
		!	Rows are described above.

		datafile = 'data_inputs/joint_assets_health.csv'
		call import_data(datafile,prob_ah)
		!			Calculate CDF and store in cdf_ah
		!			We'll need it during simulation cycle
		do i = 1,4
			call cumsum(prob_ah(i,:), cdf_ah(i,:))
		enddo
		!1.1.4.2. Means and varcov matrix of joint lognormal distribution of mexp, aime and wage, assets == 0
		!			9 columns, denote following:
		!mu(logaime) mu(logwage) mu(logmexp) Var(logaime) Cov(logaime,logwage)	Var(logwage) Cov(logaime,logmexp) Cov(logwage,logmexp)	Var(logmexp)
		!	Rows are described above.

		datafile = 'data_inputs/joint_amw_zeroas.csv'
		call import_data(datafile,joint_awma0)

		!1.1.4.3. Means and varcov matrix of joint lognormal distribution of mexp, aime and wage and assets when assets > 0
		!			14 columns, denote following:
		!mu(laime) mu(lwage) mu(lmexp) mu(las) Var(laime) Cov(laime,lwage) Var(lwage) Cov(laime,lmexp) Cov(lwage,lmexp)	Var(lmexp) Cov(las,laime) Cov(las,lwage) Cov(las,lmexp) Var(las)
		!	Rows are described above.

		datafile = 'data_inputs/joint_aamw_nzeroas.csv'
		call import_data(datafile,joint_awma)

		!2.	Construct wage and medical expenses grids
		!Import a wage profile, profile of sd and mean of parenting normal distribution (which we are going to use to draw a wage)
		!		And med exp profile, profile of sd and mean of parenting normal distribution (which we are going to use to draw a med exp), additionally conditional on health

		!1. Matrices of bin borders
		datafile = 'data_inputs/wagebin_borders_smooth'//trim(expr)//'.csv'
		call import_data(datafile,wbinborders)

		WRITE(datafile,'(a,i1,a)') "data_inputs/medbin_borders_smooth.csv"
		call import_data(datafile,mbinborders)

		do i = 1,ntypes
			!import from a proper file
			if (mode == 2 .AND. i == itype) then !import from updated file
				datafile = 'data_inputs/lwage_smooth_simul.csv'
			else !import from initial file
				WRITE(datafile,'(a,i1,a)') "data_inputs/lwage_smooth_",i,".csv"
			endif
			call import_data(datafile,logwagedata(:,:,i))
			logwage(:,i) = logwagedata(:,1,i)

			wzmean(:,i) = logwagedata(:,2,i)
			wzsd(:,i) = logwagedata(:,3,i)
			!construct a wage grid:
			!values of wage on the borders of bins
			do k = 1,nwagebins
				do t = 1,lifespan
					call truncated_normal_ab_mean (wzmean(t,i), wzsd(t,i), wbinborders(t,k), wbinborders(t,k+1), bin_mean)
					wagegrid(t,k,i) = exp(logwage(t,i)	+	bin_mean)
				enddo
			enddo


		!Same for medical expenses, additionally condition on health
			do j = 1,nhealth
				WRITE(datafile,'(a,i1,i1,a)') "data_inputs/lmexp_smooth_",i,j-1,".csv"
				call import_data(datafile,logmexpdata(:,:,j,i))
				logmexp(:,j,i) = logmexpdata(:,1,j,i)
				mzmean(:,j,i) = logmexpdata(:,2,j,i)
				mzsd(:,j,i) = logmexpdata(:,3,j,i)
				do k = 1,nmedbins
					do t = 1,lifespan
						call truncated_normal_ab_mean (mzmean(t,j,i), mzsd(t,j,i), mbinborders(t,k), mbinborders(t,k+1), bin_mean)
						mexpgrid(t,j,k,i) = exp(logmexp(t,j,i)	+	bin_mean)
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
					!print *, ms(k,:)
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
	end subroutine inputs
	
	subroutine moment_wage_calc(strpar, momnorm, dwage)
	!	Function takes a vector of structural parameters ("strpar") and spits out norm of moment conditions as a function result.
		use parameters
		use procedures
		use nrtype
		implicit none
		real (DP), dimension(:), intent(inout)	:: strpar ! ltot kappa d cmin age0 beta delta eta sigma ugamma xi nu
		real (kind = 8), dimension(lifespan), intent(out) :: dwage
		real (DP), intent(out) :: momnorm !norm of the moments
!		real (DP), dimension(:), allocatable :: stp

!		Local var
		real (kind = 8), dimension(:,:), allocatable :: multmat, weightmat
		real (kind = 8) mn(1,1), variances(4)
		integer ctype, ccoh


	!	BODY OF FUNCTION
	!	1. fill in the parameters, controlling for possible range
!		allocate(stp(size(strpar)))
		strpar		= max(strpar,0.0)
		strpar(1)	= min(strpar(1), 5.0d3) !ltot E [0,5000]	
		strpar(5) 	= min(max(strpar(5),0.9d0),0.99999d0) !beta E [0.9;0.99999]
		
		ltot 	= strpar(1)
		kappa 	= strpar(2)
		d 		= strpar(3)
		cmin 	= strpar(4)
		beta 	= strpar(5)
		delta 	= strpar(6)
		eta 	= strpar(7)
		sigma 	= strpar(8)
		cgamma 	= strpar(9)
		lgamma 	= strpar(10)
		xi 		= strpar(11)
		nu 		= strpar(12)

	!	2. Calculate value function given these parameters
		!labor choice
		lchoice = -10
		lchoice_na = -20
		!value function
		ctype = 2
		ccoh = 1
		
		call valfun(ccoh,wagegrid(:,:,itype),mexpgrid(:,:,:,itype),transmat(:,:,:,itype),surv(:,:,itype),ms, &
		vf,vf_na,pf,pf_na,lchoice,lchoice_na,app_policy)

		!open files to save results
		write(datafile,'(a,i1.1,a)') "data_output/life_assets_cd_",itype,".txt"
		open(unit=11,file=datafile,status='replace')
		write(datafile,'(a,i1.1,a)') "data_output/life_labor_cd_",itype,".txt"
		open(unit=12,file=datafile,status='replace')
		write(datafile,'(a,i1.1,a)') "data_output/life_earnings_cd_",itype,".txt"
		open(unit=13,file=datafile,status='replace')
		write(datafile,'(a,i1.1,a)') "data_output/life_wage_cd_",itype,".txt"
		open(unit=14,file=datafile,status='replace')
		write(datafile,'(a,i1.1,a)') "data_output/life_medexp_cd_",itype,".txt"
		open(unit=15,file=datafile,status='replace')
		write(datafile,'(a,i1.1,a)') "data_output/app_age_cd_",itype,".txt"
		open(unit=16,file=datafile,status='replace')
		write(datafile,'(a,i1.1,a)') "life_health_cd_",itype,".txt"
		open(unit=17,file=datafile,status='replace')

		counter = 0
		first = .TRUE.

		! initialize means

		mean_lifelabor 		= 0
		mean_lifeearnings 	= 0
		mean_lifewage 		= 0
		mean_lifemedexp 	= 0
		mean_lfp 			= 0
		mean_applied		= 0
		nalive 				= 0
		nwork				= 0
		worker_wage			= 0
		mean_workerwage 	= 0

		do i = 1,nsim !simulate "nsim" individuals of given type
			!First, make a draw to define jointly health status AND presence of assets
			!Once the state is clear, make proper random draw from corresponding distribution (multivariate lognormal)

			call random_number(rndnum)
			if (rndnum <= cdf_ah(ind_subtype,2)) then 	!assets = 0, trivariate lognormal draw
		!		counter = counter + 1
		!		if (counter == 1) then
		!			first = .TRUE.
		!		else
		!			first = .FALSE.
		!		endif
				call random_mvnorm(3, 	real(joint_awma0(ind_subtype,1:3),4), &
										real(joint_awma0(ind_subtype,4:9),4), fchol3, first, draw3, ier)
				in_aime 	= real(exp(draw3(1)),8)
				in_wage 	= real(exp(draw3(2)),8)
				in_mexp 	= real(exp(draw3(3)),8)
				in_asset	= 0.0d0
				if (rndnum <= cdf_ah(ind_subtype,1)) then
					in_health = 1 !bad health
				else
					in_health = 2 !good health
				endif
			else										!assets>0, quadrivariate lognormal draw
		!		counter = counter + 1
		!		if (counter == 1) then
		!			first = .TRUE.
		!		else
		!			first = .FALSE.
		!		endif
				call random_mvnorm(4, 	real(joint_awma(ind_subtype,1:4),4), &
										real(joint_awma(ind_subtype,5:14),4), fchol4, first, draw4, ier)
				in_aime 	= real(exp(draw4(1)),8)
				in_wage 	= real(exp(draw4(2)),8)
				in_mexp 	= real(exp(draw4(3)),8)
				in_asset	= real(exp(draw4(4)),8)
				if (rndnum <= cdf_ah(ind_subtype,3)) then
					in_health = 1 !bad health
				else
					in_health = 2 !good health
				endif
			endif

		!	print *, i
		!	if ( i == 77) then
		!		print *, 'pause'
		!	endif
			
			call lc_simulation(ccoh,in_asset,in_aime,in_wage,in_mexp,in_health, init_workyears, wagegrid(:,:,itype), &
						logwage(:,itype),wzmean(:,itype),wzsd(:,itype), wbinborders, &
						mexpgrid(:,:,:,itype),logmexp(:,:,itype),mzmean(:,:,itype),mzsd(:,:,itype), mbinborders, &
						ms, transmat(:,:,:,itype),vf,vf_na,pf,pf_na,lchoice,lchoice_na,app_policy, surv(:,:,itype), &
						life_assets,life_labor,life_earnings, life_wage, life_medexp,life_health,app_age,tn_seed)
			!construct a vector of number of ss applications at a given age
			if (app_age > 0) then
				mean_applied(app_age-49) = mean_applied(app_age-49) + 1
			endif
			!construct individual lfp given life_labor: 0 means not participating, 1 is participating, -1 means dead
			where (life_labor>0)
				life_lfp = 1
				mean_lfp = mean_lfp + life_lfp
				nwork = nwork + 1
				worker_wage = life_wage
				mean_workerwage = worker_wage + mean_workerwage
			end where
			!	find those who are alive in each age
			! 	and sum up statistics for those alive
			where (life_labor > -1)
				nalive = nalive + 1
				mean_lifeearnings = mean_lifeearnings + life_earnings
				mean_lifelabor = mean_lifelabor + life_labor
				mean_lifemedexp = mean_lifemedexp + life_medexp
				mean_lifewage = life_wage + mean_lifewage
			end where

			where (life_assets > -1) mean_lifeassets = mean_lifeassets+life_assets



			!save the simulation results to .csv files
			do j = 11,17
				call csv_write(j,i,.false.) !put the simulation number in the begining of the line, not advancing to next line
			enddo

			!write vectors and advance to next line
			call csv_write(11,life_assets,.true.)
			call csv_write(12,life_labor,.true.)
			call csv_write(13,life_earnings,.true.)
			call csv_write(14,life_wage,.true.)
			call csv_write(15,life_medexp,.true.)
			call csv_write(16,app_age,.true.)
			call csv_write(17,life_health,.true.)
		enddo

		!close opened files
		do i = 11,17
			close(i)
		enddo

		!calculate means
		mean_lfp = mean_lfp/nalive
		mean_lifelabor = mean_lifelabor/nalive
		mean_lifeassets = mean_lifeassets/nalive
		mean_applied = mean_applied/nalive
		where (nwork > 0)
			mean_lifewage = mean_lifewage/nalive
			mean_workerwage = mean_workerwage/nwork
			!wage bias correction in percentage points
			dwage = (mean_workerwage - mean_lifewage)/mean_workerwage
		elsewhere
			mean_lifewage = 1.0d-6 !black box, we know nothing about the wage of those people apart from the fact it is below threshold value
			mean_workerwage = mean_lifewage !just for consistency
			!no wage bias in this case... ??
			dwage = 0
		end where

		
		
		
		
		!calculate difference between data and model moments
		if (icoh == 1) then
			moments(:,1) = (/datalfp - mean_lfp(8:lifespan),datahours-mean_lifelabor(8:lifespan), &
			datassets - mean_lifeassets(8:lifespan), dataret - mean_applied(8:lifespan)/)
!			variances(1) = variance((/datalfp - mean_lfp(8:lifespan)/))
!			variances(2) = variance((/datahours-mean_lifelabor(8:lifespan)/))
!			variances(3) = variance((/datassets - mean_lifeassets(8:lifespan)/))
!			variances(4) = variance((/dataret - mean_applied(8:lifespan)/))
		else
			moments(:,1) = (/datalfp - mean_lfp(1:18),datahours-mean_lifelabor(1:18), &
			datassets - mean_lifeassets(1:18),dataret - mean_applied(1:18)/)
!			variances(1) = variance((/datalfp - mean_lfp(1:18)/))
!			variances(2) = variance((/datahours-mean_lifelabor(1:18)/))
!			variances(3) = variance((/datassets - mean_lifeassets(1:18)/))
!			variances(4) = variance((/dataret - mean_applied(1:18)/))
		endif
		
		!allocate moments' weight
		!allocate(multmat(3*tmom(icoh),3*tmom(icoh)))
		allocate(weightmat(4*tmom(icoh),4*tmom(icoh)))		
		!multmat = matmul(moments,transpose(moments))
		!call inverse(multmat,weightmat,3*tmom(icoh))
		weightmat = 0.0d0
		do i = 1,tmom(icoh)
			weightmat(i,i) = 1.0d0/sdlfp(i)
			weightmat(tmom(icoh)+i,tmom(icoh)+ i) = 1.0d0/sdhours(i)
			weightmat(2*tmom(icoh)+i,2*tmom(icoh)+i) = 0.0d0 !1.0d0/sdassets(i)
			weightmat(3*tmom(icoh)+i,3*tmom(icoh)+i) = 0.0d0 !1.0d0/sdret(i)
		enddo
!		weightmat(4*tmom(icoh),4*tmom(icoh)) = 1.0d0/variances(4)
!		do i = 1,3*tmom(icoh)
!			if (i<=2*tmom(icoh)) then
!				weightmat(i,i) = 1.0d0
!			else
!				weightmat(i,i) = 1.0d-4
!			endif
!		enddo
		
		mn = matmul(matmul(transpose(moments),weightmat),moments)
		momnorm = mn(1,1)
		!deallocate(multmat)
		deallocate(weightmat)
	end subroutine moment_wage_calc

	function mom_calc(strpar)
		use parameters
		use procedures
		use nrtype
		implicit none
		real (DP), dimension(:), intent(inout)	:: strpar ! ltot kappa d cmin age0 beta delta eta sigma ugamma xi nu
		real (DP) :: mom_calc !norm of the moments
		call moment_wage_calc(strpar, mom_calc, delta_wage)
	end function mom_calc

end !main program

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
    real (kind = 8) vf_na_interp(statesize), vf_interp(statesize) !containers for interpolated value function
    real (kind = 8) val_max_na(grid_asset) !technical variable, maxed vf_na BEFORE comparison vith vf
    real (kind = 8) withheld_b, bnext_erntest

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

!    if (size(wage,1) .NE. lifespan) then
!        print *, 'Wage matrix "lifespan" size is incorrect!'
!        vf = -1
!        return
!    elseif (size(wage,2) .NE. nwagebins) then
!        print *, 'Wage matrix "wage bins" size is incorrect!'
!        vf = -2
!        return
!    endif

!    if (size(medexp,1) .NE. lifespan) then
!        print *, 'Med matrix "lifespan" size is incorrect!'
!        vf = -3
!        return
!    elseif (size(medexp,2) .NE. nhealth) then
!        print *, 'Med matrix "health state" size is incorrect!'
!        vf = -4
!        return
!    elseif (size(medexp,3) .NE. nmedbins) then
!        print *, 'Med matrix "med bins" size is incorrect!'
!        vf = -5
!        return
!    endif

!    if (size(transmat,1) .NE. statesize) then
!        print *, 'Transition matrix "current state" size is incorrect!'
!        vf = -6
!        return
!    elseif (size(transmat,2) .NE. statesize) then
!        print *, 'Transition matrix "next-period state" size is incorrect!'
!        vf = -7
!        return
!    elseif (size(transmat,3) .NE. lifespan-1) then
!        print *, 'Transition matrix "lifespan-1" size is incorrect!'
!        vf = -8
!        return
!    endif

!    if (size(vf,1) .NE. grid_asset) then
!        print *, 'Value function "asset" size is incorrect!'
!        vf = -9
!        return
!    elseif (size(vf,2) .NE. statesize) then
!        print *, 'Value function "state" size is incorrect!'
!        vf = -10
!        return
!    elseif (size(vf,3) .NE. grid_ss) then
!        print *, 'Value function "social security" size is incorrect!'
!        vf = -11
!        return
!    elseif (size(vf,4) .NE. lifespan+1) then
!        print *, 'Value function "lifespan+1" size is incorrect!'
!        vf = -12
!        return
!    endif

!    if (size(surv,1) .NE. lifespan) then
!        print *, 'Survival matrix "lifespan" size is incorrect!'
!        vf = -13
!        return
!    elseif (size(surv,2) .NE. 2) then
!        print *, 'Survival matrix "health state" size is incorrect!'
!        vf = -14
!        return
!    endif

!    if (size(pf,1) .NE. grid_asset) then
!        print *, 'Policy function "asset" size is incorrect!'
!        vf = -15
!        return
!    elseif (size(pf,2) .NE. statesize) then
!        print *, 'Policy function "state" size is incorrect!'
!        vf = -16
!        return
!    elseif (size(pf,3) .NE. grid_ss) then
!        print *, 'Policy function "social security" size is incorrect!'
!        vf = -17
!    elseif (size(pf,4) .NE. lifespan) then
!        print *, 'Policy function "lifespan" size is incorrect!'
!        vf = -18
!        return
!    endif

!    if (size(lchoice,1) .NE. grid_asset) then
!        print *, 'Labchoice "asset" size is incorrect!'
!        vf = -19
!        return
!    elseif (size(lchoice,2) .NE. statesize) then
!        print *, 'Labchoice "state" size is incorrect!'
!        vf = -20
!        return
!    elseif (size(lchoice,3) .NE. grid_ss) then
!        print *, 'Labchoice "social security" size is incorrect!'
!        vf = -21
!        return
!    elseif (size(lchoice,4) .NE. lifespan) then
!        print *, 'Labchoice "lifespan" size is incorrect!'
!        vf = -22
!        return
!    endif

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
!	                    if (t == 31 .AND. k == 2 .AND. s == 41 .AND. j == 1) then
!							print *, 'stop'
!						endif
						call labchoice(opt_labor, cons, util, assets(i),assets(j),h,wage(t,w),medexp(t,h,m),ss(k),age,cohort)
	                    lab_cur_next(i,j) = opt_labor
	                    !print *, opt_labor

	                    !**** EARNINGS TEST PT.2
	                    !First part is inside "labchoice" routine
						!calculate the updated benefit level for the next period
	                    call earnings_test(cohort,age,ss(k), opt_labor*wage(t,w), withheld_b,bnext_erntest)

	                    if (k == grid_ss) then !highest possible benefits
							bnext_erntest = ss(grid_ss) !in the unlikely event that individual benefits are already on the highest level, one can't get more
						endif
	                    !**** EARNINGS TEST PT.2

	                    if (bnext_erntest > ss(k+1) .AND. k /= grid_ss) then
							print *, 'bnext_erntest is too large! Press any button to proceed:'
							read (*,*)
						endif
	                    !if updated benefit differs from current one
	                    if (bnext_erntest > ss(k)) then
	                    !we know it can't be above ss(k+1)
							!Calculate set of interpolations for next-period vf given next-period: one interpolation for every state s
							do st = 1,statesize
								vf_interp(st) = linproj(bnext_erntest,ss(k),ss(k+1),vf(j,st,k,t+1),vf(j,st,k+1,t+1)) !vf(j - next period asset, st - next-period state;k,k+1.. - indices of ss, t+1 - next period)
							enddo
						elseif (bnext_erntest < ss(k)) then
							print *, 'Next period benefit can''t decrease here! An error! Press any button to proceed:'
							read (*,*)
		                else !no need to interpolate
							vf_interp = vf(j,:,k,t+1)
						endif
	                    !unmaxed "value function"
	                    val_nomax(i,j) = util + beta*(surv(t,h)*dot_product(vf_interp,transmat(s,:,t)) &		!notice vf_interp here, it plays a role for earnings test
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
!	            if (sum(lchoice(:,s,k,t))>0) then
!	                     print *, lchoice(:,s,k,t)
!	            endif
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

	vf_na = -103.0d0
    pf_na = -104.0d0

	vf_na(:,:,:,lifespan:lifespan+1,1) = vf(:,:,:,lifespan:lifespan+1) !upon death and in the last period, the two value functions are identical (everyone applies at age 90 with certainty)
	vf_na(:,:,:,lifespan:lifespan+1,2) = vf(:,:,:,lifespan:lifespan+1)
	app_policy = 0 !initialize application policy
	app_policy(:,:,:,lifespan,:) = 1 ! We claim that in the last period of life everybody applies for ss with certainty

	do bcalc = 1,2
		!bcalc = 1: we calculate "next period" benefits as if person has less than 35 working years
		!bcalc = 2: calculate bnext as if an individual has more than 35 working years
		do t = lifespan-1,1,-1 !age 89 to 50
			age = t+49 !real age
			if (age == 58) then
				print *, 58
			endif
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
							call labchoice(opt_labor, cons, util, assets(i),assets(j),h,wage(t,w),medexp(t,h,m),0.0d0,age,cohort)
		                    lab_cur_next_na(i,j) = opt_labor

		                    !NOW, HERE ARE THE BIG DIFFERENCES BETWEEN THIS VF AND THE FIRST ONE
		                    !First, we calculate next-period B (it's going to be off the grid), and use linear interpolation to calculate
		                    !value functions between points of B. Then we have to take expectation of these interpolations.
		                    !On top of that, we need to control for application age: an individual can't apply if she's younger than 62
!		                    if (k == grid_ss) then
!								print *, 'checkpoint'
!		                    endif

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

subroutine lc_simulation(cohort,in_asset,in_aime,in_wage,in_mexp,in_health,years_aime, &
						wage,logwage,wzmean,wzsd,wbinborders,mexp,logmexp,mzmean,mzsd, &
						mbinborders,ms,transmat,vf,vf_na,pf,pf_na,lchoice,lchoice_na,app_policy, surv, &
						life_assets,life_labor,life_earnings, life_wage, life_medexp,life_health,app_age,tn_seed)
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
integer, intent(in)	:: cohort, in_health
real (kind = 8) 	:: in_asset, in_aime,in_wage,in_mexp
integer, intent(in) :: years_aime
real (kind = 8), dimension(lifespan,nwagebins), intent(in) :: wage
real (kind = 8), dimension(lifespan), intent(in) :: wzmean, wzsd,logwage
real (kind = 8), intent(in) :: wbinborders(lifespan,nwagebins+1), mbinborders(lifespan,nmedbins+1)
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
real (kind = 8), dimension(lifespan,nhealth), intent(in) :: surv
real (kind = 8), dimension(grid_asset,statesize,grid_ss,lifespan,2), intent(in) :: lchoice_na
real (kind = 8), intent(out) :: life_assets(lifespan)
real (kind = 8), dimension(lifespan), intent(out) :: life_labor, life_earnings, life_wage, life_medexp, life_health
integer, intent(out) :: app_age !the variable wil contain age of social security application
integer, intent(inout) :: tn_seed


!LOCAL VARIABLES
real (kind = 8) anext, workyears, aime, apx, lpx, acur, wcur, mcur, bcur, vx1,vx2
real (kind = 8) wage_tn, med_tn !wage and medical expenses residiual drawn from truncated normal
real (kind = 8) bin_mean !technical variable to store mean of given bin of residuals
integer iaprev, ianext, ibprev, ibnext
integer iwprev,iwnext,imprev,imnext !indices of previous and next gridpoints in (w,m) - space

!values of bin borders IN NUMBRS, NOT IN LOG DEVIATIONS (NOTE THE DIFFERENCE FROM WBINBORDERS);
!FIRST BIN NOT LIMITED FROM BELOW, UPPER BIN IS NOT LIMITED FRO ABOVE
!ONLY PURPOSE OF THIS IS TO DETERMINE INITIAL STATE
real (kind = 8) init_wage_bins(nwagebins-1), init_mexp_bins(nmedbins-1)


real (kind = 8) rndnum, statvec(statesize), cumstvec(statesize), cumtrans(statesize)
real (kind = 8) assets(grid_asset), ss(grid_ss)
real (kind = 8) bound !variable used to calculate truncation bound for truncated normal
integer t, curstate, nextstate , bcalc, age, i
integer w,m,h !wage, med, health index
integer death_age


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
real (kind = 8) :: lmin = 5.0d3, lnmin = 5.0d3
!Corresponding state values in s-space
integer state00, state01,state10,state11



!initialize vectors
life_earnings 	= -1.0d0
life_assets 	= -1.0d0
life_labor		= -1.0d0
life_medexp 	= -1.0d0
life_wage		= -1.0d0
life_health		= -1.0d0

!   create assets grid and social security grid
assets(1)	= 0
ss(1)		= 0
call logspace(asset_min,asset_max,grid_asset-1,assets(2:grid_asset))
call logspace(ss_min,ss_max,grid_ss-1,ss(2:grid_ss))


h = in_health
!ensure that inputs are within the borders
in_wage 	= min(max(in_wage,wage(1,1)),wage(1,5))
in_mexp 	= min(max(in_mexp,mexp(1,h,1)),mexp(1,h,5))
in_asset 	= min(in_asset,asset_max)
in_aime		= min(in_aime,aime_max)


!determine current state, given initial health, wage and mexp, given health status



	do i = 1,nwagebins-1
		init_wage_bins(i) = exp(logwage(1)+wbinborders(1,i+1))
	enddo
	do i = 1,nmedbins-1
		init_mexp_bins(i) = exp(logmexp(1,h)+mbinborders(1,i+1))
	enddo

	!determine wage and mexp bin
	w = locate_greater(init_wage_bins,in_wage)
	m = locate_greater(init_mexp_bins,in_mexp)
	!determine last bins
	if 		(w == 0) then
		w = nwagebins
	endif

	if 	(m == 0) then
		m = nmedbins
	endif

	nextstate = h+nhealth*(m-1)+nmedbins*nhealth*(w-1)
!END determine current state,

wcur = in_wage 	!input wage
mcur = in_mexp	!input mexp
acur = in_asset !input asset
aime = in_aime 	!input aime
workyears = years_aime !and years worked
app_age = 1000 !we first set it to very high number, just for initialization


do t = 1,lifespan
	age = t+49
	curstate = nextstate !update state

	!put current values into personal history
	life_assets(t) 	= acur
	life_wage(t) 	= wcur
	life_medexp(t) 	= mcur
	life_health(t) 	= h

	! Find nearest wage grid points
	if (life_wage(t)>=wage(t,w) .AND. w /= 5) then
		iwprev = w
		iwnext = w + 1
	else
		iwprev = w - 1
		iwnext = w
	endif

	! Nearest med grid points
	if (life_medexp(t)>=mexp(t,h,m) .AND. m /= 5) then
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
	bcur = min(bcur,ss_max)
	!given current benefit, calculate two closest locations on benefit grid
	ibnext = locate_greater(ss,bcur) !index of closest asset on the grid from above
	if (ibnext == 0) then
		ibnext = grid_ss
	endif

	if (ibnext /= 1) then
		ibprev = ibnext-1 !index of closest asset on the grid from below
	else
		ibprev = 1
		ibnext = 2
	endif
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

	!	Also, among the values of labor supply, I look for the smallest nonzero number. That is, a very cheap way to avoid interpolation error,
	!	which produces small positive annual hours. This smallest nonzero number will serve as a treshold: interpolated value either zero or larger than this 
	!	smallest number.

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
	if (l0000 < lmin .AND. l0000>0) then 
		lmin = l0000
	endif
	if (ln0000 < lnmin .AND. ln0000>0) then 
		lnmin = ln0000
	endif
		!0001
	v0001 	= vf			(iaprev,state00,ibnext,t)
	vn0001 	= vf_na			(iaprev,state00,ibnext,t,bcalc)
	a0001 	= assets(pf		(iaprev,state00,ibnext,t))
	an0001 	= assets(pf_na	(iaprev,state00,ibnext,t,bcalc))
	l0001 	= lchoice		(iaprev,state00,ibnext,t)
	ln0001 	= lchoice_na	(iaprev,state00,ibnext,t,bcalc)
	if (l0001 < lmin .AND. l0001>0) then 
		lmin = l0001
	endif
	if (ln0001 < lnmin .AND. ln0001>0) then 
		lnmin = ln0001
	endif
		!0010
	v0010	= vf			(iaprev,state01,ibprev,t)
	vn0010	= vf_na			(iaprev,state01,ibprev,t,bcalc)
	a0010 	= assets(pf		(iaprev,state01,ibprev,t))
	an0010 	= assets(pf_na	(iaprev,state01,ibprev,t,bcalc))
	l0010 	= lchoice		(iaprev,state01,ibprev,t)
	ln0010 	= lchoice_na	(iaprev,state01,ibprev,t,bcalc)
	if (l0010 < lmin .AND. l0010>0) then 
		lmin = l0010
	endif
	if (ln0010 < lnmin .AND. ln0010>0) then 
		lnmin = ln0010
	endif
		!0011
	v0011	= vf			(iaprev,state01,ibnext,t)
	vn0011	= vf_na			(iaprev,state01,ibnext,t,bcalc)
	a0011 	= assets(pf		(iaprev,state01,ibnext,t))
	an0011 	= assets(pf_na	(iaprev,state01,ibnext,t,bcalc))
	l0011 	= lchoice		(iaprev,state01,ibnext,t)
	ln0011 	= lchoice_na	(iaprev,state01,ibnext,t,bcalc)
	if (l0011 < lmin .AND. l0011>0) then 
		lmin = l0011
	endif
	if (ln0011 < lnmin .AND. ln0011>0) then 
		lnmin = ln0011
	endif
		!0100
	v0100	= vf			(iaprev,state10,ibprev,t)
	vn0100	= vf_na			(iaprev,state10,ibprev,t,bcalc)
	a0100 	= assets(pf		(iaprev,state10,ibprev,t))
	an0100 	= assets(pf_na	(iaprev,state10,ibprev,t,bcalc))
	l0100 	= lchoice		(iaprev,state10,ibprev,t)
	ln0100 	= lchoice_na	(iaprev,state10,ibprev,t,bcalc)
	if (l0100 < lmin .AND. l0100>0) then 
		lmin = l0100
	endif
	if (ln0100 < lnmin .AND. ln0100>0) then 
		lnmin = ln0100
	endif
		!0101
	v0101	= vf			(iaprev,state10,ibnext,t)
	vn0101	= vf_na			(iaprev,state10,ibnext,t,bcalc)
	a0101 	= assets(pf		(iaprev,state10,ibnext,t))
	an0101 	= assets(pf_na	(iaprev,state10,ibnext,t,bcalc))
	l0101 	= lchoice		(iaprev,state10,ibnext,t)
	ln0101 	= lchoice_na	(iaprev,state10,ibnext,t,bcalc)
	if (l0101 < lmin .AND. l0101>0) then 
		lmin = l0101
	endif
	if (ln0101 < lnmin .AND. ln0101>0) then 
		lnmin = ln0101
	endif
		!0110
	v0110	= vf			(iaprev,state11,ibprev,t)
	vn0110	= vf_na			(iaprev,state11,ibprev,t,bcalc)
	a0110 	= assets(pf		(iaprev,state11,ibprev,t))
	an0110 	= assets(pf_na	(iaprev,state11,ibprev,t,bcalc))
	l0110 	= lchoice		(iaprev,state11,ibprev,t)
	ln0110 	= lchoice_na	(iaprev,state11,ibprev,t,bcalc)
	if (l0110 < lmin .AND. l0110>0) then 
		lmin = l0110
	endif
	if (ln0110 < lnmin .AND. ln0110>0) then 
		lnmin = ln0110
	endif
		!0111
	v0111	= vf			(iaprev,state11,ibnext,t)
	vn0111	= vf_na			(iaprev,state11,ibnext,t,bcalc)
	a0111 	= assets(pf		(iaprev,state11,ibnext,t))
	an0111 	= assets(pf_na	(iaprev,state11,ibnext,t,bcalc))
	l0111 	= lchoice		(iaprev,state11,ibnext,t)
	ln0111 	= lchoice_na	(iaprev,state11,ibnext,t,bcalc)
	if (l0111 < lmin .AND. l0111>0) then 
		lmin = l0111
	endif
	if (ln0111 < lnmin .AND. ln0111>0) then 
		lnmin = ln0111
	endif
		!1000
	v1000 	= vf			(ianext,state00,ibprev,t)
	vn1000 	= vf_na			(ianext,state00,ibprev,t,bcalc)
	a1000 	= assets(pf		(ianext,state00,ibprev,t))
	an1000 	= assets(pf_na	(ianext,state00,ibprev,t,bcalc))
	l1000 	= lchoice		(ianext,state00,ibprev,t)
	ln1000 	= lchoice_na	(ianext,state00,ibprev,t,bcalc)
	if (l1000 < lmin .AND. l1000>0) then 
		lmin = l1000
	endif
	if (ln1000 < lnmin .AND. ln1000>0) then 
		lnmin = ln1000
	endif
		!1001
	v1001 	= vf			(ianext,state00,ibnext,t)
	vn1001 	= vf_na			(ianext,state00,ibnext,t,bcalc)
	a1001 	= assets(pf		(ianext,state00,ibnext,t))
	an1001 	= assets(pf_na	(ianext,state00,ibnext,t,bcalc))
	l1001 	= lchoice		(ianext,state00,ibnext,t)
	ln1001 	= lchoice_na	(ianext,state00,ibnext,t,bcalc)
	if (l1001 < lmin .AND. l1001>0) then 
		lmin = l1001
	endif
	if (ln1001 < lnmin .AND. ln1001>0) then 
		lnmin = ln1001
	endif
		!1010
	v1010	= vf			(ianext,state01,ibprev,t)
	vn1010	= vf_na			(ianext,state01,ibprev,t,bcalc)
	a1010 	= assets(pf		(ianext,state01,ibprev,t))
	an1010 	= assets(pf_na	(ianext,state01,ibprev,t,bcalc))
	l1010 	= lchoice		(ianext,state01,ibprev,t)
	ln1010 	= lchoice_na	(ianext,state01,ibprev,t,bcalc)
	if (l1010 < lmin .AND. l1010>0) then 
		lmin = l1010
	endif
	if (ln1010 < lnmin .AND. ln1010>0) then 
		lnmin = ln1010
	endif
		!1011
	v1011	= vf			(ianext,state01,ibnext,t)
	vn1011	= vf_na			(ianext,state01,ibnext,t,bcalc)
	a1011 	= assets(pf		(ianext,state01,ibnext,t))
	an1011 	= assets(pf_na	(ianext,state01,ibnext,t,bcalc))
	l1011 	= lchoice		(ianext,state01,ibnext,t)
	ln1011 	= lchoice_na	(ianext,state01,ibnext,t,bcalc)
	if (l1011 < lmin .AND. l1011>0) then 
		lmin = l1011
	endif
	if (ln1011 < lnmin .AND. ln1011>0) then 
		lnmin = ln1011
	endif
		!1100
	v1100	= vf			(ianext,state10,ibprev,t)
	vn1100	= vf_na			(ianext,state10,ibprev,t,bcalc)
	a1100 	= assets(pf		(ianext,state10,ibprev,t))
	an1100 	= assets(pf_na	(ianext,state10,ibprev,t,bcalc))
	l1100 	= lchoice		(ianext,state10,ibprev,t)
	ln1100 	= lchoice_na	(ianext,state10,ibprev,t,bcalc)
	if (l1100 < lmin .AND. l1100>0) then 
		lmin = l1100
	endif
	if (ln1100 < lnmin .AND. ln1100>0) then 
		lnmin = ln1100
	endif
		!1101
	v1101	= vf			(ianext,state10,ibnext,t)
	vn1101	= vf_na			(ianext,state10,ibnext,t,bcalc)
	a1101 	= assets(pf		(ianext,state10,ibnext,t))
	an1101 	= assets(pf_na	(ianext,state10,ibnext,t,bcalc))
	l1101 	= lchoice		(ianext,state10,ibnext,t)
	ln1101 	= lchoice_na	(ianext,state10,ibnext,t,bcalc)
	if (l1101 < lmin .AND. l1101>0) then 
		lmin = l1101
	endif
	if (ln1101 < lnmin .AND. ln1101>0) then 
		lnmin = ln1101
	endif
		!1110
	v1110	= vf			(ianext,state11,ibprev,t)
	vn1110	= vf_na			(ianext,state11,ibprev,t,bcalc)
	a1110 	= assets(pf		(ianext,state11,ibprev,t))
	an1110 	= assets(pf_na	(ianext,state11,ibprev,t,bcalc))
	l1110 	= lchoice		(ianext,state11,ibprev,t)
	ln1110 	= lchoice_na	(ianext,state11,ibprev,t,bcalc)
	if (l1110 < lmin .AND. l1110>0) then 
		lmin = l1110
	endif
	if (ln1110 < lnmin .AND. ln1110>0) then 
		lnmin = ln1110
	endif
		!1111
	v1111	= vf			(ianext,state11,ibnext,t)
	vn1111	= vf_na			(ianext,state11,ibnext,t,bcalc)
	a1111 	= assets(pf		(ianext,state11,ibnext,t))
	an1111 	= assets(pf_na	(ianext,state11,ibnext,t,bcalc))
	l1111 	= lchoice		(ianext,state11,ibnext,t)
	ln1111 	= lchoice_na	(ianext,state11,ibnext,t,bcalc)
	if (l1111 < lmin .AND. l1111>0) then 
		lmin = l1111
	endif
	if (ln1111 < lnmin .AND. ln1111>0) then 
		lnmin = ln1111
	endif

	!check whether individual is:
	!1) age<62: can't apply; use only vf_na (not vf)
	!2) age>=62 AND age<app_age: not yet applied, but we must check whether applies: compare values of two quadrilinear interpolations and use higher one (change app_age if it is the case)
	!3) age>=62 AND age>age_app: already applied, use only vf (not vf_na)

	!Before calculating interpolated NEXT PERIOD assets and CURRENT labor choice, we check if individual is eligible and didn't apply yet. If it's the case, we
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
		!Here, we only use vf_na to construct NEXT PERIOD assets and CURRENT PERIOD labor choice; quadlin
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
		!Careful here: we'd like to avoid situations where an individual works few hours a year solely because of interpolation.
		!First, let's find out if any of lnXXXX are 0. If none, then use lpx as is. If there is, then compare lpx and the smallest value of ln that is larger than 0 (lnmin).
		!If lpx is even smaller, than set it to 0.
		if (min(ln0000,ln0001,ln0010,ln0011,ln0100,ln0101,ln0110,ln0111, & 
				ln1000,ln1001,ln1010,ln1011,ln1100,ln1101,ln1110,ln1111) == 0 .AND. lpx >0) then
			if (lpx < lnmin) then
				lpx = 0.0d0
			endif
		endif

		!we have to update AIME
		if (bcalc == 1) then !individual worked less than 35 years
			aime = aime + wcur*lpx/35.0 		      !next period aime grows
			aime = min(aime,aime_max)
		else !individual aready has more than 35 working years
			aime = aime + max(0.0,(wcur*lpx-aime)/35.0) !next-period aime doesn't decline
			aime = min(aime,aime_max)
		endif

		!update number of working years (if individual worked)
		if (lpx>0) then
			workyears = workyears+1
		endif
	elseif (age>=62 .AND. age>=app_age) then !individual applied
		!calculate interpolated next-period assets and labor choice using functions of applicants
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
		!Careful here: we'd like to avoid situations where an individual works few hours a year solely because of interpolation.
		!First, let's find out if any of lXXXX are 0. If none, then use lpx as is. If there is, then compare lpx and the smallest value of lnXXXX that is larger than 0 (lmin).
		!If lpx is even smaller, than set it to 0.
		if (min(l0000,l0001,l0010,l0011,l0100,l0101,l0110,l0111, & 
				l1000,l1001,l1010,l1011,l1100,l1101,l1110,l1111) == 0 .AND. lpx>0) then
			if (lpx < lmin) then
				lpx = 0.0d0
			endif
		endif
		!In this case, we don't care about the AIME and working years anymore, it doesn't play any role
	else !Excluded case, just fishing for errors
		print *, "Wrong scenario, check SS applicaton code in lc_simulation"
		stop
	endif

	!define anext
	anext = apx	!notice, that we will put this value into the vector of lifetime assets in the beginning of next cycle
	!put labor into personal working history

	!	LABOR HISTORY UPDATE
	life_labor(t) = lpx
	!Fill in individual history

	!	EARNINGS HISTORY UPDATE
	life_earnings(t) = 	wcur*life_labor(t)

!-----------------------------------------------------------------------
!	Now we know everything for current period, including labor hoice and optimal next period asset choice.

!	We need to update state, using transition matrix
	call cumsum(transmat(curstate,:,t),cumtrans)
	call random_number(rndnum)
	nextstate = locate_greater(cumtrans,rndnum)

	!Check if an individual survives to the next period
	!notice that at age 90 an individual doesn't survive
	call random_number(rndnum)
	if (rndnum > surv(t,h)) then !a person is dead
!		life_assets(t+2:lifespan+1) = -1
!		life_labor(t+1:lifespan) 	= -1
!		life_wage(t+1:lifespan) 	= -1
!		life_earnings(t+1:lifespan) = -1
!		life_medexp(t+1:lifespan) 	= -1
		if (app_age>age) then !individual didn't apply for the benefits
			app_age = -1
		endif

		exit !exit DO loop
	endif

!	Given that we know next period state, we can determine wage and medical expenses for the next period

	!determine wage shock and med shock corresponding to the next-period state,
	!that is, which bin we're at now for wage and for med exp
	w = int(ms(nextstate,1))
	m = int(ms(nextstate,2))
	h = int(ms(nextstate,3))
	!given the bin, make two draws from truncated normal for both wage and med exp next period
	!First, draw wage residual
	if (w /= 1 .AND. w /= nwagebins) then !not first and not last bin
		call truncated_normal_ab_sample(wzmean(t+1),wzsd(t+1),wbinborders(t+1,w), &
										wbinborders(t+1,w+1),tn_seed,wage_tn)
	elseif (w == 1) then !first bin
		call truncated_normal_ab_mean (wzmean(t+1), wzsd(t+1), wbinborders(t+1,1), wbinborders(t+1,2), bin_mean) !calculate mean of the first bin
		call truncated_normal_ab_sample(wzmean(t+1),wzsd(t+1),bin_mean, & 		!draw residal from mean of the first bin to edge of the first bin
										wbinborders(t+1,2),tn_seed,wage_tn)
	elseif (w == nwagebins) then !last bin
		call truncated_normal_ab_mean (wzmean(t+1), wzsd(t+1), wbinborders(t+1,nwagebins), &
		 wbinborders(t+1,nwagebins+1), bin_mean) !calculate mean of the last bin
		call truncated_normal_ab_sample(wzmean(t+1),wzsd(t+1),wbinborders(t+1,nwagebins), &
										bin_mean,tn_seed,wage_tn)
	else
		print *, "Something's wrong with current wage state"
		read (*,*)
		stop
	endif


	! Calculate wage
	wcur = exp(logwage(t+1)+wage_tn)


	! Second draw med exp residuals
	if (m /= 1 .AND. m /= nmedbins) then
		call truncated_normal_ab_sample(mzmean(t+1,h),mzsd(t+1,h),mbinborders(t+1,m), &
										mbinborders(t+1,m+1),tn_seed,med_tn)
	elseif (m == 1) then
		call truncated_normal_ab_mean 	(mzmean(t+1,h), mzsd(t+1,h), mbinborders(t+1,1), mbinborders(t+1,2), bin_mean) !calculate mean of the first bin
		call truncated_normal_ab_sample	(mzmean(t+1,h),	mzsd(t+1,h), bin_mean, &
										mbinborders(t+1,2),tn_seed,med_tn)
	elseif (m == nmedbins) then
		call truncated_normal_ab_mean 	(mzmean(t+1,h), mzsd(t+1,h), mbinborders(t+1,nmedbins), &
		 mbinborders(t+1,nmedbins+1), bin_mean) !calculate mean of the last bin
		call truncated_normal_ab_sample	(mzmean(t+1,h),	mzsd(t+1,h), mbinborders(t+1,nmedbins), &
										bin_mean,tn_seed,med_tn)
	else
		print *, "Something's wrong with current med exp state"
		read (*,*)
		stop
	endif
	! Calculate med exp
	mcur = exp(logmexp(t+1,h)+med_tn)

	!update assets
	acur = anext

enddo

end subroutine lc_simulation




