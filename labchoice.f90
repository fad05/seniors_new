subroutine labchoice (lchoice, cons, util, acur,anext, h, wage, mexp,ss,age,cohort)
	use procedures
    use parameters
    implicit none
    real (kind = 8), intent(in) :: acur, anext, wage,mexp,ss !current assets, next period assets, current wage, mexp, social scurity
    integer, intent(in) :: h, age, cohort !health status, age and cohort
    real (kind = 8), intent(out) :: lchoice !optimal labor choice
    real (kind = 8), intent(out) :: cons !consumption under optimal choice
    real (kind = 8), intent(out) :: util !value of utility under optimal choice

    !SCALED VARIABLES
    real (kind = 8) acur_scaled, anext_scaled, mexp_scaled, ss_scaled, lchoice_scaled
    real (kind = 8) :: ltot_scaled
	real (kind = 8) :: kappa_scaled, xi_scaled
	real (kind = 8) :: cmin_scaled
	!variables necessary to use minpack and local variables
    integer (kind = 4) info, i, u_max_index, u_min_index
    integer (kind = 4) :: counter !counter for search attempts
    real (kind = 8) fvec(1), lc(1), lmax, lmin, labgrid_max, uchoice, ctemp, lguess, utemp
    real (kind = 8) u0, c0 !utility and consumption IF labor is 0, GIVEN assets, mexp, wage
    real (kind = 8) bma ! bma stands for "b medicaid"; computed in a way such that consuption floor cmin is satisfied; in this case anext MUST be 0
    real (kind = 8), dimension(100) :: labgrid, cgrid, ugrid
    real (kind = 8) a,b,c,dnew,ua,ub,uc,ud,cd ! triplet of points for bracketing method

    interface
        subroutine lmdif1 ( fcn, m, n, x, fvec, tol, info )
            implicit none
            external fcn
            integer ( kind = 4 ) m
            integer ( kind = 4 ) n
            real ( kind = 8 ) x(n)
            real ( kind = 8 ) fvec(m)
            real ( kind = 8 ) tol
            integer ( kind = 4 ) info
        end subroutine lmdif1
    end interface

    !BODY OF SUBROUTINE
    
    
    lchoice = -5 ! initialize with negative value

!=========================SCALING=========================================================================================== 
!1. WE NEED TO SCALE THE PROBLEM, UTILITY TO GIVE SENSIBLE VALUES
!	Recall: c = (1+r)*a +w*l + ss - a' - m
!	Currently c,a, m are expressed in US dollars, and l in hours, and w in 1usd/1hr.
!	We would like to express (C, A, A', M, SS) in 1000s of US dollars, (W) in 1000USD/1000hrs and L in 1000s of hours.
!	C = c/1000 (thd USD); A = a/1000 (thd USD); M = m/1000; SS = ss/1000
!	L = l/1000 (1000hrs).
!	NOTE: we don't need to scale wage!
!	W = w (since 1000USD/1000hrs = 1USD/1hr). 
!	Then:
!	c/1000 = C = (1+r)*A + WL - A' - M = (1+r)a/1000 + w*l/1000 + ss/1000 -a'/1000 -m/1000
		
	acur_scaled = acur/scale_factor
	anext_scaled = anext/scale_factor
	mexp_scaled = mexp/scale_factor
	ss_scaled = ss/scale_factor
	
	lmax = (ltot - kappa - 1)/scale_factor !maximum possible labor minus one working hour, SCALED
	lmin = hmin/scale_factor
	
	!Scale parameters of utility functon!
	ltot_scaled = ltot/scale_factor
	kappa_scaled = kappa/scale_factor
	xi_scaled = xi/scale_factor
	cmin_scaled = cmin/scale_factor
	
    !check some points on the grid of labor, looking for good starting point for solver.
    !The grid spans from small positive (nonzero) number to max possible labor (ltot-kappa) minus small number, all scaled
    !The grid has 100 equally spaced nodes.
    call linspace(lmin,lmax,100,labgrid) !GRID is ALSO SCALED, starting with ONE WORKING HOUR
!=========================END SCALING========================================================================================
	
    

    
    do i = 1,100 !on this labor grid, find corresponding utility and consumption
	    call cons_util_calc(labgrid(i), wage, acur_scaled, anext_scaled, mexp_scaled,ss_scaled, cgrid(i), ugrid(i))
    enddo

    
    u_max_index = maxloc(ugrid,1) !index of maximal utility on the grid
    labgrid_max = labgrid(u_max_index)
    lguess = labgrid_max !initial guess on optimal labor choice OFF the grid: maximum ON the grid
    

!	1.Bracketing method, slowest but sure
!	We're looking for maximum 

	if (u_max_index == 1) then
		lchoice_scaled = labgrid(1)
	else
		a = labgrid(u_max_index-1)
		ua = ugrid(u_max_index-1)
		b = labgrid(u_max_index)
		ub = ugrid(u_max_index)
		c = labgrid(u_max_index+1)
		uc = ugrid(u_max_index+1)
	
		counter = 0
		
		do while (c-a>1.0d-3) !the tolerance is ohe hour (1E-3 of 1000s of hours)
			counter = counter + 1
			if (counter > 1000) then
				print *, 'pause'
			endif
			if (b-a<=c-b) then
				dnew = (b+c)/2.0
				call cons_util_calc(dnew, wage, acur_scaled, anext_scaled, mexp_scaled, ss_scaled,cd, ud)  
			else
				dnew = (a+b)/2.0
				call cons_util_calc(dnew, wage, acur_scaled, anext_scaled, mexp_scaled, ss_scaled,cd, ud)  
			endif
			
			if (dnew<b .AND. ud<ub) then
				a = dnew
				ua = ud
			elseif (dnew<b .AND. ud>ub) then
				c = b
				uc = ub
				b = dnew
				ub = ud
			elseif (dnew>b .AND. ud>ub) then
				a = b
				ua = ub
				b = dnew
				ub = ud
			else 
				c = dnew
				uc = ud
			endif
		enddo
		lchoice_scaled = c	
	endif
		
	call cons_util_calc(lchoice_scaled, wage, acur_scaled, anext_scaled, mexp_scaled, ss_scaled,cons, util)
	lchoice = lchoice_scaled*scale_factor !solution, expressed in hours
		
    ! Utility and consumption IF labor is 0:
    call cons_util_calc(0.0d0,wage,acur_scaled,anext_scaled,mexp_scaled,ss_scaled,c0,u0)
    
    !if utility at 0 is at least as good as utility in nonzero
    if (u0 >= util) then
		lchoice = 0.0d0
	endif

return !end of subroutine

    contains
        
        subroutine cons_util_calc(lbr, wg, acrnt, anxt, mxp, ssec, cns, utl)
			real (kind = 8), intent(in) :: lbr, wg, acrnt, anxt, mxp, ssec
			real (kind = 8), intent(out) :: cns, utl
			real (kind = 8) taxes, income, taxable_income, earnings, withheld_b, ssec_received, bnext
		
		
			!BODY
			!1.Calculate consumption and utility
			!two possible cases, depending on next period assets:
			!1. anext = 0; in this case an individual with consumption lower than cmin (due to medical expenses) is entitled to a Medicaid bobus bma to compensate his expenses exactly to cmin
			!2. anext > 0; if here consumption is lower than cmin, the case is not valid: agend sacrifices consumption for savings. If in this case consumption is lover than cmin, agent is not compensated;
			!if consumption is negative, the utility is set to large negative number, in order never to be optimal.
			!Another fork is whether labor is zero or not. If it is, there is no utility penalty, expressed as KAPPA.
			
			
			!**** EARNINGS TEST PT.1
			!I have to check the earnings test
			earnings = lbr*wg !earnings
			call earnings_test(cohort,age,ssec*scale_factor, earnings*scale_factor, withheld_b, bnext) !calculate withheld amount
			withheld_b = withheld_b/scale_factor 	!scale withheld amount
			ssec_received = ssec - withheld_b 		!calculate actually received socsec benefit
			!**** EARNINGS TEST PT.1
			!The further implementation of earnings test is inside "valfun" routine
			
			income = r*acrnt + earnings !asset return + labor income
			taxable_income = taxable_amount(income*scale_factor, ssec_received*scale_factor)
			taxes = income_tax(taxable_income)/scale_factor
			
			cns = (1+r)*acrnt + lbr*wg +ssec_received- anxt - mxp - taxes !uncompensated consumption; in the case when lbr==0, labor income is ZERO
			if (lbr == 0) then !labor is zero
				if (anxt == 0) then !next-period assets are zero: agent is eligible for minimal consulption level					
					if (cns > cmin_scaled) then !consumption is hgher than minimal, no transfers
						utl = ((1+delta*(h-1))*(cns**ugamma*ltot_scaled**(1-ugamma))**(1-sigma))/(1-sigma) !NO KAPPA IN UTILITY: zero labor
					else
						cns = cmin_scaled !agent reseives a granted consumption floor
						utl = ((1+delta*(h-1))*(cns**ugamma*ltot_scaled**(1-ugamma))**(1-sigma))/(1-sigma) !NO KAPPA IN UTILITY: zero labor
					endif
				else 
					if (cns > 0) then
						utl = ((1+delta*(h-1))*(cns**ugamma*ltot_scaled**(1-ugamma))**(1-sigma))/(1-sigma) !NO KAPPA IN UTILITY: zero labor
					else
						utl = -1.0d5
					endif	
				endif
			else !labor is larger than zero; the changes are in utility functon: notice KAPPA and XI there
				if (anxt == 0) then !next-period assets are zero: agent is eligible for minimal consumption level					
					if (cns > cmin_scaled) then
						utl = ((1+delta*(h-1))*(cns**ugamma*(ltot_scaled-lbr-kappa_scaled)**(1-ugamma))**(1-sigma))/(1-sigma) 
					else
						cns = cmin_scaled !consumption floor
						utl = ((1+delta*(h-1))*(cns**ugamma*(ltot_scaled-lbr-kappa_scaled)**(1-ugamma))**(1-sigma))/(1-sigma)
					endif
				else 
					if (cns > 0) then
						utl = ((1+delta*(h-1))*(cns**ugamma*(ltot_scaled-lbr-kappa_scaled)**(1-ugamma))**(1-sigma))/(1-sigma)
					else
						utl = -1.0d5
					endif	
				endif
			endif
						
		end subroutine cons_util_calc
		
end subroutine labchoice
