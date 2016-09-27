subroutine labchoice (lchoice, cons, util, acur,anext, h, wage, mexp,ss)
	use procedures
    use parameters
    implicit none
    real (kind = 8), intent(in) :: acur, anext, wage,mexp,ss !current assets, next period assets, current wage, mexp, social scurity
    integer, intent(in) :: h !health status, 1 or 2
    real (kind = 8), intent(out) :: lchoice !optimal labor choice
    real (kind = 8), intent(out) :: cons !consumption under optimal choice
    real (kind = 8), intent(out) :: util !value of utility under optimal choice

    !SCALED VARIABLES
    real (kind = 8) acur_scaled, anext_scaled, mexp_scaled, ss_scaled
    real (kind = 8) :: ltot_scaled
	real (kind = 8) :: kappa_scaled
	real (kind = 8) :: cmin_scaled
	!variables necessary to use minpack and local variables
    integer (kind = 4) info, i, u_max_index, u_min_index
    integer (kind = 4) :: counter = 0 !counter for search attempts
    real (kind = 8) fvec(1), lc(1), lmax, labgrid_max, uchoice, ctemp, lguess, utemp
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
	
	!Scale parameters of utility functon!
	ltot_scaled = ltot/scale_factor
	kappa_scaled = kappa/scale_factor
	cmin_scaled = cmin/scale_factor
	
    !check some points on the grid of labor, looking for good starting point for solver.
    !The grid spans from small positive (nonzero) number to max possible labor (ltot-kappa) minus small number, all scaled
    !The grid has 100 equally spaced nodes.
    call linspace(1.0d-3,lmax,100,labgrid) !GRID is ALSO SCALED, starting with ONE WORKING HOUR
!=========================END SCALING========================================================================================
	
    

    
    do i = 1,100 !on this labor grid, find corresponding utility and consumption
	    call cons_util_calc(labgrid(i), wage, acur_scaled, anext_scaled, mexp_scaled,ss_scaled, cgrid(i), ugrid(i))
    enddo

    
    u_max_index = maxloc(ugrid,1) !index of maximal utility on the grid
    labgrid_max = labgrid(u_max_index)
    lguess = labgrid_max !initial guess on optimal labor choice OFF the grid: maximum ON the grid
    
    !IMPORTANT REFERENCE POINTS:
    ! Utility and consumption IF labor is 0:
    call cons_util_calc(0.0d0,wage,acur_scaled,anext_scaled,mexp_scaled,ss_scaled,c0,u0)
    ! Utility and cons at first grdpoint (1 hour of work)
    call cons_util_calc(labgrid(1),wage,acur_scaled, anext_scaled,mexp_scaled,ss_scaled,ctemp,utemp)
    
    ! IF maximum is achieved almost at 0 on gridpoint 1, we compare it to utility at 0
    if (u_max_index == 1 .AND. u0>=utemp) then 
			lchoice = 0.0d0
			util = u0
			cons = c0
	!ELSE, IF MAXIMUM is achieved elsewhere on the grid, we optimize.
	else   
!	1.Bracketing method, slowest but sure
!	We're looking for maximum 
		if (u_max_index == 1 .AND. u0<utemp) then !IF UTILITY ALMOST AT 0 is LARGER THAN UTILITY AT 0
			call cons_util_calc(0.0d0,wage,acur_scaled,anext_scaled,mexp_scaled,ss_scaled,c0,u0)
			call cons_util_calc(labgrid(1),wage,acur_scaled, anext_scaled,mexp_scaled,ss_scaled,ctemp,utemp)
			print *, 'U at 0:', u0
			print *, 'U at labgrid(1):', utemp
			print *, 'Utility at 0 is lower than utility almost at zero, something''s wrong. Press any button and check the labchoice.f95'
			read (*,*)
			print *, 'Proceeding with bracketing method'
			a = 0.0d0
			ua = u0
			b = labgrid(1)
			ub = ugrid(1)
			c = labgrid(2)
			uc = ugrid(2)
		else	
			a = labgrid(u_max_index-1)
			ua = ugrid(u_max_index-1)
			b = labgrid(u_max_index)
			ub = ugrid(u_max_index)
			c = labgrid(u_max_index+1)
			uc = ugrid(u_max_index+1)
		endif
		do while (c-a>1.0d-3) !the tolerance is ohe hour (1E-3 of 1000s of hours)
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
		
		lchoice = c	!solution	
		call cons_util_calc(lchoice, wage, acur_scaled, anext_scaled, mexp_scaled, ss_scaled,cons, util)
		
!		!scale back
!		lchoice = lchoice*scale_factor
!		cons = cons*scale_factor
!		util = util*scale_factor**(1-sigma)
		
	
		
!	    !Look for maximum utility on the interval (0,lmax] using grid maximum as starting point    
!		10  lc(1) = lguess 
		
!	    call lmdif1(lc_minpack, 1, 1, lc, fvec,tol, info) !find optimal labor choice
!	    uchoice = -fvec(1) !utility equals negative of the outcome of the procedure
!	    !print *, uchoice
!	    if (info>=1 .AND. info<=3 .AND. lc(1)>=0.0d0 .AND. lc(1) <= lmax) then !optimal solution is found and it is within feasible range
!			print *, info
!	        lchoice = lc(1)
!	        utemp = uchoice
!	        call cons_util_calc(lchoice, wage, acur_scaled, anext_scaled, mexp_scaled, cons, util)  
!	        if (real(utemp) .NE. real(util)) then
!				print *, 'U outcome of minimzer:', real(utemp)
!				print *, 'U calculated with optimal labor:', real(util)
!				print *, 'Utility outcome of minimizer not equal to utility calculated with optimal labor! Press any button:'
!				read (*,*)				
!	        endif   
!	        counter = 0 !Nullify counter
!	    elseif (info>=1 .AND. info<=3 .AND. lc(1)<0.0d0) then !optimal solution found and it's negative: something's wrong with the code
!	        print *, info
!	        print *, 'Negative hours! Check the code!'
!	        print *, char(7) !beep
!	        call sleep(120)
!	    elseif (info>=1 .AND. info<=3 .AND. lc(1)>lmax) then !optimal solution found and it's outside the range; set it to 
!			print *, info
!			print *, 'Hours too large! Check the code!'
!			print *, char(7) !beep
!			call sleep(120)
!		elseif (info >= 4) then
!			counter = counter+1
!			lguess = lguess + max((-1)**counter*(counter*100/scale_factor),0) !update initial guess
!			print *, info
!			if ((counter > 10 .AND. labgrid_max >0) .OR. (counter > 3 .AND. labgrid_max == 0)) then !if maximum utility is at zero labor, we only repeat search 3 times; if max utility on the grid achieved NOT in zero labor, we search 10 times with different initial guesses 
!				lchoice = labgrid_max !in this case, pick labor that corresponds to a maximum on the grid
!				if (lchoice /= 0.0d0) then
!					print *,  'Nonzero labor on the grid is chosen, somethings''s wrong:', lchoice
!					call sleep(120)
!				else !ZERO labor, no KAPPA in utility
!					print *, 'Zero labor is chosen from the grid!'
!					!call sleep(3)
!					call cons_util_calc(0.0d0,wage,acur_scaled, anext_scaled, mexp_scaled,cons,util)
!				endif
!				counter = 0 !Nullify counter
!			else
!				goto 10 !make a new search attempt with different guess
!			endif
!	    else 
!			print *, 'Hell knows what''s going on!'
!			print *, char(7) !beep
!	        call sleep(120)
!	    endif
	endif ! END IF maximum achieved not on the first value of the grid

    return !end of subroutine
    !function that allows to use minpack
    contains
        subroutine lc_minpack( m, n, x, fvec1, iflag )
            integer ( kind = 4 ) n, m
            real (kind = 8) fvec1(m)
            integer ( kind = 4 ) iflag
            real (kind = 8)  x(n)

            !internal variables
            real (kind = 8) lab, cl, ul

            !body
            lab = x(n)
            if (lab<0) then
				print *, 'Input labor is negative! Push any button to proceed.'
				read (*,*)
				lab = 1.0d-6
			endif 
			!if (lab < 0) then
			!	lab = 0 !no negative hours are feasible
			!endif
			
			if (lab > lmax) then !choice is infeasible
				fvec(1) = 1.0d5 !large POSITIVE utility, since we're minimizing; we can't afford violation of this constraint; thus it's never optimal to choose this amount of hours
			else !lab is between zero and lmax, thus we just calculate consumption and utility
				call cons_util_calc(lab, wage, acur_scaled, anext_scaled, mexp_scaled, ss_scaled,cl, ul)
				fvec(1) = -ul ! we're minimizing, so in order to find maximum, we need negative of utility function
			endif
			x(n) = lab
			if (isnan(x(n)) .OR. isnan(fvec(1))) then
				print *, x(n) 
				print *, fvec(1)
				print *, 'NAN. Press any button:'
				read (*,*)
			endif
			!print *, x(n)
			!print *, fvec(1)
        end subroutine lc_minpack
        
        subroutine cons_util_calc(lbr, wg, acrnt, anxt, mxp, ssec, cns, utl)
			real (kind = 8), intent(in) :: lbr, wg, acrnt, anxt, mxp, ssec
			real (kind = 8), intent(out) :: cns, utl
		
		
			!BODY
			!1.Calculate consumption and utility
			!two possible cases, depending on next period assets:
			!1. anext = 0; in this case an individual with consumption lower than cmin (due to medical expenses) is entitled to a Medicaid bobus bma to compensate his expenses exactly to cmin
			!2. anext > 0; if here consumption is lower than cmin, the case is not valid: agend sacrifices consumption for savings. If in this case consumption is lover than cmin, agent is not compensated;
			!if consumption is negative, the utility is set to large negative number, in order never to be optimal.
			!Another fork is whether labor is zero or not. If it is, there is no utility penalty, expressed as KAPPA.
			cns = (1+r)*acrnt + lbr*wg +ssec- anxt - mxp !uncompensated consumption; in the case when lbr==0, labor income is ZERO
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
			else !labor is larger than zero; the changes are in utility functon: notice KAPPA there
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
		
!		subroutine findmin(l0,epsl,epsu,lopt,uopt)
!		!Newton's method to find minimum of the utility function
!		real (kind = 8), intent(in) :: l0 !initial guess
!		real (kind = 8), intent(in) ::epsl, epsu !tolerance for l and util
!		real (kind = 8), intent(out) :: lopt, uopt !solution and utility at solution
		
!		!local vars
!		real (kind = 8) f1, f2, l
		
		
!		!===========================================================================
		
!		!Check inputs
!		if (l0<0) then
!			print *, 'Choose another guess, labor can''t be negative'
!			read (*,*)
!			return
!		elseif (l0>lmax) then
!			print *, 'Choose another guess, labor can''t be this large'
!			read (*,*)
!			return
!		endif
		
!		l = l0
!		!calculate first derivative analytically
!		!f1 = (1+delta*(h-1))*(ugamma*wage*(lmax-l)^(1-ugamma)*((1+r)*acur+l*wage-anext-mexp))
		
		
!		end subroutine findmin
end subroutine labchoice
