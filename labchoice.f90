subroutine labchoice (lchoice, cons, util, acur,anext, h, wage, mexp)
    use parameters
    implicit none
    real (kind = 8), intent(in) :: acur, anext, wage,mexp !current assets, next period assets, current wage, mexp
    integer, intent(in) :: h !health status, 1 or 2
    real (kind = 8), intent(out) :: lchoice !optimal labor choice
    real (kind = 8), intent(out) :: cons !consumption under optimal choice
    real (kind = 8), intent(out) :: util !value of utility under optimal choice

    !variables necessary to use minpack and local variables
    integer ( kind = 4 ) info, i, u_max_index
    integer (kind = 4) :: counter = 0 !counter for search attempts
    real (kind = 8) fvec(1), lc(1), lmax, labgrid_max, uchoice, ctemp, lguess, utemp
    real (kind = 8) u0, c0 !utility and consumption IF labor is 0, GIVEN assets, mexp, wage
    real (kind = 8) bma ! bma stands for "b medicaid"; computed in a way such that consuption floor cmin is satisfied; in this case anext MUST be 0
    real (kind = 8), dimension(100) :: labgrid, cgrid, ugrid

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
        subroutine  linspace(a,b,n,outvec)
            real (kind = 8), intent(in) :: a, b !starting point, endpoint
            integer, intent(in) :: n !number of points
            real (kind = 8), dimension(n), intent(out) :: outvec !output vector
        end subroutine linspace
    end interface

    !BODY OF SUBROUTINE
    

    lchoice = -5 ! initialize with negative value
    lmax = ltot - kappa - 1.0d-6 !maximum possible labor
    !check some points on the grid of labor, looking for good starting point for solver.
    !The grid spans from small positive (nonzero) number to max possible labor (ltot-kappa) minus small number.
    !The grid has 100 equally spaced nodes.
    call linspace(1.0d-6,lmax,100,labgrid)
    
    do i = 1,100 !on this labor grid, find corresponding utility and consumption
	    call cons_util_calc(labgrid(i), wage, acur, anext, mexp, cgrid(i), ugrid(i))
    enddo
    
    u_max_index = maxloc(ugrid,1) !index of maximal utility on the grid
    labgrid_max = labgrid(u_max_index)
    lguess = labgrid_max !initial guess on optimal labor choice OFF the grid: maximum ON the grid
    
    !IMPORTANT REFERENCE POINT: Utility and consumption IF labor is 0:
    call cons_util_calc(0.0d0,wage,acur,anext,mexp,c0,u0)
    
    ! IF maximum is achieved almost at 0 on gridpoint 1, we compare it to utility at 0
    if (u_max_index == 1) then 
		call cons_util_calc(labgrid(1),wage,acur, anext,mexp,ctemp,utemp)
		if (u0>utemp) then
			lchoice = 0.0d0
		else
			print*, 'Utility at 0 is lower than utility almost at zero, something''s wrong. Press any button and check the labchoice.f95'
			read (*,*)
		endif	
	!ELSE, IF MAXIMUM is achieved elsewhere on the grid, we use minpack to optimize.
	else    
	    !Look for maximum utility on the interval (0,lmax] using grid maximum as starting point    
		10  lc(1) = lguess 
		
	    call lmdif1(lc_minpack, 1, 1, lc, fvec,tol, info) !find optimal labor choice
	    uchoice = -fvec(1) !utility equals negative of the outcome of the procedure
	    print *, uchoice
	    if (info>=1 .AND. info<=3 .AND. lc(1)>=0.0d0 .AND. lc(1) <= lmax) then !optimal solution is found and it is within feasible range
			print *, info
	        lchoice = lc(1)
	        utemp = uchoice
	        call cons_util_calc(lchoice, wage, acur, anext, mexp, cons, util)  
	        if (real(utemp) .NE. real(util)) then
				print *, 'Utility outcome of minimizer not equal to utility calculated with optimal labor!vPress any button:'
				read (*,*)
				print *, utemp
				print *, util
				!call sleep(3)
	        endif   
	        counter = 0 !Nullify counter
	    elseif (info>=1 .AND. info<=3 .AND. lc(1)<0.0d0) then !optimal solution found and it's negative: something's wrong with the code
	        print *, info
	        print *, 'Negative hours! Check the code!'
	        print *, char(7) !beep
	        call sleep(120)
	    elseif (info>=1 .AND. info<=3 .AND. lc(1)>lmax) then !optimal solution found and it's outside the range; set it to 
			print *, info
			print *, 'Hours too large! Check the code!'
			print *, char(7) !beep
			call sleep(120)
		elseif (info >= 4) then
			counter = counter+1
			lguess = lguess + max((-1)**counter*counter*100,0) !update initial guess
			print *, info
			if ((counter > 10 .AND. labgrid_max >0) .OR. (counter > 3 .AND. labgrid_max == 0)) then !if maximum utility is at zero labor, we only repeat search 3 times; if max utility on the grid achieved NOT in zero labor, we search 10 times with different initial guesses 
				lchoice = labgrid_max !in this case, pick labor that corresponds to a maximum on the grid
				if (lchoice /= 0.0d0) then
					print *,  'Nonzero labor on the grid is chosen, somethings''s wrong:', lchoice
					call sleep(120)
				else !ZERO labor, no KAPPA in utility
					print *, 'Zero labor is chosen from the grid!'
					!call sleep(3)
					call cons_util_calc(0.0d0,wage,acur,anext,mexp,cons,util)
				endif
				counter = 0 !Nullify counter
			else
				goto 10 !make a new search attempt with different guess
			endif
	    else 
			print *, 'Hell knows what''s going on!'
			print *, char(7) !beep
	        call sleep(120)
	    endif
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
				call cons_util_calc(lab, wage, acur, anext, mexp, cl, ul)
				fvec(1) = -ul ! we're minimizing, so in order to find maximum, we need negative of utility function
			endif
			x(n) = lab
			if (isnan(x(n)) .OR. isnan(fvec(1))) then
				print *, x(n) 
				print *, fvec(1)
			endif
			print *, x(n)
			print *, fvec(1)
			return
        end subroutine lc_minpack
        
        subroutine cons_util_calc(lbr, wg, acrnt, anxt, mxp, cns, utl)
			real (kind = 8), intent(in) :: lbr, wg, acrnt, anxt, mxp
			real (kind = 8), intent(out) :: cns, utl
			!BODY
			!two possible cases, depending on next period assets:
			!1. anext = 0; in this case an individual with consumption lower than cmin (due to medical expenses) is entitled to a Medicaid bobus bma to compensate his expenses exactly to cmin
			!2. anext > 0; if here consumption is lower than cmin, the case is not valid: agend sacrifices consumption for savings. If in this case consumption is lover than cmin, agent is not compensated;
			!if consumption is negative, the utility is set to large negative number, in order never to be optimal.
			!Another fork is whether labor is zero or not. If it is, there is no utility penalty, expressed as KAPPA.
			cns = (1+r)*acrnt + lbr*wg - anxt - mxp !uncompensated consumption; in the case when lbr==0, labor income is ZERO
			if (lbr == 0) then !labor is zero
				if (anxt == 0) then !next-period assets are zero: agent is eligible for minimal consulption level					
					if (cns > cmin) then !consumption is hgher than minimal, no transfers
						utl = ((1+delta*(h-1))*(cns**ugamma*ltot**(1-ugamma))**(1-sigma))/(1-sigma) !NO KAPPA IN UTILITY: zero labor
					else
						cns = cmin !agent reseives a granted consumption floor
						utl = ((1+delta*(h-1))*(cns**ugamma*ltot**(1-ugamma))**(1-sigma))/(1-sigma) !NO KAPPA IN UTILITY: zero labor
					endif
				else 
					if (cns > 0) then
						cns = ((1+delta*(h-1))*(cns**ugamma*ltot**(1-ugamma))**(1-sigma))/(1-sigma) !NO KAPPA IN UTILITY: zero labor
					else
						cns = -1.0d5
					endif	
				endif
			else !labor is larger than zero; the changes are in utility functon: notice KAPPA there
				if (anxt == 0) then !next-period assets are zero: agent is eligible for minimal consumption level					
					if (cns > cmin) then
						utl = ((1+delta*(h-1))*(cns**ugamma*(ltot-lbr-kappa)**(1-ugamma))**(1-sigma))/(1-sigma) 
					else
						cns = cmin !consumption floor
						utl = ((1+delta*(h-1))*(cns**ugamma*(ltot-lbr-kappa)**(1-ugamma))**(1-sigma))/(1-sigma)
					endif
				else 
					if (cns > 0) then
						cns = ((1+delta*(h-1))*(cns**ugamma*(ltot-lbr-kappa)**(1-ugamma))**(1-sigma))/(1-sigma)
					else
						cns = -1.0d5
					endif	
				endif
			endif
			
		end subroutine cons_util_calc
end subroutine labchoice
