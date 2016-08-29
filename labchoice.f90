subroutine labchoice (lchoice, acur,anext, h, wage, mexp)
    use parameters
    implicit none
    real (kind = 8), intent(in) :: acur, anext, wage,mexp !current assets, next period assets, current wage, mexp
    integer, intent(in) :: h !health status, 1 or 2
    real (kind = 8), intent(out) :: lchoice !optimal labor choice

    !variables necessary to use minpack and local variables
    integer ( kind = 4 ) info, i, u_max_index
    integer (kind = 4) :: counter = 0 !counter for search attempts
    real ( kind = 8 ) fvec(1), lc(1), lmax, labgrid_max, uchoice, ctemp, lguess
    real (kind = 8 ) bma ! bma stands for "b medicaid"; computed in a way such that consuption floor cmin is satisfied; in this case anext MUST be 0
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
    !check some points on the grid of labor, looking for good starting point for solver
    call linspace(0.0d0,lmax,100,labgrid)
    
    !two possible cases, depending on next period assets:
    !1. anext = 0; in this case an individual with consumption lower than cmin (due to medical expenses) is entitled to a Medicaid bobus bma to compensate his expenses exactly to cmin
    !2. anext > 0; if here consumption is lower than cmin, the case is not valid: agend sacrifices consumption for savings. If in this case consumption is lover than cmin, agent is not compensated;
    !if consumption is negative, the utility is set to large negative number, in order never to be optimal.
    
    if (anext == 0) then !next-period savings are 0; check consumption floor!
		cgrid(1) = (1+r)*acur-mexp-anext !uncompensated consumption on a first gridpoint 
		if (cgrid(1) > cmin) then !consumption is larger than floor value
			ugrid(1) = ((1+delta*(h-1))*(cgrid(1)**ugamma*ltot**(1-ugamma))**(1-sigma))/(1-sigma)
		else !consumption value is too low, and no savings for next period
			cgrid(1) = cmin
			ugrid(1) = ((1+delta*(h-1))*(cgrid(1)**ugamma*ltot**(1-ugamma))**(1-sigma))/(1-sigma)
		endif
		do i = 2,100
			cgrid(i) = (1+r)*acur+wage*labgrid(i)-mexp-anext !uncompensated consumption on a ith gridpoint 
			if (cgrid(i) > cmin) then
				ugrid(i) = ((1+delta*(h-1))*(cgrid(i)**ugamma*(ltot-labgrid(i)-kappa)**(1-ugamma))**(1-sigma))/(1-sigma)
			else
				cgrid(i) = cmin
				ugrid(i) = ((1+delta*(h-1))*(cgrid(i)**ugamma*(ltot-labgrid(i)-kappa)**(1-ugamma))**(1-sigma))/(1-sigma)
			endif
		enddo
    else !anext >0; agent is not entitled to compensation
		cgrid(1) = (1+r)*acur-mexp-anext
		if (cgrid(1) > 0) then
			ugrid(1) = ((1+delta*(h-1))*(cgrid(1)**ugamma*ltot**(1-ugamma))**(1-sigma))/(1-sigma)
		else
			ugrid(1) = -1.0d5 !if cons is zero or negative, utility is highly negative
		endif
		
		do i = 2,100
			cgrid(i) = (1+r)*acur+wage*labgrid(i)-mexp-anext
			if (cgrid(i)>0) then
				ugrid(i) = ((1+delta*(h-1))*(cgrid(i)**ugamma*(ltot-labgrid(i)-kappa)**(1-ugamma))**(1-sigma))/(1-sigma)				
			else
				ugrid(i) = -1.0d5
			endif
		enddo
    endif
    u_max_index = maxloc(ugrid,1) !index of maximal utility on the grid
    labgrid_max = labgrid(u_max_index)
    lguess = labgrid_max !initial guess on optimal labor choice OFF the grid: maximum ON the grid

 10   lc(1) = lguess 
    call lmdif1(lc_minpack, 1, 1, lc, fvec,tol, info) !find optimal labor choice, unconstrained (if possible)
    uchoice = -fvec(1)
    if (info>=1 .AND. info<=3 .AND. lc(1)>=0.0d0 .AND. lc(1) <= lmax) then !optimal solution is found and it is within feasible range
		print *, info
        lchoice = lc(1)
        counter = 0
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
		lguess = lguess + max((-1)**counter*counter*100,0) !try another guess
		print *, info
		if (counter > 5) then
			lchoice = labgrid_max !in this case, pick labor that corresponds to a maximum on the grid
			print *, 'Labor on the grid is chosen:', lchoice
		else
			goto 10 !make a new search attempt with different guess
		endif
    else 
		print *, 'Hell knows what''s going on!'
		print *, char(7) !beep
        call sleep(120)
    endif

    return !end of subroutine
    !function that allows to use minpack
    contains
        subroutine lc_minpack( m, n, x, fvec1, iflag )
            integer ( kind = 4 ) n, m
            real (kind = 8) fvec1(m)
            integer ( kind = 4 ) iflag
            real (kind = 8)  x(n)

            !internal variables
            real (kind = 8) lab, cl

            !body
            lab = x(n)
			if (lab < 0) then
				lab = 0 !no negative hours are feasible
			endif
			
			if (lab > lmax) then !choice is infeasible
				fvec(1) = 1.0d5 !large POSITIVE utility, since we're minimizing; we can't afford violation of this constraint; thus it's never optimal to choose this amount of hours
			else
				cl = (1+r)*acur+wage*lab-mexp-anext
				!Now, calculate consumption depending of anext (whether or not there is a consumption floor satisfied)
				!Remember, here fvec is NEGATIVE of utility, so to find a minimum!
				if (anext == 0) then
					bma = 0
					if (cl < cmin) then
						bma = cmin - cl
						cl = cmin
					endif
					fvec1(1) = - ((1+delta*(h-1))*(cl**ugamma*(ltot-lab-kappa)**(1-ugamma))**(1-sigma))/(1-sigma)
				else
					bma = 0
					if (cl > 0) then 
						fvec1(1) = - ((1+delta*(h-1))*(cl**ugamma*(ltot-lab-kappa)**(1-ugamma))**(1-sigma))/(1-sigma)
					else
						fvec(1) = 1.0d5 !large POSITIVE utility, since we're minimizing
					endif
				endif
            endif
        end subroutine lc_minpack

end subroutine labchoice
