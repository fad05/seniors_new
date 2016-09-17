subroutine labchoice_new (lchoice, cons, util, acur,anext, h, wage, mexp)
    use parameters
    implicit none
    !inputs
    real (kind = 8), intent(in) :: acur, anext, wage,mexp !current assets, next period assets, current wage, mexp
    integer, intent(in) :: h !health status, 1 or 2
    !outputs
    real (kind = 8), intent(out) :: lchoice !optimal labor choice
    real (kind = 8), intent(out) :: cons !consumption under optimal choice
    real (kind = 8), intent(out) :: util !value of utility under optimal choice
    
    !BODY OF SUBROUTINE   
    
    lchoice = -5 ! initialize with negative value
    !1.Find max on the grid ad use it as starting value for optimization
    !2.Use newton method to solve for max, using analytical derivtives. Here cons is function of labor, thus I need 
    !derivatives wrt to both variables c and l, where c is expressed through l.
    
    return !end of subroutine
    contains
		subroutine findmin(l0,epsl,epsu,lopt,uopt)
			!Newton's method to find minimum of the utility function
			real (kind = 8), intent(in) :: l0 !initial guess
			real (kind = 8), intent(in) ::epsl, epsu !tolerance for l and util
			real (kind = 8), intent(out) :: lopt, uopt !solution and utility at solution
			
			!local vars
			real (kind = 8) f1, f2, l
			
			
			!===========================================================================
			
			!Check inputs
			if (l0<0) then
				print *, 'Choose another guess, labor can''t be negative'
				read (*,*)
				return
			elseif (l0>lmax) then
				print *, 'Choose another guess, labor can''t be this large'
				read (*,*)
				return
			end
			
			!If anext>0, no consumption floor provided => derivatives are continuous,
			!and we can use Newton's method
			if (anext >0) then !Use Newton
				l = l0
			
				!calculate first derivative analytically
				f1 = ((1+delta*(h-1))*(ugamma*wage*(ltot-l-kappa)^(1-ugamma)
			else !Use amoeba search, since we have discontinuity at cmin
			
			endif 
		
		
		end subroutine findmin
    
end subroutine labchoice_new
