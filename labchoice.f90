subroutine labchoice (lchoice, acur,anext, h, wage, mexp)
    use parameters
    implicit none
    real (kind = 8), intent(in) :: acur, anext, wage,mexp !current assets, next period assets, current wage, mexp
    integer, intent(in) :: h !health status, 1 or 2
    real (kind = 8), intent(out) :: lchoice !optimal labor choice

    !variables necessary to use minpack and local variables
    integer ( kind = 4 ) info, i, u_max_index
    real ( kind = 8 ) fvec(1), lc(1), lmax, labgrid_max, uchoice
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

    !body

    lchoice = -5 ! initialize with negative value
    lmax = ltot - kappa - 1.0d-6 !maximum possible labor
    !check some points on the grid of labor, looking for good starting point for solver
    call linspace(0.0d0,lmax,100,labgrid)
    cgrid(1) = (1+r)*acur-mexp-anext
    ugrid(1) = ((1+delta*(h-1))*(cgrid(1)**ugamma*ltot**(1-ugamma))**(1-sigma))/(1-sigma)
    do i = 2,100
        cgrid(i) = (1+r)*acur+wage*labgrid(i)-mexp-anext
        ugrid(i) = ((1+delta*(h-1))*(cgrid(i)**ugamma*(ltot-labgrid(i)-kappa)**(1-ugamma))**(1-sigma))/(1-sigma)
    enddo
    u_max_index = maxloc(ugrid,1) !index of maximal utility on the grid
    labgrid_max = labgrid(u_max_index)

    lc(1) = labgrid_max !initial guess on optimal labor choice

    call lmdif1(lc_minpack, 1, 1, lc, fvec,tol, info) !find optimal labor choice, unconstrained (if possible)
    uchoice = -fvec(1)
    if (info>=1 .AND. info<=3 .AND. lc(1)>=0.0d0) then
        lchoice = lc(1)
    else !need to check negative hours!!
        print *, info
        if (lc(1)<0) then
            lchoice = 0 !CHECK THIS!!!!!!!
        elseif (uchoice>ugrid(u_max_index)) then !THIS IS UNFINISHED
            lchoice = lc(1)
        else
            lchoice = labgrid_max
        endif
    endif

!!!!!!========================== UNFINISHED PROCEDURE
    !now, I need to check if unconstrained optimal labor choice lies between 0 and kappa, and check utiliy on borders
    !THIS IS UNFINISHED
!!!!!!========================== UNFINISHED PROCEDURE

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

            cl = (1+r)*acur+wage*lab-mexp-anext

            fvec1(1) = - ((1+delta*(h-1))*(cl**ugamma*(ltot-lab-kappa)**(1-ugamma))**(1-sigma))/(1-sigma)

        end subroutine lc_minpack

end subroutine labchoice
