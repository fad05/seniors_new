!From Numerical Recipes, Press et al
subroutine amoeba(p,y,ftol,func,iter,info)
use nrtype; use nrutil, only: assert_eq, imaxloc, iminloc, nrerror, swap
implicit none
integer(I4B), intent(out) :: iter, info
real(DP), intent(in) :: ftol
real(DP), dimension(:), intent(inout) :: y
real(DP), dimension(:,:), intent(inout) :: p
interface
    function func(x)
    use nrtype
    implicit none
    real(DP), dimension(:), intent(in) :: x
    real(DP) :: func
    end function func
end interface
integer(I4B), parameter :: ITMAX = 10000
real(DP), parameter :: TINY = 1.0e-10

integer(I4B) :: ihi, ndim
real(DP), dimension(size(p,2)) :: psum
call amoeba_private
return
contains

subroutine amoeba_private
implicit none
integer(I4B) :: i, ilo, inhi
real(DP) :: rtol, ysave, ytry, ytmp
ndim = assert_eq(size(p,2),size(p,1)-1,size(y)-1,'amoeba')
iter = 0
psum(:) = sum(p(:,:),dim=1)
do
    ilo = iminloc(y(:))
    ihi = imaxloc(y(:))
    ytmp = y(ihi)
    y(ihi) = y(ilo)
    inhi = imaxloc(y(:))
    y(ihi) = ytmp
    rtol = 2.0_sp*abs(y(ihi)-y(ilo))/(abs(y(ihi))+abs(y(ilo))+TINY)
    print *,'rtol =', rtol
    print *, 'abs error value:', abs(y(ilo))
    if (rtol<ftol) then
        call swap(y(1),y(ilo))
        call swap(p(1,:),p(ilo,:))
        info = 0 !solution found
        return
    endif
    if (iter>=ITMAX) then
    call nrerror('ITMAX exceeded in amoeba')
    info = 1 !no solution found
    return
    !read *
    endif
    ytry = amotry(-1.0_dp)
    iter = iter+1
    if (ytry<=y(ilo)) then
        ytry = amotry(2.0_dp)
        iter = iter+1
    elseif (ytry>=y(inhi)) then
        ysave = y(ihi)
        ytry = amotry(0.5_dp)
        iter = iter+1
        if (ytry>=ysave) then
            p(:,:) = 0.5_dp*(p(:,:)+spread(p(ilo,:),1,size(p,1)))
            do i = 1,ndim+1
                if (i/=ilo) y(i) = func(p(i,:))
            enddo
            iter = iter + ndim
            psum = sum(p(:,:),dim = 1)
        endif
    endif
enddo
end subroutine amoeba_private

function amotry(fac)
implicit none
real(DP), intent(in) :: fac
real(DP) :: amotry

real(DP) :: fac1, fac2, ytry
real(DP), dimension(size(p,2)) :: ptry

fac1 = (1.0_sp-fac)/ndim
fac2 =  fac1-fac
ptry(:) = psum(:)*fac1 - p(ihi,:)*fac2
ytry =  func(ptry)
if (ytry<y(ihi)) then
    y(ihi)=ytry
    psum(:) = psum(:) - p(ihi,:) + ptry(:)
    p(ihi,:) = ptry(:)
endif
amotry = ytry
end function amotry
end subroutine amoeba
