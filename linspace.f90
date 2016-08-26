subroutine  linspace(a,b,n,outvec)
real (kind = 8), intent(in) :: a, b !starting point, endpoint
integer, intent(in) :: n !number of points
real (kind = 8), dimension(n), intent(out) :: outvec !output vector
real (kind = 8) step
integer i

step = (b - a)/(n-1)
outvec(1) = a

do i = 2,n
    outvec(i) = outvec(i-1)+step
enddo

return
end subroutine linspace
