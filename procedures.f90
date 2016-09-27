module procedures
implicit none
public :: linspace, locate_greater, linproj, import_data, import_matrices, import_vector
contains

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

function locate_greater(inarray,num)
	integer locate_greater
	real (kind = 8), intent(in), dimension(:) :: inarray
	real (kind = 8), intent(in) :: num
	!local variables
	real (kind = 8), dimension(:), allocatable :: temparray
	integer loc(1)
	
	allocate(temparray(size(inarray)))	
	temparray = inarray-num
	
	where (temparray>0.0d0)
	temparray = 1.0d0
	endwhere	
	
	loc = maxloc(temparray)	
	locate_greater = loc(1)	
	
	deallocate(temparray)
end function locate_greater

function linproj(x,a,b,fa,fb)
	! a function calculates a linearly approximated value of
	! some function f(x) on a point x, a<x<b, given f(a) and f(b)
	real (kind = 8), intent(in) :: x, a, b, fa, fb
	real (kind = 8) linproj
	
	linproj =  fa+(fb-fa)*(x-a)/(b-a)
end function linproj

subroutine import_vector(datafile,outdata)
    !The subroutine to read one-dimensional vector of data from 'datafile' to a vector "outdata"
    implicit none
    character*200, intent(in) :: datafile ! address of a data file to read
    real (kind = 8), dimension(:), intent(out) :: outdata

    !local variables
    integer n !vector size
    integer i

    !body
    !print *, datafile
    n = size(outdata)
    open(unit = 1, file = datafile, action = "read")

    do i = 1,n
        read(1,*) outdata(i)
    enddo

    close(1)
    !print *, outdata(1)

end subroutine import_vector


subroutine import_matrices(datafile,outdata)
    !The subroutine to read data from 'datafile' to a variable "outdata"
    implicit none
    character*200, intent(in) :: datafile ! address of a data fie to read
    real (kind = 8), dimension(:,:,:), intent(out) :: outdata

    !local variables
    integer n,m,k !matrix sizes
    integer i,j

    !body
    n = size(outdata,1)
    m = size(outdata,2)
    k = size(outdata,3)

    open(unit = 1, file = datafile, action = "read")

    do j = 1,k
        do i = 1,n
            read(1,*) outdata(i,:,j)
        enddo
    enddo

    close(1)
end subroutine import_matrices

subroutine import_data(datafile,outdata)
    !The subroutine to read data from 'datafile' to a variable "outdata"
    implicit none
    character*200, intent(in) :: datafile ! address of a data fie to read
    real (kind = 8), dimension(:,:), intent(out) :: outdata

    !inner variables
    integer i, matsize, namesize
    character*30, dimension(:) , allocatable:: names
    !Beginning of the program
    matsize = size(outdata,1)
    namesize = size(outdata,2)
    allocate(names(namesize))
    outdata = 42

    open(unit = 1, file = datafile, action = "read")

    read(1,*) names(:)
    !print *, names
    do i = 1,matsize
        read(1,*) outdata(i,:)
    enddo
    close(1)

end subroutine import_data

end module procedures
