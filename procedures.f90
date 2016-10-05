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

subroutine aime_calc(cohort,age,benefit,aime)
!routine calculates current AIME given current benefit input,
!cohort (1 or 2) and current age, 
use parameters
integer, intent(in) :: cohort, age
real (kind = 8), intent(in) :: benefit
real (kind = 8), intent(out) :: aime

!local variables
real (kind = 8) pia

!First we have to calculate PIA (real value, as if we're at NRA)

if (age<62) then
	pia = benefit
elseif (age>=nra(cohort)) then !get full benefit plus whatever bonus (which is only up to age 70)
	pia = benefit/(1+credit_delret(cohort)*(min(age,70)-nra(cohort)))
elseif (cohort == 2 .AND. nra(cohort)-age>3) then !62 to 63 for second cohort
	pia = benefit/(0.8-penalty_long*(nra(cohort)-age-3))
else !62 to 65 for first cohort, 64 to 67 for second cohort
	pia = benefit/(1-penalty_er3*(nra(cohort)-age))
endif

!Once we've calculated PIA, we can calculate AIME
	
if (pia <= aime_bend(1,cohort)*0.9) then !lower than firs bendpoint
	aime = pia/0.9
elseif (pia > aime_bend(1,cohort)*0.9 .AND. pia <= aime_bend(1,cohort)*0.9 + 0.32*(aime_bend(2,cohort)-aime_bend(1,cohort))) then !between first and second bendpoint
	aime = (pia - aime_bend(1,cohort)*0.9)/0.32 + aime_bend(1,cohort)
else !pia > aime_bend(1,cohort)*0.9 + 0.32*(aime_bend(2,cohort)-aime_bend(1,cohort))  : larger than third bendpoint
	aime = (pia-0.9*aime_bend(1,cohort)-0.32*(aime_bend(2,cohort)-aime_bend(1,cohort)))/0.15+aime_bend(2,cohort)
endif
	
end subroutine	aime_calc 

subroutine benefit_calc(cohort,age,aime,benefit)
!A reverse of aime_calc: given current aime, age, cohort
!a procedure calculates correct benefit level.
use parameters
integer, intent(in) :: cohort, age
real (kind = 8), intent(in) :: aime 
real (kind = 8), intent(out) :: benefit

!local variables
real (kind = 8) pia, cap

!First, we have to calculate PIA from AIME

if (aime <= aime_bend(1,cohort)) then !less than first bendpoint
	pia = 0.9*aime !90% of his/her average indexed monthly earnings below the first $aime_bend(1)
elseif (aime>aime_bend(1,cohort) .AND. aime <=aime_bend(2,cohort)) then !between first and second bendpoint
	pia = 0.9*aime_bend(1,cohort) + 0.32*(aime-aime_bend(1,cohort)) !90% of the first $aime_bend(1) plus 32 percent of his/her average indexed monthly earnings over $aime_bend(1) and through $aime_bend(2)
else !aime is over the second bendpoint
	pia = 0.9*aime_bend(1,cohort) + 0.32*(aime_bend(2,cohort)-aime_bend(1,cohort)) + 0.15*(aime-aime_bend(2,cohort)) !on top of previous, 15% of aime over second bendpoint
endif

!Second, from PIA we recover benefit level given application age bonuses or penalties, and applying benefit cap
if (age<62) then
	benefit = 0 !an individual can't get benefits this early
elseif (age>=nra(cohort)) then !full benefit plus bonus (bonus is only up to age 70)
	benefit = pia*(1+credit_delret(cohort)*(min(age,70)-nra(cohort)))
	!remember that we have benefit cap (which is also adjusted by age of applicaton in the similar way)
	cap = sscap_nra(cohort)*(1+credit_delret(cohort)*(min(age,70)-nra(cohort)))
	benefit = min(benefit,cap)
elseif (cohort == 2 .AND. nra(cohort)-age>3) then !62 to 63 for second cohort
	benefit = pia*(0.8-penalty_long*(nra(cohort)-age-3))
	!benefit cap
	cap = sscap_nra(cohort)*(0.8-penalty_long*(nra(cohort)-age-3))
	benefit = min(benefit,cap)
else !62 to 65 for first cohort, 64 to 67 for second cohort
	benefit = pia*(1-penalty_er3*(nra(cohort)-age))
	!benefit cap
	cap = sscap_nra(cohort)*(1-penalty_er3*(nra(cohort)-age))
	benefit = min(benefit,cap)
endif

end subroutine benefit_calc

subroutine bequest_value(assets,beq)
use parameters
real (kind = 8), dimension(grid_asset), intent(in) :: assets
real (kind = 8), dimension(grid_asset), intent(out) :: beq

beq = eta*((assets+d)/scale_factor)**(1-sigma)/(1-sigma) !notice scale factor here

end subroutine bequest_value

end module procedures
