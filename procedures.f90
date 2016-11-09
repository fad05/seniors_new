module procedures
implicit none
!public :: linspace, logspace,locate_greater, linproj, import_data, import_matrices, import_vector
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

subroutine logspace(a,b,n,outvec)
	real (kind = 8), intent(in) :: a, b !starting point, endpoint
	integer, intent(in) :: n !number of points
	real (kind = 8), dimension(n), intent(out) :: outvec !output vector
	real (kind = 8) step
	integer i
	real (kind = 8) atemp, btemp, stpt, endpt

	btemp = b
	atemp = a
	
	if (b<a) then !inputs are incorrect
		print *, 'b should be larger than a! Input correct b:'
		read (*,*) btemp
		print *, 'a = ', atemp
		print *, 'b = ', btemp
	endif
		
	if (a<=0.0) then
		atemp = 1
		btemp = btemp+(atemp-a) !atemp-a is just shift
	endif
	
	stpt	=	log(atemp)
	endpt	=	log(btemp)
	call linspace(stpt,endpt,n,outvec)
	
	outvec = exp(outvec)
	
	if (a<=0) then
		outvec = outvec - (atemp-a) !shift backwards
	endif
end subroutine logspace

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
	
function lpsimp(x,b,f0,fb)
	! simplified version of linproj, when a = 0
	! some function f(x) on a point x, 0<x<b, given f(0) and f(b)
	real (kind = 8), intent(in) :: x, b, f0, fb
	real (kind = 8) lpsimp, fx
	
	fx = linproj(x,0.0d0,b,f0,fb)
	lpsimp = fx
	
end function lpsimp

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
	
if (pia <= aime_bend(1,cohort)*0.9) then !lower than first bendpoint
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
	benefit = pia !an individual can't get benefits this early
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

function distance2d(x1,y1,x2,y2)
!a function that calculates Euclidean metric of a distance between two ponts (x1,y1) and (x2,y2) on XY-plane
real (kind = 8), intent(in) :: x1,y1,x2,y2
real (kind = 8) distance2d
distance2d = dsqrt((x2-x1)**2+(y2-y1)**2)
end function distance2d

function income_tax(income)
use parameters
!function takes taxable income as input and
!spits out taxes to pay
real (kind = 8), intent(in) :: income 
real tax_rate, income_in_thousands
real income_tax

income_in_thousands = income/1000 !formula uses income in thousands of dollars of 2000
tax_rate = btax*(1.0-(1.0+stax*income_in_thousands**ptax)**(-1.0/ptax))
income_tax = income*tax_rate

end function income_tax

function taxable_amount(income,benefit)
!function takes whole income, benefit lvl and marital state as inputs and
!spits out taxable amount (which can be used as input to calculate income tax)
use parameters
real (kind = 8), intent(in) :: income, benefit
integer :: ms = 1
real taxable_amount, taxable_benefit

!local vars
real (kind = 8) provisional_income, first_tier, second_tier

provisional_income = income + 0.5*benefit


do while ((ms /= 1) .AND. (ms /= 2))
	print *, 'Wrong marital status'
	print *, 'Input correct marital status: 1 or 2'
	read (*,*) ms
enddo


first_tier = max(0.0d0,0.5*min(benefit,min(provisional_income-first_treshold(ms),second_treshold(ms)-first_treshold(ms))))
taxable_benefit = max(0.0d0,min(0.85*benefit,first_tier+0.85*(provisional_income-second_treshold(ms))))

taxable_amount = income+taxable_benefit !total taxable amount: income +benefits

end function taxable_amount

function quadlin(x,y,z,q,x0,y0,z0,q0,x1,y1,z1,q1, & 
v0000,v0001,v0010,v0011,v0100,v0101,v0110,v0111, &
v1000,v1001,v1010,v1011,v1100,v1101,v1110,v1111)
!4-dimensional linear interpolation
real (kind = 8), intent(in) :: 	x,y,z,q,x0,y0,z0,q0,x1,y1,z1,q1, &
								v0000,v0001,v0010,v0011,v0100,v0101,v0110,v0111, &
								v1000,v1001,v1010,v1011,v1100,v1101,v1110,v1111
								
real (kind = 8) quadlin, v
real (kind = 8) n0000,n0001,n0010,n0011,n0100,n0101,n0110,n0111,n1000,n1001,n1010,n1011,n1100,n1101,n1110,n1111
real (kind = 8) denom
							
!4D-linear interpolation
!(x,y,z,q) - coordinates of a point in 4d space
!(x0,y0,z0,q0) and (x1,y1,z1,q1) - polytope vertex nearest the origin and farthest from the origin ponts respectively
!v0000 - v1111: values of a function on a points (example: v1010 = v(x1,y0,z1,q0))

denom = (x1-x0)*(y1-y0)*(z1-z0)*(q1-q0)


n0000 = (x1-x)*(y1-y)*(z1-z)*(q1-q)/denom
n0001 = (x1-x)*(y1-y)*(z1-z)*(q-q0)/denom
n0010 = (x1-x)*(y1-y)*(z-z0)*(q1-q)/denom
n0011 = (x1-x)*(y1-y)*(z-z0)*(q-q0)/denom
n0100 = (x1-x)*(y-y0)*(z1-z)*(q1-q)/denom
n0101 = (x1-x)*(y-y0)*(z1-z)*(q-q0)/denom
n0110 = (x1-x)*(y-y0)*(z-z0)*(q1-q)/denom
n0111 = (x1-x)*(y-y0)*(z-z0)*(q-q0)/denom
n1000 = (x-x0)*(y1-y)*(z1-z)*(q1-q)/denom
n1001 = (x-x0)*(y1-y)*(z1-z)*(q-q0)/denom
n1010 = (x-x0)*(y1-y)*(z-z0)*(q1-q)/denom
n1011 = (x-x0)*(y1-y)*(z-z0)*(q-q0)/denom
n1100 = (x-x0)*(y-y0)*(z1-z)*(q1-q)/denom
n1101 = (x-x0)*(y-y0)*(z1-z)*(q-q0)/denom
n1110 = (x-x0)*(y-y0)*(z-z0)*(q1-q)/denom
n1111 = (x-x0)*(y-y0)*(z-z0)*(q-q0)/denom

v = v0000*n0000+v0001*n0001+v0010*n0010+v0011*n0011+	&
	v0100*n0100+v0101*n0101+v0110*n0110+v0111*n0111+	&
	v1000*n1000+v1001*n1001+v1010*n1010+v1011*n1011+	&
	v1100*n1100+v1101*n1101+v1110*n1110+v1111*n1111
	
quadlin = v

end function quadlin

!function trilin(x,y,z,x0,y0,z0,x1,y1,z1, & 
!v000,v001,v010,v011,v100,v101,v110,v111)
!!3-dimensional linear interpolation

!real (kind = 8), intent(in) :: 	x,y,z,q,x0,y0,z0,q0,x1,y1,z1,q1, &
!								v000,v001,v010,v011,v100,v101,v110,v111)
								
!real (kind = 8) trilin, v
!real (kind = 8) n000,n001,n010,n011,n100,n101,n110,n111
!real (kind = 8) denom
							
!!3D-linear interpolation
!!(x,y,z) - coordinates of a point in 3d space
!!(x0,y0,z0) and (x1,y1,z1) - polytope vertex nearest the origin and farthest from the origin ponts respectively
!!v000 - v111: values of a function on a points (example: v101 = v(x1,y0,z1))

!denom = (x1-x0)*(y1-y0)*(z1-z0)


!n000 = (x1-x)*(y1-y)*(z1-z)/denom
!n001 = (x1-x)*(y1-y)*(z-z0)/denom
!n010 = (x1-x)*(y-y0)*(z1-z)/denom
!n011 = (x1-x)*(y-y0)*(z-z0)/denom
!n100 = (x-x0)*(y1-y)*(z1-z)/denom
!n101 = (x-x0)*(y1-y)*(z-z0)/denom
!n110 = (x-x0)*(y-y0)*(z1-z)/denom
!n111 = (x-x0)*(y-y0)*(z-z0)/denom


!v = v000*n000+v001*n001+v010*n010+v011*n011+	&
!	v100*n100+v101*n101+v110*n110+v011*n011

	
!trilin = v

!end function trilin

subroutine stationary_dist(transmat, statvec)
!A proceduree takes stocahstic matrix as input and calculates stationary probabilistic vector
!corresponding to this matrix.
real (kind = 8), dimension(:,:), intent(in) :: transmat
real (kind = 8), dimension(:), intent(out) :: statvec

integer n,m,k !dimensions of the input
real (kind = 8) :: eps=1.0d-6, deltasum = 1.0d0
real (kind = 8), dimension(:), allocatable :: statvec_upd

	n = size(transmat,1)
	m = size(transmat,2)
	k = size(statvec)

	if (n /= m) then
		print *, "Input matrix is not square!"
		read (*,*)
		return
	elseif (n /= k) then
		print *, "Dimensions of input matrix and output vector should correspond!"
		read (*,*)
		return
	endif

	allocate(statvec_upd(k))

	statvec = 1.0d0/k !initialize

	do while (deltasum>eps)
		statvec_upd = matmul(statvec,transmat)
		deltasum = sum(abs(statvec-statvec_upd))
		statvec = statvec_upd
	enddo

	deallocate(statvec_upd)

end subroutine stationary_dist

subroutine cumsum(invec,outvec)
!A function that calculate cumulative sums of consecutive elements of a given vector invec
!and stores them into outvec

real (kind = 8), intent(in), dimension(:) :: invec
real (kind = 8), intent(out), dimension(:) :: outvec

integer n, i

n = size(invec)

outvec(1) = invec(1)

do i = 2,n
	outvec(i) = outvec(i-1)+invec(i)
enddo

end subroutine cumsum

end module procedures
