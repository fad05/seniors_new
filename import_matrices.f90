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
