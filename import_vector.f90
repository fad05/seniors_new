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
