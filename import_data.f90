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
