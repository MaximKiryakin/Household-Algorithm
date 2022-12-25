program worst_program
    implicit None

    integer :: i, j, INFO, LDA, LWORK, N

    integer, parameter :: matrix_size = 4096
    double precision, dimension(matrix_size , matrix_size) :: matrix

    double precision, dimension(matrix_size) :: D
    double precision, dimension(matrix_size - 1) :: E, TAU

    double precision, dimension(1) :: tmp
    double precision, allocatable :: WORK(:)
    real :: start, finish

    do i = 1, matrix_size
        do j = 1, matrix_size
            if (i >= j) then
                call random_number(matrix(i, j))
                matrix(j, i) = matrix(i, j)
            end if
        end do
    end do

    call cpu_time(start)

    LWORK = -1

    call DSYTRD('U', matrix_size, matrix, matrix_size, D, E, TAU, tmp, LWORK,  INFO)

    LWORK = tmp(1)

    allocate (WORK(LWORK))

    call DSYTRD('U', matrix_size, matrix, matrix_size, D, E, TAU, WORK, LWORK,  INFO)
    
    call cpu_time(finish)
    print '("Time = ",f6.3," seconds.")',finish-start

end program worst_program
