program worst_program
    implicit None

    integer :: i, INFO, LDA, LWORK, N

    integer, parameter :: matrix_size = 4
    double precision, dimension(matrix_size , matrix_size) :: matrix

    double precision, dimension(matrix_size) :: D
    double precision, dimension(matrix_size - 1) :: E, TAU

    double precision, dimension(1) :: tmp
    double precision, allocatable :: WORK(:)

    matrix = reshape((/ 3, 6, 2, 4, 6, 1, 6, 1, 2, 6, 5, 5, 4, 1, 5, 4 /), shape(matrix))

    ! сначала передаю LWORK какой попало
    LWORK = 1

    WORK = -1
    tmp(1) = 10
    call SSYTRD('U', matrix_size, matrix, matrix_size, D, E, TAU, tmp, LWORK,  INFO)

    LWORK = tmp(1)
    write(*,*) LWORK

    allocate (WORK(LWORK))
    ! call SSYTRD('U', matrix_size, matrix, matrix_size, D, E, TAU, WORK, LWORK,  INFO)

    print *, D
    print *, E


end program worst_program
