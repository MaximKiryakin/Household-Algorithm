program program_1
    implicit none

    integer, parameter :: matrix_size = 4
    double precision, dimension(matrix_size , matrix_size) :: matrix
    real :: start, finish
    integer :: i, j
    !read (*,*) (matrix(i), i = 1, matrix_size **2)
    !matrix = reshape((/ 3, 6, 2, 4, 6, 1, 6, 1, 2, 6, 5, 5, 4, 1, 5, 4 /), shape(matrix))
    do i = 1, matrix_size
        do j = 1, matrix_size
            if (i >= j) then
                call random_number(matrix(i, j))
                matrix(j, i) = matrix(i, j)
            end if
        end do
    end do

    call cpu_time(start)
    call hausholders_transformation(matrix)
    call print_matrix(matrix)
    call cpu_time(finish)
    print '("Time = ",f6.3," seconds.")',finish-start


contains


    !процедура печати матрицы


    subroutine print_matrix(matrix)
        implicit none
        double precision, dimension(matrix_size,matrix_size) :: matrix
        integer :: i, j
        do i = 1, matrix_size
            do j = 1, matrix_size
                if (j < matrix_size) then
                    write(*, 101) matrix(i, j)
                    101 format(F10.5, ' ', $)
                else
                    write(*, 102) matrix(i, j)
                    102 format(F10.5, ' ')
                end if

                !if (j < matrix_size) then
                !    write(*, 101) matrix((i-1) * matrix_size + j)
                !    101 format(F10.5, ' ', $)
                !else
                !    write(*, 102) matrix((i-1) * matrix_size + j)
                !    102 format(F10.5, ' ')
                !end if
            end do
        end do
	end subroutine print_matrix



    ! процедура выполняет перемножение матриц PAP
    subroutine pap_transformation(a, b, matrix, digit)
        implicit none
        integer, intent(in) :: digit

        double precision, dimension(:), intent(in) :: a, b
        double precision, dimension(matrix_size, matrix_size), intent(inout) :: matrix

        integer :: i, j, size1
        size1 = size(a)

        do i = 1, size1
            do j = 1, size1
                !matrix((i - 1)*size1 + j) = matrix((i - 1)*size1 + j) - a(i) * b(j) * digit
                !matrix((i - 1)*size1 + j) = matrix((i - 1)*size1 + j) - a(j) * b(i) * digit
                matrix(i, j) = matrix(i, j) - a(j) * b(i) * digit
                matrix(i, j) = matrix(i, j) - a(i) * b(j) * digit
            end do
        end do

    end subroutine pap_transformation

    !скалярное умножение векторов
    subroutine dot(a, b, ans)
        implicit none

        double precision :: ans
        double precision, dimension(:) :: a, b
        !integer :: i

        ans = 0
        do i = 1, size(a)
            ans = ans + a(i) * b(i)
        end do
	!ans = ddot(size(a), a, 1, b, 1)
    end subroutine dot

    !умножение матрицы на вектор
    subroutine dot_matrix_vector(matrix, b, ans)
        implicit none

        double precision, dimension(matrix_size, matrix_size) :: matrix
        double precision, dimension(:) :: b
        double precision, dimension(:) :: ans

        integer :: k, i

        do k = 1, size(ans)
            ans(k) = 0
        end do

        !call dsymv('U',  matrix_size, 1, matrix, matrix_size, b, 1, 0, ans, 1)
        do k = 1, size(ans)
            ans(k) = 0
            do i = 1, size(ans)
                ans(k) = ans(k) + matrix(k, i) * b(i)
            end do
        end do

    end subroutine dot_matrix_vector



    !преобразование матрицы по алгоритму хаусхолдера
    subroutine hausholders_transformation(matrix)
        double precision, dimension(matrix_size, matrix_size), intent(inout) :: matrix

        integer :: size2, i, j
        double precision :: norm
        double precision :: s, r, k
        double precision :: c
        double precision :: yyy
        double precision, dimension(matrix_size) :: w, v, q

        size2 = matrix_size


        do i = 1, size2 -2

            norm = 0
            do j = i + 1, size2
                norm = norm + matrix(j, i) * matrix(j, i)
            end do


            norm = sqrt(norm)


            yyy = 1
            if (matrix(i + 1, i) > 0) then
                s =  norm
            else
                s = (-1)* norm
            end if
            r =sqrt(2 * matrix(i + 1, i) * s + 2 * s * s)

            do k = 1, size2
                w(k) = 0
            end do

            w(i + 1) = (matrix(i + 1, i) + s) * (1 / r)

            do j = i + 2, size2
                if (j >= i + 2) then
                    w(j) = matrix(j, i) * (1/ r)
                end if

                if (j <= i) then
                    w(j) = 0
                end if
            end do

            do k = 1, size2
                v(k) = 0
            end do

            call dot_matrix_vector(matrix, w, v)


            c = 0
            call dot(v, w, c)


            do j = 1, size2
                q(j) = v(j) - c * w(j)
            end do

            call pap_transformation(w, q, matrix, 2)
        end do
    end subroutine hausholders_transformation



    !6, 5, 3, 4, 5, 2, 3, 4, 3, 3, 4, 3, 4, 4, 3, 5
    !3, 6, 2, 4, 6, 1, 6, 1, 2, 6, 5, 5, 4, 1, 5, 4



end program program_1

