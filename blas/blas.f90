PROGRAM blas
    implicit none

    DOUBLE PRECISION, dimension (:, :), allocatable :: A, D
    DOUBLE PRECISION, dimension (:), allocatable :: coss, sins
    integer :: n, i, j, cur_rots_size
    real :: start, finish

    cur_rots_size = 1
    read(*, *) n

    allocate (A(n, n), D(n, n), coss(n * n), sins(n * n))

    call random_number(A)

    do i = 1, n
        do j = i, n
            A(i, j) = A(j, i)
        end do
    end do

    !write(*, *) A

    call cpu_time(start)
    call tridiagonalization(A, D, n, coss, sins, cur_rots_size)
    call cpu_time(finish)
    if (n < 5) then
        do i = 1, n
            print *, D(i, :)
        end do
    endif
    print '("Time = ",f6.3," seconds.")',finish-start

    deallocate (A, D, coss, sins)

END PROGRAM blas

subroutine applyRotationCR(D, n, col, coss, sins, rots_size)
    implicit none
    integer, INTENT(IN) :: n, col
    DOUBLE PRECISION, INTENT(OUT), dimension(n, n) :: D
    DOUBLE PRECISION, INTENT(OUT), dimension(n * n) :: coss, sins
    integer, INTENT(OUT) :: rots_size

    integer :: cur_rot, rots_start, j, i
    DOUBLE PRECISION :: modd, c, s, aii, aij

    rots_start = rots_size

    ! Apply rotation to cols: A -> AQ^T, accamulate sin, cos
    do i  = col + 2, n
        modd = sqrt(D(col+1, col)**2 + D(i,col)**2)
        s = -D(i, col) / modd
        c = D(col, col+1) / modd
        call drot(n, D(:, col+1), 1, D(:, i), 1, c, -s)
        coss(rots_size) = c
        sins(rots_size) = s
        rots_size = rots_size + 1
    end do

    ! Apply rotation to cols AQ^T -> QAQ^T with accamulated sin, cos
    do i = col, n
        cur_rot = rots_start
        do j = col + 2, n
            c = coss(cur_rot)
            s = sins(cur_rot)
            aii = c * D(col + 1, i) - s * D(j, i)
            aij = s * D(col + 1, i) + c * D(j, i)
            D(col + 1, i) = aii
            D(j, i) = aij
            cur_rot = cur_rot + 1
        end do
    end do
end subroutine applyRotationCR

subroutine tridiagonalization(A, D, n, coss, sins, cur_rots_size)
    implicit none
    integer, intent(in) :: n
    DOUBLE PRECISION, intent(in), dimension(n, n) :: A
    DOUBLE PRECISION, intent(out), dimension(n, n) :: D
    DOUBLE PRECISION, intent(out), dimension(n * n) :: coss, sins
    integer, intent(out) :: cur_rots_size

    integer :: i

    D = A
    do i = 1, n-2
        call applyRotationCR(D, n, i, coss, sins, cur_rots_size)
    end do
end subroutine tridiagonalization