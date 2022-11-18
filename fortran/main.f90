PROGRAM main
    implicit none

    DOUBLE PRECISION, dimension (:, :), allocatable :: A, D
    DOUBLE PRECISION, dimension (:), allocatable :: coss, sins
    integer :: n, i, j, cur_rots_size

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

    call tridiagonalization(A, D, n, coss, sins, cur_rots_size)

    !write(*, *) D

    deallocate (A, D, coss, sins)

END PROGRAM main


subroutine rotationRow(A,n, i, j, c, s, start)
    implicit none
    integer :: k
    integer, INTENT(IN) :: n, i, j, start
    DOUBLE PRECISION, INTENT(IN) :: c, s
    DOUBLE PRECISION, INTENT(OUT), dimension(n, n) :: A
    DOUBLE PRECISION :: aik, ajk
    do k = start, n
        aik = c * A(i, k) - s * A(j, k)
        ajk = s * A(i, k) + c * A(j, k)
        A(i, k) = aik
        A(j, k) = ajk
    end do
end subroutine rotationRow

subroutine applyRotation(D, n, row, coss, sins, rots_size)
    implicit none
    integer, INTENT(IN) :: n, row
    DOUBLE PRECISION, INTENT(OUT), dimension(n, n) :: D
    DOUBLE PRECISION, INTENT(OUT), dimension(n * n) :: coss, sins
    integer, INTENT(OUT) :: rots_size

    integer :: cur_rot, rots_start, j, i
    DOUBLE PRECISION :: modd, c, s, aii, aij

    rots_start = rots_size

    ! Apply rotation to rows: A -> QA
    do j  = row + 2, n
        modd = sqrt(D(row+1, row)**2 + D(j,row)**2)
        s = -D(j, row) / modd
        c = D(row + 1, row) / modd
        call rotationRow(D,n, row + 1 ,j, c, s, row)
        coss(rots_size) = c
        sins(rots_size) = s
        rots_size = rots_size + 1
    end do

    ! Apply rotation to cols QA -> QAQ^T
    do i = row, n
        cur_rot = rots_start
        do j = row + 2, n
            c = coss(cur_rot)
            s = sins(cur_rot)
            aii = c * D(i, row + 1) - s * D(i, j)
            aij = s * D(i, row + 1) + c * D(i, j)
            D(i, row + 1) = aii
            D(i, j) = aij
            cur_rot = cur_rot + 1
        end do
    end do

end subroutine applyRotation

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
        call applyRotation(D, n, i, coss, sins, cur_rots_size)
    end do
end subroutine tridiagonalization