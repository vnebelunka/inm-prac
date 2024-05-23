! Created by nebil on 25.12.22.

program lapack
    implicit none

    double precision, dimension (:, :), allocatable :: A
    double precision, dimension(:), allocatable :: D, E, tau
    integer :: i, j, n, info, LDA
    character :: UPLO = 'U'
    real :: start, finish

    read(*, *) n
    LDA = n
    allocate (A(n, n), D(n), E(n-1), tau(n-1))
    call random_number(A)
    do i = 1, n
        do j = i, n
            A(i, j) = A(j, i)
        end do
    end do

    call cpu_time(start)
    call dsytd2(UPLO, n, A, LDA, D, E, tau, info)
    call cpu_time(finish)

    print '("Time = ",f6.3," seconds.")',finish-start

    deallocate(A, D, E, tau)
end program lapack