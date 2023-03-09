! Created by nebil on 09.03.23.

program main
    implicit none
    ! RCS for input matrix
    integer :: nmax, nzmax , n, nnz, ierr, iounit, i
    real*8, dimension (:), allocatable :: a, rhs
    integer, dimension (:), allocatable :: ia, ja

    ! ILU(k) space
    real*8, dimension (:), allocatable :: alu, w
    integer, dimension (:), allocatable :: jlu, ju, levs, jw
    integer :: lfil, iwk

    !PGMRES space
    real*8, dimension (:), allocatable :: sol, vv
    real*8 :: eps
    integer :: maxits, im, iout

    ! read matrix
    nmax = 5000
    nzmax = 53000
    iounit = 1

    allocate(ia(nmax+1), ja(nzmax), a(nzmax))
    open(iounit, file="../data/1.dat", status ="old")

    call readsk(nmax, nzmax, n, nnz, a, ja, ia, 1, ierr)
    write(*, *) ierr, n

    ! rhs
    allocate(rhs(n))
    do i = 1, n
        rhs(i) = sin(real(i))
    end do

    ! ILU(k)
    iwk = nzmax * 100
    allocate(jw(3 * n), w(n), alu(iwk), jlu(iwk), levs(iwk), ju(n))
    lfil = 2

    call iluk(n,a,ja,ia,lfil,alu,jlu,ju,levs,iwk,w,jw,ierr)

    write(*, *) ierr

    ! PGMRES
    eps = 1e-3 ! во сколько упадёт невязка
    im = 10 ! кол-во крыловских подпро-в (гиперпараметр)
    maxits = n * n ! max итераций в gmres
    iout = 2 ! дескриптор файла для логов GMRES

    open(iout, file="../log.txt", status ="old")
    allocate(sol(n), vv(n * (im + 1))) ! sol: вектор для решения. vv: рабочее пространство

    call pgmres(n, im, rhs, sol, vv, eps, maxits, iout, a, ja, ia, alu, jlu, ju, ierr)

    write(*, *)(ierr)

    !do i = 1, n
    !    write(*, *  )sol(i)
    !end do

    deallocate(ia, ja, a, jw, w, alu, jlu, levs, ju, sol, vv, rhs)

end program main