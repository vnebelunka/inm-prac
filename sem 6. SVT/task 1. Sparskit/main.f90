! Created by nebil on 09.03.23.

program main
    implicit none
    real :: start, finish
    ! RCS for input matrix
    integer :: nmax, nzmax , n, nnz, ierr, iounit, i
    real*8, dimension (:), allocatable :: a, rhs
    integer, dimension (:), allocatable :: ia, ja

    character(len=30) :: fdat

    ! ILU(k) space
    real*8, dimension (:), allocatable :: alu, w
    integer, dimension (:), allocatable :: jlu, ju, levs, jw
    integer :: lfil, iwk

    !PGMRES space
    real*8, dimension (:), allocatable :: sol, vv
    real*8 :: eps
    integer :: maxits, im, iout

    ! read matrix

    call get_command_argument(1,fdat)
    nmax = 264752 ! max rows in all 4 .dat files
    nzmax = 3541545 ! max nz in all 4 .dat files
    iounit = 1

    allocate(ia(nmax+1), ja(nzmax), a(nzmax))
    open(iounit, file=fdat, status ="old")

    call readsk(nmax, nzmax, n, nnz, a, ja, ia, 1, ierr)
    print '("size of matrix: ",i0)',n
    if(ierr /= 0) then
        print '("read problems")'
        stop ierr
    end if

    ! rhs
    allocate(rhs(n))
    do i = 1, n
        rhs(i) = sin(real(i))
    end do

    ! ILU(k)
    iwk = nnz * 100
    allocate(jw(3 * n), w(n), alu(iwk), jlu(iwk), levs(iwk), ju(n))
    lfil = 0

    call cpu_time(start)
    call iluk(n,a,ja,ia,lfil,alu,jlu,ju,levs,iwk,w,jw,ierr)
    call cpu_time(finish)
    print '("Time for initialization = ",f6.3," seconds.")',finish-start
    if(ierr /= 0) then
        print '("iluk problems")'
        stop ierr
    end if

    ! PGMRES
    eps = 1e-8 ! во сколько упадёт невязка
    im = 10 ! кол-во крыловских подпро-в (гиперпараметр)
    maxits = n * n ! max итераций в gmres
    iout = 2 ! дескриптор файла для логов GMRES

    open(iout, file="../log.txt", status ="old")
    allocate(sol(n), vv(n * (im + 1))) ! sol: вектор для решения. vv: рабочее пространство

    call cpu_time(start)
    call pgmres(n, im, rhs, sol, vv, eps, maxits, iout, a, ja, ia, alu, jlu, ju, ierr)
    call cpu_time(finish)
    print '("Time for solving = ",f9.3," seconds.")',finish-start
    if(ierr /= 0) then
        print '("pgmres problems")'
        stop ierr
    end if


    deallocate(ia, ja, a, jw, w, alu, jlu, levs, ju, sol, vv, rhs)

end program main