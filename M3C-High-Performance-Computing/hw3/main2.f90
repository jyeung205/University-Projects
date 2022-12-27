! This is a main program which can be used with the tribes module
! if and as you wish as you develop your codes.
! It reads the problem parameters from a text file, data.in, which must be created.
!
! The subroutine simulate2_f90 is called, and output is written to the text files,
! s.txt and fc_ave.txt which can be read in Python using np.loadtxt (see below)
!
! You should not submit this code with your assignment.
! To compile: gfortran -fopenmp -O3 -o main.exe hw3_dev.f90 hw3_main.f90
! TO RUN: gfortran -fopenmp -O3 -o main.exe hw3_dev.f90 hw3_main.f90 && ./main.exe
! gfortran -fopenmp -O3 -o main.exe hw3_dev2.f90 hw3_main.f90 && ./main.exe
program hw3_main  !FOR TESTING PURPOSES WITH RANDOM NUMBERS
  use tribes
  implicit none
  integer :: n,nt,i1,j1,m1!m!, threadID
  integer, allocatable, dimension(:) :: m
  real(kind=8), allocatable, dimension(:) :: fc_ave
  integer, allocatable, dimension(:,:,:) :: s
  real(kind=8), dimension(20) :: timing
  real(kind=8) :: start,finish,start2,finish2,timefin1, timefin2, timefin
  real(kind=8), allocatable, dimension(:,:,:,:) :: r !random numbers


! DEBUGGING TODAY
  !integer, dimension(21,21,1000) :: nb

  !Read in problem parameters from text file, data.in
  open(unit=11,file='data.in')
  read(11,*) n !n x n villages
  read(11,*) nt !number of time steps
  read(11,*) tr_b !model parameters
  read(11,*) tr_e
  read(11,*) tr_g
  !read(11,*) m !number of trials
  read(11,*) numthreads !not used below
  close(11)

  !Change This for testing
  m1 = 3 !Number of m's to test
  allocate(m(m1))
  m = [500, 750, 850]    !Add more if m1 changes



  call omp_set_num_threads(numthreads)
!UNCOMMENT FOR TIME TESTING
  !$OMP parallel
  numThreads = omp_get_num_threads()
  !$OMP end parallel
  print *, 'numThreads is=', numThreads
  allocate(fc_ave(nt+1))
do j1 = 1,m1
     allocate(s(n,n,m(j1)), r(n,n,m(j1),nt))
     call random_number(r)
        do i1 = 1,10
            call cpu_time(start)
            call simulate2_f90(n,nt,m(j1),r,s,fc_ave)
            call cpu_time(finish)
            timing(i1) = finish-start
        end do
        timefin1 = sum(timing)/10
        !print '("Average time to run non-parallelized version 10 times = ",f20.15," seconds.")',timefin1
        do i1 = 1,10
            call cpu_time(start2)
            call simulate2_omp(n,nt,m(j1),r,s,fc_ave)
            call cpu_time(finish2)
            timing(i1) = finish2-start2
        end do
        timefin2 = sum(timing)/10
        !print '("Average time to run parallelized version 10 times = ",f20.15," seconds.")',timefin2
    timefin = (timefin1-timefin2)/timefin1 * 100
    print*, 'M = ', m(j1)
    print '("Parallelized is ",f10.4,"% faster")', timefin
    deallocate(s,r)
end do
deallocate(fc_ave,m)
allocate(m(1))

!!!!!!!!!
m = 500                 !Change for individual Game
!!!!!!!!!

allocate(fc_ave(nt+1),s(n,n,m(1)),r(n,n,m(1),nt))
call random_number(r)
call simulate2_f90(n,nt,m(1),r,s,fc_ave)
print*, fc_ave((Nt-1):)
deallocate(fc_ave)
allocate(fc_ave(nt+1))
deallocate(s)
allocate(s(n,n,m(1)))
call simulate2_omp(n,nt,m(1),r,s,fc_ave)
print*, fc_ave((Nt-1):)
deallocate(fc_ave,m,s,r)
print*, 'All calculations finished'
end program hw3_main
