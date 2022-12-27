!M3C 2018 Homework 3
!This module contains four module variables and two subroutines;
!one of these routines must be developed for this assignment.
!Module variables--
! tr_b, tr_e, tr_g: the parameters b, e, and g in the tribe competition model
! numthreads: The number of threads that should be used in parallel regions within simulate2_omp
!
!Module routines---
! simulate2_f90: Simulate tribal competition over m trials. Return: all s matrices at final time
! and fc at nt+1 times averaged across the m trials.
! simulate2_omp: Same input/output functionality as simulate2.f90 but parallelized with OpenMP

module tribes
  use omp_lib
  implicit none
  integer :: numthreads
  real(kind=8) :: tr_b,tr_e,tr_g
contains


!Simulate m trials of Cooperator vs. Mercenary competition using the parameters, tr_b and tr_e.
!Input:
! n: defines n x n grid of villages
! nt: number of time steps
! m: number of trials
!Output:
! s: status matrix at final time step for all m trials
! fc_ave: fraction of cooperators at each time step (including initial condition)
! averaged across the m trials
subroutine simulate2_f90(n,nt,m,s,fc_ave)
  implicit none
  integer, intent(in) :: n,nt,m
  integer, intent(out), dimension(n,n,m) :: s
  real(kind=8), intent(out), dimension(nt+1) :: fc_ave
  integer :: i1,j1
  real(kind=8) :: n2inv
  integer, dimension(n,n,m) :: nb,nc
  integer, dimension(n+2,n+2,m) :: s2
  real(kind=8), dimension(n,n,m) :: f,p,a,pden,nbinv
  real(kind=8), dimension(n+2,n+2,m) :: f2,f2s2
  real(kind=8), allocatable, dimension(:,:,:) :: r !random numbers

  !---Problem setup----
  !Initialize arrays and problem parameters

  !initial condition
  s=1
  j1 = (n+1)/2
  s(j1,j1,:) = 0

  n2inv = 1.d0/dble(n*n)
  fc_ave(1) = sum(s)*(n2inv/m)

  s2 = 0
  f2 = 0.d0

  !Calculate number of neighbors for each point
  nb = 8
  nb(1,2:n-1,:) = 5
  nb(n,2:n-1,:) = 5
  nb(2:n-1,1,:) = 5
  nb(2:n-1,n,:) = 5
  nb(1,1,:) = 3
  nb(1,n,:) = 3
  nb(n,1,:) = 3
  nb(n,n,:) = 3

  !print *, "nb_orig=", nb

  nbinv = 1.d0/nb
  allocate(r(n,n,m))
  !---finished Problem setup---
  !print *,"1nt=",nt
  !----Time marching----
  do i1=1,nt

    call random_number(r) !Random numbers used to update s every time step

    !Set up coefficients for fitness calculation in matrix, a
    a = 1
    where(s==0)
      a=tr_b
    end where

    !print *, "a =", a

    !create s2 by adding boundary of zeros to s
    s2(2:n+1,2:n+1,:) = s

    !print *, "1s2=", s2

    !Count number of C neighbors for each point
    nc = s2(1:n,1:n,:) + s2(1:n,2:n+1,:) + s2(1:n,3:n+2,:) + &
         s2(2:n+1,1:n,:)                  + s2(2:n+1,3:n+2,:) + &
         s2(3:n+2,1:n,:)   + s2(3:n+2,2:n+1,:)   + s2(3:n+2,3:n+2,:)

    !print *,"1nc=", nc
    !Calculate fitness matrix, f----
    f = nc*a
    !print *,"shape f + (nb-nc)*tr_e=", shape(f + (nb-nc)*tr_e)
    where(s==0)
      f = f + (nb-nc)*tr_e
    end where
    f = f*nbinv
    !print *,"shape f + (nb-nc)*tr_e=", shape(f + (nb-nc)*tr_e)
    !print *,"f1" , f
    !-----------

    !Calculate probability matrix, p----
    f2(2:n+1,2:n+1,:) = f
    f2s2 = f2*s2

    !print *,"1f2s2" , f2s2

    !Total fitness of cooperators in community
    p = f2s2(1:n,1:n,:) + f2s2(1:n,2:n+1,:) + f2s2(1:n,3:n+2,:) + &
           f2s2(2:n+1,1:n,:) + f2s2(2:n+1,2:n+1,:)  + f2s2(2:n+1,3:n+2,:) + &
          f2s2(3:n+2,1:n,:)   + f2s2(3:n+2,2:n+1,:)   + f2s2(3:n+2,3:n+2,:)

          !print *,"1p=",p

    !Total fitness of all members of community
    pden = f2(1:n,1:n,:) + f2(1:n,2:n+1,:) + f2(1:n,3:n+2,:) + &
           f2(2:n+1,1:n,:) + f2(2:n+1,2:n+1,:)  + f2(2:n+1,3:n+2,:) + &
          f2(3:n+2,1:n,:)   + f2(3:n+2,2:n+1,:)   + f2(3:n+2,3:n+2,:)

          !print *,"pden1=",pden


    p = (p/pden)*tr_g + 0.5d0*(1.d0-tr_g) !probability matrix
    !print*,"p1=",p
    !----------

    !Set new affiliations based on probability matrix and random numbers stored in R
    s = 0
    where (R<=p)
        s = 1
    end where

    !print*,"s1=",s

    fc_ave(i1+1) = sum(s)*(n2inv/m)


  end do
  print*,"fc_ave1=",fc_ave
end subroutine simulate2_f90


!Simulate m trials of Cooperator vs. Mercenary competition using the parameters, tr_b and tr_e.
!Same functionality as simulate2_f90, but parallelized with OpenMP
!Parallel regions should use numthreads threads.
!Input:
! n: defines n x n grid of villages
! nt: number of time steps
! m: number of trials
!Output:
! s: status matrix at final time step for all m trials
! fc_ave: fraction of cooperators at each time step (including initial condition)
! averaged across the m trials
subroutine simulate2_omp(n,nt,m,s,fc_ave)
  implicit none
  integer, intent(in) :: n,nt,m
  integer, intent(out), dimension(n,n,m) :: s
  real(kind=8), intent(out), dimension(nt+1) :: fc_ave
  integer :: i1,j1
  real(kind=8), allocatable, dimension(:,:,:) :: r !random numbers

  !Add further variables as needed
  real(kind=8) :: n2inv
  integer, dimension(n,n,m) :: nb,nc
  integer, dimension(n+2,n+2,m) :: s2
  real(kind=8), dimension(n,n,m) :: f,p,a,pden,nbinv
  real(kind=8), dimension(n+2,n+2,m) :: f2,f2s2
  integer:: threadID, numthreads, m1,j2

  !initial condition and r allocation (does not need to be parallelized)
  s=1
  j1 = (n+1)/2
  s(j1,j1,:) = 0
  allocate(r(n,n,m))
  !------------------
  n2inv = 1.d0/dble(n*n)
  fc_ave(1) = sum(s)*(n2inv/m)
  s2 = 0
  f2 = 0.d0
  nb = 8

  do i1 = 1, m
    nb(1,2:n-1,i1) = 5
    nb(n,2:n-1,i1) = 5
    nb(2:n-1,1,i1) = 5
    nb(2:n-1,n,i1) = 5
    nb(1,1,i1) = 3
    nb(1,n,i1) = 3
    nb(n,1,i1) = 3
    nb(n,n,i1) = 3
  end do

  nbinv = 1.d0/nb
  !print *, "nb_omp=", nb
  !print *,"2nt=",nt

  do m1=1,nt

    call random_number(r)!Random numbers used to update s every time step

    a = 1
    where(s==0)
      a=tr_b
    end where
    !print*,"a=",a

    do j2 = 1,m

      !create s2 by adding boundary of zeros to s
      s2(2:n+1,2:n+1,j2) = s(:,:,j2)

      !print *, "2s2=", s2

      !Count number of C neighbors for each point
      nc(:,:,j2) = s2(1:n,1:n,j2) + s2(1:n,2:n+1,j2) + s2(1:n,3:n+2,j2) + &
                s2(2:n+1,1:n,j2)                  + s2(2:n+1,3:n+2,j2) + &
                s2(3:n+2,1:n,j2)   + s2(3:n+2,2:n+1,j2)   + s2(3:n+2,3:n+2,j2)

      !print *,"2nc=",nc

      f(:,:,j2) = nc(:,:,j2)*a(:,:,j2)

      !print *,"shape f(:,:,j2)=", shape(f(:,:,j2))
      !print *,"shape nb(:,:,j2)=", shape(nb(:,:,j2))
      !print *,"shape nc(:,:,j2)=", shape(nc(:,:,j2))
      !print *,"shape tr_e=", shape(tr_e)
      !print *, "tr_e=", tr_e
      where(s(:,:,j2)==0)
        f(:,:,j2) = f(:,:,j2) + (nb(:,:,j2)-nc(:,:,j2))*tr_e
      end where
      f(:,:,j2) = f(:,:,j2)*nbinv(:,:,j2)

      !print *,"f2" , f

      !Calculate probability matrix, p----
      f2(2:n+1,2:n+1,j2) = f(:,:,j2)
      f2s2(:,:,j2) = f2(:,:,j2)*s2(:,:,j2)
      !print *,"2f2s2" , f2s2

      !Total fitness of cooperators in community
      p(:,:,j2) = f2s2(1:n,1:n,j2) + f2s2(1:n,2:n+1,j2) + f2s2(1:n,3:n+2,j2) + &
                  f2s2(2:n+1,1:n,j2) + f2s2(2:n+1,2:n+1,j2)  + f2s2(2:n+1,3:n+2,j2) + &
                  f2s2(3:n+2,1:n,j2)   + f2s2(3:n+2,2:n+1,j2)   + f2s2(3:n+2,3:n+2,j2)

      !print *,"2p=",p

      !Total fitness of all members of community
      pden(:,:,j2) = f2(1:n,1:n,j2) + f2(1:n,2:n+1,j2) + f2(1:n,3:n+2,j2) + &
                     f2(2:n+1,1:n,j2) + f2(2:n+1,2:n+1,j2)  + f2(2:n+1,3:n+2,j2) + &
                     f2(3:n+2,1:n,j2)   + f2(3:n+2,2:n+1,j2)   + f2(3:n+2,3:n+2,j2)

      !print *,"pden2=",pden

      p(:,:,j2) = (p(:,:,j2)/pden(:,:,j2))*tr_g + 0.5d0*(1.d0-tr_g) !probability matrix
      !print*,"p2=",p

      s(:,:,j2) = 0
      where (R(:,:,j2)<=p(:,:,j2))
         s(:,:,j2) = 1
      end where
      !print*,"s2=",s
    end do
    fc_ave(m1+1) = sum(s)*(n2inv/m)

  end do

  print*,"fc_ave2=",fc_ave

  deallocate(r)
end subroutine simulate2_omp


end module tribes
