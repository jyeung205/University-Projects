!M3C 2018 Homework 3 01203261 Jimmy Yeung
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

  nbinv = 1.d0/nb
  allocate(r(n,n,m))
  !---finished Problem setup---

  !----Time marching----
  do i1=1,nt

    call random_number(r) !Random numbers used to update s every time step

    !Set up coefficients for fitness calculation in matrix, a
    a = 1
    where(s==0)
      a=tr_b
    end where

    !create s2 by adding boundary of zeros to s
    s2(2:n+1,2:n+1,:) = s

    !Count number of C neighbors for each point
    nc = s2(1:n,1:n,:) + s2(1:n,2:n+1,:) + s2(1:n,3:n+2,:) + &
         s2(2:n+1,1:n,:)                  + s2(2:n+1,3:n+2,:) + &
         s2(3:n+2,1:n,:)   + s2(3:n+2,2:n+1,:)   + s2(3:n+2,3:n+2,:)

    !Calculate fitness matrix, f----
    f = nc*a
    where(s==0)
      f = f + (nb-nc)*tr_e
    end where
    f = f*nbinv

    !Calculate probability matrix, p----
    f2(2:n+1,2:n+1,:) = f
    f2s2 = f2*s2

    !Total fitness of cooperators in community
    p = f2s2(1:n,1:n,:) + f2s2(1:n,2:n+1,:) + f2s2(1:n,3:n+2,:) + &
           f2s2(2:n+1,1:n,:) + f2s2(2:n+1,2:n+1,:)  + f2s2(2:n+1,3:n+2,:) + &
          f2s2(3:n+2,1:n,:)   + f2s2(3:n+2,2:n+1,:)   + f2s2(3:n+2,3:n+2,:)

    !Total fitness of all members of community
    pden = f2(1:n,1:n,:) + f2(1:n,2:n+1,:) + f2(1:n,3:n+2,:) + &
           f2(2:n+1,1:n,:) + f2(2:n+1,2:n+1,:)  + f2(2:n+1,3:n+2,:) + &
          f2(3:n+2,1:n,:)   + f2(3:n+2,2:n+1,:)   + f2(3:n+2,3:n+2,:)

    p = (p/pden)*tr_g + 0.5d0*(1.d0-tr_g) !probability matrix

    !Set new affiliations based on probability matrix and random numbers stored in R
    s = 0
    where (R<=p)
        s = 1
    end where

    fc_ave(i1+1) = sum(s)*(n2inv/m)
  end do
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
  !My approach to parallisation was to split the m trials over a number of threads.
  !First, I changed the vectorised code so that it calculates the m trials using a do loop.
  !Then, I used OMP parallel do to split up the do loop iterations over a number of threads - set using omp_set_num_threads.
  !It was not possible to parallelise the loop from 1 to nt because the S calculated in each year depends on the previous year.
  !So I parallelised the loop from 1 to m which is inside the loop from 1 to nt.
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
  integer:: i2,j2

  !initial condition and r allocation (does not need to be parallelized)
  s=1
  j1 = (n+1)/2
  s(j1,j1,:) = 0
  allocate(r(n,n,m))
  !------------------
  n2inv = 1.d0/dble(n*n)
  fc_ave(1) = sum(s)*(n2inv/m)

  !$ call omp_set_num_threads(numthreads)

  !$OMP parallel do
  do i1 = 1, m
    s2(:,:,i1) = 0
    f2(:,:,i1) = 0.d0
    nb(:,:,i1) = 8
    nb(1,2:n-1,i1) = 5
    nb(n,2:n-1,i1) = 5
    nb(2:n-1,1,i1) = 5
    nb(2:n-1,n,i1) = 5
    nb(1,1,i1) = 3
    nb(1,n,i1) = 3
    nb(n,1,i1) = 3
    nb(n,n,i1) = 3
    nbinv(:,:,i1) = 1.d0/nb(:,:,i1)
  end do
  !$OMP end parallel do

  do i2=1,nt

    !$OMP parallel do

    do j2 = 1,m

      call random_number(r(:,:,j2))!Random numbers used to update s every time step

      a(:,:,j2) = 1
      where(s(:,:,j2)==0)
        a(:,:,j2)=tr_b
      end where

      !create s2 by adding boundary of zeros to s
      s2(2:n+1,2:n+1,j2) = s(:,:,j2)

      !Count number of C neighbors for each point
      nc(:,:,j2) = s2(1:n,1:n,j2) + s2(1:n,2:n+1,j2) + s2(1:n,3:n+2,j2) + &
                s2(2:n+1,1:n,j2)                  + s2(2:n+1,3:n+2,j2) + &
                s2(3:n+2,1:n,j2)   + s2(3:n+2,2:n+1,j2)   + s2(3:n+2,3:n+2,j2)

      f(:,:,j2) = nc(:,:,j2)*a(:,:,j2)

      where(s(:,:,j2)==0)
        f(:,:,j2) = f(:,:,j2) + (nb(:,:,j2)-nc(:,:,j2))*tr_e
      end where
      f(:,:,j2) = f(:,:,j2)*nbinv(:,:,j2)

      !Calculate probability matrix, p----
      f2(2:n+1,2:n+1,j2) = f(:,:,j2)
      f2s2(:,:,j2) = f2(:,:,j2)*s2(:,:,j2)

      !Total fitness of cooperators in community
      p(:,:,j2) = f2s2(1:n,1:n,j2) + f2s2(1:n,2:n+1,j2) + f2s2(1:n,3:n+2,j2) + &
                  f2s2(2:n+1,1:n,j2) + f2s2(2:n+1,2:n+1,j2)  + f2s2(2:n+1,3:n+2,j2) + &
                  f2s2(3:n+2,1:n,j2)   + f2s2(3:n+2,2:n+1,j2)   + f2s2(3:n+2,3:n+2,j2)

      !Total fitness of all members of community
      pden(:,:,j2) = f2(1:n,1:n,j2) + f2(1:n,2:n+1,j2) + f2(1:n,3:n+2,j2) + &
                     f2(2:n+1,1:n,j2) + f2(2:n+1,2:n+1,j2)  + f2(2:n+1,3:n+2,j2) + &
                     f2(3:n+2,1:n,j2)   + f2(3:n+2,2:n+1,j2)   + f2(3:n+2,3:n+2,j2)

      p(:,:,j2) = (p(:,:,j2)/pden(:,:,j2))*tr_g + 0.5d0*(1.d0-tr_g) !probability matrix

      s(:,:,j2) = 0
      where (R(:,:,j2)<=p(:,:,j2))
         s(:,:,j2) = 1
      end where

    end do

    !$OMP end parallel do
    fc_ave(i2+1) = sum(s)*(n2inv/m)

  end do

  deallocate(r)
end subroutine simulate2_omp

end module tribes
