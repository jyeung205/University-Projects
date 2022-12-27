!M3C 2018 Homework 2 Jimmy Yeung 01203261
!This module contains two module variables and four subroutines;
!two of these routines must be developed for this assignment.
!Module variables--
! nm_x: training images, typically n x d with n=784 and d<=60000
! nm_y: labels for training images, d-element array containing 0s and 1s
! corresponding to images of even and odd integers, respectively.
!
!Module routines---
! data_init: allocate nm_x and nm_y using input variables n and d. Used by sgd, may be used elsewhere if needed
! sgd: Use simple stochastic descent algorithm to iteratively find fitting parameters using either snmodel (when m=0) or
! nnmodel (when m>1)
! snmodel: compute cost function and gradient using single neuron model (SNM) with nm_x and nm_y, and
! with fitting parameters provided as input
! nnmodel: compute cost function and gradient using neural network model (NNM) with 1 output neuron and
! m neurons in the inner layer, and with nm_x and nm_y, and with fitting parameters provided as input

module nmodel
  implicit none
  real(kind=8), allocatable, dimension(:,:) :: nm_x
  integer, allocatable, dimension(:) :: nm_y

contains

!---allocate nm_x and nm_y deallocating first if needed (used by sgd)---
subroutine data_init(n,d)
  implicit none
  integer, intent(in) :: n,d
  if (allocated(nm_x)) deallocate(nm_x)
  if (allocated(nm_y)) deallocate(nm_y)
  allocate(nm_x(n,d),nm_y(d))
end subroutine data_init


subroutine snforward(fvec,n,d,a)
  implicit none
  integer, intent(in) :: n,d !training data sizes
  real(kind=8), dimension(n+1), intent(in) :: fvec !fitting parameters
  real(kind=8), dimension(d) :: z
  real(kind=8), dimension(d), intent(out) :: a
  integer :: k1

  do k1 = 1,d
    z(k1) = sum(fvec(1:n) * nm_x(:,k1)) + fvec(n+1)
    a(k1) = 1.d0/(1.d0 + exp(-z(k1)))
  end do

end subroutine snforward


!Compute cost function and its gradient for single neuron model
!for d images (in nm_x) and d labels (in nm_y) along with the
!fitting parameters provided as input in fvec.
!The weight vector, w, corresponds to fvec(1:n) and
!the bias, b, is stored in fvec(n+1)
!Similarly, the elements of dc/dw should be stored in cgrad(1:n)
!and dc/db should be stored in cgrad(n+1)
!Note: nm_x and nm_y must be allocated and set before calling this subroutine.
subroutine snmodel(fvec,n,d,c,cgrad)
  implicit none
  integer, intent(in) :: n,d !training data sizes
  real(kind=8), dimension(n+1), intent(in) :: fvec !fitting parameters
  real(kind=8), intent(out) :: c !cost
  real(kind=8), dimension(n+1), intent(out) :: cgrad !gradient of cost

  !Declare other variables as needed
  real(kind=8), dimension(d) :: a
  integer :: l1

  !Add code to compute c and cgrad
  call snforward(fvec,n,d,a)

  c = (1.d0/(2.d0*real(d))) * sum((a - nm_y)**2.d0)

  cgrad(n+1) = 1/real(d) * sum( (a - nm_y) * (a*(1-a)) )

  do l1 = 1, n
    cgrad(l1) = 1/real(d) * sum( (a - nm_y) * nm_x(l1,:) * a * (1-a) )
  end do

end subroutine snmodel


subroutine nnforward(fvec,n,m,d,a,a_outer)
  implicit none
  real(kind=8), dimension(m*(n+2)+1), intent(in) :: fvec
  integer, intent(in) :: n,m,d
  real(kind=8), dimension(m,d), intent(out) :: a
  real(kind=8), dimension(m,d) :: z
  real(kind=8), dimension(d), intent(out) :: a_outer
  real(kind=8), dimension(d) :: z_outer
  integer :: k1, j1, i1

  real(kind=8), dimension(m,n) :: w_inner
  real(kind=8), dimension(m) :: b_inner,w_outer
  real(kind=8) :: b_outer

  !unpack fitting parameters (use if needed)
  do i1=1,n
    j1 = (i1-1)*m+1
    w_inner(:,i1) = fvec(j1:j1+m-1) !inner layer weight matrix
  end do
  b_inner = fvec(n*m+1:n*m+m) !inner layer bias vector
  w_outer = fvec(n*m+m+1:n*m+2*m) !output layer weight vector
  b_outer  = fvec(n*m+2*m+1) !output layer bias

  do k1 = 1, d
    do j1= 1, m
     z(j1,k1) = sum(w_inner(j1,:)*nm_x(:,k1)) + b_inner(j1)
     a(j1,k1) = 1.d0/(1.d0 + exp(-z(j1,k1)))
    end do
  end do

  do k1 = 1, d
   z_outer(k1) = sum(w_outer * a(:,k1)) + b_outer
   a_outer(k1) = 1.d0/(1.d0 + exp(-z_outer(k1)))
  end do

end subroutine nnforward


!!Compute cost function and its gradient for neural network model
!for d images (in nm_x) and d labels (in nm_y) along with the
!fitting parameters provided as input in fvec. The network consists of
!an inner layer with m neurons and an output layer with a single neuron.
!fvec contains the elements of dw_inner, b_inner, w_outer, and b_outer
! Code has been provided below to "unpack" fvec
!The elements of dc/dw_inner,dc/db_inner, dc/dw_outer,dc/db_outer should be stored in cgrad
!and should be "packed" in the same order that fvec was unpacked.
!Note: nm_x and nm_y must be allocated and set before calling this subroutine.
subroutine nnmodel(fvec,n,m,d,c,cgrad)
  implicit none
  integer, intent(in) :: n,m,d !training data and inner layer sizes
  real(kind=8), dimension(m*(n+2)+1), intent(in) :: fvec !fitting parameters
  real(kind=8), intent(out) :: c !cost
  real(kind=8), dimension(m*(n+2)+1), intent(out) :: cgrad !gradient of cost
  integer :: i1, j1, l1, k1
  real(kind=8), dimension(m,n) :: w_inner
  real(kind=8), dimension(m) :: b_inner,w_outer
  real(kind=8) :: b_outer

  !Declare other variables as needed
  real(kind=8), dimension(m,d) :: z, a, dadw_outer, dadb_inner
  real(kind=8), dimension(d) :: z_outer, a_outer, dadb_outer
  real(kind=8), dimension(m,d,n) :: dadw_inner
  real(kind=8) :: dcdb_outer
  real(kind=8), dimension(m) :: dcdw_outer, dcdb_inner
  real(kind=8), dimension(m,n) :: dcdw_inner

  !unpack fitting parameters (use if needed)
  do i1=1,n
    j1 = (i1-1)*m+1
    w_inner(:,i1) = fvec(j1:j1+m-1) !inner layer weight matrix
  end do
  b_inner = fvec(n*m+1:n*m+m) !inner layer bias vector
  w_outer = fvec(n*m+m+1:n*m+2*m) !output layer weight vector
  b_outer  = fvec(n*m+2*m+1) !output layer bias

  !Add code to compute c and cgrad
  call nnforward(fvec,n,m,d,a,a_outer)

  c = (1.d0/(2.d0*real(d))) * sum((a_outer - nm_y)**2)

  !!first set of equations

  !da/db_outer
  do k1 = 1, d
    dadb_outer(k1) = a_outer(k1) * (1.d0 - a_outer(k1))
  end do

  !da/dw_outer and da/db_inner
  do j1 = 1, m
    do k1 = 1, d
      dadw_outer(j1,k1) =   a(j1,k1) * a_outer(k1) * (1.d0 - a_outer(k1))

      dadb_inner(j1,k1) = a_outer(k1) * (1.d0 - a_outer(k1)) * w_outer(j1) * a(j1,k1) * (1.d0 - a(j1,k1))
    end do
  end do

  !da/dw_inner
  do j1 = 1, m
    do k1 = 1, d
      do l1 = 1, n
       dadw_inner(j1,k1,l1) = a_outer(k1) * (1.d0 - a_outer(k1)) * w_outer(j1) * a(j1, k1) * (1.d0 - a(j1,k1)) * nm_x(l1, k1)
      end do
    end do
  end do

  !!second set of equations
  dcdb_outer = 1.d0/real(d) * sum((a_outer - nm_y) * dadb_outer)

  do j1 = 1,m
    dcdw_outer(j1) = 1.d0/real(d) * sum((a_outer - nm_y)  * dadw_outer(j1,:))
    dcdb_inner(j1) = 1.d0/real(d) * sum((a_outer - nm_y) * dadb_inner(j1,:))
  end do

  do j1 = 1,m
    do l1 = 1,n
      dcdw_inner(j1,l1) = 1.d0/real(d) * sum((a_outer - nm_y) * dadw_inner(j1,:,l1))
    end do
  end do

  !!putting in gradient components
  do i1=1,n
    j1 = (i1-1)*m+1
    cgrad(j1:j1+m-1) = dcdw_inner(:,i1)!inner layer weight matrix
  end do
  cgrad(n*m+1:n*m+m) = dcdb_inner!inner layer bias vector
  cgrad(n*m+m+1:n*m+2*m) = dcdw_outer !output layer weight vector
  cgrad(n*m+2*m+1) = dcdb_outer   !output layer bias

end subroutine nnmodel


!Use crude implementation of stochastic gradient descent
!to move towards optimal fitting parameters using either
! snmodel or nnmodel. Iterates for 400 "epochs" and final fitting
!parameters are stored in fvec.
!Input:
!fvec_guess: initial vector of fitting parameters
!n: number of pixels in each image (should be 784)
!m: number of neurons in inner layer; snmodel is used if m=0
!d: number of training images to be used; only the 1st d images and labels stored
!in nm_x and nm_y are used in the optimization calculation
!alpha: learning rate, it is fine to keep this as alpha=0.1 for this assignment
!Output:
!fvec: fitting parameters, see comments above for snmodel and nnmodel to see how
!weights and biases are stored in the array.
!Note: nm_x and nm_y must be allocated and set before calling this subroutine.
subroutine sgd(fvec_guess,n,m,d,alpha,fvec)
  implicit none
  integer, intent(in) :: n,m,d
  real(kind=8), dimension(:), intent(in) :: fvec_guess
  real(kind=8), intent(in) :: alpha
  real(kind=8), dimension(size(fvec_guess)), intent(out) :: fvec
  integer :: i1, j1, i1max=400
  real(kind=8) :: c
  real(kind=8), dimension(size(fvec_guess)) :: cgrad
  real(kind=8), allocatable, dimension(:,:) :: xfull
  integer, allocatable, dimension(:) :: yfull
  real(kind=8), dimension(d) :: a
  real(kind=8), dimension(d+1) :: r
  integer, dimension(d+1) :: j1array

  !store full nm_x,nm_y
  allocate(xfull(size(nm_x,1),size(nm_x,2)),yfull(size(nm_y)))
  xfull = nm_x
  yfull = nm_y

  !will only use one image at a time, so need to reallocate nm_x,nm_y
  call data_init(n,1)

  fvec = fvec_guess
  do i1=1,i1max
    call random_number(r)
    j1array = floor(r*d+1.d0) !d random integers falling between 1 and d (inclusive); will compute c, cgrad for one image at a time cycling through these integers

    do j1 = 1,d
      nm_x(:,1) = xfull(:,j1array(j1))
      nm_y = yfull(j1array(j1))

      !compute cost and gradient with randomly selected image
      if (m==0) then
        call snmodel(fvec,n,1,c,cgrad)
      else
        call nnmodel(fvec,n,m,1,c,cgrad)
      end if
      fvec = fvec - alpha*cgrad !update fitting parameters using gradient descent step
    end do

    if (mod(i1,50)==0) print *, 'completed epoch # ', i1

  end do

 !reset nm_x,nm_y to intial state at beginning of subroutine
  call data_init(size(xfull,1),size(xfull,2))
  nm_x = xfull
  nm_y = yfull
  deallocate(xfull,yfull)

  end subroutine sgd


end module nmodel
