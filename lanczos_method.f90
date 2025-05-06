program lanczos_method
!Compile: gfortran -o lanczos_method lanczos_method.f90 -llapack -lblas
!Run: ./lanczos_method
! This program implements the Lanczos method for finding eigenvalues and eigenvectors
! of a symmetric matrix A. It reads the matrix from the file 'matrix.inp', performs 
! the Lanczos process, constructs the tridiagonal matrix T, computes its eigenvalues 
! and eigenvectors, and computes the Ritz vectors and residual norms.
! The matrix A is assumed to be symmetric and real-valued.

  implicit none
  integer, parameter :: dp = kind(1.0d0)
  integer :: n, m
  integer :: i, j, k, init_method
  real(dp), allocatable :: A(:,:), A_copy(:,:), V(:,:), ritz_vectors(:,:)
  real(dp), allocatable :: alpha(:), beta(:), eval_T(:), offdiag_T(:), Z(:,:), work_T(:)
  real(dp), allocatable :: v0(:), w(:)
  real(dp) :: norm_w, res_norm
  real(dp), allocatable :: T(:,:)  ! Declare T

  ! For eigenvalue computations
  real(dp), allocatable :: eval_A(:), work_A(:)
  integer :: lwork_A, info_A, info_T

  ! Declare variables for shift-and-invert
  integer :: use_shift_invert, info
  real(dp) :: sigma
  integer, allocatable :: ipiv(:)

  ! Prompt the user for the size of the matrix (n)
  print *, 'Enter the size of the matrix (n):'
  read(*, *) n
  if (n <= 0) then
     print *, 'Invalid matrix size. n must be greater than 0.'
     stop
  end if

  ! Prompt the user for the number of Lanczos iterations (m)
  print *, 'Enter the number of Lanczos iterations (m):'
  read(*, *) m
  if (m <= 0 .or. m > n) then
     print *, 'Invalid number of iterations. m must be between 1 and n.'
     stop
  end if

  ! Allocate arrays based on n and m
  allocate(A(n,n), A_copy(n,n), V(n,m), ritz_vectors(n,m))
  allocate(alpha(m), beta(m), eval_T(m), offdiag_T(m), Z(m,m), work_T(2*m))
  allocate(v0(n), w(n), eval_A(n))
  allocate(T(m,m))  ! Allocate memory for T

  ! Read matrix A from file
  open(unit=10, file='matrix.inp', status='old', action='read')
  do i = 1, n
     read(10, *) (A(i,j), j=1,n)
  end do
  close(10)

  ! Print original matrix A
  print *, 'Original symmetric matrix A:'
  do i = 1, n
     write(*,'(5F10.4)') (A(i,j), j=1,n)
  end do
  print *

  ! Initialize normalized vector v0
  print *, 'Choose initialization method for v0:'
  print *, '1 - Random initialization'
  print *, '2 - Read from file (v0.dat)'
  read(*, *) init_method

  if (init_method == 1) then
      ! Random initialization
      call random_seed()
      call random_number(v0)
      print *, 'Randomly initialized v0:'
  else if (init_method == 2) then
      ! Read from file
      open(unit=11, file='v0.dat', status='old', action='read')
      read(11, *) v0
      close(11)
      print *, 'v0 read from file:'
  else
      print *, 'Invalid choice. Exiting.'
      stop
  end if

  ! Normalize the vector v0
  v0 = v0 / sqrt(dot_product(v0, v0))
  V(:,1) = v0
  beta(1) = 0.0d0

  ! Ask the user whether to use the shift-and-invert technique
  print *, 'Do you want to use the shift-and-invert technique? (1 = Yes, 0 = No):'
  read(*, *) use_shift_invert

  if (use_shift_invert == 1) then
      print *, 'Enter the shift value (sigma):'
      read(*, *) sigma

      ! Transform A to (A - sigma*I)
      do i = 1, n
          A(i,i) = A(i,i) - sigma
      end do

      ! Allocate pivot array for solving linear systems
      allocate(ipiv(n))
  end if

  ! Lanczos process
  do k = 1, m
      if (use_shift_invert == 1) then
          ! Solve (A - sigma*I) * w = V(:,k)
          w = V(:,k)
          call dgesv(n, 1, A, n, ipiv, w, n, info)
          if (info /= 0) then
              print *, 'Error solving shifted system'
              stop
          end if
      else
          ! Standard matrix-vector multiplication
          w = matmul(A, V(:,k))
      end if

      if (k > 1) w = w - beta(k-1) * V(:,k-1)
      alpha(k) = dot_product(V(:,k), w)
      w = w - alpha(k) * V(:,k)
      if (k > 1) w = w - dot_product(V(:,k-1), w) * V(:,k-1)

      norm_w = sqrt(dot_product(w, w))
      if (k < m) then
          beta(k) = norm_w
          if (norm_w == 0.0d0) exit
          V(:,k+1) = w / norm_w
      end if
  end do

  ! Print Lanczos vectors
  print *, 'Lanczos vectors V (each column is a vector v_k):'
  do i = 1, n
     write(*,'(5F10.6)') (V(i,j), j=1,m)
  end do
  print *

  ! Print alpha and beta
  print *, 'Diagonal elements (alpha):'
  print '(F10.6)', alpha
  print *
  print *, 'Off-diagonal elements (beta):'
  print '(F10.6)', beta(1:m-1)
  print *

  ! Initialize T to zero
  T = 0.0d0

  ! Construct tridiagonal matrix T
  do i = 1, m
     T(i,i) = alpha(i)
     if (i < m) then
        T(i,i+1) = beta(i)
        T(i+1,i) = beta(i)
     end if
  end do

  ! Print tridiagonal matrix T
  print *, 'Tridiagonal matrix T:'
  do i = 1, m
     write(*,'(5F10.4)') (T(i,j), j=1,m)
  end do
  print *

  ! Compute eigenvalues of A
  A_copy = A
  lwork_A = 3*n
  allocate(work_A(lwork_A))
  call dsyev('N', 'U', n, A_copy, n, eval_A, work_A, lwork_A, info_A)
  if (info_A /= 0) then
     print *, 'Error computing eigenvalues of A'
  else
     print *, 'Eigenvalues of matrix A:'
     print '(F10.6)', eval_A
  end if
  deallocate(work_A)
  print *

  ! Compute eigenvalues and eigenvectors of T
  eval_T = alpha
  offdiag_T(1:m-1) = beta(1:m-1)
  Z = 0.0d0
  do i = 1, m
     Z(i,i) = 1.0d0
  end do
  call dsteqr('I', m, eval_T, offdiag_T, Z, m, work_T, info_T)
  if (info_T /= 0) then
     print *, 'Error computing eigenvalues of T'
  else
     print *, 'Eigenvalues of tridiagonal matrix T:'
     print '(F10.6)', eval_T
  end if
  print *

  ! Recover original eigenvalues if shift-and-invert was used
  if (use_shift_invert == 1) then
      print *, 'Eigenvalues of T before recovery:'
      print '(F10.6)', eval_T
      do i = 1, m
          if (abs(eval_T(i)) > 1.0d-12) then
              eval_T(i) = sigma + 1.0d0 / eval_T(i)
          else
              print *, 'Warning: Eigenvalue too close to zero during recovery.'
              eval_T(i) = sigma
          end if
      end do
      print *, 'Recovered eigenvalues of A:'
      print '(F10.6)', eval_T
  end if

  ! Compute Ritz vectors: V * Z
  ritz_vectors = matmul(V, Z)

  ! Recover original eigenvectors if shift-and-invert was used
  if (use_shift_invert == 1) then
      do j = 1, m
          ! Normalize the Ritz vector
          ritz_vectors(:,j) = ritz_vectors(:,j) / sqrt(dot_product(ritz_vectors(:,j), ritz_vectors(:,j)))
      end do
  end if

  ! Print Ritz vectors
  print *, 'Ritz vectors (approximate eigenvectors of A):'
  do j = 1, m
      print *, 'Ritz vector ', j
      write(*,'(5F10.6)') ritz_vectors(:,j)
  end do
  print *

  ! Compute residual norms ||A*x - lambda*x||
  print *, 'Residual norms ||A*x - lambda*x|| for Ritz pairs:'
  do j = 1, m
     w = matmul(A, ritz_vectors(:,j)) - eval_T(j) * ritz_vectors(:,j)
     res_norm = sqrt(dot_product(w, w))
     print '(I3, F10.6)', j, res_norm
  end do

  ! Open a file for output
  open(unit=20, file='output.txt', status='replace', action='write')

  ! Print original matrix A
  write(20, *) 'Original symmetric matrix A:'
  do i = 1, n
     write(20, '(5F10.4)') (A(i,j), j=1,n)
  end do
  write(20, *)

  ! Print Lanczos vectors
  write(20, *) 'Lanczos vectors V (each column is a vector v_k):'
  do i = 1, n
     write(20, '(5F10.6)') (V(i,j), j=1,m)
  end do
  write(20, *)

  ! Print alpha and beta
  write(20, *) 'Diagonal elements (alpha):'
  write(20, '(F10.6)', advance='no') alpha
  write(20, *)
  write(20, *) 'Off-diagonal elements (beta):'
  write(20, '(F10.6)', advance='no') beta(1:m-1)
  write(20, *)

  ! Print tridiagonal matrix T
  write(20, *) 'Tridiagonal matrix T:'
  do i = 1, m
     write(20, '(5F10.4)') (T(i,j), j=1,m)
  end do
  write(20, *)

  ! Print eigenvalues of A
  if (info_A /= 0) then
     write(20, *) 'Error computing eigenvalues of A'
  else
     write(20, *) 'Eigenvalues of matrix A:'
     write(20, '(F10.6)', advance='no') eval_A
  end if
  write(20, *)

  ! Print eigenvalues of T
  if (info_T /= 0) then
     write(20, *) 'Error computing eigenvalues of T'
  else
     write(20, *) 'Eigenvalues of tridiagonal matrix T:'
     write(20, '(F10.6)', advance='no') eval_T
  end if
  write(20, *)

  ! Print Ritz vectors
  write(20, *) 'Ritz vectors (approximate eigenvectors of A):'
  do j = 1, m
     write(20, *) 'Ritz vector ', j
     write(20, '(5F10.6)') ritz_vectors(:,j)
  end do
  write(20, *)

  ! Print residual norms
  write(20, *) 'Residual norms ||A*x - lambda*x|| for Ritz pairs:'
  do j = 1, m
     w = matmul(A, ritz_vectors(:,j)) - eval_T(j) * ritz_vectors(:,j)
     res_norm = sqrt(dot_product(w, w))
     write(20, '(I3, F10.6)') j, res_norm
  end do

  ! Close the output file
  close(20)

end program lanczos_method
