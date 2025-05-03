program lanczos_method
!ifx: source /opt/intel/oneapi/setvars.sh intel64
!Compile: ifx -qmkl -o lanczos_method lanczos_method.f90
  implicit none
  integer, parameter :: dp = kind(1.0d0)
  integer, parameter :: n = 5, m = 5
  integer :: i, j, k
  real(dp), dimension(n,n) :: A, A_copy
  real(dp), dimension(n,m) :: V
  real(dp), dimension(m) :: alpha, beta
  real(dp), dimension(m,m) :: T
  real(dp), dimension(n) :: v0, w
  real(dp) :: norm_w

  ! For eigenvalue computations
  real(dp), dimension(n) :: eval_A
  real(dp), allocatable :: work_A(:)
  integer :: lwork_A, info_A

  real(dp), dimension(m) :: eval_T, offdiag_T
  real(dp), dimension(m,m) :: Z
  real(dp), dimension(2*m) :: work_T
  integer :: info_T

  real(dp), dimension(n,m) :: ritz_vectors  ! Declare Ritz vectors array
  real(dp) :: res_norm  ! Declare residual norm variable

  interface
    subroutine dsyev(jobz, uplo, n, a, lda, w, work, lwork, info)
      character(len=1), intent(in) :: jobz, uplo
      integer, intent(in) :: n, lda, lwork
      double precision, intent(inout) :: a(lda,*)
      double precision, intent(out) :: w(n), work(*)
      integer, intent(out) :: info
    end subroutine dsyev

    subroutine dsteqr(compz, n, d, e, z, ldz, work, info)
      character(len=1), intent(in) :: compz
      integer, intent(in) :: n, ldz
      double precision, intent(inout) :: d(n), e(n-1), z(ldz,*), work(*)
      integer, intent(out) :: info
    end subroutine dsteqr
  end interface

  ! Read matrix A from file
  open(unit=10, file='matrix.dat', status='old', action='read')
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

  ! Initialize normalized random vector v0
  call random_seed()
  call random_number(v0)
  v0 = v0 / sqrt(dot_product(v0, v0))
  V(:,1) = v0
  beta(1) = 0.0d0

  ! Lanczos process
  do k = 1, m
     w = matmul(A, V(:,k))
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

  ! Construct tridiagonal matrix T
  T = 0.0d0
  do i = 1, m
     T(i,i) = alpha(i)
     if (i < m) then
        T(i,i+1) = beta(i)
        T(i+1,i) = beta(i)
     end if
  end do

  print *, 'Tridiagonal matrix T:'
  do i = 1, m
     write(*,'(5F10.4)') (T(i,j), j=1,m)
  end do
  print *

  ! --- Eigenvalues of matrix A ---
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

  ! Eigenvalues and eigenvectors of matrix T
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

  ! Compute Ritz vectors: V * Z
  ritz_vectors = matmul(V, Z)

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

end program lanczos_method
