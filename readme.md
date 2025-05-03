Lanczos Method for Symmetric Matrix Eigenvalue Approximation
=============================================================

Overview
--------
This Fortran program performs eigenvalue approximation using the Lanczos algorithm 
for a real symmetric matrix. The algorithm projects the input matrix onto a Krylov 
subspace and computes the eigenvalues and eigenvectors of the resulting tridiagonal 
matrix.

It also compares the results with exact eigenvalues from LAPACK routines and computes 
the residual norms for each approximate eigenpair (Ritz pair).

Input
-----
- `matrix.dat`: A file containing the symmetric matrix A of size nxn.
  Format: each row of the matrix on a new line, with n space-separated values.

  > **Note**: In this code, both `n` (matrix size) and `m` (number of Lanczos steps) 
  > are declared as parameters:
  >
  > ```fortran
  > integer, parameter :: n = 5, m = 5
  > ```

Compilation
-----------
Make sure you have the Intel Fortran compiler (`ifx`) and Intel MKL installed.

To compile the code, use:

    ifx -qmkl -o lanczos_method lanczos_method.f90
    gfortran -o lanczos_method lanczos_method.f90 -llapack

If required, source the Intel environment setup:

    source /opt/intel/oneapi/setvars.sh intel64

Execution
---------
After compilation, run the executable:

    ./lanczos_method

The program performs the following steps:
-----------------------------------------
1. Reads a symmetric nxn matrix A from `matrix.dat`.
2. Initializes a random normalized starting vector `v0`.
3. Runs the Lanczos iteration for `m=5` steps, producing vectors `V` and scalars `alpha`, `beta`.
4. Constructs the tridiagonal matrix `T` from `alpha` and `beta`.
5. Computes:
    - The exact eigenvalues of A using LAPACK's `dsyev`.
    - The eigenvalues and eigenvectors of T using LAPACK's `dsteqr`.
6. Forms approximate eigenvectors of A (Ritz vectors) by multiplying `V * Z`.
7. Computes and displays the residual norm for each Ritz pair.

Dependencies
------------
- Intel Fortran compiler (`ifx`)
- Intel Math Kernel Library (MKL), for LAPACK routines `dsyev` and `dsteqr`

Output
------
The program prints:
- The original matrix A
- Lanczos vectors
- Tridiagonal matrix T
- Exact eigenvalues of A
- Approximate eigenvalues (from T)
- Ritz vectors
- Residual norms for each Ritz pair

Author
------
Developed by Dr. Nguyen Le Anh  
Date: May 3, 2025
