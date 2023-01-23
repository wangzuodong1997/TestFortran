!===README
! We assemble two matrices using (using mtype=2 or mtype=11; does not matter).
!
! The code crashes if one uses the multirecursive option for BOTH matrices (iparm(32)=1 and iparm2(32)=1)
! and if the size of the matrices are DIFFERENT, (n/=n2).
! Here is the message for mtype =2
!/home/oschenk/Work/pardiso/pardiso-git/pardiso/src/PILS/corelib/smat.c:1652: smat_scale: Assertion `mat->m == u->n' failed
! and for mtype = 11
!/home/oschenk/Work/pardiso/pardiso-git/pardiso/src/PILS/corelib/precon_iluc.c:135: precon_iluc: Assertion `ilu->L->m == x->n' failed.
!
! The code runs perfectly well with the multirecursive option if the two matrices have the SAME size (n=n2).
! 
! The code also runs well with two diferent size if one use iparm(32)=0 and iparm2(32)=1,
! or iparm(32)=1 and iparm2(32)=0, or iparm(32)=0 and iparm2(32)=0.
PROGRAM pardiso_sym
  IMPLICIT NONE
  INTEGER, PARAMETER :: dp = KIND(1.0d0)
  TYPE PARDISO_HANDLE
     INTEGER :: DUMMY
  END TYPE PARDISO_HANDLE
  TYPE(PARDISO_HANDLE), ALLOCATABLE :: pt(:), pt2(:)
  INTEGER maxfct, mnum, mnum2, mtype, phase, n, n2, nrhs, error, msglvl, nnz, nnz2
  INTEGER, ALLOCATABLE :: iparm(:), iparm2(:)
  REAL(KIND = DP), ALLOCATABLE :: dparm(:), dparm2(:)
  INTEGER, ALLOCATABLE :: ia(:), ia2(:)
  INTEGER, ALLOCATABLE :: ja(:), ja2(:)
  REAL(KIND = DP), ALLOCATABLE :: a(:), a2(:)
  REAL(KIND = DP), ALLOCATABLE :: b(:), b2(:)
  REAL(KIND = DP), ALLOCATABLE :: x(:), x2(:)
  REAL(KIND = DP), ALLOCATABLE :: rhs(:), rhs2(:)
  INTEGER i, k, p, idum(1)
  REAL(KIND = DP) ddum(1)
  REAL(KIND = DP) :: aa(3), aa2(3)
  aa  = (/-1.d0,4.d0,-1.d0/)
  aa2 = (/-1.d0,3.d0,-1.d0/)
 
!!$  !===Non symmetric matrices
!!$  !===Matrix 1
!!$  mnum = 1
!!$  n = 100000
!!$  ALLOCATE(rhs(n))
!!$  nnz = 3*(n-2) + 2
!!$  
!!$  ALLOCATE(ia (n + 1))
!!$  ALLOCATE(ja (nnz))
!!$  ALLOCATE(a (nnz))
!!$  ia(1) = 1
!!$  ia(2) = ia(1) + 1
!!$  ja(1) = 1
!!$  a(ja(1)) = 1
!!$  DO i = 2, n-1
!!$     ia (i+1) = ia(i) + 3
!!$     DO k = 1, 3
!!$        p = ia(i) + k - 1
!!$        ja(p) = i + k - 2
!!$        a(p) = aa(k)
!!$     END DO
!!$  END DO
!!$  ia (n+1) = ia(n) + 1
!!$  ja(ia(n)) = n
!!$  a(ia(n)) = 1
!!$
!!$   !===Matrix 2
!!$  mnum2 = 2
!!$  n2 = 20000
!!$  ALLOCATE(rhs2(n2))
!!$  nnz2 = 3*(n2-2) + 2      +1
!!$  
!!$  ALLOCATE(ia2 (n2 + 1))
!!$  ALLOCATE(ja2 (nnz2))
!!$  ALLOCATE(a2 (nnz2))
!!$  ia2(1) = 1
!!$  ia2(2) = ia2(1) + 2
!!$  ja2(1) = 1
!!$  ja2(2) = 3
!!$  
!!$  a2(ja2(1)) = 1
!!$  DO i = 2, n2-1
!!$     ia2(i+1) = ia2(i) + 3
!!$     DO k = 1, 3
!!$        p = ia2(i) + k - 1
!!$        ja2(p) = i + k - 2
!!$        a2(p) = aa2(k)
!!$     END DO
!!$  END DO
!!$  ia2(n2+1) = ia2(n2) + 1
!!$  ja2(ia2(n2)) = n2
!!$  a2(ia2(n2)) = 1
!!$

  !===Symmetric matrices
  !===Matrix 1
  n = 10000
  mnum = 1
  ALLOCATE(rhs(n))
  nnz = 2*(n-2) + 3
  ALLOCATE(ia (n + 1))
  ALLOCATE(ja (nnz))
  ALLOCATE(a (nnz))
  ia(1) = 1
  ia(2) = ia(1) + 2
  ja(1) = 1
  ja(2) = 2
  a(ja(1)) = 1
  a(ja(2)) = -1
  DO i = 2, n-1
     ia (i+1) = ia(i) + 2
     DO k = 1, 2
        p = ia(i) + k - 1
        ja(p) = i + k - 1
        a(p) = aa(k+1)
     END DO
  END DO
  ia (n+1) = ia(n) + 1
  ja(ia(n)) = n
  a(ia(n)) = 1

  !===Matrix 2
  n2 = n   !===BUG for mutirecursive option if n2/= n
  mnum2 = 2
  ALLOCATE(rhs2(n2))
  nnz2 = 2*(n2-2) + 3
  ALLOCATE(ia2 (n2 + 1))
  ALLOCATE(ja2 (nnz2))
  ALLOCATE(a2 (nnz2))
  ia2(1) = 1
  ia2(2) = ia2(1) + 2
  ja2(1) = 1
  ja2(2) = 2
  a2(ja2(1)) = 1
  a2(ja2(2)) = -1
  DO i = 2, n2-1
     ia2 (i+1) = ia2(i) + 2
     DO k = 1, 2
        p = ia2(i) + k - 1
        ja2(p) = i + k - 1
        a2(p) = aa2(k+1)
     END DO
  END DO
  ia2 (n2+1) = ia2(n2) + 1
  ja2(ia2(n2)) = n2 
  a2(ia2(n2)) =  1

  !===Right hand side solution
  ALLOCATE(b (n))
  ALLOCATE(x (n))
  ALLOCATE(b2 (n2))
  ALLOCATE(x2 (n2))
 
  !===Set up PARDISO control parameter
  nrhs = 1
  maxfct = 2
  error = 0  ! initialize error flag
  msglvl = 0 ! print statistical information
  mtype = 2  !===Real symmetric !=== 11: real non symmetric
  !===iparm parameters
  ALLOCATE(iparm (64))
  ALLOCATE(iparm2 (64))
  DO i = 1, 64
     iparm(i) = 0
     iparm2(i) = 0
  END DO
  
  iparm(1) = 1 ! no solver default
  iparm(2) = 2 ! fill-in reordering from METIS
  iparm(4) = 0 ! no iterative-direct algorithm
  iparm(5) = 0 ! no user fill-in reducing permutation
  iparm(6) = 0 ! =0 solution on the first n compoments of x
  iparm(8) = 9 ! numbers of iterative refinement steps
  iparm(10) = 13 ! perturbe the pivot elements with 1E-13
  iparm(11) = 1 ! use nonsymmetric permutation and scaling MPS
  iparm(13) = 0 ! maximum weighted matching algorithm is switched-off (default for symmetric). 
  iparm(14) = 0 ! Output: number of perturbed pivots
  iparm(18) = -1 ! Output: number of nonzeros in the factor LU
  iparm(19) = -1 ! Output: Mflops for LU factorization
  iparm(20) = 0 ! Output: Numbers of CG Iterations

  iparm2(1) = 1 ! no solver default
  iparm2(2) = 2 ! fill-in reordering from METIS
  iparm2(4) = 0 ! no iterative-direct algorithm
  iparm2(5) = 0 ! no user fill-in reducing permutation
  iparm2(6) = 0 ! =0 solution on the first n compoments of x
  iparm2(8) = 9 ! numbers of iterative refinement steps
  iparm2(10) = 13 ! perturbe the pivot elements with 1E-13
  iparm2(11) = 1 ! use nonsymmetric permutation and scaling MPS
  iparm2(13) = 0 ! maximum weighted matching algorithm is switched-off (default for symmetric).
  iparm2(14) = 0 ! Output: number of perturbed pivots
  iparm2(18) = -1 ! Output: number of nonzeros in the factor LU
  iparm2(19) = -1 ! Output: Mflops for LU factorization
  iparm2(20) = 0 ! Output: Numbers of CG Iterations

  !===Solver type. 0: direct solve; 1: multirecursive
  iparm(32)  = 1 !===SOLVER TYPE HERE
  iparm2(32) = 1 !===SOLVER TYPE HERE

  !===dparm parameters
  ALLOCATE(dparm(64))
  ALLOCATE(dparm2(64))
  DO i = 1, 64
     dparm(i) = 0.d0
     dparm2(i) = 0.d0
  END DO

  dparm(1) = 300
  dparm(2) = 1.d-15
  dparm(3) = 5000
  dparm(4) = 10
  dparm(5) = 1.d-2
  dparm(6) = 5.d-3
  dparm(7) = 10
  dparm(8) = 500
  dparm(9) = 25
  
  dparm2(1) = 300
  dparm2(2) = 1.d-15
  dparm2(3) = 5000
  dparm2(4) = 10
  dparm2(5) = 1.d-2
  dparm2(6) = 5.d-3
  dparm2(7) = 10
  dparm2(8) = 500
  dparm2(9) = 25

  !===Pointers
  ALLOCATE (pt (64))
  ALLOCATE (pt2 (64))
  DO i = 1, 64
     pt(i)%DUMMY = 0
     pt2(i)%DUMMY = 0
  END DO
  
  !===Reordering and Symbolic Factorization
  phase = 11 !===only reordering and symbolic factorization
  CALL pardiso (pt,  maxfct, mnum,  mtype, phase, n,  a,  ia,  ja,  idum, nrhs, iparm,  msglvl, ddum, ddum, error, dparm)
  CALL pardiso (pt2, maxfct, mnum2, mtype, phase, n2, a2, ia2, ja2, idum, nrhs, iparm2, msglvl, ddum, ddum, error, dparm2)
  WRITE(*, *) 'Reordering completed ... '
  IF (error /= 0) THEN
     WRITE(*, *) 'The following ERROR was detected: ', error
     STOP
  END IF
  WRITE(*, *) 'Number of nonzeros in factors = ', iparm(18)
  WRITE(*, *) 'Number of factorization MFLOPS = ', iparm(19)

  !===Factorization.
  phase = 22 !===only factorization
  CALL pardiso (pt,  maxfct, mnum,  mtype, phase, n,  a,  ia,  ja, &
       idum, nrhs, iparm,  msglvl, ddum, ddum, error, dparm)
  CALL pardiso (pt2, maxfct, mnum2, mtype, phase, n2, a2, ia2, ja2, &
       idum, nrhs, iparm2, msglvl, ddum, ddum, error, dparm2)
  WRITE(*, *) 'Factorization completed ... '
  IF (error /= 0) THEN
     WRITE(*, *) 'The following ERROR was detected: ', error
     STOP
  ENDIF

  !===Back substitution and iterative refinement
  iparm(8) = 2 !===max numbers of iterative refinement steps
  phase = 33 !==only solution

  !===Set up rhight-hand sides
  b = 1.d0
  b2 = 2.d0

  !===First series
  !===First linear system
  CALL pardiso (pt, maxfct, mnum, mtype, phase, n, a, ia, ja, &
       idum, nrhs, iparm, msglvl, b, x, error, dparm)
  WRITE(*, *) 'Solve 1 completed ... '
  IF (error /= 0) THEN
     WRITE(*, *) 'The following ERROR was detected: ', error
     STOP
  ENDIF
  CALL Verify(a,ia,ja,x,rhs,mtype)
  WRITE(*, *) "Error 1", MAXVAL(ABS(rhs-b))/MAXVAL(ABS(b))
  
  !===Second linear system
  CALL pardiso (pt2, maxfct, mnum2, mtype, phase, n2, a2, ia2, ja2, &
       idum, nrhs, iparm2, msglvl, b2, x2, error, dparm2)
  WRITE(*, *) 'Solve 2 completed ... '
  IF (error /= 0) THEN
     WRITE(*, *) 'The following ERROR was detected: ', error
     STOP
  ENDIF
  CALL Verify(a2,ia2,ja2,x2,rhs2,mtype)
  WRITE(*, *) "Error 2", MAXVAL(ABS(rhs2-b2))/MAXVAL(ABS(b2))

  !===Second series
  !===First linear system
  b = x
  b2 = x2
  CALL pardiso (pt, maxfct, mnum, mtype, phase, n, a, ia, ja, &
       idum, nrhs, iparm, msglvl, b, x, error, dparm)
  WRITE(*, *) 'Second solve 1 completed ... '
  IF (error /= 0) THEN
     WRITE(*, *) 'The following ERROR was detected: ', error
     STOP
  ENDIF
  CALL Verify(a,ia,ja,x,rhs,mtype)
  WRITE(*, *) "Error 1", MAXVAL(ABS(rhs-b))/MAXVAL(ABS(b))
  
  !===Second linear system
  CALL pardiso (pt2, maxfct, mnum2, mtype, phase, n2, a2, ia2, ja2, &
       idum, nrhs, iparm2, msglvl, b2, x2, error, dparm2)
  WRITE(*, *) 'Second solve 2 completed ... '
  IF (error /= 0) THEN
     WRITE(*, *) 'The following ERROR was detected: ', error
     STOP
  ENDIF
  CALL Verify(a2,ia2,ja2,x2,rhs2,mtype)
  WRITE(*, *) "Error 2", MAXVAL(ABS(rhs2-b2))/MAXVAL(ABS(b2))
CONTAINS
  SUBROUTINE verify(a,ia,ja,x,rhs,mtype)
    IMPLICIT NONE
    REAL(KIND = DP), DIMENSION(:) :: a, x, rhs
    INTEGER, DIMENSION(:) :: ia, ja
    INTEGER :: mtype, i, j, n
    n = SIZE(ia)-1
    DO i = 1, n
       rhs(i) = SUM(a(ia(i):ia(i+1)-1)*x(ja(ia(i):ia(i+1)-1)))
    END DO
    SELECT CASE(mtype)
    CASE(1,2,-2)
       DO j = 1, n
          DO p = ia(j)+1, ia(j+1) -1
             i = ja(p)
             rhs(i) = rhs(i) + a(p)*x(j)
          END DO
       END DO
    END SELECT
  END SUBROUTINE verify
END PROGRAM pardiso_sym
