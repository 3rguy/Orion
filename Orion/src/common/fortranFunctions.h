/* Some fortran function declarations.*/

#ifndef fortranFunctions_h_
#define fortranFunctions_h_

typedef int idxtype;

/************************************************************************/
// Carlo's viscoplastic constitutive law

extern "C" void tangent_(double* strain,double* history,double* stress,
			 double* constitutiveTens,double& kappa,
			 double& G,double& D0,double& Z0,double& Z1,
			 double& AM,double& AN,double& Dt,
			 int& iterationStep,int& rank,int& ierr);


extern "C" void invda_(double* A,double* ALN);


extern "C" void dexpo_(double* A,double* exp,double* dexp,double& tol,
		       int& rank,int& ierr);


/************************************************************************/
// convert a integer number to a string

extern "C" void convertIntString_(int& number,char* string,int& size);

// open log files for fortran code
extern "C" void openfortranlogs_(int& rank);

/***********************************************************************/
/* Computation of a suitable integration point distribution over the used
   processors in parallel.*/
extern "C" void ParMETIS_PartKway(idxtype*,idxtype*,idxtype*,idxtype*,
				  idxtype*,int*,int*,int*,int*,int*,
				  idxtype*,MPI_Comm*);

/***********************************************************************
*     sgefa factors a real matrix by gaussian elimination.
*
*     sgefa is usually called by sgeco, but it can be called
*     directly with a saving in time if  rcond  is not needed.
*     (time for sgeco) = (1 + 9/n)*(time for sgefa) .
*
*     on entry
*
*        a       real(lda, n)
*                the matrix to be factored.
*
*        lda     integer
*                the leading dimension of the array  a .
*
*        n       integer
*                the order of the matrix  a .
*
*     on return
*
*        a       an upper triangular matrix and the multipliers
*                which were used to obtain it.
*                the factorization can be written  a = l*u  where
*                l  is a product of permutation and unit lower
*                triangular matrices and  u  is upper triangular.
*
*        ipvt    integer(n)
*                an integer vector of pivot indices.
*
*        info    integer
*                = 0  normal value.
*                = k  if  u(k,k) .eq. 0.0 .  this is not an error
*                     condition for this subroutine, but it does
*                     indicate that sgesl or sgedi will divide by zero
*                     if called.  use  rcond  in sgeco for a reliable
*                     indication of singularity.*/

extern "C" void sgefa_(float* A,int& lda,int& n,int* ipvt,int& info);


/***********************************************************************
*     sgedi computes the determinant and inverse of a matrix
*     using the factors computed by sgeco or sgefa.
*
*     on entry
*
*        a       real(lda, n)
*                the output from sgeco or sgefa.
*
*        lda     integer
*                the leading dimension of the array  a .
*
*        n       integer
*                the order of the matrix  a .
*
*        ipvt    integer(n)
*                the pivot vector from sgeco or sgefa.
*
*        work    real(n)
*                work vector.  contents destroyed.
*
*        job     integer
*                = 11   both determinant and inverse.
*                = 01   inverse only.
*                = 10   determinant only.
*
*     on return
*
*        a       inverse of original matrix if requested.
*                otherwise unchanged.
*
*        det     real(2)
*                determinant of original matrix if requested.
*                otherwise not referenced.
*                determinant = det(1) * 10.0**det(2)
*                with  1.0 .le. abs(det(1)) .lt. 10.0
*                or  det(1) .eq. 0.0 .
*
*     error condition
*
*        a division by zero will occur if the input factor contains
*        a zero on the diagonal and the inverse is requested.
*        it will not occur if the subroutines are called correctly
*        and if sgeco has set rcond .gt. 0.0 or sgefa has set
*        info .eq. 0 .  */

extern "C" void sgedi_(float* A,int& lda,int& n,int* ipvt,float* det,
		       float* work,int& job);


/***********************************************************************/
/***********************************************************************/
/*    SUBROUTINE DGECO(a,lda,n,ipvt,rcond,z)
 *    integer lda,n,ipvt(1)
 *    double precision a(lda,1),z(1)
 *    double precision rcond
 *
 *    dgeco factors a double precision matrix by gaussian elimination
 *    and estimates the condition of the matrix.
 *
 *    if  rcond  is not needed, dgefa is slightly faster.
 *    to solve  a*x = b , follow dgeco by dgesl.
 *    to compute  inverse(a)*c , follow dgeco by dgesl.
 *    to compute  determinant(a) , follow dgeco by dgedi.
 *    to compute  inverse(a) , follow dgeco by dgedi.
 *
 *    on entry
 *
 *       a       double precision(lda, n)
 *               the matrix to be factored.
 *
 *       lda     integer
 *               the leading dimension of the array  a .
 *
 *       n       integer
 *               the order of the matrix  a .
 *
 *    on return
 *
 *       a       an upper triangular matrix and the multipliers
 *               which were used to obtain it.
 *               the factorization can be written  a = l*u  where
 *               l  is a product of permutation and unit lower
 *               triangular matrices and  u  is upper triangular.
 *
 *       ipvt    integer(n)
 *               an integer vector of pivot indices.
 *
 *       rcond   double precision
 *               an estimate of the reciprocal condition of  a .
 *               for the system  a*x = b , relative perturbations
 *               in  a  and  b  of size  epsilon  may cause
 *               relative perturbations in  x  of size  epsilon/rcond .
 *               if  rcond  is so small that the logical expression
 *                          1.0 + rcond .eq. 1.0
 *               is true, then  a  may be singular to working
 *               precision.  in particular,  rcond  is zero  if
 *               exact singularity is detected or the estimate
 *               underflows.
 *
 *       z       double precision(n)
 *               a work vector whose contents are usually unimportant.
 *               if  a  is close to a singular matrix, then  z  is
 *               an approximate null vector in the sense that
 *               norm(a*z) = rcond*norm(a)*norm(z)
 */

extern "C" void dgeco_(double* a,int& lda,int& n,int* ipvt,double& rcond,
		       double* z);



/***********************************************************************/
/***********************************************************************/
/*     SUBROUTINE DGEDI (a,lda,n,ipvt,det,work,job)
 *      integer lda,n,ipvt(1),job
 *      double precision a(lda,1),det(2),work(1)
 *
 *    dgedi computes the determinant and inverse of a matrix
 *    using the factors computed by dgeco or dgefa.
 *
 *    on entry
 *
 *       a       double precision(lda, n)
 *               the output from dgeco or dgefa.
 *
 *       lda     integer
 *               the leading dimension of the array  a .
 *
 *       n       integer
 *               the order of the matrix  a .
 *
 *       ipvt    integer(n)
 *               the pivot vector from dgeco or dgefa.
 *
 *       work    double precision(n)
 *               work vector.  contents destroyed.
 *
 *       job     integer
 *               = 11   both determinant and inverse.
 *               = 01   inverse only.
 *               = 10   determinant only.
 *
 *    on return
 *
 *       a       inverse of original matrix if requested.
 *               otherwise unchanged.
 *
 *       det     double precision(2)
 *               determinant of original matrix if requested.
 *               otherwise not referenced.
 *               determinant = det(1) * 10.0**det(2)
 *               with  1.0 .le. dabs(det(1)) .lt. 10.0
 *               or  det(1) .eq. 0.0 .
 *
 *    error condition
 *
 *       a division by zero will occur if the input factor contains
 *       a zero on the diagonal and the inverse is requested.
 *       it will not occur if the subroutines are called correctly
 *       and if dgeco has set rcond .gt. 0.0 or dgefa has set
 *       info .eq. 0 .
 *
 *    linpack. this version dated 08/14/78 .
 *    cleve moler, university of new mexico, argonne national lab.
 *
 *    subroutines and functions
 *
 *    blas daxpy,dscal,dswap
 *    fortran dabs,mod
 */

extern "C" void dgedi_(double* a,int& lda,int& n,int* ipvt,double* det,
		       double* work,int& job);


/***********************************************************************/
/***********************************************************************/
/*     SUBROUTINE DGETRF( M, N, A, LDA, IPIV, INFO )
 *
 *
 *     .. Scalar Arguments ..
 *     INTEGER            INFO, LDA, M, N
 *     ..
 *     .. Array Arguments ..
 *     INTEGER            IPIV( * )
 *     DOUBLE PRECISION   A( LDA, * )
 *     ..
 *
 *  Purpose
 *  =======
 *
 *  DGETRF computes an LU factorization of a general m-by-n matrix A
 *  using partial pivoting with row interchanges.
 *
 *  The factorization has the form
 *     A = P * L * U
 *  where P is a permutation matrix, L is lower triangular with unit
 *  diagonal elements (lower trapezoidal if m > n), and U is upper
 *  triangular (upper trapezoidal if m < n).
 *
 *  This is the right-looking Level 3 BLAS version of the algorithm.
 *
 *  Arguments
 *  =========
 *
 *  M       (input) INTEGER
 *          The number of rows of the matrix A.  M >= 0.
 *
 *  N       (input) INTEGER
 *          The number of columns of the matrix A.  N >= 0.
 *
 *  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
 *          On entry, the m by n matrix to be factored.
 *          On exit, the factors L and U from the factorization
 *          A = P*L*U; the unit diagonal elements of L are not stored.
 *
 *  LDA     (input) INTEGER
 *          The leading dimension of the array A.  LDA >= max(1,M).
 *
 *  IPIV    (output) INTEGER array, dimension (min(M,N))
 *          The pivot indices; for 1 <= i <= min(M,N), row i of the
 *          matrix was interchanged with row IPIV(i).
 *
 *  INFO    (output) INTEGER
 *          = 0: successful exit
 *          < 0: if INFO = -k, the k-th argument had an illegal value
 *          > 0: if INFO = k, U(k,k) is exactly zero. The factorization
 *               has been completed, but the factor U is exactly
 *               singular, and division by zero will occur if it is used
 *               to solve a system of equations.
 */

extern "C" void dgetrf_(int& M,int& N,double* A,int& LDA,int* IPIV,
			int& INFO);

/***********************************************************************/
/*      SUBROUTINE DGETRI( N, A, LDA, IPIV, WORK, LWORK, INFO )
 *
 *     .. Scalar Arguments ..
 *     INTEGER            INFO, LDA, LWORK, N
 *     ..
 *     .. Array Arguments ..
 *     INTEGER            IPIV( * )
 *     DOUBLE PRECISION   A( LDA, * ), WORK( LWORK )
 *     ..
 *
 *  Purpose
 *  =======
 *
 *  DGETRI computes the inverse of a matrix using the LU factorization
 *  computed by DGETRF.
 *
 *  This method inverts U and then computes inv(A) by solving the system
 *  inv(A)*L = inv(U) for inv(A).
 *
 *  Arguments
 *  =========
 *
 *  N       (input) INTEGER
 *          The order of the matrix A.  N >= 0.
 *
 *  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
 *          On entry, the factors L and U from the factorization
 *          A = P*L*U as computed by DGETRF.
 *          On exit, if INFO = 0, the inverse of the original matrix A.
 *
 *  LDA     (input) INTEGER
 *          The leading dimension of the array A.  LDA >= max(1,N).
 *
 *  IPIV    (input) INTEGER array, dimension (N)
 *          The pivot indices from DGETRF; for 1<=i<=N, row i of the
 *          matrix was interchanged with row IPIV(i).
 *
 *  WORK    (workspace/output) DOUBLE PRECISION array, dimension (LWORK)
 *          On exit, if INFO=0, then WORK(1) returns the optimal LWORK.
 *
 *  LWORK   (input) INTEGER
 *          The dimension of the array WORK.  LWORK >= max(1,N).
 *          For optimal performance LWORK >= N*NB, where NB is
 *          the optimal blocksize returned by ILAENV.
 *
 *  INFO    (output) INTEGER
 *          = 0:  successful exit
 *          < 0:  if INFO = -i, the i-th argument had an illegal value
 *          > 0:  if INFO = i, U(i,i) is exactly zero; the matrix is
 *                singular and its inverse could not be computed.
 *
 */
extern "C" void dgetri_(int& N,double* A,int& LDA,int* IPIV,double* WORK,
			int& LWORK,int& INFO);


/***********************************************************************/
/*      SUBROUTINE DGEEV( JOBVL, JOBVR, N, A, LDA, WR, WI, VL, LDVL, VR,
 *    $                  LDVR, WORK, LWORK, INFO )
 *
 *  -- LAPACK driver routine (version 2.0) --
 *     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
 *     Courant Institute, Argonne National Lab, and Rice University
 *     September 30, 1994
 *
 *     .. Scalar Arguments ..
 *     CHARACTER          JOBVL, JOBVR
 *     INTEGER            INFO, LDA, LDVL, LDVR, LWORK, N
 *     ..
 *     .. Array Arguments ..
 *     DOUBLE PRECISION   A( LDA, * ), VL( LDVL, * ), VR( LDVR, * ),
 *    $                   WI( * ), WORK( * ), WR( * )
 *     ..
 *
 *  Purpose
 *  =======
 *
 *  DGEEV computes for an N-by-N real nonsymmetric matrix A, the
 *  eigenvalues and, optionally, the left and/or right eigenvectors.
 *
 *  The right eigenvector v(j) of A satisfies
 *                   A * v(j) = lambda(j) * v(j)
 *  where lambda(j) is its eigenvalue.
 *  The left eigenvector u(j) of A satisfies
 *                u(j)**H * A = lambda(j) * u(j)**H
 *  where u(j)**H denotes the conjugate transpose of u(j).
 *
 *  The computed eigenvectors are normalized to have Euclidean norm
 *  equal to 1 and largest component real.
 *
 *  Arguments
 *  =========
 *
 *  JOBVL   (input) CHARACTER*1
 *          = 'N': left eigenvectors of A are not computed;
 *          = 'V': left eigenvectors of A are computed.
 *
 *  JOBVR   (input) CHARACTER*1
 *          = 'N': right eigenvectors of A are not computed;
 *          = 'V': right eigenvectors of A are computed.
 *
 *  N       (input) INTEGER
 *          The order of the matrix A. N >= 0.
 *
 *  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
 *          On entry, the N-by-N matrix A.
 *          On exit, A has been overwritten.
 *
 *  LDA     (input) INTEGER
 *          The leading dimension of the array A.  LDA >= max(1,N).
 *
 *  WR      (output) DOUBLE PRECISION array, dimension (N)
 *  WI      (output) DOUBLE PRECISION array, dimension (N)
 *          WR and WI contain the real and imaginary parts,
 *          respectively, of the computed eigenvalues.  Complex
 *          conjugate pairs of eigenvalues appear consecutively
 *          with the eigenvalue having the positive imaginary part
 *          first.
 *
 *  VL      (output) DOUBLE PRECISION array, dimension (LDVL,N)
 *          If JOBVL = 'V', the left eigenvectors u(j) are stored one
 *          after another in the columns of VL, in the same order
 *          as their eigenvalues.
 *          If JOBVL = 'N', VL is not referenced.
 *          If the j-th eigenvalue is real, then u(j) = VL(:,j),
 *          the j-th column of VL.
 *          If the j-th and (j+1)-st eigenvalues form a complex
 *          conjugate pair, then u(j) = VL(:,j) + i*VL(:,j+1) and
 *          u(j+1) = VL(:,j) - i*VL(:,j+1).
 *
 *  LDVL    (input) INTEGER
 *          The leading dimension of the array VL.  LDVL >= 1; if
 *          JOBVL = 'V', LDVL >= N.
 *
 *  VR      (output) DOUBLE PRECISION array, dimension (LDVR,N)
 *          If JOBVR = 'V', the right eigenvectors v(j) are stored one
 *          after another in the columns of VR, in the same order
 *          as their eigenvalues.
 *          If JOBVR = 'N', VR is not referenced.
 *          If the j-th eigenvalue is real, then v(j) = VR(:,j),
 *          the j-th column of VR.
 *          If the j-th and (j+1)-st eigenvalues form a complex
 *          conjugate pair, then v(j) = VR(:,j) + i*VR(:,j+1) and
 *          v(j+1) = VR(:,j) - i*VR(:,j+1).
 *
 *  LDVR    (input) INTEGER
 *          The leading dimension of the array VR.  LDVR >= 1; if
 *          JOBVR = 'V', LDVR >= N.
 *
 *  WORK    (workspace/output) DOUBLE PRECISION array, dimension (LWORK)
 *          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
 *
 *  LWORK   (input) INTEGER
 *          The dimension of the array WORK.  LWORK >= max(1,3*N), and
 *          if JOBVL = 'V' or JOBVR = 'V', LWORK >= 4*N.  For good
 *          performance, LWORK must generally be larger.
 *
 *  INFO    (output) INTEGER
 *          = 0:  successful exit
 *          < 0:  if INFO = -i, the i-th argument had an illegal value.
 *          > 0:  if INFO = i, the QR algorithm failed to compute all the
 *                eigenvalues, and no eigenvectors have been computed;
 *                elements i+1:N of WR and WI contain eigenvalues which
 *                have converged.
 */

extern "C" void dgeev_(char& JOBVL,char& JOBVR,int& N,double* A,int& LDA,
		       double* WR,double* WI,double* VL,int& LDVL,
		       double* VR,int& LDVR,double* WORK,int& LWORK,
		       int& INFO);


/***********************************************************************/
/*      SUBROUTINE DGGEV( JOBVL, JOBVR, N, A, LDA, B, LDB, ALPHAR, ALPHAI,
 *    $                  BETA, VL, LDVL, VR, LDVR, WORK, LWORK, INFO )
 *
 *  -- LAPACK driver routine (version 3.0) --
 *     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
 *     Courant Institute, Argonne National Lab, and Rice University
 *     June 30, 1999
 *
 *     .. Scalar Arguments ..
 *     CHARACTER          JOBVL, JOBVR
 *     INTEGER            INFO, LDA, LDB, LDVL, LDVR, LWORK, N
 *     ..
 *     .. Array Arguments ..
 *     DOUBLE PRECISION   A( LDA, * ), ALPHAI( * ), ALPHAR( * ),
 *    $                   B( LDB, * ), BETA( * ), VL( LDVL, * ),
 *    $                   VR( LDVR, * ), WORK( * )
 *     ..
 *
 *  Purpose
 *  =======
 *
 *  DGGEV computes for a pair of N-by-N real nonsymmetric matrices (A,B)
 *  the generalized eigenvalues, and optionally, the left and/or right
 *  generalized eigenvectors.
 *
 *  A generalized eigenvalue for a pair of matrices (A,B) is a scalar
 *  lambda or a ratio alpha/beta = lambda, such that A - lambda*B is
 *  singular. It is usually represented as the pair (alpha,beta), as
 *  there is a reasonable interpretation for beta=0, and even for both
 *  being zero.
 *
 *  The right eigenvector v(j) corresponding to the eigenvalue lambda(j)
 *  of (A,B) satisfies
 *
 *                   A * v(j) = lambda(j) * B * v(j).
 *
 *  The left eigenvector u(j) corresponding to the eigenvalue lambda(j)
 *  of (A,B) satisfies
 *
 *                   u(j)**H * A  = lambda(j) * u(j)**H * B .
 *
 *  where u(j)**H is the conjugate-transpose of u(j).
 *
 *
 *  Arguments
 *  =========
 *
 *  JOBVL   (input) CHARACTER*1
 *          = 'N':  do not compute the left generalized eigenvectors;
 *          = 'V':  compute the left generalized eigenvectors.
 *
 *  JOBVR   (input) CHARACTER*1
 *          = 'N':  do not compute the right generalized eigenvectors;
 *          = 'V':  compute the right generalized eigenvectors.
 *
 *  N       (input) INTEGER
 *          The order of the matrices A, B, VL, and VR.  N >= 0.
 *
 *  A       (input/output) DOUBLE PRECISION array, dimension (LDA, N)
 *          On entry, the matrix A in the pair (A,B).
 *          On exit, A has been overwritten.
 *
 *  LDA     (input) INTEGER
 *          The leading dimension of A.  LDA >= max(1,N).
 *
 *  B       (input/output) DOUBLE PRECISION array, dimension (LDB, N)
 *          On entry, the matrix B in the pair (A,B).
 *          On exit, B has been overwritten.
 *
 *  LDB     (input) INTEGER
 *          The leading dimension of B.  LDB >= max(1,N).
 *
 *  ALPHAR  (output) DOUBLE PRECISION array, dimension (N)
 *  ALPHAI  (output) DOUBLE PRECISION array, dimension (N)
 *  BETA    (output) DOUBLE PRECISION array, dimension (N)
 *          On exit, (ALPHAR(j) + ALPHAI(j)*i)/BETA(j), j=1,...,N, will
 *          be the generalized eigenvalues.  If ALPHAI(j) is zero, then
 *          the j-th eigenvalue is real; if positive, then the j-th and
 *          (j+1)-st eigenvalues are a complex conjugate pair, with
 *          ALPHAI(j+1) negative.
 *
 *          Note: the quotients ALPHAR(j)/BETA(j) and ALPHAI(j)/BETA(j)
 *          may easily over- or underflow, and BETA(j) may even be zero.
 *          Thus, the user should avoid naively computing the ratio
 *          alpha/beta.  However, ALPHAR and ALPHAI will be always less
 *          than and usually comparable with norm(A) in magnitude, and
 *          BETA always less than and usually comparable with norm(B).
 *
 *  VL      (output) DOUBLE PRECISION array, dimension (LDVL,N)
 *          If JOBVL = 'V', the left eigenvectors u(j) are stored one
 *          after another in the columns of VL, in the same order as
 *          their eigenvalues. If the j-th eigenvalue is real, then
 *          u(j) = VL(:,j), the j-th column of VL. If the j-th and
 *          (j+1)-th eigenvalues form a complex conjugate pair, then
 *          u(j) = VL(:,j)+i*VL(:,j+1) and u(j+1) = VL(:,j)-i*VL(:,j+1).
 *          Each eigenvector will be scaled so the largest component have
 *          abs(real part)+abs(imag. part)=1.
 *          Not referenced if JOBVL = 'N'.
 *
 *  LDVL    (input) INTEGER
 *          The leading dimension of the matrix VL. LDVL >= 1, and
 *          if JOBVL = 'V', LDVL >= N.
 *
 *  VR      (output) DOUBLE PRECISION array, dimension (LDVR,N)
 *          If JOBVR = 'V', the right eigenvectors v(j) are stored one
 *          after another in the columns of VR, in the same order as
 *          their eigenvalues. If the j-th eigenvalue is real, then
 *          v(j) = VR(:,j), the j-th column of VR. If the j-th and
 *          (j+1)-th eigenvalues form a complex conjugate pair, then
 *          v(j) = VR(:,j)+i*VR(:,j+1) and v(j+1) = VR(:,j)-i*VR(:,j+1).
 *          Each eigenvector will be scaled so the largest component have
 *          abs(real part)+abs(imag. part)=1.
 *          Not referenced if JOBVR = 'N'.
 *
 *  LDVR    (input) INTEGER
 *          The leading dimension of the matrix VR. LDVR >= 1, and
 *          if JOBVR = 'V', LDVR >= N.
 *
 *  WORK    (workspace/output) DOUBLE PRECISION array, dimension (LWORK)
 *          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
 *
 *  LWORK   (input) INTEGER
 *          The dimension of the array WORK.  LWORK >= max(1,8*N).
 *          For good performance, LWORK must generally be larger.
 *
 *          If LWORK = -1, then a workspace query is assumed; the routine
 *          only calculates the optimal size of the WORK array, returns
 *          this value as the first entry of the WORK array, and no error
 *          message related to LWORK is issued by XERBLA.
 *
 *  INFO    (output) INTEGER
 *          = 0:  successful exit
 *          < 0:  if INFO = -i, the i-th argument had an illegal value.
 *          = 1,...,N:
 *                The QZ iteration failed.  No eigenvectors have been
 *                calculated, but ALPHAR(j), ALPHAI(j), and BETA(j)
 *                should be correct for j=INFO+1,...,N.
 *          > N:  =N+1: other than QZ iteration failed in DHGEQZ.
 *                =N+2: error return from DTGEVC.
 *
 */

extern "C" void dggev_(char& JOBVL,char& JOBVR,int& N,double* A,int& LDA,
		       double* B,int& LDB,double* ALPHAR,double* ALPHAI,
		       double* BETA,double* VL,int& LDVL,double* VR,int& LDVR,
		       double* WORK,int& LWORK,int& INFO);

/************************************************************************/
/*     SUBROUTINE DSYEVD( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, IWORK,
 *     $                   LIWORK, INFO )
 *
 *  -- LAPACK driver routine (version 3.1) --
 *     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
 *     November 2006
 *
 *     .. Scalar Arguments ..
 *      CHARACTER          JOBZ, UPLO
 *      INTEGER            INFO, LDA, LIWORK, LWORK, N
 *     ..
 *     .. Array Arguments ..
 *      INTEGER            IWORK( * )
 *      DOUBLE PRECISION   A( LDA, * ), W( * ), WORK( * )
 *     ..
 *
 *  Purpose
 *  =======
 *
 *  DSYEVD computes all eigenvalues and, optionally, eigenvectors of a
 *  real symmetric matrix A. If eigenvectors are desired, it uses a
 *  divide and conquer algorithm.
 *
 *  The divide and conquer algorithm makes very mild assumptions about
 *  floating point arithmetic. It will work on machines with a guard
 *  digit in add/subtract, or on those binary machines without guard
 *  digits which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or
 *  Cray-2. It could conceivably fail on hexadecimal or decimal machines
 *  without guard digits, but we know of none.
 *
 *  Because of large use of BLAS of level 3, DSYEVD needs N**2 more
 *  workspace than DSYEVX.
 *
 *  Arguments
 *  =========
 *
 *  JOBZ    (input) CHARACTER*1
 *          = 'N':  Compute eigenvalues only;
 *          = 'V':  Compute eigenvalues and eigenvectors.
 *
 *  UPLO    (input) CHARACTER*1
 *          = 'U':  Upper triangle of A is stored;
 *          = 'L':  Lower triangle of A is stored.
 *
 *  N       (input) INTEGER
 *          The order of the matrix A.  N >= 0.
 *
 *  A       (input/output) DOUBLE PRECISION array, dimension (LDA, N)
 *          On entry, the symmetric matrix A.  If UPLO = 'U', the
 *          leading N-by-N upper triangular part of A contains the
 *          upper triangular part of the matrix A.  If UPLO = 'L',
 *          the leading N-by-N lower triangular part of A contains
 *          the lower triangular part of the matrix A.
 *          On exit, if JOBZ = 'V', then if INFO = 0, A contains the
 *          orthonormal eigenvectors of the matrix A.
 *          If JOBZ = 'N', then on exit the lower triangle (if UPLO='L')
 *          or the upper triangle (if UPLO='U') of A, including the
 *          diagonal, is destroyed.
 *
 *  LDA     (input) INTEGER
 *          The leading dimension of the array A.  LDA >= max(1,N).
 *
 *  W       (output) DOUBLE PRECISION array, dimension (N)
 *          If INFO = 0, the eigenvalues in ascending order.
 *
 *  WORK    (workspace/output) DOUBLE PRECISION array,
 *                                         dimension (LWORK)
 *          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
 *
 *  LWORK   (input) INTEGER
 *          The dimension of the array WORK.
 *          If N <= 1,               LWORK must be at least 1.
 *          If JOBZ = 'N' and N > 1, LWORK must be at least 2*N+1.
 *          If JOBZ = 'V' and N > 1, LWORK must be at least
 *                                                1 + 6*N + 2*N**2.
 *
 *          If LWORK = -1, then a workspace query is assumed; the routine
 *          only calculates the optimal sizes of the WORK and IWORK
 *          arrays, returns these values as the first entries of the WORK
 *          and IWORK arrays, and no error message related to LWORK or
 *          LIWORK is issued by XERBLA.
 *
 *  IWORK   (workspace/output) INTEGER array, dimension (MAX(1,LIWORK))
 *          On exit, if INFO = 0, IWORK(1) returns the optimal LIWORK.
 *
 *  LIWORK  (input) INTEGER
 *          The dimension of the array IWORK.
 *          If N <= 1,                LIWORK must be at least 1.
 *          If JOBZ  = 'N' and N > 1, LIWORK must be at least 1.
 *          If JOBZ  = 'V' and N > 1, LIWORK must be at least 3 + 5*N.
 *
 *          If LIWORK = -1, then a workspace query is assumed; the
 *          routine only calculates the optimal sizes of the WORK and
 *          IWORK arrays, returns these values as the first entries of
 *          the WORK and IWORK arrays, and no error message related to
 *          LWORK or LIWORK is issued by XERBLA.
 *
 *  INFO    (output) INTEGER
 *          = 0:  successful exit
 *          < 0:  if INFO = -i, the i-th argument had an illegal value
 *          > 0:  if INFO = i and JOBZ = 'N', then the algorithm failed
 *                to converge; i off-diagonal elements of an intermediate
 *                tridiagonal form did not converge to zero;
 *                if INFO = i and JOBZ = 'V', then the algorithm failed
 *                to compute an eigenvalue while working on the submatrix
 *                lying in rows and columns INFO/(N+1) through
 *                mod(INFO,N+1).
 */

extern "C" void dsyevd_(char& JOBZ,char& UPLO,int& N,double* A,int& LDA,
			double* W,double* WORK,int& LWORK,int* IWORK,
			int& LIWORK,int& INFO);


/************************************************************************/
/* SUBROUTINE DGESVD( JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT,
 *                    WORK, LWORK, INFO )
 *
 *
 *  -- LAPACK driver routine (version 3.1) --
 *     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
 *     November 2006
 *
 *     .. Scalar Arguments ..
 *     CHARACTER          JOBU, JOBVT
 *     INTEGER            INFO, LDA, LDU, LDVT, LWORK, M, N
 *     ..
 *     .. Array Arguments ..
 *     REAL               A( LDA, * ), S( * ), U( LDU, * ),
 *                        VT( LDVT, * ), WORK( * )
 *     ..
 *
 *  Purpose
 *  =======
 *
 *  SGESVD computes the singular value decomposition (SVD) of a real
 *  M-by-N matrix A, optionally computing the left and/or right singular
 *  vectors. The SVD is written
 *
 *       A = U * SIGMA * transpose(V)
 *
 *  where SIGMA is an M-by-N matrix which is zero except for its
 *  min(m,n) diagonal elements, U is an M-by-M orthogonal matrix, and
 *  V is an N-by-N orthogonal matrix.  The diagonal elements of SIGMA
 *  are the singular values of A; they are real and non-negative, and
 *  are returned in descending order.  The first min(m,n) columns of
 *  U and V are the left and right singular vectors of A.
 *
 *  Note that the routine returns V**T, not V.
 *
 *  Arguments
 *  =========
 *
 *  JOBU    (input) CHARACTER*1
 *          Specifies options for computing all or part of the matrix U:
 *          = 'A':  all M columns of U are returned in array U:
 *          = 'S':  the first min(m,n) columns of U (the left singular
 *                  vectors) are returned in the array U;
 *          = 'O':  the first min(m,n) columns of U (the left singular
 *                  vectors) are overwritten on the array A;
 *          = 'N':  no columns of U (no left singular vectors) are
 *                  computed.
 *
 *  JOBVT   (input) CHARACTER*1
 *          Specifies options for computing all or part of the matrix
 *          V**T:
 *          = 'A':  all N rows of V**T are returned in the array VT;
 *          = 'S':  the first min(m,n) rows of V**T (the right singular
 *                  vectors) are returned in the array VT;
 *          = 'O':  the first min(m,n) rows of V**T (the right singular
 *                  vectors) are overwritten on the array A;
 *          = 'N':  no rows of V**T (no right singular vectors) are
 *                  computed.
 *
 *          JOBVT and JOBU cannot both be 'O'.
 *
 *  M       (input) INTEGER
 *          The number of rows of the input matrix A.  M >= 0.
 *
 *  N       (input) INTEGER
 *          The number of columns of the input matrix A.  N >= 0.
 *
 *  A       (input/output) REAL array, dimension (LDA,N)
 *          On entry, the M-by-N matrix A.
 *          On exit,
 *          if JOBU = 'O',  A is overwritten with the first min(m,n)
 *                          columns of U (the left singular vectors,
 *                          stored columnwise);
 *          if JOBVT = 'O', A is overwritten with the first min(m,n)
 *                          rows of V**T (the right singular vectors,
 *                          stored rowwise);
 *          if JOBU .ne. 'O' and JOBVT .ne. 'O', the contents of A
 *                          are destroyed.
 *
 *  LDA     (input) INTEGER
 *          The leading dimension of the array A.  LDA >= max(1,M).
 *
 *  S       (output) REAL array, dimension (min(M,N))
 *          The singular values of A, sorted so that S(i) >= S(i+1).
 *
 *  U       (output) REAL array, dimension (LDU,UCOL)
 *          (LDU,M) if JOBU = 'A' or (LDU,min(M,N)) if JOBU = 'S'.
 *          If JOBU = 'A', U contains the M-by-M orthogonal matrix U;
 *          if JOBU = 'S', U contains the first min(m,n) columns of U
 *          (the left singular vectors, stored columnwise);
 *          if JOBU = 'N' or 'O', U is not referenced.
 *
 *  LDU     (input) INTEGER
 *          The leading dimension of the array U.  LDU >= 1; if
 *          JOBU = 'S' or 'A', LDU >= M.
 *
 *  VT      (output) REAL array, dimension (LDVT,N)
 *          If JOBVT = 'A', VT contains the N-by-N orthogonal matrix
 *          V**T;
 *          if JOBVT = 'S', VT contains the first min(m,n) rows of
 *          V**T (the right singular vectors, stored rowwise);
 *          if JOBVT = 'N' or 'O', VT is not referenced.
 *
 *  LDVT    (input) INTEGER
 *          The leading dimension of the array VT.  LDVT >= 1; if
 *          JOBVT = 'A', LDVT >= N; if JOBVT = 'S', LDVT >= min(M,N).
 *
 *  WORK    (workspace/output) REAL array, dimension (MAX(1,LWORK))
 *          On exit, if INFO = 0, WORK(1) returns the optimal LWORK;
 *          if INFO > 0, WORK(2:MIN(M,N)) contains the unconverged
 *          superdiagonal elements of an upper bidiagonal matrix B
 *          whose diagonal is in S (not necessarily sorted). B
 *          satisfies A = U * B * VT, so it has the same singular values
 *          as A, and singular vectors related by U and VT.
 *
 *  LWORK   (input) INTEGER
 *          The dimension of the array WORK.
 *          LWORK >= MAX(1,3*MIN(M,N)+MAX(M,N),5*MIN(M,N)).
 *          For good performance, LWORK should generally be larger.
 *
 *          If LWORK = -1, then a workspace query is assumed; the routine
 *          only calculates the optimal size of the WORK array, returns
 *          this value as the first entry of the WORK array, and no error
 *          message related to LWORK is issued by XERBLA.
 *
 *  INFO    (output) INTEGER
 *          = 0:  successful exit.
 *          < 0:  if INFO = -i, the i-th argument had an illegal value.
 *          > 0:  if SBDSQR did not converge, INFO specifies how many
 *                superdiagonals of an intermediate bidiagonal form B
 *                did not converge to zero. See the description of WORK
 *                above for details.
 *
 */

extern "C" void dgesvd_(char& JOBU,char& JOBVT,int& M,int& N,double* A,
			int& LDA,double* S,double* U,int& LDU,
			double* VT,int& LDVT,double* WORK,int& LWORK,
			int& INFO);

#endif
