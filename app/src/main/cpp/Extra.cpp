//
// Created by w1 on 6/6/19.
//



#ifndef FINTEGER
#define FINTEGER long int
#endif
#include "f2c.h"

extern "C" {
/*
typedef long int integer;
typedef unsigned long int uinteger;
typedef char *address;
typedef short int shortint;
typedef float real;
typedef double doublereal;
typedef struct { real r, i; } complex;
typedef struct { doublereal r, i; } doublecomplex;
typedef long int logical;
typedef short int shortlogical;
typedef char logical1;
typedef char integer1;

typedef long int flag;
typedef long int ftnlen;
typedef long int ftnint;

#define abs(x) ((x) >= 0 ? (x) : -(x))
#define dabs(x) (doublereal)abs(x)
#define min(a,b) ((a) <= (b) ? (a) : (b))
#define max(a,b) ((a) >= (b) ? (a) : (b))
#define dmin(a,b) (doublereal)min(a,b)
#define dmax(a,b) (doublereal)max(a,b)
#define bit_test(a,b)	((a) >> (b) & 1)
#define bit_clear(a,b)	((a) & ~((uinteger)1 << (b)))
#define bit_set(a,b)	((a) |  ((uinteger)1 << (b)))
*/

extern /* Subroutine */ int sorg2r_(integer *, integer *, integer *, real
*, integer *, real *, real *, integer *);

extern /* Subroutine */ int sgeqr2_(integer *, integer *, real *, integer
*, real *, real *, integer *);
extern /* Subroutine */ int slarfb_(char *, char *, char *, char *,
                                    integer *, integer *, integer *, real *, integer *, real *,
                                    integer *, real *, integer *, real *, integer *), xerbla_(char *, integer *);
extern integer ilaenv_(integer *, char *, char *, integer *, integer *,
                       integer *, integer *, ftnlen, ftnlen);
extern /* Subroutine */ int slarft_(char *, char *, integer *, integer *,
                                    real *, integer *, real *, real *, integer *);

extern /* Subroutine */ int xerbla_(char *, integer *);



int sgeqrf_ (FINTEGER *m, FINTEGER *n, float *a, FINTEGER *lda,
             float *tau, float *work, FINTEGER *lwork, FINTEGER *info)
{
/*  -- LAPACK routine (version 3.0) --
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
       Courant Institute, Argonne National Lab, and Rice University
       June 30, 1999


    Purpose
    =======

    SGEQRF computes a QR factorization of a real M-by-N matrix A:
    A = Q * R.

    Arguments
    =========

    M       (input) INTEGER
            The number of rows of the matrix A.  M >= 0.

    N       (input) INTEGER
            The number of columns of the matrix A.  N >= 0.

    A       (input/output) REAL array, dimension (LDA,N)
            On entry, the M-by-N matrix A.
            On exit, the elements on and above the diagonal of the array
            contain the min(M,N)-by-N upper trapezoidal matrix R (R is
            upper triangular if m >= n); the elements below the diagonal,
            with the array TAU, represent the orthogonal matrix Q as a
            product of min(m,n) elementary reflectors (see Further
            Details).

    LDA     (input) INTEGER
            The leading dimension of the array A.  LDA >= max(1,M).

    TAU     (output) REAL array, dimension (min(M,N))
            The scalar factors of the elementary reflectors (see Further
            Details).

    WORK    (workspace/output) REAL array, dimension (LWORK)
            On exit, if INFO = 0, WORK(1) returns the optimal LWORK.

    LWORK   (input) INTEGER
            The dimension of the array WORK.  LWORK >= max(1,N).
            For optimum performance LWORK >= N*NB, where NB is
            the optimal blocksize.

            If LWORK = -1, then a workspace query is assumed; the routine
            only calculates the optimal size of the WORK array, returns
            this value as the first entry of the WORK array, and no error
            message related to LWORK is issued by XERBLA.

    INFO    (output) INTEGER
            = 0:  successful exit
            < 0:  if INFO = -i, the i-th argument had an illegal value

    Further Details
    ===============

    The matrix Q is represented as a product of elementary reflectors

       Q = H(1) H(2) . . . H(k), where k = min(m,n).

    Each H(i) has the form

       H(i) = I - tau * v * v'

    where tau is a real scalar, and v is a real vector with
    v(1:i-1) = 0 and v(i) = 1; v(i+1:m) is stored on exit in A(i+1:m,i),
    and tau in TAU(i).

    =====================================================================


       Test the input arguments

       Parameter adjustments */
/* Table of constant values */
static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__3 = 3;
static integer c__2 = 2;

/* System generated locals */
integer a_dim1, a_offset, i__1, i__2, i__3, i__4;
/* Local variables */
static integer i__, k, nbmin, iinfo;
    static integer ib, nb, nx;

static integer ldwork, lwkopt;
static logical lquery;
static integer iws;
#define a_ref(a_1,a_2) a[(a_2)*a_dim1 + a_1]


a_dim1 = *lda;
a_offset = 1 + a_dim1 * 1;
a -= a_offset;
--tau;
--work;

/* Function Body */
*info = 0;
nb = ilaenv_(&c__1, "SGEQRF", " ", m, n, &c_n1, &c_n1, (ftnlen)6, (ftnlen)
        1);
lwkopt = *n * nb;
work[1] = (real) lwkopt;
lquery = *lwork == -1;
if (*m < 0) {
*info = -1;
} else if (*n < 0) {
*info = -2;
} else if (*lda < max(1,*m)) {
*info = -4;
} else if (*lwork < max(1,*n) && ! lquery) {
*info = -7;
}
if (*info != 0) {
i__1 = -(*info);
xerbla_("SGEQRF", &i__1);
return 0;
} else if (lquery) {
return 0;
}

/*     Quick return if possible */

k = min(*m,*n);
if (k == 0) {
work[1] = 1.f;
return 0;
}

nbmin = 2;
nx = 0;
iws = *n;
if (nb > 1 && nb < k) {

/*        Determine when to cross over from blocked to unblocked code.

   Computing MAX */
i__1 = 0, i__2 = ilaenv_(&c__3, "SGEQRF", " ", m, n, &c_n1, &c_n1, (
        ftnlen)6, (ftnlen)1);
nx = max(i__1,i__2);
if (nx < k) {

/*           Determine if workspace is large enough for blocked code. */

ldwork = *n;
iws = ldwork * nb;
if (*lwork < iws) {

/*              Not enough workspace to use optimal NB:  reduce NB and
                determine the minimum value of NB. */

nb = *lwork / ldwork;
/* Computing MAX */
i__1 = 2, i__2 = ilaenv_(&c__2, "SGEQRF", " ", m, n, &c_n1, &
        c_n1, (ftnlen)6, (ftnlen)1);
nbmin = max(i__1,i__2);
}
}
}

if (nb >= nbmin && nb < k && nx < k) {

/*        Use blocked code initially */

i__1 = k - nx;
i__2 = nb;
for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
/* Computing MIN */
i__3 = k - i__ + 1;
ib = min(i__3,nb);

/*           Compute the QR factorization of the current block
             A(i:m,i:i+ib-1) */

i__3 = *m - i__ + 1;
sgeqr2_(&i__3, &ib, &a_ref(i__, i__), lda, &tau[i__], &work[1], &
        iinfo);
if (i__ + ib <= *n) {

/*              Form the triangular factor of the block reflector
                H = H(i) H(i+1) . . . H(i+ib-1) */

i__3 = *m - i__ + 1;
slarft_("Forward", "Columnwise", &i__3, &ib, &a_ref(i__, i__),
lda, &tau[i__], &work[1], &ldwork);

/*              Apply H' to A(i:m,i+ib:n) from the left */

i__3 = *m - i__ + 1;
i__4 = *n - i__ - ib + 1;
slarfb_("Left", "Transpose", "Forward", "Columnwise", &i__3, &
i__4, &ib, &a_ref(i__, i__), lda, &work[1], &ldwork, &
a_ref(i__, i__ + ib), lda, &work[ib + 1], &ldwork);
}
/* L10: */
}
} else {
i__ = 1;
}

/*     Use unblocked code to factor the last or only block. */

if (i__ <= k) {
i__2 = *m - i__ + 1;
i__1 = *n - i__ + 1;
sgeqr2_(&i__2, &i__1, &a_ref(i__, i__), lda, &tau[i__], &work[1], &
        iinfo);
}

work[1] = (real) iws;
return 0;

/*     End of SGEQRF */

}






/* Subroutine */ int sorgqr_(integer *m, integer *n, integer *k, real *a,
                             integer *lda, real *tau, real *work, integer *lwork, integer *info)
{
/*  -- LAPACK routine (version 3.0) --
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
       Courant Institute, Argonne National Lab, and Rice University
       June 30, 1999


    Purpose
    =======

    SORGQR generates an M-by-N real matrix Q with orthonormal columns,
    which is defined as the first N columns of a product of K elementary
    reflectors of order M

          Q  =  H(1) H(2) . . . H(k)

    as returned by SGEQRF.

    Arguments
    =========

    M       (input) INTEGER
            The number of rows of the matrix Q. M >= 0.

    N       (input) INTEGER
            The number of columns of the matrix Q. M >= N >= 0.

    K       (input) INTEGER
            The number of elementary reflectors whose product defines the
            matrix Q. N >= K >= 0.

    A       (input/output) REAL array, dimension (LDA,N)
            On entry, the i-th column must contain the vector which
            defines the elementary reflector H(i), for i = 1,2,...,k, as
            returned by SGEQRF in the first k columns of its array
            argument A.
            On exit, the M-by-N matrix Q.

    LDA     (input) INTEGER
            The first dimension of the array A. LDA >= max(1,M).

    TAU     (input) REAL array, dimension (K)
            TAU(i) must contain the scalar factor of the elementary
            reflector H(i), as returned by SGEQRF.

    WORK    (workspace/output) REAL array, dimension (LWORK)
            On exit, if INFO = 0, WORK(1) returns the optimal LWORK.

    LWORK   (input) INTEGER
            The dimension of the array WORK. LWORK >= max(1,N).
            For optimum performance LWORK >= N*NB, where NB is the
            optimal blocksize.

            If LWORK = -1, then a workspace query is assumed; the routine
            only calculates the optimal size of the WORK array, returns
            this value as the first entry of the WORK array, and no error
            message related to LWORK is issued by XERBLA.

    INFO    (output) INTEGER
            = 0:  successful exit
            < 0:  if INFO = -i, the i-th argument has an illegal value

    =====================================================================


       Test the input arguments

       Parameter adjustments */
    /* Table of constant values */
    static integer c__1 = 1;
    static integer c_n1 = -1;
    static integer c__3 = 3;
    static integer c__2 = 2;

    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;
    /* Local variables */
    static integer i__, j, l, nbmin, iinfo, ib;

    static integer nb, ki, kk, nx;

    static integer ldwork, lwkopt;
    static logical lquery;
    static integer iws;
#define a_ref(a_1,a_2) a[(a_2)*a_dim1 + a_1]


    a_dim1 = *lda;
    a_offset = 1 + a_dim1 * 1;
    a -= a_offset;
    --tau;
    --work;

    /* Function Body */
    *info = 0;
    nb = ilaenv_(&c__1, "SORGQR", " ", m, n, k, &c_n1, (ftnlen)6, (ftnlen)1);
    lwkopt = max(1,*n) * nb;
    work[1] = (real) lwkopt;
    lquery = *lwork == -1;
    if (*m < 0) {
        *info = -1;
    } else if (*n < 0 || *n > *m) {
        *info = -2;
    } else if (*k < 0 || *k > *n) {
        *info = -3;
    } else if (*lda < max(1,*m)) {
        *info = -5;
    } else if (*lwork < max(1,*n) && ! lquery) {
        *info = -8;
    }
    if (*info != 0) {
        i__1 = -(*info);
        xerbla_("SORGQR", &i__1);
        return 0;
    } else if (lquery) {
        return 0;
    }

/*     Quick return if possible */

    if (*n <= 0) {
        work[1] = 1.f;
        return 0;
    }

    nbmin = 2;
    nx = 0;
    iws = *n;
    if (nb > 1 && nb < *k) {

/*        Determine when to cross over from blocked to unblocked code.

   Computing MAX */
        i__1 = 0, i__2 = ilaenv_(&c__3, "SORGQR", " ", m, n, k, &c_n1, (
                ftnlen)6, (ftnlen)1);
        nx = max(i__1,i__2);
        if (nx < *k) {

/*           Determine if workspace is large enough for blocked code. */

            ldwork = *n;
            iws = ldwork * nb;
            if (*lwork < iws) {

/*              Not enough workspace to use optimal NB:  reduce NB and
                determine the minimum value of NB. */

                nb = *lwork / ldwork;
/* Computing MAX */
                i__1 = 2, i__2 = ilaenv_(&c__2, "SORGQR", " ", m, n, k, &c_n1,
                                         (ftnlen)6, (ftnlen)1);
                nbmin = max(i__1,i__2);
            }
        }
    }

    if (nb >= nbmin && nb < *k && nx < *k) {

/*        Use blocked code after the last block.
          The first kk columns are handled by the block method. */

        ki = (*k - nx - 1) / nb * nb;
/* Computing MIN */
        i__1 = *k, i__2 = ki + nb;
        kk = min(i__1,i__2);

/*        Set A(1:kk,kk+1:n) to zero. */

        i__1 = *n;
        for (j = kk + 1; j <= i__1; ++j) {
            i__2 = kk;
            for (i__ = 1; i__ <= i__2; ++i__) {
                a_ref(i__, j) = 0.f;
/* L10: */
            }
/* L20: */
        }
    } else {
        kk = 0;
    }

/*     Use unblocked code for the last or only block. */

    if (kk < *n) {
        i__1 = *m - kk;
        i__2 = *n - kk;
        i__3 = *k - kk;
        sorg2r_(&i__1, &i__2, &i__3, &a_ref(kk + 1, kk + 1), lda, &tau[kk + 1]
                , &work[1], &iinfo);
    }

    if (kk > 0) {

/*        Use blocked code */

        i__1 = -nb;
        for (i__ = ki + 1; i__1 < 0 ? i__ >= 1 : i__ <= 1; i__ += i__1) {
/* Computing MIN */
            i__2 = nb, i__3 = *k - i__ + 1;
            ib = min(i__2,i__3);
            if (i__ + ib <= *n) {

/*              Form the triangular factor of the block reflector
                H = H(i) H(i+1) . . . H(i+ib-1) */

                i__2 = *m - i__ + 1;
                slarft_("Forward", "Columnwise", &i__2, &ib, &a_ref(i__, i__),
                        lda, &tau[i__], &work[1], &ldwork);

/*              Apply H to A(i:m,i+ib:n) from the left */

                i__2 = *m - i__ + 1;
                i__3 = *n - i__ - ib + 1;
                slarfb_("Left", "No transpose", "Forward", "Columnwise", &
                        i__2, &i__3, &ib, &a_ref(i__, i__), lda, &work[1], &
                                ldwork, &a_ref(i__, i__ + ib), lda, &work[ib + 1], &
                                ldwork);
            }

/*           Apply H to rows i:m of current block */

            i__2 = *m - i__ + 1;
            sorg2r_(&i__2, &ib, &ib, &a_ref(i__, i__), lda, &tau[i__], &work[
                    1], &iinfo);

/*           Set rows 1:i-1 of current block to zero */

            i__2 = i__ + ib - 1;
            for (j = i__; j <= i__2; ++j) {
                i__3 = i__ - 1;
                for (l = 1; l <= i__3; ++l) {
                    a_ref(l, j) = 0.f;
/* L30: */
                }
/* L40: */
            }
/* L50: */
        }
    }

    work[1] = (real) iws;
    return 0;

/*     End of SORGQR */

} /* sorgqr_ */

#undef a_ref


/* Subroutine */ int ssyrk_(char *uplo, char *trans, integer *n, integer *k,
                            real *alpha, real *a, integer *lda, real *beta, real *c, integer *ldc)
{


    /* System generated locals */
    integer a_dim1, a_offset, c_dim1, c_offset, i__1, i__2, i__3;

    /* Local variables */
    static integer info;
    static real temp;
    static integer i, j, l;
    extern logical lsame_(char *, char *);
    static integer nrowa;
    static logical upper;



/*  Purpose
    =======

    SSYRK  performs one of the symmetric rank k operations

       C := alpha*A*A' + beta*C,

    or

       C := alpha*A'*A + beta*C,

    where  alpha and beta  are scalars, C is an  n by n  symmetric matrix

    and  A  is an  n by k  matrix in the first case and a  k by n  matrix

    in the second case.

    Parameters
    ==========

    UPLO   - CHARACTER*1.
             On  entry,   UPLO  specifies  whether  the  upper  or  lower

             triangular  part  of the  array  C  is to be  referenced  as

             follows:

                UPLO = 'U' or 'u'   Only the  upper triangular part of  C

                                    is to be referenced.

                UPLO = 'L' or 'l'   Only the  lower triangular part of  C

                                    is to be referenced.

             Unchanged on exit.

    TRANS  - CHARACTER*1.
             On entry,  TRANS  specifies the operation to be performed as

             follows:

                TRANS = 'N' or 'n'   C := alpha*A*A' + beta*C.

                TRANS = 'T' or 't'   C := alpha*A'*A + beta*C.

                TRANS = 'C' or 'c'   C := alpha*A'*A + beta*C.

             Unchanged on exit.

    N      - INTEGER.
             On entry,  N specifies the order of the matrix C.  N must be

             at least zero.
             Unchanged on exit.

    K      - INTEGER.
             On entry with  TRANS = 'N' or 'n',  K  specifies  the number

             of  columns   of  the   matrix   A,   and  on   entry   with

             TRANS = 'T' or 't' or 'C' or 'c',  K  specifies  the  number

             of rows of the matrix  A.  K must be at least zero.
             Unchanged on exit.

    ALPHA  - REAL            .
             On entry, ALPHA specifies the scalar alpha.
             Unchanged on exit.

    A      - REAL             array of DIMENSION ( LDA, ka ), where ka is

             k  when  TRANS = 'N' or 'n',  and is  n  otherwise.
             Before entry with  TRANS = 'N' or 'n',  the  leading  n by k

             part of the array  A  must contain the matrix  A,  otherwise

             the leading  k by n  part of the array  A  must contain  the

             matrix A.
             Unchanged on exit.

    LDA    - INTEGER.
             On entry, LDA specifies the first dimension of A as declared

             in  the  calling  (sub)  program.   When  TRANS = 'N' or 'n'

             then  LDA must be at least  max( 1, n ), otherwise  LDA must

             be at least  max( 1, k ).
             Unchanged on exit.

    BETA   - REAL            .
             On entry, BETA specifies the scalar beta.
             Unchanged on exit.

    C      - REAL             array of DIMENSION ( LDC, n ).
             Before entry  with  UPLO = 'U' or 'u',  the leading  n by n

             upper triangular part of the array C must contain the upper

             triangular part  of the  symmetric matrix  and the strictly

             lower triangular part of C is not referenced.  On exit, the

             upper triangular part of the array  C is overwritten by the

             upper triangular part of the updated matrix.
             Before entry  with  UPLO = 'L' or 'l',  the leading  n by n

             lower triangular part of the array C must contain the lower

             triangular part  of the  symmetric matrix  and the strictly

             upper triangular part of C is not referenced.  On exit, the

             lower triangular part of the array  C is overwritten by the

             lower triangular part of the updated matrix.

    LDC    - INTEGER.
             On entry, LDC specifies the first dimension of C as declared

             in  the  calling  (sub)  program.   LDC  must  be  at  least

             max( 1, n ).
             Unchanged on exit.


    Level 3 Blas routine.

    -- Written on 8-February-1989.
       Jack Dongarra, Argonne National Laboratory.
       Iain Duff, AERE Harwell.
       Jeremy Du Croz, Numerical Algorithms Group Ltd.
       Sven Hammarling, Numerical Algorithms Group Ltd.



       Test the input parameters.


   Parameter adjustments
       Function Body */

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]
#define C(I,J) c[(I)-1 + ((J)-1)* ( *ldc)]

    if (lsame_(trans, "N")) {
        nrowa = *n;
    } else {
        nrowa = *k;
    }
    upper = lsame_(uplo, "U");

    info = 0;
    if (! upper && ! lsame_(uplo, "L")) {
        info = 1;
    } else if (! lsame_(trans, "N") && ! lsame_(trans, "T") &&
               ! lsame_(trans, "C")) {
        info = 2;
    } else if (*n < 0) {
        info = 3;
    } else if (*k < 0) {
        info = 4;
    } else if (*lda < max(1,nrowa)) {
        info = 7;
    } else if (*ldc < max(1,*n)) {
        info = 10;
    }
    if (info != 0) {
        xerbla_("SSYRK ", &info);
        return 0;
    }

/*     Quick return if possible. */

    if (*n == 0 || (*alpha == 0.f || *k == 0) && *beta == 1.f) {
        return 0;
    }

/*     And when  alpha.eq.zero. */

    if (*alpha == 0.f) {
        if (upper) {
            if (*beta == 0.f) {
                i__1 = *n;
                for (j = 1; j <= *n; ++j) {
                    i__2 = j;
                    for (i = 1; i <= j; ++i) {
                        C(i,j) = 0.f;
/* L10: */
                    }
/* L20: */
                }
            } else {
                i__1 = *n;
                for (j = 1; j <= *n; ++j) {
                    i__2 = j;
                    for (i = 1; i <= j; ++i) {
                        C(i,j) = *beta * C(i,j);
/* L30: */
                    }
/* L40: */
                }
            }
        } else {
            if (*beta == 0.f) {
                i__1 = *n;
                for (j = 1; j <= *n; ++j) {
                    i__2 = *n;
                    for (i = j; i <= *n; ++i) {
                        C(i,j) = 0.f;
/* L50: */
                    }
/* L60: */
                }
            } else {
                i__1 = *n;
                for (j = 1; j <= *n; ++j) {
                    i__2 = *n;
                    for (i = j; i <= *n; ++i) {
                        C(i,j) = *beta * C(i,j);
/* L70: */
                    }
/* L80: */
                }
            }
        }
        return 0;
    }

/*     Start the operations. */

    if (lsame_(trans, "N")) {

/*        Form  C := alpha*A*A' + beta*C. */

        if (upper) {
            i__1 = *n;
            for (j = 1; j <= *n; ++j) {
                if (*beta == 0.f) {
                    i__2 = j;
                    for (i = 1; i <= j; ++i) {
                        C(i,j) = 0.f;
/* L90: */
                    }
                } else if (*beta != 1.f) {
                    i__2 = j;
                    for (i = 1; i <= j; ++i) {
                        C(i,j) = *beta * C(i,j);
/* L100: */
                    }
                }
                i__2 = *k;
                for (l = 1; l <= *k; ++l) {
                    if (A(j,l) != 0.f) {
                        temp = *alpha * A(j,l);
                        i__3 = j;
                        for (i = 1; i <= j; ++i) {
                            C(i,j) += temp * A(i,l);
/* L110: */
                        }
                    }
/* L120: */
                }
/* L130: */
            }
        } else {
            i__1 = *n;
            for (j = 1; j <= *n; ++j) {
                if (*beta == 0.f) {
                    i__2 = *n;
                    for (i = j; i <= *n; ++i) {
                        C(i,j) = 0.f;
/* L140: */
                    }
                } else if (*beta != 1.f) {
                    i__2 = *n;
                    for (i = j; i <= *n; ++i) {
                        C(i,j) = *beta * C(i,j);
/* L150: */
                    }
                }
                i__2 = *k;
                for (l = 1; l <= *k; ++l) {
                    if (A(j,l) != 0.f) {
                        temp = *alpha * A(j,l);
                        i__3 = *n;
                        for (i = j; i <= *n; ++i) {
                            C(i,j) += temp * A(i,l);
/* L160: */
                        }
                    }
/* L170: */
                }
/* L180: */
            }
        }
    } else {

/*        Form  C := alpha*A'*A + beta*C. */

        if (upper) {
            i__1 = *n;
            for (j = 1; j <= *n; ++j) {
                i__2 = j;
                for (i = 1; i <= j; ++i) {
                    temp = 0.f;
                    i__3 = *k;
                    for (l = 1; l <= *k; ++l) {
                        temp += A(l,i) * A(l,j);
/* L190: */
                    }
                    if (*beta == 0.f) {
                        C(i,j) = *alpha * temp;
                    } else {
                        C(i,j) = *alpha * temp + *beta * C(i,j);
                    }
/* L200: */
                }
/* L210: */
            }
        } else {
            i__1 = *n;
            for (j = 1; j <= *n; ++j) {
                i__2 = *n;
                for (i = j; i <= *n; ++i) {
                    temp = 0.f;
                    i__3 = *k;
                    for (l = 1; l <= *k; ++l) {
                        temp += A(l,i) * A(l,j);
/* L220: */
                    }
                    if (*beta == 0.f) {
                        C(i,j) = *alpha * temp;
                    } else {
                        C(i,j) = *alpha * temp + *beta * C(i,j);
                    }
/* L230: */
                }
/* L240: */
            }
        }
    }

    return 0;

/*     End of SSYRK . */

} /* ssyrk_ */

}