/* ch.f -- translated by f2c (version 20000121).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static doublereal c_b43 = 1.;

/* Subroutine */ int ch_(nm, n, ar, ai, w, matz, zr, zi, fv1, fv2, fm1, ierr)
integer *nm, *n;
doublereal *ar, *ai, *w;
integer *matz;
doublereal *zr, *zi, *fv1, *fv2, *fm1;
integer *ierr;
{
    /* System generated locals */
    integer ar_dim1, ar_offset, ai_dim1, ai_offset, zr_dim1, zr_offset, 
	    zi_dim1, zi_offset, i__1, i__2;

    /* Local variables */
    static integer i__, j;
    extern /* Subroutine */ int htridi_(), htribk_(), tqlrat_(), tql2_();



/*     this subroutine calls the recommended sequence of */
/*     subroutines from the eigensystem subroutine package (eispack) */
/*     to find the eigenvalues and eigenvectors (if desired) */
/*     of a complex hermitian matrix. */

/*     on input */

/*        nm  must be set to the row dimension of the two-dimensional */
/*        array parameters as declared in the calling program */
/*        dimension statement. */

/*        n  is the order of the matrix  a=(ar,ai). */

/*        ar  and  ai  contain the real and imaginary parts, */
/*        respectively, of the complex hermitian matrix. */

/*        matz  is an integer variable set equal to zero if */
/*        only eigenvalues are desired.  otherwise it is set to */
/*        any non-zero integer for both eigenvalues and eigenvectors. */

/*     on output */

/*        w  contains the eigenvalues in ascending order. */

/*        zr  and  zi  contain the real and imaginary parts, */
/*        respectively, of the eigenvectors if matz is not zero. */

/*        ierr  is an integer output variable set equal to an error */
/*           completion code described in the documentation for tqlrat */
/*           and tql2.  the normal completion code is zero. */

/*        fv1, fv2, and  fm1  are temporary storage arrays. */

/*     questions and comments should be directed to burton s. garbow, */
/*     mathematics and computer science div, argonne national laboratory */

/*     this version dated august 1983. */

/*     ------------------------------------------------------------------ */

    /* Parameter adjustments */
    fm1 -= 3;
    --fv2;
    --fv1;
    zi_dim1 = *nm;
    zi_offset = 1 + zi_dim1 * 1;
    zi -= zi_offset;
    zr_dim1 = *nm;
    zr_offset = 1 + zr_dim1 * 1;
    zr -= zr_offset;
    --w;
    ai_dim1 = *nm;
    ai_offset = 1 + ai_dim1 * 1;
    ai -= ai_offset;
    ar_dim1 = *nm;
    ar_offset = 1 + ar_dim1 * 1;
    ar -= ar_offset;

    /* Function Body */
    if (*n <= *nm) {
	goto L10;
    }
    *ierr = *n * 10;
    goto L50;

L10:
    htridi_(nm, n, &ar[ar_offset], &ai[ai_offset], &w[1], &fv1[1], &fv2[1], &
	    fm1[3]);
    if (*matz != 0) {
	goto L20;
    }
/*     .......... find eigenvalues only .......... */
    tqlrat_(n, &w[1], &fv2[1], ierr);
    goto L50;
/*     .......... find both eigenvalues and eigenvectors .......... */
L20:
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {

	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
	    zr[j + i__ * zr_dim1] = 0.;
/* L30: */
	}

	zr[i__ + i__ * zr_dim1] = 1.;
/* L40: */
    }

    tql2_(nm, n, &w[1], &fv1[1], &zr[zr_offset], ierr);
    if (*ierr != 0) {
	goto L50;
    }
    htribk_(nm, n, &ar[ar_offset], &ai[ai_offset], &fm1[3], n, &zr[zr_offset],
	     &zi[zi_offset]);
L50:
    return 0;
} /* ch_ */

doublereal epslon_(x)
doublereal *x;
{
    /* System generated locals */
    doublereal ret_val, d__1;

    /* Local variables */
    static doublereal a, b, c__, eps;


/*     estimate unit roundoff in quantities of size x. */


/*     this program should function properly on all systems */
/*     satisfying the following two assumptions, */
/*        1.  the base used in representing floating point */
/*            numbers is not a power of three. */
/*        2.  the quantity  a  in statement 10 is represented to */
/*            the accuracy used in floating point variables */
/*            that are stored in memory. */
/*     the statement number 10 and the go to 10 are intended to */
/*     force optimizing compilers to generate code satisfying */
/*     assumption 2. */
/*     under these assumptions, it should be true that, */
/*            a  is not exactly equal to four-thirds, */
/*            b  has a zero for its last bit or digit, */
/*            c  is not exactly equal to one, */
/*            eps  measures the separation of 1.0 from */
/*                 the next larger floating point number. */
/*     the developers of eispack would appreciate being informed */
/*     about any systems where these assumptions do not hold. */

/*     this version dated 4/6/83. */

    a = 1.3333333333333333;
L10:
    b = a - 1.;
    c__ = b + b + b;
    eps = (d__1 = c__ - 1., abs(d__1));
    if (eps == 0.) {
	goto L10;
    }
    ret_val = eps * abs(*x);
    return ret_val;
} /* epslon_ */

/* Subroutine */ int htribk_(nm, n, ar, ai, tau, m, zr, zi)
integer *nm, *n;
doublereal *ar, *ai, *tau;
integer *m;
doublereal *zr, *zi;
{
    /* System generated locals */
    integer ar_dim1, ar_offset, ai_dim1, ai_offset, zr_dim1, zr_offset, 
	    zi_dim1, zi_offset, i__1, i__2, i__3;

    /* Local variables */
    static doublereal h__;
    static integer i__, j, k, l;
    static doublereal s, si;



/*     this subroutine is a translation of a complex analogue of */
/*     the algol procedure trbak1, num. math. 11, 181-195(1968) */
/*     by martin, reinsch, and wilkinson. */
/*     handbook for auto. comp., vol.ii-linear algebra, 212-226(1971). */

/*     this subroutine forms the eigenvectors of a complex hermitian */
/*     matrix by back transforming those of the corresponding */
/*     real symmetric tridiagonal matrix determined by  htridi. */

/*     on input */

/*        nm must be set to the row dimension of two-dimensional */
/*          array parameters as declared in the calling program */
/*          dimension statement. */

/*        n is the order of the matrix. */

/*        ar and ai contain information about the unitary trans- */
/*          formations used in the reduction by  htridi  in their */
/*          full lower triangles except for the diagonal of ar. */

/*        tau contains further information about the transformations. */

/*        m is the number of eigenvectors to be back transformed. */

/*        zr contains the eigenvectors to be back transformed */
/*          in its first m columns. */

/*     on output */

/*        zr and zi contain the real and imaginary parts, */
/*          respectively, of the transformed eigenvectors */
/*          in their first m columns. */

/*     note that the last component of each returned vector */
/*     is real and that vector euclidean norms are preserved. */

/*     questions and comments should be directed to burton s. garbow, */
/*     mathematics and computer science div, argonne national laboratory */

/*     this version dated august 1983. */

/*     ------------------------------------------------------------------ */

    /* Parameter adjustments */
    tau -= 3;
    ai_dim1 = *nm;
    ai_offset = 1 + ai_dim1 * 1;
    ai -= ai_offset;
    ar_dim1 = *nm;
    ar_offset = 1 + ar_dim1 * 1;
    ar -= ar_offset;
    zi_dim1 = *nm;
    zi_offset = 1 + zi_dim1 * 1;
    zi -= zi_offset;
    zr_dim1 = *nm;
    zr_offset = 1 + zr_dim1 * 1;
    zr -= zr_offset;

    /* Function Body */
    if (*m == 0) {
	goto L200;
    }
/*     .......... transform the eigenvectors of the real symmetric */
/*                tridiagonal matrix to those of the hermitian */
/*                tridiagonal matrix. .......... */
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {

	i__2 = *m;
	for (j = 1; j <= i__2; ++j) {
	    zi[k + j * zi_dim1] = -zr[k + j * zr_dim1] * tau[(k << 1) + 2];
	    zr[k + j * zr_dim1] *= tau[(k << 1) + 1];
/* L50: */
	}
    }

    if (*n == 1) {
	goto L200;
    }
/*     .......... recover and apply the householder matrices .......... */
    i__2 = *n;
    for (i__ = 2; i__ <= i__2; ++i__) {
	l = i__ - 1;
	h__ = ai[i__ + i__ * ai_dim1];
	if (h__ == 0.) {
	    goto L140;
	}

	i__1 = *m;
	for (j = 1; j <= i__1; ++j) {
	    s = 0.;
	    si = 0.;

	    i__3 = l;
	    for (k = 1; k <= i__3; ++k) {
		s = s + ar[i__ + k * ar_dim1] * zr[k + j * zr_dim1] - ai[i__ 
			+ k * ai_dim1] * zi[k + j * zi_dim1];
		si = si + ar[i__ + k * ar_dim1] * zi[k + j * zi_dim1] + ai[
			i__ + k * ai_dim1] * zr[k + j * zr_dim1];
/* L110: */
	    }
/*     .......... double divisions avoid possible underflow .......... */
	    s = s / h__ / h__;
	    si = si / h__ / h__;

	    i__3 = l;
	    for (k = 1; k <= i__3; ++k) {
		zr[k + j * zr_dim1] = zr[k + j * zr_dim1] - s * ar[i__ + k * 
			ar_dim1] - si * ai[i__ + k * ai_dim1];
		zi[k + j * zi_dim1] = zi[k + j * zi_dim1] - si * ar[i__ + k * 
			ar_dim1] + s * ai[i__ + k * ai_dim1];
/* L120: */
	    }

/* L130: */
	}

L140:
	;
    }

L200:
    return 0;
} /* htribk_ */

/* Subroutine */ int htridi_(nm, n, ar, ai, d__, e, e2, tau)
integer *nm, *n;
doublereal *ar, *ai, *d__, *e, *e2, *tau;
{
    /* System generated locals */
    integer ar_dim1, ar_offset, ai_dim1, ai_offset, i__1, i__2, i__3;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt();

    /* Local variables */
    static doublereal f, g, h__;
    static integer i__, j, k, l;
    static doublereal scale, fi, gi, hh;
    static integer ii;
    static doublereal si;
    extern doublereal pythag_();
    static integer jp1;



/*     this subroutine is a translation of a complex analogue of */
/*     the algol procedure tred1, num. math. 11, 181-195(1968) */
/*     by martin, reinsch, and wilkinson. */
/*     handbook for auto. comp., vol.ii-linear algebra, 212-226(1971). */

/*     this subroutine reduces a complex hermitian matrix */
/*     to a real symmetric tridiagonal matrix using */
/*     unitary similarity transformations. */

/*     on input */

/*        nm must be set to the row dimension of two-dimensional */
/*          array parameters as declared in the calling program */
/*          dimension statement. */

/*        n is the order of the matrix. */

/*        ar and ai contain the real and imaginary parts, */
/*          respectively, of the complex hermitian input matrix. */
/*          only the lower triangle of the matrix need be supplied. */

/*     on output */

/*        ar and ai contain information about the unitary trans- */
/*          formations used in the reduction in their full lower */
/*          triangles.  their strict upper triangles and the */
/*          diagonal of ar are unaltered. */

/*        d contains the diagonal elements of the the tridiagonal matrix. */

/*        e contains the subdiagonal elements of the tridiagonal */
/*          matrix in its last n-1 positions.  e(1) is set to zero. */

/*        e2 contains the squares of the corresponding elements of e. */
/*          e2 may coincide with e if the squares are not needed. */

/*        tau contains further information about the transformations. */

/*     calls pythag for  dsqrt(a*a + b*b) . */

/*     questions and comments should be directed to burton s. garbow, */
/*     mathematics and computer science div, argonne national laboratory */

/*     this version dated august 1983. */

/*     ------------------------------------------------------------------ */

    /* Parameter adjustments */
    tau -= 3;
    --e2;
    --e;
    --d__;
    ai_dim1 = *nm;
    ai_offset = 1 + ai_dim1 * 1;
    ai -= ai_offset;
    ar_dim1 = *nm;
    ar_offset = 1 + ar_dim1 * 1;
    ar -= ar_offset;

    /* Function Body */
    tau[(*n << 1) + 1] = 1.;
    tau[(*n << 1) + 2] = 0.;

    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L100: */
	d__[i__] = ar[i__ + i__ * ar_dim1];
    }
/*     .......... for i=n step -1 until 1 do -- .......... */
    i__1 = *n;
    for (ii = 1; ii <= i__1; ++ii) {
	i__ = *n + 1 - ii;
	l = i__ - 1;
	h__ = 0.;
	scale = 0.;
	if (l < 1) {
	    goto L130;
	}
/*     .......... scale row (algol tol then not needed) .......... */
	i__2 = l;
	for (k = 1; k <= i__2; ++k) {
/* L120: */
	    scale = scale + (d__1 = ar[i__ + k * ar_dim1], abs(d__1)) + (d__2 
		    = ai[i__ + k * ai_dim1], abs(d__2));
	}

	if (scale != 0.) {
	    goto L140;
	}
	tau[(l << 1) + 1] = 1.;
	tau[(l << 1) + 2] = 0.;
L130:
	e[i__] = 0.;
	e2[i__] = 0.;
	goto L290;

L140:
	i__2 = l;
	for (k = 1; k <= i__2; ++k) {
	    ar[i__ + k * ar_dim1] /= scale;
	    ai[i__ + k * ai_dim1] /= scale;
	    h__ = h__ + ar[i__ + k * ar_dim1] * ar[i__ + k * ar_dim1] + ai[
		    i__ + k * ai_dim1] * ai[i__ + k * ai_dim1];
/* L150: */
	}

	e2[i__] = scale * scale * h__;
	g = sqrt(h__);
	e[i__] = scale * g;
	f = pythag_(&ar[i__ + l * ar_dim1], &ai[i__ + l * ai_dim1]);
/*     .......... form next diagonal element of matrix t .......... */
	if (f == 0.) {
	    goto L160;
	}
	tau[(l << 1) + 1] = (ai[i__ + l * ai_dim1] * tau[(i__ << 1) + 2] - ar[
		i__ + l * ar_dim1] * tau[(i__ << 1) + 1]) / f;
	si = (ar[i__ + l * ar_dim1] * tau[(i__ << 1) + 2] + ai[i__ + l * 
		ai_dim1] * tau[(i__ << 1) + 1]) / f;
	h__ += f * g;
	g = g / f + 1.;
	ar[i__ + l * ar_dim1] = g * ar[i__ + l * ar_dim1];
	ai[i__ + l * ai_dim1] = g * ai[i__ + l * ai_dim1];
	if (l == 1) {
	    goto L270;
	}
	goto L170;
L160:
	tau[(l << 1) + 1] = -tau[(i__ << 1) + 1];
	si = tau[(i__ << 1) + 2];
	ar[i__ + l * ar_dim1] = g;
L170:
	f = 0.;

	i__2 = l;
	for (j = 1; j <= i__2; ++j) {
	    g = 0.;
	    gi = 0.;
/*     .......... form element of a*u .......... */
	    i__3 = j;
	    for (k = 1; k <= i__3; ++k) {
		g = g + ar[j + k * ar_dim1] * ar[i__ + k * ar_dim1] + ai[j + 
			k * ai_dim1] * ai[i__ + k * ai_dim1];
		gi = gi - ar[j + k * ar_dim1] * ai[i__ + k * ai_dim1] + ai[j 
			+ k * ai_dim1] * ar[i__ + k * ar_dim1];
/* L180: */
	    }

	    jp1 = j + 1;
	    if (l < jp1) {
		goto L220;
	    }

	    i__3 = l;
	    for (k = jp1; k <= i__3; ++k) {
		g = g + ar[k + j * ar_dim1] * ar[i__ + k * ar_dim1] - ai[k + 
			j * ai_dim1] * ai[i__ + k * ai_dim1];
		gi = gi - ar[k + j * ar_dim1] * ai[i__ + k * ai_dim1] - ai[k 
			+ j * ai_dim1] * ar[i__ + k * ar_dim1];
/* L200: */
	    }
/*     .......... form element of p .......... */
L220:
	    e[j] = g / h__;
	    tau[(j << 1) + 2] = gi / h__;
	    f = f + e[j] * ar[i__ + j * ar_dim1] - tau[(j << 1) + 2] * ai[i__ 
		    + j * ai_dim1];
/* L240: */
	}

	hh = f / (h__ + h__);
/*     .......... form reduced a .......... */
	i__2 = l;
	for (j = 1; j <= i__2; ++j) {
	    f = ar[i__ + j * ar_dim1];
	    g = e[j] - hh * f;
	    e[j] = g;
	    fi = -ai[i__ + j * ai_dim1];
	    gi = tau[(j << 1) + 2] - hh * fi;
	    tau[(j << 1) + 2] = -gi;

	    i__3 = j;
	    for (k = 1; k <= i__3; ++k) {
		ar[j + k * ar_dim1] = ar[j + k * ar_dim1] - f * e[k] - g * ar[
			i__ + k * ar_dim1] + fi * tau[(k << 1) + 2] + gi * ai[
			i__ + k * ai_dim1];
		ai[j + k * ai_dim1] = ai[j + k * ai_dim1] - f * tau[(k << 1) 
			+ 2] - g * ai[i__ + k * ai_dim1] - fi * e[k] - gi * 
			ar[i__ + k * ar_dim1];
/* L260: */
	    }
	}

L270:
	i__3 = l;
	for (k = 1; k <= i__3; ++k) {
	    ar[i__ + k * ar_dim1] = scale * ar[i__ + k * ar_dim1];
	    ai[i__ + k * ai_dim1] = scale * ai[i__ + k * ai_dim1];
/* L280: */
	}

	tau[(l << 1) + 2] = -si;
L290:
	hh = d__[i__];
	d__[i__] = ar[i__ + i__ * ar_dim1];
	ar[i__ + i__ * ar_dim1] = hh;
	ai[i__ + i__ * ai_dim1] = scale * sqrt(h__);
/* L300: */
    }

    return 0;
} /* htridi_ */

/* Subroutine */ int tqlrat_(n, d__, e2, ierr)
integer *n;
doublereal *d__, *e2;
integer *ierr;
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(), d_sign();

    /* Local variables */
    static doublereal b, c__, f, g, h__;
    static integer i__, j, l, m;
    static doublereal p, r__, s, t;
    static integer l1, ii;
    extern doublereal pythag_(), epslon_();
    static integer mml;



/*     this subroutine is a translation of the algol procedure tqlrat, */
/*     algorithm 464, comm. acm 16, 689(1973) by reinsch. */

/*     this subroutine finds the eigenvalues of a symmetric */
/*     tridiagonal matrix by the rational ql method. */

/*     on input */

/*        n is the order of the matrix. */

/*        d contains the diagonal elements of the input matrix. */

/*        e2 contains the squares of the subdiagonal elements of the */
/*          input matrix in its last n-1 positions.  e2(1) is arbitrary. */

/*      on output */

/*        d contains the eigenvalues in ascending order.  if an */
/*          error exit is made, the eigenvalues are correct and */
/*          ordered for indices 1,2,...ierr-1, but may not be */
/*          the smallest eigenvalues. */

/*        e2 has been destroyed. */

/*        ierr is set to */
/*          zero       for normal return, */
/*          j          if the j-th eigenvalue has not been */
/*                     determined after 30 iterations. */

/*     calls pythag for  dsqrt(a*a + b*b) . */

/*     questions and comments should be directed to burton s. garbow, */
/*     mathematics and computer science div, argonne national laboratory */

/*     this version dated august 1983. */

/*     ------------------------------------------------------------------ */

    /* Parameter adjustments */
    --e2;
    --d__;

    /* Function Body */
    *ierr = 0;
    if (*n == 1) {
	goto L1001;
    }

    i__1 = *n;
    for (i__ = 2; i__ <= i__1; ++i__) {
/* L100: */
	e2[i__ - 1] = e2[i__];
    }

    f = 0.;
    t = 0.;
    e2[*n] = 0.;

    i__1 = *n;
    for (l = 1; l <= i__1; ++l) {
	j = 0;
	h__ = (d__1 = d__[l], abs(d__1)) + sqrt(e2[l]);
	if (t > h__) {
	    goto L105;
	}
	t = h__;
	b = epslon_(&t);
	c__ = b * b;
/*     .......... look for small squared sub-diagonal element .......... */
L105:
	i__2 = *n;
	for (m = l; m <= i__2; ++m) {
	    if (e2[m] <= c__) {
		goto L120;
	    }
/*     .......... e2(n) is always zero, so there is no exit */
/*                through the bottom of the loop .......... */
/* L110: */
	}

L120:
	if (m == l) {
	    goto L210;
	}
L130:
	if (j == 30) {
	    goto L1000;
	}
	++j;
/*     .......... form shift .......... */
	l1 = l + 1;
	s = sqrt(e2[l]);
	g = d__[l];
	p = (d__[l1] - g) / (s * 2.);
	r__ = pythag_(&p, &c_b43);
	d__[l] = s / (p + d_sign(&r__, &p));
	h__ = g - d__[l];

	i__2 = *n;
	for (i__ = l1; i__ <= i__2; ++i__) {
/* L140: */
	    d__[i__] -= h__;
	}

	f += h__;
/*     .......... rational ql transformation .......... */
	g = d__[m];
	if (g == 0.) {
	    g = b;
	}
	h__ = g;
	s = 0.;
	mml = m - l;
/*     .......... for i=m-1 step -1 until l do -- .......... */
	i__2 = mml;
	for (ii = 1; ii <= i__2; ++ii) {
	    i__ = m - ii;
	    p = g * h__;
	    r__ = p + e2[i__];
	    e2[i__ + 1] = s * r__;
	    s = e2[i__] / r__;
	    d__[i__ + 1] = h__ + s * (h__ + d__[i__]);
	    g = d__[i__] - e2[i__] / g;
	    if (g == 0.) {
		g = b;
	    }
	    h__ = g * p / r__;
/* L200: */
	}

	e2[l] = s * g;
	d__[l] = h__;
/*     .......... guard against underflow in convergence test .......... */
	if (h__ == 0.) {
	    goto L210;
	}
	if ((d__1 = e2[l], abs(d__1)) <= (d__2 = c__ / h__, abs(d__2))) {
	    goto L210;
	}
	e2[l] = h__ * e2[l];
	if (e2[l] != 0.) {
	    goto L130;
	}
L210:
	p = d__[l] + f;
/*     .......... order eigenvalues .......... */
	if (l == 1) {
	    goto L250;
	}
/*     .......... for i=l step -1 until 2 do -- .......... */
	i__2 = l;
	for (ii = 2; ii <= i__2; ++ii) {
	    i__ = l + 2 - ii;
	    if (p >= d__[i__ - 1]) {
		goto L270;
	    }
	    d__[i__] = d__[i__ - 1];
/* L230: */
	}

L250:
	i__ = 1;
L270:
	d__[i__] = p;
/* L290: */
    }

    goto L1001;
/*     .......... set error -- no convergence to an */
/*                eigenvalue after 30 iterations .......... */
L1000:
    *ierr = l;
L1001:
    return 0;
} /* tqlrat_ */

doublereal pythag_(a, b)
doublereal *a, *b;
{
    /* System generated locals */
    doublereal ret_val, d__1, d__2, d__3;

    /* Local variables */
    static doublereal p, r__, s, t, u;


/*     finds dsqrt(a**2+b**2) without overflow or destructive underflow */

/* Computing MAX */
    d__1 = abs(*a), d__2 = abs(*b);
    p = max(d__1,d__2);
    if (p == 0.) {
	goto L20;
    }
/* Computing MIN */
    d__2 = abs(*a), d__3 = abs(*b);
/* Computing 2nd power */
    d__1 = min(d__2,d__3) / p;
    r__ = d__1 * d__1;
L10:
    t = r__ + 4.;
    if (t == 4.) {
	goto L20;
    }
    s = r__ / t;
    u = s * 2. + 1.;
    p = u * p;
/* Computing 2nd power */
    d__1 = s / u;
    r__ = d__1 * d__1 * r__;
    goto L10;
L20:
    ret_val = p;
    return ret_val;
} /* pythag_ */

/* Subroutine */ int tql2_(nm, n, d__, e, z__, ierr)
integer *nm, *n;
doublereal *d__, *e, *z__;
integer *ierr;
{
    /* System generated locals */
    integer z_dim1, z_offset, i__1, i__2, i__3;
    doublereal d__1, d__2;

    /* Builtin functions */
    double d_sign();

    /* Local variables */
    static doublereal c__, f, g, h__;
    static integer i__, j, k, l, m;
    static doublereal p, r__, s, c2, c3;
    static integer l1, l2;
    static doublereal s2;
    static integer ii;
    extern doublereal pythag_();
    static doublereal dl1, el1;
    static integer mml;
    static doublereal tst1, tst2;



/*     this subroutine is a translation of the algol procedure tql2, */
/*     num. math. 11, 293-306(1968) by bowdler, martin, reinsch, and */
/*     wilkinson. */
/*     handbook for auto. comp., vol.ii-linear algebra, 227-240(1971). */

/*     this subroutine finds the eigenvalues and eigenvectors */
/*     of a symmetric tridiagonal matrix by the ql method. */
/*     the eigenvectors of a full symmetric matrix can also */
/*     be found if  tred2  has been used to reduce this */
/*     full matrix to tridiagonal form. */

/*     on input */

/*        nm must be set to the row dimension of two-dimensional */
/*          array parameters as declared in the calling program */
/*          dimension statement. */

/*        n is the order of the matrix. */

/*        d contains the diagonal elements of the input matrix. */

/*        e contains the subdiagonal elements of the input matrix */
/*          in its last n-1 positions.  e(1) is arbitrary. */

/*        z contains the transformation matrix produced in the */
/*          reduction by  tred2, if performed.  if the eigenvectors */
/*          of the tridiagonal matrix are desired, z must contain */
/*          the identity matrix. */

/*      on output */

/*        d contains the eigenvalues in ascending order.  if an */
/*          error exit is made, the eigenvalues are correct but */
/*          unordered for indices 1,2,...,ierr-1. */

/*        e has been destroyed. */

/*        z contains orthonormal eigenvectors of the symmetric */
/*          tridiagonal (or full) matrix.  if an error exit is made, */
/*          z contains the eigenvectors associated with the stored */
/*          eigenvalues. */

/*        ierr is set to */
/*          zero       for normal return, */
/*          j          if the j-th eigenvalue has not been */
/*                     determined after 30 iterations. */

/*     calls pythag for  dsqrt(a*a + b*b) . */

/*     questions and comments should be directed to burton s. garbow, */
/*     mathematics and computer science div, argonne national laboratory */

/*     this version dated august 1983. */

/*     ------------------------------------------------------------------ */

    /* Parameter adjustments */
    z_dim1 = *nm;
    z_offset = 1 + z_dim1 * 1;
    z__ -= z_offset;
    --e;
    --d__;

    /* Function Body */
    *ierr = 0;
    if (*n == 1) {
	goto L1001;
    }

    i__1 = *n;
    for (i__ = 2; i__ <= i__1; ++i__) {
/* L100: */
	e[i__ - 1] = e[i__];
    }

    f = 0.;
    tst1 = 0.;
    e[*n] = 0.;

    i__1 = *n;
    for (l = 1; l <= i__1; ++l) {
	j = 0;
	h__ = (d__1 = d__[l], abs(d__1)) + (d__2 = e[l], abs(d__2));
	if (tst1 < h__) {
	    tst1 = h__;
	}
/*     .......... look for small sub-diagonal element .......... */
	i__2 = *n;
	for (m = l; m <= i__2; ++m) {
	    tst2 = tst1 + (d__1 = e[m], abs(d__1));
	    if (tst2 == tst1) {
		goto L120;
	    }
/*     .......... e(n) is always zero, so there is no exit */
/*                through the bottom of the loop .......... */
/* L110: */
	}

L120:
	if (m == l) {
	    goto L220;
	}
L130:
	if (j == 30) {
	    goto L1000;
	}
	++j;
/*     .......... form shift .......... */
	l1 = l + 1;
	l2 = l1 + 1;
	g = d__[l];
	p = (d__[l1] - g) / (e[l] * 2.);
	r__ = pythag_(&p, &c_b43);
	d__[l] = e[l] / (p + d_sign(&r__, &p));
	d__[l1] = e[l] * (p + d_sign(&r__, &p));
	dl1 = d__[l1];
	h__ = g - d__[l];
	if (l2 > *n) {
	    goto L145;
	}

	i__2 = *n;
	for (i__ = l2; i__ <= i__2; ++i__) {
/* L140: */
	    d__[i__] -= h__;
	}

L145:
	f += h__;
/*     .......... ql transformation .......... */
	p = d__[m];
	c__ = 1.;
	c2 = c__;
	el1 = e[l1];
	s = 0.;
	mml = m - l;
/*     .......... for i=m-1 step -1 until l do -- .......... */
	i__2 = mml;
	for (ii = 1; ii <= i__2; ++ii) {
	    c3 = c2;
	    c2 = c__;
	    s2 = s;
	    i__ = m - ii;
	    g = c__ * e[i__];
	    h__ = c__ * p;
	    r__ = pythag_(&p, &e[i__]);
	    e[i__ + 1] = s * r__;
	    s = e[i__] / r__;
	    c__ = p / r__;
	    p = c__ * d__[i__] - s * g;
	    d__[i__ + 1] = h__ + s * (c__ * g + s * d__[i__]);
/*     .......... form vector .......... */
	    i__3 = *n;
	    for (k = 1; k <= i__3; ++k) {
		h__ = z__[k + (i__ + 1) * z_dim1];
		z__[k + (i__ + 1) * z_dim1] = s * z__[k + i__ * z_dim1] + c__ 
			* h__;
		z__[k + i__ * z_dim1] = c__ * z__[k + i__ * z_dim1] - s * h__;
/* L180: */
	    }

/* L200: */
	}

	p = -s * s2 * c3 * el1 * e[l] / dl1;
	e[l] = s * p;
	d__[l] = c__ * p;
	tst2 = tst1 + (d__1 = e[l], abs(d__1));
	if (tst2 > tst1) {
	    goto L130;
	}
L220:
	d__[l] += f;
/* L240: */
    }
/*     .......... order eigenvalues and eigenvectors .......... */
    i__1 = *n;
    for (ii = 2; ii <= i__1; ++ii) {
	i__ = ii - 1;
	k = i__;
	p = d__[i__];

	i__2 = *n;
	for (j = ii; j <= i__2; ++j) {
	    if (d__[j] >= p) {
		goto L260;
	    }
	    k = j;
	    p = d__[j];
L260:
	    ;
	}

	if (k == i__) {
	    goto L300;
	}
	d__[k] = d__[i__];
	d__[i__] = p;

	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
	    p = z__[j + i__ * z_dim1];
	    z__[j + i__ * z_dim1] = z__[j + k * z_dim1];
	    z__[j + k * z_dim1] = p;
/* L280: */
	}

L300:
	;
    }

    goto L1001;
/*     .......... set error -- no convergence to an */
/*                eigenvalue after 30 iterations .......... */
L1000:
    *ierr = l;
L1001:
    return 0;
} /* tql2_ */

