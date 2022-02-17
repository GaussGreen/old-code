/*******************************************************************
 *      cumchn.c:
 * Routines translated by f2c from
 *           DSTATTAB
 * Double Precision Statistical Tables
 * For cumchn(2197  ,10  ,2617) other implementations
 * return 0. This one returns  .603e-5.
 *(see http://odin.mdacc.tmc.edu/anonftp/page_2.html)
 * Using f2c translator for translating fortran to C.
 * (see http://www.netlib.org/)
 * Author (?) :Cyril Godart
 * The whole library has been translated and the makefile as well
 * as the fortran library have been put under source safe on
 * Ldgs000007\SourceSafe\Sourcesafe\first\stat
 * Database : FirstSwaps\dstattab
 * Date: 14/12/98
 ********************************************************************/
#include "UTALLHDR.H>
#include "math.h"
#include "num_h_f2c.h"

double d_sign(doublereal *a, doublereal *b) {
  double x;
  x = (*a >= 0 ? *a : -*a);
  return (*b >= 0 ? x : -x);
}

double pow_di(doublereal *ap, integer *bp) {
  double pow, x;
  integer n;
  unsigned long u;

  pow = 1;
  x = *ap;
  n = *bp;

  if (n != 0) {
    if (n < 0) {
      n = -n;
      x = 1 / x;
    }
    for (u = n;;) {
      if (u & 01)
        pow *= x;
      if (u >>= 1)
        x *= x;
      else
        break;
    }
  }
  return (pow);
}

/* IPMPAR.F -- translated by f2c (version 19980913).
   You must link the resulting object file with the libraries:
        -lf2c -lm   (in that order)
*/

integer ipmpar_(integer *i__) {
  /* Initialized data */

  static integer imach[10] = {2,    31,  2147483647, 2,     24,
                              -125, 128, 53,         -1021, 1024};

  /* System generated locals */
  integer ret_val;

  /* ----------------------------------------------------------------------- */

  /*     IPMPAR PROVIDES THE INTEGER MACHINE CONSTANTS FOR THE COMPUTER */
  /*     THAT IS USED. IT IS ASSUMED THAT THE ARGUMENT I IS AN INTEGER */
  /*     HAVING ONE OF THE VALUES 1-10. IPMPAR(I) HAS THE VALUE ... */

  /*  INTEGERS. */

  /*     ASSUME INTEGERS ARE REPRESENTED IN THE N-DIGIT  , BASE-A FORM */

  /*               SIGN ( X(N-1)*A**(N-1) + ... + X(1)*A + X(0) ) */

  /*               WHERE 0 .LE. X(I) .LT. A FOR I=0  ,...  ,N-1. */

  /*     IPMPAR(1) = A  , THE BASE. */

  /*     IPMPAR(2) = N  , THE NUMBER OF BASE-A DIGITS. */

  /*     IPMPAR(3) = A**N - 1  , THE LARGEST MAGNITUDE. */

  /*  FLOATING-POINT NUMBERS. */

  /*     IT IS ASSUMED THAT THE SINGLE AND DOUBLE PRECISION FLOATING */
  /*     POINT ARITHMETICS HAVE THE SAME BASE  , SAY B  , AND THAT THE */
  /*     NONZERO NUMBERS ARE REPRESENTED IN THE FORM */

  /*               SIGN (B**E) * (X(1)/B + ... + X(M)/B**M) */

  /*               WHERE X(I) = 0  ,1  ,...  ,B-1 FOR I=1  ,...  ,M  , */
  /*               X(1) .GE. 1  , AND EMIN .LE. E .LE. EMAX. */

  /*     IPMPAR(4) = B  , THE BASE. */

  /*  SINGLE-PRECISION */

  /*     IPMPAR(5) = M  , THE NUMBER OF BASE-B DIGITS. */

  /*     IPMPAR(6) = EMIN  , THE SMALLEST EXPONENT E. */

  /*     IPMPAR(7) = EMAX  , THE LARGEST EXPONENT E. */

  /*  DOUBLE-PRECISION */

  /*     IPMPAR(8) = M  , THE NUMBER OF BASE-B DIGITS. */

  /*     IPMPAR(9) = EMIN  , THE SMALLEST EXPONENT E. */

  /*     IPMPAR(10) = EMAX  , THE LARGEST EXPONENT E. */

  /* ----------------------------------------------------------------------- */

  /*     TO DEFINE THIS FUNCTION FOR THE COMPUTER BEING USED  , ACTIVATE */
  /*     THE DATA STATMENTS FOR THE COMPUTER BY REMOVING THE C FROM */
  /*     COLUMN 1. (ALL THE OTHER DATA STATEMENTS SHOULD HAVE C IN */
  /*     COLUMN 1.) */

  /* ----------------------------------------------------------------------- */

  /*     IPMPAR IS AN ADAPTATION OF THE FUNCTION I1MACH  , WRITTEN BY */
  /*     P.A. FOX  , A.D. HALL  , AND N.L. SCHRYER (BELL LABORATORIES). */
  /*     IPMPAR WAS FORMED BY A.H. MORRIS (NSWC). THE CONSTANTS ARE */
  /*     FROM BELL LABORATORIES  , NSWC  , AND OTHER SOURCES. */

  /* ----------------------------------------------------------------------- */
  /*     .. Scalar Arguments .. */
  /*     .. */
  /*     .. Local Arrays .. */
  /*     .. */
  /*     .. Data statements .. */

  /*     DATA IMACH(10) /  127 / */

  ret_val = imach[*i__ - 1];
  return ret_val;
} /* ipmpar_ */

/* EXPARG.F -- translated by f2c (version 19980913).
   You must link the resulting object file with the libraries:
        -lf2c -lm   (in that order)
*/

/* Table of constant values */

static integer c__4 = 4;
static integer c__8 = 8;
static integer c__9 = 9;
static integer c__10 = 10;

doublereal exparg_(integer *l) {
  /* System generated locals */
  doublereal ret_val;

  /* Builtin functions */
  double log(doublereal);

  /* Local variables */
  static integer b, m;
  extern integer ipmpar_(integer *);
  static doublereal lnb;

  /* -------------------------------------------------------------------- */
  /*     IF L = 0 THEN  EXPARG(L) = THE LARGEST POSITIVE W FOR WHICH */
  /*     EXP(W) CAN BE COMPUTED. */

  /*     IF L IS NONZERO THEN  EXPARG(L) = THE LARGEST NEGATIVE W FOR */
  /*     WHICH THE COMPUTED VALUE OF EXP(W) IS NONZERO. */

  /*     NOTE... ONLY AN APPROXIMATE VALUE FOR EXPARG(L) IS NEEDED. */
  /* -------------------------------------------------------------------- */
  /*     .. Scalar Arguments .. */
  /*     .. */
  /*     .. Local Scalars .. */
  /*     .. */
  /*     .. External Functions .. */
  /*     .. */
  /*     .. Intrinsic Functions .. */
  /*     .. */
  /*     .. Executable Statements .. */

  b = ipmpar_(&c__4);
  if (b != 2) {
    goto L10;
  }
  lnb = .69314718055995;
  goto L40;
L10:
  if (b != 8) {
    goto L20;
  }
  lnb = 2.0794415416798;
  goto L40;
L20:
  if (b != 16) {
    goto L30;
  }
  lnb = 2.7725887222398;
  goto L40;
L30:
  lnb = log((doublereal)b);

L40:
  if (*l == 0) {
    goto L50;
  }
  m = ipmpar_(&c__9) - 1;
  ret_val = m * lnb * .99999;
  return ret_val;
L50:
  m = ipmpar_(&c__10);
  ret_val = m * lnb * .99999;
  return ret_val;
} /* exparg_ */

/* DEVLPL.F -- translated by f2c (version 19980913).
   You must link the resulting object file with the libraries:
        -lf2c -lm   (in that order)
*/

doublereal devlpl_(doublereal *a, integer *n, doublereal *x) {
  /* System generated locals */
  doublereal ret_val;

  /* Local variables */
  static doublereal term;
  static integer i__;

  /* ********************************************************************** */

  /*     DOUBLE PRECISION FUNCTION DEVLPL(A  ,N  ,X) */
  /*              Double precision EVALuate a PoLynomial at X */

  /*                              Function */

  /*     returns */
  /*          A(1) + A(2)*X + ... + A(N)*X**(N-1) */

  /*                              Arguments */

  /*     A --> Array of coefficients of the polynomial. */
  /*                                        A is DOUBLE PRECISION(N) */

  /*     N --> Length of A  , also degree of polynomial - 1. */
  /*                                        N is INTEGER */

  /*     X --> Point at which the polynomial is to be evaluated. */
  /*                                        X is DOUBLE PRECISION */

  /* ********************************************************************** */

  /*     .. Scalar Arguments .. */
  /*     .. */
  /*     .. Array Arguments .. */
  /*     .. */
  /*     .. Local Scalars .. */
  /*     .. */
  /*     .. Executable Statements .. */
  /* Parameter adjustments */
  --a;

  /* Function Body */
  term = a[*n];
  for (i__ = *n - 1; i__ >= 1; --i__) {
    term = a[i__] + term * *x;
    /* L10: */
  }
  ret_val = term;
  return ret_val;
} /* devlpl_ */

/* ALNGAM.F -- translated by f2c (version 19980913).
   You must link the resulting object file with the libraries:
        -lf2c -lm   (in that order)
*/

/* Table of constant values */

static integer c__5 = 5;

doublereal alngam_(doublereal *x) {
  /* Initialized data */

  static doublereal scoefn[9] = {
      62.003838007127258804, 36.036772530024836321, 20.782472531792126786,
      6.338067999387272343,  2.15994312846059073,   .3980671310203570498,
      .1093115956710439502,  .0092381945590275995,  .0029737866448101651};
  static doublereal scoefd[4] = {62.003838007126989331, 9.822521104713994894,
                                 -8.906016659497461257, 1.};
  static doublereal coef[5] = {.083333333333333023564, -.0027777777768818808,
                               7.9365006754279e-4, -5.94997310889e-4,
                               8.065880899e-4};

  /* System generated locals */
  integer i__1;
  doublereal ret_val, d__1, d__2;

  /* Builtin functions */
  double log(doublereal);

  /* Local variables */
  static doublereal prod;
  static integer i__, n;
  static doublereal xx, offset;
  extern doublereal devlpl_(doublereal *, integer *, doublereal *);

  /* ********************************************************************** */

  /*     DOUBLE PRECISION FUNCTION ALNGAM(X) */
  /*                 double precision LN of the GAMma function */

  /*                              Function */

  /*     Returns the natural logarithm of GAMMA(X). */

  /*                              Arguments */

  /*     X --> value at which scaled log gamma is to be returned */
  /*                    X is DOUBLE PRECISION */

  /*                              Method */

  /*     If X .le. 6.0  , then use recursion to get X below 3 */
  /*     then apply rational approximation number 5236 of */
  /*     Hart et al  , Computer Approximations  , John Wiley and */
  /*     Sons  , NY  , 1968. */

  /*     If X .gt. 6.0  , then use recursion to get X to at least 12 and */
  /*     then use formula 5423 of the same source. */

  /* ********************************************************************** */

  /*     .. Parameters .. */
  /*     .. */
  /*     .. Scalar Arguments .. */
  /*     .. */
  /*     .. Local Scalars .. */
  /*     .. */
  /*     .. Local Arrays .. */
  /*     .. */
  /*     .. External Functions .. */
  /*     .. */
  /*     .. Intrinsic Functions .. */
  /*     .. */
  /*     .. Data statements .. */
  /*     .. */
  /*     .. Executable Statements .. */
  if (!(*x <= 6.)) {
    goto L70;
  }
  prod = 1.;
  xx = *x;
  if (!(*x > 3.)) {
    goto L30;
  }
L10:
  if (!(xx > 3.)) {
    goto L20;
  }
  xx += -1.;
  prod *= xx;
  goto L10;
L20:
L30:
  if (!(*x < 2.)) {
    goto L60;
  }
L40:
  if (!(xx < 2.)) {
    goto L50;
  }
  prod /= xx;
  xx += 1.;
  goto L40;
L50:
L60:
  d__1 = xx - 2.;
  d__2 = xx - 2.;
  ret_val = devlpl_(scoefn, &c__9, &d__1) / devlpl_(scoefd, &c__4, &d__2);

  /*     COMPUTE RATIONAL APPROXIMATION TO GAMMA(X) */

  ret_val *= prod;
  ret_val = log(ret_val);
  goto L110;
L70:
  offset = .91893853320467274178;

  /*     IF NECESSARY MAKE X AT LEAST 12 AND CARRY CORRECTION IN OFFSET */

  n = (integer)(12. - *x);
  if (!(n > 0)) {
    goto L90;
  }
  prod = 1.;
  i__1 = n;
  for (i__ = 1; i__ <= i__1; ++i__) {
    prod *= *x + (doublereal)(i__ - 1);
    /* L80: */
  }
  offset -= log(prod);
  xx = *x + (doublereal)n;
  goto L100;
L90:
  xx = *x;

  /*     COMPUTE POWER SERIES */

L100:
  /* Computing 2nd power */
  d__2 = xx;
  d__1 = 1. / (d__2 * d__2);
  ret_val = devlpl_(coef, &c__5, &d__1) / xx;
  ret_val = ret_val + offset + (xx - .5) * log(xx) - xx;
L110:
  return ret_val;
} /* alngam_ */

/* GAM1.F -- translated by f2c (version 19980913).
   You must link the resulting object file with the libraries:
        -lf2c -lm   (in that order)
*/

doublereal gam1_(doublereal *a) {
  /* Initialized data */

  static doublereal p[7] = {.577215664901533,   -.409078193005776,
                            -.230975380857675,  .0597275330452234,
                            .0076696818164949,  -.00514889771323592,
                            5.89597428611429e-4};
  static doublereal q[5] = {1., .427569613095214, .158451672430138,
                            .0261132021441447, .00423244297896961};
  static doublereal r__[9] = {
      -.422784335098468,  -.771330383816272,   -.244757765222226,
      .118378989872749,   9.30357293360349e-4, -.0118290993445146,
      .00223047661158249, 2.66505979058923e-4, -1.32674909766242e-4};
  static doublereal s1 = .273076135303957;
  static doublereal s2 = .0559398236957378;

  /* System generated locals */
  doublereal ret_val;

  /* Local variables */
  static doublereal d__, t, w, bot, top;

  /*     ------------------------------------------------------------------ */
  /*     COMPUTATION OF 1/GAMMA(A+1) - 1  FOR -0.5 .LE. A .LE. 1.5 */
  /*     ------------------------------------------------------------------ */
  /*     .. Scalar Arguments .. */
  /*     .. */
  /*     .. Local Scalars .. */
  /*     .. */
  /*     .. Local Arrays .. */
  /*     .. */
  /*     .. Data statements .. */
  /*     ------------------- */
  /*     ------------------- */
  /*     ------------------- */
  /*     ------------------- */
  /*     .. */
  /*     .. Executable Statements .. */
  /*     ------------------- */
  t = *a;
  d__ = *a - .5;
  if (d__ > 0.) {
    t = d__ - .5;
  }
  if (t < 0.) {
    goto L40;
  } else if (t == 0) {
    goto L10;
  } else {
    goto L20;
  }

L10:
  ret_val = 0.;
  return ret_val;

L20:
  top = (((((p[6] * t + p[5]) * t + p[4]) * t + p[3]) * t + p[2]) * t + p[1]) *
            t +
        p[0];
  bot = (((q[4] * t + q[3]) * t + q[2]) * t + q[1]) * t + 1.;
  w = top / bot;
  if (d__ > 0.) {
    goto L30;
  }
  ret_val = *a * w;
  return ret_val;
L30:
  ret_val = t / *a * (w - .5 - .5);
  return ret_val;

L40:
  top = (((((((r__[8] * t + r__[7]) * t + r__[6]) * t + r__[5]) * t + r__[4]) *
               t +
           r__[3]) *
              t +
          r__[2]) *
             t +
         r__[1]) *
            t +
        r__[0];
  bot = (s2 * t + s1) * t + 1.;
  w = top / bot;
  if (d__ > 0.) {
    goto L50;
  }
  ret_val = *a * (w + .5 + .5);
  return ret_val;
L50:
  ret_val = t * w / *a;
  return ret_val;
} /* gam1_ */

/* ERF.F -- translated by f2c (version 19980913).
   You must link the resulting object file with the libraries:
        -lf2c -lm   (in that order)
*/

/* Table of constant values */

static doublereal c_b5 = 1.;

doublereal erf_(doublereal *x) {
  /* Initialized data */

  static doublereal c__ = .564189583547756;
  static doublereal a[5] = {7.7105849500132e-5, -.00133733772997339,
                            .0323076579225834, .0479137145607681,
                            .128379167095513};
  static doublereal b[3] = {.00301048631703895, .0538971687740286,
                            .375795757275549};
  static doublereal p[8] = {-1.36864857382717e-7, .564195517478974,
                            7.21175825088309,     43.1622272220567,
                            152.98928504694,      339.320816734344,
                            451.918953711873,     300.459261020162};
  static doublereal q[8] = {1.,
                            12.7827273196294,
                            77.0001529352295,
                            277.585444743988,
                            638.980264465631,
                            931.35409485061,
                            790.950925327898,
                            300.459260956983};
  static doublereal r__[5] = {2.10144126479064, 26.2370141675169,
                              21.3688200555087, 4.6580782871847,
                              .282094791773523};
  static doublereal s[4] = {94.153775055546, 187.11481179959, 99.0191814623914,
                            18.0124575948747};

  /* System generated locals */
  doublereal ret_val;

  /* Builtin functions */
  double exp(doublereal), d_sign(doublereal *, doublereal *);

  /* Local variables */
  static doublereal t, x2, ax, bot, top;

  /* ----------------------------------------------------------------------- */
  /*             EVALUATION OF THE REAL ERROR FUNCTION */
  /* ----------------------------------------------------------------------- */
  /*     .. Scalar Arguments .. */
  /*     .. */
  /*     .. Local Scalars .. */
  /*     .. */
  /*     .. Local Arrays .. */
  /*     .. */
  /*     .. Intrinsic Functions .. */
  /*     .. */
  /*     .. Data statements .. */
  /* ------------------------- */
  /* ------------------------- */
  /* ------------------------- */
  /* ------------------------- */
  /*     .. */
  /*     .. Executable Statements .. */
  /* ------------------------- */
  ax = fabs(*x);
  if (ax > .5) {
    goto L10;
  }
  t = *x * *x;
  top = (((a[0] * t + a[1]) * t + a[2]) * t + a[3]) * t + a[4] + 1.;
  bot = ((b[0] * t + b[1]) * t + b[2]) * t + 1.;
  ret_val = *x * (top / bot);
  return ret_val;

L10:
  if (ax > 4.) {
    goto L20;
  }
  top = ((((((p[0] * ax + p[1]) * ax + p[2]) * ax + p[3]) * ax + p[4]) * ax +
          p[5]) *
             ax +
         p[6]) *
            ax +
        p[7];
  bot = ((((((q[0] * ax + q[1]) * ax + q[2]) * ax + q[3]) * ax + q[4]) * ax +
          q[5]) *
             ax +
         q[6]) *
            ax +
        q[7];
  ret_val = .5 - exp(-(*x) * *x) * top / bot + .5;
  if (*x < 0.) {
    ret_val = -ret_val;
  }
  return ret_val;

L20:
  if (ax >= 5.8) {
    goto L30;
  }
  x2 = *x * *x;
  t = 1. / x2;
  top = (((r__[0] * t + r__[1]) * t + r__[2]) * t + r__[3]) * t + r__[4];
  bot = (((s[0] * t + s[1]) * t + s[2]) * t + s[3]) * t + 1.;
  ret_val = (c__ - top / (x2 * bot)) / ax;
  ret_val = .5 - exp(-x2) * ret_val + .5;
  if (*x < 0.) {
    ret_val = -ret_val;
  }
  return ret_val;

L30:
  ret_val = d_sign(&c_b5, x);
  return ret_val;
} /* erf_ */

/* SPMPAR.F -- translated by f2c (version 19980913).
   You must link the resulting object file with the libraries:
        -lf2c -lm   (in that order)
*/

/* Table of constant values */

doublereal spmpar_(integer *i__) {
  /* System generated locals */
  integer i__1;
  doublereal ret_val;

  /* Builtin functions */
  double pow_di(doublereal *, integer *);

  /* Local variables */
  static integer emin, emax;
  static doublereal binv, b;
  static integer m, ibeta;
  static doublereal w, z__;
  extern integer ipmpar_(integer *);
  static doublereal bm1, one;

  /* ----------------------------------------------------------------------- */

  /*     SPMPAR PROVIDES THE SINGLE PRECISION MACHINE CONSTANTS FOR */
  /*     THE COMPUTER BEING USED. IT IS ASSUMED THAT THE ARGUMENT */
  /*     I IS AN INTEGER HAVING ONE OF THE VALUES 1  , 2  , OR 3. IF THE */
  /*     SINGLE PRECISION ARITHMETIC BEING USED HAS M BASE B DIGITS AND */
  /*     ITS SMALLEST AND LARGEST EXPONENTS ARE EMIN AND EMAX  , THEN */

  /*        SPMPAR(1) = B**(1 - M)  , THE MACHINE PRECISION  , */

  /*        SPMPAR(2) = B**(EMIN - 1)  , THE SMALLEST MAGNITUDE  , */

  /*        SPMPAR(3) = B**EMAX*(1 - B**(-M))  , THE LARGEST MAGNITUDE. */

  /* ----------------------------------------------------------------------- */
  /*     WRITTEN BY */
  /*        ALFRED H. MORRIS  , JR. */
  /*        NAVAL SURFACE WARFARE CENTER */
  /*        DAHLGREN VIRGINIA */
  /* ----------------------------------------------------------------------- */
  /* ----------------------------------------------------------------------- */
  /*     MODIFIED BY BARRY W. BROWN TO RETURN DOUBLE PRECISION MACHINE */
  /*     CONSTANTS FOR THE COMPUTER BEING USED.  THIS MODIFICATION WAS */
  /*     MADE AS PART OF CONVERTING BRATIO TO DOUBLE PRECISION */
  /* ----------------------------------------------------------------------- */
  /*     .. Scalar Arguments .. */
  /*     .. */
  /*     .. Local Scalars .. */
  /*     .. */
  /*     .. External Functions .. */
  /*     .. */
  /*     .. Intrinsic Functions .. */
  /*     .. */
  /*     .. Executable Statements .. */

  if (*i__ > 1) {
    goto L10;
  }
  b = (doublereal)ipmpar_(&c__4);
  m = ipmpar_(&c__8);
  i__1 = 1 - m;
  ret_val = pow_di(&b, &i__1);
  return ret_val;

L10:
  if (*i__ > 2) {
    goto L20;
  }
  b = (doublereal)ipmpar_(&c__4);
  emin = ipmpar_(&c__9);
  one = 1.;
  binv = one / b;
  i__1 = emin + 2;
  w = pow_di(&b, &i__1);
  ret_val = w * binv * binv * binv;
  return ret_val;

L20:
  ibeta = ipmpar_(&c__4);
  m = ipmpar_(&c__8);
  emax = ipmpar_(&c__10);

  b = (doublereal)ibeta;
  bm1 = (doublereal)(ibeta - 1);
  one = 1.;
  i__1 = m - 1;
  z__ = pow_di(&b, &i__1);
  w = ((z__ - one) * b + bm1) / (b * z__);

  i__1 = emax - 2;
  z__ = pow_di(&b, &i__1);
  ret_val = w * z__ * b * b;
  return ret_val;
} /* spmpar_ */

/* GAMMA.F -- translated by f2c (version 19980913).
   You must link the resulting object file with the libraries:
        -lf2c -lm   (in that order)
*/
/* Table of constant values */

static integer c__3 = 3;
static integer c__0 = 0;

doublereal gamma_(doublereal *a) {
  /* Initialized data */

  static doublereal pi = 3.1415926535898;
  static doublereal d__ = .41893853320467274178;
  static doublereal p[7] = {5.39637273585445e-4,
                            .0026193926004269,
                            .020449366759492,
                            .0730981088720487,
                            .279648642639792,
                            .553413866010467,
                            1.};
  static doublereal q[7] = {-8.32979206704073e-4,
                            .00470059485860584,
                            .022521113103534,
                            -.17045896931336,
                            -.056790276197494,
                            1.13062953091122,
                            1.};
  static doublereal r1 = 8.20756370353826e-4;
  static doublereal r2 = -5.95156336428591e-4;
  static doublereal r3 = 7.93650663183693e-4;
  static doublereal r4 = -.00277777777770481;
  static doublereal r5 = .0833333333333333;

  /* System generated locals */
  integer i__1;
  doublereal ret_val;

  /* Builtin functions */
  double sin(doublereal), log(doublereal), exp(doublereal);

  /* Local variables */
  static doublereal g;
  static integer i__, j, m, n;
  static doublereal s, t, w, x, z__;
  extern doublereal exparg_(integer *), spmpar_(integer *);
  static doublereal bot, lnx, top;

  /* ----------------------------------------------------------------------- */

  /*         EVALUATION OF THE GAMMA FUNCTION FOR REAL ARGUMENTS */

  /*                           ----------- */

  /*     GAMMA(A) IS ASSIGNED THE VALUE 0 WHEN THE GAMMA FUNCTION CANNOT */
  /*     BE COMPUTED. */

  /* ----------------------------------------------------------------------- */
  /*     WRITTEN BY ALFRED H. MORRIS  , JR. */
  /*          NAVAL SURFACE WEAPONS CENTER */
  /*          DAHLGREN  , VIRGINIA */
  /* ----------------------------------------------------------------------- */
  /*     .. Scalar Arguments .. */
  /*     .. */
  /*     .. Local Scalars .. */
  /*     .. */
  /*     .. Local Arrays .. */
  /*     .. */
  /*     .. External Functions .. */
  /*     .. */
  /*     .. Intrinsic Functions .. */
  /*     .. */
  /*     .. Data statements .. */
  /* -------------------------- */
  /*     D = 0.5*(LN(2*PI) - 1) */
  /* -------------------------- */
  /* -------------------------- */
  /* -------------------------- */
  /*     .. */
  /*     .. Executable Statements .. */
  /* -------------------------- */
  ret_val = 0.;
  x = *a;
  if (fabs(*a) >= 15.) {
    goto L110;
  }
  /* ----------------------------------------------------------------------- */
  /*            EVALUATION OF GAMMA(A) FOR ABS(A) .LT. 15 */
  /* ----------------------------------------------------------------------- */
  t = 1.;
  m = (integer)(*a) - 1;

  /*     LET T BE THE PRODUCT OF A-J WHEN A .GE. 2 */

  if (m < 0) {
    goto L40;
  } else if (m == 0) {
    goto L30;
  } else {
    goto L10;
  }
L10:
  i__1 = m;
  for (j = 1; j <= i__1; ++j) {
    x += -1.;
    t = x * t;
    /* L20: */
  }
L30:
  x += -1.;
  goto L80;

  /*     LET T BE THE PRODUCT OF A+J WHEN A .LT. 1 */

L40:
  t = *a;
  if (*a > 0.) {
    goto L70;
  }
  m = -m - 1;
  if (m == 0) {
    goto L60;
  }
  i__1 = m;
  for (j = 1; j <= i__1; ++j) {
    x += 1.;
    t = x * t;
    /* L50: */
  }
L60:
  x = x + .5 + .5;
  t = x * t;
  if (t == 0.) {
    return ret_val;
  }

L70:

  /*     THE FOLLOWING CODE CHECKS IF 1/T CAN OVERFLOW. THIS */
  /*     CODE MAY BE OMITTED IF DESIRED. */

  if (fabs(t) >= 1e-30) {
    goto L80;
  }
  if (fabs(t) * spmpar_(&c__3) <= 1.0001) {
    return ret_val;
  }
  ret_val = 1. / t;
  return ret_val;

  /*     COMPUTE GAMMA(1 + X) FOR  0 .LE. X .LT. 1 */

L80:
  top = p[0];
  bot = q[0];
  for (i__ = 2; i__ <= 7; ++i__) {
    top = p[i__ - 1] + x * top;
    bot = q[i__ - 1] + x * bot;
    /* L90: */
  }
  ret_val = top / bot;

  /*     TERMINATION */

  if (*a < 1.) {
    goto L100;
  }
  ret_val *= t;
  return ret_val;
L100:
  ret_val /= t;
  return ret_val;
/* ----------------------------------------------------------------------- */
/*            EVALUATION OF GAMMA(A) FOR ABS(A) .GE. 15 */
/* ----------------------------------------------------------------------- */
L110:
  if (fabs(*a) >= 1e3) {
    return ret_val;
  }
  if (*a > 0.) {
    goto L120;
  }
  x = -(*a);
  n = (integer)x;
  t = x - n;
  if (t > .9) {
    t = 1. - t;
  }
  s = sin(pi * t) / pi;
  if (n % 2 == 0) {
    s = -s;
  }
  if (s == 0.) {
    return ret_val;
  }

  /*     COMPUTE THE MODIFIED ASYMPTOTIC SUM */

L120:
  t = 1. / (x * x);
  g = ((((r1 * t + r2) * t + r3) * t + r4) * t + r5) / x;

  /*     ONE MAY REPLACE THE NEXT STATEMENT WITH  LNX = ALOG(X) */
  /*     BUT LESS ACCURACY WILL NORMALLY BE OBTAINED. */

  lnx = log(x);

  /*     FINAL ASSEMBLY */

  z__ = x;
  g = d__ + g + (z__ - .5) * (lnx - 1.);
  w = g;
  t = g - w;
  if (w > exparg_(&c__0) * .99999) {
    return ret_val;
  }
  ret_val = exp(w) * (t + 1.);
  if (*a < 0.) {
    ret_val = 1. / (ret_val * s) / x;
  }
  return ret_val;
} /* gamma_ */

/* REXP.F -- translated by f2c (version 19980913).
   You must link the resulting object file with the libraries:
        -lf2c -lm   (in that order)
*/

doublereal rexp_(doublereal *x) {
  /* Initialized data */

  static doublereal p1 = 9.14041914819518e-10;
  static doublereal p2 = .0238082361044469;
  static doublereal q1 = -.499999999085958;
  static doublereal q2 = .107141568980644;
  static doublereal q3 = -.0119041179760821;
  static doublereal q4 = 5.95130811860248e-4;

  /* System generated locals */
  doublereal ret_val;

  /* Builtin functions */
  double exp(doublereal);

  /* Local variables */
  static doublereal w;

  /* ----------------------------------------------------------------------- */
  /*            EVALUATION OF THE FUNCTION EXP(X) - 1 */
  /* ----------------------------------------------------------------------- */
  /*     .. Scalar Arguments .. */
  /*     .. */
  /*     .. Local Scalars .. */
  /*     .. */
  /*     .. Intrinsic Functions .. */
  /*     .. */
  /*     .. Data statements .. */
  /*     .. */
  /*     .. Executable Statements .. */
  /* ----------------------- */
  if (fabs(*x) > .15) {
    goto L10;
  }
  ret_val = *x * (((p2 * *x + p1) * *x + 1.) /
                  ((((q4 * *x + q3) * *x + q2) * *x + q1) * *x + 1.));
  return ret_val;

L10:
  w = exp(*x);
  if (*x > 0.) {
    goto L20;
  }
  ret_val = w - .5 - .5;
  return ret_val;
L20:
  ret_val = w * (.5 - 1. / w + .5);
  return ret_val;
} /* rexp_ */

/* RLOG.F -- translated by f2c (version 19980913).
   You must link the resulting object file with the libraries:
        -lf2c -lm   (in that order)
*/

doublereal rlog_(doublereal *x) {
  /* Initialized data */

  static doublereal a = .0566749439387324;
  static doublereal b = .0456512608815524;
  static doublereal p0 = .333333333333333;
  static doublereal p1 = -.224696413112536;
  static doublereal p2 = .00620886815375787;
  static doublereal q1 = -1.27408923933623;
  static doublereal q2 = .354508718369557;

  /* System generated locals */
  doublereal ret_val;

  /* Builtin functions */
  double log(doublereal);

  /* Local variables */
  static doublereal r__, t, u, w, w1;

  /*     ------------------- */
  /*     COMPUTATION OF  X - 1 - LN(X) */
  /*     ------------------- */
  /*     .. Scalar Arguments .. */
  /*     .. */
  /*     .. Local Scalars .. */
  /*     .. */
  /*     .. Intrinsic Functions .. */
  /*     .. */
  /*     .. Data statements .. */
  /*     ------------------- */
  /*     .. */
  /*     .. Executable Statements .. */
  /*     ------------------- */
  if (*x < .61 || *x > 1.57) {
    goto L40;
  }
  if (*x < .82) {
    goto L10;
  }
  if (*x > 1.18) {
    goto L20;
  }

  /*              ARGUMENT REDUCTION */

  u = *x - .5 - .5;
  w1 = 0.;
  goto L30;

L10:
  u = *x - .7;
  u /= .7;
  w1 = a - u * .3;
  goto L30;

L20:
  u = *x * .75 - 1.;
  w1 = b + u / 3.;

  /*               SERIES EXPANSION */

L30:
  r__ = u / (u + 2.);
  t = r__ * r__;
  w = ((p2 * t + p1) * t + p0) / ((q2 * t + q1) * t + 1.);
  ret_val = t * 2. * (1. / (1. - r__) - r__ * w) + w1;
  return ret_val;

L40:
  r__ = *x - .5 - .5;
  ret_val = r__ - log(*x);
  return ret_val;
} /* rlog_ */

/* ERFC1.F -- translated by f2c (version 19980913).
   You must link the resulting object file with the libraries:
        -lf2c -lm   (in that order)
*/

/* Table of constant values */

static integer c__1 = 1;

doublereal erfc1_(integer *ind, doublereal *x) {
  /* Initialized data */

  static doublereal c__ = .564189583547756;
  static doublereal a[5] = {7.7105849500132e-5, -.00133733772997339,
                            .0323076579225834, .0479137145607681,
                            .128379167095513};
  static doublereal b[3] = {.00301048631703895, .0538971687740286,
                            .375795757275549};
  static doublereal p[8] = {-1.36864857382717e-7, .564195517478974,
                            7.21175825088309,     43.1622272220567,
                            152.98928504694,      339.320816734344,
                            451.918953711873,     300.459261020162};
  static doublereal q[8] = {1.,
                            12.7827273196294,
                            77.0001529352295,
                            277.585444743988,
                            638.980264465631,
                            931.35409485061,
                            790.950925327898,
                            300.459260956983};
  static doublereal r__[5] = {2.10144126479064, 26.2370141675169,
                              21.3688200555087, 4.6580782871847,
                              .282094791773523};
  static doublereal s[4] = {94.153775055546, 187.11481179959, 99.0191814623914,
                            18.0124575948747};

  /* System generated locals */
  doublereal ret_val, d__1;

  /* Builtin functions */
  double exp(doublereal);

  /* Local variables */
  static doublereal e, t, w, ax;
  extern doublereal exparg_(integer *);
  static doublereal bot, top;

  /* ----------------------------------------------------------------------- */
  /*         EVALUATION OF THE COMPLEMENTARY ERROR FUNCTION */

  /*          ERFC1(IND  ,X) = ERFC(X)            IF IND = 0 */
  /*          ERFC1(IND  ,X) = EXP(X*X)*ERFC(X)   OTHERWISE */
  /* ----------------------------------------------------------------------- */
  /*     .. Scalar Arguments .. */
  /*     .. */
  /*     .. Local Scalars .. */
  /*     .. */
  /*     .. Local Arrays .. */
  /*     .. */
  /*     .. External Functions .. */
  /*     .. */
  /*     .. Intrinsic Functions .. */
  /*     .. */
  /*     .. Data statements .. */
  /* ------------------------- */
  /* ------------------------- */
  /* ------------------------- */
  /* ------------------------- */
  /*     .. */
  /*     .. Executable Statements .. */
  /* ------------------------- */

  /*                     ABS(X) .LE. 0.5 */

  ax = fabs(*x);
  if (ax > .5) {
    goto L10;
  }
  t = *x * *x;
  top = (((a[0] * t + a[1]) * t + a[2]) * t + a[3]) * t + a[4] + 1.;
  bot = ((b[0] * t + b[1]) * t + b[2]) * t + 1.;
  ret_val = .5 - *x * (top / bot) + .5;
  if (*ind != 0) {
    ret_val = exp(t) * ret_val;
  }
  return ret_val;

  /*                  0.5 .LT. ABS(X) .LE. 4 */

L10:
  if (ax > 4.) {
    goto L20;
  }
  top = ((((((p[0] * ax + p[1]) * ax + p[2]) * ax + p[3]) * ax + p[4]) * ax +
          p[5]) *
             ax +
         p[6]) *
            ax +
        p[7];
  bot = ((((((q[0] * ax + q[1]) * ax + q[2]) * ax + q[3]) * ax + q[4]) * ax +
          q[5]) *
             ax +
         q[6]) *
            ax +
        q[7];
  ret_val = top / bot;
  goto L40;

  /*                      ABS(X) .GT. 4 */

L20:
  if (*x <= -5.6) {
    goto L60;
  }
  if (*ind != 0) {
    goto L30;
  }
  if (*x > 100.) {
    goto L70;
  }
  if (*x * *x > -exparg_(&c__1)) {
    goto L70;
  }

L30:
  /* Computing 2nd power */
  d__1 = 1. / *x;
  t = d__1 * d__1;
  top = (((r__[0] * t + r__[1]) * t + r__[2]) * t + r__[3]) * t + r__[4];
  bot = (((s[0] * t + s[1]) * t + s[2]) * t + s[3]) * t + 1.;
  ret_val = (c__ - t * top / bot) / ax;

  /*                      FINAL ASSEMBLY */

L40:
  if (*ind == 0) {
    goto L50;
  }
  if (*x < 0.) {
    ret_val = exp(*x * *x) * 2. - ret_val;
  }
  return ret_val;
L50:
  w = *x * *x;
  t = w;
  e = w - t;
  ret_val = (.5 - e + .5) * exp(-t) * ret_val;
  if (*x < 0.) {
    ret_val = 2. - ret_val;
  }
  return ret_val;

  /*             LIMIT VALUE FOR LARGE NEGATIVE X */

L60:
  ret_val = 2.;
  if (*ind != 0) {
    ret_val = exp(*x * *x) * 2.;
  }
  return ret_val;

  /*             LIMIT VALUE FOR LARGE POSITIVE X */
  /*                       WHEN IND = 0 */

L70:
  ret_val = 0.;
  return ret_val;
} /* erfc1_ */

/* GRATIO.F -- translated by f2c (version 19980913).
   You must link the resulting object file with the libraries:
        -lf2c -lm   (in that order)
*/

/* Subroutine */ int gratio_(doublereal *a, doublereal *x, doublereal *ans,
                             doublereal *qans, integer *ind) {
  /* Initialized data */

  static doublereal acc0[3] = {5e-15, 5e-7, 5e-4};
  static doublereal d10 = -.00185185185185185;
  static doublereal d1[12] = {
      -.00347222222222222, .00264550264550265,   -9.9022633744856e-4,
      2.05761316872428e-4, -4.01877572016461e-7, -1.809855033449e-5,
      7.64916091608111e-6, -1.61209008945634e-6, 4.64712780280743e-9,
      1.37863344691572e-7, -5.7525456035177e-8,  1.19516285997781e-8};
  static doublereal d20 = .00413359788359788;
  static doublereal d2[10] = {-.00268132716049383,  7.71604938271605e-4,
                              2.0093878600823e-6,   -1.07366532263652e-4,
                              5.29234488291201e-5,  -1.27606351886187e-5,
                              3.42357873409614e-8,  1.37219573090629e-6,
                              -6.29899213838006e-7, 1.42806142060642e-7};
  static doublereal d30 = 6.49434156378601e-4;
  static doublereal d3[8] = {2.29472093621399e-4, -4.69189494395256e-4,
                             2.67720632062839e-4, -7.56180167188398e-5,
                             -2.3965051138673e-7, 1.10826541153473e-5,
                             -5.6749528269916e-6, 1.42309007324359e-6};
  static doublereal d40 = -8.61888290916712e-4;
  static doublereal d4[6] = {7.84039221720067e-4,  -2.9907248030319e-4,
                             -1.46384525788434e-6, 6.64149821546512e-5,
                             -3.96836504717943e-5, 1.13757269706784e-5};
  static doublereal d50 = -3.36798553366358e-4;
  static doublereal d5[4] = {-6.97281375836586e-5, 2.77275324495939e-4,
                             -1.99325705161888e-4, 6.79778047793721e-5};
  static doublereal big[3] = {20., 14., 10.};
  static doublereal d60 = 5.31307936463992e-4;
  static doublereal d6[2] = {-5.92166437353694e-4, 2.70878209671804e-4};
  static doublereal d70 = 3.44367606892378e-4;
  static doublereal e00[3] = {2.5e-4, .025, .14};
  static doublereal x00[3] = {31., 17., 9.7};
  static doublereal alog10 = 2.30258509299405;
  static doublereal rt2pin = .398942280401433;
  static doublereal rtpi = 1.77245385090552;
  static doublereal third = .333333333333333;
  static doublereal d0[13] = {
      .0833333333333333,    -.0148148148148148,   .00115740740740741,
      3.52733686067019e-4,  -1.78755144032922e-4, 3.91926317852244e-5,
      -2.18544851067999e-6, -1.85406221071516e-6, 8.29671134095309e-7,
      -1.76659527368261e-7, 6.7078535434015e-9,   1.02618097842403e-8,
      -4.38203601845335e-9};

  /* System generated locals */
  integer i__1;
  doublereal d__1;

  /* Builtin functions */
  double log(doublereal), exp(doublereal), sqrt(doublereal);

  /* Local variables */
  static doublereal a2nm1, b2nm1;
  extern doublereal rlog_(doublereal *);
  static doublereal twoa;
  extern doublereal rexp_(doublereal *), erfc1_(integer *, doublereal *);
  static doublereal c__, e, g, h__;
  static integer i__;
  static doublereal j, l;
  static integer m, n;
  static doublereal r__, s, t, u, w, y, z__;
  extern doublereal gamma_(doublereal *);
  static doublereal c0, c1, c2, c3, c4, c5, c6, e0, t1, x0, an, wk[20], am0,
      an0, a2n, b2n;
  extern doublereal spmpar_(integer *);
  static doublereal acc, cma, amn;
  extern doublereal erf_(doublereal *);
  static doublereal apn;
  static integer max__;
  static doublereal rta;
  static integer iop;
  static doublereal tol, sum, rtx;
  extern doublereal gam1_(doublereal *);

  /* ---------------------------------------------------------------------- */
  /*        EVALUATION OF THE INCOMPLETE GAMMA RATIO FUNCTIONS */
  /*                      P(A  ,X) AND Q(A  ,X) */

  /*                        ---------- */

  /*     IT IS ASSUMED THAT A AND X ARE NONNEGATIVE  , WHERE A AND X */
  /*     ARE NOT BOTH 0. */

  /*     ANS AND QANS ARE VARIABLES. GRATIO ASSIGNS ANS THE VALUE */
  /*     P(A  ,X) AND QANS THE VALUE Q(A  ,X). IND MAY BE ANY INTEGER. */
  /*     IF IND = 0 THEN THE USER IS REQUESTING AS MUCH ACCURACY AS */
  /*     POSSIBLE (UP TO 14 SIGNIFICANT DIGITS). OTHERWISE  , IF */
  /*     IND = 1 THEN ACCURACY IS REQUESTED TO WITHIN 1 UNIT OF THE */
  /*     6-TH SIGNIFICANT DIGIT  , AND IF IND .NE. 0  ,1 THEN ACCURACY */
  /*     IS REQUESTED TO WITHIN 1 UNIT OF THE 3RD SIGNIFICANT DIGIT. */

  /*     ERROR RETURN ... */
  /*        ANS IS ASSIGNED THE VALUE 2 WHEN A OR X IS NEGATIVE  , */
  /*     WHEN A*X = 0  , OR WHEN P(A  ,X) AND Q(A  ,X) ARE INDETERMINANT. */
  /*     P(A  ,X) AND Q(A  ,X) ARE COMPUTATIONALLY INDETERMINANT WHEN */
  /*     X IS EXCEEDINGLY CLOSE TO A AND A IS EXTREMELY LARGE. */
  /* ---------------------------------------------------------------------- */
  /*     WRITTEN BY ALFRED H. MORRIS  , JR. */
  /*        NAVAL SURFACE WEAPONS CENTER */
  /*        DAHLGREN  , VIRGINIA */
  /*     -------------------- */
  /*     .. Scalar Arguments .. */
  /*     .. */
  /*     .. Local Scalars .. */
  /*     .. */
  /*     .. Local Arrays .. */
  /*     .. */
  /*     .. External Functions .. */
  /*     .. */
  /*     .. Intrinsic Functions .. */
  /*     .. */
  /*     .. Data statements .. */
  /*     -------------------- */
  /*     -------------------- */
  /*     ALOG10 = LN(10) */
  /*     RT2PIN = 1/SQRT(2*PI) */
  /*     RTPI   = SQRT(PI) */
  /*     -------------------- */
  /*     -------------------- */
  /*     -------------------- */
  /*     -------------------- */
  /*     -------------------- */
  /*     -------------------- */
  /*     -------------------- */
  /*     -------------------- */
  /*     -------------------- */
  /*     .. */
  /*     .. Executable Statements .. */
  /*     -------------------- */
  /*     ****** E IS A MACHINE DEPENDENT CONSTANT. E IS THE SMALLEST */
  /*            FLOATING POINT NUMBER FOR WHICH 1.0 + E .GT. 1.0 . */

  e = spmpar_(&c__1);

  /*     -------------------- */
  if (*a < 0. || *x < 0.) {
    goto L430;
  }
  if (*a == 0. && *x == 0.) {
    goto L430;
  }
  if (*a * *x == 0.) {
    goto L420;
  }

  iop = *ind + 1;
  if (iop != 1 && iop != 2) {
    iop = 3;
  }
  /* Computing MAX */
  d__1 = acc0[iop - 1];
  acc = DMAX(d__1, e);
  e0 = e00[iop - 1];
  x0 = x00[iop - 1];

  /*            SELECT THE APPROPRIATE ALGORITHM */

  if (*a >= 1.) {
    goto L10;
  }
  if (*a == .5) {
    goto L390;
  }
  if (*x < 1.1) {
    goto L160;
  }
  t1 = *a * log(*x) - *x;
  u = *a * exp(t1);
  if (u == 0.) {
    goto L380;
  }
  r__ = u * (gam1_(a) + 1.);
  goto L250;

L10:
  if (*a >= big[iop - 1]) {
    goto L30;
  }
  if (*a > *x || *x >= x0) {
    goto L20;
  }
  twoa = *a + *a;
  m = (integer)twoa;
  if (twoa != (doublereal)m) {
    goto L20;
  }
  i__ = m / 2;
  if (*a == (doublereal)i__) {
    goto L210;
  }
  goto L220;
L20:
  t1 = *a * log(*x) - *x;
  r__ = exp(t1) / gamma_(a);
  goto L40;

L30:
  l = *x / *a;
  if (l == 0.) {
    goto L370;
  }
  s = .5 - l + .5;
  z__ = rlog_(&l);
  if (z__ >= 700. / *a) {
    goto L410;
  }
  y = *a * z__;
  rta = sqrt(*a);
  if (fabs(s) <= e0 / rta) {
    goto L330;
  }
  if (fabs(s) <= .4) {
    goto L270;
  }

  /* Computing 2nd power */
  d__1 = 1. / *a;
  t = d__1 * d__1;
  t1 = (((t * .75 - 1.) * t + 3.5) * t - 105.) / (*a * 1260.);
  t1 -= y;
  r__ = rt2pin * rta * exp(t1);

L40:
  if (r__ == 0.) {
    goto L420;
  }
  if (*x <= DMAX(*a, alog10)) {
    goto L50;
  }
  if (*x < x0) {
    goto L250;
  }
  goto L100;

  /*                 TAYLOR SERIES FOR P/R */

L50:
  apn = *a + 1.;
  t = *x / apn;
  wk[0] = t;
  for (n = 2; n <= 20; ++n) {
    apn += 1.;
    t *= *x / apn;
    if (t <= .001) {
      goto L70;
    }
    wk[n - 1] = t;
    /* L60: */
  }
  n = 20;

L70:
  sum = t;
  tol = acc * .5;
L80:
  apn += 1.;
  t *= *x / apn;
  sum += t;
  if (t > tol) {
    goto L80;
  }

  max__ = n - 1;
  i__1 = max__;
  for (m = 1; m <= i__1; ++m) {
    --n;
    sum += wk[n - 1];
    /* L90: */
  }
  *ans = r__ / *a * (sum + 1.);
  *qans = .5 - *ans + .5;
  return 0;

  /*                 ASYMPTOTIC EXPANSION */

L100:
  amn = *a - 1.;
  t = amn / *x;
  wk[0] = t;
  for (n = 2; n <= 20; ++n) {
    amn += -1.;
    t *= amn / *x;
    if (fabs(t) <= .001) {
      goto L120;
    }
    wk[n - 1] = t;
    /* L110: */
  }
  n = 20;

L120:
  sum = t;
L130:
  if (fabs(t) <= acc) {
    goto L140;
  }
  amn += -1.;
  t *= amn / *x;
  sum += t;
  goto L130;

L140:
  max__ = n - 1;
  i__1 = max__;
  for (m = 1; m <= i__1; ++m) {
    --n;
    sum += wk[n - 1];
    /* L150: */
  }
  *qans = r__ / *x * (sum + 1.);
  *ans = .5 - *qans + .5;
  return 0;

  /*             TAYLOR SERIES FOR P(A  ,X)/X**A */

L160:
  an = 3.;
  c__ = *x;
  sum = *x / (*a + 3.);
  tol = acc * 3. / (*a + 1.);
L170:
  an += 1.;
  c__ = -c__ * (*x / an);
  t = c__ / (*a + an);
  sum += t;
  if (fabs(t) > tol) {
    goto L170;
  }
  j = *a * *x * ((sum / 6. - .5 / (*a + 2.)) * *x + 1. / (*a + 1.));

  z__ = *a * log(*x);
  h__ = gam1_(a);
  g = h__ + 1.;
  if (*x < .25) {
    goto L180;
  }
  if (*a < *x / 2.59) {
    goto L200;
  }
  goto L190;
L180:
  if (z__ > -.13394) {
    goto L200;
  }

L190:
  w = exp(z__);
  *ans = w * g * (.5 - j + .5);
  *qans = .5 - *ans + .5;
  return 0;

L200:
  l = rexp_(&z__);
  w = l + .5 + .5;
  *qans = (w * j - l) * g - h__;
  if (*qans < 0.) {
    goto L380;
  }
  *ans = .5 - *qans + .5;
  return 0;

  /*             FINITE SUMS FOR Q WHEN A .GE. 1 */
  /*                 AND 2*A IS AN INTEGER */

L210:
  sum = exp(-(*x));
  t = sum;
  n = 1;
  c__ = 0.;
  goto L230;

L220:
  rtx = sqrt(*x);
  sum = erfc1_(&c__0, &rtx);
  t = exp(-(*x)) / (rtpi * rtx);
  n = 0;
  c__ = -.5;

L230:
  if (n == i__) {
    goto L240;
  }
  ++n;
  c__ += 1.;
  t = *x * t / c__;
  sum += t;
  goto L230;
L240:
  *qans = sum;
  *ans = .5 - *qans + .5;
  return 0;

  /*              CONTINUED FRACTION EXPANSION */

L250:
  /* Computing MAX */
  d__1 = e * 5.;
  tol = DMAX(d__1, acc);
  a2nm1 = 1.;
  a2n = 1.;
  b2nm1 = *x;
  b2n = *x + (1. - *a);
  c__ = 1.;
L260:
  a2nm1 = *x * a2n + c__ * a2nm1;
  b2nm1 = *x * b2n + c__ * b2nm1;
  am0 = a2nm1 / b2nm1;
  c__ += 1.;
  cma = c__ - *a;
  a2n = a2nm1 + cma * a2n;
  b2n = b2nm1 + cma * b2n;
  an0 = a2n / b2n;
  if ((d__1 = an0 - am0, fabs(d__1)) >= tol * an0) {
    goto L260;
  }

  *qans = r__ * an0;
  *ans = .5 - *qans + .5;
  return 0;

  /*                GENERAL TEMME EXPANSION */

L270:
  if (fabs(s) <= e * 2. && *a * e * e > .00328) {
    goto L430;
  }
  c__ = exp(-y);
  d__1 = sqrt(y);
  w = erfc1_(&c__1, &d__1) * .5;
  u = 1. / *a;
  z__ = sqrt(z__ + z__);
  if (l < 1.) {
    z__ = -z__;
  }
  if ((i__1 = iop - 2) < 0) {
    goto L280;
  } else if (i__1 == 0) {
    goto L290;
  } else {
    goto L300;
  }

L280:
  if (fabs(s) <= .001) {
    goto L340;
  }
  c0 = ((((((((((((d0[12] * z__ + d0[11]) * z__ + d0[10]) * z__ + d0[9]) * z__ +
                d0[8]) *
                   z__ +
               d0[7]) *
                  z__ +
              d0[6]) *
                 z__ +
             d0[5]) *
                z__ +
            d0[4]) *
               z__ +
           d0[3]) *
              z__ +
          d0[2]) *
             z__ +
         d0[1]) *
            z__ +
        d0[0]) *
           z__ -
       third;
  c1 = (((((((((((d1[11] * z__ + d1[10]) * z__ + d1[9]) * z__ + d1[8]) * z__ +
               d1[7]) *
                  z__ +
              d1[6]) *
                 z__ +
             d1[5]) *
                z__ +
            d1[4]) *
               z__ +
           d1[3]) *
              z__ +
          d1[2]) *
             z__ +
         d1[1]) *
            z__ +
        d1[0]) *
           z__ +
       d10;
  c2 = (((((((((d2[9] * z__ + d2[8]) * z__ + d2[7]) * z__ + d2[6]) * z__ +
             d2[5]) *
                z__ +
            d2[4]) *
               z__ +
           d2[3]) *
              z__ +
          d2[2]) *
             z__ +
         d2[1]) *
            z__ +
        d2[0]) *
           z__ +
       d20;
  c3 =
      (((((((d3[7] * z__ + d3[6]) * z__ + d3[5]) * z__ + d3[4]) * z__ + d3[3]) *
             z__ +
         d3[2]) *
            z__ +
        d3[1]) *
           z__ +
       d3[0]) *
          z__ +
      d30;
  c4 = (((((d4[5] * z__ + d4[4]) * z__ + d4[3]) * z__ + d4[2]) * z__ + d4[1]) *
            z__ +
        d4[0]) *
           z__ +
       d40;
  c5 = (((d5[3] * z__ + d5[2]) * z__ + d5[1]) * z__ + d5[0]) * z__ + d50;
  c6 = (d6[1] * z__ + d6[0]) * z__ + d60;
  t = ((((((d70 * u + c6) * u + c5) * u + c4) * u + c3) * u + c2) * u + c1) *
          u +
      c0;
  goto L310;

L290:
  c0 = (((((d0[5] * z__ + d0[4]) * z__ + d0[3]) * z__ + d0[2]) * z__ + d0[1]) *
            z__ +
        d0[0]) *
           z__ -
       third;
  c1 = (((d1[3] * z__ + d1[2]) * z__ + d1[1]) * z__ + d1[0]) * z__ + d10;
  c2 = d2[0] * z__ + d20;
  t = (c2 * u + c1) * u + c0;
  goto L310;

L300:
  t = ((d0[2] * z__ + d0[1]) * z__ + d0[0]) * z__ - third;

L310:
  if (l < 1.) {
    goto L320;
  }
  *qans = c__ * (w + rt2pin * t / rta);
  *ans = .5 - *qans + .5;
  return 0;
L320:
  *ans = c__ * (w - rt2pin * t / rta);
  *qans = .5 - *ans + .5;
  return 0;

  /*               TEMME EXPANSION FOR L = 1 */

L330:
  if (*a * e * e > .00328) {
    goto L430;
  }
  c__ = .5 - y + .5;
  w = (.5 - sqrt(y) * (.5 - y / 3. + .5) / rtpi) / c__;
  u = 1. / *a;
  z__ = sqrt(z__ + z__);
  if (l < 1.) {
    z__ = -z__;
  }
  if ((i__1 = iop - 2) < 0) {
    goto L340;
  } else if (i__1 == 0) {
    goto L350;
  } else {
    goto L360;
  }

L340:
  c0 = ((((((d0[6] * z__ + d0[5]) * z__ + d0[4]) * z__ + d0[3]) * z__ + d0[2]) *
             z__ +
         d0[1]) *
            z__ +
        d0[0]) *
           z__ -
       third;
  c1 = (((((d1[5] * z__ + d1[4]) * z__ + d1[3]) * z__ + d1[2]) * z__ + d1[1]) *
            z__ +
        d1[0]) *
           z__ +
       d10;
  c2 = ((((d2[4] * z__ + d2[3]) * z__ + d2[2]) * z__ + d2[1]) * z__ + d2[0]) *
           z__ +
       d20;
  c3 = (((d3[3] * z__ + d3[2]) * z__ + d3[1]) * z__ + d3[0]) * z__ + d30;
  c4 = (d4[1] * z__ + d4[0]) * z__ + d40;
  c5 = (d5[1] * z__ + d5[0]) * z__ + d50;
  c6 = d6[0] * z__ + d60;
  t = ((((((d70 * u + c6) * u + c5) * u + c4) * u + c3) * u + c2) * u + c1) *
          u +
      c0;
  goto L310;

L350:
  c0 = (d0[1] * z__ + d0[0]) * z__ - third;
  c1 = d1[0] * z__ + d10;
  t = (d20 * u + c1) * u + c0;
  goto L310;

L360:
  t = d0[0] * z__ - third;
  goto L310;

  /*                     SPECIAL CASES */

L370:
  *ans = 0.;
  *qans = 1.;
  return 0;

L380:
  *ans = 1.;
  *qans = 0.;
  return 0;

L390:
  if (*x >= .25) {
    goto L400;
  }
  d__1 = sqrt(*x);
  *ans = erf_(&d__1);
  *qans = .5 - *ans + .5;
  return 0;
L400:
  d__1 = sqrt(*x);
  *qans = erfc1_(&c__0, &d__1);
  *ans = .5 - *qans + .5;
  return 0;

L410:
  if (fabs(s) <= e * 2.) {
    goto L430;
  }
L420:
  if (*x <= *a) {
    goto L370;
  }
  goto L380;

  /*                     ERROR RETURN */

L430:
  *ans = 2.;
  return 0;
} /* gratio_ */

/* CUMGAM.F -- translated by f2c (version 19980913).
   You must link the resulting object file with the libraries:
        -lf2c -lm   (in that order)
*/

/* Table of constant values */

/* Subroutine */ int cumgam_(doublereal *x, doublereal *a, doublereal *cum,
                             doublereal *ccum) {
  extern /* Subroutine */ int gratio_(doublereal *, doublereal *, doublereal *,
                                      doublereal *, integer *);

  /* ********************************************************************** */

  /*     SUBROUTINE CUMGAM(X  ,A  ,CUM  ,CCUM) */
  /*           Double precision cUMulative incomplete GAMma distribution */

  /*                              Function */

  /*     Computes   the  cumulative        of    the     incomplete   gamma */
  /*     distribution  , i.e.  , the integral from 0 to X of */
  /*          (1/GAM(A))*EXP(-T)*T**(A-1) DT */
  /*     where GAM(A) is the complete gamma function of A  , i.e.  , */
  /*          GAM(A) = integral from 0 to infinity of */
  /*                    EXP(-T)*T**(A-1) DT */

  /*                              Arguments */

  /*     X --> The upper limit of integration of the incomplete gamma. */
  /*                                                X is DOUBLE PRECISION */

  /*     A --> The shape parameter of the incomplete gamma. */
  /*                                                A is DOUBLE PRECISION */

  /*     CUM <-- Cumulative incomplete gamma distribution. */
  /*                                        CUM is DOUBLE PRECISION */

  /*     CCUM <-- Compliment of Cumulative incomplete gamma distribution. */
  /*                                                CCUM is DOUBLE PRECISIO */

  /*                              Method */

  /*     Calls the routine GRATIO. */

  /* ********************************************************************** */

  /*     .. */
  /*     .. Scalar Arguments .. */
  /*     .. */
  /*     .. External Routines .. */
  /*     .. */
  /*     .. Executable Statements .. */
  if (!(*x <= 0.)) {
    goto L10;
  }
  *cum = 0.;
  *ccum = 1.;
  return 0;
L10:
  gratio_(a, x, cum, ccum, &c__0);
  /*     Call gratio routine */
  return 0;
} /* cumgam_ */

/* CUMCHI.F -- translated by f2c (version 19980913).
   You must link the resulting object file with the libraries:
        -lf2c -lm   (in that order)
*/

/* Subroutine */
int cumchi_(doublereal *x, doublereal *df, doublereal *cum, doublereal *ccum) {
  static doublereal a;
  extern /* Subroutine */ int cumgam_(doublereal *, doublereal *, doublereal *,
                                      doublereal *);
  static doublereal xx;

  /* ********************************************************************** */

  /*     SUBROUTINE FUNCTION CUMCHI(X  ,DF  ,CUM  ,CCUM) */
  /*             CUMulative of the CHi-square distribution */

  /*                              Function */

  /*     Calculates the cumulative chi-square distribution. */

  /*                              Arguments */

  /*     X       --> Upper limit of integration of the */
  /*                 chi-square distribution. */
  /*                                                 X is DOUBLE PRECISION */

  /*     DF      --> Degrees of freedom of the */
  /*                 chi-square distribution. */
  /*                                                 DF is DOUBLE PRECISION */

  /*     CUM <-- Cumulative chi-square distribution. */
  /*                                                 CUM is DOUBLE PRECISIO */

  /*     CCUM <-- Compliment of Cumulative chi-square distribution. */
  /*                                                 CCUM is DOUBLE PRECISI */

  /*                              Method */

  /*     Calls incomplete gamma function (CUMGAM) */

  /* ********************************************************************** */
  /*     .. Scalar Arguments .. */
  /*     .. */
  /*     .. Local Scalars .. */
  /*     .. */
  /*     .. External Subroutines .. */
  /*     .. */
  /*     .. Executable Statements .. */
  a = *df * .5;
  xx = *x * .5;
  cumgam_(&xx, &a, cum, ccum);
  return 0;
} /* cumchi_ */

/* CUMCHN.F -- translated by f2c (version 19980913).
   You must link the resulting object file with the libraries:
        -lf2c -lm   (in that order)
*/

/* Subroutine */ int cumchn_(doublereal *x, doublereal *df, doublereal *pnonc,
                             doublereal *cum, doublereal *ccum) {
  /* Initialized data */

  static doublereal eps = 1e-8;

  /* System generated locals */
  doublereal d__1;

  /* Builtin functions */
  double log(doublereal), exp(doublereal);

  /* Local variables */
  static doublereal term, chid2;
  static integer i__;
  static doublereal lfact;
  static integer icent;
  static doublereal pcent, xnonc, pterm;
  extern doublereal alngam_(doublereal *);
  static doublereal centaj;
  extern /* Subroutine */ int cumchi_(doublereal *, doublereal *, doublereal *,
                                      doublereal *);
  static doublereal lcntaj, wt, sumadj, centwt, lcntwt, adj, sum, dfd2;

  /* *********************************************************************** */

  /*     SUBROUTINE CUMCHN(X  ,DF  ,PNONC  ,CUM  ,CCUM) */
  /*             CUMulative of the Non-central CHi-square distribution */

  /*                               Function */

  /*     Calculates     the       cumulative      non-central    chi-square */
  /*     distribution  , i.e.  ,  the probability   that  a   random   variable
   */
  /*     which    follows  the  non-central chi-square  distribution  ,  with */
  /*     non-centrality  parameter    PNONC  and   continuous  degrees   of */
  /*     freedom DF  , is less than or equal to X. */

  /*                              Arguments */

  /*     X       --> Upper limit of integration of the non-central */
  /*                 chi-square distribution. */
  /*                                                 X is DOUBLE PRECISION */

  /*     DF      --> Degrees of freedom of the non-central */
  /*                 chi-square distribution. */
  /*                                                 DF is DOUBLE PRECISION */

  /*     PNONC   --> Non-centrality parameter of the non-central */
  /*                 chi-square distribution. */
  /*                                                 PNONC is DOUBLE PRECIS */

  /*     CUM <-- Cumulative non-central chi-square distribution. */
  /*                                                 CUM is DOUBLE PRECISIO */

  /*     CCUM <-- Compliment of Cumulative non-central chi-square distribut */
  /*                                                 CCUM is DOUBLE PRECISI */

  /*                                Method */

  /*     Uses  formula  26.4.25   of  Abramowitz  and  Stegun  , Handbook  of */
  /*     Mathematical    Functions  ,  US   NBS   (1966)    to calculate  the */
  /*     non-central chi-square. */

  /*                                Variables */

  /*     EPS     --- Convergence criterion.  The sum stops when a */
  /*                 term is less than EPS*SUM. */
  /*                                                 EPS is DOUBLE PRECISIO */

  /*     CCUM <-- Compliment of Cumulative non-central */
  /*              chi-square distribution. */
  /*                                                 CCUM is DOUBLE PRECISI */

  /* *********************************************************************** */

  /*     .. Scalar Arguments .. */
  /*     .. */
  /*     .. Local Scalars .. */
  /*     .. */
  /*     .. External Functions .. */
  /*     .. */
  /*     .. External Subroutines .. */
  /*     .. */
  /*     .. Intrinsic Functions .. */
  /*     .. */
  /*     .. Statement Functions .. */
  /*     .. */
  /*     .. Data statements .. */
  /*     .. */
  /*     .. Statement Function definitions .. */
  /*     .. */

  if (!(*x <= 0.)) {
    goto L10;
  }
  *cum = 0.;
  *ccum = 1.;
  return 0;
L10:
  if (!(*pnonc <= 1e-10)) {
    goto L20;
  }

  /*     When non-centrality parameter is (essentially) zero  , */
  /*     use cumulative chi-square distribution */

  cumchi_(x, df, cum, ccum);
  return 0;
L20:
  xnonc = *pnonc / 2.;
  /* *********************************************************************** */

  /*     The following code calcualtes the weight  , chi-square  , and */
  /*     adjustment term for the central term in the infinite series. */
  /*     The central term is the one in which the poisson weight is */
  /*     greatest.  The adjustment term is the amount that must */
  /*     be subtracted from the chi-square to move up two degrees */
  /*     of freedom. */

  /* *********************************************************************** */
  icent = (integer)xnonc;
  if (icent == 0) {
    icent = 1;
  }
  chid2 = *x / 2.;

  /*     Calculate central weight term */

  d__1 = (doublereal)(icent + 1);
  lfact = alngam_(&d__1);
  lcntwt = -xnonc + icent * log(xnonc) - lfact;
  centwt = exp(lcntwt);

  /*     Calculate central chi-square */

  d__1 = *df + (doublereal)icent * 2.;
  cumchi_(x, &d__1, &pcent, ccum);

  /*     Calculate central adjustment term */

  dfd2 = (*df + (doublereal)icent * 2.) / 2.;
  d__1 = dfd2 + 1.;
  lfact = alngam_(&d__1);
  lcntaj = dfd2 * log(chid2) - chid2 - lfact;
  centaj = exp(lcntaj);
  sum = centwt * pcent;
  /* *********************************************************************** */

  /*     Sum backwards from the central term towards zero. */
  /*     Quit whenever either */
  /*     (1) the zero term is reached  , or */
  /*     (2) the term gets small relative to the sum  , or */

  /* *********************************************************************** */
  sumadj = 0.;
  adj = centaj;
  wt = centwt;
  i__ = icent;

  goto L40;
L30:
  if (sum < 1e-20 || term < eps * sum || i__ == 0) {
    goto L50;
  }
L40:
  dfd2 = (*df + (doublereal)i__ * 2.) / 2.;

  /*     Adjust chi-square for two fewer degrees of freedom. */
  /*     The adjusted value ends up in PTERM. */

  adj = adj * dfd2 / chid2;
  sumadj += adj;
  pterm = pcent + sumadj;

  /*     Adjust poisson weight for J decreased by one */

  wt *= i__ / xnonc;
  term = wt * pterm;
  sum += term;
  --i__;
  goto L30;
L50:
  sumadj = centaj;
  /* *********************************************************************** */

  /*     Now sum forward from the central term towards infinity. */
  /*     Quit when either */
  /*     (1) the term gets small relative to the sum  , or */

  /* *********************************************************************** */
  adj = centaj;
  wt = centwt;
  i__ = icent;

  goto L70;
L60:
  if (sum < 1e-20 || term < eps * sum) {
    goto L80;
  }

  /*     Update weights for next higher J */

L70:
  wt *= xnonc / (i__ + 1);

  /*     Calculate PTERM and add term to sum */

  pterm = pcent - sumadj;
  term = wt * pterm;
  sum += term;

  /*     Update adjustment term for DF for next iteration */

  ++i__;
  dfd2 = (*df + (doublereal)i__ * 2.) / 2.;
  adj = adj * chid2 / dfd2;
  sumadj += adj;
  goto L60;
L80:
  *cum = sum;
  *ccum = .5 - *cum + .5;

  return 0;
} /* cumchn_ */

/************************************************************
 * A nicer interface for griffin
 * cumchn
 ************************************************************/

Err cumchn(double x, double df, double pnonc, double *cum, double *ccum) {
  cumchn_(&x, &df, &pnonc, cum, ccum);
  return NULL;
}

Err dchi2nc(double x, double df, double pnonc, double *cum) {
  double *ccum;
  ccum = (double *)calloc(1, sizeof(double));
  cumchn_(&x, &df, &pnonc, cum, ccum);
  return NULL;
}