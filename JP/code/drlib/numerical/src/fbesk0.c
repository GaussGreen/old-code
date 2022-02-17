#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif

static Mfloat PROTO(l_dbsk0e,(Mfloat x));

	/*
	 * SERIES FOR BK0 ON THE INTERVAL 0.0  TO  4.00000D+00 WITH WEIGHTED
	 * ERROR        3.08D-33 LOG WEIGHTED ERROR        32.51 SIGNIFICANT
	 * FIGURES REQD. 32.05 DECIMAL PLACES REQUIRED   33.11
	 */
static Mfloat lv_bk0cs[] = {
	 -.353273932339027687201140060063e-1,
	 .344289899924628486886344927529e0,
	 .359799365153615016265721303687e-1,
	 .126461541144692592338479508674e-2,
	 .228621210311945178608269830298e-4,
	 .253479107902614945730790013428e-6,
	 .190451637722020885897214059381e-8,
	 .103496952576336245851008317853e-10,
	 .42598161427910825765244532717e-13,
	 .13744654358807508969423832544e-15,
	 .357089652850837359099688597333e-18,
	 .763164366011643737667498666667e-21,
	 .136542498844078185908053333333e-23,
	 .20752752669066680832e-26,
	 .27128142180729856e-29,
	 .308259388791466666666666666667e-32
};
	/*
	 * SERIES FOR AK0 ON THE INTERVAL 1.25000D-01 TO  5.00000D-01 WITH
	 * WEIGHTED ERROR        2.85D-32 LOG WEIGHTED ERROR        31.54
	 * SIGNIFICANT FIGURES REQD. 30.19 DECIMAL PLACES REQUIRED   32.33
	 */
static Mfloat lv_ak0cs[] = {
	 -.764394790332794142408297827009e-1,
	 -.223565260569981905202309555079e-1,
	 .773418115469385823530061817405e-3,
	 -.428100668888609946445214643542e-4,
	 .308170017386297474365001482666e-5,
	 -.263936722200966497406744889272e-6,
	 .256371303640346920629408826574e-7,
	 -.274270554990020126385721191524e-8,
	 .31694296580974995920808328734e-9,
	 -.390235328696218414160106571796e-10,
	 .506804069818857540205009212729e-11,
	 -.688957474100787067954171355798e-12,
	 .974497849782591769138820133683e-13,
	 -.142733284188454850538985534012e-13,
	 .215641257102146303955806297653e-14,
	 -.334965425514956277218878205853e-15,
	 .53352602169529116921452803926e-16,
	 -.869366998089075380763962237884e-17,
	 .144640434786221222788776344235e-17,
	 -.245288982550012968240467875157e-18,
	 .42337545262321715728217063424e-19,
	 -.742794652645446419569534129493e-20,
	 .13231505293926668662779674624e-20,
	 -.23905871647396494513359814656e-21,
	 .437682758592322614016571255467e-22,
	 -.811370060734511805933901141333e-23,
	 .152181991383217295831037815467e-23,
	 -.288604194148339777023595861333e-24,
	 .553062066705471797999261013333e-25,
	 -.107037732924989872859163306667e-25,
	 .209108689314238430029632853333e-26,
	 -.412171372364620382741026133333e-27,
	 .819348397112130764013568e-28,
	 -.164200027545929772678075733333e-28,
	 .331614328148022719589034666667e-29,
	 -.674686364414529594108586666667e-30,
	 .138242914631842467763541333333e-30,
	 -.285187416735983257081173333333e-31
};
	/*
	 * SERIES FOR AK02 ON THE INTERVAL 0.0  TO  1.25000D-01 WITH WEIGHTED
	 * ERROR        2.30D-32 LOG WEIGHTED ERROR        31.64 SIGNIFICANT
	 * FIGURES REQD. 29.68 DECIMAL PLACES REQUIRED   32.40
	 */
static Mfloat lv_ak02cs[] = {
	 -.120186982630759223983934621245e-1,
	 -.917485269102569531065256107571e-2,
	 .144455093177500582104884387806e-3,
	 -.401361417543570972867102107788e-5,
	 .156783181085231067259034899033e-6,
	 -.777011043852173771031579975446e-8,
	 .461118257617971788253313052959e-9,
	 -.315859299786056577052666580331e-10,
	 .243501803936504112783588781433e-11,
	 -.207433138739834789770985337351e-12,
	 .192578728058991708474273650469e-13,
	 -.192755480583895610360034718222e-14,
	 .206219802919781827828523786964e-15,
	 -.234168511757924240260364019507e-16,
	 .280590281064304224681517882846e-17,
	 -.353050763116180794581548246357e-18,
	 .464529542293510826742421633707e-19,
	 -.636862594134426647392205346133e-20,
	 .90695213109865155676223488e-21,
	 -.1337974785423690739845005312e-21,
	 .203983602185995231552208896e-22,
	 -.320702748136784050006086997333e-23,
	 .518974441366230996362635946667e-24,
	 -.8629501497540572192964608e-25,
	 .14721611831025598552080384e-25,
	 -.2573069023867011283812352e-26,
	 .460177408664351658737664e-27,
	 -.841155532420109373713066666667e-28,
	 .156980630663536893930154666667e-28,
	 -.29882264530057577889792e-29,
	 .579683137521683652061866666667e-30,
	 -.114503599434768133215573333333e-30,
	 .230126659424968280200533333333e-31
};
		
/* Structured by FOR_STRUCT, v0.2, on 08/19/90 at 12:15:18
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  BSK0/DBSK0 (Single/float precision version)

    Computer:   FORC/float

    Revised:    November 29, 1988

    Purpose:    Evaluate the modified Bessel function of the third kind
                of order zero.

    Usage:      BSK0(X)

    Arguments:
       X      - Argument for which the function value is desired.
                (Input)
       BSK0   - Function value.  (Output)

    Remark:
       Informational error:
       Type Code
         2   1  The function underflows because X is too large.

    GAMS:       C10b1

    Chapter:    SFUN/LIBRARY Bessel Functions

    Copyright:  1984 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
Mfloat imsl_f_bessel_K0(Mfloat x)
#else
Mfloat imsl_f_bessel_K0(x)
	Mfloat          x;
#endif
{
	Mfloat          dbsk0_v, y;
	static Mint     ntk0 = 0;
	static Mfloat   xsml = 0.0e0;
	static Mfloat   xmax = 0.0e0;

	E1PSH("imsl_f_bessel_K0", "imsl_d_bessel_K0");
	dbsk0_v = imsl_amach(6);

	if (ntk0 == 0) {
		ntk0 = imsl_inits(lv_bk0cs, 16, 0.1*imsl_amach(3));
		xsml = sqrt(F_FOUR * imsl_amach(3));
		xmax = -log(imsl_amach(1));
		xmax += -F_HALF * xmax * log(xmax) / (xmax + F_HALF) - 0.01e0;
	}
	if (x <= F_ZERO) {
                /* X = %(r1) must be greater than zero. */
		imsl_e1str(1, x);
		imsl_ermes(IMSL_TERMINAL, IMSL_NON_POSITIVE_ARGUMENT);
	} else if (x <= F_TWO) {
		y = F_ZERO;
		if (x > xsml)
			y = x * x;
		dbsk0_v = -log(F_HALF * x) * imsl_f_bessel_I0(x) - .25e0 + imsl_csevl(F_HALF *
						    y - F_ONE, lv_bk0cs, ntk0);
	} else {
		if (x <= xmax) {
			dbsk0_v = exp(-x) * l_dbsk0e(x);
		} else {
			dbsk0_v = F_ZERO;
                            /* The function underflows because X = %(r1) is
                               greater than %(r2).  The result is set to zero. */
			imsl_e1str(1, x);
			imsl_e1str(2, xmax);
			imsl_ermes(IMSL_ALERT, IMSL_LARGE_ARG_UNDERFLOW);
		}
	}
	E1POP("imsl_f_bessel_K0", "imsl_d_bessel_K0");
	return (dbsk0_v);
}				/* end of function */



/* Structured by FOR_STRUCT, v0.2, on 08/20/90 at 15:32:27
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  BSK0E/DBSK0E (Single/float precision version)
    Computer:   FORC/float

    Revised:    November 29, 1988

    Purpose:    Evaluate the exponentially scaled modified Bessel
                function of the third kind of order zero.

    Usage:      BSK0E(X)

    Arguments:
       X      - Argument for which the function value is desired.
                (Input)
       BSK0E  - Function value.  (Output)

    GAMS:       C10b1

    Chapter:    SFUN/LIBRARY Bessel Functions

    Copyright:  1984 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static Mfloat l_dbsk0e(Mfloat x)
#else
static Mfloat l_dbsk0e(x)
	Mfloat          x;
#endif
{
	Mfloat           eta;
	Mfloat           dbsk0e_v, y;
	static Mint      ntk0 = 0;
	static Mint      ntak0 = 0;
	static Mint      ntak02 = 0;
	static Mfloat    xsml = 0.0e0;

	imsl_e1psh("l_dbsk0e");
	dbsk0e_v = imsl_amach(6);

	if (ntk0 == 0) {
		eta = 0.1 * (Mfloat) (imsl_amach(3));
		ntk0 = imsl_inits(lv_bk0cs, 16, eta);
		ntak0 = imsl_inits(lv_ak0cs, 38, eta);
		ntak02 = imsl_inits(lv_ak02cs, 33, eta);
		xsml = sqrt(F_FOUR * imsl_amach(3));
	}
	if (x <= F_ZERO) {
                /* X = %(r1) must be greater than zero. */
		imsl_e1str(1, x);
		imsl_ermes(IMSL_TERMINAL, IMSL_NON_POSITIVE_ARGUMENT);
	} else if (x <= F_TWO) {
		y = F_ZERO;
		if (x > xsml)
			y = x * x;
		dbsk0e_v = exp(x) * (-log(F_HALF * x) * imsl_f_bessel_I0(x) - .25e0 +
				imsl_csevl(F_HALF * y - F_ONE, lv_bk0cs, ntk0));
	} else if (x <= F_EIGHT) {
		dbsk0e_v = (1.25e0 + imsl_csevl((16.0e0 / x - F_FIVE) / F_THREE, lv_ak0cs,
						 ntak0)) / sqrt(x);
	} else {
		dbsk0e_v = (1.25e0 + imsl_csevl(16.0e0 / x - F_ONE, lv_ak02cs, ntak02)) /
			sqrt(x);
	}
	imsl_e1pop("l_dbsk0e");
	return (dbsk0e_v);
}
