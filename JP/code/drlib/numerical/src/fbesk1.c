#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif

static Mfloat PROTO(l_dbsk1e,(Mfloat x));

	/*
	 * SERIES FOR BK1 ON THE INTERVAL 0.0  TO  4.00000D+00 WITH WEIGHTED
	 * ERROR        9.16D-32 LOG WEIGHTED ERROR        31.04 SIGNIFICANT
	 * FIGURES REQD. 30.61 DECIMAL PLACES REQUIRED   31.64
	 */
static Mfloat lv_bk1cs[] = {
	 .253002273389477705325311208685e-1,
	 -.353155960776544875667238316918e0,
	 -.1226111808226571482347906793e0,
	 -.697572385963986435018129202961e-2,
	 -.17302889575130520630176507369e-3,
	 -.243340614156596823496007350302e-5,
	 -.221338763073472585583152525451e-7,
	 -.141148839263352776109583302126e-9,
	 -.666690169419932900608537512644e-12,
	 -.242744985051936593392631968649e-14,
	 -.7023863479386287597178379712e-17,
	 -.165432751551009946754910293333e-19,
	 -.323383474599444919918933333333e-22,
	 -.533127505292652749994666666667e-25,
	 -.751304071621572266666666666667e-28,
	 -.915508571765418666666666666667e-31
};
	/*
	 * SERIES FOR AK1 ON THE INTERVAL 1.25000D-01 TO  5.00000D-01 WITH
	 * WEIGHTED ERROR        3.07D-32 LOG WEIGHTED ERROR        31.51
	 * SIGNIFICANT FIGURES REQD. 30.71 DECIMAL PLACES REQUIRED   32.30
	 */
static Mfloat lv_ak1cs[] = {
	 .274431340697388296952576662273e0,
	 .757198995319936781708923781493e-1,
	 -.144105155647540612298531161756e-2,
	 .66501169551257479394251385477e-4,
	 -.436998470952014076605808450892e-5,
	 .354027749976305267994171390085e-6,
	 -.331116377929329202089826882457e-7,
	 .34459775819010534532311499771e-8,
	 -.389893234747542710489819374928e-9,
	 .47208197504658356400947449339e-10,
	 -.604783566287535623453735915629e-11,
	 .812849487486587478881938379857e-12,
	 -.11386945747147891428923915951e-12,
	 .165403584084622823259729482051e-13,
	 -.248090256770688482215160104405e-14,
	 .382923789070240969484292272992e-15,
	 -.606473410400124181877682103774e-16,
	 .983242562326486160381940046507e-17,
	 -.162841687382843800356666201156e-17,
	 .275015364967526237182841203371e-18,
	 -.47289666463953250924281069568e-19,
	 .826815000281099327223920503467e-20,
	 -.146814051366249563371939648853e-20,
	 .264476392692082459780858948267e-21,
	 -.482901575648563878979698688e-22,
	 .892930207436101301806563328e-23,
	 -.167083971689725171769977514667e-23,
	 .316164560340406949313686186667e-24,
	 -.604620553122749891065064106667e-25,
	 .116787989420427327007184213333e-25,
	 -.2277374158265399623286784e-26,
	 .448110973007736757953058133333e-27,
	 -.88932884769020194062336e-28,
	 .17794680018850275131392e-28,
	 -.358845559673290958219946666667e-29,
	 .7290629049269425799168e-30,
	 -.14918449845546227073024e-30,
	 .307365738729342763008e-31
};
	/*
	 * SERIES FOR AK12 ON THE INTERVAL 0.0  TO  1.25000D-01 WITH WEIGHTED
	 * ERROR        2.41D-32 LOG WEIGHTED ERROR        31.62 SIGNIFICANT
	 * FIGURES REQD. 30.25 DECIMAL PLACES REQUIRED   32.38
	 */
static Mfloat lv_ak12cs[] = {
	 .63793083437390010366004885341e-1,
	 .283288781304972093583503028471e-1,
	 -.247537067390525034541454556673e-3,
	 .577197245160724882047097662576e-5,
	 -.206893921953654830274553319655e-6,
	 .973998344138180418030921309789e-8,
	 -.558533614038062498468889551113e-9,
	 .373299663404618524022121285473e-10,
	 -.282505196102322544513506575493e-11,
	 .237201900248414417364349695549e-12,
	 -.217667738799175397926830166794e-13,
	 .215791416161603245393956268971e-14,
	 -.229019693071826927599155133815e-15,
	 .258288572982327496191993956523e-16,
	 -.307675264126846318762109817344e-17,
	 .38514877212804915970948968448e-18,
	 -.50447948976415289771172825088e-19,
	 .6888673850418544237018292224e-20,
	 -.977504154195011830300213248e-21,
	 .143741621852383646100165973333e-21,
	 -.218505949734434737349973333333e-22,
	 .34262456218092206316453888e-23,
	 -.5531064394246408232501248e-24,
	 .917660150568599540378282666667e-25,
	 -.156228720361802491144874666667e-25,
	 .272541937548433313234944e-26,
	 -.486567491007482799237802666667e-27,
	 .887938855272350258735786666667e-28,
	 -.165458591803925754893653333333e-28,
	 .3145111321357848674304e-29,
	 -.6092998312193127612416e-30,
	 .1202021939369815834624e-30,
	 -.241293080145940884138666666667e-31
};

/* Structured by FOR_STRUCT, v0.2, on 08/19/90 at 12:16:12
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  BSK1/DBSK1 (Single/float precision version)

    Computer:   FORC/float

    Revised:    November 29, 1988

    Purpose:    Evaluate the modified Bessel function of the third kind
                of order one.

    Usage:      BSK1(X)

    Arguments:
       X      - Argument for which the function value is desired.
                (Input)
       BSK1   - Function value.  (Output)

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
Mfloat imsl_f_bessel_K1(Mfloat x)
#else
Mfloat imsl_f_bessel_K1(x)
	Mfloat          x;
#endif
{
	Mfloat          dbsk1_v, y;
	static Mint     ntk1 = 0;
	static Mfloat   xmin = 0.0e0;
	static Mfloat   xsml = 0.0e0;
	static Mfloat   xmax = 0.0e0;

	E1PSH("imsl_f_bessel_K1", "imsl_d_bessel_K1");
	dbsk1_v = imsl_amach(6);

	if (ntk1 == 0) {
		ntk1 = imsl_inits(lv_bk1cs, 16, (Mfloat) (0.1 * (Mfloat) (imsl_amach(3))));
		xmin = exp(imsl_f_max(log(imsl_amach(1)), -log(imsl_amach(2))) +
				0.01e0);
		xsml = sqrt(F_FOUR * imsl_amach(3));
		xmax = -log(imsl_amach(1));
		xmax += -F_HALF * xmax * log(xmax) / (xmax + F_HALF) - 0.01e0;
	}
	if (x <= F_ZERO) {
                /* X = %(r1) must be greater than zero. */
		imsl_e1str(1, x);
		imsl_ermes(IMSL_TERMINAL, IMSL_NON_POSITIVE_ARGUMENT);
	} else if (x <= F_TWO) {
		if (x >= xmin) {
			y = F_ZERO;
			if (x > xsml)
				y = x * x;
			dbsk1_v = log(F_HALF * x) * imsl_f_bessel_I1(x) + (0.75e0 + imsl_csevl(F_HALF *
					       y - F_ONE, lv_bk1cs, ntk1)) / x;
		} else {
                        /* The function overflows because X = %(r1) is
                          less than %(R2). */
			imsl_e1str(1, x);
			imsl_e1str(2, xmin);
			imsl_ermes(IMSL_FATAL, IMSL_SMALL_ARG_OVERFLOW);
		}
	} else {
		if (x <= xmax) {
			dbsk1_v = exp(-x) * l_dbsk1e(x);
		} else {
			dbsk1_v = F_ZERO;
                            /* The function underflows because X = %(r1) is
                               greater than %(r2).  The result is set to zero. */
			imsl_e1str(1, x);
			imsl_e1str(2, xmax);
			imsl_ermes(IMSL_ALERT, IMSL_LARGE_ARG_UNDERFLOW);
		}
	}
	E1POP("imsl_f_bessel_K1", "imsl_d_bessel_K1");
	return (dbsk1_v);
}				/* end of function */



/* Structured by FOR_STRUCT, v0.2, on 08/20/90 at 16:14:37
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  BSK1E/DBSK1E (Single/float precision version)

    Computer:   FORC/float

    Revised:    November 29, 1988

    Purpose:    Evaluate the exponentially scaled modified Bessel
                function of the third kind of order one.

    Usage:      BSK1E(X)

    Arguments:
       X      - Argument for which the function value is desired.
                (Input)
       BSK1E  - Function value.  (Output)

    GAMS:       C10b1

    Chapter:    SFUN/LIBRARY Bessel Functions

    Copyright:  1984 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static Mfloat l_dbsk1e(Mfloat x)
#else
static Mfloat l_dbsk1e(x)
	Mfloat          x;
#endif
{
	Mfloat          eta;
	Mfloat          dbsk1e_v, y;
	static Mint     ntk1 = 0;
	static Mint     ntak1 = 0;
	static Mint     ntak12 = 0;
	static Mfloat   xmin = 0.0e0;
	static Mfloat   xsml = 0.0e0;

	imsl_e1psh("l_dbsk1e");
	dbsk1e_v = imsl_amach(6);

	if (ntk1 == 0) {
		eta = 0.1 * (Mfloat) (imsl_amach(3));
		ntk1 = imsl_inits(lv_bk1cs, 16, eta);
		ntak1 = imsl_inits(lv_ak1cs, 38, eta);
		ntak12 = imsl_inits(lv_ak12cs, 33, eta);
		xmin = exp(imsl_f_max(log(imsl_amach(1)), -log(imsl_amach(2))) +
				0.01e0);
		xsml = sqrt(F_FOUR * imsl_amach(3));
	}
	if (x <= F_ZERO) {
                /* X = %(r1) must be greater than zero. */
		imsl_e1str(1, x);
		imsl_ermes(IMSL_TERMINAL, IMSL_NON_POSITIVE_ARGUMENT);
	} else if (x <= F_TWO) {
		if (x >= xmin) {
			y = F_ZERO;
			if (x > xsml)
				y = x * x;
			dbsk1e_v = exp(x) * (log(F_HALF * x) * imsl_f_bessel_I1(x) + (0.75e0 +
			   imsl_csevl(F_HALF * y - F_ONE, lv_bk1cs, ntk1)) / x);
		} else {
                        /* The function overflows because X = %(r1) is
                          less than %(R2). */
			imsl_e1str(1, x);
			imsl_e1str(2, xmin);
			imsl_ermes(IMSL_FATAL, IMSL_SMALL_ARG_OVERFLOW);
		}
	} else if (x <= F_EIGHT) {
		dbsk1e_v = (1.25e0 + imsl_csevl((16.0e0 / x - F_FIVE) / F_THREE, lv_ak1cs,
						 ntak1)) / sqrt(x);
	} else {
		dbsk1e_v = (1.25e0 + imsl_csevl(16.0e0 / x - F_ONE, lv_ak12cs, ntak12)) /
			sqrt(x);
	}
	imsl_e1pop("l_dbsk1e");
	return (dbsk1e_v);
}				/* end of function */
