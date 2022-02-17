#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif

static Mfloat PROTO(l_erfce,(Mfloat *x));

#define COUNT(OS)   (sizeof(OS)/sizeof(Mfloat))

	/*
	 * SERIES FOR ERFI ON THE INTERVAL 0.00000E-01  TO  5.62500E-01 WITH
	 * WEIGHTED ERROR        6.81E-17 LOG WEIGHTED ERROR        16.17
	 * SIGNIFICANT FIGURES REQD. 15.17 DECIMAL PLACES REQUIRED   16.84
	 */
static Mfloat lv_erfics[] = {
     -.0510281607605825074013733865176e0,
     .0976712136542976168051122672031e0,
     .0106798426748877000197021227858e0,
     .00145839261934235994760845092566e0,
     .000221395617935346452980889958716e0,
     .0000357546426953869466455481344805e0,
     .00000601305243584652659027746712369e0,
     .00000104056052819195311829805188061e0,
     .00000018393099988652242913597575563e0,
     .00000003304859580256689297257310837e0,
     .00000000601592296492422089765598882e0,
     .00000000110675054607426935593741602e0,
     .00000000020540500062873634187190224e0,
     .00000000003840485495835378474193636e0,
     .00000000000722613514301305057424678e0,
     .00000000000136709063800064215760026e0,
     .00000000000025987080985380550245075e0,
     .00000000000004960649711523937004868e0,
     .00000000000000950461236358905433088e0,
     .0000000000000018271407438754928837e0,
     .00000000000000035229350975383006662e0,
     .00000000000000006810937482440024885e0,
     .00000000000000001319990081950303068e0,
     .00000000000000000256390641344054319e0,
     .00000000000000000049902175821710612e0,
     .000000000000000000097308471941199e0,
     .00000000000000000001900781230046419e0,
     .00000000000000000000371884031926145e0,
     .00000000000000000000072866273739281e0,
     .00000000000000000000014296976677763e0,
     .00000000000000000000002808794800238e0,
     .0000000000000000000000055248115833e0,
     .00000000000000000000000108793649641e0,
     .00000000000000000000000021446095381e0,
     .0000000000000000000000000423179096e0,
     .00000000000000000000000000835807968e0,
     .00000000000000000000000000165223675e0,
     .00000000000000000000000000032688917e0,
     .00000000000000000000000000006472499e0,
     .00000000000000000000000000001282534e0,
     .00000000000000000000000000000254316e0,
     .00000000000000000000000000000050463e0,
     .0000000000000000000000000000001002e0,
     .00000000000000000000000000000001991e0
};

/* Structured by FOR_STRUCT, v0.2, on 06/05/90 at 17:20:08
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  ERFCI/DERFCI (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    November 28, 1988

    Purpose:    Evaluate the inverse complementary error function.

    Usage:      ERFCI(X)

    Arguments:
       X      - Argument for which the function value is desired.
                (Input)
       ERFCI  - Function value.  (Output)

    Remark:
       Informational error:
       Type Code
         3   2  Result of ERFCI(X) is accurate to less than one half
                precision because the argument is too close to 2.0.

    GAMS:       C8a

    Chapter:    SFUN/LIBRARY Error Function and Related Functions

    Copyright:  1984 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
Mfloat imsl_f_erfc_inverse(Mfloat x)
#else
Mfloat imsl_f_erfc_inverse(x)
	Mfloat           x;
#endif
{
	Mfloat           alneta, alnz, corr, erfc_inverse_v, eta, etas, y, yabs, z;
	static Mint      nerfi = 0;
	static Mfloat    xbig = 0.0e0;
	static Mfloat    eps = 0.0e0;
	static Mfloat    sqrtpi = 1.77245385090551602729816748334e0;

	E1PSH("imsl_f_erfc_inverse", "imsl_d_erfc_inverse");

	erfc_inverse_v = imsl_amach(6);

	if (nerfi == 0) {
		etas = 0.1 * imsl_amach(3);
		nerfi = imsl_inits(lv_erfics, COUNT(lv_erfics), etas);
		eps = sqrt(imsl_amach(4) / (F_FOUR * 3.1415));
		xbig = -eps / log(1.77 * eps);
		xbig += -(xbig * log(1.77 * xbig) + eps) / (F_ONE + log(1.77 * xbig));
		xbig = F_TWO - xbig;
		eps = F_TEN * imsl_amach(4);
	}
	if (x <= F_ZERO || x >= F_TWO ) {
                    /* %(L1) must be between %(R1) and $(R2),
                       but %(L1) = %(R3). */
                imsl_e1stl(1, "x");
                imsl_e1str(1, F_ZERO);
                imsl_e1str(2, F_TWO);
		imsl_e1str(3, x);
		imsl_ermes(IMSL_FATAL, IMSL_REAL_OUT_OF_RANGE);
		goto L_9000;
	} else if (x > xbig) {
                    /* The result is accurate to less than one half
                       precision because %(R1) is greater than %(R2). */
		imsl_e1str(1, x);
                imsl_e1str(2, xbig);
		imsl_ermes(IMSL_WARNING, IMSL_LARGE_ARG_WARN);
	}
	y = F_ONE - x;
	yabs = fabs(y);

	if (yabs >= 0.75) {
                Mint    iter;
		if (y >= F_ZERO) {
			z = x;
		} else {
			z = F_ONE + (F_ONE - x);
		}
		alnz = log(z);
		eta = -log(sqrtpi * z);
		alneta = log(eta);
		erfc_inverse_v = sqrt(eta - F_HALF * alneta + (0.25 * alneta - F_HALF) / eta);

		for (iter = 1; iter <= 100; iter++) {
			corr = F_HALF * sqrtpi * (l_erfce(&erfc_inverse_v) - exp(alnz + erfc_inverse_v *
								  erfc_inverse_v));
			erfc_inverse_v += corr;
			if (fabs(corr) < eps * erfc_inverse_v)
				goto L_20;
		}
                    /* Sufficiently accurate result was not achieved.
                       An algorithmic problem exists. */
		imsl_ermes(IMSL_FATAL, IMSL_ERF_ALGORITHM);
		goto L_9000;
L_20:
		if (y < F_ZERO)
			erfc_inverse_v = -erfc_inverse_v;
	} else {
		erfc_inverse_v = y * (F_ONE + imsl_csevl(32.0 * y * y / F_NINE - F_ONE, lv_erfics, nerfi));
	}

L_9000:
	E1POP("imsl_f_erfc_inverse", "imsl_d_erfc_inverse");
	return (erfc_inverse_v);
}


	/*
	 * SERIES FOR ERF ON THE INTERVAL 0.0  TO  1.00000E+00 WITH WEIGHTED
	 * ERROR        7.10E-18 LOG WEIGHTED ERROR        17.15 SIGNIFICANT
	 * FIGURES REQD. 16.31 DECIMAL PLACES REQUIRED   17.71
	 */
static Mfloat lv_erfcs[] = {
     -.490461212346918080399845440334e-1,
     -.142261205103713642378247418996e0,
     .100355821875997955757546767129e-1,
     -.576876469976748476508270255092e-3,
     .274199312521960610344221607915e-4,
     -.110431755073445076041353812959e-5,
     .384887554203450369499613114982e-7,
     -.118085825338754669696317518016e-8,
     .323342158260509096464029309534e-10,
     -.799101594700454875816073747086e-12,
     .179907251139614556119672454866e-13,
     -.371863548781869263823168282095e-15,
     .710359900371425297116899083947e-17,
     -.126124551191552258324954248533e-18,
     .209164069417692943691705002667e-20,
     -.3253973102931407298236416e-22,
     .476686720979767483323733333333e-24,
     -.659801207828513431552e-26,
     .865501146996376261973333333333e-28,
     -.107889251774980642133333333333e-29,
     .128118839930170026666666666667e-31,
};
	/*
	 * SERIES FOR ERC2 ON THE INTERVAL 2.50000E-01  TO  1.00000E+00 WITH
	 * WEIGHTED ERROR        5.22E-17 LOG WEIGHTED ERROR        16.28
	 * APPROX SIGNIFICANT FIGURES REQUIRED                  15.0 DECIMAL
	 * PLACES REQUIRED   16.96
	 */
static Mfloat lv_erc2cs[] = {
     -.69601346602309501127391508262e-1,
     -.411013393626208934898221208467e-1,
     .391449586668962688156114370524e-2,
     -.490639565054897916128093545077e-3,
     .715747900137703638076089414183e-4,
     -.115307163413123283380823284791e-4,
     .199467059020199763505231486771e-5,
     -.364266647159922287393611843071e-6,
     .694437261000501258993127721463e-7,
     -.137122090210436601953460514121e-7,
     .278838966100713713196386034809e-8,
     -.581416472433116155186479105032e-9,
     .123892049175275318118016881795e-9,
     -.269063914530674343239042493789e-10,
     .594261435084791098244470968384e-11,
     -.133238673575811957928775442057e-11,
     .30280468061771320171736972433e-12,
     -.696664881494103258879586758895e-13,
     .162085454105392296981289322763e-13,
     -.380993446525049199987691305773e-14,
     .904048781597883114936897101298e-15,
     -.2164006195089607347809812047e-15,
     .522210223399585498460798024417e-16,
     -.126972960236455533637241552778e-16,
     .310914550427619758383622741295e-17,
     -.766376292032038552400956671481e-18,
     .190081925136274520253692973329e-18,
     -.474220727906903954522565599997e-19,
     .118964920007652838288068307845e-19,
     -.300003559032578025684527131307e-20,
     .76029934530432461730193852771e-21,
     -.193590944760687288156981104913e-21,
     .495139912477333788100004238677e-22,
     -.127180748133637187960862198989e-22,
     .328004960046951304331584165205e-23,
     -.84923201768228965689247924224e-24,
     .22069178928075602235198799872e-24,
     -.57556172456965284983128195072e-25,
     .15061915336392342503541440512e-25,
     -.3954502959018796953104285696e-26,
     .104152970415150097998464505173e-26,
     -.275148779527876507945017890133e-27,
     .729005820549755740899770368e-28,
     -.193693964591594780407750109867e-28,
     .516035711205148729837005482667e-29,
     -.13784193221930940993896448e-29,
     .369132679310706904225109333333e-30,
     -.990938959062436542065322666667e-31,
     .266649170519538841332394666667e-31,
};
	/*
	 * SERIES FOR ERFC ON THE INTERVAL 0.0  TO  2.50000E-01 WITH WEIGHTED
	 * ERROR        4.81E-17 LOG WEIGHTED ERROR        16.32 APPROX
	 * SIGNIFICANT FIGURES REQUIRED                  15.0 DECIMAL PLACES
	 * REQUIRED   17.01
	 */
static Mfloat lv_erfccs[] = {
     .715179310202924774503697709496e-1,
     -.265324343376067157558893386681e-1,
     .171115397792085588332699194606e-2,
     -.163751663458517884163746404749e-3,
     .198712935005520364995974806758e-4,
     -.284371241276655508750175183152e-5,
     .460616130896313036969379968464e-6,
     -.822775302587920842057766536366e-7,
     .159214187277090112989358340826e-7,
     -.329507136225284321486631665072e-8,
     .72234397604005554658126115389e-9,
     -.166485581339872959344695966886e-9,
     .401039258823766482077671768814e-10,
     -.100481621442573113272170176283e-10,
     .260827591330033380859341009439e-11,
     -.699111056040402486557697812476e-12,
     .192949233326170708624205749803e-12,
     -.547013118875433106490125085271e-13,
     .158966330976269744839084032762e-13,
     -.47268939801975548392036958429e-14,
     .14358733767849847867287399784e-14,
     -.444951056181735839417250062829e-15,
     .140481088476823343737305537466e-15,
     -.451381838776421089625963281623e-16,
     .147452154104513307787018713262e-16,
     -.489262140694577615436841552532e-17,
     .164761214141064673895301522827e-17,
     -.562681717632940809299928521323e-18,
     .194744338223207851429197867821e-18,
     -.682630564294842072956664144723e-19,
     .242198888729864924018301125438e-19,
     -.869341413350307042563800861857e-20,
     .315518034622808557122363401262e-20,
     -.115737232404960874261239486742e-20,
     .428894716160565394623737097442e-21,
     -.160503074205761685005737770964e-21,
     .606329875745380264495069923027e-22,
     -.231140425169795849098840801367e-22,
     .888877854066188552554702955697e-23,
     -.344726057665137652230718495566e-23,
     .134786546020696506827582774181e-23,
     -.531179407112502173645873201807e-24,
     .210934105861978316828954734537e-24,
     -.843836558792378911598133256738e-25,
     .339998252494520890627359576337e-25,
     -.13794523880732420900223837711e-25,
     .563449031183325261513392634811e-26,
     -.2316490434477065448234277527e-26,
     .958446284460181015263158381226e-27,
     -.399072288033010972624224850193e-27,
     .167212922594447736017228709669e-27,
     -.704599152276601385638803782587e-28,
     .297976840286420635412357989444e-28,
     -.126252246646061929722422632994e-28,
     .539543870454248793985299653154e-29,
     -.238099288253145918675346190062e-29,
     .10990528301027615735972668375e-29,
     -.486771374164496572732518677435e-30,
     .152587726411035756763200828211e-30
};

/* Structured by FOR_STRUCT, v0.2, on 06/05/90 at 17:18:14
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  ERFCE/DERFCE (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    November 28, 1988

    Purpose:    Evaluate the exponentially scaled complementary
                error function.

    Usage:      ERFCE(X)

    Arguments:
       X      - Argument for which the function value is desired.
                (Input)
       ERFCE  - Function value.  (Output)

    Remark:
       Informational error:,
       Type Code
         2   1  The function underflows because X is too large.

    GAMS:       C8a

    Chapter:    SFUN/LIBRARY Error Function and Related Functions

    Copyright:  1984 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static Mfloat l_erfce(Mfloat *x)
#else
static Mfloat l_erfce(x)
	Mfloat          *x;
#endif
{
	Mfloat           erfce_v, eta, y;
	static Mfloat    sqrtpi = 1.77245385090551602729816748334e0;
	static Mint      nterf = 0;
	static Mint      nterfc = 0;
	static Mint      nterc2 = 0;
	static Mfloat    xmin = 0.0;
	static Mfloat    xsml = 0.0;
	static Mfloat    xmax = 0.0;
	static Mfloat    sqeps = 0.0;

	imsl_e1psh("l_erfce");
	erfce_v = imsl_amach(6);

	if (nterf == 0) {
		eta = 0.1 * imsl_amach(3);
		nterf = imsl_inits(lv_erfcs, COUNT(lv_erfcs), eta);
		nterfc = imsl_inits(lv_erfccs, COUNT(lv_erfccs), eta);
		nterc2 = imsl_inits(lv_erc2cs, COUNT(lv_erc2cs), eta);
		xmin = -sqrt(log(imsl_amach(2) / F_TWO)) + 0.01;
		xsml = -sqrt(-log(sqrtpi * imsl_amach(3)));
		xmax = exp(imsl_f_min(-log(sqrtpi * imsl_amach(1)), log(imsl_amach(2))) -
			   .01);
		sqeps = sqrt(F_TWO * imsl_amach(3));
	}
	if (*x < xmin) {
                    /* The function overflows because X = %(r1) is
                       less than %(R2). */
		imsl_e1str(1, *x);
		imsl_e1str(2, xmin);
		imsl_ermes(IMSL_FATAL, IMSL_SMALL_ARG_OVERFLOW);
	} else if (*x <= xsml) {
		/*
		 * ERFC(X) = 1.0 - ERF(X) FOR X .LT. XSML
		 */
		erfce_v = F_TWO * exp(*x ** x);
	} else if (*x <= xmax) {
		y = fabs(*x);
		/*
		 * ERFC(X) = 1.0 - ERF(X) FOR -1.0 .LE. X .LE. 1.0
		 */
		if (y <= F_ONE) {
			if (y < sqeps) {
				erfce_v = F_ONE - F_TWO ** x / sqrtpi;
			} else {
				erfce_v = exp(*x ** x) * (F_ONE - *x * (F_ONE + imsl_csevl(F_TWO *
					     *x ** x - F_ONE, lv_erfcs, nterf)));
			}
			/*
			 * ERFC(X) = 1.0 - ERF(X) FOR 1.0 .LT. ABS(X) .LE.
			 * XMAX
			 */
		} else {
			y *= y;
			if (y <= F_FOUR) {
				erfce_v = (F_HALF + imsl_csevl((F_EIGHT / y - F_FIVE) / F_THREE, lv_erc2cs,
							 nterc2)) / fabs(*x);
			} else {
				erfce_v = (F_HALF + imsl_csevl(F_EIGHT / y - F_ONE, lv_erfccs, nterfc)) /
					fabs(*x);
			}
			if (*x < F_ZERO)
				erfce_v = F_TWO * exp(y) - erfce_v;
		}
	} else {
                    /* The function underflows because X = %(r1) is
                       greater than %(r2).  The result is set to zero. */
		imsl_e1str(1, *x);
		imsl_e1str(2, xmax);
		imsl_ermes(IMSL_ALERT, IMSL_LARGE_ARG_UNDERFLOW);
		erfce_v = F_ZERO;
	}
	imsl_e1pop("l_erfce");
	return (erfce_v);
}				/* end of function */
