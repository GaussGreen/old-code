#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif

#define COUNT(OS)   (sizeof(OS)/sizeoffloat)

static Mfloat PROTO(l_dbsi1e,(Mfloat));

     	/*
	 * SERIES FOR BI1 ON THE INTERVAL 0.0  TO  9.00000D+00 WITH WEIGHTED
	 * ERROR        1.44D-32 LOG WEIGHTED ERROR        31.84 SIGNIFICANT
	 * FIGURES REQD. 31.45 DECIMAL PLACES REQUIRED   32.46
	 */
static Mfloat lv_bi1cs[] = {
     -.197171326109985973161385032182e-2,
     .40734887667546480608155393652e0,
     .348389942999594558662450377838e-1,
     .154539455630012360385984010585e-2,
     .418885210983777841294588320041e-4,
     .764902676483621147419597039661e-6,
     .100424939247411786891798080372e-7,
     .993220779192381064813712980549e-10,
     .766380179184476372752001716814e-12,
     .474141892381673949803880919482e-14,
     .24041144040745181799863172032e-16,
     .101715050070937136491211008e-18,
     .364509356578669494584917333333e-21,
     .112057495025620393448106666667e-23,
     .29875441934468088832e-26,
     .697323109391947093333333333333e-29,
     .143679482206208e-31,
};
     	/*
	 * SERIES FOR AI1 ON THE INTERVAL 1.25000D-01 TO  3.33333D-01 WITH
	 * WEIGHTED ERROR        2.81D-32 LOG WEIGHTED ERROR        31.55
	 * SIGNIFICANT FIGURES REQD. 29.93 DECIMAL PLACES REQUIRED   32.38
	 */
static Mfloat lv_ai1cs[] = {
     -.284674418188147867410037246831e-1,
     -.192295323144322065104444877498e-1,
     -.611518585794378898225624991779e-3,
     -.206997125335022770888282377798e-4,
     .858561914581072556553694467314e-5,
     .104949824671159086251745399786e-5,
     -.29183389184479022020934323267e-6,
     -.155937814663173900016068096908e-7,
     .131801236714494470552530287391e-7,
     -.144842341818307831763913446782e-8,
     -.290851224399314209482504099301e-9,
     .12663889178753823873111596904e-9,
     -.166494777291922067062417839858e-10,
     -.1666653644609432976095937155e-11,
     .124260241429076826523216847202e-11,
     -.273154937967243239725146142863e-12,
     .202394788164580378070026268898e-13,
     .730795001811688363619869812612e-14,
     -.333290563440467494381377861713e-14,
     .717534655851295374354225466567e-15,
     -.698253032479625635585062922366e-16,
     -.129994420156276076006044608059e-16,
     .812094286424279889205467834286e-17,
     -.219401620741073689815626664378e-17,
     .363051617002965484827986093233e-18,
     -.16951397724391041663068667904e-19,
     -.128818482989790780711688253822e-19,
     .569442860496705278010999107311e-20,
     -.145959700909048005654550990029e-20,
     .251454601067571731408469133449e-21,
     -.184475888313912481816040002901e-22,
     -.6339760596227948641928609792e-23,
     .346144110203101111110814662656e-23,
     -.101706233537139354759654102357e-23,
     .214987714709043144596250077867e-24,
     -.304525242523867640174620617387e-25,
     .523808214472128598217763498667e-27,
     .1443583107089382446416789504e-26,
     -.612130207489004273320067072e-27,
     .170001111746781841834918980267e-27,
     -.359658910798424415853521578667e-28,
     .544817857894841857665051306667e-29,
     -.273183178968908498916256426667e-30,
     -.1858905021708600715771904e-30,
     .921268297451393344112776533333e-31,
     -.281383515565356110637083306667e-31,
};
     	/*
	 * SERIES FOR AI12 ON THE INTERVAL 0.0 TO  1.25000D-01 WITH WEIGHTED
	 * ERROR        1.83D-32 LOG WEIGHTED ERROR        31.74 SIGNIFICANT
	 * FIGURES REQD. 29.97 DECIMAL PLACES REQUIRED   32.66
	 */
static Mfloat lv_ai12cs[] = {
     .285762350182801204744984594847e-1,
     -.97610974913614684077651644573e-2,
     -.110588938762623716291256921278e-3,
     -.388256480887769039345654477627e-5,
     -.251223623787020892529452002212e-6,
     -.263146884688951950683705236523e-7,
     -.383538038596423702204500678797e-8,
     -.558974346219658380686811252223e-9,
     -.189749581235054123449892503324e-10,
     .325260358301548823855508067995e-10,
     .141258074366137813316336633285e-10,
     .203562854414708950722452613684e-11,
     -.719855177624590851209258989045e-12,
     -.408355111109219731822849963969e-12,
     -.210154184277266431301984572746e-13,
     .4272440016711951354297788337e-13,
     .104202769841288027641741449995e-13,
     -.38144030724370078047670725354e-14,
     -.188035477551078244851273453396e-14,
     .330820231092092828273190335241e-15,
     .296262899764595013906854654205e-15,
     -.320952592199342395877837353289e-16,
     -.465030536848935832557128281898e-16,
     .441434832307170794994611375964e-17,
     .75172963108421048054254580803e-17,
     -.931417886732688337568484784516e-18,
     -.12421932751948909561167844887e-17,
     .241427671945484846900515390218e-18,
     .202694438405328517897192286069e-18,
     -.639426718826909778704391988681e-19,
     -.304981245237309589608488450357e-19,
     .161284185165148022513462230769e-19,
     .356091396430992505451027090462e-20,
     -.375201794793643907966682800325e-20,
     -.578703742707479934595198231074e-22,
     .775999751164816196198236963209e-21,
     -.145279089720223339406445987409e-21,
     -.131822528673903670212192275337e-21,
     .611665486290307070187999133172e-22,
     .137627976242712642773024338363e-22,
     -.169083768995934788491983938231e-22,
     .143059608859543315398720108539e-23,
     .34095578280905940204053677299e-23,
     -.130945766627076022784573872642e-23,
     -.394070641124025743609352141756e-24,
     .427713742698087658080616679735e-24,
     -.442463483098260688190028312303e-25,
     -.873411319623071497211530978875e-25,
     .404540133568353339214340414243e-25,
     .706710065809468946565160771781e-26,
     -.124946334456510522300286451861e-25,
     .286739224440343703297948339143e-26,
     .204429289250429267028177957421e-26,
     -.151863663382046256837134680291e-26,
     .811018109818757588613227910704e-28,
     .358037935477358609112717370327e-27,
     -.169292901892790250959305717545e-27,
     -.222290249970242763906775852777e-28,
     .542453512714596965504860040113e-28,
     -.17870684015780186887649129933e-28,
     -.656547906872281493882392943788e-29,
     .780701316506114528092206770684e-29,
     -.181659526066897971737933315222e-29,
     -.128770495266008482037687559896e-29,
     .111454817298816454741370927369e-29,
     -.180834314503933693915936887669e-30,
     -.223167771820377195223244822894e-30,
     .161902959608034151061790980361e-30,
     -.183407990880494141390130843921e-31,
};
     

/* Structured by FOR_STRUCT, v0.2, on 08/19/90 at 12:12:45
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  BSI1/DBSI1 (Single/float precision version)

    Computer:   FORC/float

    Revised:    November 29, 1988

    Purpose:    Evaluate the modified Bessel function of the first
                kind of order one.

    Usage:      BSI1(X)

    Arguments:
       X      - Argument for which the function value is desired.
                (Input)
       BSI1   - Function value.  (Output)

    Remark:
       Informational error:
       Type Code
         2   1  The function underflows because the absolute value
                of X is too small.

    GAMS:       C10b1

    Chapter:    SFUN/LIBRARY Bessel Functions

    Copyright:  1984 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
Mfloat imsl_f_bessel_I1(Mfloat x)
#else
Mfloat imsl_f_bessel_I1(x)
	Mfloat          x;
#endif
{
	Mfloat          dbsi1_v, y;
	static Mint      nti1 = 0;
	static Mfloat   xmin = 0.0e0;
	static Mfloat   xsml = 0.0e0;
	static Mfloat   xmax = 0.0e0;

	E1PSH("imsl_f_bessel_I1", "imsl_d_bessel_I1");
	dbsi1_v = imsl_amach(6);

	if (nti1 == 0) {
		nti1 = imsl_inits(lv_bi1cs, 17, 0.1*imsl_amach(3));
		xmin = F_TWO * imsl_amach(1);
		xsml = sqrt(F_EIGHT * imsl_amach(3));
		xmax = log(imsl_amach(2));
	}
	y = fabs(x);

	if (y <= F_THREE) {
		if (y == F_ZERO) {
			dbsi1_v = F_ZERO;
		} else if (y > xsml) {
			dbsi1_v = x * (.875e0 + imsl_csevl(y * y / 4.5e0 - F_ONE, lv_bi1cs,
							    nti1));
		} else if (y > xmin) {
			dbsi1_v = F_HALF * x;
		} else {
                        /* The function underflows because ABS(%(d1)) is less
                           than %(d2). The result is set to zero. */
			imsl_e1str(1, x);
			imsl_e1str(2, xmin);
			imsl_ermes(IMSL_ALERT, IMSL_SMALL_ABS_ARG_UNDERFLOW);
			dbsi1_v = F_ZERO;
		}
	} else {
		if (y <= xmax) {
			dbsi1_v = exp(y) * l_dbsi1e(x);
		} else {
			imsl_e1str(1, x);
			imsl_e1str(2, xmax);
	                imsl_ermes(IMSL_FATAL, IMSL_LARGE_ABS_ARG_FATAL);
		}
	}
	E1POP("imsl_f_bessel_I1", "imsl_d_bessel_I1");
	return (dbsi1_v);
}				/* end of function */


/* Structured by FOR_STRUCT, v0.2, on 08/20/90 at 10:00:30
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  BSI1E/DBSI1E (Single/float precision version)

    Computer:   FORC/float

    Revised:    November 29, 1988

    Purpose:    Evaluate the exponentially scaled modified Bessel
                function of the first kind of order one.

    Usage:      BSI1E(X)

    Arguments:
       X      - Argument for which the function value is desired.
                (Input)
       BSI1E  - Function value.  (Output)

    Remark:
       Informational error:
       Type Code
         2   1  The function underflows because the absolute value
                of X is too small.

    GAMS:       C10b1

    Chapter:    SFUN/LIBRARY Bessel Functions

    Copyright:  1984 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static Mfloat l_dbsi1e(Mfloat x)
#else
static Mfloat l_dbsi1e(x)
	Mfloat          x;
#endif
{
	Mfloat          eta;
	Mfloat          dbsi1e_v, y;
	static Mint     nti1 = 0;
	static Mint     ntai1 = 0;
	static Mint     ntai12 = 0;
	static Mfloat   xmin = 0.0e0;
	static Mfloat   xsml = 0.0e0;

	imsl_e1psh("l_dbsi1e");
	dbsi1e_v = imsl_amach(6);

	if (nti1 == 0) {
		eta = 0.1 * imsl_amach(3);
		nti1 = imsl_inits(lv_bi1cs, 17, eta);
		ntai1 = imsl_inits(lv_ai1cs, 46, eta);
		ntai12 = imsl_inits(lv_ai12cs, 69,  eta);
		xmin = F_TWO * imsl_amach(1);
		xsml = sqrt(F_EIGHT * imsl_amach(3));
	}
	y = fabs(x);

	if (y <= F_THREE) {
		if (y == F_ZERO) {
			dbsi1e_v = F_ZERO;
		} else if (y > xsml) {
			dbsi1e_v = x * (.875e0 + imsl_csevl(y * y / 4.5e0 - F_ONE, lv_bi1cs,
							     nti1));
			dbsi1e_v *= exp(-y);
		} else if (y > xmin) {
			dbsi1e_v = F_HALF * x;
			dbsi1e_v *= exp(-y);
		} else {
                        /* The function underflows because ABS(%(d1)) is less
                           than %(d2). The result is set to zero. */
			imsl_e1str(1, x);
			imsl_e1str(2, xmin);
			imsl_ermes(IMSL_ALERT, IMSL_SMALL_ABS_ARG_UNDERFLOW);
		}
	} else {
		if (y <= F_EIGHT) {
			dbsi1e_v = (.375e0 + imsl_csevl((48.0e0 / y - 11.0e0) / F_FIVE,
					      lv_ai1cs, ntai1)) / sqrt(y);
		} else {
			dbsi1e_v = (.375e0 + imsl_csevl(16.0e0 / y - F_ONE, lv_ai12cs,
						    ntai12)) / sqrt(y);
		}
		dbsi1e_v = sign(dbsi1e_v, x);
	}

	imsl_e1pop("l_dbsi1e");
	return (dbsi1e_v);
}				/* end of function */
