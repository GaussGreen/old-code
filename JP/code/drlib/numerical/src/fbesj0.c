#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif

#define COUNT(OS)   (sizeof(OS)/sizeof(Mfloat))
static void PROTO(l_d9b0mp,(Mfloat *x, Mfloat *ampl, Mfloat *theta));

	/*
	 * SERIES FOR BJ0 ON THE INTERVAL 0.0  TO  1.60000D+01 WITH WEIGHTED
	 * ERROR        4.39D-32 LOG WEIGHTED ERROR        31.36 SIGNIFICANT
	 * FIGURES REQD. 31.21 DECIMAL PLACES REQUIRED   32.00
	 */
static Mfloat lv_bj0cs[] = {
     .100254161968939137010731272641e0,
     -.665223007764405131776787578311e0,
     .248983703498281313704604687267e0,
     -.332527231700357696538843415039e-1,
     .231141793046940154629049241177e-2,
     -.991127741995080923390485193366e-4,
     .289167086439988088847339037471e-5,
     -.612108586630326350578184074815e-7,
     .983865079385678413247687486364e-9,
     -.124235515973017651455158970068e-10,
     .126543363025590457979158272104e-12,
     -.10619456495287244546914817513e-14,
     .74706210758024567437098915584e-17,
     -.44697032274412780547627008e-19,
     .230242815843374362005230933333e-21,
     .406081782748733227008e-26,
     -.1414383600524091392e-28,
    .4391090549669888e-31
};
	/*
	 * SERIES FOR BY0 ON THE INTERVAL 0.0  TO  1.60000D+01 WITH WEIGHTED
	 * ERROR        8.14D-32 LOG WEIGHTED ERROR        31.09 SIGNIFICANT
	 * FIGURES REQD. 30.31 DECIMAL PLACES REQUIRED   31.73
	 */
static Mfloat lv_by0cs[] = {
     -.112778393928655732179398054603e-1,
     -.128345237560420346048088453184e0,
     -.104378847997942493658176227662e0,
     .236627491839696954092415926461e-1,
     -.209039164770048623919622395034e-2,
     .103975453939057252099924657638e-3,
     -.336974716242397209671877534504e-5,
     .772938426767066715852136721637e-7,
     -.132497677266425959144347606896e-8,
     .176482326154045279210038936316e-10,
     -.188105507158019620060282301207e-12,
     .164186548536614950279223718575e-14,
     -.119565943860460608574599100672e-16,
     .737729629744018584249411242667e-19,
     -.390684347671043733074090666667e-21,
     .179550366443615794982912e-23,
     -.722962712544801047893333333333e-26,
     .257172793163516859733333333333e-28,
     -.814126881416369493333333333333e-31
};
	/*
	 * SERIES FOR BM0 ON THE INTERVAL 1.56250D-02 TO  6.25000D-02 WITH
	 * WEIGHTED ERROR        4.40D-32 LOG WEIGHTED ERROR        31.36
	 * SIGNIFICANT FIGURES REQD. 30.02 DECIMAL PLACES REQUIRED   32.14
	 */
static Mfloat lv_bm0cs[] = {
     .921165624682774271257376773018e-1,
     -.105059099727190510248071637176e-2,
     .147015984076875975405639285095e-4,
     -.50585576060385542233479293277e-6,
     .278725453863244417663035613788e-7,
     -.206236361178091480261884101897e-8,
     .187021431313887967513817259626e-9,
     -.196933097113563620024173077783e-10,
     .232597379399927544401250881805e-11,
     -.300952034493825027285122473448e-12,
     .419452133385066918147120676865e-13,
     -.621944931218844582597326742956e-14,
     .971826041133606846960176588527e-15,
     -.158847858570107520736663596694e-15,
     .270007219367130889008621732446e-16,
     -.475009236523400899247750478677e-17,
     .861512816260437087319170374656e-18,
     -.160560868695614481574560270336e-18,
     .30665139873144829751885398016e-19,
     -.598776422319395643069650561707e-20,
     .119297125374824830648906984107e-20,
     -.242096914204480548948468258133e-21,
     .499675176051061645337100288e-22,
     -.1047493639351158510095040512e-22,
     .222778684379746810104818346667e-23,
     -.480181323939816286237054293333e-24,
     .104796272347095995647699626667e-24,
     -.23138581656786153251012608e-25,
     .51648230884626742116352e-26,
     -.11646911918500653895254016e-26,
     .2651788486043319282958336e-27,
     -.609255950382572849769130666667e-28,
     .141180468614425930803882666667e-28,
     -.329809496123173724575061333333e-29,
     .776393114307406503171413333333e-30,
     -.184103134366145847842133333333e-30,
     .43958801385943107371008e-31
};
	/*
	 * SERIES FOR BTH0 ON THE INTERVAL 0.0  TO  1.56250D-02 WITH WEIGHTED
	 * ERROR        2.66D-32 LOG WEIGHTED ERROR        31.57 SIGNIFICANT
	 * FIGURES REQD. 30.67 DECIMAL PLACES REQUIRED   32.40
	 */
static Mfloat lv_bth0cs[] = {
     -.2490178086212893671770979379e0,
     .485502996096237492410486155355e-3,
     -.545118373450172049506562735635e-5,
     .135586730594059640543774459299e-6,
     -.556913989022276262275832184149e-8,
     .326090318249943353040042057195e-9,
     -.24918807862461341125237903878e-10,
     .234493774208825205543524135649e-11,
     -.260965344443103877621775747661e-12,
     .333531404200973951058699550149e-13,
     -.478900004405726846467507705574e-14,
     .759561784361922159726425685453e-15,
     -.131315560168914403827733974876e-15,
     .244836183452408574954268207384e-16,
     -.488057298106187776832567619183e-17,
     .103272850297863161492237563612e-17,
     -.23057633815057217157004744527e-18,
     .540444430018926939930171084838e-19,
     -.132406951943665727241550328824e-19,
     .337807956213719702034247921247e-20,
     -.894576291571117790030269262923e-21,
     .245199068892193170908999086514e-21,
     -.693884228768663186801399331577e-22,
     .202282787148901383929463033378e-22,
     -.606285000023354831057941953718e-23,
     .186497489640376353818237883963e-23,
     -.587837323848498945602450365309e-24,
     .189585914479995634855311795035e-24,
     -.624819793722588589592916207286e-25,
     .210179016845510246866386335291e-25,
     -.720843009352092536908139339925e-26,
     .251813638924742408671564059768e-26,
     -.895180422587857788061439459536e-27,
     .323572374797622985332562358686e-27,
     -.118830105198553536570471441138e-27,
     .443062869073581048205792319417e-28,
     -.167610096488348294957920101357e-28,
     .642929469212074669725323939661e-29,
     -.249922611669786524212072136828e-29,
     .983997942995219556728282603553e-30,
     -.392203752424080163979891316262e-30,
     .158181070300565221385906188457e-30,
     -.645255061448907159443440983654e-31,
     .266111113691993561371770183464e-31
};
	/*
	 * SERIES FOR BM02 ON THE INTERVAL 0.0  TO  1.56250D-02 WITH WEIGHTED
	 * ERROR        4.72D-32 LOG WEIGHTED ERROR        31.33 SIGNIFICANT
	 * FIGURES REQD. 30.00 DECIMAL PLACES REQUIRED   32.13
	 */
static Mfloat lv_bm02cs[] = {
     .950041514522838136933086133556e-1,
     -.380186468236567099174808156685e-3,
     .225833930103148119295182992722e-5,
     -.389572580237222876473062141261e-7,
     .124688641651208169793099052973e-8,
     -.606594902210250377980383505839e-10,
     .400846165142174699101527597105e-11,
     -.335099818339809421846729879457e-12,
     .3377119716517417367063264342e-13,
     -.396458590163501270056935629582e-14,
     .528611150388385721738793974474e-15,
     -.785251908345085231365464024349e-16,
     .128030057338668220101163407345e-16,
     -.226399629639142977628709924488e-17,
     .430049692965679038864641029048e-18,
     -.870574980513258707974753545146e-19,
     .186586271396209514118144277205e-19,
     -.42104824860930654573450869723e-20,
     .995667696422840099158162741784e-21,
     -.245735744280531335960592147855e-21,
     .630769216076203156808735370706e-22,
     -.167877369144074014269333117239e-22,
     .462025906467390443377087813609e-23,
     -.13117822668603087322376934025e-23,
     .383408756411630282774792244028e-24,
     -.115145932407774127107261329358e-24,
     .354721000752333852307697134521e-25,
     -.111921838581500464626435594218e-25,
     .361187942762983783169840499426e-26,
     -.119068776591333315009264176246e-26,
     .400509405940396813180247644954e-27,
     -.137316942245221239059519391602e-27,
     .479419908874253158599649152644e-28,
     -.170296562762410958400699447645e-28,
     .614951242893633007150357516132e-29,
     -.225576689658182834994430023724e-29,
     .83997075092942994860616583532e-30,
     -.317299759556260235556742393615e-30,
     .121520529888129855458333302651e-30,
     -.471585274975443869301321056805e-31
};
	/*
	 * SERIES FOR BT02 ON THE INTERVAL 1.56250D-02 TO  6.25000D-02 WITH
	 * WEIGHTED ERROR        2.99D-32 LOG WEIGHTED ERROR        31.52
	 * SIGNIFICANT FIGURES REQD. 30.61 DECIMAL PLACES REQUIRED   32.32
	 */
static Mfloat lv_bt02cs[] = {
     -.245482952134245974620504672493e0,
     .125441210390846157807853317783e-2,
     -.312539504148715228549734467096e-4,
     .147097782499408311644534269693e-5,
     -.995434889379500336434688503512e-7,
     .854931667332030412475787113978e-8,
     -.869897595265543345579855121792e-9,
     .100520995335597910845401010822e-9,
     -.128282306017088929034836236855e-10,
     .17731700781805131705655750451e-11,
     -.261745745694855774886362841809e-12,
     .408283513899720596219664812211e-13,
     -.667516682397427200546067495543e-14,
     .1136576139307162944839246955e-14,
     -.200511896206471602505592664121e-15,
     .364979787947662696357205914641e-16,
     -.683096375645823031693558437888e-17,
     .13107583145670756620057104268e-17,
     -.257233631018506077787571306496e-18,
     .515216574418639599252677809493e-19,
     -.105130175637588026379407414613e-19,
     .218203819911948138473010845013e-20,
     -.460047012103621605772259054933e-21,
     .984070069254668185209536512e-22,
     -.213340380357283758447359863467e-22,
     .468310364239733652960662869333e-23,
     -.104002136919857472365133824e-23,
     .233491056773015100517777408e-24,
     -.529568253233186157880497493333e-25,
     .12126341952959756829196288e-25,
     -.280188970822894287602756266667e-26,
     .652926789870128733425937066667e-27,
     -.153379800618733464278357333333e-27,
     .363058843063645366823594666667e-28,
     -.865607557136291224791722666667e-29,
     .207799099725362845712384e-29,
     -.502111702214172216743253333333e-30,
     .12208360279441714184192e-30,
     -.298600562670399134542506666667e-31
};

/* Structured by FOR_STRUCT, v0.2, on 08/19/90 at 12:10:23
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  BSJ0/DBSJ0 (Single/float precision version)

    Computer:   FORC/float

    Revised:    November 29, 1988

    Purpose:    Evaluate the Bessel function of the first kind of
                order zero.

    Usage:      BSJ0(X)

    Arguments:
       X      - Argument for which the function value is desired.
                (Input)
       BSJ0   - Function value.  (Output)

    GAMS:       C10a1

    Chapter:    SFUN/LIBRARY Bessel Functions

    Copyright:  1984 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
Mfloat imsl_f_bessel_J0(Mfloat x)
#else
Mfloat imsl_f_bessel_J0(x)
	Mfloat          x;
#endif
{
	Mfloat           ampl, bsj0_v, theta, y, y_float;
	static Mint      ntj0 = 0;
	static Mfloat    xsml = 0.0e0;
	static Mfloat    xwarn = 0.0e0;

	E1PSH("imsl_f_bessel_J0", "imsl_d_bessel_J0");
	bsj0_v = imsl_amach(6);

	if (ntj0 == 0) {
		ntj0 = imsl_inits(lv_bj0cs, COUNT(lv_bj0cs), 0.1*imsl_amach(3));
		xsml = sqrt(F_FOUR * imsl_amach(3));
		xwarn = F_ONE/sqrt(imsl_amach(3));
	}
	y = fabs(x);

	if (y <= F_FOUR) {
		if (y <= xsml) {
			bsj0_v = F_ONE;
		} else {
			bsj0_v = imsl_csevl(.125e0*y*y-F_ONE, lv_bj0cs, ntj0);
		}
	} else {
		y_float = y;
		l_d9b0mp(&y_float, &ampl, &theta);
		if ((imsl_n1rty(0) == IMSL_TERMINAL)||(imsl_n1rty(0) == IMSL_FATAL))
			goto L_9000;
		bsj0_v = ampl * cos(theta);
                if (y > xwarn) {
                    /* The result is accurate to less than one half
                       precision because abs(%(R1)) is greater
                       than %(R2). */
                    imsl_e1str(1, x);
                    imsl_e1str(2, xwarn);
                    imsl_ermes(IMSL_WARNING, IMSL_LARGE_ABS_ARG_WARN);
                }
	}
L_9000:
	E1POP("imsl_f_bessel_J0", "imsl_d_bessel_J0");
	return (bsj0_v);
}				/* end of function */


/* Structured by FOR_STRUCT, v0.2, on 08/19/90 at 12:13:34
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  BSY0/DBSY0 (Single/float precision version)

    Computer:   FORC/float

    Revised:    November 29, 1988

    Purpose:    Evaluate the Bessel function of the second kind of
                order zero.

    Usage:      BSY0(X)

    Arguments:
       X      - Argument for which the function value is desired.
                (Input)
       BSY0   - Function value.  (Output)

    GAMS:       C10c

    Chapter:    SFUN/LIBRARY Bessel Functions

    Copyright:  1984 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
Mfloat imsl_f_bessel_Y0(Mfloat x)
#else
Mfloat imsl_f_bessel_Y0(x)
	Mfloat          x;
#endif
{
	Mfloat          ampl, dbsy0_v, theta, y, x_float;
	/* TWODPI = 2.0/PI */
	static Mfloat   twodpi = 0.63661977236758134307553505349e0;
	static Mfloat   alnhaf = -0.693147180559945309417232121458e0;
	static Mint     nty0 = 0;
	static Mfloat   xsml = 0.0e0;
        static Mfloat   xwarn = 0.0e0;

	E1PSH("imsl_f_bessel_Y0", "imsl_d_bessel_Y0");
	dbsy0_v = imsl_amach(6);

	if (nty0 == 0) {
		nty0 = imsl_inits(lv_by0cs, COUNT(lv_by0cs), 0.1*imsl_amach(3));
		xsml = sqrt(F_FOUR * imsl_amach(3));
                xwarn = F_ONE/sqrt(imsl_amach(3));
	}
	if (x <= F_ZERO) {
                /* X = %(r1) must be greater than zero. */
		imsl_e1str(1, x);
		imsl_ermes(IMSL_TERMINAL, IMSL_NON_POSITIVE_ARGUMENT);
	} else if (x <= F_FOUR) {
		if (x <= xsml) {
			y = F_ZERO;
		} else {
			y = x * x;
		}
		dbsy0_v = twodpi * (alnhaf + log(x)) * imsl_f_bessel_J0(x) + .375e0 +
			imsl_csevl(.125e0 * y - F_ONE, lv_by0cs, nty0);
	} else {
		x_float = x;
		l_d9b0mp(&x_float, &ampl, &theta);
		if ((imsl_n1rty(0) == IMSL_TERMINAL)||(imsl_n1rty(0) == IMSL_FATAL))
			goto L_9000;
		dbsy0_v = ampl * sin(theta);
                if (x > xwarn) {
                    /* The result is accurate to less than one half
                       precision because abs(%(R1)) is greater
                       than %(R2). */
                    imsl_e1str(1, x);
                    imsl_e1str(2, xwarn);
                    imsl_ermes(IMSL_WARNING, IMSL_LARGE_ABS_ARG_WARN);
                }
	}
L_9000:
	E1POP("imsl_f_bessel_Y0", "imsl_d_bessel_Y0");
	return (dbsy0_v);
}				/* end of function */


/* Structured by FOR_STRUCT, v0.2, on 08/19/90 at 12:50:00
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  D9B0MP

    Computer:   FORC/float

    Revised:    January 1, 1984

    Purpose:    Evaluate the modulus and phase for the J0 and Y0
                Bessel functions.

    Usage:      CALL D9B0MP (X, AMPL, THETA)

    Arguments:
       X      - float precision argument for which the modulus and
                phase are desired. (Input)
       AMPL   - float precision modulus for J0 and Y0. (Output)
       THETA  - float precision phase for J0 and Y0. (Output)

    Copyright:  1984 by IMSL, Inc. All rights reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_d9b0mp(Mfloat *x, Mfloat *ampl, Mfloat *theta)
#else
static void l_d9b0mp(x, ampl, theta)
	Mfloat         *x, *ampl, *theta;
#endif
{
	Mfloat          eta;
	Mfloat          z;
	static Mfloat   pi4 = 0.78539816339744830961566084582e0;
	static Mint     nbm0 = 0;
	static Mint     nbt02 = 0;
	static Mint     nbm02 = 0;
	static Mint     nbth0 = 0;
	static Mfloat   xmax = 0.0e0;

	imsl_e1psh("l_d9b0mp");
	*ampl = imsl_amach(6);
	*theta = imsl_amach(6);

	if (nbm0 == 0) {
		eta = 0.1 * imsl_amach(3);
		nbm0 = imsl_inits(lv_bm0cs, COUNT(lv_bm0cs),  eta);
		nbt02 = imsl_inits(lv_bt02cs, COUNT(lv_bt02cs), eta);
		nbm02 = imsl_inits(lv_bm02cs, COUNT(lv_bm02cs), eta);
		nbth0 = imsl_inits(lv_bth0cs, COUNT(lv_bth0cs), eta);
		xmax = F_ONE / imsl_amach(4);
	}
	if (*x <= F_EIGHT) {
		z = (128.0e0 / (*x ** x) - F_FIVE) / F_THREE;
		*ampl = (.75e0 + imsl_csevl(z, lv_bm0cs, nbm0)) / sqrt(*x);
		*theta = *x - pi4 + imsl_csevl(z, lv_bt02cs, nbt02) / *x;
	} else if (*x <= xmax) {
		z = 128.0e0 / (*x ** x) - F_ONE;
		*ampl = (.75e0 + imsl_csevl(z, lv_bm02cs, nbm02)) / sqrt(*x);
		*theta = *x - pi4 + imsl_csevl(z, lv_bth0cs, nbth0) / *x;
	} else {
		imsl_e1str(1, *x);
		imsl_e1str(2, xmax);
                imsl_ermes(IMSL_FATAL, IMSL_LARGE_ABS_ARG_FATAL);
	}

	imsl_e1pop("l_d9b0mp");
	return;
}				/* end of function */
