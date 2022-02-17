#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif

#define COUNT(OS)   (sizeof(OS)/sizeof(Mfloat))
static void PROTO(l_d9b1mp,(Mfloat *x, Mfloat *ampl, Mfloat *theta));
	/*
	 * SERIES FOR BJ1 ON THE INTERVAL 0.0  TO  1.60000e+01 WITH WEIGHTED
	 * ERROR        1.16e-33 LOG WEIGHTED ERROR        32.93 SIGNIFICANT
	 * FIGURES REQD. 32.36 DECIMAL PLACES REQUIRED   33.57
	 */
static Mfloat lv_bj1cs[19] = {
	-.117261415133327865606240574524e+0,
	-.253615218307906395623030884555e+0,
	.501270809844695685053656363204e-1,
	-.46315148096250819184261972879e-2,
	.247996229415914024539124064592e-3,
	-.867894868627882584521246435176e-5,
	.214293917143793691502766250991e-6,
	-.393609307918317979229322764073e-8,
	.559118231794688004018248059864e-10,
	-.632761640466139302477695274015e-12,
	.584099161085724700326945563268e-14,
	-.4482533818701258190391350592e-16,
	.290538449262502466306018688e-18,
	-.161173219784144165412118186667e-20,
	.773947881939274637298346666667e-23,
	-.324869378211199841143466666667e-25,
	.12022376772274102272e-27,
	-.395201221265134933333333333333e-30,
	.116167808226645333333333333333e-32
};
	/*
	 * SERIES FOR BY1 ON THE INTERVAL 0.0  TO  1.60000e+01 WITH WEIGHTED
	 * ERROR        8.65e-33 LOG WEIGHTED ERROR        32.06 SIGNIFICANT
	 * FIGURES REQD. 32.17 DECIMAL PLACES REQUIRED   32.71
	 */
static Mfloat lv_by1cs[] = {
    .320804710061190862932352018628e-1,
     .126270789743350044953431726e1,
     .649996189992317500097490637314e-2,
     -.89361645288605041165314416001e-1,
     .13250881221757095451237551037e-1,
     -.897905911964835237753039508298e-3,
     .364736148795830678242287368165e-4,
     -.100137438166600055549075523845e-5,
     .199453965739017397031159372421e-7,
     -.302306560180338167284799332521e-9,
     .360987815694781196116252914243e-11,
     -.348748829728758242414552947409e-13,
     .278387897155917665813507698517e-15,
     -.186787096861948768766825352533e-17,
     .106853153391168259757070336e-19,
     -.527472195668448228943872e-22,
     .227019940315566414370133333333e-24,
     -.859539035394523108693333333333e-27,
     .288540437983379456e-29,
     -.864754113893717333333333333333e-32
};
	/*
	 * SERIES FOR BM1 ON THE INTERVAL 1.56250e-02 TO 6.25000e-02 WITH
	 * WEIGHTED ERROR        4.91e-32 LOG WEIGHTED ERROR        31.31
	 * SIGNIFICANT FIGURES REQD. 30.04 DECIMAL PLACES REQUIRED   32.09
	 */
static Mfloat lv_bm1cs[37] = {
	.106984545261806301496998530854e+0,
	.327491503971596490072905514345e-2,
	-.298778326683169859203044577794e-4,
	.833123717799197453139322266902e-6,
	-.41126656903020073048963817255e-7,
	.285534422878921522071975766316e-8,
	-.248540830541562387806002659606e-9,
	.254339333807258244274248439717e-10,
	-.294104577282296752348975082791e-11,
	.374339202549390330926505615363e-12,
	-.514911829382116721872054824353e-13,
	.75525359498651439080340407642e-14,
	-.116940970682884644416629062246e-14,
	.189656244943479157172182460506e-15,
	-.320195536869328642066477531639e-16,
	.559954839931620411448416990549e-17,
	-.101021589473043244311939044454e-17,
	.187384498572756298330204271957e-18,
	-.356353747032858021927430144e-19,
	.693128381997123833042276352e-20,
	-.137605945340650015225140893013e-20,
	.2783430784107080220599779328e-21,
	-.572759536432056168934866944e-22,
	.11973614459188926725357568e-22,
	-.253992850989187197664144042667e-23,
	.54613782896572959730696192e-24,
	-.118921134177332028898628949333e-24,
	.2620150977340081594957824e-25,
	-.583681077425568590192093866667e-26,
	.1313743500080595773423616e-26,
	-.298581462251038035533277866667e-27,
	.68483904713346049376256e-28,
	-.158440156822247672119296e-28,
	.369564100657093805430101333333e-29,
	-.868711592114466824301226666667e-30,
	.205708084615876346292906666667e-30,
	-.490522576111622551852373333333e-31
};
	/*
	 * SERIES FOR BT12 ON THE INTERVAL 1.56250e-02 TO  6.25000e-02 WITH
	 * WEIGHTED ERROR        3.33e-32 LOG WEIGHTED ERROR        31.48
	 * SIGNIFICANT FIGURES REQD. 31.05 DECIMAL PLACES REQUIRED   32.27
	 */
static Mfloat lv_bt12cs[39] = {
	.738238601287429746626208397928e+0,
	-.333611131744839063844701476812e-2,
	.614634548880469646985148994202e-4,
	-.240245851616023742649776354696e-5,
	.146635555775097461532105919972e-6,
	-.11841917305589180567005147505e-7,
	.115741989639191970521254663031e-8,
	-.130011611294391874493660077946e-9,
	.162453911413617319377421662737e-10,
	-.220896368214031887521554417701e-11,
	.321803042585531770904743586538e-12,
	-.496531479327684807855520211354e-13,
	.804389004328478259855588826393e-14,
	-.135891213101612913846947126823e-14,
	.23810504397147214869676529606e-15,
	-.430814663638491067244712414208e-16,
	.802025440327710024349935125504e-17,
	-.153163106424623118642300274688e-17,
	.299286063527155689240730405547e-18,
	-.597099646580854433938156366507e-19,
	.121402896694151850241608526507e-19,
	-.251151146966129489010069777067e-20,
	.527905671703287448507383808e-21,
	-.112605092275504983243611613867e-21,
	.243482773595763266596634624e-22,
	-.533172612369318001300384426667e-23,
	.118136150597071210392059904e-23,
	-.264653682833535235148567893333e-24,
	.599033940413615039455778133333e-25,
	-.13690854630829503109136384e-25,
	.315767901543802283264136533333e-26,
	-.734579150820843564914005333333e-27,
	.1722808148072274793070592e-27,
	-.407169079612865079410688e-28,
	.969347451367796227003733333333e-29,
	-.232376363377657167653546666667e-29,
	.560745106735220294068906666667e-30,
	-.136164653915390058605226666667e-30,
	.332631092338946543889066666667e-31
};
	/*
	 * SERIES FOR BM12 ON THE INTERVAL 0.0  TO  1.56250e-02 WITH WEIGHTED
	 * ERROR        5.01e-32 LOG WEIGHTED ERROR        31.30 SIGNIFICANT
	 * FIGURES REQD. 29.99 DECIMAL PLACES REQUIRED   32.10
	 */
static Mfloat lv_bm12cs[40] = {
	.980797915623305002727209354694e-1,
	.11509611895046853061754834846e-2,
	-.431248216433820540988935809773e-5,
	.595183961008881630781302980183e-7,
	-.170484401982690985740070158648e-8,
	.77982654136111095086581738274e-10,
	-.495898612676641580949175495187e-11,
	.403843241642114151683820226514e-12,
	-.399304616372517544576548384665e-13,
	.461988618311896649431334243278e-14,
	-.608920801909538330134547261933e-15,
	.896093091643387648215704804125e-16,
	-.144962942394202312291651891893e-16,
	.254646315853777605616514964807e-17,
	-.480947287464783644425926371862e-18,
	.968768466829259904908727583912e-19,
	-.206721337227796602324503811755e-19,
	.464665155915038473180276780959e-20,
	-.109496612884833413824135132834e-20,
	.269389279728868286090570761279e-21,
	-.689499291093037447781897002686e-22,
	.183026826275206290989066855474e-22,
	-.502506424635191642815611355322e-23,
	.142354519445480603963169363419e-23,
	-.41521912036164503880688867698e-24,
	.124460920150397932588233007655e-24,
	-.382733637056930429943191866129e-25,
	.120559135781561753537472398184e-25,
	-.388453624637648807643185936112e-26,
	.127868952872040972190489528346e-26,
	-.429514668944794627206193691591e-27,
	.147068911782907088645680270798e-27,
	-.51283156651060731281803740178e-28,
	.181950958547116938548143737329e-28,
	-.656303131484198086761863505037e-29,
	.240489897691996065319891487583e-29,
	-.894596674469061247323495824298e-30,
	.337608516065723102663714897824e-30,
	-.129179145462065636091309991697e-30,
	.500863446295881052068495150125e-31,
};
	/*
	 * SERIES FOR BTH1 ON THE INTERVAL 0.0  TO  1.56250e-02 WITH WEIGHTED
	 * ERROR        2.82e-32 LOG WEIGHTED ERROR        31.55 SIGNIFICANT
	 * FIGURES REQD. 31.12 DECIMAL PLACES REQUIRED   32.37
	 */
static Mfloat lv_bth1cs[44] = {
	.747499572035872760554434839697e+0,
	-.124007771446517112525457775414e-2,
	.992524424044245273766414976896e-5,
	-.203036907371597110524193753756e-6,
	.753596177056908857121840175836e-8,
	-.416616127153435501076300238562e-9,
	.307016180708348904812451020912e-10,
	-.281784996376052139923240088839e-11,
	.307906967390402954760281468217e-12,
	-.388033002628034341127873475548e-13,
	.550960396086309049345617262086e-14,
	-.86590060768383779940103398954e-15,
	.148560491415367490034236890607e-15,
	-.27519529815904085805371212125e-16,
	.545507960904810896250362236409e-17,
	-.114865345019836427495436310272e-17,
	.255352133779739002231990525335e-18,
	-.596214901974134503957682879079e-19,
	.145566229023727186202883020058e-19,
	-.370221854224505382015797760196e-20,
	.977630741253453576641684345179e-21,
	-.267268216396684884687237753931e-21,
	.754533003849832717940381906558e-22,
	-.219478999198027448978923833717e-22,
	.656483946239552621789069998175e-23,
	-.201556042983702075707840768695e-23,
	.634177685567761434921446671857e-24,
	-.204192778853378956348137699556e-24,
	.671914642207205674866589800186e-25,
	-.225690791102075735957090036873e-25,
	.772977198929897063709269598719e-26,
	-.269674445122946409132114240809e-26,
	.957493445185026980722955219336e-27,
	-.345691684488901130001756808276e-27,
	.126812348173984365042119862384e-27,
	-.472325366307226398604649937135e-28,
	.178500084781863761778586197964e-28,
	-.684043610045103954062152235668e-29,
	.265660286717204193582934226722e-29,
	-.104504025279144529177141614847e-29,
	.416182908253771443068619171971e-30,
	-.167716392036437148565013478829e-30,
	.683619977766643891735359280285e-31,
	-.281722478612336411667395746228e-31
};


/* Structured by FOR_STRUCT, v0.2, on 06/05/90 at 17:05:34
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  BSJ1/DBSJ1 (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    November 29, 1988

    Purpose:    Evaluate the Bessel function of the first kind of
                order one.

    Usage:      BSJ1(X)

    Arguments:
       X      - Argument for which the function value is desired.
                (Input)
       BSJ1   - Function value.  (Output)

    Remark:
       Informational error:
       Type Code
         2   1  The function underflows because the absolute value
                of X is too small.

    GAMS:       C10a1

    Chapter:    SFUN/LIBRARY Bessel Functions

    Copyright:  1984 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
Mfloat imsl_f_bessel_J1(Mfloat x)
#else
Mfloat imsl_f_bessel_J1(x)
	Mfloat           x;
#endif
{
	Mfloat          ampl, bsj1_v, theta, y;
	static Mint     ntj1 = 0;
	static Mfloat   xsml = 0.0e0;
	static Mfloat   xmin = 0.0e0;
	static Mfloat   xwarn = 0.0e0;

	E1PSH("imsl_f_bessel_J1", "imsl_d_bessel_J1");
	bsj1_v = imsl_amach(6);

	if (ntj1 == 0) {
		ntj1 = imsl_inits(lv_bj1cs, COUNT(lv_bj1cs), 0.1*imsl_amach(3));
		xsml = sqrt(F_FOUR * imsl_amach(3));
		xmin = F_TWO * imsl_amach(1);
		xwarn = F_ONE/sqrt(imsl_amach(3));
	}
	y = fabs(x);

	if (y <= F_FOUR) {
		if (y == F_ZERO) {
			bsj1_v = F_ZERO;
		} else if (y > xsml) {
			bsj1_v = x * (.25e0 + imsl_csevl(.125e0*y*y-1.,
                                    lv_bj1cs, ntj1));
		} else if (y > xmin) {
			bsj1_v = F_HALF * x;
		} else {
                        /* The function underflows because abs(%(r1)) 
                           is less than %(r2).  The result is set to
                           zero. */
			imsl_e1str(1, x);
			imsl_e1str(2, xmin);
                        imsl_ermes(IMSL_ALERT, IMSL_SMALL_ABS_ARG_UNDERFLOW);
		}
	} else {
		l_d9b1mp(&y, &ampl, &theta);
		if ((imsl_n1rty(0) == IMSL_TERMINAL)||(imsl_n1rty(0) == IMSL_FATAL)) goto L_9000;
		(x < 0) ? (bsj1_v = -ampl*cos(theta)) : (bsj1_v = ampl*cos(theta));
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
	E1POP("imsl_f_bessel_J1", "imsl_d_bessel_J1");
	return (bsj1_v);
}				/* end of function */



/* Structured by FOR_STRUCT, v0.2, on 08/19/90 at 12:14:25
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  BSY1/DBSY1 (Single/Double precision version)

    Computer:   FORC/DOUBLE

    Revised:    November 29, 1988

    Purpose:    Evaluate the Bessel function of the second kind of
                order one.

    Usage:      BSY1(X)

    Arguments:
       X      - Argument for which the function value is desired.
                (Input)
       BSY1   - Function value.  (Output)

    GAMS:       C10c

    Chapter:    SFUN/LIBRARY Bessel Functions

    Copyright:  1984 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
Mfloat imsl_f_bessel_Y1(Mfloat temp_x)
#else
Mfloat imsl_f_bessel_Y1(temp_x)
	Mfloat          temp_x;
#endif
{
	Mfloat          ampl, dbsy1_v, theta, x, y;
	/* TWODPI = 2.0/PI */
	static Mfloat   twodpi = 0.63661977236758134307553505349e0;
	static Mint      nty1 = 0;
	static Mfloat   xmin = 0.0e0;
	static Mfloat   xsml = 0.0e0;
        static Mfloat   xwarn = 0.0e0;

	E1PSH("imsl_f_bessel_Y1", "imsl_d_bessel_Y1");
	dbsy1_v = imsl_amach(6);
	x = (Mfloat)temp_x;

	if (nty1 == 0) {
		nty1 = imsl_inits(lv_by1cs, COUNT(lv_by1cs), 0.1*imsl_amach(3));
		xmin = 1.571e0 * exp(imsl_f_max(log(imsl_amach(1)), -log(imsl_amach(2))) +
					  0.01e0);
		xsml = sqrt(F_FOUR * imsl_amach(3));
                xwarn = F_ONE/sqrt(imsl_amach(3));
	}
	if (x <= F_ZERO) {
                /* X = %(r1) must be greater than zero. */
		imsl_e1str(1, x);
		imsl_ermes(IMSL_TERMINAL, IMSL_NON_POSITIVE_ARGUMENT);
	} else if (x <= F_FOUR) {
		if (x >= xmin) {
			if (x <= xsml) {
				y = F_ZERO;
			} else {
				y = x * x;
			}
			dbsy1_v = twodpi * log(F_HALF * x) * imsl_f_bessel_J1(x) + (F_HALF +
			  imsl_csevl(.125e0 * y - F_ONE, lv_by1cs, nty1)) / x;
		} else {
                        /* The function overflows because X = %(r1) is
                          less than %(R2). */
			imsl_e1str(1, x);
			imsl_e1str(2, xmin);
			imsl_ermes(IMSL_FATAL, IMSL_SMALL_ARG_OVERFLOW);
		}
	} else {
		l_d9b1mp(&x, &ampl, &theta);
		if ((imsl_n1rty(0) == IMSL_TERMINAL)||(imsl_n1rty(0) == IMSL_FATAL))
			goto L_9000;
		dbsy1_v = ampl * sin(theta);
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
	E1POP("imsl_f_bessel_Y1", "imsl_d_bessel_Y1");
	return (dbsy1_v);
}				/* end of function */

/* Structured by FOR_STRUCT, v0.2, on 08/19/90 at 13:40:36
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  D9B1MP

    Computer:   FORC/float

    Revised:    January 1, 1984

    Purpose:    Evaluate the modulus and phase for the J1 and Y1
                Bessel functions.

    Usage:      CALL D9B1MP (X, AMPL, THETA)

    Arguments:
       X      - float precision argument for which the modulus and
                phase are desired. (Input)
       AMPL   - float precision modulus for J1 and Y1. (Output)
       THETA  - float precision phase for J1 and Y1. (Output)

    Copyright:  1984 by IMSL, Inc. All rights reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_d9b1mp(Mfloat *x, Mfloat *ampl, Mfloat *theta)
#else
static void l_d9b1mp(x, ampl, theta)
	Mfloat         *x, *ampl, *theta;
#endif
{
	Mfloat          eta;
	Mfloat          z;
	static Mfloat   pi4 = 0.78539816339744830961566084582e0;
	static Mint     nbm1 = 0;
	static Mint     nbt12 = 0;
	static Mint     nbm12 = 0;
	static Mint     nbth1 = 0;
	static Mfloat   xmax = 0.0e0;

	imsl_e1psh("l_d9b1mp");
	*ampl = imsl_amach(6);
	*theta = imsl_amach(6);

	if (nbm1 == 0) {
		eta = 0.1 * imsl_amach(3);
		nbm1 = imsl_inits(lv_bm1cs, COUNT(lv_bm1cs), eta);
		nbt12 = imsl_inits(lv_bt12cs, COUNT(lv_bt12cs), eta);
		nbm12 = imsl_inits(lv_bm12cs, COUNT(lv_bm12cs), eta);
		nbth1 = imsl_inits(lv_bth1cs, COUNT(lv_bth1cs), eta);
		xmax = F_ONE / imsl_amach(4);
	}
	if (*x <= F_EIGHT) {
		z = (128.0e0 / (*x ** x) - F_FIVE) / F_THREE;

		*ampl = (0.75e0 + imsl_csevl(z, lv_bm1cs, nbm1)) / sqrt(*x);
		*theta = *x - F_THREE * pi4 + imsl_csevl(z, lv_bt12cs, nbt12) / *x;

	} else if (*x <= xmax) {
		z = 128.0e0 / (*x ** x) - F_ONE;
		*ampl = (0.75e0 + imsl_csevl(z, lv_bm12cs, nbm12)) / sqrt(*x);
		*theta = *x - F_THREE * pi4 + imsl_csevl(z, lv_bth1cs, nbth1) / *x;

	} else {
		imsl_e1str(1, *x);
		imsl_e1str(2, xmax);
		imsl_ermes(IMSL_FATAL, IMSL_LARGE_ABS_ARG_FATAL);
	}

	imsl_e1pop("l_d9b1mp");
	return;
}				/* end of function */
