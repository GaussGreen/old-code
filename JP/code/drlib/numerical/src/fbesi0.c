#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif

static Mfloat PROTO(l_dbsi0e,(Mfloat));

#define COUNT(OS)   (sizeof(OS)/sizeof(Mfloat))

#if 0
	/*
	 * SERIES FOR BI0 ON THE INTERVAL 0.0  TO  9.00000D+00 WITH WEIGHTED
	 * ERROR        9.51D-34 LOG WEIGHTED ERROR        33.02 SIGNIFICANT
	 * FIGURES REQD. 33.31 DECIMAL PLACES REQUIRED   33.65
	 */
static Mfloat lv_bi0cs[] = {
     -.766054725283914495108189497624e-1,
     .192733795399380826995240875088e1,
     .228264458692030133893702929233e0,
     .130489146670729042807933421069e-1,
     .434427090081648745137868268103e-3,
     .942265768600193466392317174412e-5,
     .143400628951069107996209187818e-6,
     .161384906966174906991541972e-8,
     .139665004453566969949509270814e-10,
     .957945172550544534462752317189e-13,
     .5333981859862502131015107744e-15,
     .245871608843747077469678592e-17,
     .953568089024877002694434133333e-20,
     .315438203972142733678933333333e-22,
     .900456410109463743146666666667e-25,
     .2240647369123670016e-27,
     .490303460324283733333333333333e-30,
     .950817260612266666666666666667e-33,
};
#endif
	/*
	 * SERIES FOR BI0 ON THE INTERVAL 0.0  TO  9.00000D+00 WITH WEIGHTED
	 * ERROR        9.51D-34 LOG WEIGHTED ERROR        33.02 SIGNIFICANT
	 * FIGURES REQD. 33.31 DECIMAL PLACES REQUIRED   33.65
	 */
static Mfloat lv_bi0cs[] = {
     -.766054725283914495108189497624e-1,
     .192733795399380826995240875088e1,
     .228264458692030133893702929233e0,
     .130489146670729042807933421069e-1,
     .434427090081648745137868268103e-3,
     .942265768600193466392317174412e-5,
     .143400628951069107996209187818e-6,
     .161384906966174906991541972e-8,
     .139665004453566969949509270814e-10,
     .957945172550544534462752317189e-13,
     .5333981859862502131015107744e-15,
     .245871608843747077469678592e-17,
     .953568089024877002694434133333e-20,
     .315438203972142733678933333333e-22,
     .900456410109463743146666666667e-25,
     .2240647369123670016e-27,
     .490303460324283733333333333333e-30,
     .950817260612266666666666666667e-33,
};
     	/*
	 * SERIES FOR AI0 ON THE INTERVAL 1.25000D-01 TO  3.33333D-01 WITH
	 * WEIGHTED ERROR        2.74D-32 LOG WEIGHTED ERROR        31.56
	 * SIGNIFICANT FIGURES REQD. 30.15 DECIMAL PLACES REQUIRED   32.39
	 */
static Mfloat lv_ai0cs[] = {
     .757599449402379594272987203744e-1,
     .75913808108233455072929787332e-2,
     .415313133892375050186319749138e-3,
     .107007646343907307358242970217e-4,
     -.790117997921289466075031948573e-5,
     -.782614350143875226978898980691e-6,
     .278384994294887080638118538986e-6,
     .82524726006120271919668291332e-8,
     -.12044639455201991790549608911e-7,
     .155964859850607644361228752793e-8,
     .229255636710331654347725480286e-9,
     -.119162288427906460367777423448e-9,
     .175785491603240983021833124774e-10,
     .112822446321890051714441135682e-11,
     -.114684862592729887772963387698e-11,
     .271559205480366287264365192161e-12,
     -.241587466656268783844247572028e-13,
     -.608446988825512506460609963922e-14,
     .31457050771754772937083602673e-14,
     -.717221292487118771796217505918e-15,
     .787449340345410339608390960333e-16,
     .100480275300946240234524457184e-16,
     -.756689536535053485342843588881e-17,
     .215038010687611988781205128785e-17,
     -.375485834183087442915158445261e-18,
     .235406584222699257690075710532e-19,
     .111466761204792853022637335511e-19,
     -.539889188439699037869677932271e-20,
     .143959879224075267704285840452e-20,
     -.259191636011109340646081840196e-21,
     .223813318399858390743409229824e-22,
     .5250672575364771172772216832e-23,
     -.324990413853323078417343228587e-23,
     .99242141032050379278572847104e-24,
     -.216499225424466952314655429973e-24,
     .3233609471943594083973332992e-25,
     -.118462020739674248982473386667e-26,
     -.1281671853950498650548338688e-26,
     .582701518227939051160556885333e-27,
     -.1668222326026109719364501504e-27,
     .36253095105415699757006848e-28,
     -.57336279990557135899459584e-29,
     .373679672206309822964258133333e-30,
     .160207398315685196336551253333e-30,
     -.8700424864057229884522496e-31,
     .274132093793748114560341333333e-31,
};
     	/*
	 * SERIES FOR AI02 ON THE INTERVAL 0.0  TO  1.25000D-01 WITH WEIGHTED
	 * ERROR        1.97D-32 LOG WEIGHTED ERROR        31.71 SIGNIFICANT
	 * FIGURES REQD. 30.15 DECIMAL PLACES REQUIRED   32.63
	 */
static Mfloat lv_ai02cs[] = {
     .544904110141088316078960962268e-1,
     .33691164782556940898978566298e-2,
     .688975834691682398426263914301e-4,
     .289137052083475648296692402323e-5,
     .204891858946906374182760534093e-6,
     .226666899049817806459327743136e-7,
     .339623202570838634515084396952e-8,
     .494060238822496958910482449784e-9,
     .118891471078464383424084525196e-10,
     -.314991652796324136453864862962e-10,
     -.132158118404477131187540739927e-10,
     -.179417853150680611777943574027e-11,
     .718012445138366623367106429347e-12,
     .385277838274214270114089801778e-12,
     .15400862175214098269132582334e-13,
     -.415056934728722208662689972016e-13,
     -.955484669882830764870214494313e-14,
     .381168066935262242074605535512e-14,
     .177256013305652638360493266676e-14,
     -.342548561967721913461924790328e-15,
     -.282762398051658348494205593759e-15,
     .346122286769746109309706250813e-16,
     .446562142029675999901042054284e-16,
     -.483050448594418207125525403795e-17,
     -.723318048787475395456227240925e-17,
     .992147541217369859888046093981e-18,
     .119365089084598208550439949924e-17,
     -.24887098371508072357205449166e-18,
     -.193842645416090592898469781133e-18,
     .644465669737344386878301949395e-19,
     .288605159628922432648171383073e-19,
     -.160195490717497180706167156201e-19,
     -.327081501059231472089193567486e-20,
     .368693228382640918114600723939e-20,
     .126829764803095015301359529711e-22,
     -.75498250193772739076963666441e-21,
     .150213357137783534963712789053e-21,
     .126519588350964853493208799248e-21,
     -.6100998370083680708629408916e-22,
     -.126880962926012826436872095924e-22,
     .166101609989074145784038487491e-22,
     -.158519433576588557937970504881e-23,
     -.330264540596821780095381766756e-23,
     .131358090283923978174039623117e-23,
     .36890402466711567933142563728e-24,
     -.42101419104616891492197824725e-24,
     .479195459108286578063171401373e-25,
     .845947039022182179529971707412e-25,
     -.403980094087283249314607937181e-25,
     -.64347146536504313473010085047e-26,
     .122574339887566599034464736991e-25,
     -.293439131602570892319879821175e-26,
     -.196131130919498292620371205729e-26,
     .15035203748221934241622990031e-26,
     -.958872051574482655203386388207e-28,
     -.348333938081704548639441108511e-27,
     .169090361026304367306244960726e-27,
     .198286653873560304389400115719e-28,
     -.531749808149181621457583002528e-28,
     .18033066298883929462350145039e-28,
     .621309334145489317588405311242e-29,
     -.769218929277216186320072806673e-29,
     .185825282611170254262556016596e-29,
     .123758514228139572489927154554e-29,
     -.110225912040922380321779478779e-29,
     .188628711803970449007787447943e-30,
     .216019687224365891314903141406e-30,
     -.160545412491974320058446594966e-30,
     .196535298459429060393884807332e-31,
};     

/* Structured by FOR_STRUCT, v0.2, on 08/19/90 at 12:12:00
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  BSI0/DBSI0 (Single/float precision version)

    Computer:   FORC/float

    Revised:    November 29, 1988

    Purpose:    Evaluate the modified Bessel function of the first kind
                of order zero.

    Usage:      BSI0(X)

    Arguments:
       X      - Argument for which the function value is desired.
                (Input)
       BSI0   - Function value.  (Output)

    GAMS:       C10b1

    Chapter:    SFUN/LIBRARY Bessel Functions

    Copyright:  1984 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
Mfloat imsl_f_bessel_I0(Mfloat x)
#else
Mfloat imsl_f_bessel_I0(x)
	Mfloat          x;
#endif
{
	Mfloat          dbsi0_v, y;
	static Mint     nti0 = 0;
	static Mfloat   xsml = 0.0e0;
	static Mfloat   xmax = 0.0e0;

	E1PSH("imsl_f_bessel_I0", "imsl_d_bessel_I0");
	dbsi0_v = imsl_amach(6);

	if (nti0 == 0) {
		nti0 = imsl_inits(lv_bi0cs, COUNT(lv_bi0cs), 0.1*imsl_amach(3));
		xsml = sqrt(F_EIGHT * imsl_amach(3));
		xmax = log(imsl_amach(2));
	}
	y = fabs(x);

	if (y <= F_THREE) {
		if (y <= xsml) {
			dbsi0_v = F_ONE;
		} else {
			dbsi0_v = 2.75e0 + imsl_csevl(y * y / 4.5e0 - F_ONE, lv_bi0cs, nti0);
		}

	} else {
		if (y <= xmax) {
			dbsi0_v = exp(y) * l_dbsi0e(x);
		} else {
			imsl_e1str(1, x);
			imsl_e1str(2, xmax);
	                imsl_ermes(IMSL_FATAL, IMSL_LARGE_ABS_ARG_FATAL);
		}
	}

	E1POP("imsl_f_bessel_I0", "imsl_d_bessel_I0");
	return (dbsi0_v);
}				/* end of function */


/* Structured by FOR_STRUCT, v0.2, on 08/20/90 at 09:13:12
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  BSI0E/DBSI0E (Single/float precision version)

    Computer:   FORC/float

    Revised:    November 29, 1988

    Purpose:    Evaluate the exponentially scaled modified Bessel
                function of the first kind of order zero.

    Usage:      BSI0E(X)

    Arguments:
       X      - Argument for which the function value is desired.
                (Input)
       BSI0E  - Function value.  (Output)

    GAMS:       C10b1

    Chapter:    SFUN/LIBRARY Bessel Functions

    Copyright:  1984 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static Mfloat l_dbsi0e(Mfloat x)
#else
static Mfloat l_dbsi0e(x)
	Mfloat          x;
#endif
{
	Mfloat          eta;
	Mfloat          dbsi0e_v, y;
	static Mint     nti0 = 0;
	static Mint     ntai0 = 0;
	static Mint     ntai02 = 0;
	static Mfloat   xsml = 0.0e0;

	imsl_e1psh("l_dbsi0e");
	dbsi0e_v = imsl_amach(6);

	if (nti0 == 0) {
		eta = 0.1*imsl_amach(3);
		nti0 = imsl_inits(lv_bi0cs, COUNT(lv_bi0cs), eta);
		ntai0 = imsl_inits(lv_ai0cs, COUNT(lv_ai0cs), eta);
		ntai02 = imsl_inits(lv_ai02cs, COUNT(lv_ai02cs), eta);
		xsml = sqrt(F_EIGHT * imsl_amach(3));
	}
	y = fabs(x);

	if (y <= F_THREE) {
		if (y <= xsml) {
			dbsi0e_v = F_ONE;
		} else {
			dbsi0e_v = exp(-y) * (2.75e0 +
                            imsl_csevl(y * y / 4.5e0 - F_ONE, lv_bi0cs, nti0));
		}

	} else {
		if (y <= F_EIGHT) {
			dbsi0e_v = (.375e0 + imsl_csevl((48.0e0 / y - 11.0e0) 
                                        / F_FIVE, lv_ai0cs, ntai0)) / sqrt(y);
		} else {
			dbsi0e_v = (.375e0 + imsl_csevl(16.0e0 / y - F_ONE,
                                        lv_ai02cs, ntai02)) / sqrt(y);
		}
	}

	imsl_e1pop("l_dbsi0e");
	return (dbsi0e_v);
}				/* end of function */
