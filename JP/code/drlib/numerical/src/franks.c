#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif

#ifdef ANSI
static void	l_r2nks(Mint, Mfloat[], Mfloat, Mint, Mint, Mfloat[], Mint[]);
static void	l_r3nks(Mint, Mint, Mint, Mint*, Mint[], Mfloat[]);
static Mfloat	l_enos(Mint, Mint);
static Mfloat	l_rnunf(void);
static VA_LIST_HACK	l_ranks(Mint, Mfloat[], va_list);
#else
static void	l_r2nks();
static void	l_r3nks();
static Mfloat	l_enos();
static Mfloat	l_rnunf();
static VA_LIST_HACK	l_ranks();
#endif

static Mfloat	*lv_ranks = NULL;

#ifdef ANSI
Mfloat *imsl_f_ranks(Mint n_observations, Mfloat x[], ...)
#else
Mfloat *imsl_f_ranks(n_observations, x, va_alist)
    Mint	n_observations;
    Mfloat	x[];
    va_dcl
#endif
{
    va_list     argptr;

    VA_START(argptr,x);

    E1PSH("imsl_f_ranks","imsl_d_ranks");

    lv_ranks = NULL;
    IMSL_CALL(l_ranks(n_observations, x, argptr));
    va_end(argptr);

    E1POP("imsl_f_ranks","imsl_d_ranks");

    return lv_ranks;
}

#ifdef ANSI
static VA_LIST_HACK l_ranks(Mint n_observations, Mfloat x[], va_list argptr)
#else
static VA_LIST_HACK l_ranks(n_observations, x, argptr)
    Mint	n_observations;
    Mfloat	x[];
    va_list	argptr;
#endif
{
    Mint	    tie_option		= 0;
    Mfloat	    fuzz_value		= F_ZERO;
    Mint	    score_option	= 0;

    Mint	        code = 1, arg_number = 2, itie=0, iscr = 0;
    Mint	user_ranks = 0;
    Mint	*iwk = NULL;

    while (code > 0){
	code = va_arg(argptr, Mint);
        ++arg_number;
	switch (code){
	    case IMSL_AVERAGE_TIE:
                if (!itie) {
                   itie = 1;
                   tie_option = 0;
                   break;
                } else {
                   imsl_ermes(IMSL_TERMINAL, IMSL_MUT_EXC_TIE_OPTION);
                   goto RETURN;
                }
	    case IMSL_HIGHEST:
                if (!itie) {
                   itie = 1;
                   tie_option = 1;
                   break;
                } else {
                   imsl_ermes(IMSL_TERMINAL, IMSL_MUT_EXC_TIE_OPTION);
                   goto RETURN;
                }
	    case IMSL_LOWEST:
                if (!itie) {
                   itie = 1;
                   tie_option = 2;
                   break;
                } else {
                   imsl_ermes(IMSL_TERMINAL, IMSL_MUT_EXC_TIE_OPTION);
                   goto RETURN;
                }
	    case IMSL_RANDOM_SPLIT:
                if (!itie) {
                   itie = 1;
                   tie_option = 3;
                   break;
                } else {
                   imsl_ermes(IMSL_TERMINAL, IMSL_MUT_EXC_TIE_OPTION);
                   goto RETURN;
                }
	    case IMSL_RANKS:
                if (!iscr) {
                   iscr = 1;
                   score_option = 0;
                   break;
                } else {
                   imsl_ermes(IMSL_TERMINAL, IMSL_MUT_EXC_SCORE_OPTION);
                   goto RETURN;
                }
	    case IMSL_BLOM_SCORES:
                if (!iscr) {
                   iscr = 1;
                   score_option = 1;
                   break;
                } else {
                   imsl_ermes(IMSL_TERMINAL, IMSL_MUT_EXC_SCORE_OPTION);
                   goto RETURN;
                }
	    case IMSL_TUKEY_SCORES:
                if (!iscr) {
                   iscr = 1;
                   score_option = 2;
                   break;
                } else {
                   imsl_ermes(IMSL_TERMINAL, IMSL_MUT_EXC_SCORE_OPTION);
                   goto RETURN;
                }
	    case IMSL_VAN_DER_WAERDEN_SCORES:
                if (!iscr) {
                   iscr = 1;
                   score_option = 3;
                   break;
                } else {
                   imsl_ermes(IMSL_TERMINAL, IMSL_MUT_EXC_SCORE_OPTION);
                   goto RETURN;
                }
	    case IMSL_EXPECTED_NORMAL_SCORES:
                if (!iscr) {
                   iscr = 1;
                   score_option = 4;
                   break;
                } else {
                   imsl_ermes(IMSL_TERMINAL, IMSL_MUT_EXC_SCORE_OPTION);
                   goto RETURN;
                }
	    case IMSL_SAVAGE_SCORES:
                if (!iscr) {
                   iscr = 1;
                   score_option = 5;
                   break;
                } else {
                   imsl_ermes(IMSL_TERMINAL, IMSL_MUT_EXC_SCORE_OPTION);
                   goto RETURN;
                }
	    case IMSL_FUZZ:
		fuzz_value = (Mfloat) va_arg(argptr, Mdouble);
		++arg_number;
		break;
	    case IMSL_FUZZ_ADR:
		fuzz_value = *(va_arg(argptr, Mfloat *));
		++arg_number;
		break;
	    case IMSL_RETURN_USER:
		lv_ranks = va_arg(argptr, Mfloat*);
		user_ranks = 1;
		++arg_number;
		break;
	    case 0: 
		break;
	    default:
                imsl_e1sti (1, code);
                imsl_e1sti (2, arg_number);
                imsl_ermes (IMSL_TERMINAL, IMSL_ILLEGAL_OPT_ARG);
                break;
	}
    }
    if (imsl_n1rty(0)) goto RETURN;

    iwk = (Mint *) imsl_malloc (n_observations*sizeof(Mint));
    if (!iwk){
        imsl_e1sti (1,n_observations);
        imsl_e1stl (1,"n_observations");
	imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY_1);
        goto FREE_SPACE;
    }
    if (!lv_ranks) {
	lv_ranks = (Mfloat *) imsl_malloc (n_observations*sizeof(Mfloat));
	if (!lv_ranks) {
            imsl_e1sti (1,n_observations);
            imsl_e1stl (1,"n_observations");
	    imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY_1);
            goto FREE_SPACE;	    
	}
      }

    l_r2nks (n_observations, x, fuzz_value, tie_option, score_option, lv_ranks, iwk);

    if ((code=imsl_n1rty(1))>3 && code!=6 && !user_ranks){
	if (lv_ranks) imsl_free(lv_ranks);
	lv_ranks = NULL;
    }

FREE_SPACE:

    if (iwk) imsl_free (iwk);

RETURN:

    return (argptr);
}
/*
  -----------------------------------------------------------------------
    IMSL Name:  R2NKS/DR2NKS (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    December 31, 1984

    Purpose:    Compute the ranks, normal scores, or exponential scores
                for a vector of observations.

    Usage:      CALL R2NKS (NOBS, X, FUZZ, ITIE, ISCORE, SCORE, IWK)

    Arguments:  See RANKS.

    Chapter:    STAT/LIBRARY Basic Statistics

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_r2nks(Mint nobs, Mfloat x[], Mfloat fuzz, Mint itie, Mint iscore, Mfloat score[], Mint iwk[])
#else
static void l_r2nks(nobs, x, fuzz, itie, iscore, score, iwk)
	Mint            nobs;
	Mfloat          x[], fuzz;
	Mint            itie, iscore;
	Mfloat          score[];
        Mint            iwk[];
#endif
{
	Mint            i, k, err;
	Mfloat          t0, t1, t;


	imsl_e1psh("l_r2nks");

        err = 0;
	if (nobs < 2) {
		imsl_e1sti(1, nobs);

/*		imsl_ermes(5, 1, "NOBS = %(i1).  The number of observations, NOBS, must be greater than or equal to 2.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_NOBS_MUST_BE_GE_TWO);
                ++err;
	}
	if (fuzz < 0) {
		imsl_e1str(1, fuzz);

/*		imsl_ermes(5, 2, "FUZZ = %(r1).  The value used to determine ties, FUZZ, must be greater than or equal to 0.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_NEED_LARGER_FUZZ_VALUE);
                ++err;
	}
	if (itie < 0 || itie > 3) {
		imsl_e1sti(1, itie);

/*		imsl_ermes(5, 4, "ITIE = %(i1).  The option used to assign a rank to a tied observation, ITIE, must be equal to 0, 1, 2, or 3.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_WRONG_ITIE_OPTION_USED);
                ++err;
	}
	if (iscore < 0 || iscore > 5) {
		imsl_e1sti(1, iscore);

/*		imsl_ermes(5, 5, "ISCORE = %(i1).  The option used to give the type of values returned in SCORE, ISCORE, must be equal to 0, 1, 2, 3, 4, or 5.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_WRONG_ISCORE_OPTION_USED);
                ++err;
	}
	if (imsl_n1rcd(0))
		goto L_9000;

	for (i = 1; i <= nobs; i++) {
		iwk[i - 1] = i;
	}
	/* Sort X into SCORE */
	imsl_svrgp(nobs, x, score, iwk);
	/* Compute the rank score */
	if (iscore == 1) {
		t0 = F_THREE / F_EIGHT;
		t1 = nobs + 0.25;
	} else if (iscore == 2) {
		t0 = F_ONE / F_THREE;
		t1 = nobs + F_ONE / F_THREE;
	} else if (iscore == 3) {
		t0 = F_ZERO;
		t1 = nobs + F_ONE;
	}
	for (i = 1; i <= nobs; i++) {
		if (iscore == 0) {
			score[iwk[i - 1] - 1] = i;
		} else if ((iscore == 1 || iscore == 2) || iscore == 3) {
                        t = (i-t0)/t1;
			score[iwk[i - 1] - 1] = imsl_f_normal_inverse_cdf(t);
		} else if (iscore == 4) {
			score[iwk[i - 1] - 1] = l_enos(i, nobs);
		} else if (iscore == 5) {
			if (i > 1) {
				score[iwk[i - 1] - 1] = score[iwk[i - 2] - 1] + F_ONE /
					(nobs - i + 1);
			} else {
				score[iwk[i - 1] - 1] = F_ONE / nobs;
			}
		}
	}
	/* Adjust for ties according to ITIE */
	k = 1;
	for (i = 2; i <= nobs; i++) {
		if (x[iwk[i - 1] - 1] - x[iwk[i - 2] - 1] <= fuzz) {
			k += 1;
		} else {
			l_r3nks(nobs, itie, i, &k, iwk, score);
		}
	}
	l_r3nks(nobs, itie, nobs + 1, &k, iwk, score);
L_9000:
	imsl_e1pop("l_r2nks");
	return;
}				/* end of function */
/*
  -----------------------------------------------------------------------
    IMSL Name:  R3NKS/DR3NKS (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    December 31, 1985

    Purpose:    Adjust ranks according to ties for routine RANKS.

    Usage:      CALL R3NKS (NOBS, ITIE, I, K, IWK, SCORE)

    Arguments:
       NOBS   - Number of elements in SCORE.  (Input)
       ITIE   - Option parameter.  (Input)
                ITIE gives the method to be used for assigning a score
                to tied observations.  These are:
                ITIE  Method
                 0    The average of the scores of the tied observations
                      is used.
                 1    The highest score in the group of ties is used.
                 2    The lowest score in the group of ties is used.
                 3    The tied observations are to be randomly untied.
       I      - Next highest untied rank.  (Input)
                SCORE(IWK(I-K)) through SCORE(IWK(I-1)) are tied.
       K      - Number of tied elements of SCORE.  (Input/ouput)
       IWK    - Vector of length NOBS containing the permutations used
                to rank SCORE.  (Input)
       SCORE  - Vector of length NOBS containing the rank of each
                observation, or a transformation of that rank.
                (Input/output)

    Chapter:    STAT/LIBRARY Basic Statistics

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
/* NOBS is not used here, but leave the calling sequence intact. */
#ifdef ANSI
static void l_r3nks(Mint nobs, Mint itie, Mint i, Mint *k, Mint iwk[], Mfloat score[])
#else
static void l_r3nks(nobs, itie, i, k, iwk, score)
	Mint            nobs, itie, i, *k, iwk[];
        Mfloat          score[];
#endif
{
	Mint            i1, j, j1, j2, l;
	Mfloat          rnk, temp;


	if (*k > 1) {
		j1 = i - *k;
		j2 = i - 1;
		if (itie == 3) {
			i1 = i - *k;
			l = *k;
	L_10:
			j = (int) ((float) (l) * l_rnunf()) + i1;
			temp = score[iwk[i1 + l - 2] - 1];
			score[iwk[i1 + l - 2] - 1] = score[iwk[j - 1] - 1];
			score[iwk[j - 1] - 1] = temp;
			l -= 1;
			if (l <= 1) {
				*k = 1;
			} else {
				goto L_10;
			}
		} else {
			if (itie == 2) {
				rnk = score[iwk[i - *k - 1] - 1];
			} else if (itie == 1) {
				rnk = score[iwk[i - 2] - 1];
			} else if (itie == 0) {
				rnk = F_ZERO;
				for (j = j1; j <= j2; j++) {
					rnk += score[iwk[j - 1] - 1];
				}
				rnk /= *k;
			}
			for (j = j1; j <= j2; j++) {
				score[iwk[j - 1] - 1] = rnk;
			}
		}
	}
	*k = 1;
	return;
}				/* end of function */
/*
  -----------------------------------------------------------------------
    IMSL Name:  ENOS/DENOS (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    February 14, 1985

    Purpose:    Evaluate the expected value of a normal order statistic.

    Usage:      ENOS(I, N)

    Arguments:
       I      - Rank of the order statistic.  (Input)
       N      - Sample size.  (Input)
       ENOS   - Function value, the expected value of the I-th order
                statistic in a sample of size N from the standard normal
                distribution.  (Output)

    Remark:
       Informational errors
       Type Code
         3   1  The rank of the order statistic is less than 1.  A rank
                of 1 is assumed.
         3   2  The rank of the order statistic is greater than sample
                size (N).  A rank of N is assumed.

    Keywords:   Ranks; Normal scores

    GAMS:       L5a2n

    Chapter:    STAT/LIBRARY Mathematical Support

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static Mfloat l_enos(Mint i, Mint n)
#else
static Mfloat l_enos(i, n)
	Mint            i, n;
#endif
{
	Mint            ii, ioe, j, m;
	Mfloat          alogp5, enos_v, f, fff, fi, fn, fr, p, ri, rn, seta, smexe;
	static Mfloat q[] = {
		-.2995732273553875e01,
		-.2302585092993988e01,
		-.1897119984885843e01,
		-.1609437912434071e01,
		-.1386294361119867e01,
		-.1203972804325917e01,
		-.1049822124498661e01,
		-.9162907318741406e00,
		-.7985076962177587e00,
		-.6931471805599337e00,
		-.5978370007556099e00,
		-.510825623765981e00,
		-.4307829160924453e00,
		-.356674943938724e00,
		-.2876820724517732e00,
		-.2231435513142025e00,
		-.162518929497768e00,
		-.1053605156578198e00,
		-.5129329438754439e-01,
		.5773159728050793e-14,
		.4879016416943751e-01,
		.9531017980433021e-01,
		.1397619423751639e00,
		.1823215567939596e00,
		.2231435513142145e00,
		.2623642644674957e00,
		.3001045924503426e00,
		.3364722366212173e00,
		.3715635564324873e00,
		.4054651081081685e00,
		.4382549309311592e00,
		.4700036292457395e00,
		.5007752879124931e00,
		.5306282510621742e00,
		.5596157879354263e00,
		.5877866649021226e00,
		.615185639090237e00,
		.6418538861723983e00,
		.6678293725756588e00,
		.6931471805599486e00,
		.7178397931503201e00,
		.7419373447293805e00,
		.7654678421395746e00,
		.7884573603642733e00,
		.8109302162163319e00,
		.8329091229351071e00,
		.8544153281560706e00,
		.8754687373539028e00,
		.8960880245566385e00,
		.9162907318741579e00,
		.9360933591703376e00,
		.9555114450274392e00,
		.9745596399981336e00,
		.9932517730102861e00,
		.1011600911678482e01,
		.1029619417181161e01,
		.1047318994280562e01,
		.1064710736992431e01,
		.1081805170351731e01,
		.1098612288668112e01,
		.1115141590619323e01,
		.1131402111491103e01,
		.1147402452837544e01,
		.1163150809805683e01,
		.1178654996341648e01,
		.1193922468472437e01,
		.1208960345836977e01,
		.1223775431622118e01,
		.1238374231043271e01,
		.125276296849537e01,
		.1266947603487327e01,
		.1280933845462066e01,
		.1294727167594402e01,
		.1308332819650181e01,
		.1321755839982321e01,
		.1335001066732342e01,
		.1348073148299695e01,
		.1360976553135603e01,
		.1373715578913033e01,
		.1386294361119893e01,
		.139871688111845e01,
		.1410986973710264e01,
		.1423108334242609e01,
		.1435084525289324e01,
		.1446918982936327e01,
		.1458615022699519e01,
		.1470175845100595e01,
		.1481604540924217e01,
		.1492904096178151e01,
		.1504077396776276e01,
		.1515127232962861e01,
		.1526056303495051e01,
		.1536867219599267e01,
		.1547562508716015e01,
		.1558144618046552e01,
		.1568615917913847e01,
		.1578978704949394e01,
		.1589235205116583e01,
		.1599387576580601e01,
		.1609437912434102e01,
		.161938824328727e01,
		.1629240539730282e01,
		.1638996714675647e01,
		.1648658625587383e01,
		.1658228076603534e01,
		.1667706820558078e01,
		.1677096560907917e01,
		.168639895357023e01,
		.1695615608675154e01,
		.1704748092238427e01,
		.1713797927758345e01,
		.1722766597741105e01,
		.1731655545158351e01,
		.1740466174840506e01,
		.1749199854809261e01,
		.1757857917552375e01,
		.1766441661243767e01,
		.1774952350911676e01,
		.178339121955754e01,
		.1791759469228057e01,
		.1800058272042752e01,
		.1808288771179267e01,
		.1816452081818428e01,
		.1824549292051047e01,
		.1832581463748312e01,
		.1840549633397489e01,
		.1848454812904602e01,
		.1856297990365628e01,
		.1864080130807683e01,
		.1871802176901593e01,
		.1879465049647162e01,
		.1887069649032381e01,
		.1894616854667764e01,
		.1902107526396922e01,
		.190954250488444e01,
		.1916922612182063e01,
		.1924248652274135e01,
		.1931521411603215e01,
		.1938741659576702e01,
		.1945910149055315e01,
		.1953027616824179e01,
		.1960094784047271e01,
		.1967112356705918e01,
		.1974081026022011e01,
		.1981001468866585e01,
		.1987874348154347e01,
		.1994700313224747e01,
		.2001480000210126e01,
		.200821403239147e01,
		.2014903020542266e01,
		.2021547563260935e01,
		.2028148247292287e01};
	static Mfloat r[] = {
		-.1476596622751468e-13,
		-.2176037128265329e-13,
		-.3186340080674248e-13,
		-.4662936703425763e-13,
		-.6805667140952435e-13,
		-.9914291609903131e-13,
		-.1438849039914305e-12,
		-.2083888617221634e-12,
		-.3010924842783875e-12,
		-.4338751580236049e-12,
		-.6238343175370696e-12,
		-.8946177132433506e-12,
		-.1279865102788698e-11,
		-.1826427897812512e-11,
		-.2600142323675495e-11,
		-.3692490757607623e-11,
		-.5230926802837565e-11,
		-.739230898719146e-11,
		-.1042099739839586e-10,
		-.1465461085825204e-10,
		-.2055788872489286e-10,
		-.287685431037107e-10,
		-.401599864482697e-10,
		-.5592504237640184e-10,
		-.776885222849778e-10,
		-.1076574385252767e-09,
		-.1488228429491196e-09,
		-.2052263914575559e-09,
		-.2823158374142204e-09,
		-.3874147670086474e-09,
		-.5303423257515197e-09,
		-.7242291213807226e-09,
		-.9865877009111562e-09,
		-.134071243086345e-08,
		-.181750781257116e-08,
		-.2457865025640902e-08,
		-.3315746016587054e-08,
		-.4462172419302919e-08,
		-.5990371439069619e-08,
		-.8022391894107592e-08,
		-.1071759031715677e-07,
		-.1428347994534444e-07,
		-.1898956265833513e-07,
		-.2518491036346734e-07,
		-.3332044906657968e-07,
		-.4397711691582782e-07,
		-.5790134206524873e-07,
		-.7604960806059773e-07,
		-.9964426813848227e-07,
		-.130243238016349e-06,
		-.1698267550907938e-06,
		-.2209050566660918e-06,
		-.2866516130081046e-06,
		-.3710674767571234e-06,
		-.4791833914221242e-06,
		-.6173075625315513e-06,
		-.7933284666339352e-06,
		-.1017083759745548e-05,
		-.1300808299914155e-05,
		-.1659676521578845e-05,
		-.2112456933732011e-05,
		-.2682299376960864e-05,
		-.3397678896843112e-05,
		-.4293523687114815e-05,
		-.5412558555603229e-05,
		-.6806899766264155e-05,
		-.8539941936205909e-05,
		-.1068858289760318e-04,
		-.1334583807120353e-04,
		-.1662390190597559e-04,
		-.2065772028172582e-04,
		-.2560914438540419e-04,
		-.3167174337748929e-04,
		-.3907636006878421e-04,
		-.4809750068383132e-04,
		-.5906065646521043e-04,
		-.7235066117105132e-04,
		-.8842119423939519e-04,
		-.1078055442862845e-03,
		-.1311287514194654e-03,
		-.1591212492720841e-03,
		-.1926341283980515e-03,
		-.2326561413768194e-03,
		-.2803325663186314e-03,
		-.3369860390946651e-03,
		-.4041394552135121e-03,
		-.4835410295067237e-03,
		-.577191585409949e-03,
		-.6873741253905021e-03,
		-.8166857098364631e-03,
		-.9680716434015077e-03,
		-.1144861935423052e-02,
		-.1350809964748202e-02,
		-.1590133239375676e-02,
		-.1867556098180926e-02,
		-.2188354156192187e-02,
		-.2558400247160518e-02,
		-.2984211568374077e-02,
		-.3472997683837415e-02,
		-.4032708994200948e-02,
		-.4672085236439271e-02,
		-.5400703534551409e-02,
		-.622902548585991e-02,
		-.7168442737160194e-02,
		-.823132048231313e-02,
		-.9431038299075455e-02,
		-.1078202773905185e-01,
		-.1229980609144914e-01,
		-.1400100575938433e-01,
		-.1590339871712135e-01,
		-.180259155577274e-01,
		-.2038865869285799e-01,
		-.2301290932896313e-01,
		-.2592112791607553e-01,
		-.2913694784510533e-01,
		-.3268516225556619e-01,
		-.3659170390600892e-01,
		-.4088361815209435e-01,
		-.455890291700683e-01,
		-.5073709965425603e-01,
		-.56357984303983e-01,
		-.6248277749609433e-01,
		-.6914345561223305e-01,
		-.7637281455373084e-01,
		-.8420440303017174e-01,
		-.926724522495226e-01,
		-.1018118026676539e00,
		-.1116578284729239e00,
		-.1222463604874165e00,
		-.133613608160869e00,
		-.1457960813170513e00,
		-.1588305122863204e00,
		-.1727537790234482e00,
		-.1876028297678836e00,
		-.203414609755745e00,
		-.2202259904404449e00,
		-.2380737016233259e00,
		-.2569942668383629e00,
		-.2770239422771288e00,
		-.2981986594829503e00,
		-.3205539719875162e00,
		-.3441250060099932e00,
		-.3689464152886534e00,
		-.3950523400687012e00,
		-.4224763702277728e00,
		-.4512515124827797e00,
		-.4814101615884777e00,
		-.5129840754094268e00,
		-.5460043537227704e00,
		-.5805014205893657e00,
		-.616505010115022e00,
		-.6540441554116743e00,
		-.7338416953693647e00,
		-.7761545927302784e00,
		-.8201120483516769e00,
		-.8657395226816005e00,
		-.9130617648111405e00,
		-.9621028181688563e00,
		-.1012886027819567e01,
		-.1065434049189583e01,
		-.111976885804932e01,
		-.1175911761593625e01,
		-.1233883410469901e01,
		-.1293703811614035e01,
		-.1355392341764083e01,
		-.1418967761531539e01,
		-.1484448229919664e01,
		-.1551851319187785e01,
		-.1621194029969489e01,
		-.1692492806561346e01,
		-.1765763552306998e01,
		-.1841021645009272e01,
		-.1918281952310279e01,
		-.1997558846986307e01,
		-.2078866222110689e01,
		-.216221750604375e01,
		-.2247625677214329e01,
		-.2335103278662453e01,
		-.2424662432317242e01,
		-.2516314852988339e01,
		-.2610071862052941e01,
		-.2705944400823903e01,
		-.280394304358751e01,
		-.2904078010302258e01,
		-.300635917895243e01,
		-.3110796097552495e01,
		-.3217397995800287e01,
		-.3326173796378607e01,
		-.3437132125906423e01,
		-.3550281325542148e01,
		-.3665629461242573e01,
		-.3783184333682046e01,
		-.3902953487837305e01,
		-.4024944222244003e01,
		-.4149163597931659e01,
		-.427561844704416e01,
		-.4404315381153292e01,
		-.4535260799273172e01,
		-.4668460895583632e01,
		-.4803921666870696e01,
		-.4941648919692438e01,
		-.5081648277278704e01,
		-.5223925186172927e01,
		-.5368484922624293e01,
		-.5515332598738667e01,
		-.5664473168396457e01,
		-.5815911432945259e01,
		-.5969652046675234e01,
		-.6125699522085338e01,
		-.6284058234947415e01,
		-.6444732429175935e01,
		-.6607726221510343e01,
		-.6773043606017728e01,
		-.6940688458421271e01,
		-.7110664540262647e01,
		-.7282975502903613e01,
		-.7457624891374113e01,
		-.7634616148071119e01,
		-.7813952616316348e01,
		-.799563754377673e01,
		-.8179674085752888e01,
		-.8366065308344028e01,
		-.8554814191487645e01,
		-.874592363188728e01,
		-.8939396445825637e01,
		-.913523537187268e01,
		-.9333443073489079e01,
		-.9534022141533047e01,
		-.9736975096667792e01,
		-.9942304391691628e01,
		-.1015001241375612e02,
		-.1036010148652729e02,
		-.1057257387225042e02,
		-.1078743177374731e02,
		-.1100467733631163e02,
		-.1122431264960134e02,
		-.1144633974936864e02,
		-.1167076061919396e02,
		-.1189757719215191e02,
		-.1212679135238504e02,
		-.123584049366392e02,
		-.1259241973571053e02,
		-.1282883749596439e02,
		-.1306765992055022e02,
		-.1330888867094168e02,
		-.1355252536795449e02,
		-.1379857159319925e02,
		-.1404702889012706e02,
		-.1429789876529278e02,
		-.1455118268930635e02,
		-.1480688209835897e02,
		-.1506499839383404e02,
		-.1532553294603481e02,
		-.1558848709213367e02,
		-.1585386213820371e02,
		-.1612165936169902e02,
		-.1639188000998599e02,
		-.1666452530256692e02,
		-.1693959643039488e02,
		-.1721709455901047e02,
		-.1749702082947138e02,
		-.1777937635198566e02,
		-.1806416222321361e02,
		-.1835137949577754e02,
		-.1864102922238369e02,
		-.189331124198754e02,
		-.1922763010220532e02,
		-.1952458319757374e02,
		-.19823972740809e02,
		-.2012579960891278e02,
		-.2043006470011414e02,
		-.2073676889383496e02,
		-.2104591330797427e02,
		-.2135749842050528e02,
		-.2167152524762946e02,
		-.2198799468102285e02,
		-.2230690739766362e02,
		-.2262826449094309e02,
		-.2295206679539271e02,
		-.2327831358784527e02,
		-.2360700885084431e02,
		-.239381499780087e02,
		-.2427173857908781e02,
		-.2460777636910266e02,
		-.2494626609563959e02,
		-.2528719836457364e02,
		-.2563058098225117e02,
		-.2597643264459843e02,
		-.263247198835277e02,
		-.2667545493252539e02,
		-.2702865902505768e02,
		-.2738426643199774e02,
		-.2774237990392581e02,
		-.281028915785649e02,
		-.2846601955651443e02,
		-.2883135892061681e02,
		-.2919937053309059e02,
		-.2956976269276488e02,
		-.299422139888006e02,
		-.3031843563374089e02,
		-.3069654585839969e02,
		-.3107731835391748e02,
		-.3145868591044658e02,
		-.3184645144145535e02};
	static Mfloat s[] = {
		-.2979893853320468e02,
		-.2942018853320467e02,
		-.2904393853320467e02,
		-.2867018853320467e02,
		-.2829893853320467e02,
		-.2793018853320467e02,
		-.2756393853320467e02,
		-.2720018853320467e02,
		-.2683893853320467e02,
		-.2648018853320467e02,
		-.2612393853320467e02,
		-.2577018853320467e02,
		-.2541893853320467e02,
		-.2507018853320467e02,
		-.2472393853320467e02,
		-.2438018853320467e02,
		-.2403893853320467e02,
		-.2370018853320467e02,
		-.2336393853320467e02,
		-.2303018853320467e02,
		-.2269893853320467e02,
		-.2237018853320467e02,
		-.2204393853320467e02,
		-.2172018853320467e02,
		-.2139893853320467e02,
		-.2108018853320467e02,
		-.2076393853320467e02,
		-.2045018853320467e02,
		-.2013893853320467e02,
		-.1983018853320467e02,
		-.1952393853320466e02,
		-.1922018853320467e02,
		-.1891893853320467e02,
		-.1862018853320466e02,
		-.1832393853320467e02,
		-.1803018853320466e02,
		-.1773893853320467e02,
		-.1745018853320467e02,
		-.1716393853320466e02,
		-.1688018853320466e02,
		-.1659893853320466e02,
		-.1632018853320466e02,
		-.1604393853320466e02,
		-.1577018853320466e02,
		-.1549893853320466e02,
		-.1523018853320466e02,
		-.1496393853320466e02,
		-.1470018853320466e02,
		-.1443893853320466e02,
		-.1418018853320466e02,
		-.1392393853320466e02,
		-.1367018853320466e02,
		-.1341893853320466e02,
		-.1317018853320466e02,
		-.1292393853320466e02,
		-.1268018853320466e02,
		-.1243893853320466e02,
		-.1220018853320466e02,
		-.1196393853320466e02,
		-.1173018853320466e02,
		-.1149893853320466e02,
		-.1127018853320466e02,
		-.1104393853320466e02,
		-.1082018853320466e02,
		-.1059893853320466e02,
		-.1038018853320466e02,
		-.1016393853320466e02,
		-.995018853320466e01,
		-.973893853320466e01,
		-.953018853320466e01,
		-.932393853320466e01,
		-.912018853320466e01,
		-.891893853320466e01,
		-.872018853320466e01,
		-.852393853320466e01,
		-.833018853320466e01,
		-.813893853320466e01,
		-.795018853320466e01,
		-.776393853320466e01,
		-.758018853320466e01,
		-.739893853320466e01,
		-.722018853320466e01,
		-.704393853320466e01,
		-.687018853320466e01,
		-.669893853320466e01,
		-.653018853320466e01,
		-.636393853320466e01,
		-.620018853320466e01,
		-.603893853320466e01,
		-.588018853320466e01,
		-.572393853320466e01,
		-.5570188533204661e01,
		-.5418938533204661e01,
		-.5270188533204661e01,
		-.5123938533204661e01,
		-.4980188533204661e01,
		-.4838938533204661e01,
		-.4700188533204661e01,
		-.4563938533204661e01,
		-.4430188533204661e01,
		-.4298938533204661e01,
		-.4170188533204661e01,
		-.4043938533204662e01,
		-.3920188533204662e01,
		-.3798938533204662e01,
		-.3680188533204662e01,
		-.3563938533204662e01,
		-.3450188533204662e01,
		-.3338938533204662e01,
		-.3230188533204662e01,
		-.3123938533204662e01,
		-.3020188533204663e01,
		-.2918938533204663e01,
		-.2820188533204663e01,
		-.2723938533204663e01,
		-.2630188533204663e01,
		-.2538938533204663e01,
		-.2450188533204664e01,
		-.2363938533204664e01,
		-.2280188533204664e01,
		-.2198938533204664e01,
		-.2120188533204665e01,
		-.2043938533204665e01,
		-.1970188533204665e01,
		-.1898938533204665e01,
		-.1830188533204665e01,
		-.1763938533204665e01,
		-.1700188533204666e01,
		-.1638938533204666e01,
		-.1580188533204666e01,
		-.1523938533204666e01,
		-.1470188533204667e01,
		-.1418938533204667e01,
		-.1370188533204667e01,
		-.1323938533204668e01,
		-.1280188533204668e01,
		-.1238938533204668e01,
		-.1200188533204668e01,
		-.1163938533204669e01,
		-.1130188533204669e01,
		-.1098938533204669e01,
		-.1070188533204669e01,
		-.104393853320467e01,
		-.102018853320467e01,
		-.9989385332046704e00,
		-.9801885332046707e00,
		-.9639385332046709e00,
		-.9501885332046712e00,
		-.9389385332046715e00,
		-.9301885332046718e00,
		-.9239385332046721e00,
		-.9201885332046724e00,
		-.920188533204673e00,
		-.9239385332046733e00,
		-.9301885332046736e00,
		-.9389385332046738e00,
		-.9501885332046741e00,
		-.9639385332046744e00,
		-.9801885332046747e00,
		-.998938533204675e00,
		-.1020188533204675e01,
		-.1043938533204676e01,
		-.1070188533204676e01,
		-.1098938533204676e01,
		-.1130188533204676e01,
		-.1163938533204677e01,
		-.1200188533204677e01,
		-.1238938533204677e01,
		-.1280188533204678e01,
		-.1323938533204678e01,
		-.1370188533204678e01,
		-.1418938533204678e01,
		-.1470188533204679e01,
		-.1523938533204679e01,
		-.1580188533204679e01,
		-.163893853320468e01,
		-.170018853320468e01,
		-.1763938533204681e01,
		-.1830188533204681e01,
		-.1898938533204681e01,
		-.1970188533204682e01,
		-.2043938533204682e01,
		-.2120188533204682e01,
		-.2198938533204683e01,
		-.2280188533204683e01,
		-.2363938533204684e01,
		-.2450188533204684e01,
		-.2538938533204684e01,
		-.2630188533204685e01,
		-.2723938533204685e01,
		-.2820188533204686e01,
		-.2918938533204686e01,
		-.3020188533204686e01,
		-.3123938533204687e01,
		-.3230188533204687e01,
		-.3338938533204688e01,
		-.3450188533204688e01,
		-.3563938533204689e01,
		-.3680188533204689e01,
		-.3798938533204689e01,
		-.392018853320469e01,
		-.404393853320469e01,
		-.4170188533204691e01,
		-.4298938533204691e01,
		-.4430188533204692e01,
		-.4563938533204692e01,
		-.4700188533204693e01,
		-.4838938533204693e01,
		-.4980188533204694e01,
		-.5123938533204694e01,
		-.5270188533204695e01,
		-.5418938533204695e01,
		-.5570188533204696e01,
		-.5723938533204696e01,
		-.5880188533204697e01,
		-.6038938533204697e01,
		-.6200188533204698e01,
		-.6363938533204698e01,
		-.6530188533204699e01,
		-.66989385332047e01,
		-.68701885332047e01,
		-.7043938533204701e01,
		-.7220188533204701e01,
		-.7398938533204702e01,
		-.7580188533204702e01,
		-.7763938533204703e01,
		-.7950188533204703e01,
		-.8138938533204704e01,
		-.8330188533204704e01,
		-.8523938533204705e01,
		-.8720188533204706e01,
		-.8918938533204706e01,
		-.9120188533204706e01,
		-.9323938533204708e01,
		-.9530188533204708e01,
		-.9738938533204708e01,
		-.995018853320471e01,
		-.1016393853320471e02,
		-.1038018853320471e02,
		-.1059893853320471e02,
		-.1082018853320471e02,
		-.1104393853320471e02,
		-.1127018853320471e02,
		-.1149893853320471e02,
		-.1173018853320471e02,
		-.1196393853320471e02,
		-.1220018853320471e02,
		-.1243893853320472e02,
		-.1268018853320472e02,
		-.1292393853320472e02,
		-.1317018853320472e02,
		-.1341893853320472e02,
		-.1367018853320472e02,
		-.1392393853320472e02,
		-.1418018853320472e02,
		-.1443893853320472e02,
		-.1470018853320472e02,
		-.1496393853320472e02,
		-.1523018853320472e02,
		-.1549893853320472e02,
		-.1577018853320472e02,
		-.1604393853320473e02,
		-.1632018853320473e02,
		-.1659893853320473e02,
		-.1688018853320473e02,
		-.1716393853320473e02,
		-.1745018853320473e02,
		-.1773893853320473e02,
		-.1803018853320473e02,
		-.1832393853320473e02,
		-.1862018853320473e02,
		-.1891893853320473e02,
		-.1922018853320473e02,
		-.1952393853320473e02,
		-.1983018853320474e02,
		-.2013893853320474e02,
		-.2045018853320474e02,
		-.2076393853320474e02,
		-.2108018853320474e02,
		-.2139893853320474e02,
		-.2172018853320474e02,
		-.2204393853320474e02,
		-.2237018853320474e02,
		-.2269893853320475e02,
		-.2303018853320474e02,
		-.2336393853320474e02,
		-.2370018853320474e02,
		-.2403893853320475e02,
		-.2438018853320475e02,
		-.2472393853320475e02,
		-.2507018853320475e02,
		-.2541893853320475e02,
		-.2577018853320475e02,
		-.2612393853320475e02,
		-.2648018853320475e02,
		-.2683893853320475e02,
		-.2720018853320476e02,
		-.2756393853320475e02,
		-.2793018853320476e02,
		-.2829893853320476e02,
		-.2867018853320476e02,
		-.2904393853320476e02,
		-.2942018853320476e02,
		-.2979893853320476e02};

	imsl_e1psh("l_enos");
	/* Compute Integral */
	ii = i;
	enos_v = F_ZERO;
	if (n < 1) {
		imsl_e1sti(1, n);

/*		imsl_ermes(5, 1, "The sample size must be at least 1 while N = %(i1) is given.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_NEED_LARGER_SAMPLE_SIZE);
		goto L_9000;
	}
	if (i < 1) {
		imsl_e1sti(1, i);

/*		imsl_ermes(3, 1, "The rank must be at least 1 while  I = %(i1) is given.  I = 1 is used.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_RANK_MUST_BE_AT_LEAST_ONE);
		ii = 1;
	}
	if (i > n) {
		imsl_e1sti(1, i);
		imsl_e1sti(2, n);

/*		imsl_ermes(3, 2, "The rank must be less than or equal to the sample size while I = %(i1) and N = %(i2)  are given.  I = %(i2) is used.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_NEED_RANK_LE_SAMPLE_SIZE);
		ii = n;
	}
	ioe = mod(n, 2);
	m = (n + 1) / 2;
	if (ii - m == 0 && ioe != 0) {
		enos_v = F_ZERO;
		goto L_9000;
	}
	seta = imsl_amach(1);
	smexe = log(seta);
	alogp5 = log(F_HALF);
	f = F_ZERO;
	fn = (float) (n + 1);
	fr = (float) (n + 1 - ii);
	fi = (float) (ii);
	fn = imsl_f_log_gamma(fn);
	fr = imsl_f_log_gamma(fr);
	fi = imsl_f_log_gamma(fi);
	fff = fn - fr - fi;
	ri = (float) (ii - 1);
	rn = (float) (n - ii);
	p = fff + q[151] + ri * r[0] + rn * r[303] + s[0];
	if (p >= smexe + alogp5) {
		f = -F_HALF * exp(p);
	}
	for (j = 2; j <= 152; j++) {
		p = fff + q[-j + 152] + ri * r[j - 1] + rn * r[-j + 304] + s[j - 1];
		if (p >= smexe) {
			f -= exp(p);
		}
	}
	for (j = 153; j <= 303; j++) {
		p = fff + q[j - 153] + ri * r[j - 1] + rn * r[-j + 304] + s[j - 1];
		if (p >= smexe) {
			f += exp(p);
		}
	}
	p = fff + q[151] + ri * r[303] + rn * r[0] + s[303];
	if (p >= smexe + alogp5)
		f += F_HALF * exp(p);
	enos_v = -.05e0 * f;
L_9000:
	imsl_e1pop("l_enos");
	return (enos_v);
}				/* end of function */

#ifdef ANSI
static Mfloat l_rnunf(void)
#else
static Mfloat l_rnunf()
#endif
{
    Mfloat r[1];

    imsl_f_random_uniform (1, IMSL_RETURN_USER, r, 0);
    return r[0];
}

