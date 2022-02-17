/* =================================================================================

   FILENAME:	srt_f_types.c

   PURPOSE:     all the names and their equivalent types used all over the 
                library, to represent;
                    - model names
					- curve types
					- underlyings types
					- monte carlo methods
                    - calibration methods

   DESCRIPTION: basically, the functions translate between Strings and 
                enumerated types.
	            Note the way the enumerated types in this are strung together;
	            eventually all srt enumerated types that are translated into 
	            strings should be strung together in this way.
  ================================================================================ */

#include "utallhdr.h"
#include "srt_h_types.h"

/* REMEMBER: defined in srt_h_nametypes.h 
typedef struct{
	long ival;			
	long mval;    		
	String svals[8];   	
	}SrtTrnLstEle;
*/

/* BOOLEAN */
static SrtTrnLstEle boolean_SrtTrnLst[] = {
{	(long)SRT_NO,			0,
	{"NO",0}
},
{	(long)SRT_YES,			
	0,
	{"YES",0}
},
{	(long)0,			0,
	{""}
}
};

/* ----------------------- Curve Type ------------------------------------ */
static SrtTrnLstEle curve_SrtTrnLst[] = {
{	(long)YIELD_CURVE,		0,
	{"YC","YIELD","YIELD_CURVE","YIELDCURVE","YLDCRV",0}
},
{	(long)CMT_CURVE,		0,
	{"CMT_CURVE","CMT","TREASURY_CURVE","TREASURYCURVE",0}
},
{	(long)DVD_CURVE,		0,
	{"DVD_CURVE","DVD", "DIVIDEND_CURVE","DIVIDEND","DVD_CRV",0}
},
{	(long)REPO_CURVE,		0,
	{"REPO_CURVE","REPO", "REPO_CRV",0}
},
{	(long)BOND_CURVE,		0,
	{"BOND_CURVE","BONDCURVE","BC", "BNDCRV",0}
},
{	(long)UNKNOWN_CURVE,		0,
	{"UNKNOWN_CURVE","UNKNOWN",0}
},
{	(long)0,			0,
	{""}
}
};

/* ----------------------- Underlying Type ------------------------------- */
static SrtTrnLstEle underlying_SrtTrnLst[] = {
{	(long)INTEREST_RATE_UND,	0,
	{"IR_UND","INTEREST_RATE_UND","INTEREST_RATE","INTEREST_RATE_UNDERLYING",0}
},
{	(long)FOREX_UND,		0,
	{"FX_UND","FOREX","FOREX_UND","FOREX_UNDERLYING",0}
},
{	(long)EQUITY_UND,		0,
	{"EQ_UND","EQUITY","EQUITY_UND","EQD_UND","EQUITY_UNDERLYING",0}
},
{	(long)UNKNOWN_UND,		0,
	{"UNKNOWN_UND","UNKNOWN",0}
},
{	(long)0,			0,
	{""}
}
};

/* ----------------------- Yield Curve Type ------------------------------- */
static SrtTrnLstEle yieldcurve_SrtTrnLst[] = {
{	(long)SRT_YC,			0,
	{"SRT_YC","SORT_YC","SRT_CURVE","SRT_YIELD_CURVE","ZC_OBJ_CURVE","ZC_OBJ",0}
},
{	(long)MAP_YC,			0,
	{"MAP_YC","SPG_YC","SPG_CURVE","SPG_YIELD_CURVE","MAPINT_CURVE","MAP_CURVE",0}
},
{	(long)WES_YC,			0,
	{"WES_YC","WESTMINSTER_YC","WES_CURVE","WETMINSTER_CURVE",0}
},
{	(long)BND_YC,			0,
	{"BND_YC","BOND_YC","BND_CURVE","BOND_CURVE",0}
},
{	(long)0,			0,
	{""}
}
};

/* -------------------------- Model Type ------------------------------- */
static SrtTrnLstEle model_SrtTrnLst[] = {
{	(long)CHEY,			(long)ONE_FAC,
	{"CHEY1F", "CHEY","LOGNORMAL","CHE","CHEY2D_R_TRUNCATE","CHEYETTE","CHEYETTE1F", 0}
},
{	(long)CHEY_BETA,		(long)ONE_FAC,
	{"CHEYBETA","CHEYBETA1F", "CHEY_BETA","CHEY_BETA1F", "CHEY_BETA_1F","CHEYETTEBETA","CHEYETTE_BETA",0}
},
{	(long)LGM,			(long)ONE_FAC,
	{"LGM1F", "LGM","NORMAL","GAUSS_MARKOV","LINEAR_GAUSS_MARKOV","LGM_1F", 0}
},
{	(long)NEWLGM,			(long)ONE_FAC,
	{"NEWLGM1F", "NEWLGM","NEWNORMAL","NEW_LGM_1F", 0}
},
{	(long)NEWCHEYBETA,			(long)ONE_FAC,
	{"NEWCHEYBETA1F", "NEWCHEYBETA","NEWBETA","NEW_CHEYBETA_1F", 0}
},
{	(long)CHEY1D,		(long)ONE_FAC,
	{"CHEY1D","CHEY_1D","CHE1D","CHE_1D",0}
},
{	(long)ETABETA,			(long)ONE_FAC,
	{"BETAETA", "BE", "BETA_ETA", "SMILE", "EB", "ETA_BETA", "ETABETA", "PAT","LE_MODEL", "THE_ONE", 0}
},
{	(long)BDT,			(long)ONE_FAC,
	{"BDT", "BLACK_DERMAN_TOY","BLACK_KARASINSKI","BK","B_K",0}
},
{	(long)LGM,			(long)TWO_FAC,
	{"LGM2FAC","LGM_2FAC","LGM_2F","LGM2F",0}
},
{	(long)CHEY,		(long)TWO_FAC,
	{"CHE2FAC","CHE_2FAC","CHEY2FAC","CHEY_2FAC","CHEY2F","CHE2F","CHEY_2F",0}
},
{	(long)CHEY_BETA,		(long)TWO_FAC,
	{"CHEYBETA2F","CHEY_BETA2F","CHEY_BETA_2F","CHEYETTEBETA2F","CHEYETTE_BETA2F",0}
},
{	(long)MIXED_BETA,		(long)TWO_FAC,
	{"MIXEDBETA2F","MIXED_BETA2F","MIXBETA2F","MIXED_BETA_2F","MIX_BETA_2F", 0}
},
{	(long)CHEY_STOCH_VOL,			(long)ONE_FAC,
	{"CHEY_STOCH_VOL","CHEY_STOCHVOL1F","CHEYSTOCHVOL1F", "CHEYSTOCHVOL","CHE_STOCHASTIC_VOL","STOCH_VOL_CHE","CSV",0}
},
{	(long)LGM_STOCH_VOL,			(long)ONE_FAC,
	{"LGM_STOCH_VOL","LGM_STOCHASTIC_VOL","STOCH_VOL_LGM","LGMSV",0}
},
{	(long)CHEY_BETA_STOCH_VOL,			(long)ONE_FAC,
	{"CHEY_BETA_STOCH_VOL","CHEYBETA_STOCHVOL1F","CHEYBETASTOCHVOL", "CHEYBETA_STOCH_VOL","CBSV",0}
},
{	(long)BDT,		(long)TWO_FAC,
	{"BDT2FAC","BDT_2FAC","BDT_2F","BDT2F","BK2F","BK_2F","BK_2FAC",0}
	},
{	(long)VASICEK,	(long)ONE_FAC,
	{"VASICEK",0}
},
{	(long)BLACK_SCHOLES,	(long)ONE_FAC,
	{"BLACK_SCHOLES","BLACKSCHOLES","BS","BLACK_SCHOLES","EQD_BS","B_S",0}
},
{	(long)FX_STOCH_RATES,	(long)ONE_FAC,
	{"BS_STOCH_RATES","STOCH_RATES_FX","FX_STOCH_RATES","STOCH_RATES", "STOCHRATES","FX_STOCHRATES",0}
},
{	(long)EQ_STOCH_RATES,	(long)ONE_FAC,
	{"BS_STOCH_RATES","STOCH_RATES_EQ","EQ_STOCH_RATES","STOCH_RATES", "STOCHRATES","EQ_STOCHRATES",0}
},
{	(long)EQ_STOCH_RATES_SRVGS,	(long)ONE_FAC,
	{"EQ_STOCH_RATES_SRVGS",0}
},
{	(long)NORMAL_BS,	(long)ONE_FAC,
	{"NORMAL_BLACK_SCHOLES","NORMAL_BLACK_SCHOLES","NBS","NORMAL_BS","NORMALBS", 0}
},
{	(long)DRIFTED_NORMAL,	(long)ONE_FAC,
	{"DRIFT","DRIFTED_NORMAL","DN","DRIFT_NORMAL","DNORMAL", 0}
},
{	(long)CLEAN_TRIN,			(long)ONE_FAC,
	{"CLEAN_TRIN","CLEAN","BOND_TRIN",0} 
},
{	(long)NONE,			(long)NO_FAC,
	{"NONE","DETERMINISTIC",0}
},
{	(long)0,			0,
	{""}
}
};

/* ---------------- Monte Carlo Discretisation Method ----------------------- */
static SrtTrnLstEle diffscheme_SrtTrnLst[] = {
{	(long)EULER,			0,
	{"EULER","EUL",0}
},
{	(long)MILSHTEIN,		0,
	{"MILSHTEIN","MIL",0}
},
{	(long)0,			0,
	{""}
}
};

/* ----------------------- Random Number Generation Method --------------------- */
static SrtTrnLstEle mcsample_SrtTrnLst[] = {
{	(long)RANDOM_GAUSS,    0,
	{"RG", "RANDOM_GAUSS", 0}
},
{	(long)ANTITHETIC,	0,
	{"ANTI", "ANTITHETIC", 0}
},
{	(long)RANDOM_GAUSS_DYNAMIC,    0,
	{"RANDSAM", "RANDOM_GAUSS_DYNAMIC", "RGD",0}
},
{	(long)ANTITHETIC_DYNAMIC,	0,
	{"ANTISAM", "ANTITHETIC_DYNAMIC", "AD",0}
},
{	(long)BAL_SAMPLING,		0,
	{"BALSAM", "BALANCED", 0}
},
{	(long)BAL_ANTI,		0,
	{"BALANTISAM", "BALANCED_ANTITHETIC", "ANTIBAL","BALANTI",0}
},
{	(long)SOBOL,			0,
	{"SOBOL","SOB",0}
},
{	(long)UWAN,			0,
	{"UWAN","U1",0}
},
{	(long)ABS,			    0,
	{"ABS",0}
},
{	(long)SPECTRUNC,			    0,
	{"SPECTRUNC","SPEC_TRUNC","SPECTRAL_TRUNCATION", "SPECTRAL",
	0}
},
{	(long)0,			0,
	{""}
}

};

/* ----------------------- Calibration Algorithm Name --------------------- */

static SrtTrnLstEle algorithm_SrtTrnLst[] = {
{	(long)LEVENBERG_MARQUARDT,			0,
	{"LM","LEVENBERG","LEVENBERG_MARQUARDT",0}
},
{	(long)SIMPLEX,			0,
	{"SIMPLEX","SIMP",0}
},
{	(long)SIMULATED_ANNEALING,			0,
	{"SA","SIMULATED_ANNEALING","S_A","ANNEALING",0}
},
{	(long)SOBENBERG,			0,
	{"SG","SOBENBERG","SOBNBG","SOBOL_LEVENBERG",0}
},
{	(long)FIXED_POINT,			0,
	{"FP","FIXEDPOINT","FIXED_POINT","FIXPOINT","FIX_POINT","QUICK",0}
},
{	(long)BOOTSTRAP,			0,
	{"BOOTSTRAP","BOOT_STRAP","CALBOOT","BOOT",0}
},
{	(long)BROYDN,			0,
	{"BROYDN","BROY",0}
},
{	(long)0,			0,
	{""}
}

};

/* ---------------------------- Calibration Method ------------------------ */

static SrtTrnLstEle calibration_SrtTrnLst[] = {
{	(long)GLOBAL_CALIB,			0,
	{"DEFAULT","GLOBAL","GLOB","CORRELATION", "CALIBRATE", "OPTIMIZE", "OPTIMISE","FIT",0}
},
{	(long)FIXED_CALIB,			0,
	{"FIXED","FIXED_CORR","FIXED_CORRELATION","FIXEDCORR","FIX_CORR","FIXED","FREEZE",0}
},
{	(long)0,			0,
	{""}
}

};

/* ---------------------------- Computation Measure ------------------------ */

static SrtTrnLstEle measure_SrtTrnLst[] = {
{	(long)SPOT_MEAS,			0,
	{"DEFAULT","SPOT","SPOT_MEASURE","SPOTMEAS", "SPOT_MEAS", "SPOTMEASURE",0}
},
{	(long)FINAL_MEAS,			0,
	{"FINAL","FINAL_MEASURE","FINAL_MEAS","FINALMEASURE","FINALMEAS",0}
},
{	(long)0,			0,
	{""}
}

};
/* ======================================================================== */

static Err srt_f_strtotypes( String s, long *ival, long *mval, SrtTrnLst lst)  
/* chng string to enum */
{
	int i;

	strupper(s);
	strip_white_space(s);

/* Last String is set to NULL: stopping condition */
	while(lst->svals[0][0])
	{
		i=0;
		while(lst->svals[i])
		{
			if(!strcmp(lst->svals[i],s))
			{
				*ival = lst->ival;
				*mval = lst->mval;
			/* Resets the string to the real default value of the string */
				strcpy(s,lst->svals[0]);
				return NULL;
			}
			i++;
		}
		lst++;
	}
   
	return serror("Don't know %s",s);
}

/* ------------------------------------------------------------------------ */
static Err srt_f_typetostr(String s, long  ival)
/* chng enum to string */
{
	static int len = 9;

	static SrtTrnLst lst[10]={0};

	int i;

	lst[0] = boolean_SrtTrnLst;
	lst[1] = curve_SrtTrnLst;
	lst[2] = underlying_SrtTrnLst;
	lst[3] = yieldcurve_SrtTrnLst;
	lst[4] = model_SrtTrnLst;
	lst[5] = diffscheme_SrtTrnLst;
	lst[6] = mcsample_SrtTrnLst;
	lst[7] = algorithm_SrtTrnLst;
	lst[9] = measure_SrtTrnLst;

	for (i=0;i<len;i++)
	{
		while(lst[i]->svals != NULL)
		{
			if(ival == lst[i]->ival)
			{
				strcpy(s, lst[i]->svals[0]);
				return NULL;
			}
			lst[i]++;
		}
	}
	return serror("Bad enumerated type:%d",ival);
	
}		

/* ======================================================================== */


Err srt_f_interp_mcdiffscheme(String sch_str, SrtMCDfSchType *val)
{
	Err err = NULL;
	long tmpLng;
	long tmpLong;
	char tmpStr[SRTSTRBUFSZ];
	SrtTrnLst lst; 
	
	lst = diffscheme_SrtTrnLst;

	strcpy(tmpStr, sch_str);
	strupper(tmpStr);
	strip_white_space(tmpStr);

	err = srt_f_strtotypes(tmpStr,&tmpLng, &tmpLong, lst);
	if (err)
		return err;
	if ((tmpLng<= FIRSTSRTMCDFSCHTYPE) || (tmpLng >=LASTSRTMCDFSCHTYPE)) 	
	{
		return serror("Unknown scheme name: %s",sch_str);
	}
	else
	{
		*val = (SrtMCDfSchType) tmpLng;
	}
	return err;
}


/* ======================================================================== */

Err srt_f_interp_model(String mdl_str, SrtMdlType *typ, SrtMdlDim *dim)
{
	Err err = NULL;
	long tmpLng;
	long tmpLong;
	char tmpStr[SRTSTRBUFSZ];
	SrtTrnLst lst; 

	lst = model_SrtTrnLst;

	strcpy(tmpStr, mdl_str);
	strupper(tmpStr);
	strip_white_space(tmpStr);

	err = srt_f_strtotypes(tmpStr,&tmpLng, &tmpLong, lst);
	if (err)
		return err;
	if ((tmpLng<= FIRSTSRTMDLTYPE) || (tmpLng >=LASTSRTMDLTYPE)) 	
	{
		return serror("Unknown model name: %s",mdl_str);
	}
	else
	{
		*typ = (SrtMdlType) tmpLng;
		*dim = (SrtMdlDim) tmpLong;
	}
	return err;
}


/* ======================================================================== */

Err srt_f_interp_mcsample(String mcs_str, SrtMCSamType *val)
{
	Err err = NULL;
	long tmpLng;
	long tmpLong;
	char tmpStr[SRTSTRBUFSZ];
	SrtTrnLst lst; 

	lst = mcsample_SrtTrnLst;

	strcpy(tmpStr, mcs_str);
	strupper(tmpStr);
	strip_white_space(tmpStr);
	
	err = srt_f_strtotypes(tmpStr,&tmpLng,&tmpLong, lst);
	if (err)
		return err;
	if ((tmpLng<= FIRSTSRTMCSAMTYPE) || (tmpLng >=LASTSRTMCSAMTYPE)) 	
	{
		return serror("Unknown MC sample name: %s",mcs_str);
	}
	else
	{
		*val = (SrtMCSamType)tmpLng;
	}
	return err;
}

/* ======================================================================== */

Err srt_f_interp_curve(String curve_str, SrtCurveType *val)
{
	Err err = NULL;
	long tmpLng;
	long tmpLong;
	char tmpStr[SRTSTRBUFSZ];
	SrtTrnLst lst; 

	lst = curve_SrtTrnLst;

	strcpy(tmpStr, curve_str);
	strupper(tmpStr);
	strip_white_space(tmpStr);
	
	err = srt_f_strtotypes(tmpStr,&tmpLng,&tmpLong, lst);
	if (err)
		return err;
	if ((tmpLng<= FIRSTSRTCURVETYPE) || 
		(tmpLng >=LASTSRTCURVETYPE)) 	
	{
		return serror("Unknown Curve type: %s",curve_str);
	}
	else
	{
		*val = (SrtCurveType)tmpLng;
	}
	return err;
}
/* ======================================================================== */

Err srt_f_interp_under(String und_str, SrtUnderlyingType *val)
{
	Err err = NULL;
	long tmpLng;
	long tmpLong;
	char tmpStr[SRTSTRBUFSZ];
	SrtTrnLst lst; 

	lst = underlying_SrtTrnLst;

	strcpy(tmpStr, und_str);
	strupper(tmpStr);
	strip_white_space(tmpStr);
	
	err = srt_f_strtotypes(tmpStr,&tmpLng,&tmpLong, lst);
	if (err)
		return err;
	if ((tmpLng<= FIRSTSRTUNDERLYINGTYPE) || 
		(tmpLng >=LASTSRTUNDERLYINGTYPE)) 	
	{
		return serror("Unknown Underlying type: %s",und_str);
	}
	else
	{
		*val = (SrtUnderlyingType)tmpLng;
	}
	return err;
}

/* ======================================================================== */

Err srt_f_interp_yc(String crv_str, SrtYieldCurveType *val)
{
	Err err = NULL;
	long tmpLng;
	long tmpLong;
	char tmpStr[SRTSTRBUFSZ];

	SrtTrnLst lst; 

	lst = yieldcurve_SrtTrnLst;

	strcpy(tmpStr, crv_str);
	strupper(tmpStr);
	strip_white_space(tmpStr);
	
	err = srt_f_strtotypes(tmpStr,&tmpLng,&tmpLong, lst);
	if (err)
		return err;
	if ((tmpLng<= FIRSTSRTYIELDCURVETYPE) || 
		(tmpLng >=LASTSRTYIELDCURVETYPE)) 	
	{
		return serror("Unknown Yield Curve Type : %s",crv_str);
	}
	else
	{
		*val = (SrtCurveType)tmpLng;
	}
	return err;
}
/* ======================================================================== */


Err srt_f_interp_calibalgo(String algo_str, SrtCalibAlgoType *val)
{
	Err err = NULL;
	long tmpLng;
	long tmpLong;
	char tmpStr[SRTSTRBUFSZ];
	SrtTrnLst lst; 
	
	lst = algorithm_SrtTrnLst;

	strcpy(tmpStr, algo_str);
	strupper(tmpStr);
	strip_white_space(tmpStr);

	err = srt_f_strtotypes(tmpStr,&tmpLng, &tmpLong, lst);
	if (err)
		return err;
	if ((tmpLng<= FIRSTSRTCALIBALGOTYPE) || (tmpLng >=LASTSRTCALIBALGOTYPE)) 	
	{
		return serror("Unknown calibration algorithm name: %s",algo_str);
	}
	else
	{
		*val = (SrtCalibAlgoType) tmpLng;
	}
	return err;
}

/* ======================================================================== */

Err srt_f_interp_calibtype(String algo_str, SrtCalibType *val)
{
	Err err = NULL;
	long tmpLng;
	long tmpLong;
	char tmpStr[SRTSTRBUFSZ];
	SrtTrnLst lst; 
	
	lst = calibration_SrtTrnLst;

	strcpy(tmpStr, algo_str);
	strupper(tmpStr);
	strip_white_space(tmpStr);

	err = srt_f_strtotypes(tmpStr,&tmpLng, &tmpLong, lst);
	if (err)
		return err;
	if ((tmpLng<= FIRSTSRTCALIBTYPE) || (tmpLng >=LASTSRTCALIBTYPE)) 	
	{
		return serror("Unknown calibration algorithm name: %s",algo_str);
	}
	else
	{
		*val = (SrtCalibType) tmpLng;
	}
	return err;
}

/* ======================================================================== */

Err srt_f_interp_measuretype(String meas_str, SrtMeasureType *val)
{
	Err err = NULL;
	long tmpLng;
	long tmpLong;
	char tmpStr[SRTSTRBUFSZ];
	SrtTrnLst lst; 
	
	lst = measure_SrtTrnLst;

	strcpy(tmpStr, meas_str);
	strupper(tmpStr);
	strip_white_space(tmpStr);

	err = srt_f_strtotypes(tmpStr,&tmpLng, &tmpLong, lst);
	if (err)
		return err;
	if ((tmpLng<= FIRSTSRTMEASURETYPE) || (tmpLng >=LASTSRTMEASURETYPE)) 	
	{
		return serror("Unknown measure name: %s",meas_str);
	}
	else
	{
		*val = (SrtMeasureType) tmpLng;
	}
	return err;
}

/* ================================================================ */

Err srt_f_translate_model(SrtMdlType type, SrtMdlDim dim, String mdl_name)
{
	Err err = NULL;
	SrtTrnLst lst;

	lst = model_SrtTrnLst;

	while(lst->ival != 0)
	{
		if ((type == lst->ival) && (dim == lst->mval))
		{
			strcpy(mdl_name, lst->svals[0]);
			return err;
		}
		lst++;
	}

	return NULL;
}

/* ======================================================================== */


