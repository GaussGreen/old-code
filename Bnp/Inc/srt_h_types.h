/* =================================================================================

   FILENAME:	srt_h_types.h

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
#ifndef SRT_H_NAMETYPES_H
#define SRT_H_NAMETYPES_H

#include "utallhdr.h"
#define 	SRTSTRBUFSZ		128

typedef struct{
	long ival;			/* int val by which is known 		*/
	long mval;    		/* other enumerated type if known (mdl_dim...) 	*/
	String svals[15];   /* array of possible strings, last String is NULL */
	}SrtTrnLstEle;

typedef SrtTrnLstEle *SrtTrnLst; /* last elm in SrtTrnLst is assumed to 
				    							have sval of NULL */

/* ----------------------- Curve Type ------------------------------------ */
typedef enum SrtCurveType{
	FIRSTSRTCURVETYPE = 0,
	YIELD_CURVE,
	CMT_CURVE,
	
	DVD_CURVE,
	REPO_CURVE,
	
	BOND_CURVE,
	UNKNOWN_CURVE,
	LASTSRTCURVETYPE
	} SrtCurveType;

/* ----------------------- Underlying Type ------------------------------- */
typedef enum SrtUnderlyingType{
	FIRSTSRTUNDERLYINGTYPE = LASTSRTCURVETYPE,
	INTEREST_RATE_UND,
	FOREX_UND, 
	EQUITY_UND, 
	BOND_UND,
	UNKNOWN_UND,
	ATMVOL_UND,
	CVM_UND,
	GENMIDAT_UND,
	LASTSRTUNDERLYINGTYPE
	} SrtUnderlyingType;

/* ----------------------- Yield Curve Type ------------------------------- */
typedef enum SrtYieldCurveType{
	FIRSTSRTYIELDCURVETYPE = LASTSRTUNDERLYINGTYPE,
	SRT_YC,
	MAP_YC,
	WES_YC,
	BND_YC,
	XXX_YC,
	LASTSRTYIELDCURVETYPE
	}SrtYieldCurveType;

/* ------------------- Discretisation Method --------------------------- */
typedef enum SrtImpType{
	FIRSTSRTIMPTYPE = LASTSRTCURVETYPE,
	IMPMONTECARLO,
	IMPBACKWARDLATTICE,
	IMPNONE,
	LASTSRTIMPTYPE
	}SrtImpType;
 
/* -------------------------- Model Type ------------------------------- */
typedef enum SrtMdlType{
	FIRSTSRTMDLTYPE=LASTSRTIMPTYPE, 
	LGM,
	NEWLGM,
	NEWCHEYBETA,
	CHEY,
	CHEY_BETA,
	MIXED_BETA,
	CHEY1D, 
	LGM_STOCH_VOL, 
	QUAD_GAUSS,
	CHEY_STOCH_VOL,
	CHEY_BETA_STOCH_VOL,
	VASICEK,
	BLACK_SCHOLES,
	EQ_STOCH_RATES,
	EQ_STOCH_RATES_SRVGS,
	NORMAL_BS,
	FX_STOCH_RATES,
	FX_LGMSV,
	FX_BETADLM,
	DRIFTED_NORMAL,
	ETABETA,
	BDT,
	CLEAN_TRIN,
	NONE,
	DETERMINISTIC,
	BGM,
	BGMSABR,
	HLX_GENERIC,	/* Helix underlyings (Stanley Mrose: 16.10.2002) */
	LGM_IR_CR,
	LASTSRTMDLTYPE
	}SrtMdlType, SrtModelType;

/* -------------------------- Model Dimension ------------------------------- */
typedef enum SrtMdlDim{
	FIRSTSRTMDLDIM = LASTSRTMDLTYPE,
	NO_FAC,
	ONE_FAC,
	TWO_FAC,
	MULTI_FAC,		// multi-factor; added for new BGM (Stanley Mrose: 16.10.2002)
	LASTSRTMDLDIM
	}SrtMdlDim;

/* -------------------------- Clsoed Form Type------------------------------- */
typedef enum SrtClsdFrmType{
	FIRSTSRTCLSDFRMTYPE=LASTSRTMDLDIM, 
	REAL_CLSDFRM,
	QUICK_CLSDFRM,
	LASTSRTCLSDFRMTYPE
	}SrtClsdFrmType, SrtClosedFormType;

/* ----------------------- Monte Carlo Discretisation ----------------------- */
typedef enum SrtMCDfSchType{
	FIRSTSRTMCDFSCHTYPE = LASTSRTCLSDFRMTYPE,
	EULER,
	MILSHTEIN,
	LASTSRTMCDFSCHTYPE
	}SrtMCDfSchType;

/* ----------------------- Random Number Generation Method --------------------- */
typedef enum SrtMCSamType{
	FIRSTSRTMCSAMTYPE = LASTSRTMCDFSCHTYPE,
	RANDOM_GAUSS,
	BAL_SAMPLING,  
	ANTITHETIC, 
	BAL_ANTI, 
	RANDOM_GAUSS_DYNAMIC,
	ANTITHETIC_DYNAMIC,
	SOBOL,
	UWAN,
	ABS,
	SPECTRUNC,
	LASTSRTMCSAMTYPE,
	SOBOLBM
	} SrtMCSamType;

/* ----------------------- Calibration Algorithm Name --------------------- */
typedef enum SrtCalibAlgoType{
	FIRSTSRTCALIBALGOTYPE = LASTSRTMCSAMTYPE,
	LEVENBERG_MARQUARDT,
	SIMPLEX,
	SIMULATED_ANNEALING,
	SOBENBERG,
	FIXED_POINT,
	BOOTSTRAP,
	BROYDN,
	LASTSRTCALIBALGOTYPE
	} SrtCalibAlgoType;

/* ---------------------------- Calibration Method ------------------------ */
typedef enum SrtCalibType{
	FIRSTSRTCALIBTYPE = LASTSRTCALIBALGOTYPE,
	GLOBAL_CALIB,
	BLIND_CALIB,
	FIXED_CALIB,
	LASTSRTCALIBTYPE
	} SrtCalibType;

/* ---------------------------- Measure of Computation ------------------------ */
typedef enum SrtMeasureType{
	FIRSTSRTMEASURETYPE = LASTSRTCALIBTYPE,
	SPOT_MEAS,
	FINAL_MEAS,
	LASTSRTMEASURETYPE
	} SrtMeasureType;

/* ------------------------------------------------------------------------ */

/* ---------------------------------------------------------------------- */
/* Main functions to establish the link between an enum type and a string:
   look at the list of possible strings in the srt_f_strtrnlst.c file */

Err srt_f_interp_mcdiffscheme(String sch_str, SrtMCDfSchType *val);
Err srt_f_interp_model(String mdl_str, SrtMdlType *typ, SrtMdlDim *dim);
Err srt_f_interp_mcsample(String mcs_str, SrtMCSamType *val);
Err srt_f_interp_under(String und_str, SrtUnderlyingType *val);
Err srt_f_interp_curve(String und_str, SrtCurveType *val);
Err srt_f_interp_yc(String crv_str, SrtYieldCurveType *val);
Err srt_f_interp_calibalgo(String algo_str, SrtCalibAlgoType *val);
Err srt_f_interp_calibtype(String calib_str, SrtCalibType *val);
Err srt_f_interp_measuretype(String meas_str, SrtMeasureType *val);
Err srt_f_translate_model(SrtMdlType type, SrtMdlDim dim, String mdl_name);
#endif
