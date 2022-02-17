#ifndef __CRIF_PROD_STRUCT_H
#define __CRIF_PROD_STRUCT_H

#include "srt_h_all.h"
#include "LGMSVCalibApprox.h"
#include "CMSSpreadOptions.h"
#include "MCEBOptimisation.h"
#include "CMSCTSProdStruct.h"


/* ------------------------------------------------------------------------------------------------------------------ */
/* Different models */

/* Deal Structures */
/* *************** */

/* Exotic Coupons */
typedef enum CRIF_COUPON_TYPE_
{
	REGULAR_CRIF = 0,
	MODIFIED_CRIF
} CRIF_COUPON_TYPE;

Err crif_interp_coupon_type(	char				*cCouponType,
								CRIF_COUPON_TYPE	*eCouponType);
typedef struct
{
	int						iStartIndex;
	int						iEndIndex;
	int						iNbUsedCpn;

	int						iNumStartIndex;
	int						iNbUsedNumCpn;
	
	int						*iFloatingType; /* 0: no floating, 1: floating */

	/* Coupon Payoff = dAlpha * PreviousCoupon  + dBeta * CMSS + dGamma  */
	double					*dAlpha;
	double					*dBeta;
	double					*dGamma;

	/* For IV control Variate  */
	double					*dCVStrike;

	/* Floor / Cap for Time Swap Coupon */
	int						*iHasFloorCap;

	/* Market Value */
	double					*dMarketValue;

	/* Fixed PV*/
	double					dFixedPV;
	
} crif_exotic_aux, *CRIF_EXOTIC_AUX;

typedef struct
{
	/* Schedule */
	int						iNbCpn;	
	long					*lStartDate;
	long					*lEndDate;
	long					*lPayDate;
	double					*dCoverage;
	double					*dNot;	
							
	/* Coupons */			
	CRIF_COUPON_TYPE		*eCouponType;	

	/* Ratchet gearing */
	double					*dPrvCpngearing;

	/* Paid CMS */			
	int						iFixingLag;
	long					*lPaidFixingDate;
	char					**cPaidTenor;
	char					cPaidFreq[20];
	char					cPaidBasis[20];
	char					cPaidRef[20];
	double					*dPaidGearing;	
	double					*dPaidMargin;
	double					*dFloor;
	double					*dCap;

	double					*dPastPaidCpn;
							

	/* calculation aux */
	CRIF_EXOTIC_AUX		sAux;

} crif_exotic, *CRIF_EXOTIC;

void crif_free_exotic(CRIF_EXOTIC sExotic);

/* Right to Call */
typedef struct
{
	int				iStartIndex;
	int				iNbUsedCall;
	
	/* Calling Details */
	int				*iStartCpnIdx;
	int				*iStartFundIdx;
	
	/* Exercised Details */
	int				iIndexExercised;


} crif_call_aux, *CRIF_CALL_AUX;

typedef struct
{	
	double			dPayRec;	/* -1: Pay, +1: Rec */

	int				iNbCall;
	long			*lExeDate;
	long			*lSettlDate;
	double			*dFee;

	/* Exercised */
	int				iIsExercised;
	long			lExercisedDate;
	int				iIsCallable;

	/* aux for calculation */
	CRIF_CALL_AUX	sAux;

} crif_call, *CRIF_CALL;

void crif_free_call(CRIF_CALL	sCall);


Err crif_interp_calib_strat(const char *constStr, CMSCTS_CALIB_STRAT *val);

/* Calibration / Model Parameters */

typedef struct
{
	long			lStartDate;
	long			lTheoEndDate;
	double			dPayRec;

	/* Funding */
	CMSCTS_FUNDING	sFunding;

	/* Exotic */
	CRIF_EXOTIC	sExotic;

	/* Call */
	CRIF_CALL		sCall;	

} crif_deal, *CRIF_DEAL;


Err	crif_allocate_deal(CRIF_DEAL sDeal);
void crif_free_deal(CRIF_DEAL sDeal);

/* ------------------------------------------------------------------------------------------------------------------ */
/* global variable for CRIF Deal pricing*/
typedef struct 
{
	double			m_dCouponNotional;
	double			m_dFundingNotional;
	double*			m_PreviousCoupon;
	double*			tempvalue;
} CRIF_Global_Structure;

void resetCRIF_Global_Structure();
void freeCRIF_Global_Structure();



/* All Pricing Params */
/* ****************** */

/* For Coupon Pricing */
typedef struct
{
	/* Call Params */
	double		dMCEBCallSpread;
} crif_cpn_numparams, *CRIF_CPN_NUMPARAMS;


typedef struct
{
	CRIF_CPN_NUMPARAMS		sCouponPricingParams;
	CMSCTS_SIMUL_PARAMS		sSimulParams;
	CMSCTS_RESERVE_PARAMS		sReserveParams;

} crif_pricing_params, *CRIF_PRICING_PARAMS;

Err	crif_allocate_pricing_params(CRIF_PRICING_PARAMS sPricingParams);
void crif_free_pricing_params(CRIF_PRICING_PARAMS sPricingParams);

/* Function to check and fill the structures */
Err crif_fill_check_all_struct(	CRIF_DEAL			sDeal,
								CMSCTS_MARKET			sMarket,
								CRIF_CPN_NUMPARAMS	sNumerParams);

typedef struct
{
	long				lNbCall;
	long				*lExeDates;
	double				*dMarketIV;
	double				*dModelIV;
	double				*dExeProba;
	double				*dNewModelIV;
	double				*dFees;
	double				*dOTCall;
	double				*dOTCallPartial;

} crif_fwdiv, *CRIF_FWDIV;


Err	crif_allocate_fwdiv(	CRIF_FWDIV		sFwdIVInfos,
							long			lNbCall);

void crif_free_fwdiv(CRIF_FWDIV	sFwdIVInfos);

typedef struct
{
	double				dFunding;
	double				dExotic;
	double				dCall;
	double				dStd;

	double				dNewIV;
	double				dCallInit;
	double				dSwitchAdjust;
	double				dGMAInitCall;

	CPD_CALIB_INST_DATA sInstDatas;	
	LGMSV_MODEL			sLGMSVModel;

	int					iCalcExeProbas;	

	CRIF_FWDIV		sFwdIVInfos;
	CPD_CALIB_INST_DATA sSwitchInstDatas;	

	CPD_CALIB_INST_DATA sLGM2F_FixedTau_InstDatas;	
	CPD_CALIB_INST_DATA sLGM2F_TauTs_InstDatas;	
} crif_outputs, *CRIF_OUTPUTS;

Err	crif_allocate_outputs(CRIF_OUTPUTS	sOutputs);
void crif_free_outputs(CRIF_OUTPUTS sOutputs);


#endif