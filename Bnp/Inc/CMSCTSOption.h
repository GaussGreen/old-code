#ifndef _CMSCTS_OPTION_H_
#define _CMSCTS_OPTION_H_

#include "CMSCTSProdStruct.h"
#include "CMSCTSCaller.h"

#define	CMSCTS_MAX_FUND_CPN		240
#define	CMSCTS_MAX_CMS_CPN		120
#define	CMSCTS_MAX_ALLCMS_CPN	240
#define	CMSCTS_MAX_CMS			2
#define CMSCTS_MAX_STRIKES	5

typedef struct
{
	int			iIndexCall;	

	/* for funding */
	int			iNbFundCoupon;
	int			iNbFundCouponPartial;
	double		tpay_tstar_alpha[CMSCTS_MAX_FUND_CPN];
	double		tpay_tstar_beta[CMSCTS_MAX_FUND_CPN];
	double		tpay_tstar_gamma[CMSCTS_MAX_FUND_CPN];	

	double		tpay_tstar_beta_2[CMSCTS_MAX_FUND_CPN];
	double		tpay_tstar_gamma_2[CMSCTS_MAX_FUND_CPN];
	double		tpay_tstar_gamma_12[CMSCTS_MAX_FUND_CPN];

	/* for settlement */
	double		tset_tstar_alpha;
	double		tset_tstar_beta;
	double		tset_tstar_gamma;

	double		tset_tstar_beta_2;
	double		tset_tstar_gamma_2;
	double		tset_tstar_gamma_12;

} cmscts_callarg, *CMSCTS_CALLARG;

typedef struct
{
	int			iIndexKO;	

	/* for funding */	
	int			iNbFundCouponPartial;
	double		tpay_tstar_alpha[CMSCTS_MAX_FUND_CPN];
	double		tpay_tstar_beta[CMSCTS_MAX_FUND_CPN];
	double		tpay_tstar_gamma[CMSCTS_MAX_FUND_CPN];	

	double		tpay_tstar_beta_2[CMSCTS_MAX_FUND_CPN];
	double		tpay_tstar_gamma_2[CMSCTS_MAX_FUND_CPN];
	double		tpay_tstar_gamma_12[CMSCTS_MAX_FUND_CPN];

	/* All DF needed by CMS */
	int			iNbDF;
	double		tdf_tstar_alpha[CMSCTS_MAX_ALLCMS_CPN];
	double		tdf_tstar_beta[CMSCTS_MAX_ALLCMS_CPN];
	double		tdf_tstar_gamma[CMSCTS_MAX_ALLCMS_CPN];

	double		tdf_tstar_beta_2[CMSCTS_MAX_ALLCMS_CPN];
	double		tdf_tstar_gamma_2[CMSCTS_MAX_ALLCMS_CPN];
	double		tdf_tstar_gamma_12[CMSCTS_MAX_ALLCMS_CPN];

	double		dAllDF[CMSCTS_MAX_ALLCMS_CPN];

	/* Each CMS reconstruction */
	int			iNbCMSDF[CMSCTS_MAX_CMS];
	int			iIndexCMSDF[CMSCTS_MAX_CMS][CMSCTS_MAX_CMS_CPN];
	double		dCMSCoverage[CMSCTS_MAX_CMS][CMSCTS_MAX_CMS_CPN];
	double		dSpread[CMSCTS_MAX_CMS];

} cmscts_koarg, *CMSCTS_KOARG;

typedef struct
{
	int			iIndexCoupon;

	int			iSavedCoupon;

	/* For payment */
	double		tpay_tstar_alpha;
	double		tpay_tstar_beta;
	double		tpay_tstar_gamma;

	double		tpay_tstar_beta_2;
	double		tpay_tstar_gamma_2;
	double		tpay_tstar_gamma_12;
	
	/* All DF needed by CMS */
	int			iNbDF;
	double		tdf_tstar_alpha[CMSCTS_MAX_ALLCMS_CPN];
	double		tdf_tstar_beta[CMSCTS_MAX_ALLCMS_CPN];
	double		tdf_tstar_gamma[CMSCTS_MAX_ALLCMS_CPN];

	double		tdf_tstar_beta_2[CMSCTS_MAX_ALLCMS_CPN];
	double		tdf_tstar_gamma_2[CMSCTS_MAX_ALLCMS_CPN];
	double		tdf_tstar_gamma_12[CMSCTS_MAX_ALLCMS_CPN];

	double		dAllDF[CMSCTS_MAX_ALLCMS_CPN];

	/* Each CMS reconstruction */
	int			iNbCMSDF[CMSCTS_MAX_CMS];
	int			iIndexCMSDF[CMSCTS_MAX_CMS][CMSCTS_MAX_CMS_CPN];
	double		dCMSCoverage[CMSCTS_MAX_CMS][CMSCTS_MAX_CMS_CPN];
	double		dSpread[CMSCTS_MAX_CMS];

	/* For Payoff */
	double		dAlpha;
	double		dBeta;
	double		dNbOpt[CMSCTS_MAX_STRIKES];

	double		dFloor;
	double		dCap;

} cmscts_couponarg, *CMSCTS_COUPONARG;

typedef struct
{
	int			iIndexFixing;
	int			iIndexCoupon;

	int			iSavedCoupon;

	/* For payment */
	double		tpay_tstar_alpha;
	double		tpay_tstar_beta;
	double		tpay_tstar_gamma;

	double		tpay_tstar_beta_2;
	double		tpay_tstar_gamma_2;
	double		tpay_tstar_gamma_12;

	/* All DF needed by CMS */
	int			iNbDF;
	double		tdf_tstar_alpha[CMSCTS_MAX_ALLCMS_CPN];
	double		tdf_tstar_beta[CMSCTS_MAX_ALLCMS_CPN];
	double		tdf_tstar_gamma[CMSCTS_MAX_ALLCMS_CPN];

	double		tdf_tstar_beta_2[CMSCTS_MAX_ALLCMS_CPN];
	double		tdf_tstar_gamma_2[CMSCTS_MAX_ALLCMS_CPN];
	double		tdf_tstar_gamma_12[CMSCTS_MAX_ALLCMS_CPN];

	double		dAllDF[CMSCTS_MAX_ALLCMS_CPN];

	/* Each CMS reconstruction */
	int			iNbCMSDF[CMSCTS_MAX_CMS];
	int			iIndexCMSDF[CMSCTS_MAX_CMS][CMSCTS_MAX_CMS_CPN];
	double		dCMSCoverage[CMSCTS_MAX_CMS][CMSCTS_MAX_CMS_CPN];
	double		dSpread[CMSCTS_MAX_CMS];

	/* For Payoff */
	double		dAlpha;
	double		dBeta;
	double		dNbOpt[CMSCTS_MAX_STRIKES];

	double		dFloor;
	double		dCap;

} cmscts_fixingarg, *CMSCTS_FIXINGARG;

typedef struct
{
	/* Disctretisation */
	int			iNbTimes;
	double		*dTimes;
	double		*dDates;
	double		*dSigma;
	double		*dAlphaSV;
	double		*dLambdaSV;
	double		*dLevelSV;
	double		*dRhoSV;
	double		*dRho2SV;
	double		*dLGMAlpha;
	double		*dLGMRho;

	/* Payoff infos */
	int			iNbEvent;
	int			*iEvalEvent;
	int			*iIndexEvent;
	void		**vPayoffParams;
	double		*dModelValue;

	/* Path Dependent Infos */
	double		*dPathInfos;
	int			iHasFloorCap;
	double		*dSavedCoupon;
	double		*dSavedDF;

	int			iHasTargetKO;
	double		*dKONotional;
	double		*dKOIndex;

	int			iHasCallFee;
		
	/* MCEB Parameters */
	int			iAdjustIV;
	MCEBPARAMS	sMCEBParams;
	int			*iOptimise;
	double		*dMarketValue;	

} cmscts_mcarg, *CMSCTS_MCARG;

Err	cmscts_allocate_mc_arg(CMSCTS_MCARG	sMCArg);
void cmscts_free_mc_arg(CMSCTS_MCARG sMCArg);

typedef struct
{
	/* Model used for pricing */
	LGMSV_MODEL		sModel;
	
	/* All the arguments for MC */
	CMSCTS_MCARG	sMCArg;

} cmscts_simularg, *CMSCTS_SIMULARG;

Err cmscts_allocate_simularg(CMSCTS_SIMULARG sSimulArg);
void cmscts_free_simularg(CMSCTS_SIMULARG sSimulArg);

typedef struct
{	
	int					iIsCoupon;
	int					iIsFixing;
	int					iIsCall;
	int					iIsKO;

	CMSCTS_DEAL			sDeal;
	CMSCTS_SIMULARG		sSimulArg;
	int					iNumeraireType;
	
	CMSCTS_COUPONARG	sCouponArg;
	CMSCTS_FIXINGARG	sFixingArg;
	CMSCTS_CALLARG		sCallArg;
	CMSCTS_KOARG		sKOArg;

	/* All Path dependent infos */
	double				*dPathInfos;
	double				*dSavedCoupon;
	double				*dSavedDF;	
	double				*dKONotional;
	double				*dKOIndex;

	int					iAdjustPayoff;
	
} cmscts_payarg, *CMSCTS_PAYARG;

Err	cmscts_calibrate_model(	CMSCTS_MARKET	sMarket,
							CMSCTS_DEAL		sDeal,
							CMSCTS_SIMULARG	sSimulArg,
							CMSCTS_CALIB	sCalibration,
							CMSCTS_OUTPUTS	sOutputs);

Err	cmscts_calculate_option(CMSCTS_MARKET			sMarket,
							CMSCTS_DEAL				sDeal,
							CMSCTS_SIMULARG			sSimulArg,						
							CMSCTS_PRICING_PARAMS	sPricingParams,
							CMSCTS_OUTPUTS			sOutputs);

Err	cmscts_calculate_switch_adjust(	CMSCTS_MARKET			sMarket,
									CMSCTS_DEAL				sDeal,
									CMSCTS_SIMULARG			sSimulArg,
									CMSCTS_CALIB			sCalibration,
									CMSCTS_PRICING_PARAMS	sPricingParams,
									CMSCTS_OUTPUTS			sOutputs);

Err cmscts_fill_outputs(CMSCTS_DEAL		sDeal,
						CMSCTS_SIMULARG sSimulArg,
						CMSCTS_OUTPUTS	sOutputs);

#endif