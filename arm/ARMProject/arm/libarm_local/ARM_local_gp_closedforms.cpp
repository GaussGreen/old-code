/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * $Log: ARM_local_gp_closedforms.cpp,v $
 * Revision 1.1  2004/05/02 15:08:43  ocroissant
 * Initial version
 *
 */

#include "firstToBeIncluded.h"

#include "ARM_local_wrapper.h"
#include "ARM_local_gp_closedforms.h"

#include <GP_ClosedForms\gpclosedforms\argconvdefault_cf.h>
#include <GP_ClosedForms\gpclosedforms\vanilla_normal.h>
#include <GP_ClosedForms\gpclosedforms\spreadoption_lognormal_interface.h>
#include <GP_ClosedForms\gpclosedforms\spreadoption_normal_interface.h>
#include <GP_ClosedForms\gpclosedforms\spreadoption_sabr_interface.h>
#include <GP_ClosedForms\gpclosedforms\spreadoption_shiftedlognormal_interface.h>
#include <GP_ClosedForms\gpclosedforms\merton_interface.h>
#include <GP_ClosedForms\gpclosedforms\vanille_bs_interface.h>
#include <GP_ClosedForms\gpclosedforms\vanille_normal_interface.h>
#include <GP_ClosedForms\gpclosedforms\extended_sabr_interface.h>
#include <GP_ClosedForms\gpclosedforms\barriere_bs_interface.h>
#include <GP_ClosedForms\gpclosedforms\heston_interface.h>
#include <GP_ClosedForms\gpclosedforms\gaussian_integrals.h>
#include <GP_ClosedForms\gpclosedforms\tri_spreadoption_lognormal_interface.h>
#include <GP_ClosedForms\gpclosedforms\CEV_interface.h>
#include <GP_ClosedForms\gpclosedforms\asian_lognormal_interface.h>
#include <GP_ClosedForms\gpclosedforms\numerics_interface.h>
#include <GP_ClosedForms\gpclosedforms\gaussian_integrals.h>
#include <GP_ClosedForms\gpclosedforms\erf.h>
#include <GP_ClosedForms\gpclosedforms\gamma.h>
#include <GP_ClosedForms\gpclosedforms\incompletebeta.h>
#include <GP_ClosedForms\gpclosedforms\hypergeometric.h>
#include <GP_ClosedForms\gpclosedforms\convertor.h>
#include <GP_ClosedForms\gpclosedforms\sabr_calibration.h>
#include <GP_ClosedForms\gpclosedforms\qgm_skewcalibration.h>
#include <GP_ClosedForms\gpclosedforms\bisabr_calibration.h>
#include <GP_ClosedForms\gpclosedforms\cppi_options.h>
#include <GP_ClosedForms\gpclosedforms\heston_calibration.h>
#include <GP_ClosedForms\gpclosedforms\StochasticVol_LN.h>
#include <GP_ClosedForms\gpclosedforms\stochasticvol_ln_interface.h>
#include <GP_ClosedForms\gpclosedforms\stochasticvol_N.h>
#include <GP_ClosedForms\gpclosedforms\glambda_calibration.h>
#include <GP_ClosedForms\gpclosedforms\spreadoption_glambda_interface.h>
#include <GP_ClosedForms\gpclosedforms\normal.h>
#include <GP_ClosedForms\gpclosedforms\numerics_interface.h>
#include <GP_ClosedForms\gpclosedforms\vanille_normal_interface.h>
#include <GP_ClosedForms\gpclosedforms\CIRBond.h>
#include <GP_ClosedForms\gpclosedforms\mepi_interface.h>
#include <GP_ClosedForms\gpclosedforms\tridiagonalsolve.h>
#include <GP_ClosedForms\gpclosedforms\spreadoption_bisabr_interface.h>
#include <GP_ClosedForms\gpclosedforms\bisabr_interface.h>
#include <GP_ClosedForms\gpclosedforms\tarnProxy.h>
#include <GP_ClosedForms\gpclosedforms\VolBondMinMaxProxy.h>
#include <GP_ClosedForms\gpclosedforms\Berm2DatesProxy.h>
#include <GP_ClosedForms\gpclosedforms\SpreadVolBondProxy.h>
#include <GP_ClosedForms\gpclosedforms\smile_calibration.h>
#include <GP_ClosedForms\gpclosedforms\trisabr_interface.h>
#include <GP_ClosedForms\gpclosedforms\spreadoption_nonparametric_interface.h>
#include <GP_ClosedForms\gpclosedforms\distribution_interface.h>
#include <GP_ClosedForms\gpclosedforms\biheston_interface.h>
#include <GP_ClosedForms\gpclosedforms\merton.h>
#include <GP_ClosedForms\gpclosedforms\ExpRiccati.h>

#include <GP_Calib\gpcalib\typedef.h>
#include <GP_Calib\gpcalib\densityfunctors.h>
#include <GP_Models\gpmodels\sabr_eq.h>
#include "gpnumlib/gaussiananalytics.h"
#include "gpbase/datemanip.h"
#include "gpinfra/modelparam.h"

#include "ARM_local_class.h"
#include "ARM_local_persistent.h"
#include "ARM_result.h"
#include "ARM_local_glob.h"
#include <glob\expt.h>
#include <ARM\local_xlarm\ARM_local_interglob.h>
#include <GP_Base\gpbase\gplinalgconvert.h>
#include <GP_Base\gpbase\gpvector.h>
#include <vector>
#include <GP_Base\gpbase\gpmatrix.h>


using std::vector;

using ARM::ConvertXLDateToJulian;
using ARM::ConvertJulianToXLDate;


using ARM::Export_normal_VanillaOption;
using ARM::Export_LogNormal_SpreadOption;
using ARM::Export_Smiled_LogNormal_SpreadOption;
using ARM::Export_LogNormal_SpreadOption_Calibrate_Correlation;
using ARM::Export_Smiled_LogNormal_SpreadOption_Calibrate_Correlation;
using ARM::Export_Normal_SpreadOption;
using ARM::Export_Smiled_Normal_SpreadOption;
using ARM::Export_SABR_ImplicitVol;
using ARM::Export_SABR2B_ImplicitVol;
using ARM::Export_Gaussian_SABR_Power_SpreadOption;
using ARM::Export_Gaussian_SABR_Power_Digital_SpreadOption;
using ARM::Export_Gaussian_SABR_Power_SpreadOption_Certitude;
using ARM::Export_BlackSholes;
using ARM::Export_BlackSholes_ImplicitVol;
using ARM::Export_Gaussian_ShiftedLN_Power_SpreadOption;
using ARM::Export_Gaussian_ShiftedLN_Power_SpreadOption_Certitude;
using ARM::Export_Merton_JumpDiffusion;
using ARM::Export_SABR_VanillaOption;
using ARM::Export_SABR_FromGaussianToDistribution;
using ARM::Export_SABR_FromDistributionToGaussian;
using ARM::Export_SABR_From_Sigma_To_Alpha;
using ARM::Export_BS_EuroBarriere;
using ARM::Export_BS_EuroDoubleBarriere;
using ARM::Export_RealPart_ComplexErf;
using ARM::Export_ImaginaryPart_ComplexErf;
using ARM::Export_GHeston_VanillaOption;
using ARM::Export_LogNormal_TriSpreadOption;
using ARM::Export_LogNormal_TriSpreadDigitalOption;
using ARM::Export_LogNormal_TriSpreadDigitalOption2;
using ARM::Export_CEV_VanillaOption;
using ARM::Export_CEV_DoubleBarrierOption;
using ARM::Export_CEV_BarrierOption;
using ARM::Export_BS_PartialTime_Start_SingleBarrier;
using ARM::Export_BS_PartialTime_End_SingleBarrier;
using ARM::Export_BS_SingleBarrier_2Asset;
using ARM::Export_Bivariate;
using ARM::Export_Gamma;
using ARM::Export_IncompleteBeta;
using ARM::Export_InverseIncompleteBeta;
using ARM::Export_Hypergeometric2F1;
using ARM::Export_Lognormal_Asian_VanillaOption;
using ARM::Export_GaussLegendre_Coefficients;
using ARM::GaussLegendre_Coefficients;
using ARM::Export_GaussHermite_Coefficients;
using ARM::GaussHermite_Coefficients;
using ARM::Export_GaussLaguerre_Coefficients;
using ARM::GaussLaguerre_Coefficients;
using ARM::ConvertFXOptionPriceToStrike;
using ARM::CreateARMGPVectorFromVECTOR;
using ARM::std::vector<double>;
using ARM::Export_GHeston_VanillaOption_ModelVector;
using ARM::Export_GHeston_Implicit_Volatility_ModelVector;
using ARM::ARM_GaussianAnalytics;
using ARM::Export_SABR_Model_Calibrate;
using ARM::SABR_ParameterSet;
using ARM::Export_GeneralizedHeston_Model_Calibrate;
using ARM::Export_GeneralizedHeston_Model_Calibrate_NoJump;
using ARM::GeneralizedHeston_ParameterSet;
using ARM::Export_Student_SABR_Power_SpreadOption;
using ARM::Export_Student_SABR_Power_Digital_SpreadOption;
using ARM::mepi_VanillaOption_STOBS;
using ARM::mepi_VanillaOption_STOBS_delta;
using ARM::mepi_VanillaOption_SABR;
using ARM::mepi_VanillaOption_SABR_delta;
using ARM::Export_Heston2b_OptionPrice;

using ARM::ARM_SABR_Eq;


using ARM::Export_Fund_VanillaOption;
using ARM::Export_StochasticVol_LN_Geometric_VanillaOption_with_Reset;
using ARM::Export_StochasticVol_LN_Arithmetic_VanillaOption_with_Reset;
using ARM::Export_StochasticVol_N_VanillaOption;
using ARM::GLambda_ParameterSet;
using ARM::GLambda_CalibrateFromSABR;
using ARM::Export_Student_GLambda_Power_SpreadOption;
using ARM::Export_Lambert_Function;
using ARM::Export_BlackSholesTimeValue;
using ARM::Export_BlackSholesTimeValue_ImplicitVol;
using ARM::Export_ShiftedLogNormal_Quantile;
using ARM::Export_SABR_Quantile;
using ARM::Export_Student_GLambda_Power_Digital_SpreadOption;
using ARM::Export_Student_GLambda_Power_Index1Digital_SpreadOption;
using ARM::Export_Student_GLambda_Power_Index2Digital_SpreadOption;
using ARM::Export_Student_GLambda_SpreadOption;
using ARM::NormalCDF;
using ARM::NormalCDFInverse;
using ARM::Export_GLambda_Distribution;
using ARM::Export_GLambda_Quantile;
using ARM::Export_Student_Quantile;
using ARM::Export_Student_Distribution;
using ARM::Export_ImcompleteBeta_Inverse;
using ARM::Export_Student_QIntegral;
using ARM::Export_Normal_ImpliedVol;
using ARM::Export_Normal_Digital_ImpliedVol;
using ARM::Export_ShiftedLogNormal_Distribution;
using ARM::Export_SABR_Distribution;
using ARM::Export_Bessel_Y;
using ARM::Export_Bessel_I;
using ARM::Export_Bessel_J;
using ARM::Export_Bessel_K;
using ARM::Export_Hypergeometric_Whittaker_M;
using ARM::Export_Hypergeometric_Whittaker_W;
using ARM::Export_Shifted_Heston_VanillaOption;
using ARM::Export_SABR_Heston_VanillaOption;
using ARM::Export_GLambda_CompleteSpreadoption;
using ARM::Export_Mepi_EDP_VanillaOption;
using ARM::Export_Util_Trigonal_Solve;
using ARM::Export_BiSABR_SpreadOption;
using ARM::Export_Hypergeometric_Appell;
using ARM::Export_Normal_DoubleDigital;
using ARM::Export_Eigenvalues4;
using ARM::Export_Eigenvalues3;
using ARM::Export_BiSABR_CalibrateToSmile;
using ARM::BiSABR_ParameterSet;
using ARM::Export_BlackScholesDigitalOption;
using ARM::Export_BiSABR_CorrelationEvolution;
using ARM::Export_LN_RatioOption;
using ARM::Export_LN_ProductOption;
using ARM::Export_BS_EuroBarriere_ImpliedVol;
using ARM::Export_BS_EuroDoubleBarriere_ImpliedVol;
using ARM::Export_GaussianSABRDigitalCall;
using ARM::Export_GaussianSABRDigitalCallPayingS1;
using ARM::Export_GaussianSABRDigitalCallPayingS2;
using ARM::Export_GaussianSABRDigitalCallPayingS3;

using ARM::ARM_TarnProxy;
using ARM::ARM_VBMinMaxProxy;
using ARM::ARM_Berm2DatesProxy;

using ARM::ARM_SpreadVBProxy;

using ARM::ARM_DensityFunctor;
using ARM::ARM_DensityFunctorPtr;
using ARM::ARM_DensityFunctorPtrVector;
using ARM::Export_BiSABR_Digital_SpreadOption;
using ARM::Export_BiSABR_Digital_SpreadOption_PayS1;
using ARM::Export_BiSABR_Digital_SpreadOption_PayS2;
using ARM::Export_BiSABR_Digital_SpreadOption_PayS3;
using ARM::Export_BiSABR_Distribution;
using ARM::Export_BiSABR_Quantile;
using ARM::Export_Shifted2LogNormal_Distribution;
using ARM::Export_Shifted2LogNormal_Quantile;
using ARM::Export_SABR_BetaEqualZero_Option;
using ARM::Export_BiSABR_SpreadOption;
using ARM::Export_BiSABR_S3_SpreadOption;
using ARM::Export_Heston_OptionPrice;
using ARM::Export_MixteHeston_OptionPrice;
using ARM::Export_Normal_Heston_VanillaOption;
using ARM::Export_TriSABR_VanillaOption;
using ARM::Export_TriSABR_Eigenvalues;
using ARM::Export_Nonparametric_CompleteSpreadoption;
using ARM::Export_NonParametric_LogVolatility;
using ARM::Export_NonParametric_NormalVolatility;
using ARM::Export_NonParametric_N_Distribution;
using ARM::Export_NonParametric_LN_Distribution;
using ARM::Export_NonParametric_N_Quantile;
using ARM::Export_NonParametric_LN_Quantile;


using ARM::SABR_CalibrateToSmileBetaFixedToOne;
using ARM::ARM_SABRToHestonSmileCalibration;
using ARM::ARM_CIRModelParams;
using ARM::ARM_ModelParam;
using ARM::ARM_HestonOptionPricer;

using ARM::ARM_SmileCalibration_Params;
using ARM::ARM_SmileCalibration_Params_SABR;
using ARM::ARM_SmileCalibration_Params_SABR2beta;
using ARM::ARM_SmileCalibration_Params_Heston;
using ARM::ARM_SmileCalibration_Params_Heston2b;
using ARM::ARM_SmileCalibration_Params_BiSABR;
using ARM::ARM_SmileCalibration_Params_Spread2Heston;
using ARM::ARM_SmileCalibration_Params_Merton;
using ARM::ARM_SmileCalibration;
using ARM::ARM_SmileCalibration_SABR;
using ARM::ARM_SmileCalibration_SABR2beta;
using ARM::ARM_SmileCalibration_Heston;
using ARM::ARM_SmileCalibration_Heston2b;
using ARM::ARM_SmileCalibration_NormalHeston;
using ARM::ARM_SmileCalibration_BiSABR;
using ARM::ARM_SmileCalibration_Spread2Heston;
using ARM::ARM_SmileCalibration_Merton;

using ARM::Spread2HestonVanilla;
using ARM::ARM_ArgConv_SABR_ImplicitVol_Formula_Extended_Flag;
using ARM::ARM_Spread2Heston_TOTEMCalibration;
using ARM::MertonOption;
using ARM::ARM_ImpliedVolBS;
using ARM::SuperNormalHeston;

using ARM::ARM_GP_Matrix;

using ARM::Export_BiShiftedHeston_VanillaOption;

using ARM::DistribSet;

using ARM::ARM_ExpRiccati;

#include "gpbase/ostringstream.h"
#include "gpbase/stringmanip.h"
using ARM::stringTrim;
using ARM::stringToUpper;

#include <deque>

/////////////////////////////////////////////////////////////////////////////////////////////////////////
///
///        Vanilla Option Normal
///
/////////////////////////////////////////////////////////////////////////////////////////////////////////

extern long ARMLOCAL_VanillaOption_Normal(
	double underling,
	double volatility,
	double  strike,
	double maturity,
	int callput,
	ARM_result& result
	)
{
    	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	try
	{
        double dResult;
            dResult = Export_normal_VanillaOption(underling,
                strike,
                volatility,
                maturity,
                callput);
		result.setDouble(dResult);
		return ARM_OK;
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}
}

extern long ARMLOCAL_VanillaOption_Normal_Der(
	double i,
	double underling,
	double volatility,
	double strike,
	double maturity,
	int callput,
	ARM_result& result
	)
{
    	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	try
	{
        double dResult;
            dResult = Export_normal_VanillaOption(i,underling,
				strike,
				volatility,
                maturity,
                callput);
		result.setDouble(dResult);
		return ARM_OK;
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}
}

extern long ARMLOCAL_VanillaOption_Normal_Der2(
	double i,
	double j,
	double underling,
	double volatility,
	double strike ,
	double maturity,
	int callput,
	ARM_result& result
	)
{
    	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	try
	{
        double dResult;
            dResult = Export_normal_VanillaOption(i,j,underling,
                strike,
				volatility,
                maturity,
                callput);

		result.setDouble(dResult);
		return ARM_OK;
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}
}


/// double digitale
extern long ARMLOCAL_DoubleDigital_Normal(
	   double fwd1, 
	   double fwd2,
	   double maturity,
	   double K1, double spread1,
	   double K2, double spread2,
	   double vol1plus, double vol1minus,
	   double vol2plus, double vol2minus,
	   double correl,
	   int callorput1,
	   int callorput2,
	   ARM_result& result)
{
    	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	try
	{
        double dResult;
            dResult = Export_Normal_DoubleDigital(
									fwd1, 
									fwd2, 
									maturity, 
									K1, spread1, 
									K2, spread2, 
									vol1plus, vol1minus, 
									vol2plus, vol2minus, 
									correl, 
									callorput1, 
									callorput2);
		result.setDouble(dResult);
		return ARM_OK;
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}
}



/////////////////////////////////////////////////////////////////////////////////////////////////////////
///
///        Lognormal Spreadoption
///
/////////////////////////////////////////////////////////////////////////////////////////////////////////

extern long ARMLOCAL_LogNormal_SpreadOption(
	double S1,
	double S2,
	double sig1,
	double sig2,
	double rho,
	double k,
	double t,
	int callput,
	int optiontype,
	int n,
	ARM_result& result
	)
{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	try
	{
		double dResult = Export_LogNormal_SpreadOption(S1,S2,sig1,sig2,rho,k,t,callput,optiontype,n);
		result.setDouble(dResult);
		return ARM_OK;
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}
}

extern long ARMLOCAL_Smiled_LogNormal_SpreadOption(
	double S1,
	double S2,
	double sig1,
	double sig2,
	double rho,
	double k,
	double t,
	double slope1,
	double slope2,
	int callput,
	int optiontype,
	int n,
	ARM_result& result
	)
{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	try
	{
		double dResult = Export_Smiled_LogNormal_SpreadOption(S1,S2,sig1,sig2,rho,k,t,
												slope1,slope2,callput,optiontype,n);
		result.setDouble(dResult);
		return ARM_OK;
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}
}

extern long ARMLOCAL_LogNormal_SpreadOption_Calibrate_Correlation(
	double S1,
	double S2,
	double sig1,
	double sig2,
	double rho,
	double k,
	double t,
	int callput,
	int optiontype,
	int n,
	ARM_result& result
	)
{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	try
	{
		double dResult = Export_LogNormal_SpreadOption_Calibrate_Correlation(S1,S2,sig1,sig2,rho,k,t,callput,optiontype,n);
		result.setDouble(dResult);
		return ARM_OK;
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}
}

extern long ARMLOCAL_Smiled_LogNormal_SpreadOption_Calibrate_Correlation(
	double S1,
	double S2,
	double sig1,
	double sig2,
	double rho,
	double k,
	double t,
	double slope1,
	double slope2,
	int callput,
	int optiontype,
	int n,
	ARM_result& result
	)
{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	try
	{
		double dResult = Export_Smiled_LogNormal_SpreadOption_Calibrate_Correlation(S1,S2,sig1,sig2,rho,k,t,
												slope1,slope2,callput,optiontype,n);
		result.setDouble(dResult);
		return ARM_OK;
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}
}



extern long ARMLOCAL_LogNormal_SpreadOption_Der(
	int i,
	double S1,
	double S2,
	double sig1,
	double sig2,
	double rho,
	double k,
	double t,
	int callput,
	int optiontype,
	int n,
	ARM_result& result
	)
{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	try
	{
		double dResult = Export_LogNormal_SpreadOption(i,S1,S2,sig1,sig2,rho,k,t,
												callput,optiontype,n);
		result.setDouble(dResult);
		return ARM_OK;
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}
}

extern long ARMLOCAL_LogNormal_SpreadOption_Der2(
	int i,
	int j,
	double S1,
	double S2,
	double sig1,
	double sig2,
	double rho,
	double k,
	double t,
	int callput,
	int optiontype,
	int n,
	ARM_result& result
	)
{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	try
	{
		double dResult = Export_LogNormal_SpreadOption(i,j,S1,S2,sig1,sig2,rho,k,t,
												callput,optiontype,n);
		result.setDouble(dResult);
		return ARM_OK;
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}
}


extern long ARMLOCAL_Smiled_LogNormal_SpreadOption_Der(
	int i,
	double S1,
	double S2,
	double sig1,
	double sig2,
	double rho,
	double k,
	double t,
	double slope1,
	double slope2,
	int callput,
	int optiontype,
	int n,
	ARM_result& result
	)
{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	try
	{
		double dResult = Export_Smiled_LogNormal_SpreadOption(i,S1,S2,sig1,sig2,rho,k,t,
												slope1,slope2,callput,optiontype,n);
		result.setDouble(dResult);
		return ARM_OK;
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}
}

extern long ARMLOCAL_Smiled_LogNormal_SpreadOption_Der2(
	int i,
	int j,
	double S1,
	double S2,
	double sig1,
	double sig2,
	double rho,
	double k,
	double t,
	double slope1,
	double slope2,
	int callput,
	int optiontype,
	int n,
	ARM_result& result
	)
{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	try
	{
		double dResult = Export_Smiled_LogNormal_SpreadOption(i,j,S1,S2,sig1,sig2,rho,k,t,
												slope1,slope2,callput,optiontype,n);
		result.setDouble(dResult);
		return ARM_OK;
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
///
///        Normal Spreadoption
///
/////////////////////////////////////////////////////////////////////////////////////////////////////////

extern long ARMLOCAL_Normal_SpreadOption(
	double S1,
	double S2,
	double sig1,
	double sig2,
	double rho,
	double k,
	double t,
	int callput,
	int optiontype,
	ARM_result& result
	)
	{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	try
	{
		double dResult = Export_Normal_SpreadOption(S1,S2,sig1,sig2,rho,k,t,
												callput,optiontype);
		result.setDouble(dResult);
		return ARM_OK;
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}
}

extern long ARMLOCAL_Smiled_Normal_SpreadOption(
	double S1,
	double S2,
	double sig1,
	double sig2,
	double rho,
	double k,
	double t,
	double slope1,
	double slope2,
	int callput,
	int optiontype,
	ARM_result& result
	)
	{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	try
	{
		double dResult = Export_Smiled_Normal_SpreadOption(S1,S2,sig1,sig2,rho,k,t,
												slope1,slope2,callput,optiontype);
		result.setDouble(dResult);
		return ARM_OK;
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}
}


extern long ARMLOCAL_Normal_SpreadOption_Der(
	int i,
	double S1,
	double S2,
	double sig1,
	double sig2,
	double rho,
	double k,
	double t,
	int callput,
	int optiontype,
	ARM_result& result
	)
	{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	try
	{
		double dResult = Export_Normal_SpreadOption(i,S1,S2,sig1,sig2,rho,k,t,
											callput,optiontype);
		result.setDouble(dResult);
		return ARM_OK;
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}
}

extern long ARMLOCAL_Smiled_Normal_SpreadOption_Der(
	int i,
	double S1,
	double S2,
	double sig1,
	double sig2,
	double rho,
	double k,
	double t,
	double slope1,
	double slope2,
	int callput,
	int optiontype,
	ARM_result& result
	)
	{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	try
	{
		double dResult = Export_Smiled_Normal_SpreadOption(i,S1,S2,sig1,sig2,rho,k,t,
												slope1,slope2,callput,optiontype);
		result.setDouble(dResult);
		return ARM_OK;
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}
}

extern long ARMLOCAL_Normal_SpreadOption_Der2(
	int i,
	int j,
	double S1,
	double S2,
	double sig1,
	double sig2,
	double rho,
	double k,
	double t,
	int callput,
	int optiontype,
	ARM_result& result
	)
	{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	try
	{
		double dResult = Export_Normal_SpreadOption(i,j,S1,S2,sig1,sig2,rho,k,t,
												callput,optiontype);
		result.setDouble(dResult);
		return ARM_OK;
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}
}

extern long ARMLOCAL_Smiled_Normal_SpreadOption_Der2(
	int i,
	int j,
	double S1,
	double S2,
	double sig1,
	double sig2,
	double rho,
	double k,
	double t,
	double slope1,
	double slope2,
	int callput,
	int optiontype,
	ARM_result& result
	)
	{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	try
	{
		double dResult = Export_Smiled_Normal_SpreadOption(i,j,S1,S2,sig1,sig2,rho,k,t,
												slope1,slope2,callput,optiontype);
		result.setDouble(dResult);
		return ARM_OK;
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
///
///        SABR Power Spreadoption
///
/////////////////////////////////////////////////////////////////////////////////////////////////////////


extern long ARMLOCAL_Gaussian_SABR_Power_SpreadOption(
	double S1,
	double S2,
	double alpha1,
	double beta1,
	double rho1,
	double nu1,
	double alpha2,
	double beta2, 
	double rho2, 
	double nu2,
	double copula_corr,
	double t,
	int flag,
	double a10,
	double b10,
	double k10,
	double a20,
	double b20,
	double k20,
	const vector<double>&  C_Parameters_vec,
	ARM_result& result
	)
	{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	std::vector<double>&		Parameters_vec		= NULL;	
	try
	{
		Parameters_vec		= CreateARMGPVectorFromVECTOR( C_Parameters_vec	);
		const vector<double> params_Vec=Parameters_vec->GetValues();
		int n=params_Vec[0];
		double alpha_exp=params_Vec[1];
		double alpha_tanh=params_Vec[2];
		double kb_tanh=params_Vec[3];

		double dResult = Export_Gaussian_SABR_Power_SpreadOption(S1,S2,
													alpha1,beta1,rho1,nu1,
													alpha2,beta2,rho2,nu2,
													copula_corr,t,flag,
												a10,b10,k10,a20,b20,k20,n,alpha_exp,alpha_tanh,kb_tanh);
		result.setDouble(dResult);
		return ARM_OK;
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}
}

extern long ARMLOCAL_Gaussian_SABR_Power_Digital_SpreadOption(
	double S1,
	double S2,
	double alpha1,
	double beta1,
	double rho1,
	double nu1,
	double alpha2,
	double beta2, 
	double rho2, 
	double nu2,
	double copula_corr,
	double t,
	int flag,
	double k10,
	double a20,
	double b20,
	double k20,
	const vector<double>&  C_Parameters_vec,
	ARM_result& result
	)
	{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	std::vector<double>&		Parameters_vec		= NULL;	
	try
	{
		Parameters_vec		= CreateARMGPVectorFromVECTOR( C_Parameters_vec	);
		const vector<double> params_Vec=Parameters_vec->GetValues();
		int n=params_Vec[0];
		double alpha_exp=params_Vec[1];
		double alpha_tanh=params_Vec[2];
		double kb_tanh=params_Vec[3];
		double dResult = Export_Gaussian_SABR_Power_Digital_SpreadOption(S1,S2,
													alpha1,beta1,rho1,nu1,
													alpha2,beta2,rho2,nu2,
													copula_corr,t,flag,
												k10,a20,b20,k20,n,alpha_exp,alpha_tanh,kb_tanh);
		result.setDouble(dResult);
		return ARM_OK;
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}
}



extern long ARMLOCAL_Gaussian_SABR_Power_SpreadOption_Der(
	int i,
	double S1,
	double S2,
	double alpha1,
	double beta1,
	double rho1,
	double nu1,
	double alpha2,
	double beta2, 
	double rho2, 
	double nu2,
	double copula_corr,
	double t,
	int flag,
	double a10,
	double b10,
	double k10,
	double a20,
	double b20,
	double k20,
	const vector<double>&  C_Parameters_vec,
	ARM_result& result
	)
		{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	std::vector<double>&		Parameters_vec		= NULL;	
	try
	{
		Parameters_vec		= CreateARMGPVectorFromVECTOR( C_Parameters_vec	);
		const vector<double> params_Vec=Parameters_vec->GetValues();
		int n=params_Vec[0];
		double alpha_exp=params_Vec[1];
		double alpha_tanh=params_Vec[2];
		double kb_tanh=params_Vec[3];
		double dResult = Export_Gaussian_SABR_Power_SpreadOption(i,S1,S2,
													alpha1,beta1,rho1,nu1,
													alpha2,beta2,rho2,nu2,
													copula_corr,t,flag,
												a10,b10,k10,a20,b20,k20,n,alpha_exp,alpha_tanh,kb_tanh);
		result.setDouble(dResult);
		return ARM_OK;
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}
}

extern long ARMLOCAL_Gaussian_SABR_Power_Digital_SpreadOption_Der(
	int i,
	double S1,
	double S2,
	double alpha1,
	double beta1,
	double rho1,
	double nu1,
	double alpha2,
	double beta2, 
	double rho2, 
	double nu2,
	double copula_corr,
	double t,
	int flag,
	double k10,
	double a20,
	double b20,
	double k20,
	const vector<double>&  C_Parameters_vec,
	ARM_result& result
	)
		{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	std::vector<double>&		Parameters_vec		= NULL;	
	try
	{
		Parameters_vec		= CreateARMGPVectorFromVECTOR( C_Parameters_vec	);
		const vector<double> params_Vec=Parameters_vec->GetValues();
		int n=params_Vec[0];
		double alpha_exp=params_Vec[1];
		double alpha_tanh=params_Vec[2];
		double kb_tanh=params_Vec[3];
		double dResult = Export_Gaussian_SABR_Power_Digital_SpreadOption(i,S1,S2,
													alpha1,beta1,rho1,nu1,
													alpha2,beta2,rho2,nu2,
													copula_corr,t,flag,
												k10,a20,b20,k20,n,alpha_exp,alpha_tanh,kb_tanh);
		result.setDouble(dResult);
		return ARM_OK;
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}
}


extern long ARMLOCAL_Gaussian_SABR_Power_SpreadOption_Der2(
	int i,
	int j,
	double S1,
	double S2,
	double alpha1,
	double beta1,
	double rho1,
	double nu1,
	double alpha2,
	double beta2, 
	double rho2, 
	double nu2,
	double copula_corr,
	double t,
	int flag,
	double a10,
	double b10,
	double k10,
	double a20,
	double b20,
	double k20,
	const vector<double>&  C_Parameters_vec,
	ARM_result& result
	)
		{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	std::vector<double>&		Parameters_vec		= NULL;	
	try
	{
		Parameters_vec		= CreateARMGPVectorFromVECTOR( C_Parameters_vec	);
		const vector<double> params_Vec=Parameters_vec->GetValues();
		int n=params_Vec[0];
		double alpha_exp=params_Vec[1];
		double alpha_tanh=params_Vec[2];
		double kb_tanh=params_Vec[3];
		double dResult = Export_Gaussian_SABR_Power_SpreadOption(i,j,S1,S2,
													alpha1,beta1,rho1,nu1,
													alpha2,beta2,rho2,nu2,
													copula_corr,t,flag,
												a10,b10,k10,a20,b20,k20,n,alpha_exp,alpha_tanh,kb_tanh);
		result.setDouble(dResult);
		return ARM_OK;
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}
}

extern long ARMLOCAL_Gaussian_SABR_Power_Digital_SpreadOption_Der2(
	int i,
	int j,
	double S1,
	double S2,
	double alpha1,
	double beta1,
	double rho1,
	double nu1,
	double alpha2,
	double beta2, 
	double rho2, 
	double nu2,
	double copula_corr,
	double t,
	int flag,
	double k10,
	double a20,
	double b20,
	double k20,
	const vector<double>&  C_Parameters_vec,
	ARM_result& result
	)
		{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	std::vector<double>&		Parameters_vec		= NULL;	
	try
	{
		Parameters_vec		= CreateARMGPVectorFromVECTOR( C_Parameters_vec	);
		const vector<double> params_Vec=Parameters_vec->GetValues();
		int n=params_Vec[0];
		double alpha_exp=params_Vec[1];
		double alpha_tanh=params_Vec[2];
		double kb_tanh=params_Vec[3];
		double dResult = Export_Gaussian_SABR_Power_Digital_SpreadOption(i,j,S1,S2,
													alpha1,beta1,rho1,nu1,
													alpha2,beta2,rho2,nu2,
													copula_corr,t,flag,
												k10,a20,b20,k20,n,alpha_exp,alpha_tanh,kb_tanh);
		result.setDouble(dResult);
		return ARM_OK;
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}
}


extern long ARMLOCAL_Gaussian_SABR_Power_SpreadOption_Certitude(
	double copula_corr,
	double t,
	int n,
	ARM_result& result
	)
		{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	try
	{
		double dResult = Export_Gaussian_SABR_Power_SpreadOption_Certitude(copula_corr,t,n);
		result.setDouble(dResult);
		return ARM_OK;
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}
}



////////////////////////////////////////////////////////////////////////////////////////
///
///			Black and Sholes Formula
///
////////////////////////////////////////////////////////////////////////////////////////


extern long ARMLOCAL_CF_BlackSholes(
								 double forward,
								 double totalvolatility,
								 double bondprice,
								 double strike,
								 double CallPut,
								 ARM_result& result
	)
		{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	try
	{
		double dResult = Export_BlackSholes(forward,totalvolatility,bondprice,strike,CallPut);
		result.setDouble(dResult);
		return ARM_OK;
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}
}

extern long ARMLOCAL_CF_BlackSholes_Der(
									 int i,
									 double forward,
									 double totalvolatility,
									 double bondprice,
									 double strike,
									 double CallPut,
									 ARM_result& result
									 )
{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	try
	{
		double dResult = Export_BlackSholes(i,forward,totalvolatility,bondprice,strike,CallPut);
		result.setDouble(dResult);
		return ARM_OK;
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}
}


extern long ARMLOCAL_CF_BlackSholes_Der2(
									  int i,
									  int j,
									  double forward,
									  double totalvolatility,
									  double bondprice,
									  double strike,
									  double CallPut,
									  ARM_result& result
									  )
		{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	try
	{
		double dResult = Export_BlackSholes(i,j,forward,totalvolatility,bondprice,strike,CallPut);
		result.setDouble(dResult);
		return ARM_OK;
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}
}


extern long ARMLOCAL_CF_BlackSholes_ImplicitVol(
												double forward,
												double bondprice,
												double strike,
												double CallPut,
												double optprice,
												int algo,
												ARM_result& result
												)
{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	
	try
	{
		double dResult = Export_BlackSholes_ImplicitVol(forward,bondprice,strike,CallPut,optprice,algo);
		result.setDouble(dResult);
		return ARM_OK;
	}
	
	catch(Exception& x)
	{
		x.DebugPrint();
		
		ARM_RESULT();
	}
}

extern long ARMLOCAL_CF_BlackSholes_ImplicitVolatility(
												double forward,
												double bondprice,
												double strike,
												double maturity,
												double CallPut,
												double optprice,
												int algo,
												ARM_result& result
												)
{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	
	try
	{
		double dResult = Export_BlackSholes_ImplicitVol(forward,bondprice,strike,CallPut,optprice,algo);
		result.setDouble(dResult/sqrt(maturity));
		return ARM_OK;
	}
	
	catch(Exception& x)
	{
		x.DebugPrint();
		
		ARM_RESULT();
	}
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////
///
///        ShiftedLognormal Power Spreadoption
///
/////////////////////////////////////////////////////////////////////////////////////////////////////////



extern long ARMLOCAL_Gaussian_ShiftedLN_Power_SpreadOption(
	double S1,
	double S2,
	double sigma1,
	double alpha1,
	double sigma2,
	double alpha2,
	double copula_corr,
	double t,
	double a10,
	double b10,
	double k10,
	double a20,
	double b20,
	double k20,
	int n,
	ARM_result& result
	)
	{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	try
	{
		double dResult = Export_Gaussian_ShiftedLN_Power_SpreadOption(S1,S2,
													sigma1,alpha1,
													sigma2,alpha2,
													copula_corr,t,
												a10,b10,k10,a20,b20,k20,n);
		result.setDouble(dResult);
		return ARM_OK;
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}
}




extern long ARMLOCAL_Gaussian_ShiftedLN_Power_SpreadOption_Der(
	int i,
	double S1,
	double S2,
	double sigma1,
	double alpha1,
	double sigma2,
	double alpha2,
	double copula_corr,
	double t,
	double a10,
	double b10,
	double k10,
	double a20,
	double b20,
	double k20,
	int n,
	ARM_result& result
	)
		{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	try
	{
		double dResult = Export_Gaussian_ShiftedLN_Power_SpreadOption(i,S1,S2,
													sigma1,alpha1,
													sigma2,alpha2,
													copula_corr,t,
												a10,b10,k10,a20,b20,k20,n);
		result.setDouble(dResult);
		return ARM_OK;
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}
}


extern long ARMLOCAL_Gaussian_ShiftedLN_Power_SpreadOption_Der2(
	int i,
	int j,
	double S1,
	double S2,
	double sigma1,
	double alpha1,
	double sigma2,
	double alpha2,
	double copula_corr,
	double t,
	double a10,
	double b10,
	double k10,
	double a20,
	double b20,
	double k20,
	int n,
	ARM_result& result
	)
		{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	try
	{
		double dResult = Export_Gaussian_ShiftedLN_Power_SpreadOption(i,j,S1,S2,
													sigma1,alpha1,
													sigma2,alpha2,
													copula_corr,t,
												a10,b10,k10,a20,b20,k20,n);
		result.setDouble(dResult);
		return ARM_OK;
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}
}


extern long ARMLOCAL_Gaussian_ShiftedLN_Power_SpreadOption_Certitude(
	double copula_corr,
	double t,
	int n,
	ARM_result& result
	)
{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	
	try
	{
		double dResult = Export_Gaussian_ShiftedLN_Power_SpreadOption_Certitude(copula_corr,t,n);
		result.setDouble(dResult);
		return ARM_OK;
	}
	
	catch(Exception& x)
	{
		x.DebugPrint();
		
		ARM_RESULT();
	}
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
///
///        Merton Formula
///
/////////////////////////////////////////////////////////////////////////////////////////////////////////



extern long ARMLOCAL_Merton_JumpDiffusion(
		double F,
		double K,
		double t,
		double sigma,
		double lambda,
		double muJ, 
		double sigmaJ,
		int callorput,int nb,
		ARM_result& result
		)
{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	
	try
	{
		double dResult = Export_Merton_JumpDiffusion(F,K,t,sigma,lambda,muJ,sigmaJ,callorput, nb);
		result.setDouble(dResult);
		return ARM_OK;
	}
	
	catch(Exception& x)
	{
		x.DebugPrint();
		
		ARM_RESULT();
	}
}




extern long ARMLOCAL_Merton_JumpDiffusion_Der(
		int i,
		double F,
		double K,
		double t,
		double sigma,
		double lambda,
		double muJ, 
		double sigmaJ,
		int callorput,
		int nb,
		ARM_result& result
	)
	{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	
	try
	{
		double dResult = Export_Merton_JumpDiffusion(i,F,K,t,sigma,lambda,muJ,sigmaJ,callorput, nb);
		result.setDouble(dResult);
		return ARM_OK;
	}
	
	catch(Exception& x)
	{
		x.DebugPrint();
		
		ARM_RESULT();
	}
}



extern long ARMLOCAL_Merton_JumpDiffusion_Der2(
		int i,
		int j,
		double F,
		double K,
		double t,
		double sigma,
		double lambda,
		double muJ, 
		double sigmaJ,
		int callorput,
		int nb,
		ARM_result& result
	)
	{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	
	try
	{
		double dResult = Export_Merton_JumpDiffusion(i,j,F,K,t,sigma,lambda,muJ,sigmaJ,callorput, nb);
		result.setDouble(dResult);
		return ARM_OK;
	}
	
	catch(Exception& x)
	{
		x.DebugPrint();
		
		ARM_RESULT();
	}
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////
///
///        SABR Implicit Vol Formula
///
/////////////////////////////////////////////////////////////////////////////////////////////////////////



extern long ARMLOCAL_SABR_ImplicitVol(
									  double f, 
									  double K, 
									  double tex,
									  double alpha, 
									  double beta, 
									  double rho, 
									  double nu,
									  int flag,
									  int nbsteps,
									  double C_Alpha_Exp,
									  double C_Alpha_Tanh,
									  double C_Kb_Tanh,
									  ARM_result& result
									  )
{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	
	try
	{
		double dResult = Export_SABR_ImplicitVol(f,K,tex,alpha,beta,rho,nu,flag, nbsteps,C_Alpha_Exp,C_Alpha_Tanh,C_Kb_Tanh);
		result.setDouble(dResult);
		return ARM_OK;
	}
	
	catch(Exception& x)
	{
		x.DebugPrint();
		
		ARM_RESULT();
	}
}




extern long ARMLOCAL_SABR_ImplicitVol_Der(
										  int i,
										  double f, 
										  double K, 
										  double tex,
										  double alpha, 
										  double beta, 
										  double rho, 
										  double nu,
										  int flag,
										  int nbsteps,
									  double C_Alpha_Exp,
									  double C_Alpha_Tanh,
									  double C_Kb_Tanh,
										  ARM_result& result
										  )
{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	
	try
	{
		double dResult = Export_SABR_ImplicitVol(i,f,K,tex,alpha,beta,rho,nu,flag, nbsteps,C_Alpha_Exp,C_Alpha_Tanh,C_Kb_Tanh);
		result.setDouble(dResult);
		return ARM_OK;
	}
	
	catch(Exception& x)
	{
		x.DebugPrint();
		
		ARM_RESULT();
	}
}



extern long ARMLOCAL_SABR_ImplicitVol_Der2(
										   int i,
										   int j,
										   double f, 
										   double K, 
										   double tex,
										   double alpha, 
										   double beta, 
										   double rho, 
										   double nu,
										   int flag,
										   int nbsteps,
									  double C_Alpha_Exp,
									  double C_Alpha_Tanh,
									  double C_Kb_Tanh,
										   ARM_result& result
										   )
{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	
	try
	{
		double dResult = Export_SABR_ImplicitVol(i,j,f,K,tex,alpha,beta,rho,nu,flag, nbsteps,C_Alpha_Exp,C_Alpha_Tanh,C_Kb_Tanh);
		result.setDouble(dResult);
		return ARM_OK;
	}
	
	catch(Exception& x)
	{
		x.DebugPrint();
		
		ARM_RESULT();
	}
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////
///
///        SABR Vanilla option Formula
///
/////////////////////////////////////////////////////////////////////////////////////////////////////////



extern long ARMLOCAL_SABR_VanillaOption(
									  double f, 
									  double K, 
									  double tex,
									  double alpha, 
									  double beta, 
									  double rho, 
									  double nu,
									  int callput,
									  int flag,
									  int nbsteps,
									  double C_Alpha_Exp,
									  double C_Alpha_Tanh,
									  double C_Kb_Tanh,
									  ARM_result& result
									  )
{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	
	try
	{
		double dResult = Export_SABR_VanillaOption(f,K,tex,alpha,beta,rho,nu,callput,flag, nbsteps,C_Alpha_Exp,C_Alpha_Tanh,C_Kb_Tanh);
		result.setDouble(dResult);
		return ARM_OK;
	}
	
	catch(Exception& x)
	{
		x.DebugPrint();
		
		ARM_RESULT();
	}
}

extern long ARMLOCAL_SABR_FromGaussianToDistribution(
									  double f, 
									  double K, 
									  double tex,
									  double alpha, 
									  double beta, 
									  double rho, 
									  double nu,
									  int flag,
									  int nbsteps,
									  double aexp,
									  double atanh,
									  double ktanh,
									  ARM_result& result
									  )
{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	
	try
	{
		double dResult = Export_SABR_FromGaussianToDistribution(f,K,tex,alpha,beta,rho,nu,flag, nbsteps,aexp,atanh,ktanh);
		result.setDouble(dResult);
		return ARM_OK;
	}
	
	catch(Exception& x)
	{
		x.DebugPrint();
		
		ARM_RESULT();
	}
}

extern long ARMLOCAL_SABR_FromDistributionToGaussian(
									  double f, 
									  double K, 
									  double tex,
									  double alpha, 
									  double beta, 
									  double rho, 
									  double nu,
									  int flag,
									  int nbsteps,
									  double aexp,
									  double atanh,
									  double ktanh,
									  ARM_result& result
									  )
{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	
	try
	{
		double dResult = Export_SABR_FromDistributionToGaussian(f,K,tex,alpha,beta,rho,nu,flag, nbsteps,aexp,atanh,ktanh);
		result.setDouble(dResult);
		return ARM_OK;
	}
	
	catch(Exception& x)
	{
		x.DebugPrint();
		
		ARM_RESULT();
	}
}



extern long ARMLOCAL_SABR_VanillaOption_Der(
										  int i,
										  double f, 
										  double K, 
										  double tex,
										  double alpha, 
										  double beta, 
										  double rho, 
										  double nu,
										  int callput,
										  int flag,
										  int nbsteps,
									  double C_Alpha_Exp,
									  double C_Alpha_Tanh,
									  double C_Kb_Tanh,
										  ARM_result& result
										  )
{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	
	try
	{
		double dResult = Export_SABR_VanillaOption(i,f,K,tex,alpha,beta,rho,nu,callput,flag, nbsteps,C_Alpha_Exp,C_Alpha_Tanh,C_Kb_Tanh);
		result.setDouble(dResult);
		return ARM_OK;
	}
	
	catch(Exception& x)
	{
		x.DebugPrint();
		
		ARM_RESULT();
	}
}



extern long ARMLOCAL_SABR_VanillaOption_Der2(
										   int i,
										   int j,
										   double f, 
										   double K, 
										   double tex,
										   double alpha, 
										   double beta, 
										   double rho, 
										   double nu,
										   int callput,
										   int flag,
										   int nbsteps,
									  double C_Alpha_Exp,
									  double C_Alpha_Tanh,
									  double C_Kb_Tanh,
										   ARM_result& result
										   )
{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	
	try
	{
		double dResult = Export_SABR_VanillaOption(i,j,f,K,tex,alpha,beta,rho,nu,callput,flag, nbsteps,C_Alpha_Exp,C_Alpha_Tanh,C_Kb_Tanh);
		result.setDouble(dResult);
		return ARM_OK;
	}
	
	catch(Exception& x)
	{
		x.DebugPrint();
		
		ARM_RESULT();
	}
}


extern long ARMLOCAL_SABR_From_Sigma_To_Alpha (
							   double f, 
							   double K, 
							   double tex,
							   double sigma, 
							   double beta, 
							   double rho, 
							   double nu,
							   int flag,
							   int nbsteps,
							   ARM_result& result)
							   {
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	
	try
	{
		double dResult =Export_SABR_From_Sigma_To_Alpha(f,K, tex, sigma,beta, rho, nu,flag,nbsteps);
		result.setDouble(dResult);
		return ARM_OK;
	}
	
	catch(Exception& x)
	{
		x.DebugPrint();
		
		ARM_RESULT();
	}
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
///
///        Barriere Option in BS Model Formula (Single and double barriere)
///
/////////////////////////////////////////////////////////////////////////////////////////////////////////



extern long ARMLOCAL_BS_EuroBarriere(
									  double f, 
									  double k, 
									  double b,
									  double r, 
									  double v, 
									  double t, 
									  double discount,
									  int callput,
									  int inout,
									  int updown,
									  ARM_result& result
									  )
{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	
	try
	{
		double rate;
		if(t>0)  rate=-log(discount)/t;else rate= 0.0;
		double dResult = Export_BS_EuroBarriere(f,k,b,r,v,t,rate,callput,inout, updown);
		result.setDouble(dResult);
		return ARM_OK;
	}
	
	catch(Exception& x)
	{
		x.DebugPrint();
		
		ARM_RESULT();
	}
}


extern long ARMLOCAL_BS_EuroDoubleBarriere(
									  double f, 
									  double k, 
									  double bup,
									  double bdown, 
									  double v, 
									  double t, 
									  double r,
									  double b,
									  int callput,
									  ARM_result& result
									  )
{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	
	try
	{
		double dResult = Export_BS_EuroDoubleBarriere(f,k,bup,bdown,v,t,r,b,callput);
		result.setDouble(dResult);
		return ARM_OK;
	}
	
	catch(Exception& x)
	{
		x.DebugPrint();
		
		ARM_RESULT();
	}
}

extern long ARMLOCAL_BS_EuroDoubleBarriere_ImpliedVol(
									  double f, 
									  double k, 
									  double bup,
									  double bdown, 
									  double opt, 
									  double t, 
									  double r,
									  double b,
									  int callput,
									  ARM_result& result
									  )
{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	
	try
	{
		double dResult = Export_BS_EuroDoubleBarriere_ImpliedVol(f,k,bup,bdown,opt,t,r,b,callput);
		result.setDouble(dResult);
		return ARM_OK;
	}
	
	catch(Exception& x)
	{
		x.DebugPrint();
		
		ARM_RESULT();
	}
}


extern long ARMLOCAL_BS_EuroBarriere_ImpliedVol(
									  double f, 
									  double k, 
									  double b,
									  double r, 
									  double op, 
									  double t, 
									  double discount,
									  int callput,
									  int inout,
									  int updown,
									  ARM_result& result
									  )
{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	
	try
	{
		double rate;
		if(t>0)  rate=-log(discount)/t;else rate= 0.0;
		double dResult = Export_BS_EuroBarriere_ImpliedVol(f,k,b,r,op,t,rate,callput,inout, updown);
		result.setDouble(dResult);
		return ARM_OK;
	}
	
	catch(Exception& x)
	{
		x.DebugPrint();
		
		ARM_RESULT();
	}
}

extern long ARMLOCAL_BS_PartialTime_Start_SingleBarrier(
									  double f,
									  double k, 
									  double barrier, 
									  double rebate,
									  double v,
									  double bendtime,
									  double t,
									  int callput,
									  int optype,
									  ARM_result& result
									  )
{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	
	try
	{
		double dResult = Export_BS_PartialTime_Start_SingleBarrier( f,k,barrier,rebate,v,bendtime,t,callput,optype);
		result.setDouble(dResult);
		return ARM_OK;
	}
	
	catch(Exception& x)
	{
		x.DebugPrint();
		
		ARM_RESULT();
	}
}

extern long ARMLOCAL_BS_PartialTime_End_SingleBarrier(
									  double f,
									  double k, 
									  double barrier, 
									  double rebate,
									  double v,
									  double bstarttime,
									  double t,
									  int callput,
									  int optype,
									  ARM_result& result
									  )
{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	
	try
	{
		double dResult = Export_BS_PartialTime_End_SingleBarrier( f,k,barrier,rebate,v,bstarttime,t,callput,optype);
		result.setDouble(dResult);
		return ARM_OK;
	}
	
	catch(Exception& x)
	{
		x.DebugPrint();
		
		ARM_RESULT();
	}
}

extern long ARMLOCAL_BS_SingleBarrier_2Asset(
									  double f1,
									  double k1,
									  double f2,
									  double k2, 
									  double v1, 
									  double v2,
									  double corr,
									  double t,
									  int    callput, 
									  int    optype,
									  ARM_result& result
									  )
{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	
	try
	{
		double dResult = Export_BS_SingleBarrier_2Asset(f1,k1,f2,k2,v1,v2,corr,t,callput,optype);
		result.setDouble(dResult);
		return ARM_OK;
	}
	
	catch(Exception& x)
	{
		x.DebugPrint();
		
		ARM_RESULT();
	}
}




extern long ARMLOCAL_BS_EuroBarriere_Der(
										 int i,
										 double f, 
										 double k, 
										 double b,
										 double r, 
										 double v, 
										 double t, 
										 double discount,
										 int callput,
										 int inout,
										 int updown,
										 ARM_result& result
										 )
{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	
	try
	{
		double rate;
		if(t>0)  rate=-log(discount)/t;else rate= 0.0;
		double dResult = Export_BS_EuroBarriere(i,f,k,b,r,v,t,rate,callput,inout, updown);
		result.setDouble(dResult);
		return ARM_OK;
	}
	
	catch(Exception& x)
	{
		x.DebugPrint();
		
		ARM_RESULT();
	}
}

extern long ARMLOCAL_BS_EuroDoubleBarriere_Der(
										 int i,
										 double f, 
										 double k, 
										 double bup,
										 double bdown, 
										 double v, 
										 double t, 
										 int callput,
										 ARM_result& result
										 )
{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	
	try
	{
		double dResult = Export_BS_EuroDoubleBarriere(i,f,k,bup,bdown,v,t,callput);
		result.setDouble(dResult);
		return ARM_OK;
	}
	
	catch(Exception& x)
	{
		x.DebugPrint();
		
		ARM_RESULT();
	}
}

extern long ARMLOCAL_BS_PartialTime_Start_SingleBarrier_Der(
									  int i,
									  double f,
									  double k, 
									  double barrier, 
									  double rebate,
									  double v,
									  double bendtime,
									  double t,
									  int callput,
									  int optype,
									  ARM_result& result
									  )
{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	
	try
	{
		double dResult = Export_BS_PartialTime_Start_SingleBarrier(i, f,k,barrier,rebate,v,bendtime,t,callput,optype);
		result.setDouble(dResult);
		return ARM_OK;
	}
	
	catch(Exception& x)
	{
		x.DebugPrint();
		
		ARM_RESULT();
	}
}

extern long ARMLOCAL_BS_PartialTime_End_SingleBarrier_Der(
									  int i,
									  double f,
									  double k, 
									  double barrier, 
									  double rebate,
									  double v,
									  double bstarttime,
									  double t,
									  int callput,
									  int optype,
									  ARM_result& result
									  )
{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	
	try
	{
		double dResult = Export_BS_PartialTime_End_SingleBarrier(i, f,k,barrier,rebate,v,bstarttime,t,callput,optype);
		result.setDouble(dResult);
		return ARM_OK;
	}
	
	catch(Exception& x)
	{
		x.DebugPrint();
		
		ARM_RESULT();
	}
}

extern long ARMLOCAL_BS_SingleBarrier_2Asset_Der(
									  int i,
									  double f1,
									  double k1,
									  double f2,
									  double k2, 
									  double v1, 
									  double v2,
									  double corr,
									  double t,
									  int    callput, 
									  int    optype,
									  ARM_result& result
									  )
{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	
	try
	{
		double dResult = Export_BS_SingleBarrier_2Asset(i,f1,k1,f2,k2,v1,v2,corr,t,callput,optype);
		result.setDouble(dResult);
		return ARM_OK;
	}
	
	catch(Exception& x)
	{
		x.DebugPrint();
		
		ARM_RESULT();
	}
}




extern long ARMLOCAL_BS_EuroBarriere_Der2(
										 int i,
										 int j,
										 double f, 
										 double k, 
										 double b,
										 double r, 
										 double v, 
										 double t, 
										 double discount,
										 int callput,
										 int inout,
										 int updown,
										 ARM_result& result
										 )
{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	
	try
	{
		double rate;
		if(t>0)  rate=-log(discount)/t;else rate= 0.0;
		double dResult = Export_BS_EuroBarriere(i,j,f,k,b,r,v,t,rate,callput,inout, updown);
		result.setDouble(dResult);
		return ARM_OK;
	}
	
	catch(Exception& x)
	{
		x.DebugPrint();
		
		ARM_RESULT();
	}
}

extern long ARMLOCAL_BS_EuroDoubleBarriere_Der2(
										 int i,
										 int j,
										 double f, 
										 double k, 
										 double bup,
										 double bdown, 
										 double v, 
										 double t, 
										 int callput,
										 ARM_result& result
										 )
{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	
	try
	{
		double dResult = Export_BS_EuroDoubleBarriere(i,j,f,k,bup,bdown,v,t,callput);
		result.setDouble(dResult);
		return ARM_OK;
	}
	
	catch(Exception& x)
	{
		x.DebugPrint();
		
		ARM_RESULT();
	}
}


extern long ARMLOCAL_BS_PartialTime_Start_SingleBarrier_Der2(
									  int i,
									  int j,
									  double f,
									  double k, 
									  double barrier, 
									  double rebate,
									  double v,
									  double bendtime,
									  double t,
									  int callput,
									  int optype,
									  ARM_result& result
									  )
{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	
	try
	{
		double dResult = Export_BS_PartialTime_Start_SingleBarrier(i,j, f,k,barrier,rebate,v,bendtime,t,callput,optype);
		result.setDouble(dResult);
		return ARM_OK;
	}
	
	catch(Exception& x)
	{
		x.DebugPrint();
		
		ARM_RESULT();
	}
}

extern long ARMLOCAL_BS_PartialTime_End_SingleBarrier_Der2(
									  int i,
									  int j,
									  double f,
									  double k, 
									  double barrier, 
									  double rebate,
									  double v,
									  double bstarttime,
									  double t,
									  int callput,
									  int optype,
									  ARM_result& result
									  )
{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	
	try
	{
		double dResult = Export_BS_PartialTime_End_SingleBarrier(i,j, f,k,barrier,rebate,v,bstarttime,t,callput,optype);
		result.setDouble(dResult);
		return ARM_OK;
	}
	
	catch(Exception& x)
	{
		x.DebugPrint();
		
		ARM_RESULT();
	}
}

extern long ARMLOCAL_BS_SingleBarrier_2Asset_Der2(
									  int i,
									  int j,
									  double f1,
									  double k1,
									  double f2,
									  double k2, 
									  double v1, 
									  double v2,
									  double corr,
									  double t,
									  int    callput, 
									  int    optype,
									  ARM_result& result
									  )
{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	
	try
	{
		double dResult = Export_BS_SingleBarrier_2Asset(i,j,f1,k1,f2,k2,v1,v2,corr,t,callput,optype);
		result.setDouble(dResult);
		return ARM_OK;
	}
	
	catch(Exception& x)
	{
		x.DebugPrint();
		
		ARM_RESULT();
	}
}

////////////////////////////////////////////////////////////////////////////////////////////
///
///			Computation of the bivariate function 
///					
////////////////////////////////////////////////////////////////////////////////////////////


extern long ARMLOCAL_Bivariate(
										 double x,
										 double y,
										 double rho,
										 int p,
										 int q,
										 ARM_result& result
										 )
{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	
	try
	{
		double dResult = Export_Bivariate(x,y,rho,p,q);
		result.setDouble(dResult);
		return ARM_OK;
	}
	
	catch(Exception& x)
	{
		x.DebugPrint();
		
		ARM_RESULT();
	}
}
////////////////////////////////////////////////////////////////////////////////////////////
///
///			Computation of the gamma function 
///					
////////////////////////////////////////////////////////////////////////////////////////////


extern long ARMLOCAL_Gamma(
										 double x,
										 ARM_result& result
										 )
{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	
	try
	{
		double dResult = Export_Gamma(x);
		result.setDouble(dResult);
		return ARM_OK;
	}
	
	catch(Exception& x)
	{
		x.DebugPrint();
		
		ARM_RESULT();
	}
}

////////////////////////////////////////////////////////////////////////////////////////////
///
///			Computation of the incomplete beta function 
///					
////////////////////////////////////////////////////////////////////////////////////////////


extern long ARMLOCAL_ImcompleteBeta(	 double x,double y,double z,
										 ARM_result& result
										 )
{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	
	try
	{
		double dResult = Export_IncompleteBeta(x,y,z);
		result.setDouble(dResult);
		return ARM_OK;
	}
	
	catch(Exception& x)
	{
		x.DebugPrint();
		
		ARM_RESULT();
	}
}

extern long ARMLOCAL_InverseImcompleteBeta(	 double x,double y,double z0,double z,
										 ARM_result& result
										 )
{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	
	try
	{
		double dResult = Export_InverseIncompleteBeta(x,y,z0,z);
		result.setDouble(dResult);
		return ARM_OK;
	}
	
	catch(Exception& x)
	{
		x.DebugPrint();
		
		ARM_RESULT();
	}
}

extern long ARMLOCAL_Hypergeometric2F1(	 double x,double y,double z0,double z,
										 ARM_result& result
										 )
{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	
	try
	{
		double dResult = Export_Hypergeometric2F1(x,y,z0,z);
		result.setDouble(dResult);
		return ARM_OK;
	}
	
	catch(Exception& x)
	{
		x.DebugPrint();
		
		ARM_RESULT();
	}
}
////////////////////////////////////////////////////////////////////////////////////////////
///
///			Computation of complex error function erf(z) 
///					
////////////////////////////////////////////////////////////////////////////////////////////


extern long ARMLOCAL_RealPart_ComplexErf(
										 double x,
										 double y,
										 int n,
										 ARM_result& result
										 )
{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	
	try
	{
		double dResult = Export_RealPart_ComplexErf(x,y,n);
		result.setDouble(dResult);
		return ARM_OK;
	}
	
	catch(Exception& x)
	{
		x.DebugPrint();
		
		ARM_RESULT();
	}
}

extern long ARMLOCAL_ImaginaryPart_ComplexErf(
										 double x,
										 double y,
										 int n,
										 ARM_result& result
										 )
{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	
	try
	{
		double dResult = Export_ImaginaryPart_ComplexErf(x,y,n);
		result.setDouble(dResult);
		return ARM_OK;
	}
	
	catch(Exception& x)
	{
		x.DebugPrint();
		
		ARM_RESULT();
	}
}


extern long ARMLOCAL_cdfNormal_Inv(
										 double x,
										 ARM_result& result
										 )
{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	
	try
	{
		double dResult = NormalCDFInverse(x);
		result.setDouble(dResult);
		return ARM_OK;
	}
	
	catch(Exception& x)
	{
		x.DebugPrint();
		
		ARM_RESULT();
	}
}

extern long ARMLOCAL_cdfNormal(
										 double x,
										 ARM_result& result
										 )
{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	
	try
	{
		double dResult = NormalCDF(x);
		result.setDouble(dResult);
		return ARM_OK;
	}
	
	catch(Exception& x)
	{
		x.DebugPrint();
		
		ARM_RESULT();
	}
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
///
///        Generalized Heston Vanilla Option Formula
///
/////////////////////////////////////////////////////////////////////////////////////////////////////////



extern long ARMLOCAL_GHeston_VanillaOption(
										   double F,
										   double K,
										   double sig,
										   double t,
										   double longtermV,
										   double theta,
										   double ksi,
										   double rho,
										   double lambda,
										   double muJ, 
										   double sigmaJ,
										   int callorput,
										   int nb,
										   ARM_result& result
										   )
{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	
	try
	{
		double dResult = Export_GHeston_VanillaOption(F,K,sig,t,longtermV,theta,ksi,rho,lambda, muJ,sigmaJ,callorput,nb);
		result.setDouble(dResult);
		return ARM_OK;
	}
	
	catch(Exception& x)
	{
		x.DebugPrint();
		
		ARM_RESULT();
	}
}




extern long ARMLOCAL_GHeston_VanillaOption_Der(
											   int i,
											   double F,
											   double K,
											   double sig,
											   double t,
											   double longtermV,
											   double theta,
											   double ksi,
											   double rho,
											   double lambda,
											   double muJ, 
											   double sigmaJ,
											   int callorput,
											   int nb,
											   ARM_result& result
											   )
{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	
	try
	{
		double dResult = Export_GHeston_VanillaOption(i,F,K,sig,t,longtermV,theta,ksi,rho,lambda, muJ,sigmaJ,callorput,nb);
		result.setDouble(dResult);
		return ARM_OK;
	}
	
	catch(Exception& x)
	{
		x.DebugPrint();
		
		ARM_RESULT();
	}
}



extern long ARMLOCAL_GHeston_VanillaOption_Der2(
												int i,
												int j,
												double F,
												double sig,
												double K,
												double t,
												double longtermV,
												double theta,
												double ksi,
												double rho,
												double lambda,
												double muJ, 
												double sigmaJ,
												int callorput,
												int nb,
												ARM_result& result
												)
{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	
	try
	{
		double dResult = Export_GHeston_VanillaOption(i,j,F,K,sig,t,longtermV,theta,ksi,rho,lambda,muJ,sigmaJ,callorput,nb);
		result.setDouble(dResult);
		return ARM_OK;
	}
	
	catch(Exception& x)
	{
		x.DebugPrint();
		
		ARM_RESULT();
	}
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////
///
///       Tri Spread Option Formula
///
/////////////////////////////////////////////////////////////////////////////////////////////////////////



extern long ARMLOCAL_LogNormal_TriSpreadOption(
										   double S1,double S2,double S3,double sig1,double sig2,double sig3,
										double rho12,double rho13,double rho23,double mu1,double mu2,double mu3,
										double a0,double a1,double a2,double a3,double t,int n,
										   ARM_result& result
										   )
{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	
	try
	{
		double dResult = Export_LogNormal_TriSpreadOption(S1,S2,S3,sig1,sig2,sig3,
										 rho12,rho13,rho23,mu1,mu2,mu3,a0,a1,a2,a3,t,n);
		result.setDouble(dResult);
		return ARM_OK;
	}
	
	catch(Exception& x)
	{
		x.DebugPrint();
		
		ARM_RESULT();
	}
}

extern long ARMLOCAL_LogNormal_TriSpreadOption_Der(int i,
										   double S1,double S2,double S3,double sig1,double sig2,double sig3,
										double rho12,double rho13,double rho23,double mu1,double mu2,double mu3,
										double a0,double a1,double a2,double a3,double t,int n,
										   ARM_result& result
										   )
{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	
	try
	{
		double dResult = Export_LogNormal_TriSpreadOption(i,S1,S2,S3,sig1,sig2,sig3,
										 rho12,rho13,rho23,mu1,mu2,mu3,a0,a1,a2,a3,t,n);
		result.setDouble(dResult);
		return ARM_OK;
	}
	
	catch(Exception& x)
	{
		x.DebugPrint();
		
		ARM_RESULT();
	}
}


extern long ARMLOCAL_LogNormal_TriSpreadOption_Der2(int i,int j,
										   double S1,double S2,double S3,double sig1,double sig2,double sig3,
										double rho12,double rho13,double rho23,double mu1,double mu2,double mu3,
										double a0,double a1,double a2,double a3,double t,int n,
										   ARM_result& result
										   )
{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	
	try
	{
		double dResult = Export_LogNormal_TriSpreadOption(i,j,S1,S2,S3,sig1,sig2,sig3,
										 rho12,rho13,rho23,mu1,mu2,mu3,a0,a1,a2,a3,t,n);
		result.setDouble(dResult);
		return ARM_OK;
	}
	
	catch(Exception& x)
	{
		x.DebugPrint();
		
		ARM_RESULT();
	}
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////
///
///       Tri Spread Digital Option Formula
///
/////////////////////////////////////////////////////////////////////////////////////////////////////////



extern long ARMLOCAL_LogNormal_TriSpreadDigitalOption(
										   double S1,double S2,double S3,double sig1,double sig2,double sig3,
										double rho12,double rho13,double rho23,double mu1,double mu2,double mu3,
										double a0,double a1,double a2,double a3,double t,int callput,int n,
										   ARM_result& result
										   )
{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	
	try
	{
		double dResult = Export_LogNormal_TriSpreadDigitalOption(S1,S2,S3,sig1,sig2,sig3,
										 rho12,rho13,rho23,mu1,mu2,mu3,a0,a1,a2,a3,t,callput,n);
		result.setDouble(dResult);
		return ARM_OK;
	}
	
	catch(Exception& x)
	{
		x.DebugPrint();
		
		ARM_RESULT();
	}
}

extern long ARMLOCAL_LogNormal_TriSpreadDigitalOption_Der(int i,
										   double S1,double S2,double S3,double sig1,double sig2,double sig3,
										double rho12,double rho13,double rho23,double mu1,double mu2,double mu3,
										double a0,double a1,double a2,double a3,double t,int callput,int n,
										   ARM_result& result
										   )
{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	
	try
	{
		double dResult = Export_LogNormal_TriSpreadDigitalOption(i,S1,S2,S3,sig1,sig2,sig3,
										 rho12,rho13,rho23,mu1,mu2,mu3,a0,a1,a2,a3,t,callput,n);
		result.setDouble(dResult);
		return ARM_OK;
	}
	
	catch(Exception& x)
	{
		x.DebugPrint();
		
		ARM_RESULT();
	}
}


extern long ARMLOCAL_LogNormal_TriSpreadDigitalOption_Der2(int i,int j,
										   double S1,double S2,double S3,double sig1,double sig2,double sig3,
										double rho12,double rho13,double rho23,double mu1,double mu2,double mu3,
										double a0,double a1,double a2,double a3,double t,int callput,int n,
										   ARM_result& result
										   )
{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	
	try
	{
		double dResult = Export_LogNormal_TriSpreadDigitalOption(i,j,S1,S2,S3,sig1,sig2,sig3,
										 rho12,rho13,rho23,mu1,mu2,mu3,a0,a1,a2,a3,t,callput,n);
		result.setDouble(dResult);
		return ARM_OK;
	}
	
	catch(Exception& x)
	{
		x.DebugPrint();
		
		ARM_RESULT();
	}
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
///
///       Tri Spread Digital Option2 Formula
///
/////////////////////////////////////////////////////////////////////////////////////////////////////////



extern long ARMLOCAL_LogNormal_TriSpreadDigitalOption2(
										   double S1,double S2,double S3,double sig1,double sig2,double sig3,
										double rho12,double rho13,double rho23,double mu1,double mu2,double mu3,
										double a0,double a2,double a3,double b0,double b1,double t,int n,
										   ARM_result& result
										   )
{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	
	try
	{
		double dResult = Export_LogNormal_TriSpreadDigitalOption2(S1,S2,S3,sig1,sig2,sig3,
										 rho12,rho13,rho23,mu1,mu2,mu3,a0,a2,a3,b0,b1,t,n);
		result.setDouble(dResult);
		return ARM_OK;
	}
	
	catch(Exception& x)
	{
		x.DebugPrint();
		
		ARM_RESULT();
	}
}

extern long ARMLOCAL_LogNormal_TriSpreadDigitalOption2_Der(int i,
										   double S1,double S2,double S3,double sig1,double sig2,double sig3,
										double rho12,double rho13,double rho23,double mu1,double mu2,double mu3,
										double a0,double a2,double a3,double b0,double b1,double t,int n,
										   ARM_result& result
										   )
{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	
	try
	{
		double dResult = Export_LogNormal_TriSpreadDigitalOption2(i,S1,S2,S3,sig1,sig2,sig3,
										 rho12,rho13,rho23,mu1,mu2,mu3,a0,a2,a3,b0,b1,t,n);
		result.setDouble(dResult);
		return ARM_OK;
	}
	
	catch(Exception& x)
	{
		x.DebugPrint();
		
		ARM_RESULT();
	}
}


extern long ARMLOCAL_LogNormal_TriSpreadDigitalOption2_Der2(int i,int j,
										   double S1,double S2,double S3,double sig1,double sig2,double sig3,
										double rho12,double rho13,double rho23,double mu1,double mu2,double mu3,
										double a0,double a2,double a3,double b0,double b1,double t,int n,
										   ARM_result& result
										   )
{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	
	try
	{
		double dResult = Export_LogNormal_TriSpreadDigitalOption2(i,j,S1,S2,S3,sig1,sig2,sig3,
										 rho12,rho13,rho23,mu1,mu2,mu3,a0,a2,a3,b0,b1,t,n);
		result.setDouble(dResult);
		return ARM_OK;
	}
	
	catch(Exception& x)
	{
		x.DebugPrint();
		
		ARM_RESULT();
	}
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
///
///       Forward value  CEV VanillaOption
///
/////////////////////////////////////////////////////////////////////////////////////////////////////////



extern long ARMLOCAL_CEV_VanillaOption(double f, double K, double T,double drift, double sig, double beta,int callput, int nbsteps,
										   ARM_result& result
										   )
{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	
	try
	{
		double dResult = Export_CEV_VanillaOption(f,K,T,drift,sig,beta,callput,nbsteps);
		result.setDouble(dResult);
		return ARM_OK;
	}
	
	catch(Exception& x)
	{
		x.DebugPrint();
		
		ARM_RESULT();
	}
}

extern long ARMLOCAL_CEV_VanillaOption_Der(int i,
										double f, double K, double T,double drift, double sig, double beta,int callput, int nbsteps,
										   ARM_result& result
										   )
{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	
	try
	{
		double dResult = Export_CEV_VanillaOption(i,f,K,T,drift,sig,beta,callput,nbsteps);
		result.setDouble(dResult);
		return ARM_OK;
	}
	
	catch(Exception& x)
	{
		x.DebugPrint();
		
		ARM_RESULT();
	}
}


extern long ARMLOCAL_CEV_VanillaOption_Der2(int i,int j,
										double f, double K, double T,double drift, double sig, double beta,int callput, int nbsteps,
										   ARM_result& result
										   )
{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	
	try
	{
		double dResult = Export_CEV_VanillaOption(i,j,f,K,T,drift,sig,beta,callput,nbsteps);
		result.setDouble(dResult);
		return ARM_OK;
	}
	
	catch(Exception& x)
	{
		x.DebugPrint();
		
		ARM_RESULT();
	}
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
///
///       Forward value  CEV DoubleBarrierOption
///
/////////////////////////////////////////////////////////////////////////////////////////////////////////



extern long ARMLOCAL_CEV_DoubleBarrierOption(double f, double K, double T,double barrdown,double barrup,double drift, double sig, double beta,int callput, int nbsteps,
										   ARM_result& result
										   )
{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	
	try
	{
		double dResult = Export_CEV_DoubleBarrierOption(f,K,T,barrdown,barrup,drift,sig,beta,callput,nbsteps);
		result.setDouble(dResult);
		return ARM_OK;
	}
	
	catch(Exception& x)
	{
		x.DebugPrint();
		
		ARM_RESULT();
	}
}

extern long ARMLOCAL_CEV_DoubleBarrierOption_Der(int i,
										double f, double K, double T,double barrdown,double barrup,double drift, double sig, double beta,int callput, int nbsteps,
										   ARM_result& result
										   )
{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	
	try
	{
		double dResult = Export_CEV_DoubleBarrierOption(i,f,K,T,barrdown,barrup,drift,sig,beta,callput,nbsteps);
		result.setDouble(dResult);
		return ARM_OK;
	}
	
	catch(Exception& x)
	{
		x.DebugPrint();
		
		ARM_RESULT();
	}
}


extern long ARMLOCAL_CEV_DoubleBarrierOption_Der2(int i,int j,
										double f, double K, double T,double barrdown,double barrup,double drift, double sig, double beta,int callput, int nbsteps,
										   ARM_result& result
										   )
{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	
	try
	{
		double dResult = Export_CEV_DoubleBarrierOption(i,j,f,K,T,barrdown,barrup,drift,sig,beta,callput,nbsteps);
		result.setDouble(dResult);
		return ARM_OK;
	}
	
	catch(Exception& x)
	{
		x.DebugPrint();
		
		ARM_RESULT();
	}
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
///
///       Forward value  CEV BarrierOption
///
/////////////////////////////////////////////////////////////////////////////////////////////////////////



extern long ARMLOCAL_CEV_BarrierOption(double f, double K, double T,double barrier,double drift, double sig, double beta,int optype,int callput, int nbsteps,
										   ARM_result& result
										   )
{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	
	try
	{
		double dResult = Export_CEV_BarrierOption(f,K,T,barrier,drift,sig,beta,optype,callput,nbsteps);
		result.setDouble(dResult);
		return ARM_OK;
	}
	
	catch(Exception& x)
	{
		x.DebugPrint();
		
		ARM_RESULT();
	}
}

extern long ARMLOCAL_CEV_BarrierOption_Der(int i,
										double f, double K, double T,double barrier,double drift, double sig, double beta,int optype,int callput, int nbsteps,
										   ARM_result& result
										   )
{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	
	try
	{
		double dResult = Export_CEV_BarrierOption(i,f,K,T,barrier,drift,sig,beta,optype,callput,nbsteps);
		result.setDouble(dResult);
		return ARM_OK;
	}
	
	catch(Exception& x)
	{
		x.DebugPrint();
		
		ARM_RESULT();
	}
}


extern long ARMLOCAL_CEV_BarrierOption_Der2(int i,int j,
										double f, double K, double T,double barrier,double drift, double sig, double beta,int optype,int callput, int nbsteps,
										   ARM_result& result
										   )
{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	
	try
	{
		double dResult = Export_CEV_BarrierOption(i,j,f,K,T,barrier,drift,sig,beta,optype,callput,nbsteps);
		result.setDouble(dResult);
		return ARM_OK;
	}
	
	catch(Exception& x)
	{
		x.DebugPrint();
		
		ARM_RESULT();
	}
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////
///
///       Asian Lognormal vanilla option (geman yor formula)
///
/////////////////////////////////////////////////////////////////////////////////////////////////////////



extern long ARMLOCAL_Lognormal_Asian_VanillaOption(double f, double K, double T, double r, double sig,int alpha,int callput, int nbsteps,
										   ARM_result& result
										   )
{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	
	try
	{
		double dResult = Export_Lognormal_Asian_VanillaOption(f,K,T,r,sig,alpha,callput,nbsteps);
		result.setDouble(dResult);
		return ARM_OK;
	}
	
	catch(Exception& x)
	{
		x.DebugPrint();
		
		ARM_RESULT();
	}
}

extern long ARMLOCAL_Lognormal_Asian_VanillaOption_Der(int i,
										double f, double K, double T, double r, double sig,int alpha,int callput, int nbsteps,
										   ARM_result& result
										   )
{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	
	try
	{
		double dResult = Export_Lognormal_Asian_VanillaOption(i,f,K,T,r,sig,alpha,callput,nbsteps);
		result.setDouble(dResult);
		return ARM_OK;
	}
	
	catch(Exception& x)
	{
		x.DebugPrint();
		
		ARM_RESULT();
	}
}


extern long ARMLOCAL_Lognormal_Asian_VanillaOption_Der2(int i,int j,
										double f, double K, double T, double r, double sig,int alpha,int callput, int nbsteps,
										   ARM_result& result
										   )
{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	
	try
	{
		double dResult = Export_Lognormal_Asian_VanillaOption(i,j,f,K,T,r,sig,alpha,callput,nbsteps);
		result.setDouble(dResult);
		return ARM_OK;
	}
	
	catch(Exception& x)
	{
		x.DebugPrint();
		
		ARM_RESULT();
	}
}

extern long ARMLOCAL_GaussianIntegrals_Legendre_Coeffs(
													   double a,
													   double b,
													   int n,ARM_result& result
										   )
{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	int i;
	
	try
	{
		GaussLegendre_Coefficients* gauss_ptr = Export_GaussLegendre_Coefficients(a,b,n);
		int nbpts=gauss_ptr->get_order();
		result.setLong(2*nbpts);

		for ( i=0;i<nbpts;i++)
		{
			result.setArray(gauss_ptr->get_point(i),2*i);
			result.setArray(gauss_ptr->get_weight(i),2*i+1);
		}
		delete gauss_ptr;
		return ARM_OK;
	}
	
	catch(Exception& x)
	{
		x.DebugPrint();
		
		ARM_RESULT();
	}
}

extern long ARMLOCAL_GaussianIntegrals_Hermite_Coeffs(int n,ARM_result& result
										   )
{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	int i;
	
	try
	{
		GaussHermite_Coefficients* gauss_ptr = Export_GaussHermite_Coefficients(n);
		int nbpts=gauss_ptr->get_order();
		result.setLong(2*nbpts);

		for ( i=0;i<nbpts;i++)
		{
			result.setArray(gauss_ptr->get_point(i),2*i);
			result.setArray(gauss_ptr->get_weight(i),2*i+1);
		}
		delete gauss_ptr;
		return ARM_OK;
	}
	
	catch(Exception& x)
	{
		x.DebugPrint();
		
		ARM_RESULT();
	}
}


extern long ARMLOCAL_GaussianIntegrals_Laguerre_Coeffs(double a,int n,ARM_result& result
										   )
{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	int i;
	
	try
	{
		GaussLaguerre_Coefficients* gauss_ptr = Export_GaussLaguerre_Coefficients(a,n);
		int nbpts=gauss_ptr->get_order();
		result.setLong(2*nbpts);

		for ( i=0;i<nbpts;i++)
		{
			result.setArray(gauss_ptr->get_point(i),2*i);
			result.setArray(gauss_ptr->get_weight(i),2*i+1);
		}
		delete gauss_ptr;
		return ARM_OK;
	}
	
	catch(Exception& x)
	{
		x.DebugPrint();
		
		ARM_RESULT();
	}
}




long ARMLOCAL_ConvertFXOptionToStrike( double C_TargetDelta, double C_Fwd, double C_TotalVol, int C_callPut, double C_initValue, ARM_result& result )
{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	try
	{
		double dResult = ConvertFXOptionPriceToStrike(C_TargetDelta, C_Fwd, C_TotalVol, C_callPut, C_initValue );
		result.setDouble(dResult);
		return ARM_OK;
	}
	
	catch(Exception& x)
	{
		x.DebugPrint();
		
		ARM_RESULT();
	}
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
///
///       Generalized Heston Pricer with a vector of models
///
/////////////////////////////////////////////////////////////////////////////////////////////////////////
long ARMLOCAL_GHeston_VanillaOption_ModelVector(
		double K,
		double T,
		int callput,
		int Interpolation_Method,
		const vector<double>& C_Maturity_Vec,
		const vector<double>& C_F_Vec,
		const vector<double>& C_InitialVol_Vec,
		const vector<double>& C_longtermV_Vec,
		const vector<double>& C_theta_Vec,
		const vector<double>& C_ksi_Vec,
		const vector<double>& C_rho_Vec,
		const vector<double>& C_lambda_Vec,
		const vector<double>& C_muJ_Vec,
		const vector<double>& C_sigmaJ_Vec,
		int nbsteps,
		ARM_result&  result)
{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	std::vector<double>& Maturity_Vec	= NULL;	
	std::vector<double>& F_Vec	= NULL;	
	std::vector<double>& InitialVol_Vec = NULL;
	std::vector<double>& longtermV_Vec = NULL;
	std::vector<double>& theta_Vec = NULL;
	std::vector<double>& ksi_Vec = NULL;
	std::vector<double>& rho_Vec = NULL;
	std::vector<double>& lambda_Vec = NULL;
	std::vector<double>& muJ_Vec = NULL;
	std::vector<double>& sigmaJ_Vec = NULL;

	try
	{

		Maturity_Vec	= CreateARMGPVectorFromVECTOR( C_Maturity_Vec );
		F_Vec			= CreateARMGPVectorFromVECTOR( C_F_Vec );
		InitialVol_Vec	= CreateARMGPVectorFromVECTOR( C_InitialVol_Vec);
		longtermV_Vec	= CreateARMGPVectorFromVECTOR( C_longtermV_Vec);
		theta_Vec		= CreateARMGPVectorFromVECTOR( C_theta_Vec);
		ksi_Vec			= CreateARMGPVectorFromVECTOR( C_ksi_Vec);
		rho_Vec			= CreateARMGPVectorFromVECTOR( C_rho_Vec);
		lambda_Vec		= CreateARMGPVectorFromVECTOR( C_lambda_Vec);
		muJ_Vec			= CreateARMGPVectorFromVECTOR( C_muJ_Vec);
		sigmaJ_Vec		= CreateARMGPVectorFromVECTOR( C_sigmaJ_Vec);
		
		double dResult = Export_GHeston_VanillaOption_ModelVector(K,T,callput,Interpolation_Method,
			Maturity_Vec,
			F_Vec,
			InitialVol_Vec,
			longtermV_Vec,
			theta_Vec,
			ksi_Vec,
			rho_Vec,
			lambda_Vec,
			muJ_Vec,
			sigmaJ_Vec,

			nbsteps);

		result.setDouble(dResult);

		delete Maturity_Vec;
		delete F_Vec;
		delete InitialVol_Vec;
		delete longtermV_Vec;
		delete theta_Vec;
		delete ksi_Vec;
		delete rho_Vec;
		delete lambda_Vec;
		delete muJ_Vec;
		delete sigmaJ_Vec;
		
		return ARM_OK;
	}
	
	catch(Exception& x)
	{
		x.DebugPrint();
		
		ARM_RESULT();
	}
}

long ARMLOCAL_GHeston_Implicit_Volatility_ModelVector(
		double K,
		double T,
		int Interpolation_Method,
		const vector<double>& C_Maturity_Vec,
		const vector<double>& C_F_Vec,
		const vector<double>& C_InitialVol_Vec,
		const vector<double>& C_longtermV_Vec,
		const vector<double>& C_theta_Vec,
		const vector<double>& C_ksi_Vec,
		const vector<double>& C_rho_Vec,
		const vector<double>& C_lambda_Vec,
		const vector<double>& C_muJ_Vec,
		const vector<double>& C_sigmaJ_Vec,
		int nbsteps,
		ARM_result&  result)
{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	std::vector<double>& Maturity_Vec	= NULL;	
	std::vector<double>& F_Vec	= NULL;	
	std::vector<double>& InitialVol_Vec = NULL;
	std::vector<double>& longtermV_Vec = NULL;
	std::vector<double>& theta_Vec = NULL;
	std::vector<double>& ksi_Vec = NULL;
	std::vector<double>& rho_Vec = NULL;
	std::vector<double>& lambda_Vec = NULL;
	std::vector<double>& muJ_Vec = NULL;
	std::vector<double>& sigmaJ_Vec = NULL;

	try
	{

		Maturity_Vec	= CreateARMGPVectorFromVECTOR( C_Maturity_Vec );
		F_Vec			= CreateARMGPVectorFromVECTOR( C_F_Vec );
		InitialVol_Vec	= CreateARMGPVectorFromVECTOR( C_InitialVol_Vec);
		longtermV_Vec	= CreateARMGPVectorFromVECTOR( C_longtermV_Vec);
		theta_Vec		= CreateARMGPVectorFromVECTOR( C_theta_Vec);
		ksi_Vec			= CreateARMGPVectorFromVECTOR( C_ksi_Vec);
		rho_Vec			= CreateARMGPVectorFromVECTOR( C_rho_Vec);
		lambda_Vec		= CreateARMGPVectorFromVECTOR( C_lambda_Vec);
		muJ_Vec			= CreateARMGPVectorFromVECTOR( C_muJ_Vec);
		sigmaJ_Vec		= CreateARMGPVectorFromVECTOR( C_sigmaJ_Vec);
		
		double dResult = Export_GHeston_Implicit_Volatility_ModelVector(K,T,Interpolation_Method,
			Maturity_Vec,
			F_Vec,
			InitialVol_Vec,
			longtermV_Vec,
			theta_Vec,
			ksi_Vec,
			rho_Vec,
			lambda_Vec,
			muJ_Vec,
			sigmaJ_Vec,

			nbsteps);

		result.setDouble(dResult);

		delete Maturity_Vec;
		delete F_Vec;
		delete InitialVol_Vec;
		delete longtermV_Vec;
		delete theta_Vec;
		delete ksi_Vec;
		delete rho_Vec;
		delete lambda_Vec;
		delete muJ_Vec;
		delete sigmaJ_Vec;
		
		return ARM_OK;
	}
	
	catch(Exception& x)
	{
		x.DebugPrint();
		
		ARM_RESULT();
	}
}

/////////////////////////////////////////////////////////////////////
///
///     Calibration of SABR Models
///
/////////////////////////////////////////////////////////////////////

long ARMLOCAL_SABR_Model_Calibrate(
		const vector<double>& C_K_Vec,
		const vector<double>& C_ImpVol_Vec,
		const vector<double>& C_Weigth_Vec,
		double f,double t,
		int flag ,int nbsteps,int algorithm,
		double alpha0,double beta0,double rho0,double nu0,
		ARM_result&  result)
{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	std::vector<double>& K_Vec	= NULL;	
	std::vector<double>& ImpVol_Vec	= NULL;	
	std::vector<double>& Weigth_Vec	= NULL;	
	
	try
	{
		K_Vec				= CreateARMGPVectorFromVECTOR( C_K_Vec );
		ImpVol_Vec			= CreateARMGPVectorFromVECTOR( C_ImpVol_Vec );
		Weigth_Vec			= CreateARMGPVectorFromVECTOR( C_Weigth_Vec );
		
		SABR_ParameterSet* setptr = Export_SABR_Model_Calibrate(K_Vec,ImpVol_Vec,Weigth_Vec,f,t,flag,
			nbsteps,algorithm,alpha0,beta0,rho0,nu0);

		result.setLong(6);
		result.setArray(setptr->get_alpha(),0);
		result.setArray(setptr->get_beta(),1);
		result.setArray(setptr->get_rho(),2);
		result.setArray(setptr->get_nu(),3);
		result.setArray(f,4);
		result.setArray(setptr->get_objective(),5);
		delete K_Vec;
		delete ImpVol_Vec;
		delete Weigth_Vec;
		delete setptr;
		
		return ARM_OK;
	}
	catch(Exception& x)
	{
		x.DebugPrint();
		
		ARM_RESULT();
	}
}

long ARMLOCAL_SABR_Model_Calibrate_WForward(
		const vector<double>& C_K_Vec,
		const vector<double>& C_ImpVol_Vec,
		const vector<double>& C_Weigth_Vec,
		double t,
		int flag ,int nbsteps,int algorithm,
		double alpha0,double beta0,double rho0,double nu0,double f0,
		ARM_result&  result)
{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	std::vector<double>& K_Vec	= NULL;	
	std::vector<double>& ImpVol_Vec	= NULL;	
	std::vector<double>& Weigth_Vec	= NULL;	
	
	try
	{
		K_Vec				= CreateARMGPVectorFromVECTOR( C_K_Vec );
		ImpVol_Vec			= CreateARMGPVectorFromVECTOR( C_ImpVol_Vec );
		Weigth_Vec			= CreateARMGPVectorFromVECTOR( C_Weigth_Vec );
		
		SABR_ParameterSet* setptr = Export_SABR_Model_Calibrate(K_Vec,ImpVol_Vec,Weigth_Vec,t,flag,
			nbsteps,algorithm,alpha0,beta0,rho0,nu0,f0);

		result.setLong(6);
		result.setArray(setptr->get_alpha(),0);
		result.setArray(setptr->get_beta(),1);
		result.setArray(setptr->get_rho(),2);
		result.setArray(setptr->get_nu(),3);
		result.setArray(setptr->get_forward(),4);
		result.setArray(setptr->get_objective(),5);
		delete K_Vec;
		delete ImpVol_Vec;
		delete Weigth_Vec;
		delete setptr;
		
		return ARM_OK;
	}
	catch(Exception& x)
	{
		x.DebugPrint();
		
		ARM_RESULT();
	}
}

long ARMLOCAL_SABR_Model_Calibrate_FixedBeta(
		const vector<double>& C_K_Vec,
		const vector<double>& C_ImpVol_Vec,
		const vector<double>& C_Weigth_Vec,
		double f,
		double beta,
		double t,
		int flag ,double Alpha_Exp,double Alpha_Tanh,double Kb_Tanh,
		int nbsteps,int algorithm,
		double alpha0,double rho0,double nu0,
		ARM_result&  result)
{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	std::vector<double>& K_Vec	= NULL;	
	std::vector<double>& ImpVol_Vec	= NULL;	
	std::vector<double>& Weigth_Vec	= NULL;	
	
	try
	{
		K_Vec				= CreateARMGPVectorFromVECTOR( C_K_Vec );
		ImpVol_Vec			= CreateARMGPVectorFromVECTOR( C_ImpVol_Vec );
		Weigth_Vec			= CreateARMGPVectorFromVECTOR( C_Weigth_Vec );
		
		SABR_ParameterSet* setptr = Export_SABR_Model_Calibrate(K_Vec,ImpVol_Vec,Weigth_Vec,f,beta,t,
			flag,Alpha_Exp,Alpha_Tanh,Kb_Tanh,nbsteps,algorithm,alpha0,rho0,nu0);

		result.setLong(6);
		result.setArray(setptr->get_alpha(),0);
		result.setArray(beta,1);
		result.setArray(setptr->get_rho(),2);
		result.setArray(setptr->get_nu(),3);
		result.setArray(f,4);
		result.setArray(setptr->get_objective(),5);
		delete K_Vec;
		delete ImpVol_Vec;
		delete Weigth_Vec;
		delete setptr;
		
		return ARM_OK;
	}
	catch(Exception& x)
	{
		x.DebugPrint();
		
		ARM_RESULT();
	}
}


long ARMLOCAL_SABR_Model_Calibrate_FixedBeta_Linked(
		const vector<double>& C_K_Vec,
		const vector<double>& C_ImpVol_Vec,
		const vector<double>& C_Weigth_Vec,
		double f,
		double beta,
		double t,
		int flag ,int nbsteps,int algorithm,
		double alpha0,double rho0,double nu0,
		double alphap,double rhop, double nup,double rweight_alpha,double rweight_rho,double rweight_nu,
		ARM_result&  result)
{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	std::vector<double>& K_Vec	= NULL;	
	std::vector<double>& ImpVol_Vec	= NULL;	
	std::vector<double>& Weigth_Vec	= NULL;	
	
	try
	{
		K_Vec				= CreateARMGPVectorFromVECTOR( C_K_Vec );
		ImpVol_Vec			= CreateARMGPVectorFromVECTOR( C_ImpVol_Vec );
		Weigth_Vec			= CreateARMGPVectorFromVECTOR( C_Weigth_Vec );
		
		SABR_ParameterSet* setptr = Export_SABR_Model_Calibrate(K_Vec,ImpVol_Vec,Weigth_Vec,f,beta,t,flag,nbsteps,
											algorithm,alpha0,rho0,nu0,alphap,rhop,nup,rweight_alpha,rweight_rho,rweight_nu);

		result.setLong(6);
		result.setArray(setptr->get_alpha(),0);
		result.setArray(beta,1);
		result.setArray(setptr->get_rho(),2);
		result.setArray(setptr->get_nu(),3);
		result.setArray(f,4);
		result.setArray(setptr->get_objective(),5);
		delete K_Vec;
		delete ImpVol_Vec;
		delete Weigth_Vec;
		delete setptr;
		
		return ARM_OK;
	}
	catch(Exception& x)
	{
		x.DebugPrint();
		
		ARM_RESULT();
	}
}

long ARMLOCAL_SABR_Calibrate_BetaFixedToOne(
											 const vector<double>& C_K_Vec,
											 const vector<double>& C_ImpVol_Vec,
											 const vector<double>& C_Weigth_Vec,
											 double f,
											 double t,
											 double atmvol,
											 ARM_result&  result)
{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	std::vector<double>& K_Vec	= NULL;	
	std::vector<double>& ImpVol_Vec	= NULL;	
	std::vector<double>& Weigth_Vec	= NULL;	
	
	try
	{
		K_Vec				= CreateARMGPVectorFromVECTOR( C_K_Vec );
		ImpVol_Vec			= CreateARMGPVectorFromVECTOR( C_ImpVol_Vec );
		Weigth_Vec			= CreateARMGPVectorFromVECTOR( C_Weigth_Vec );
		
		SABR_ParameterSet* setptr = SABR_CalibrateToSmileBetaFixedToOne(K_Vec,ImpVol_Vec,Weigth_Vec,f,t,atmvol);

		result.setLong(6);
		result.setArray(setptr->get_alpha(),0);
		result.setArray(1.,1);
		result.setArray(setptr->get_rho(),2);
		result.setArray(setptr->get_nu(),3);
		result.setArray(f,4);
		result.setArray(setptr->get_objective(),5);
		delete K_Vec;
		delete ImpVol_Vec;
		delete Weigth_Vec;
		delete setptr;
		
		return ARM_OK;
	}
	catch(Exception& x)
	{
		x.DebugPrint();
		
		ARM_RESULT();
	}
}

/////////////////////////////////////////////////////////////////////////////////////
///
///		Calibration of GHeston
///
/////////////////////////////////////////////////////////////////////////////////////

long ARMLOCAL_GHeston_Model_Calibrate_Total(
			const vector<double>& C_K_Vec,
			const vector<double>& C_ImpVol_Vec,
			double C_f,
			double C_t,
			int C_nb,
			int C_algorithm,
			double C_V0,
			double C_omega,
			double C_theta,
			double C_ksi,
			double C_rho,
			double C_muJ,
			double C_sigmaJ,
			double C_lambda,
			ARM_result&  result)
{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	std::vector<double>& K_Vec	= NULL;	
	std::vector<double>& ImpVol_Vec	= NULL;	
	
	try
	{
		K_Vec				= CreateARMGPVectorFromVECTOR( C_K_Vec );
		ImpVol_Vec			= CreateARMGPVectorFromVECTOR( C_ImpVol_Vec );
		
		GeneralizedHeston_ParameterSet* setptr = Export_GeneralizedHeston_Model_Calibrate(K_Vec,ImpVol_Vec,C_f,C_t,C_nb,C_algorithm,C_V0,C_omega,C_theta,C_ksi,C_rho,C_muJ,C_sigmaJ,C_lambda);

		result.setLong(9);
		result.setArray(setptr->get_V0(),0);
		result.setArray(setptr->get_omega(),1);
		result.setArray(setptr->get_theta(),2);
		result.setArray(setptr->get_ksi(),3);
		result.setArray(setptr->get_rho(),4);
		result.setArray(setptr->get_muJ(),5);
		result.setArray(setptr->get_sigmaJ(),6);
		result.setArray(setptr->get_lambda(),7);
		result.setArray(setptr->get_objective(),8);
		delete K_Vec;
		delete ImpVol_Vec;
		delete setptr;
		
		return ARM_OK;
	}
	catch(Exception& x)
	{
		x.DebugPrint();
		
		ARM_RESULT();
	}

}

long ARMLOCAL_GHeston_Model_Calibrate_NoJump(
			const vector<double>& C_K_Vec,
			const vector<double>& C_ImpVol_Vec,
			double C_f,
			double C_t,
			int C_nb,
			int C_algorithm,
			double C_V0,
			double C_omega,
			double C_theta,
			double C_ksi,
			double C_rho,
			double C_muJ0,
			double C_sigmaJ0,
			double C_lambda0,
			ARM_result&  result)
{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	std::vector<double>& K_Vec	= NULL;	
	std::vector<double>& ImpVol_Vec	= NULL;	
	
	try
	{
		K_Vec				= CreateARMGPVectorFromVECTOR( C_K_Vec );
		ImpVol_Vec			= CreateARMGPVectorFromVECTOR( C_ImpVol_Vec );
		
		GeneralizedHeston_ParameterSet* setptr = Export_GeneralizedHeston_Model_Calibrate_NoJump(K_Vec,ImpVol_Vec,C_f,C_t,C_muJ0,C_sigmaJ0,C_lambda0,C_nb,C_algorithm,C_V0,C_omega,C_theta,C_ksi,C_rho);

		result.setLong(9);
		result.setArray(setptr->get_V0(),0);
		result.setArray(setptr->get_omega(),1);
		result.setArray(setptr->get_theta(),2);
		result.setArray(setptr->get_ksi(),3);
		result.setArray(setptr->get_rho(),4);
		result.setArray(setptr->get_muJ(),5);
		result.setArray(setptr->get_sigmaJ(),6);
		result.setArray(setptr->get_lambda(),7);
		result.setArray(setptr->get_objective(),8);
		delete K_Vec;
		delete ImpVol_Vec;
		delete setptr;
		
		return ARM_OK;
	}
	catch(Exception& x)
	{
		x.DebugPrint();
		
		ARM_RESULT();
	}

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
///
///        SABR Student Power Spreadoption
///
/////////////////////////////////////////////////////////////////////////////////////////////////////////


extern long ARMLOCAL_Student_SABR_Power_SpreadOption(
	double S1,
	double S2,
	double alpha1,
	double beta1,
	double rho1,
	double nu1,
	double alpha2,
	double beta2, 
	double rho2, 
	double nu2,
	const vector<double>&  C_Copula_vec,
	double t,
	int flag,
	double a10,
	double b10,
	double k10,
	double a20,
	double b20,
	double k20,
	const vector<double>&  C_Parameters_vec,
	ARM_result& result
	)
	{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
		std::vector<double>&		Parameters_vec		= NULL;	
		std::vector<double>&		Copula_vec		= NULL;	

	try
	{	
		Copula_vec		= CreateARMGPVectorFromVECTOR( C_Copula_vec	);
		const vector<double> copula_Vec=Copula_vec->GetValues();
		double copula_corr=copula_Vec[0];
		double copula_degre=copula_Vec[1];
		
		Parameters_vec		= CreateARMGPVectorFromVECTOR( C_Parameters_vec	);
		const vector<double> params_Vec=Parameters_vec->GetValues();
		int n=params_Vec[0];
		double alpha_exp=params_Vec[1];
		double alpha_tanh=params_Vec[2];
		double kb_tanh=params_Vec[3];
		double dResult = Export_Student_SABR_Power_SpreadOption(S1,S2,
													alpha1,beta1,rho1,nu1,
													alpha2,beta2,rho2,nu2,
													copula_corr,copula_degre,t,flag,
												a10,b10,k10,a20,b20,k20,n,alpha_exp,alpha_tanh,kb_tanh);
		result.setDouble(dResult);
		return ARM_OK;
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}
}

extern long ARMLOCAL_Student_SABR_Power_Digital_SpreadOption(
	double S1,
	double S2,
	double alpha1,
	double beta1,
	double rho1,
	double nu1,
	double alpha2,
	double beta2, 
	double rho2, 
	double nu2,
	double copula_corr,
	double copula_degre,
	double t,
	int flag,
	double k10,
	double a20,
	double b20,
	double k20,
	int n,double alpha_exp,double alpha_tanh,double kb_tanh,
	ARM_result& result
	)
	{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	try
	{
		double dResult = Export_Student_SABR_Power_Digital_SpreadOption(S1,S2,
													alpha1,beta1,rho1,nu1,
													alpha2,beta2,rho2,nu2,
													copula_corr,copula_degre,t,flag,
												k10,a20,b20,k20,n,alpha_exp,alpha_tanh,kb_tanh);
		result.setDouble(dResult);
		return ARM_OK;
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}
}



/////////////////////////////////////////////////////////////////////////////////////////////////////////
///
///        Mepi Vanilla Option 
///
/////////////////////////////////////////////////////////////////////////////////////////////////////////



extern long ARMLOCAL_Mepi_VanillaOption_STOBS(
	double date,
	double P0,
	double K,
	double T,
	int zero_curve_id,
	double b,
	double YearlyFees,
	double cashspread,
	int sabr_model_id,
	double minExp,
	double maxExp,
	double riskFac,
	double g0,
	double g,
	int AveragingPeriodNb,
	double asianreset,
	int	   callOrPut,
	int	   n,
	int nz,
	int nh,
	ARM_result& result
	)
{
    	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	try
	{
        double dResult;
		double minport=0.0;
		double maxport=2.0*P0;

		
		ARM_ZeroCurve* zero_curve_ptr = NULL;	
		if( !GetObjectFromId( &zero_curve_ptr, zero_curve_id, ARM_ZERO_CURVE ) )	
		{
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"interest rate fwd curve is not of a good type:");	
		}
		ARM_SABR_Eq* sabr_model_ptr = (ARM_SABR_Eq*)(LOCAL_PERSISTENT_OBJECTS->GetObject(sabr_model_id));	

		dResult=	mepi_VanillaOption_STOBS(date,P0,K, T, zero_curve_ptr,"EUR",
						   b,YearlyFees, cashspread,sabr_model_ptr,"NTHG",
						  minExp,maxExp,riskFac,g0,g,
						  minport,maxport,AveragingPeriodNb,asianreset,
						  callOrPut,n,nz,nh);
            
		result.setDouble(dResult);
		return ARM_OK;
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}
}

extern long ARMLOCAL_Mepi_VanillaOption_STOBS_delta(
	double date,
	double P0,
	double K,
	double T,
	int zero_curve_id,
	double b,
	double YearlyFees,
	double cashspread,
	int sabr_model_id,
	double minExp,
	double maxExp,
	double riskFac,
	double g0,
	double g,
	int AveragingPeriodNb,
	double asianreset,
	int	   callOrPut,
	int	   n,
	int nz,
	int nh,
	ARM_result& result
	)
{
    	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	try
	{
        double dResult;
		double minport=0.0;
		double maxport=2.0*P0;

		
		ARM_ZeroCurve* zero_curve_ptr = NULL;	
		if( !GetObjectFromId( &zero_curve_ptr, zero_curve_id, ARM_ZERO_CURVE ) )	
		{
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"interest rate fwd curve is not of a good type:");	
		}
		ARM_SABR_Eq* sabr_model_ptr = (ARM_SABR_Eq*)(LOCAL_PERSISTENT_OBJECTS->GetObject(sabr_model_id));	

		dResult=	mepi_VanillaOption_STOBS_delta(date,P0,K, T, zero_curve_ptr,"",
						   b,YearlyFees,cashspread, sabr_model_ptr,"",
						  minExp,maxExp,riskFac,g0,g,
						  minport,maxport,AveragingPeriodNb,asianreset,
						  callOrPut,n,nz,nh);
            
		result.setDouble(dResult);
		return ARM_OK;
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}
}


extern long ARMLOCAL_Mepi_VanillaOption_SABR(
	double date,
	double P0,
	double K,
	double T,
	int zero_curve_id,
	double b,
	double YearlyFees,
	double cashspread,
	int sabr_model_id,
	CCString C_SpotName,
	double minExp,
	double maxExp,
	double riskFac,
	double g0,
	double g,
	int AveragingPeriodNb,
	double asianreset,
	int	   callOrPut,
	int	   n,
	int nz,
	int nh,
	ARM_result& result
	)
{
    	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	try
	{
        double dResult;
		double minport=0.0;
		double maxport=2.0*P0;

		
		ARM_ZeroCurve* zero_curve_ptr = NULL;	
		if( !GetObjectFromId( &zero_curve_ptr, zero_curve_id, ARM_ZERO_CURVE ) )	
		{
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"interest rate fwd curve is not of a good type:");	
		}

		ARM_ZeroCurve*  zerocurve=&(*zero_curve_ptr);
		string C_CurveName=zerocurve->GetStrCurrency();
		ARM_SABR_Eq* sabr_model_ptr = (ARM_SABR_Eq*)(LOCAL_PERSISTENT_OBJECTS->GetObject(sabr_model_id));	

		dResult=	mepi_VanillaOption_SABR(ConvertXLDateToJulian(date),P0,K, ConvertXLDateToJulian(T), zero_curve_ptr,C_CurveName.c_str(),
						   b,YearlyFees,cashspread, sabr_model_ptr,C_SpotName.c_str(),
						  minExp,maxExp,riskFac,g0,g,
						  minport,maxport,AveragingPeriodNb,asianreset,n);
            
		result.setDouble(dResult);
		return ARM_OK;
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}
}

extern long ARMLOCAL_Mepi_VanillaOption_SABR_delta(
	double date,
	double P0,
	double K,
	double T,
	int zero_curve_id,
	double b,
	double YearlyFees,
	double cashspread,
	int sabr_model_id,
	CCString C_SpotName,
	double minExp,
	double maxExp,
	double riskFac,
	double g0,
	double g,
	int AveragingPeriodNb,
	double asianreset,
	int	   callOrPut,
	int	   n,
	int nz,
	int nh,
	ARM_result& result
	)
{
    	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	try
	{
        double dResult;
		double minport=0.0;
		double maxport=2.0*P0;

		
		ARM_ZeroCurve* zero_curve_ptr = NULL;	
		if( !GetObjectFromId( &zero_curve_ptr, zero_curve_id, ARM_ZERO_CURVE ) )	
		{
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"interest rate fwd curve is not of a good type:");	
		}
		ARM_ZeroCurve*  zerocurve=&(*zero_curve_ptr);
		string C_CurveName=zerocurve->GetStrCurrency();
	

		ARM_SABR_Eq* sabr_model_ptr = (ARM_SABR_Eq*)(LOCAL_PERSISTENT_OBJECTS->GetObject(sabr_model_id));	

		dResult=	mepi_VanillaOption_SABR_delta(ConvertXLDateToJulian(date),P0,K, ConvertXLDateToJulian(T), zero_curve_ptr,C_CurveName.c_str(),
						   b,YearlyFees,cashspread, sabr_model_ptr,C_SpotName.c_str(),
						  minExp,maxExp,riskFac,g0,g,
						  minport,maxport,AveragingPeriodNb,asianreset,n,0.01);
            
		result.setDouble(dResult);
		return ARM_OK;
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////
///
///       Present value  of a VanillaOption on an underlying with stochastic volatility
///
/////////////////////////////////////////////////////////////////////////////////////////////////////////





extern long ARMLOCAL_StochasticVol_LN_VanillaOption(double f, double K, double T,double drift, double sig, double VolDrift,double VolVol,	double averaging,double reset,int callput, int nbsteps,
										   ARM_result& result
										   )
{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	
	try
	{
		double dResult = Export_StochasticVol_LN_Geometric_VanillaOption_with_Reset(f,K,T,drift,sig,VolDrift,VolVol,averaging,reset,callput,nbsteps);
		result.setDouble(dResult);
		return ARM_OK;
	}
	
	catch(Exception& x)
	{
		x.DebugPrint();
		
		ARM_RESULT();
	}
}

extern long ARMLOCAL_StochasticVol_LN_VanillaOption_Der(int i,double f, double K, double T,double drift, double sig, double VolDrift,double VolVol,	double averaging,double reset,int callput, int nbsteps,
										   ARM_result& result
										   )
{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	
	try
	{
		double dResult = Export_StochasticVol_LN_Geometric_VanillaOption_with_Reset(i,f,K,T,drift,sig,VolDrift,VolVol,averaging,reset,callput,nbsteps);
		result.setDouble(dResult);
		return ARM_OK;
	}
	
	catch(Exception& x)
	{
		x.DebugPrint();
		
		ARM_RESULT();
	}
}


extern long ARMLOCAL_StochasticVol_LN_VanillaOption_Der2(int i,int j,double f, double K, double T,double drift, double sig, double VolDrift,double VolVol,	double averaging,double reset,int callput, int nbsteps,
										   ARM_result& result
										   )
{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	
	try
	{
		double dResult = Export_StochasticVol_LN_Geometric_VanillaOption_with_Reset(i,j,f,K,T,drift,sig,VolDrift,VolVol,averaging,reset,callput,nbsteps);
		result.setDouble(dResult);
		return ARM_OK;
	}
	
	catch(Exception& x)
	{
		x.DebugPrint();
		
		ARM_RESULT();
	}
}

extern long ARMLOCAL_StochasticVol_LN_Ari_VanillaOption(double f, double K, double T,double drift, double sig, double VolDrift,double VolVol,	double averaging,double reset,int callput, int nbsteps,
										   ARM_result& result
										   )
{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	
	try
	{
		double dResult = Export_StochasticVol_LN_Arithmetic_VanillaOption_with_Reset(f,K,T,drift,sig,VolDrift,VolVol,averaging,reset,callput,nbsteps);
		result.setDouble(dResult);
		return ARM_OK;
	}
	
	catch(Exception& x)
	{
		x.DebugPrint();
		
		ARM_RESULT();
	}
}

extern long ARMLOCAL_StochasticVol_LN_Ari_VanillaOption_Der(int i,double f, double K, double T,double drift, double sig, double VolDrift,double VolVol,	double averaging,double reset,int callput, int nbsteps,
										   ARM_result& result
										   )
{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	
	try
	{
		double dResult = Export_StochasticVol_LN_Arithmetic_VanillaOption_with_Reset(i,f,K,T,drift,sig,VolDrift,VolVol,averaging,reset,callput,nbsteps);
		result.setDouble(dResult);
		return ARM_OK;
	}
	
	catch(Exception& x)
	{
		x.DebugPrint();
		
		ARM_RESULT();
	}
}


extern long ARMLOCAL_StochasticVol_LN_Ari_VanillaOption_Der2(int i,int j,double f, double K, double T,double drift, double sig, double VolDrift,double VolVol,	double averaging,double reset,int callput, int nbsteps,
										   ARM_result& result
										   )
{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	
	try
	{
		double dResult = Export_StochasticVol_LN_Arithmetic_VanillaOption_with_Reset(i,j,f,K,T,drift,sig,VolDrift,VolVol,averaging,reset,callput,nbsteps);
		result.setDouble(dResult);
		return ARM_OK;
	}
	
	catch(Exception& x)
	{
		x.DebugPrint();
		
		ARM_RESULT();
	}
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
///
///       Optimize a Vector of Mean Skew
///
/////////////////////////////////////////////////////////////////////////////////////////////////////////

extern long ARMLOCAL_OptimizeSkewVector(
	const vector<double>&	C_tagetVect,
	const vector<double>&   C_weights,
    const vector<double>&   C_presicions,
	const vector<double>&   C_InitVector,
	const vector<double>&   C_LBoundVector,
	const vector<double>&   C_UBoundVector,
	vector<double>&   C_DataResult,
	const long&				C_algo,
	ARM_result&				result)
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	std::vector<double>& optimizedVector = NULL;
	ARM::QGM_ParameterSet* qgmSet= NULL;

	std::vector<double>& tagetVect = new std::vector<double>(C_tagetVect);
	std::vector<double>& weights = new std::vector<double>(C_weights);
    std::vector<double>& precisions = new std::vector<double>(C_presicions);
	std::vector<double>& InitVector = new std::vector<double>(C_InitVector);
	std::vector<double>& LBoundVector = new std::vector<double>(C_LBoundVector);
	std::vector<double>& UBoundVector = new std::vector<double>(C_UBoundVector);

	try
	{

		qgmSet = ARM::QGM_CalibrateSkew(tagetVect,
            weights,
            precisions,
            InitVector,
            LBoundVector,
            UBoundVector,
            C_algo);
		optimizedVector  = qgmSet->get_X();
		C_DataResult.resize(optimizedVector->size());
		for (int i=0;i<optimizedVector->size();i++)
				C_DataResult[i] = (*optimizedVector )[i];
		return ARM_OK;
	}
	
	catch(Exception& x)
	{
		delete tagetVect;
		x.DebugPrint();
		ARM_RESULT();
	}
}




/////////////////////////////////////////////////////////////////////////////////////////////////////////
///
///       Present value  of a Normal VanillaOption on an underlying with stochastic volatility
///
/////////////////////////////////////////////////////////////////////////////////////////////////////////
extern long ARMLOCAL_StochasticVol_N_VanillaOption(double f, double K, double T,double drift, double sig, double VolDrift,double VolVol,int callput, int nbsteps,
										   ARM_result& result
										   )
{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	
	try
	{
		double dResult = Export_StochasticVol_N_VanillaOption(f,K,T,drift,sig,VolDrift,VolVol,callput,nbsteps);
		result.setDouble(dResult);
		return ARM_OK;
	}
	
	catch(Exception& x)
	{
		x.DebugPrint();
		
		ARM_RESULT();
	}
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
///
///       GLambda Distribution
///
/////////////////////////////////////////////////////////////////////////////////////////////////////////

long ARMLOCAL_GLambda_From_SABR_Calibrate(

	double  C_forward,
	double	C_alpha,
	double	C_beta,
	double	C_rho,
	double	C_nu,
	double	C_T,
	double	C_SabrType,
	double	C_Scope,
	double	C_IniL1,
	double	C_IniL2,
	double	C_IniL3,
	double	C_IniL4,
	double	C_IniL5,
	double	C_IniL6,
	double	C_Algo,
	ARM_result&  result)
{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	std::vector<double>& K_Vec	= NULL;	
	std::vector<double>& ImpVol_Vec	= NULL;	
	
	try
	{	
		GLambda_ParameterSet* setptr = GLambda_CalibrateFromSABR(C_forward,
												C_alpha,
												C_beta,
												C_rho,
												C_nu, 
												C_T,
												C_SabrType,
												120,
												C_Scope,	/// scope=0.2 means we use strikes between f*0.8 and f*1.2
												C_IniL1,
												C_IniL2,
												C_IniL3,
												C_IniL4,
												C_IniL5,
												C_IniL6,
												C_Algo);

		result.setLong(7);
		result.setArray(setptr->get_l1(),0);
		result.setArray(setptr->get_l2(),1);
		result.setArray(setptr->get_l3(),2);
		result.setArray(setptr->get_l4(),3);
		result.setArray(setptr->get_l5(),4);
		result.setArray(setptr->get_l6(),5);
		result.setArray(setptr->get_objective(),6);
		delete setptr;
		
		return ARM_OK;
	}
	catch(Exception& x)
	{
		x.DebugPrint();
		
		ARM_RESULT();
	}
}



extern long ARMLOCAL_Student_GLambda_Power_SpreadOption(
	double l1a,
	double l2a,
	double l3a,
	double l4a,
	double l5a,
	double l6a,

	double l1b,
	double l2b,
	double l3b,
	double l4b,
	double l5b,
	double l6b,

	double copula_corr,
	double copula_degre,

	double a10,
	double b10,
	double k10,
	double a20,
	double b20,
	double k20,
	int n,
	ARM_result& result
	)
	{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	try
	{
		double dResult = Export_Student_GLambda_Power_SpreadOption(	l1a,l2a,l3a,l4a,l5a,l6a,
																	l1b,l2b,l3b,l4b,l5b,l6b,
																	copula_corr,copula_degre,
																	a10,b10,k10,a20,b20,k20,n);
		result.setDouble(dResult);
		return ARM_OK;
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}
}

extern long ARMLOCAL_Student_GLambda_Power_Digital_SpreadOption(
	double l1a,
	double l2a,
	double l3a,
	double l4a,
	double l5a,
	double l6a,

	double l1b,
	double l2b,
	double l3b,
	double l4b,
	double l5b,
	double l6b,

	double copula_corr,
	double copula_degre,

	double a20,
	double b20,
	double k20,
	int n,
	ARM_result& result
	)
	{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	try
	{
		double dResult = Export_Student_GLambda_Power_Digital_SpreadOption(	l1a,l2a,l3a,l4a,l5a,l6a,
																	l1b,l2b,l3b,l4b,l5b,l6b,
																	copula_corr,copula_degre,
																	a20,b20,k20,n);
		result.setDouble(dResult);
		return ARM_OK;
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}
}

extern long ARMLOCAL_Student_GLambda_Power_Index1Digital_SpreadOption(
	double l1a,
	double l2a,
	double l3a,
	double l4a,
	double l5a,
	double l6a,

	double l1b,
	double l2b,
	double l3b,
	double l4b,
	double l5b,
	double l6b,

	double copula_corr,
	double copula_degre,

	double a20,
	double b20,
	double k20,
	int n,
	ARM_result& result
	)
	{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	try
	{
		double dResult = Export_Student_GLambda_Power_Index1Digital_SpreadOption(	l1a,l2a,l3a,l4a,l5a,l6a,
																	l1b,l2b,l3b,l4b,l5b,l6b,
																	copula_corr,copula_degre,
																	a20,b20,k20,n);
		result.setDouble(dResult);
		return ARM_OK;
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}
}

extern long ARMLOCAL_Student_GLambda_Power_Index2Digital_SpreadOption(
	double l1a,
	double l2a,
	double l3a,
	double l4a,
	double l5a,
	double l6a,

	double l1b,
	double l2b,
	double l3b,
	double l4b,
	double l5b,
	double l6b,

	double copula_corr,
	double copula_degre,

	double a20,
	double b20,
	double k20,
	int n,
	ARM_result& result
	)
	{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	try
	{
		double dResult = Export_Student_GLambda_Power_Index2Digital_SpreadOption(	l1a,l2a,l3a,l4a,l5a,l6a,
																	l1b,l2b,l3b,l4b,l5b,l6b,
																	copula_corr,copula_degre,
																	a20,b20,k20,n);
		result.setDouble(dResult);
		return ARM_OK;
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}
}

extern long ARMLOCAL_Student_GLambda_SpreadOption(
	double l1a,
	double l2a,
	double l3a,
	double l4a,
	double l5a,
	double l6a,

	double l1b,
	double l2b,
	double l3b,
	double l4b,
	double l5b,
	double l6b,

	double copula_corr,
	double copula_degre,

	double k20,
	int n,
	ARM_result& result
	)
	{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	try
	{
		double dResult = Export_Student_GLambda_SpreadOption(	l1a,l2a,l3a,l4a,l5a,l6a,
																	l1b,l2b,l3b,l4b,l5b,l6b,
																	copula_corr,copula_degre,
																k20,n);
		result.setDouble(dResult);
		return ARM_OK;
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}
}
////////////////////////////////////////////////////////////////////////////////////////////
///
///			Computation of the Lambert function 
///					
////////////////////////////////////////////////////////////////////////////////////////////


extern long ARMLOCAL_Lambert_Function(
										 double x,
										 ARM_result& result
										 )
{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	
	try
	{
		double dResult = Export_Lambert_Function(x);
		result.setDouble(dResult);
		return ARM_OK;
	}
	
	catch(Exception& x)
	{
		x.DebugPrint();
		
		ARM_RESULT();
	}
}

////////////////////////////////////////////////////////////////////////////////////////
///
///			Black and Sholes Asymptotic time value Formula
///
////////////////////////////////////////////////////////////////////////////////////////


extern long ARMLOCAL_CF_BlackSholesTimeValue(
								 double forward,
								 double totalvolatility,
								 double bondprice,
								 double strike,
								 double CallPut,
								 ARM_result& result
	)
		{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	try
	{
		double dResult = Export_BlackSholesTimeValue(forward,totalvolatility,bondprice,strike,CallPut);
		result.setDouble(dResult);
		return ARM_OK;
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}
}

extern long ARMLOCAL_CF_BlackSholesTimeValue_ImplicitVol(
												double forward,
												double bondprice,
												double strike,
												double CallPut,
												double optprice,
												ARM_result& result
												)
{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	
	try
	{
		double dResult = Export_BlackSholesTimeValue_ImplicitVol(forward,bondprice,strike,CallPut,optprice);
		result.setDouble(dResult);
		return ARM_OK;
	}
	
	catch(Exception& x)
	{
		x.DebugPrint();
		
		ARM_RESULT();
	}
}

extern long ARMLOCAL_CF_ShiftedLogNormal_Quantile(
												double f,
												double k,
												double t,
												double sigma,
												double alpha,
												ARM_result& result
												)
{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	
	try
	{
		double dResult = Export_ShiftedLogNormal_Quantile(f,k,t,sigma,alpha);
		result.setDouble(dResult);
		return ARM_OK;
	}
	
	catch(Exception& x)
	{
		x.DebugPrint();
		
		ARM_RESULT();
	}
}

extern long ARMLOCAL_CF_ShiftedLogNormal_Distribution(
												double f,
												double k,
												double t,
												double sigma,
												double alpha,
												ARM_result& result
												)
{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	
	try
	{
		double dResult = Export_ShiftedLogNormal_Distribution(f,k,t,sigma,alpha);
		result.setDouble(dResult);
		return ARM_OK;
	}
	
	catch(Exception& x)
	{
		x.DebugPrint();
		
		ARM_RESULT();
	}
}

extern long ARMLOCAL_CF_SABR_Quantile(
												double f,
												double k,
												double t,
												double alpha,
												double beta,
												double rho,
												double nu,
												double Sabr_Type,
												double n,
												double alpha_exp,
												double alpha_tanh,
												double kb_tanh,
												ARM_result& result
												)
{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	
	try
	{
		double dResult = Export_SABR_Quantile(f,k,t,alpha,beta,rho,nu,Sabr_Type,n,alpha_exp,alpha_tanh,kb_tanh);
		result.setDouble(dResult);
		return ARM_OK;
	}
	
	catch(Exception& x)
	{
		x.DebugPrint();
		
		ARM_RESULT();
	}
}

extern long ARMLOCAL_CF_SABR_Distribution(
												double f,
												double k,
												double t,
												double alpha,
												double beta,
												double rho,
												double nu,
												double Sabr_Type,
												double n,
												double alpha_exp,
												double alpha_tanh,
												double kb_tanh,
												ARM_result& result
												)
{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	
	try
	{
		double dResult = Export_SABR_Distribution(f,k,t,alpha,beta,rho,nu,Sabr_Type,n,alpha_exp,alpha_tanh,kb_tanh);
		result.setDouble(dResult);
		return ARM_OK;
	}
	
	catch(Exception& x)
	{
		x.DebugPrint();
		
		ARM_RESULT();
	}
}



extern long ARMLOCAL_GLambda_Distribution(
										  double l1,
										  double l2,
										  double l3,
										  double l4,
										  double l5,
										  double l6,
										  double x,
										  ARM_result& result
										  )
{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	
	try
	{
		double dResult = Export_GLambda_Distribution(l1,l2,l3,l4,l5,l6,x);
		result.setDouble(dResult);
		return ARM_OK;
	}
	
	catch(Exception& x)
	{
		x.DebugPrint();
		
		ARM_RESULT();
	}
}

extern long ARMLOCAL_GLambda_Quantile(
										  double l1,
										  double l2,
										  double l3,
										  double l4,
										  double l5,
										  double l6,
										  double x,
										  ARM_result& result
										  )
{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	
	try
	{
		double dResult = Export_GLambda_Quantile(l1,l2,l3,l4,l5,l6,x);
		result.setDouble(dResult);
		return ARM_OK;
	}
	
	catch(Exception& x)
	{
		x.DebugPrint();
		
		ARM_RESULT();
	}
}

extern long ARMLOCAL_Student_Quantile(
										  double deg,
										  double x,
										  ARM_result& result
										  )
{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	
	try
	{
		double dResult = Export_Student_Quantile(deg,x);
		result.setDouble(dResult);
		return ARM_OK;
	}
	
	catch(Exception& x)
	{
		x.DebugPrint();
		
		ARM_RESULT();
	}
}

extern long ARMLOCAL_Student_Distribution(
										  double deg,
										  double x,
										  ARM_result& result
										  )
{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	
	try
	{
		double dResult = Export_Student_Distribution(deg,x);
		result.setDouble(dResult);
		return ARM_OK;
	}
	
	catch(Exception& x)
	{
		x.DebugPrint();
		
		ARM_RESULT();
	}
}

extern long ARMLOCAL_ImcompleteBeta_Inverse(
										  double a,
										  double b,
										  double x,
										  ARM_result& result
										  )
{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	
	try
	{
		double dResult = Export_ImcompleteBeta_Inverse(a,b,x);
		result.setDouble(dResult);
		return ARM_OK;
	}
	
	catch(Exception& x)
	{
		x.DebugPrint();
		
		ARM_RESULT();
	}
}

extern long ARMLOCAL_Student_QIntegral(
										  double a,
										  double b,
										  double x,
										  ARM_result& result
										  )
{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	
	try
	{
		double dResult = Export_Student_QIntegral(a,b,x);
		result.setDouble(dResult);
		return ARM_OK;
	}
	
	catch(Exception& x)
	{
		x.DebugPrint();
		
		ARM_RESULT();
	}
}

extern long ARMLOCAL_Normal_ImpliedVol(
										  double a,
										  double b,
										  double x,
										  double callput,
										  ARM_result& result
										  )
{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	
	try
	{
		double dResult = Export_Normal_ImpliedVol(a,b,x,callput);
		result.setDouble(dResult);
		return ARM_OK;
	}
	
	catch(Exception& x)
	{
		x.DebugPrint();
		
		ARM_RESULT();
	}
}

extern long ARMLOCAL_Normal_Digital_ImpliedVol(
										  double a,
										  double b,
										  double x,
										  double callput,
										  ARM_result& result
										  )
{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	
	try
	{
		double dResult = Export_Normal_Digital_ImpliedVol(a,b,x,callput);
		result.setDouble(dResult);
		return ARM_OK;
	}
	
	catch(Exception& x)
	{
		x.DebugPrint();
		
		ARM_RESULT();
	}
}


extern long ARMLOCAL_Hypergeometric_Whittaker_M(
										  double a,
										  double b,
										  double x,
										  ARM_result& result
										  )
{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	
	try
	{
		double dResult = Export_Hypergeometric_Whittaker_M(a,b,x);
		result.setDouble(dResult);
		return ARM_OK;
	}
	
	catch(Exception& x)
	{
		x.DebugPrint();
		
		ARM_RESULT();
	}
}

extern long ARMLOCAL_Hypergeometric_Whittaker_W(
										  double a,
										  double b,
										  double x,
										  ARM_result& result
										  )
{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	
	try
	{
		double dResult = Export_Hypergeometric_Whittaker_W(a,b,x);
		result.setDouble(dResult);
		return ARM_OK;
	}
	
	catch(Exception& x)
	{
		x.DebugPrint();
		
		ARM_RESULT();
	}
}


extern long ARMLOCAL_Bessel_Y(
										  double a,
										  double x,
										  ARM_result& result
										  )
{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	
	try
	{
		double dResult = Export_Bessel_Y(a,x);
		result.setDouble(dResult);
		return ARM_OK;
	}
	
	catch(Exception& x)
	{
		x.DebugPrint();
		
		ARM_RESULT();
	}
}

extern long ARMLOCAL_Bessel_I(
										  double a,
										  double x,
										  ARM_result& result
										  )
{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	
	try
	{
		double dResult = Export_Bessel_I(a,x);
		result.setDouble(dResult);
		return ARM_OK;
	}
	
	catch(Exception& x)
	{
		x.DebugPrint();
		
		ARM_RESULT();
	}
}

extern long ARMLOCAL_Bessel_J(
										  double a,
										  double x,
										  ARM_result& result
										  )
{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	
	try
	{
		double dResult = Export_Bessel_J(a,x);
		result.setDouble(dResult);
		return ARM_OK;
	}
	
	catch(Exception& x)
	{
		x.DebugPrint();
		
		ARM_RESULT();
	}
}
extern long ARMLOCAL_Bessel_K(
										  double a,
										  double x,
										  ARM_result& result
										  )
{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	
	try
	{
		double dResult = Export_Bessel_K(a,x);
		result.setDouble(dResult);
		return ARM_OK;
	}
	
	catch(Exception& x)
	{
		x.DebugPrint();
		
		ARM_RESULT();
	}
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
///
///        Shifted Heston Vanilla Option Formula
///
/////////////////////////////////////////////////////////////////////////////////////////////////////////



extern long ARMLOCAL_Shifted_Heston_VanillaOption(
										   double F,
										   double K,
										   double sig,
										   double t,
										   double kappa,
										   double theta,
										   double ksi,
										   double rho,
										   double shift,
										   int callorput,
										   int nb1,int nb,int nbS, int nbO,double prec,
										   ARM_result& result
										   )
{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	
	try
	{
		double dResult = Export_Shifted_Heston_VanillaOption(F,K,sig,t,kappa,theta,ksi,rho,shift,callorput,nb1,nb,nbS,nbO,prec);
		result.setDouble(dResult);
		return ARM_OK;
	}
	
	catch(Exception& x)
	{
		x.DebugPrint();
		
		ARM_RESULT();
	}
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
///
///        SABR Heston Vanilla Option Formula
///
/////////////////////////////////////////////////////////////////////////////////////////////////////////



extern long ARMLOCAL_SABR_Heston_VanillaOption(
										   double F,
										   double K,
										   double sig,
										   double t,
										   double kappa,
										   double theta,
										   double ksi,
										   double rho,
										   double beta,
										   int callorput,
										   int nb1,int nb,int nbS, int nbO,double prec,
										   ARM_result& result
										   )
{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	
	try
	{
		double dResult = Export_SABR_Heston_VanillaOption(F,K,sig,t,kappa,theta,ksi,rho,beta,callorput,nb1,nb,nbS,nbO,prec);
		result.setDouble(dResult);
		return ARM_OK;
	}
	
	catch(Exception& x)
	{
		x.DebugPrint();
		
		ARM_RESULT();
	}
}


long ARMLOCAL_GLambda_CompleteSpreadoption(
										   const vector<double>& C_l1a_Vec,
										   const vector<double>& C_l2a_Vec,
										   const vector<double>& C_l3a_Vec,
										   const vector<double>& C_l4a_Vec,
										   const vector<double>& C_l5a_Vec,
										   const vector<double>& C_l6a_Vec,
										   const vector<double>& C_l1b_Vec,
										   const vector<double>& C_l2b_Vec,
										   const vector<double>& C_l3b_Vec,
										   const vector<double>& C_l4b_Vec,
										   const vector<double>& C_l5b_Vec,
										   const vector<double>& C_l6b_Vec,
										   const vector<double>& C_Discount_Vec,
										   double C_copula_corr,
										   double C_copula_degre,
										   double C_k,
										   double C_n,
										   ARM_result&  result)
{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	std::vector<double>&  C_l1a_VVec	= NULL;	
	std::vector<double>&  C_l2a_VVec	= NULL;
	std::vector<double>&  C_l3a_VVec	= NULL;
	std::vector<double>&  C_l4a_VVec	= NULL;
	std::vector<double>&  C_l5a_VVec	= NULL;
	std::vector<double>&  C_l6a_VVec	= NULL;
	std::vector<double>&  C_l1b_VVec	= NULL;
	std::vector<double>&  C_l2b_VVec	= NULL;
	std::vector<double>&  C_l3b_VVec	= NULL;
	std::vector<double>&  C_l4b_VVec	= NULL;
	std::vector<double>&  C_l5b_VVec	= NULL;
	std::vector<double>&  C_l6b_VVec	= NULL;
	std::vector<double>&  C_Dis_VVec	= NULL;	
	try
	{
		C_l1a_VVec	= CreateARMGPVectorFromVECTOR(C_l1a_Vec 		 );
		C_l2a_VVec	= CreateARMGPVectorFromVECTOR(C_l2a_Vec 		 );
		C_l3a_VVec	= CreateARMGPVectorFromVECTOR(C_l3a_Vec 		 );
		C_l4a_VVec	= CreateARMGPVectorFromVECTOR(C_l4a_Vec 		 );
		C_l5a_VVec	= CreateARMGPVectorFromVECTOR(C_l5a_Vec 		 );
		C_l6a_VVec	= CreateARMGPVectorFromVECTOR(C_l6a_Vec 		 );
		C_l1b_VVec	= CreateARMGPVectorFromVECTOR(C_l1b_Vec 		 );
		C_l2b_VVec	= CreateARMGPVectorFromVECTOR(C_l2b_Vec 		 );
		C_l3b_VVec	= CreateARMGPVectorFromVECTOR(C_l3b_Vec 		 );
		C_l4b_VVec	= CreateARMGPVectorFromVECTOR(C_l4b_Vec 		 );
		C_l5b_VVec	= CreateARMGPVectorFromVECTOR(C_l5b_Vec 		 );
		C_l6b_VVec	= CreateARMGPVectorFromVECTOR(C_l6b_Vec 		 );
		C_Dis_VVec	= CreateARMGPVectorFromVECTOR(C_Discount_Vec 	 );		
		
		CCString msg ("");
		
		double dResult = Export_GLambda_CompleteSpreadoption(
			
			C_l1a_VVec,
			C_l2a_VVec,
			C_l3a_VVec,
			C_l4a_VVec,
			C_l5a_VVec,
			C_l6a_VVec,
			C_l1b_VVec,
			C_l2b_VVec,
			C_l3b_VVec,
			C_l4b_VVec,
			C_l5b_VVec,
			C_l6b_VVec,
			C_Dis_VVec,
			C_copula_corr,
			C_copula_degre,
			C_k,
			C_n
			);
		result.setDouble(dResult);
		return ARM_OK;
	}
	
	catch(Exception& x)
	{
		x.DebugPrint();
		
		ARM_RESULT();
	}
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////
///
///        Jump Diffusion Mepi Vanilla Option Formula
///
/////////////////////////////////////////////////////////////////////////////////////////////////////////



extern long ARMLOCAL_JumpDiffusion_Mepi_Call(
											   double   C_P0,
											   double	C_f0,
											   double	C_T,
											   double	C_K,
											   double	C_R,
											   double	C_Emin,
											   double	C_Lmax,
											   double	C_gamma0,
											   double	C_gamma1,
											   double	C_sig,
											   double	C_lambda,
											   double	C_sigJ,
											   double	C_r,
											   double	C_s,
											   double	C_mu,
											   double	C_fees,
											   double	C_volDrift,
											   double	C_volVol,
											   int	C_CallPut,
											   const vector<double>& C_params,
											   ARM_result& result)
{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	std::vector<double>& A_params	= NULL;	
	try
	{
		A_params				= CreateARMGPVectorFromVECTOR( C_params );

		vector<double>* setptr = Export_Mepi_EDP_VanillaOption(
			C_P0,
			C_f0,
			C_T,
			C_K,
			C_R,
			C_Emin,
			C_Lmax,
			C_gamma0,
			C_gamma1,
			C_sig,
			C_lambda,
			C_sigJ,
			C_r,
			C_s,
			C_mu,
			C_fees,
			C_volDrift,
			C_volVol,
			C_CallPut,
			A_params
			);
		int n=setptr->size();
		result.setLong(n);
		for(int i=0;i<n;i++)
		{
			result.setArray((*setptr)[i],i);
		}
		delete setptr;



		return ARM_OK;
	}


	catch(Exception& x)
	{
		x.DebugPrint();
		
		ARM_RESULT();
	}
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
///
///        Trigonal solve
///
/////////////////////////////////////////////////////////////////////////////////////////////////////////



extern long ARMLOCAL_Util_TrigonalSolve(
										const vector<double>& C_A_Vec,
										const vector<double>& C_B_Vec,
										const vector<double>& C_C_Vec,
										const vector<double>& C_R_Vec,
										ARM_result& result)
{
		/// used in the MACRO ARM_RESULT
	CCString msg ("");
	std::vector<double>& A_Vec	= NULL;	
	std::vector<double>& B_Vec	= NULL;	
	std::vector<double>& C_Vec	= NULL;	
	std::vector<double>& R_Vec	= NULL;	
	
	try
	{
		A_Vec				= CreateARMGPVectorFromVECTOR( C_A_Vec );
		B_Vec				= CreateARMGPVectorFromVECTOR( C_B_Vec );
		C_Vec				= CreateARMGPVectorFromVECTOR( C_C_Vec );
		R_Vec				= CreateARMGPVectorFromVECTOR( C_R_Vec );
		
		vector<double>* setptr = Export_Util_Trigonal_Solve(A_Vec,B_Vec,C_Vec,R_Vec);
		int n=setptr->size();
		result.setLong(n);
		for(int i=0;i<n;i++)
		{
			result.setArray((*setptr)[i],i);
		}
		delete A_Vec;
		delete B_Vec;
		delete C_Vec;
		delete R_Vec;
		delete setptr;
		
		return ARM_OK;
	}
	catch(Exception& x)
	{
		x.DebugPrint();
		
		ARM_RESULT();
	}
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
///
///        BiSABR Spread Option Formula
///
/////////////////////////////////////////////////////////////////////////////////////////////////////////



extern long ARMLOCAL_BiSABR_SpreadOption(
											   double   C_F1,
											   double	C_alpha1,
											   double	C_beta1,
											   double	C_rho1,
											   double	C_nu1,
											   double	C_F2,
											   double	C_alpha2,
											   double	C_beta2,
											   double	C_rho2,
											   double	C_nu2,
											   double	C_K,
											   double	C_T,
											   double	C_CallPut,
											   double	C_rhos,
											   double	C_rhov,
											   double	C_rhoc12,
											   double	C_rhoc21,
											   double	C_flag,
											  
											   ARM_result& result)
{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	
	try
	{
		double dResult = Export_BiSABR_SpreadOption(
			C_F1,
			C_alpha1,
			C_beta1,
			C_rho1,
			C_nu1,
			C_F2,
			C_alpha2,
			C_beta2,
			C_rho2,
			C_nu2,
			C_rhos,
			C_rhov,
			C_rhoc12,
			C_rhoc21,
			C_K,
			C_T,
			C_CallPut,
			C_flag
			);
		result.setDouble(dResult);
		return ARM_OK;
	}
	
	catch(Exception& x)
	{
		x.DebugPrint();
		
		ARM_RESULT();
	}
}


extern long ARMLOCAL_Hypergeometric_Appell(
										   double   C_a,
										   double	C_b1,
										   double	C_b2,
										   double	C_c,
										   double	C_x,
										   double	C_xim,
										   double	C_y,
										   double	C_yim,
										   double	C_nb,
										   
										   ARM_result& result)
{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	
	try
	{
		double dResult = Export_Hypergeometric_Appell(
			C_a,
			C_b1,
			C_b2,
			C_c,
			C_x,
			C_xim,
			C_y,
			C_yim,
			C_nb
			);
		result.setDouble(dResult);
		return ARM_OK;
	}
	
	catch(Exception& x)
	{
		x.DebugPrint();
		
		ARM_RESULT();
	}
}

extern long ARMLOCAL_Util_Eigenvalues3(
										const double rho12,
										const double rho13,
										const double rho23,
										ARM_result& result)
{
		/// used in the MACRO ARM_RESULT
	CCString msg ("");
	
	try
	{
		double e1,e2,e3;
		Export_Eigenvalues3(rho12,rho13,rho23,&e1,&e2,&e3);
		result.setLong(3);
		result.setArray(e1,0);
		result.setArray(e2,1);
		result.setArray(e3,2);		
		return ARM_OK;
	}
	catch(Exception& x)
	{
		x.DebugPrint();
		
		ARM_RESULT();
	}
}


extern long ARMLOCAL_Util_Eigenvalues4(
										const double rho12,
										const double rho13,
										const double rho14,
										const double rho23,
										const double rho24,
										const double rho34,
										ARM_result& result)
{
		/// used in the MACRO ARM_RESULT
	CCString msg ("");
	
	try
	{
		double e1,e2,e3,e4;
		Export_Eigenvalues4(rho12,rho13,rho14,rho23,rho24,rho34,&e1,&e2,&e3,&e4);
		result.setLong(4);
		result.setArray(e1,0);
		result.setArray(e2,1);
		result.setArray(e3,2);
		result.setArray(e4,3);

		
		return ARM_OK;
	}
	catch(Exception& x)
	{
		x.DebugPrint();
		
		ARM_RESULT();
	}
}

extern long ARMLOCAL_Util_BiSABR_CorrelationEvolution(
										const double rho1,
										const double rho2,
										const double rhos,
										const double rhov,
										const double rhoc12,
										const double rhoc21,
										const double newrho1,
										const double newrho2,
										ARM_result& result)
{
		/// used in the MACRO ARM_RESULT
	CCString msg ("");
	
	try
	{
		double newrhov,newrhoc12,newrhoc21;
		Export_BiSABR_CorrelationEvolution(rho1,rho2,rhos,
			rhov,rhoc12,rhoc21,
			newrho1,newrho2,
			&newrhov,&newrhoc12,&newrhoc21);
		result.setLong(3);
		result.setArray(newrhov,0);
		result.setArray(newrhoc12,1);
		result.setArray(newrhoc21,2);

		
		return ARM_OK;
	}
	catch(Exception& x)
	{
		x.DebugPrint();
		
		ARM_RESULT();
	}
}


extern long ARMLOCAL_BiSABR_Calibrate(
										const vector<double>& C_F1_vec		, 
										const vector<double>& C_alpha1_vec	,
										const vector<double>& C_beta1_vec	,
										const vector<double>& C_rho1_vec	,
										const vector<double>& C_nu1_vec		,
										const vector<double>& C_F2_vec		,
										const vector<double>& C_alpha2_vec	,
										const vector<double>& C_beta2_vec	,
										const vector<double>& C_rho2_vec	,
										const vector<double>& C_nu2_vec		,
										const vector<double>& C_strike_vec	,
										const vector<double>& C_maturity_vec,
										const vector<double>& C_callput_vec	,
										const vector<double>& C_price_vec	,
										const vector<double>& C_weight_vec	,
										const vector<double>& C_initialparams_vec,
										ARM_result& result)
{
		/// used in the MACRO ARM_RESULT
	CCString msg ("");
	std::vector<double>&		F1_vec					= NULL;	
	std::vector<double>&		alpha1_vec				= NULL;	
	std::vector<double>&		beta1_vec				= NULL;	
	std::vector<double>&		rho1_vec				= NULL;	
	std::vector<double>&		nu1_vec					= NULL;	
	std::vector<double>&		F2_vec					= NULL;	
	std::vector<double>&		alpha2_vec				= NULL;	
	std::vector<double>&		beta2_vec				= NULL;	
	std::vector<double>&		rho2_vec				= NULL;	
	std::vector<double>&		nu2_vec					= NULL;	
	std::vector<double>&		strike_vec				= NULL;	
	std::vector<double>&		maturity_vec			= NULL;	
	std::vector<double>&		callput_vec				= NULL;	
	std::vector<double>&		price_vec				= NULL;	
	std::vector<double>&		weight_vec				= NULL;	
	std::vector<double>&		initialparams_vec		= NULL;	


	
	try
	{
		F1_vec						= CreateARMGPVectorFromVECTOR( C_F1_vec					);
		alpha1_vec					= CreateARMGPVectorFromVECTOR( C_alpha1_vec				);
		beta1_vec					= CreateARMGPVectorFromVECTOR( C_beta1_vec				);
		rho1_vec					= CreateARMGPVectorFromVECTOR( C_rho1_vec				);
		nu1_vec						= CreateARMGPVectorFromVECTOR( C_nu1_vec				);
		F2_vec						= CreateARMGPVectorFromVECTOR( C_F2_vec					);
		alpha2_vec					= CreateARMGPVectorFromVECTOR( C_alpha2_vec				);
		beta2_vec					= CreateARMGPVectorFromVECTOR( C_beta2_vec				);
		rho2_vec					= CreateARMGPVectorFromVECTOR( C_rho2_vec				);
		nu2_vec						= CreateARMGPVectorFromVECTOR( C_nu2_vec				);
		strike_vec					= CreateARMGPVectorFromVECTOR( C_strike_vec				);
		maturity_vec				= CreateARMGPVectorFromVECTOR( C_maturity_vec			);
		callput_vec					= CreateARMGPVectorFromVECTOR( C_callput_vec			);
		price_vec					= CreateARMGPVectorFromVECTOR( C_price_vec				);
		weight_vec					= CreateARMGPVectorFromVECTOR( C_weight_vec				);
		initialparams_vec			= CreateARMGPVectorFromVECTOR( C_initialparams_vec		);

		const vector<double> initialparams_Vec=initialparams_vec->GetValues();
		double rhos_0=initialparams_Vec[0];
		double rhov_0=initialparams_Vec[1];
		double rhoc12_0=initialparams_Vec[2];
		double rhoc21_0=initialparams_Vec[3];
		double converg_prec=initialparams_Vec[4];
		double nbIter_max=initialparams_Vec[5];
		double first_step_max=initialparams_Vec[6];
		double rhos_flag=initialparams_Vec[7];
		int flag=initialparams_Vec[8];
			
		BiSABR_ParameterSet* setptr = Export_BiSABR_CalibrateToSmile(
			F1_vec,
			alpha1_vec,
			beta1_vec,
			rho1_vec,
			nu1_vec,
			F2_vec,
			alpha2_vec,
			beta2_vec,
			rho2_vec,
			nu2_vec,
			strike_vec,
			maturity_vec,
			price_vec,
			weight_vec,
			rhos_0,rhov_0,rhoc12_0,rhoc21_0,
			converg_prec, nbIter_max, first_step_max, rhos_flag,flag);

		
		result.setLong(5);
		result.setArray(setptr->get_rhos(),0);
		result.setArray(setptr->get_rhov(),1);
		result.setArray(setptr->get_rhoc12(),2);
		result.setArray(setptr->get_rhoc21(),3);
		result.setArray(setptr->get_objective(),4);

		delete F1_vec			;		
		delete alpha1_vec		;
		delete beta1_vec		;
		delete rho1_vec			;
		delete nu1_vec			;
		delete F2_vec			;
		delete alpha2_vec		;
		delete beta2_vec		;
		delete rho2_vec			;
		delete nu2_vec			;
		delete strike_vec		;
		delete maturity_vec		;
		delete callput_vec		;
		delete price_vec		;
		delete initialparams_vec;
		delete setptr;
		
		return ARM_OK;
	}
	catch(Exception& x)
	{
		x.DebugPrint();
		
		ARM_RESULT();
	}
}

extern long ARMLOCAL_LN_DigitalOption(
											   double   C_forward,
											   double	C_strike,
											   double	C_maturity,
											   double	C_CallPut,
											   double	C_volatility,
											   ARM_result& result)
{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	
	try
	{
		double dResult = Export_BlackScholesDigitalOption(
			C_forward,
			C_strike,
			C_maturity,
			C_CallPut,
			C_volatility
			);
		result.setDouble(dResult);
		return ARM_OK;
	}
	
	catch(Exception& x)
	{
		x.DebugPrint();
		
		ARM_RESULT();
	}
}

extern long ARMLOCAL_LN_RatioOption(
											   double   C_S1,
											   double	C_Mu1,
											   double	C_Sigma1,
											   double   C_S2,
											   double	C_Mu2,
											   double	C_Sigma2,
											   double	C_Rho,
											   double	C_K,
											   double	C_T,
											   double	C_CallPut,
											   ARM_result& result)
{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	
	try
	{
		double dResult = Export_LN_RatioOption(
			C_S1,
			C_Mu1,
			C_Sigma1,
			C_S2,
			C_Mu2,
			C_Sigma2,
			C_Rho,
			C_K,
			C_T,
			C_CallPut
			);
		result.setDouble(dResult);
		return ARM_OK;
	}
	
	catch(Exception& x)
	{
		x.DebugPrint();
		
		ARM_RESULT();
	}
}

extern long ARMLOCAL_LN_ProductOption(
											   double   C_S1,
											   double	C_Mu1,
											   double	C_Sigma1,
											   double   C_S2,
											   double	C_Mu2,
											   double	C_Sigma2,
											   double	C_Rho,
											   double	C_K,
											   double	C_T,
											   double	C_CallPut,
											   ARM_result& result)
{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	
	try
	{
		double dResult = Export_LN_ProductOption(
			C_S1,
			C_Mu1,
			C_Sigma1,
			C_S2,
			C_Mu2,
			C_Sigma2,
			C_Rho,
			C_K,
			C_T,
			C_CallPut
			);
		result.setDouble(dResult);
		return ARM_OK;
	}
	
	catch(Exception& x)
	{
		x.DebugPrint();
		
		ARM_RESULT();
	}
}

extern long ARMLOCAL_SABR_GaussianSABRDigitalCall(
		double C_f1,
		double C_alpha1,
		double C_beta1,
		double C_rho1,
		double C_nu1,
		double C_flag1,
		double C_f2,
		double C_alpha2,
		double C_beta2,
		double C_rho2,
		double C_nu2,
		double C_flag2,
		double C_rho,
		double C_K,
		double C_T,
		double C_LegendreNb,double C_alpha_exp,double C_alpha_tanh,double C_kb_tanh,
		ARM_result& result)
{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	
	try
	{
		double dResult = Export_GaussianSABRDigitalCall(
			C_f1,
			C_alpha1,
			C_beta1,
			C_rho1,
			C_nu1,
			C_flag1,
			C_f2,
			C_alpha2,
			C_beta2,
			C_rho2,
			C_nu2,
			C_flag2,
			C_rho,
			C_K,
			C_T,
			C_LegendreNb, C_alpha_exp, C_alpha_tanh, C_kb_tanh
			);
		result.setDouble(dResult);
		return ARM_OK;
	}
	
	catch(Exception& x)
	{
		x.DebugPrint();
		
		ARM_RESULT();
	}
}

extern long ARMLOCAL_SABR_GaussianSABRDigitalCallPayingS1(
		double C_f1,
		double C_alpha1,
		double C_beta1,
		double C_rho1,
		double C_nu1,
		double C_flag1,
		double C_f2,
		double C_alpha2,
		double C_beta2,
		double C_rho2,
		double C_nu2,
		double C_flag2,
		double C_rho,
		double C_K,
		double C_T,
		double C_LegendreNb,double C_alpha_exp,double  C_alpha_tanh,double  C_kb_tanh,
		ARM_result& result)
{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	
	try
	{
		double dResult = Export_GaussianSABRDigitalCallPayingS1(
			C_f1,
			C_alpha1,
			C_beta1,
			C_rho1,
			C_nu1,
			C_flag1,
			C_f2,
			C_alpha2,
			C_beta2,
			C_rho2,
			C_nu2,
			C_flag2,
			C_rho,
			C_K,
			C_T,
			C_LegendreNb , C_alpha_exp, C_alpha_tanh, C_kb_tanh
			);
		result.setDouble(dResult);
		return ARM_OK;
	}
	
	catch(Exception& x)
	{
		x.DebugPrint();
		
		ARM_RESULT();
	}
}


extern long ARMLOCAL_SABR_GaussianSABRDigitalCallPayingS2(
		double C_f1,
		double C_alpha1,
		double C_beta1,
		double C_rho1,
		double C_nu1,
		double C_flag1,
		double C_f2,
		double C_alpha2,
		double C_beta2,
		double C_rho2,
		double C_nu2,
		double C_flag2,
		double C_rho,
		double C_K,
		double C_T,
		double C_LegendreNb,double C_alpha_exp,double  C_alpha_tanh,double  C_kb_tanh,
		ARM_result& result)
{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	
	try
	{
		double dResult = Export_GaussianSABRDigitalCallPayingS2(
			C_f1,
			C_alpha1,
			C_beta1,
			C_rho1,
			C_nu1,
			C_flag1,
			C_f2,
			C_alpha2,
			C_beta2,
			C_rho2,
			C_nu2,
			C_flag2,
			C_rho,
			C_K,
			C_T,
			C_LegendreNb, C_alpha_exp, C_alpha_tanh, C_kb_tanh
			);
		result.setDouble(dResult);
		return ARM_OK;
	}
	
	catch(Exception& x)
	{
		x.DebugPrint();
		
		ARM_RESULT();
	}
}


extern long ARMLOCAL_SABR_GaussianSABRDigitalCallPayingS3(
		double C_f1,
		double C_alpha1,
		double C_beta1,
		double C_rho1,
		double C_nu1,
		double C_flag1,
		double C_f2,
		double C_alpha2,
		double C_beta2,
		double C_rho2,
		double C_nu2,
		double C_flag2,
		double C_f3,
		double C_sigma3,
		const vector<double>& C_Correlations_Vec,
		double C_K,
		double C_T,
		const vector<double>& C_Params_Vec,
		ARM_result& result)
{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
		std::vector<double>&		Correlations_vec		= NULL;	
		std::vector<double>&		Params_vec		= NULL;	
	
	try
	{
		Correlations_vec		= CreateARMGPVectorFromVECTOR( C_Correlations_Vec	);
		const vector<double> Correlations_Vec=Correlations_vec->GetValues();
		double C_rho12=Correlations_Vec[0];
		double C_rho13=Correlations_Vec[1];
		double C_rho23=Correlations_Vec[2];
	

		Params_vec		= CreateARMGPVectorFromVECTOR( C_Params_Vec	);
		const vector<double> params_Vec=Params_vec->GetValues();
		int C_LegendreNb=params_Vec[0];
		double C_alpha_exp=params_Vec[1];
		double C_alpha_tanh=params_Vec[2];
		double C_kb_tanh=params_Vec[3];


		double dResult = Export_GaussianSABRDigitalCallPayingS3(
			C_f1,
			C_alpha1,
			C_beta1,
			C_rho1,
			C_nu1,
			C_flag1,
			C_f2,
			C_alpha2,
			C_beta2,
			C_rho2,
			C_nu2,
			C_flag2,
			C_f3,
			C_sigma3,
			C_rho12,
			C_rho13,
			C_rho23,
			C_K,
			C_T,
			C_LegendreNb, C_alpha_exp, C_alpha_tanh, C_kb_tanh
			);
		result.setDouble(dResult);
		return ARM_OK;
	}
	
	catch(Exception& x)
	{
		x.DebugPrint();
		
		ARM_RESULT();
	}
}

////////////////////////////////////////////
//// Function to create a Tarn Proxy
////////////////////////////////////////////
extern long ARMLOCAL_TarnProxy_Create(
	const vector<double>&	C_ResetDates,
	const vector<double>&	C_Fwds,
	const VECTOR<long>&		C_DensityFunctorId,
	const vector<double>&	C_Df,
	const vector<double>&	C_LevPrec,
	const vector<double>&	C_Lev,
	const vector<double>&	C_Fix,
	const vector<double>&	C_Cap,
	const vector<double>&	C_Floor,
	const vector<double>&	C_Fees,
	const vector<double>&	C_dcf,
	const double&			C_Target,
	const bool&				C_Globalcap,
	const bool&				C_Globalfloor,
	const double&			C_CorrelInput,
	const int&				C_NbSimul,
	ARM_result&		result, 
	long			objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");
		

	try
	{
		std::vector<double> rd(C_ResetDates);
		std::vector<double> fwd(C_Fwds);
		std::vector<double> df(C_Df);
		std::vector<double> levprec(C_LevPrec);
		std::vector<double> lev(C_Lev);
		std::vector<double> fix(C_Fix);
		std::vector<double> cap(C_Cap);
		std::vector<double> floor(C_Floor);
		std::vector<double> fees(C_Fees);
		std::vector<double> dcf(C_dcf);

		vector<ARM_DensityFunctor*> densityFunctor( C_DensityFunctorId.size() );
		ARM_DensityFunctor* functor;

		for(size_t i=0; i<C_DensityFunctorId.size(); ++i )
		{
			if (C_DensityFunctorId[i] != ARM_NULL_OBJECT)
			{
				functor = dynamic_cast<ARM_DensityFunctor*>(LOCAL_PERSISTENT_OBJECTS->GetObject(C_DensityFunctorId[i]));

			    if (!functor)
				{
					result.setMsg ("ARM_ERR: functor is not of good type!");
					return ARM_KO;
				}
				
				densityFunctor[i] = static_cast<ARM_DensityFunctor*> (functor->Clone()) ;
			}

		}

		ARM_TarnProxy* proxy=new ARM_TarnProxy(rd,df,levprec,lev,fix,cap,floor,fees,dcf,C_Target,C_Globalcap,C_Globalfloor);
		proxy->Build(fwd,densityFunctor);
		proxy->Price(C_CorrelInput,C_NbSimul);

		/// assign object
		if( !assignObject( proxy, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }

	}

	catch(Exception& x)
	{
		x.DebugPrint();
		ARM_RESULT();
	}

}

long ARMLOCAL_TarnProxy_GetPrice(
	const long&	C_TarnProxyId, 
	const double& C_NbPayoff, 
	ARM_result&	result)
{
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	try
	{
		ARM_Object* armObj =  LOCAL_PERSISTENT_OBJECTS->GetObject(C_TarnProxyId);
		ARM_TarnProxy* proxy = NULL;
		if( !(proxy = dynamic_cast< ARM_TarnProxy* >( armObj )) )
		{
			result.setMsg ("ARM_ERR: proxy is not of a good type");
			return ARM_KO;
		};
		
		double price = proxy->GetPrice(C_NbPayoff);
		result.setDouble(price);
		return ARM_OK;
	}

	catch(Exception& x)
	{
		x.DebugPrint();
		ARM_RESULT();
	}
}


long ARMLOCAL_VBMinMaxProxy_Create(
	const double&			C_AsOf,
	const vector<double>&	C_resetDates,
	const vector<double>&	C_fwdRates,
	const vector<double>&	C_totalVol,
	const vector<double>&	C_leftVol,
	const vector<double>&	C_rightVol,
	const vector<double>&	C_nu,
	const vector<double>&	C_rho,
	const int&				C_nbSimul,
	const bool&				C_sabrDiff,
	const int&				C_typeprice,
	const int&				C_type1sens,
	const int&				C_type2sens,
	const vector<double>&	C_rate1Lev,
	const vector<double>&	C_rate1Add,
	const vector<double>&	C_capRateLev,
	const vector<double>&	C_capRateAdd,
	const vector<double>&	C_vbLev,
	const int&				C_maxChoice,
	const int&				C_minChoice,
	const double&			C_minmaxFreq,
	ARM_result&				result,
	long					objId)
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");
		

	try
	{
		std::vector<double> rd(C_resetDates);
		std::vector<double> fwd(C_fwdRates);
		std::vector<double> totalvol(C_totalVol);
		std::vector<double> leftvol(C_leftVol);
		std::vector<double> rightvol(C_rightVol);
		std::vector<double> rho(C_rho);
		std::vector<double> nu(C_nu);
		std::vector<double> rate1Lev(C_rate1Lev);
		std::vector<double> rate1Add(C_rate1Add);
		std::vector<double> caplev(C_capRateLev);
		std::vector<double> capadd(C_capRateAdd);
		std::vector<double> vblev(C_vbLev);

		ARM_VBMinMaxProxy * proxy = new ARM_VBMinMaxProxy(C_AsOf, rd, fwd, totalvol, leftvol, rightvol, nu, rho, 
			C_nbSimul, C_sabrDiff, C_typeprice, C_type1sens, C_type2sens, rate1Lev, rate1Add, caplev, capadd, vblev,
			C_maxChoice,C_minChoice, C_minmaxFreq);

		/// assign object
		if( !assignObject( proxy, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }

	}

	catch(Exception& x)
	{
		x.DebugPrint();
		ARM_RESULT();
	}
}

long ARMLOCAL_VBMinMaxProxy_GetInfo(
	const long& VBMinMaxProxyId,
	const int& info,
	vector<double>& vresult,
	ARM_result& result)
{
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	try
	{
		ARM_Object* armObj =  LOCAL_PERSISTENT_OBJECTS->GetObject(VBMinMaxProxyId);
		ARM_VBMinMaxProxy* proxy = NULL;
		if( !(proxy = dynamic_cast< ARM_VBMinMaxProxy* >( armObj )) )
		{
			result.setMsg ("ARM_ERR: proxy is not of a good type");
			return ARM_KO;
		};
		 
		if(info == 1) // STD VB
		{
			vresult.resize(proxy->GetStdVB1().size());
			for(int k = 0; k < vresult.size(); k++) vresult[k] = proxy->GetStdVB1(k);
		}
		else if(info == 2) 
		{
			vresult.resize(proxy->GetStdVB2().size());
			for(int k = 0; k < vresult.size(); k++) vresult[k] = proxy->GetStdVB2(k);
		}
		else if(info == 3) 
		{
			vresult.resize(proxy->GetStdVB3().size());
			for(int k = 0; k < vresult.size(); k++) vresult[k] = proxy->GetStdVB3(k);
		}
		else if(info == 31) // MAXRATE
		{
			vresult.resize(proxy->GetMaxRate().size());
			for(int k = 0; k < vresult.size(); k++) vresult[k] = proxy->GetMaxRate()[k];
		}
		else if(info == 32) // MINRATE
		{
			vresult.resize(proxy->GetMinRate().size());
			for(int k = 0; k < vresult.size(); k++) vresult[k] = proxy->GetMinRate()[k];
		}
		else if(info == 4)
		{
			vresult.resize(proxy->GetCoupon().size());
			for(int k = 0; k < vresult.size(); k++) vresult[k] = proxy->GetCoupon()[k];
		}
		else
		{
			result.setMsg("ARM_ERR : info should price, maxrate or minrate");
			return ARM_KO;
		}

		return ARM_OK;
	}

	catch(Exception& x)
	{
		x.DebugPrint();
		ARM_RESULT();
	}
}

extern long ARMLOCAL_Berm2DatesProxy_Create(
	const double&			C_AsOf,
	const vector<double>&	C_resetDates,
	const vector<double>&	C_fwdRates,
	const double&			C_strike,
	const vector<double>&	C_DFs,
	const vector<double>&	C_fwdRatesVol,
	const vector<double>&	C_fwdRatesPartVol,
	const vector<double>&	C_vols,
	const vector<double>&	C_volvols,
	const vector<double>&	C_rho,
	const double&			C_Rho12,
	const int&				nbSimul,
	const int&				C_TypeDiff,
	ARM_result&				result,
	long					objId)
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");
		

	try
	{
		std::vector<double> rd(C_resetDates);
		std::vector<double> fwd(C_fwdRates);
		std::vector<double> dfs(C_DFs);
		std::vector<double> fwdvol(C_fwdRatesVol);
		std::vector<double> vols(C_vols);
		std::vector<double> vvols(C_volvols);
		std::vector<double> rho(C_rho);
		std::vector<double> fwdpartvol(C_fwdRatesPartVol);

		ARM_Berm2DatesProxy * proxy = new ARM_Berm2DatesProxy(C_AsOf, rd, fwd, C_strike, dfs, fwdvol, fwdpartvol, vols, vvols, rho, C_Rho12, nbSimul, C_TypeDiff);

		/// assign object
		if( !assignObject( proxy, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }

	}

	catch(Exception& x)
	{
		x.DebugPrint();
		ARM_RESULT();
	}
}

extern long ARMLOCAL_Berm2DatesProxy_GetPrice(
	const long& Berm2DatesProxyId,
	const int& info,
	ARM_result& result)
{
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	try
	{
		ARM_Object* armObj =  LOCAL_PERSISTENT_OBJECTS->GetObject(Berm2DatesProxyId);
		ARM_Berm2DatesProxy * proxy = NULL;
		if( !(proxy = dynamic_cast< ARM_Berm2DatesProxy* >( armObj )) )
		{
			result.setMsg ("ARM_ERR: proxy is not of a good type");
			return ARM_KO;
		};
		
		double price;
		if(info == 0)
			price = proxy->GetPrice();
		else if(info == 1)
			price = proxy->GetEuro(0);
		else if(info == 2)
			price = proxy->GetEuro(1);

		result.setDouble(price);
		return ARM_OK;
	}

	catch(Exception& x)
	{
		x.DebugPrint();
		ARM_RESULT();
	}
}


	


extern long ARMLOCAL_SpreadVBProxy_Create(
	const double&			C_AsOf,
	const vector<double>&	C_resetDates,
	const vector<double>&	C_fwdRates1,
	const vector<double>&	C_fwdRatesVol1,
	const vector<double>&	C_fwdRatesPartVol1,
	const vector<double>&	C_vols1,
	const vector<double>&	C_nu1,
	const vector<double>&	C_rho1,
	const vector<double>&	C_fwdRates2,
	const vector<double>&	C_fwdRatesVol2,
	const vector<double>&	C_fwdRatesPartVol2,
	const vector<double>&	C_vols2,
	const vector<double>&	C_nu2,
	const vector<double>&	C_rho2,
	const vector<double>&	C_Rate1Rate2Correl,
	const int&				nbSimul,
	const bool&				sabrDiff,
	const vector<double>&	C_fwdLev,
	const vector<double>&	C_fwdStrikes,
	const int&				C_typeprice,
	const int&				C_type1sens,
	const int&				C_type2sens,
	const double&			C_CorrMeanRev,
	const double&			C_CorrVol,
	const vector<double>&	C_Levier,
	const vector<double>&	C_Fixed,
	ARM_result&				result,
	long					objId)
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");
		

	try
	{
		std::vector<double> rd(C_resetDates);
		std::vector<double> fwd1(C_fwdRates1);
		std::vector<double> fwdvol1(C_fwdRatesVol1);
		std::vector<double> vols1(C_vols1);
		std::vector<double> nu1(C_nu1);
		std::vector<double> rho1(C_rho1);
		std::vector<double> fwdpartvol1(C_fwdRatesPartVol1);
		std::vector<double> fwd2(C_fwdRates2);
		std::vector<double> fwdvol2(C_fwdRatesVol2);
		std::vector<double> vols2(C_vols2);
		std::vector<double> nu2(C_nu2);
		std::vector<double> rho2(C_rho2);
		std::vector<double> fwdpartvol2(C_fwdRatesPartVol2);
		std::vector<double> rcorrel(C_Rate1Rate2Correl);
		std::vector<double> fwdLev(C_fwdLev);
		std::vector<double> strikes(C_fwdStrikes);
		std::vector<double> levier(C_Levier);
		std::vector<double> fix(C_Fixed);

		ARM_SpreadVBProxy * proxy = new ARM_SpreadVBProxy(C_AsOf, rd,fwd1,fwd2,
			fwdvol1,fwdpartvol1,vols1,nu1,rho1,fwdvol2,fwdpartvol2,vols2,nu2,rho2,rcorrel,
			nbSimul,fwdLev, strikes,C_typeprice,C_type1sens,C_type2sens,sabrDiff,C_CorrMeanRev,C_CorrVol,
			levier,fix);

		/// assign object
		if( !assignObject( proxy, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }

	}

	catch(Exception& x)
	{
		x.DebugPrint();
		ARM_RESULT();
	}
}

extern long ARMLOCAL_SpreadVBProxy_GetInfo(
	const long& SpreadVBProxyId,
	const int& info,
	vector<double>&	vresult,
	ARM_result& result
	)
{
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	try
	{
		ARM_Object* armObj =  LOCAL_PERSISTENT_OBJECTS->GetObject(SpreadVBProxyId);
		ARM_SpreadVBProxy * proxy = NULL;
		if( !(proxy = dynamic_cast< ARM_SpreadVBProxy* >( armObj )) )
		{
			result.setMsg ("ARM_ERR: proxy is not of a good type");
			return ARM_KO;
		};
		 
		if(info == 1) // STD VB
		{
			vresult.resize(proxy->GetStdVB1().size());
			for(int k = 0; k < vresult.size(); k++) vresult[k] = proxy->GetStdVB1(k);
		}
		else if(info == 2) 
		{
			vresult.resize(proxy->GetStdVB2().size());
			for(int k = 0; k < vresult.size(); k++) vresult[k] = proxy->GetStdVB2(k);
		}
		else if(info == 3)
		{
			vresult.resize(proxy->GetCoupon().size());
			for(int k = 0; k < vresult.size(); k++) vresult[k] = proxy->GetCoupon(k);
		}
		else
		{
			result.setMsg("ARM_ERR : info should stdvb1 or stdvb2");
			return ARM_KO;
		}

		return ARM_OK;
	}

	catch(Exception& x)
	{
		x.DebugPrint();
		ARM_RESULT();
	}
}





/////////////////////////////////////////////////////////////////////////////////////////////////////////
///
///        BiSABR Digital Spread Option Formula
///
/////////////////////////////////////////////////////////////////////////////////////////////////////////



extern long ARMLOCAL_BiSABR_Digital_SpreadOption(
											   double   C_F1,
											   double	C_alpha1,
											   double	C_beta1,
											   double	C_rho1,
											   double	C_nu1,
											   double	C_F2,
											   double	C_alpha2,
											   double	C_beta2,
											   double	C_rho2,
											   double	C_nu2,
											   double	C_K,
											   double	C_T,
											   double	C_CallPut,
											   double	C_rhos,
											   double	C_rhov,
											   double	C_rhoc12,
											   double	C_rhoc21,
											   double	C_flag,
											  
											   ARM_result& result)
{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	
	try
	{
		double dResult = Export_BiSABR_Digital_SpreadOption(
			C_F1,
			C_alpha1,
			C_beta1,
			C_rho1,
			C_nu1,
			C_F2,
			C_alpha2,
			C_beta2,
			C_rho2,
			C_nu2,
			C_K,
			C_T,
			C_CallPut,
			C_rhos,
			C_rhov,
			C_rhoc12,
			C_rhoc21,
			C_flag
			);
		result.setDouble(dResult);
		return ARM_OK;
	}
	
	catch(Exception& x)
	{
		x.DebugPrint();
		
		ARM_RESULT();
	}
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
///
///        BiSABR Digital Pays S1 Spread Option Formula
///
/////////////////////////////////////////////////////////////////////////////////////////////////////////



extern long ARMLOCAL_BiSABR_DigitalPaysS1_SpreadOption(
											   double   C_F1,
											   double	C_alpha1,
											   double	C_beta1,
											   double	C_rho1,
											   double	C_nu1,
											   double	C_F2,
											   double	C_alpha2,
											   double	C_beta2,
											   double	C_rho2,
											   double	C_nu2,
											   double	C_K,
											   double	C_T,
											   double	C_CallPut,
											   double	C_rhos,
											   double	C_rhov,
											   double	C_rhoc12,
											   double	C_rhoc21,
											   double	C_flag,
											  
											   ARM_result& result)
{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	
	try
	{
		double dResult = Export_BiSABR_Digital_SpreadOption_PayS1(
			C_F1,
			C_alpha1,
			C_beta1,
			C_rho1,
			C_nu1,
			C_F2,
			C_alpha2,
			C_beta2,
			C_rho2,
			C_nu2,
			C_K,
			C_T,
			C_CallPut,
			C_rhos,
			C_rhov,
			C_rhoc12,
			C_rhoc21,
			C_flag
			);
		result.setDouble(dResult);
		return ARM_OK;
	}
	
	catch(Exception& x)
	{
		x.DebugPrint();
		
		ARM_RESULT();
	}
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
///
///        BiSABR Digital Pays S2 Spread Option Formula
///
/////////////////////////////////////////////////////////////////////////////////////////////////////////



extern long ARMLOCAL_BiSABR_DigitalPaysS2_SpreadOption(
											   double   C_F1,
											   double	C_alpha1,
											   double	C_beta1,
											   double	C_rho1,
											   double	C_nu1,
											   double	C_F2,
											   double	C_alpha2,
											   double	C_beta2,
											   double	C_rho2,
											   double	C_nu2,
											   double	C_K,
											   double	C_T,
											   double	C_CallPut,
											   double	C_rhos,
											   double	C_rhov,
											   double	C_rhoc12,
											   double	C_rhoc21,
											   double	C_flag,
											  
											   ARM_result& result)
{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	
	try
	{
		double dResult = Export_BiSABR_Digital_SpreadOption_PayS2(
			C_F1,
			C_alpha1,
			C_beta1,
			C_rho1,
			C_nu1,
			C_F2,
			C_alpha2,
			C_beta2,
			C_rho2,
			C_nu2,
			C_K,
			C_T,
			C_CallPut,
			C_rhos,
			C_rhov,
			C_rhoc12,
			C_rhoc21,
			C_flag
			);
		result.setDouble(dResult);
		return ARM_OK;
	}
	
	catch(Exception& x)
	{
		x.DebugPrint();
		
		ARM_RESULT();
	}
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
///
///        BiSABR Digital Pays S3 Spread Option Formula
///
/////////////////////////////////////////////////////////////////////////////////////////////////////////





extern long ARMLOCAL_BiSABR_DigitalPaysS3_SpreadOption(
											   double   C_F1,
											   double	C_alpha1,
											   double	C_beta1,
											   double	C_rho1,
											   double	C_nu1,
											   double	C_F2,
											   double	C_alpha2,
											   double	C_beta2,
											   double	C_rho2,
											   double	C_nu2,
											   double	C_rhos,
											   double	C_rhov,
											   double	C_rhoc12,
											   double	C_rhoc21,
											   const vector<double>& C_S3params_vec,
											   double	C_K,
											   double	C_T,
											   double	C_CallPut,
											   double	C_flag,
											  
											   ARM_result& result)
{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	std::vector<double>&		S3params_vec		= NULL;	
	try
	{
		S3params_vec			= CreateARMGPVectorFromVECTOR( C_S3params_vec		);
		const vector<double> S3params_Vec=S3params_vec->GetValues();
		double S3=S3params_Vec[0];
		double sigma3=S3params_Vec[1];
		double rho13=S3params_Vec[2];
		double rho23=S3params_Vec[3];
		double dResult = Export_BiSABR_Digital_SpreadOption_PayS3(
			C_F1,
			C_alpha1,
			C_beta1,
			C_rho1,
			C_nu1,
			C_F2,
			C_alpha2,
			C_beta2,
			C_rho2,
			C_nu2,
			C_rhos,
			C_rhov,
			C_rhoc12,
			C_rhoc21,
			S3,
			sigma3,
			rho13,
			rho23,
			C_K,
			C_T,
			C_CallPut,
			C_flag
			);
		result.setDouble(dResult);
		return ARM_OK;
	}
	
	catch(Exception& x)
	{
		x.DebugPrint();
		
		ARM_RESULT();
	}
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
///
///        BiSABR Quantile and distribution
///
/////////////////////////////////////////////////////////////////////////////////////////////////////////





extern long ARMLOCAL_BiSABR_Distribution(
											   double   C_F1,
											   double	C_alpha1,
											   double	C_beta1,
											   double	C_rho1,
											   double	C_nu1,
											   double	C_F2,
											   double	C_alpha2,
											   double	C_beta2,
											   double	C_rho2,
											   double	C_nu2,
											   double	C_rhos,
											   double	C_rhov,
											   double	C_rhoc12,
											   double	C_rhoc21,
											   double	C_K,
											   double	C_T,
											   double	C_flag,
											  
											   ARM_result& result)
{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	try
	{
		double dResult = Export_BiSABR_Distribution(
			C_F1,
			C_alpha1,
			C_beta1,
			C_rho1,
			C_nu1,
			C_F2,
			C_alpha2,
			C_beta2,
			C_rho2,
			C_nu2,
			C_rhos,
			C_rhov,
			C_rhoc12,
			C_rhoc21,
			C_K,
			C_T,
			C_flag
			);
		result.setDouble(dResult);
		return ARM_OK;
	}
	
	catch(Exception& x)
	{
		x.DebugPrint();
		
		ARM_RESULT();
	}
}

extern long ARMLOCAL_BiSABR_Quantile(
											   double   C_F1,
											   double	C_alpha1,
											   double	C_beta1,
											   double	C_rho1,
											   double	C_nu1,
											   double	C_F2,
											   double	C_alpha2,
											   double	C_beta2,
											   double	C_rho2,
											   double	C_nu2,
											   double	C_rhos,
											   double	C_rhov,
											   double	C_rhoc12,
											   double	C_rhoc21,
											   double	C_K,
											   double	C_T,
											   double	C_flag,
											  
											   ARM_result& result)
{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	try
	{
		double dResult = Export_BiSABR_Quantile(
			C_F1,
			C_alpha1,
			C_beta1,
			C_rho1,
			C_nu1,
			C_F2,
			C_alpha2,
			C_beta2,
			C_rho2,
			C_nu2,
			C_rhos,
			C_rhov,
			C_rhoc12,
			C_rhoc21,
			C_K,
			C_T,
			C_flag
			);
		result.setDouble(dResult);
		return ARM_OK;
	}
	
	catch(Exception& x)
	{
		x.DebugPrint();
		
		ARM_RESULT();
	}
}

extern long ARMLOCAL_Shifted2LogNormal_Distribution(
											   double   C_F1,
											   double	C_sigma1,
											   double	C_F2,
											   double	C_sigma2,
											   double	C_alpha,
											   double	C_rho,
											   double	C_K,
											   double	C_T,
											   double	C_n,
											   ARM_result& result)
{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	try
	{
		double dResult = Export_Shifted2LogNormal_Distribution(
			C_F1,
			C_sigma1,
			C_F2,
			C_sigma2,
			C_alpha,
			C_rho,
			C_K,
			C_T,
			C_n
			);
		result.setDouble(dResult);
		return ARM_OK;
	}
	
	catch(Exception& x)
	{
		x.DebugPrint();
		
		ARM_RESULT();
	}
}

extern long ARMLOCAL_Shifted2LogNormal_Quantile(
											   double   C_F1,
											   double	C_sigma1,
											   double	C_F2,
											   double	C_sigma2,
											   double	C_alpha,
											   double	C_rho,
											   double	C_K,
											   double	C_T,
											   double	C_n,
											   ARM_result& result)
{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	try
	{
		double dResult = Export_Shifted2LogNormal_Quantile(
			C_F1,
			C_sigma1,
			C_F2,
			C_sigma2,
			C_alpha,
			C_rho,
			C_K,
			C_T,
			C_n
			);
		result.setDouble(dResult);
		return ARM_OK;
	}
	catch(Exception& x)
	{
		x.DebugPrint();
		
		ARM_RESULT();
	}
}

extern long ARMLOCAL_BetaEqualZeroSABR(
											   double   C_F,
											   double	C_K,
											   double	C_T,
											   double	C_mu,
											   double	C_alpha,
											   double	C_rho,
											   double	C_nu,
											   double	C_CallPut,
											   ARM_result& result)
{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	try
	{
		double dResult = Export_SABR_BetaEqualZero_Option(
			C_F,
			C_K,
			C_T,
			C_mu,
			C_alpha,
			C_rho,
			C_nu,
			C_CallPut
			);
		result.setDouble(dResult);
		return ARM_OK;
	}
	
	catch(Exception& x)
	{
		x.DebugPrint();
		
		ARM_RESULT();
	}
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////
///
///        BiSABR : Local_BiSABR_S3_SpreadOption  pays A1*(S1-S2) +B1*S3 -K1 if A2*(S1-S2) +B2*S3 -K2 >0
///
/////////////////////////////////////////////////////////////////////////////////////////////////////////





extern long ARMLOCAL_BiSABR_S3_SpreadOption(
											const vector<double>& C_S1params_vec,
											const vector<double>& C_S2params_vec,
											const vector<double>& C_S3params_vec,
											double	C_rhos,
											double	C_rhov,
											double	C_rhoc12,
											double	C_rhoc21,
											double	C_Correlation,
											double	C_T,
											double	C_A1,
											double	C_B1,
											double	C_K1,
											double	C_A2,
											double	C_B2,
											double	C_K2,
											double	C_flag,
											double  C_nbsteps,
											ARM_result& result)
{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	std::vector<double>&		S1params_vec		= NULL;	
	std::vector<double>&		S2params_vec		= NULL;	
	std::vector<double>&		S3params_vec		= NULL;	
	
	try
	{
		S1params_vec			= CreateARMGPVectorFromVECTOR( C_S1params_vec		);
		const vector<double> S1params_Vec=S1params_vec->GetValues();
		double F1=S1params_Vec[0];
		double alpha1=S1params_Vec[1];
		double beta1=S1params_Vec[2];
		double rho1=S1params_Vec[3];
		double nu1=S1params_Vec[4];
		
		S2params_vec			= CreateARMGPVectorFromVECTOR( C_S2params_vec		);
		const vector<double> S2params_Vec=S2params_vec->GetValues();
		double F2=S2params_Vec[0];
		double alpha2=S2params_Vec[1];
		double beta2=S2params_Vec[2];
		double rho2=S2params_Vec[3];
		double nu2=S2params_Vec[4];
		
		S3params_vec			= CreateARMGPVectorFromVECTOR( C_S3params_vec		);
		const vector<double> S3params_Vec=S3params_vec->GetValues();
		double F3=S3params_Vec[0];
		double alpha3=S3params_Vec[1];
		double beta3=S3params_Vec[2];
		double rho3=S3params_Vec[3];
		double nu3=S3params_Vec[4];

		double dResult = Export_BiSABR_S3_SpreadOption(
			 F1,
			 alpha1,
			 beta1,
			 rho1,
			 nu1,
								 
			 F2,
			 alpha2,
			 beta2,
			 rho2,
			 nu2,
								 
			 F3,
			 alpha3,
			 beta3,
			 rho3,
			 nu3,
			 
			 C_rhos,
			 C_rhov,
			 C_rhoc12,
			 C_rhoc21,
			 C_Correlation,
			 C_T,
			 C_A1,
			 C_B1,
			 C_K1,
			 C_A2,
			 C_B2,
			 C_K2,
			 C_flag,
			 C_nbsteps);
		result.setDouble(dResult);
		return ARM_OK;
	}
	
	catch(Exception& x)
	{
		x.DebugPrint();
		
		ARM_RESULT();
	}
}


extern long ARMLOCAL_Heston_OptionPrice(
												double	C_AsOfDate,
												double	C_ResetDate,
												double	C_Forward,
												double	C_Strike,
												int		C_CallOrPut,
												double	C_V0,
												double	C_Kappa,
												double	C_Rho,
												double	C_Theta,
												double	C_VVol,
												double	C_Shift,
												const vector<double>& C_Times,
												const vector<double>& C_Levels,
												ARM_result& result )
{
	CCString msg("");

	try
	{
		std::vector<double> times(C_Times);
		std::vector<double> levels(C_Levels);
		times -= C_AsOfDate;
		times /= 365.;

		double price = Export_Heston_OptionPrice((C_ResetDate - C_AsOfDate) / 365.,
												 C_Forward,
												 C_Strike,
												 C_CallOrPut,
												 C_V0,
												 C_Kappa,
												 C_Theta,
												 C_VVol,
												 C_Rho,
												 C_Shift,
												 times,
												 levels);

		result.setDouble(price);

		return ARM_OK;
	}
	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}
}

extern long ARMLOCAL_Heston2b_OptionPrice(
												double	C_AsOfDate,
												double	C_ResetDate,
												double	C_Forward,
												double	C_Strike,
												int		C_CallOrPut,
												double	C_V01,
												double	C_Kappa1,
												double	C_Rho1,
												double	C_Theta1,
												double	C_VVol1,
												double	C_V02,
												double	C_Kappa2,
												double	C_Rho2,
												double	C_Theta2,
												double	C_VVol2,
												double	C_Shift,
												const vector<double>& C_Times,
												const vector<double>& C_Levels,
												ARM_result& result )
{
	CCString msg("");

	try
	{
		std::vector<double> times(C_Times);
		std::vector<double> levels(C_Levels);
		times -= C_AsOfDate;
		times /= 365.;

		double price = Export_Heston2b_OptionPrice((C_ResetDate - C_AsOfDate) / 365.,
												 C_Forward,
												 C_Strike,
												 C_CallOrPut,
												 C_V01,
												 C_Kappa1,
												 C_Theta1,
												 C_VVol1,
												 C_Rho1,
												 C_V02,
												 C_Kappa2,
												 C_Theta2,
												 C_VVol2,
												 C_Rho2,
												 C_Shift,
												 times,
												 levels);

		result.setDouble(price);

		return ARM_OK;
	}
	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}
}

extern long ARMLOCAL_MixteHeston_OptionPrice(
	const double&	C_asOf,
	const double&	C_resetDate,
	const double&	C_forward,
	const double&	C_strike,
	const int&		C_callPut,
	const double&	C_sigma,
	const double&	C_v0,
	const double&	C_kappa,
	const double&	C_rho,
	const double&	C_theta,
	const double&	C_vvol,
	const double&	C_shift,
	const vector<double>& C_times,
	const vector<double>& C_levels,
	ARM_result& result )
{
	CCString msg("");

	try
	{
		std::vector<double> times(C_times);
		std::vector<double> levels(C_levels);
		times -= C_asOf;
		times /= 365.;

		double price = Export_MixteHeston_OptionPrice((C_resetDate - C_asOf) / 365.,
												 C_forward,
												 C_strike,
												 C_callPut,
												 C_sigma,
												 C_v0,
												 C_kappa,
												 C_theta,
												 C_vvol,
												 C_rho,
												 C_shift,
												 times,
												 levels);

		result.setDouble(price);

		return ARM_OK;
	}
	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}
}


extern long ARMLOCAL_Normal_Heston_VanillaCall(
											   double   C_rho,
											   double	C_lambdaV,
											   double	C_thetaV,
											   double	C_kappaV,
											   double	C_V0,
											   double	C_S0,
											   double	C_k,
											   double	C_T,
											   double	C_lambdaB,
											   double	C_callput,
											   double	C_nbfirst,
											   double	C_nb,
											   double	C_NbStage,
											   double	C_NbOscill,
											   double	C_prec,
											   ARM_result& result)
{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	try
	{
		double dResult = Export_Normal_Heston_VanillaOption(
			C_rho,
			C_lambdaV,
			C_thetaV,
			C_kappaV,
			C_V0,
			C_S0,
			C_k,
			C_T,
			C_lambdaB,
			C_callput,
			C_nbfirst,
			C_nb,
			C_NbStage,
			C_NbOscill,
			C_prec);
		result.setDouble(dResult);
		return ARM_OK;
	}
	
	catch(Exception& x)
	{
		x.DebugPrint();
		
		ARM_RESULT();
	}
}

extern long ARMLOCAL_SuperNormal_Heston_VanillaCall(
	const double& C_rho1,
	const double& C_lambda1,
	const double& C_theta1,
	const double& C_kappa1,
	const double& C_V01,
	const double& C_rho2,
	const double& C_lambda2,
	const double& C_theta2,
	const double& C_kappa2,
	const double& C_V02,
	const double& C_S0,
	const double& C_k,
	const double& C_T,
	const double& C_CallPut,
	const double& C_nb,
	ARM_result& result)
{
	CCString msg("");

	try
	{
		double res = SuperNormalHeston(C_rho1,C_lambda1,C_theta1,C_kappa1,C_V01,
									   C_rho2,C_lambda2,C_theta2,C_kappa2,C_V02,
									   C_S0,C_k,C_T,C_CallPut,C_nb);

		result.setDouble(res);
		return ARM_OK;
	}
	catch(Exception& x)
	{
		x.DebugPrint();
		
		ARM_RESULT();
	}

}


extern long ARMLOCAL_SABR_To_Heston_SmileCalibration_GetValue(
	const long& CalibrationId,
	const int& info,
	ARM_result& result
	)
{
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	try
	{
		ARM_Object* armObj =  LOCAL_PERSISTENT_OBJECTS->GetObject(CalibrationId);
		ARM_SABRToHestonSmileCalibration * object = NULL;
		if( !(object = dynamic_cast< ARM_SABRToHestonSmileCalibration* >( armObj )) )
		{
			result.setMsg ("ARM_ERR: object is not of a good type");
			return ARM_KO;
		};
		
		if(info == 0)
			result.setDouble(object->GetV0());
		else if(info == 1)
			result.setDouble(object->GetKappa());
		else if(info == 2)
			result.setDouble(object->GetTheta());
		else if(info == 3)
			result.setDouble(object->GetVVol());
		else if(info == 4)
			result.setDouble(object->GetRho());
		else if(info == 5)
			result.setDouble(object->GetShift());
		else if(info == 6)
			result.setDouble(object->GetLevel());
		else if(info == 7)
			result.setDouble(object->GetCalibError());

		return ARM_OK;
	}

	catch(Exception& x)
	{
		x.DebugPrint();
		ARM_RESULT();
	}
}

long ARM_Heston_CalibrateFunctor::operator()( ARM_result& result, long objId )
{
	/// input checks
	if ( !GlobalPersistanceOk(result) )
	   return(ARM_KO);

	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	ARM_GenericParams* genericParams = GetGenericParams();
	try
	{
		double asOfDate = genericParams->GetParamValue("AsOfDate").GetDouble();

		ResizeRetValues(2,2);
		SetValue(0,0,1);
		SetValue(0,1,"OK");
		SetValue(1,0,"OK");
		SetValue(1,1,1);

		return ARM_OK; 
	}

	catch(Exception& x)
	{
		x.DebugPrint();
		ARM_RESULT();
	}
}

long ARM_SABR_SmileCalibration_ParamFunctor::operator () (ARM_result& result, long objId)
{
	if ( !GlobalPersistanceOk(result) )
	   return(ARM_KO);

	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	ARM_GenericParams* genericParams = GetGenericParams();
	try
	{
		double alpha = genericParams->GetParamValue("alpha").GetDouble();
		double beta = genericParams->GetParamValue("beta").GetDouble();
		double rho = genericParams->GetParamValue("rho").GetDouble();
		double nu = genericParams->GetParamValue("nu").GetDouble();
		double shift = genericParams->GetParamValue("shift").GetDouble();
		int flag = ARM_ArgConv_SABR_ImplicitVol_Formula_Extended_Flag.GetNumber(genericParams->GetParamValue("sabr flag").GetString());

		bool calibalpha = fabs(alpha + 999.) < K_DOUBLE_TOL ? true : false;
		bool calibbeta = fabs(beta + 999.) < K_DOUBLE_TOL ? true : false;
		bool calibrho = fabs(rho + 999.) < K_DOUBLE_TOL ? true : false;
		bool calibnu  = fabs(nu + 999.) < K_DOUBLE_TOL ? true : false;
		bool calibshift = fabs(shift) < K_DOUBLE_TOL ? true : false;
		
		string yesNo = genericParams->GetParamValue("localRhoCalib").GetString();
		stringToUpper(yesNo);
		bool calibeachrho = yesNo == "Y" ? true : false;

		ARM_SmileCalibration_Params_SABR * obj = new ARM_SmileCalibration_Params_SABR(alpha, beta, rho, nu, flag, shift,
			calibalpha, calibbeta, calibrho, calibnu, calibshift, calibeachrho);

		/// assign object
		if( !assignObject( obj, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }

		return ARM_OK; 
	}

	catch(Exception& x)
	{
		x.DebugPrint();
		ARM_RESULT();
	}
}

long ARM_SABR2B_SmileCalibration_ParamFunctor::operator () (ARM_result& result, long objId)
{
	if ( !GlobalPersistanceOk(result) )
	   return(ARM_KO);

	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	ARM_GenericParams* genericParams = GetGenericParams();
	try
	{
		double alpha = genericParams->GetParamValue("alpha").GetDouble();
		double beta1 = genericParams->GetParamValue("beta1").GetDouble();
		double beta2 = genericParams->GetParamValue("beta2").GetDouble();
		double rho = genericParams->GetParamValue("rho").GetDouble();
		double nu = genericParams->GetParamValue("nu").GetDouble();
		double zero = genericParams->GetParamValue("zero").GetDouble();
		double lambda = genericParams->GetParamValue("lambda").GetDouble();
		
		bool calibalpha = fabs(alpha + 999.) < K_DOUBLE_TOL ? true : false;
		bool calibbeta1 = fabs(beta1 + 999.) < K_DOUBLE_TOL ? true : false;
		bool calibbeta2 = fabs(beta2 + 999.) < K_DOUBLE_TOL ? true : false;
		bool calibrho = fabs(rho + 999.) < K_DOUBLE_TOL ? true : false;
		bool calibnu  = fabs(nu + 999.) < K_DOUBLE_TOL ? true : false;
		bool calibzero = fabs(zero + 999.) < K_DOUBLE_TOL ? true : false;
		bool caliblambda = fabs(lambda + 999.) < K_DOUBLE_TOL ? true : false;

		ARM_SmileCalibration_Params_SABR2beta * obj = new ARM_SmileCalibration_Params_SABR2beta(
			alpha, beta1, beta2, rho, nu, zero, lambda,
			calibalpha, calibbeta1, calibbeta2, calibrho, calibnu, calibzero, caliblambda);

		/// assign object
		if( !assignObject( obj, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }

		return ARM_OK; 
	}

	catch(Exception& x)
	{
		x.DebugPrint();
		ARM_RESULT();
	}
}

long ARM_Heston_SmileCalibration_ParamFunctor::operator () (ARM_result& result, long objId)
{
	if ( !GlobalPersistanceOk(result) )
	   return(ARM_KO);

	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	ARM_GenericParams* genericParams = GetGenericParams();
	try
	{
		double v0 = genericParams->GetParamValue("v0").GetDouble();
		double kappa = genericParams->GetParamValue("kappa").GetDouble();
		double theta = genericParams->GetParamValue("theta").GetDouble();
		double rho = genericParams->GetParamValue("rho").GetDouble();
		double vvol = genericParams->GetParamValue("vvol").GetDouble();
		double shift = genericParams->GetParamValue("shift").GetDouble();
		double level = genericParams->GetParamValue("level").GetDouble();
		long sigmaId = genericParams->GetParamValue("sigma").GetObjectId();
		
		std::vector<double> * sigma = NULL;
		if(sigmaId > -1)
		{
			if( !GetObjectFromId( &sigma, sigmaId, ARM_GP_VECTOR) )
			{
				result.setMsg ("ARM_ERR: sigma is not of a good type");
				return ARM_KO;
			}			
		}
		else
		{
			sigma = new std::vector<double>(1,0.);
		}

		bool calibkappa = fabs(kappa + 999.) < K_DOUBLE_TOL ? true : false;
		bool calibtheta = fabs(theta + 999.) < K_DOUBLE_TOL ? true : false;
		bool calibrho = fabs(rho + 999.) < K_DOUBLE_TOL ? true : false;
		bool calibvvol = fabs(vvol + 999.) < K_DOUBLE_TOL ? true : false;

		bool calibshift = fabs(shift) < K_DOUBLE_TOL ? true : false;

		string yesNo = genericParams->GetParamValue("localRhoCalib").GetString();
		stringToUpper(yesNo);
		bool calibeachrho = yesNo == "Y" ? true : false;

		yesNo = genericParams->GetParamValue("normalheston").GetString();
		stringToUpper(yesNo);
		bool normheston = yesNo == "Y" ? true : false;
		
		yesNo = genericParams->GetParamValue("bootstraplevel").GetString();
		stringToUpper(yesNo);
		bool bootstrap = yesNo == "Y" ? true : false;

		ARM_SmileCalibration_Params_Heston * obj = new ARM_SmileCalibration_Params_Heston(v0, kappa, theta, rho, vvol, shift, 
			level, *sigma, calibkappa, calibtheta, calibrho, calibvvol, calibshift, calibeachrho, bootstrap, normheston);

		if(sigmaId == -1)
		{
			delete sigma;
		}

		/// assign object
		if( !assignObject( obj, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }

		return ARM_OK; 
	}

	catch(Exception& x)
	{
		x.DebugPrint();
		ARM_RESULT();
	}
}

long ARM_Heston2b_SmileCalibration_ParamFunctor::operator () (ARM_result& result, long objId)
{
	if ( !GlobalPersistanceOk(result) )
	   return(ARM_KO);

	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	ARM_GenericParams* genericParams = GetGenericParams();
	try
	{
		double v01 = genericParams->GetParamValue("v01").GetDouble();
		double kappa1 = genericParams->GetParamValue("kappa1").GetDouble();
		double theta1 = genericParams->GetParamValue("theta1").GetDouble();
		double rho1 = genericParams->GetParamValue("rho1").GetDouble();
		double vvol1 = genericParams->GetParamValue("vvol1").GetDouble();
		double v02 = genericParams->GetParamValue("v02").GetDouble();
		double kappa2 = genericParams->GetParamValue("kappa2").GetDouble();
		double theta2 = genericParams->GetParamValue("theta2").GetDouble();
		double rho2 = genericParams->GetParamValue("rho2").GetDouble();
		double vvol2 = genericParams->GetParamValue("vvol2").GetDouble();
		double shift = genericParams->GetParamValue("shift").GetDouble();
		double level = genericParams->GetParamValue("level").GetDouble();

		bool calibkappa1 = fabs(kappa1 + 999.) < K_DOUBLE_TOL ? true : false;
		bool calibtheta1 = fabs(theta1 + 999.) < K_DOUBLE_TOL ? true : false;
		bool calibrho1 = fabs(rho1 + 999.) < K_DOUBLE_TOL ? true : false;
		bool calibvvol1 = fabs(vvol1 + 999.) < K_DOUBLE_TOL ? true : false;

		bool calibkappa2 = fabs(kappa2 + 999.) < K_DOUBLE_TOL ? true : false;
		bool calibtheta2 = fabs(theta2 + 999.) < K_DOUBLE_TOL ? true : false;
		bool calibrho2 = fabs(rho2 + 999.) < K_DOUBLE_TOL ? true : false;
		bool calibvvol2 = fabs(vvol2 + 999.) < K_DOUBLE_TOL ? true : false;

		bool calibshift = fabs(shift) < K_DOUBLE_TOL ? true : false;

		string yesNo = genericParams->GetParamValue("localRhoCalib1").GetString();
		stringToUpper(yesNo);
		bool calibeachrho1 = yesNo == "Y" ? true : false;
	
		yesNo = genericParams->GetParamValue("localRhoCalib2").GetString();
		stringToUpper(yesNo);
		bool calibeachrho2 = yesNo == "Y" ? true : false;

		yesNo = genericParams->GetParamValue("bootstraplevel").GetString();
		stringToUpper(yesNo);
		bool bootstrap = yesNo == "Y" ? true : false;

		ARM_SmileCalibration_Params_Heston2b * obj = new ARM_SmileCalibration_Params_Heston2b(v01, kappa1, theta1, rho1, vvol1,
			v02, kappa2, theta2, rho2, vvol2, shift, level, 
			calibkappa1, calibtheta1, calibrho1, calibvvol1,
			calibkappa2, calibtheta2, calibrho2, calibvvol2,
			calibshift, calibeachrho1, calibeachrho2, bootstrap);

		/// assign object
		if( !assignObject( obj, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }

		return ARM_OK; 
	}

	catch(Exception& x)
	{
		x.DebugPrint();
		ARM_RESULT();
	}
}

long ARM_BiSABR_SmileCalibration_ParamFunctor::operator () (ARM_result& result, long objId)
{
	if ( !GlobalPersistanceOk(result) )
	   return(ARM_KO);

	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	ARM_GenericParams* genericParams = GetGenericParams();
	try
	{
		double fwd1 = genericParams->GetParamValue("forward1").GetDouble();
		double alpha1 = genericParams->GetParamValue("alpha1").GetDouble();
		double rho1 = genericParams->GetParamValue("rho1").GetDouble();
		double nu1 = genericParams->GetParamValue("nu1").GetDouble();
		double beta1 = genericParams->GetParamValue("beta1").GetDouble();

		double fwd2 = genericParams->GetParamValue("forward2").GetDouble();
		double alpha2 = genericParams->GetParamValue("alpha2").GetDouble();
		double rho2 = genericParams->GetParamValue("rho2").GetDouble();
		double nu2 = genericParams->GetParamValue("nu2").GetDouble();
		double beta2 = genericParams->GetParamValue("beta2").GetDouble();

		double rhos1s2 = genericParams->GetParamValue("rhos1s2").GetDouble();
		double rhos1v2 = genericParams->GetParamValue("rhos1v2").GetDouble();
		double rhos2v1 = genericParams->GetParamValue("rhos2v1").GetDouble();
		double rhov1v2 = genericParams->GetParamValue("rhov1v2").GetDouble();

		bool calibrhoss = fabs(rhos1s2 + 999.) < K_DOUBLE_TOL ? true : false;
		bool calibrhosv = fabs(rhos1v2 + 999.) < K_DOUBLE_TOL || fabs(rhos2v1 + 999.) < K_DOUBLE_TOL ? true : false;
		bool calibrhovv = fabs(rhov1v2 + 999.) < K_DOUBLE_TOL ? true : false;

		ARM_SmileCalibration_Params_BiSABR * obj = new ARM_SmileCalibration_Params_BiSABR(fwd1, alpha1, beta1, rho1, nu1,
			fwd2, alpha2, beta2, rho2, nu2, rhos1s2, rhos1v2, rhos2v1, rhov1v2, 
			calibrhoss, calibrhosv, calibrhovv);

		/// assign object
		if( !assignObject( obj, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }

		return ARM_OK; 
	}

	catch(Exception& x)
	{
		x.DebugPrint();
		ARM_RESULT();
	}
}

long ARM_Merton_SmileCalibration_ParamFunctor::operator () (ARM_result& result, long objId)
{
	if ( !GlobalPersistanceOk(result) )
	   return(ARM_KO);

	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	ARM_GenericParams* genericParams = GetGenericParams();
	try
	{
		double sigma = genericParams->GetParamValue("sigma").GetDouble();
		double lambda1 = genericParams->GetParamValue("lambda1").GetDouble();
		double U1 = genericParams->GetParamValue("U1").GetDouble();
		double lambda2 = genericParams->GetParamValue("lambda2").GetDouble();
		double U2 = genericParams->GetParamValue("U2").GetDouble();

		bool calibsigma		= sigma < K_DOUBLE_TOL ? true : false;
		bool caliblambda1	= lambda1 < K_DOUBLE_TOL ? true : false;
		bool calibU1		= U1 < K_DOUBLE_TOL ? true : false;
		bool caliblambda2	= lambda2 < K_DOUBLE_TOL ? true : false;
		bool calibU2		= U2 < K_DOUBLE_TOL ? true : false;

		ARM_SmileCalibration_Params_Merton * obj = new ARM_SmileCalibration_Params_Merton(sigma, lambda1, U1, 
			lambda2, U2, calibsigma, caliblambda1, calibU1, caliblambda2, calibU2);

		/// assign object
		if( !assignObject( obj, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }

		return ARM_OK; 
	}

	catch(Exception& x)
	{
		x.DebugPrint();
		ARM_RESULT();
	}
}


extern long ARMLOCAL_SmileCalibration(
	double C_AsOfDate,
	const vector<double>& C_CalibTimes,
	const vector<double>& C_Forwards,
	const vector< vector<double> >& C_MktVols,
	const vector< vector<double> >& C_Strikes,
	const long& C_CalibParamId,
	const vector<double>& C_ConstraintStrikes,
	const vector<double>& C_ConstraintVols,
	const bool& C_CalibConstraint,
	const vector<double>& C_Weights,
	int& rowsResult,
	int& colsResult,
	vector<double>& outResult,
	//FIXMEFRED : no vector<bool> anymore
	std::deque<bool>& outResultBool,
	ARM_result& result )
{
	CCString msg ("");

	try
	{
		std::vector<double> calibtimes(C_CalibTimes);
		calibtimes -= C_AsOfDate;
		calibtimes /= 365.;

		std::vector<double> forwards(C_Forwards);
		std::vector<double> weights(C_Weights);
		std::vector<double> constraintstrikes(C_ConstraintStrikes);
		std::vector<double> constraintvols(C_ConstraintVols);

		vector< std::vector<double>& > mktvols, strikes;

		if(C_MktVols.size() > 0 && C_MktVols[0].size() == 1)
		{
			mktvols.resize(1);
			mktvols[0] = new std::vector<double>;
			mktvols[0]->resize(C_MktVols.size());
			for(int i = 0; i < C_MktVols.size(); i++) (*(mktvols[0]))[i] = C_MktVols[i][0];
		}
		else
		{
			mktvols.resize(C_MktVols.size());
			for(int k = 0; k < mktvols.size(); k++)
			{
				mktvols[k] = new std::vector<double>;
				mktvols[k]->resize(C_MktVols[k].size());
				for(int i = 0; i < mktvols[k]->size(); i++) (*(mktvols[k]))[i] = C_MktVols[k][i];
			}
		}
		
		if(C_Strikes.size() > 0 && C_Strikes[0].size() == 1)
		{
			strikes.resize(1);
			strikes[0] = new std::vector<double>;
			strikes[0]->resize(C_Strikes.size());
			for(int i = 0; i < C_Strikes.size(); i++) (*(strikes[0]))[i] = C_Strikes[i][0];
		}
		else
		{
			strikes.resize(C_Strikes.size());
			for(int k = 0; k < strikes.size(); k++)
			{
				strikes[k] = new std::vector<double>;
				strikes[k]->resize(C_Strikes[k].size());
				for(int i = 0; i < strikes[k]->size(); i++) (*(strikes[k]))[i] = C_Strikes[k][i];
			}
		}

		ARM_Object* armObj =  LOCAL_PERSISTENT_OBJECTS->GetObject(C_CalibParamId);

		ARM_SmileCalibration_Params * params;
		ARM_SmileCalibration_Params_SABR * sabrparamObj = dynamic_cast<ARM_SmileCalibration_Params_SABR *>(armObj);
		ARM_SmileCalibration_Params_SABR2beta * sabr2bparamObj = dynamic_cast<ARM_SmileCalibration_Params_SABR2beta *>(armObj);
		ARM_SmileCalibration_Params_Heston * hestonparamObj = dynamic_cast<ARM_SmileCalibration_Params_Heston *>(armObj);;
		ARM_SmileCalibration_Params_Heston2b * heston2bparamObj = dynamic_cast<ARM_SmileCalibration_Params_Heston2b *>(armObj);;
		ARM_SmileCalibration_Params_BiSABR * bisabrparamObj = dynamic_cast<ARM_SmileCalibration_Params_BiSABR *>(armObj);;
		ARM_SmileCalibration_Params_Merton * mertonObj = dynamic_cast<ARM_SmileCalibration_Params_Merton *>(armObj);

		ARM_SmileCalibration * calibtool = NULL;

		int type = 0;

		if((params = sabrparamObj) != NULL)
		{
			calibtool = new ARM_SmileCalibration_SABR;
			type = 1;
		}
		else if((params = hestonparamObj) != NULL)
		{
			if(hestonparamObj->normalHeston() == false)
				calibtool = new ARM_SmileCalibration_Heston;
			else
				calibtool = new ARM_SmileCalibration_NormalHeston;

			type = 2;
		}
		else if((params = bisabrparamObj) != NULL)
		{
			calibtool = new ARM_SmileCalibration_BiSABR;
			type = 3;
		}
		else if((params = sabr2bparamObj) != NULL)
		{
			calibtool = new ARM_SmileCalibration_SABR2beta;
			type = 4;
		}
		else if((params = mertonObj) != NULL)
		{
			calibtool = new ARM_SmileCalibration_Merton;
			type = 5;
		}
		else if((params = heston2bparamObj) != NULL)
		{
			calibtool = new ARM_SmileCalibration_Heston2b;
			type = 6;
		}
		else
		{
			result.setMsg ("ARM_ERR: object is not of a good type");
			return ARM_KO;
		};
		
		if(calibtimes.size() == 1)
		{
			calibtool->Init(calibtimes[0], forwards[0], *(mktvols[0]), *(strikes[0]), 
				C_CalibConstraint, constraintstrikes[0], constraintvols[0], params);
		}
		else
		{
			calibtool->Init(calibtimes, forwards, mktvols, strikes, 
				C_CalibConstraint, constraintstrikes, constraintvols, weights, params);
		}

		calibtool->Calibrate();

		if(type == 1)
		{
			int cols = 4 + (sabrparamObj->shifts().size() > 0) + (sabrparamObj->levels().size() > 0);
			int rows = 1;
			if(sabrparamObj->rhos().size() > rows) rows = sabrparamObj->rhos().size();
			if(sabrparamObj->shifts().size() > rows) rows = sabrparamObj->shifts().size();
			if(sabrparamObj->levels().size() > rows) rows = sabrparamObj->levels().size();

			outResult.resize(rows*cols);
			outResultBool.resize(rows*cols,false);
			
			outResult[0 * rows] = sabrparamObj->alpha();
			outResult[1 * rows] = sabrparamObj->beta();
			outResult[3 * rows] = sabrparamObj->nu();

			outResultBool[0 * rows] = true;
			outResultBool[1 * rows] = true;
			outResultBool[3 * rows] = true;
            int k;
			int col = 2;
			for(k = 0; k < sabrparamObj->rhos().size(); k++) 
			{
				outResult[k + col * rows] = sabrparamObj->rhos()[k];
				outResultBool[k + col * rows] = true;
			}

			col = 4;
			for(k = 0; k < sabrparamObj->shifts().size(); k++) 
			{
				outResult[k + col * rows] = sabrparamObj->shifts()[k];
				outResultBool[k + col * rows] = true;
			}
			if(sabrparamObj->shifts().size() > 0) col ++;

			for(k = 0; k < sabrparamObj->levels().size(); k++) 
			{
				outResult[k + col * rows] = sabrparamObj->levels()[k];
				outResultBool[k + col * rows] = true;
			}

			rowsResult = rows;
			colsResult = cols;
		}
		else if(type == 2)
		{
			int cols = 5 + (hestonparamObj->shifts().size() > 0) + (hestonparamObj->levels().size() > 0);
			int rows = 1;
			if(hestonparamObj->rhos().size() > rows) rows = hestonparamObj->rhos().size();
			if(hestonparamObj->shifts().size() > rows) rows = hestonparamObj->shifts().size();
			if(hestonparamObj->levels().size() > rows) rows = hestonparamObj->levels().size();

			outResult.resize(rows*cols);
			outResultBool.resize(rows*cols,false);
			
			outResult[0 * rows] = hestonparamObj->v0();
			outResult[1 * rows] = hestonparamObj->kappa();
			outResult[2 * rows] = hestonparamObj->theta();
			outResult[4 * rows] = hestonparamObj->nu();

			outResultBool[0 * rows] = true;
			outResultBool[1 * rows] = true;
			outResultBool[2 * rows] = true;
			outResultBool[4 * rows] = true;

			int col = 3;
			int k;
			for(k = 0; k < hestonparamObj->rhos().size(); k++) 
			{
				outResult[k + col * rows] = hestonparamObj->rhos()[k];
				outResultBool[k + col * rows] = true;
			}
			col = 5;
			for(k = 0; k < hestonparamObj->shifts().size(); k++) 
			{
				outResult[k + col * rows] = hestonparamObj->shifts()[k];
				outResultBool[k + col * rows] = true;
			}
			if(hestonparamObj->shifts().size() > 0) col ++;

			for(k = 0; k < hestonparamObj->levels().size(); k++) 
			{
				outResult[k + col * rows] = hestonparamObj->levels()[k];
				outResultBool[k + col * rows] = true;
			}

			rowsResult = rows;
			colsResult = cols;
		}
		else if(type == 3)
		{
			outResult.resize(4);
			outResultBool.resize(4,true);
			outResult[0] = bisabrparamObj->rhoS1S2();
			outResult[1] = bisabrparamObj->rhoS1V2();
			outResult[2] = bisabrparamObj->rhoS2V1();
			outResult[3] = bisabrparamObj->rhoV1V2();

			rowsResult = 1;
			colsResult = 4;
		}
		else if(type == 4)
		{
			outResult.resize(7);
			outResultBool.resize(7,true);
			
			outResult[0] = sabr2bparamObj->alpha();
			outResult[1] = sabr2bparamObj->beta1();
			outResult[2] = sabr2bparamObj->beta2();
			outResult[3] = sabr2bparamObj->rho();
			outResult[4] = sabr2bparamObj->nu();
			outResult[5] = sabr2bparamObj->zero();
			outResult[6] = sabr2bparamObj->lambda();

			rowsResult = 1;
			colsResult = 7;
		}
		else if(type == 5)
		{
			outResult.resize(5);
			outResultBool.resize(5,true);

			outResult[0] = mertonObj->sigma();
			outResult[1] = mertonObj->lambda1();
			outResult[2] = mertonObj->U1();
			outResult[3] = mertonObj->lambda2();
			outResult[4] = mertonObj->U2();

			rowsResult = 1;
			colsResult = 5;
		}
		else if(type == 6)
		{
			int cols = 10 + (heston2bparamObj->shifts().size() > 0) + (heston2bparamObj->levels().size() > 0);
			int rows = 1;
			if(heston2bparamObj->rhos1().size() > rows) rows = heston2bparamObj->rhos1().size();
			if(heston2bparamObj->rhos2().size() > rows) rows = heston2bparamObj->rhos2().size();
			if(heston2bparamObj->shifts().size() > rows) rows = heston2bparamObj->shifts().size();
			if(heston2bparamObj->levels().size() > rows) rows = heston2bparamObj->levels().size();

			outResult.resize(rows*cols);
			outResultBool.resize(rows*cols,false);
			
			outResult[0 * rows] = heston2bparamObj->v01();
			outResult[1 * rows] = heston2bparamObj->kappa1();
			outResult[2 * rows] = heston2bparamObj->theta1();
			outResult[4 * rows] = heston2bparamObj->nu1();
			outResult[5 * rows] = heston2bparamObj->v02();
			outResult[6 * rows] = heston2bparamObj->kappa2();
			outResult[7 * rows] = heston2bparamObj->theta2();
			outResult[9 * rows] = heston2bparamObj->nu2();

			outResultBool[0 * rows] = true;
			outResultBool[1 * rows] = true;
			outResultBool[2 * rows] = true;
			outResultBool[4 * rows] = true;
			outResultBool[5 * rows] = true;
			outResultBool[6 * rows] = true;
			outResultBool[7 * rows] = true;
			outResultBool[9 * rows] = true;

			int col = 3;
			for(int k = 0; k < heston2bparamObj->rhos1().size(); k++) 
			{
				outResult[k + col * rows] = heston2bparamObj->rhos1()[k];
				outResultBool[k + col * rows] = true;
			}
			col = 8;
			for(k = 0; k < heston2bparamObj->rhos2().size(); k++) 
			{
				outResult[k + col * rows] = heston2bparamObj->rhos2()[k];
				outResultBool[k + col * rows] = true;
			}
			col = 10;
			for(k = 0; k < heston2bparamObj->shifts().size(); k++) 
			{
				outResult[k + col * rows] = heston2bparamObj->shifts()[k];
				outResultBool[k + col * rows] = true;
			}
			if(heston2bparamObj->shifts().size() > 0) col ++;

			for(k = 0; k < heston2bparamObj->levels().size(); k++) 
			{
				outResult[k + col * rows] = heston2bparamObj->levels()[k];
				outResultBool[k + col * rows] = true;
			}

			rowsResult = rows;
			colsResult = cols;
		}

		delete calibtool;
		
		return ARM_OK;
	}
	catch(Exception& x)
	{
		x.DebugPrint();
		
		ARM_RESULT();
	}
}

extern long ARMLOCAL_SpreadSmileCalibration2Heston(
	double C_ResetTime,
	double C_Fwd1,
	const vector<double>& C_MktVols1,
	const vector<double>& C_Strikes1,
	double C_ConstrVol1,
	double C_ConstrK1,
	double C_Fwd2,
	const vector<double>& C_MktVols2,
	const vector<double>& C_Strikes2,
	double C_ConstrVol2,
	double C_ConstrK2,
	const vector<double>& C_MktVolsSpread,
	const vector<double>& C_StrikesSpread,
	double C_ConstrVolSpread,
	double C_ConstrKSpread,
	const double& C_v0,
	const double& C_kappa,
	const double& C_theta,
	const bool& C_calibTheta,
	int& rowsResult,
	int& colsResult,
	vector<double>& outResult,
	//FIXMEFRED : no vector<bool> anymore
	std::deque<bool>& outResultBool,
	ARM_result& result)
{
	CCString msg ("");

	try
	{
		std::vector<double> mktvols1(C_MktVols1);
		std::vector<double> strikes1(C_Strikes1);
		std::vector<double> mktvols2(C_MktVols2);
		std::vector<double> strikes2(C_Strikes2);
		std::vector<double> mktvolspread(C_MktVolsSpread);
		std::vector<double> strikespread(C_StrikesSpread);

		ARM_SmileCalibration_Params_Spread2Heston params(C_v0, C_kappa, C_theta, C_calibTheta);

		ARM_SmileCalibration_Spread2Heston calibtool;
		
		calibtool.Init(C_ResetTime, C_Fwd1, mktvols1, strikes1, C_ConstrVol1, C_ConstrK1, 
			C_Fwd2, mktvols2, strikes2, C_ConstrVol2, C_ConstrK2, 
			mktvolspread, strikespread, C_ConstrVolSpread, C_ConstrKSpread, &params);

		calibtool.Calibrate();

		outResult.resize(11);
		outResultBool.resize(11,true);

		outResult[0] = params.v0();
		outResult[1] = params.kappa();
		outResult[2] = params.theta();
		outResult[3] = params.nu();
		outResult[4] = params.rho1();
		outResult[5] = params.shift1();
		outResult[6] = params.level1();
		outResult[7] = params.rho2();
		outResult[8] = params.shift2();
		outResult[9] = params.level2();
		outResult[10] = params.correl();

		rowsResult = 1;
		colsResult = 11;

		return ARM_OK;
	}
	catch(Exception& x)
	{
		x.DebugPrint();
		
		ARM_RESULT();
	}
}

extern long ARMLOCAL_Spread2Heston_TOTEMCalibration(
	const vector<double>& TOTEMMat,
	const vector<double>& TOTEMStrikes,
	const vector<double>& TOTEMPrices,
	const vector<double>& FullScheduleReset,
	const vector<double>& FullScheduleAnnuity,
	const vector<double>& FullScheduleFwd1,
	const vector<double>& FullScheduleFwd2,
	const vector<double>& FwdCalibReset,
	const vector<double>& LongFwds,
	const long& LongVolsId,
	const long& LongStrikesId,
	const vector<double>& ShortFwds,
	const long& ShortVolsId,
	const long& ShortStrikesId,
	const double& v0,
	const double& kappa,
	const double& theta,
	const bool& constrCorrel,
	int& nbRowsRes,
	int& nbColsRes,
	vector<double>& Res,
	//FIXMEFRED : no vector<bool> anymore
	std::deque<bool>& outResBool,
	ARM_result& result
	)
{
	CCString msg ("");

	try
	{
		std::vector<double> TotemMat(TOTEMMat), TotemStrikes(TOTEMStrikes), TotemPrices(TOTEMPrices);
		std::vector<double> FullSchedReset(FullScheduleReset), FullSchedAnnuity(FullScheduleAnnuity), FullSchedFwd1(FullScheduleFwd1), FullSchedFwd2(FullScheduleFwd2);
		std::vector<double> longFwds(LongFwds), shortFwds(ShortFwds);
		std::vector<double> fwdCalibReset(FwdCalibReset);
		ARM_GP_Matrix * longVols, * longK, * shortVols, * shortK;
		std::vector<double> Mat, Theta, Nu, LongRho, LongShift, LongLevel, ShortRho, ShortShift, ShortLevel, Correl;
		
		ARM_Object * obj;

		obj = LOCAL_PERSISTENT_OBJECTS->GetObject(LongVolsId);
		if((longVols = dynamic_cast<ARM_GP_Matrix*>(obj)) == NULL)
		{
			result.setMsg("ARM_ERR : Long CMS Vol is not of good type, expected gp matrix object");
			return ARM_KO;
		}

		obj = LOCAL_PERSISTENT_OBJECTS->GetObject(LongStrikesId);
		if((longK = dynamic_cast<ARM_GP_Matrix*>(obj)) == NULL)
		{
			result.setMsg("ARM_ERR : Long CMS Strikes is not of good type, expected gp matrix object");
			return ARM_KO;
		}

		obj = LOCAL_PERSISTENT_OBJECTS->GetObject(ShortVolsId);
		if((shortVols = dynamic_cast<ARM_GP_Matrix*>(obj)) == NULL)
		{
			result.setMsg("ARM_ERR : Short CMS Vol is not of good type, expected gp matrix object");
			return ARM_KO;
		}

		obj = LOCAL_PERSISTENT_OBJECTS->GetObject(ShortStrikesId);
		if((shortK = dynamic_cast<ARM_GP_Matrix*>(obj)) == NULL)
		{
			result.setMsg("ARM_ERR : Short CMS Strikes is not of good type, expected gp matrix object");
			return ARM_KO;
		}

		ARM_Spread2Heston_TOTEMCalibration(TotemMat, TotemStrikes, TotemPrices, 
				FullSchedReset, FullSchedAnnuity, FullSchedFwd1, FullSchedFwd2,
				fwdCalibReset,
				longFwds, *longVols, *longK, shortFwds, *shortVols, *shortK,
				v0, kappa, theta, constrCorrel, 
				Mat, Theta, Nu, LongRho, LongShift, LongLevel, ShortRho, ShortShift, ShortLevel, Correl);


		int k, corrsize = Mat.size(), fwdsize = fwdCalibReset.size();
		

		nbRowsRes = corrsize > fwdsize ? corrsize : fwdsize;
		nbColsRes = 11;

		Res.resize(nbRowsRes * nbColsRes);
		outResBool.resize(nbRowsRes * nbColsRes, true);

		int col = 0;
		for(k = 0; k < corrsize; k++) Res[k + col * nbRowsRes]  = Mat[k];
		for(k = corrsize; k < nbRowsRes; k++) outResBool[k + col * nbRowsRes] = false;
		col ++;
		for(k = 0; k < corrsize; k++) Res[k + col * nbRowsRes]  = Correl[k];
		for(k = corrsize; k < nbRowsRes; k++) outResBool[k + col * nbRowsRes] = false;

		col ++;
		for(k = 0; k < fwdsize; k++) Res[k + col * nbRowsRes]  = fwdCalibReset[k];
		for(k = fwdsize; k < nbRowsRes; k++) outResBool[k + col * nbRowsRes] = false;
		
		col ++;
		for(k = 0; k < fwdsize; k++) Res[k + col * nbRowsRes]  = Theta[k];
		for(k = fwdsize; k < nbRowsRes; k++) outResBool[k + col * nbRowsRes] = false;
		
		col ++;
		for(k = 0; k < fwdsize; k++) Res[k + col * nbRowsRes]  = Nu[k];
		for(k = fwdsize; k < nbRowsRes; k++) outResBool[k + col * nbRowsRes] = false;
		
		col ++;
		for(k = 0; k < fwdsize; k++) Res[k + col * nbRowsRes]  = LongRho[k];
		for(k = fwdsize; k < nbRowsRes; k++) outResBool[k + col * nbRowsRes] = false;
		
		col ++;
		for(k = 0; k < fwdsize; k++) Res[k + col * nbRowsRes]  = LongShift[k];
		for(k = fwdsize; k < nbRowsRes; k++) outResBool[k + col * nbRowsRes] = false;
		
		col ++;
		for(k = 0; k < fwdsize; k++) Res[k + col * nbRowsRes]  = LongLevel[k];
		for(k = fwdsize; k < nbRowsRes; k++) outResBool[k + col * nbRowsRes] = false;
		
		col ++;
		for(k = 0; k < fwdsize; k++) Res[k + col * nbRowsRes]  = ShortRho[k];
		for(k = fwdsize; k < nbRowsRes; k++) outResBool[k + col * nbRowsRes] = false;
		
		col ++;
		for(k = 0; k < fwdsize; k++) Res[k + col * nbRowsRes]  = ShortShift[k];
		for(k = fwdsize; k < nbRowsRes; k++) outResBool[k + col * nbRowsRes] = false;
		
		col ++;
		for(k = 0; k < fwdsize; k++) Res[k + col * nbRowsRes]  = ShortLevel[k];
		for(k = fwdsize; k < nbRowsRes; k++) outResBool[k + col * nbRowsRes] = false;
		
		return ARM_OK;
		
	}
	catch(Exception& x)
	{
		x.DebugPrint();
		
		ARM_RESULT();
	}
}

extern long ARMLOCAL_Spread2HestonVanilla(
	double C_ResetTime,
	double C_Fwd1,
	double C_Fwd2,
	double C_Strike,
	double C_CallPut,
	double C_V0,
	double C_Kappa,
	double C_Theta,
	double C_Nu,
	double C_Rho1,
	double C_Rho2,
	double C_Shift1,
	double C_Shift2,
	double C_Level1,
	double C_Level2,
	double C_Correl,
	double C_Index1Lev,
	double C_Index2Lev,
	ARM_result& result)
{
	CCString msg ("");
	
	try
	{
		double price = Spread2HestonVanilla(C_ResetTime, C_Fwd1, C_Fwd2, C_Strike, (int)C_CallPut, 
							C_V0, C_Kappa, C_Theta, C_Nu, C_Rho1, C_Rho2, C_Shift1, C_Shift2, 
							C_Level1, C_Level2, C_Correl, C_Index1Lev, C_Index2Lev);

		result.setDouble(price);

		return ARM_OK;
	}	
	catch(Exception& x)
	{
		x.DebugPrint();
		
		ARM_RESULT();
	}

}

extern long ARMLOCAL_TRiSABR_VanillaOption(
											const vector<double>& C_S1params_vec,
											const vector<double>& C_S2params_vec,
											const vector<double>& C_S3params_vec,
											double	C_rhos12,
											double	C_rhos23,
											double	C_rhos13,
											double	C_rhov12,
											double	C_rhov23,
											double	C_rhov13,
											double	C_rhoc12,
											double	C_rhoc21,
											double	C_rhoc23,
											double	C_rhoc32,
											double	C_rhoc13,
											double	C_rhoc31,
											double	C_K,
											double	C_T,	
											double	C_callput,
											double	C_flag,
								
											ARM_result& result)
{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	std::vector<double>&		S1params_vec		= NULL;	
	std::vector<double>&		S2params_vec		= NULL;	
	std::vector<double>&		S3params_vec		= NULL;	
	
	try
	{
		S1params_vec			= CreateARMGPVectorFromVECTOR( C_S1params_vec		);
		const vector<double> S1params_Vec=S1params_vec->GetValues();
		double F1=S1params_Vec[0];
		double alpha1=S1params_Vec[1];
		double beta1=S1params_Vec[2];
		double rho1=S1params_Vec[3];
		double nu1=S1params_Vec[4];
		
		S2params_vec			= CreateARMGPVectorFromVECTOR( C_S2params_vec		);
		const vector<double> S2params_Vec=S2params_vec->GetValues();
		double F2=S2params_Vec[0];
		double alpha2=S2params_Vec[1];
		double beta2=S2params_Vec[2];
		double rho2=S2params_Vec[3];
		double nu2=S2params_Vec[4];
		
		S3params_vec			= CreateARMGPVectorFromVECTOR( C_S3params_vec		);
		const vector<double> S3params_Vec=S3params_vec->GetValues();
		double F3=S3params_Vec[0];
		double alpha3=S3params_Vec[1];
		double beta3=S3params_Vec[2];
		double rho3=S3params_Vec[3];
		double nu3=S3params_Vec[4];

		double	C_nbsteps=120;  ///  default option for a numerical integration of the SABR

		double dResult = Export_TriSABR_VanillaOption(
			 F1,
			 alpha1,
			 beta1,
			 rho1,
			 nu1,
								 
			 F2,
			 alpha2,
			 beta2,
			 rho2,
			 nu2,
								 
			 F3,
			 alpha3,
			 beta3,
			 rho3,
			 nu3,
			 
			 C_rhos12,
			 C_rhos23,
			 C_rhos13,
			 C_rhov12,
			 C_rhov23,
			 C_rhov13,
			 C_rhoc12,
			 C_rhoc21,
			 C_rhoc23,
			 C_rhoc32,
			 C_rhoc13,
			 C_rhoc31,
			 C_K,
			 C_T,
			 C_callput,
			 C_flag,
			 C_nbsteps);
		result.setDouble(dResult);
		return ARM_OK;
	}
	
	catch(Exception& x)
	{
		x.DebugPrint();
		
		ARM_RESULT();
	}
}

extern long ARMLOCAL_SABR_To_Heston_SmileCalibration_Create(
	const double&  C_resetTime,
	const double&  C_fwdRate,
	const double&  C_ATMVol,
	const double&  C_Alpha,
	const double&  C_Beta,
	const double&  C_RhoSABR,
	const double&  C_Nu,
	const double&  C_Sabr_Type,
	const double&  C_InitialVar,
	const double&  C_Kappa,
	const double&  C_Theta,
	const double&  C_Rho,
	const double&  C_Shift,
	const bool& C_CalibTheta,
	const bool& C_CalibRho,
	const bool& C_CalibShift,
	ARM_result&				result,
	long					objId)
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");
		

	try
	{
		ARM_SABRToHestonSmileCalibration * obj = new ARM_SABRToHestonSmileCalibration(C_resetTime, C_fwdRate, C_ATMVol,
															C_Alpha, C_Beta, C_RhoSABR, C_Nu, (int)C_Sabr_Type, 
															C_InitialVar, C_Kappa, C_Theta, C_Rho, C_Shift, 
															C_CalibTheta, C_CalibRho, C_CalibShift);

		/// assign object
		if( !assignObject( obj, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }

	}

	catch(Exception& x)
	{
		x.DebugPrint();
		ARM_RESULT();
	}
}



////////////////////////////////////////////
//// Function to create an CIR_ParamModel
////////////////////////////////////////////

extern long ARMLOCAL_CIR_ModelParamsCreate(
	const vector<long>&		modelParamVec,
	ARM_result&				result, 
	long					objId )
{
	if( !GlobalPersistanceOk( result ) )	return ARM_KO;
	CCString msg ("");
	ARM_CIRModelParams*	modParam=NULL;

	try	{
		
		CC_STL_VECTOR(ARM_ModelParam* )	vecModParam;
		for( size_t i=0; i<modelParamVec.size(); ++i ){
			ARM_Object* armObj	=  LOCAL_PERSISTENT_OBJECTS->GetObject(modelParamVec[i] );
			ARM_ModelParam* modelParam = dynamic_cast<ARM_ModelParam*>(armObj);

			if( !modelParam ){
				result.setMsg ("ARM_ERR: model parameter is not of a good type");
				return ARM_KO;
			}
			vecModParam.push_back(modelParam);
		}

		modParam	= new ARM_CIRModelParams( vecModParam );

		if( !assignObject( modParam, result, objId ) )	return ARM_KO; 
		else											return ARM_OK; 
	}
	
	catch(Exception& x){
		delete modParam;
		x.DebugPrint();
		ARM_RESULT();
	}
}

////////////////////////////////////////////
//// Function to price an ARM_CIRBondPrice
////////////////////////////////////////////

extern long ARMLOCAL_CIRBondPrice(
	const long&				modelParamsId,
	const double&			time,
	const double&			r0,
	const double&			phi,
	ARM_result&				result )
{
	if( !GlobalPersistanceOk( result ) )	return ARM_KO;
	CCString msg ("");
	
	try	{

		ARM_Object* obj	=  LOCAL_PERSISTENT_OBJECTS->GetObject(modelParamsId );
		ARM_CIRModelParams* modParam = dynamic_cast<ARM_CIRModelParams*>(obj);
		if( !modParam ){
			result.setMsg ("ARM_ERR: model parameter is not of a good type");
			return ARM_KO;
		}

		double price = modParam->CptExp_Bond(time,r0,phi);

		result.setDouble(price);
	}
	catch(Exception& x){
		x.DebugPrint();
		ARM_RESULT();
	}

	catch (...){
		result.setMsg("ARM_ERR: unrecognized failure in ARMLOCAL_CIRBondPrice");
        return(ARM_KO);
	}
	return ARM_OK;
}

////////////////////////////////////////////
//// Function to price an ARM_CIRBondDensity
////////////////////////////////////////////

extern long ARMLOCAL_CIRBondDensity(
	const long&				modelParamsId,
	const double&			time,
	const double&			r0,
	const double&			r,
	const double&			frequency,
	ARM_result&				result )
{
	if( !GlobalPersistanceOk( result ) )	return ARM_KO;
	CCString msg ("");
	
	try	{

		ARM_Object* obj	=  LOCAL_PERSISTENT_OBJECTS->GetObject(modelParamsId );
		ARM_CIRModelParams* modParam = dynamic_cast<ARM_CIRModelParams*>(obj);
		if( !modParam ){
			result.setMsg ("ARM_ERR: model parameter is not of a good type");
			return ARM_KO;
		}

		double density= modParam->CptDensity_Bond(time,r0,r );

		result.setDouble(density);
	}
	catch(Exception& x){
		x.DebugPrint();
		ARM_RESULT();
	}

	catch (...){
		result.setMsg("ARM_ERR: unrecognized failure in ARMLOCAL_CIRBondDensity");
        return(ARM_KO);
	}
	return ARM_OK;
}



extern long ARMLOCAL_CIRBondDistribution(
	const long&				modelParamsId,
	const double&			time,
	const double&			r0,
	const double&			nbDiscr,
	const double&			factor,
	ARM_result&				result){
	if( !GlobalPersistanceOk( result ) )	return ARM_KO;
	CCString msg ("");
	
	try	{

		ARM_Object* obj	=  LOCAL_PERSISTENT_OBJECTS->GetObject(modelParamsId );
		ARM_CIRModelParams* modParam = dynamic_cast<ARM_CIRModelParams*>(obj);
		if( !modParam ){
			result.setMsg ("ARM_ERR: model parameter is not of a good type");
			return ARM_KO;
		}

		DistribSet dist= modParam->CptDistribSet(time, r0, nbDiscr, factor);

		result.setLong(3*nbDiscr);
		for (int i=0; i<nbDiscr; i++){
			result.setArray(dist.itsPoint[i],3*i);
			result.setArray(dist.itsWeight[i],3*i+1);
			result.setArray(dist.itsDensity[i],3*i+2);
		}

	}
	catch(Exception& x){
		x.DebugPrint();
		ARM_RESULT();
	}

	catch (...){
		result.setMsg("ARM_ERR: unrecognized failure in ARMLOCAL_CIRBondDensity");
        return(ARM_KO);
	}
	return ARM_OK;;
}

////////////////////////////////////////////
//// Function to price an ARM_CIRBondWeight
////////////////////////////////////////////

extern long ARMLOCAL_CIRBondWeight(
	const long&				modelParamsId,
	ARM_result&				result ){

	if( !GlobalPersistanceOk( result ) )	return ARM_KO;
	CCString msg ("");
	
	try	{

		ARM_Object* obj	=  LOCAL_PERSISTENT_OBJECTS->GetObject(modelParamsId );
		ARM_CIRModelParams* modParam = dynamic_cast<ARM_CIRModelParams*>(obj);
		if( !modParam ){
			result.setMsg ("ARM_ERR: model parameter is not of a good type");
			return ARM_KO;
		}

		double res = modParam->CptWeight_Bond();

		result.setDouble(res);
	}
	catch(Exception& x){
		x.DebugPrint();
		ARM_RESULT();
	}

	catch (...){
		result.setMsg("ARM_ERR: unrecognized failure in ARMLOCAL_CIRBondWeight");
        return(ARM_KO);
	}
	return ARM_OK;
}

////////////////////////////////////////////
//// Function to price an ARM_CIRBondExpectation
////////////////////////////////////////////

extern long ARMLOCAL_CIRBondExpectation(
	const long&				modelParamsId,
	ARM_result&				result ){

	if( !GlobalPersistanceOk( result ) )	return ARM_KO;
	CCString msg ("");
	
	try	{

		ARM_Object* obj	=  LOCAL_PERSISTENT_OBJECTS->GetObject(modelParamsId );
		ARM_CIRModelParams* modParam = dynamic_cast<ARM_CIRModelParams*>(obj);
		if( !modParam ){
			result.setMsg ("ARM_ERR: model parameter is not of a good type");
			return ARM_KO;
		}

		double res = modParam->CptExpectation_Bond();

		result.setDouble(res);
	}
	catch(Exception& x){
		x.DebugPrint();
		ARM_RESULT();
	}

	catch (...){
		result.setMsg("ARM_ERR: unrecognized failure in ARMLOCAL_CIRBondExpectation");
        return(ARM_KO);
	}
	return ARM_OK;
}

////////////////////////////////////////////
//// Function to price an ARM_CIRBondVariance
////////////////////////////////////////////

extern long ARMLOCAL_CIRBondVariance(
	const long&				modelParamsId,
	ARM_result&				result ){

	if( !GlobalPersistanceOk( result ) )	return ARM_KO;
	CCString msg ("");
	
	try	{

		ARM_Object* obj	=  LOCAL_PERSISTENT_OBJECTS->GetObject(modelParamsId );
		ARM_CIRModelParams* modParam = dynamic_cast<ARM_CIRModelParams*>(obj);
		if( !modParam ){
			result.setMsg ("ARM_ERR: model parameter is not of a good type");
			return ARM_KO;
		}

		double res = modParam->CptVariance_Bond();

		result.setDouble(res);
	}
	catch(Exception& x){
		x.DebugPrint();
		ARM_RESULT();
	}

	catch (...){
		result.setMsg("ARM_ERR: unrecognized failure in ARMLOCAL_CIRBondVariance");
        return(ARM_KO);
	}
	return ARM_OK;
}


extern long ARMLOCAL_EXP_RICCATI_Create(
	const double&			alpha,
	const double&			beta,
	const double&			delta,
	const double&			lambda,
	const double&			x0,
	const double&			x1,
	const double&			x2,
	const double&			y0,
	const double&			t0,
	ARM_result&				result, 
	long					objId ){

	if( !GlobalPersistanceOk( result ) )	return ARM_KO;
	CCString msg ("");
	
	ARM_ExpRiccati * Eq = NULL;

	try	{

		Eq	= new ARM_ExpRiccati(alpha,beta,delta,lambda,x0,x1,x2,y0,t0);

		if( !assignObject( Eq, result, objId ) )	return ARM_KO; 
		else										return ARM_OK; 

	}
	catch(Exception& x){
		x.DebugPrint();
		ARM_RESULT();
	}

	catch (...){
		result.setMsg("ARM_ERR: unrecognized failure in ARMLOCAL_EXP_RICCATI_Create");
        return(ARM_KO);
	}
	return ARM_OK;
}

extern long ARMLOCAL_EXP_RICCATI_Price(
	const long&				RiccatiEq,
	const double&			t,
	ARM_result&				result ){

	if( !GlobalPersistanceOk( result ) )	return ARM_KO;
	CCString msg ("");
	
	try	{

		ARM_Object* obj	=  LOCAL_PERSISTENT_OBJECTS->GetObject(RiccatiEq );
		ARM_ExpRiccati* Eq = dynamic_cast<ARM_ExpRiccati*>(obj);
		if( !Eq ){
			result.setMsg ("ARM_ERR: this is not a Riccati equation");
			return ARM_KO;
		}

		double price= std::real( (*Eq)[t] );
		result.setDouble(price);
	}
	catch(Exception& x){
		x.DebugPrint();
		ARM_RESULT();
	}

	catch (...){
		result.setMsg("ARM_ERR: unrecognized failure in ARMLOCAL_CIRBondDensity");
        return(ARM_KO);
	}
	return ARM_OK;
}



extern long ARMLOCAL_EXP_RICCATI_Int(
	const long&				RiccatiEq,
	const double&			t,
	const double&			T,
	ARM_result&				result ){
	
	if( !GlobalPersistanceOk( result ) )	return ARM_KO;
	CCString msg ("");
	
	try	{

		ARM_Object* obj	=  LOCAL_PERSISTENT_OBJECTS->GetObject(RiccatiEq );
		ARM_ExpRiccati* Eq = dynamic_cast<ARM_ExpRiccati*>(obj);
		if( !Eq ){
			result.setMsg ("ARM_ERR: this is not a Riccati equation");
			return ARM_KO;
		}


		double price= std::real((*Eq)(t,T));
		result.setDouble(price);
		}

	catch(Exception& x){
		x.DebugPrint();
		ARM_RESULT();
	}

	catch (...){
		result.setMsg("ARM_ERR: unrecognized failure in ARMLOCAL_CIRBondDensity");
        return(ARM_KO);
	}
	return ARM_OK;
}




///
extern long ARMLOCAL_Util_TriSABR_Eigenvalues(
			double rho1,   double rho2,   double rho3,
			double rhos12, double rhos23, double rhos13,
			double rhov12, double rhov23, double rhov13,
			double rhoc12, double rhoc13,
			double rhoc21, double rhoc23,
			double rhoc31, double rhoc32,
			ARM_result& result)
{
		/// used in the MACRO ARM_RESULT
	CCString msg ("");
	
	try
	{
		double e1,e2,e3,e4,e5,e6;
		Export_TriSABR_Eigenvalues(
			 rho1,	  rho2,	   rho3,
			 rhos12,  rhos23,  rhos13,
			 rhov12,  rhov23,  rhov13,
			 rhoc12,  rhoc13,
			 rhoc21,  rhoc23,
			 rhoc31,  rhoc32,
			&e1,&e2,&e3,&e4,&e5,&e6);
		result.setLong(4);
		result.setArray(e1,0);
		result.setArray(e2,1);
		result.setArray(e3,2);
		result.setArray(e4,3);
		result.setArray(e5,4);
		result.setArray(e6,5);

		
		return ARM_OK;
	}
	catch(Exception& x)
	{
		x.DebugPrint();
		
		ARM_RESULT();
	}
}

extern long ARMLOCAL_Nonparametric_CompleteOption(
											const vector<double>& C_Strike1_vec,
											const vector<double>& C_Vol1_vec,
											const vector<double>& C_Strike2_vec,
											const vector<double>& C_Vol2_vec,
											const vector<double>& C_S1params_vec,
											const vector<double>& C_S2params_vec,
											double	C_correlation,
											double	C_maturity,
											double	C_a1,
											double	C_b1,
											double	C_k1,
											double	C_a2,
											double	C_b2,
											double	C_k2,
											double	C_nbsteps,
											double	C_algo,
											double  C_smiletype,
								
											ARM_result& result)
{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	std::vector<double>&		strike1_vec			= NULL;	
	std::vector<double>&		vol1_vec			= NULL;	
	std::vector<double>&		strike2_vec			= NULL;	
	std::vector<double>&		vol2_vec			= NULL;	
	std::vector<double>&		S1params_vec		= NULL;	
	std::vector<double>&		S2params_vec		= NULL;	
	
	try
	{
		strike1_vec			= CreateARMGPVectorFromVECTOR( C_Strike1_vec	);		
		vol1_vec			= CreateARMGPVectorFromVECTOR( C_Vol1_vec		);		
		strike2_vec			= CreateARMGPVectorFromVECTOR( C_Strike2_vec	);
		vol2_vec			= CreateARMGPVectorFromVECTOR( C_Vol2_vec		);
			
		S1params_vec		= CreateARMGPVectorFromVECTOR( C_S1params_vec	);
		const vector<double> S1params_Vec=S1params_vec->GetValues();
		double index_begin1=S1params_Vec[0];
		double index_end1=S1params_Vec[1];
		double flag_begin1=S1params_Vec[2];
		double flag_end1=S1params_Vec[3];
		double S1=S1params_Vec[4];
	
		
		S2params_vec		= CreateARMGPVectorFromVECTOR( C_S2params_vec	);
		const vector<double> S2params_Vec=S2params_vec->GetValues();
		double index_begin2=S2params_Vec[0];
		double index_end2=S2params_Vec[1];
		double flag_begin2=S2params_Vec[2];
		double flag_end2=S2params_Vec[3];
		double S2=S2params_Vec[4];

		double dResult =  Export_Nonparametric_CompleteSpreadoption(
		 strike1_vec,
		 vol1_vec,
		 strike2_vec,
		 vol2_vec,S1,S2,
		 index_begin1, index_end1, flag_begin1, flag_end1,
		 index_begin2, index_end2, flag_begin2, flag_end2,
		 C_correlation,
		 C_maturity,
		 C_a1,
		 C_b1,
		 C_k1,
		 C_a2,
		 C_b2,
		 C_k2,
		 C_nbsteps,
		 C_algo,
		 C_smiletype
		 );

		result.setDouble(dResult);
		return ARM_OK;
	}
	
	catch(Exception& x)
	{
		x.DebugPrint();
		
		ARM_RESULT();
	}
}


extern long ARMLOCAL_Nonparametric_LogVolatility(
											const vector<double>& C_Strike1_vec,
											const vector<double>& C_Vol1_vec,
											const vector<double>& C_S1params_vec,
											double	C_k,
											ARM_result& result)
{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	std::vector<double>&		strike1_vec			= NULL;	
	std::vector<double>&		vol1_vec			= NULL;	
	std::vector<double>&		S1params_vec		= NULL;	
	
	try
	{
		strike1_vec			= CreateARMGPVectorFromVECTOR( C_Strike1_vec	);		
		vol1_vec			= CreateARMGPVectorFromVECTOR( C_Vol1_vec		);		
			
		S1params_vec		= CreateARMGPVectorFromVECTOR( C_S1params_vec	);
		const vector<double> S1params_Vec=S1params_vec->GetValues();
		double index_begin1=S1params_Vec[0];
		double index_end1=S1params_Vec[1];
		double flag_begin1=S1params_Vec[2];
		double flag_end1=S1params_Vec[3];
	
		

		double dResult =  Export_NonParametric_LogVolatility(
		 strike1_vec,vol1_vec,index_begin1, index_end1, flag_begin1, flag_end1,C_k
		 );

		result.setDouble(dResult);
		return ARM_OK;
	}
	
	catch(Exception& x)
	{
		x.DebugPrint();
		
		ARM_RESULT();
	}
}

extern long ARMLOCAL_Nonparametric_NormalVolatility(
											const vector<double>& C_Strike1_vec,
											const vector<double>& C_Vol1_vec,
											const vector<double>& C_S1params_vec,
											double	C_k,
											ARM_result& result)
{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	std::vector<double>&		strike1_vec			= NULL;	
	std::vector<double>&		vol1_vec			= NULL;	
	std::vector<double>&		S1params_vec		= NULL;	
	
	try
	{
		strike1_vec			= CreateARMGPVectorFromVECTOR( C_Strike1_vec	);		
		vol1_vec			= CreateARMGPVectorFromVECTOR( C_Vol1_vec		);		
			
		S1params_vec		= CreateARMGPVectorFromVECTOR( C_S1params_vec	);
		const vector<double> S1params_Vec=S1params_vec->GetValues();
		double index_begin1=S1params_Vec[0];
		double index_end1=S1params_Vec[1];
		double flag_begin1=S1params_Vec[2];
		double flag_end1=S1params_Vec[3];
	
		

		double dResult =  Export_NonParametric_NormalVolatility(
		 strike1_vec,vol1_vec,index_begin1, index_end1, flag_begin1, flag_end1,C_k
		 );

		result.setDouble(dResult);
		return ARM_OK;
	}
	
	catch(Exception& x)
	{
		x.DebugPrint();
		
		ARM_RESULT();
	}
}

extern long ARMLOCAL_Nonparametric_NormalDistribution(
											const vector<double>& C_Strike1_vec,
											const vector<double>& C_Vol1_vec,
											const vector<double>& C_S1params_vec,
											 double C_S, double C_T,double	C_k,
											ARM_result& result)
{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	std::vector<double>&		strike1_vec			= NULL;	
	std::vector<double>&		vol1_vec			= NULL;	
	std::vector<double>&		S1params_vec		= NULL;	
	
	try
	{
		strike1_vec			= CreateARMGPVectorFromVECTOR( C_Strike1_vec	);		
		vol1_vec			= CreateARMGPVectorFromVECTOR( C_Vol1_vec		);		
			
		S1params_vec		= CreateARMGPVectorFromVECTOR( C_S1params_vec	);
		const vector<double> S1params_Vec=S1params_vec->GetValues();
		double index_begin1=S1params_Vec[0];
		double index_end1=S1params_Vec[1];
		double flag_begin1=S1params_Vec[2];
		double flag_end1=S1params_Vec[3];
		double dResult =  Export_NonParametric_N_Distribution(
		 strike1_vec,vol1_vec,index_begin1, index_end1, flag_begin1, flag_end1,C_S,C_T,C_k
		 );

		result.setDouble(dResult);
		return ARM_OK;
	}
	catch(Exception& x)
	{
		x.DebugPrint();
		
		ARM_RESULT();
	}
}


extern long ARMLOCAL_Nonparametric_NormalQuantile(
											const vector<double>& C_Strike1_vec,
											const vector<double>& C_Vol1_vec,
											const vector<double>& C_S1params_vec,
											 double C_S, double C_T,double	C_k,
											ARM_result& result)
{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	std::vector<double>&		strike1_vec			= NULL;	
	std::vector<double>&		vol1_vec			= NULL;	
	std::vector<double>&		S1params_vec		= NULL;	
	
	try
	{
		strike1_vec			= CreateARMGPVectorFromVECTOR( C_Strike1_vec	);		
		vol1_vec			= CreateARMGPVectorFromVECTOR( C_Vol1_vec		);		
			
		S1params_vec		= CreateARMGPVectorFromVECTOR( C_S1params_vec	);
		const vector<double> S1params_Vec=S1params_vec->GetValues();
		double index_begin1=S1params_Vec[0];
		double index_end1=S1params_Vec[1];
		double flag_begin1=S1params_Vec[2];
		double flag_end1=S1params_Vec[3];
	
		

		double dResult =  Export_NonParametric_N_Quantile(
		 strike1_vec,vol1_vec,index_begin1, index_end1, flag_begin1, flag_end1,C_S,C_T,C_k
		 );

		result.setDouble(dResult);
		return ARM_OK;
	}
	
	catch(Exception& x)
	{
		x.DebugPrint();
		
		ARM_RESULT();
	}
}

extern long ARMLOCAL_Nonparametric_LogNormalDistribution(
											const vector<double>& C_Strike1_vec,
											const vector<double>& C_Vol1_vec,
											const vector<double>& C_S1params_vec,
												 double C_S, double C_T,double	C_k,
											ARM_result& result)
{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	std::vector<double>&		strike1_vec			= NULL;	
	std::vector<double>&		vol1_vec			= NULL;	
	std::vector<double>&		S1params_vec		= NULL;	
	
	try
	{
		strike1_vec			= CreateARMGPVectorFromVECTOR( C_Strike1_vec	);		
		vol1_vec			= CreateARMGPVectorFromVECTOR( C_Vol1_vec		);		
			
		S1params_vec		= CreateARMGPVectorFromVECTOR( C_S1params_vec	);
		const vector<double> S1params_Vec=S1params_vec->GetValues();
		double index_begin1=S1params_Vec[0];
		double index_end1=S1params_Vec[1];
		double flag_begin1=S1params_Vec[2];
		double flag_end1=S1params_Vec[3];
	
		

		double dResult =  Export_NonParametric_LN_Distribution(
		 strike1_vec,vol1_vec,index_begin1, index_end1, flag_begin1, flag_end1,C_S,C_T,C_k
		 );

		result.setDouble(dResult);
		return ARM_OK;
	}
	catch(Exception& x)
	{
		x.DebugPrint();
		
		ARM_RESULT();
	}
}



extern long ARMLOCAL_Nonparametric_LogNormalQuantile(
											const vector<double>& C_Strike1_vec,
											const vector<double>& C_Vol1_vec,
											const vector<double>& C_S1params_vec,
												 double C_S, double C_T,double	C_k,
											ARM_result& result)
{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	std::vector<double>&		strike1_vec			= NULL;	
	std::vector<double>&		vol1_vec			= NULL;	
	std::vector<double>&		S1params_vec		= NULL;	
	
	try
	{
		strike1_vec			= CreateARMGPVectorFromVECTOR( C_Strike1_vec	);		
		vol1_vec			= CreateARMGPVectorFromVECTOR( C_Vol1_vec		);		
			
		S1params_vec		= CreateARMGPVectorFromVECTOR( C_S1params_vec	);
		const vector<double> S1params_Vec=S1params_vec->GetValues();
		double index_begin1=S1params_Vec[0];
		double index_end1=S1params_Vec[1];
		double flag_begin1=S1params_Vec[2];
		double flag_end1=S1params_Vec[3];
	
		

		double dResult =  Export_NonParametric_LN_Quantile(
		 strike1_vec,vol1_vec,index_begin1, index_end1, flag_begin1, flag_end1,C_S,C_T,C_k
		 );

		result.setDouble(dResult);
		return ARM_OK;
	}
	
	catch(Exception& x)
	{
		x.DebugPrint();
		
		ARM_RESULT();
	}
}
extern long ARMLOCAL_SABR2B_ImplicitVol(
									  double f, 
									  double K, 
									  double tex,
									  double alpha, 
									  double beta1, 
									  double beta2, 
									  double rho,
									  double nu,
									  double zero,
									  double lambda,
									  ARM_result& result
									  )
{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	
	try
	{
		double dResult = Export_SABR2B_ImplicitVol(f,K,tex,alpha,beta1,beta2,rho,nu,zero,lambda);
		result.setDouble(dResult);
		return ARM_OK;
	}
	
	catch(Exception& x)
	{
		x.DebugPrint();
		
		ARM_RESULT();
	}
}

extern long ARMLOCAL_BiShiftedHeston_VanillaOption(
										const double C_F1,
										const double C_V1,
										const double C_Vinfini1,
										const double C_lambda1,
										const double C_nu1,
										const double C_rho1,
										const double C_gamma1,
										const double C_F2,
										const double C_V2,
										const double C_Vinfini2,
										const double C_lambda2,
										const double C_nu2,
										const double C_rho2,
										const double C_gamma2,
										const vector<double>&  C_Correlations_vec,
										const double C_k,
										const double C_T,
										const double C_CallPut,
										const double C_LambdaB,
										const double C_Flag,
										ARM_result& result)
{
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	std::vector<double>&		Correlations_vec			= NULL;	

	
	try
	{
		Correlations_vec			= CreateARMGPVectorFromVECTOR( C_Correlations_vec	);		
		const vector<double> Correlations_Vec0=Correlations_vec->GetValues();
		double rhos=Correlations_Vec0[0];
		double rhov=Correlations_Vec0[1];
		double rhoc12=Correlations_Vec0[2];
		double rhoc21=Correlations_Vec0[3];
		double	C_nbfirst=Correlations_Vec0[4];
		double	C_nb=Correlations_Vec0[5];
		double	C_NbStage=Correlations_Vec0[6];
		double	C_NbOscill=Correlations_Vec0[7];
		double	C_prec=Correlations_Vec0[8];
			
		
		
		double dResult =  Export_BiShiftedHeston_VanillaOption(
			C_F1,
			C_V1,
			C_Vinfini1,
			C_lambda1,
			C_nu1,
			C_rho1,
			C_gamma1,
			C_F2,
			C_V2,
			C_Vinfini2,
			C_lambda2,
			C_nu2,
			C_rho2,
			C_gamma2,
			rhos,
			rhov,
			rhoc12,
			rhoc21,
			C_k,
			C_T,
			C_CallPut,
			C_LambdaB,
			C_Flag,
			C_nbfirst,
			C_nb,
			C_NbStage,
			C_NbOscill,
			C_prec
			
		 );

		result.setDouble(dResult);
		return ARM_OK;
	}
	
	catch(Exception& x)
	{
		x.DebugPrint();
		
		ARM_RESULT();
	}
}

extern long ARMLOCAL_MertonOption(
	const double& Fwd,
	const double& Strike,
	const double& t,
	const double& CallPut,
	const double& sigma,
	const double& lambda1,
	const double& U1,
	const double& lambda2,
	const double& U2,
	const int& N,
	ARM_result& result)
{
	CCString msg("");

	try
	{
		double price = MertonOption(Fwd, Strike, (int)CallPut, t, sigma, lambda1, U1, lambda2, U2, N);

		result.setDouble(price);

		return ARM_OK;
	}
	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}

}

extern long ARMLOCAL_BSImpliedVol(
	const double& F,
	const double& K,
	const double& T,
	const double& CallPut,
	const double& Target,
	ARM_result& result)
{
	CCString msg("");

	try
	{
		ARM_ImpliedVolBS func(F, K, T, (int)CallPut);

		result.setDouble(func.vol(Target));

		return ARM_OK;
	}
	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}
}
