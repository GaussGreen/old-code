/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file argconvdefault.cpp
 *
 *  \brief
 *
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date October 2003
 */


/// gpbase
#include "gpbase/curve.h"

/// gpinfra
#include "gpinfra/argconvdefault.h"
#include "gpinfra/modelparamtype.h"
#include "gpinfra/numeraire.h"
#include "gpinfra/correlmatparam.h"
#include "gpinfra/mktdatacst.h"
#include "gpinfra/gramfunctorcf.h"
#include "gpinfra/gramfunctormepicf.h"



CC_BEGIN_NAMESPACE( ARM )

ARGConvTable DayCountTable[] = 
{
	/// name		/// number
	"ACTUAL",		KACTUAL_ACTUAL,
	"A365",			KACTUAL_365,
	"A360",			KACTUAL_360,
	"30/360",		K30_360,
	"ACTREAL",		KACTUAL_REAL,
	"ACT29",		KACTUAL_FEB29,
	"ISMA",			KACTUAL_ISMA,
	"30E",          K30_360E,
    "NOBASE",		KNOBASE,

	/// very important as it tells that this is the end
	ENDOFLINE_CHAR	
};

const ARM_ArgConv ARM_ArgConv_DayCount( DayCountTable, "DayCount" );
const ARM_ArgConvReverse ARM_ArgConvReverse_DayCount( DayCountTable, "DayCount" );


ARGConvTable ModelParamTable[] = 
{
	/// name					/// number
	"Volatility",				ARM_ModelParamType::Volatility,
	"MeanReversion",			ARM_ModelParamType::MeanReversion,
	"VolatilityRatio",			ARM_ModelParamType::VolatilityRatio,
	"MeanReversionSpread",		ARM_ModelParamType::MeanReversionSpread,
	"Shift",					ARM_ModelParamType::Shift,
    "Alpha",				    ARM_ModelParamType::Alpha, 
    "Beta",						ARM_ModelParamType::Beta,
   	"Correlation",				ARM_ModelParamType::Correlation,
    "VolOfVol",					ARM_ModelParamType::VolOfVol,
	"VolDrift",					ARM_ModelParamType::VolDrift,
	"Dividend",					ARM_ModelParamType::Dividend,
	"Multiplier",				ARM_ModelParamType::Multiplier,
	"Skew",						ARM_ModelParamType::Skew,
	"InitialVol",				ARM_ModelParamType::InitialVol,
	"LongTermVol",				ARM_ModelParamType::LongTermVol,
	"VolMeanReversion",			ARM_ModelParamType::VolMeanReversion,
	"JumpProba",				ARM_ModelParamType::JumpProba,
	"JumpSize",					ARM_ModelParamType::JumpSize,
	"JumpVol",					ARM_ModelParamType::JumpVol,
	"Q",						ARM_ModelParamType::QParameter,
	"LNVol",					ARM_ModelParamType::LNVol,
	"NVol",						ARM_ModelParamType::NVol,
	"QVol",						ARM_ModelParamType::QVol,
	"Drift",					ARM_ModelParamType::Drift,
	"Smile",					ARM_ModelParamType::Smile,
	"BrownianCorrelation",		ARM_ModelParamType::BrownianCorrelation,
	"ForwardAdjustment",		ARM_ModelParamType::ForwardAdjustment,
	"StrikeAdjustment",			ARM_ModelParamType::StrikeAdjustment,
	"Hump",						ARM_ModelParamType::Hump,
	"BetaCorrelation",			ARM_ModelParamType::BetaCorrelation,
	"CrossFactor",				ARM_ModelParamType::CrossFactor,
	"ReCorrelation",			ARM_ModelParamType::ReCorrelation,
	"Sigma",			        ARM_ModelParamType::Sigma,
	"ScalingVol",		        ARM_ModelParamType::ScalingVol,
	"CompoundVol",		        ARM_ModelParamType::CompoundVol,
	"Unknown",					ARM_ModelParamType::Unknown,	

	/// very important as it tells that this is the end
	ENDOFLINE_CHAR	
};

const ARM_ArgConv ARM_ArgConv_ModelParam( ModelParamTable, "ModelParam" );
const ARM_ArgConvReverse ARM_ArgConvReverse_ModelParam( ModelParamTable, "ModelParam" );


ARGConvTable ModelParamDataTypeTable[] = 
{
	/// name		    /// number
	"Values",		    ARM_ModelParamType::Values,
	"BreakPointTimes",	ARM_ModelParamType::BreakPointTimes,
	"Tenors",			ARM_ModelParamType::Tenors,
	"Strikes",			ARM_ModelParamType::Strikes,

	/// very important as it tells that this is the end
	ENDOFLINE_CHAR	
};


const ARM_ArgConv ARM_ArgConv_ModelParamDataType( ModelParamDataTypeTable, "Model Param Data Type" );
const ARM_ArgConvReverse ARM_ArgConvReverse_ModelParamDataType( ModelParamDataTypeTable, "Model Param Data Type" );


ARGConvTable NumeraireTable[] = 
{
	/// name				/// number
	"Cash",					ARM_Numeraire::Cash,
	"TerminalZc",			ARM_Numeraire::TerminalZc,
	"TerminalEventZc",		ARM_Numeraire::TerminalEventZc,
	"RollingPayment",		ARM_Numeraire::RollingPayment,
	"RollingEvent",			ARM_Numeraire::RollingEvent,
	"RollingCash",			ARM_Numeraire::RollingCash,
	"Unknown",				ARM_Numeraire::Unknown,	

	/// very important as it tells that this is the end
	ENDOFLINE_CHAR	
};

const ARM_ArgConv ARM_ArgConv_Numeraire( NumeraireTable, "Numeraire" );
const ARM_ArgConvReverse ARM_ArgConvReverse_Numeraire( NumeraireTable, "Numeraire" );


ARGConvTable CapFloorTable[] = 
{
	/// name				/// number
	"Cap",					K_CAP,
	"Floor",				K_FLOOR,
	"C",					K_CAP,
	"F",					K_FLOOR,

	/// very important as it tells that this is the end
	ENDOFLINE_CHAR	
};

const ARM_ArgConv ARM_ArgConv_CapFloor( CapFloorTable, "CapFloor" );
/// no argconv reverse because it is not a one to one mapping!


ARGConvTable CallPutTable[] = 
{
	/// name				/// number
	"Call",					K_CALL,
	"Put",					K_PUT,
	"C",					K_CAP,
	"P",					K_PUT,

	/// very important as it tells that this is the end
	ENDOFLINE_CHAR	
};

const ARM_ArgConv ARM_ArgConv_CallPut( CallPutTable, "CallPut" );


ARGConvTable CallPutTable2[] = 
{
	/// name				/// number
	"Call",					K_CALL,
	"Put",					K_PUT,

	/// very important as it tells that this is the end
	ENDOFLINE_CHAR	
};

const ARM_ArgConv ARM_ArgConv_CallPut2( CallPutTable2, "CallPut2" );
const ARM_ArgConvReverse ARM_ArgConvReverse_CallPut2( CallPutTable2, "CallPut2" );

ARGConvTable InOutTable[] = 
{
	/// name				/// number
	"IN",					K_IN,
	"OUT",					K_OUT,

	/// very important as it tells that this is the end
	ENDOFLINE_CHAR
};

const ARM_ArgConv ARM_ArgConv_InOut( InOutTable, "InOut" );
const ARM_ArgConvReverse ARM_ArgConvReverse_InOut( InOutTable, "InOut" );


ARGConvTable PayRecTable[] = 
{
	/// name				/// number
	"Pay",					K_PAY,
	"Rec",					K_RCV,
	"P",					K_PAY,
	"R",					K_RCV,

	/// very important as it tells that this is the end
	ENDOFLINE_CHAR	
};

ARGConvTable BijectivePayRecTable[] = 
{
	/// name				/// number
	"P",					K_PAY,
	"R",					K_RCV,

	/// very important as it tells that this is the end
	ENDOFLINE_CHAR	
};
const ARM_ArgConv ARM_ArgConv_PayRec( PayRecTable, "Pay/Rec" );
const ARM_ArgConvReverse ARM_ArgConvReverse_PayRec( BijectivePayRecTable, "P/R" );
/// no argconv reverse because it is not a one to one mapping!

ARGConvTable StdFrequencyTable[] = 
{
	/// name		/// number
	"A",		    K_ANNUAL,
	"S",			K_SEMIANNUAL,
	"Q",			K_QUARTERLY,
	"B",		    K_BIMONTHLY,
	"M",		    K_MONTHLY,
	"W",		    K_WEEKLY,
	"D",		    K_DAILY,
	"Z",			K_ZEROCOUPON,

	/// very important as it tells that this is the end
	ENDOFLINE_CHAR	
};

const ARM_ArgConv ARM_ArgConv_StdFrequency( StdFrequencyTable, "Std Frequency" );
const ARM_ArgConvReverse ARM_ArgConvReverse_StdFrequency( StdFrequencyTable, "Std Frequency" );

ARGConvTable MatFrequencyTable[] = 
{
	/// name		/// number
	"12M",		    K_ANNUAL,
	"1Y",		    K_ANNUAL,
	"6M",			K_SEMIANNUAL,
	"3M",			K_QUARTERLY,
	"2M",		    K_BIMONTHLY,
	"1M",		    K_MONTHLY,
	"1W",		    K_WEEKLY,
	"1D",		    K_DAILY,

	/// very important as it tells that this is the end
	ENDOFLINE_CHAR	
};


ARGConvTable MatFrequencyReverseTable[] = 
{
	/// name		/// number
	"1Y",		    K_ANNUAL,
	"6M",			K_SEMIANNUAL,
	"3M",			K_QUARTERLY,
	"2M",		    K_BIMONTHLY,
	"1M",		    K_MONTHLY,
	"1W",		    K_WEEKLY,
	"1D",		    K_DAILY,

	/// very important as it tells that this is the end
	ENDOFLINE_CHAR	
};

const ARM_ArgConv ARM_ArgConv_MatFrequency( MatFrequencyTable, "Mat Frequency" );
const ARM_ArgConvReverse ARM_ArgConvReverse_MatFrequency( MatFrequencyReverseTable, "Mat Frequency Reverse" );

ARGConvTable IndexTypeTable[] = 
{
	/// name		    /// number
	"LIBOR1M",		    LIBOR1M,
	"LIBOR2M",			LIBOR2M,
	"LIBOR3M",			LIBOR3M,
	"LIBOR6M",		    LIBOR6M,
	"LIBOR1Y",		    LIBOR1Y,
	"EURIBOR1M",		EURIBOR1M,
	"EURIBOR2M",        EURIBOR2M,
	"EURIBOR3M",        EURIBOR3M,
	"EURIBOR6M",        EURIBOR6M,
	"EURIBOR1Y",        EURIBOR1Y,

	/// very important as it tells that this is the end
	ENDOFLINE_CHAR	
};

const ARM_ArgConv ARM_ArgConv_IndexType( IndexTypeTable, "Index Type" );
const ARM_ArgConvReverse ARM_ArgConvReverse_IndexType( IndexTypeTable, "Index Type" );



ARGConvTable TimingTable[] = 
{
	/// name		/// number
	"ADV",		    K_ADVANCE,
	"ARR",			K_ARREARS,

	/// very important as it tells that this is the end
	ENDOFLINE_CHAR	
};

const ARM_ArgConv ARM_ArgConv_Timing( TimingTable, "Timing" );
const ARM_ArgConvReverse ARM_ArgConvReverse_Timing( TimingTable, "Timing" );



ARGConvTable volTypeTable[] = 
{
	/// name		/// number
	"ATM",			ARM_MarketData::MKT_ATM_VOL,
	"ATMANDSMILE",	ARM_MarketData::MKT_ATM_AND_SMILE_VOL,

	/// very important as it tells that this is the end
	ENDOFLINE_CHAR
};

const ARM_ArgConv ARM_ArgConv_VolType( volTypeTable, "VolType" );
const ARM_ArgConvReverse ARM_ArgConvReverse_VolType( volTypeTable, "VolType" );



ARGConvTable volMktTypeTable[] = 
{
	/// name		/// number
	"IRG",			ARM_MarketData::MKT_CAPORCAPLET_VOL,
	"SWOPT",		ARM_MarketData::MKT_SWAPTION_VOL,

	/// very important as it tells that this is the end
	ENDOFLINE_CHAR
};

const ARM_ArgConv ARM_ArgConv_VolMktType( volMktTypeTable, "VolMktType" );
const ARM_ArgConvReverse ARM_ArgConvReverse_VolMktType( volMktTypeTable, "VolMktType" );

ARGConvTable interpolCurveTypeTable[] = 
{
	/// name		/// number
	"LINEAR",			K_LININTERPOL_REF,
	"STEPUPRIGHT",		K_STEPUP_RIGHT,
	"STEPUPLEFT",		K_STEPUP_LEFT,

	/// very important as it tells that this is the end
	ENDOFLINE_CHAR
};

const ARM_ArgConv ARM_ArgConv_interpolCurveType( interpolCurveTypeTable, "interpolCurveType" );
const ARM_ArgConvReverse ARM_ArgConvReverse_interpolCurveType( interpolCurveTypeTable, "interpolCurveType" );


ARGConvTable cfDistTypeTable[] = 
{
	/// name		/// number
	"LN",			ARM_CFDispatcher::K_LN_DIST,
	"N",			ARM_CFDispatcher::K_Normal_DIST,

	/// very important as it tells that this is the end
	ENDOFLINE_CHAR
};

const ARM_ArgConv ARM_ArgConv_cfDistType( cfDistTypeTable, "cfDistTypeTable" );
const ARM_ArgConvReverse ARM_ArgConvReverse_cfDistType( cfDistTypeTable, "cfDistTypeTable" );


ARGConvTable cfMepiDistTypeTable[] = 
{
	/// name		/// number
	"STOBS",			ARM_CF_MepiDispatcher::K_STOCHASTIC_BLACKSCHOLES_DIST,
	"STOBS_GEO",		ARM_CF_MepiDispatcher::K_STOCHASTIC_BLACKSCHOLES_GEOMETRIC_DIST,
	"STOBS_ARI",		ARM_CF_MepiDispatcher::K_STOCHASTIC_BLACKSCHOLES_ARITHMETIC_DIST,
	"SABR",				ARM_CF_MepiDispatcher::K_SABR_DIST,
	"SABR_MC",			ARM_CF_MepiDispatcher::K_SABR_MC_DIST,

	/// very important as it tells that this is the end
	ENDOFLINE_CHAR
};

const ARM_ArgConv ARM_ArgConv_cfMepiDistType( cfMepiDistTypeTable, "cfMepiDistTypeTable" );
const ARM_ArgConvReverse ARM_ArgConvReverse_cfMepiDistType( cfMepiDistTypeTable, "cfMepiDistTypeTable" );


ARGConvTable cfGreekTypeTable[] = 
{
	/// name		/// number
	"Delta",			ARM_CFDispatcher::K_Delta,
	"Vega",				ARM_CFDispatcher::K_Vega,
	"Gamma",			ARM_CFDispatcher::K_Gamma,

	/// very important as it tells that this is the end
	ENDOFLINE_CHAR
};

const ARM_ArgConv ARM_ArgConv_cfGreekType( cfGreekTypeTable, "cfGreekTypeTable" );
const ARM_ArgConvReverse ARM_ArgConvReverse_cfGreekType( cfGreekTypeTable, "cfGreekTypeTable" );

ARGConvTable IndexClassTable[] = 
{
	/// name		/// number
	"FIXED",			K_FIXED,
	"LIBOR",			K_LIBOR,
	"CMS",				K_CMS,

	/// very important as it tells that this is the end
	ENDOFLINE_CHAR
};

const ARM_ArgConv ARM_ArgConv_IndexClass( IndexClassTable, "IndexClassTable" );
const ARM_ArgConvReverse ARM_ArgConvReverse_IndexClass( IndexClassTable, "IndexClassTable" );

ARGConvTable IntRuleTable[] = 
{
	/// name		/// number
	"ADJ",			K_ADJUSTED,
	"UNADJ",		K_UNADJUSTED,

	/// very important as it tells that this is the end
	ENDOFLINE_CHAR
};

const ARM_ArgConv ARM_ArgConv_IntRule( IntRuleTable, "IntRuleTable" );
const ARM_ArgConvReverse ARM_ArgConvReverse_IntRule( IntRuleTable, "IntRuleTable" );

ARGConvTable InterpolInfTypeTable[] = 
{
	/// name		/// number
	"CPILINEAR",		K_CPILINEAR,	
	"CPISTEPWISE",		K_CPISTEPWISE,
	"ZCLINEA",			K_ZCLINEAR,
	"ZCCTFWD",			K_ZCCTFWD,
	"CPISTEPWISESTART",	K_CPISTEPWISESTART,
	"CPISTEPWISEEND",	K_CPISTEPWISEEND,
	"CPISTEPWISEMIDDLE",K_CPISTEPWISEMIDDLE,

	/// very important as it tells that this is the end
	ENDOFLINE_CHAR
};

const ARM_ArgConv ARM_ArgConv_InterpolInfType( InterpolInfTypeTable, "InterpolInfTypeTable" );
const ARM_ArgConvReverse ARM_ArgConvReverse_InterpolInfType( InterpolInfTypeTable, "InterpolInfTypeTable" );


ARGConvTable YesNoTable[] = 
{
	/// name		/// number
	"Y",			K_ADJUSTED,
	"N",			K_UNADJUSTED,

	/// very important as it tells that this is the end
	ENDOFLINE_CHAR
};

const ARM_ArgConv ARM_ArgConv_YesNo( YesNoTable, "IntRuleTable" );
const ARM_ArgConvReverse ARM_ArgConvReverse_YesNo( YesNoTable, "IntRuleTable" );

ARGConvTable BumpParamTable[] = 
{
	/// name		/// number
	"ISCUMULATIVE",			ARM_ModelParamBump::isCumulative,
	"ISPERTURBATIVE",		ARM_ModelParamBump::isPerturbative,
	"ISNOTHING",			ARM_ModelParamBump::isNothing,
		

	/// very important as it tells that this is the end
	ENDOFLINE_CHAR
};

const ARM_ArgConv ARM_ArgConv_BumpParam( BumpParamTable, "BumpParamTable" );
const ARM_ArgConvReverse ARM_ArgConvReverse_BumpParam( BumpParamTable, "BumpParamTable" );



ARGConvTable FwdRuleTable[]=
{
	/// methodFlag				methodName
		"P",		K_PREVIOUS		, 
		"MP",		K_MOD_PREVIOUS	, 
		"F",		K_FOLLOWING		, 
		"MF",		K_MOD_FOLLOWING	, 
		
	/// very important as it tells that this is the end
	ENDOFLINE_CHAR
};

const ARM_ArgConv ARM_ArgConv_FwdRule( FwdRuleTable, "Fwd Rule" );
const ARM_ArgConvReverse ARM_ArgConvReverse_FwdRule( FwdRuleTable, "Fwd Rule" );

ARGConvTable DigitTypeTable[]=
{
	/// methodFlag				methodName
		"ANALYTIC",		ARM_FXDigitType::analytic, 
		"CENTRED",		ARM_FXDigitType::centred, 
		"BACKWARD",		ARM_FXDigitType::backward, 
		"FORWARD",		ARM_FXDigitType::forward, 
		
	/// very important as it tells that this is the end
	ENDOFLINE_CHAR
};

const ARM_ArgConv ARM_ArgConv_DigitType( DigitTypeTable, "Digital Type" );
const ARM_ArgConvReverse ARM_ArgConvReverse_DigitType( DigitTypeTable, "Digital Type" );

CC_END_NAMESPACE()


/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

