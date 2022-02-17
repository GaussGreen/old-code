#ifndef ICM_LOCAL_SWAP_H
#define ICM_LOCAL_SWAP_H


#include "ICMKernel/glob/icm_enums.h"
class ARM_result;



extern long ICMLOCAL_CDS   (double FixedRate,
					 double	EffectiveDate,
					 double EndDate,
					 double ReferenceDate,
					 double FirstCpnEffDate,
					 int	Frequency,
					 int	DayCount,
					 double FixedPayerAmount,
					 double FloatingPayerAmount,
					 int	Stubrule,
					 CCString Currency,
					 int	intRule,
					 int	CreditLag,
					 bool	IncludeMaturity,
					 double ProtectionStartDate,
					 double ProtectionEndDate,
					 CCString name,
					 double binary,
					 int	intStartAdj,
					 qPAYMENT_PREMIUM_LEG q_accuredOnDef,
					long l_NotionalEch_Type, 
					long l_NotionalEchange,
					 ARM_result& result,
					 long objId = -1);

extern long ICMLOCAL_CMCDS  (double PartRate,
							 double	EffectiveDate,
							 double EndDate,
							 double FirstCpnEffDate,
							 int	IndexId,
							 double ReferenceDate,
							 int	Frequency,
							 int	DayCount,
							 double FixedPayerAmount,
							 double FloatingPayerAmount,
							 int	Stubrule,
							 CCString Currency,
							 int	intRule,
							 int	CreditLag,
							 bool	IncludeMaturity,
							 double ProtectionStartDate,
							 double ProtectionEndDate,
							 ARM_result& result,
							 long objId = -1);

extern long ICMLOCAL_CAPFLOORCMCDS  (double PartRate,
									 double	EffectiveDate,
									 double EndDate,
									 double CapLevel,
									 double FloorLevel,
									 int	IndexId,
									 double ReferenceDate,
									 int	Frequency,
									 int	DayCount,
									 double FixedPayerAmount,
									 double FloatingPayerAmount,
									 int	Stubrule,
									 CCString Currency,
									 int	intRule,
									 int	CreditLag,
									 bool	IncludeMaturity,
									 double ProtectionStartDate,
									 double ProtectionEndDate,
									 ARM_result& result,
									 long objId = -1);


extern long ICMLOCAL_CDSGEN(int FeelegId,
					 int	DeflegId,
					 double RcvFee,
					 double TradedNotional,
					 ARM_result& result,
					 long	objId = -1);

extern long ICMLOCAL_FTD   (double FixedRate,
					 double	EffectiveDate,
					 double EndDate,
					 double ReferenceDate,
					 double FirstCpnEffDate,
					 VECTOR<CCString> labels,
					 int	Frequency,
					 int	DayCount,
					 double IssuerNotional,
					 qPAYMENT_PREMIUM_LEG	AccruedOnDefault,
					 CCString Currency,
					 int	CreditLag,
					 int	stub,
					 int intRule,
					 int startAdj,
					 ARM_result& result,
					 long	objId = -1);

extern long ICMLOCAL_NTHTD   (double FixedRate,
					 double	EffectiveDate,
					 double EndDate,
					 double ReferenceDate,
					 double FirstCpnEffDate,
					 int	FirstNumDefault,
					 int	LastNumDefault,
					 VECTOR<CCString> labels,
					 int	Frequency,
					 int	DayCount,
					 double IssuerNotional,
					 qPAYMENT_PREMIUM_LEG AccruedOnDefault,
					 CCString Currency,
					 int	CreditLag,
					 int	stub,
					 int	freqdefleg,
					 double Binary,
					 CCString PayCal,
					 double RcvFee,
					 double TradedNotional,
					 bool	IncludeMaturity,
					 int intRule,
					 int startAdj,
					 ARM_result& result,
					 long objId = -1);

extern long ICMLOCAL_MEZZANINE (double FixedRate,
					 double	EffectiveDate,
					 double EndDate,
					 double ReferenceDate,
					 double FirstCpnEffDate,
					 double MezzAmount,
					 double SubAmount,
					 VECTOR<CCString> labels,
					 VECTOR<double> notionals,
					 int	FreqFeeLeg,
					 int	FreqDefLeg,
					 int	DayCount,
					 double FixedPayerAmount,
					 double FloatingPayerAmount,
					 qPAYMENT_PREMIUM_LEG AccruedOnDefault,
					 CCString Currency,
					 int	CreditLag,
					 int	stub,
					 double Binary,
					 CCString PayCal,
					 double RcvFee,
					 double TradedNotional,
					 long TypeFeeLeg,
					 long TypeDefLeg,
					 bool IncludeMaturity,
					 int intRule,
					 int adjStartDate,
					 ARM_result& result,
					 long objId = -1);

extern long ICMLOCAL_CMTranche(
					 double	EffectiveDate,
					 double EndDate,
					 double ReferenceDate,
					 double MezzAmount,
					 double SubAmount,
					 VECTOR<CCString> labels,
					 VECTOR<double> notionals,
					 int	IndexId,
					 double PartRate,
					 int	FreqFeeLeg,
					 int	FreqDefLeg,
					 int	DayCount,
					 double FixedPayerAmount,
					 double FloatingPayerAmount,
					 qPAYMENT_PREMIUM_LEG AccruedOnDefault,
					 CCString Currency,
					 int CreditLag,
					 int stub,
					 double Binary,
					 CCString PayCal,
					 double RcvFee,
					 double TradedNotional,
					 double FwdFixedDate,
					 bool	IncludeMaturity,
					 ARM_result& result,
					 long objId = -1);

extern long ICMLOCAL_CDO2   (double Spread,
					 double	EffectiveDate,
					 double EndDate,
					 double ReferenceDate,
					 double FirstCpnEffDate,
					 double MezzAmount,
					 double SubAmount,
					 int	FreqFeeLeg,
					 int	DayCount,
					 qPAYMENT_PREMIUM_LEG AccruedOnDefault,
					 CCString Currency,
					 long	PtfId, 	
					 int	CreditLag,
					 int	stub,
					 int	FreqDefLeg,
					 double Binary,
					 CCString PayCal,
					 double RcvFee,
					 double TradedNotional,
					 bool   CrossSubordination,
					 bool	IncludeMaturity,
					 ARM_result& result,
					 long	objId = -1);

extern long ICMLOCAL_CDSIndex(double FixedRate,
							  int		IndexId,
							  double	EffectiveDate,
							  double	EndDate,
							  double	ReferenceDate,
							  int		Frequency,
							  int		DayCount,
							  double	FixedPayerAmount,
							  double	FloatingPayerAmount,
							  int		Stubrule,
							  CCString	Currency,
							  int		intRule,
							  int		CreditLag,
							  bool		IncludeMaturity,
							  double	ProtectionStartDate,
							  double	ProtectionEndDate,
							  ARM_result& result,
							  long objId = -1);

extern long ICMLOCAL_Option(string UnderlyingMaturity,
							ARM_Date UnderlyingMaturityDate,
							double Expiry,
							string ccy,
							qCDS_ADJ cds_adj,
							bool endAdj,
							double Strike,
							int OptionType,
							qDEF_MAT KO,
							double Notional,
							ARM_result& result,
							long objId=-1);

extern long ICMLOCAL_SpreadOption(long idCDS,
							double Strike,
							bool isCall,
							int,//qKoStyle,
							int, //qAccelerationStyle ,
							const std::vector<double>& exerciseDates,
							int exerciseStyle,
							int exerciseFrequency,
							int,//qUnderlying_Maturity_Style matuStyle,
							ARM_result& result,
							long objId) ;





extern long ICMLOCAL_CPPI  (double	StartDate,
							double MaturityDate,
							long idSecurity,
							CCString Currency,
							vector<double> Min,
							vector<double> Max,
							vector<double> ValueMin,
							vector<double> ValueMax,
							double Notional,
							double ProtectedAmount,
							double ManagementCost,
							double AdditionalLeverage,
							double DesactivateCushion,
							CCString CorrelName,
							ARM_result& result,
							long objId = -1);


extern long ICMLOCAL_NTDGEN(const int& CdsId,
					 const int& FirstNumDefault,
					 const int& LastNumDefault,
					 const int& CollateralId,
					 const double& Binary,
					 const double& RcvFee,
					 ARM_result& result,
					 long	objId = -1);

extern long ICMLOCAL_CDOGEN(const int& CdsId,
					 const double& SubAmountId,
// FIXMEFRED: mig.vc8 (30/05/2007 18:00:45):int defautl used
					 const int strike_type,
					 const int& CollateralId,
					 const double& Binary,
					 const double& RcvFee,
					 ARM_result& result,
					 long	objId = -1);

extern long ICMLOCAL_getLastFixingDate(long instId,double asofDate,ARM_result&result) ;
extern long ICMLOCAL_setPastFixing(long instId,double resetDate,double fixingValue,ARM_result&result) ;

extern long ICMLOCAL_SetPricerForRatesComputation(long legId,
									 long pricerId,
									 ARM_result& result);

extern long ICMLOCAL_GetBounds(long CdoId,
							   double& low,
							   double& up,
							   ARM_result& result);

extern long ICMLOCAL_Customized_CDO(
							const VECTOR<CCString>&	Labels,
							const VECTOR<double>&	notionals,
							const CCString			Currency,
							const long&				CP_DefaultLegId,
							const long&				CP_PremiumLegId,
							const long&				CP_ProductParametersId,
							ARM_result&			result,
							long				objId = -1);

/**
extern long ICMLOCAL_Option_Index_Gen(
							const long&		ScheduleParametersId,
							const long&		DataParametersId,
							ARM_result&		result,
							long			objId = -1);


extern long ICMLOCAL_Corridor_Index_Gen(
							const long&		ScheduleParametersId,
							const long&		DataParametersId,
							const long&		SubScheduleParametersId,
							ARM_result&		result,
							long			objId = -1);
							*/ 


/*
extern long ICMLOCAL_Get_Instrument_Data(
							long				InstrumentId,
							VECTOR<double*>&	OutputMatrix,
							VECTOR<CCString>&	OutputLabels,
							int&				OutputNbRows,
							int&				OutputNbCols,
							ARM_result&			result);
*/


extern long ICMLOCAL_LssGapOption(long MezzId,
							long nblines_spread_triggers,
							long nbcols_spread_triggers,
							vector<double> spread_triggers,
							long nblines_default_triggers,
							long nbcols_default_triggers,
							vector<double> default_triggers,
							long nblines_mtm_triggers,
							long nbcols_mtm_triggers,
							vector<double> mtm_triggers,
							double	mtm_single_cond,
							ARM_result& result,
							long objId=-1);

extern long ICMLOCAL_CPDO (const long&			RiskyLegId,
					 const long&		RollLegId,
					 const long&		NoRiskyLegId,
					 const double&		CPDOMaturity,
					 const double&		InitialValo,
					 const double&		Target,
					 const CCString&	CouponType,
					 const double&		UFFees,
					 const double&		RunningFees,
					 const double&		VExpo,
					 const double&		V0Expo,
					 const double&		Alpha,
					 const double&		Beta,
					 const double&		Desactivation,
					 const double&		NbAssets,	
					 ARM_result& result,
					 long objId=-1);

extern long ICMLOCAL_RESTRIKABLE_CDO(const ARM_Date* TriggerStartDate,
								const ARM_Date* Expiry ,
								double Strike,
								double InitSpread,
								int OptionType,
								long idSecurity,
								double Rehauss,
								double Frequency,
								int DiffCDO,
								int IsCMSpread,
								double CMSpreadMatu,
								long objId = -1);

#endif	// ARM_LOCAL_SWAP_H

/*----End Of File ----*/
// EOF %M%
