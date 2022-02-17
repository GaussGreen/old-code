#ifndef ICM_LOCAL_GLOB_H
#define ICM_LOCAL_GLOB_H

#include <ARM\libarm\ARM_result.h>
#include "ICMKernel\util\icm_matrix.h"
#include "ICMKernel\glob\icm_enums.h"

extern std::string  ICMLOCAL_Version();

extern long ICMLOCAL_Spread(long pricerId, double MtM, ARM_result& result);

extern long ICMLOCAL_Price (long pricerId, double AsOfDate, ARM_result& result);

// extern long ICMLOCAL_NPV(long pricerId, qCMPMETH measure, ARM_result& result);
double ICMLOCAL_NPV(long pricerId, qCMPMETH measure); 

extern double ICMLOCAL_GetBeta (long pricerId, CCString label, ARM_result& result);

extern double ICMLOCAL_GetDefProbTranche (long pricerId, double yearterm, ARM_result& result);

extern double ICMLOCAL_VirtualCdsSpread (long pricerId, double MaturityDate, ARM_result& result);

extern double ICMLOCAL_GetDuration (long pricerId, ARM_result& result);

extern long ICMLOCAL_GetLabel (long curveId, ARM_result& result);

extern long ICMLOCAL_SetLabel (long curveId, CCString label, ARM_result& result);

extern long ICMLOCAL_Sensitivity (long pricerId ,qSENSITIVITY_TYPE curvtype, const std::string& plot, const std::string& label, double epsilon, double epsilonGamma, ARM_result& result);

extern long ICMLOCAL_DTR (long pricerId , int NbDefaults,CCString S_or_L , VECTOR<CCString>& Labels, VECTOR<double> RecoveryRates , ARM_result& result );

// extern double ICMLOCAL_CorrelationSmile (long pricerId, int smiletype, double mktdata, double seed, double Upfront,int datatype, ARM_result& result);

// extern double ICMLOCAL_Get_MixCopula_Factor (long pricerId, int SmileType, double BC1, double BC2, double Seed1, double Seed2, double C_Accuracy, int FactorType, ARM_result& result);

// extern double ICMLOCAL_Get_ReducedMixCopula_Factor (long pricerId, int SmileType, double BC1, double BC2, double IndepSeed, double BetaSeed, int FactorType, ARM_result& result);

// extern double ICMLOCAL_GetBaseCorrelation (long pricerId, double EquityAmount,VECTOR<double>& MktSpreads, double smiletype,VECTOR<double>& Seeds,VECTOR<double>& Upfronts, ARM_result& result);

// extern double ICMLOCAL_GetImpliedCorrelation (long pricerId, double Range,VECTOR<double>& MktSpreads, double smiletype,VECTOR<double>& Seeds, VECTOR<double>& Upfronts, ARM_result& result);

// extern long ICM_ConvCurvType (const CCString& aOptType, ARM_result& result);

extern long ARM_Credit_DisplayScheduleDates (long instId, long datesType, long recId, ARM_result& result);

extern long ARM_Credit_DisplayScheduleValues (long pricerId, long valuesType, long recId, ARM_result& result);

extern long ICMLOCAL_DisplayMatrix (long pricerId, VECTOR<CCString>& VectorOut,int& NbCol,VECTOR<CCString>& ColNames, ARM_result& result);

/**extern long ICMLOCAL_SensitivityFull (long pricerId, VECTOR<CCString>& VectorOut, int& NbCol, VECTOR<CCString>& VectorName, ARM_result& result);
**/ 

double ICMLOCAL_RiskyDuration (long DefCurveId, const ARM_Date& date);
double ICMLOCAL_RiskyDuration (long DefCurveId, const std::string& tenor);

extern long ICMLOCAL_RiskyDuration (long ModelId,CCString issuer, double date, CCString Tenor, ARM_result& result);

// extern long ICMLOCAL_RiskyPV01 (long DefCurveId, double date1,double date2 ,CCString& tenor, ARM_result& result) ;
double ICMLOCAL_RiskyPV01 (long DefCurveId, const ARM_Date* date1,const ARM_Date& date2 ); 
double ICMLOCAL_RiskyPV01 (long DefCurveId, const ARM_Date* date1,const std::string& tenor ); 
double ICMLOCAL_RiskyPV01AsSensitivity (long DefCurveId, const std::string& tenor ); 


extern long ICMLOCAL_Delivery (double AsOfDate,CCString maturity,ARM_result& result);

extern long ICMLOCAL_CORRMATRIX (const ARM_Date& AsOfDate,
								 const std::string& Name,
								 const std::vector<std::string>& labels,
								 const VECTOR<double>& coefs, 
								 ARM_result& result, 
								 long objId = -1);

extern long ICMLOCAL_EXTRACTCORRMATRIX (long CorrMatrixId,
										const std::vector<std::string> &labels,
								 VECTOR<double>& coefsout, 
								 ARM_result& result);

extern long ICMLOCAL_GenSchedule (double	EffectiveDate,
								double EndDate,
								double ReferenceDate,
								int	Frequency,
								int	DayCount,
								CCString Currency,
								long datesType,
								int modfoll,
								int gapcredit,
								ARM_result& result);

/**
extern long ICMLOCAL_ZCPrice(double ValoSpread,
							 double Maturity,
							 double Recovery,
							 double NotionalAmount,
							 ARM_result& result);

extern long ICMLOCAL_ZCFDPrice(VECTOR<double>& ValoSpreads,
							   double Maturity,
							   double Beta,
							   double Recovery,
							   double NotionalAmount,
							   ARM_result& result);

extern long ICMLOCAL_ZCNDPrice(VECTOR<double>& ValoSpreads,
							   double Maturity,
							   double Beta,
							   int Counter,
							   double Recovery,
							   double NotionalAmount,
							   ARM_result& result);

extern long ICMLOCAL_ZCNDDelta(VECTOR<double>& ValoSpreads,
							   double Maturity,
							   double Beta,
							   int Counter,
							   int Underl,
							   double Recovery,
							   ARM_result& result);

extern long ICMLOCAL_ZCHdgNDGamma(VECTOR<double>& ValoSpreads,
							   VECTOR<double>& SpreadsAnnVol, 
							   double SpreadBeta,
							   double Maturity,
							   double Beta,
							   int Counter,
							   double Recovery,
							   double NotionalAmount,
							   ARM_result& result);

extern long ICMLOCAL_ZCNDTimeShift(VECTOR<double>& ValoSpreads, 
							   double Maturity,
							   double Beta,
							   int Counter,
							   double Recovery,
							   double NotionalAmount,
							   ARM_result& result);


**/ 
extern long ICMLOCAL_GetCleanPriceFromSummit (double date,
									   const CCString& BondSummitId,
									   const CCString& SummitId,
									   const CCString& currency,
									   const CCString& SummitType,
									   ARM_result& result);


extern long ICMLOCAL_GetASWFromSummit (double date,
									   const CCString& BondSummitId,
									   const CCString& SummitId,
									   const CCString& currency,
									   const CCString& SummitType,
									   ARM_result& result);

extern long ICMLOCAL_GetEquityPriceFromSummit (double date,
									   const CCString& EquityId,
									   const CCString& SummitId,
									   const CCString& currency,
									   ARM_result& result);

extern long ICMLOCAL_GetEquityVolFromSummit (double date,
									   const CCString& EquityId,
									   const CCString& SummitId,
									   const CCString& currency,
									   ARM_result& result);

extern long ICMLOCAL_SetCorrelationMatrix (long ModelMultiCurvesId,
									long CorrMatrixId,
									ARM_result& result);

extern long ICMLOCAL_CloneCorrMatrixBary(int CorrMatrixId,
										  double Beta,
										  int UpOrDown,	
										  ARM_result& result,
										  long objId = -1);

extern long ICMLOCAL_SetInterpolationType (long VolCurveId,
										   long InterpolType,
										   ARM_result& result);

extern long ICMLOCAL_BSGreeks (long pricerId,
							   long greektype,
							   ARM_result& result);

extern long ICMLOCAL_Greeks (long pricerId,
					         long greektype,
							 ARM_result& result);

extern long ICMLOCAL_ImpliedVol(long pricerId,
							    double Price,
								ARM_result& result);

extern long ICMLOCAL_FwdSpreadPricer(long pricerId,
									 const ARM_Date& matu1,
									 const ARM_Date& matu2,
									 ARM_result& result);

long ICMLOCAL_CORRELATION_STRIKE(const ARM_Date &AsOfDate,
								 const std::string& Name,
								 const std::vector<std::string>& labels,
									const std::vector<long>& volcurves, 
									const std::vector<double>& proportions,
									const std::vector<double>& smilestrikelow,
									const std::vector<double>& smilestrikehight,
									const ICM_QMatrix<double>& fullStrikeLow,
									const ICM_QMatrix<double>& fullStrikeUp,
									const std::vector<long>& IndexIds,
									ARM_result& result, 
									long objId = -1);
long ICMLOCAL_CORRELATION_SMILE_STRIKE(const double& AsOfDate,
									const std::string& Name,
									const std::vector<string>& labels,
									const std::vector<long>& volcurves, 
									const ARM_Vector& proportions,
									const ARM_Vector& smilestrikelow,
									const ARM_Vector& smilestrikehight,
									const ICM_QMatrix<double>& fullStrikeLow,
									const ICM_QMatrix<double>& fullStrikeUp,
									const std::vector<long>& IndexIds,
									long objId = -1);

extern long ICMLOCAL_BETA_CORRELATION(const ARM_Date& AsOfDate,
									  const std::string& Name,
									std::vector<std::string>& labels,
									VECTOR<double>& betas,
									long idIndex1,
									long idIndex2,
									ARM_result& result, 
									long objId = -1);

long ICMLOCAL_FLAT_CORRELATION(const ARM_Date&AsOf,
							   const std::string& structName,
							   double correlValue,
							   long idIndex1,
							   long idIndex2,
							   long prevId); 

extern long ICMLOCAL_SetCorrelation (long ModelMultiCurvesId,
							  long CorrelationId,
							  ARM_result& result);


extern long ICMLOCAL_GetEqStrikeDown (long CorrelId,
							   CCString indexname,
							  ARM_result& result);

extern long ICMLOCAL_GetEqStrikeUp (long CorrelId,
							 CCString indexname,
							  ARM_result& result);

extern long ICMLOCAL_GetCorrelStrikeDown (long CorrelId,
							  double yf_maturity,
							  ARM_result& result);

extern long ICMLOCAL_GetCorrelStrikeUp (long CorrelId,
							  double yf_maturity,
							  ARM_result& result);

extern long ICMLOCAL_GetCorrelation (long ModelId,
							  ARM_result& result,
							  long objId = -1);

extern long ICMLOCAL_GetExpectedLoss(long pricerId,
									 double yearterm,
									 ARM_result& result);

extern long ICMLOCAL_FwdSpreadAsIndex (const long& defCurve,
								const double& StartDate,
								const double& EndDate,	
								ARM_result& result);

extern long ICMLOCAL_SetProportionsInfos (const long& CorrelId,
								   const CCString& IndexName,
								   const double& proportion,
								   const double& forcedstrikelow,
								   const double& forcedstrikehigh,
								   ARM_result& result);

extern long ICMLOCAL_ComputeImplicitCurveForCDO2(const long& PricerId, 
								   const CCString& Name,
								   const CCString& Tenor,
								   ARM_result& result);

extern long ICMLOCAL_Credit_AddPeriod(double AsOfDate , 
									  CCString Maturity, 
									  CCString Currency , 
									  bool Adj, 
									  qCDS_ADJ AdjCDS, 
									  ARM_result& result);

extern long ICMLOCAL_GetBaseCorrelFromSummit(const ARM_Date& AsOfDate , 
									  CCString Index, 
									  CCString CurveType, 
									  VECTOR<CCString> Currency, 
									  VECTOR<CCString> CvIssuerName, 
									  VECTOR<double> Proportions,
									  VECTOR<double> smilestrikelow,
									  VECTOR<double> smilestrikehight,
									  CCString CorrelName,
									  ARM_result& result,
									  int ObjectId = -1);

extern long ICMLOCAL_INDEX_CORRELATION(const double& AsOfDate,
								const CCString& Name,
								const long& CalMethod,
								const long& CreditIndexId,
								VECTOR<double>& StrikeLow,
								VECTOR<double>& StrikeHigh,
								VECTOR<double>& MktBid,
								VECTOR<double>& MktAsk,
								VECTOR<double>& UpfBid,
								VECTOR<double>& UpfAsk,
								VECTOR<double>& InitialCorrelation,
								double RFLBeta0,
								const long ParamId,
								const long MMCId,
								ARM_result& result, 
								long objId = -1);

extern double ICMLOCAL_GetDataFromLabel (long pricerId, 
										CCString Label,
										ARM_result& result);


extern long ICMLOCAL_CPT_BASE_CORRELATION(const double& AsOfDate,
								//const CCString& Name,
								const long& CalMethod,
								const long& CreditIndexId,
								VECTOR<double>& StrikeLow,
								VECTOR<double>& StrikeHigh,
								VECTOR<double>& MktBid,
								VECTOR<double>& MktAsk,
								VECTOR<double>& UpfBid,
								VECTOR<double>& UpfAsk,
								VECTOR<double>& InitialCorrelation,
								VECTOR<double>& Leverages,
								VECTOR<double>& basesCorrelation,
								const long& mmcId,
								const int& intstep,
								const int& startlag,
								const int& creditlag,
								VECTOR<long>& PrevCreditIndexId,
								VECTOR<double>& MatrixPrevBC,
								const double step,
								const qOPTIMIZE_TYPE method,
								const long ParamId,
								ARM_result& result);
								

extern long ICMLOCAL_GetEqStrike (long CorrelId,
							CCString indexname,
							int isUp,
							std::vector<double>& matu,
							std::vector<double>& strikes,
							ARM_result& result);


extern long ICMLOCAL_CptLeverageLevels(long pricerId,
									 long ParametersId,
									 long	TriggerCorrelationId,
									 VECTOR<double>		MatrixMultiplesInput,
									 VECTOR<double>		MatrixFlagsInput,
									 VECTOR<double>& VectOfLosses,
									 VECTOR<double>& VectOfMaturitiesInYF,
									 ICM_QMatrix<double>*&	MatrixMultiples,
									 ARM_result& result);

extern long ICMLOCAL_SetMatuLabel (long CurveId, 
								  VECTOR<CCString>& MatuLabels,
								  ARM_result& result);

extern long ICMLOCAL_SetRecovCoef (long SecId, 
								  const double& RecovCoef,
								  ARM_result& result);

extern long ICMLOCAL_SetFees(long securityId,
							 long RefvalueId,
							 ARM_result& result);

extern long ICMLOCAL_SECTORIAL_CORRELATION(
										   const ARM_Date& AsOf,
										   const std::string& structName,
											qTWO_FACTORS_CORRELATION_TYPE	Sectorial_Correlation_Id,
											const std::vector<std::string>&	Labels,
											const VECTOR<int>&		Sector_Membership,
											double				Intra_Sector_Correlation,
											double				Inter_Sector_Correlation,
											const VECTOR<double>&		Sector_Betas,
											const VECTOR<double>&		Sector_Lambdas,
											const VECTOR<double>&		Sector_Betas_Down,
											const VECTOR<double>&		Sector_Lambdas_Down,
											ARM_result&			result, 
											long				objId = -1);


extern long ICMLOCAL_GenTsCorrfromBaseCorr (long irCurveId, 
											long defCurveId, 
											long volcurveId, 
											int creditlag,
											ICM_QMatrix<double>*& Correls,
											ARM_result& result);

extern long ICMLOCAL_Get_IR_Curve_Moved_In_Time (long IRCurveId,
												 double MoveDate,
												ARM_result&			result, 
												long				objId = -1);

extern long ICMLOCAL_Math_Bivariate_normale(double x,double y, double rho, ARM_result& result);
extern long ICMLOCAL_Math_random_uniform(int seed, ARM_result& result);
extern long ICMLOCAL_Math_random_normal(double a,double b,int seed, ARM_result& result);
extern long ICMLOCAL_Math_Interpol(VECTOR<double>& a,
								   VECTOR<double>& b,
								   double value,
								   int type,
								   double smooth,
								   VECTOR<double>& weights,
								   int modespline,
								   int withC1condition,
								   double leftSlope,
								   double rightSlope,
								   ARM_result& result);

extern long ICMLOCAL_QMatrix(const ICM_QMatrix<double>& aQMatrix, long objId = -1);

extern long ICMLOCAL_QuickELoss(double pdefault,
							   double recovery,
							   int nbnames,
							   double strikedw,
							   double strikeup,
							   double correldw,
							   double correlup,
							   int lossesno,
							   int intstep,
							   bool lhp,
							   ARM_result& result);

extern long ICMLOCAL_QuickCDOPV(long discid,
						 long pdefid,
						 double startdate,
						 double enddate,
						 int frequency,
						 double rate,
						 double strikedw,
					     double strikeup,
					     double correldw,
					     double correlup,
					     double notfeeleg,
					     double notdefleg,
					     double recovery,
					     int nbnames,
					     int intstep,
					     bool LHP,
						 int pvtype,
						 long correlid,
						 ARM_result& result);

extern long ICMLOCAL_SCHEDULE_INFO(const double	EffectiveDate,
									const double MaturityDate,			
									const int	payFrequency,
									const int ResetFreq ,
									const int	DayCount,
									const int	Stubrule,
									const int	intRule,
									const std::string& payCalName,
									const int PayTiming,
									const int ResetTiming,
									const int fwdRule,
									const bool	IncludeMaturity,
									const int adj,
									const int	intStartAdj,
									const int AccDayCount,
									const double ReferenceDate,
									const double FirstCpnEffDate,
									const int CDSAdj,
									long objId = -1);

void ICMLOCAL_PriceVector(long pricerId, 
								 const std::string& measure, 
								 ARM_Vector& output);

extern void ICMLOCAL_GenPrice(long pricerId, 
							  const std::string& measure, 
							  double& price);

extern long ICMLOCAL_Math_CF_SpreadOption(double& yearterm,
										double& strike,
										double& correlation,
										VECTOR<double>& coefs,
										VECTOR<double>& spots,
										VECTOR<double>& vol,
										int intstep,
										ARM_result& result);

extern long ICMLOCAL_Register(long address,ARM_result& result);


extern long ICMLOCAL_Calibrator(VECTOR<long>& security_vector,
								VECTOR<double>& price_vectorBid,
								VECTOR<double>& price_vectorAsk,
								VECTOR<CCString>& parameters_vector,
								VECTOR<double>& tsparams_vector,
								long l_model,
								long pricerType,
								long l_parameters,
								long l_parameters_inf,
								long l_parameters_sup,
								int iPricingType,
								long l_Optimparameters,
								long prevId = -1);

extern long ICMLOCAL_RandomGenerator(const qRAN_GEN& RandomType,
									 long parameterId, long prevId=-1 );

void ICMLOCAL_GenerateOneRandom(long RandomGenId, double& Random);

void ICMLOCAL_GenerateRandoms(long RandomGenId, ARM_Vector& RandomVector);

void ICMLOCAL_ResetRandom(long RandomGenId);

#endif	// ARM_LOCAL_GLOB_H

// EOF %M%
