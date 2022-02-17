#ifndef ICM_LOCAL_CURVE_H
#define ICM_LOCAL_CURVE_H

class ARM_result;


extern long ICMLOCAL_SurvivalProba (long idCurve, 
									double matu,
									ARM_result& result);


extern long ICMLOCAL_DefaultIntensity (long idCurve, 
									   double matu, 
									   long meth,
									   ARM_result& result);


extern long ICMLOCAL_DefaultProba (long idCurve, 
								   double matu,
								   ARM_result& result);

double ICMLOCAL_CptImplicitSpreadInterpol (long defcurveId, 
										   const std::string& plot); 
double ICMLOCAL_CptImplicitSpreadInterpol (long defcurveId, 
										   const ARM_Date& date); 
long ICMLOCAL_createFlatCurve(long defcurveId,const ARM_Date&date,long prevId); 
long ICMLOCAL_createFlatCurve(long defcurveId,const std::string& ,long prevId); 
long ICMLOCAL_createDefCurveFromBase(long defcurveIdCDS, long defcurveIdIndex ,const ARM_Vector& vBase,long prevId);

extern long ICMLOCAL_CptImplicitSpreadInterpolOLD (long defcurveId, 
												CCString plot,
												double Date,
												double slope,
												double  Exactdate,
												VECTOR<double>& output,
												ARM_result& result);
long
ICMLOCAL_GetDPFromCalypso (const ARM_Date& AsOf,
							const std::string& issuer,
							const std::string& Seniority,
							const std::string& Currency,
							const std::string& PricingEnv,
							const std::string& ForceCurveName,
							const long& PWCcurveId,	
							const std::string& label,
							const std::string& xmlFileName,
							long objId=-1) ;

extern long ICMLOCAL_GetDPFromSummit (const ARM_Date&AsOf,
							   const CCString& issuer,
							   const CCString& CurveName,
							   const long& PWCcurveId,	
							   const CCString& label,
							   ARM_result& result,
							   long objId = -1);


extern long ICMLOCAL_DPMktDataFromSummit (double AsOf,
										  const CCString& issuer,
										  const CCString& CurveName,
										  VECTOR<CCString>&	Matu,
										  VECTOR<double>& Spread,
										  VECTOR<double>& Recovery,
										  CCString& Currency,
										  CCString& IndexName,
										  long& AIMMADJ, 
										  ARM_result& result);

void  ICMLOCAL_DPMktDataFromCalypso(const ARM_Date& AsOf,
									const std::string & pricingEnv,
									const std::string & issuer,
									const std::string & seniority,
									const std::string & ccy,
									const std::string & forceCurveName,
									const std::string & xmlFile,
									std::vector<std::string>& Matu,
									std::vector<double>& Spread,
									double & Recovery,
									std::string & Currency,
									std::string & IndexName,
									qCDS_ADJ & AIMMADJ ) ; 

extern long ICMLOCAL_CorrFromSummit(double C_AsOfDate,
									const CCString& C_IssuerId,
									const CCString& C_CurveName,
									const CCString& C_Ccy,	
									ARM_result& C_result);




extern	long ICMLOCAL_LinearDefaultCurve (const ARM_Date& AsOfDate, 
							   const std::vector<std::string>&  Tenors, 
							   VECTOR<double>& Rates,
							   double recovery,
							   long ircurveid,
							   const CCString& ccy, 
							   const CCString& Label, 
							   qCDS_ADJ AdjCalType,
							   bool IsSummitCurve,
							   ARM_result& result, 
							   long objId = -1);

long ICMLOCAL_ConstantDefaultCurve (const ARM_Date& AsOfDate, 
								const std::vector<std::string>& Tenors, 
							   VECTOR<double>& Rates,
							   double recovery,
							   long ircurveid,
							   int intRule,
							   int adjStartRule,
							   const std::string& ccy, 
							   const std::string& Label,
							   qCDS_ADJ AdjCalType,
							   bool IsSummitCurve,
								qDEFCURVE_CALIB_ALGO calibrationAlgo,
							   long VolCurveId,
							   const std::string& calibrationMethod,
							   int lag,
							   long paramId,
							   long objId);

/**
long ICMLOCAL_DefaultCurveIndexPWC (const ARM_Date& AsOfDate, 
							   const std::vector<std::string>&  Tenors, 
							   VECTOR<double>& Rates,
							   VECTOR<double>& RefRates,
							   double recovery,
							   long ircurveid,
							   const CCString& ccy, 
							   const CCString& Label,
							   qCDS_ADJ AdjCalType,
							   bool IsSummitCurve,
							   qDEFCURVE_CALIB_ALGO calibAlgo,
							   long VolCurveId,
							   const std::string& calibrationMethod,
							   int lag,
							   // ARM_result& result, 
							   long objId = -1);
**/ 
extern	long ICMLOCAL_FwdSpread (long DefCurveId, 
								 double matu1,
								 double matu2,
								 double fwdstart,
								 double fwdend,
								 long VolId, 
								 ARM_result& result);

extern long ICMLOCAL_GetZC_DP_FromSummit (const CCString& Issuer,
								   const CCString& currency,
								   const CCString& cvName,
								   double aSdate,
								   long ircurveid,	
								   const CCString& Label,
								   ARM_result& result,
								   long objId = -1);

long ICMLOCAL_GetZC_DP_FromCalypso(const ARM_Date& AsOfDate,
								   const std::string& pricingEnv,
								   const std::string& issuer,
								   const std::string& seniority,
								   const std::string& ccy,
								   const std::string& forceCurveName,
								   const std::string& xmlFilename,
								   long ircurveid,
								   const std::string& forceLabel, 
								   long objId=-1) ; 

/**
long ICMLOCAL_CstDefCurve (const ARM_Date& AsOfDate, 
						   const std::vector<double>& Yearfractions, 
								  const std::vector<double>& Inputs,
								  double Type ,
								  double recovery,
								  long ircurveid,
								  int intRule,
								  int adjStartDateRule,
								  const std::string & ccy, 
								  const std::string & Label,
								  long VolCurveId,
								  const std::string& calibrationMethod,
								  int lag,
								  long objId = -1);
								  **/ 

long ICMLOCAL_CstDefCurve_Dates (const ARM_Date& AsOfDate, 
										VECTOR<double>& Dates, 
										VECTOR<double>& Inputs,
										double Type ,
										double recovery,
										long ircurveid,
										const CCString& ccy, 
										const CCString& Label,
										long VolCurveId,
										const std::string& calibrationMethod,
										int lag,
										long paramId,
										// ARM_result& result, 
										long prevId );


extern long ICMLOCAL_InputDefCurve_Dates (const ARM_Date&  AsOfDate, 
										  const std::vector<ARM_Date>&   Dates,
										const ARM_Vector&   Inputs,
										double recovery,
										long ircurveid,
										const std::string & ccy, 
										const std::string & Label,
										// qINTERPOL_TYPE InterpolType, 
										// ARM_result& result, 
										long prevId= -1);

extern long ICMLOCAL_ImpliedLossTree(double AsOfDate, 
							  long DefCurveId,
							  long VolAsCorrelId,
							  long correltype,
							  bool adjusted,
							  int step,
							  ARM_result& result, 
							  long objId = -1);

extern long ICMLOCAL_ModifyLossTree (long TreeId,
									 long indexno,
							  const VECTOR<double>& yearterms,
							  const VECTOR<double>& losses,
							  const VECTOR<double>& transprobas,
							  const VECTOR<double>& statelosses);


extern long ICMLOCAL_DefProbInverse (long idCurve, 
							 double DefaultProba,
							 ARM_result& result);


extern long ICMLOCAL_QuickDefaultCurve (double spread,
								 double recovery,
								 CCString label,
							     ARM_result& result, 
							     long objId=-1);

extern long ICMLOCAL_Fixing_Curve (const std::vector<ARM_Date>& vDates,
								   const std::vector<double>& vValues,
								   const ARM_Date& asOf,
								   const string& indexName,
								   long indexNameID,
							       ARM_result& result, 
							       long objId=-1);

long ICMLOCAL_ABS_PWCDefaultCurve (const ARM_Date& AsOfDate, 
								const std::vector<std::string>& Tenors, 
							   VECTOR<double>& Rates,
							   double recovery,
							   long ircurveid,
							   const std::string& ccy, 
							   const std::string& Label,
							   qCDS_ADJ AdjCalType,
							   bool IsSummitCurve,
							   // bool IsBrentCalibration,
								qDEFCURVE_CALIB_ALGO calibrationAlgo,
							   long VolCurveId,
							   const std::string& calibrationMethod,
							   int lag,
							   VECTOR<double>& Upfront, 
							   VECTOR<long>& RefValId, 
							   long objId = -1);


#endif	// ARM_LOCAL_ZCCURVE_H

// EOF %M%
