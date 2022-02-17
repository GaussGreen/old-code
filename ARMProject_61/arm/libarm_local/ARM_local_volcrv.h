#ifndef ARM_LOCAL_VOLCRV_H
#define ARM_LOCAL_VOLCRV_H



class ARM_result;
class ARM_VolLInterpol;

#include <ARM\local_xlarm\ARM_local_interglob.h>



extern long ARMLOCAL_volcurv(const VECTOR<double>& matu,
							 const VECTOR<double>& strikes,
							 const VECTOR<double>& vols,
							 double date,
							 long strikeType,
							 long volType,
                             long interpolType,
                             const CCString& ccy,
							 const CCString& indexName,
							 long indexId,
							 ARM_result& result,
							 long objId = -1);

extern long ARMLOCAL_FxVolcurv(double date,
                               const VECTOR<double>& matu,
					           const VECTOR<double>& fxFwds,
                               const VECTOR<double>& pivotVols,
                               const VECTOR<double>& pivotTypes,
                               const VECTOR<double>& deltaPuts,
                               const VECTOR<double>& deltaCalls,
					           const VECTOR<double>& volsPuts,
                               const VECTOR<double>& volsCalls,
                               const VECTOR<double>& interpolTypes,
							   long whatIsInterpolated,
							   long correctSplineWithLinear,
                               double fxSpot,
                               long domZcCurveId,
                               long forZcCurveId,
                               int inRRSTR,
                               int isATM,
					           ARM_result& result,
					           long objId = -1);

extern long ARMLOCAL_NewComputeFxVolatility(long idCurve,
								            double moneyness, 
                                            double maturity,
								            ARM_result& result);

long ARMLOCAL_GetInfoFromFxVolatility(long idCurve, double*& res, 
                                      int &nbColumns, int &nbRows,
								      ARM_result& result);

extern long ARMLOCAL_volflat(double volFlat,
							 double date,
                             const CCString& ccy,
							 ARM_result& result,
							 long objId = -1);

extern long ARMLOCAL_ComputeVolatility(long idCurve,
									   double matu,
									   double strike,
									   double tenor,
									   ARM_result& result);

extern long ARMLOCAL_ComputeModelVolatility(long idModel,
											double matu,
											double tenor,
											double fwd,
											double strike,
											ARM_result& result,
											int useSabr = K_NO);

extern long ARMLOCAL_VolCube(long ATMVolId,
							 const VECTOR<long>& volCurveIds,
							 const VECTOR<double>& tenors,
							 long volType,
							 long checkCcy,
                             ARM_result& result,
							 long objId = -1);


extern long ARMLOCAL_GetVolFromCalypso(const string& index,
									  const string& currency,
									  const string& cvName,
									  ARM_Date date,
									  const string& vtype,
									  const string& matuIndex,
									  const string& impOrHist,
									  long indexId,
									  ARM_result& result,
									  long objId = -1);

extern long ARMLOCAL_GetVolCubeFromCalypso(const string& index,
										  const string& currency,
										  const string& cvName,
										  ARM_Date date,
										  const string& vtype,
                                          const string& type,
                                          const string suffix,
										  VECTOR<string>& tenors,
										  const string& smileOrNot,
										  long indexId,
										  ARM_result& result,
										  long objId = -1);

extern long ARMLOCAL_GetVolFromSummit(const CCString& index,
									  const CCString& currency,
									  const CCString& cvName,
									  double date,
									  const CCString& vtype,
									  const CCString& matuIndex,
									  const CCString& impOrHist,
									  long indexId,
									  ARM_result& result,
									  long objId = -1);

extern long ARMLOCAL_GetVolCubeFromSummit(const CCString& index,
										  const CCString& currency,
										  const CCString& cvName,
										  double date,
										  const CCString& vtype,
										  VECTOR<CCString>& tenors,
										  const CCString& smileOrNot,
										  long indexId,
										  int smileType,
										  ARM_result& result,
										  long objId = -1);

extern long ARMLOCAL_ARM_BumpVolatility(long VolId,
										double valueToBump,
										long nthLine,
										long nthCol,
										long cumulId,
										long absoluteId,
										ARM_result& result,
										long objId = -1);

extern long ARMLOCAL_ARM_BumpSmile(	long VolId,
									double valueToBump,
									double tenor,
									long nthLine,
									long nthCol,
									long cumulId,
									long absoluteId,
									ARM_result& result,
									long objId = -1);

extern long ARMLOCAL_ARM_BumpHyperCubeSmile(long VolId,
											double valueToBump,
											double cubeTenor,
											double smileTenor,
											long nthLine,
											long nthCol,
											long cumulId,
											long absoluteId,
											ARM_result& result,
											long objId = -1);

extern long ARMLOCAL_ARM_FXBumpRRorSTR(	long VolId,
										double valueToBump,
										long nthLine,
										long nthCol,
										double spotFX,
										long isCumul,
										long isAbsolute,
										long isRR,
										ARM_result& result,
										long objId = -1);

extern long ARMLOCAL_GetFXVolFromSummit(const CCString& ccy1,
										const CCString& ccy2,
										double date,
										const CCString& cvName,
										const CCString& impOrHist,
										const CCString& volType,
										ARM_result& result,
										long objId = -1);

extern long ARMLOCAL_GetNewFXVolFromSummit(const CCString& ccy1,
										   const CCString& ccy2,
										   double date,
										   const CCString& cvName,
										   long domZcId,
										   long forZcId,
										   double fxSpot,
										   const VECTOR<double>& forwards,
										   long WhatIsInterpolated,
										   long correctSplineWithLinear,
										   long isATM,
										   CCString& curClass,
										   ARM_result& result,
										   long objId = -1);

extern long ARMLOCAL_GetInitialFXVolFromSummit(const CCString& ccy1,
											   const CCString& ccy2,
											   double date,
											   const CCString& cvName,
											   const CCString& impOrHist,
											   const CCString& volType,
											   VECTOR<CCString> *maturities,
											   VECTOR<double> *tenors,
											   VECTOR<double> *vol,
											   ARM_result& result);

extern long ARMLOCAL_GetInitialVolFromSummit (const CCString& index,
											  const CCString& currency,
											  const CCString& cvName,
											  double date,
											  const CCString& vtype,
											  const CCString& matuIndex,
											  VECTOR<CCString> *maturities,
											  VECTOR<CCString> *tenors,
											  VECTOR<double> *vol,
											  ARM_result& result);

extern long ARMLOCAL_GetFXCorrelFromSummit(const CCString& ccy1,
										   const CCString& index,
										   const CCString& ccy2,
										   double date,
										   const CCString& cvName,
										   const VECTOR<CCString>& tenors,
										   ARM_result& result,
										   long objId = -1);

extern long ARMLOCAL_GetCorrelFromSummit(const CCString& ccy1,
										 const CCString& index1,
										 const CCString& ccy2,
										 const CCString& index2,
										 double date,
										 const CCString& cvName,
										 ARM_result& result,
										 long objId = -1);

extern long ARMLOCAL_GetNthMaturity (long idCurve,
									 long nLine,
									 ARM_result& result);


extern long ARMLOCAL_ARM_CONVERTFROMBSTONORMALVOL(long volId,
												  long zcId,
												  long isSwoptVol,
												  long inPct,
												  ARM_result& result,
												  long objId = -1);

extern long ARMLOCAL_ARM_CONVERTFROMNORMALTOBSVOL(long volId,
												  long zcId,
												  long isSwoptVol,
												  long inPct,
												  long outPct,
												  ARM_result& result,
												  long objId = -1);

extern long ARMLOCAL_ARM_CONV3FFROMSPOTTOFWDVOL(long tree3fId,
										        long PrcsId,
										        const VECTOR<double>& dForwardVolDates,
                                                ARM_result& result,
										        long objId = -1);

extern long ARMLOCAL_InterpolInStrikeFwdTime (long idCurve,
											  double forward,
											  double strike,		
											  double matu,
											  double precision,
											  double sigmaATMF,
											  long y2NULL,
											  ARM_result& result);

extern long ARMLOCAL_ComputeFxVol (	long idCurve,
									double asof,
									double matu,
									double calcmatu,
									double fxspot,
									double strike,
									long idDiscCrv,
									long idDivCrv,
									ARM_result& result);

extern long ARMLOCAL_ARM_CONV3FFROMSPOTTOFWDVOL(double asOf,
										        long dZcId,
                                                long fZcId,
                                                long dBSZcId,
				                                long fBSZcId,
                                                long volSwopBaseId,
						                        long volSwopForeignId,
                                                long fxVolId,
                                                double dMeanReversionBase,
						                        double dMeanReversionForeign,
                                                long dFxRdCorrId,
			                                    long dFxRfCorrId,
                                                long dRdRfCorrId,
                                                double dCutOff,
                                                double dVolLongTerm,
                                                long calibBasisIncluded,
                                                const VECTOR<double>& dForwardVolDates,
										        ARM_result& result,
										        long objId = -1);


extern long ARMLOCAL_ARM_CONV3FFROMFWDTOSPOTVOL(long tree3fId,
                                                long PrcsId,
                                                ARM_result& result,
                                                long objId = -1);

long ARMLOCAL_ARM_CONV3FFROMFWDTOSPOTVOL(   double asOf,
										    long dZcId,
                                            long fZcId,
                                            long dBSZcId,
				                            long fBSZcId,
                                            long volSwopBaseId,
						                    long volSwopForeignId,
                                            long fxVolId,
                                            double dMeanReversionBase,
						                    double dMeanReversionForeign,
                                            long dFxRdCorrId,
			                                long dFxRfCorrId,
                                            long dRdRfCorrId,
                                            double dCutOff,
                                            double dVolLongTerm,
                                            long calibBasisIncluded,
                                            ARM_result& result,
										    long objId = - 1);


ARM_VolLInterpol* ARMLOCAL_GetVolATMFromSummit(const CCString& index,
											   const CCString& currency,
											   const CCString& cvName,
											   ARM_Date date,
											   const CCString& vtype,
				   							   ARM_result& result);

ARM_VolLInterpol* ARMLOCAL_GetVolATMFromCalypso(const string& index,
											   const string& currency,
											   const string& cvName,
											   ARM_Date date,
											   const string& vtype,
				   							   ARM_result& result);

extern double ARMLOCAL_SetVolCurveName(long idCurve,
									   CCString name,
									   ARM_result& result);

// HyperCube
extern long ARMLOCAL_HyperCube(vector<long>& aVolCubeIds, 
							   vector<CCString>& aKeys, 
							   ARM_result& result, 
							   long objId = -1);

extern long ARMLOCAL_CreateCorrelCubeByExpiry(long C_hyperCubeId,
											  vector<CCString>& aTenorList,
											  vector<CCString>& aExpiryList,
											  CCString& aIntersurfaceInterpol,
											  ARM_result& result,
											  long objId = -1);

extern long ARMLOCAL_ComputeCorrelByExpiry(long C_correlCubeId,
										   double C_Expiry,
										   double C_Tenor1,
										   double C_Tenor2,
										   ARM_result& result);

extern long ARMLOCAL_ComputeHyperCorrel(long C_correlCubeId,
										   double C_Tenor1,
										   double C_Tenor2,
										   double C_Expiry,
										   double C_Moneyness,
										   ARM_result& result);

// IndexIndexVolCube
extern long ARMLOCAL_IndexIndexCorrelCube(vector<long>& aVolCurveIds,
										  vector<CCString>& aTenors1List,
										  vector<CCString>& aTenors2List,
										  CCString& aIntersurfaceInterpol,
										  ARM_result& result,
										  long objId = -1);

extern long ARMLOCAL_ComputeIndexIndexCorrel(long C_correlCubeId,
											 double C_Tenor1,
											 double C_Tenor2,
											 double C_Expiry1,
											 double C_Expiry2,
											 ARM_result& result);

extern long	ARMLOCAL_CreateGenCorrelatorManager(VECTOR<CCString>& mktTags,
											     vector<long>& IndexIndexVolIds, 
												 vector<long>& HyperDiagVolIds,
												 vector<long>& CorrelVolIds,
												 vector<long>& IndexVolIds,
												 vector<long>& IRVolHyperCubeIds,
												 vector<long>& VolVolHyperCubeIds,
												 vector<long>& FXVolHyperCubeIds,
												 ARM_result& result, 
												 long objId = ARM_NULL_OBJECT_ID);

extern long ARMLOCAL_ComputeIdIdCorrelFromCorrelatorManager(const string& ccy,
															const string& tenor1,
															const string& tenor2,
															double expiry1,
															double expiry2,
															long correlManagerId,
															ARM_result& result );

extern long ARMLOCAL_ComputeHyperCorrelFromCorrelatorManager(const string& ccy,
															 const string& tenor1,
															 const string& tenor2,
															 double expiry,
															 long correlManagerId,
															 ARM_result& result,
															 bool isByExpiry = false);

extern long ARMLOCAL_ComputeCorrSimplFromCorrelatorManager(const string& ccy,
													 const string& tenor,
													 double expiry,
													 long correlManagerId,
													 ARM_result& result,
													 bool isCorr = true);

extern long ARMLOCAL_ComputeFXCorrelFromCorrelatorManager(const string& ccy,
														  const string& ccy2,
														  const string& tenor,
														  double expiry,
														  long correlManagerId,
														  ARM_result& result);


extern long ARMLOCAL_BumpVolatilityCorrelManager( long correlManagerId,
												  string& C_ccy,
												  long TypeCorrel,
												  double value,
												  long nthLine,
												  long nthCol,
												  long isCumul,
												  long isAbsolute,
												  long isToClone,
												  ARM_result& result,
												  long ObjId = -1);

extern long ARMLOCAL_BumpVolatilityCorrelManager( long correlManagerId,
												  vector<string> mktTag,
												  vector<string> intraMktTag,
												  double value,
												  long nthLine,
												  long nthCol,
												  long isCumul,
												  long isAbsolute,
												  long isToClone,
												  ARM_result& result,
												  long ObjId = -1);

extern long ARMLOCAL_GetCorrelDiag(long IdIdCorrelId,
								   string& C_Tenor1,
								   vector<CCString>& C_Tenor,
								   ARM_result& result,
								   long ObjId = -1);

extern long ARMLOCAL_GetMixtureParamsFromSummit ( const CCString& index,
												  const CCString& currency,
												  const CCString& cvName,
												  double date,
												  const CCString& interpolMethod,
												  ARM_result& result,
												  long objId = -1 );

// ARM_SABRVol Functions
extern long ARMLOCAL_SABRVol(long SigmaOrAlphaId,
							 long RhoId,
							 long BetaId,
							 long NuId,
							 long SigmaOrAlphaFlag,
							 long ModelType,
							 double Weight,
							 ARM_result& result,
							 long objId = -1);

extern long ARMLOCAL_GetSABRVolFromSummit(const CCString& index,
										  const CCString& currency,
										  const CCString& cvName,
										  double date,
										  const CCString& vtype,
										  const CCString& matuIndex,
										  const CCString& impOrHist,
										  long indexId,
										  long SigmaOrAlphaFlag,
										  long ModelType,
										  double Weight,
										  ARM_result& result,
										  long objId = -1);

extern long ARMLOCAL_OldVolCurve(
	const long& volId,
	const double& asOf,
	ARM_result& result,
	long objId);


extern long ARMLOCAL_ConvIndexInYearTerm(const CCString& index,
									double AsOf, 
									const CCString& Ccy,
									ARM_result& result);

#endif 
/* ARM_LOCAL_VOLCRV_H */
