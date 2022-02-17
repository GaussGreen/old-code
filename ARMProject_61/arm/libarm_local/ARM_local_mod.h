#ifndef ARM_LOCAL_MOD_H
#define ARM_LOCAL_MOD_H

class ARM_Date;

#include <ARM\libarm\ARM_result.h>

///////////////////////////////////
/// this function creates either 
/// a YCMODEL or a Y2CMODEL
/// depending whether we specicy
/// or not a discount Curve
/// use ARM_NULL_OBJECT for
/// idDiscCurve to have a simple
/// YCMOD
///////////////////////////////////
extern long ARMLOCAL_ycmod(long idZeroCurve,
						   long idDiscCurve,
						   ARM_result& result, 
						   long objId = -1 );

extern long ARMLOCAL_bsmodel(double date,
							 double spot,
							 long dividend_type,
							 double dividend,
							 long discrate_type,
							 double discrate,
							 long volat_type,
							 double volat,
							 long typstk,
							 ARM_result& result,
							 long objId = -1);

long ARMLOCAL_NewBSModel(long yieldCurveId,
						 long volatilityId,
						 long capletVolatId,
						 long correlmanagerId,
						 long convadjmanagerId,
						 long spreadLockId,
						 long discountCurveId,
						 ARM_result& result,
						 long objId = -1);


extern long ARMLOCAL_BSMODELGEN(long yieldCurveId,
						 long spreadLockId,
						 long convadjustVolatId, 
						 long volatilityId,
						 long correlmanagerId,
						 long convadjmanagerId,//convxMngrId
						 long discountCurveId,
						 long correlId,
						 long cashVolId,
						 long spreadVolId,
						 long modelTypeId,
						 long spreadVolTypeId,
						 long sabrModId,
                         bool isLnVol,
						 long numSteps,
						 ARM_result& result,
						 long objId = -1);

extern long ARMLOCAL_BSPricingModel(long yieldCurveId,
									long convexityManagerId,
									long capModelId, 
									long swoptModelId,
									long correlManagerId,
									long modelTypeId,
									long spreadVolCurveId,
									long discountCurveId,
                                    long adjConvVolCurveId,
									const VECTOR<double>& calibInfos,
									const VECTOR<double>& numInfos,
									ARM_result& result,
									long objId = -1);

extern long ARMLOCAL_GTWOYC(double dMeanRevSpeed,
							double dSigma,
							long dZcId,
							long fZcId,
							double ratesCorr,
							double fMeanRevSpeed,
							double fSigma,
							ARM_result& result,
							long objId = -1);

extern long ARMLOCAL_GYCMODEL(long zcId,
							  double a,
							  double sigma,
							  ARM_result& result,
							  long objId = -1);

extern long ARMLOCAL_HWFNMONTECARLOSV(long zcId,
									  double horizon,
									  double a,
									  const VECTOR<double>& sigmaDate,
									  const VECTOR<double>& sigmaVal,
									  long nbTraj,
									  ARM_result& result,
									  long objId = -1);

extern long ARMLOCAL_HWFNMONTECARLO(long zcId,
									double horizon,
									double a,
									double sigma,
									long nbTraj,
									ARM_result& result,
									long objId = -1);

extern long ARMLOCAL_DFHWSIGVARTREE(double startDate,
									double horizon,
									long numSteps,
									long dZcId,
									long fZcId,
									double dMeanRevSpeed,
									const VECTOR<double>& dDate,
									const VECTOR<double>& dSigma,
									double fMeanRevSpeed,
									const VECTOR<double>& fDate,
									const VECTOR<double>& fSigma,
									long dFxCorrId,
									long fFxCorrId,
									long fxVolId,
									double ratesCorr,
									double fxSpotRate,
									ARM_result& result,
									long objId = -1);

extern long ARMLOCAL_DFHWBASISTREE (double startDate,
									double horizon,
									long numSteps,
									long dZcId,
									long fZcId,
									double dMeanRevSpeed,
									const VECTOR<double>& dDate,
									const VECTOR<double>& dSigma,
									double fMeanRevSpeed,
									const VECTOR<double>& fDate,
									const VECTOR<double>& fSigma,
									long dFxCorrId,
									long fFxCorrId,
									long fxVolId,
									double ratesCorr,
									double fxSpotRate,
									long dNonBasisZcId,
									long fNonBasisZcId,
									long SwaptionBSVolId,
									long SwaptionBSSmileId,
									long calageType,
									long pfType,
									double amin,
									double amax,
									double volmin,
									double volmax,
									ARM_result& result,
									long objId = -1);

extern long ARMLOCAL_GetParameter(long modId,
								  long paraId,
								  ARM_result& result);

extern long ARMLOCAL_FRMTREE(long FrmAnaId,
							 double horizon,
							 long fineMonth,
							 const VECTOR<double>& corrMatu,
							 const VECTOR<double>& corrMatrix,
							 const VECTOR<double>& corr2Matu,
							 ARM_result& result,
							 long objId = -1);

extern long ARMLOCAL_FRMTREE_AUTO(long zcId,
								  long volId,
								  long smileId,
								  long autoModeId,
								  double horizon,
								  long fineMonthId,
								  double shapeDecay,
								  double shapeSlope,
								  double shapeAsymptote,
								  long nbFactor,
								  const VECTOR<double>& corrMatu,
								  const VECTOR<double>& corrMatrix,
								  const VECTOR<double>& corr2Matu,
								  ARM_result& result,
								  long objId = -1);

extern long ARMLOCAL_FRMTREE_AUTO_G(long zcId,
									long zdId,
									long volId,
									long smileId,
									long irgvolId,
									long irgsmileId,
									long autoModeId,
									double horizon,
									long fineMonthId,
									double shapeDecay,
									double shapeSlope,
									double shapeAsymptote,
									long nbFactor,
									const VECTOR<double>& corrMatu,
									const VECTOR<double>& corrMatrix,
									const VECTOR<double>& corr2Matu,
									ARM_result& result,
									long objId = -1);

extern long ARMLOCAL_FRMANA(long zcId, 
							const VECTOR<CCString>& resetDates,
							const VECTOR<double>& spotVols,
							long shapeType,
							double shapeDecay,
							double shapeSlope,
							double shapeAsymptote,
							long nbFactor,
							const VECTOR<double>& corrMatu,
							const VECTOR<double>& corrMatrix,
							const VECTOR<double>& corr2Matu,
							ARM_result& result,
							long objId = -1);

extern long ARMLOCAL_FRMANA_PORT(long zcId,
								 long portId,
								 const VECTOR<CCString>& resetDates,
								 double precision,
								 double min_paras,
								 double max_paras,
								 long max_iters,
								 long shapeType,
								 double shapeDecay,
								 double shapeSlope,
								 double shapeAsymptote,
								 long nbFactor,
								 const VECTOR<double>& corrMatu,
								 const VECTOR<double>& corrMatrix,
								 const VECTOR<double>& corr2Matu,
								 ARM_result& result,
								 long objId = -1);

extern long ARMLOCAL_FRMLSMC(long FrmAnaId,
							 double horizon,
							 long nbTraj,
							 long fineMonth,
							 long mcMethod,
							 const VECTOR<double>& corrMatu,
							 const VECTOR<double>& corrMatrix,
							 const VECTOR<double>& corr2Matu,
							 ARM_result& result,
							 long objId = -1);

extern long ARMLOCAL_FRMLSMC_AUTO(long zcId,
								  long volId,
								  long smileId,
								  long intAutoMode,
								  double horizon,
								  long fineMonth,
								  long nbTraj,
								  long mcMethod,
								  long noControl,
								  double shapeDecay,
								  double shapeSlope,
								  double shapeAsymptote,
								  long nbFactor,
								  const VECTOR<double>& corrMatu,
								  const VECTOR<double>& corrMatrix,
								  const VECTOR<double>& corr2Matu,
								  ARM_result& result,
								  long objId = -1);

extern long ARMLOCAL_FRMLSMC_AUTO_2(long zcId, 
									long swoptVolId, 
									long swoptSmileId,
									long irgVolId, 
									long irgSmileId,
									long intAutoMode,
									double horizon,
									long fineMonth,
									long nbTraj,
									long mcMethod,
									long noControl,
									double shapeDecay,
									double shapeSlope,
									double shapeAsymptote,
									long nbFactor,
									const VECTOR<double>& corrMatu,
									const VECTOR<double>& corrMatrix,
									const VECTOR<double>& corr2Matu,
									ARM_result& result, long objId = -1);

extern long ARMLOCAL_FRMLSMC_AUTO_G(long zcId,
									long zdId,
									long swoptVolId, 
									long swoptSmileId,
									long irgVolId, 
									long irgSmileId,
									long intAutoMode,
									double horizon,
									long fineMonth,
									long nbTraj,
									long mcMethod,
									long noControl,
									double shapeDecay,
									double shapeSlope,
									double shapeAsymptote,
									long nbFactor,
									const VECTOR<double>& corrMatu,
									const VECTOR<double>& corrMatrix,
									const VECTOR<double>& corr2Matu,
									ARM_result& result, long objId = -1);

extern long ARMLOCAL_bsslmodel(double startDate,
							   long zc_type,
							   double zc,
							   long volSpreadLock_type,
							   double volSpreadLock,
							   long capVol_type,
							   double capVol,
							   long indexVol_type,
							   double indexVol,
							   ARM_result& result,
							   long objId = -1);

extern long ARMLOCAL_BKIRTREE (long zcId,
							   double startDate,
							   double endDate,
							   long pas,
							   double a,
							   double sigma,
							   ARM_result& result,
							   long objId = -1);

extern long ARMLOCAL_LDCMC_FROM_ANA(long anaModId,
									double horizon,
									long nbTraj,
									long mcmethod,
									long pricerTypeId,
									ARM_result& result,
									long objId = -1);

extern long ARMLOCAL_HWTREE(long zcId,
							double begDate,
							double endDate,
							long nbSteps,
							double a,
							double sigma,
							ARM_result& result,
							long objId = -1);

extern long ARMLOCAL_HWSIGCONST(long zcId,
								double begDate,
								double endDate,
								long nbSteps,
								double a,
								double sigma,
								ARM_result& result,
								long objId = -1);

extern long ARMLOCAL_HWSIGVAR(long zcId,
							  long begDate,
							  long endDate,
							  long nbSteps,
							  double a,
							  long real_size,
							  const VECTOR<CCString>& dates,
							  const VECTOR<double>& sigmas,
							  ARM_result& result,
							  long objId = -1);

extern long ARMLOCAL_DFGYC (double dMeanRevSpeed,
							double fMeanRevSpeed,
							double dSigma,
							double fSigma,
							double fxCorr,
							double fxVol,
							double ratesCorr,
							long dZcId,
							long fZcId,
							ARM_result& result,
							long objId = -1);

extern long ARMLOCAL_IDTREEHW(double startDate,
							  double horizon,
							  long nbSteps,
							  double dMeanRevSpeed,
							  double fMeanRevSpeed,
							  double dSigma,
							  double fSigma,
							  double prtyCorr,
							  double prtyVol,
							  double ratesCorr,
							  long dZcId,
							  long fZcId,
							  ARM_result& result,
							  long objId = -1);

extern long ARMLOCAL_HWANALYTICSIGVAR(long zcId,
									  double a,
									  long real_size,
									  const VECTOR<CCString>& dates,
									  const VECTOR<double>& sigmas,
									  ARM_result& result,
									  long objId = -1);

extern long ARMLOCAL_BASISMCFRM(long zcId,
								long baZcId,
								long volId,
								long smileId,
								long productType,
								double horizon,
								long nbTraj,
								long MCGeneratorType,
								double shapeDecay,
								double shapeSlope,
								double shapeAsymptote,
								long nbFactor,
								const VECTOR<double>& indexes,
								const VECTOR<double>& correlatedIndex,
								const VECTOR<double>& corrMatrix,
								long control,
								long seed,
								ARM_result& result,
								long objId = -1);

extern long ARMLOCAL_BASISMCFRM2CR(long zcId,
								   long baZcId,
								   long swoptvolId,
								   long swoptsmileId,
								   long irgvolId,
								   long irgsmileId,
								   long productType,
								   double horizon,
								   long nbTraj,
								   long MCGeneratorType,
								   double shapeDecay,
								   double shapeSlope,
								   double shapeAsymptote,
								   long nbFactor,
								   const VECTOR<double>& indexes,
								   const VECTOR<double>& correlatedIndex,
								   const VECTOR<double>& corrMatrix,
								   long control,
								   long seed,
								   ARM_result& result,
								   long objId = -1);

extern long ARMLOCAL_GetCalibrationOutputPFSize(long modId, ARM_result& result);
extern long ARMLOCAL_GetCalibrationOutputHWVolSize(long modId, ARM_result& result);
extern long ARMLOCAL_GetCalibrationOutputMeanRev(long modId, ARM_result& result);
extern long ARMLOCAL_GetCalibrationOutputDate(long modId, long i, ARM_result& result);
extern long ARMLOCAL_GetCalibrationVolSched(long modId, long i, ARM_result& result);
extern long ARMLOCAL_GetCalibrationOutputSecExeDate(long modId, long i, ARM_result& result);
extern long ARMLOCAL_GetCalibrationOutputSecStartDate(long modId, long i, ARM_result& result);
extern long ARMLOCAL_GetCalibrationOutputSecEndDate(long modId, long i, ARM_result& result);
extern long ARMLOCAL_GetCalibrationOutputSecStrike(long modId, long i, ARM_result& result);
extern long ARMLOCAL_GetCalibrationOutputInputVol(long modId, long i, ARM_result& result);
extern long ARMLOCAL_GetCalibrationOutputOutputVol(long modId, long i, ARM_result& result);
extern long ARMLOCAL_GetCalibrationOutputErrorVol(long modId, long i, ARM_result& result);
extern long ARMLOCAL_GetCalibrationOutputInputPrice(long modId, long i, ARM_result& result);
extern long ARMLOCAL_GetCalibrationOutputOutputPrice(long modId, long i, ARM_result& result);
extern long ARMLOCAL_GetCalibrationOutputErrorPrice(long modId, long i, ARM_result& result);

extern long ARMLOCAL_SMILEDLDCANA(long anaModId,
								  double dVolSpotUp,
								  double Asymp_RowUp,
								  double Asymp_ColUp,
								  double pUp,
								  double dVolSpotDo,
								  double Asymp_RowDo,
								  double Asymp_ColDo,
								  double pDo,
								  ARM_result& result,
								  long objId = -1);

extern long ARMLOCAL_SMILEDLDCFROMANA(long anaModId,
									  double horizon,
									  long nbTraj,
									  long mcmethod,
									  ARM_result& result,
									  long objId = -1);

extern long ARMLOCAL_bssmiledmodel(double date,
								   double spot,
								   long dividend_type,
								   double dividend,
								   long discrate_type,
								   double discrate,
								   long volat_type,
								   double volat,
								   long typstk,
								   const VECTOR<double>& matu,
                                   long rhoType,
                                   long rhoObj,
								   const VECTOR<double>& rho,
                                   long nuType,
                                   long nuObj,
								   const VECTOR<double>& nu,
								   long isSABR,
								   long beta_type,
								   double beta,
								   const VECTOR<double>& beta_vect,
								   double weight,
                                   long sigmaOrAlphaFlag,
								   long correlManagerId,
								   long rdiscrate_type,
								   double rdiscrate,
								   long adjCvxVolId,
								   long convToAdjVolWhenAlpha,
                                   ARM_result& result,
								   long objId = -1);

extern long ARMLOCAL_bssmiledmodel_SABR(double date,
										double spot,
										long dividend_type,
										double dividend,
										long discrate_type,
										double discrate,
										long SABRVolId,
										long typstk,
										long correlManagerId,
										long rdiscrate_type,
										double rdiscrate,
										long adjCvxVolId,
										long convToAdjVolWhenAlpha,
										ARM_result& result,
										long objId = -1);

extern long ARMLOCAL_CRRTREE(double begDate,
							 double endDate,
							 long nbSteps,
							 double spot,
							 double dividend,
							 double discrate,
							 double volat,
							 ARM_result& result,
							 long objId = -1);

extern long ARMLOCAL_BSCORRMODEL(double startDate,
								 long zcId,
								 long spreadlockId,
								 long capirgvolId,
								 long capcashvolId,
								 long indexadjvolId,
								 long spreadvolId,
								 long correlationsId,
								 long modelTypeId,
								 long volTypeId,
								 ARM_result& result,
								 long objId = -1);

long ARMLOCAL_FRMTREE_AUTO_B(long zcId,
							 long zdId, 
							 long volId, 
							 long smileId,
							 long autoModeId,
							 double horizon,
							 long fineMonthId,
							 double shapeDecay,
							 double shapeSlope,
							 double shapeAsymptote,
							 long nbFactor,
							 const VECTOR<double>& corrMatu,
							 const VECTOR<double>& corrMatrix,
							 const VECTOR<double>& corr2Matu,
							 ARM_result& result,
							 long objId = -1);

extern long ARMLOCAL_DFFXBS (long dVolId,
							 long fVolId,
							 long dZcId,
							 long fZcId,
							 long dFxCorrId,
							 long fFxCorrId,
							 long fxVolId,
							 double ratesCorr,
							 long dBSZcId,
							 long fBSZcId,
							 double spot,
                             long FundZcId,
                             long FundBSZcId,
                             double fundFxSpot,
							 long fxSmileId,
                             bool isLnVols,
                             long RhoDomId,
                             long NuDomId,
                             long BetaDomId,
                             long RhoForId,
                             long NuForId,
                             long BetaForId,
							 long fxVolModelId,
							 ARM_result& result,
							 long objId = -1);

extern long ARMLOCAL_DFFXBS (long dBSId,
							 long fBSId,
							 long dFxCorrId,
							 long fFxCorrId,
							 long fxVolId,
							 double ratesCorr,
							 long dBSZcId,
							 long fBSZcId,
							 double spot,
                             long FundZcId,
                             long FundBSZcId,
                             double fundFxSpot,
							 long fxSmileId,
							 long fxVolModelId,
							 ARM_result& result,
							 long objId = -1);

extern long ARMLOCAL_GetComputedSigmaSABR(long modId,
										  ARM_result& result);

extern long ARMLOCAL_CALIBRATIONHWSV(double asof,
									 long zcId,
									 long volId,
									 long secId,
									 double amin,
									 double amax,
									 double volmin,
									 double volmax,
									 const VECTOR<double>& dates,
									 long pfId,
									 long RhoSwoptId,
									 long NuSwoptId,
									 long BetaSwoptId,
									 ARM_result& result,
									 long objId = -1);

extern long ARMLOCAL_CALIBRATE (long calibId,
								long calibType,
								ARM_result& result,
								long objId = -1);

extern long ARMLOCAL_HWSIGVARFROMANA(long modId,
									 long endDate,
									 long nbSteps,
									 ARM_result& result,
									 long objId = -1);

extern long ARMLOCAL_BSSMILEDCALIBRATE(long modelId,
									   long pfId,
									   long calVolOrNotId,
									   long calRhoOrNotId,
									   long calNuOrNotId,
									   long calBetaOrNotId,
                                       double volTenor,
									   double timeStep,
									   double minSig,
									   double maxSig,
									   double minRho,
									   double maxRho,
									   double minNu,
									   double maxNu,
                                       double minBeta,
                                       double maxBeta,
									   long interpMethodId,
									   double tol,
									   long maxIter,
									   long gradCalcId,
									   double lambda,
									   long globOrBootstrapId,
									   double beginSmoothMatu,
									   ARM_result& result,
									   long objId = -1);

extern long ARMLOCAL_GETCALIBRATED_SIGRHONU(long modelId,
											long paramId,
											ARM_result& result,
											long objId = -1);

extern long ARMLOCAL_CalibrationFRMModel(long zcId,
										 long pf1Id,
                                         long pf2Id,
										 long nbFactors,
										 long liborTypeId,
										 const VECTOR<double>& powers,
										 const VECTOR<double>& smoothParams,
										 const VECTOR<double>& ACalibrationSchedule,
										 const VECTOR<double>& KCalibrationSchedule,
										 long initA,
										 long initK,
										 const VECTOR<double>& meanRevA,
										 const VECTOR<double>& meanRevK,
                                         const VECTOR<double>& SchedPrice,
                                         const VECTOR<double>& mdec,
										 const VECTOR<double>& Bounds,
										 const VECTOR<double>& optimizerParams,
                                         long calibrate,
										 ARM_result& result,
										 long objId = -1);

extern long ARMLOCAL_GETFROMFRMMODEL(long modelId,
									 const CCString& param,
									 ARM_result& result);


extern long ARMLOCAL_BootstCalibFRMModel( long zcId,
										  long pf1Id,
										  long pf2Id,
                                          long pf3Id,
                                          long marketmodelId,
										  long liborTypeId,
										  const VECTOR<double>& CalibParams,
										  const VECTOR<double>& initCurve,
										  const VECTOR<double>& mdec,
										  double meanRev,
                                          long   nbfactor,
										  long   nbrows,
										  long   nbcolumns,
                                          const VECTOR<double>& correlmatrix,
                                          long  VolType,
                                          long calibrate,
                                          long PreInitialise,
                                          double Presicion,
                                          double VegaLevel,
										  ARM_result& result,
										  long objId = -1);

extern long ARMLOCAL_BootstCalibFRMModelMixture(
										int N,
										long zcId,
										long pf1Id,
										long pf2Id,
										int nbProduct,
										long liborTypeId,
										const VECTOR<double>& lambdas,
										double meanRev,
										long   nbfactor,
										long   nbrows,
										long nbcolumns,
										const VECTOR<double>& correlmatrix,
										long VolType,
										const VECTOR<double>& lowerBound,
										const VECTOR<double>& upperBound,
										long calibrate,
										ARM_result& result,
										long objId = -1);

extern long ARMLOCAL_GETFROMFRMMODELMIXTURE(long modelId,
											int modelNumber,
											const CCString& param,
											ARM_result& result);

extern long ARMLOCAL_MCFRMMODEL( long FRMModId,
                    	         long nbTraj,
						         long nbStepIn,
                                 long pricerTypeId,
                                 long mcMethod,
						         ARM_result& result,
						         long objId = -1);

extern long ARMLOCAL_MCFRMMODELMIXTURE( long FRMModMixId,
                    	         long nbTraj,
						         long nbStepIn,
                                 long pricerTypeId,
						         ARM_result& result,
						         long objId = -1);

extern long ARMLOCAL_FRMMARKOVTREE(double startDate,
                                   double horizon,
								   long zcId,
                                   long PathNumber,
                                   const VECTOR<double>& Params,
                                   long FRMModelId,
								   ARM_result& result,
								   long objId = -1);

extern long ARMLOCAL_LOGDECANA(long zcId,
							   long frequency,
							   const VECTOR<double>& resetMatu,
							   const VECTOR<double>& shifts,
							   const VECTOR<double>& fwdVols,
							   long isWeighted,
							   ARM_result& result,
							   long objId = -1);

extern long ARMLOCAL_QMODEL(double date,
							double spot,
							long dividend_type,
							double dividend,
							long discrate_type,
							double discrate,
							long volat_type,
							double volat,
							double q0,
							double q1,
							long typstk,
							ARM_result& result,
							long objId = -1);

extern long ARMLOCAL_CROSSMODEL(double date,
								long DbsmodId,
								long FbsmodId,
								long fxvolId,
								long dfxcorrelId,
								long ffxcorrelId,
								long dfcorrelId,
								long discountZcId,
								double correlforAdj,
								double adjustFlag,
								double slopeFlag,
								ARM_result& result,
								long objId = -1);


extern long ARMLOCAL_GLOBDFBS(long DomBSId,
							  long DomCurrId,
							  long FrgBSId,
							  long FrgCurrId,
							  long fxVolCrvId,
							  long FFxCorrId,
							  long RatesCorrId,
							  long fxVolModelId,
							  ARM_result& result,
							  long objId = -1);

extern long ARMLOCAL_BSConvAdjust_Create(int SUMMITFormulaeUsed, int UseSabrCMS, ARM_result& result, long objId = -1);

extern long ARMLOCAL_ReplicConvAdjust_Create(int Payoff_ReplicMode,
									  double Payoff_ReplicPrecision,
									  int Payoff_StopMode,
									  double Payoff_StopThreshold,
									  int Sensi_ReplicMode,
									  double Sensi_ReplicPrecision,
									  int Sensi_StopMode,
									  double Sensi_StopThreshold,
									  long UsedModelId,
									  double StrikeMinReplic,
									  double StrikeMaxReplic,
									  ARM_result& result, 
									  long objId = -1);

extern long ARMLOCAL_BSConvAdjustRep_Create(long UsedModelId, long swoptVolCurveId, 
											const VECTOR<double>& stddev, 
											int NbPtsForRepliq,
											long MR,
											bool FullRepliq, 
											double upperProba,
											double lowerProba,
											ARM_result& result,
											long objId = -1);

extern long ARMLOCAL_MapConvAdjust_Create(long LiborArrearAdjId,
								   long NaturalCMSAdjId,
								   long PaymentLagAdjId,
								   ARM_result& result,
								   long objId);


long ARMLOCAL_ReplicModel_Create(
							int ReplicMode,
							double StepOrReplicPrecision,
							int StopMode,
							double StopThreshold,
                            int SensiReplicMode,
							double SensiStepOrReplicPrecision,
                            int SensiStopMode,
							double SensiStopThreshold,
							long UsedModelId,
                            int SUMMITFormulaeUsed,
							double StrikeMinReplic,
							double StrikeMaxReplic,
							ARM_result& result,
							long objId);

long ARMLOCAL_SetReplicDebugMode(
                    int DebugMode,
                    ARM_result& result);

extern long ARMLOCAL_TriBSModel (long Model1Id,
								 long Model2Id,
								 long DiscModelId,
								 long Fx1DiscVolId,
								 long Fx2DiscVolId,
								 long Idx1Idx2CorrId,
								 long Idx1DiscIdxCorrId,
								 long Idx2DiscIdxCorrId,
								 long Idx1FxCorrId,
                                 long Idx2FxCorrId,
								 int quantoadjflag,
								 ARM_result& result,
								 long objId = -1);

extern long ARMLOCAL_SetDiscountPricingMode(long modelId,
											long discPriceMode,
											ARM_result& result);

extern long ARMLOCAL_CALIBRATORSFRM(long secId,
									long capModelId,
									long swoptModelId,
									double meanRev,
									const VECTOR<double>& calibParams,
									long preInitFlagId,
									const VECTOR<double>& initSigmaCrv,
									const VECTOR<double>& initBetaOrShift,
									long sigmaPfId,
									long betaPfId,
									long meanRevPfId,
									long voltypeId,
									long correlId,
                                    const VECTOR<double>& SecurityParams,
                                    const CCString tocalswaptATM,
									ARM_result& result,
									long objId = -1);

extern long ARMLOCAL_SFRMCALIBRATE (long calibratorSFRMId,
                                    long mktcapmodelId,
                                    long mktswaptmodelId,
                                    CCString tocalibrateBeta,
                                    CCString tocalibrateMR,
                                    long kerneltoGP,
                                    ARM_result& result,
                                    long objId = -1);

extern long ARMLOCAL_TriBSDualModel (long Model1Id,
								     long Model2Id,
								     long DiscModelId,
								     long Fx1DiscVolId,
								     long Fx2DiscVolId,
								     long Idx1Idx2CorrId,
								     long Idx1DiscIdxCorrId,
								     long Idx2DiscIdxCorrId,
								     long Idx1FxCorrId,
                                     long Idx2FxCorrId,
								     int quantoadjflag,
                                     double correlforadj,
                                     int withslopeflag,
								     ARM_result& result,
								     long objId = -1);

extern long ARMLOCAL_CRACALCULATOR(long secId,
								   long zcId,
								   long capVolATMId,
								   long rhoId,
								   long nuId,
								   long swptVolId,
								   long betaSABRId,
								   double meanRev,
								   const VECTOR<double>& calibParams,
								   long preInitFlagId,
								   const VECTOR<double>& initSigmaCrv,
								   const VECTOR<double>& initBetaOrShift,
								   long voltypeId,
                                   const VECTOR<double>& SecurityParams,
								   double horizon,
								   long PathNumber,
								   const VECTOR<double>& TreeParams,
								   const VECTOR<double>& ReCalibFlags,
								   CCString CalswaptATM,
								   long rhoSwoptId,
								   long nuSwoptId,
								   long betaSwoptId,
                                   long SABRSigmaOrAlphaFlag,
								   ARM_result& result,
								   long objId = -1);

extern long ARMLOCAL_BUMPCRACALCULATOR(long cracalcId,
									   long zcId,
									   long capVolATMId,
									   long rhoId,
									   long nuId,
									   long swptVolId,
									   long betaSABRId,
									   long rhoSwoptId,
									   long nuSwoptId,
									   long betaSwoptId,
                                       long SABRSigmaOrAlphaFlag,
									   ARM_result& result,
									   long objId = -1);

extern long ARMLOCAL_COMPUTEFWDBSDELTA(	double fwd, 
									    double strike, 
										double vol,
										double T, 
										int CallPut, 
										ARM_result& result);

extern long ARMLOCAL_ComputeSplinedSigmaATMF(VECTOR<double>& deltas, 
											 VECTOR<double>& sigmas, 
											 double matu, 
											 double SigmaZDS, 
											 double Precision,
                                             double FX_SPOT,
											 ARM_result& result);

extern long ARMLOCAL_ComputeDeltaFwdFromDeltaWP(double AsOf, 									
												double matu,
												double sigma,
												double fxSpot,
												double deltaWithPremium,
												long domCrvId,
												long foreignCrvId,
												ARM_result& result);

extern long ARMLOCAL_CalculateImpliedStrikeFromDeltaWithPremium (double AsOf, 									
																 double matu,
																 double sigma,
																 double fxSpot,
																 double deltaWithPremium,
																 long domCrvId,
																 long foreignCrvId,
																 ARM_result& result);


extern long ARMLOCAL_CalcFwdFXSpot(double AsofDate,
								   double Spot,
								   double aFwdDate,
								   long NumDiscountCurveID,
								   long UndDiscountCurveID,
								   ARM_result& result);


extern long ARMLOCAL_GetQModelQVolatility(long QModelId,ARM_result& result);

extern long ARMLOCAL_FromSigmaToAlpha( long bssimilemodId,
									   long sigmaCurveId,
									   double strike,
									   ARM_result& result,
									   long objId = -1);

extern long ARMLOCAL_FromAlphaToSigma( long bssimilemodId,
									   long alphaCurveId,
									   double strike,
									   ARM_result& result,
									   long objId = -1);

extern long ARMLOCAL_TriXBSModel (	long Model1Id,
									long Model2Id,
									long Model3Id,
									long Idx1Idx2CorrCurveId,
									long Idx1Idx3CorrCurveId,
									long Idx2Idx3CorrCurveId,
									ARM_result& result,
									long objId = -1);
									
extern long ARMLOCAL_GetAsOf(long modelId,
							 ARM_result& result);

#endif	// ARM_LOCAL_MOD_H

// EOF %M%
