#ifndef ICM_LOCAL_STR_H
#define ICM_LOCAL_STR_H

#include <vector>
#include <ARM\libarm\ARM_result.h>
#include "ARMKernel\glob\linalg.h"
#include "CCstring.h"
// using namespace std;

extern long ICMLOCAL_STR_Cppi (vector<double> CrbDefaut,
							   long nbRowCrbDef,
							   long nbColCrbDef,
							   double Recov,
							   vector<double> correl,
							   long nbRowCorrel,
							   long nbColCorrel,
							   vector<double> DataCppi,
							   vector<double> riskyFactor,
							   long nbRowRiskyFactor,
							   long nbColRiskyFactor,
							   vector<double> DataPlacement,
							   vector<double> DataSoult,
							   vector<double> DataVasicek,
							   vector<double> DataModel,
							   vector<double> DataCoupon,
							   vector<double> DataWidener,
							   vector<double> DataSimul,
							   vector<double> DataBarriere,
							   vector<double> taux,
							   long nbRowtaux,
							   long nbColtaux,
							   CCString Fonction,
							   vector<double>& Resultat,
							   ARM_result& result);

/*extern long ICMLOCAL_STR_GenAlea (int Seed,
								  int NbTir,
								  vector<double>& Resultat,
								  ARM_result& result);
*/

// CDS
extern long ICMLOCAL_STR_CDS (double SpreadCDS,
							  double Maturity,
							  double PasPaiement,
							  vector<double> VRecov,
							  double Nominal,
							  vector<double> SpreadCourbe,	
							  vector<double> PlotCourbe,
							  double taux,
							  double DateShift,
							  double ValShift,
							  CCString C_Stripping,
							  CCString C_Function,
							  vector<double>& Resultat,
							  ARM_result& result);
// TRANCHES
extern long ICMLOCAL_STR_SpreadForwardCurves(vector<double> SpreadInit,
											 vector<double> Plot,
											 double Recov,
											 double Taux,
											 double TimeInit,
											 double Duree,
											 CCString Stripping,
											 ARM_result& result);

// Discount Factor --------------------------------------
/*extern long ICMLOCAL_STR_DPrice ( double  YF_Maturity,
								  double  Taux,
								  ARM_result& result);*/
// Proba défaut ------------------------------------------

/*extern long ICMLOCAL_STR_PDefaultCurves (double YF_Maturity,
										 double IssuersRecovery,
										 vector<double> Spread,
										 vector<double> YF_plot,
										 double taux,
										 CCString Stripping,
										 ARM_result& result);
*/
/*
extern long ICMLOCAL_STR_ProbDefTranche (double YF_Maturity,
										   vector<double> Attachement,
										   double dNb_Issuers,
										   vector<double> Recovery,
										   vector<double> MatIssuersSpread,
										   vector<double> YF_Plot,		  
										   vector<double> VBeta,
										   double taux,
										   CCString Stripping,
										   ARM_result& result);
										   */
// Densité loi jointes -------------------------------------
/*
extern long ICMLOCAL_STR_DensityCurves (double YF_Maturity,
										double dNb_Issuers,
										double IssuersRecovery,
		  								vector<double> VIssuersSpread,
										vector<double> YF_Plot,
										char** Label,
										vector<double> VBetas,
										double taux,
										vector<double>& Density,
										ARM_result& result);
										*/
// CDO2 SMILE (code ARM)
extern long ICMLOCAL_STR_CdO2Smile(vector<double> DataMez,		
								   vector<double> NotionelTranche,
								   vector<CCString>& StrikeCdo,
								   int nbColStrikeCdo,
								   vector<CCString> DataDefaultCurve,
								   int nbColDefaultCurve,
								   int nbRowDefaultCurve,
								   vector<CCString> IssuersLabel,
								   vector<double> SpreadTraxx,
								   vector<CCString>& MatriceOverlap,
								   int nbColMatriceOverlap,
								   int nbRowMatriceOverlap,
								   vector<double> DataTraxx,
								   vector<CCString> Recov,
								   int nbColRecov,
								   int nbRowRecov,
								   vector<CCString> ForcedCorrel,
								   int nbColForcedCorrel,
								   int nbRowForcedCorrel,
								   vector<CCString>& CorrelParStrike,
								   int nbRowCorrel,
								   int nbColCorrel,
								   CCString LabelToShift,
								   double DateToShift,
								   double ValShift,
								   vector<CCString> NbStepNumMethode,
								   CCString C_Function,
								   ARM_result& result);

// -----------------------------------------------------------
// INTERFACAGE PROTO - CODE DAMIEN
// -----------------------------------------------------------
// ----------------------------------------------------------- 
// Paniers
extern long ICMLOCAL_STR_Basket(double Maturity,
							    vector<double> IndiceDefault,
								double Nb_Issuers,
								double IssuersNotionals,
								vector<double> CalibRecovery,
								vector<double> LossRecovery,
								vector<double> YF_plot,
								vector<double> MatIssuersSpread,								
								vector<CCString> IssuersLabel,
								double FrequencyFeeLeg,
								double FixedRate,
								vector<double> IssuersBeta,
								double taux,
								CCString LabelToShift,
								double DateToShift,
								double ValShift,
								CCString C_Stripping,
								CCString C_Function,
								ARM_result& result);


// CDO
extern long ICMLOCAL_STR_CDO(int MezTtpe,
							 vector<double> DataMez,
							 vector<double> Notional,
							 vector<double> CalibRecovery,
							 vector<double> LossRecovery,
							 vector<double> YF_plot,
							 vector<double> MatIssuersSpread,								
							 vector<CCString> IssuersLabel,
							 vector<double> IssuersBeta,
							 double taux,
							 CCString LabelToShift,
							 double DateToShift,
							 double ValShift,
							 vector<double> SpreadCrbFwd,
							 double VolAjustConvexity,
							 vector<CCString> NumMethod,
							 CCString C_Function,
							 ARM_result& result);

/**
extern long ICMLOCAL_STR_CDOSMILE ( int MezType, 
								    vector<double> DataMez,
								    vector<double> Recovery,
									vector<double> PlotPtf,
									vector<double> MatIssuersSpread,								
									vector<CCString> IssuersLabel,
									vector<double> PlotTraxx,
									vector<double> SpreadTraxx,
									vector<double> MaturityBaseCorrel,
									vector<CCString>& CorrelParStrike,
									int nbRowCorrel,
									int nbColCorrel,
									vector<double> DataTraxx,
									double taux,
									CCString LabelToShift,
									double DateToShift,
									double ValShift,
									vector<double> SpreadCrbFwd,
									double VolAjustConvexity,
									vector<CCString> NumMethod,
									CCString C_Function,
									vector<double>& Resultat,
									ARM_result& result);
									**/ 


// Expected Loss
extern long ICMLOCAL_STR_GetExpectedLoss(double YearTerm,
										 vector<double> Attachement,
										 double Nb_Issuers,
										 vector<double> CalibRecovery,
										 vector<double> LossRecovery,
										 vector<double> YF_plot,
										 vector<double> MatIssuersSpread,								
										 vector<double> IssuersBeta,
										 double taux,
										 vector<CCString> NumMethod,
										 ARM_result& result);

extern long ICMLOCAL_STR_ExpectedLossSMILE( vector<double> DataMez,
								 		    vector<double> CalibRecovery,
										    vector<double> LossRecovery,
											vector<double> YF_plot,
											vector<double> MatIssuersSpread,								
											vector<double> SpreadTraxx,
											vector<double> MaturityCorrel,
											vector<CCString>& CorrelParStrike,
										    int nbRowCorrel,
										    int nbColCorrel,
											vector<double> DataTraxx,
								    		vector<CCString> NumMethod,
								    		vector<double>& Resultat,
											ARM_result& result);

// Interpolatino lineaire (utile pour retrouver la correl interpolée : X = Strike et Y = Correl) 
/** 
extern long ICMLOCAL_STR_LinearInterpol(vector<double> XRef,
										vector<double> YRef,
										vector<double> XTarget,
										vector<double>& YTarget,
										ARM_result& result);
										**/ 
						
/*
long ICMLOCAL_Str_CorrelInterTranche(long NbSimul,
											int NbTranche,
											int NbName,
											double BetaIssuers,
											vector<double> CalibRecovery,
											double LossRecovery,
											double Maturite,
											vector<CCString>& StrikeDownUp,
											int nbColStrikeDownUp,
											vector<CCString>& Appartenance,
											int nbColIssuersInCDO,
											vector<double> YF_plot,
											vector<double> MatIssuersSpread,
											double taux, 
											ARM_Vector & ProbaDefautTranche,
											ARM_Vector & VarianceTranche,
											vector<double*>& MatriceEvtDefaut,
											vector<double*>& MatriceTpsDefaut,
											ARM_result& result);

*/
/*
extern long ICMLOCAL_Str_PriceAndSensiForCDSOptions(CCString PricingType,
									  double Spread ,
									  double Spread5Y, 
									  int Frequency, 
									  double NbIssuers , 
									  double Notional , 
									  double TradeDate , 
									  double OptionMaturity , 
									  double CDSMaturity , 
									  double Strike ,
									  int OptionType, 
									  double Vol, 
									  int KOType, 
									  CCString UnderlyingTp,  
									  double Recovery,
									  double Rate, 
									  VECTOR<CCString> Maturities, 
									  VECTOR<double> Spreads,
									  double MktPrice,
									  VECTOR<double>& Res, 
									  ARM_result& result);


*/

#endif	// ICMLOCAL_STR_TEST