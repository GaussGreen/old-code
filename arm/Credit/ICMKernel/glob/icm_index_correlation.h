/********************************************************************************/
/*! \file icm_index_correlation.h
 * 
 *  \brief Describes an Index Correlation Matrix Object
 *  \author 
 *	\version 1.0
 *	\date   June 2005 */
/*
 *********************************************************************************/

#ifndef _ICM_INDEX_CORRELATION_H_
#define _ICM_INDEX_CORRELATION_H_

/*********************************************************************************/
/*! \class  ICM_Index_Correlation icm_index_correlation.h "icm_index_correlation.h"
 *  \author D Pouponneau
 *	\version 1.0
 *	\date   April 2004
 *	\brief Creates an Index correlation object */
/***********************************************************************************/

#include "ICMKernel/glob/icm_correlation.h"
//#include "ICMKernel/inst/icm_credit_index.h"
//#include "ARMKernel/util/interpol.h"
//#include "ICMKernel/util/icm_functors.h"
//#include "ICMKernel/pricer/icm_pricer_cds.h"
//#include "ICMKernel/mod/modelmulticurves.h"
//#include "ICMKernel/glob/icm_smile_correlation.h"

double PSW_objfun_RF_Allfactors_AllQuotes(int n, double *y, double *lb, double *ub);

class ICM_Parameters;
class ICM_Mez;
class ICM_Smile_Correlation; 
class ICM_Credit_Index; 
class ICM_ModelMultiCurves; 

class ICM_Index_Correlation : public ICM_Correlation
{
public :
	ARM_Vector* itsStrikeLow;			//vector of Strikes low (in %)
	ARM_Vector* itsStrikeHigh;			//vector of Strikes high (in %)

	//Internal Use Only
	ICM_Pricer* its_cal_Pricer;
	ICM_Pricer* its_cal_PrevPricer;

	double		its_cal_CurrentPrice;
	double		its_cal_K1;
	double		its_cal_K2;
	double		its_cal_CorrelDown;
	std::vector<ARM_Security*> its_cal_securitites;
	std::vector<double> its_cal_CurrentPrices;
	ICM_Parameters* its_cal_Param;
	ARM_Model*	its_cal_model;

	double		its_cal_valueparam;
	double		its_cal_fixedparam;
	double		its_cal_thres;
	double		its_cal_Hedge;
	double		its_cal_cdo_Notional;
	double		its_RFL_Beta0;

	std::vector< std::vector<double> >	itsImpliedParameters;	//for correlation parametrisation


private:

	qCAL_INDEX_CORR_TYPE itsCalMethod;	//Calibration method Base Correlation...
	ICM_Credit_Index*    itsCreditIndex;//Credit Index

	vector<ICM_Credit_Index*>    itsPrevCreditIndex;// Vector of Credit Index for term structure BC
	ICM_QMatrix<double>			 itsPrevBaseCorrels;// Vector of previous BC for term structure BC

	ARM_Vector* itsMktBid;				//Market Bid (in %)
	ARM_Vector* itsMktAsk;				//Market Ask (in %)
	ARM_Vector* itsUpfBid;				//Market Bid (in %)
	ARM_Vector* itsUpfAsk;				//Market Ask (in %)
	ARM_Vector* itsInitialCorrelation;  //Initial Correlation for optimisation
	ARM_Vector* itsImpliedCorrelation;  //for base correlations
	ARM_Vector* itsLeverages;			//Leverages for delta hedges
	ICM_ModelMultiCurves* itsModelMultiCurves; //modelmulticurves for dispertion
	bool		itsIsLinearInterpol;	//interpolation type 

	int		itsIntegrationStep;	//integration step for calibration
	int		itsLagForStartDate;	//Lag for start dates
	int		itsCreditLag;		//credit lag

	ARM_Vector	its_cal_yts;
	ARM_Vector	its_cal_strikes;
	ARM_Matrix	its_cal_bc;

	double its_ts_prev_term;
	double its_ts_last_EL_prevMaturity;

	double itsStep;

	qOPTIMIZE_TYPE itsCalMeth;

	ICM_Smile_Correlation* itsTemplateBaseCorrel;
	ICM_Parameters*	itsInParams;

	void Init(void) ;

public: 

	ICM_Index_Correlation() { Init();}

	ICM_Index_Correlation(const ARM_Date& AsOf,
			 const string& name,
			 const qCAL_INDEX_CORR_TYPE& CalMethod,
			 ICM_Credit_Index* CreditIndex,
			 ARM_Vector* StrikeLow,
			 ARM_Vector* StrikeHigh,
			 ARM_Vector* MktBid,
			 ARM_Vector* MktAsk,
			 ARM_Vector* UpfBid,
			 ARM_Vector* UpfAsk,
			 ARM_Vector* Leverage,
			 ARM_Vector* InitialCorrelation = NULL,
			 ICM_ModelMultiCurves* mmc = NULL,
			 const int IntegrationStep = 60,
			 const int LagStartDate = 1,
			 const int CreditLag = DEFAULT_CREDIT_LAG_INDX,
			 vector<ICM_Credit_Index*> PrevCreditIndex = vector<ICM_Credit_Index*>(),
			 ICM_QMatrix<double> PrevBCMatrix = ICM_QMatrix<double>(),
			 const double step = 0.01,
			 const qOPTIMIZE_TYPE CalMeth = qNEWTON) 
	{
		Init();
		Set(AsOf,name,CalMethod,CreditIndex,StrikeLow,StrikeHigh,MktBid,MktAsk,UpfBid,UpfAsk,Leverage,InitialCorrelation,mmc,IntegrationStep,LagStartDate,CreditLag,PrevCreditIndex,PrevBCMatrix,step,CalMeth);
	}	

	void Set(const ARM_Date& AsOf,
			 const string& name,
			 const qCAL_INDEX_CORR_TYPE& CalMethod,
			 ICM_Credit_Index* CreditIndex,
			 ARM_Vector* StrikeLow,
			 ARM_Vector* StrikeHigh,
			 ARM_Vector* MktBid,
			 ARM_Vector* MktAsk,
			 ARM_Vector* UpfBid,
			 ARM_Vector* UpfAsk,
			 ARM_Vector* Leverage,
			 ARM_Vector* InitialCorrelation = NULL,
			 ICM_ModelMultiCurves* mmc = NULL,
			 const int IntegrationStep = 60,
			 const int LagStartDate = 1,
			 const int CreditLag = DEFAULT_CREDIT_LAG_INDX,
			 vector<ICM_Credit_Index*> PrevCreditIndex = vector<ICM_Credit_Index*>(),
			 ICM_QMatrix<double> PrevBCMatrix = ICM_QMatrix<double>(),
			 const double step = 0.01,
			 const qOPTIMIZE_TYPE CalMeth = qNEWTON) ; 
 
	void GlobalCalibrate() ;


	void SetCalMethod(const qCAL_INDEX_CORR_TYPE& CalMethod) {itsCalMethod = CalMethod;}
	qCAL_INDEX_CORR_TYPE GetCalMethod() { return itsCalMethod;}

	void SetCreditIndex(ICM_Credit_Index* index) ;
	ICM_Credit_Index* GetIndex(void) { return itsCreditIndex;}

	void SetUpfBid(ARM_Vector* UpfBid) 
	{
		if (itsUpfBid)	delete itsUpfBid;
		if (UpfBid) itsUpfBid = (ARM_Vector*) UpfBid->Clone();
	}
	ARM_Vector* GetUpfBid(void) { return itsUpfBid;}

	void SetUpfAsk(ARM_Vector* UpfAsk) 
	{
		if (itsUpfAsk)	delete itsUpfAsk;
		if (UpfAsk) itsUpfAsk = (ARM_Vector*) UpfAsk->Clone();
	}
	ARM_Vector* GetUpfAsk(void) { return itsUpfAsk;}

	void SetStrikeLow(ARM_Vector* StrikeLow) 
	{
		if (itsStrikeLow)	delete itsStrikeLow;
		if (StrikeLow) itsStrikeLow = (ARM_Vector*) StrikeLow->Clone();
	}
	ARM_Vector* GetStrikeLow(void) { return itsStrikeLow;}

	void SetStrikeHigh(ARM_Vector* StrikeHigh) 
	{
		if (itsStrikeHigh)	delete itsStrikeHigh;
		if (StrikeHigh) itsStrikeHigh = (ARM_Vector*) StrikeHigh->Clone();
	}
	ARM_Vector* GetStrikeHigh(void) { return itsStrikeHigh;}

	void SetMktBid(ARM_Vector* MktBid) 
	{
		if (itsMktBid)	delete itsMktBid;
		if (MktBid) itsMktBid = (ARM_Vector*) MktBid->Clone();
	}
	ARM_Vector* GetMktBid(void) { return itsMktBid;}

	void SetMktAsk(ARM_Vector* MktAsk) 
	{
		if (itsMktAsk)	delete itsMktAsk;
		if (MktAsk) itsMktAsk = (ARM_Vector*) MktAsk->Clone();
	}
	ARM_Vector* GetMktAsk(void) { return itsMktAsk;}

	void SetInitialCorrelation(ARM_Vector* InitialCorrelation) 
	{
		if (itsInitialCorrelation)	delete itsInitialCorrelation;
		if (InitialCorrelation) itsInitialCorrelation = (ARM_Vector*) InitialCorrelation->Clone();
	}
	ARM_Vector* GetInitialCorrelation(void) { return itsInitialCorrelation;}

	void SetLeverages(ARM_Vector* Leverages) 
	{
		if (itsLeverages)	delete itsLeverages;
		if (Leverages) itsLeverages = (ARM_Vector*) Leverages->Clone();
	}
	ARM_Vector* GetLeverages(void) { return itsLeverages;}

	void SetImpliedCorrelation(ARM_VolCurve* ImpliedCorrelation) 
	{
		if (itsImpliedCorrelation)	delete itsImpliedCorrelation;
		if (ImpliedCorrelation) itsImpliedCorrelation = (ARM_Vector*) ImpliedCorrelation->Clone();
	}
	ARM_Vector* GetImpliedCorrelation(void) { return itsImpliedCorrelation;}
	
	void SetModelMultiCurves(ICM_ModelMultiCurves* ModelMultiCurves) ;

	ICM_ModelMultiCurves* GetModelMultiCurves() { return itsModelMultiCurves;}

	virtual ~ICM_Index_Correlation() ;

	// ----------------------------
	//	Copy of members data
	// ----------------------------
	void BitwiseCopy(const ARM_Object* src) ;

	// -------------
	//	Copy Method 
	// -------------
	void Copy(const ARM_Object* src)
	{
		ICM_Correlation::Copy(src);
		BitwiseCopy(src);
	}

	// --------------
	//	Clone Method
	// --------------
	ARM_Object* Clone(void)
	{
     ICM_Index_Correlation* theClone = new ICM_Index_Correlation();
     theClone->Copy(this);
     return(theClone);
	}

	void View(char* id, FILE* ficOut);

	void ValuationRF(double* x,double*y);
	void ValuationRF_2factors_AllQuotes(double* x,double* y);
	void ValuationRF_3factors_AllQuotes(double* x,double* y);
	double ValuationBaseCorrelation(const double& x);
	double ValuationBaseCorrelation_TermStructure(const double& x);

	void CalibrateBC(ARM_CLASS_NAME name,ICM_Parameters& Params,qCAL_INDEX_CORR_TYPE caltype);
	void CalibrateBC_TermStructure(ARM_CLASS_NAME name,ICM_Parameters& Params,qCAL_INDEX_CORR_TYPE caltype);
	void CalibrateRF(ARM_CLASS_NAME name,ICM_Parameters& Params,qCAL_INDEX_CORR_TYPE caltype);
	void CalibrateRF_Single(ARM_CLASS_NAME name,ICM_Parameters& Params,qCAL_INDEX_CORR_TYPE caltype);

	virtual double GetCorrelation(const std::string&  issuer1,
								  const std::string&  issuer2,
								  double maturity = CREDIT_DEFAULT_VALUE,
								  double strike = CREDIT_DEFAULT_VALUE) ;
 

	//JLA: not the same interface as base ?? 
	virtual double GetBeta(const std::string&issuer,
						   double maturity = CREDIT_DEFAULT_VALUE,
						   double strike = CREDIT_DEFAULT_VALUE)
	{
		double correlation = GetCorrelation(issuer,issuer,maturity,strike);
		correlation = sqrt(correlation);

		return (correlation);
	}

	virtual ARM_Vector* ComputeBetas(int nbissuers,const std::vector<std::string>& labels,const ARM_Vector&nominals, const ARM_Date &Maturity)
	{
		return NULL;
	}

	  void GetFactors(ICM_Parameters& p,ARM_Vector& v);
	 

	double HedgeImpact(ARM_Date& Start,ARM_Date& Maturity,const double& Leverage,const double& RefSpread) ;

	// function for gen implyCorrel in a file !
	void GenImplyCorrelation(char* id, FILE* ficOut = NULL);

	inline ICM_Parameters* GetParams() {return its_cal_Param;}
	inline void SetRFLBeta0(double& value) {its_RFL_Beta0=value;}
	inline double SetRFLBeta0() {return its_RFL_Beta0;}

	inline void SetParams(ICM_Parameters* value) {itsInParams=value;}

};

ICM_Mez* CdoIndexDefinition(const ARM_Date& Start,const ARM_Date& Maturity, const std::vector<std::string>&Labels,
							  // const int& size,
							  const double& K1,const double& K2,
							  const double& spread,double& cdoamount,
							  bool incmatu /*= INCLUDE_MATURITY*/,bool adjstart/* = false*/,
							  int creditlag /*= DEFAULT_CREDIT_LAG_INDX*/, const string ccy);

#endif