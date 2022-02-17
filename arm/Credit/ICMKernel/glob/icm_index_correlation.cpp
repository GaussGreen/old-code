#include "ARMKernel/glob/firsttoinc.h" 
#include <nag.h>
#include "ICMKernel\glob\icm_index_correlation.h"
#include "ICMKernel\inst\icm_mez.h"
#include "ICMKernel\pricer\icm_pricer_adviser.h"
#include "ICMKernel\crv\icm_defaultcurve.h"
#include "ICMKernel\glob\icm_smile_correlation.h"
#include "ICMKernel\mod\modelmulticurves.h"
#include "ICMKernel\util\icm_RootFinderND.h"
//#include "ICMKernel/util/icm_macro.h" 
//#include "ICMKernel\crv\icm_distribloss.h"
#include "ICMKernel\pricer\ICM_Pricer_homogeneous_smile_rf.h"
#include "ICMKernel\glob\icm_betas_correlation.h"
#include "ICMKernel\util\icm_rootfinder1D.h"
#include "ICMKernel\pricer\icm_pricer_cds.h"
//#include "ICMKernel\mod\modelmulticurves.h"
#include "ICMKernel\inst\icm_credit_index.h"
#include "ARMKernel\util\interpol.h"
#include "ICMKernel\pricer\icm_pricer_index.h"
#include "ICMKernel\glob\icm_calibrator.h"

#include <vector>
#include <math.h>
#include <stdio.h>


#include <nage04.h>

static void __stdcall objfun_RF_Allfactors_AllQuotes(Integer n, double x[], double *objf,
                   double g[], Nag_Comm *comm);

static void __stdcall objfun_RF_Betafactors_AllQuotes(Integer n, double x[], double *objf,
                   double g[], Nag_Comm *comm);
//-----------------------------------------------------------------------------
void ICM_Index_Correlation::Init(void)
{
	ICM_Correlation::Init(); 
	SetName(ICM_INDEX_CORRELATION);
	itsCalMethod = qCAL_BASE_CORRELATION;
	itsStrikeLow=NULL;
	itsStrikeHigh=NULL;
	its_cal_Pricer = NULL;
	its_cal_Param=NULL;
	its_cal_model=NULL;

	itsCreditIndex = NULL;
	itsPrevCreditIndex.clear();
	itsMktBid=NULL;
	itsMktAsk=NULL;
	itsUpfBid=NULL;
	itsUpfAsk=NULL;
	itsLeverages=NULL;

	itsInitialCorrelation=NULL;
	itsImpliedCorrelation = NULL;
	
	its_cal_CurrentPrice = 0.;
	its_cal_cdo_Notional = 0.;
			
	itsImpliedParameters.clear();
	its_cal_securitites.clear();
	its_cal_CurrentPrices.clear();
	its_cal_thres=0.;
	its_cal_Hedge=0.;

	itsModelMultiCurves=NULL;

	itsIntegrationStep=60;	
	itsLagForStartDate=0;	
	itsCreditLag=DEFAULT_CREDIT_LAG_INDX;		
	its_cal_PrevPricer = NULL;
	itsIsLinearInterpol = false;

	its_ts_prev_term = 0.;
	//its_ts_prev_correl_dw= 0.;
	//its_ts_prev_correl_up= 0.;
	its_ts_last_EL_prevMaturity = 0.;

	itsTemplateBaseCorrel = NULL;
	itsStep=0.;
	itsCalMeth = qNEWTON;
	itsInParams = NULL;
}
ICM_Index_Correlation::~ICM_Index_Correlation() 
	{
	 if (itsCreditIndex) delete itsCreditIndex;
	 if (itsStrikeLow) delete itsStrikeLow;
	 if (itsStrikeHigh) delete itsStrikeHigh;

	 if (itsMktBid) delete itsMktBid;
	 if (itsMktAsk) delete itsMktAsk;
	 if (itsInitialCorrelation) delete itsInitialCorrelation;
	 if (itsImpliedCorrelation)	delete itsImpliedCorrelation;

	if (its_cal_Pricer) its_cal_Pricer=NULL;
	if (its_cal_PrevPricer) its_cal_PrevPricer=NULL;
	
	if (its_cal_securitites.size()>0)
	{
		for (int i=0; i<its_cal_securitites.size();i++)
		{if (its_cal_securitites[i]) delete its_cal_securitites[i];}
	}
	if (its_cal_Param) delete its_cal_Param;
	if (its_cal_model) delete its_cal_model;
	if (itsModelMultiCurves) delete itsModelMultiCurves;
	
	for (int il=0; il<itsPrevCreditIndex.size(); il++)
	{if (itsPrevCreditIndex[il]) delete itsPrevCreditIndex[il];itsPrevCreditIndex[il]=NULL;}

	if (itsTemplateBaseCorrel) delete itsTemplateBaseCorrel;

	 if (itsUpfBid) delete itsUpfBid;
	 if (itsUpfAsk) delete itsUpfAsk;
	 if (itsLeverages) delete itsLeverages;

	};
//	---------------------------------------------------------------------------------------------
void ICM_Index_Correlation::Set(const ARM_Date& AsOf,
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
			 ARM_Vector* InitialCorrelation ,
			 ICM_ModelMultiCurves* mmc ,
			 const int IntegrationStep,
			 const int LagStartDate,
			 const int CreditLag ,
			 vector<ICM_Credit_Index*> PrevCreditIndex ,
			 ICM_QMatrix<double> PrevBCMatrix,
			 const double step ,
			 const qOPTIMIZE_TYPE CalMeth ) 
	{
		std::vector<std::string> labels( CreditIndex->GetNbCurves() ); 
		for(int i=0;i<labels.size();i++) labels[i]=CreditIndex->GetLabels()[i] ;


		ICM_Correlation::Set(AsOf,labels,name,(ARM_IRIndex*)0,(ARM_IRIndex*)0); 

		itsCalMethod = CalMethod;
		SetCreditIndex(CreditIndex);
		SetStrikeLow(StrikeLow);
		SetStrikeHigh(StrikeHigh);
		SetMktBid(MktBid);
		SetMktAsk(MktAsk);
		SetInitialCorrelation(InitialCorrelation);
		SetUpfBid(UpfBid);
		SetUpfAsk(UpfAsk);
		SetLeverages(Leverage);
		SetModelMultiCurves(mmc);
		itsIntegrationStep=IntegrationStep;	
		itsLagForStartDate=LagStartDate;	
		itsCreditLag=CreditLag;

		itsPrevCreditIndex.resize(PrevCreditIndex.size());
		for ( i=0;i<PrevCreditIndex.size();i++)
		{if (PrevCreditIndex[i]) itsPrevCreditIndex[i]=(ICM_Credit_Index*)PrevCreditIndex[i]->Clone();}
		
		itsPrevBaseCorrels = PrevBCMatrix;
		itsStep = step;
		itsCalMeth = CalMeth;
	}

void ICM_Index_Correlation::GlobalCalibrate()
{
	ICM_Parameters Params;
	ARM_CLASS_NAME cname;
	double TS = (double)qTermStructure;

	if (itsInParams)
	{Params = *itsInParams;}

	switch (itsCalMethod)
	{
	case qCAL_PWL_CORREL:	
	case qCAL_PWC_CORREL:
	{
		ARM_Vector PIntegrationStep(1,51.);
		ARM_Vector PCopula(1,1.);
		ARM_Vector PIntegrationStep2(1,0.);
		ARM_Vector PFreedomDeg(1,0.);
		Params.Push(&PIntegrationStep,"INTEGRATION_STEP_1");
		Params.Push(&PCopula,"COPULA");
		Params.Push(&PIntegrationStep2,"INTEGRATION_STEP_2");
		Params.Push(&PFreedomDeg,"FREEDOM_DEGREE");
		ARM_Vector PparamsModel(GetInitialCorrelation()->GetSize() ,0.5);
		ICM_QMatrix<string> PModel(1,1);
		Params.Push(&PparamsModel,"RF_PARAMS");			
		ARM_Vector RFLParams(1 ,its_RFL_Beta0);
		Params.Push(&RFLParams,"RFL_BETA0");			
		(PModel)(0,0)="PWC";
		Params.Push(&PModel,"STR_RF_MODEL");			
		cname = ICM_PRICER_HOMOGENEOUS_SMILE_RF;
		CalibrateRF(cname,Params,itsCalMethod);
		break;
	}
	case qCAL_BASE_CORRELATION_TSR :
		TS = (double)qTermStructureR;
	case qCAL_BASE_CORRELATION_TS:
	{
		ARM_Vector PIntegrationStep(1,itsIntegrationStep);
		ARM_Vector PCopula(1,1.);
		ARM_Vector PIntegrationStep2(1,0.);
		ARM_Vector PFreedomDeg(1,0.);
		ARM_Vector TERMresc(1,TS);
		Params.Push(&PIntegrationStep,"INTEGRATION_STEP_1");
		Params.Push(&PCopula,"COPULA");
		Params.Push(&PIntegrationStep2,"INTEGRATION_STEP_2");
		Params.Push(&PFreedomDeg,"FREEDOM_DEGREE");
		Params.Push(&TERMresc,"TERMS_RESCALING");
		cname = ICM_PRICER_HOMOGENEOUS_SMILE;
		CalibrateBC_TermStructure(cname,Params,itsCalMethod);
		break;
	}
	case qCAL_BASE_CORRELATION_LINEAR_MATURITY:
		itsIsLinearInterpol = true;
	case qCAL_BASE_CORRELATION:
	default:
	{
		ARM_Vector PIntegrationStep(1,itsIntegrationStep);
		ARM_Vector PCopula(1,1.);
		ARM_Vector PIntegrationStep2(1,0.);
		ARM_Vector PFreedomDeg(1,0.);
		Params.Push(&PIntegrationStep,"INTEGRATION_STEP_1");
		Params.Push(&PCopula,"COPULA");
		Params.Push(&PIntegrationStep2,"INTEGRATION_STEP_2");
		Params.Push(&PFreedomDeg,"FREEDOM_DEGREE");
		cname = ICM_PRICER_HOMOGENEOUS_SMILE;
		CalibrateBC(cname,Params,itsCalMethod);
		break;
	}
	}

}
//-----------------------------------------------------------------------------
//Method View
//-----------------------------------------------------------------------------
void ICM_Index_Correlation::View(char* id, FILE* ficOut)
{
	FILE* fOut;
	char  fOutName[200];
	int i = 0;

	 if ( ficOut == NULL )
	{
		ARM_GetViewFile(ARM_VIEW_FILE, id, fOutName);
		(void) unlink(fOutName);
		fOut = fopen(fOutName, "w"); 
	}
	else
	{
		fOut = ficOut;
	} 

	int size = GetSize();

	fprintf(fOut, "------------------------------------------------\n");
	fprintf(fOut, "---                Index Correlation         ---\n");
	fprintf(fOut, "------------------------------------------------\n\n");

	fprintf(fOut, "Calibration Type:");
	switch (itsCalMethod)
	{
	case qCAL_BASE_CORRELATION :
		fprintf(fOut, "Base Correlation Type\n");
		break;
	case qCAL_PWC_CORREL :
		fprintf(fOut, "PWC Correlation Type\n");
		break;
	}

	fprintf(fOut, "\n");
	fprintf(fOut,"StrLow\t\tStrHigh\t\tMktBid\t\tMktAsk\t\t");

	if (itsInitialCorrelation) fprintf(fOut,"InitCor\t\t");

	switch (itsCalMethod)
	{
	case qCAL_BASE_CORRELATION :
		if (itsImpliedCorrelation) fprintf(fOut,"ImplCorr");
		break;
	case qCAL_PWC_CORREL :
	case qCAL_PWL_CORREL :
		for (int il=0;il<itsImpliedParameters[0].size();il++) {fprintf(fOut, "Param%d\t\t",il);}
		break;
	}

	fprintf(fOut,"\n");

	for (i=0;i<itsStrikeLow->GetSize();i++)
	{
		fprintf(fOut,"%.2lf\t\t",itsStrikeLow->Elt(i));
		fprintf(fOut,"%.2lf\t\t",itsStrikeHigh->Elt(i));
		fprintf(fOut,"%.4lf\t\t",itsMktBid->Elt(i));
		fprintf(fOut,"%.4lf\t\t",itsMktAsk->Elt(i));
		if (itsInitialCorrelation){if (itsInitialCorrelation->GetSize()>i) fprintf(fOut,"%.4lf\t\t",itsInitialCorrelation->Elt(i)); else fprintf(fOut,"\t\t");}

		switch (itsCalMethod)
		{	
		case qCAL_BASE_CORRELATION :
		if (itsImpliedCorrelation){ fprintf(fOut,"%.4lf",itsImpliedCorrelation->Elt(i));}
		break;
		case qCAL_PWC_CORREL :
		case qCAL_PWL_CORREL :
		if (i<itsImpliedParameters.size())
		{for (int il=0;il<itsImpliedParameters[i].size();il++)	
		{ fprintf(fOut,"%.4lf\t\t",itsImpliedParameters[i][il]);}}
		break;
		}

		fprintf(fOut,"\n");
	}

	ICM_Correlation::View(id,fOut);
	itsCreditIndex->View(id,fOut);

	if ( ficOut == NULL )
	{
		fclose(fOut);
	}

}

//-----------------------------------------------------------------------------
//Generation of Index Cdo
//-----------------------------------------------------------------------------

ICM_Mez* CdoIndexDefinition(const ARM_Date& Start,
						   const ARM_Date& Maturity,
						   // char** Labels_,
						   const std::vector<std::string>& Labels,
						   // const int& size,
						   const double& K1,
						   const double& K2,
						   const double& spread,
						   double& cdoamount,
						   bool incmatu,
						   bool adjstart,
						   int creditlag,
						   const string ccy)
{
	double IndivNot = 10000000.;
	std::vector<double> IssuersNotionals(Labels.size(),IndivNot) ; 
	// double*	IssuersNotionals = new double[size];
	// for (int i=0; i<size; i++) {IssuersNotionals[i] = IndivNot;}
	double SubAmount = Labels.size()*IndivNot*K1;
	double CDOAmount = Labels.size()*IndivNot*(K2-K1);
	cdoamount = CDOAmount;

	// std::vector<std::string> Labels(size); 
	// for(i=0;i<size;i++) Labels[i]=Labels_[i] ;
	ICM_Mez* cdo = new ICM_Mez(Start,
						   Maturity,
						   (ARM_Date*)0,
						   (ARM_Date*)0,
						   spread,
						   K_MATUNADJUSTED,		// intRule
						   adjstart,	// adjStartDate
						   SubAmount,
						   CDOAmount,
						   Labels,
						   IssuersNotionals,
						   // size,
						   // "NULL",
						   K_QUARTERLY,
						   KACTUAL_360,
						   CDOAmount,
						   qACCRUED_SETTLED,
						   ccy, /*ARM_DEFAULT_COUNTRY,*/
						   CDOAmount,
						   K_SHORTSTART,
						   creditlag,
						   K_QUARTERLY,	//freq defleg
						   CREDIT_DEFAULT_VALUE, //binary
						   "EUR",
						   qRunning_Leg,
						   qStandart_Recovery_Leg,
  						   incmatu);

	cdo->SetMezzAmount(CDOAmount);
	//if (IssuersNotionals) delete[] IssuersNotionals;
	return (cdo);
}


//-----------------------------------------------------------------------------
//Valuation Price for base correlation
//-----------------------------------------------------------------------------
double ICM_Index_Correlation::ValuationBaseCorrelation(const double& x)
{
	const ICM_DefaultCurve* defcurve = itsCreditIndex->GetForcedCurve();
	ARM_Date StartDate = defcurve->GetAsOfDate();

	ICM_Smile_Correlation* Correlation = FixedBaseCorrelation(StartDate,
															  &ARM_Currency(defcurve->GetCurrency().c_str()),
															  its_cal_K1,
															  its_cal_K2,
															  its_cal_CorrelDown,
															  x,
															  "NAME_CORR");
	ICM_ModelMultiCurves* model = (ICM_ModelMultiCurves*) its_cal_Pricer->GetModel();
	model->SetCorrelation(Correlation);

	delete Correlation;

	its_cal_Pricer->ResetPricer();
	double value = -its_cal_Pricer->Price(qCMPPRICE);
	value-= (its_cal_CurrentPrice-its_cal_Hedge);

	return (value);
}


//-----------------------------------------------------------------------------
// calibration function for Base Correlation
//-----------------------------------------------------------------------------
void ICM_Index_Correlation::CalibrateBC(ARM_CLASS_NAME name,
										ICM_Parameters& Params,
										qCAL_INDEX_CORR_TYPE caltype)
{
	double IndivNot = 10000000.;
	const ICM_DefaultCurve* defcurve = itsCreditIndex->GetForcedCurve();
	ARM_Date Maturity = itsCreditIndex->GetMaturity();
	// char** Collateral = itsCreditIndex->GetLabels();
	const std::vector<std::string>& Collateral = itsCreditIndex->GetLabels();
	ARM_Date StartDate = defcurve->GetAsOfDate();
	StartDate.AddDays(itsLagForStartDate);

	int sizeindx = itsCreditIndex->GetNbCurves();
	ICM_Mez* cdo = NULL;
	double _inf=1.e-4,_sup=0.75,slope=0.,K1=0.,K2=0.,MidPrice=0.,result=0.,spread=0.,guess=1.e-4;
	int IntegrationStep = 60;
	ARM_Vector* VParams=NULL;
	double step = itsStep;
	bool autoguess = false;

	if (its_cal_Param) delete its_cal_Param;
	its_cal_Param = (ICM_Parameters*) Params.Clone();
	
	if (itsInitialCorrelation->GetSize()>0) 
	{if (itsInitialCorrelation->Elt(0)>0.)
		{guess=itsInitialCorrelation->Elt(0);}}

	ICM_Beta_Correlation Correlation(StartDate,1.,"NONE",(ARM_IRIndex*)0,(ARM_IRIndex*)0); 
	
	if (itsImpliedCorrelation) delete itsImpliedCorrelation;
	itsImpliedCorrelation = new ARM_Vector(itsStrikeLow->GetSize(),0.);

	
	// cas homogène
	std::vector<const ICM_DefaultCurve*> DefaultCurves(1); 
	// const ICM_DefaultCurve* DefaultCurves[1]; 
	DefaultCurves[0] = defcurve; // from index
	ICM_ModelMultiCurves mmc(DefaultCurves,defcurve->GetZeroCurve(),NULL,&Correlation);
	
	its_cal_CorrelDown=0.; 
	its_cal_CurrentPrices.clear();
	its_cal_securitites.clear();
	for (int i=0; i<itsStrikeLow->GetSize();i++)
	{
		//cas correlation down invalide
		if (its_cal_CorrelDown==CREDIT_DEFAULT_VALUE)
		{   its_cal_CorrelDown=itsImpliedCorrelation->Elt(i)=CREDIT_DEFAULT_VALUE;
		continue; }

		its_cal_K1=K1=itsStrikeLow->Elt(i);
		its_cal_K2=K2=itsStrikeHigh->Elt(i);
		
		spread=(itsMktBid->Elt(i)+itsMktAsk->Elt(i))/2.;
		its_cal_CurrentPrice=(itsUpfBid->Elt(i)+itsUpfAsk->Elt(i))/2.*(K2-K1)*sizeindx*IndivNot;
		its_cal_CurrentPrices.push_back(its_cal_CurrentPrice);
		if (!cdo) cdo = CdoIndexDefinition(StartDate,Maturity,Collateral,K1,K2,spread,its_cal_cdo_Notional,INCLUDE_MATURITY,false,itsCreditLag, defcurve->GetCurrency());
		if (its_cal_Pricer) delete its_cal_Pricer;its_cal_Pricer=NULL;

		ICM_Pricer_Advisor Advisor;
		if (itsModelMultiCurves) // cas non homogène
			its_cal_Pricer = Advisor.GeneratePricer(cdo,itsModelMultiCurves,name,CREDIT_DEFAULT_VALUE,&Params,GetAsOfDate());
		else
			its_cal_Pricer = Advisor.GeneratePricer(cdo,&mmc,name,CREDIT_DEFAULT_VALUE,&Params,GetAsOfDate());

		its_cal_Hedge=0.;
		if (itsLeverages->GetSize()>i)
		{
			its_cal_Hedge = HedgeImpact(StartDate,Maturity,itsLeverages->Elt(i),itsCreditIndex->GetSpread());
		}

		_inf = guess;
		_sup = MIN(guess+0.1,0.9);

		bool output = RootFinder1D(ff1::mem_call(&ICM_Index_Correlation::ValuationBaseCorrelation,(*this))).ZeroBracketMax(_inf,_sup,guess,0.999);	

		if (output==false)
		{	
			output = RootFinder1D(ff1::mem_call(&ICM_Index_Correlation::ValuationBaseCorrelation,(*this))).ZeroBracketMax(_inf,_sup,1.e-4,0.999);	

			if (output==false)
			{
			its_cal_CorrelDown=itsImpliedCorrelation->Elt(i)=CREDIT_DEFAULT_VALUE;
			if (its_cal_Pricer) delete its_cal_Pricer;its_cal_Pricer=NULL;
			if (cdo) delete cdo;cdo=NULL;
			continue; 
			}
		}

		guess = (_inf+_sup)/2.;

		if (itsCalMeth==qNEWTON)
		{result = RootFinder1D(ff1::mem_call(&ICM_Index_Correlation::ValuationBaseCorrelation,(*this))).NewtonRaphsonWithBisection(_inf,_sup,100,0.01,1.E-5,0.1);}
		else
		{result = RootFinder1D(ff1::mem_call(&ICM_Index_Correlation::ValuationBaseCorrelation,(*this))).Dichotomy(_inf,_sup,100,1.E-5,1.E-5,false);}
		
		guess=its_cal_CorrelDown=itsImpliedCorrelation->Elt(i)=result;

		if (itsInitialCorrelation)
		{	
		 if (itsInitialCorrelation->GetSize()>i+1)
			guess=itsInitialCorrelation->Elt(i+1);
		}	

		if (its_cal_Pricer) delete its_cal_Pricer;its_cal_Pricer=NULL;
		if (cdo) delete cdo;cdo=NULL;
	}

}


//-----------------------------------------------------------------------------
//Valuation Price for base correlation
//-----------------------------------------------------------------------------
double ICM_Index_Correlation::ValuationBaseCorrelation_TermStructure(const double& x)
{
	const ICM_DefaultCurve* defcurve = itsCreditIndex->GetForcedCurve();
	ARM_Date StartDate = defcurve->GetAsOfDate();

	ICM_Smile_Correlation* Correlation = NULL;
	ARM_Date Maturity = itsCreditIndex->GetMaturity();
	double yfmaturity = (Maturity-StartDate)/K_YEAR_LEN;

	bool succeed = false;
	Correlation = itsTemplateBaseCorrel;
	UpdateSingleCorrelation(Correlation,its_cal_K1,its_cal_K2,yfmaturity,x,succeed);
	if (succeed==false)
		{ICMLOG("ICM_Index_Correlation::ValuationBaseCorrelation_TermStructure failed for "<<" K1="<<its_cal_K1<<" K2="<<its_cal_K2<<" Maturity ="<<yfmaturity);}

	ICM_ModelMultiCurves* model = (ICM_ModelMultiCurves*) its_cal_Pricer->GetModel();
	model->SetCorrelation(Correlation);

	its_cal_Pricer->ResetPricer();
	double value = -its_cal_Pricer->Price(qCMPPRICE);
	value-= (its_cal_CurrentPrice-its_cal_Hedge);

	return (value);
}

//-----------------------------------------------------------------------------
//Valuation Price for base correlation
//-----------------------------------------------------------------------------

void ICM_Index_Correlation::ValuationRF(double* x,double* y)
{ //to modified
   ICM_ModelMultiCurves* model = (ICM_ModelMultiCurves*) its_cal_Pricer->GetModel();

	ICM_Pricer_Advisor Advisor;
	ICM_Parameters Params =  its_cal_Pricer->GetParameters();
	ARM_Vector* MParams = Params.GetColVect("RF_PARAMS");

	//Beta1 & Beta2 must be in [-1,1]
	if (x[0]<=0.) x[0]=0.;
	if (x[0]>=1.) x[0] = 0.9999;
	if (x[1]<=0.) x[1] = 0.;
	if (x[1]>=1.) x[1] = 0.9999;
   
	for (int i=0;i<MParams->GetSize();i++) {MParams->Elt(i)=x[i];}

	ICM_Pricer* Pricer = Advisor.GeneratePricer(its_cal_Pricer->GetSecurity(),
									model,
									ICM_PRICER_HOMOGENEOUS_SMILE_RF,
									CREDIT_DEFAULT_VALUE,
									&Params,
									GetAsOfDate());
	
	double value = -Pricer->Price(qCMPPRICE);
	value-=its_cal_CurrentPrice;
	y[0] =y[1] =y[2] = value;

	if (Pricer) delete Pricer;
}

//-----------------------------------------------------------------------------
//Valuation 2 factors/All quotes Price for base correlation
//-----------------------------------------------------------------------------

void ICM_Index_Correlation::ValuationRF_2factors_AllQuotes(double* x,double* y)
{  //to modified
	ICM_Parameters Params = (*its_cal_Param);
	ARM_Vector* MParams = Params.GetColVect("RF_PARAMS");
	double value = 0.,diff=0.;
	double IndivNot = 10000000.;
	int sizeindx = itsCreditIndex->GetNbCurves();


	//Beta1 & Beta2 must be in [-1,1]
	if (x[0]<=-0.) x[0]=0.;
	if (x[0]>=1.) x[0] = 0.9999;
	if (x[1]<=-0.) x[1] = 0.;
	if (x[1]>=1.) x[1] = 0.9999;
   
	for (int i=0;i<2;i++) {MParams->Elt(i)=x[i];}
	MParams->Elt(2) = its_cal_thres;

	for (int il=0;il<its_cal_securitites.size();il++)
	{
		ICM_Pricer_Distrib_RandFactors Pricer; Pricer.Set(its_cal_securitites[il],
										its_cal_model,
										Params,GetAsOfDate());

		diff = -Pricer.Price(qCMPPRICE) - its_cal_CurrentPrices[il];
		diff /= (itsStrikeHigh->Elt(il)-itsStrikeLow->Elt(il)) * sizeindx * IndivNot;
		value += diff*diff;
	}
	
	y[0] =y[1] = value;
}


//-----------------------------------------------------------------------------
//Valuation 3 factors/All quotes Price for base correlation
//-----------------------------------------------------------------------------

void ICM_Index_Correlation::ValuationRF_3factors_AllQuotes(double* x,double* y)
{  //to modified
	ICM_Parameters Params = (*its_cal_Param);
	ARM_Vector* MParams = Params.GetColVect("RF_PARAMS");
	double value = 0.,diff=0.;
	double IndivNot = 10000000.;
	int sizeindx = itsCreditIndex->GetNbCurves();

	//Beta1 & Beta2 must be in [-1,1]
	if (x[0]<=-0.) x[0]=0.;
	if (x[0]>=1.) x[0] = 0.9999;
	if (x[1]<=-0.) x[1] = 0.;
	if (x[1]>=1.) x[1] = 0.9999;
   
	for (int i=0;i<MParams->GetSize();i++) {MParams->Elt(i)=x[i];}

	for (int il=0;il<its_cal_securitites.size();il++)
	{
		ICM_Pricer_Distrib_RandFactors Pricer; Pricer.Set(its_cal_securitites[il],
										its_cal_model,
										Params,
										GetAsOfDate());

		diff = -Pricer.Price(qCMPPRICE) - its_cal_CurrentPrices[il];
		diff /= (itsStrikeHigh->Elt(il)-itsStrikeLow->Elt(il)) * sizeindx * IndivNot;
		value += diff*diff;
	}
	
	y[0] =y[1] =y[2] = value;
}

//-----------------------------------------------------------------------------
// calibration function for Base Correlation
//-----------------------------------------------------------------------------

void ICM_Index_Correlation::CalibrateBC_TermStructure(ARM_CLASS_NAME name,
										ICM_Parameters& Params,
										qCAL_INDEX_CORR_TYPE caltype)
{
	double IndivNot = 10000000.;
	const ICM_DefaultCurve* defcurve = itsCreditIndex->GetForcedCurve();
	// char** Collateral = itsCreditIndex->GetLabels();
	const std::vector<std::string>& Collateral = itsCreditIndex->GetLabels();

	ARM_Date AsOf = defcurve->GetAsOfDate();
	ARM_Date StartDate = defcurve->GetAsOfDate();StartDate.AddDays(itsLagForStartDate);
	ARM_Date Maturity = itsCreditIndex->GetMaturity();
	ARM_Date PreviousMaturity = Maturity;

	//Gestion du calcul de la Base Correlation par term
	//if (itsPrevCreditIndex.size()>0) {PreviousMaturity=itsPrevCreditIndex[0]->GetMaturity();}

	int sizeindx = itsCreditIndex->GetNbCurves();
	ICM_Mez* cdo = NULL;
	double _inf=1.e-4,_sup=0.75,slope=0.,K1=0.,K2=0.,MidPrice=0.,result=0.,spread=0.,guess=0.0001;
	int IntegrationStep = 60;
	ARM_Vector* VParams=NULL;
	double step = itsStep;
	bool autoguess = false;

	if (its_cal_Param) delete its_cal_Param;
	its_cal_Param = (ICM_Parameters*) Params.Clone();
	
	if (itsInitialCorrelation->GetSize()>0) 
	{if (itsInitialCorrelation->Elt(0)>0.)
		{guess=itsInitialCorrelation->Elt(0);}}

	std::vector<const ICM_DefaultCurve*> DefaultCurves(1); 
	DefaultCurves[0] = defcurve;
	ICM_Beta_Correlation Correlation(StartDate,1.,"NONE",(ARM_IRIndex*)0,(ARM_IRIndex*)0); 
	
	if (itsImpliedCorrelation) delete itsImpliedCorrelation;
	itsImpliedCorrelation = new ARM_Vector(itsStrikeLow->GetSize(),0.);

	ICM_ModelMultiCurves mmc(DefaultCurves,defcurve->GetZeroCurve(),NULL,&Correlation);
	its_cal_CorrelDown=0.; 
	its_cal_CurrentPrices.clear();
	its_cal_securitites.clear();

	int r=0,r2=0;

	its_cal_yts.Resize(itsPrevCreditIndex.size());
	if (itsPrevCreditIndex.size()>0){
	for (r=0;r<itsPrevCreditIndex.size();r++)
	{	double yearterm = (itsPrevCreditIndex[r]->GetMaturity()-AsOf)/K_YEAR_LEN;
	its_cal_yts.Elt(r)=yearterm;}}

	its_cal_yts.push_back((Maturity-AsOf)/K_YEAR_LEN);

	its_cal_strikes.Resize(itsStrikeLow->GetSize());
	for (r=0;r<its_cal_strikes.size();r++)
	{	its_cal_strikes.Elt(r)=itsStrikeLow->Elt(r);}

	its_cal_strikes.push_back(itsStrikeHigh->Elt(itsStrikeHigh->size()-1));

	its_cal_bc.Resize(its_cal_yts.size(),its_cal_strikes.size());
	for (r=0;r<its_cal_yts.size()-1;r++)
		for (r2=1;r2<its_cal_strikes.size();r2++)
		{	its_cal_bc.Elt(r,r2)=itsPrevBaseCorrels(r,r2-1);}

	if (itsTemplateBaseCorrel) delete itsTemplateBaseCorrel;

	itsTemplateBaseCorrel = FixedBaseCorrelationEmpty(AsOf,
												&ARM_Currency(defcurve->GetCurrency().c_str()),
												its_cal_strikes,
												its_cal_yts,
												its_cal_bc);

	for (int i=0; i<itsStrikeLow->GetSize();i++)
	{
		its_cal_K1=K1=itsStrikeLow->Elt(i);
		its_cal_K2=K2=itsStrikeHigh->Elt(i);
		
		spread=(itsMktBid->Elt(i)+itsMktAsk->Elt(i))/2.;
		its_cal_CurrentPrice=(itsUpfBid->Elt(i)+itsUpfAsk->Elt(i))/2.*(K2-K1)*sizeindx*IndivNot;

		if (cdo==NULL)
		{ cdo = CdoIndexDefinition(StartDate,Maturity,Collateral,K1,K2,spread,its_cal_cdo_Notional,INCLUDE_MATURITY,false,itsCreditLag, defcurve->GetCurrency()); }

		if (its_cal_Pricer) delete its_cal_Pricer;its_cal_Pricer=NULL;

		ICM_Pricer_Advisor Advisor;
		if (itsModelMultiCurves){
			its_cal_Pricer = Advisor.GeneratePricer(cdo,itsModelMultiCurves,name,CREDIT_DEFAULT_VALUE,&Params,GetAsOfDate());}
		else{
			its_cal_Pricer = Advisor.GeneratePricer(cdo,&mmc,name,CREDIT_DEFAULT_VALUE,&Params,GetAsOfDate());}

		its_cal_Hedge=0.;
		if (itsLeverages->GetSize()>i)
		{its_cal_Hedge = HedgeImpact(StartDate,Maturity,itsLeverages->Elt(i),itsCreditIndex->GetSpread());}


		_inf = guess;
		_sup = MIN(guess+0.1,0.9);

		bool output = RootFinder1D(ff1::mem_call(&ICM_Index_Correlation::ValuationBaseCorrelation_TermStructure,(*this))).ZeroBracketMax(_inf,_sup,guess,0.999);	

		if (output==false)
		{	
			output = RootFinder1D(ff1::mem_call(&ICM_Index_Correlation::ValuationBaseCorrelation_TermStructure,(*this))).ZeroBracketMax(_inf,_sup,1.e-4,0.999);	

			if (output==false)
			{
			its_cal_CorrelDown=itsImpliedCorrelation->Elt(i)=CREDIT_DEFAULT_VALUE;
			if (its_cal_Pricer) delete its_cal_Pricer;its_cal_Pricer=NULL;
			if (cdo) delete cdo;cdo=NULL;
			continue; 
			}
		}

		guess = (_inf+_sup)/2.;

		if (itsCalMeth==qNEWTON)
		{result = RootFinder1D(ff1::mem_call(&ICM_Index_Correlation::ValuationBaseCorrelation_TermStructure,(*this))).NewtonRaphsonWithBisection(_inf,_sup,100,0.01,1.E-5,0.1);}
		else
		{result = RootFinder1D(ff1::mem_call(&ICM_Index_Correlation::ValuationBaseCorrelation_TermStructure,(*this))).Dichotomy(_inf,_sup,100,1.E-5,1.E-5,false);}
		
		guess=its_cal_CorrelDown=itsImpliedCorrelation->Elt(i)=result;
		
		if (result == CREDIT_DEFAULT_VALUE)
			UpdateSingleCorrelation(itsTemplateBaseCorrel,its_cal_K1,its_cal_K2,(Maturity-AsOf)/K_YEAR_LEN,_inf,output);
		else
			UpdateSingleCorrelation(itsTemplateBaseCorrel,its_cal_K1,its_cal_K2,(Maturity-AsOf)/K_YEAR_LEN,result,output);

		if (output==false)
		{ICMLOG("ICM_Index_Correlation::CalibrateBC_TermStructure failed for "<<" K1="<<its_cal_K1<<" K2="<<its_cal_K2<<" Maturity ="<<(Maturity-AsOf)/K_YEAR_LEN);}
		
		if (result == CREDIT_DEFAULT_VALUE)
			its_cal_CorrelDown=itsImpliedCorrelation->Elt(i)=_inf;
		else
			its_cal_CorrelDown=itsImpliedCorrelation->Elt(i)=result;

		if (itsInitialCorrelation)
		{	
		 if (itsInitialCorrelation->GetSize()>i+1)
			guess=itsInitialCorrelation->Elt(i+1);
		}	

		if (its_cal_Pricer) delete its_cal_Pricer;its_cal_Pricer=NULL;
		if (cdo) delete cdo;cdo=NULL;
	}

}


//-----------------------------------------------------------------------------
// calibration function for Range Factor (single)
//-----------------------------------------------------------------------------

void ICM_Index_Correlation::CalibrateRF_Single(ARM_CLASS_NAME name,
											ICM_Parameters& Params,
											qCAL_INDEX_CORR_TYPE caltype)
{
	double IndivNot = 10000000.;
	const ICM_DefaultCurve* defcurve = itsCreditIndex->GetForcedCurve();
	ARM_Date Maturity = itsCreditIndex->GetMaturity();
	// char** Collateral = itsCreditIndex->GetLabels();
	const std::vector<std::string>& Collateral = itsCreditIndex->GetLabels();

	ARM_Date StartDate = defcurve->GetAsOfDate();
	int sizeindx = itsCreditIndex->GetNbCurves();
	ICM_Mez* cdo = NULL;
	double _inf=1.e-4,_sup=0.75,slope=0.,K1=0.,K2=0.,MidPrice=0.,result=0.,spread=0.,guess=0.15;
	int IntegrationStep = 60;
	ARM_Vector* VParams=NULL;

	if (its_cal_Param) delete its_cal_Param;
	its_cal_Param = (ICM_Parameters*) Params.Clone();
	
	if (itsInitialCorrelation){
		VParams = its_cal_Param->GetColVect("RF_PARAMS");
		for (int il=0;il<VParams->GetSize();il++) {VParams->Elt(il)=itsInitialCorrelation->Elt(il);}}

	std::vector<const ICM_DefaultCurve*> DefaultCurves(1); 
	DefaultCurves[0] = defcurve;
	ICM_Correlation* Correlation = new ICM_Beta_Correlation(StartDate,1.,"NONE",(ARM_IRIndex*)0,(ARM_IRIndex*)0); 
	
	if (itsImpliedCorrelation) delete itsImpliedCorrelation;
	itsImpliedCorrelation = new ARM_Vector(itsStrikeLow->GetSize(),0.);

	ICM_ModelMultiCurves mmc(DefaultCurves,defcurve->GetZeroCurve(),NULL,Correlation);
	its_cal_CorrelDown=0.; 
	its_cal_CurrentPrices.clear();
	its_cal_securitites.clear();
	if (its_cal_model) delete its_cal_model;

	if (itsModelMultiCurves) // cas non homogène
		its_cal_model = (ARM_Model*) itsModelMultiCurves->Clone();
	else
		its_cal_model = (ARM_Model*) mmc.Clone();


	for (int i=0; i<itsStrikeLow->GetSize();i++)
	{
		its_cal_K1=K1=itsStrikeLow->Elt(i);
		its_cal_K2=K2=itsStrikeHigh->Elt(i);
		
		spread=(itsMktBid->Elt(i)+itsMktAsk->Elt(i))/2.;
		its_cal_CurrentPrice=(itsUpfBid->Elt(i)+itsUpfAsk->Elt(i))/2.*(K2-K1)*sizeindx*IndivNot;
		its_cal_CurrentPrices.push_back(its_cal_CurrentPrice);
		cdo = CdoIndexDefinition(StartDate,Maturity,Collateral,K1,K2,spread,its_cal_cdo_Notional,INCLUDE_MATURITY,false,itsCreditLag, defcurve->GetCurrency());

		ICM_Pricer_Advisor Advisor;
		its_cal_Pricer = Advisor.GeneratePricer(cdo,&mmc,name,CREDIT_DEFAULT_VALUE,&Params,GetAsOfDate());

		switch (caltype)
		{
		case qCAL_PWC_CORREL:
		case qCAL_PWL_CORREL:
		default:
			{
			bool check=false;
			ARM_Vector Cparams = (*VParams);
			try{result = RootFinderND(ff1::mem_call(&ICM_Index_Correlation::ValuationRF,(*this))).NewtonMultiDim(Cparams.GetElt(),Cparams.GetSize(),100,1.e-5,1.);}
			catch (...)
			{ICMTHROW(ERR_INVALID_ARGUMENT,"unable to calibrate Indx("<<its_cal_K1<<","<<its_cal_K2<<")");}
			vector<double> vect;
			for (int il=0;il<Cparams.GetSize();il++) {vect.push_back(Cparams.Elt(il));}
			itsImpliedParameters.push_back(vect);
			break;
			}
		}

		if (its_cal_Pricer) delete its_cal_Pricer;its_cal_Pricer=NULL;
		if (cdo) delete cdo;
	}

	if (its_cal_model) delete its_cal_model;its_cal_model=NULL;
}



//-----------------------------------------------------------------------------
// calibration function for Range Factor
//-----------------------------------------------------------------------------

void ICM_Index_Correlation::CalibrateRF(ARM_CLASS_NAME name,
										ICM_Parameters& Params,
										qCAL_INDEX_CORR_TYPE caltype)
{

	qOPTIMIZE_TYPE opttype = qOPT_PSWARM;
	int nbiter_=50;
	
	if (itsInParams)
	{
		ARM_Vector* OPTTYPE = itsInParams->GetColVect("OPTTYPE");
		if (OPTTYPE) opttype = (qOPTIMIZE_TYPE)(int) OPTTYPE->Elt(0);

		ARM_Vector* NBITER = itsInParams->GetColVect("NBITER");
		if (NBITER) nbiter_ = (int) NBITER->Elt(0);
	}


	double IndivNot = 10000000.;
	const ICM_DefaultCurve* defcurve = itsCreditIndex->GetForcedCurve();
	ARM_Date Maturity = itsCreditIndex->GetMaturity();
	// char** Collateral = itsCreditIndex->GetLabels();
	const std::vector<std::string>& Collateral = itsCreditIndex->GetLabels();
	ARM_Date StartDate = defcurve->GetAsOfDate();
	int sizeindx = itsCreditIndex->GetNbCurves();
	ICM_Mez* cdo = NULL;
	double _inf=1.e-4,_sup=0.75,slope=0.,K1=0.,K2=0.,MidPrice=0.,result=0.,spread=0.,guess=0.15;
	int i=0,il=0,il1=0,il2=0;
	ARM_Vector VParams;
    
	//NAG parameters
    double objf;
    Nag_BoundType bound;
    
	NagError fail ; 
	INIT_FAIL(fail); // 

    // Initialise options structure 
	Nag_E04_Opt options ;
	ICM_NAGSILENT(options) ;
	options.init_state=Nag_Init_None; //JLA: ?? 
	options.optim_tol=1.e-5;
	options.max_iter=50;

	Nag_Comm comm;
	comm.p = (Pointer)(this) ;

	if (its_cal_Param) delete its_cal_Param;its_cal_Param=NULL;
	its_cal_Param = (ICM_Parameters*) Params.Clone();
	
	//Estimation de la taille du vecteur de paramètres (taille de Param et si pas de valeur : 3)
	int nb_params_defaut = 3;
	if (itsInitialCorrelation)
	{
		VParams = *its_cal_Param->GetColVect("RF_PARAMS");
		nb_params_defaut = VParams.GetSize();
		for (int il=0;il<VParams.GetSize();il++) {VParams.Elt(il)=itsInitialCorrelation->Elt(il);}
	}

	//Variable : fonction de calibration
	vector<double> xvals(nb_params_defaut),gvals(nb_params_defaut),bl(nb_params_defaut),bu(nb_params_defaut);
	
	std::vector<const ICM_DefaultCurve*> DefaultCurves(1); 
	DefaultCurves[0] = defcurve;
	ICM_Beta_Correlation Correlation(StartDate,1.,"NONE",(ARM_IRIndex*)0,(ARM_IRIndex*)0); 
	
	if (itsImpliedCorrelation) delete itsImpliedCorrelation;
	itsImpliedCorrelation = new ARM_Vector(itsStrikeLow->GetSize(),0.);

	ICM_ModelMultiCurves mmc(DefaultCurves,defcurve->GetZeroCurve(),NULL,&Correlation);

	its_cal_CorrelDown=0.; 
	its_cal_CurrentPrices.clear();
	its_cal_securitites.clear();
	if (its_cal_model) delete its_cal_model;

	if (itsModelMultiCurves) // cas non homogène
		its_cal_model = (ARM_Model*) itsModelMultiCurves->Clone();
	else
		its_cal_model = (ARM_Model*) mmc.Clone();

	//its_cal_model = (ARM_Model*) mmc.Clone();

	//Construction des tranches d'indices
	for (i=0; i<itsStrikeLow->GetSize();i++)
	{
		its_cal_K1=K1=itsStrikeLow->Elt(i);
		its_cal_K2=K2=itsStrikeHigh->Elt(i);
		
		spread=(itsMktBid->Elt(i)+itsMktAsk->Elt(i))/2.;
		its_cal_CurrentPrice=(itsUpfBid->Elt(i)+itsUpfAsk->Elt(i))/2.*(K2-K1)*sizeindx*IndivNot;
		its_cal_CurrentPrices.push_back(its_cal_CurrentPrice);
		cdo = CdoIndexDefinition(StartDate,Maturity,Collateral,K1,K2,spread,its_cal_cdo_Notional,INCLUDE_MATURITY,false,itsCreditLag, defcurve->GetCurrency());
		// its_cal_securitites.push_back((ARM_Security*) cdo);
		its_cal_securitites.push_back(dynamic_cast<ARM_Security*>(cdo));
	}

	bound = Nag_Bounds;
	//MaJ des bornes sur chacuns des paramètres de type Beta et seuil
	int idx = 0,mult = 0;

	if (nb_params_defaut%2)
	{

		for (i=0; i<nb_params_defaut ; i++)
		{
			idx = i % 3;

			if (idx<2)
			{
				if ((i==0) && (its_RFL_Beta0>=0.))
				{
					bl[i] = its_RFL_Beta0 - 1.E-5;
					bu[i] = its_RFL_Beta0 + 1.E-5;
				}
				else
				{
					bl[i] = 0.;
					bu[i] = 0.999;
				}
			}
			else if (idx == 2)
			{
				bl[i] = -3.;
				bu[i] = 3.;
			}
		}	
	}
	else 
	{
		int nbfac = nb_params_defaut/2;

		for (i=0; i<nb_params_defaut ; i++)
		{
			idx = i % nb_params_defaut;

			if (idx<nbfac)
			{
				if ((i==0) && (its_RFL_Beta0>=0.))
				{
					bl[i] = its_RFL_Beta0 - 1.E-5;
					bu[i] = its_RFL_Beta0 + 1.E-5;
				}
				else
				{
					bl[i] = 0.;
					bu[i] = 0.999;
				}
			}
			else if ((idx>=nbfac) && (idx<(2*nbfac-1)))
			{
				bl[i] = -3.;
				bu[i] = 3.;
			}
			else if (idx>=(2*nbfac-1))
			{
				bl[i] = 0.0001;
				bu[i] = 0.9999;
			}
		}	
	}

	double valeur_min = 1.e12;
	
	//Recherche du pint de départ (donné ou non)
	if (itsInitialCorrelation)
	{
		for (il=0;il<nb_params_defaut;il++) {xvals[il]=VParams.Elt(il);}
	}
	else
	{
		//Recherche d'un point de départ X0 valide : modèle a deux états
		int discretization = 5;	
// FIXMEFRED: mig.vc8 (28/05/2007 10:24:16):cast
		int nb_grid_x = pow(static_cast<double>(discretization),nb_params_defaut);

		vector<vector<double> > Param0(nb_grid_x);
		vector<double> min_param(nb_params_defaut), max_param(nb_params_defaut);
		vector<int> base(nb_params_defaut);
		int ind_xmin=0;
		double valeur=0.;

		//Init Bound
		for (i=0;i<nb_params_defaut;i++)
		{
			//Bounds
			if (i<(nb_params_defaut+1)/2)
			{
				min_param[i] = 0.;
				max_param[i] = 0.999;
			}
			else
			{
				min_param[i] = -3.;
				max_param[i] = 3.;
			}
		}

		//Init Grid
		for (int i=0;i<nb_grid_x ;i++)
		{
			Param0[i].resize(nb_params_defaut);
		
			//Mise a jour du jeu de param
			int base_init = i;
			for (int d=0;d<nb_params_defaut;d++)
			{
// FIXMEFRED: mig.vc8 (28/05/2007 10:24:40):cast
				base[d] = (int) (base_init / (int) (pow(static_cast<double>(nb_params_defaut),nb_params_defaut-1-d)));
				base_init -=base[d] * pow(static_cast<double>(nb_params_defaut),nb_params_defaut-1-d);
				Param0[i][d] = min_param[d] + base[d]*(max_param[d]-min_param[d])/(double)discretization;
			}
				
		}

		//Opt on the grid
		for (il=0; il<nb_grid_x; il++)
		{
			//Init Xvals
			xvals = Param0[il];
			
			//Nb Iter max
			options.max_iter=5;

			//Opt routine
// FIXMEFRED: mig.vc8 (28/05/2007 10:25:19):cast
			e04jbc(nb_params_defaut, objfun_RF_Allfactors_AllQuotes, bound, &(*bl.begin()), &(*bu.begin()), &(*xvals.begin()), &objf,
			        &(*gvals.begin()), &options, &comm, &fail);
			
			//Min atteint
			if (objf<valeur_min) 
			{
				ind_xmin = il;

				valeur_min=objf;
			}
		}
		
		//Return the best init point
		xvals = Param0[ind_xmin];
	}

	//Recherche d'une solution optimale
	options.max_iter=nbiter_;

	_context_cal ctxt;
	ctxt.m_evrythg = this;
	ctxt.m_Parameters = its_cal_Param;
	OptimTools::GansoLibContext = &ctxt;

	if (opttype == qOPT_PSWARM)
	{
	objf = OptimTools::MultiOptimisator(opttype,(void*)&ctxt,(void*)PSW_objfun_RF_Allfactors_AllQuotes,xvals,
		bl,bu,itsUpfBid->size(),1.e-3,nbiter_,0.,0.);

//	objf = OptimTools::AleaOptimSearch(opttype,(void*)&ctxt,(void*)PSW_objfun_RF_Allfactors_AllQuotes,xvals,
//		bl,bu,itsUpfBid->size(),1e20,720,1.e-5,nbiter_,1.e-5,0.5,0.);
	}
	else{
	// Call optimization routine 
	e04jbc(nb_params_defaut, objfun_RF_Allfactors_AllQuotes, bound, bl.begin(), bu.begin(), xvals.begin(), &objf,
		gvals.begin(), &options, &comm, &fail);}

	//Stockage du résultat pour l'affichage
	vector<double> res = xvals;
	res.push_back(objf);
	itsImpliedParameters.push_back(res);

	ICM_NAGFREE(options); 

	if (its_cal_model) delete its_cal_model;its_cal_model=NULL;

/*
    if (fail.code != NE_NOERROR && fail.code != NW_COND_MIN) 
		ICMTHROW(ERR_INVALID_ARGUMENT,"exit(EXIT_FAILURE);"); 
*/
}



static void __stdcall objfun_RF_Betafactors_AllQuotes(Integer n, double x[], double *objf,
                   double g[], Nag_Comm *comm)
{
  // Routine to evaluate objective function. 

  ICM_Index_Correlation* ptr = (ICM_Index_Correlation*)comm->p ;
	
  ICM_Parameters Params = (*ptr->its_cal_Param);
  ARM_Vector* MParams = Params.GetColVect("RF_PARAMS");
  double value = 0.,diff=0.;
  double IndivNot = 10000000.;
  int sizeindx = ptr->GetIndex()->GetNbCurves();


  if (ptr->its_cal_fixedparam==0)	
  {
	MParams->Elt(0)=ptr->its_cal_valueparam;
	MParams->Elt(1)=x[0];
	MParams->Elt(2)=x[1];
  }else if (ptr->its_cal_fixedparam==1)	
  {
	MParams->Elt(0)=x[0];
	MParams->Elt(1)=ptr->its_cal_valueparam;
	MParams->Elt(2)=x[1];
  }else			 
  {
	MParams->Elt(0)=x[0];
	MParams->Elt(1)=x[1];
	MParams->Elt(2)=ptr->its_cal_valueparam;
  }

  for (int il=0;il<ptr->its_cal_securitites.size();il++)
	{
		ICM_Pricer_Distrib_RandFactors Pricer; Pricer.Set(ptr->its_cal_securitites[il],
										ptr->its_cal_model,
										Params,ptr->GetAsOfDate());

		diff = -Pricer.Price(qCMPPRICE) - ptr->its_cal_CurrentPrices[il];
		diff /= (ptr->itsStrikeHigh->Elt(il)-ptr->itsStrikeLow->Elt(il)) * sizeindx * IndivNot;
		value += diff*diff;
	}
	
    value =sqrt(value);
	*objf = value;
}                              
 
// ----------------------------------------------------------------------------
// NAG RFL
// ----------------------------------------------------------------------------
static void __stdcall objfun_RF_Allfactors_AllQuotes(Integer n, double x[], double *objf,
                   double g[], Nag_Comm *comm)
{
  /* Routine to evaluate objective function. */

  ICM_Index_Correlation* ptr = (ICM_Index_Correlation*)comm->p ;
	
  ICM_Parameters Params = (*ptr->its_cal_Param);
  ARM_Vector* MParams = Params.GetColVect("RF_PARAMS");
  double value = 0.,diff=0.;
  double IndivNot = 10000000.;
  int sizeindx = ptr->GetIndex()->GetNbCurves();

  int idx=0;
  double ecart =0.;

  //Init Param Value
  for (int i=0; i<n; i++)
  {
	  if (n%2)
	  {
		idx = i % 3;
		MParams->Elt(i)=x[i];
		//Bound Test : [0,1] beta et [-3,3] thresholds(size_param+1)/2;
		if ((idx<2) && (x[i]<=0.) ) MParams->Elt(i)=0.0;
		if ((idx<2) && (x[i]>=1.) ) MParams->Elt(i)=0.9999;
		if ((idx==2) && (x[i]<=-3.) ) MParams->Elt(i)=-3.;
		if ((idx==2) && (x[i]>=3.) ) MParams->Elt(i)=3.;
	  }
	  else
	  {
		idx = i % 4;
		MParams->Elt(i)=x[i];
		//Bound Test : [0,1] beta et [-3,3] thresholds(size_param+1)/2;
		if ((idx<2) && (x[i]<=0.) ) MParams->Elt(i)=0.0;
		if ((idx<2) && (x[i]>=1.) ) MParams->Elt(i)=0.9999;
		if (idx==1) 
		{if (x[i-1]>=x[i]) 
			{ecart=x[i-1]-x[i];}
		}
		if ((idx==2) && (x[i]>=3.) ) MParams->Elt(i)=3.;
		if ((idx==2) && (x[i]<=-3.) ) MParams->Elt(i)=-3.;
		if ((idx==3) && (x[i]<0.) ) MParams->Elt(i)=0.001;
		if ((idx==3) && (x[i]>1.) ) MParams->Elt(i)=0.9999;
	  }
  }

  //Price Securities
  for (int il=0;il<ptr->its_cal_securitites.size();il++)
	{
		ICM_Pricer_Distrib_RandFactors Pricer; Pricer.Set(ptr->its_cal_securitites[il],
										ptr->its_cal_model,
										Params,ptr->GetAsOfDate());

		diff = -Pricer.Price(qCMPPRICE) - ptr->its_cal_CurrentPrices[il];
		diff /= (ptr->itsStrikeHigh->Elt(il)-ptr->itsStrikeLow->Elt(il)) * sizeindx * IndivNot;
		value += diff*diff;
	}
	
    ecart = 0.;
    value =sqrt(value)*exp(ecart);
	*objf = value;
} 

// --------------------------------------------------------------------------
// RFL PSWARM
// --------------------------------------------------------------------------
extern "C" struct Stats stats;

double PSW_objfun_RF_Allfactors_AllQuotes(int n, double *y, double *lb, double *ub)
{
	int i=0;

	double fx;
	vector<double> x;x.resize(n);
	_context_cal* p = (_context_cal*)OptimTools::GansoLibContext ;
	double* scal = OptimTools::scale_ = p->m_scale;
	ICM_Index_Correlation* ptr = (ICM_Index_Correlation*) p->m_evrythg;

	stats.objfunctions++;

	for(i=0;i<n;i++)
	{
	  if(y[i]<lb[i] || y[i]>ub[i])
	  {return 1e20;}
	  x[i]=y[i]*scal[i];
	}

  ICM_Parameters Params = (*ptr->its_cal_Param);
  ARM_Vector* MParams = Params.GetColVect("RF_PARAMS");
  double value = 0.,diff=0.;
  double IndivNot = 10000000.;
  int sizeindx = ptr->GetIndex()->GetNbCurves();

  int idx=0;
  double ecart =0.;

  //Init Param Value
  for (i=0; i<n; i++)
  {
	  if (n%2)
	  {
		idx = i % 3;
		MParams->Elt(i)=x[i];
		//Bound Test : [0,1] beta et [-3,3] thresholds(size_param+1)/2;
/*		if ((idx<2) && (x[i]<=0.) ) MParams->Elt(i)=0.0;
		if ((idx<2) && (x[i]>=1.) ) MParams->Elt(i)=0.9999;
		if ((idx==2) && (x[i]<=-3.) ) MParams->Elt(i)=-3.;
		if ((idx==2) && (x[i]>=3.) ) MParams->Elt(i)=3.;
*/	  }
	  else
	  {
		idx = i % 4;
		MParams->Elt(i)=x[i];
		//Bound Test : [0,1] beta et [-3,3] thresholds(size_param+1)/2;
/*		if ((idx<2) && (x[i]<=0.) ) MParams->Elt(i)=0.0;
		if ((idx<2) && (x[i]>=1.) ) MParams->Elt(i)=0.9999;
		if (idx==1) 
		{if (x[i-1]>=x[i]) 
			{ecart=x[i-1]-x[i];}
		}
		if ((idx==2) && (x[i]>=3.) ) MParams->Elt(i)=3.;
		if ((idx==2) && (x[i]<=-3.) ) MParams->Elt(i)=-3.;
		if ((idx==3) && (x[i]<0.) ) MParams->Elt(i)=0.001;
		if ((idx==3) && (x[i]>1.) ) MParams->Elt(i)=0.9999;
*/	  }
  }

  //Price Securities
  for (int il=0;il<ptr->its_cal_securitites.size();il++)
	{
		ICM_Pricer_Distrib_RandFactors Pricer; Pricer.Set(ptr->its_cal_securitites[il],
										ptr->its_cal_model,
										Params,ptr->GetAsOfDate());

		diff = -Pricer.Price(qCMPPRICE) - ptr->its_cal_CurrentPrices[il];
		diff /= (ptr->itsStrikeHigh->Elt(il)-ptr->itsStrikeLow->Elt(il)) * sizeindx * IndivNot;
		value += diff*diff;
	}
	
    ecart = 0.;
    value =sqrt(value)*exp(ecart);
	fx = value;

	return fx;
}                              

void ICM_Index_Correlation::GenImplyCorrelation( char* id, FILE* ficOut)
{
    FILE* fOut;
    char fOutName[200];

    int vectSize;
	int prev = 0;

    if ( ficOut == NULL )
    {
       ARM_GetViewFile(ARM_VIEW_FILE, id, fOutName);

       (void) unlink(fOutName);

       fOut = fopen(fOutName, "w");
    }
    else
    {
       fOut = ficOut;
    }

    
    if ( itsImpliedCorrelation == NULL )
    {
       vectSize = 0;
    }
    else
    {
       vectSize = itsImpliedCorrelation->GetSize();
    }

    // Now Print out Vect values in date format
    fprintf(fOut,"%d\n", vectSize);
	double temp = 0.;
    for (int i = 0;  i< vectSize; i++)
    {
        temp = itsImpliedCorrelation->Elt(i);
		fprintf(fOut,"%lf\n",temp);
    }

    if ( ficOut == NULL )
    {
       fclose(fOut);
    }
}

double ICM_Index_Correlation::HedgeImpact(ARM_Date& Start,ARM_Date& Maturity,const double& Leverage,const double& RefSpread)
	{
		ARM_Date date;
		//definition for the flat curve	
		const ICM_DefaultCurve* DefCurve = itsCreditIndex->GetForcedCurve();
		//ICM_DefaultCurveModel model(DefCurve,DefCurve->GetZeroCurve(),NULL,false);
		std::vector<const ICM_DefaultCurve*> DefaultCurves(1); 
		DefaultCurves[0] = DefCurve; // from index
		ICM_ModelMultiCurves model(DefaultCurves, DefCurve->GetZeroCurve(), NULL,false);
		

		//FILE* pFile = NULL;
		//pFile = fopen("C:\\temp\\CDSHedgeImpact.txt","w");
		ICM_Cds* cds = (ICM_Cds*)DefCurve->stdCDS(Start,Maturity,RefSpread/1.e4,its_cal_cdo_Notional*Leverage ,true);
		//cds->View("",pFile);
		// ICM_Pricer_Cds sltpricer(cds,&model);
		ICM_Parameters param;
		param.setParam("FLAT_COMPUTE",(bool) 1);
		//ICM_Pricer_Cds sltpricer; sltpricer.Set(cds,&model,ICM_Parameters(),model.GetAsOfDate());

		ICM_Pricer_CDSIndex sltpricer; sltpricer.Set(cds,&model,param,model.GetAsOfDate());

		double PV = sltpricer.Price(qCMPPRICE);//feelegpv
		//sltpricer.View("",pFile);
		//if(pFile) fclose(pFile);
		if (cds) delete cds;
		return (PV);
	}
void ICM_Index_Correlation::SetModelMultiCurves(ICM_ModelMultiCurves* ModelMultiCurves) 
	{
		if (itsModelMultiCurves) delete itsModelMultiCurves;
		if (ModelMultiCurves) itsModelMultiCurves = (ICM_ModelMultiCurves*) ModelMultiCurves->Clone();
	}
void ICM_Index_Correlation::SetCreditIndex(ICM_Credit_Index* index) 
	{
		if (itsCreditIndex)	delete itsCreditIndex;
		if (index) itsCreditIndex = (ICM_Credit_Index*) index->Clone();
	}
	// ----------------------------
	//	Copy of members data
	// ----------------------------
void ICM_Index_Correlation::BitwiseCopy(const ARM_Object* src)
	{
	    ICM_Index_Correlation* Correl = (ICM_Index_Correlation *) src;

		itsCalMethod = Correl->itsCalMethod;
		if (Correl->itsCreditIndex) itsCreditIndex = (ICM_Credit_Index*)Correl->itsCreditIndex->Clone();
		if (Correl->itsStrikeLow) itsStrikeLow = (ARM_Vector*)Correl->itsStrikeLow->Clone();
		if (Correl->itsStrikeHigh) itsStrikeHigh = (ARM_Vector*)Correl->itsStrikeHigh->Clone();
		if (Correl->itsMktBid) itsMktBid = (ARM_Vector*)Correl->itsMktBid->Clone();
		if (Correl->itsMktAsk) itsMktAsk = (ARM_Vector*)Correl->itsMktAsk->Clone();
		if (Correl->itsInitialCorrelation) itsInitialCorrelation = (ARM_Vector*)Correl->itsInitialCorrelation->Clone();
		if (Correl->itsImpliedCorrelation) itsImpliedCorrelation = (ARM_Vector*) Correl->itsImpliedCorrelation->Clone();
	}
// virtual 
double ICM_Index_Correlation::GetCorrelation(const std::string&  issuer1,
								  const std::string&  issuer2,
								  double maturity ,
								  double strike )
	{
		double Correlation = 0.;

		switch (itsCalMethod)
		{
		case qCAL_PWC_CORREL:
		case qCAL_PWL_CORREL:
			{
/*			ICM_CorrelFunction CF(1,, qCAL_INDEX_CORR_TYPE type, vector<double> &param)
			ICM_CorrelFunction.BetaCond
			std::vector< std::vector<double> >	itsImpliedParameters
*/			break;
			}
		case qCAL_BASE_CORRELATION:
		default:
			{Correlation = linInterpol(itsStrikeHigh,strike,itsImpliedCorrelation);}
		}
		
		return (Correlation);
	}
// inline 
void ICM_Index_Correlation::GetFactors(ICM_Parameters& p,ARM_Vector& v)
	{
		int i=0;
		v.Resize(0);

		switch (itsCalMethod)
		{
		case qCAL_PWL_CORREL:	
		case qCAL_PWC_CORREL:
			{
				ARM_Vector* vector = p.GetColVect("RF_PARAMS");
				v.Resize(vector->GetSize());
				for (i=0;i<vector->GetSize();i++) {v.Elt(i)=vector->Elt(i);}
				break;
			}
		default:;
		}
	}