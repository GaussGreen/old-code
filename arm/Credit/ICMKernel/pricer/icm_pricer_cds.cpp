#include "firsttoinc.h"
#include "ICMKernel/pricer/icm_pricer_cds.h"
#include "ICMKernel/glob/icm_constants.h"
#include "ICMKernel/mod/icm_defcurvemodel.h"
#include "ICMKernel/glob/icm_mktdatamng.h"
#include "ICMKernel/util/icm_polynoms.h"
#include "ICMKernel/inst/icm_credit_index.h"

#include "ICMKernel/util/icm_integrator.h"

// previously member of defcurve.. 
#define STEP_INT_FEE_LEG 7
#define STEP_INT_DEF_LEG 30

static const std::string CDSPricingMethod ("CDSPricingMethod"); 
static const std::string CDSIntegrationMethod ("INTEGRATION_METHOD"); 
static const std::string CDSIntegrationStep ("INTEGRATION_STEP"); 
static const double openDaysCoeff =1.4; 
static const int defIntegrationStep = 3; 
static const qIntegratorChoice  defIntegratorChoice = qGAUSS_LEGENDRE ;
static const int defPricingMethod = -1 ; // -1: default, 0: OldNonSummit, 1: Summit, 2: NonSummit

// integrated function for def+recovery 
//	friend function of ICM_Pricer_Cds
void objfun_cds(void* param_,double x,double&res)
{
	AddressVector* param = (AddressVector*)param_; 
	ICM_Pricer_Cds* pricer = (ICM_Pricer_Cds*)param->Get(0); 
	double * reco =  (double*)param->Get(1); 
	int * lagNotio =  (int*)param->Get(2); 
	int * lagReco=  (int*)param->Get(3); 
	double lambda = pricer->itsDefaultCurve->GetPWCIntensity(x) ;
	double sp = pricer->itsDefaultCurve->SurvivalProba(x) ;
	double dfNotio = 	pricer->itsDiscountCurve->DiscountPrice(x+*lagNotio/365.);
	double dfReco = 	pricer->itsDiscountCurve->DiscountPrice(x+*lagReco/365.);
	res = lambda*sp*(dfNotio-*reco*dfReco) / exp(-x*x*0.5) * SQRT2PI; 
}; 

// integrated function for accrued 
//	friend function of ICM_Pricer_Cds
void objfun_cds_acc(void* param_,double x,double&res)
{
	AddressVector* param = (AddressVector*)param_; 
	ICM_Pricer_Cds* pricer = (ICM_Pricer_Cds*)param->Get(0); 
	// ICM_Security* sec = ((ICM_Cds*)pricer->GetSecurity())->GetFeeLeg()->GetCreditInfos(); 
	// ARM_Date xDate(pricer->GetAsOfDate().GetJulian() + x*365.); 
	int * lagAcc =  (int*)param->Get(1); 
	double * yfStart=  (double*)param->Get(2); 
	double lambda = pricer->itsDefaultCurve->GetPWCIntensity(x) ;
	double sp = pricer->itsDefaultCurve->SurvivalProba(x) ;
	double dfNotio = 	pricer->itsDiscountCurve->DiscountPrice(x+*lagAcc/365.);
	res = (x-*yfStart)*lambda*sp*dfNotio / exp(-x*x*0.5) * SQRT2PI; 
}; 

void ICM_Pricer_Cds::Set(ARM_Security*sec,ARM_Object*mod,const ICM_Parameters&params,const ARM_Date&asof) 
{	
	ICM_Pricer::Set(sec, mod,params,&asof);
	LoadMarketData(mod,GetAsOfDate()); 
}

void ICM_Pricer_Cds::LoadMarketData(ARM_Object* object,const ARM_Date&AsOf)
{
	ICM_Cds* cds = (ICM_Cds*) GetSecurity();
	ICM_Security* _feeleg = (ICM_Security*) cds->GetFeeLeg()->GetCreditInfos();
	ICM_Leg* feefleg = (ICM_Leg*) cds->GetFeeLeg();
	ICM_Leg* defleg = (ICM_Leg*) cds->GetDefLeg();

	if (object->GetName()==ICM_MKTDATAMNG)
	{
		ICM_MktDataMng* Mng = dynamic_cast<ICM_MktDataMng*>(object);
		string zc_name = GetZeroCurveName(_feeleg->GetCcy(),AsOf);
		string cc_name = GetZeroCurveName(_feeleg->GetCcy(),AsOf);
		string dc_name = GetDefaultCurveName(feefleg->GetSingleName(),AsOf);
		itsDiscountCurve = dynamic_cast<ARM_ZeroCurve*>(Mng->find(zc_name));
		itsDefaultCurve = dynamic_cast<ICM_DefaultCurve*>(Mng->find(dc_name));
		itsCouponCurve = dynamic_cast<ARM_ZeroCurve*>(Mng->find(cc_name));
	}
	else
	{
		ICM_DefaultCurveModel* Model = dynamic_cast<ICM_DefaultCurveModel*>(object);
		itsDiscountCurve = dynamic_cast<ARM_ZeroCurve*>(Model->GetZeroCurve());
		itsDefaultCurve =  Model->GetDefaultCurve() ;
		itsCouponCurve = dynamic_cast<ARM_ZeroCurve*>(Model->GetZeroCurve());
	}
}
// ------------------------------------------------------------------------
// Calcul de l'accrued
// ------------------------------------------------------------------------
double ICM_Pricer_Cds::Accrued( )
{
	ICM_Cds* cds = (ICM_Cds*) GetSecurity();
	ARM_Date AsOf = GetAsOfDate();

	ICM_Security* security = cds->GetFeeLeg()->GetCreditInfos();
	security->ComputeYF(AsOf);

	ARM_ZeroCurve* ircurve = itsDiscountCurve;
	const ICM_DefaultCurve* defaultcurve = itsDefaultCurve;

	// 14514 int NumFlows = security->GetNumFlows();
	const ARM_Vector& AccStartDates = security->GetAccStartDates();
	const ARM_Vector& AccEndDates = security->GetAccEndDates();
	const ARM_Vector& PayDates = security->GetPayDates();
	int NumFlows = AccStartDates.size() ; 

	double Accrued = 0.,Coupon =0.,DP=0.,SP=0.;
	bool AccruedFlag = false;

	int iPeriodIndex = security->PeriodIndex(AsOf); 
	if (iPeriodIndex!=-1) 
	{
		Accrued =  security->AccruedCoupon(AsOf) ;
 
	}
	
	SetAccrued(Accrued); 
	return (Accrued);
}


//----------------------------------------------------------------------
// Calcul de la feeleg PV
//----------------------------------------------------------------------
double ICM_Pricer_Cds::FeeLegPV ()
{
	// ExecutionDate unsued. 
	if (GetFeeLegPriceFlg()) return GetFeeLegPrice();
	
	{
		// FeeLegPV_Old(); 
		// FeeLegPV_Summit();  
		// FeeLegPV_NonSummit();  
	} 

	long cdsPricingMethod (defPricingMethod); 
	if (!getParam(CDSPricingMethod,cdsPricingMethod,false)) cdsPricingMethod=defPricingMethod; // default case
	switch (cdsPricingMethod)
	{
		case -1: break; 
		case 0: return FeeLegPV_Old(); 
		case 1: return FeeLegPV_Summit();  
		case 2: return FeeLegPV_NonSummit();  
		default:
			ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_Pricer_Cds::FeeLegPV: Can't handle CDSPricingMethod "<<cdsPricingMethod); 
	}

	const ICM_DefaultCurve* defaultcurve = itsDefaultCurve;
	bool IsSummitPricer = defaultcurve->GetIsSummitCurve();
	if (IsSummitPricer) return FeeLegPV_Summit();
	return FeeLegPV_Old(); 
/** 
	if (GetFeeLegPriceFlg()) return GetFeeLegPrice();

	double PayLag = 0.,FeePV = 0.,FeeLegPV = 0.;
	double Coupon =0.,DP =0.,SP=0.,TmpFlow = 0.;
	int NumFlows = 0,indice = 0,int_size = 0;
	double t0 = 0., t1=0., tn=0.,Accrued =0.,ti=0.,ti_next=0.,irrate =0.;
	bool GoodPeriod = true;
	double Duration=0,Intermediate_Duration=0; 
	
	ICM_Cds* cds = (ICM_Cds*) GetSecurity();
	ARM_Date AsOf = GetAsOfDate();
	double fixedSpread = cds->GetCreditSpread(); 

	ICM_Leg* Feeleg = cds->GetFeeLeg();
	Feeleg->GetCreditInfos()->ComputeYF(AsOf);
	ICM_Security* security = Feeleg->GetCreditInfos(); 

	ARM_ZeroCurve* ircurve = itsDiscountCurve;

	// 14514 NumFlows = security->GetNumFlows();
	const ARM_Vector& AccStartDates = security->GetAccStartDates();
	const ARM_Vector& AccEndDates = security->GetAccEndDates();
	const ARM_Vector&Notionals = security->GetNotionals(); 
	NumFlows = AccStartDates.size() ;
	ARM_Vector* AccSch = NULL;
	MergeDates(&AccSch,AccStartDates,AccEndDates);

	const ARM_Vector& PayDates = security->GetPayDates();

	ARM_Vector* PeriodFeeLegPV = NULL;
	if (!GetFaster()) PeriodFeeLegPV = new ARM_Vector(NumFlows,0.);

	tn = AccEndDates.Elt(NumFlows-1);
	//ARM_Vector* Integration_Schedule = GenerateIntSch(MAX(AccStartDates.Elt(0),AsOf.GetJulian()),tn,STEP_INT_FEE_LEG );
	int IntFeeLeg = ( (cds->GetFeeLeg()->GetPaymentFreq()==K_ZEROCOUPON) ? (50000) : STEP_INT_FEE_LEG );
	ARM_Vector* Integration_Schedule = GenerateIntSch(MAX(AccStartDates.Elt(0),AsOf.GetJulian()),tn,IntFeeLeg);
	t0 = Integration_Schedule->Elt(0);
	t1 = Integration_Schedule->Elt(1);

	if (Feeleg->GetCreditIndex()) {Feeleg->GetCreditIndex()->ResetDefCurve();}
	Feeleg->CptCoupons(GetModel(),GetAsOfDate());

	// ----------------------------------------------------------------------------------------------
	// Full Period Pricing
	// ----------------------------------------------------------------------------------------------
	for (int i=0; i<NumFlows; i++)
	{
		DP =0.;SP=0.,Accrued =0.;GoodPeriod = true;Intermediate_Duration=0.;
		Coupon = security->FullCoupon((ARM_Date) AccStartDates.Elt(i)); //Full Coupon

		if ((AsOf.GetJulian()>AccStartDates.Elt(i)) &&  (AsOf.GetJulian()<AccEndDates.Elt(i))) 
		{
			Accrued = security->AccruedCoupon(AsOf); //Accrued coupon
			Coupon -= Accrued; //broken period
		}
		else if ((AsOf.GetJulian()>AccStartDates.Elt(i)) &&  (AsOf.GetJulian()>=AccEndDates.Elt(i))) 
		{ Coupon = 0.; GoodPeriod= false; }//Empty Period 

		if (GoodPeriod)
		{
		DP = ircurve->DiscountPrice((ARM_Date)PayDates.Elt(i));
		SP = defaultcurve->SurvivalProba((ARM_Date)AccEndDates.Elt(i));
		}

		Intermediate_Duration = DP * SP * Coupon/(Notionals.Elt(i)*fixedSpread/100.);
		Duration += Intermediate_Duration;
		FeeLegPV = Coupon * DP * SP + Accrued;
		if (!GetFaster()) PeriodFeeLegPV->Elt(i)=FeeLegPV;
		FeePV +=FeeLegPV; 
	}

	// Gestion de la premiere periode d'integration------------------------------------------------------
	SP = defaultcurve->SurvivalProba((ARM_Date)t0);
	// indice = GetIndice(AccSch,t1);  
	indice = locateIndex(AccSch,t1);

	TmpFlow = security->DiffAccrued((ARM_Date)AccSch->Elt(indice),(ARM_Date)t1)*
					ircurve->DiscountPrice((ARM_Date)(t1+PayLag))* SP;
	FeePV+=TmpFlow;
	int	ii;
	ii = security->PaymentIndex(AsOf);
	Duration+=TmpFlow/(Notionals.Elt(ii)*fixedSpread/100.) ;

	if (!GetFaster()) PeriodFeeLegPV->Elt(0)=TmpFlow + PeriodFeeLegPV->Elt(0);
	// End Gestion de la premiere periode d'integration--------------------------------------------------

	// Gestion de la derniere periode d'integration------------------------------------------------------
	// indice = GetIndice(AccSch,tn);
	indice = locateIndex(AccSch,tn);

	SP = defaultcurve->SurvivalProba((ARM_Date)tn);
	TmpFlow = security->DiffAccrued((ARM_Date)AccSch->Elt(indice),(ARM_Date)tn)*
					ircurve->DiscountPrice((ARM_Date)(tn+PayLag))*SP;
	FeePV-=TmpFlow;
	Duration -= TmpFlow/(Notionals.Elt(NumFlows-1)*fixedSpread/100.); 
	if (!GetFaster()) PeriodFeeLegPV->Elt(NumFlows-1)=-TmpFlow + PeriodFeeLegPV->Elt(NumFlows-1);
	// End Gestion de la derniere periode d'integration--------------------------------------------------

	int_size = Integration_Schedule->GetSize();
	ARM_Vector irrates(Integration_Schedule->GetSize()+(int)PayLag,0.);

	for (i=1; i<Integration_Schedule->GetSize()+(int)PayLag; i++)
	{
		if (i<Integration_Schedule->GetSize())
			irrate = ircurve->DiscountPrice((ARM_Date)Integration_Schedule->Elt(i));
		else
			irrate = ircurve->DiscountPrice((ARM_Date)(Integration_Schedule->Elt(int_size-1)+ (double)(i-int_size+1)));

		irrates.Elt(i)=irrate;
	}

	for (i=1; i<Integration_Schedule->GetSize()-1; i++)
	{
		ti_next = Integration_Schedule->Elt(i+1);
		// indice = GetIndice(AccSch,ti_next);
		indice = locateIndex(AccSch,ti_next);

		DP = irrates.Elt(i+1+(int)PayLag);

		TmpFlow = security->DiffAccrued((ARM_Date)AccSch->Elt(indice),(ARM_Date)ti_next)*DP;

		ti = Integration_Schedule->Elt(i);
		// indice = GetIndice(AccSch,ti);
		indice = locateIndex(AccSch,ti);

		DP = irrates.Elt(i+(int)PayLag);

		TmpFlow -= security->DiffAccrued((ARM_Date)AccSch->Elt(indice),(ARM_Date)ti)*DP;
		
		SP = defaultcurve->SurvivalProba((ARM_Date)ti);
		TmpFlow *= SP;

		if (!GetFaster()) 
			PeriodFeeLegPV->Elt(indice)=TmpFlow + PeriodFeeLegPV->Elt(indice);

		FeePV +=TmpFlow; 
		Duration += TmpFlow / (Notionals.Elt(indice)*fixedSpread/100.) ;
	}

	if (!GetFaster()) 
		security->SetPeriodFeeLegPV(*PeriodFeeLegPV);

	SetFeeLegPrice(FeePV);
	SetDuration(fabs(Duration)) ;

	if (Integration_Schedule) delete Integration_Schedule;
	if (AccSch) delete AccSch;

	return (FeePV);
	**/ 
}


//----------------------------------------------------------------------
// Calcul de la feeleg PV / Cas Non Summit
//----------------------------------------------------------------------
double ICM_Pricer_Cds::FeeLegPV_NonSummit ()
{
	// ExecutionDate unsued. 
	// if (GetFeeLegPriceFlg()) return GetFeeLegPrice();
	ICMQUANTIFYER("FeeLegPV_NonSummit "); 
	// GETTING SOME DATA

	// Security
	ICM_Cds* cds = (ICM_Cds*) GetSecurity();
	ARM_Date AsOf = GetAsOfDate();
	ICM_Leg* Feeleg = cds->GetFeeLeg();
	// Compute schedule Year Fractions
	Feeleg->GetCreditInfos()->ComputeYF(AsOf);
	ICM_Security* security = Feeleg->GetCreditInfos(); 

	// Zero Curve
	ARM_ZeroCurve* ircurve = itsDiscountCurve;
	const ICM_DefaultCurve* defaultcurve = itsDefaultCurve;

	// The Default Lag
	double DefPayLag = cds->GetCreditLag();

	// Flows
	const ARM_Vector& AccStartDates = security->GetAccStartDates();
	const ARM_Vector& AccEndDates = security->GetAccEndDates();
	int	NumFlows = AccStartDates.size() ;
	// ---------------------------------------------------------------
	// PAYMENT TYPE
	qPAYMENT_PREMIUM_LEG	ThePremiumLegPaymentType;
	ThePremiumLegPaymentType	=	cds->GetFeeLeg()->GetAccruedOnDefault();
	// ---------------------------------------------------------------

	
	// ---------------------------------------------------------------
	// Start and End dates
	const ARM_Vector& YFAccStartDates = security->GetYFAccStartDates();
	const ARM_Vector& YFAccEndDates = security->GetYFAccEndDates();
	const ARM_Vector& YFPayDates = security->GetYFPayDates();
	const ARM_Vector& Notionals = security->GetNotionals();
	const ARM_Vector& YFInterestDays = security->GetYFInterestDays();
	ARM_Vector DefProb ;  
	ARM_Vector Discount; 

	if (!GetFaster()) 
	{	DefProb.Resize(Notionals.GetSize());
		Discount.Resize(Notionals.GetSize());	
	}

	double YFFirstAccStart	= MAX(YFAccStartDates.Elt(0),0.);
	double YFAccEndDate		= YFAccEndDates[NumFlows-1];
	// ---------------------------------------------------------------

	ICM_Poly	Poly;

	double t1=0.,t2=0.,Lambda=0.,integral=0.;
	double	DP, SP, Coupon, FeeLegPV;
	double	CumPremium = 0.0;
	double Duration = 0, Im_Duration = 0,DeltaYF=0.;; 
	double fixedSpread  = cds->GetCreditSpread(); 

	
	if (Feeleg->GetCreditIndex()) {Feeleg->GetCreditIndex()->ResetDefCurve();}
	Feeleg->CptCoupons(GetModel(),GetAsOfDate());

	//  We iterate over all the flows that are to be paid. 
	// 
	double CumPremium2=0; 
	double Duration2=0; 
	int iPayment ; // 
	int iPayment0 = security->PaymentIndex(AsOf) ;	// index of first flow to be paid (-1 if no flows) 
	int iAccrued = security->PeriodIndex(AsOf) ;	// index of the period including AsOf (-1 if no period) 
	if (iPayment0==-1) 
	{
		ICMLOG("WARNING: Leg has no payments... "); 
		SetFeeLegPrice(0);
		SetDuration(0); 
		return 0;
	}
	for (iPayment=iPayment0;iPayment<NumFlows;iPayment++) 
	{
		Coupon = security->FullCoupon(iPayment); 
		DP = ircurve->DiscountPrice(YFPayDates[iPayment]);
		switch (ThePremiumLegPaymentType)
		{
			case qCONTINUE_TO_MATURITY:
				SP	=	1.0;
				break;

			case qACCRUED_NOT_SETTLED:
			case qACCRUED_SETTLED:
				SP = defaultcurve->SurvivalProba(YFAccEndDates.Elt(iPayment));
				break;

			case qCOMPLETE_CURRENT_PERIOD:
				SP = defaultcurve->SurvivalProba(YFAccStartDates.Elt(iPayment));
				break;
		}
		double DeltaYF= YFInterestDays.Elt(iPayment) ;
		
		Im_Duration = DeltaYF * DP * SP;
		Duration2 += Im_Duration;
		//Duration based on NetPrice
		if (iPayment!=iAccrued) 
			FeeLegPV = Coupon  / DeltaYF * Im_Duration  ;
		else 
		{
			double accrued = security->AccruedCoupon(AsOf) ; 
			double DeltaAcc = DeltaYF - CountYears(security->GetAccrualBasis(),AsOf.GetJulian(),AccEndDates.Elt(iPayment));
			if (ThePremiumLegPaymentType == qACCRUED_SETTLED)
			{
				FeeLegPV = (Coupon-accrued) / DeltaYF * Im_Duration  + accrued * DP;
				Duration2 += DeltaAcc * (DP-SP*DP);	
			}
			else
				FeeLegPV = Coupon  / DeltaYF * Im_Duration ;
			
			Duration2 += - DeltaAcc;	
		}
		CumPremium2 += FeeLegPV; 
		// view only 
		if (!GetFaster())
		{
			DefProb.Elt(iPayment)= 1. - SP;
			Discount.Elt(iPayment)= DP;
		}
	}
	CumPremium=CumPremium2; 
	Duration=Duration2 ;
	
	if (ThePremiumLegPaymentType == qCONTINUE_TO_MATURITY) 
	{
		CumPremium	*=	defaultcurve->SurvivalProba(YFFirstAccStart);
		Duration	*=	defaultcurve->SurvivalProba(YFFirstAccStart);
	}


	// ---------------------------------------------------------------
	// ACCRUED
	
	double	CumPremiumAccrued, FeePV;
	CumPremiumAccrued	=	0.0;
	FeePV=CumPremium ;
	// ---------------------------------------------------------------
	if (ThePremiumLegPaymentType == qACCRUED_SETTLED)
	{
		ICM_Integrator integrator; 
		//
		qIntegratorChoice integratorChoice(defIntegratorChoice); 
		long lIntegrationChoice; 
		if (getParam(CDSIntegrationMethod,lIntegrationChoice,false)==true) 
			integratorChoice = (qIntegratorChoice)lIntegrationChoice ;
		integrator.SetIntegrationType(integratorChoice); 
		long integrationStep; 
		if (!getParam(CDSIntegrationStep,integrationStep,false)) 
			integrationStep= defIntegrationStep; 
		integrator.SetIntegrationStep(integrationStep); 
		AddressVector param(4); 

		int lagAcc = cds->GetCreditLag() * openDaysCoeff; 
		param.Set(0,this); 
		param.Set(1,&lagAcc); 
		double yfStart; 
		param.Set(2,&yfStart); 
		int nointegcoeff; 
		param.Set(3,&nointegcoeff); 
		for(int iPayment=iPayment0;iPayment<security->GetYFAccStartDates().size();iPayment++)
		{
			yfStart = security->GetYFAccStartDates().Elt(iPayment); 
			double yfEnd = security->GetYFAccEndDates().Elt(iPayment);  
			double res=0; 
			integrator.Integrate(yfStart,yfEnd,&objfun_cds_acc,&param,res); 
			// double notional = security->GetNotionals().Elt(iPayment) ; 
			double coupon = security->FullCoupon(iPayment); 
			CumPremiumAccrued+=res*coupon/(yfEnd-yfStart) ; // notional is already inside the coupon
		}
	}
	FeePV	+= CumPremiumAccrued;


	//on applique ce correctif pour le momment
	if (fixedSpread) Duration = 100.*(FeePV - security->AccruedCoupon(AsOf))/(Notionals.Elt(0)*fixedSpread);

	SetFeeLegPrice(FeePV);
	SetDuration(fabs(Duration)); 

	if (!GetFaster()) 
	{
		security->SetDefaultProbability(DefProb) ; 
		security->SetDiscountRate(Discount); 
	}

	return (FeePV);
}


//----------------------------------------------------------------------
// Calcul de la feeleg PV / Cas Non Summit
//----------------------------------------------------------------------
double ICM_Pricer_Cds::FeeLegPV_Old()
{
	
	ICMQUANTIFYER("FeeLegPV_Old"); 
	//if (GetFeeLegPriceFlg()) return GetFeeLegPrice();
	const ICM_DefaultCurve* defaultcurve = itsDefaultCurve;

	double PayLag = 0.,FeePV = 0.,FeeLegPV = 0.;
	double Coupon =0.,DP =0.,SP=0.,TmpFlow = 0.;
	int NumFlows = 0,indice = 0,int_size = 0;
	double t0 = 0., t1=0., tn=0.,Accrued =0.,ti=0.,ti_next=0.,irrate =0.;
	bool GoodPeriod = true;
	double Duration=0,Intermediate_Duration=0; 
	
	ICM_Cds* cds = (ICM_Cds*) GetSecurity();
	ARM_Date AsOf = GetAsOfDate();
	double fixedSpread = cds->GetCreditSpread(); 

	ICM_Leg* Feeleg = cds->GetFeeLeg();
	Feeleg->GetCreditInfos()->ComputeYF(AsOf);
	ICM_Security* security = Feeleg->GetCreditInfos(); 

	ARM_ZeroCurve* ircurve = itsDiscountCurve;

	// 14514 NumFlows = security->GetNumFlows();
	const ARM_Vector& AccStartDates = security->GetAccStartDates();
	const ARM_Vector& AccEndDates = security->GetAccEndDates();
	const ARM_Vector&Notionals = security->GetNotionals(); 
	NumFlows = AccStartDates.size() ;
	ARM_Vector* AccSch = NULL;
	MergeDates(&AccSch,AccStartDates,AccEndDates);

	const ARM_Vector& PayDates = security->GetPayDates();

	ARM_Vector* PeriodFeeLegPV = NULL;
	if (!GetFaster()) PeriodFeeLegPV = new ARM_Vector(NumFlows,0.);

	tn = AccEndDates.Elt(NumFlows-1);
	//ARM_Vector* Integration_Schedule = GenerateIntSch(MAX(AccStartDates.Elt(0),AsOf.GetJulian()),tn,STEP_INT_FEE_LEG );
	int IntFeeLeg = ( (cds->GetFeeLeg()->GetPaymentFreq()==K_ZEROCOUPON) ? (50000) : STEP_INT_FEE_LEG );
	ARM_Vector* Integration_Schedule = GenerateIntSch(MAX(AccStartDates.Elt(0),AsOf.GetJulian()),tn,IntFeeLeg);
	t0 = Integration_Schedule->Elt(0);
	t1 = Integration_Schedule->Elt(1);

	if (Feeleg->GetCreditIndex()) {Feeleg->GetCreditIndex()->ResetDefCurve();}
	Feeleg->CptCoupons(GetModel(),GetAsOfDate());

	// ----------------------------------------------------------------------------------------------
	// Full Period Pricing
	// ----------------------------------------------------------------------------------------------
	for (int i=0; i<NumFlows; i++)
	{
		DP =0.;SP=0.,Accrued =0.;GoodPeriod = true;Intermediate_Duration=0.;
		Coupon = security->FullCoupon((ARM_Date) AccStartDates.Elt(i)); //Full Coupon

		if ((AsOf.GetJulian()>AccStartDates.Elt(i)) &&  (AsOf.GetJulian()<AccEndDates.Elt(i))) 
		{
			Accrued = security->AccruedCoupon(AsOf); //Accrued coupon
			Coupon -= Accrued; //broken period
		}
		else if ((AsOf.GetJulian()>AccStartDates.Elt(i)) &&  (AsOf.GetJulian()>=AccEndDates.Elt(i))) 
		{ Coupon = 0.; GoodPeriod= false; }//Empty Period 

		if (GoodPeriod)
		{
		DP = ircurve->DiscountPrice((ARM_Date)PayDates.Elt(i));
		SP = defaultcurve->SurvivalProba((ARM_Date)AccEndDates.Elt(i));
		}

		Intermediate_Duration = DP * SP * Coupon/(Notionals.Elt(i)*fixedSpread/100.);
		Duration += Intermediate_Duration;
		FeeLegPV = Coupon * DP * SP + Accrued;
		if (!GetFaster()) PeriodFeeLegPV->Elt(i)=FeeLegPV;
		FeePV +=FeeLegPV; 
	}

	// Gestion de la premiere periode d'integration------------------------------------------------------
	SP = defaultcurve->SurvivalProba((ARM_Date)t0);
	/** indice = GetIndice(AccSch,t1); **/
	indice = locateIndex(AccSch,t1);

	TmpFlow = security->DiffAccrued((ARM_Date)AccSch->Elt(indice),(ARM_Date)t1)*
					ircurve->DiscountPrice((ARM_Date)(t1+PayLag))* SP;
	FeePV+=TmpFlow;
	int	ii;
	ii = security->PaymentIndex(AsOf);
	Duration+=TmpFlow/(Notionals.Elt(ii)*fixedSpread/100.) ;

	if (!GetFaster()) PeriodFeeLegPV->Elt(0)=TmpFlow + PeriodFeeLegPV->Elt(0);
	// End Gestion de la premiere periode d'integration--------------------------------------------------

	// Gestion de la derniere periode d'integration------------------------------------------------------
	// indice = GetIndice(AccSch,tn);
	indice = locateIndex(AccSch,tn);

	SP = defaultcurve->SurvivalProba((ARM_Date)tn);
	TmpFlow = security->DiffAccrued((ARM_Date)AccSch->Elt(indice),(ARM_Date)tn)*
					ircurve->DiscountPrice((ARM_Date)(tn+PayLag))*SP;
	FeePV-=TmpFlow;
	Duration -= TmpFlow/(Notionals.Elt(NumFlows-1)*fixedSpread/100.); 
	if (!GetFaster()) PeriodFeeLegPV->Elt(NumFlows-1)=-TmpFlow + PeriodFeeLegPV->Elt(NumFlows-1);
	// End Gestion de la derniere periode d'integration--------------------------------------------------

	int_size = Integration_Schedule->GetSize();
	ARM_Vector irrates(Integration_Schedule->GetSize()+(int)PayLag,0.);

	for (i=1; i<Integration_Schedule->GetSize()+(int)PayLag; i++)
	{
		if (i<Integration_Schedule->GetSize())
			irrate = ircurve->DiscountPrice((ARM_Date)Integration_Schedule->Elt(i));
		else
			irrate = ircurve->DiscountPrice((ARM_Date)(Integration_Schedule->Elt(int_size-1)+ (double)(i-int_size+1)));

		irrates.Elt(i)=irrate;
	}

	for (i=1; i<Integration_Schedule->GetSize()-1; i++)
	{
		ti_next = Integration_Schedule->Elt(i+1);
		// indice = GetIndice(AccSch,ti_next);
		indice = locateIndex(AccSch,ti_next);

		DP = irrates.Elt(i+1+(int)PayLag);

		TmpFlow = security->DiffAccrued((ARM_Date)AccSch->Elt(indice),(ARM_Date)ti_next)*DP;

		ti = Integration_Schedule->Elt(i);
		// indice = GetIndice(AccSch,ti);
		indice = locateIndex(AccSch,ti);

		DP = irrates.Elt(i+(int)PayLag);

		TmpFlow -= security->DiffAccrued((ARM_Date)AccSch->Elt(indice),(ARM_Date)ti)*DP;
		
		SP = defaultcurve->SurvivalProba((ARM_Date)ti);
		TmpFlow *= SP;

		if (!GetFaster()) 
			PeriodFeeLegPV->Elt(indice)=TmpFlow + PeriodFeeLegPV->Elt(indice);

		FeePV +=TmpFlow; 
		Duration += TmpFlow / (Notionals.Elt(indice)*fixedSpread/100.) ;
	}

	if (!GetFaster()) 
		security->SetPeriodFeeLegPV(*PeriodFeeLegPV);

	SetFeeLegPrice(FeePV);
	SetDuration(fabs(Duration)) ;

	if (Integration_Schedule) delete Integration_Schedule;
	if (AccSch) delete AccSch;

	return (FeePV);
}

//----------------------------------------------------------------------
// Calcul de la defleg PV
//----------------------------------------------------------------------
double ICM_Pricer_Cds::DefLegPV ()
{
	if (GetDefLegPriceFlg()) return GetDefLegPrice();

	// ExecutionDate unsued. 
	long cdsPricingMethod (defPricingMethod); 
	if (!getParam(CDSPricingMethod,cdsPricingMethod,false)) cdsPricingMethod=defPricingMethod; // default case
	
	{
		// DefLegPV_Old(); 
		// DefLegPV_Summit();  
		// DefLegPV_NonSummit();  
	}
	
	
	switch (cdsPricingMethod)
	{
		case -1: break; 
		case 0: return DefLegPV_Old(); 
		case 1: return DefLegPV_Summit();  
		case 2: return DefLegPV_NonSummit();  
		default:
			ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_Pricer_Cds::DefLegPV: Can't handle CDSPricingMethod "<<cdsPricingMethod); 
	}

	const ICM_DefaultCurve* defaultcurve = itsDefaultCurve;
	bool IsSummitPricer = defaultcurve->GetIsSummitCurve();
	if (IsSummitPricer) return DefLegPV_Summit();
	return DefLegPV_Old(); 
/** 
	// ExecutionDate unsued. 
	const ICM_DefaultCurve* defaultcurve = itsDefaultCurve;
	bool IsSummitPricer = defaultcurve->GetIsSummitCurve();

	if (IsSummitPricer) return DefLegPV_Summit();
	if (GetDefLegPriceFlg()) return GetDefLegPrice();

	ICM_Cds* cds = (ICM_Cds*) GetSecurity();
	ARM_Date AsOf = GetAsOfDate();
	
	double DefPayLag = cds->GetCreditLag();
 
	ICM_Leg* Defleg = cds->GetDefLeg();
	Defleg->GetCreditInfos()->ComputeYF(AsOf);
	ICM_Security* security = Defleg->GetCreditInfos(); 

	bool isConstantNotional(true) ; 
	double constantNotionalValue(0); 
	if ( Defleg->GetAmount() && Defleg->GetAmount()->size()>1) isConstantNotional=false; 
	else {
		int i = security->PeriodIndex(AsOf); 
		if (i<0) i=0; 
		constantNotionalValue= security->GetNotionals().Elt(i) ; 
	}

	ARM_ZeroCurve* ircurve = itsDiscountCurve;

	
	const ARM_Vector& AccStartDates = security->GetAccStartDates();
	const ARM_Vector& AccEndDates = security->GetAccEndDates();
	// 14514 int NumFlows = security->GetNumFlows();
	int NumFlows = AccStartDates.size(); 

	ARM_Vector* AccSch = NULL;
	MergeDates(&AccSch,AccStartDates,AccEndDates);

	ARM_Vector* PeriodDefLegPV = NULL;
	if (!GetFaster()) PeriodDefLegPV = new ARM_Vector(NumFlows,0.);

	double tn = AccEndDates.Elt(NumFlows-1);
	//ARM_Vector* Integration_Schedule = GenerateIntSch(MAX(AccStartDates.Elt(0),AsOf.GetJulian()),tn,STEP_INT_DEF_LEG );
	int IntDefLeg = ( (cds->GetDefLeg()->GetPaymentFreq()==K_ZEROCOUPON) ? (50000) : STEP_INT_DEF_LEG );
	ARM_Vector* Integration_Schedule = GenerateIntSch(MAX(AccStartDates.Elt(0),AsOf.GetJulian()),tn,IntDefLeg);
	double t1 = Integration_Schedule->Elt(1);
	double t0 = Integration_Schedule->Elt(0);

	double DefPV = 0.,LR = 1.;
	double DP1 =0.,DP2 =0.,SP1 =0.,SP2 =0.;
	double TmpFlow = 0., irrate = 0.;// notional = 0.
	int indice = 0,i = 0,int_size =0;
	double summitCoef = 1.4;

	if (Defleg->GetBinaryFlg())
		LR= 1.- Defleg->GetBinary();  
	else
		LR= 1.- defaultcurve->GetRecovery();  

	const ARM_Vector&notionals = security->GetNotionals();
	int_size = Integration_Schedule->GetSize();

	ARM_Vector irrates(Integration_Schedule->GetSize(),0.);

	for (i=1; i<Integration_Schedule->GetSize(); i++)
	{
		irrate = ircurve->DiscountPrice((ARM_Date)(Integration_Schedule->Elt(i)+summitCoef*DefPayLag));
		irrates.Elt(i)=irrate;
	}

	// notional = notionals->Elt(0);

	for (i=1; i<Integration_Schedule->GetSize()-1; i++)
	{
		DP1 = irrates.Elt(i);
		DP2 = irrates.Elt(i+1);
		SP1 = defaultcurve->SurvivalProba((ARM_Date)Integration_Schedule->Elt(i));

		if ( !GetFaster() || !isConstantNotional ) 
		{
		TmpFlow = LR * (DP2 - DP1)* SP1;

		// indice = GetIndice(AccSch,Integration_Schedule->Elt(i));
		indice = locateIndex(AccSch,Integration_Schedule->Elt(i));

		if (indice<0) indice = 0;
		// if (indice<0) indice = security->PeriodIndex(AsOf) ;

		TmpFlow *= 	notionals.Elt(indice);

		DefPV +=TmpFlow; 
		PeriodDefLegPV->Elt(indice)= TmpFlow+PeriodDefLegPV->Elt(indice);
		}
		else
		{
		TmpFlow = LR * (DP2 - DP1)* SP1 * constantNotionalValue ; 
		DefPV +=TmpFlow; 
		}
	}

	DP1 = ircurve->DiscountPrice((ARM_Date)(Integration_Schedule->Elt(1)+summitCoef*DefPayLag));
	DP2 = ircurve->DiscountPrice((ARM_Date)(tn + summitCoef*DefPayLag));
	SP1 = defaultcurve->SurvivalProba((ARM_Date)Integration_Schedule->Elt(0));
	SP2 = defaultcurve->SurvivalProba((ARM_Date)tn);

	indice = security->PeriodIndex((ARM_Date)Integration_Schedule->Elt(1));
	if (indice<0) indice = 0;

	TmpFlow = LR*DP1*SP1*notionals.Elt(indice);
	if (!GetFaster()) PeriodDefLegPV->Elt(indice)=TmpFlow+PeriodDefLegPV->Elt(indice);

	DefPV += TmpFlow;

	indice = security->PeriodIndex((ARM_Date)tn);
	TmpFlow = -LR*DP2*SP2*notionals.Elt(indice);
	if (!GetFaster()) PeriodDefLegPV->Elt(NumFlows-1)=TmpFlow+PeriodDefLegPV->Elt(NumFlows-1);

	DefPV += TmpFlow;

	if (!GetFaster()) security->SetPeriodDefLegPV(*PeriodDefLegPV);

	if (defaultcurve->GetIsNameInDefault()) DefPV	= LR * constantNotionalValue ; 


	SetDefLegPrice(DefPV);

	if (Integration_Schedule) delete Integration_Schedule;
	if (AccSch) delete AccSch;

	return (DefPV);
	**/ 
}

//----------------------------------------------------------------------
// Calcul de la defleg PV / CAS NON SUMMIT. 
//----------------------------------------------------------------------
double ICM_Pricer_Cds::DefLegPV_NonSummit ()
{
	ICMQUANTIFYER("DefLegPV_NonSummit "); 
	// 
	ICM_Cds* cds = (ICM_Cds*) GetSecurity();
	ARM_Date AsOf = GetAsOfDate();
	//Default Leg
	ICM_Leg* Defleg = cds->GetDefLeg();
	Defleg->GetCreditInfos()->ComputeYF(AsOf);
	ICM_Security* security = Defleg->GetCreditInfos(); 
	// 
	int iPayment0 = security->PaymentIndex(AsOf) ;	// index of first flow to be paid (-1 if no flows) 
	// int iAccrued = security->PeriodIndex(AsOf) ;	// index of the period including AsOf (-1 if no period) 
	if (iPayment0==-1) 
	{
		ICMLOG("WARNING: Leg has no payments... "); 
		SetDefLegPrice(1);
		return 0; 
	}
	double pv=0; 
	// 
	ICM_Integrator integrator; 
	//
	qIntegratorChoice integratorChoice(qGAUSS_LEGENDRE); long lIntegrationChoice=0;
	if (getParam(CDSIntegrationMethod,lIntegrationChoice,false)==true) 
		integratorChoice = (qIntegratorChoice)lIntegrationChoice ;
	integrator.SetIntegrationType(integratorChoice); 
	long integrationStep ;
	if (!getParam(CDSIntegrationStep,integrationStep,false)) 
		integrationStep= defIntegrationStep; 
	integrator.SetIntegrationStep(integrationStep); 
	AddressVector param(5); 

	double reco=0 ;
	if (Defleg->GetBinaryFlg()) reco = Defleg->GetBinary();  
	else reco = itsDefaultCurve->GetRecovery();  

	int lagNotio = cds->GetCreditLag() * openDaysCoeff; 
	int lagReco = 3 * openDaysCoeff;
	param.Set(0,this); 
	param.Set(1,&reco); 
	param.Set(2,&lagNotio); 
	param.Set(3,&lagReco); 
	int nointegcoeff; 
	param.Set(4,&nointegcoeff); 
	for(int iPayment=iPayment0;iPayment<security->GetYFAccStartDates().size();iPayment++)
	{
		double yfStart = security->GetYFAccStartDates().Elt(iPayment); 
		double yfEnd = security->GetYFAccEndDates().Elt(iPayment);  
		double res=0; 
		integrator.Integrate(yfStart,yfEnd,&objfun_cds,&param,res); 
		double notional = security->GetNotionals().Elt(iPayment) ; 
		pv+=res*notional; 
	}
	// 
	SetDefLegPrice(pv);
	return pv ;
}

//----------------------------------------------------------------------
// Calcul de la defleg PV / CAS NON SUMMIT. 
//----------------------------------------------------------------------
double ICM_Pricer_Cds::DefLegPV_Old ()
{
	ICMQUANTIFYER("DefLegPV_Old"); 
	// ExecutionDate unsued. 
	
	// bool IsSummitPricer = defaultcurve->GetIsSummitCurve();

	// if (IsSummitPricer) return DefLegPV_Summit();
	

	const ICM_DefaultCurve* defaultcurve = itsDefaultCurve;
	ICM_Cds* cds = (ICM_Cds*) GetSecurity();
	ARM_Date AsOf = GetAsOfDate();
	
	double DefPayLag = cds->GetCreditLag();
 
	ICM_Leg* Defleg = cds->GetDefLeg();
	Defleg->GetCreditInfos()->ComputeYF(AsOf);
	ICM_Security* security = Defleg->GetCreditInfos(); 

	bool isConstantNotional(true) ; 
	double constantNotionalValue(0); 
	if ( Defleg->GetAmount() && Defleg->GetAmount()->size()>1) isConstantNotional=false; 
	else {
		int i = security->PeriodIndex(AsOf); 
		if (i<0) i=0; 
		constantNotionalValue= security->GetNotionals().Elt(i) ; 
	}

	ARM_ZeroCurve* ircurve = itsDiscountCurve;

	
	const ARM_Vector& AccStartDates = security->GetAccStartDates();
	const ARM_Vector& AccEndDates = security->GetAccEndDates();
	// 14514 int NumFlows = security->GetNumFlows();
	int NumFlows = AccStartDates.size(); 

	ARM_Vector* AccSch = NULL;
	MergeDates(&AccSch,AccStartDates,AccEndDates);

	ARM_Vector* PeriodDefLegPV = NULL;
	if (!GetFaster()) PeriodDefLegPV = new ARM_Vector(NumFlows,0.);

	double tn = AccEndDates.Elt(NumFlows-1);
	//ARM_Vector* Integration_Schedule = GenerateIntSch(MAX(AccStartDates.Elt(0),AsOf.GetJulian()),tn,STEP_INT_DEF_LEG );
	int IntDefLeg = ( (cds->GetDefLeg()->GetPaymentFreq()==K_ZEROCOUPON) ? (50000) : STEP_INT_DEF_LEG );
	ARM_Vector* Integration_Schedule = GenerateIntSch(MAX(AccStartDates.Elt(0),AsOf.GetJulian()),tn,IntDefLeg);
	double t1 = Integration_Schedule->Elt(1);
	double t0 = Integration_Schedule->Elt(0);

	double DefPV = 0.,LR = 1.;
	double DP1 =0.,DP2 =0.,SP1 =0.,SP2 =0.;
	double TmpFlow = 0., irrate = 0.;// notional = 0.
	int indice = 0,i = 0,int_size =0;
	double summitCoef = 1.4;

	if (Defleg->GetBinaryFlg())
		LR= 1.- Defleg->GetBinary();  
	else
		LR= 1.- defaultcurve->GetRecovery();  

	const ARM_Vector&notionals = security->GetNotionals();
	int_size = Integration_Schedule->GetSize();

	ARM_Vector irrates(Integration_Schedule->GetSize(),0.);

	for (i=1; i<Integration_Schedule->GetSize(); i++)
	{
		irrate = ircurve->DiscountPrice((ARM_Date)(Integration_Schedule->Elt(i)+summitCoef*DefPayLag));
		irrates.Elt(i)=irrate;
	}

	// notional = notionals->Elt(0);

	for (i=1; i<Integration_Schedule->GetSize()-1; i++)
	{
		DP1 = irrates.Elt(i);
		DP2 = irrates.Elt(i+1);
		SP1 = defaultcurve->SurvivalProba((ARM_Date)Integration_Schedule->Elt(i));

		if ( !GetFaster() || !isConstantNotional ) 
		{
		TmpFlow = LR * (DP2 - DP1)* SP1;

		// indice = GetIndice(AccSch,Integration_Schedule->Elt(i));
		indice = locateIndex(AccSch,Integration_Schedule->Elt(i));

		if (indice<0) indice = 0;
		// if (indice<0) indice = security->PeriodIndex(AsOf) ;

		TmpFlow *= 	notionals.Elt(indice);

		DefPV +=TmpFlow; 
		PeriodDefLegPV->Elt(indice)= TmpFlow+PeriodDefLegPV->Elt(indice);
		}
		else
		{
		TmpFlow = LR * (DP2 - DP1)* SP1 * constantNotionalValue ; 
		DefPV +=TmpFlow; 
		}
	}

	DP1 = ircurve->DiscountPrice((ARM_Date)(Integration_Schedule->Elt(1)+summitCoef*DefPayLag));
	DP2 = ircurve->DiscountPrice((ARM_Date)(tn + summitCoef*DefPayLag));
	SP1 = defaultcurve->SurvivalProba((ARM_Date)Integration_Schedule->Elt(0));
	SP2 = defaultcurve->SurvivalProba((ARM_Date)tn);

	indice = security->PeriodIndex((ARM_Date)Integration_Schedule->Elt(1));
	if (indice<0) indice = 0;

	TmpFlow = LR*DP1*SP1*notionals.Elt(indice);
	if (!GetFaster()) PeriodDefLegPV->Elt(indice)=TmpFlow+PeriodDefLegPV->Elt(indice);

	DefPV += TmpFlow;

	indice = security->PeriodIndex((ARM_Date)tn);
	TmpFlow = -LR*DP2*SP2*notionals.Elt(indice);
	if (!GetFaster()) PeriodDefLegPV->Elt(NumFlows-1)=TmpFlow+PeriodDefLegPV->Elt(NumFlows-1);

	DefPV += TmpFlow;

	if (!GetFaster()) security->SetPeriodDefLegPV(*PeriodDefLegPV);

	if (defaultcurve->GetIsNameInDefault()) DefPV	= LR * constantNotionalValue ; 


	SetDefLegPrice(DefPV);

	if (Integration_Schedule) delete Integration_Schedule;
	if (AccSch) delete AccSch;

	return (DefPV);
}

//----------------------------------------------------------------------------
//	  Compute Implied Credit Spread
//----------------------------------------------------------------------------
double ICM_Pricer_Cds::ComputeSpread(const double& MtM )
{
	double Result=0.,initialnotional=0.;

	if (GetSpreadFlg() && MtM == 0 ) return GetSpread();

	ICM_Cds* cds = dynamic_cast<ICM_Cds*>(GetSecurity());
	
	if (cds)
	{ if (cds->GetSecurityType() == qCM_CDS) return ComputeImpliedPartRate(); }

	ComputePrice(qCMPPRICE);

	if (cds)
	{ 
		// 14514 initialnotional = cds->GetFeeLeg()->GetCreditInfos()->GetInitialNotional();
		initialnotional = cds->GetFeeLeg()->GetCreditInfos()->GetNotionals().Elt(0);
	}
	else
	{	
		ICM_Leg* leg = dynamic_cast<ICM_Leg*>(GetSecurity());
		// 14514 initialnotional = leg->GetCreditInfos()->GetInitialNotional();
		initialnotional = leg->GetCreditInfos()->GetNotionals().Elt(0);
	}

	Result=(MtM + Price(qCMPDEFLEGPV))/(Price(qCMPDURATION)*initialnotional);

	Result = fabs(Result*10000.);
	if (MtM == 0)
		SetSpread(Result); // set ds le cache

	return (Result);
}


// -------------------------------------------------------------
// Computes Risky Duration
// -------------------------------------------------------------

double ICM_Pricer_Cds::ComputeDuration(void)
{
	this->FeeLegPV(); 
	// return GetDuration(); 
	return getValue(qCMPDURATION); 
}

// --------------------------------------------------------------
// Recherche de l'indice k tq max(=0..N/Tk<u)
// --------------------------------------------------------------
/** 
int ICM_Pricer_Cds::GetIndice(const ARM_Vector& vecteur,double t)
{
	bool flag = true;
	int size = vecteur.GetSize();

	for (int i=0; i<size; i++)
	{
		if (fabs(vecteur.Elt(i)-t)<1E-4)
		{
			flag = false;
			break;
		}
		else if (vecteur.Elt(i)>t)
		{
			i--;
			flag = false;
			break;
		}
	}
 
	if (flag) return -1;

	return (i);
}
**/ 

// ------------------------------------------------------------
// Compute Sensitivity Method
// ------------------------------------------------------------
double ICM_Pricer_Cds::ComputeSensitivity(const qSENSITIVITY_TYPE& typesensi, 
										  const std::string& plot, 
										  const std::string&  label, 
										  double  epsvalue, double epsilonGamma // useless
										  )
{
    double sensitivity =0.;

	ICM_DefaultCurveModel* DefModel = (ICM_DefaultCurveModel*) GetModel();

	ARM_Date ExecutionDate = GetModel()->GetStartDate();

	double result = 0.0;
	double initialprice = 0.0;
	double modifiedprice = 0.0;

    
    {

		if (GetInitialPriceFlg())
			initialprice = GetInitialPrice();
		else
			initialprice = Price(qCMPPRICE);

		ResetPricer();

		switch (typesensi)
		{
			case ICMRECOVERY_TYPE :
			case ICMIRCURVE_TYPE :
			case ICMSPREAD_TYPE :
			case ICMSPRELSHIFT_TYPE :
			case ICM_DTR_TYPE :
			{
				ICM_DefaultCurveModel* ModelDef = DefModel->GenerateShiftModel(typesensi,
																			   plot, 
																			   epsvalue);

				SetModel(ModelDef);	


				modifiedprice = ComputePrice(qCMPPRICE);

				if (ModelDef)
					delete ModelDef;
				ModelDef = NULL;

				SetModel(DefModel); //On reset le model initial
			}
			break;
			default :
			result = -99999999.0;
		}

	ResetPricer();

	if (!result)
		result=modifiedprice - initialprice;

	}
   

	return (result);

}

//----------------------------------------------------------------------
// Calcul de la defleg PV Summit compliant
//----------------------------------------------------------------------
double ICM_Pricer_Cds::DefLegPV_Summit()
{

	ICMQUANTIFYER("DefLegPV_Summit"); 
	// ExecutionDate unsued. 
	// if (GetDefLegPriceFlg()) return GetDefLegPrice();

	ICM_Cds* cds = (ICM_Cds*) GetSecurity();
	ARM_Date AsOf = GetAsOfDate();

	//Default Leg
	ICM_Leg* Defleg = cds->GetDefLeg();

	Defleg->GetCreditInfos()->ComputeYF(AsOf);
	ICM_Security* security = Defleg->GetCreditInfos(); 

	const ICM_DefaultCurve* defaultcurve = itsDefaultCurve;

	double DefPayLag = cds->GetCreditLag();
	int NumFlows = security->GetAccEndDates().GetSize();

	bool isConstantNotional(true) ; 
	double constantNotionalValue(0); 
	// the AMOUNT is not reproduced in the CreditInfo object
	if ( Defleg->GetAmount() && Defleg->GetAmount()->size()>1) isConstantNotional=false; 
	else 
	{
		int i = security->PeriodIndex(AsOf);
		if (i<0) i= 0;
		constantNotionalValue= security->GetNotionals().Elt(i) ; 
	}

	
	ARM_Vector PeriodDefLegPV ; 
	ARM_Vector PeriodRefDefLegPV ; 
	ARM_Vector PeriodRecDefLegPV ; 
	if (!GetFaster())
	{
		PeriodDefLegPV.Resize(NumFlows); 
		PeriodRefDefLegPV.Resize(NumFlows); 
		PeriodRecDefLegPV.Resize(NumFlows); 
	}

	const ARM_Vector& YFAccStartDates = security->GetYFAccStartDates();
	double YFFirstAccStart = MAX(YFAccStartDates.Elt(0),0.);
	double YFAccEndDate = security->GetYFAccEndDates().Elt(NumFlows-1);

	double DefPV = 0.,TmpFlow = 0. ; // JLA , notional = 0.;
	double ProtectNot = 0.,ProtectRecovery = 0.;
	double CumProtectNot = 0.,CumProtectRecovery = 0.;

	double Recovery = defaultcurve->GetRecovery();  
	double Ref = itsRefNotional;

	if (Defleg->GetBinaryFlg()) {Recovery= Defleg->GetBinary();}  

	int indice = 0;
	const ARM_Vector& notionals = security->GetNotionals();
	double t1=0.,t2=0.,Lambda=0.,zrates=0.,integral=0.;
	double SP = 0.,ZP = 0.,slope = 0.;

	int nxchg = Defleg->GetNxFlag();

	//JLA notional = notionals->Elt(0);

	// ---------------------------------------------------------------
	// Modified by LJ --> entity
	// ---------------------------------------------------------------
	double	the_lag;
	double	adjusted_lag;
	// ---------------------------------------------------------------
	// ENTITY SCHEDULE: lag is tsl for Entity Leg
	double	ssl_lag	=	DefPayLag;
	the_lag	=	ssl_lag;
	adjusted_lag	=	1.4 *	the_lag;
	
	ARM_Vector Integration_Schedule_Entity ; 
	// JLA-Chaching GenerateMergeSchedule_Summit(qStyle_Recovery_Leg, adjusted_lag,Integration_Schedule_Entity);
	getSchedule_Summit(qStyle_Recovery_Leg,adjusted_lag,Integration_Schedule_Entity); 
 
	// ---------------------------------------------------------------

	int		i;
	double	Start;
	double	Maturity;
	double	Coef, alpha_poly;
	ICM_Poly	Poly;

	if (Ref)
	{
	for (i=0; i<Integration_Schedule_Entity.GetSize()-1; i++)
	{
		t1 = Integration_Schedule_Entity.Elt(i);
		t2 = Integration_Schedule_Entity.Elt(i+1);

		if (((t1 >= YFFirstAccStart) && (t1 < YFAccEndDate))
				|| ((t1<YFFirstAccStart) && (t2 > YFFirstAccStart)))
		{
			Start	= MAX(t1,YFFirstAccStart);
			Maturity = MIN(t2,YFAccEndDate);
			
			SP = defaultcurve->SurvivalProba(Start);
			// JLA Lambda = Intensity->Elt(i);
			Lambda = defaultcurve->GetPWCIntensity( t1 ); 
			//FIX: is curve in default Coef is likely to bo INF, but flow should be 0
			if (defaultcurve->GetIsNameInDefault()) Coef=0; 
			else Coef = exp(Lambda*Start);

			alpha_poly = 0.;
			// Laurent
			CptPolyIrCurve(Poly,alpha_poly,t1,t2,Lambda,adjusted_lag);
			integral = Poly.EvaluateIntegraleExp(alpha_poly,Start,Maturity);
			// ProtectNot = Coef*SP*Ref*Lambda*integral* notional;
			
			if (!isConstantNotional || !GetFaster() ) 
			{
				// indice = GetIndice(security->GetYFAccStartDates(),t1);
				indice = locateIndex(&unconst(security->GetYFAccStartDates()),t1);

				if (indice<0) indice = security->PeriodIndex(AsOf);
				if (indice<0) indice = 0;
			}
			if ( isConstantNotional ) 
				ProtectNot = Coef*SP*Ref*Lambda*integral* constantNotionalValue;
			else 
			{
				ProtectNot = Coef*SP*Ref*Lambda*integral* notionals.Elt(indice);
			}
			
			CumProtectNot += ProtectNot;
	
			if (!GetFaster()) 
			{
				if (indice<0) indice = 0;
				PeriodRefDefLegPV.Elt(indice)=ProtectNot +PeriodRefDefLegPV.Elt(indice); 
			}
		}
	}
	}

	// Restore IR Curve
	// useless now : RestoreIRCurve(adjusted_lag); 

	// ---------------------------------------------------------------
	// REFERENCE SCHEDULE: lag is ssl for Default Leg
	double	tsl_lag	=	3.0;
	the_lag	=	tsl_lag;
	adjusted_lag	=	1.4 * the_lag;
	ARM_Vector Integration_Schedule_Reference ; 
	// GenerateMergeSchedule_Summit(qStyle_Recovery_Leg, adjusted_lag,Integration_Schedule_Reference );
	getSchedule_Summit(qStyle_Recovery_Leg, adjusted_lag,Integration_Schedule_Reference );

	// ---------------------------------------------------------------

	for (i=0; i<Integration_Schedule_Reference.GetSize()-1; i++)
	{
		t1 = Integration_Schedule_Reference.Elt(i);
		t2 = Integration_Schedule_Reference.Elt(i+1);

		if (((t1 >= YFFirstAccStart) && (t1 < YFAccEndDate))
				|| ((t1<YFFirstAccStart) && (t2 > YFFirstAccStart)))
		{
			Start	= MAX(t1,YFFirstAccStart);
			Maturity = MIN(t2,YFAccEndDate);
			
			//JLA Lambda = Intensity->Elt(i);
			Lambda =defaultcurve->GetPWCIntensity( t1 ) ;
			SP = defaultcurve->SurvivalProba(Start);
			if (defaultcurve->GetIsNameInDefault()) Coef=0; 
			else Coef = exp(Lambda*Start);
		
			alpha_poly = 0.;
			// Laurent
			CptPolyIrCurve(Poly,alpha_poly,t1,t2,Lambda,adjusted_lag);
			integral = Poly.EvaluateIntegraleExp(alpha_poly,Start,Maturity);
			
			//JLA ProtectRecovery = Coef*SP*Recovery*Lambda*integral*notional;
			if ( !isConstantNotional || !GetFaster() ) 
			{
				// indice = GetIndice(security->GetYFAccStartDates(),t1);
				indice = locateIndex(&unconst(security->GetYFAccStartDates()),t1);

				if (indice<0) indice = security->PeriodIndex(AsOf);
				if (indice<0) indice = 0;
			}
			if (isConstantNotional) 
				ProtectRecovery = Coef*SP*Recovery*Lambda*integral*constantNotionalValue;
			else 
				ProtectRecovery = Coef*SP*Recovery*Lambda*integral*notionals.Elt(indice);
			
			CumProtectRecovery += ProtectRecovery;

			if (!GetFaster()) 
			{
			if (indice<0) indice = 0;
			PeriodRecDefLegPV.Elt(indice)=ProtectRecovery +PeriodRecDefLegPV.Elt(indice); 
			}

		}

	}

	// Restore IR Curve
	// useless now : RestoreIRCurve(adjusted_lag);

	if (defaultcurve->GetIsNameInDefault()) 
		DefPV	= (1.-Recovery) * constantNotionalValue ; 
	else 
		DefPV	= (CumProtectNot - CumProtectRecovery); 


	//for notional exchange at maturity (CLN case)
	if ((nxchg==K_NX_END)||(nxchg==K_NX_BOTH))
	{	
		ARM_Date Maturity = security->GetEndDateNA();
		DefPV -= constantNotionalValue *
				defaultcurve->SurvivalProba(Maturity)*
				itsDiscountCurve->DiscountPrice(Maturity);}

	if (!GetFaster())
	{	
		security->SetPeriodRefDefLegPV(PeriodRefDefLegPV);
		security->SetPeriodRecDefLegPV(PeriodRecDefLegPV);

		for (i=0;i<PeriodRefDefLegPV.GetSize();i++)
		{PeriodDefLegPV.Elt(i)=PeriodRecDefLegPV.Elt(i)-PeriodRefDefLegPV.Elt(i);}

		security->SetPeriodDefLegPV(PeriodDefLegPV);
	}

	SetDefLegPrice(DefPV);

	return (DefPV);
}

//----------------------------------------------------------------------
// Calcul de la Fee Leg PV Summit compliant
//----------------------------------------------------------------------

double ICM_Pricer_Cds::FeeLegPV_Summit()
{
	ICMQUANTIFYER("FeeLegPV_Summit"); 
	// ExecutionDate unsued. 
	//if (GetFeeLegPriceFlg()) return GetFeeLegPrice();

	// GETTING SOME DATA

	// Security
	ICM_Cds* cds = (ICM_Cds*) GetSecurity();
	
	// Model
	//ICM_DefaultCurveModel* model = (ICM_DefaultCurveModel*) GetModel();
	// Date
	ARM_Date AsOf = GetAsOfDate();

	// Fee Leg
	ICM_Leg* Feeleg = cds->GetFeeLeg();
	// Compute schedule Year Fractions
	Feeleg->GetCreditInfos()->ComputeYF(AsOf);
	ICM_Security* security = Feeleg->GetCreditInfos(); 

	// Zero Curve
	ARM_ZeroCurve* ircurve = itsDiscountCurve;
	const ICM_DefaultCurve* defaultcurve = itsDefaultCurve;

	// The Default Lag
	double DefPayLag = cds->GetCreditLag();

	// Flows
	
	
	const ARM_Vector& AccStartDates = security->GetAccStartDates();
	const ARM_Vector& AccEndDates = security->GetAccEndDates();
	// 14514 int	NumFlows = security->GetNumFlows();
	int	NumFlows = AccStartDates.size() ;
	// ---------------------------------------------------------------
	// PAYMENT TYPE
	qPAYMENT_PREMIUM_LEG	ThePremiumLegPaymentType;
	ThePremiumLegPaymentType	=	cds->GetFeeLeg()->GetAccruedOnDefault();
	// ---------------------------------------------------------------

	// ---------------------------------------------------------------
	// Modified by LJ --> entity
	// ---------------------------------------------------------------
	
	// ---------------------------------------------------------------
	// Start and End dates
	const ARM_Vector& YFAccStartDates = security->GetYFAccStartDates();
	const ARM_Vector& YFAccEndDates = security->GetYFAccEndDates();
	const ARM_Vector& YFPayDates = security->GetYFPayDates();
	const ARM_Vector& Notionals = security->GetNotionals();
	const ARM_Vector& YFInterestDays = security->GetYFInterestDays();
	ARM_Vector DefProb ; // = security->GetDefaultProbability();
	ARM_Vector Discount; // = security->GetDiscountRate();

	if (!GetFaster()) 
	{	DefProb.Resize(Notionals.GetSize());
		Discount.Resize(Notionals.GetSize());	
	}

	double YFFirstAccStart	= MAX(YFAccStartDates.Elt(0),0.);
	double YFAccEndDate		= YFAccEndDates[NumFlows-1];
	// ---------------------------------------------------------------

	double	Start;
	double	Maturity;
	double	Coef, alpha_poly;
	ICM_Poly	Poly;

	double t1=0.,t2=0.,Lambda=0.,integral=0.;
	double	DP, SP, Coupon, Accrued, FeeLegPV;
	double	CumPremium = 0.0;
// 	bool	GoodPeriod;
	double Duration = 0, Im_Duration = 0,DeltaYF=0.;; 
	double fixedSpread  = cds->GetCreditSpread(); 

	double	tStart, tEnd;
	if (Feeleg->GetCreditIndex()) {Feeleg->GetCreditIndex()->ResetDefCurve();}
	Feeleg->CptCoupons(GetModel(),GetAsOfDate());

	// ----------------------------------------------------------------------------------------------
	// Full Period Pricing
	// ----------------------------------------------------------------------------------------------

	{
		// JLA. Rewriting. We iterate over all the flows that are to be paid. 
		// 
		double CumPremium2=0; 
		double Duration2=0; 
		int iPayment ; // 
		int iPayment0 = security->PaymentIndex(AsOf) ;	// index of first flow to be paid (-1 if no flows) 
		int iAccrued = security->PeriodIndex(AsOf) ;	// index of the period including AsOf (-1 if no period) 
		if (iPayment0==-1) 
		{
			ICMLOG("WARNING: Leg has no payments... "); 
			SetFeeLegPrice(0);
			SetDuration(0); 
			return 0;
		}
		// if (iAccrued<iPayment0) ICMTHROW(ERR_INVALID_ARGUMENT,"Accrued unhandled ... ") ; 
		for (iPayment=iPayment0;iPayment<NumFlows;iPayment++) 
		{
			Coupon = security->FullCoupon(iPayment); 
			DP = ircurve->DiscountPrice(YFPayDates[iPayment]);
			switch (ThePremiumLegPaymentType)
			{
				case qCONTINUE_TO_MATURITY:
					SP	=	1.0;
					break;

				case qACCRUED_NOT_SETTLED:
				case qACCRUED_SETTLED:
					SP = defaultcurve->SurvivalProba(YFAccEndDates.Elt(iPayment));
					break;

				case qCOMPLETE_CURRENT_PERIOD:
					SP = defaultcurve->SurvivalProba(YFAccStartDates.Elt(iPayment));
					break;
			}
			double DeltaYF= YFInterestDays.Elt(iPayment) ;
			
			Im_Duration = DeltaYF * DP * SP;
			Duration2 += Im_Duration;
			//Duration based on NetPrice
			if (iPayment!=iAccrued) 
				FeeLegPV = Coupon  / DeltaYF * Im_Duration  ;
			else 
			{
				double accrued = security->AccruedCoupon(AsOf) ; 
				double DeltaAcc = DeltaYF - CountYears(security->GetAccrualBasis(),AsOf.GetJulian(),AccEndDates.Elt(iPayment));
				if (ThePremiumLegPaymentType == qACCRUED_SETTLED)
				{
					FeeLegPV = (Coupon-accrued) / DeltaYF * Im_Duration  + accrued * DP;
					Duration2 += DeltaAcc * (DP-SP*DP);	
				}
				else
					FeeLegPV = Coupon  / DeltaYF * Im_Duration ;
				
				Duration2 += - DeltaAcc;	
			}
			CumPremium2 += FeeLegPV; 
			// ICMLOG("j ="<<iPayement<<" DP="<<DP<<" SP="<<SP<<" YF="<<DeltaYF<<" PV="<<FeeLegPV); 
			// view only 
			if (!GetFaster())
			{
				DefProb.Elt(iPayment)= 1. - SP;
				Discount.Elt(iPayment)= DP;
			}
		}
		CumPremium=CumPremium2; 
		Duration=Duration2 ;
	}
	

	if (ThePremiumLegPaymentType == qCONTINUE_TO_MATURITY) 
	{
		CumPremium	*=	defaultcurve->SurvivalProba(YFFirstAccStart);
		Duration	*=	defaultcurve->SurvivalProba(YFFirstAccStart);
	}


	// ---------------------------------------------------------------
	// ACCRUED
	
	double	CumPremiumAccrued, PremiumAccrued, FeePV;
	CumPremiumAccrued	=	0.0;
	
	// ---------------------------------------------------------------
	if (ThePremiumLegPaymentType == qACCRUED_SETTLED)
	{
		// Get Merged Schedule  Summit like
		// ARM_Vector* Integration_Schedule_Entity = 0; 
		// ARM_Vector* Intensity	=	0; 

		double the_lag	=	DefPayLag;
		double adjusted_lag	=	1.4 *	the_lag;

		// Get Merged Schedule  Summit like
		ARM_Vector  Integration_Schedule_Entity ; 
		// GenerateMergeSchedule_Summit(qStyle_Premium_Leg, adjusted_lag,Integration_Schedule_Entity );
		getSchedule_Summit(qStyle_Premium_Leg, adjusted_lag,Integration_Schedule_Entity );

		int indice = 0;
	
		// ---------------------------------------------------------------
		// PREMIUM SCHEDULE: lag is tsl for Premium Leg (Accrued cases)	
		// ---------------------------------------------------------------

		//JLA useless int	Poly_order;

		// otherwise, play with AccSch and indice
		int i=0;
		while (YFAccStartDates[i] <= YFFirstAccStart)
		{
			i++;
			if (i  == NumFlows)
				break;
		}

		if ((i == NumFlows) && (YFAccStartDates[NumFlows-1] < YFFirstAccStart))
			// pb. no flows
			return -1.e25;

		int next_index;
		double	next_fee_schedule_date;
		double	last_fee_schedule_date;
		double	ratio_schedule_date;

		if (i != NumFlows)
		{
			next_index	=	i;
			next_fee_schedule_date	=	YFAccStartDates[i];
			last_fee_schedule_date	=	YFAccStartDates[i-1];
			ratio_schedule_date	= (next_fee_schedule_date - last_fee_schedule_date);
		}
		else
		{
			next_index	=	i-1;
			next_fee_schedule_date	=	YFAccEndDates[i-1];
			last_fee_schedule_date	=	YFAccStartDates[i-1];
			ratio_schedule_date	= (next_fee_schedule_date - last_fee_schedule_date);
		}


		for (i=0; i<Integration_Schedule_Entity.GetSize()-1; i++)
		{
			t1 = Integration_Schedule_Entity.Elt(i);
			t2 = Integration_Schedule_Entity.Elt(i+1);
			Im_Duration = 0.;

			if (((t1 >= YFFirstAccStart) && (t1 < YFAccEndDate))
					|| ((t1<YFFirstAccStart) && (t2 > YFFirstAccStart)))
			{
				Start	= MAX(t1,YFFirstAccStart);
				Maturity = MIN(t2,YFAccEndDate);
				
				SP = defaultcurve->SurvivalProba(Start);

				// Lambda = Intensity->Elt(i);
				Lambda = defaultcurve->GetPWCIntensity(t1) ;
				//FIX: si nom en défaut 
				if (defaultcurve->GetIsNameInDefault()) Coef=0 ;
				else Coef = exp(Lambda*Start); 

				alpha_poly = 0.;
				// Laurent
				CptPolyIrCurve(Poly,alpha_poly,t1,t2,Lambda,adjusted_lag);

				// ----------------------------------------------------------
				// Adjust the polynom
				// JLA useless Poly_order	=	Poly.order();
				// vector of coeffecients -> one order more

				if (t1 >= next_fee_schedule_date)
				{
					last_fee_schedule_date	=	next_fee_schedule_date;
					next_fee_schedule_date	=	YFAccEndDates[next_index];
					next_index++;
					if (! ratio_schedule_date)
						// error
						return -1e25;
					ratio_schedule_date	=	(next_fee_schedule_date - last_fee_schedule_date);
				}

				ICM_Poly	PolyRebate ; 
				PolyRebate.ResetOrder(1); 
				PolyRebate(0)=last_fee_schedule_date; 
				PolyRebate(1)= -1.0; 
				PolyRebate.multiplication(PolyRebate, Poly);

				// ----------------------------------------------------------

				integral = PolyRebate.EvaluateIntegraleExp(alpha_poly,Start,Maturity);
				
				// ----------------------------------------------------------
				if (next_index)
				{
					tEnd = AccEndDates.Elt(next_index-1);
					tStart = AccStartDates.Elt(next_index-1);
					DeltaYF = YFInterestDays.Elt(next_index-1);
				}
				else
				{
					// only one period
					tEnd = AccEndDates.Elt(0);
					tStart = AccStartDates.Elt(0);
					DeltaYF = YFInterestDays.Elt(0);
				}

				Coupon = security->FullCoupon((ARM_Date) tStart); //Full Coupon

				if ((AsOf.GetJulian()>tStart) &&  (AsOf.GetJulian()<tEnd)) 
				{
					Accrued = security->AccruedCoupon(AsOf); //Accrued coupon
					Coupon -= Accrued; //broken period
					DeltaYF += -CountYears(security->GetAccrualBasis(),AsOf.GetJulian(),tEnd);
				}
				else if /*(*/((AsOf.GetJulian()>tStart) &&  (AsOf.GetJulian()>=tEnd)) /*|| (AsOf.GetJulian()<tStart)) */
					Coupon = 0.; //Empty Period 
				// ----------------------------------------------------------

				Im_Duration = DeltaYF*Coef*SP*Lambda*integral/ratio_schedule_date;
				//Duration += Im_Duration;
				PremiumAccrued = Im_Duration*Coupon/DeltaYF;

				CumPremiumAccrued += PremiumAccrued;
				Duration -= Im_Duration;
			}
		}

		// Restore IR Curve
		// useless now : RestoreIRCurve(adjusted_lag);

		FeePV	= (CumPremium - CumPremiumAccrued); 
	}
	else
		FeePV	= CumPremium;

	//on applique ce correctif pour le momment
	if (fixedSpread) Duration = 100.*(FeePV - security->AccruedCoupon(AsOf))/(Notionals.Elt(0)*fixedSpread);

	SetFeeLegPrice(FeePV);
	SetDuration(fabs(Duration)); 

	if (!GetFaster()) 
	{
		security->SetDefaultProbability(DefProb) ; 
		security->SetDiscountRate(Discount); 
	}

	return (FeePV);
}


// ----------------------------------------------------------------
//Création de l'échéancier mergé
// ----------------------------------------------------------------
void ICM_Pricer_Cds::GenerateMergeSchedule_Summit(const qCredit_Leg_Style& Leg_Type, 
												  const double& lag,ARM_Vector& res)
{
	// ICMQUANTIFYER("ICM_Pricer_Cds::GenerateMergeSchedule_Summit"); 

	// Merge of Interest Rates and Default Intensity Schedule without any lag 
	// always recomputes...

	int ind = 0,i = 0;

	ICM_Cds* cds = (ICM_Cds*) GetSecurity();
	ARM_Date AsOf = GetAsOfDate();

	ARM_ZeroCurve* ircurve = itsDiscountCurve;
	const ICM_DefaultCurve* defaultcurve = itsDefaultCurve;

	ICM_Security* security	=	NULL;
	if (Leg_Type == qStyle_Premium_Leg)
	{
		// I will just add, the first begin date
		ICM_Leg* Feeleg = cds->GetFeeLeg();
		Feeleg->GetCreditInfos()->ComputeYF(AsOf);
		security = Feeleg->GetCreditInfos(); 
	}
	else 
	{
		ICM_Leg* Defleg = cds->GetDefLeg();
		Defleg->GetCreditInfos()->ComputeYF(AsOf);
		security = Defleg->GetCreditInfos(); 
	}

	int NumFlows = security->GetAccEndDates().GetSize();
	double Maturity = security->GetYFAccEndDates().Elt(NumFlows-1);

	// ------------------------------------------------
	// INTEREST RATES SCHEDULE
	// Adjusted with the lag
	ARM_Vector YFIRsch (* ircurve->GetYearTerms()) ;
	YFIRsch -=lag/ K_YEAR_LEN; 

	
	// -------------------------------------------------
	// Adjust the 0.0 Date
	const ARM_Vector* YFDPsch = & defaultcurve->GetInterpolYF(); 
	ARM_Vector  YFDPsch_Copy(defaultcurve->GetInterpolYF()) ; 
	YFDPsch_Copy.Elt(0)	=	INFINITY_MINUS;

	ARM_Vector* Tmp_Vector = NULL;
	ARM_Vector* Tmp_Vector2 = NULL;

	// First Merge
	MergeDates(&Tmp_Vector,&YFIRsch,&YFDPsch_Copy);



	// ------------------------------------------------
	// SCHEDULE
	if (Leg_Type == qStyle_Premium_Leg)
	{
		// I will just add, the first begin date
		// ICM_Leg* Feeleg = cds->GetFeeLeg();
		// Feeleg->GetCreditInfos()->ComputeYF(AsOf);
		// security = Feeleg->GetCreditInfos(); 

		const ARM_Vector& YFAccStartDates = security->GetYFAccStartDates();
		const ARM_Vector& YFAccEndDates = security->GetYFAccEndDates();
		
		ind	=	YFAccStartDates.GetSize();

		// with no lag, Begin(i+1) = End(i) / i=1..NFlow-1
		ARM_Vector YFBeginEndDates(ind + 1,0.);
			YFBeginEndDates[0]	=	YFAccStartDates[0];
		for (i=0;i<ind;i++)
			YFBeginEndDates[i+1]	=	YFAccEndDates[i];


		// ------------------------------------------------

		// Second Merge, with Schedule
		MergeDates(&Tmp_Vector2,Tmp_Vector,&YFBeginEndDates);
	}
	else
	{
		// ICM_Leg* Defleg = cds->GetDefLeg();
		// Defleg->GetCreditInfos()->ComputeYF(AsOf);
		// security = Defleg->GetCreditInfos(); 

		Tmp_Vector2	=	Tmp_Vector;
	}
	// ------------------------------------------------

	// Finally, do not keep not relevant periods
	ind = 0;
	for (i=0;i<Tmp_Vector2->GetSize();i++)
	{
		if ((Tmp_Vector2->Elt(i)>Maturity) && (fabs(Tmp_Vector2->Elt(i)-Maturity)>DB_TOL))
			break;
	}

	// ind = MAX(Tmp_Vector2->GetSize(),i+1); 
	ind = MIN (Tmp_Vector2->GetSize(),i+1); 
	
 
	// Integration_Schedule = new ARM_Vector(ind,0.);
	// ICM_Leg* leg=NULL;
	// if (Leg_Type == qStyle_Premium_Leg) {leg = cds->GetFeeLeg();}
	// else {leg = cds->GetDefLeg();}

	res.Resize(ind); 
	for (i=0;i<ind;i++)	{res.Elt(i)=Tmp_Vector2->Elt(i);}


	if (Tmp_Vector) delete Tmp_Vector;
	if (Tmp_Vector2 && (Leg_Type == qStyle_Premium_Leg))
		delete Tmp_Vector2;
}
 
 

//------------------------------------------------------//
//Calcul du Polynome pour la IrCurve sur l'intervalle [yf1,yf2]
//------------------------------------------------------//
void ICM_Pricer_Cds::CptPolyIrCurve(ICM_Poly& poly,double& alpha, double yf1,double yf2,double Lambda,double lag)
{
	int i = 0, size=0;
	double yf_ti=0,yf_ti_next=0,rate=0,rate_next=0,slope=0;
	double Debut=0, a0=0, a1=0, a2 = 0;
	double ssl=0., tsl=0.;
	// useless ARM_Date date1,date2; 

	ARM_ZeroCurve* ircurve = itsDiscountCurve;

	ARM_Date Asof =  GetAsOfDate();

	ARM_Vector* YFIRsch = ircurve->GetYearTerms(); // JLA : curve items are not lagged. 
	ARM_Vector* IRrates = ircurve->GetZeroRates();

	int	NbIRPlots	=	YFIRsch->GetSize();

	double	*coefs = NULL;

	alpha=0.;
	// Deal with cases outside the IR range
	// if (yf1 < YFIRsch->Elt(0))
	if (yf1 < YFIRsch->Elt(0) - lag/K_YEAR_LEN )
	{
		// the ZC rate will be kept constant, linear in the Exp
		// rate = ircurve->DiscountYield(0)/100.;
		rate = IRrates->Elt(0)/100.; 
		a0	=	- rate * lag / 365.0;
		a1	=	- (rate + Lambda);
		alpha	=	- a1;
		poly.ResetOrder(0); 
		poly(0)=exp(-a0) ;
		return;
	}
	// else if (yf2 > YFIRsch->Elt(NbIRPlots-1))
	else if (yf2 > YFIRsch->Elt(NbIRPlots-1) - lag/K_YEAR_LEN )
	{
		// the ZC rate will be kept constant, linear in the Exp
		// rate = ircurve->DiscountYield(NbIRPlots-1)/100.;
		rate = IRrates->Elt(NbIRPlots-1)/100.; 
		a0	=	- rate * lag / 365.0;
		a1	=	- (rate + Lambda);
		alpha	=	- a1;
		poly.ResetOrder(0); 
		poly(0)=exp(-a0);
		return;
	}

	/** for (i=0;i<NbIRPlots-1;i++)	
	{
		if( leq(YFIRsch->Elt(i)- lag/K_YEAR_LEN,yf1) 
			&& geq(YFIRsch->Elt(i+1)- lag/K_YEAR_LEN,yf2) ) break; 
	}
	**/ 
	
	i=locateIndex(YFIRsch,yf1 + lag/K_YEAR_LEN); 
	if (i==-1) 
		ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_Pricer_Cds::CptPolyIrCurve: can't locateIndex: "<<yf1 + lag/K_YEAR_LEN); 

	yf_ti = YFIRsch->Elt(i);
	yf_ti_next = YFIRsch->Elt(MIN(i+1,NbIRPlots-1));

	if (fabs(yf_ti-yf_ti_next)<DB_TOL)
	{
		ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_Pricer_Cds::CptPolyIrCurve: pb with Schedule when computing SUMMIT ExpPolyoms "); 
		// throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
        //  "Parameters :  pb with Schedule when computing SUMMIT ExpPolyoms");
	}

	// not required to interpol : we are on pivot
	rate = IRrates->Elt(i)/100. ;
	rate_next = IRrates->Elt(MIN(i+1,NbIRPlots-1))/100. ;

	// coeff a
	slope = (rate_next-rate)/(yf_ti_next-yf_ti);	// /365.0;
	// coeff b
	// double constant	=	-  yf_ti * slope + rate; // * slope * 365.0 + rate;
	double constant	=	- ( yf_ti - lag/K_YEAR_LEN  )* slope + rate; // * slope * 365.0 + rate;

	// coeffecients for P polynom in YF
	a0	=	- constant * lag / 365.0;
	a1	=	- (constant + slope * lag / 365.0 + Lambda);
	a2	=	- slope;

	// remark
//	double	tc;
//	tc	=	(yf2 + yf1) / 2.0;
//	alpha	=	- (2 * a2 * tc + a1);

	ICM_Poly::TaylorExpansion(a0,a1,a2,yf1,yf2,alpha,poly);
}
 
// -- JLA Useless

//*******************************************************************************************
// Compute Fwd Spread
//*******************************************************************************************

double ICM_Pricer_Cds::Compute_Fwd_Spread(const ARM_Date &Mty, const ARM_Date & CDS_ExpiryDate, double& FlatRbp_Fwd)
{
	auto_ptr<ICM_DefaultCurve> pdefaultcurve(dynamic_cast<ICM_DefaultCurve*>(itsDefaultCurve->Clone()));
	if (pdefaultcurve.get())
		return pdefaultcurve->FwdSpread(Mty, CDS_ExpiryDate) ;
	else
		return CREDIT_DEFAULT_VALUE;

}	



// ---------------------------------------------------------------
// Full Taylor Expansion for expression Exp(a0+a1*X+a2*X²) on [t1,t2] returns alpha & Poly_x 
// for Exp(alpha*X)*Poly_x
// ---------------------------------------------------------------
/** 

void ICM_Pricer_Cds::TaylorExpansion(const double& a0,
											const double& a1,
											const double& a2,
											const double& t1,
											const double& t2,
											double& alpha,
											ICM_Poly& Poly_x)
	{
		double Px = 0.;
		// double* coefs = NULL;
		int degree = 1;
		double epsilon = 1.E-6;
		int maxloop = 10;
		int deg = 0;
		double tmp = 0.;
		double tc = (t1+t2)/2.;

		//definition du Polynome Q d'ordre 2
		
		ICM_Poly Q ; 
		Q.ResetOrder(2); 
		Q(0) = a0;
		Q(1) = a1;
		Q(2) = a2;

		alpha = Q.EvaluateDerivatives(tc);

		//definition du Polynome P d'ordre 1
		double	Qtc	=	Q.Evaluate(tc);
		tmp = Qtc-tc*alpha;
		Px = exp(tmp);
		ICM_Poly P ; 
		P.ResetOrder(0); 
		P(0)= Px; 
			
		ICM_Poly Q2 ; 
		Q2.ResetOrder(2); 
		Q2(0)=Q(0) - tmp;
		Q2(1)=Q(1) - alpha;
		Q2(2)=a2;

		ICM_Poly R ; 
		R.ResetOrder(0) ;
		R(0)=1; 

		// JLA ICM_Poly Qd;
		// JLA Qd = R;
		ICM_Poly Qd(R); 
		deg = 1;

		double exp_q2_t1 = exp(Q2.Evaluate(t1)); 
		double exp_q2_t2 = exp(Q2.Evaluate(t2)); 
		while (((fabs((exp_q2_t1  -  R.Evaluate(t1))/exp_q2_t1)>epsilon) ||
			 (fabs((exp_q2_t2 -  R.Evaluate(t2))/exp_q2_t2)>epsilon)) &&
			 degree <maxloop)

		{
			Qd.multiplication(Qd,Q2);
			deg *= degree;

			ICM_Poly Qd2(Qd); 

			for (int i =0; i<=Qd.order(); i++)
					// Qd2.itsPoly[i] /= (double)deg;
					Qd2(i) /= (double)deg;

			R.somme(R,Qd2);

			degree++;
		}

		
		Poly_x.multiplication(R,P); 

		return;
	}
	**/ 
// virtual 
ICM_Pricer* ICM_Pricer_Cds::CloneOnlyParams(ARM_Security* sec,ARM_Model* mod)
	{
	ICM_Pricer_Cds* theClone = new ICM_Pricer_Cds; 
	theClone->Set(sec,mod,ICM_Parameters(),GetAsOfDate());
    return(theClone);
	}


void ICM_Pricer_Cds::Init()
{
	SetName(ICM_PRICER_CDS);
	itsRefNotional = 1.;
	itsDefaultCurve = NULL;
	itsDiscountCurve = NULL;
	itsCouponCurve = NULL;
	itsIntegrationSchedules.clear(); 
}
void 
ICM_Pricer_Cds::getSchedule_Summit(const qCredit_Leg_Style& Leg_Type, const double& lag , ARM_Vector& res)
{
	//ICMQUANTIFYER("ICM_Pricer_Cds::getSchedule_Summit"); 
	ICM_Pricer_Cds::disc_index key; 
	key.lag=lag; 
	key.legStyle=Leg_Type; 
	discs_t::iterator it = itsIntegrationSchedules.find(key) ; 
	if (it!=itsIntegrationSchedules.end()) res=it->second; 
	else 
	{
		GenerateMergeSchedule_Summit(Leg_Type,lag,res); 
		itsIntegrationSchedules[key]=res; 
	}
}
