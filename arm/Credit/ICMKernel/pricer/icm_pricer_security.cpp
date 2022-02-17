#include "firsttoinc.h"
#include "ICMKernel/pricer/icm_pricer_security.h"
#include "ICMKernel/mod/icm_defcurvemodel.h"
#include "ICMKernel/inst/icm_ftd.h"
#include "ICMKernel/inst/icm_mez.h"
#include "ICMKernel/inst/icm_credit_index.h"
#include "ICMKernel/glob/icm_smile_correlation.h"
#include "ICMKernel/inst/icm_corridorleg.h"
#include "ICMKernel/crv/icm_distriblosscalculator.h" 
#include "ICMKernel/mod/modelmulticurves.h" 


void ICM_Pricer_Security::Set(ARM_Security *sec, ARM_Model *mod,const ICM_Parameters&params,const ARM_Date&asof)
{
	if ((mod->GetName()!=ICM_MODELMULTICURVES))			
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "ERROR : pricer not compatible with this security or this model ");

	((ICM_Mez*)sec)->RebuildAfterDefault((ICM_ModelMultiCurves *)mod);

	ICM_Pricer::Set(sec, mod,params,&asof);
}

// ------------------------------------------------------------------------
// Calcul de l'expected loss tranche
// ------------------------------------------------------------------------
void ICM_Pricer_Security::CptExpectedLossTranche()	
{
	if (!itsDistribLoss.empty()) return; 

	ICM_Cds* cds = (ICM_Cds*) GetSecurity();
	ICM_ModelMultiCurves* model = (ICM_ModelMultiCurves*) GetModel();
	ICM_Correlation* correl = model->GetCorrelation();
	ARM_Date AsOf = model->GetStartDate();
	ARM_Vector*	Schedule = NULL;

	vector<double> losses;
	itsScheduleForEL.clear();
	itsDetailLosses.clear();


	const ICM_Parameters& parameters = GetParameters();
	// if (parameters)
	if (!parameters.empty())
	{	ARM_Vector* VResc = parameters.GetColVect("TERMS_RESCALING");
		if (VResc) 
		{
			qTERM_STRUCTURE type;
			bool ret;
			ICM_EnumsCnv::cnv(VResc->Elt(0),type , ret);
			SetTermStructurePricing(type);	
		}
	}

	//frozen maturity case
	ICM_Cds* oldsec = NULL;
	ICM_Cds* newsec = NULL;
	bool frozenmatu = cds->GetFeeLeg()->GetCreditInfos()->IsFrozenMaturity();		//frozen maturity case
	double YFfrozenmatu = 0.;														//frozen maturity case	
	if (frozenmatu)																	//frozen maturity case	
	{																				//frozen maturity case
		oldsec = cds;																//frozen maturity case				
		newsec = (ICM_Cds*) cds->Clone();											//frozen maturity case
		ARM_Date FrozMat = cds->GetFeeLeg()->GetCreditInfos()->GetFrozenMaturity();	//frozen maturity case
		YFfrozenmatu = MAX((FrozMat-AsOf)/365.,0.);									//frozen maturity case		
		newsec->CptCashFlowDatesCredit(FrozMat);									//frozen maturity case
		SetSecurity(newsec);														//frozen maturity case
		cds = newsec;																//frozen maturity case
	}																				//frozen maturity case

	//Estimation de la Loss Distribution pour chaque pas de temps
	//Si utilisation de la Term Structure :
	//Il faut merger les dates des bourbes de Base Correl par terme
	//Il faut faire les estimations deux fois pour ne pas prendre en compte
	//les sauts aux points de discontinuité

	//Génération du Schedule
	cds->GenerateSchedule(AsOf,itsStepForDistribLoss);
	Schedule = cds->GetSchedule();
	if (Schedule == NULL) return;
	int ScheduleSize	=	Schedule->GetSize();
	itsTermStructureSched.Resize(0);

	if ((Schedule != NULL)&&
		(GetTermStructurePricing()==qTermStructure))
	//test si correl existe (différence rescaling ou pas)
	{
		double YFStart = (cds->GetStartDateNA()-AsOf)/365.;

		int i=0,l=0;
		double checkvalue = 0.;

		//Merge avec les dates des bases correl
		ARM_Vector BaseCorrelSchedule;
		ARM_Vector* NewSchedule=NULL;
		
		correl->GetCorrelationTerms(BaseCorrelSchedule);

		double StartSchedule = Schedule->Elt(1);
		//cas AsOf == start Date
		if (Schedule->GetSize()==2) {StartSchedule = Schedule->Elt(0);}

		double MatSchedule = Schedule->Elt(Schedule->size()-1);
		int maxidx = BaseCorrelSchedule.GetSize();
		int minid = 0;
		
		int intermeth = correl->GetInterpType();

		// On supprime les points de la BC avant la start date & apres la maturity
		for (i=0;i<BaseCorrelSchedule.GetSize();i++)
		{
		 if ((BaseCorrelSchedule.Elt(i)<StartSchedule) &&
			!CHECK_EQUAL(BaseCorrelSchedule.Elt(i),StartSchedule)) 
			{minid=i+1;}
		 else if ((BaseCorrelSchedule.Elt(i)>MatSchedule) &&
			!CHECK_EQUAL(BaseCorrelSchedule.Elt(i),MatSchedule)) 
			{maxidx=i;break;}	
		}

		ARM_Vector NewScheduleBC(maxidx-minid,0.);
		for (i=minid;i<maxidx;i++){NewScheduleBC.Elt(i-minid)=BaseCorrelSchedule.Elt(i);}

		MergeDates(&NewSchedule,Schedule,&NewScheduleBC);
		itsTermStructureSched.Resize(NewSchedule->size());
		itsTermStructureSched = *NewSchedule; 
				
		itsDistribLoss.SetUseStdEL(false);
		//Modif du type : continue à gauche (sur les courbes EUR et US)
		//correl->SetInterpType(K_ICM_STEPUP_RIGHT_MATURITY);
		correl->SetInterpType(K_LINEAR);
		//correl->SetInterpType(K_ICM_STEPUP_LEFT_LINEAR);
		
		//Estimation de la distrib a gauche
		for (i=0; i<NewScheduleBC.size(); i++)
		{checkvalue = itsDistribLoss[NewScheduleBC.Elt(i)] = ExpectedLossTranche(NewScheduleBC.Elt(i),losses);}

		itsDistribLoss.SetUseStdEL(true);
		//Modif du type : continue à droite (sur les courbes EUR et US)
		//correl->SetInterpType(K_ICM_STEPUP_LEFT_MATURITY);
		correl->SetInterpType(K_LINEAR);
		//correl->SetInterpType(K_ICM_STEPUP_LEFT_LINEAR);

		//Estimation de la distrib a droite
		for (i=0; i<NewSchedule->size(); i++)
		{
			losses.clear();
			checkvalue = itsDistribLoss[NewSchedule->Elt(i)] = ExpectedLossTranche(NewSchedule->Elt(i),losses);
			itsScheduleForEL.push_back(NewSchedule->Elt(i));	
			itsDetailLosses.push_back(losses);
		}

		//on reset la méthode d'interpolation initiale
		correl->SetInterpType(intermeth);

		if (NewSchedule) delete NewSchedule;NewSchedule=NULL;
	}
	else
	{
		for (int i=0; i<ScheduleSize; i++)
		{
			losses.clear();
			itsDistribLoss[Schedule->Elt(i)] = ExpectedLossTranche(Schedule->Elt(i),losses);
			itsScheduleForEL.push_back(Schedule->Elt(i));	
			itsDetailLosses.push_back(losses);
		}
	}

	

	if (frozenmatu)																	//frozen maturity case
	{																				//frozen maturity case
		itsDistribLoss[1000.] = itsDistribLoss.InterpolEL(YFfrozenmatu);			//frozen maturity case
		if (newsec) delete newsec;newsec=NULL;										//frozen maturity case
		SetSecurity(oldsec);cds = oldsec;											//frozen maturity case
	}																				//frozen maturity case
}

double ICM_Pricer_Security::ExpectedLossTranche(const double& yearterm,vector<double>& losses)
{return cpt_elt_default(yearterm,*this,*GetModel(),*GetSecurity(),losses) ;}

// ------------------------------------------------------------------------
// computation of the riksy notional for the feeleg between 2 dates
// ------------------------------------------------------------------------
double ICM_Pricer_Security::RiskyNotional(const int& index)
{
	ARM_Date AsOf = GetModel()->GetStartDate();
	ICM_Cds* cds = (ICM_Cds*) GetSecurity();
	ICM_Security* security = cds->GetFeeLeg()->GetCreditInfos();
	ICM_Security* securitydef = cds->GetDefLeg()->GetCreditInfos();
	security->ComputeYF(AsOf);

	const ARM_Vector& StartDates = security->GetAccStartDates();
	const ARM_Vector& EndDates = security->GetAccEndDates();
	const ARM_Vector& YFStartDates = security->GetYFAccStartDates();
	const ARM_Vector& YFEndDates = security->GetYFAccEndDates();

	//old version
	//ARM_Vector* YFAccEndDates = security->GetYFAccEndDates();
	//return security->GetNotionals()->Elt(index)*
	//		(1.-getDistribLoss().InterpolEL(YFAccEndDates->Elt(index)));

	ARM_Vector* sch = GetSchedule();
	int begin=-1,end=-1,i=0;
	double avgrn=0.;

	//Cas FeeLeg Zero Coupon
	if ((security->GetPaymentFreq()==0) &&(securitydef->GetPaymentFreq()==0))
	{	double EL = getDistribLoss().InterpolEL(YFEndDates.Elt(index)); 
		avgrn =(1.-EL)*ExpectedAmortTranche(YFEndDates.Elt(index));
		avgrn *= security->GetNotionals().Elt(index);
		return avgrn; }

	for (i=0; i<sch->GetSize();i++)
	{
		if CHECK_EQUAL(sch->Elt(i),YFStartDates.Elt(index)) {begin = i;}
		if CHECK_EQUAL(sch->Elt(i),YFEndDates.Elt(index)) {end = i;break;}
	}

	if (security->GetPaymentFreq()!=0)
	{
		for (i=begin; i<=end;i++)
		{double EL = getDistribLoss().InterpolEL(sch->Elt(i)); 
		avgrn+=(1.-EL)*ExpectedAmortTranche(sch->Elt(i));}

		avgrn /= (double)(end-begin+1);
	}
	else //cas feeleg zero coupon only
	{	for (i=begin; i<=end;i++)
		{double EL = getDistribLoss().InterpolEL(sch->Elt(i)); 
		 double YearFraction = 0.;

		 if (i>begin){YearFraction = (365./360.)*(sch->Elt(i)-sch->Elt(i-1));}
		 else {YearFraction = (365./360.)*sch->Elt(i);}

		avgrn+=YearFraction*(1.-EL)*ExpectedAmortTranche(sch->Elt(i));}
		avgrn/=(365./360.)*(sch->Elt(end)-sch->Elt(begin));
	}

	avgrn *= security->GetNotionals().Elt(index);

	return (avgrn);
}

// ------------------------------------------------------------------------
// Calcul de l'accrued
// ------------------------------------------------------------------------
double ICM_Pricer_Security::Accrued()
{
	ICM_Cds* cds = (ICM_Cds*) GetSecurity();
	ICM_ModelMultiCurves* model = (ICM_ModelMultiCurves*) GetModel();
	ARM_Date AsOf = model->GetStartDate();

	ARM_ZeroCurve* ircurve = model->GetZeroCurve();

	ICM_Security* security = NULL;
	security = cds->GetFeeLeg()->GetCreditInfos();
	security->ComputeYF(AsOf);
	// JLA : should be done when computing Accrued first for CMCDO.
//	if (cds->GetFeeLeg()->GetCreditIndex()) cds->GetFeeLeg()->GetCreditIndex()->ResetDefCurve();
	cds->GetFeeLeg()->CptCoupons(model,AsOf);

	
	const ARM_Vector& YFRiskEndDates = security->GetYFEndRiskDates();
	const ARM_Vector& YFPayDates = security->GetYFPayDates();
	// 14514 int NumFlows = security->GetNumFlows();
	int NumFlows = YFRiskEndDates.size() ;

	double Accrued = 0.;
	double coupon =0.,DP =0., EL=0.,notional=0.;
	bool AccruedFlag = false;


/** 	int iPeriodIndex = security->PeriodIndex(AsOf); 
	if (iPeriodIndex!=-1) 
	{
		Coupon = security->AccruedCoupon(AsOf) ;
		switch (cds->GetFeeLeg()->GetAccruedOnDefault())
		{
			case qACC_In_Risk : 
				DP = ircurve->DiscountPrice((ARM_Date)PayDates->Elt(iPeriodIndex));
				SP = defaultcurve->SurvivalProba((ARM_Date)AccEndDates->Elt(iPeriodIndex));
				Accrued = Coupon * DP * SP;
				break;
			case qNOS :
				Accrued = 0.;
				break;
			case qACCRUED_SETTLED : 
			default : 
				Accrued = Coupon;
		}
	}
	
	return (Accrued);
	*/ 

	int iPeriodIndex = security->PeriodIndex(AsOf); 
	if (iPeriodIndex!=-1) 
	{
		Accrued = cds->GetTradedCoef() * security->AccruedCoupon(AsOf) ;
		/** 
		coupon=0 ;
		switch (cds->GetFeeLeg()->GetAccruedOnDefault())
		{
		case qACC_In_Risk : 
			DP = ircurve->DiscountPrice(YFPayDates.Elt(iPeriodIndex));
			EL = getDistribLoss().InterpolEL(YFRiskEndDates.Elt(iPeriodIndex));
			notional = RiskyNotional(iPeriodIndex);
			coupon = security->AccruedCoupon(AsOf,false);
			Accrued = notional * coupon * DP * (1. - EL);
			break;
		case qNOS :
			Accrued = 0.;
			break;
		case qACCRUED_SETTLED : 
		default : 
			DP = ircurve->DiscountPrice(YFPayDates.Elt(iPeriodIndex));
			coupon = security->AccruedCoupon(AsOf);
			Accrued = coupon * DP;
		}
		Accrued *= cds->GetTradedCoef();
		**/ 
	}
	SetAccrued(Accrued);
	return (Accrued);
}



// ------------------------------------------------------------------------
// Calcul de la FeeLeg PV
// ------------------------------------------------------------------------
double ICM_Pricer_Security::FeeLegPV ()
{

	if (GetFeeLegPriceFlg()) return GetFeeLegPrice();

	ICM_Cds* cds = (ICM_Cds*) GetSecurity();
	ICM_ModelMultiCurves* model = (ICM_ModelMultiCurves*) GetModel();
	ARM_Date AsOf = model->GetStartDate();
	ICM_Security* security = NULL;
	ICM_Leg* FeeLeg = cds->GetFeeLeg();

	if (FeeLeg->GetCreditLegType() == qNone_Leg) { SetFeeLegPrice(0.); return 0.;}

	security = FeeLeg->GetCreditInfos();
	security->ComputeYF(AsOf);
	ARM_ZeroCurve* ircurve = model->GetZeroCurve();

	
	// (unused) ARM_Vector* YFAccrualDates = security->GetYFAccrualDates();
	const ARM_Vector& AccStartDates = security->GetAccStartDates();
	const ARM_Vector& AccEndDates = security->GetAccEndDates();
	const ARM_Vector& YFAccEndDates = security->GetYFAccEndDates();
	const ARM_Vector& YFAccStartDates = security->GetYFAccStartDates();
	const ARM_Vector& YFPayDates = security->GetYFPayDates();
	int NumFlows = AccStartDates.size(); 
	ARM_Vector  DefProb ; // = security->GetDefaultProbability();
	ARM_Vector Discount ; // = security->GetDiscountRate();
	const ARM_Vector& YFRiskEndDates = security->GetYFEndRiskDates();

	if (!GetFaster()) 
	{	
		DefProb.Resize(AccStartDates.GetSize());
		Discount.Resize(AccStartDates.GetSize());	
	}

	ARM_Vector PeriodFeeLegPV(NumFlows,0.);
	double  notional=0.;
	ARM_Vector* InitSpreads = NULL;
	double FeePV = 0.;

//	if (FeeLeg->GetCreditIndex()) FeeLeg->GetCreditIndex()->ResetDefCurve();
	FeeLeg->CptCoupons(model,GetAsOfDate());

	ARM_ReferenceValue* NotXchange = cds->GetFeeLeg()->GetNotExchange();

	if (NotXchange)
	{	
		FeePV = NotXchange->CptReferenceValue(AccStartDates.Elt(0));
		double DP = ircurve->DiscountPrice(AccStartDates.Elt(0));
		FeePV *= DP;
	}

	int iPayment ; // 
	int iPayment0 = security->PaymentIndex(AsOf) ;	// index of first flow to be paid (-1 if no flows) 
	int iAccrued = security->PeriodIndex(AsOf) ;	// index of the period including AsOf (-1 if no period) 
	if (iPayment0==-1) 
	{
		ICMLOG("WARNING: Leg has no payments... "); 
		SetFeeLegPrice(0);
		return 0;
	}
	for(iPayment=iPayment0;iPayment<NumFlows;iPayment++)
	{
		double flowPV=0; 
		double Coupon = security->FullCouponRate(iPayment) * RiskyNotional(iPayment); 
		double DP = ircurve->DiscountPrice(YFPayDates[iPayment]);

		//View Only
		if (!GetFaster()) Discount.Elt(iPayment)= DP ;
		if (!GetFaster()) DefProb.Elt(iPayment)= 1. - getDistribLoss().InterpolEL(YFAccEndDates[iPayment]);

		if (NotXchange) Coupon += NotXchange->CptReferenceValue(AccEndDates.Elt(iPayment));
		
		//Estimation du coupon en fonction de la gestion du couru en cas de défaut
		switch (cds->GetFeeLeg()->GetAccruedOnDefault()) 
		{
			//Premium leg sans risque (toujours sur le nominal initial)
			case qCONTINUE_TO_MATURITY:
				flowPV = security->FullCouponRate(iPayment) * DP * security->GetNotionals().Elt(iPayment);
				break;
			//Pas de paiement de couru en cas de défaut (sur le nominal en fin de période)
			case qACCRUED_NOT_SETTLED : 
				{
					flowPV = 1 - getDistribLoss().InterpolEL(YFAccEndDates.Elt(iPayment)); 
					flowPV *= ExpectedAmortTranche(YFAccEndDates.Elt(iPayment));
					flowPV *= security->GetNotionals().Elt(iPayment);
					flowPV *= security->FullCouponRate(iPayment) * DP;
					break;
				}
			
			//Paiement du coupon complet sur la période (nominal en début de période)
			case qCOMPLETE_CURRENT_PERIOD :
				{
					flowPV = 1 - getDistribLoss().InterpolEL(YFAccStartDates.Elt(iPayment)); 
					flowPV *= ExpectedAmortTranche(YFAccStartDates.Elt(iPayment));
					flowPV *= security->GetNotionals().Elt(iPayment);
					flowPV *= security->FullCouponRate(iPayment) * DP;
					break;
				}
			//Paiement du couru sur la période considérée en cas de défaut (nominal moyen sur la période)
			case qACCRUED_SETTLED : 
			default : 
				flowPV = Coupon *DP ;
		}
		PeriodFeeLegPV.Elt(iPayment)=flowPV;
		FeePV += flowPV ; 
	}

	security->SetPeriodFeeLegPV(PeriodFeeLegPV);

	FeePV *= cds->GetTradedCoef();
	if (!GetFaster()) 
	{	
		security->SetDefaultProbability(DefProb); 
		security->SetDiscountRate(Discount); 
	}
	SetFeeLegPrice(FeePV);
	return (FeePV);
}


// ------------------------------------------------------------------------
// Calcul de la DefaultLeg PV
// ------------------------------------------------------------------------
double ICM_Pricer_Security::DefLegPV ()
{

	if (GetDefLegPriceFlg()) return GetDefLegPrice();

	ICM_Cds* cds = (ICM_Cds*) GetSecurity();
	ICM_ModelMultiCurves* model = (ICM_ModelMultiCurves*) GetModel();
	ARM_Date AsOf = model->GetStartDate();
	ICM_Security* security = NULL;
	ARM_Vector sch;

	if (cds->GetDefLeg()->GetCreditLegType() == qNone_Leg) { SetDefLegPrice(0.); return 0.;}

	security = cds->GetDefLeg()->GetCreditInfos();
	security->ComputeYF(AsOf);
	ARM_ZeroCurve* ircurve = model->GetZeroCurve();

	
	const ARM_Vector& YFAccStartDates = security->GetYFAccStartDates();
	const ARM_Vector& YFAccEndDates = security->GetYFAccEndDates();
	int NumFlows = YFAccStartDates.size();
	//Cas Default Zero Coupon
	if (security->GetPaymentFreq()==0) 
	{sch.Resize(2);
	 sch.Elt(0)=YFAccStartDates.Elt(0);;
	 sch.Elt(1)=YFAccEndDates.Elt(NumFlows-1);}
	else 
	{sch = *GetRiskSchedule();}

	//Special for pricing in term structure mode
	//We want in all cases use correlations maturities
	double cptEl = getDistribLoss().InterpolEL(sch.Elt(0));
	if (itsIsTermStructurePricing==qTermStructure) 
		{sch = itsTermStructureSched;}

	double Tranche = GetTranche(AsOf);
	double DefLegPV = 0.,LR = 1.,Notional =0.,DP =0.,EL1 =0.,EL2 =0.;
	double FlagGoodPeriod = 0.,StartPeriod =0.,DefPv = 0.;
	double yflag = ((double)cds->GetCreditLag())/365.*1.4;

	Notional = security->GetNotionals().Elt(0);
	int nxchg = cds->GetDefLeg()->GetNxFlag();

	if (AsOf <= security->GetStartDateNA())
	{	if (nxchg==K_NX_START) {DefPv+=Notional;}}

	//if (cds->GetNbIssuers() == 1) LR= 1. - model->GetRecoveryRate((char*)cds->GetSingleName().c_str());  

	for (int i=0; i<sch.GetSize()-1; i++)
	{
		FlagGoodPeriod = 1.;

		if ((sch.Elt(i)<YFAccStartDates.Elt(0)) && 
			(!CHECK_EQUAL(sch.Elt(i),YFAccStartDates.Elt(0)))) {FlagGoodPeriod = 0.;}

		if ((sch.Elt(i)>YFAccEndDates.Elt(NumFlows-1)) || 
			(CHECK_EQUAL(sch.Elt(i),YFAccEndDates.Elt(NumFlows-1)))) {break;}
		if (cds->GetDefLeg()->GetPaymentFreq()==K_ZEROCOUPON)
			DP = ircurve->DiscountPrice(sch.Elt(i+1)+yflag);
		else 
			DP = (ircurve->DiscountPrice(sch.Elt(i+1)+yflag) + ircurve->DiscountPrice(sch.Elt(i)+yflag))/2;
		EL1 = getDistribLoss().InterpolEL(sch.Elt(i));
		EL2 = getDistribLoss().InterpolEL(sch.Elt(i+1));

		if ((nxchg==K_NX_START) && (cds->GetDefLeg()->GetPaymentFreq()==K_ZEROCOUPON))
			DefLegPV = FlagGoodPeriod * Notional * LR * DP * Tranche* (1. - EL2);
		else
			DefLegPV = FlagGoodPeriod * Notional * LR * DP * Tranche* (EL2 - EL1);

		DefPv+= DefLegPV;
	}

	DefPv *= cds->GetTradedCoef();
	
	SetDefLegPrice(DefPv);
	return (DefPv);
}

/*----------------------------------------------------------------------------*
  Compute Breakeven Spread
*----------------------------------------------------------------------------*/ 
double ICM_Pricer_Security::ComputeSpread(const double& MtM )
{
	if (GetSpreadFlg()) return GetSpread(); //Already computed
	
	ICM_Cds* OriginalCds = (ICM_Cds*) GetSecurity();
	
	// if matured : returns 0 
	ICM_Security* security = OriginalCds->GetFeeLeg()->GetCreditInfos();
	if (security->PaymentIndex(GetModel()->GetStartDate())==-1) 
	{
		ICMLOG("ICM_Pricer_Security::ComputeSpread: matured"); 
		SetSpread(0) ;
		return 0; 
	}

	if (OriginalCds->GetSecurityType() == qCM_TRANCHE) return ComputeImpliedPartRate(); //CM-Cds spread

	ICM_Cds* cds = (ICM_Cds*) OriginalCds->Clone();
	ICM_Leg* Feeleg = (ICM_Leg*) cds->GetFeeLeg();
	qCredit_Leg_Type feelegtype = Feeleg->GetCreditLegType();

	SetSecurity(cds);

	double Result = 0.,NPV = 0.,bp = 0.01,InitialSpread=0.,NewNPV=0.;

	if ((Feeleg->GetLegType() == K_FIXED_LEG) && (Feeleg->GetName() != ICM_CORRIDORLEG))
	{
		InitialSpread = Feeleg->GetCreditSpread();
		//JLA useless, done by Price.. if (!GetPriceFlg()) 
		NPV = Price(qCMPPRICE);

		if  (((GetSecurity()->GetName() == ICM_MEZ ) || 
			  (GetSecurity()->GetName() == ICM_NTD ) ||
			  (GetSecurity()->GetName() == ICM_CDO2 )) 
			  &&  (InitialSpread) 
			  &&  (feelegtype != qStyle_None_Leg))
		{	// Result = InitialSpread * GetDefLegPrice()/(GetFeeLegPrice()-GetAccrued());
			Result = InitialSpread * Price(qCMPDEFLEGPV)/(Price(qCMPFEELEGPV)-Price(qCMPACCRUED));
			Result = fabs(Result*100);	
			SetSpread(Result);

			if (cds) delete cds; cds=NULL;
			SetSecurity(OriginalCds);
			return (Result);
		}

		ResetRootPricer();

		Feeleg->SetCreditSpread(bp);
		NewNPV = Price(qCMPPRICE);
		Result = bp * Price(qCMPDEFLEGPV)/(Price(qCMPFEELEGPV)-Price(qCMPACCRUED));
		Feeleg->SetCreditSpread(InitialSpread);

	} else if (Feeleg->GetName() == ICM_CORRIDORLEG){
		ResetRootPricer();
		double bp1 = 0.01;
		((ICM_CorridorLeg*)Feeleg)->SetLeverageFloatIdx(0.);
		ARM_Vector * vDate = new ARM_Vector(1);
		vDate->Elt(0) = this->GetAsOfDate().GetJulian();
		ARM_Vector * vValue = new ARM_Vector(1);
		vValue->Elt(0) = 0.01;
		ARM_ReferenceValue* newSpread = new ARM_ReferenceValue();
		newSpread->SetDiscreteDates(vDate);
		newSpread->SetDiscreteValues(vValue);
		((ICM_CorridorLeg*)Feeleg)->SetSpreads(newSpread);		
		if (newSpread) delete newSpread;
		if (vDate) delete vDate;
		NewNPV = Price(qCMPPRICE);
		Result = bp1 * Price(qCMPDEFLEGPV)/(Price(qCMPFEELEGPV)-Price(qCMPACCRUED));
		Feeleg->SetCreditSpread(InitialSpread);
		
	}else
	{
		ResetRootPricer();

		Feeleg->SetCreditSpread(bp);
		NewNPV = Price(qCMPPRICE);
		Result = bp * Price(qCMPDEFLEGPV)/(Price(qCMPFEELEGPV)-Price(qCMPACCRUED));
		Feeleg->SetCreditSpread(InitialSpread);
	}

	ResetRootPricer();
	SetSecurity(OriginalCds);
	if (cds) delete cds; cds=NULL;

	Result = fabs(Result*100.);	// result in BP ?
	SetSpread(Result);

	return (Result);
}

// *************************************************************
// Computing the Tranche's Duration
// *************************************************************

double ICM_Pricer_Security::ComputeDuration(void)
{
	// if (GetDurationFlg()) return GetDuration();

	ICM_Cds* cds = (ICM_Cds*) GetSecurity();
	ICM_ModelMultiCurves* model = (ICM_ModelMultiCurves*) GetModel();

	double Duration =0.,bp = 0.01,Fair_Spread=0.,spread=0.;
	double FeeLeg = 0.,Accrued = 0.,cds_NPV =0.,init_NPV=0.;

	// 14514 double notional = cds->GetInitialNotional()*cds->GetTradedCoef();
	double notional = cds->GetFeeLeg()->GetCreditInfosRef().GetNotionals().Elt(0)*cds->GetTradedCoef();

	if (cds->GetFeeLeg()->GetLegType() != K_FIXED_LEG || cds->GetFeeLeg()->GetName() == ICM_CORRIDORLEG)
	//if ((!cds->GetFeeLeg()->IsFixedLeg()) && cds->GetFeeLeg()->GetVariableSpread()) //variable leg
	{
		std::auto_ptr<ICM_Cds> Fixed_Cds ( (ICM_Cds*)cds->Clone()) ;
		// ICM_Cds* Fixed_Cds = (ICM_Cds*) cds->Clone();
		Fixed_Cds->GetFeeLeg()->SetCreditLegType(qRunning_Leg); //Force Fix Leg
		Fixed_Cds->GetFeeLeg()->SetRatesPricer(NULL);			//Force Pricer at NULL
		Fixed_Cds->SetCreditSpread(bp);							//Force Coupon Rate != 0
		Fixed_Cds->GetFeeLeg()->SetRefPartRate(&ARM_ReferenceValue(1)); 
		Fixed_Cds->GetFeeLeg()->GetCreditInfos()->SetFwdType(qCPN_FIXED); //Force Fixed coupon
		// Fixed_Cds->GetFeeLeg()->SetRefSpread(1.);				//Force Ref Spread if exist (CM)
		if (cds->GetFeeLeg()->GetName() == ICM_CORRIDORLEG) {
			Fixed_Cds->GetFeeLeg()->SetCreditLegType(qCorridor_Leg); //Force 		
			((ICM_CorridorLeg*)	Fixed_Cds->GetFeeLeg())->SetLeverageFloatIdx(0.);
			ARM_Date date = this->GetModel()->GetStartDate();
			ARM_Vector * vDate = new ARM_Vector(1);
			vDate->Elt(0) = date.GetJulian();
			ARM_Vector * vValue = new ARM_Vector(1);
			vValue->Elt(0) = 0.01;
			ARM_ReferenceValue* newSpread = new ARM_ReferenceValue();
			newSpread->SetDiscreteDates(vDate);
			newSpread->SetDiscreteValues(vValue);
			((ICM_CorridorLeg*)	Fixed_Cds->GetFeeLeg())->SetSpreads(newSpread);		
			if (newSpread) delete newSpread;
			if (vDate) delete vDate;
		}
		
		try 
		{
			SetSecurity(Fixed_Cds.get());							//New Security
			ResetRootPricer();										//Reset intermediate pricing parameters		
			SetInitialPriceFlg(false);								//Reset Initial NPV
			cds_NPV = Price(qCMPPRICE);										//Pricing 
			FeeLeg = Price(qCMPFEELEGPV);								//Feeleg Price
			Accrued = Price(qCMPACCRUED);								//accrued pv
			Duration = 100.*(FeeLeg-Accrued)/(notional*bp);				//Duration	
			SetSecurity(cds);										//Reset Initial Cds
			ResetRootPricer();										//Reset intermediate pricing parameters		
			SetInitialPriceFlg(false);								//Reset Initial NPV
		}
		catch(...)
		{
			SetSecurity(cds);										//Reset Initial Cds
			ResetRootPricer();										//Reset intermediate pricing parameters		
			SetInitialPriceFlg(false);								//Reset Initial NPV
			throw ;													// forward the exception
		}
		// if (Fixed_Cds) delete Fixed_Cds;Fixed_Cds=NULL;			//Free Memory
	}
	else if(cds->GetFeeLeg()->GetLegType() == K_FIXED_LEG)//Fixed Leg
	{
		spread = cds->GetCreditSpread();						//Get Initial Spread
		cds_NPV = Price(qCMPPRICE);										//Pricing 
		FeeLeg = Price(qCMPFEELEGPV);								//Feeleg Price
		Accrued = Price(qCMPACCRUED);								//accrued pv
		Fair_Spread = ComputeSpread(0.);						//comopute breakeven spread

		if(FeeLeg) Duration = (FeeLeg-Accrued)/(notional*(100.*spread)/10000.);
		else Duration = cds_NPV/(notional*(100.*spread - Fair_Spread)/10000.);
	}
	else if(cds->GetFeeLeg()->GetLegType() == K_FIXED_LEG)//Fixed Leg

	SetDuration(fabs(Duration));
	return fabs(Duration) ;
}


void ICM_Pricer_Security::View(char* id, FILE* ficOut)
{	
	FILE* fOut;
	char  fOutName[200];

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

	int size =0;

	if (itsIsTermStructurePricing) fprintf(fOut, "\n\tTermStructurePricing : YES \n\n"); else fprintf(fOut, "\n\tTermStructurePricing : NO \n\n");

	itsDistribLoss.View(id,fOut); 
	fprintf(fOut, "\n");
	
	ICM_Pricer::View(id, fOut);

	if (itsScheduleForEL.size()>0)
	{
	std::stringstream sstr ;
	sstr<<"\t\t\t ----------------- Losses /Yearterms /Unit ----------------- \n\n"<<std::endl ;
	sstr<<"\n"<<std::endl ;
	sstr<<"\t"<<std::endl ;

	for(int i1=0;i1<itsDetailLosses[0].size();i1++)
	{sstr<<i1<<"\t";}
	sstr<<"\n"<<std::endl ;

	for(int i2=0;i2<itsDetailLosses.size();i2++)
	{
		sstr<<"\t"<<((double)(int)(1.e3*itsScheduleForEL[i2]))/1.e3<<"\t";
		for(int i3=0;i3<itsDetailLosses[i2].size();i3++)
		{sstr<<((double)(int)(1.e3*itsDetailLosses[i2][i3]))/1.e3<<"\t";}		
		sstr<<std::endl ;
	}
	sstr<<"\n\n"<<std::endl ;
	fprintf(fOut,"%s\n",sstr.str().c_str()); 
	}

	if ( ficOut == NULL )
	{
		fclose(fOut);
	}
}

// --------------------------------------------
// Compute Forward Coupons for CM Tranche
// --------------------------------------------
double ICM_Pricer_Security::ComputeImpliedPartRate()
{
	if (GetSpreadFlg()) return GetSpread();
	ICM_Cds* OldCds = (ICM_Cds*) GetSecurity();
	std::auto_ptr<ICM_Cds> cds( dynamic_cast<ICM_Cds*> ( OldCds->Clone() ) ) ;
	ICM_Leg* Feeleg = cds->GetFeeLeg();
	// SetFeeLegPriceFlg(false);
	unsetFlg(qCMPFEELEGPV); 
	try {
		SetSecurity(cds.get());	
		// Feeleg->SetRefSpread(1.);
		ResetRootPricer(); 
		Feeleg->SetRefPartRate(&ARM_ReferenceValue(1.));
		double Result = DefLegPV()/FeeLegPV();	
		Result = fabs(Result);
		// SetFeeLegPriceFlg(false);	
		unsetFlg(qCMPFEELEGPV) ;
		SetSecurity(OldCds);	
		ResetRootPricer();
		SetSpread(Result);
		return (Result);
	}
	// in case of exception we restore the state of the pricer.. 
	catch(...)
	{
		// SetFeeLegPriceFlg(false);	
		unsetFlg(qCMPFEELEGPV); 
		SetSecurity(OldCds);	
		ResetRootPricer();
		throw; // and we forward the exception to the upper level. 
	}
	ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_Pricer_Security::ComputeImpliedPartRate: Should not occur."); 
}


// --------------------------------------------
double ICM_Pricer_Security::Compute_Fwd_Spread(const ARM_Date &Mty, const ARM_Date &CDS_ExpiryDate, double& dur)
{
	ICMTHROW(ERR_INVALID_ARGUMENT,"Not implemented"); 
}
// --------------------------------------------
double ICM_Pricer_Security::ComputeImpliedVol(const double& Price)
{
	ICMTHROW(ERR_INVALID_ARGUMENT,"Not implemented"); 
}
