#include "ARMKernel/glob/firsttoinc.h" 
#include "icm_distriblosscalculator.h" 

#include "ICMKernel\glob\icm_correlation_sector.h"
#include "ICMKernel\mod\icm_Homogene_FastLossDistrib_Gaussian_2F.h"
#include "ICMKernel\mod\icm_Homogene_FastLossDistrib_Gaussian.h"
#include "ICMKernel\mod\modelmulticurves.h"
#include "ICMKernel\pricer\icm_pricer_homogeneous_smile_collat_fwd.h"

#include "ICMKernel\mod\icm_Homogene_FastLossDistrib_Gaussian_MF.h"

#include "ICMKernel\inst\icm_mez.h"
#include "ICMKernel\pricer\icm_pricer_cds.h"
#include "ICMKernel\util\icm_rootfinder1D.h"
#include "ICMKernel\util\icm_RootFinderND.h"
//#include "ICMKernel\util\icm_brentsolver.h"
#include "ICMKernel\inst\icm_collateral.h"
#include "ICMKernel\mod\icm_lossunits.h"
#include "ICMKernel\inst\icm_credit_index.h"



// ---------------------------------------------------------------
// Fast Hedge loss distribution with beta correlation computation 
// ---------------------------------------------------------------
/**
17783 
double 
cpt_elt_pricer_distrib_fast(const double& yf,ICM_Pricer& pricer_, ARM_Model& model_,ARM_Security& sec_)
// double ICM_Pricer_Distrib::ExpectedLossTrancheForFastSpreadHedge(const double& yearterm)
{
	ICM_Pricer_Distrib& pricer = dynamic_cast<ICM_Pricer_Distrib&>(pricer_) ; 
	ICM_ModelMultiCurves& model = dynamic_cast<ICM_ModelMultiCurves&>(model_) ; 
	ICM_Ftd& ftd = dynamic_cast<ICM_Ftd&>(sec_);
	ICM_Collateral*collat =ftd.GetCollateral(); 
const ICM_LossUnits& lossUnits= pricer.getLossUnits(); 

// 17783 	if (!pricer.IsActivateShift())
// 17783 		ICMTHROW(ERR_INVALID_ARGUMENT,"ERROR : Mode Hedge expected in 'ExpectedLossTrancheForFastSpreadHedge'"); 

	//ICM_DefaultCurve* defaultcurve = model.GetDefaultCurves()[0]; 

	ICM_Security* security = ftd.GetDefLeg()->GetCreditInfos();
	
	//JLA 14514 int numflows = security->GetNumFlows();
	int nbnames = collat->GetNbIssuers();
	
	ARM_Vector* Betas = pricer.GetBetas();
	// double* beta = new double[nbnames];
	// memcpy(beta,Betas->GetElt(),sizeof(double)*nbnames);
	// double* Pdefault = NULL;

	int i, j, k;

	ICM_QMatrix<double>* ShiftMatrix = NULL;

	ARM_Date date = (ARM_Date)(model.GetStartDate().GetJulian() + K_YEAR_LEN*yf);
	double tranche_down = pricer.GetTranche_Down(date);
	double tranche_up = pricer.GetTranche_Up(date);

	double ELT=0.0;

	
	j	=	pricer.GetTenorShift();
	i	=	pricer.GetIssuerShift();
	// For Fast Spread Sensisitivity
	// first computation for this date yearterm
// 17783 		if (pricer.GetFirstComputationForHedgeFlg())
// 17783 		{
// 17783 			// get Perturb Probas from SensiManager
// 17783 			ShiftMatrix = pricer.PerturbProbasFromSensiManager(yf);
// 17783 	
// 17783 		}
// 17783 		else
	{
		// already computed

		// which index from yearterm in order to get the correct ShiftMatrix in the array 'LossTShifts'
		k	=	ftd.GetYearTermSearchInSchedule(yf);
		ELT	=	pricer.GetShifts(k)->Getvalue(i,j);
		if (ELT<=0)
			ELT =0.;
		else
			ELT /= (tranche_up-tranche_down);
		
		return ELT;
	}


	ARM_Vector Pdefault(nbnames) ; 
	const std::vector<int> &collatRank = lossUnits.getCollatRank(date); 
	for (i=0; i<nbnames; i++) 
		Pdefault[i] = model.GetDefaultCurve(collat->GetIssuersLabels(collatRank[i]))->DefaultProba(yf);

	ICM_Gauss1FLossDistrib Distrib(nbnames,Pdefault,*Betas,lossUnits.getLossRates(date));
	// if (Pdefault) delete[] Pdefault;

	Distrib.setPercentIndep(model.GetIndependantPart());
	Distrib.setPercentFullCorrel(model.GetFullCorrelPart());

	tranche_down = ABS(round(tranche_down));
	tranche_up = ABS(round(tranche_up));

	Distrib.compute_expectedlosstranche_fast_spread_hedge(tranche_up,
													tranche_down,
													lossUnits.getLossUnit(date),
													ShiftMatrix
													);
	// -----------------------------------
	// push the ShiftMatrix in the vector (Clone???)
	pricer.AppendShifts((ICM_QMatrix<double> *)Distrib.GetShifts()->Clone());
	// -----------------------------------

	// if (beta) delete[] beta;

	if (pricer.GetFirstComputationForHedgeFlg())
	{
		// which index from yearterm in order to get the correct ShiftMatrix in the array 'LossTShifts'
		k	=	ftd.GetYearTermSearchInSchedule(yf);

		// remark: i=j=0
		ELT	=	pricer.GetShifts(k)->Getvalue(0, 0);
		if (ELT<=0)
			ELT =0.;
		else
			ELT /= (tranche_up-tranche_down);
	}

	return (ELT);
}
// 17783 **/ 

// ---------------------------------------------------------------
// Fast Hedge loss distribution with base correlation computation 
// ---------------------------------------------------------------
/**  17783 
double 
cpt_elt_pricer_distrib_smile_fast(const double& yf,ICM_Pricer& pricer_, ARM_Model& model_,ARM_Security& sec_)
// double ICM_Pricer_Distrib_Smile::ExpectedLossTrancheForFastSpreadHedge(const double& yearterm)
{
	ICM_Pricer_Distrib_Smile& pricer = dynamic_cast<ICM_Pricer_Distrib_Smile&>(pricer_) ; 
	ICM_ModelMultiCurves& model = dynamic_cast<ICM_ModelMultiCurves&>(model_) ; 
	ICM_Ftd& ftd = dynamic_cast<ICM_Ftd&>(sec_);
	ICM_Collateral*collat =ftd.GetCollateral(); 
	const ICM_LossUnits& lossUnits= pricer.getLossUnits(); 

// 17783 	if (!pricer.IsActivateShift())
// 17783 		ICMTHROW(ERR_INVALID_ARGUMENT,"ERROR : Mode Hedge expected in 'ExpectedLossTrancheForFastSpreadHedge'"); 

	ARM_Date Maturity = ftd.GetEndDateNA();
	double YF_Maturity = (Maturity-model.GetStartDate())/365.;

	int nbnames = collat->GetNbIssuers();
	
	//ARM_Vector* Betas = GetBetas();
	// double* beta = new double[nbnames];
	// memset(beta,'\0',sizeof(double)*nbnames);
	// double* Pdefault = NULL;
	ARM_Vector beta(nbnames,0.); 

	int i, j, k;
 
	ICM_QMatrix<double>* ShiftMatrix = NULL;

	ARM_Date date = (ARM_Date)(model.GetStartDate().GetJulian() + K_YEAR_LEN*yf);
	double tranche_down = pricer.GetTranche_Down(date);
	double tranche_up = pricer.GetTranche_Up(date);

	tranche_down = ABS(round(tranche_down));
	tranche_up = ABS(round(tranche_up));

	double ELT=0.0;

	// if (!IsActivateShift())
	// {
	// 	// error
	// 	throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
    //         "ERROR : Mode Hedge expected in 'ExpectedLossTrancheForFastSpreadHedge'");
	// }
	
	j	=	pricer.GetTenorShift();
	i	=	pricer.GetIssuerShift();
	// For Fast Spread Sensisitivity
	// first computation for this date yearterm
// 17783 		if (pricer.GetFirstComputationForHedgeFlg())
// 17783 		{
// 17783 			// get Perturb Probas from SensiManager
// 17783 			ShiftMatrix = pricer.PerturbProbasFromSensiManager(yf);
// 17783 		}
// 17783 		else
	{
		// already computed

		// which index from yearterm in order to get the correct ShiftMatrix in the array 'LossTShifts'
		k	=	ftd.GetYearTermSearchInSchedule(yf) - 1;
		k	=	ftd.GetYearTermSearchInSchedule(yf)-1;
		// index 0 is skipped

		ELT	=	pricer.GetShifts(k)->Getvalue(i,j);
		if (ELT<=0)
			ELT =0.;
		else
			ELT /= (tranche_up-tranche_down);
		
		return ELT;
	}

	qIntegratorChoice	TheIntegratorType;
	// if (! itsIntegrationStep1) 
	if (!pricer.GetIntegrationStep1()) 
	{
		TheIntegratorType =qGAUSS_HERMITE;
		pricer.SetIntegrationStep1(40) ; 
		// itsIntegrationStep1 = 40;
	}
	else	
		TheIntegratorType	=	pricer.GetGaussLegendreMethodFromIntegrationStep(pricer.GetIntegrationStep1());

	const std::vector<int> &collatRank = lossUnits.getCollatRank(date); 
	ARM_Vector Pdefault (nbnames); 
	for (i=0; i<nbnames; i++) 
		Pdefault[i] = model.GetDefaultCurve(collat->GetIssuersLabels(collatRank[i]))->DefaultProba(yf);

	ICM_Gauss1FLossDistrib_LJ Distrib(nbnames,Pdefault,beta,
		lossUnits.getLossRates(date), // its_LossRate,
		pricer.GetIntegrationStep1(), // itsIntegrationStep1,
		pricer.GetCopulaType(), // itsCopulaType,
		TheIntegratorType,
		pricer.GetIntegrationStep1() // itsIntegrationStep1
		);

	//if (Pdefault) delete[] Pdefault;

  	// Smile Parameters
	ICM_Correlation* correlation = model.GetCorrelation();
	if (!pricer.GetFaster()) correlation->ComputeStrikesEq(&pricer,pricer.GetRescalType());

	ARM_Vector corr_low(nbnames);
	ARM_Vector corr_Hight(nbnames);

	//Reconstruction du vecteur des betas
	for (i=0; i<nbnames; i++) 
	{ 
		correlation->SetForcedStrikeType(qStrike_LOW);
		corr_low[i] = correlation->GetBeta(collat->GetIssuersLabels(collatRank[i]),YF_Maturity,CREDIT_DEFAULT_VALUE,yf);
		correlation->SetForcedStrikeType(qStrike_UP);
		corr_Hight[i] = correlation->GetBeta(collat->GetIssuersLabels(collatRank[i]),YF_Maturity,CREDIT_DEFAULT_VALUE,yf);
	}

	correlation->SetForcedStrikeType(qStrike_NONE);
	Distrib.UpdateCorrelation(corr_low,corr_Hight);

	// if (beta) delete[] beta;

	// New...
	Distrib.SetCopulaType(qGAUSSIAN);
	Distrib.SetIntegratorType(TheIntegratorType,pricer.GetIntegrationStep1());
	Distrib.SetIntegrationStep(pricer.GetIntegrationStep1());
		
	double Notional = collat->SumNotionals(ftd.GetStartDateNA());

	if (Notional)
	{

		Distrib.compute_expectedlosstranche_fast_spread_hedge(tranche_up,
														tranche_down,
														lossUnits.getLossUnit(date),
														ShiftMatrix
														);

		// -----------------------------------
		// push the ShiftMatrix in the vector (Clone???)
		pricer.AppendShifts((ICM_QMatrix<double> *)Distrib.GetShifts()->Clone());
		// -----------------------------------

	}
	else
	{	int	nbscenario	=	ShiftMatrix->Getnbrows();

		ICM_QMatrix<double>* DummyShift = new ICM_QMatrix<double>(nbnames, nbscenario, 0.0);

		// -----------------------------------
		// push the ShiftMatrix in the vector (Clone???)
		pricer.AppendShifts((ICM_QMatrix<double> *)DummyShift->Clone());
		// -----------------------------------
	}

	if (pricer.GetFirstComputationForHedgeFlg())
	{
		// which index from yearterm in order to get the correct ShiftMatrix in the array 'LossTShifts'
		k	=	ftd.GetYearTermSearchInSchedule(yf) - 1;
		// -1, because 0.0 is skipped
		k	=	ftd.GetYearTermSearchInSchedule(yf)-1;

		if (pricer.GetLossTShifts().size()<=k)
			ICMTHROW(ERR_INVALID_MODEL,"Parameters :  FAST HEDGE PROBLEM"); 
// 			throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
// 				"Parameters :  FAST HEDGE PROBLEM");

		// remark: i=j=0
		ELT	=	pricer.GetShifts(k)->Getvalue(0, 0);
		if (ELT<=0)
			ELT =0.;
		else
			ELT /= (tranche_up-tranche_down);
	}

	return (ELT);
}
17783 
**/ 

// --------------------------------------------------------------------
// tranche expected loss computation
// --------------------------------------------------------------------
double cpt_quick_el_full_homog(double pdefault, double recovery, int nbnames, double strikedw,
							   double strikeup, double correldw, double correlup, vector<double>& losses,
							   int intstep, bool LHP)
{

	ICM_LightDistrib Distrib(nbnames);
	double lossunit = 1.-recovery; 
	
	double notunit = 10.e6;

	double ELT = 0.;

	if (!LHP)
	{
	if (strikedw>strikeup) return 0.;

	ELT =Distrib.ComputeEL_FullHomog(lossunit*notunit,strikedw*notunit*nbnames,
									 strikeup*notunit*nbnames,sqrt(correldw), sqrt(correlup),
									 pdefault,intstep,losses);

	ELT = ELT/((strikeup-strikedw)*notunit*nbnames);

	}
	else
	{
	ELT =Distrib.ComputeEL_LHP(strikedw, strikeup, sqrt(correldw), sqrt(correlup), pdefault,
							   recovery, losses);
	}

	if (ELT<0.) return 0.;

	return ELT;

}

// --------------------------------------------------------------------
// Quick CDO pricing with dates vectors (for single correlation)
// --------------------------------------------------------------------
double FastCDOPricing(ARM_Vector* YFstartdates, ARM_Vector* YFenddates,  double rate,
					  double strikedw, double strikeup, double correldw, double correlup,
					  double notfeeleg,double notdefleg, vector<double>& pdef,
					  vector<double>& disc,double recovery,int nbnames,int intstep,
					  double& feepv,double& defpv, bool LHP)
{
	vector<double> el;
	int nbperiods = YFstartdates->GetSize();
	el.resize(nbperiods+1);

	int i=0;
	vector<double> losses;

	feepv = defpv = 0.;
	double indivnot = 10.e6;
	double notcdo = (strikeup-strikedw)*nbnames*indivnot;

	for (i=0;i<nbperiods;i++)
	{
		el[i]=0.;

		if (YFstartdates->Elt(i)>0.)
		{
		el[i] = cpt_quick_el_full_homog(pdef[i], recovery, nbnames, strikedw, strikeup,
							   correldw, correlup, losses, intstep, LHP);
		}
	}

	el[nbperiods] = cpt_quick_el_full_homog(pdef[nbperiods], recovery, nbnames,
							   strikedw, strikeup, correldw,correlup, losses, intstep, LHP);

	for (i=0;i<nbperiods;i++)
	{
		double begin = YFstartdates->Elt(i);
		double end = YFenddates->Elt(i);

		if (end<=0.)
		{	continue; }
		else if ((end>0.) && (begin<0.))
		{	begin=0.; }

		double deltaT = (365./360.)*(end-begin);
		feepv += deltaT*rate*notfeeleg*(1.-el[i+1])*disc[i+1];
		defpv += notdefleg*(-el[i]+el[i+1])*disc[i+1];
	}

	return (feepv-defpv);

}


// --------------------------------------------------------------------
// Quick CDO pricing 
// --------------------------------------------------------------------
double FastCDOPricing(ARM_Vector& schedule, const int& begin1,const int& end1,  double rate,
					  const vector<double>& el,double notfeeleg,double notdefleg, 
					  ARM_ZeroCurve* zc, double& feepv,double& defpv)
{
	int beg_ = MAX(0,begin1);
	int end_ = MIN(schedule.GetSize()-1,end1);

	feepv = defpv = 0.;

	if (begin1==end1)
	{return 0.;}

	for (int i=beg_;i<end_;i++)
	{
		double begin = schedule.Elt(i);
		double end = schedule.Elt(i+1);

		if (end<=0.)
		{	continue; }
		else if ((end>0.) && (begin<0.))
		{	begin=0.; }

		double disc = zc->ForwardPrice(schedule.Elt(beg_),end);

		double deltaT = (365./360.)*(end-begin);
		feepv += deltaT*rate*notfeeleg*(1.-el[i+1])*disc;
		defpv += notdefleg*(-el[i]+el[i+1])*disc;
	}

	return (feepv-defpv);

}

// --------------------------------------------------------------------
// Quick CDO pricing with start date & end date (for single correlation)
// --------------------------------------------------------------------
double FastCDOPricing(ARM_Date& AsOf, ARM_Date& startdate, ARM_Date& enddate,int frequency,
					  double rate, double strikedw, double strikeup,  double correldw,
					  double correlup, double notfeeleg, double notdefleg,  ICM_DefaultCurve* pdef,
					  ARM_ZeroCurve* disc, double recovery, int nbnames, int intstep,
					  double& feepv, double& defpv, bool LHP)
{
	vector<double> Vpdef;
	vector<double> Vdisc;

	ARM_Vector YFstartdates=NULL;
	ARM_Vector YFenddates=NULL;

	ARM_Vector* startdates=NULL;
	ARM_Vector* enddates=NULL;

	int i=0;

	startdates = CptStartDates(startdate,enddate,frequency,
                                            K_MOD_FOLLOWING, K_LONGSTART, K_MATUNADJUSTED,
                                            "EUR",0);

	enddates = new ARM_Vector(startdates);

	YFstartdates.Resize(startdates->size());
	YFenddates.Resize(startdates->size());

	for (i=0;i<enddates->size()-1;i++)
	{enddates->Elt(i)=enddates->Elt(i+1);}

	enddates->Elt(enddates->size()-1)=enddate.GetJulian();

	for (i=0;i<enddates->size();i++)
	{
		YFstartdates.Elt(i)=(startdates->Elt(i)-AsOf.GetJulian())/365.;
		if (YFstartdates.Elt(i)<0.) {YFstartdates.Elt(i)=0.;}
		YFenddates.Elt(i)=(enddates->Elt(i)-AsOf.GetJulian())/365.;
		if (YFenddates.Elt(i)<0.) {YFenddates.Elt(i)=0.;}
	}

	for (i=0;i<YFstartdates.GetSize();i++)
	{
		Vdisc.push_back(disc->DiscountPrice(YFstartdates.Elt(i)));
		Vpdef.push_back(pdef->DefaultProba(YFstartdates.Elt(i)));
	}
	Vdisc.push_back(disc->DiscountPrice(YFenddates.Elt(YFstartdates.GetSize()-1)));
	Vpdef.push_back(pdef->DefaultProba(YFenddates.Elt(YFstartdates.GetSize()-1)));


	double result = FastCDOPricing(&YFstartdates,&YFenddates,rate,strikedw,strikeup,correldw,correlup,
				   notfeeleg,notdefleg,Vpdef,Vdisc,recovery,nbnames,intstep,feepv,defpv,LHP);

	if (startdates) delete startdates;
	if (enddates) delete enddates;

	return result;
}


// --------------------------------------------------------------------
// Portfolio Expected loss computation
// --------------------------------------------------------------------
double CptPtf_ELoss(ICM_Pricer* pricer,int payfreq)
{
	ICM_Mez* cdo = (ICM_Mez*) pricer->GetSecurity();
	ICM_Collateral * collat = cdo->GetCollateral(); 
	double eloss_bespoke = 0.,eloss_indice = 0.;
	ARM_Date AsOf = pricer->GetModel()->GetStartDate();
	//ARM_Date Start = AsOf; Start.AddDays(1);
	ARM_Date TrancheStart = ((ARM_SwapLeg*)(cdo->GetFeeLeg()))->GetStartDate();
	ARM_Date Maturity = cdo->GetFeeLeg()->GetMaturity();
	int bespoke_nbnames = collat->GetNbIssuers();
	double recovery_ratio = 0., poids = 0., nominal_total =0.;

	//Eloss du portefeuille : sommes des eloss ind
	ICM_ModelMultiCurves* mod = (ICM_ModelMultiCurves*) pricer->GetModel();
	//Nominal total du portefeuille
	nominal_total = cdo->GetCollateral()->SumNotionals(cdo->GetStartDateNA());

	//const std::vector<int> &collatRank = lossUnits.getCollatRank(date); 
	for (int i=0; i<bespoke_nbnames; i++) 
	{
		//Creation du modèle
		ICM_DefaultCurveModel DefModel; 
		const ICM_DefaultCurve* CDSDefCurve =  mod->GetDefaultCurve(collat->GetIssuersLabels(i));
		// bool calibtype = CDSDefCurve->GetIsSummitCurve();
		ARM_ZeroCurve* CDSRateCurve =	mod->GetDiscountCurve();
		if (payfreq==K_ZEROCOUPON)
		{
			//JLA added:
			ICM_DefaultCurve* summitCurve = dyn_clone(CDSDefCurve); 
			summitCurve->SetIsSummitCurve(false);
			DefModel.Set(summitCurve, CDSRateCurve); 
			delete summitCurve; 
		}
		else 
		{
			DefModel.Set(CDSDefCurve, CDSRateCurve); 
		}
		
		
		string ccy = CDSDefCurve->GetCurrency();
		//Creation du CDS (fonction du pricer)
		switch (pricer->GetName())
		{
		case ICM_PRICER_HOMOGENEOUS_SMILE:
			{
				ICM_Cds CDS(TrancheStart,Maturity,(ARM_Date*)0,// stub 
					0,// first cpn 
					TrancheStart,Maturity,
					0., // rate 
					ARM_ReferenceValue(100),
					ARM_ReferenceValue(100),// 0,0, // notionals 
					payfreq, KACTUAL_360, qACCRUED_SETTLED, ccy/*ARM_DEFAULT_COUNTRY*/,
					cdo->GetFeeLeg()->GetStubMeth(),// const int& stubrule  
					cdo->GetDefLeg()->GetCreditLag(), // const double& CreditLag 
					payfreq, // const int& FrequencyDefLeg 
					cdo->GetFeeLeg()->GetIRIndex()->GetIntRule(), // const int& intRule  
					cdo->GetFeeLeg()->GetIncludeMaturity(), // const bool& includematurity  
					cdo->GetDefLeg()->GetAdjStartDateFlag(), // const int& adjStartDate 
					std::string(), // const std::string& payCalName  
					qRunning_Leg, // const qCredit_Leg_Type& TypeFeeLeg 
					qStandart_Recovery_Leg, // const qCredit_Leg_Type& TypeDefLeg 
					ISSUER_UNDEFINE, // const string& name  
					CREDIT_DEFAULT_VALUE // const double& Binary  
					);
				//Mis à jour du recovery si Recov Fixe
				if (cdo->GetBinaryFlg() == true)
					recovery_ratio = (1-cdo->GetBinary())/(1-CDSDefCurve->GetRecovery());
				else if (cdo->GetCollateral()->GetRecovCoefFlg())
					recovery_ratio = (1-MIN(CDSDefCurve->GetRecovery()*cdo->GetCollateral()->GetRecovCoef(),1))/(1-CDSDefCurve->GetRecovery());
				else
					recovery_ratio = 1.;

				//Pricer
				// ICM_Pricer_Cds PricerCds(&CDS,&DefModel);
				ICM_Pricer_Cds PricerCds; PricerCds.Set(&CDS,&DefModel,ICM_Parameters(),AsOf);
				
				//ELoss (exprimée en %)
				// poids = collat->GetIssuersNotional(collat->GetIssuersLabels(i)) / nominal_total;
				poids = collat->GetIssuersNotional(i) / nominal_total;
				// eloss_bespoke += PricerCds.DefLegPV() * poids * recovery_ratio / 100;
				eloss_bespoke += PricerCds.Price(qCMPDEFLEGPV) * poids * recovery_ratio / 100;
				break;
			}
		case ICM_PRICER_HOMOGENEOUS_SMILE_COLLAT_FWD:
			{
				ICM_Cds CDS(TrancheStart,Maturity,(ARM_Date*)0, // stub 
					0, // fst coupon 
					TrancheStart,Maturity,
					0., // rate , 
					ARM_ReferenceValue(100),
					ARM_ReferenceValue(100),// 0,0,// notionals 
					payfreq, KACTUAL_360, qACCRUED_SETTLED, ccy/*ARM_DEFAULT_COUNTRY*/ ,
						cdo->GetFeeLeg()->GetStubMeth(),// const int& stubrule  
						cdo->GetDefLeg()->GetCreditLag(), // const double& CreditLag 
						payfreq, // const int& FrequencyDefLeg 
						cdo->GetFeeLeg()->GetIRIndex()->GetIntRule(), // const int& intRule  
						cdo->GetFeeLeg()->GetIncludeMaturity(), // const bool& includematurity  
						cdo->GetDefLeg()->GetAdjStartDateFlag(), // const int& adjStartDate 
						std::string(), // const std::string& payCalName  
						qRunning_Leg, // const qCredit_Leg_Type& TypeFeeLeg 
						qStandart_Recovery_Leg, // const qCredit_Leg_Type& TypeDefLeg 
						ISSUER_UNDEFINE, // const string& name  
						CREDIT_DEFAULT_VALUE // const double& Binary  
			);
				//Mis à jour du recovery si Recov Fixe
				if (cdo->GetBinaryFlg() == true)
					recovery_ratio = (1-cdo->GetBinary())/(1-CDSDefCurve->GetRecovery());
				else if (cdo->GetCollateral()->GetRecovCoefFlg())
					recovery_ratio = (1-MIN(CDSDefCurve->GetRecovery()*cdo->GetCollateral()->GetRecovCoef(),1))/(1-CDSDefCurve->GetRecovery());
				else
					recovery_ratio = 1.;

				//Pricer
				// ICM_Pricer_Cds PricerCds(&CDS,&DefModel);
				ICM_Pricer_Cds PricerCds; PricerCds.Set(&CDS,&DefModel,ICM_Parameters(),AsOf);
				
				//ELoss (exprimée en %)
				// poids = collat->GetIssuersNotional(collat->GetIssuersLabels(i)) / nominal_total;
				poids = collat->GetIssuersNotional(i) / nominal_total;
				// eloss_bespoke += PricerCds.DefLegPV() * poids * recovery_ratio / 100;
				eloss_bespoke += PricerCds.Price(qCMPDEFLEGPV) * poids * recovery_ratio / 100;
				break;
			}
		default :
			{
				ICM_Cds CDS(TrancheStart,Maturity,NULL, // stub 
					0, // fist cpn 
					TrancheStart,Maturity,
					0., // rate 
					ARM_ReferenceValue(100),
					ARM_ReferenceValue(100),
					// ,0,0, // notionals 
					payfreq, KACTUAL_360,  qACCRUED_SETTLED, ccy/*ARM_DEFAULT_COUNTRY*/ ,
						cdo->GetFeeLeg()->GetStubMeth(),// const int& stubrule  
						cdo->GetDefLeg()->GetCreditLag(), // const double& CreditLag 
						payfreq, // const int& FrequencyDefLeg 
						cdo->GetFeeLeg()->GetIRIndex()->GetIntRule(), // const int& intRule  
						cdo->GetFeeLeg()->GetIncludeMaturity(), // const bool& includematurity  
						cdo->GetDefLeg()->GetAdjStartDateFlag(), // const int& adjStartDate 
						std::string(), // const std::string& payCalName  
						qRunning_Leg, // const qCredit_Leg_Type& TypeFeeLeg 
						qStandart_Recovery_Leg, // const qCredit_Leg_Type& TypeDefLeg 
						ISSUER_UNDEFINE, // const string& name  
						CREDIT_DEFAULT_VALUE // const double& Binary  
					);
				//Mis à jour du recovery si Recov Fixe
				if (cdo->GetBinaryFlg() == true)
					recovery_ratio = (1-cdo->GetBinary())/(1-CDSDefCurve->GetRecovery());
				else if (cdo->GetCollateral()->GetRecovCoefFlg())
					recovery_ratio = (1-MIN(CDSDefCurve->GetRecovery()*cdo->GetCollateral()->GetRecovCoef(),1))/(1-CDSDefCurve->GetRecovery());
				else
					recovery_ratio = 1.;

				//Pricer
				// ICM_Pricer_Cds PricerCds(&CDS,&DefModel);
				ICM_Pricer_Cds PricerCds; PricerCds.Set(&CDS,&DefModel,ICM_Parameters(),AsOf);
				
				//ELoss (exprimée en %)
				// poids = collat->GetIssuersNotional(collat->GetIssuersLabels(i)) / nominal_total;
				poids = collat->GetIssuersNotional(i) / nominal_total;
				// eloss_bespoke += PricerCds.DefLegPV() * poids * recovery_ratio / 100;
				eloss_bespoke += PricerCds.Price(qCMPDEFLEGPV) * poids * recovery_ratio / 100;
				break;
			}
		}	
		
		// if (payfreq==K_ZEROCOUPON) delete CDSDefCurve; 
		// if (payfreq==K_ZEROCOUPON)
		// {CDSDefCurve->SetIsSummitCurve(calibtype);}
	}

	//if (cdo) delete cdo;
	return (eloss_bespoke);
}

// --------------------------------------------------------------------
// Portfolio Expected loss computation
// --------------------------------------------------------------------
double CptPtf_ELoss_Index(ICM_ModelMultiCurves* mod,
					ICM_Credit_Index* index ,
					ARM_Date& Maturity,
					int payfreq)
{
	double eloss_bespoke = 0.;
	const vector<string> indexlabels=index->GetLabels();
	int idx_nbnames = indexlabels.size();

	if (index->IsHomogeneous())
	{ return CptPtf_ELoss_Index(mod,indexlabels[0],Maturity,payfreq);}

	double poids = 0., nominal_total =0.;

	//Nominal total du portefeuille
	nominal_total = idx_nbnames;

	for (int i=0; i<idx_nbnames; i++) 
	{
		//Creation du modèle
		const ICM_DefaultCurve* CDSDefCurve =  mod->GetDefaultCurve(indexlabels[i]);
		
		/** 
		//  JLA: useless here
		//	done in CptPtf_ELoss_Index
		// 
		bool calibtype = CDSDefCurve->GetIsSummitCurve();
		if (payfreq==K_ZEROCOUPON)
		{CDSDefCurve->SetIsSummitCurve(false);}
		**/ 

		//Creation du CDS (fonction du pricer)
		double defpv = CptPtf_ELoss_Index(mod,indexlabels[i],Maturity,payfreq);				

		//ELoss (exprimée en %)
		poids = 1. / nominal_total;
		eloss_bespoke += defpv * poids ;
	}

	//if (cdo) delete cdo;
	return (eloss_bespoke);
}

// --------------------------------------------------------------------
// Portfolio Expected loss computation for index
// --------------------------------------------------------------------
double CptPtf_ELoss_Index(ICM_ModelMultiCurves* mod,
						  const std::string& indexlabel,
						ARM_Date& Maturity,
						int payfreq)
{

	double eloss_indice = 0.;
	ARM_Date AsOf = mod->GetStartDate();
	ARM_Date Start = AsOf; Start.AddDays(1);
	double recovery_ratio=1.;
	
	const ICM_DefaultCurve* CDSDefCurveIndex =		mod->GetDefaultCurve(indexlabel);
	ARM_ZeroCurve* CDSRateCurveIndex =		CDSDefCurveIndex->GetZeroCurve();
	ICM_DefaultCurveModel DefModelIndex ;
	//	(CDSDefCurveIndex,CDSRateCurveIndex);
	// bool calibtype = CDSDefCurveIndex->GetIsSummitCurve();
	if (payfreq==K_ZEROCOUPON)
	{
		ICM_DefaultCurve * summitCurve= dyn_clone(CDSDefCurveIndex); 
		summitCurve->SetIsSummitCurve(false);
		// CDSDefCurveIndex->SetIsSummitCurve(false);
		DefModelIndex.Set(summitCurve,CDSRateCurveIndex); 
		delete summitCurve; 
	}
	else
	{
		DefModelIndex.Set(CDSDefCurveIndex,CDSRateCurveIndex);
	}

	string ccy = CDSDefCurveIndex->GetCurrency();
	//Creation du cds sur indice (base 100)
	ICM_Cds CDSIndex(Start,Maturity,(ARM_Date*)0, // stub 
					0, // fst cpn 
					Start,Maturity,
					0., // rate 
					ARM_ReferenceValue(100),
					ARM_ReferenceValue(100),// 0,0, // notionals 
					payfreq, KACTUAL_360, qACCRUED_SETTLED, ccy /*ARM_DEFAULT_COUNTRY*/,
					K_SHORTSTART,// const int& stubrule  
					DEFAULT_CREDIT_LAG_INDX, // const double& CreditLag 
					payfreq, // const int& FrequencyDefLeg 
					K_MATUNADJUSTED, // const int& intRule  
					INCLUDE_MATURITY, // const bool& includematurity  
					K_UNADJUSTED, // const int& adjStartDate 
					std::string(), // const std::string& payCalName  
					qRunning_Leg, // const qCredit_Leg_Type& TypeFeeLeg 
					qStandart_Recovery_Leg, // const qCredit_Leg_Type& TypeDefLeg 
					ISSUER_UNDEFINE, // const string& name  
					CREDIT_DEFAULT_VALUE // const double& Binary  
					);
	//Creation du modèle					
	// ICM_DefaultCurveModel DefModelIndex(CDSDefCurveIndex,CDSRateCurveIndex);
				
	//Pricer
	// ICM_Pricer_Cds PricerCdsIndex(&CDSIndex,&DefModelIndex);
	ICM_Pricer_Cds PricerCdsIndex; PricerCdsIndex.Set(&CDSIndex,&DefModelIndex,ICM_Parameters(),AsOf);

	//Estimation de l'eloss sur l'indice (exprimée en %)	
	// eloss_indice = PricerCdsIndex.DefLegPV()/100.;
	eloss_indice = PricerCdsIndex.Price(qCMPDEFLEGPV)/100.;

	// if (payfreq==K_ZEROCOUPON)
	// {CDSDefCurveIndex->SetIsSummitCurve(calibtype);}

	return (eloss_indice);			
}


// --------------------------------------------------------------------
// tranche expected loss computation Term Structure
// --------------------------------------------------------------------
double cpt_quick_el_full_homog_TSR(double T1,double T2,double pdefault_T1,double pdefault_T2, double recovery, int nbnames, double strikedw, double strikeup, 
								   double correldw_T1,double correldw_T2, double correlup_T1, double correlup_T2,
								   vector<double>& losses,int intstep, bool LHP)
{

	ICM_LightDistrib Distrib(nbnames);
	double lossunit = 1.-recovery; 
	
	double notunit = 10.e6;

	double ELT = 0.;

	if (!LHP)
	{
	if (strikedw>strikeup) return 0.;

	ELT =Distrib.ComputeEL_FullHomog_TSR(lossunit*notunit,strikedw*notunit*nbnames,
									 strikeup*notunit*nbnames,sqrt(correldw_T1),sqrt(correldw_T2), 
									 sqrt(correlup_T1),sqrt(correlup_T2),pdefault_T1,pdefault_T2,
									 /* to modify */ pdefault_T1,pdefault_T2,
									 T1,T2,intstep,losses);

	ELT = ELT/((strikeup-strikedw)*notunit*nbnames);

	}
	else
	{
	ELT =Distrib.ComputeEL_LHP(strikedw, strikeup, sqrt(correldw_T1), sqrt(correlup_T1), pdefault_T1,
							   recovery, losses);
	}

	if (ELT<0.) return 0.;

	return ELT;

}

// --------------------------------------------------------------------
// Quick CDO pricing with dates vectors (for single correlation) Term Structure
// --------------------------------------------------------------------
double FastCDOPricing_TSR(ARM_Vector* YFstartdates, ARM_Vector* YFenddates,  double rate,
					  double strikedw, double strikeup,ARM_VolCurve* correl,
					  double notfeeleg,double notdefleg,vector<double>& pdef_T1,vector<double>& pdef_T2,
					  vector<double>& disc,double recovery,int nbnames,int intstep,
					  double& feepv,double& defpv, bool LHP)
{
	vector<double> el;
	int nbperiods = YFstartdates->GetSize();
	el.resize(nbperiods+1);

	int i=0;
	vector<double> losses;

	vector<double> correldw_T1;
	vector<double> correlup_T1;
	vector<double> correldw_T2;
	vector<double> correlup_T2;

	correldw_T1.resize(nbperiods+1);
	correlup_T1.resize(nbperiods+1);
	correldw_T2.resize(nbperiods+1);
	correlup_T2.resize(nbperiods+1);

	feepv = defpv = 0.;
	double indivnot = 10.e6;
	double notcdo = (strikeup-strikedw)*nbnames*indivnot;

	//year fractions base correlation
	ARM_Vector BaseCorrelSchedule = *correl->GetExpiryTerms();

	double YF_start=0.;
	double YF_end=YFstartdates->Elt(0);
	double yf=0.,yf1=0.,yf2=0.;

	for (i=0;i<nbperiods;i++)
	{
		el[i]=0.;

		if (i){
		YF_start=YFstartdates->Elt(i-1);
		YF_end=YFstartdates->Elt(i);}

		yf = ((YF_start<=BaseCorrelSchedule.Elt(0)) ? FlatVectorInterpol(BaseCorrelSchedule,BaseCorrelSchedule,YF_start,false) : FlatVectorInterpol(BaseCorrelSchedule,BaseCorrelSchedule,YF_start,true));
		correldw_T1[i]=correl->ComputeVolatility(yf,100.*strikedw)/100.;
		correlup_T1[i]=correl->ComputeVolatility(yf,100.*strikeup)/100.;

		yf = ((YF_end<=BaseCorrelSchedule.Elt(0)) ? 0. : FlatVectorInterpol(BaseCorrelSchedule,BaseCorrelSchedule,YF_end,false));
		correldw_T2[i]=correl->ComputeVolatility(yf,100.*strikedw)/100.;
		correlup_T2[i]=correl->ComputeVolatility(yf,100.*strikeup)/100.;

		yf1 = ((YF_start<=BaseCorrelSchedule.Elt(0)) ? YF_start : FlatVectorInterpol(BaseCorrelSchedule,BaseCorrelSchedule,YF_start,true));
		yf2 = ((YF_end<=BaseCorrelSchedule.Elt(0)) ? 0 : YF_end);

		if (YFstartdates->Elt(i)>0.)
		{
		el[i] = cpt_quick_el_full_homog_TSR(yf1,yf2,pdef_T1[i],pdef_T2[i], recovery, nbnames, strikedw, strikeup,
							   correldw_T1[i],correldw_T2[i], correlup_T1[i],correlup_T2[i], losses, intstep, LHP);
		}
	}

	YF_start=YFstartdates->Elt(nbperiods-1);
	YF_end=YFenddates->Elt(nbperiods-1);

	yf = ((YF_start<=BaseCorrelSchedule.Elt(0)) ? FlatVectorInterpol(BaseCorrelSchedule,BaseCorrelSchedule,YF_start,false) : FlatVectorInterpol(BaseCorrelSchedule,BaseCorrelSchedule,YF_start,true));
	correldw_T1[nbperiods]=correl->ComputeVolatility(yf,100.*strikedw)/100.;
	correlup_T1[nbperiods]=correl->ComputeVolatility(yf,100.*strikeup)/100.;

	yf = ((YF_end<=BaseCorrelSchedule.Elt(0)) ? 0. : FlatVectorInterpol(BaseCorrelSchedule,BaseCorrelSchedule,YF_end,false));
	correldw_T2[nbperiods]=correl->ComputeVolatility(YF_end,100.*strikedw)/100.;
	correlup_T2[nbperiods]=correl->ComputeVolatility(YF_end,100.*strikeup)/100.;

	yf1 = ((YF_start<=BaseCorrelSchedule.Elt(0)) ? YF_start : FlatVectorInterpol(BaseCorrelSchedule,BaseCorrelSchedule,YF_start,true));
	yf2 = ((YF_end<=BaseCorrelSchedule.Elt(0)) ? 0 : YF_end);

	el[nbperiods] = cpt_quick_el_full_homog_TSR(yf1,yf2,pdef_T1[nbperiods],pdef_T2[nbperiods], recovery, nbnames,
							   strikedw, strikeup, correldw_T1[nbperiods],correldw_T2[nbperiods],correlup_T1[nbperiods],correlup_T2[nbperiods], losses, intstep, LHP);

	for (i=0;i<nbperiods;i++)
	{
		double begin = YFstartdates->Elt(i);
		double end = YFenddates->Elt(i);

		if (end<=0.)
		{	continue; }
		else if ((end>0.) && (begin<0.))
		{	begin=0.; }

		double deltaT = (365./360.)*(end-begin);
		feepv += deltaT*rate*notfeeleg*(1.-el[i+1])*disc[i+1];
		defpv += notdefleg*(-el[i]+el[i+1])*disc[i+1];
	}

	return (feepv-defpv);

}

// --------------------------------------------------------------------
// Quick CDO pricing with start date & end date (for single correlation)  Term Structure
// --------------------------------------------------------------------
double FastCDOPricing_TSR(ARM_Date& AsOf, ARM_Date& startdate, ARM_Date& enddate,int frequency,
					  double rate, double strikedw, double strikeup, double notfeeleg, double notdefleg,  ICM_DefaultCurve* pdef,
					  ARM_ZeroCurve* disc, double recovery, int nbnames,ARM_VolCurve* correl, int intstep,
					  double& feepv, double& defpv, bool LHP)
{
	vector<double> Vpdef_T1;
	vector<double> Vpdef_T2;
	vector<double> Vdisc;

	ARM_Vector YFstartdates=NULL;
	ARM_Vector YFenddates=NULL;

	ARM_Vector* startdates=NULL;
	ARM_Vector* enddates=NULL;

	int i=0;

	double correldw=0.;
	double correlup=0.;

	startdates = CptStartDates(startdate,enddate,frequency,
                                            K_MOD_FOLLOWING, K_LONGSTART, K_MATUNADJUSTED,
                                            "EUR",0);

	enddates = new ARM_Vector(startdates);

	YFstartdates.Resize(startdates->size());
	YFenddates.Resize(startdates->size());

	for (i=0;i<enddates->size()-1;i++)
	{enddates->Elt(i)=enddates->Elt(i+1);}

	enddates->Elt(enddates->size()-1)=enddate.GetJulian();

	//year fractions base correlation
	ARM_Vector BaseCorrelSchedule = *correl->GetExpiryTerms();
	double yf=0.;

	for (i=0;i<enddates->size();i++)
	{
		YFstartdates.Elt(i)=(startdates->Elt(i)-AsOf.GetJulian())/365.;
		if (YFstartdates.Elt(i)<0.) {YFstartdates.Elt(i)=0.;}
		YFenddates.Elt(i)=(enddates->Elt(i)-AsOf.GetJulian())/365.;
		if (YFenddates.Elt(i)<0.) {YFenddates.Elt(i)=0.;}
	}

	for (i=0;i<YFstartdates.GetSize();i++)
	{
		Vdisc.push_back(disc->DiscountPrice(YFstartdates.Elt(i)));

		if (YFstartdates.Elt(i)>BaseCorrelSchedule.Elt(0))
			{yf=FlatVectorInterpol(BaseCorrelSchedule,BaseCorrelSchedule,YFstartdates.Elt(i));}
		else
			{yf=YFstartdates.Elt(i);}

		Vpdef_T1.push_back(pdef->DefaultProba(yf));
		Vpdef_T2.push_back(pdef->DefaultProba(YFstartdates.Elt(i)));
	}
	Vdisc.push_back(disc->DiscountPrice(YFenddates.Elt(YFstartdates.GetSize()-1)));

	if (YFenddates.Elt(YFstartdates.GetSize()-1)>BaseCorrelSchedule.Elt(0))
		{yf=FlatVectorInterpol(BaseCorrelSchedule,BaseCorrelSchedule,YFenddates.Elt(YFstartdates.GetSize()-1));}
	else
		{yf=YFenddates.Elt(YFstartdates.GetSize()-1);}

	Vpdef_T1.push_back(pdef->DefaultProba(yf));
	Vpdef_T2.push_back(pdef->DefaultProba(YFenddates.Elt(YFstartdates.GetSize()-1)));


	double result = FastCDOPricing_TSR(&YFstartdates,&YFenddates,rate,strikedw,strikeup,correl,notfeeleg,notdefleg,
								Vpdef_T1,Vpdef_T2,Vdisc,recovery,nbnames,intstep,feepv,defpv,LHP);

	if (startdates) delete startdates;
	if (enddates) delete enddates;

	return result;
}