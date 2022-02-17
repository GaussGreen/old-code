#include "firsttoinc.h"
#include "ICMKernel/util/icm_utils.h"
#include "icm_pricer_tranche_restrikable.h"
#include "ICMKernel/inst/icm_option_tranche.h"
#include "ICMKernel/inst/icm_collateral.h"
#include "ICMKernel/inst/icm_mez.h"
#include "ICMKernel\crv\icm_defaultcurve.h"
#include "ICMKernel\crv\icm_constant_piecewise.h"
#include "ICMKernel\mod\modelmulticurves.h"
#include "ICMKernel\glob\icm_smile_correlation.h"
#include "ICMKernel\pricer\ICM_Pricer_homogeneous_smile.h"

void ICM_Pricer_Tranche_Restrikable :: Init(void)
	{
		ICM_Pricer::Init();
		
		itsNbSimul = 1;
		itsTimeStep=0.0001;

		itsSigma=0.0;
		itsMRS = 0.0;

		itsBeta0 = 0.;
		itsBeta1 = 0.;
		itsBeta2 = 0.;
		itsBeta3 = 0.;

		itsIsHetero = 0;

		itsRecov = 0.4;
		//itsIR = 0.0;
		
		itsDatesYF = NULL;
		itsprobas = NULL;

		itsInitialSpread = 0.0;
		itsInitialSpreadFlag = false;
		itsELCoeff=1. ;
		itsLossCoeff = 1. ;
		itsAvgSpreadFreq = 0.001 ; 
		itsUsePeriodIsEL = 1 ;
		itsFeeLegUnit = 0.0 ;
		itsFeeLegFlg = false ;

		SetName(ICM_PRICER_TRANCHE_RESTRIKABLE);
	}

void ICM_Pricer_Tranche_Restrikable :: Set(ARM_Security* sec , ARM_Model* mod ,const ICM_Parameters& params, const ARM_Date& asof)
{
	
	ICM_Pricer::Set(sec,(ARM_Object*)mod,params,&asof);
	
	// Simul params 

	ARM_Vector* Simuls  = params.GetColVect("NbSimulations");
	if(Simuls)
		itsNbSimul = Simuls->Elt(0);

	ARM_Vector* TimeStep = params.GetColVect("TimeStep");
	if(TimeStep)
		itsTimeStep = TimeStep->Elt(0);


	ARM_Vector* MRS = params.GetColVect("MRS");
	if ( MRS)
		itsMRS = MRS->Elt(0);

	ARM_Vector* Hetero = params.GetColVect("IsHeterogene");
	if ( Hetero)
		itsIsHetero = Hetero->Elt(0);

	ARM_Vector* ELCoeff = params.GetColVect("EL_Coeff");
	if ( ELCoeff)
		itsELCoeff = ELCoeff->Elt(0);

	ARM_Vector* LossCoeff = params.GetColVect("Loss_Coeff");
	if ( LossCoeff)
		itsLossCoeff = LossCoeff->Elt(0);

	ARM_Vector* AvrgFreq = params.GetColVect("AverageSpread_Freq");
	if(AvrgFreq)
		itsAvgSpreadFreq = AvrgFreq->Elt(0);

	ARM_Vector* UsePeriodIsEL = params.GetColVect("UsePeriodInEL");
	if(UsePeriodIsEL)
		itsUsePeriodIsEL = UsePeriodIsEL->Elt(0);

	/*****************************************************/
	
	ARM_Vector* Beta0 = params.GetColVect("Beta0");
	if ( Beta0)
		itsBeta0 =sqrt( Beta0->Elt(0));

	ARM_Vector* Beta1 = params.GetColVect("Beta1");
	if ( Beta1)
		itsBeta1 =sqrt( Beta1->Elt(0));

	ARM_Vector* Beta2 = params.GetColVect("Beta2");
	if ( Beta2)
		itsBeta2 =sqrt( Beta2->Elt(0));

	ARM_Vector* Beta3 = params.GetColVect("Beta3");
	if ( Beta3)
		itsBeta3 =sqrt( Beta3->Elt(0));


	


}

ICM_Pricer_Tranche_Restrikable::~ICM_Pricer_Tranche_Restrikable(void)
{
	if (itsDatesYF)
		delete[] itsDatesYF;
	itsDatesYF = NULL;
	
	if ( itsprobas)
		delete[] itsprobas;
	itsprobas = NULL;

}
void ICM_Pricer_Tranche_Restrikable::SetSigmaFromModel()
{
	ICM_ModelMultiCurves* model = (ICM_ModelMultiCurves *) GetModel() ;

	ARM_VolCurve* volcurve = model ->GetVolCurve();

	if(volcurve)
	{
		itsSigma = volcurve->GetVolatilities()->Elt(0,0);
	}
	else
	{
		itsSigma=0.0;
	}

}
void ICM_Pricer_Tranche_Restrikable::ComputeBetas()
{

	ICM_ModelMultiCurves* model =(ICM_ModelMultiCurves *) ((ICM_ModelMultiCurves *) GetModel())->Clone() ;
	
	ICM_Option_Tranche* TriggerOption = (ICM_Option_Tranche*) GetSecurity();

	ICM_Mez* MezzRestrikable = (ICM_Mez*) ((ICM_Mez*) TriggerOption ->GetUnderlying())->Clone();
	
	
		// Need to price the CDO
	ICM_Pricer_Distrib_Smile* tmpPricer = new ICM_Pricer_Distrib_Smile() ;
	tmpPricer->Set(MezzRestrikable,model,GetParameters(),GetAsOfDate());
	// Pricing
	tmpPricer->Price(qCMPPRICE);

	ICM_Smile_Correlation* correl=(ICM_Smile_Correlation*) model ->GetCorrelation();


	double yfMatu =( MezzRestrikable->GetEndDateNA() - model->GetAsOfDate()) /365. ;
	itsBeta0 = correl->GetCorrelStrikeDown(yfMatu);
	itsBeta1 = correl->GetCorrelStrikeUp(yfMatu);

	// computing betas with rehauss: 


	
	((ICM_Smile_Correlation*) model ->GetCorrelation())->SetAlready_rescal(false);
	
	ICM_Pricer_Distrib_Smile* tmpPricer2 = new ICM_Pricer_Distrib_Smile() ;
	tmpPricer2->Set(MezzRestrikable,model,GetParameters(),GetAsOfDate());
	
	double subAmount = 0. ;//, MezzAmount = 0. ;

	subAmount = MezzRestrikable->GetSubAmount(MezzRestrikable->GetFeeLeg()->GetStartDate()) + (TriggerOption->GetRehauss() * MezzRestrikable->GetCollateral()->SumNotionals(MezzRestrikable->GetFeeLeg()->GetStartDate()) );
	
	
	MezzRestrikable->SetSubAmount(subAmount);
	tmpPricer2->Price(qCMPPRICE);

	itsBeta2 = correl->GetCorrelStrikeDown(yfMatu);
	itsBeta3 = correl->GetCorrelStrikeUp(yfMatu);

	if (model)
		delete model;
	model = NULL;

	if (MezzRestrikable)
		delete MezzRestrikable;
	MezzRestrikable = NULL;

	if (tmpPricer)
		delete tmpPricer;
	tmpPricer = NULL;


	if (tmpPricer2)
		delete tmpPricer2;
	tmpPricer2 = NULL;



}
void ICM_Pricer_Tranche_Restrikable::ComputeInitialProbas()
{
	// fonction générant les probas de défaut à diffuser
	// Calcul de probas à partir des spreads moyens du portefeuille 
	ARM_Vector* AvgSpreads ;

	ICM_DefaultCurve* CurrentDefcurve ;
	
	ICM_ModelMultiCurves* model = (ICM_ModelMultiCurves *) GetModel() ;

	ICM_Collateral* collat = ((ICM_Mez*) ((ICM_Option_Tranche*) GetSecurity()) ->GetUnderlying())->GetCollateral();
	
	CurrentDefcurve =(ICM_DefaultCurve*) model->GetDefaultCurve(collat->GetIssuersLabels(0));
	
	AvgSpreads = new ARM_Vector(CurrentDefcurve->GetYearTerms()->GetSize()-1,0.);
	double AvgRecovery = 0;
	// Moyenne pondéré par les nominaux
	double Notional = 0;
	double sumNotio = 0;
	for(int i=0;i< collat->GetNbIssuers();i++)
	{
		CurrentDefcurve =(ICM_DefaultCurve*) model->GetDefaultCurve(collat->GetIssuersLabels(i));
		Notional = collat->GetIssuersNotional(i);
		sumNotio +=Notional;
		for ( int j=0;j<CurrentDefcurve->GetRates()->GetSize()-1;j++)
		{
			AvgSpreads->Elt(j) +=CurrentDefcurve->GetRate(j+1) * Notional ;
		}
		AvgRecovery+=(CurrentDefcurve->GetRecovery()* Notional);
	
	}
	for ( int j=0;j<AvgSpreads->GetSize();j++)
		{
			AvgSpreads->Elt(j) = AvgSpreads->Elt(j)/sumNotio;
		}
	AvgRecovery=AvgRecovery/sumNotio;
	
	// create defcurve with avfspreads
	ICM_DefaultCurve* AvgDefcurve = new ICM_Constant_Piecewise(CurrentDefcurve->GetAsOfDate(),
																CurrentDefcurve->GetTerms(),
																AvgSpreads,
																AvgRecovery,
																CurrentDefcurve->GetZeroCurve(),
																K_ADJUSTED,
																K_ADJUSTED,
																CurrentDefcurve ->GetCdsAdj(),
																CurrentDefcurve->GetCurrency(),
																"", //label
																CurrentDefcurve->GetIsSummitCurve(),
																CurrentDefcurve ->GetVolCurve(),
																CurrentDefcurve ->GetSTDCDS_Frequency(),
																CurrentDefcurve ->GetCalibrationAlgo(),
																CurrentDefcurve ->GetCalibrationData(),
																CurrentDefcurve ->GetLagStartDate(),
																CurrentDefcurve->GetParameters()
																);

	
	itsDatesYF =new  ARM_Vector(AvgDefcurve ->GetInterpolYF());
	itsprobas =new ARM_Vector(AvgDefcurve ->GetInterpolSP());

	itsRecov  = AvgRecovery; // TODO case Heterogene: vector of recov

	#ifdef _DEBUG
				FILE* pFile = NULL;
				pFile = fopen("C:\\temp\\AvgDefCurve.txt","w");
				AvgDefcurve ->View("",pFile);
				if ( pFile) fclose(pFile);
	#endif

	if( AvgSpreads)
		delete[] AvgSpreads;
	AvgSpreads = NULL;

	if (AvgDefcurve)
		delete AvgDefcurve;
	AvgDefcurve = NULL;
											

}

double ICM_Pricer_Tranche_Restrikable::FeeLegPV () 
{
	if (GetFeeLegPriceFlg()) return GetFeeLegPrice();

	double Spread = GetInitialSpread();
	double ret =GetPrice(1,Spread);
	SetFeeLegPrice(ret);

	return ret; 
}

double ICM_Pricer_Tranche_Restrikable::GetFeeLegUnit()
{
	if ( itsFeeLegFlg) return itsFeeLegUnit ;
	
	double spread = GetInitialSpread();
	return itsFeeLegUnit ;

	
}
double ICM_Pricer_Tranche_Restrikable::DefLegPV () 
{
	if (GetDefLegPriceFlg()) return GetDefLegPrice();

	double Spread = GetInitialSpread();
	double ret = -GetPrice(2,Spread);
	SetDefLegPrice(ret);

	return ret; 
}

double ICM_Pricer_Tranche_Restrikable::ComputeSpread() 
{
	ICMTHROW(ERR_INVALID_ARGUMENT,"Not implemented"); 
}

void ICM_Pricer_Tranche_Restrikable::ResetPricer(void)
{
	ICM_Pricer::ResetPricer();
	itsInitialSpread = 0.0;
	itsInitialSpreadFlag = false;
	itsFeeLegUnit = 0.0 ;
	itsFeeLegFlg = false ;

}


double ICM_Pricer_Tranche_Restrikable::Accrued()
{
	SetAccrued(0.0) ;
	return 0.0 ; // TODO

	//ICMTHROW(ERR_INVALID_ARGUMENT,"Not implemented"); 
}


double ICM_Pricer_Tranche_Restrikable::ComputeDuration(void)
{
	this->FeeLegPV(); 
	// return GetDuration(); 
	return getValue(qCMPDURATION); 
}

double ICM_Pricer_Tranche_Restrikable::CptDuration(void) 
{
	/*double spread = GetInitialSpread();
	return getValue(qCMPDURATION); */

	double FeeLeg0 ;
	FeeLeg0=0.;


	int initialNbSimul = itsNbSimul;
	double initialSigma = itsSigma;
	


	ICM_Option_Tranche* TriggerOption = (ICM_Option_Tranche*) GetSecurity();
	ICM_Mez* MezzRestrikable = (ICM_Mez*) TriggerOption ->GetUnderlying();

	double TrancheSpread  =MezzRestrikable ->GetCreditSpread(); 
	double InitialTriggerFreq =TriggerOption->GetTriggerFreq() ;
	int initialdiffcdo =  TriggerOption->GetDiffCDO();
	double InitialStrike = TriggerOption->GetStrike();

	double initialBeta2 = itsBeta2;
	double initialBeta3 = itsBeta3;
	
	itsBeta2 = itsBeta0;
	itsBeta3 = itsBeta1;

	double InitialRehauss  = TriggerOption ->GetRehauss();
	
	// DefLeg(t=0) : Temporary Modifs
	itsNbSimul=1 ; 
	TriggerOption->SetTriggerFreq(52.0) ; // Weekly Observation

	double k1 = MezzRestrikable ->GetPercentLow(MezzRestrikable->GetFeeLeg()->GetStartDate());
	TriggerOption->SetDiffCDO(0);
	TriggerOption->SetStrike(k1);
	itsSigma = 0;
	
	Compute() ;
	FeeLeg0 = itsFeeLegUnit * 0.0001 ;

	//FeeLeg0 = GetPrice(1,0.0001);
	

	double duration=0;
	duration = - FeeLeg0*10000. ;
	SetDuration(fabs(duration));
	

	// Remettre paramètres initiaux
	TriggerOption->SetStrike(InitialStrike);
	TriggerOption ->SetTriggerFreq(InitialTriggerFreq);
	TriggerOption->SetDiffCDO(initialdiffcdo);
	MezzRestrikable->SetCreditSpread(TrancheSpread);
	itsNbSimul = initialNbSimul;
	itsSigma = initialSigma;
	itsBeta2 = initialBeta2;
	itsBeta3 = initialBeta3;

	return fabs(duration) ; 
}

double ICM_Pricer_Tranche_Restrikable::ComputeSpread(const double& MtM ) 
{
	if (GetSpreadFlg()) return GetSpread();

	double Spread =Price(qCMPPRICE)*10000./getValue(qCMPDURATION) ; //GetInitialSpread();

	SetSpread(Spread);
	return Spread ;
}

double ICM_Pricer_Tranche_Restrikable::ComputeImpliedVol(const double& Price)
{
	ICMTHROW(ERR_INVALID_ARGUMENT,"Not implemented"); 
}


// 	virtual 
double ICM_Pricer_Tranche_Restrikable::ComputeSensitivity(const qSENSITIVITY_TYPE& typesensi, 
		const std::string& plot,
		const std::string& label,
		double epsilon , double epsilonGamma) 
{
	//return ICM_Pricer::ComputeSensitivity(typesensi,plot,label,epsilon, epsilonGamma); 

	double sensitivity = 0. ;

	ICM_ModelMultiCurves* DefModel = (ICM_ModelMultiCurves*) GetModel();
	//ICM_Option_Tranche* Restrikable_CDO = (ICM_Option_Tranche*) GetSecurity();
	ICM_ModelMultiCurves* ModelDef2 = NULL;
	int nocurve = -1;

	double initialprice=0. ,modifiedprice = 0.;

	initialprice = Price(qCMPPRICE) ;
	ResetPricer();

	switch(typesensi)
	{
		case ICMSPREAD_TYPE : 
		{
		
			ModelDef2 = DefModel->GenerateShiftModel(typesensi,
													plot, 
													label,
													nocurve,
													epsilon);
			SetModel(ModelDef2);
			
			modifiedprice = Price(qCMPPRICE);
			
			if (ModelDef2)
				delete ModelDef2;

			SetModel(DefModel);
			break;
		}
			
		case ICM_GREEK_VEGA_TYPE:
			
			ARM_VolCurve* volcurve = (ARM_VolCurve*)(DefModel->GetVolCurve()->Clone());
			ModelDef2 =(ICM_ModelMultiCurves*) DefModel->Clone() ;

			volcurve ->BumpVolatility(epsilon*100);

			ModelDef2->SetVolCurve(volcurve);

			if (volcurve)
				delete volcurve ;

			SetModel(ModelDef2); // recalcul de initialspread et duration non necessaire
			
			modifiedprice = Price(qCMPPRICE);
			
			if (ModelDef2)
				delete ModelDef2;

			SetModel(DefModel);
			break;

	}
	
	ResetPricer();
	sensitivity = modifiedprice - initialprice ;

	return sensitivity ;

	ICMTHROW(ERR_INVALID_ARGUMENT,"Not implemented"); 

}


double ICM_Pricer_Tranche_Restrikable::GetInitialSpread()
{
	// Initial Spread
	if ( itsInitialSpreadFlag ) 
		return itsInitialSpread;

	
	ComputeInitialProbas();
	//ComputeBetas();
	SetSigmaFromModel();

	// compute duration
	double duration = 0 ;

	duration=CptDuration(); 
	ResetPricer() ;

	SetDuration(fabs(duration));
	// calcul FeeLeg et DefLeg
	
		// Use Initial spread from security
	ICM_Option_Tranche* TriggerOption = (ICM_Option_Tranche*) GetSecurity();

	itsInitialSpread = TriggerOption->GetInitSpread() ;
	Compute() ;

	itsInitialSpreadFlag = true ;
	return itsInitialSpread ;
}

// Nicolas copy

 double ICM_Pricer_Tranche_Restrikable::GetPrice(int typeLeg, double SpreadInit)
{
	int flOK, dlOK;
	double price = 0.0 ;
	switch (typeLeg )
	{
	case 0 : // feeleg - defleg
		flOK = 1;
		dlOK = 1;
		break;
	case 1 : // feeleg
		flOK = 1;
		dlOK = 0;
		break;
	case 2 : // feeleg
		flOK = 0;
		dlOK = 1;
		break;
	}
	price = flOK*SpreadInit*GetFeeLegUnit() -dlOK*DefLegPV();
	return price;
}

 void ICM_Pricer_Tranche_Restrikable::Compute()
{
	 
	double K1, K2;

	ICM_Option_Tranche* TriggerOption = (ICM_Option_Tranche*) GetSecurity();
	ICM_Mez* MezzRestrikable = (ICM_Mez*) TriggerOption ->GetUnderlying();
	ICM_ModelMultiCurves* model = (ICM_ModelMultiCurves *) GetModel() ;

	ARM_ZeroCurve* ircurve = model->GetZeroCurve();

	K1 = MezzRestrikable ->GetPercentLow(MezzRestrikable->GetFeeLeg()->GetStartDate());
	K2 =( MezzRestrikable ->GetPercentHight(MezzRestrikable->GetFeeLeg()->GetStartDate()) +MezzRestrikable ->GetPercentLow(MezzRestrikable->GetFeeLeg()->GetStartDate())) ;
	
	double dealMatu =(MezzRestrikable->GetEndDateNA()- GetAsOfDate())/365.;
	
	int nbNoms = MezzRestrikable->GetCollateral()->GetNbIssuers();
	double triggerFreq = TriggerOption ->GetTriggerFreq();

	double triggerMatu =( TriggerOption ->GetExpiry()- GetAsOfDate())/365.;
	double TriggerStartDate = ( TriggerOption->GetStartDateNA() - GetAsOfDate())/365. ;
	
	double seuilOC = TriggerOption ->GetStrike();
	
	int diffCDO = TriggerOption ->GetDiffCDO();

	double SumSpread = 0.0;
	int nbSpread =0 ;
	double  matuSpread ;
	//double SpreadInit = MezzRestrikable ->GetCreditSpread();
	
	int i, k, h;
	time_t now = time (NULL);
	double seed = int(now);


	double *probaNomt;
	double probaPFt;
	if (itsIsHetero)
	{
		probaNomt = new double[nbNoms];
	}
	else
	{
		probaNomt = NULL;
	}


	double loss, EL, OC;
	double sommeFeeLeg = 0.0, feeLeg, sommeDefLeg = 0.0, defLeg;
	double sommeFeeLegR = 0.0, feeLegR, sommeDefLegR = 0.0, defLegR; 
	double feeLegH1, feeLegH2, feeLegH1R, feeLegH2R; 
	double defLegH1, defLegH2, defLegH1R, defLegH2R; 
	
	//double price, sommePrice = 0.0;
	double spread, dur;
	double K1R = K1+ (TriggerOption->GetRehauss()), K2R = K2+(TriggerOption->GetRehauss());
	double H1, H2, H1R, H2R , H0=0.0;
	int NbNomAvecDef;

	double dY, xt, zt, Dxt, intxt, eintxt, n01;
	
	bool continuer, isResetDate;
	double resetDate;
	
	int triggerOC = 0;
	
	//double spreadTranche;

	double sqrtdt = sqrt(itsTimeStep);
	int kMax = int(dealMatu/itsTimeStep);
	double t, tolt = 10e-9;
	 

	double nbCond = 0.;		
	double probaNTD, probaNTDPrec, lambdaNTD;	
	
	
	int degreeHerm = 12; // valeur possible 40, 20 , 12
	GaussHermite gaussHermite(degreeHerm);
	double *vt = gaussHermite.getx(); // vecteur variable systemique ( taille = degreeHerm)
	double *w  =  gaussHermite.getw(); 
	
	DiffStochIntensityHetero diffusion(itsMRS,  itsSigma, itsRecov, ircurve, nbNoms, itsDatesYF, itsprobas,
										degreeHerm,itsTimeStep,dealMatu,triggerMatu,1.0/triggerFreq,
										itsBeta0,itsBeta1,itsBeta2,itsBeta3,K1,K2,K1R,K2R,8); // degree legendre = 8??? 

	// loss
	int nbDefauts;
	double lgd = 1.0/nbNoms*(1-itsRecov);
	//double *probloss_v_l;

	// Alloc : matrice probloss_v P(L=l | v) 
	double ** probloss_v = new double*[nbNoms+1];
	for (k=0; k <= nbNoms; k++)
	{
		probloss_v[k] = new double[degreeHerm];
	}

	// Pour chaque nom P(def |v)
	double ** probanoms_v;
	double *probapf_v;
	if (itsIsHetero)
	{
		probapf_v = NULL;
		probanoms_v = new double*[nbNoms];
		for (k=0; k < nbNoms; k++)
		{
			probanoms_v[k] = new double[degreeHerm];
		}
	}
	else
	{
		probanoms_v = NULL;
		probapf_v = new double[degreeHerm];
	}


	for (i=0; i < itsNbSimul; i++)
	{
		t = 0.0;
		xt = 0.0;
		intxt = 0.0;
		nbDefauts = 0;
		probaNTDPrec = 0.0;
		//price = 0.0;
		feeLeg=0.0 ;
		defLeg = 0.0;
		feeLegR = 0.0 ;
		defLegR = 0.0 ;

		triggerOC = 0;

		SumSpread= 0.0 ;
		nbSpread = 0 ;

		continuer = true;				
		resetDate = 1/triggerFreq;

		for (k=1; k < kMax && continuer; k++)
		{
			t += itsTimeStep;

			if ( TriggerOption->IsCMSpread())
				matuSpread = t + TriggerOption->GetCMSpreadMatu() ;
			else
				matuSpread = dealMatu;

			n01 = Normal(seed);
			Dxt = - itsMRS * xt * itsTimeStep + itsSigma * sqrtdt * n01;	// discretisation de dX = vol*dW -aXdt
			xt +=  Dxt;		
			//if (xt < 0.0) xt = 0.0;
			intxt += itsTimeStep*Dxt;
			eintxt = exp(-intxt);
			

			// calcul des proba pour chaque valeur de V
			if (itsIsHetero)
			{
				for (h=0; h < nbNoms; h++)
					probaNomt[h]  = diffusion.ProbaSurvInitNom(h, t)*eintxt;

				CptProbaLossHeterog(probloss_v, probanoms_v, probaNomt, vt, degreeHerm, itsBeta0,  nbNoms, nbDefauts+2);
			}
			else
			{
				probaPFt = diffusion.ProbaSurvInitPF(t)*eintxt;
				if (probaPFt >= 1.0)
					probaPFt = 1.0;
				//CptProbaLossHomog(probloss_v, probapf_v, probaPFt, vt, degreeHerm, itsBeta0,  nbNoms, nbDefauts+2);
				diffusion.CptProbaLossHomog(probloss_v, probapf_v, probaPFt, vt, degreeHerm, itsBeta0,  
											nbNoms, nbDefauts+2, t, dealMatu);
			}

			// integration sur V 
			probaNTD = gaussHermite.GaussHermiteInt(probloss_v[nbDefauts+1]);
			
			lambdaNTD = -log(1-probaNTD)/t;
			

			//dY = GenerateDefaut(lambdaNTD, dt, seed);
			dY = GenerateDefaut2(probaNTD-probaNTDPrec , seed);

			if (dY > 0.0) 
				probaNTDPrec = gaussHermite.GaussHermiteInt(probloss_v[nbDefauts+2]);
			else
				probaNTDPrec = probaNTD;
			
			nbDefauts += dY;
			loss = nbDefauts*lgd;
			
			if (  t + tolt  >= resetDate - itsAvgSpreadFreq)
			{
				SumSpread+= diffusion.SpreadPF(t, matuSpread, xt, dur);
				nbSpread++;

			}
			// gestion resetdate			
			if (  t + tolt  >= resetDate)
			{
				isResetDate = true;
				resetDate += 1/triggerFreq;
				spread =SumSpread/nbSpread ; //diffusion.SpreadPF(t, dealMatu, xt, dur);
				nbSpread=0;
				SumSpread = 0.0 ;
			}
			else
			{
				isResetDate = false;
			}

			if ((t - tolt <= triggerMatu) && (TriggerStartDate < t-tolt))
			{													
				// on ne calcule le trigger que si Loss < Kd
				if (( triggerOC == 0) && (loss < K1 || K1 ==0) && (t < triggerMatu))
				{
					if (isResetDate)
					{
						
						EL =itsELCoeff* spread * ( itsUsePeriodIsEL<1 ? 1 : (dealMatu - t));
						OC = K1 - (itsLossCoeff *loss) - EL;						
						triggerOC = ( OC < seuilOC ? 1 : 0);
					}
				}
			}
			
			// Calcul du prix tranche
			if (triggerOC)
			//if ((t < triggerMatu+tolt) && (t >= triggerMatu-tolt))
			{
				H1 = K1-loss;
				H2 = K2 - loss;
				H1R = K1R-loss;
				H2R = K2R - loss;
				NbNomAvecDef = nbNoms-nbDefauts;


				zt =  xt;
				PriceTranche(feeLegH1,defLegH1,H0,H1,t,zt,itsBeta0,diffusion,gaussHermite,NbNomAvecDef,lgd,probloss_v,0);
				
				PriceTranche(feeLegH2,defLegH2,H0,H2,t,zt,itsBeta1,diffusion,gaussHermite,NbNomAvecDef,lgd,probloss_v,1);
							
				defLeg = (H2*defLegH2-H1*defLegH1)/(H2-H1);
				feeLeg = (H2*feeLegH2-H1*feeLegH1)/(H2-H1);
					
				//PriceTranche(feeLeg,defLeg,H1,H2,t,zt,itsBeta1,diffusion,gaussHermite,NbNomAvecDef,lgd,probloss_v);

				if (diffCDO > 0)
				{
				
					PriceTranche(feeLegH1R,defLegH1R,H0,H1R,t,zt,itsBeta2,diffusion,gaussHermite,NbNomAvecDef,lgd,probloss_v,2);
					PriceTranche(feeLegH2R,defLegH2R,H0,H2R,t,zt,itsBeta3,diffusion,gaussHermite,NbNomAvecDef,lgd,probloss_v,3);
				
					
					defLegR = (H2R*defLegH2R-H1R*defLegH1R)/(H2R-H1R);
					feeLegR = (H2R*feeLegH2R-H1R*feeLegH1R)/(H2R-H1R);

					//PriceTranche(feeLegR,defLegR,H1R,H2R,t,zt,itsBeta2,diffusion,gaussHermite,NbNomAvecDef,lgd,probloss_v);
				}
				else
				{
					feeLegR = defLegR = 0.0;
				}
				
				
				//price = flOK*SpreadInit*(feeLegR-feeLeg) -dlOK*(defLegR- defLeg);
				continuer = false;
			}

			if (t > triggerMatu)
				continuer = false;

		} // fin for k
		
		sommeFeeLeg += (feeLegR-feeLeg) ;
		sommeDefLeg += (defLegR- defLeg) ;
		//sommePrice += price;

		//if (price < 10e-6)
		//{
		//	price = 0.0;
		//}

		#ifdef _DEBUG
			//fprintf(stream_hw_tree,"\n");
		#endif
		

	} // fin for i


	
	// Delete
	if (probaNomt) 
		delete probaNomt;
	probaNomt = NULL;

	
	for (k=0; k <= nbNoms; k++)
	{
		if (probloss_v[k]) delete probloss_v[k];
	}
	if (probloss_v) 
		delete probloss_v;
	probloss_v = NULL;


	for (k=0; k < nbNoms; k++)
	{
		if (probanoms_v)
			if (probanoms_v[k])	delete probanoms_v[k];
	}
	if (probanoms_v) 
		delete probanoms_v;
	probanoms_v = NULL;
	
	if (probapf_v)
		delete[] probapf_v;
	probapf_v = NULL;

	//return sommePrice/i;
	itsFeeLegUnit = sommeFeeLeg / i ;
	itsFeeLegFlg = true ;

	SetDefLegPrice(sommeDefLeg/i) ;



}


void ICM_Pricer_Tranche_Restrikable::PriceTranche(double &fl, double &dl,double &K1, double &K2, double &t ,
												  double &xt ,double &beta, DiffStochIntensityHetero &diffusion , 
												  GaussHermite &gaussHermite ,int &nbNoms, double &lgd, 
												  double** probloss_v , int numK, int adj)
{
	double probaPFti;
	double un_lossTR; // 1 - loss sur la tranche
	double ELTR_ti_1, ELTR_ti; // 1 - Expextedloss sur la tranche en ti et ti-1

	double T ;
	// getting percentLow and percentHigh of the tranche
	

	ICM_Option_Tranche* TriggerOption = (ICM_Option_Tranche*) GetSecurity();
	ICM_Mez* MezzRestrikable = (ICM_Mez*) TriggerOption ->GetUnderlying();

	ICM_ModelMultiCurves* model = (ICM_ModelMultiCurves *) GetModel() ;

	ARM_ZeroCurve* ircurve = model->GetZeroCurve();

	
	T = (MezzRestrikable->GetEndDateNA()- GetAsOfDate())/365.;

	double lgdTR = lgd/(K2-K1); // impact d'un defaut sur la tranche
	double *vt = gaussHermite.getx();
	int degreeHerm = gaussHermite.getDeg();

	double * probapf_v = new double[degreeHerm];

	int lMax = floor(K2/lgd);
	int lMin = floor(K1/lgd);
	double *probloss = new double[lMax+1];
	int k; //


	double ti = T;
	double df_ti;

	fl = 0.0;
	dl = 0.0;
	while (ti > t)
	{
		probaPFti = diffusion.ProbaSurvPF(t, ti, xt);

		//CptProbaLossHomog(probloss_v, probapf_v, probaPFti, vt, degreeHerm, beta,  nbNoms, lMax);
		diffusion.CptProbaLossHomog(probloss_v,probapf_v,probaPFti,vt,degreeHerm,beta,nbNoms,lMax,t,ti,adj);

		// integration sur V 
		ELTR_ti = 0.0;
		un_lossTR;
		for (k= 0; k <= lMax; k++)
		{
			un_lossTR = 1.0 - (DMAX(k*lgd-K1, 0) - DMAX(k*lgd-K2, 0))/(K2-K1); 
			probloss[k] = gaussHermite.GaussHermiteInt(probloss_v[k]);

			ELTR_ti += un_lossTR*probloss[k];			
		}
		
		
		ELTR_ti = 1.0 - ELTR_ti;

		
		//df_ti = exp(-ti*itsIR);
		df_ti = ircurve->DiscountPrice(ti);

		fl += df_ti * (1.0 - ELTR_ti)*0.25;
		if ( ti < T) 
			dl += df_ti*(ELTR_ti_1 - ELTR_ti);

		ELTR_ti_1 = ELTR_ti;
		ti -= 0.25;
	}

	df_ti = ircurve->DiscountPrice(ti);
	dl += df_ti * ELTR_ti;
	// 

	fl *= ircurve->DiscountPrice(t);
	dl *= ircurve->DiscountPrice(t);

	int jt = t/diffusion.getDtReset()-1;
	dl  -=diffusion.GetCorrET(numK,jt);
}

void ICM_Pricer_Tranche_Restrikable::View(char* id, FILE* ficOut) 
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

	fprintf(fOut, "-----------------------------------------------------------\n");
    fprintf(fOut, "---                      Restrikeable MC Pricer         ---\n");
	fprintf(fOut, "-----------------------------------------------------------\n");
	
	fprintf(fOut, "\n");
	fprintf(fOut, " Beta 0-Kdown :  %f\n",itsBeta0*100);
	fprintf(fOut, " Beta 0-Kup : %f\n",itsBeta1*100);
	fprintf(fOut, " Beta 0-Kdown+Rehauss :  %f\n",itsBeta2*100);
	fprintf(fOut, " Beta 0-Kup+Rehauss :  %f\n",itsBeta3*100);

	
	fprintf(fOut, "\n");
	fprintf(fOut, " Is Heterogene :  %i\n",itsIsHetero);
	fprintf(fOut,"\n");

	fprintf(fOut, "    ----------------------- Spread diffusion --------------------            \n");
	fprintf(fOut, "\n");
	fprintf(fOut, " Mean Reversion :  %f\n",itsMRS);
	fprintf(fOut, " Volatility :  %f\n",itsSigma);
	fprintf(fOut, "\n");

	int size = 0;
	size = itsDatesYF->GetSize();
	fprintf(fOut, "     -------------   Initial Probas   -------------         \n");
	fprintf(fOut, "\n");
	for ( int i=0; i<size;i++)
	{
		fprintf(fOut, "\t%f\t",itsDatesYF->Elt(i));
		fprintf(fOut, "%f\t", itsprobas ->Elt(i));
		fprintf(fOut, "\n");

	}
	fprintf(fOut, "\n");

	ICM_Pricer::View("",fOut);
	
	if ( ficOut == NULL )
    {
       fclose(fOut);
    }

}

/* __________________________________________________________________________________________________________*/
 

double Uniform(double &seed)
{
   
    double A = seed;
    double test = A * 16807 / 2147483647;
	int enttest = int(test);
	double reste;

    if (test > 0 )
        reste = (A * 16807) - enttest * 2147483647.;
    else
        reste = (A * 16807) - (enttest + 1) * 2147483647.;
    
    
    double u = reste / 2147483647.;
    seed = reste;
    
    return u;
}

double Normal(double &seed)
{
	return NormInv(Uniform(seed));
}

double GenerateDefaut(double lambda, double dt, double &seed)
{
	// Exponentielle
	double dY = 0.0;
	double prob = 1. -exp(-lambda*dt);

	if (Uniform(seed) < prob)
		dY = 1.0;
	
	return dY;
}

double GenerateDefaut2(double prodadef, double &seed)
{
	double dY = 0.0;

	if (Uniform(seed) < prodadef)
		dY = 1.;
	
	return dY;
}

void GenerateVar(double lambda, double dt, double &seed, double &W1, double& W2, double& dY)
{
    double pi = 3.1415927;
    int h = 0;
    
    
    double U1 = Uniform(seed);
    double U2 = Uniform(seed);

    
    W1 = sqrt(-2 * log(U1)) * cos(2 * pi * U2);         //Normale centrée réduite
    //W2 = sqrt(-2 * log(U1)) * sin(2 * pi * U2);         //Normale centrée réduite
   
	

	// Exponentielle
	
	if (Uniform(seed) < 1 -exp(-125*lambda*dt))
		dY = 1;
	else
		dY = 0;

	// Poisson
    /*
	double  X = 1.0;
    while (X >= exp(-pdef * dt))
	{
		X = X * Uniform(seed);
        h = h + 1;
    }

    
    dY = h - 1;
/* 
	// dY = F-1(U2)
	double F_1 = 0, proba_i = 0.0;
	int i = 0, k, cni, n = 125;
	int ifact;
	double p = 0.1; //lambda*dt;
	
	while (U2 > F_1)
	{
		// calcul CNI
		ifact = 1;
		for (k = 1; k <= i; k++)
			ifact *=k;

		cni = 1;
		for (k = n-i+1; k <= n ; k++)
			cni *=k;

		cni /= ifact;
		proba_i = cni*pow(1-p,n-i)*pow(p,i);
		F_1 += proba_i;
		i++;
	}
	dY = i-1;
*/
}


// prob : proba de  defaut en t
// calcule la probabilte de loss sur un portefeuille pour chacune des valeurs du vect vt (variable systemique)
// input prob_noms
//		 vsyst
// output  probaloss_v : proba loss sachant v (nbnom+1 x Degree)
//		   proba_noms_v : proba def nom sachant v (nbnom x Degree)
// util 

void CptProbaLossHeterog(double ** &probaloss_v, double ** &proba_noms_v, double *&prob_noms, double* &vsyst, 
				  int &degreeHerm, double &beta, int &nbNoms, int nbDefmax)
{
	int k, j, h;
	double phi_1, beta2;
	beta2 = beta*beta;

	// calcul des proba pour chaque valeur de v	

	for (k=0; k < nbNoms; k++)
	{
		phi_1 = NormInv(prob_noms[k]);

		for (h = 0; h < degreeHerm; h++)
		{
			proba_noms_v[k][h] = probaSachV( phi_1,  vsyst[h], beta, beta2);
		}
	}
	
	

	// init proba si aucun nom
	for (h = 0; h < degreeHerm; h++)
		probaloss_v[0][h] = 1.0;

	for (k=1; k <= nbDefmax; k++)
	{
		for (h = 0; h < degreeHerm; h++)
			probaloss_v[k][h] = 0.0;
	}

	// Boucle jusqu'au nbre de noms
	for (j=0; j < nbNoms; j++)
	{
		for (k=nbDefmax; k >0 ; k--)
		{
			for (h = 0; h < degreeHerm; h++)
			{
				probaloss_v[k][h] = probaloss_v[k][h]*proba_noms_v[j][h] + probaloss_v[k-1][h]*(1-proba_noms_v[j][h]);
			}
		}
	
		for (h = 0; h < degreeHerm; h++)
		{
			probaloss_v[0][h] = probaloss_v[0][h]*proba_noms_v[j][h];
		}
	}

}


void CptProbaLossHomog(double ** &probaloss_v, double * &proba_pf_v, double &prob_noms, double* &vsyst, 
				  int &degreeHerm, double &beta, int &nbNoms, int nbDefmax)
{
	int k, j, h;
	double phi_1, beta2;
	beta2 = beta*beta;

	// calcul des proba pour chaque valeur de v		
	phi_1 = NormInv(1.0-prob_noms);

	for (h = 0; h < degreeHerm; h++)
	{
			proba_pf_v[h] = 1.0- probaSachV( phi_1,  vsyst[h], beta, beta2);
	}
		

	// init proba si aucun nom
	for (h = 0; h < degreeHerm; h++)
		probaloss_v[0][h] = 1.0;

	for (k=1; k <= nbDefmax; k++)
	{
		for (h = 0; h < degreeHerm; h++)
			probaloss_v[k][h] = 0.0;
	}

	// Boucle jusqu'au nbre de noms
	for (j=0; j < nbNoms; j++)
	{
		for (k=nbDefmax; k >0 ; k--)
		{
			for (h = 0; h < degreeHerm; h++)
			{
				probaloss_v[k][h] = probaloss_v[k][h]*proba_pf_v[h] + probaloss_v[k-1][h]*(1-proba_pf_v[h]);
			}
		}
	
		for (h = 0; h < degreeHerm; h++)
		{
			probaloss_v[0][h] = probaloss_v[0][h]*proba_pf_v[h];
		}
	}

}
/*_________________________________________________________________________________________________________________________*/

double DiffStochIntensity::ProbaSurvInitPF(double &t)
{
	// t doit etre en year term

	if (t <= itsdateProba->Elt(0))
		return itsprobaInitPF->Elt(0);

	if (t >= itsdateProba->Elt(itssizeProbaInit-1))
		return itsprobaInitPF->Elt(itssizeProbaInit-1);

	int k = 0;

	while (itsdateProba->Elt(k) < t)
		k++;

	double prob = itsprobaInitPF->Elt(k-1) + (itsprobaInitPF->Elt(k)-itsprobaInitPF->Elt(k-1)) *(t-itsdateProba->Elt(k-1))/(itsdateProba->Elt(k)-itsdateProba->Elt(k-1));
	return prob;
}


// probabilité qu un defaut survienne apres T sachant pas de défaut en t
// P(to> T)
double DiffStochIntensity::ProbaSurvPF(double &t, double &T, double &xt)
{
	double Beta;
    double Phi0;
    double Phi1;
    double Alpha;
    
    Beta = (1 - exp(- itsMRS * (T - t))) / itsMRS;
    Phi0 = itssigma * itssigma * (1 - exp(-2 * itsMRS * t)) / (2 * itsMRS);
    Phi1 = (1 - exp(-itsMRS * t)) - (1 - exp(-2 * itsMRS * t)) / 2;
    Phi1 = Phi1 * itssigma * itssigma / (itsMRS * itsMRS);
    
    Alpha = -Beta * (Beta * Phi0 / 2 + Phi1);
    
    double proba = ProbaSurvInitPF(T)/ProbaSurvInitPF(t) * exp(Alpha - Beta * xt);
    
    return proba;	
}

double DiffStochIntensity::MgE(double &t, double &T, double &xt)
{
	double Beta;
    double Phi0;
    double Phi1;
    double Alpha;
    
    Beta = (1 - exp(- itsMRS * (T - t))) / itsMRS;
    Phi0 = itssigma * itssigma * (1 - exp(-2 * itsMRS * t)) / (2 * itsMRS);
    Phi1 = (1 - exp(-itsMRS * t)) - (1 - exp(-2 * itsMRS * t)) / 2;
    Phi1 = Phi1 * itssigma * itssigma / (itsMRS * itsMRS);
    
    Alpha = -Beta * (Beta * Phi0 / 2 + Phi1);
    
    double proba =  exp(Alpha - Beta * xt);
    
    return proba;	
}

double DiffStochIntensity::SpreadPF(double &t, double &T, double &xt, double &duration)
{
	double ti_1, ti = T;
	double feeLeg = 0.0;
	double defLeg = 0.0;
	double theta = itsbaseCDS; 

	double Pi, Pi_1 = ProbaSurvPF(t, T, xt);
	
	while (ti > t)
	{
        Pi = Pi_1;
        ti_1 = ti - itsbaseCDS;

        if (ti_1 <= t) 
		{
            Pi_1 = 1;
            theta = ti - t;
		}
        else
		{
            Pi_1 = ProbaSurvPF(t, ti_1, xt);
            theta = itsbaseCDS;
		}
        
        //feeLeg += theta * exp(-itsIR * (ti - t)) * Pi;
        //defLeg += (1 - itsrecovery) * (Pi_1 - Pi) * exp(-itsIR * (ti - t));
		feeLeg+= theta *Pi *itsZeroCurve->DiscountPrice((ti-t)) ;
		defLeg += (1 - itsrecovery) * (Pi_1 - Pi) * itsZeroCurve->DiscountPrice((ti-t)) ;
        ti = ti_1;
    }
    
    double spreadCDS = defLeg / feeLeg;

	duration = feeLeg;

	return spreadCDS;
}



double DiffStochIntensity::lambdaY(double &t)
{
	double t_eps = t + EPS;
	double DlogP = (ProbaSurvInitPF(t)-ProbaSurvInitPF(t_eps))/EPS;

	return DlogP;
}


void DiffStochIntensity::CptProbaLossHomog(double ** &probaloss_v, double * &proba_pf_v, double &prob_noms, double* &vsyst, 
				  int &degreeHerm, double &beta, int &nbNoms, int nbDefmax, double t, double T, int adj)
{
	
	double phi_1, beta2;
	beta2 = beta*beta;

	// calcul des proba pour chaque valeur de v		
	phi_1 = NormInv(1.0-prob_noms);
	
	int h;
	int j = t/itsdt;
	int k = T/itsbaseCDS;


	for (h = 0; h < degreeHerm; h++)
	{
		proba_pf_v[h] = 1.0- probaSachV( phi_1,  vsyst[h], beta, beta2);

		// Ajustement _corrEPV = E[p|V] - p0|V
		if (adj == 1)
			proba_pf_v[h] -= itsCorrEPV[h][j][k];
		
	}

	


	// init proba si aucun nom
	for (h = 0; h < degreeHerm; h++)
		probaloss_v[0][h] = 1.0;

	for (k=1; k <= nbDefmax; k++)
	{
		for (h = 0; h < degreeHerm; h++)
			probaloss_v[k][h] = 0.0;
	}

	// Boucle jusqu'au nombre de noms
	for (j=0; j < nbNoms; j++)
	{
		for (k=nbDefmax; k >0 ; k--)
		{
			for (h = 0; h < degreeHerm; h++)
			{
				probaloss_v[k][h] = probaloss_v[k][h]*proba_pf_v[h] + probaloss_v[k-1][h]*(1-proba_pf_v[h]);
			}
		}
	
		for (h = 0; h < degreeHerm; h++)
		{
			probaloss_v[0][h] = probaloss_v[0][h]*proba_pf_v[h];
		}
	}

}

void DiffStochIntensity::cptCorrEpv(double sigma)
{
	
	int i, j, k, h;
	double beta2 = itsBetas[0]*itsBetas[0], p0, sigmax, ds;

	
	// borne d'integration correpondant à x
	double a , b;
	int imax = itsDegreeHerm;
	int jmax = itsDiffmatu/itsdt + 1;
	int kmax = itsDealMatu/itsbaseCDS + 1;

	itsCorrEPV = new double **[imax];
	for (i=0; i < imax; i++)
	{
		itsCorrEPV[i] = new double*[jmax];

		for (j=0; j < jmax; j++)
		{
			itsCorrEPV[i][j] = new double[kmax];
			for (k=0; k < kmax; k++)
			{
				itsCorrEPV[i][j][k] = 0.0;
			}
		}
	}


	double *fx = new double[itsDegreeLegendre]; // valeur des proba = f(x) 

	double *v = itsGaussHermite.getx(); // vecteur variable syst ( taille = degreeHerm)
	double *vx = itsGaussLegendre.getx(); // vecteur  pour integration ( taille = degreeHerm)
	double *w  =  itsGaussLegendre.getw(); // poids

	double x, t , proba, phi_1;
	double T = itsbaseCDS;
	double V;

	if (sigma == 0.0)
		return;

	
	for (i=0; i < imax; i++) // boucle sur variable systemique
	{
		V =  v[i];
		t = itsdt;
		for (j=0; j < jmax; j++) // boucle sur t
		{
			sigmax = sigma*sqrt((1.0-exp(-2.0*0.025*t))/(2*0.025));
			a = -6.0*sigmax;
			b= -a;

			T = itsDealMatu;
			for (k=kmax-1; (k >=0)&& (t < T); k--) // boucle sur T
			{
				
				// calcul de E[p|V] grace gaussLegendre
				for (h = 0; h < itsDegreeLegendre ; h++)
				{
					x = a+ (b-a)* vx[h];
					proba = 1.0 - ProbaSurvPF(t, T, x);
					phi_1 = NormInv(proba);

					//if (sigma > 0.0)
						ds = exp(-x*x/(2.0*sigmax*sigmax));
					//else
					//	ds = 1.0;
					fx[h] = (1.0 - probaSachV( phi_1,  V, itsBetas[0], beta2))*ds;
					//fx[h] = exp(-x*x/(2.0*sigmax*sigmax));
					
				}

				
				p0= 1.0-ProbaSurvInitPF(T)/ProbaSurvInitPF(t);
				phi_1 = NormInv(p0);
				p0 = 1.0-probaSachV( phi_1,  V, itsBetas[0], beta2);

				itsCorrEPV[i][j][k] = itsGaussLegendre.GaussIntegration(fx, a, b);
				itsCorrEPV[i][j][k] /= SQRT2PI*sigmax;
				itsCorrEPV[i][j][k] -= p0;

				T -= itsbaseCDS;

			}
			t += itsdt;
		}
	}

}
// **************************************************************************

void DiffStochIntensityHetero::cptCorrET(void)
{
	int i, j, h;
	double t = 0.0;

	int jmax = itsDiffmatu/itsdtReset;
	//int jmax = _diffMatu/_dt + 1;

	double *vx = itsGaussLegendre.getx();
	double *fx = new double[itsDegreeLegendre]; // valeur des proba = f(x)
	double sigmax, a, b, x, ds;

	double feeLeg, defLeg, prix0;
	double K0 = 0.0;
	
	for (i=0; i < 4; i++)
	{
		itsCorrET[i] = new double[jmax];
		for (j=0; j < jmax; j++) // boucle sur t
		{
			itsCorrET[i][j] = 0.0;
		}
	}

	if( itssigma == 0.0)
		return;

	for (i=0; i < 4; i++)
	{
		itsCorrET[i] = new double[jmax];
		t = 0.0;

		for (j=0; j < jmax; j++) // boucle sur t
		{
			t += itsdtReset;
			sigmax = itssigma*sqrt((1.0-exp(-2.0*itsMRS*t))/(2*itsMRS));
			a = -6.0*sigmax;
			b = -a;

			// calcul de E[TR|V] grace gaussLegendre
			for (h = 0; h < itsDegreeLegendre ; h++)
			{
				x = a+ (b-a)* vx[h];
				
				PrixTranche(feeLeg, defLeg, K0, itsStrikes[i], t, itsDealMatu, x, itsBetas[i], 0);

				//if (sigma > 0.0)
					ds = exp(-x*x/(2.0*sigmax*sigmax));
				//else
				//	ds = 1.0;

				//fx[h] = (feeLeg*0.0001-defLeg)*ds;
				fx[h] = (defLeg)*ds;
										
			} // fin for h

			itsCorrET[i][j] = itsGaussLegendre.GaussIntegration(fx, a, b);
			itsCorrET[i][j] /= SQRT2PI*sigmax;

			PrixTranche0(feeLeg, defLeg, K0, itsStrikes[i], t, itsDealMatu, itsBetas[i], 0);
			prix0 = defLeg;
			itsCorrET[i][j] -= prix0;

		} // fin for j
	} // fin for i
}

double DiffStochIntensityHetero::ProbaSurvInitNom(int &nom, double &t)
{
	// t doit etre en year term
	ARM_Vector* probaInit =new ARM_Vector(itsprobaInitNom,0,nom);

	if (t < itsdateProba->Elt(0))
		return probaInit->Elt(0);

	if (t >= itsdateProba->Elt(itssizeProbaInit-1))
		return probaInit->Elt(itssizeProbaInit-1);

	int k = 0;

	while (itsdateProba->Elt(k) < t)
		k++;

	double prob = probaInit->Elt(k-1) + (probaInit->Elt(k)- probaInit->Elt(k-1)) *(t-itsdateProba->Elt(k-1))/(itsdateProba->Elt(k)-itsdateProba->Elt(k-1));
	return prob;
}






// probabilité qu'un defaut survienne apres T sachant pas de défaut en t
double DiffStochIntensityHetero::ProbaSurvNom(int &nom, double &t, double &T, double &xt)
{
	double Beta;
    double Phi0;
    double Phi1;
    double Alpha;
    
    Beta = (1 - exp(- itsMRS * (T - t))) / itsMRS;
    Phi0 = itssigma * itssigma * (1 - exp(-2 * itsMRS * t)) / (2 * itsMRS);
    Phi1 = (1 - exp(-itsMRS * t)) - (1 - exp(-2 * itsMRS * t)) / 2;
    Phi1 = Phi1 * itssigma * itssigma / (itsMRS * itsMRS);
    
    Alpha = -Beta * (Beta * Phi0 / 2 + Phi1);
    
    double proba = ProbaSurvInitNom(nom, T)/ProbaSurvInitNom(nom, t) * exp(Alpha - Beta * xt);
    
    return proba;	
}


void DiffStochIntensityHetero::PrixTranche(double &fl, double &dl, double &K1, double &K2, double &t, double &T, 
				 double &xt, double &beta, int adj)
{
	
	double probaPFti;
	double un_lossTR; // 1 - loss sur la tranche
	double ELTR_ti_1, ELTR_ti; // 1 - Expextedloss sur la tranche en ti et ti-1
	double lgdTR = itsLgd/(K2-K1); // impact d'un defaut sur la tranche
	double *vt = itsGaussHermite.getx();

	double * probapf_v = new double[itsDegreeHerm];

	int lMax = floor(K2/itsLgd);
	int lMin = floor(K1/itsLgd);
	double *probloss = new double[lMax+1];
	int k; //


	double ti = T;
	double df_ti;
	fl = 0.0;
	dl = 0.0;
	while (ti > t)
	{
		probaPFti = ProbaSurvPF(t, ti, xt);

		CptProbaLossHomog(itsprobloss_v, probapf_v, probaPFti, vt, itsDegreeHerm, beta,  itsnbNoms, lMax,
								t, ti, adj);

		// integration sur V 
		ELTR_ti = 0.0;
		un_lossTR;
		for (k= 0; k <= lMax; k++)
		{
			un_lossTR = 1.0 - (DMAX(k*itsLgd-K1, 0) - DMAX(k*itsLgd-K2, 0))/(K2-K1); 
			probloss[k] = itsGaussHermite.GaussHermiteInt(itsprobloss_v[k]);

			ELTR_ti += un_lossTR*probloss[k];			
		}
		
		
		ELTR_ti = 1.0 - ELTR_ti;


		df_ti =itsZeroCurve ->DiscountPrice(ti);// exp(-ti * itsIR);

		fl += df_ti * (1.0 - ELTR_ti)*0.25;
		if ( ti < T) 
			dl += df_ti*(ELTR_ti_1 - ELTR_ti);

		ELTR_ti_1 = ELTR_ti;
		ti -= 0.25;
	}

	fl *=itsZeroCurve->DiscountPrice(t);// exp(-t * itsIR);
	dl *= itsZeroCurve->DiscountPrice(t);// exp(-t * itsIR);

	delete probloss;
	delete probapf_v;
	
}



// fl feeleg
// dl defleg
void DiffStochIntensityHetero::PrixTranche0(double &fl, double &dl, double &K1, double &K2, double &t, double &T, 
				 double &beta, int adj)
{
	
	double probaPFti;
	double un_lossTR; // 1 - loss sur la tranche
	double ELTR_ti_1, ELTR_ti; // 1 - Expextedloss sur la tranche en ti et ti-1
	double lgdTR = itsLgd/(K2-K1); // impact d'un defaut sur la tranche
	double *vt = itsGaussHermite.getx();
	

	double * probapf_v = new double[itsDegreeHerm];

	int lMax = floor(K2/itsLgd);
	int lMin = floor(K1/itsLgd);
	double *probloss = new double[lMax+1];
	int k; //


	double ti = T;
	double df_ti;
	fl = 0.0;
	dl = 0.0;
	while (ti > t)
	{
		probaPFti = ProbaSurvInitPF(ti)/ProbaSurvInitPF(t);

		CptProbaLossHomog(itsprobloss_v, probapf_v, probaPFti, vt, itsDegreeHerm, beta,  itsnbNoms, lMax,
								t, ti, adj);

		// integration sur V 
		ELTR_ti = 0.0;
		un_lossTR;
		for (k= 0; k <= lMax; k++)
		{
			un_lossTR = 1.0 - (DMAX(k*itsLgd-K1, 0) - DMAX(k*itsLgd-K2, 0))/(K2-K1); 
			probloss[k] = itsGaussHermite.GaussHermiteInt(itsprobloss_v[k]);

			ELTR_ti += un_lossTR*probloss[k];			
		}
		
		
		ELTR_ti = 1.0 - ELTR_ti;


		df_ti =itsZeroCurve->DiscountPrice(ti);//  exp(-ti * itsIR);

		fl += df_ti * (1.0 - ELTR_ti)*0.25;
		if ( ti < T) 
			dl += df_ti*(ELTR_ti_1 - ELTR_ti);

		ELTR_ti_1 = ELTR_ti;
		ti -= 0.25;
	}

	fl *= itsZeroCurve->DiscountPrice(t);// exp(-t * itsIR);
	dl *= itsZeroCurve->DiscountPrice(t);// exp(-t * itsIR);

	delete probloss;
	delete probapf_v;
	
}
