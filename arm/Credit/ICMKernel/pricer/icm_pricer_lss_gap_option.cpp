#include "firsttoinc.h"
#include "ICMKernel\pricer\icm_pricer_lss_gap_option.h"
#include "ICMKernel\inst\icm_lss_gap_option.h"

#include "ICMKernel\pricer\icm_pricer_adviser.h"
#include "ICMKernel\glob\icm_maths.h"
#include "ICMKernel\glob\icm_smile_correlation.h"
#include "ICMKernel\mod\modelmulticurves.h"
#include "ICMKernel\crv\icm_defaultcurveSimple.h"
#include "ICMKernel\crv\icm_distriblosscalculator.h"
#include "ICMKernel\util\icm_utils.h"
#include "ICMKernel\inst\icm_collateral.h"
//-----------------------------------------------------------------------------
// Pricing du Cdo sous-jacent
//-----------------------------------------------------------------------------
void ICM_Pricer_LssGapOption::PriceUnderlying()
{
	if (its_SS_Pricer) delete its_SS_Pricer;

	ICM_Pricer_Advisor Advisor;

	ICM_Mez* underlying_cdo = ((ICM_Lss_Gap_Option*) GetSecurity())->GetUnderlyingCdo();

	ICM_ModelMultiCurves* mmc = (ICM_ModelMultiCurves*)GetModel();
	ICM_Correlation* correl = mmc->GetCorrelation();

	correl->ResetBetas();


	its_SS_Pricer = Advisor.GeneratePricer(underlying_cdo,mmc,
										   ICM_PRICER_HOMOGENEOUS_SMILE,
										   CREDIT_DEFAULT_VALUE,
										   &GetParameters(),GetAsOfDate());

	its_SS_Pricer->Price(qCMPPRICE);

}


//-----------------------------------------------------------------------------
// Pricing du Cdo sous-jacent
//-----------------------------------------------------------------------------
void ICM_Pricer_LssGapOption::LoadingParameters(long& seed,long& nbdef,
	double& initial_price,double& Matu_Spread,double& loss,double& Simply,double& Time2Matu,
	double& Leverage,double& Fast,double& volspreads, double& volcorrel,double& VarReduction)
{

	// if (GetParameters() != NULL)
	if (!GetParameters().empty())
	{
		//seed for monte-carlo
		if (GetParameters().GetColVect("SEED"))
		{seed = (long) GetParameters().GetColVect("SEED")->Elt(0);}

		//exact number of defaults (random choose)
		if (GetParameters().GetColVect("NBDEF"))
		{nbdef = (long) GetParameters().GetColVect("NBDEF")->Elt(0);}

		//initial pv for the CDO
		if (GetParameters().GetColVect("INITIALPV"))
		{initial_price = GetParameters().GetColVect("INITIALPV")->Elt(0);}

		//forced maturity for initial spreads
		if (GetParameters().GetColVect("MATU_SPREAD"))
		{Matu_Spread = GetParameters().GetColVect("MATU_SPREAD")->Elt(0);}

		//imposed loss in %
		if (GetParameters().GetColVect("LOSS"))
		{loss = GetParameters().GetColVect("LOSS")->Elt(0);}

		//simplify case for curve building
		if (GetParameters().GetColVect("SIMPLIFY"))
		{Simply = GetParameters().GetColVect("SIMPLIFY")->Elt(0);}

		//time to maturity
		if (GetParameters().GetColVect("TIME2MATU"))
		{Time2Matu = GetParameters().GetColVect("TIME2MATU")->Elt(0);}

		//Levrage for Super Senior
		if (GetParameters().GetColVect("LEVERAGE"))
		{Leverage = GetParameters().GetColVect("LEVERAGE")->Elt(0);}

		//Levrage for Super Senior
		if (GetParameters().GetColVect("FAST"))
		{Fast = GetParameters().GetColVect("FAST")->Elt(0);}

		//volatility spreads
		if (GetParameters().GetColVect("VOL_SPREADS"))
		{volspreads = GetParameters().GetColVect("VOL_SPREADS")->Elt(0);}

		//volatility correl
		if (GetParameters().GetColVect("VOL_CORREL"))
		{volcorrel = GetParameters().GetColVect("VOL_CORREL")->Elt(0);}

		//volatility correl
		if (GetParameters().GetColVect("VAR_REDUCTION"))
		{VarReduction = GetParameters().GetColVect("VAR_REDUCTION")->Elt(0);}
	}

}

//-----------------------------------------------------------------------------
// Pricing de la Gap Option en Monte-Carlo
//-----------------------------------------------------------------------------
double ICM_Pricer_LssGapOption::ComputePrice(qCMPMETH mode)
{

	if (getFlg(mode)) return getValue(mode);

	vector<double> InitialSpreads;
	vector<double> recovery;
	vector<double> triggers_yt;
	vector<double> triggers_sp;
	vector<double> pdef;
	vector<double> disc;

	long seed = 0,nbdef = -1,row=0,col=0,i=0,nbtrigger = 0,nbloss = 0,nbnames=0;
	double initpv = -1.,initial_price=0.,Matu_Spread=-1.,loss = -1.,Simply = -1;
	double Time2Matu = -1,Leverage = 1.,feepv = 0.,defpv = 0.,Fast = 1.,volspreads=-1.,volcorrel=-1.;
	double spreadgap=0.,Correldw0=0.,CorrelUp0=0.,avgspread=0.,price=0.,	m_strike_dw=-1.,m_strike_up=-1.;
	double sumpv = 0.,pv=0.,rate=0.,VarReduction=-999.,LossPtf=0.,Corr00=0.;
	double pvcdo =0.;

	ICM_ModelMultiCurves* Modelo = NULL;
	ICM_Pricer*	Pricer=NULL;
	ICM_Mez* cdo_modified = NULL;

	LoadingParameters(seed,nbdef,initial_price,Matu_Spread,loss,
						Simply,Time2Matu,Leverage,Fast,volspreads,volcorrel,VarReduction);

	if (seed>=0) {NAG_random_init_repeatable(seed);}
	//if (volspreads>=0.) {its_vol_spreads=volspreads;}
	//if (volcorrel>=0.) {its_vol_correl=volcorrel;}

	its_vol_spreads=volspreads;
	its_vol_correl=volcorrel;

	ICM_Lss_Gap_Option* GapOption= (ICM_Lss_Gap_Option*) GetSecurity();
	ICM_Mez* underlying_cdo = GapOption->GetUnderlyingCdo();

	if (!initial_price)
	{
	PriceUnderlying();
	feepv = its_SS_Pricer->Price(qCMPFEELEGPV);
	defpv = its_SS_Pricer->Price(qCMPDEFLEGPV);
	initial_price = -feepv + Leverage*defpv;
	if (its_SS_Pricer) delete its_SS_Pricer;
	its_SS_Pricer=NULL;
	}

	// recuperation du notional initial
	// double notional0 = underlying_cdo->GetInitialNotional();
	double notional0 = underlying_cdo->GetFeeLeg()->GetCreditInfosRef().GetNotionals().Elt(0) ;
	notional0 *= underlying_cdo->GetTradedCoef();

	ICM_ModelMultiCurves* mod = (ICM_ModelMultiCurves*) GetModel();
	ARM_ZeroCurve* ircurve = mod->GetZeroCurve();
	
	ARM_Date StartDate = underlying_cdo->GetStartDateNA();
	ARM_Date EndDate = underlying_cdo->GetEndDateNA();
	long nbdays = EndDate-StartDate;
	long frequency = underlying_cdo->GetFeeLeg()->GetPaymentFreq();

	double YTMaturityCDO = (underlying_cdo->GetEndDateNA()-GetAsOfDate())/365.;
	double MaturitySimuls = (YTMaturityCDO-Time2Matu);
	if (MaturitySimuls<0.) return 0.;

	// recuperation de la matrice de triggers sur les spreads
	ICM_QMatrix<double>* Triggers_spreads = GapOption->GetSpreadTrigger();

	// recuperation des recovery
	//recovery.resize(underlying_cdo->GetCollateral()->GetNbIssuers());
	recovery.resize(1);
	for (int il=0;il<recovery.size();il++)
	{recovery[il]=mod->GetRecoveryRate(underlying_cdo->GetCollateral()->GetIssuersLabels(il));}

	//recuperation du spread trigger à la maturité
	spreadgap = GetSpreadLevel(Time2Matu,loss,row,col);

	//recuperation des triggers spreads jusqu'a maturité
	double TimeToTrigg =0.;
	for (i=1;i<=col;i++)
	{	TimeToTrigg = YTMaturityCDO-(*Triggers_spreads)(0,i);
		if (TimeToTrigg>=0.)
		{triggers_yt.push_back(TimeToTrigg);
		triggers_sp.push_back((*Triggers_spreads)(row,i));}}

	//generation du cdo impacté des losses
	if (nbdef>0) //cas nbdef défauts tirés aléatoirement
	{cdo_modified = GenerateCDOWithNNamesInDefault(underlying_cdo,nbdef,m_strike_dw,m_strike_up);
	}else //cas on impose loss% de défaut
	{cdo_modified = GenerateCDOImpactInLoss(underlying_cdo,loss,m_strike_dw,m_strike_up);}

	//give asof values for spreads & correlation
	FillInitalValues(cdo_modified,YTMaturityCDO,InitialSpreads,Correldw0,CorrelUp0,Simply);

	//fill initial correlation
	double corr00 = CptCorrForCDO(underlying_cdo,(*Triggers_spreads)(row,0));

	//calcul de la loss du portefeuille
	const ICM_DefaultCurve* dc = mod->GetDefaultCurve(underlying_cdo->GetCollateral()->GetIssuersLabels(0));
	double _pdef = dc->DefaultProba(MaturitySimuls);
	if ((MaturitySimuls<1.e-3) && (loss==0.))
		{LossPtf=1.;}
	else if (MaturitySimuls<1.e-3)
		{LossPtf=0.;}
	else
		{LossPtf =1. - ComputeEL_LHP((*Triggers_spreads)(row,0),sqrt(corr00),_pdef,recovery[0]);}

	//recuperation du notional apres impact des losses
	// 14514 double NotAfterLosses = cdo_modified->GetInitialNotional();
	double NotAfterLosses = cdo_modified->GetFeeLeg()->GetCreditInfosRef().GetNotionals().Elt(0) ;
	NotAfterLosses *= cdo_modified->GetTradedCoef();

	//fill initial correlation according to impacted CDO
	//FillInitalCorrels(underlying_cdo,Correldw00,CorrelUp00);

	//use only for fast case
	ICM_Security* security = cdo_modified->GetFeeLeg()->GetCreditInfos();
	rate = security->GetCouponRates().Elt(0)/100.;
	nbnames = cdo_modified->GetCollateral()->GetNbIssuers();

	double PSMT_S=0.;
	if (VarReduction==-999.) 
		{VarReduction=InitialSpreads[0];}
	else if (VarReduction==-998.) 
		{VarReduction=InitialSpreads[0];}

	if (VarReduction) {
	PSMT_S = 1. - NAG_cumul_normal(
		(1./(its_vol_spreads*sqrt(MaturitySimuls)))* 
		(log(VarReduction/InitialSpreads[0])+0.5*its_vol_spreads*its_vol_spreads*MaturitySimuls)
		);}

	ARM_Vector YFstartdates;
	ARM_Vector YFenddates;

	int stopflag = 0;

	for (int nbsim=0;nbsim<itsnbsimuls;nbsim++)
	{
		if (nbsim==stopflag)
		{stopflag=stopflag;}

		cdo_modified->ResetSchedules();
		double AsOfTriggerd = MaturitySimuls;

		vector<double> GenSpreads;
		double correls_dw =0.,correls_up =0.;
		bool istriggered=false;

		//build space of randoms values
		//if (Simply<0)
		//{// generation des markets data - toutes les courbes de defauts sont générées
		//GenerateValues(InitialSpreads,MaturitySimuls,its_drift_spreads,its_vol_spreads,
		//			GenSpreads,Correldw0,CorrelUp0,its_drift_correl,its_vol_correl,
		//			correls_dw,correls_up,avgspread);
		//}
		//else 
		//{// generation des markets data - 1 courbe de défaut est générée
		GenerateValuesLight(InitialSpreads,triggers_yt,triggers_sp,its_drift_spreads,its_vol_spreads,
					GenSpreads,Correldw0,CorrelUp0,its_drift_correl,its_vol_correl,
					correls_dw,correls_up,avgspread,AsOfTriggerd,istriggered,VarReduction);
		//}
		

		if (!Fast) //normal pricing - default mode
		{
			// construction d'un model rapide avec les default curves naïves 
			Modelo = BuildFastModel(ircurve,GenSpreads,recovery,m_strike_dw,
									m_strike_up,correls_dw,correls_up,AsOfTriggerd,
									Simply);

			ICM_Pricer_Advisor Advisor;
			Pricer = Advisor.GeneratePricer(cdo_modified, Modelo, ICM_PRICER_HOMOGENEOUS_SMILE,
										   CREDIT_DEFAULT_VALUE,&GetParameters(),GetAsOfDate()
										   );

			feepv = Pricer->Price(qCMPFEELEGPV);
			defpv = Pricer->Price(qCMPDEFLEGPV);
			pv = -feepv + Leverage*defpv;
		}
		else	// Homogeneous case
		{
			pdef.clear();disc.clear();

			// calcul des year fractions du schedule pour AsOfTriggered
			ARM_Date AsOfTrig = (ARM_Date)floor(GetAsOfDate().GetJulian() + AsOfTriggerd*365.);
			/*cdo_modified->GetFeeLeg()->ComputeYF(AsOfTrig);
			ICM_Security* security = cdo_modified->GetFeeLeg()->GetCreditInfos();
			ARM_Vector YFstartdates=security->GetYFAccStartDates();
			ARM_Vector YFenddates=security->GetYFAccEndDates();
			*/

			ARM_Date _start = StartDate;
			//ARM_Date _start = (ARM_Date)floor(GetAsOfDate().GetJulian() + MaturitySimuls*365.);
			ARM_Date _end = EndDate;

			GenerateScheduleYF(AsOfTrig, _start, _end,frequency,YFstartdates,YFenddates);

			// construction de la default curve naïve
			ICM_DefCurveSimple dc(GenSpreads[0],recovery[0],underlying_cdo->GetCollateral()->GetIssuersLabels(0));

			// determination pdef & discount price for year fractions
			for (i=0;i<YFstartdates.GetSize();i++)
			{pdef.push_back(dc.DefaultProba(YFstartdates.Elt(i)));
			disc.push_back(ircurve->DiscountPrice(YFstartdates.Elt(i)));}
			pdef.push_back(dc.DefaultProba(YFenddates.Elt(YFstartdates.GetSize()-1)));
			disc.push_back(ircurve->DiscountPrice(YFenddates.Elt(YFstartdates.GetSize()-1)));

			if CHECK_EQUAL(Fast,1.) // pricing en mode LHP 
				pv = -FastCDOPricing(&YFstartdates,&YFenddates,rate,m_strike_dw,m_strike_up,correls_dw,
								correls_up,NotAfterLosses,NotAfterLosses*Leverage,pdef,disc,recovery[0],
								nbnames,40,feepv,defpv,true);
			else					// pricing par la formule de récurence 
				pv = -FastCDOPricing(&YFstartdates,&YFenddates,rate,m_strike_dw,m_strike_up,correls_dw,
								correls_up,NotAfterLosses,NotAfterLosses*Leverage,pdef,disc,recovery[0],
								nbnames,40,feepv,defpv,false);

		}

		double DP = ircurve->DiscountPrice(AsOfTriggerd); // on actualise la pv en Asof
		double pv_delta_mtm_cdo = MAX(pv*DP -initial_price,0.);

		if ((VarReduction)&&(!CHECK_EQUAL(AsOfTriggerd,MaturitySimuls)))
		{	pv_delta_mtm_cdo = pv_delta_mtm_cdo;  }
		else if (VarReduction)
		{	pv_delta_mtm_cdo = pv_delta_mtm_cdo*PSMT_S; }
 
		pvcdo += pv*DP; 

		//application des triggers & description du payoff
		if ((spreadgap<avgspread) || (istriggered))
		{nbtrigger++;}

		if (pv_delta_mtm_cdo > notional0)
		{	nbloss++;
			price += (pv_delta_mtm_cdo-notional0);}


		if (Modelo) delete Modelo;
		if (Pricer) delete Pricer;

	}

	if (cdo_modified) delete cdo_modified;

	//quelle est la proba que le trigger se soit levé ? 
	its_proba_trigger=(((double)nbtrigger)/(double)(itsnbsimuls));
	setValue(qCMPFEELEGPV,its_proba_trigger);

	//quelle est la proba qu'il y ait une perte
	its_proba_loss=(((double)nbloss)/(double)(itsnbsimuls));
	if (VarReduction) {its_proba_loss = its_proba_loss*PSMT_S;}
	setValue(qCMPDEFLEGPV,its_proba_loss);

	double result = 0.;
	result = price/itsnbsimuls;
	pvcdo = pvcdo/itsnbsimuls;;

	//if (VarReduction) {result = result*PSMT_S;}

	setValue(qCMPPRICE,result);
	setValue(qCMPACCRUED,LossPtf);
	setValue(qCMPPREMIUM,pvcdo);

	return getValue(mode);
}

//-----------------------------------------------------------------------------
// Récupération des conditions de marché initiales
//-----------------------------------------------------------------------------
void ICM_Pricer_LssGapOption::FillInitalValues(ICM_Mez* cdo,
											   double matu_spread,
											   vector<double>& InitialSpreads,
											   double& Correldw0,
											   double& CorrelUp0,
											   double simplify)
{
	InitialSpreads.clear();

	ICM_ModelMultiCurves* mmc = (ICM_ModelMultiCurves*)GetModel();
	ICM_Correlation* correl = mmc->GetCorrelation();

	correl->ResetBetas();

	double strike_dw = cdo->GetPercentLow(cdo->GetFeeLeg()->GetStartDate());
	double strike_up = strike_dw + cdo->GetPercentHight(cdo->GetFeeLeg()->GetStartDate());

	double YF_Maturity = (cdo->GetEndDateNA()- GetAsOfDate())/365.;

	ICM_Pricer_Advisor Advisor;

	ICM_Pricer* pricer = Advisor.GeneratePricer(cdo,
										   GetModel(),
										   ICM_PRICER_HOMOGENEOUS_SMILE,
										   CREDIT_DEFAULT_VALUE,
										   &GetParameters(),GetAsOfDate());

	pricer->Price(qCMPPRICE);
	

	correl->SetForcedStrikeType(qStrike_LOW);
	Correldw0=	correl->GetBeta(cdo->GetCollateral()->GetIssuersLabels(0),YF_Maturity,CREDIT_DEFAULT_VALUE,YF_Maturity);
	Correldw0 *=Correldw0;

	correl->SetForcedStrikeType(qStrike_UP);
	CorrelUp0=	correl->GetBeta(cdo->GetCollateral()->GetIssuersLabels(0),YF_Maturity,CREDIT_DEFAULT_VALUE,YF_Maturity);
	CorrelUp0 *=CorrelUp0;

	InitialSpreads.resize(1);

	for (int i=0;i<InitialSpreads.size();i++)
	{
		if ((simplify>=2.)&&(i>=1))
		{InitialSpreads[i] = InitialSpreads[0];}
		else
		{
		const ICM_DefaultCurve* dc = mmc->GetDefaultCurve(cdo->GetCollateral()->GetIssuersLabels(i));
		InitialSpreads[i] = dc->ImpliedSpreadInterpol(cdo->GetEndDateNA());
		}
	}

	if (pricer) delete pricer;
}

//-----------------------------------------------------------------------------
// Récupération des conditions de marché initiales
//-----------------------------------------------------------------------------
void ICM_Pricer_LssGapOption::FillInitalCorrels(ICM_Mez* mez,
												double& Correldw0,
												double& CorrelUp0)
{

	ICM_ModelMultiCurves* mmc = (ICM_ModelMultiCurves*)GetModel();
	ICM_Correlation* correl = mmc->GetCorrelation();

	ICM_Mez* underlying_cdo = mez;

	correl->ResetBetas();

	double strike_dw = underlying_cdo->GetPercentLow(underlying_cdo->GetFeeLeg()->GetStartDate());
	double strike_up = strike_dw + underlying_cdo->GetPercentHight(underlying_cdo->GetFeeLeg()->GetStartDate());

	double YF_Maturity = (underlying_cdo->GetEndDateNA()- GetAsOfDate())/365.;
	

	correl->SetForcedStrikeType(qStrike_LOW);
	Correldw0=	correl->GetBeta(underlying_cdo->GetCollateral()->GetIssuersLabels(0),YF_Maturity,CREDIT_DEFAULT_VALUE,YF_Maturity);
	Correldw0 *=Correldw0;

	correl->SetForcedStrikeType(qStrike_UP);
	CorrelUp0=	correl->GetBeta(underlying_cdo->GetCollateral()->GetIssuersLabels(0),YF_Maturity,CREDIT_DEFAULT_VALUE,YF_Maturity);
	CorrelUp0 *=CorrelUp0;
}	

//-----------------------------------------------------------------------------
// Generation des spreads & des correlations suivant un brownien géométrique
//-----------------------------------------------------------------------------
void ICM_Pricer_LssGapOption::GenerateValues(vector<double> InitialSpreads,
											double Maturity,
											double drift_spreads,
											double vol_spreads,
											vector<double>& spreads,
											double Correldw0,
											double CorrelUp0,
											double drift_correls,
											double vol_correls,
											double& correls_dw,
											double& correls_up,
											double& avgspreads)
{
	/*
	spreads.clear();
	double spread=0.;

	spreads.resize(InitialSpreads.size());

	correls_dw=0.;
	correls_up=0.;

	avgspreads = 0.;

	//on tire les spreads
	for (int k=0;k<InitialSpreads.size();k++)
	{
		spreads[k] = InitialSpreads[k] + 
					InitialSpreads[k]*drift_spreads*Maturity + 
					InitialSpreads[k]*vol_spreads*sqrt(Maturity)*NAG_random_normal(0,1);

		spreads[k] = MAX(spreads[k],0);

		avgspreads += spreads[k];
	}

	avgspreads /= (double) InitialSpreads.size();

	//on tire la correlation down
	correls_dw +=	Correldw0 +	
					Correldw0*drift_spreads*Maturity + 
					Correldw0*vol_spreads*sqrt(Maturity)*NAG_random_normal(0,1);

	//on tire la correlation up
	correls_up +=	CorrelUp0 +	
					CorrelUp0*drift_spreads*Maturity + 
					CorrelUp0*vol_spreads*sqrt(Maturity)*NAG_random_normal(0,1);
	
	correls_dw = MAX(correls_dw,0);
	correls_up = MAX(correls_up,0);
	*/
}

//-----------------------------------------------------------------------------
// Generation d'un cdo ayant N issuers en défauts
//-----------------------------------------------------------------------------
ICM_Mez* ICM_Pricer_LssGapOption::GenerateCDOWithNNamesInDefault(ICM_Mez* initialCDO,
																 int nbdefaults,
																 double& strike_dw,
																 double& strike_up)
{
	ICM_Collateral* collat = initialCDO->GetCollateral();
	int nbnames = collat->GetNbIssuers();
	double random = 0.;
	vector<string> defaults;
	double already = false;

	for (int i=0;i<nbdefaults;i++)
	{
		random=(int)(NAG_random_continuous_uniform()*nbnames);
		string name = collat->GetIssuersLabels(i);
		already = false;

		for (int j=0;j<defaults.size();j++)
		{
			if (defaults[j]==name)
			{
				i--;
				already = true;
				break;
			}
		}

		if (!already)
		{
			defaults.push_back(name);
			nbnames--;
		}
	}

	ICM_Mez* CDO= (ICM_Mez*)initialCDO->Clone();

	for (int k=0;k<defaults.size();k++)
	{CDO->ExcludeIssuer(defaults[k]);}

	strike_dw = CDO->GetPercentLow(CDO->GetFeeLeg()->GetStartDate());
	strike_up = CDO->GetPercentLow(CDO->GetFeeLeg()->GetStartDate()) + CDO->GetPercentHight(CDO->GetFeeLeg()->GetStartDate());

	return (CDO);
}

//-----------------------------------------------------------------------------
// Build Quckly Model Multi Curves
//-----------------------------------------------------------------------------
ICM_ModelMultiCurves* ICM_Pricer_LssGapOption::BuildFastModel(ARM_ZeroCurve* ircurve,
														vector<double> spreads,
														vector<double> recovery,
														double Strike_Dw,
														double Strike_Up,
														double Correl_Dw,
														double Correl_Up,
														double AsOfTriggered,
														double quick)
{
	int nbcurves = spreads.size();	

	if (quick>=0)
	{ nbcurves = 1; }

	std::vector<const ICM_DefaultCurve*> DefaultCurves ( nbcurves );
	ARM_Vector RecoveryRates (nbcurves);

	ICM_Lss_Gap_Option* GapOption= (ICM_Lss_Gap_Option*) GetSecurity();
	ICM_Mez* underlying_cdo = GapOption->GetUnderlyingCdo();

	ARM_Date AsOf = (ARM_Date)floor(GetAsOfDate().GetJulian() + AsOfTriggered*365.);

	int i=0;

	for (i=0; i<nbcurves;i++)
	{RecoveryRates[i]=recovery[i];}		

	for (int j=0; j<nbcurves;j++)
	{DefaultCurves[j]=new ICM_DefCurveSimple(spreads[j],
											recovery[j],
											underlying_cdo->GetCollateral()->GetIssuersLabels(j));}		

	ICM_Smile_Correlation* correl = FixedBaseCorrelation(AsOf,
														underlying_cdo->GetCurrencyUnit(),
														Strike_Dw,
														Strike_Up,
														Correl_Dw,
														Correl_Up,
														"TRAX");

	ARM_ZeroCurve* ir = (ARM_ZeroCurve*) ircurve->Clone();
	ir->SetAsOfDate(AsOf);

	ICM_ModelMultiCurves*  MMC = new ICM_ModelMultiCurves(/** nbcurves, **/ 
														DefaultCurves,
														ir,
														RecoveryRates,
														correl);
	// if (RecoveryRates) delete[] RecoveryRates;
	if (ir) delete ir;
	if (correl) delete correl;


	// if (DefaultCurves)
	// {
	for (i=0; i<nbcurves;i++)
	{delete[] DefaultCurves[i];}		
	
		// delete[] DefaultCurves;
	// }

	return MMC;

}

//-----------------------------------------------------------------------------
// Generation d'un cdo ayant subi N% losses
//-----------------------------------------------------------------------------
ICM_Mez* ICM_Pricer_LssGapOption::GenerateCDOImpactInLoss(ICM_Mez* initialCDO,
														double loss,
														double& strike_dw,
														double& strike_up)
{
	ICM_Mez* CDO= (ICM_Mez*)initialCDO->Clone();

	double percent_dw= CDO->GetPercentLow(CDO->GetFeeLeg()->GetStartDate());
	double percent_up= percent_dw + CDO->GetPercentHight(CDO->GetFeeLeg()->GetStartDate());

	if (percent_dw>0)
	{ percent_dw=MAX(percent_dw-loss,0.);}

	if (percent_up>0)
	{ percent_up=MAX(percent_up-loss,0.);}

	strike_dw = percent_dw;
	strike_up = percent_up;

	double total = CDO->GetCollateral()->SumNotionals(CDO->GetStartDateNA());
	double Subamount = percent_dw*total;
	double Mezzamount = (percent_up-percent_dw)*total;

	CDO->SetSubAmount(Subamount);
	CDO->SetMezzAmount(Mezzamount);

	return (CDO);
}

//-----------------------------------------------------------------------------
// Generation d'un cdo ayant subi N% losses
//-----------------------------------------------------------------------------
double ICM_Pricer_LssGapOption::CptCorrForCDO(ICM_Mez* initialCDO,double strike)
{
	ICM_Mez* CDO= (ICM_Mez*)initialCDO->Clone();

	double total = CDO->GetCollateral()->SumNotionals(CDO->GetStartDateNA());
	double Subamount = 0.;
	double Mezzamount = strike*total;

	CDO->SetSubAmount(Subamount);
	CDO->SetMezzAmount(Mezzamount);

	ICM_ModelMultiCurves* mmc = (ICM_ModelMultiCurves*)GetModel();
	ICM_Correlation* correl = mmc->GetCorrelation();

	correl->ResetBetas();

	ICM_Pricer_Advisor Advisor;

	ICM_Pricer* pricer = Advisor.GeneratePricer(CDO,
										   GetModel(),
										   ICM_PRICER_HOMOGENEOUS_SMILE,
										   CREDIT_DEFAULT_VALUE,
										   &GetParameters(),GetAsOfDate());

	pricer->Price(qCMPPRICE);

	double YF_Maturity = (CDO->GetEndDateNA()- GetAsOfDate())/365.;
	
	correl->SetForcedStrikeType(qStrike_UP);
	double corr=	correl->GetBeta(CDO->GetCollateral()->GetIssuersLabels(0),YF_Maturity,CREDIT_DEFAULT_VALUE,YF_Maturity);
	corr *=corr;

	if (CDO) delete CDO;
	if (pricer) delete pricer;

	return (corr);
}

//-----------------------------------------------------------------------------
// Generation des spreads & des correlations suivant un brownien géométrique
//-----------------------------------------------------------------------------
void ICM_Pricer_LssGapOption::GenerateValuesLight(vector<double>& InitialSpreads,
											vector<double>& Maturity,
											vector<double>& Triggers_spreads,
											double drift_spreads,
											double vol_spreads,
											vector<double>& spreads,
											double Correldw0,
											double CorrelUp0,
											double drift_correls,
											double vol_correls,
											double& correls_dw,
											double& correls_up,
											double& avgspreads,
											double& MatuTriggered,
											bool& istriggered,
											double VarReduction)
{
	spreads.clear();
	double spread=0.;

	spreads.resize(InitialSpreads.size());

	correls_dw=0.;
	correls_up=0.;
	int i=0;

	vector<double> norm;
	for (i=0;i<Maturity.size();i++)
	{norm.push_back(NAG_random_normal(0,1));}

	istriggered = false;
	int _max=Maturity.size();

	MatuTriggered = Maturity[Maturity.size()-1];

	//on tire les spreads
	for (int k=0;k<InitialSpreads.size();k++)
	{
		for (i=0;i<_max;i++)
		{
			if (istriggered==false)
			{
			if (i==0) 
				{spreads[k] = St(drift_spreads,vol_spreads,InitialSpreads[k],0.,Maturity[i],norm[i]);}
			else
				{spreads[k] = St(drift_spreads,vol_spreads,spreads[k],Maturity[i-1],Maturity[i],norm[i]);}
			}

			if ((spreads[k]>Triggers_spreads[i])&&
				(istriggered==false)&&
				(i<_max-1))
			{MatuTriggered=Maturity[i];istriggered = true;break;}

			if ((i==_max-1)&&(istriggered==false))
			{ 
			if (spreads[k]>Triggers_spreads[i])
				{istriggered = true;}

			if ((VarReduction)&&(spreads[k]<VarReduction)) 
				{while (spreads[k]<VarReduction)
					{spreads[k] = St(drift_spreads,vol_spreads,InitialSpreads[k],0.,Maturity[i],NAG_random_normal(0,1));}}
			}
		}
	}

	avgspreads = spreads[0];

	//on tire la pente de la correl suivant la même vol que la correlation
	double ecart_absolu = (CorrelUp0-Correldw0);
		//on rejete les pentes<0
	if (ecart_absolu<0.) { ecart_absolu = 0.; }

	ecart_absolu = St(drift_correls,vol_correls,ecart_absolu,0.,MatuTriggered,NAG_random_normal(0,1));
	correls_dw = St(drift_correls,vol_correls,Correldw0,0.,MatuTriggered,NAG_random_normal(0,1));

	//on tire la correlation up
	correls_up =	correls_dw + ecart_absolu;
	
	correls_dw = MIN(MAX(correls_dw,0),0.99);
	correls_up = MIN(MAX(correls_up,0),0.99);
}


//-----------------------------------------------------------------------------
// return fixed spread level for tim2matu & loslevel
//-----------------------------------------------------------------------------
double ICM_Pricer_LssGapOption::GetSpreadLevel(double time2matu,double loslevel,long& row, long& col)
{
	double TOL = 7./365.;

	ICM_Lss_Gap_Option* GapOption= (ICM_Lss_Gap_Option*) GetSecurity();
	ICM_QMatrix<double>* Triggers_spreads = GapOption->GetSpreadTrigger();

	row=-1;	col=-1;

	for (int i=1;i<Triggers_spreads->Getnbrows();i++)
	{double  value = (*Triggers_spreads)(i,0);

		if (fabs(value-loslevel)<=1.e-3)
		{row = i;break;}}

	for (int j=1;j<Triggers_spreads->Getnbcols();j++)
	{double  value = (*Triggers_spreads)(0,j);
		if (fabs(value-time2matu)<=TOL)
		{col = j;break;}}

	if ((row>=0)&&(col>=0))
		return (*Triggers_spreads)(row,col);
	else
		return CREDIT_DEFAULT_VALUE;
}