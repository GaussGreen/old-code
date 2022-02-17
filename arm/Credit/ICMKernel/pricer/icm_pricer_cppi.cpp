#include "firsttoinc.h"
#include "icm_pricer_cppi.h"
#include "ICMKernel/inst/icm_cppi.h"
#include "ICMKernel/util/icm_diffusion.h"
#include <numeric>
#include "nag.h"
#include <nags.h>
#include <nagg01.h>
#include "nagg05.h"
#include <ICMKernel\glob\icm_smile_correlation.h>
#include <ICMKernel\pricer\icm_pricer_homogeneous_smile.h>
#include <ICMKernel\inst\icm_mez.h>
#include <ARMKernel\crv\volflat.h>
#include <ARMKernel\crv\volint.h>
#include "ICMKernel\crv\icm_default_str.h"
#include "ICMKernel/glob/icm_mktdatamng.h"
#include "ICMKernel\inst\icm_collateral.h"
#include "ICMKernel/mod/modelmulticurves.h"


//------------------------------
#ifndef PI
#define PI		3.14159265359
#endif

void ICM_Pricer_Cppi::Init()
{
	SetName(ICM_PRICER_CPPI);

	itsFunding = 0.; 
	itsFinancement = 0.;
	itsRunning = 0.; 
	itsLambda = 0.;
	itsBidMid = 0.;
	itsSwitch = 0.;
	itsR0 = 0.;
	itsIRVol = 0.;
	itsA = 0.;
	itsNbSimul = 10;
	itsStep = 1./365.;
	itsIndexIndice = 0;
	itsVolImplied = 0.;
	itsMeanImplied = 1.;
	itsSpeedReversionImplied = 0.;
	itsImplied.clear();
	itsE.clear();
	itsf.clear();
	itsN.clear();
	itsV.clear();
	itsP.clear();
	itsAlphaMin.clear();
	itsAlphaMax.clear();
	itsAlpha.clear();
	itsNbRebal.clear();
	itsZeroCoupon.clear();
	itsDur.clear();
	itsSoult.clear();
	itsR.clear();
	itsS.clear();
	itsTimeToRoll.clear();
	itsJump.clear();
	itsSaut.clear();
	itsSigmaV.clear();
	itsMeanV.clear();
	itsWidener.clear();
	itsNbCouponPayed.clear();
	itsBarriereDownTime.clear();
	itsDistribTRI = NULL;

	itsSpreadCoupon = 0.;
	itsPrimeUpFront = 0.;
	itsFactorCoupon = 0.;
	itsFixedEurib = 0;
	itsStopMaturity =0.;

	itsInitWidener = 0.;
	itsMuWidener = 0.;
	itsSigmaWidener = 0.;
	itsPartWidener = 0.;

	itsBarriereDown = 0.0;
	itsLockedTime = 0.0;
}

void ICM_Pricer_Cppi::Set(ARM_Security *sec, ARM_Object *mod,const ICM_Parameters& parameters,const ARM_Date&asof)
{
	ICM_Pricer::Set(sec,mod,parameters,&asof);

	// if (parameters != NULL)
	if (!parameters.empty())
	{
		ARM_Vector* Funding = parameters.GetColVect("FUNDING");
		if (Funding)		{itsFunding		= Funding->Elt(0);}

		ARM_Vector* Financement = parameters.GetColVect("FINANCEMENT");
		if (Financement)	{itsFinancement		= Financement->Elt(0);}

		ARM_Vector* Running = parameters.GetColVect("RUNNING");
		if (Running)		{itsRunning		= Running->Elt(0);}

		ARM_Vector* IntensitePoisson = parameters.GetColVect("INTENSITE_POISSON");
		if (IntensitePoisson) {itsLambda = IntensitePoisson->Elt(0);}

		ARM_Vector* BidMid = parameters.GetColVect("SOULT_BIDMID");
		if (BidMid) {itsBidMid = BidMid->Elt(0);}
		
		ARM_Vector* Switch = parameters.GetColVect("SOULT_SWITCH");
		if (Switch) {itsSwitch = Switch->Elt(0);}
		
		ARM_Vector* Taux = parameters.GetColVect("TAUX");
		if (Taux)	{itsR0 = Taux->Elt(0);}

		ARM_Vector* Vasicek_param = parameters.GetColVect("VASICEK_PARAM");
		if (Vasicek_param)	{itsA = Vasicek_param->Elt(0);}
	
		ARM_Vector* Vasicek_Vol = parameters.GetColVect("VASICEK_VOL");
		if (Vasicek_Vol)	{itsIRVol = Vasicek_Vol->Elt(0);}

		ARM_Vector* MC_Pas = parameters.GetColVect("MC_PAS");
		if (MC_Pas)	{itsStep = MC_Pas->Elt(0);}

		ARM_Vector* MC_NbSimul = parameters.GetColVect("MC_NB_SIMUL");
		if (MC_NbSimul)	{itsNbSimul = MC_NbSimul->Elt(0);}

		ARM_Vector* Indice = parameters.GetColVect("INDICE");
		if (Indice)	{itsIndexIndice = Indice->Elt(0);}

		ARM_Vector* VolImplied = parameters.GetColVect("VOL_IMPLIED");
		if (VolImplied)	{itsVolImplied = VolImplied->Elt(0);}

		ARM_Vector* MeanImplied = parameters.GetColVect("MEAN_IMPLIED");
		if (MeanImplied)	{itsMeanImplied = MeanImplied->Elt(0);}

		ARM_Vector* SpeedReversionmplied = parameters.GetColVect("SPEED_REVERSION_IMPLIED");
		if (SpeedReversionmplied)	{itsSpeedReversionImplied = SpeedReversionmplied->Elt(0);}

		ARM_Vector* ZCoupon = parameters.GetColVect("CRB_TAUX");
		if (ZCoupon)
		{
			int size = ZCoupon->GetSize();
			itsZeroCoupon.resize(size);
			for (int i = 0; i< size;i++)
			{itsZeroCoupon[i] = ZCoupon->Elt(i);}
		}

		ARM_Vector* SpreadCoupon = parameters.GetColVect("SPREAD_COUPON");
		if (SpreadCoupon)	{itsSpreadCoupon = SpreadCoupon->Elt(0);}
		
		ARM_Vector* FactorCoupon = parameters.GetColVect("FACTOR_COUPON");
		if (FactorCoupon)	{itsFactorCoupon = FactorCoupon->Elt(0);}

		ARM_Vector* PrimeUpFront = parameters.GetColVect("PRIME_UPFRONT");
		if (PrimeUpFront)	{itsPrimeUpFront = PrimeUpFront->Elt(0);}

		ARM_Vector* StopMaturity = parameters.GetColVect("STOP_MAT");
		if (StopMaturity)	{itsStopMaturity = StopMaturity->Elt(0);}

		ARM_Vector* FixedEurib = parameters.GetColVect("FIXED_EURIB");
		if (FixedEurib)	{itsFixedEurib = FixedEurib->Elt(0);}

		ARM_Vector* InitWidener = parameters.GetColVect("INIT_WIDENER");
		if (InitWidener)	{itsInitWidener = InitWidener->Elt(0);}

		ARM_Vector* MuWidener = parameters.GetColVect("MU_WIDENER");
		if (MuWidener)	{itsMuWidener = MuWidener->Elt(0);}
		
		ARM_Vector* SigmaWidener = parameters.GetColVect("SIGMA_WIDENER");
		if (SigmaWidener)	{itsSigmaWidener = SigmaWidener->Elt(0);}

		ARM_Vector* PartWidener = parameters.GetColVect("PART_WIDENER");
		if (PartWidener)	{itsPartWidener = PartWidener->Elt(0);}

		ARM_Vector* BarriereDown = parameters.GetColVect("BARRIERE_DOWN");
		if (BarriereDown)	{itsBarriereDown = BarriereDown->Elt(0);}

		ARM_Vector* LockedTime = parameters.GetColVect("LOCKED_TIME");
		if (LockedTime)	{itsLockedTime = LockedTime->Elt(0);}
	}
}




double ICM_Pricer_Cppi::ComputePrice(qCMPMETH measure )
{	
	ICM_Cppi* cppi = (ICM_Cppi*)GetSecurity();
	ICM_Mez* mez = (ICM_Mez*)cppi->GetSecurity();

	// Acquisition des paramètres propre au cppi
	ARM_Date StartDate = cppi->GetStartDate();
	double T = (cppi->GetEndDate().GetJulian() - StartDate.GetJulian())/365.;
	double V0 = cppi->GetNotional();
	double Garanti		 = cppi->GetProtectedAmount();
	double ManagementCost= cppi->GetManagementCost();
	double AddLeverage   = cppi->GetAdditionalLeverage();
	double DCushion		 = cppi->GetDesactivateCushion();
	ICM_RangeFactor* rfactor = cppi->GetRangeFactor();
	vector<double> MinSpread = rfactor->GetInf();
	vector<double> MaxSpread = rfactor->GetMax();;
	vector<double> MidValue  = rfactor->GetValueMid();
	vector<double> MinValue  = rfactor->GetValueMin();
	vector<double> MaxValue  = rfactor->GetValueMax();
	ARM_Currency* CCy= cppi->GetCurency();
	
	// Acquisition des parametres propres au pricer
	double Funding		= itsFunding;
	double Financement	= itsFinancement;
	double Running		= itsRunning;
	double Lambda		= itsLambda;
	double BidMid		= itsBidMid;
	double Switch		= itsSwitch;
	double r0			= itsR0;
	double IRVol		= itsIRVol;
	double a			= itsA;
	double MeanImplied  = log(itsMeanImplied);
	double Speed		= itsSpeedReversionImplied;
	double VolImplied	= itsVolImplied;
	double BarriereTaux	= itsBarriereDown;
	int	   TimeLocked	= (int)itsLockedTime;
	
	double SpreadCoupon = itsSpreadCoupon;
	double PrimeUpFront = itsPrimeUpFront;
	double FactorCoupon = itsFactorCoupon;
	int StopMat		= (int)itsStopMaturity;
	if (abs(StopMat)<1e-6)
		StopMat = (int)T;

	int FixedEurib		= itsFixedEurib;

	double InitWidener  = itsInitWidener;
	double MuWidener	= itsMuWidener;
	double SigmaWidener = itsSigmaWidener;
	double PartWidener	= itsPartWidener;
	// --------------------------------------------------------------------------------------	
	// Maturité de la courbe
	ARM_Date EndDate = mez->GetEndDateNA();
	// --------------------------------------------------------------------------------------	
	// Date de Roll (choix entre Septembre et Mars. Même année que StartDate);
	int yearDate = StartDate.GetYear();

	ARM_Date DateRollMarch = StartDate;
	DateRollMarch.ChangeDate(20,3,yearDate);
	
	ARM_Date DateRollSeptember = StartDate;	
	DateRollSeptember.ChangeDate(20,9,yearDate);
	
	ARM_Date rollDate = DateRollMarch;
	if (StartDate>DateRollMarch)
		rollDate = DateRollSeptember;
	// --------------------------------------------------------------------------------------	
	// On recupere les infos présentent dans le mktdatamanager
	ICM_MktDataMng* MktDataMng = GetMktDataMng();
	ICM_DefaultCurve* DefCurveTraxx = MktDataMng->GetDefaultCurve(mez->GetCollateral()->GetIssuersLabels(0),StartDate);
	ICM_DefaultCurve* DefCurveTraxxTemp = NULL;
	ARM_ZeroCurve* Zerocurve = MktDataMng->GetZeroCurve(CCy->GetCcyName(),StartDate);
	// --------------------------------------------------------------------------------------
	// On construit le vecteur Zero Coupon utile pour le Vasicek
	double taux = 0.;
	int sizeZeroCoupon =  (int)4*T; 
	if (itsZeroCoupon.size() == 0)
	{
		itsZeroCoupon.resize(sizeZeroCoupon);
		for (int iZc=0;iZc<sizeZeroCoupon;iZc++)
		{
			taux = Zerocurve->DiscountYield((iZc+1.0)*0.25);
			itsZeroCoupon[iZc] = exp(-taux * (iZc+1.0)*0.25);
		}
	}
	// --------------------------------------------------------------------------------------
	// Objet Correlation
	ICM_Correlation* Corr = MktDataMng->GetCorrelation(cppi->GetCorrelName(),StartDate);
	// --------------------------------------------------------------------------------------
	// Courbe Spread
	int NbIssuers = mez->GetCollateral()->GetNbIssuers();
	int sizeRange = MaxSpread.size();
	vector<double> MaxImplied; 	MaxImplied.resize(sizeRange);
	vector<double> MinImplied;	MinImplied.resize(sizeRange);

	int sizeSpread = DefCurveTraxx->GetRates()->GetSize() - 1; // On ne prend pas le plot 0
	ARM_Vector NewSpread = ARM_Vector(sizeSpread,1.);
	ARM_Vector Plot = ARM_Vector(sizeSpread,1.);
	for (int ind2 = 0 ; ind2<sizeSpread;ind2++)
	{
		NewSpread.Elt(ind2) = DefCurveTraxx->GetRate(ind2 + 1);
		Plot.Elt(ind2) = DefCurveTraxx->GetYearTerm(ind2 +1 );
	}

	ARM_Vector NewSpreadLow = ARM_Vector(NewSpread);
	ARM_Vector NewSpreadHigh = ARM_Vector(NewSpreadLow);
	
	// char* Label    = DefCurveTraxx->GetLabel();
	std::vector<const ICM_DefaultCurve*> TabCrbDefaut (1); 
	// ICM_DefaultCurve** TabCrbDefaut = new ICM_DefaultCurve*[1];
	ARM_Vector RecovTraxx (1) ;RecovTraxx[0] = DefCurveTraxx->GetRecovery();

	// Pour le pricing de l'implied Init
	ICM_ModelMultiCurves* ModelInit = NULL;;
	// ICM_Pricer_Distrib_Smile* Pricer_MezInit = NULL;
	TabCrbDefaut[0] = DefCurveTraxx;
	ModelInit     = new ICM_ModelMultiCurves(TabCrbDefaut ,/** 1,TabCrbDefaut,**/ Zerocurve,RecovTraxx,Corr);
	ICM_Pricer_Distrib_Smile Pricer_MezInit; 
	Pricer_MezInit.Set(mez,ModelInit,ICM_Parameters(),GetAsOfDate());
	double ImpliedInit = Pricer_MezInit.ComputeSpread(0)/10000.0;

	if (ModelInit)
		delete ModelInit;
	ModelInit = NULL;
	// if (Pricer_MezInit)
	// 	delete Pricer_MezInit;
	// Pricer_MezInit = NULL;
	// ------------------------------------------------
	ICM_ModelMultiCurves* Model = NULL;;
	ICM_Pricer_Distrib_Smile* Pricer_Mez = NULL;;

	double deltaLowSpread = 0.,deltaHighSpread = 0.;
	for (int ival =0; ival<sizeRange;ival ++)
	{
		// On considere en dur le fait que l'utilisateur entre toujours des plots 1 3 5 7 10
		// On recupere l'ecart entre le spread 5 ans du traxx actuel et le spread 5 ans inputter dans les colonnes MinSpread et MaxSpread
		deltaLowSpread  = MinSpread[ival] - NewSpread.Elt(2);
		deltaHighSpread = MaxSpread[ival] - NewSpread.Elt(2);
		for (int ind = 0;ind < sizeSpread; ind ++)
		{
			NewSpreadLow.Elt(ind) = NewSpread.Elt(ind) + deltaLowSpread;
			if (NewSpreadLow.Elt(ind)<0)
				NewSpreadLow.Elt(ind) = 0.0001;

			NewSpreadHigh.Elt(ind) = NewSpread.Elt(ind) + deltaHighSpread;
			if (NewSpreadHigh.Elt(ind)<0)
				NewSpreadHigh.Elt(ind) = 0.0001;
		}

		// Implied Min
		DefCurveTraxxTemp = new ICM_DefCurvStr(&NewSpreadLow,&Plot,RecovTraxx[0],Zerocurve,DefCurveTraxx->GetLabel());
		TabCrbDefaut[0] = DefCurveTraxxTemp;
		Model     = new ICM_ModelMultiCurves(TabCrbDefaut ,Zerocurve,RecovTraxx,Corr);
		ICM_Pricer_Distrib_Smile Pricer_Mez; Pricer_Mez.Set(mez,Model,ICM_Parameters(),GetAsOfDate());
		MinImplied[ival] = Pricer_Mez.ComputeSpread(0)/10000.;
		if (DefCurveTraxxTemp)
			delete DefCurveTraxxTemp;
		DefCurveTraxxTemp = NULL;
		if (Model)
			delete Model;
		Model = NULL;
		// if (Pricer_Mez)
		// 	delete Pricer_Mez;
		// Pricer_Mez = NULL;
	
		// Implied Max
		DefCurveTraxxTemp = new ICM_DefCurvStr(&NewSpreadHigh,&Plot,RecovTraxx[0],Zerocurve,DefCurveTraxx->GetLabel());
		TabCrbDefaut[0] = DefCurveTraxxTemp;
		Model     = new ICM_ModelMultiCurves(TabCrbDefaut,Zerocurve,RecovTraxx,Corr);
		ICM_Pricer_Distrib_Smile Pricer_Mez2; Pricer_Mez2.Set(mez,Model,ICM_Parameters(),GetAsOfDate());
		MaxImplied[ival] = Pricer_Mez2.ComputeSpread(0)/10000.;

		if (DefCurveTraxxTemp)
			delete DefCurveTraxxTemp;
		DefCurveTraxxTemp = NULL;
		if (Model)
			delete Model;
		Model = NULL;
		// if (Pricer_Mez)
		// 	delete Pricer_Mez;
		// Pricer_Mez = NULL;
	}
	// --------------------------------------------------------------------------------------
	// Range Factor contenant les implied de tranche 
	ICM_RangeFactor* riskFactorImplied = new ICM_RangeFactor(MinImplied,MaxImplied,MinValue,MaxValue);
	// --------------------------------------------------------------------------------------
	// Actif
	int index = 0.;	
	double IssuerNotional = mez->GetCollateral()->GetIssuersNotional(0);
	double JumpInit = - (1 - RecovTraxx[0])* IssuerNotional/mez->GetMezzAmount(mez->GetFeeLeg()->GetStartDate());
	double Jump = JumpInit;
	double RollPeriodInit = (rollDate - StartDate)/365.0, RollPeriod = 0.5;
	double RefMat = 0.;
	double mat = (EndDate.GetJulian() - StartDate.GetJulian())/365.;
	if (fabs(mat - 5)<0.1)
		RefMat = 4.75;
	else if (fabs(mat - 7)<0.1)
		RefMat = 6.75;
	else if (fabs(mat - 10)<0.1)
		RefMat = 9.75;

	double Implied = 0., ImpliedVasicek = 0.,logImpliedVasicek = 0.;
	// --------------------------------------------------------------	
	// Taux : On reconstruit un vecteur de Zéro coupon en prenant en consideration le plot T = 0 
	// et un plot supplementaire pour les calculs de Vasicek
	vector<double> ZCCurve; 	ZCCurve.resize(itsZeroCoupon.size() + 2);
	
	ZCCurve[0] = 1.;
	for (int il2 = 0; il2<itsZeroCoupon.size(); il2++)
		ZCCurve[il2+1] = itsZeroCoupon[il2];

	ZCCurve[itsZeroCoupon.size() + 1] = ZCCurve[itsZeroCoupon.size()];
	// --------------------------------------------------------------
	// nbSimul
	double dt = itsStep;
	int nbSimul = itsNbSimul;
	int nbstep = (int)(StopMat/dt);
	// --------------------------------------------------------------
	// Variables divers
	bool liquidate = false; // Indique si le fond est liquidé ou pas 
	int NbLoss = 0, CountRebal = 0., NbCouponPayed = 0.;
	double AlphaTheo = 0.,Alpha = 0.;
	double Loss = 0., Tdesactivation = 0. , VFinal = 0.;
	double Nmin = 0.;
	double Rebal = 0., NRebalOffer = 0., NRebalBid = 0.;
	double r = 0., s = 0.,  ZCR = 0., ZCRstop = 0.;
	double TimeToRoll = 0.,Dur = 0.;
	double Soult = 0., SoultNew = 0., SoultBid = 0., SoultOffer = 0.,SoultNewBid = 0.; 
	double PrixBid = 0., PrixOffer = 0.;
	double N = 0.,Etheo = 0., E = 0., f = 0., V = 0., P = 0.;	
	double Eurib1Y = 0.;
	double Eonia = 0.;
	double Widener = 0.; // permet de simuler un nom qui s'ecarte
	bool Barriere = false;
	double BarriereDown = 0.0;
	vector<double> TauxCoupon; // Taux euribor en vigueur lorsque l'on détache le coupon
	TauxCoupon.resize(15);
	double MeanTaux = 0.;
	double ProductTaux = 0.;
	double Vtot = 0.;
	vector<double> C; C.resize(itsZeroCoupon.size());
	double ValoCoupon = 0.0;

	double rdt = 0.;
	double TRI = 0.;
	// --------------------------------------------------------------
	// Loi normale centrée réduite
	double W = 0., ZImplied = 0., h = 0.0, x = 0.0, dY = 0., barriere = 0.;
	vector<double> vZimplied;	vZimplied.resize(nbstep);
	vector<double> vdY;			vdY.resize(nbstep);
	vector<double> vW;			vW.resize(nbstep);
	vector<double> vWidener;	vWidener.resize(nbstep);

	ICM_QMatrix<double>* ValTri =  new ICM_QMatrix<double>(nbSimul,StopMat,0.);
	// --------------------------------------------------------------
	// Si l'utilisateur veut visualiser la vie d'une trajectoire
	bool Traj = false;
	itsNbCouponPayed.resize(StopMat + 2);
	if (itsNbSimul == 1)
	{
		itsZImplied.resize(nbstep);		itsImplied.resize(nbstep);
		itsE.resize(nbstep);			itsf.resize(nbstep);
		itsN.resize(nbstep);			itsV.resize(nbstep);
		itsP.resize(nbstep);			itsAlphaMin.resize(nbstep);
		itsAlphaMax.resize(nbstep);		itsAlpha.resize(nbstep);
		itsNbRebal.resize(nbstep);		itsDur.resize(nbstep);
		itsR.resize(nbstep);			itsS.resize(nbstep);
		itsSoult.resize(nbstep);		itsTimeToRoll.resize(nbstep);
		itsJump.resize(nbstep);			itsSaut.resize(nbstep);
		itsWidener.resize(nbstep);		
		Traj =true;
	}
	else 
	{
		itsNbRebal.resize(3);
		itsNbRebal[0] = 0;		// on stocke la moyenne
		itsNbRebal[1] = 0;		// le max
		itsNbRebal[2] = 1e6;	// le min
		itsSigmaV.resize(StopMat);
		itsMeanV.resize(StopMat);
		if (itsDistribTRI)
			delete itsDistribTRI;
		itsDistribTRI = new ICM_QMatrix<double>(nbstep,StopMat);
		itsBarriereDownTime.resize((StopMat-(int)TimeLocked)*4);  // On sépare la durée en trimestre
	}
	// --------------------------------------------------------------
	// --------------------------------------------------------------
	// Simulation
	for (int i = 0;i<nbSimul;i++)
		{
			// --------------------------------------------------------------
			// Initialisation
			CountRebal = 0;
			index = itsIndexIndice;
			ImpliedVasicek = ImpliedInit;
			double ImpliedVasicek2 = (MinImplied[index] + MaxImplied[index])/2.;
			MeanImplied = log(itsMeanImplied);
			Widener = InitWidener;
			Implied = ImpliedVasicek + PartWidener * Widener;
			AlphaTheo = MidValue[index];
			MeanTaux = 0.;
			ProductTaux = 1.;
			Vtot = 0.;
			ValoCoupon = 0.;
			for (int indiceCoup=0;indiceCoup<C.size();indiceCoup++)
				C[indiceCoup] = 1.0;

			rdt = 0.;

			Barriere = false;
			
			r = r0;
			s = 0.;
			Jump = JumpInit;
			V = V0;
			ZCR = TauxVasicek(0,T,a,s,IRVol,ZCCurve);
			P = Garanti * exp(-ZCR *T);
			Nmin = 1./100. * V0;
			Etheo = MAX(0,MIN(AlphaTheo * (V-P),V + AddLeverage*V0));
			TimeToRoll = RollPeriodInit;
			Dur = ( 1 - exp(-(r+Implied) * (RefMat + TimeToRoll)))/(r + Implied);
			Soult = (Implied - Running) * Dur;
			SoultBid = Soult - BidMid;
			SoultOffer = Soult + BidMid;
			PrixBid = 1 - SoultOffer;
			PrixOffer = 1 - SoultBid;


			NbCouponPayed = 0;
			liquidate = false;
			// --------------------------------------------------------------
			// Mise en place du fond
			N = Etheo/PrixOffer;
			E = Etheo;
			f = V0 - E;
			f -= PrimeUpFront*Garanti;
			Alpha = E/(V-P);
			// --------------------------------------------------------------
			// Simulation des variables
			for (int k = 0; k<nbstep;k++)
			{
				vZimplied[k] = sqrt(-2*log(g05cac())) * cos(2*PI*g05cac());		
				vW[k] = sqrt(-2*log(g05cac())) * sin(2*PI*g05cac());
				vWidener[k] = sqrt(-2*log(g05cac())) * sin(2*PI*g05cac());

				vdY[k] = 0.;
				x = g05dac(0,1);

				barriere = exp(-Lambda*dt);
				while (x >= barriere)
				{
					x *= g05dac(0,1);
					vdY[k] ++;
				}
			}

			// Data trajectoire
			if (Traj)
			{ 
				itsZImplied = vZimplied;
				itsSaut = vdY;
			}
			// --------------------------------------------------------------
			// TRAJECTOIRES
			for (k = 0; k<nbstep;k++)
			{	
				ZCR = TauxVasicek(k * dt, T, a, s, IRVol,ZCCurve);
				P = Garanti * exp(-(ZCR) * (T - k * dt));
				ValoCoupon = 0.;
				// --------------------------------------------------------------
				// Market Moves
				if (liquidate == false)
				{	
					f +=(f+N)*(r-Funding)*dt+MIN(f,0)*(Financement+Funding)*dt-ManagementCost*Garanti*dt+Running*N*dt;

					logImpliedVasicek = log(ImpliedVasicek);
					logImpliedVasicek = logImpliedVasicek*(1 - Speed*dt) + MeanImplied*Speed*dt + sqrt(dt)*VolImplied*vZimplied[k];
					ImpliedVasicek = exp(logImpliedVasicek);
	
					Widener *= (1 + MuWidener*dt + SigmaWidener*sqrt(dt)*vWidener[k]);

					Implied = ImpliedVasicek + PartWidener * Widener;

					N *= (1+Jump*vdY[k]);
						
					while (vdY[k] > 0)
					{
						Jump *= 1./(1.+Jump);
						vdY[k] -- ;
					}
					
					Dur = (1 - exp(-(r + Implied) * (RefMat + TimeToRoll))) / (r + Implied);
					Soult = (Implied - Running) * Dur;
					E = N * (1-Soult);
					V = E + f;
					Alpha = E/(V-P);

					s = s - a*s*dt + IRVol*sqrt(dt)*vW[k];
					r = Forward((k+1)*dt,dt,ZCCurve) + s;
					TimeToRoll -= dt;
				}
				else
				{V *= exp(ZCRstop * dt);}
				// --------------------------------------------------------------
				// Swtich
				if ( (TimeToRoll<1e-8) && ( liquidate == false))
				{
					Jump = JumpInit;
					TimeToRoll = RollPeriod;
					Dur = (1 - exp(-(r + Implied) * (RefMat + TimeToRoll))) / (r + Implied);
					SoultNew = (Implied - Running) * Dur;
					SoultNewBid = SoultNew - BidMid;
					N *= (1-SoultOffer)/(1-SoultNewBid);
					Soult = SoultNew;
					E = N*(1-Soult);
					V = E + f;
					Widener = InitWidener;
				}
				// --------------------------------------------------------------
				// On calcule le Fwd à tous les pas de temps
				Eurib1Y = Forward((k+1)*dt,365*dt,ZCCurve) + s;
				Eonia = Forward((k+1)*dt,dt,ZCCurve) + s;
				ProductTaux *= 1.+Eonia;
				MeanTaux = pow(ProductTaux,1./(k+1)) - 1.0;
				// --------------------------------------------------------------
				// Paiement du coupon
				if (((k+1) % 365) == 0)		
				{	
					if ( ((V - P)/V0 >= FactorCoupon * (SpreadCoupon + FixedEurib*Eurib1Y)) && ((k+1) <= StopMat * 365) )
					{
						f -= Garanti * (SpreadCoupon + FixedEurib * Eurib1Y);
						TauxCoupon[NbCouponPayed] = Eurib1Y;
						NbCouponPayed ++;
					}
				}
				// --------------------------------------------------------------
				// Rebalancement
				if (liquidate == false)
				{
					index = riskFactorImplied->FindIndex(Implied);
				
					AlphaTheo = MidValue[index];
					Rebal = MAX(0,MIN(AlphaTheo * (V - P),V + AddLeverage * V0)) - E;

					SoultBid = Soult - BidMid;
					PrixOffer = 1 - SoultBid;
					NRebalOffer = Rebal / PrixOffer;

					SoultOffer = Soult + BidMid;
					PrixBid = 1 - SoultOffer;
					NRebalBid = Rebal / PrixBid;

					if ((Alpha < MinValue[index]) && (fabs(NRebalOffer) > Nmin))
                    {
						CountRebal ++;
						E += Rebal;
						N += NRebalOffer;
						f = V - E;
						Alpha = E/(V-P);
					}
               
                    if ((Alpha > MaxValue[index]) && (fabs(NRebalBid) > Nmin))
                    {
						CountRebal ++;
						E += Rebal;
						N += NRebalBid;
						f = V - E;
						Alpha = E/(V-P);
					}

					if (V < P + DCushion*Garanti)
                    {
						ZCRstop = ZCR;
						liquidate = true;
						Tdesactivation += k*dt;
						NbLoss ++;
					}
                }	// fin du if
				
				if (Traj)
				{
					itsImplied[k]	= Implied; 			itsE[k]			= E;
					itsf[k]			= f;				itsN[k]			= N;
					itsP[k]			= P;				itsV[k]			= V;
					itsAlphaMin[k]	= MinValue[index];	itsAlphaMax[k]	= MaxValue[index];
					itsAlpha[k]		= Alpha;			itsNbRebal[k]	= CountRebal;
					itsDur[k]		= Dur;				itsSoult[k]		= Soult;
					itsR[k]			= r;				itsS[k]			= s;
					itsTimeToRoll[k]= TimeToRoll;		itsJump[k]		= Jump;
					itsWidener[k]	= Widener;
				}

				// On stocke la valeur de V dans une matrice
				if ((Traj == false) && (((k+1) % 365) == 0))
				{
					int time = (k+1)/365;
					itsMeanV[time-1] +=V;
					// Le TRI est calculé en prenant en compte les coupons
					TRI = exp(log(Vtot/V0)/time) - 1;
					ValTri->SetValue(i,time-1,TRI);
				}

				// On verifie si on est passé en dessous de la barriere down apres la durée lockée
				//if (((int)((k+1)/365) >= TimeLocked) && (V < BarriereDown) && (Barriere == false) && (Traj == false))
				if (((int)((k+1)/365) >= TimeLocked) && (rdt < BarriereTaux) && (Barriere == false) && (Traj == false))
				{
					double tpsSupTimeLocked  = (k+1)/365. - TimeLocked;
					int nbtrimestreSupTimeLocked = (int)(tpsSupTimeLocked/0.25);
					int indiceBarDTime = nbtrimestreSupTimeLocked;
					itsBarriereDownTime[indiceBarDTime] ++;
					Barriere = true;
				}
				
				// Valeur du fond en tenant compte des coupons versés
				for (int indice = 0; indice<NbCouponPayed; indice ++)
				{
					C[indice] *= pow((1+Eonia),dt);
					ValoCoupon += C[indice];
				}
					//C +=  Garanti * (SpreadCoupon + FixedEurib * TauxCoupon[indice])*exp(TauxCoupon[indice]*(k+1-(indice+1)*365)*dt);

				Vtot = V+ValoCoupon;

				//  Calcule du rendement
				rdt = log(Vtot/V0)/((k+1)*dt)-MeanTaux;
			} // fin du for k // pas de temps

			if (Traj == false)
			{
				itsNbRebal[0] += CountRebal;
				if (CountRebal>itsNbRebal[1])
					itsNbRebal[1] = CountRebal;
				if (CountRebal<itsNbRebal[2])
					itsNbRebal[2] = CountRebal;
			}

			// Dans le cas ou l'on n'a pas été liquidé, il faut racheter la protection sur la tranche
			if (liquidate == false)
				V = f + N*(1-SoultOffer);

			VFinal += V;
			if (V < Garanti)
			{
				Loss += V - Garanti;
				V = Garanti;
			}
			if (liquidate == false)
				Tdesactivation += StopMat;

			itsNbCouponPayed[0] += NbCouponPayed;
			itsNbCouponPayed[NbCouponPayed+1] ++;
		}// fin for i (simul)

	if (NbLoss>0)
		Loss *= ZCCurve[ZCCurve.size() - 1]/NbLoss;	
	
	VFinal *= 1./nbSimul;
	Tdesactivation *= 1./nbSimul;
	itsNbCouponPayed[0] *= 1./nbSimul;

	((ICM_Cppi*)this->GetSecurity())->SetAvgLossValue(Loss);
	((ICM_Cppi*)this->GetSecurity())->SetDefTime(Tdesactivation);
	((ICM_Cppi*)this->GetSecurity())->SetNbLoss(NbLoss);
	((ICM_Cppi*)this->GetSecurity())->SetAvgFinalValue(VFinal);

	if (Traj == false)
	{
		// Pour recuperer la moyenne du nombre de rebalancement on divise par le nombre de simul
		itsNbRebal[0] *= 1.0/nbSimul;

		// On recupere la moyenne, la variance de la valeur V chaque année et la distribution
		for (int  z = 0; z<StopMat;z++)
		{
			// On trie la z ieme colonne de la matrice de valeur par ordre croissant
			ValTri->Sort(z);

			// Moyenne/Sigma
			vector<double> Col = ValTri->ColAsStdVector(z);
			vector<double> Vval; Vval.resize(nbSimul);

			for (int n = 0 ; n<nbSimul;n++)
				Vval[n] = V0 * exp((z+1)*log(Col[n]+1));

			itsMeanV[z] /= nbSimul;
			itsSigmaV[z] = Sigma(Vval,itsMeanV[z]);
			itsSigmaV[z] *= 1./itsMeanV[z];


			// Distribution
			int indDistrib = 0;
			for (int y = 0; y<nbSimul; y++)
			{
				while (Col[y] > 0.5/100 * indDistrib)
					indDistrib ++;

				(*itsDistribTRI)(indDistrib,z) ++;
			}
		}
	}

	double prix = Loss + ManagementCost*Garanti*Tdesactivation*ZCCurve[ZCCurve.size() - 1];

	// --------------------------
	// Delete
	// if (RecovTraxx)
	// 	delete[] RecovTraxx ;
	// RecovTraxx = NULL;

	if (riskFactorImplied)
		delete riskFactorImplied;
	riskFactorImplied = NULL;

	if (ValTri)
		delete ValTri;
	ValTri = NULL;
	// --------------------------

	return (prix);
}


double ICM_Pricer_Cppi::Accrued()
{
	ICMTHROW(ERR_INVALID_ARGUMENT,"Not implemented"); 
}
double ICM_Pricer_Cppi::FeeLegPV () 
{
	ICMTHROW(ERR_INVALID_ARGUMENT,"Not implemented"); 
}

double ICM_Pricer_Cppi::DefLegPV () 
{
	ICMTHROW(ERR_INVALID_ARGUMENT,"Not implemented"); 
}

double ICM_Pricer_Cppi::ComputeDuration(void) 
{
	ICMTHROW(ERR_INVALID_ARGUMENT,"Not implemented"); 
}

double ICM_Pricer_Cppi::ComputeSpread(const double& MtM ) 
{
	ICMTHROW(ERR_INVALID_ARGUMENT,"Not implemented"); 
}

double ICM_Pricer_Cppi::ComputeImpliedVol(const double& Price)
{
	ICMTHROW(ERR_INVALID_ARGUMENT,"Not implemented"); 
}

double ICM_Pricer_Cppi::Compute_Fwd_Spread(const ARM_Date &Mty, const ARM_Date &CDS_ExpiryDate, double& dur) 
{
	ICMTHROW(ERR_INVALID_ARGUMENT,"Not implemented"); 
}


// 	virtual 
double ICM_Pricer_Cppi::ComputeSensitivity(const qSENSITIVITY_TYPE& typesensi, 
		const std::string& plot,
		const std::string& label,
		double epsilon , double epsilonGamma) 
{
	return ICM_Pricer::ComputeSensitivity(typesensi,plot,label,epsilon, epsilonGamma); 
}

