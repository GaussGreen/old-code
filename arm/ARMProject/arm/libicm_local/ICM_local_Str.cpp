


#include "firstToBeIncluded.h"

#include <math.h>
#include <time.h>
#include <sys/timeb.h>
#include <ICMKernel\glob\icm_enums.h>
#include <ICMKernel\glob\icm_constants.h>

#include <ICMKernel\str\icm_str_models.h>
#include <ICMKernel\str\icm_util_str.h>
#include <ICMKernel\util\icm_utils.h>
#include <ARMKernel\crv\volflat.h>
#include <ARMKernel\crv\volint.h>

#include <ARM\libarm_local\ARM_local_class.h>
#include <ARM\libarm_local\ARM_local_glob.h>
#include <ARM\libarm_local\ARM_local_persistent.h>
#include <ARM\local_xlarm\ARM_local_interglob.h>
#include <ICMKernel\inst\icm_nthtd.h>
#include <ICMKernel\inst\icm_mez.h>
#include <ICMKernel\inst\icm_cppi.h>
#include <ICMKernel\inst\icm_pf.h>
#include <ICMKernel\inst\icm_cdo2.h>
#include <ICMKernel\pricer\icm_pricer_cppi.h>
#include "ICM_local_Str.h"
#include <ICMKernel\pricer\icm_pricer_homogeneous.h>
#include <ICMKernel\pricer\icm_pricer_homogeneous_smile.h>
#include <ICMKernel\pricer\icm_pricer_analytic_cdo2_smile.h>
#include <ICMKernel\mod\icm_Homogene_FastLossDistrib.h>
#include <ICMKernel\glob\icm_smile_correlation.h>
#include <ICMKernel\glob\icm_betas_correlation.h>
#include <ICMKernel\glob\icm_corrmatrix.h>
#include <ICMKernel\pricer\icm_pricer_adviser.h>
#include <ICMKernel\inst\icm_option.h>
#include "ICMKernel/pricer/icm_pricer_option.h"
#include "ICMKernel/glob/icm_mktdatamng.h"
#include "ICMKernel/inst/icm_credit_index.h"

//------------------------------
// Include utiles pour la macro
#include <sstream>
#include <iomanip>
#include <iostream>
#include <windows.h> 

using namespace std;
//------------------------------
#ifndef unix 
#define ITOA _itoa
#else
#define ITOA intToStr
#endif

static void FreePointerTabChar(char**& Tab,const int& size)
{
	if (Tab == NULL) return;

	for (int i=0; i<size; i++)
	{
		if (Tab[i])	delete[] Tab[i];
		Tab[i] = NULL;
	}

	if (Tab) delete[] Tab;

	Tab = NULL;
}

void MyDeleteTab(const std::vector<const ICM_DefaultCurve*>& items)
{
	unsigned int i=0; 
	for(i=0;i<items.size();i++) delete items[i] ; 
}

//------------------------------------------------------------

// ----------------------------------------------------------
// Fonction de test
// ----------------------------------------------------------
extern long ICMLOCAL_STR_Cppi (vector<double> CrbDefaut,
							   long nbRowCrbDef,
							   long nbColCrbDef,
							   double Recov,
							   vector<double> Correl,
							   long nbRowCorrel,
							   long nbColCorrel,
							   vector<double> DataCppi,
							   vector<double> rFactor,
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
							   ARM_result& result)
{
	ICMMSG(WARN,"Using ICMLOCAL_STR_Cppi" ); 
	CCString msg ("");
	int i = 0;
	double Basis = 365.;
	DWORD temps_debut = GetTickCount();
	DWORD temps_execution = 0.;
try
	{

	vector<double> Zc;
	Zc.resize(nbRowtaux);
	
	for (i = 0;i<nbRowtaux;i++)
		Zc[i] = taux[i*nbColtaux + 1];

	double Taux = -log(Zc[19])/5.;
	
	ARM_Date AsOfDate;
	AsOfDate.Today();
	
	// -------------------------------------------------------------------
	// CrbTaux
	ARM_ZeroFlat CrbTaux = ARM_ZeroFlat(AsOfDate, Taux);
	// -------------------------------------------------------------------
	// Courbes de défaut
	vector<double> YForDates;
	YForDates.resize(nbColCrbDef);
	vector<double> Spread;
	Spread.resize(nbColCrbDef);
	for (i=0;i<nbColCrbDef;i++)
	{
		YForDates[i] = CrbDefaut[i]; 
		Spread[i] = CrbDefaut[i + nbColCrbDef];
	}
	bool IsDate = false;
	
	if (YForDates[0] > 100)
		IsDate = true;

	ICM_DefCurvStr CrbTraxx (Spread,YForDates,Recov,&CrbTaux,"Traxx",IsDate);

	ICM_DefaultCurve** TabCrbDefaut = new ICM_DefaultCurve*[1];  
	TabCrbDefaut[0] = &CrbTraxx;
	// -------------------------------------------------------------------
	// Correlation
	int NbNomTraxx = 125;
	ARM_Vector MaturityTraxx = ARM_Vector(1,5);
	ARM_Vector Proportion = ARM_Vector(1,1);
	
	std::vector<std::string> LabelTraxx(1); 
	LabelTraxx[0]="Traxx";

	std::vector<std::string> LabelTraxxEur(NbNomTraxx); 
	for (i=0;i<NbNomTraxx;i++)
		LabelTraxxEur[i] = "Traxx";

	ARM_Vector StrikEUR = ARM_Vector(nbColCorrel,0.);
	for (i=0;i<nbColCorrel;i++)
		StrikEUR.Elt(i) = 100 * Correl[i];

	ARM_Vector CorrelEUR = ARM_Vector(nbColCorrel,0.);
	for (i=0;i<nbColCorrel;i++)
		CorrelEUR.Elt(i) = 100 * Correl[i + nbColCorrel];

	ARM_Matrix* CorrelationEur = new ARM_Matrix(MaturityTraxx.GetSize(),CorrelEUR.GetSize());
	for (i=0;i<CorrelEUR.GetSize();i++)
		CorrelationEur->Elt(0,i) = CorrelEUR.Elt(i);

	vector<const ARM_VolCurve*> Vector_VolCurve(1);
	Vector_VolCurve[0] = new ARM_VolLInterpol(AsOfDate,
											 (ARM_Vector*)MaturityTraxx.Clone(),
											 (ARM_Vector*)StrikEUR.Clone(),				// Strike en bp
											 (ARM_Matrix*)CorrelationEur->Clone(),		// Valeur de la correl en %
											 1,
	 										 K_SMILE_VOL);

	ARM_Vector Vsmilestrikelow(1,CorrelEUR.Elt(0));
	ARM_Vector Vsmilestrikehigh(1,CorrelEUR.Elt(1));

// FIXMEFRED: mig.vc8 (28/05/2007 14:20:19):cast
	ICM_Smile_Correlation Corr(AsOfDate,
									 "Traxx",
									 &(*Vector_VolCurve.begin()), 
									 LabelTraxx,
									 &Proportion,
									 &Vsmilestrikelow,
									 &Vsmilestrikehigh);
	// ---------------------------------------------------------------------------
	// Objet Mezz
	double MezzAmount = 0.03 * 125 * 1e7;
	double SubAmount = 0.;
	double FloatingPayerAmount = MezzAmount;
	int stubrule = 1;
	int FrequencyDL = 4;
	std::vector<std::string> Label(NbNomTraxx);
	// double * pIssuersNotionals = new double[NbNomTraxx];
	std::vector<double> pIssuersNotionals (NbNomTraxx); 
	for (i=0;i<NbNomTraxx;i++)
	{
		Label[i] = "Traxx";
		pIssuersNotionals[i] = 1e7;
	}
	ARM_Date EndDate = AsOfDate;
	double addMatu = 5.0;
//	double addMatu = 7.0;

	EndDate.AddYears(addMatu);

ICM_Mez*	mez = new ICM_Mez(AsOfDate,
					  EndDate,
					  (ARM_Date*)0,
					  (ARM_Date*)0,
					  0.01,			// FixedRate
					  K_ADJUSTED,	// intRule
					  K_ADJUSTED,	// adjStartDate
					  SubAmount,
					  MezzAmount,
					  Label,
					  pIssuersNotionals,
					  // NbNomTraxx,
					  K_QUARTERLY,			// Frequency
					  KACTUAL_360,
					  MezzAmount,
					  qACCRUED_SETTLED,
					  ARM_DEFAULT_COUNTRY,
					  FloatingPayerAmount,
					  stubrule,
					  0,
					  FrequencyDL,
					CREDIT_DEFAULT_VALUE, // const double& Binary /* = CREDIT_DEFAULT_VALUE*/,
					std::string(), // const std::string&  payCalName /* = NULL*/,
					qRunning_Leg, // const qCredit_Leg_Type& TypeFeeLeg/* =qRunning_Leg*/,
					qStandart_Recovery_Leg, // const qCredit_Leg_Type& TypeDefLeg/* =qStandart_Recovery_Leg*/,
					INCLUDE_MATURITY // const bool& IncludeMaturity /* = INCLUDE_MATURITY*/);
					  );
	// ---------------------------------------------------------------------------
	// Objet Risky Factor
	vector<double> inf;inf.resize(nbRowRiskyFactor);
	vector<double> max;max.resize(nbRowRiskyFactor);
	vector<double> valuemin;valuemin.resize(nbRowRiskyFactor);
	vector<double> valuemax;valuemax.resize(nbRowRiskyFactor);

	for (i=0;i<nbRowRiskyFactor;i++)
	{
		inf[i] = rFactor[i*nbColRiskyFactor];
		max[i] = rFactor[1 + i*nbColRiskyFactor];
		valuemin[i] = rFactor[2 + i*nbColRiskyFactor];
		valuemax[i] = rFactor[3 + i*nbColRiskyFactor];
	}

	ICM_RangeFactor* riskyFactor = new ICM_RangeFactor(inf,max,valuemin,valuemax);
	// ---------------------------------------------------------------------------
	// Objet CPPI
	double Maturity = DataCppi[0];
	double AsOfExcelDate = AsOfDate.GetJulian() - 2415019;
	ARM_Date CppiEndDate;
	CppiEndDate.Today();
	char* pEndDate =  new char[11];
	if (Maturity<100)
	{
		Local_XLDATE2ARMDATE(AsOfExcelDate + Basis*Maturity,pEndDate);
		CppiEndDate = (ARM_Date) pEndDate;
	}
	else
	{
		Local_XLDATE2ARMDATE(Maturity,pEndDate);
		CppiEndDate = (ARM_Date) pEndDate;
	}	

	double notional = DataCppi[1];
	double ProtectedNotional = DataCppi[2];
	double AdditionalLeverage = DataCppi[3];
	double DesactivateCushion = DataCppi[4];
	double ManagementCost = DataCppi[5];

	ICM_Cppi* cppi = new ICM_Cppi(AsOfDate,
						   	 CppiEndDate,
							 mez,
							 notional,
							 ProtectedNotional,
							 AdditionalLeverage,
							 riskyFactor,
							 DesactivateCushion,
							 ManagementCost,
							 ARM_DEFAULT_CURRENCY,
							 "Traxx");
	// ---------------------------------------------------------------------------
	// Market Data Manager
	ICM_MktDataMng datamkt ; // = ICM_MktDataMng();
	datamkt.insert((ARM_Object*)CrbTaux.Clone());
	datamkt.insert((ARM_Object*)CrbTraxx.Clone());
	datamkt.insert((ARM_Object*)Corr.Clone());
	// ---------------------------------------------------------------------------
	// Parameters
	ICM_Parameters parameters = ICM_Parameters();
	char* fund = "FUNDING";
	ARM_Vector funding =  ARM_Vector(1,DataPlacement[0]);
	parameters.Push((ARM_Vector*)funding.Clone(),fund);

	char* fin = "FINANCEMENT";
	ARM_Vector financement =  ARM_Vector(1,DataPlacement[1]);
	parameters.Push((ARM_Vector*)financement.Clone(),fin);

	char* run = "RUNNING";
	ARM_Vector running =  ARM_Vector(1,DataPlacement[2]);
	parameters.Push((ARM_Vector*)running.Clone(),run);

	char* lamb = "INTENSITE_POISSON";
	ARM_Vector lambda =  ARM_Vector(1,DataModel[0]);
	parameters.Push((ARM_Vector*)lambda.Clone(),lamb);
	
	char* vol = "VOL_IMPLIED";
	ARM_Vector vol_implied =  ARM_Vector(1,DataModel[1]);
	parameters.Push((ARM_Vector*)vol_implied.Clone(),vol);

	char* mean = "MEAN_IMPLIED";
	ARM_Vector mean_implied =  ARM_Vector(1,DataModel[2]);
	parameters.Push((ARM_Vector*)mean_implied.Clone(),mean);

	char* speed = "SPEED_REVERSION_IMPLIED";
	ARM_Vector SReversion_implied =  ARM_Vector(1,DataModel[3]);
	parameters.Push((ARM_Vector*)SReversion_implied.Clone(),speed);

	char* soultbid = "SOULT_BIDMID";
	ARM_Vector S_bid =  ARM_Vector(1,DataSoult[0]);
	parameters.Push((ARM_Vector*)S_bid.Clone(),soultbid);

	char* soultswitch = "SOULT_SWITCH";
	ARM_Vector S_Switch =  ARM_Vector(1,DataSoult[1]);
	parameters.Push((ARM_Vector*)S_Switch.Clone(),soultswitch);

	char* vasi_taux = "TAUX";
	ARM_Vector V_taux =  ARM_Vector(1,DataVasicek[0]);
	parameters.Push((ARM_Vector*)V_taux.Clone(),vasi_taux);

	char* vasi_param = "VASICEK_PARAM";
	ARM_Vector V_param =  ARM_Vector(1,DataVasicek[1]);
	parameters.Push((ARM_Vector*)V_param.Clone(),vasi_param);

	char* vasi_vol = "VASICEK_VOL";
	ARM_Vector V_vol =  ARM_Vector(1,DataVasicek[2]);
	parameters.Push((ARM_Vector*)V_vol.Clone(),vasi_vol);

	char* nbsim = "MC_NB_SIMUL";
	if (Fonction == "T")
	{	DataSimul[0] = 1; }

	ARM_Vector nbS =  ARM_Vector(1,DataSimul[0]);
	parameters.Push((ARM_Vector*)nbS.Clone(),nbsim);

	char* step = "MC_PAS";
	ARM_Vector Step =  ARM_Vector(1,DataSimul[1]);
	parameters.Push((ARM_Vector*)Step.Clone(),step);

	char* indice = "INDICE";
	ARM_Vector Ind =  ARM_Vector(1,DataSimul[2]);
	parameters.Push((ARM_Vector*)Ind.Clone(),indice);

	char* crbtaux = "CRB_TAUX";
	ARM_Vector ZC =  ARM_Vector(Zc);
	parameters.Push((ARM_Vector*)ZC.Clone(),crbtaux);

	char* PrimeUpFront = "PRIME_UPFRONT";
	ARM_Vector PUpFront =  ARM_Vector(1,DataCoupon[0]);
	parameters.Push((ARM_Vector*)PUpFront.Clone(),PrimeUpFront);
 
	char* SCoupon = "SPREAD_COUPON";
	ARM_Vector SpreadCoupon =  ARM_Vector(1,DataCoupon[1]);
	parameters.Push((ARM_Vector*)SpreadCoupon.Clone(),SCoupon);

	char* FCoupon = "FACTOR_COUPON";
	ARM_Vector FactorCoupon =  ARM_Vector(1,DataCoupon[2]);
	parameters.Push((ARM_Vector*)FactorCoupon.Clone(),FCoupon);

	char* SMat = "STOP_MAT";
	ARM_Vector StopMat =  ARM_Vector(1,DataCoupon[3]);
	parameters.Push((ARM_Vector*)StopMat.Clone(),SMat);

	char* FEurib = "FIXED_EURIB";
	ARM_Vector FixedEurib =  ARM_Vector(1,DataCoupon[4]);
	parameters.Push((ARM_Vector*)FixedEurib.Clone(),FEurib);

	char* InitWide = "INIT_WIDENER";
	ARM_Vector InitWidener  =  ARM_Vector(1,DataWidener[0]);
	parameters.Push((ARM_Vector*)InitWidener.Clone(),InitWide);
	
	char* MuWide = "MU_WIDENER";
	ARM_Vector MuWidener  =  ARM_Vector(1,DataWidener[1]);
	parameters.Push((ARM_Vector*)MuWidener.Clone(),MuWide);

	char* SigmaWide = "SIGMA_WIDENER";
	ARM_Vector SigmaWidener  =  ARM_Vector(1,DataWidener[2]);
	parameters.Push((ARM_Vector*)SigmaWidener.Clone(),SigmaWide);

	char* PartWide = "PART_WIDENER";
	ARM_Vector PartWidener  =  ARM_Vector(1,DataWidener[3]);
	parameters.Push((ARM_Vector*)PartWidener.Clone(),PartWide);

	char* BarriereDown = "BARRIERE_DOWN";
	ARM_Vector BarDown  =  ARM_Vector(1,DataBarriere[0]);
	parameters.Push((ARM_Vector*)BarDown.Clone(),BarriereDown);

	char* LockedTime = "LOCKED_TIME";
	ARM_Vector LTime  =  ARM_Vector(1,DataBarriere[1]);
	parameters.Push((ARM_Vector*)LTime.Clone(),LockedTime);

	// ---------------------------------------------------------------------------
	// Pricer CPPI
	// ICM_Pricer_Cppi pricerCppi = ICM_Pricer_Cppi(cppi,&datamkt,&parameters);
	ICM_Pricer_Cppi pricerCppi; pricerCppi.Set(cppi,&datamkt,parameters,AsOfDate);
	// ---------------------------------------------------------------------------
	// Pricing
	double prix = pricerCppi.ComputePrice(qCMPPRICE);

	temps_execution += GetTickCount() - temps_debut;

	int nbcolRes = 20;
	int nbstep = (int)((CppiEndDate.GetJulian() - AsOfDate.GetJulian())/365./DataSimul[1]);
	
Resultat.resize(nbcolRes*nbstep);

// On garde les resultats
Resultat[0] = prix;
Resultat[nbcolRes] = ((ICM_Cppi*)pricerCppi.GetSecurity())->GetAvgFinalValue();
Resultat[2*nbcolRes] = ((ICM_Cppi*)pricerCppi.GetSecurity())->GetNbLoss();
Resultat[3*nbcolRes] = ((ICM_Cppi*)pricerCppi.GetSecurity())->GetAvgLossValue();
Resultat[4*nbcolRes] = ((ICM_Cppi*)pricerCppi.GetSecurity())->GetDefTime();
Resultat[5*nbcolRes] = temps_execution;
Resultat[6*nbcolRes] = pricerCppi.GetNbRebal(0);
Resultat[7*nbcolRes] = pricerCppi.GetNbRebal(1);
Resultat[8*nbcolRes] = pricerCppi.GetNbRebal(2);

for (i = 0; i<pricerCppi.GetNbCouponPayed().size() ;i ++)
	Resultat[(9+i)*nbcolRes] = pricerCppi.GetNbCouponPayed(i);
	  
if (Fonction == "T")	
{
	for (int j = 0;j<nbstep;j++)
	{
		Resultat[1 + j*nbcolRes] = pricerCppi.GetZImplied(j);
		Resultat[2 + j*nbcolRes] = pricerCppi.GetImplied(j);
		Resultat[3 + j*nbcolRes] = pricerCppi.GetE(j);
		Resultat[4 + j*nbcolRes] = pricerCppi.Getf(j);
		Resultat[5 + j*nbcolRes] = pricerCppi.GetN(j);
		Resultat[6 + j*nbcolRes] = pricerCppi.GetV(j);
		Resultat[7 + j*nbcolRes] = pricerCppi.GetP(j);
		Resultat[8 + j*nbcolRes] = pricerCppi.GetAlphaMin(j);
		Resultat[9 + j*nbcolRes] = pricerCppi.GetAlphaMax(j);
		Resultat[10 + j*nbcolRes] = pricerCppi.GetAlpha(j);
		Resultat[11 + j*nbcolRes] = pricerCppi.GetNbRebal(j);
		Resultat[12 + j*nbcolRes] = pricerCppi.GetDur(j);
		Resultat[13 + j*nbcolRes] = pricerCppi.GetSoult(j);
		Resultat[14 + j*nbcolRes] = pricerCppi.GetR(j);
		Resultat[15 + j*nbcolRes] = pricerCppi.GetS(j);
		Resultat[16 + j*nbcolRes] = pricerCppi.GetTimeToRoll(j);
		Resultat[17 + j*nbcolRes] = pricerCppi.GetJump(j);
		Resultat[18 + j*nbcolRes] = pricerCppi.GetSaut(j);
		Resultat[19 + j*nbcolRes] = pricerCppi.GetWidener(j);
	}
}
else
{
	// On recupere le sigma et la moyenne chaque année
	for (int indSigma = 0; indSigma<pricerCppi.GetSigmaV().size(); indSigma ++)
	{
		Resultat[1 + indSigma*nbcolRes] = pricerCppi.GetMeanV(indSigma);
		Resultat[2 + indSigma*nbcolRes] = pricerCppi.GetSigmaV(indSigma);
	}

	// Distribution de TRI
	int ColDistrib = pricerCppi.GetColDistrib();
	int RowDistrib = pricerCppi.GetRowDistrib();

	for (int k = 0; k<ColDistrib;k++)
	{
		for (int j = 0;j<RowDistrib;j++)
			Resultat[3 + k + j*nbcolRes] = pricerCppi.GetDistribTRI(j,k);
	}

	// BarrierDownTime
	for (k = 0;k<pricerCppi.GetBarriereDownTime().size();k++)
		Resultat[3 + ColDistrib + k*nbcolRes] = pricerCppi.GetBarriereDownTime(k); 
}
	// ----------------------------------------------------------------
	// Delete
	MyDeleteTab(pEndDate);
	//MyDeleteTab(pIssuersNotionals);
	// MyDeleteTab(Label);
	/* MyDeleteTab(Vector_VolCurve,1); */ delete Vector_VolCurve[0];
	// MyDeleteTab(LabelTraxxEur);
	// MyDeleteTab(LabelTraxx);
	MyDeleteTab(TabCrbDefaut);
	MyDelete(cppi);
	MyDelete(riskyFactor);
	MyDelete(mez);
	// ----------------------------------------------------------------

	result.setDouble(temps_execution);				
	return ARM_OK;
}
	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}

	catch (...)
	{
		result.setMsg ("ARM_ERR: Cppi : unrecognized failure");
		return ARM_KO;
	}

};
// ----------------------------------------------------------------------
// Generateur Alea (Methode de Lehmer)
// ----------------------------------------------------------------------
/*extern long ICMLOCAL_STR_GenAlea (int Seed,
								  int NbTir,
								  vector<double>& Resultat,
								  ARM_result& result)
{
	ICMMSG(WARN,"Using ICMLOCAL_STR_GenAlea" ); 
	CCString msg ("");
try
{
  _timeb tstruct;
  _ftime( &tstruct );

	if (Seed == CREDIT_DEFAULT_VALUE)
			Seed = (int)(tstruct.millitm);
	// Resultat
	Resultat.resize(NbTir);
	for (int i =0;i<NbTir;i++)
		 Resultat[i] = GenerateurAlea(Seed);
	
	double retour = 1.;
	result.setDouble(retour);				
	return ARM_OK;
}
 	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}

	catch (...)
	{
		result.setMsg ("ARM_ERR: Default Curve : unrecognized failure");
		return ARM_KO;
	}

};
*/
// ----------------------------------------------------------------------
// FONCTIONS TRANCHE
// ----------------------------------------------------------------------

// ------------------------------------------------------------------------
// Spread Forward Curves
// ------------------------------------------------------------------------
extern long ICMLOCAL_STR_SpreadForwardCurves(vector<double> Spread,
											 vector<double> Plot,
											 double Recov,
											 double taux,
											 double TimeInit,
											 double Duree,
											 CCString C_Stripping,
											 ARM_result& result)
{
	ICMMSG(WARN,"Using ICMLOCAL_STR_SpreadForwardCurves " << CCSTringToSTLString(C_Stripping)); 
	CCString msg ("");
	int i =0;
	double	Taux = taux * 100.;
	try
{
	double Basis = 365.;
	// -----------------------------------------
	// Date
	ARM_Date AsOfDate;
	AsOfDate.Today();
	// -----------------------------------------
	// Plots
	int sizePlot = Plot.size();
	ARM_Vector YF_plot = ARM_Vector(Plot);
	if (Plot[0]>100)
	{
		// On met la date du jour au format Excel
		double ExcelAsOf = AsOfDate.GetJulian() - 2415019;
		for (i=0;i<sizePlot;i++)
			YF_plot.Elt(i) = (Plot[i] - ExcelAsOf)/Basis;
	}
	vector<string> Term;
	for (i=0;i<sizePlot;i++)
	{
		char * t = NULL;
		ITOA(YF_plot.Elt(i),t,10);
		string s = string(t) + string("Y");
		Term.push_back(s);
		if (t) delete t;
	}
	// -----------------------------------------
	// Spreads
	int sizespread = Spread.size();
	if (sizespread == 1)
	{
		for (i = 1;i<sizePlot;i++)
			Spread.push_back(Spread[0]);
	}
	ARM_Vector Sp = ARM_Vector(Spread);
	// -----------------------------------------
	// Courbes
	ARM_ZeroFlat* CrbTaux = new ARM_ZeroFlat(AsOfDate, Taux);
	ICM_DefaultCurve* DCurve = NULL;
	// -----------------------------------------
	// Stripping
	   string Stripping = C_Stripping.c_str();
	   bool FastStrip = FastStripping(Stripping);
	// -----------------------------------------
	if (FastStrip)
	{	DCurve = new ICM_DefCurvStr(&Sp,&YF_plot,Recov,(ARM_ZeroCurve*)CrbTaux,"CrbFwd");}
	else
	{	DCurve = new ICM_Constant_Piecewise(AsOfDate,
											Term,
											(ARM_Vector*)Sp.Clone(),
											Recov,
											(ARM_ZeroCurve*)CrbTaux,
												K_ADJUSTED,	// intRule
												K_ADJUSTED,	// adjStartDate
											qCredit_Adjust20,
											ARM_DEFAULT_COUNTRY,
											"CrbFwd",
											false,
											NULL,
											//2 NULL,
											K_QUARTERLY,qDEFCURVE_DICHO,"STD", 
											ARM_Currency(ARM_DEFAULT_COUNTRY).GetCreditStartDateLag(),
											ICM_Parameters()
											);
	}
	// -----------------------------------------
	// Resultat
	double SForward = SpreadForward(DCurve,TimeInit,Duree);

	result.setDouble(SForward);
	// -----------------------------------------
	// DELETE
	MyDelete(DCurve);
	MyDelete(CrbTaux);
		
	return ARM_OK;
}
	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}
	catch (...)
	{
		result.setMsg ("ARM_ERR: Spread Forward Curves : unrecognized failure");
		return ARM_KO;
	}
};

// ----------------------------------------------------------
// Vecteur de Discount Price
// ----------------------------------------------------------
extern long ICMLOCAL_STR_DPrice (double YF_Maturity,
								 double taux,
								 ARM_result& result)
{
	ICMMSG(WARN,"Using ICMLOCAL_STR_DPrice" ); 

	CCString msg ("");
try
{
	double Taux = taux * 100.;
	ARM_Date AsOfDate;
	AsOfDate.Today();
	ARM_ZeroFlat* CrbTaux = new ARM_ZeroFlat(AsOfDate, Taux);

	double Maturity = YF_Maturity;
	MaturityYF(Maturity);
	
	double DP = CrbTaux->DiscountPrice(Maturity);
	
	result.setDouble(DP);

	// DELETE
	MyDelete(CrbTaux);

	return ARM_OK;
}
	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}
	catch (...)
	{
		result.setMsg ("ARM_ERR: Discount Price : unrecognized failure");
		return ARM_KO;
	}

};

// ---------------------------------------------------------------------------------------------------
// Vecteur de Proba de defaut (1 nom à tous les pas) : P_defaut calculée avec Courbe de défaut entière
// ---------------------------------------------------------------------------------------------------
extern long ICMLOCAL_STR_PDefaultCurves (double YF_Maturity,
										 double IssuersRecovery,
										 vector<double> Spread,
										 vector<double> YF_plot,
										 double taux,
										 CCString C_Stripping,
										 ARM_result& result)
{
	ICMMSG(WARN,"Using ICMLOCAL_STR_PDefaultCurves " << C_Stripping); 
	double Basis = 365.;
	CCString msg ("");
	int i = 0;
	double Taux = taux * 100.;
try
{
	ARM_Date AsOfDate;
	AsOfDate.Today();
	double AsOfJulian = AsOfDate.GetJulian();
	double ExcelAsOf = AsOfJulian - 2415019;

	double Maturity = YF_Maturity;
	MaturityYF(Maturity,AsOfJulian);

	ARM_ZeroFlat CrbTaux = ARM_ZeroFlat(AsOfDate, Taux);

	int sizeS = Spread.size();
	int sizeP = YF_plot.size();

	if (sizeS == 1)
	{
		for (i=1;i<sizeP;i++)
			Spread[i] = Spread[0];
	}
	
	ARM_Vector Sp = ARM_Vector(Spread);
	ARM_Vector Plot = ARM_Vector(YF_plot);
	if (YF_plot[0]>100.)
	{
		for (i=0;i<sizeP;i++)
			Plot.Elt(i) = (YF_plot[i] - ExcelAsOf)/Basis;
	}

	const int sizePlot = Plot.GetSize();

	vector<string> Term;
	for (i=0;i<sizePlot;i++)
	{
		char * t = NULL;
		ITOA(YF_plot[i],t,10);
		string s = string(t) + string("Y");
		Term.push_back(s);
		if (t) delete t;
	}
	// -----------------------------------------
	// Stripping
	   string Stripping = C_Stripping.c_str();
	   bool FastStrip = FastStripping(Stripping);
	// -----------------------------------------
	ICM_DefaultCurve* DCurve = NULL;
	if (FastStrip)
	{ DCurve = new ICM_DefCurvStr(&Sp,&Plot,IssuersRecovery,(ARM_ZeroCurve*)(&CrbTaux),"CrbDefaut");}
	else
	{ DCurve = new ICM_Constant_Piecewise(AsOfDate,
										  Term,
										  (ARM_Vector*)Sp.Clone(),
										  IssuersRecovery,
										  &CrbTaux,
												K_ADJUSTED,	// intRule
												K_ADJUSTED,	// adjStartDate
										  qCredit_Adjust20,
										  ARM_DEFAULT_COUNTRY,
										  "CrbDefaut",
										  false,
										  NULL,//2 NULL,
										  K_QUARTERLY,qDEFCURVE_DICHO,"STD", 
										  ARM_Currency(ARM_DEFAULT_COUNTRY).GetCreditStartDateLag(),
										  ICM_Parameters());
	}

	double P_defaut =  DCurve->DefaultProba(Maturity);
	result.setDouble(P_defaut);
	

	return ARM_OK;
}
	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}
	catch (...)
	{
		result.setMsg ("ARM_ERR: Proba de défaut : unrecognized failure");
		return ARM_KO;
	}
};
// ----------------------------------------------------------
// Proba de défaut d'une tranche ZC
// ----------------------------------------------------------
/** 
extern long ICMLOCAL_STR_ProbDefTranche (double YF_Maturity,
										 vector<double> Attachement,
										 double dNb_Issuers,
										 vector<double> Recovery,
										 vector<double> MatIssuersSpread,
										 vector<double> YF_Plot,		  
										 vector<double> VBeta,
										 double taux,
										 CCString C_Stripping,
										 ARM_result& result)
{
	ICMMSG(WARN,"Using ICMLOCAL_STR_ProbDefTranche " << CCSTringToSTLString(C_Stripping)); 
	CCString msg ("");

	int i,j = 0;
	int NbIssuers = (int)dNb_Issuers;
	double Taux = 100.* taux;
try
{
	double Basis = 365.;
	ARM_Date AsOfDate;
	AsOfDate.Today();
	double ExcelAsOf = AsOfDate.GetJulian() - 2415019;

	double Maturity = YF_Maturity;
	MaturityYF(Maturity);

	// Recovery
	double CalibRecovery = Recovery[0];
	double LossRecovery = 0.;
	if (Recovery.size() == 2)
		LossRecovery = Recovery[1];
	else
		LossRecovery = CalibRecovery;

	// Spread Issuers
	int NbRow = NbIssuers;
	int NbCol = YF_Plot.size();
	ICM_QMatrix<double>* IssuersSpreads = new ICM_QMatrix<double>(NbRow,NbCol,CREDIT_DEFAULT_VALUE);
	if (MatIssuersSpread.size() != NbCol*NbRow)
	{
		for (i=0;i<NbRow;i++)
			for (j=0;j<NbCol;j++)
				IssuersSpreads->SetValue(i,j,MatIssuersSpread[j]);
	}
	else
	{
		for (i=0;i<NbRow;i++)
			for (j=0;j<NbCol;j++)
				IssuersSpreads->SetValue(i,j,MatIssuersSpread[i*NbCol+j]);
	}

	// Beta
	int SizeBeta = VBeta.size();
	vector<double> Betas;
	Betas.resize(NbIssuers);
	if ( (SizeBeta<NbIssuers) && (SizeBeta == 1))
	{
		for (i=0;i<NbIssuers;i++)
			Betas[i] = VBeta[0];
	}
	else if (SizeBeta == NbIssuers)
	{
		for (i=0;i<NbIssuers;i++)
			Betas[i] = VBeta[i];
	}
	else 
	{
		MyDelete(IssuersSpreads);
			
		result.setMsg ("ARM_ERR: Vector Betas : Size Error");
		return ARM_KO;
	}	
	
	double LossUnit = (1-LossRecovery)*1e7;				// Pour ne pas avoir de probleme avec l'algo de calcul de loss unit
	double TdownValue = Attachement[0] * NbIssuers*1e7;
	double TupValue = Attachement[1] * NbIssuers*1e7;
// ----------------------------------------------------------------------
	ARM_ZeroFlat* CrbTaux = new ARM_ZeroFlat(AsOfDate, Taux);
	ARM_Vector Plot = ARM_Vector(YF_Plot);
	if (YF_Plot[0]>100)
	{
		for (i=0;i<NbCol;i++)
			Plot.Elt(i) = (YF_Plot[i] - ExcelAsOf)/Basis;
	}
	
	const int sizePlot = Plot.GetSize();

	vector<string> Term;
	for (i=0;i<sizePlot;i++)
	{
		char * t = NULL;
		ITOA(Plot.Elt(i),t,10);
		string s = string(t) + string("Y");
		Term.push_back(s);
		if (t) delete t;
	}

	
	// -----------------------------------------
	// Stripping
	   string Stripping = C_Stripping.c_str();
	   bool FastStrip = FastStripping(Stripping);
	// -----------------------------------------
	char* Label = new char[20];
	ICM_DefaultCurve** DefaultCurv = NULL;
	if (FastStrip)
	{
		for (i = 0;i<NbIssuers;i++)
		{
			sprintf(Label,"ISSUER_%i",i);
			ARM_Vector* Spread = IssuersSpreads->RowAsVector(i); 
			DefaultCurv[i] = new ICM_DefCurvStr(Spread,&Plot,CalibRecovery,CrbTaux,Label);

			MyDelete(Spread);
		}
	}
	else
	{
		for (i = 0;i<NbIssuers;i++)
		{
			sprintf(Label,"ISSUER_%i",i);
			ARM_Vector* Spread = IssuersSpreads->RowAsVector(i); 
			
			DefaultCurv[i] = new ICM_Constant_Piecewise(AsOfDate,
												        Term,
														(ARM_Vector*)Spread->Clone(),
														CalibRecovery,
					                                    CrbTaux,
												K_ADJUSTED,	// intRule
												K_ADJUSTED,	// adjStartDate
												qCredit_Adjust20,
												ARM_DEFAULT_COUNTRY,
												Label,
												false,
												NULL,//2 NULL,
												K_QUARTERLY,qDEFCURVE_DICHO,"STD", 
												ARM_Currency(ARM_DEFAULT_COUNTRY).GetCreditStartDateLag(),
												ICM_Parameters());

			MyDelete(Spread);
		}

	}

	double* P_defaut = new double[NbIssuers];
	memset(P_defaut,'\0',sizeof(double)*NbIssuers);
	
	for (int j=0;j<NbIssuers;j++)
		P_defaut[j] =  (DefaultCurv[j])->DefaultProba(Maturity);

	ICM_FFTDISTRIB* Distrib = new ICM_FFTDISTRIB(NbIssuers,P_defaut,Betas);

	double ExpectedLoss = Distrib->compute_expectedlosstranche(TupValue,TdownValue,LossUnit);

	double ProbDefaut = ExpectedLoss/(TupValue-TdownValue);

	result.setDouble(ProbDefaut);

	// --------------------------
	// Delete
	MyDelete(Distrib);
	MyDeleteTab(P_defaut);
	MyDelete(DefaultCurv,NbIssuers);
	MyDelete(Label);
	MyDelete(CrbTaux);
	MyDelete(IssuersSpreads);
	
	return ARM_OK;
}
	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}
	catch (...)
	{
		result.setMsg ("ARM_ERR: Proba defaut tranche : unrecognized failure");
		return ARM_KO;
	}
};
**/ 
// ------------------------  For Curves --------------------------------------
/** 
extern long ICMLOCAL_STR_DensityCurves (double YF_Maturity,
										double dNb_Issuers,
										double IssuersRecovery,
		  								vector<double> VIssuersSpreads,
										vector<double> YF_Plot,
										char** Label,
										vector<double> VBetas,
										double taux,
										vector<double>& Density,
										ARM_result& result)
{
	ICMMSG(WARN,"Using ICMLOCAL_STR_DensityCurves" ); 
	CCString msg ("");
	int i = 0;
	double Taux = taux * 100.;
try
{
	double Basis = 365.;
	ARM_Date AsOfDate;
	AsOfDate.Today();
	ARM_ZeroFlat* CrbTaux = new ARM_ZeroFlat(AsOfDate, Taux);
	double AsOfExcelDate = AsOfDate.GetJulian() - 2415019;
	double Maturity = YF_Maturity;
	if (YF_Maturity>100)
		Maturity = ( YF_Maturity - AsOfExcelDate)/Basis;

	int nbissuers = (int)dNb_Issuers;
	int sizePlot = YF_Plot.size();
	int MatSpreadSize = nbissuers*sizePlot;
	vector<double> VIssuersSpread;
	VIssuersSpread.resize(MatSpreadSize);
	if (VIssuersSpreads.size() == 1)
	{
		for (int i=0;i<nbissuers;i++)
		{
			for (int j=0;j<sizePlot;j++)
				VIssuersSpread[i*sizePlot + j] = VIssuersSpreads[0];
		}
	}
	if (VIssuersSpreads.size() == sizePlot)
	{
		for (int i=0;i<nbissuers;i++)
		{
			for (int j=0;j<sizePlot;j++)
				VIssuersSpread[i*sizePlot + j] = VIssuersSpreads[j];
		}
	}
	if (VIssuersSpreads.size() == MatSpreadSize)
	{
		for (int i=0;i<MatSpreadSize;i++)
			VIssuersSpread[i] = VIssuersSpreads[i];
	}
	
	vector<double> Betas;
	Betas.resize(nbissuers);
	if (VBetas.size() < nbissuers)
	{	
		for (int j=0;j<nbissuers;j++)
			Betas[j] = VBetas[0];
	}
	else
	{
		for (int j=0;j<nbissuers;j++)
			Betas[j] = VBetas[j];
	}
	ARM_Vector Plot = ARM_Vector(YF_Plot);
	if (YF_Plot[0]>100)
	{
		for (i=0;i<YF_Plot.size();i++)
			Plot.Elt(i) = (YF_Plot[i] - AsOfExcelDate)/Basis;
	}

	double* P_defaut = new double[nbissuers];
	memset(P_defaut,'\0',sizeof(double)*nbissuers);
	// On construit les courbes de défaut de chaque nom et on calcule la proba de défaut
	for (int i=0;i<nbissuers;i++)
	{
		ARM_Vector Sp = ARM_Vector(sizePlot,0.);
		for (int j=0;j<sizePlot;j++)
			Sp.Elt(j) = VIssuersSpread[i*sizePlot + j];

		ICM_DefCurvStr* DefCurve = new ICM_DefCurvStr(&Sp,&Plot,IssuersRecovery,(ARM_ZeroCurve*)CrbTaux,Label[i]);

		P_defaut[i] = 1 - DefCurve->SurvivalFunction(Maturity);

		MyDelete(DefCurve);
	}

	ICM_FFTDISTRIB* Distrib = new ICM_FFTDISTRIB(nbissuers, P_defaut, Betas);

	Density.resize(nbissuers+1);
	for (i=0;i<=nbissuers;i++)
		Density[i] = Distrib->getDensity()[i];

	// DELETE
	MyDeleteTab(P_defaut);
	MyDelete(Distrib);

	return ARM_OK;
}
 
	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}

	catch (...)
	{
		result.setMsg ("ARM_ERR: Density Curve : unrecognized failure");
		return ARM_KO;
	}

};
**/ 
// --------------------------------------------------------------------------------------------------------------
// CDO 2 SMILE
// --------------------------------------------------------------------------------------------------------------
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
								   double valShift,
								   vector<CCString> NumMethod,
								   CCString C_Function,
								   ARM_result& result)
{
ICMMSG(WARN,"Using ICMLOCAL_STR_CdO2Smile" ); 
	// DataMez[0] = Maturity;
	// DataMez[1] = K down CDO2;
	// DataMez[2] = K Up CDO2;
	// DataMez[3] = Pas paiement;
	// DataMez[4] = Coupon;
	// DataMez[5] = Nb Nom total;
	// DataMez[6] = Nb tranche fille;
	// DataMez[7] = taux en %

	CCString msg ("");



try
{
	int i = 0,j = 0;
	double temp = 0.;
	ARM_Date AsOf;
	AsOf.Today();
double Basis = 365.;
	// -------- Messages d'erreur --------------------
	if (!((C_Function == "NPV")||(C_Function == "BE")||(C_Function == "DL")||(C_Function == "FL")||(C_Function == "DUR")||(C_Function == "SENSI")||(C_Function == "ALL")||(C_Function == "STRIKE")))
	{
		result.setMsg ("ARM_ERR: Choix de Fonction Invalide");
		return ARM_KO;
	}
	if (DataMez.size() != 8)
	{
		result.setMsg ("ARM_ERR: Le champ DataMez est mal renséigné");
		return ARM_KO;
	}
	// ------------------------------------------------
	// Recuperation des données présentes dans DataMez
	double Maturity			= DataMez[0];
	double StrikeDownCDO2	= DataMez[1];
	double StrikeUpCDO2		= DataMez[2];
	double Pas				= DataMez[3];
	double Coupon			= DataMez[4];
	int NbName				= (int) DataMez[5];
	int NbTranche			= (int) DataMez[6];
	double Taux				= DataMez[7] * 100;			// DataMez[7] est en %; on doit passer des bp dans les instruments ARM
	// -----------------------------------------------
	// Dates : On récupère des données au format Excel date ou au format year fraction
		double AsOfExcelDate = AsOf.GetJulian() - 2415019;
		if (DataMez[0]>100)
			Maturity = (DataMez[0] - AsOfExcelDate)/Basis;

		ARM_Date EndDate;
		char* pEndDate =  new char[11];
		if (Maturity<100)
		{
			Local_XLDATE2ARMDATE(AsOfExcelDate + Basis*Maturity,pEndDate);
			EndDate = (ARM_Date) pEndDate;
		}
		else
		{
			Local_XLDATE2ARMDATE(Maturity,pEndDate);
			EndDate = (ARM_Date) pEndDate;
		}
	// ------------------------------------------------
	// Notionel Tranche
		vector<double> NotTranche;
		NotTranche.resize(NbTranche);
		if (NotionelTranche.size() == 0)
		{
			result.setMsg ("ARM_ERR: Le champ NominalTranche est mal renséigné");
			return ARM_KO;
		}
		else if ( (NotionelTranche.size() < NbTranche) && (NotionelTranche.size() > 0))
		{
			for (i=0;i<NbTranche;i++)
				NotTranche[i] = NotionelTranche[0];
		}
		else if (NotionelTranche.size() >= NbTranche)
		{
			for (i=0;i<NbTranche;i++)
				NotTranche[i] = NotionelTranche[i];
		}
	// ------------------------------------------------
	// Strike CDO
		ICM_QMatrix<double>* TdownTup = new ICM_QMatrix<double>(2,NbTranche);
		if (nbColStrikeCdo == 0)
		{
			result.setMsg ("ARM_ERR: Matrice de strike des tranches filles invalide");
			return ARM_KO;
		}
		else if ((nbColStrikeCdo>0) && (nbColStrikeCdo<NbTranche))
		{
			for (i=0;i<NbTranche;i++)
			{
				(*TdownTup)(0,i) = atof((char*)StrikeCdo[0]);					// On recupere le strike Down
				(*TdownTup)(1,i) = atof((char*)StrikeCdo[0+nbColStrikeCdo]);	// On recupere le strike Up
			}
		}
		else if (nbColStrikeCdo>=NbTranche)
		{
			for (i=0;i<NbTranche;i++)
			{
				(*TdownTup)(0,i) = atof((char*)StrikeCdo[i]);					// On recupere le strike Down
				(*TdownTup)(1,i) = atof((char*)StrikeCdo[i+nbColStrikeCdo]);	// On recupere le strike Up
			}
		}
	// ------------------------------------------------
	// DataSpread (la premiere ligne correspond aux plots, les autres lignes correspondent aux courbes de défaut de chaque nom
	if ((DataDefaultCurve.size() == 0) || (nbRowDefaultCurve <2))
	{
		result.setMsg ("ARM_ERR: DataDefaultCurve invalide");
		return ARM_KO;
	}
		
		// On recupere les plots
		double Data_i = 0.;
		ARM_Vector PlotNum = ARM_Vector(nbColDefaultCurve,1);
		for (i=0;i<nbColDefaultCurve;i++)
		{
			Data_i = atof((char*)DataDefaultCurve[i]);
			if (Data_i>100)
				PlotNum.Elt(i) = (Data_i - AsOfExcelDate)/Basis;
			else
				PlotNum.Elt(i) = Data_i;
		}


		const int sizePlot = PlotNum.GetSize();

		vector<string> Term;
		for (i=0;i<sizePlot;i++)
		{
			char * t = NULL;
			ITOA(PlotNum[i],t,10);
			string s = string(t) + string("Y");
			Term.push_back(s);
			if (t) delete t;
		}

		// On recupere les spreads
		ICM_QMatrix<double>* NameSpreads = new ICM_QMatrix<double>(NbName,nbColDefaultCurve,CREDIT_DEFAULT_VALUE);
		for (i=0;i<NbName;i++)
		{
			if (nbRowDefaultCurve < NbName)
			{
				for (j=0;j<nbColDefaultCurve;j++)
					(*NameSpreads)(i,j) = atof((char*)DataDefaultCurve[j + nbColDefaultCurve]);
			}
			if (nbRowDefaultCurve >= NbName)
			{
				for (j=0;j<nbColDefaultCurve;j++)
					(*NameSpreads)(i,j) = atof((char*)DataDefaultCurve[j + (i+1)*nbColDefaultCurve]);
			}
		}
	// ---------------------------------------------------------------------------
	// Label
	if (IssuersLabel.size() < NbName)
	{
		result.setMsg ("ARM_ERR: Le vecteur de Label n'est pas de la bonne taille");
		return ARM_KO;
	}
	
	char** Label = new char*[NbName];
	for (i =0; i<NbName;i++)
		Label[i] = IssuersLabel[i];
	// ---------------------------------------------------------------------------
	// Spread Traxx
	ARM_Vector SpreadTraxxEur = ARM_Vector(nbColDefaultCurve,1.);
	ARM_Vector SpreadTraxxUS =  ARM_Vector(nbColDefaultCurve,1.);
	
	if(SpreadTraxx.size() == 1)
	{
		for (i=0;i<nbColDefaultCurve;i++)
		{
			SpreadTraxxEur.Elt(i) = SpreadTraxx[0];
			SpreadTraxxUS.Elt(i)  = SpreadTraxx[0];
		}
	}
	else if(SpreadTraxx.size() == 2)
	{
		
		for (i=0;i<nbColDefaultCurve;i++)
			SpreadTraxxEur.Elt(i) = SpreadTraxx[0];
		
		for (i=0;i<nbColDefaultCurve;i++)
			SpreadTraxxUS.Elt(i) = SpreadTraxx[1];
	}
	else if(SpreadTraxx.size() == nbColDefaultCurve)
	{
		for (i=0;i<nbColDefaultCurve;i++)
		{
			SpreadTraxxEur.Elt(i) = SpreadTraxx[i];
			SpreadTraxxUS.Elt(i) = SpreadTraxx[i];	 
		}
	}
	else if(SpreadTraxx.size() == 2*nbColDefaultCurve) 
	{	
		for (i=0;i<nbColDefaultCurve;i++)
		{
			SpreadTraxxEur.Elt(i) = SpreadTraxx[i];
			SpreadTraxxUS.Elt(i) = SpreadTraxx[nbColDefaultCurve+i];	 
		}
	}
	else
	{
		result.setMsg ("ARM_ERR: Vector Spread Traxx : Size Error");
		return ARM_KO;
	}
	// ------------------------------------------------
	// Matrice des overlaps
		ICM_QMatrix<int>* Overlap = new ICM_QMatrix<int>(NbName,NbTranche);
		if ( (nbColMatriceOverlap < NbTranche) || (nbRowMatriceOverlap < NbName))
		{
			result.setMsg ("ARM_ERR: Matrice Overlap : Size Error");
			return ARM_KO;
		}

		// On remplit la matrice
		for (i=0;i<NbName;i++)
		{
			// On recupere le nom j de la tranche i
			for (j=0;j<NbTranche;j++)
				(*Overlap)(i,j) = atoi((char*)MatriceOverlap[j + nbColMatriceOverlap * i]);		
		}
	// ---------------------------------------------------------------------------
	// DataTraxx
		double MaturiteTraxx = 5.;
		// Le premier element de ce vecteur correspond à la recov utilisée pour le Traxx
		if (DataTraxx.size() <2)
		{
			result.setMsg ("ARM_Err : DataTraxx : Size Error");
			return ARM_KO;
		}
		// Dans ce cas, l'utilisateur a renseigné la recov Traxx et une seule proportion
		if (DataTraxx.size() < NbTranche + 1 )
		{
			// On verifie que la proportion est comprise entre 0 et 1
			if ( (DataTraxx[1] < 0) || (DataTraxx[1] > 1))
			{
				result.setMsg ("ARM_Err : DataTraxx : Porportion < 0 ou > 100% !");
				return ARM_KO;
			}	
			for (i=0;i<NbTranche - 1;i++)
				DataTraxx.push_back(DataTraxx[1]);
		}

		double RecovTraxx = DataTraxx[0];
	// ---------------------------------------------------------------------------
	// Recov
		double* RecovCalib = new double[NbName + 2];
		memset(RecovCalib,'\0',sizeof(double) * (NbName+2));
		// double* RecovLoss = new double[NbName + 2];
		// memset(RecovLoss,'\0',sizeof(double) * (NbName + 2) );
		ARM_Vector RecovLoss(NbName+2); 
		
		if ( (nbColRecov == 0) || (nbColRecov > 2))
		{
			result.setMsg ("ARM_ERR: Recovery : Size Error");
			return ARM_KO;
		}
		if (nbRowRecov < NbName)
		{
			for (i = 0;i<NbName;i++)
				RecovCalib[i] = atof((char*)Recov[0]);
		}
		else
			for (i = 0;i<NbName;i++)
				RecovCalib[i] = atof((char*)Recov[i*nbColRecov]);

		if (nbColRecov == 1)
		{
			for (i=0;i<NbName;i++)
				RecovLoss[i] = RecovCalib[i];
		}
		else 
		{
			if (nbRowRecov < NbName)
			{
				for (i = 0;i<NbName;i++)
					RecovLoss[i] = atof((char*)Recov[1]);
			}
			else
				for (i = 0;i<NbName;i++)
					RecovLoss[i] = atof((char*)Recov[i*nbColRecov + 1]);
		}

		RecovCalib[NbName] = RecovTraxx;
		RecovCalib[NbName+1] = RecovTraxx;
		
		RecovLoss[NbName] = RecovTraxx;
		RecovLoss[NbName+1] = RecovTraxx;

	// ------------------------------------------------
	// Forced Beta : On construit une matrice de correl par strike mais on force la correl aux niveaux souhaitse
	bool TestForcedCorrelCDO2 = false;
	bool TestForcedCorrelCDO = false;
	if ( (nbRowForcedCorrel> 1) && (nbRowForcedCorrel < NbTranche + 1))
	{
		result.setMsg("ARM_ERR : ForcedCorrel: Size Error");
		return ARM_KO;
	}
	
	// On recupere les infos contenues dans le vecteur de CCString et on les place dans un vecteur de double
	vector<double> CorrelForced;
	CorrelForced.resize(ForcedCorrel.size());
	for (i=0;i<ForcedCorrel.size();i++)
		CorrelForced[i] = atof((char*)ForcedCorrel[i]);

	if (CorrelForced[0] > 0)
		TestForcedCorrelCDO2 = true;

	// L'utilisateur a entré des betas forcés pour les tranches filles
	vector<double> ForcedCorrelCDODown;
	ForcedCorrelCDODown.resize(NbTranche);
	vector<double> ForcedCorrelCDOUp;
	ForcedCorrelCDOUp.resize(NbTranche);
	if (nbRowForcedCorrel >= NbTranche + 1)  
	{
		if (CorrelForced[2] > 0)
		{
			TestForcedCorrelCDO = true;
			
			if (nbColForcedCorrel >= 2)			
			{
				for (i=0;i<NbTranche;i++)
				{
					ForcedCorrelCDODown[i]	= CorrelForced[2*(i+1)];
					ForcedCorrelCDOUp[i]	= CorrelForced[2*(i+1) + 1];
				}
			}
			if (nbColForcedCorrel == 1)
			{
				for (i=0;i<NbTranche;i++)
				{
					ForcedCorrelCDODown[i]	= CorrelForced[i+1];
					ForcedCorrelCDOUp[i]	= 	ForcedCorrelCDODown[i];
				}
			}

		}
	}
	// On construit une ICM_Matrix que l'on utilisera pour construire l'objet portfolio
	ICM_Parameters* matrice = new ICM_Parameters();
	char* traxEur = "ITRX_EUR_IG_S2";
	// Proportion des tranches filles a projete sur le traxx Eur et sur le traxx US
	ARM_Vector PropTraxxEuro =  ARM_Vector(NbTranche,1);
	for (i=0;i<NbTranche;i++)
		PropTraxxEuro.Elt(i) = DataTraxx[i+1];
	matrice->Push((ARM_Vector*)PropTraxxEuro.Clone(),traxEur);
		
	char* traxUS = "CDX_NA_IG_S3";
	for (i=0;i<NbTranche;i++)
		PropTraxxEuro.Elt(i) = 1 - DataTraxx[i+1];
	matrice->Push((ARM_Vector*)PropTraxxEuro.Clone(),traxUS);

	
	if (TestForcedCorrelCDO)
	{
		char* CorrelDownCDO = "FORCED_CORREL_DOWN_UNDERLYING";
		ARM_Vector ForcedCorrelDownUnderlying = ARM_Vector(ForcedCorrelCDODown);
		matrice->Push((ARM_Vector*)ForcedCorrelDownUnderlying.Clone(),CorrelDownCDO);

		char* CorrelUpCDO = "FORCED_CORREL_UP_UNDERLYING";
		ARM_Vector ForcedCorrelUpUnderlying = ARM_Vector(ForcedCorrelCDOUp);
		matrice->Push((ARM_Vector*)ForcedCorrelUpUnderlying.Clone(),CorrelUpCDO);
		
		char* StrikeDownCDO = "FORCED_STRIKE_DOWN_UNDERLYING";
		ARM_Vector* ForcedStrikeDownUnderlying = TdownTup->RowAsVector(0);	// on recupere les strikes down des CDOs
		matrice->Push((ARM_Vector*)ForcedStrikeDownUnderlying->Clone(),StrikeDownCDO);

		char* StrikeUpCDO = "FORCED_STRIKE_UP_UNDERLYING";
		ARM_Vector* ForcedStrikeUpUnderlying = TdownTup->RowAsVector(1);	// On recupere les strike up des CDOs
		matrice->Push((ARM_Vector*)ForcedStrikeUpUnderlying->Clone(),StrikeUpCDO);
	
		MyDelete(ForcedStrikeDownUnderlying);
		MyDelete(ForcedStrikeUpUnderlying);
	}

	// ------------------------------------------------
	// Correl par strike
	// Le nombre de ligne est soit : 2 soit 4. Sinon il y a une erreur
	vector<double> CorrelStrike;
	CorrelStrike.resize(4*nbColCorrel);
	for (i=0;i<CorrelParStrike.size();i++)
		CorrelStrike[i] = atof((char*)CorrelParStrike[i]);
	
	if ((nbRowCorrel != 2) && (nbRowCorrel != 4))
	{
		result.setMsg ("ARM_ERR: Correl par strike : Size Error");
		return ARM_KO;
	}
	else if ((nbRowCorrel == 2)) // Cas ou l'utilisateur ne renseigne que le traxx euro. On duplique l'info
	{
		for (i=0;i<nbColCorrel;i++)
		{
			CorrelStrike[2*nbColCorrel+i] = CorrelStrike[i];
			CorrelStrike[3*nbColCorrel+i] = CorrelStrike[nbColCorrel + i];
		}
	}
	// ------------------------------------------------
	// Data Shift
		// Date
		string DateToShiftString = "NONE";
		if (DateToShift != CREDIT_DEFAULT_VALUE)
		{
			if (DateToShift>100)
			{
				DateToShift -= AsOfExcelDate;
				DateToShift /= Basis;
			}
			_gcvt(DateToShift,NBDECIMAL,(char*)DateToShiftString.c_str());
			strcat((char*)DateToShiftString.c_str(),"Y");
		}
		
		double ValShift = 100.*valShift; // On met la valeur en %
	// ------------------------------------------------
	// Methode NUMERIQUE
	// 0 - qGAUSSLEGENDRE; 	// 1 - qGAUSSHERMITE; 	// 2 - qTRAPEZE
	int NbStep = atoi((char*)NumMethod[0]);
	bool FastStrip = true;
	
	double MethodNum = 1.;
	if (NumMethod.size() > 1)
	{
		if ((NumMethod[1] == "l") || (NumMethod[1] == "L"))
			MethodNum = 0.;
		else if ((NumMethod[1] == "t") || (NumMethod[1] == "T"))
			MethodNum = 2.;

		if (NumMethod.size() > 2)
		{
			string Stripping = NumMethod[2].c_str();
			FastStrip = FastStripping(Stripping);
		}	
	}
// --------------------------------------------------
// --------------------------------------------------
//			CONSTRUCTION DES OBJETS
// --------------------------------------------------
// --------------------------------------------------
	ARM_Security** securities = new ARM_Security*[NbTranche];
// --------------------------------------------------	
// CDO filles
	// --------------------------------------------------------------------------
	// Frequency
	int  Frequency = 0;
	if ( (fabs(Pas)>1E-6) )
		Frequency = FindFrequency(1./Pas,Maturity);

	// Donnees permettant d'utiliser les constructeurs de Mez
	double SubAmountCDO = 0.;
	double MezzAmountCDO = 0.;
	double FloatingPayerAmount = 0.;
	int stubrule = 1;
	int FrequencyDL = Frequency;

	for (i=0;i<NbTranche;i++)
	{
		int NbNomTemp = 0;
		// On recupere les labels et le nombre de nom qui se trouve dans la tranche i
		vector<string> VectLabelTemp;
		for (int j=0;j<NbName;j++)
		{
			if (Overlap->Getvalue(j,i) == 1)
			{
				NbNomTemp++;
				VectLabelTemp.push_back(Label[j]);
			}
		}

		vector<double> IssuersNotional;
		IssuersNotional.resize(NbNomTemp);
		for (int k = 0;k<NbNomTemp;k++)
			IssuersNotional[k] = 1e7;	// Valeur en dur assurant qu'il n'y ait pas de probleme avec l'algo pgcd

		SubAmountCDO = TdownTup->Getvalue(0,i)*NbNomTemp*1e7;			
		MezzAmountCDO= (TdownTup->Getvalue(1,i) - TdownTup->Getvalue(0,i))*NbNomTemp*1e7;
		FloatingPayerAmount = MezzAmountCDO;
		
		// Construction tranche fille 
		ICM_Mez* mez = new ICM_Mez(AsOf,
								   EndDate,	
								   (ARM_Date*)0,
			  					   (ARM_Date*)0,
								   1.,			// Fixed Rate : Inutile ici puisqu'on ne calcule que des implied spread
								   K_ADJUSTED, // intRule
								   K_ADJUSTED, // adjStartDate
								   SubAmountCDO,
								   MezzAmountCDO,
								   VectLabelTemp,
								   IssuersNotional,
								   // NbNomTemp,
								   Frequency,
								   KACTUAL_360,
								   MezzAmountCDO,
								   qACCRUED_SETTLED,
								   ARM_DEFAULT_COUNTRY,
								   FloatingPayerAmount,
								   stubrule,
								   0,
								   FrequencyDL,
								CREDIT_DEFAULT_VALUE, // const double& Binary /* = CREDIT_DEFAULT_VALUE*/ ,
								std::string(), //  const std::string&  payCalName /* = NULL*/ ,
					qRunning_Leg, // const qCredit_Leg_Type& TypeFeeLeg/* =qRunning_Leg*/,
					qStandart_Recovery_Leg, // const qCredit_Leg_Type& TypeDefLeg/* =qStandart_Recovery_Leg*/,
								INCLUDE_MATURITY // const bool& IncludeMaturity /* = INCLUDE_MATURITY*/ ); 
								   );

		securities[i] = (ARM_Security*)mez->Clone();


		//DELETE
		MyDelete(mez);
	}	
// --------------------------------------------------
// Portefeuille
	ICM_Portfolio* ptf = new ICM_Portfolio(securities,NbTranche,matrice);
// --------------------------------------------------
	// On commence par calculer l'amout total du portefeuille constitué des tranches filles.
	double NominalPortfolio = 0.;
	for (i=0;i<NbTranche;i++)
		NominalPortfolio += NotTranche[i];

	double SubAmountCDO2 = StrikeDownCDO2*NominalPortfolio;
	double MezzAmountCDO2 = (StrikeUpCDO2 - StrikeDownCDO2)*NominalPortfolio;	
// --------------------------------------------------
// CDO2
	ICM_Cdo2* cdo2 = new ICM_Cdo2(AsOf,
								  EndDate,
								  (ARM_Date*)0,
								  (ARM_Date*)0,
								  Coupon,
								  K_ADJUSTED,	// intRule
								  K_ADJUSTED,	// adjStartDate
								  SubAmountCDO2,
								  MezzAmountCDO2,
								  Frequency,
								  KACTUAL_360,
								  qACCRUED_SETTLED,
								  ARM_DEFAULT_COUNTRY,
								  stubrule,
								  ptf,
							 	  CREDIT_DEFAULT_VALUE,
								  FrequencyDL,
								CREDIT_DEFAULT_VALUE, // const double& Binary /* = CREDIT_DEFAULT_VALUE*/ ,
								std::string(), // const std::string &  payCalName /* = NULL*/ ,
								false, // const bool& CrossSub /* = false*/ ,
								INCLUDE_MATURITY  // const bool& IncludeMaturity /* = INCLUDE_MATURITY */ );
								  );
// --------------------------------------------------
// Courbe de Taux
	ARM_ZeroFlat* CrbTaux = new ARM_ZeroFlat(AsOf, Taux);
// --------------------------------------------------
// Courbes de défaut.
	std::vector<const ICM_DefaultCurve*> TabCrbDefaut (  NbName+2 );  // Le rajout de 2 est lié aux courbes de traxx
	

	if (FastStrip)
	{
		for (i = 0;i<NbName;i++)
		{
			ARM_Vector* Spread = NameSpreads->TruncRowAsVector(i,nbColDefaultCurve); 

			TabCrbDefaut[i] = new ICM_DefCurvStr(Spread,
												 &PlotNum,
												 RecovCalib[i],
												 (ARM_ZeroCurve*)CrbTaux,
												 Label[i]);
			// Delete
			MyDelete(Spread);
		}

		// On cree les 2 courbes de Traxx
		TabCrbDefaut[NbName] = new ICM_DefCurvStr(&SpreadTraxxEur,
												 &PlotNum,
												 RecovTraxx,
												 (ARM_ZeroCurve*)CrbTaux,
												 "ITRX_EUR_IG_S2");

		TabCrbDefaut[NbName+1] = new ICM_DefCurvStr(&SpreadTraxxUS,
												   &PlotNum,
												   RecovTraxx,
												   (ARM_ZeroCurve*)CrbTaux,
												   "CDX_NA_IG_S3");
	}
	else
	{
		for (i = 0;i<NbName;i++)
		{
			ARM_Vector* Spread = NameSpreads->TruncRowAsVector(i,nbColDefaultCurve); 


			TabCrbDefaut[i] = new ICM_Constant_Piecewise(AsOf,
													     Term,
														 (ARM_Vector*)Spread->Clone(),
														 RecovCalib[i],
														 CrbTaux,
												K_ADJUSTED,	// intRule
												K_ADJUSTED,	// adjStartDate
														 qCredit_Adjust20,
														 ARM_DEFAULT_COUNTRY,
														 Label[i],
														 false,
														 NULL,//2 NULL,
														 K_QUARTERLY,qDEFCURVE_DICHO,"STD", ARM_Currency(ARM_DEFAULT_COUNTRY).GetCreditStartDateLag(),ICM_Parameters());
			// Delete
			MyDelete(Spread);
		}

		// On cree les 2 courbes de Traxx
		TabCrbDefaut[NbName] = new ICM_Constant_Piecewise(AsOf,
													     Term,
														 (ARM_Vector*)SpreadTraxxEur.Clone(),
														 RecovTraxx,
														 CrbTaux,
												K_ADJUSTED,	// intRule
												K_ADJUSTED,	// adjStartDate
														 qCredit_Adjust20,
														 ARM_DEFAULT_COUNTRY,
														 "ITRX_EUR_IG_S2",
														 false,
														 NULL,//2 NULL,
														 K_QUARTERLY,qDEFCURVE_DICHO,"STD", ARM_Currency(ARM_DEFAULT_COUNTRY).GetCreditStartDateLag(),
														 ICM_Parameters());

		TabCrbDefaut[NbName+1] = new ICM_Constant_Piecewise(AsOf,
													     Term,
														 (ARM_Vector*)SpreadTraxxUS.Clone(),
														 RecovTraxx,
														 CrbTaux,
												K_ADJUSTED,	// intRule
												K_ADJUSTED,	// adjStartDate
														 qCredit_Adjust20,
														 ARM_DEFAULT_COUNTRY,
														 "CDX_NA_IG_S3",
														 false,
														 NULL,//2 NULL,
														 K_QUARTERLY,qDEFCURVE_DICHO,"STD", 
														 ARM_Currency(ARM_DEFAULT_COUNTRY).GetCreditStartDateLag()
														 ,ICM_Parameters());
	}
// --------------------------------------------------
// Correlation
	int NbNomTraxx = 125;
	ARM_Vector MaturityTraxx = ARM_Vector(1,MaturiteTraxx);
	
	ARM_Vector proportion =  ARM_Vector(2,1);
	proportion.Elt(0) = 1.;	 // Proportion Traxx EUR
	proportion.Elt(1) = 0.;	 // Proportion Traxx US

	std::vector<std::string> LabelTraxx(2); 
	LabelTraxx[0] = "ITRX_EUR_IG_S2";
	LabelTraxx[1] = "CDX_NA_IG_S3";


	std::vector<std::string> LabelTraxxEur(NbNomTraxx); 
	// char** LabelTraxxEur = new char*[NbNomTraxx];
	for (i=0;i<NbNomTraxx;i++)
		LabelTraxxEur[i] = "ITRX_EUR_IG_S2";
		
	std::vector<std::string> LabelTraxxUS(NbNomTraxx); 
	// char** LabelTraxxUS = new char*[NbNomTraxx];
	for (i=0;i<NbNomTraxx;i++)
		LabelTraxxUS[i] = "CDX_NA_IG_S3";

	ARM_Vector StrikEUR = ARM_Vector(nbColCorrel,0.);
	for (i=0;i<nbColCorrel;i++)
		StrikEUR.Elt(i) = 100 * CorrelStrike[i];

	ARM_Matrix* CorrelationEur = new ARM_Matrix(MaturityTraxx.GetSize(),nbColCorrel);
	for (i=0;i<nbColCorrel;i++)
		CorrelationEur->Elt(0,i) = 100 * CorrelStrike[i + nbColCorrel];

	ARM_Vector StrikUS =  ARM_Vector(nbColCorrel,0.);
	for (i=0;i<nbColCorrel;i++)
		StrikUS.Elt(i) = 100 * CorrelStrike[2*nbColCorrel+ i];

	ARM_Matrix* CorrelationUS = new ARM_Matrix(MaturityTraxx.GetSize(),nbColCorrel);
	for (i=0;i<nbColCorrel;i++)
		CorrelationUS->Elt(0,i) =  100 * CorrelStrike[i + 3*nbColCorrel];


	// Objet Vol Curve : on doit en faire un pour chaque Traxx
	vector<const ARM_VolCurve*> Vector_VolCurve(2);
	Vector_VolCurve[0] = new ARM_VolLInterpol(AsOf,
											 (ARM_Vector*)MaturityTraxx.Clone(),
											  (ARM_Vector*)StrikEUR.Clone(),				// Strike en bp
											  (ARM_Matrix*)CorrelationEur->Clone(),		// Valeur de la correl en %
											  1,
	 										  K_SMILE_VOL);

	Vector_VolCurve[1] = new ARM_VolLInterpol(AsOf,
											  (ARM_Vector*)MaturityTraxx.Clone(),
											  (ARM_Vector*)StrikUS.Clone(),				// Strike en bp
											  (ARM_Matrix*)CorrelationUS->Clone(),		// Valeur de la correl en %
											  1,
											  K_SMILE_VOL);


	// On cree un index pour chaque traxx
	vector<const ICM_Credit_Index*> VIndex(2);
	ARM_Vector YT(1); YT.Elt(0) = 5.0;
	ARM_Vector spread(1); spread.Elt(0) = 0.;
	VIndex[0] = new ICM_Credit_Index(KACTUAL_360, 
									  4, 
									  4, 
									  YT,
									  ARM_DEFAULT_COUNTRY,
									  LabelTraxxEur[0],
									  LabelTraxxEur,
									  qAVERAGE,spread,NULL,0,0,K_ARREARS,0,K_ARREARS,0,qCredit_Adjust20,5,2
									  );

	VIndex[1]  = new ICM_Credit_Index(KACTUAL_360, 
									  4, 
									  4, 
									  YT,
									  ARM_DEFAULT_COUNTRY,
									  LabelTraxxUS[0],
									  LabelTraxxUS,
									  qAVERAGE,spread,NULL,0,0,K_ARREARS,0,K_ARREARS,0,qCredit_Adjust20,5,2);

// FIXMEFRED: mig.vc8 (28/05/2007 14:20:33):cast
	ICM_Smile_Correlation* Corr = new ICM_Smile_Correlation(AsOf,
														"CORRSTR",
														&(*Vector_VolCurve.begin()),
														LabelTraxx,
														&proportion,
														&(*VIndex.begin()));
// ---------------------------------------------------------------------------
// Model
ICM_ModelMultiCurves* Model = new ICM_ModelMultiCurves(
													   TabCrbDefaut,
													   CrbTaux, 
													   RecovLoss,
													   Corr);
// ---------------------------------------------------------------------------
// Pricer
	// Construction de la matrix
	ICM_Parameters* MatricePricer = new ICM_Parameters();
	
	// Methode underlying
	// On fixe en dur le nombre de pas a utiliser pour les cdo underlying
	char* MethodeUnderlying = "INTEGRATION_METHOD_UNDERLYING";
	ARM_Vector MethUnderlying =  ARM_Vector(1.,1);		// 1 -> qGauss_HERMITE
	MatricePricer->Push((ARM_Vector*)MethUnderlying.Clone(),MethodeUnderlying);

	// On fixe en dur le nombre de pas a utiliser pour les cdo underlying
	char* StepUnderlying = "INTEGRATION_STEP_UNDERLYING";
	ARM_Vector PasUnderlying =  ARM_Vector(1.,40);
	MatricePricer->Push((ARM_Vector*)PasUnderlying.Clone(),StepUnderlying);
	
	// FORCED_STRIKE_DOWN_CDO2
	char* ForcedKDownCDO2 = "FORCED_STRIKE_DOWN_CDO2";
	ARM_Vector ForcedKD_CDO2 =  ARM_Vector(1.,StrikeDownCDO2);
	MatricePricer->Push((ARM_Vector*)ForcedKD_CDO2.Clone(),ForcedKDownCDO2);
	
	// FORCED_STRIKE_UP_CDO2
	char* ForcedKUpCDO2 = "FORCED_STRIKE_UP_CDO2";
	ARM_Vector ForcedKU_CDO2 =  ARM_Vector(1.,StrikeUpCDO2);
	MatricePricer->Push((ARM_Vector*)ForcedKU_CDO2.Clone(),ForcedKUpCDO2);

	// FORCED_CORREL_DOWN_CDO2
	char* ForcedCorrelDownCDO2 = "FORCED_CORREL_DOWN_CDO2";
	ARM_Vector ForcedRhoD_CDO2 =  ARM_Vector(1.,CorrelForced[0]);
	MatricePricer->Push((ARM_Vector*)ForcedRhoD_CDO2.Clone(),ForcedCorrelDownCDO2);

	// FORCED_CORREL_UP_CDO2
	char* ForcedCorrelUpCDO2 = "FORCED_CORREL_UP_CDO2";
	ARM_Vector ForcedRhoU_CDO2 =  ARM_Vector(1.,CorrelForced[1]);
	MatricePricer->Push((ARM_Vector*)ForcedRhoU_CDO2.Clone(),ForcedCorrelUpCDO2);
	
	// BETA
	char* BETA = "BETA";
	ARM_Vector BetaForced =  ARM_Vector(1.,-999.);
	MatricePricer->Push((ARM_Vector*)BetaForced.Clone(),BETA);
	
	// NUM_SIMUL_CORREL
	char* SimulCorrel = "NUM_SIMUL_CORREL";
	ARM_Vector SimulationCorrel =  ARM_Vector(1.,10.);
	MatricePricer->Push((ARM_Vector*)SimulationCorrel.Clone(),SimulCorrel);

	// UNDER_RESCAL
	char* RescalUnderlying = "UNDER_RESCAL";
	ARM_Vector Rescal =  ARM_Vector(1.,0.);
	if (!TestForcedCorrelCDO)
		Rescal.Elt(0) = 1.;
	MatricePricer->Push((ARM_Vector*)Rescal.Clone(),RescalUnderlying);

	// COPULA
	char* Copula = "COPULA";
	ARM_Vector Copul = ARM_Vector(1,1);
	MatricePricer->Push((ARM_Vector*)Copul.Clone(), (char*)Copula);

	// INTEGRATION_STEP_1
	char* Step = "INTEGRATION_STEP_1";
	ARM_Vector StepIntegration = ARM_Vector(1,NbStep);
	MatricePricer->Push((ARM_Vector*)StepIntegration.Clone(), (char*)Step);
	
	// FREEDOM
	const char* freedom = "FREEDOM_DEGREE";
	ARM_Vector FreedomVector = ARM_Vector(1,1);
	MatricePricer->Push((ARM_Vector*)FreedomVector.Clone(), (char*)freedom);
	
	// INTEGRATION METHOD
	const char* tmp3 = "INTEGRATION_METHOD";
	// 0 -> qGAUSSLEGENDRE; // 1 -> qGAUSSHERMITE;	// 2 -> qTRAPEZE
	ARM_Vector V3 = ARM_Vector(1,MethodNum );
	MatricePricer->Push((ARM_Vector*)V3.Clone(), (char*)tmp3);


	ICM_Pricer_Analytic_Cdo2_Smile  Pricer ; Pricer .Set(cdo2,Model,*MatricePricer,AsOf);
// ---------------------------------------------------------------------------
// Fonction
	double Prix = 0.;
	if (C_Function == "NPV")
		Prix = - Pricer.Price(qCMPPRICE);
	else if (C_Function == "DL")
		Prix = Pricer.Price(qCMPDEFLEGPV);
	else if (C_Function == "FL")
		Prix = Pricer.Price(qCMPFEELEGPV);
	else if (C_Function == "BE")
		Prix = Pricer.ComputeSpread()/10000.;
	else if (C_Function == "DUR")
		Prix = Pricer.Price(qCMPDURATION);
	else if (C_Function == "SENSI") //  Sensi au spread
		Prix = - Pricer.Hedge(ICMSPREAD_TYPE,(char*)DateToShiftString.c_str(),CCSTringToSTLString(LabelToShift),ValShift);

	// Resultat
	result.setDouble(Prix);
		
	// DELETE	
	// MyDelete(Pricer);
	MyDelete(MatricePricer);
	MyDelete(Model);
	MyDelete(Corr);
	/*MyDelete(VIndex,2); */  delete VIndex[0]; delete VIndex[1];
	/*MyDelete(Vector_VolCurve,2);*/ delete Vector_VolCurve[0]; delete Vector_VolCurve[1];
	MyDelete(CorrelationUS);
	MyDelete(CorrelationEur);
	// MyDeleteTab(LabelTraxxUS);
	// MyDeleteTab(LabelTraxxEur);
	// MyDeleteTab(LabelTraxx);
	MyDeleteTab(TabCrbDefaut);
	MyDelete(cdo2);
	MyDelete(ptf);
	MyDelete(securities,NbTranche);
	MyDelete(matrice);
	MyDeleteTab(RecovCalib);
	// MyDeleteTab(RecovLoss);
	MyDelete(Overlap);

	if (Label)
		FreePointerTabChar(Label, NbName);

	MyDelete(TdownTup);
	MyDeleteTab(pEndDate);

	return ARM_OK;
	}
 
	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}

	catch (...)
	{
		result.setMsg ("ARM_ERR: CDO2 SMILE : unrecognized failure");
		return ARM_KO;
	}



}

// -----------------------------------------------------------------------------------
// BASKETS : Interface Proto / Code ARM Damien
extern long ICMLOCAL_STR_Basket(double YF_Maturity,
							    vector<double> Defaults,
								double Nb_Issuers,
								double IssuerNotional,
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
								ARM_result& result)
{
	ICMMSG(WARN,"Using ICMLOCAL_STR_Basket" ); 
	CCString msg ("");

	int i,j = 0;
	double Basis = 365.;
	double Taux = taux*100.;
	double CteARMDateToExcelDate = 2415019.;
	ARM_Date AsOf; 
	AsOf.Today();
	double AsOfExcelDate = AsOf.GetJulian() - 2415019;
	double Maturity = YF_Maturity;
	if (YF_Maturity>100)
		Maturity = ( YF_Maturity - AsOfExcelDate)/Basis;

	ARM_Date EndDate;
	char* pEndDate =  new char[11];
	if (Maturity<100)
	{
		Local_XLDATE2ARMDATE(AsOfExcelDate + Basis*Maturity,pEndDate);
		EndDate = (ARM_Date) pEndDate;
	}
	else
	{
		Local_XLDATE2ARMDATE(Maturity,pEndDate);
		EndDate = (ARM_Date) pEndDate;
	}
	

	double price = 0.;
	
	try
	{
	// -------- Messages d'erreur --------------------
	if (IssuerNotional<10000)
	{
		result.setMsg ("ARM_ERR: le notional d'un nom doit etre supérieur à 10000");
		return ARM_KO;
	}		
	// -----------------------------------------------
	// Issuers Data
	int NbIssuers = (int)Nb_Issuers;
	// -----------------------------------------------
	// Recovery	
	double* CalibRecov =  new double[NbIssuers];
	memset(CalibRecov,'\0',sizeof(double)*NbIssuers);

	if (CalibRecovery.size()<NbIssuers)
		for (i=0;i<NbIssuers;i++)
			CalibRecov[i] = CalibRecovery[0];
	else
		for (i=0;i<NbIssuers;i++)
			CalibRecov[i] = CalibRecovery[i];


	ARM_Vector LossRecov (NbIssuers);
	// memset(LossRecov,'\0',sizeof(double)*NbIssuers);

	if (LossRecovery.size()<NbIssuers)
		for (i=0;i<NbIssuers;i++)
			LossRecov[i] = LossRecovery[0];
	else
		for (i=0;i<NbIssuers;i++)
			LossRecov[i] = LossRecovery[i];
	// -----------------------------------------------
	// Matrice de Spread
	int nb_plot = YF_plot.size();
	int NbRow = NbIssuers;
	int NbCol = nb_plot; 
	ICM_QMatrix<double>* NameSpreads = new ICM_QMatrix<double>(NbRow,NbCol,CREDIT_DEFAULT_VALUE);
	for (i=0;i<NbRow;i++)
	{
		// Cas Vraie Matrice
		if (MatIssuersSpread.size()==NbCol*NbRow)
		{
			for (j=0;j<NbCol;j++)
				NameSpreads->SetValue(i,j,MatIssuersSpread[i*NbCol+j]);
		}
		// Cas d'un input = une ligne
		else if (MatIssuersSpread.size() == NbCol)
		{
			for (j=0;j<NbCol;j++)
				NameSpreads->SetValue(i,j,MatIssuersSpread[j]);
		}
		// Cas d'une colonne
		else if (MatIssuersSpread.size() == NbIssuers)
		{
			for (j=0;j<NbCol;j++)
				NameSpreads->SetValue(i,j,MatIssuersSpread[i]);
		}
		// Cas une seule valeur
		else
		{
			for (j=0;j<NbCol;j++)
				NameSpreads->SetValue(i,j,MatIssuersSpread[0]);
		}

	}
	// -----------------------------------------------
	// Vecteur de plot
	ARM_Vector PlotNum = ARM_Vector(YF_plot);
	if (YF_plot[0]>100)
	{
		for (i=0;i<YF_plot.size();i++)
			PlotNum.Elt(i) = (YF_plot[i] - AsOfExcelDate)/Basis;
	}
	// Term
	int sizePlot = PlotNum.GetSize();

	vector<string> Term;
	for (i=0;i<sizePlot;i++)
	{
		char * t = NULL;
		ITOA(YF_plot[i],t,10);
		string s = string(t) + string("Y");
		Term.push_back(s);
		if (t) delete t;
	}
	// -----------------------------------------------
	//  Date To Shift
	char* DateToShiftString = new char[10];
	strcpy(DateToShiftString,"NONE");

	if (DateToShift != CREDIT_DEFAULT_VALUE)
	{
		ITOA(DateToShift,DateToShiftString,10);
		strcat(DateToShiftString,"Y");
	}
	// -----------------------------------------------
	// Label
	std::vector<std::string> Label (NbIssuers) ; 
	// char** Label = new char*[NbIssuers];
	for (int i =0; i<NbIssuers;i++)
		Label[i] = IssuersLabel[i];
	// -----------------------------------------------
	// Betas	
	int SizeBeta = IssuersBeta.size();
	ARM_Vector BetasIssuers = ARM_Vector(NbIssuers,1.);;
	if ( (SizeBeta<NbIssuers) && (SizeBeta == 1))
	{
		for (i=0;i<NbIssuers;i++)
			BetasIssuers.Elt(i) = IssuersBeta[0];
	}
	else if (SizeBeta >= NbIssuers)
	{
		for (i=0;i<NbIssuers;i++)
			BetasIssuers.Elt(i) = IssuersBeta[i];
	}
	else 
	{
		result.setMsg ("ARM_ERR: Vector Betas : Size Error");
		return ARM_KO;
	}	
	// -----------------------------------------------
	// Frequency
	int Frequency = FindFrequency(FrequencyFeeLeg,Maturity);
	// -----------------------------------------------
	// Vecteur d'Issuers Notionals
	// double* pIssuersNotionals = new double[NbIssuers];
	std::vector<double> IssuersNotionals (NbIssuers,IssuerNotional); 
	// for (int il=0; il<NbIssuers; il++) 
	// 	pIssuersNotionals[il]=IssuerNotional;
	// -----------------------------------------------
	// Noms sujets a la protection
	if (Defaults.size() == 1)
		Defaults.push_back(Defaults[0] + 1.);

	int FirstNumDefault = (int) Defaults[0];
	int LastNumDefault  = (int) Defaults[1];
	// -----------------------------------------------
	// ValShift
	ValShift *= 100.;
	// --------------------------------------------------	
	// Stripping
	string Stripping = C_Stripping.c_str();
	bool FastStrip = FastStripping(Stripping);
	// Courbe de Taux
	ARM_ZeroFlat* CrbTaux = new ARM_ZeroFlat(AsOf, Taux);
	// -----------------------------------------------
	// Courbes de défaut.
	std::vector<const ICM_DefaultCurve*> TabCrbDefaut (NbIssuers);
	
if (FastStrip)
{
	for (i = 0;i<NbIssuers;i++)
	{
		ARM_Vector* Spread = NameSpreads->TruncRowAsVector(i,nb_plot); 

		TabCrbDefaut[i] = new ICM_DefCurvStr(Spread,
											 &PlotNum,
											 CalibRecov[i],
											 CrbTaux,
											 Label[i]);
		// Delete
		MyDelete(Spread);
	}
}
else
{
	for (i = 0;i<NbIssuers;i++)
	{
		ARM_Vector* Spread = NameSpreads->TruncRowAsVector(i,nb_plot); 

		TabCrbDefaut[i] = new ICM_Constant_Piecewise(AsOf,
												     Term,
													(ARM_Vector*)Spread->Clone(),
													CalibRecov[i],
													CrbTaux,
												K_ADJUSTED,	// intRule
												K_ADJUSTED,	// adjStartDate
													qCredit_Adjust20,
													ARM_DEFAULT_COUNTRY,	
													Label[i],
													false,
													NULL,//2 NULL,
													K_QUARTERLY,qDEFCURVE_DICHO,"STD", 
													ARM_Currency(ARM_DEFAULT_COUNTRY).GetCreditStartDateLag()
													,ICM_Parameters());
		// Delete
		MyDelete(Spread);
	}

}
// Correlation
ICM_Beta_Correlation* Betas = new ICM_Beta_Correlation(AsOf,"STRCORR",BetasIssuers,Label,(ARM_IRIndex*)0,(ARM_IRIndex*)0);
	

// Objet NthTD
int stubrule = 1;
ICM_Nthtd* NthTD = new ICM_Nthtd(AsOf,
								 EndDate,
								 (ARM_Date*)0,
			     			     (ARM_Date*)0,
								 FixedRate,
								 K_ADJUSTED,	// intRule,
								 K_ADJUSTED,	// adjStartDate
								 FirstNumDefault,
								 LastNumDefault,
								 //NbIssuers,
								 Label,
								 IssuersNotionals,
								 Frequency,
								 KACTUAL_360,
								 IssuerNotional,
								 qACCRUED_SETTLED,
								 ARM_DEFAULT_COUNTRY,
								 IssuerNotional,
								 stubrule,
								 CREDIT_DEFAULT_VALUE,
								 Frequency,
								CREDIT_DEFAULT_VALUE, // const double& Binary /* = CREDIT_DEFAULT_VALUE*/,
								std::string(), // const std::string&  Paycal /* = NULL*/,
								INCLUDE_MATURITY // const bool& IncludeMaturity /* = INCLUDE_MATURITY*/) ;
								 );
// Model
ICM_ModelMultiCurves* Model = new ICM_ModelMultiCurves(// NbIssuers,
													   TabCrbDefaut,
													   (ARM_ZeroCurve*) CrbTaux, 
													   LossRecov,
													   Betas);
// Pricer
//ICM_Pricer_Distrib* Pricer_NTD = new ICM_Pricer_Distrib(NthTD,Model);
ICM_Pricer_Distrib_Smile Pricer_NTD ; Pricer_NTD .Set(NthTD,Model,ICM_Parameters(),AsOf);

// Resu
	if (C_Function == "NPV")
		price = - Pricer_NTD.Price(qCMPPRICE);
	if (C_Function == "DL")
		price = Pricer_NTD.Price(qCMPDEFLEGPV);
	if (C_Function == "FL")
		price = Pricer_NTD.Price(qCMPFEELEGPV);
	if (C_Function == "BE")
		price = Pricer_NTD.ComputeSpread(0)/10000.;
	if (C_Function == "DUR")
		price = Pricer_NTD.Price(qCMPDURATION);
	if (C_Function == "SENSI") //  Sensi au spread
		price = - Pricer_NTD.Hedge(ICMSPREAD_TYPE,DateToShiftString,CCSTringToSTLString(LabelToShift),ValShift);


	result.setDouble(price);
		
	// DELETE
	// MyDelete(Pricer_NTD);
	MyDelete(Model);
	MyDelete(NthTD);
	MyDelete(Betas);
	MyDeleteTab(TabCrbDefaut);
	MyDeleteTab(DateToShiftString);
	// MyDeleteTab(pIssuersNotionals);

	// if (Label)
	// 	FreePointerTabChar(Label,NbIssuers);
	
	MyDelete(NameSpreads);

	MyDeleteTab(CalibRecov);
	// MyDeleteTab(LossRecov);
	MyDeleteTab(pEndDate);

	return ARM_OK;
	}
 
	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}

	catch (...)
	{
		result.setMsg ("ARM_ERR: Price : unrecognized failure");
		return ARM_KO;
	}
};


// Pricing de Tranche
// -----------------------------------------------------------------------------------
// CDO : Interface Proto / Code ARM Damien
extern long ICMLOCAL_STR_CDO(int MezType,
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
							 ARM_result& result)
{

	ICMMSG(WARN,"Using ICMLOCAL_STR_CDO" ); 
	CCString msg ("");

	int i,j = 0;
	double price = 0.;
	
	try
	{
		double Basis = 365.;
	// --------------------------------------------------------------------------
	// Message d'erreur
	if (Notional[0]<10000)
	{
		result.setMsg ("ARM_ERR: le notional d'un nom doit etre supérieur à 10000");
		return ARM_KO;
	}
	if (DataMez.size() != 6)
	{
		result.setMsg ("ARM_ERR: Le vecteur de Data Mez est mal renseigné");
		return ARM_KO;
	}
// --------------------------------------------------------------------------
// ARM Compliant
	double Taux = taux*100.;
	double ValBumpSensi = ValShift * 100.; 
// --------------------------------------------------------------------------
// Date
	double CteARMDateToExcelDate = 2415019.;
	ARM_Date AsOf; 
	AsOf.Today();
	double AsofJulian = AsOf.GetJulian();
	double AsOfExcelDate = AsofJulian - 2415019;

	double Maturity = DataMez[0];
	MaturityYF(Maturity,AsofJulian);

	ARM_Date EndDate(AsofJulian + Maturity*Basis);;
// --------------------------------------------------------------------------
// Issuers Data
	int NbIssuers = (int)DataMez[5];
// --------------------------------------------------------------------------
// Recovery	
	CompleteVector(CalibRecovery,NbIssuers);
	CompleteVector(LossRecovery,NbIssuers);

	double* CalibRecov =  new double[NbIssuers];
	memset(CalibRecov,'\0',sizeof(double)*NbIssuers);

	ARM_Vector LossRecov (NbIssuers);
	// memset(LossRecov,'\0',sizeof(double)*NbIssuers);

	for (i=0;i<NbIssuers;i++)
	{
		CalibRecov[i] = CalibRecovery[i];
		LossRecov[i] = LossRecovery[i];
	}
// --------------------------------------------------------------------------
// Matrice de Spread
	int nb_plot = YF_plot.size();
	int NbRow = NbIssuers;
	int NbCol = nb_plot; 
	ICM_QMatrix<double>* NameSpreads = new ICM_QMatrix<double>(NbRow,NbCol,CREDIT_DEFAULT_VALUE);
	for (i=0;i<NbRow;i++)
	{
		// Cas Vraie Matrice
		if (MatIssuersSpread.size()==NbCol*NbRow)
		{
			for (j=0;j<NbCol;j++)
				NameSpreads->SetValue(i,j,MatIssuersSpread[i*NbCol+j]);
		}
		// Cas d'un input = une ligne
		else if (MatIssuersSpread.size() == NbCol)
		{
			for (j=0;j<NbCol;j++)
				NameSpreads->SetValue(i,j,MatIssuersSpread[j]);
		}
		// Cas d'une colonne
		else if (MatIssuersSpread.size() == NbIssuers)
		{
			for (j=0;j<NbCol;j++)
				NameSpreads->SetValue(i,j,MatIssuersSpread[i]);
		}
		// Cas une seule valeur
		else if (MatIssuersSpread.size() == 1)
		{
			for (j=0;j<NbCol;j++)
				NameSpreads->SetValue(i,j,MatIssuersSpread[0]);
		}
		else
		{
			result.setMsg ("ARM_ERR: Issuers Spread : Size Error");
			return ARM_KO;
		}	

	}
// --------------------------------------------------------------------------
// Vecteur de plot
	ARM_Vector PlotNum = ARM_Vector(YF_plot);
	if (YF_plot[0]>100)
	{
		for (i=0;i<YF_plot.size();i++)
			PlotNum.Elt(i) = (YF_plot[i] - AsOfExcelDate)/Basis;
	}

	int sizePlot = YF_plot.size();
	
	vector<string> Term;
	for (i=0;i<sizePlot;i++)
	{
		char * t = NULL;
		ITOA(YF_plot[i],t,10);
		string s = string(t) + string("Y");
		Term.push_back(s);
		if (t) delete t;
	}
// --------------------------------------------------------------------------
//  Date To Shift
	string DateToShiftString = "NONE";
	if (DateToShift != CREDIT_DEFAULT_VALUE)
	{
		if (DateToShift > 100)
		{
			DateToShift -= AsOfExcelDate;
			DateToShift *= 1/Basis;
		}
		_gcvt(DateToShift,NBDECIMAL,(char*)DateToShiftString.c_str());
		strcat((char*)DateToShiftString.c_str(),"Y");
	}

	ValShift *= 100.;
// --------------------------------------------------------------------------
// Label
	if (IssuersLabel.size() < NbIssuers)
	{
		result.setMsg ("ARM_ERR: Le vecteur de Label n'est pas de la bonne taille");
		return ARM_KO;
	}
	std::vector<std::string> Label (NbIssuers); 
	// char** Label = new char*[NbIssuers];
	vector<string> LabelVect;
	LabelVect.resize(NbIssuers);
	for (int i =0; i<NbIssuers;i++)
	{
		Label[i] = IssuersLabel[i];
		LabelVect[i] = IssuersLabel[i];
	}
// --------------------------------------------------------------------------
// Betas
ARM_Vector BetasIssuers = ARM_Vector(NbIssuers,1.);
	CompleteVector(IssuersBeta,NbIssuers);
	for (i=0;i<NbIssuers;i++)
		BetasIssuers.Elt(i) = IssuersBeta[i];
// --------------------------------------------------------------------------
// Frequency
	int  Frequency = 0;
	// DataMez[3] : Pas 	// DataMez[0] : Maturity
	if (fabs(DataMez[3])>1E-6)
		Frequency = FindFrequency(1./DataMez[3],Maturity);
// --------------------------------------------------------------------------
// Vecteur d'Issuers Notionals
	std::vector<double> pIssuersNotionals (NbIssuers); 
	// double* pIssuersNotionals = new double[NbIssuers];
	// memset(pIssuersNotionals,'\0',sizeof(double)*NbIssuers);
	vector<double> VIssuersNotionals;
	VIssuersNotionals.resize(NbIssuers);
	for (int il=0; il<NbIssuers; il++) 
	{
		VIssuersNotionals[il] = Notional[0];
		pIssuersNotionals[il] = Notional[0];
	}
// --------------------------------------------------------------------------
// Points d'attachement 
	double CDOAmount = 0.;
	for (i=0;i<NbIssuers;i++)
		CDOAmount += VIssuersNotionals[i];

	double MezzAmount = (DataMez[2] - DataMez[1])*CDOAmount;
	double SubAmount  = DataMez[1] * CDOAmount;

	if (Notional.size() == 1)
		Notional.push_back(MezzAmount);
// --------------------------------------------------------------------------
// Methode NUMERIQUE
	// 0 - qGAUSSLEGENDRE; 	// 1 - qGAUSSHERMITE; 	// 2 - qTRAPEZE
	int NbStep = atoi((char*)NumMethod[0]);
	bool FastStrip = true;

	double MethodNum = 1.;
	if (NumMethod.size() > 1)
	{
		if ((NumMethod[1] == "l") || (NumMethod[1] == "L"))
			MethodNum = 0.;
		else if ((NumMethod[1] == "t") || (NumMethod[1] == "T"))
			MethodNum = 2.;

		if (NumMethod.size() == 3)
		{
			string Stripping = NumMethod[2].c_str();
			FastStrip = FastStripping(Stripping);
		}
	}
// -----------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------
// CONSTRUCTION D'OBJETS

// Courbe de Taux
	ARM_ZeroFlat* CrbTaux = new ARM_ZeroFlat(AsOf, Taux);
// -----------------------------------------------------------------------------------
// Courbes de défaut.
	std::vector<const ICM_DefaultCurve*> TabCrbDefaut (NbIssuers);
	
	if (FastStrip)
	{
		for (i = 0;i<NbIssuers;i++)
		{
			ARM_Vector* Spread = NameSpreads->TruncRowAsVector(i,nb_plot); 

			TabCrbDefaut[i] = new ICM_DefCurvStr(Spread,
												 &PlotNum,
												 CalibRecov[i],
												 CrbTaux,
												 Label[i]);
			// Delete
			MyDelete(Spread);
		}
	}
	else
	{
		for (i = 0;i<NbIssuers;i++)
		{
			ARM_Vector* Spread = NameSpreads->TruncRowAsVector(i,nb_plot) ;

			TabCrbDefaut[i] = new ICM_Constant_Piecewise(AsOf,
													Term,
													(ARM_Vector*)Spread->Clone(),
													CalibRecov[i],
													CrbTaux,
												K_ADJUSTED,	// intRule
												K_ADJUSTED,	// adjStartDate
													qCredit_Adjust20,
													ARM_DEFAULT_COUNTRY,
													Label[i],
													false,
													NULL,//2 NULL,
													K_QUARTERLY,qDEFCURVE_DICHO,"STD", 
													ARM_Currency(ARM_DEFAULT_COUNTRY).GetCreditStartDateLag()
													,ICM_Parameters());
			// Delete
			MyDelete(Spread);
		}
	}
// -----------------------------------------------------------------------------------
// Correlation
ICM_Beta_Correlation* Betas = new ICM_Beta_Correlation(AsOf,"STRCORR",BetasIssuers,Label,(ARM_IRIndex*)0,(ARM_IRIndex*)0);
// -----------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------
// Objet Vol
ARM_VolFlat* VolIndex = new ARM_VolFlat(AsOf,VolAjustConvexity*100.);	// Multiplication parce qu'on passe des valeurs en bp.
// -----------------------------------------------------------------------------------

// Model
ICM_ModelMultiCurves* Model = new ICM_ModelMultiCurves(// NbIssuers,
													   TabCrbDefaut,
													   (ARM_ZeroCurve*) CrbTaux, 
													   LossRecov,
													   Betas,VolIndex);
// Objet Index
ICM_Credit_Index* Index = NULL;
ICM_DefaultCurve* ForcedCurve = NULL;
ICM_DefCurvStr* CourbeTemp = NULL;
ARM_Vector SpreadFwd = ARM_Vector(SpreadCrbFwd);
double yearterm = 5.; // Maturité du Forward
if (MezType == 4) // Cas CMTranche
{
	if (SpreadCrbFwd[0] != CREDIT_DEFAULT_VALUE)
	{
		CourbeTemp = new ICM_DefCurvStr(&SpreadFwd,
										&PlotNum,
										CalibRecov[0],
										(ARM_ZeroCurve*)CrbTaux,
										"CrbFwd");

		ForcedCurve = (ICM_DefaultCurve*) (CourbeTemp->ConvertStrToARM())->Clone();

		MyDelete(CourbeTemp);
	}
	else 
		ForcedCurve = (ICM_DefaultCurve*) AverageCurveStr(TabCrbDefaut)->ConvertStrToARM()->Clone();

	ARM_Vector YT(1); YT.Elt(0) = yearterm;
	ARM_Vector spread(1); spread.Elt(0) = 0.;
	Index = new ICM_Credit_Index(KACTUAL_360, 
								 Frequency, 
								 Frequency, 
								 YT,
								 ARM_DEFAULT_COUNTRY,
								 Label[0],
								 Label,
								 qAVERAGE,
								 spread,
								 //NULL,
								 ForcedCurve,
								 0,0,K_ARREARS,0,K_ARREARS,0,qCredit_Adjust20,5,2);
}
// -----------------------------------------------------------------------------------
// Objet Mezz
ICM_Mez* mez = NULL;
double FloatingPayerAmount = MezzAmount;
int stubrule = 1;
int FrequencyDL = Frequency;
double FixedRate = DataMez[4];
string pRefDate = "NULL";

if (MezType != 4)
{
	if (MezType == 2) // Cas Infine
			Frequency = K_ZEROCOUPON;

	mez = new ICM_Mez(AsOf,
					  EndDate,
					  (ARM_Date*)0,
					  (ARM_Date*)0,
					  FixedRate,
					  K_ADJUSTED,	// intRule
					  K_ADJUSTED,	// adjStartDate
					  SubAmount,
					  MezzAmount,
					  LabelVect,
					  VIssuersNotionals,
					  //NbIssuers,
					  // "NULL",//pRefDate,	
					  Frequency,
					  KACTUAL_360,
					  MezzAmount,
					  qACCRUED_SETTLED,
					  ARM_DEFAULT_COUNTRY,
					  FloatingPayerAmount,
					  stubrule,
					  0,
					  FrequencyDL,
					CREDIT_DEFAULT_VALUE, // const double& Binary /* = CREDIT_DEFAULT_VALUE*/ ,
					std::string(), // const std::string&  payCalName /* = NULL*/ ,
					qRunning_Leg, // const qCredit_Leg_Type& TypeFeeLeg/* =qRunning_Leg*/,
					qStandart_Recovery_Leg, // const qCredit_Leg_Type& TypeDefLeg/* =qStandart_Recovery_Leg*/,
					INCLUDE_MATURITY // const bool& IncludeMaturity /* = INCLUDE_MATURITY*/ ); 
					  );
}
else // On est dans le cas CMTranche
{
	mez = new ICM_Mez(AsOf,
					  EndDate,
					  (ARM_Date*)0,
					  (ARM_Date*)0,
					  K_ADJUSTED,	// intRule
					  K_ADJUSTED,	// adjStartDate
					  SubAmount,
					  MezzAmount,
					  Label,
					  pIssuersNotionals,
					  //NbIssuers,
				      Index, 
					  FixedRate,
					  Frequency,
					  KACTUAL_360,
					  MezzAmount,
					  qACCRUED_SETTLED,
					  ARM_DEFAULT_COUNTRY,
					  FloatingPayerAmount,
					  stubrule,
					  0,
					  FrequencyDL,
						CREDIT_DEFAULT_VALUE, // const double& Binary /*  CREDIT_DEFAULT_VALUE*/ ,
						std::string(), // const std::string& /* char* payCalName = NULL*/ ,
						0, // const ARM_Date* FwdFixedDate /* = "NULL"*/ ,
						INCLUDE_MATURITY // const bool& IncludeMaturity /* = INCLUDE_MATURITY*/ );
					  );
}
// -----------------------------------------------------------------------------------
// Construction de l'objet ICM_Matrix qui va nous permettre de passer les infos concernant la méthode numérique.
ICM_Parameters* Matrix = new ICM_Parameters();

const char* tmp0 = "COPULA";
	ARM_Vector V0 = ARM_Vector(1,1);
	Matrix->Push((ARM_Vector*)V0.Clone(), (char*)tmp0);
	
	const char* tmp1 = "INTEGRATION_STEP_1";
	ARM_Vector V1 = ARM_Vector(1,NbStep);
	Matrix->Push((ARM_Vector*)V1.Clone(), (char*)tmp1);
	
	const char* tmp2 = "FREEDOM_DEGREE";
	ARM_Vector V2 (1,1);
	Matrix->Push((ARM_Vector*)V2.Clone(), (char*)tmp2);

	const char* tmp3 = "INTEGRATION_METHOD";
	// 0 -> qGAUSSLEGENDRE; // 1 -> qGAUSSHERMITE; // 2 -> qTRAPEZE
	ARM_Vector V3 (1,MethodNum );
	Matrix->Push((ARM_Vector*)V3.Clone(), (char*)tmp3);
// -----------------------------------------------------------------------------------
// Pricer
ICM_Pricer_Distrib_Smile  Pricer_Mez;Pricer_Mez.Set(mez,Model,*Matrix,AsOf);
// -----------------------------------------------------------------------------------
// Resu
	if (C_Function == "NPV")
		price = - Pricer_Mez.Price(qCMPPRICE);
	else if (C_Function == "DL")
		price = Pricer_Mez.Price(qCMPDEFLEGPV);
	else if (C_Function == "FL")
		price = Pricer_Mez.Price(qCMPFEELEGPV);
	else if (C_Function == "BE")
	{
		if (MezType != 4)
			price = Pricer_Mez.ComputeSpread(0)/10000.;
		else
			price = Pricer_Mez.ComputeSpread(0);
	}
	else if (C_Function == "DUR")
		price = Pricer_Mez.Price(qCMPDURATION);
	else if (C_Function == "SENSI") //  Sensi au spread
		price = - Pricer_Mez.Hedge(ICMSPREAD_TYPE,(char*)DateToShiftString.c_str(),CCSTringToSTLString(LabelToShift),ValBumpSensi);
	

	double Prix = price;
	if ( (fabs(Notional[1] - MezzAmount)>1E-6) && (fabs(MezzAmount)> 1e-10))
		if (!((C_Function == "BE") || (C_Function == "DUR")))
			Prix = (Notional[1]/MezzAmount) * price;

	result.setDouble(Prix);
// -----------------------------------------------------------------------------------		
	// DELETE
	// MyDelete(Pricer_Mez);;
	MyDelete(Matrix);
	MyDelete(Model);
	MyDelete(mez);
	MyDelete(Betas);
	MyDeleteTab(TabCrbDefaut);
	//MyDeleteTab(pIssuersNotionals);

	// if (Label)
	// 	FreePointerTabChar(Label,NbIssuers);

	MyDelete(NameSpreads);
	
	MyDeleteTab(CalibRecov);
//	MyDeleteTab(LossRecov);
	MyDelete(Index);
	MyDelete(ForcedCurve);
	MyDelete(VolIndex);
	MyDelete(CrbTaux);
	
	return ARM_OK;
	}
 
	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}

	catch (...)
	{
		result.setMsg ("ARM_ERR: Price : unrecognized failure");
		return ARM_KO;
	}
};


// Pricing de Tranche
// -----------------------------------------------------------------------------------
// CDO : Expected Loss
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
										 ARM_result& result)
{
	ICMMSG(WARN,"Using ICMLOCAL_STR_GetExpectedLoss" ); 
	CCString msg ("");

	int i,j = 0;
	double Taux = taux*100.;
	double Basis = 365.;
	ARM_Date AsOf; 
	AsOf.Today();
	AsOf.AddDays(-1);
	double AsOfExcelDate = AsOf.GetJulian() - 2415019;

	double YearTermYF = YearTerm;
	if (YearTerm>100)
		YearTerm = ( YearTerm - AsOfExcelDate)/Basis;

	ARM_Date EndDate;
	char* pEndDate =  new char[11];

	Local_XLDATE2ARMDATE(AsOfExcelDate + Basis*YearTermYF,pEndDate);
	EndDate = (ARM_Date) pEndDate;

	
	MyDeleteTab(pEndDate);
	
	double EL = 0.;
// ---------------------------------------------------------------------------------	
	try
	{
// ---------------------------------------------------------------------------------
// Issuers Data
	int NbIssuers = (int)Nb_Issuers;
// ---------------------------------------------------------------------------------
// Recovery	
	double* CalibRecov =  new double[NbIssuers];
	memset(CalibRecov,'\0',sizeof(double)*NbIssuers);

	if (CalibRecovery.size()<NbIssuers)
		for (i=0;i<NbIssuers;i++)
			CalibRecov[i] = CalibRecovery[0];
	else
		for (i=0;i<NbIssuers;i++)
			CalibRecov[i] = CalibRecovery[i];


	ARM_Vector LossRecov (NbIssuers);
	// memset(LossRecov,'\0',sizeof(double)*NbIssuers);

	if (LossRecovery.size()<NbIssuers)
		for (i=0;i<NbIssuers;i++)
			LossRecov[i] = LossRecovery[0];
	else
		for (i=0;i<NbIssuers;i++)
			LossRecov[i] = LossRecovery[i];
// ---------------------------------------------------------------------------------
// Matrice de Spread
	int nb_plot = YF_plot.size();
	int NbRow = NbIssuers;
	int NbCol = nb_plot; 
	ICM_QMatrix<double>* NameSpreads = new ICM_QMatrix<double>(NbRow,NbCol,CREDIT_DEFAULT_VALUE);
	for (i=0;i<NbRow;i++)
	{
		// Cas Vraie Matrice
		if (MatIssuersSpread.size()==NbCol*NbRow)
		{
			for (j=0;j<NbCol;j++)
				NameSpreads->SetValue(i,j,MatIssuersSpread[i*NbCol+j]);
		}
		// Cas d'un input = une ligne
		else if (MatIssuersSpread.size() == NbCol)
		{
			for (j=0;j<NbCol;j++)
				NameSpreads->SetValue(i,j,MatIssuersSpread[j]);
		}
		// Cas d'une colonne
		else if (MatIssuersSpread.size() == NbIssuers)
		{
			for (j=0;j<NbCol;j++)
				NameSpreads->SetValue(i,j,MatIssuersSpread[i]);
		}
		// Cas une seule valeur
		else if (MatIssuersSpread.size() == 1)
		{
			for (j=0;j<NbCol;j++)
				NameSpreads->SetValue(i,j,MatIssuersSpread[0]);
		}
		else
		{
			result.setMsg ("ARM_ERR: Issuers Spread : Size Error");
			return ARM_KO;
		}	

	}
// ---------------------------------------------------------------------------------
// Vecteur de plot
	ARM_Vector PlotNum (YF_plot);
	if (YF_plot[0]>100)
	{
		for (i=0;i<YF_plot.size();i++)
			PlotNum.Elt(i) = (YF_plot[i] - AsOfExcelDate)/Basis;
	}

	vector<string> Term;
	for (i=0;i<YF_plot.size();i++)
	{
		char * t = NULL;
		ITOA(PlotNum.Elt(i),t,10);
		string s = string(t) + string("Y");
		Term.push_back(s);
		if (t) delete t;
	}
// ---------------------------------------------------------------------------------
// Betas	
	int SizeBeta = IssuersBeta.size();
	ARM_Vector BetasIssuers = ARM_Vector(NbIssuers,1.);;
	if ( (SizeBeta<NbIssuers) && (SizeBeta == 1))
	{
		for (i=0;i<NbIssuers;i++)
			BetasIssuers.Elt(i) = IssuersBeta[0];
	}
	else if (SizeBeta >= NbIssuers)
	{
		for (i=0;i<NbIssuers;i++)
			BetasIssuers.Elt(i) = IssuersBeta[i];
	}
	else 
	{
		result.setMsg ("ARM_ERR: Vector Betas : Size Error");
		return ARM_KO;
	}	
// ---------------------------------------------------------------------------------
// Labels
	std::vector<std::string> Label(NbIssuers); 
	// char** Label = new char*[NbIssuers];
	for (int i=0;i<NbIssuers;i++)
	{
		std::stringstream sstr; sstr<<"ISSUER_"<<i; 
		sstr>>Label[i]; 
		// Label[i] = new char[60];
		// sprintf(Label[i],"ISSUER_%i",i);
	}
// ---------------------------------------------------------------------------------
// Vecteur d'Issuers Notionals
	// double* pIssuersNotionals = new double[NbIssuers];
	std::vector<double> pIssuersNotionals (NbIssuers); 
	for (int il=0; il<NbIssuers; il++) 
		pIssuersNotionals[il]=10000000.;
// ---------------------------------------------------------------------------------
// Points d'attachement 
	double CDOAmount  = NbIssuers * 10000000.;
	double MezzAmount = (Attachement[1] - Attachement[0])*CDOAmount;
	double SubAmount  = Attachement[0] * CDOAmount;
// ---------------------------------------------------------------------------------
// Methode NUMERIQUE
	// 0 - qGAUSSLEGENDRE; 	// 1 - qGAUSSHERMITE; 	// 2 - qTRAPEZE
	int NbStep = atoi((char*)NumMethod[0]);
	bool FastStrip = true;

	double MethodNum = 1.;
	if (NumMethod.size() > 1)
	{
		if ((NumMethod[1] == "l") || (NumMethod[1] == "L"))
			MethodNum = 0.;
		else if ((NumMethod[1] == "t") || (NumMethod[1] == "T"))
			MethodNum = 2.;

		if (NumMethod.size() == 3)
		{
			string Stripping = NumMethod[2].c_str();
			FastStrip = FastStripping(Stripping);
		}
	}
// ---------------------------------------------------------------------------------
// Courbe de Taux
	ARM_ZeroFlat* CrbTaux = new ARM_ZeroFlat(AsOf, Taux);
// ---------------------------------------------------------------------------------
// Courbes de défaut.
	std::vector<const ICM_DefaultCurve*> TabCrbDefaut (NbIssuers);
	
	if (FastStrip)
	{
		for (i = 0;i<NbIssuers;i++)
		{
			ARM_Vector* Spread = NameSpreads->TruncRowAsVector(i,nb_plot); 

			TabCrbDefaut[i] = new ICM_DefCurvStr(Spread,
												 &PlotNum,
												 CalibRecov[i],
												 CrbTaux,
												 Label[i]);
			// Delete
			MyDelete(Spread);
		}
	}
	else
	{
		for (i = 0;i<NbIssuers;i++)
		{
			ARM_Vector* Spread = NameSpreads->TruncRowAsVector(i,nb_plot); 

			TabCrbDefaut[i] = new ICM_Constant_Piecewise(AsOf,
													     Term,
														 (ARM_Vector*)Spread->Clone(),
														 CalibRecov[i],
														 CrbTaux,
												K_ADJUSTED,	// intRule
												K_ADJUSTED,	// adjStartDate
														 qCredit_Adjust20,
														 ARM_DEFAULT_COUNTRY,
														 Label[i],
														 false,
														 NULL,//2 NULL,
														 K_QUARTERLY,qDEFCURVE_DICHO,"STD", 
														 ARM_Currency(ARM_DEFAULT_COUNTRY).GetCreditStartDateLag()
														 ,ICM_Parameters());
			// Delete
			MyDelete(Spread);
		}
	}
// ---------------------------------------------------------------------------------
// Correlation
ICM_Beta_Correlation* Betas = new ICM_Beta_Correlation(AsOf,"STRCORR",BetasIssuers,Label,(ARM_IRIndex*)0,(ARM_IRIndex*)0);
// ---------------------------------------------------------------------------------
// Model
ICM_ModelMultiCurves* Model = new ICM_ModelMultiCurves(// NbIssuers,
													   TabCrbDefaut,
													   (ARM_ZeroCurve*) CrbTaux, 
													   LossRecov,
													   Betas);
// ---------------------------------------------------------------------------------
// Objet Mezz
ICM_Mez* mez = NULL;
double FixedRate = 0.01;
double FloatingPayerAmount = MezzAmount;
int stubrule = 1;

mez = new ICM_Mez(AsOf,
				  EndDate,
				 (ARM_Date*)0,
				 (ARM_Date*)0,
				  FixedRate,
				  K_ADJUSTED,	// intRule
				  K_ADJUSTED,	// adjStartDate
				  SubAmount,
				  MezzAmount,
				  Label,
				  pIssuersNotionals,
				  //NbIssuers,
				  K_QUARTERLY,
				  KACTUAL_360,
				  MezzAmount,
				  qACCRUED_SETTLED,
				  ARM_DEFAULT_COUNTRY,
				  FloatingPayerAmount,
				  stubrule,
				  0,
				  K_QUARTERLY,
				CREDIT_DEFAULT_VALUE, // const double& Binary /* = CREDIT_DEFAULT_VALUE*/,
				std::string(), // const std::string&  payCalName /* = NULL*/,
				qRunning_Leg, // const qCredit_Leg_Type& TypeFeeLeg/* =qRunning_Leg*/,
				qStandart_Recovery_Leg, // const qCredit_Leg_Type& TypeDefLeg/* =qStandart_Recovery_Leg*/,
				INCLUDE_MATURITY // const bool& IncludeMaturity /* = INCLUDE_MATURITY*/);

				  );
// ---------------------------------------------------------------------------------
// Construction de l'objet ICM_Matrix qui va nous permettre de passer les infos concernant la méthode numérique.
ICM_Parameters* Matrix = new ICM_Parameters();

const char* tmp0 = "COPULA";
	ARM_Vector V0 = ARM_Vector(1,1);
	Matrix->Push((ARM_Vector*)V0.Clone(), (char*)tmp0);
	
	const char* tmp1 = "INTEGRATION_STEP_1";
	ARM_Vector V1 = ARM_Vector(1,NbStep);
	Matrix->Push((ARM_Vector*)V1.Clone(), (char*)tmp1);
	
	const char* tmp2 = "FREEDOM_DEGREE";
	ARM_Vector V2 = ARM_Vector(1,1);
	Matrix->Push((ARM_Vector*)V2.Clone(), (char*)tmp2);

	const char* tmp3 = "INTEGRATION_METHOD";
	// 0 -> qGAUSSLEGENDRE;	// 1 -> qGAUSSHERMITE; // 2 -> qTRAPEZE
	ARM_Vector V3 = ARM_Vector(1,MethodNum);
	Matrix->Push((ARM_Vector*)V3.Clone(), (char*)tmp3);
// ---------------------------------------------------------------------------------
// Pricer
	ICM_Pricer_Distrib_Smile  Pricer_Mez ; Pricer_Mez.Set(mez,Model,*Matrix,AsOf);
// ---------------------------------------------------------------------------------
// Resu
	vector<double> losses;
	EL = Pricer_Mez.ExpectedLossTranche(YearTerm,losses);
	
	result.setDouble(EL);
// ---------------------------------------------------------------------------------	
	// DELETE
	// MyDelete(Pricer_Mez);
	MyDelete(Matrix);
	MyDelete(Model);
	MyDelete(mez);
	MyDelete(Betas);
	MyDeleteTab(TabCrbDefaut);
	//MyDeleteTab(pIssuersNotionals);	
	
	// if (Label)
	//	FreePointerTabChar(Label,NbIssuers);

	MyDelete(NameSpreads);	
	MyDeleteTab(CalibRecov);	
//	MyDeleteTab(LossRecov);	
	MyDelete(CrbTaux);	

	return ARM_OK;
	}
 
	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}

	catch (...)
	{
		result.setMsg ("ARM_ERR: Expected Loss : unrecognized failure");
		return ARM_KO;
	}
};

// -----------------------------------------------------------------------------------
// CDO Smile Shocked : Interface Proto / Code ARM Damien
/* 
extern long ICMLOCAL_STR_CDOSMILE (int MezType, 
								   vector<double> DataMez,
								   vector<double> Recovery,
								   vector<double> PlotPtf,
								   vector<double> MatIssuersSpread,								
								   vector<CCString> IssuersLabel,
								   vector<double> PlotTraxx,
								   vector<double> SpreadTraxx,
								   vector<double> MatBaseCorrel,
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
								   ARM_result& result)
{
	ICMMSG(WARN,"Using ICMLOCAL_STR_CDOSMILE" ); 
	CCString msg ("");

	try
	{

	int i =0;
	int j =0;

	double price = 0.;
	double Basis = 365.;
	// --------------------------------------------------------------------------
	// Message d'erreur
	if (DataMez.size() != 7)
	{
		result.setMsg ("ARM_ERR: Le vecteur de Data Mez est mal renseigné");
		return ARM_KO;
	}
	// --------------------------------------------------------------------------
	// Notional
	vector<double> Notional;
	Notional.resize(2);
	Notional[0] = 1000000.0;
	Notional[1] = DataMez[6];
	// --------------------------------------------------------------------------
	// Date
	double Taux  = taux*100.;
	double YF_Maturity  = DataMez[0];
	double FixedRate = DataMez[4];
	ARM_Date AsOf; 
	AsOf.Today();

	double AsOfExcelDate = AsOf.GetJulian() - 2415019;
	double Maturity = YF_Maturity;
	if (YF_Maturity>100)
		Maturity = ( YF_Maturity - AsOfExcelDate)/Basis;

	double AsOfJuDate = AsOf.GetJulian();

	ARM_Date EndDate;
	EndDate.Today();
	char* pEndDate =  new char[11];

	Local_XLDATE2ARMDATE(AsOfExcelDate + Basis*Maturity,pEndDate);
	EndDate = (ARM_Date) pEndDate;

	MyDeleteTab(pEndDate);
// --------------------------------------------------------------------------
   // Maturity Base Correl
	vector<double> MaturityBaseCorrel;
	MaturityBaseCorrel.resize(MatBaseCorrel.size());

	// On distingue le cas Year Fraction du cas Excel Date
	if (MatBaseCorrel[0]<100)
		for (i = 0; i<MatBaseCorrel.size();i++)
			MaturityBaseCorrel[i] = MatBaseCorrel[i];
	else    
		for (i = 0; i<MatBaseCorrel.size();i++)
			MaturityBaseCorrel[i] = (MatBaseCorrel[i]  - AsOfExcelDate)/Basis;
// ---------------------------------------------------------------------------
// Issuers Data
	int NbIssuers = (int)DataMez[5];
// ---------------------------------------------------------------------------	
// Recovery	
	double* CalibRecov =  new double[NbIssuers+2];
	memset(CalibRecov,'\0',sizeof(double)*(NbIssuers+2));

	ARM_Vector LossRecov (NbIssuers+2);
	// memset(LossRecov,'\0',sizeof(double)*NbIssuers+2);

	if (Recovery.size() == 1)
	{
		for (i=0;i<NbIssuers;i++)
		{
			CalibRecov[i] = Recovery[0];
			LossRecov[i] = Recovery[0];
		}
	}
	else if (Recovery.size() == 2)
	{
		for (i=0;i<NbIssuers;i++)
		{
			CalibRecov[i] = Recovery[0];
			LossRecov[i] = Recovery[1];
		}
	}
	else if (Recovery.size() < NbIssuers)
	{
		for (i=0;i<NbIssuers;i++)
		{
			CalibRecov[i] = Recovery[0];
			LossRecov[i] = Recovery[0];
		}
	}
	else if (Recovery.size() == NbIssuers)
	{
		for (i=0;i<NbIssuers;i++)
		{
			CalibRecov[i] = Recovery[i];
			LossRecov[i] = Recovery[i];
		}
	}
	else if (Recovery.size() == 2*NbIssuers)
	{
		for (i=0;i<NbIssuers;i++)
		{
			CalibRecov[i] = Recovery[2*i];
			LossRecov[i] = Recovery[2*i+1];
		}
	}
// ---------------------------------------------------------------------------
// Matrice de Spread
	int nb_plot = PlotPtf.size();
	int NbRow = NbIssuers;
	int NbCol = nb_plot; 
	ICM_QMatrix<double>* NameSpreads = new ICM_QMatrix<double>(NbRow,NbCol,CREDIT_DEFAULT_VALUE);
	for (i=0;i<NbRow;i++)
	{
		// Cas Vraie Matrice
		if (MatIssuersSpread.size()==NbCol*NbRow)
		{
			for (j=0;j<NbCol;j++)
				NameSpreads->SetValue(i,j,MatIssuersSpread[i*NbCol+j]);
		}
		// Cas d'un input = une ligne
		else if (MatIssuersSpread.size() == NbCol)
		{
			for (j=0;j<NbCol;j++)
				NameSpreads->SetValue(i,j,MatIssuersSpread[j]);
		}
		// Cas d'une colonne
		else if (MatIssuersSpread.size() == NbIssuers)
		{
			for (j=0;j<NbCol;j++)
				NameSpreads->SetValue(i,j,MatIssuersSpread[i]);
		}
		// Cas une seule valeur
		else if (MatIssuersSpread.size() == 1)
		{
			for (j=0;j<NbCol;j++)
				NameSpreads->SetValue(i,j,MatIssuersSpread[0]);
		}
		else
		{
			result.setMsg ("ARM_ERR: Issuers Spread : Size Error");
			return ARM_KO;
		}	
	}
// ---------------------------------------------------------------------------
// Vecteur de PlotPtf
	const int sizePlotPtf = PlotPtf.size();
	ARM_Vector PlotNumPtf = ARM_Vector(PlotPtf);
	if (PlotPtf[0]>100)
	{
		for (i=0;i<sizePlotPtf;i++)
			PlotNumPtf.Elt(i) = (PlotPtf[i] - AsOfExcelDate)/Basis;
	}

	vector<string> TermPtf;
	for (i=0;i<sizePlotPtf;i++)
	{
			char * t = new char[5];
			ITOA(PlotNumPtf.Elt(i),t,10);
			string s = string(t) + string("Y");
			TermPtf.push_back(s);
			if (t) delete t;
	}
// ---------------------------------------------------------------------------
// Vecteur de PlotTraxx
	const int sizePlotTraxx = PlotTraxx.size();
	ARM_Vector PlotNumTraxx = ARM_Vector(PlotTraxx);
	if (PlotTraxx[0]>100)
	{
		for (i=0;i<sizePlotTraxx;i++)
			PlotNumTraxx.Elt(i) = (PlotTraxx[i] - AsOfExcelDate)/Basis;
	}

	vector<string> TermTraxx;
	for (i=0;i<sizePlotTraxx;i++)
	{
			char * t = new char[5];
			ITOA(PlotNumTraxx.Elt(i),t,10);
			string s = string(t) + string("Y");
			TermTraxx.push_back(s);
			if (t) delete t;
	}

// ---------------------------------------------------------------------------
//  Date To Shift
	string DateToShiftString = "NONE";
	if (DateToShift != CREDIT_DEFAULT_VALUE)
	{
		if (DateToShift > 100)
		{
			DateToShift -= AsOfExcelDate;
			DateToShift *= 1/Basis;
		}
		_gcvt(DateToShift,NBDECIMAL,(char*)DateToShiftString.c_str());
		strcat((char*)DateToShiftString.c_str(),"Y");
	}
// ---------------------------------------------------------------------------
// ValShift
	ValShift *= 100.;
// ---------------------------------------------------------------------------
// Label
	if (IssuersLabel.size() < NbIssuers)
	{
		result.setMsg ("ARM_ERR: Le vecteur de Label n'est pas de la bonne taille");
		return ARM_KO;
	}
	std::vector<std::string> Label(NbIssuers);
	for (i =0; i<NbIssuers;i++)
		Label[i] = IssuersLabel[i];
// ---------------------------------------------------------------------------
// Frequency
	int Frequency = FindFrequency(DataMez[3],DataMez[0]);
// ---------------------------------------------------------------------------
// Vecteur d'Issuers Notionals
	// double* pIssuersNotionals = new double[NbIssuers];
	std::vector<double> pIssuersNotionals (NbIssuers); 
	for (int il=0; il<NbIssuers; il++) 
		pIssuersNotionals[il]=Notional[0];
// ---------------------------------------------------------------------------
// Points d'attachement 
	double CDOAmount = 0.;
	for (i=0;i<NbIssuers;i++)
		CDOAmount += pIssuersNotionals[i];

	double MezzAmount = (DataMez[2] - DataMez[1])*CDOAmount;
	double SubAmount  = DataMez[1] * CDOAmount;
// ---------------------------------------------------------------------------
// Traxx Data
	double MaturiteTraxx		= Maturity;
	double ProportionTraxxEuro	= DataTraxx[0];
	double RecovTraxx			= DataTraxx[1];
	LossRecov[NbIssuers]		= RecovTraxx;
	LossRecov[NbIssuers+1]		= RecovTraxx;

	CalibRecov[NbIssuers]		= RecovTraxx;
	CalibRecov[NbIssuers+1]		= RecovTraxx;
// ---------------------------------------------------------------------------
// Correlation
	int sizeMatBaseCorrel = MaturityBaseCorrel.size();

	if (nbRowCorrel == sizeMatBaseCorrel + 1)	// On duplique les informations de l'indice Euro
	{
		for (i=0;i<(sizeMatBaseCorrel + 1)*nbColCorrel;i++)
			CorrelParStrike.push_back(CorrelParStrike[i]);
		
		ProportionTraxxEuro = 1.;				// On force la projection sur le traxx euro
	}

	vector<double> Correl;
	Correl.resize(CorrelParStrike.size());
	for (i=0;i<CorrelParStrike.size();i++)
		Correl[i] = atof((char*)CorrelParStrike[i]);
// ---------------------------------------------------------------------------
// spread
	ARM_Vector SpreadTraxxEur = ARM_Vector(sizePlotTraxx,1.);
	ARM_Vector SpreadTraxxUS = ARM_Vector(sizePlotTraxx,1.);

	bool Rescal = true;
	if (SpreadTraxx[0] < 0.)
	{	
		Rescal = false;
		SpreadTraxx[0] = 0.01/100.;
	}

	if(SpreadTraxx.size() == 1)
	{
		for (i=0;i<sizePlotTraxx;i++)
		{
			SpreadTraxxEur.Elt(i) = SpreadTraxx[0];
			SpreadTraxxUS.Elt(i) = SpreadTraxx[0];
		}
	}
	else if(SpreadTraxx.size() == 2)
	{
		
		for (i=0;i<sizePlotTraxx;i++)
		{
			SpreadTraxxEur.Elt(i) = SpreadTraxx[0];
			SpreadTraxxUS.Elt(i) = SpreadTraxx[1];
		}
	}
	else if(SpreadTraxx.size() == sizePlotTraxx)
	{
		for (i=0;i<sizePlotTraxx;i++)
		{
			SpreadTraxxEur.Elt(i) = SpreadTraxx[i];
			SpreadTraxxUS.Elt(i) = SpreadTraxx[i];	 
		}
	}
	else if(SpreadTraxx.size() == 2*sizePlotTraxx) 
	{	
		for (i=0;i<sizePlotTraxx;i++)
		{
			SpreadTraxxEur.Elt(i) = SpreadTraxx[i];
			SpreadTraxxUS.Elt(i) = SpreadTraxx[sizePlotTraxx+i];	 
		}
	}
	else
	{
		result.setMsg ("ARM_ERR: Vector Spread Traxx : Size Error");
		return ARM_KO;
	}
// ---------------------------------------------------------------------------
// Methode NUMERIQUE
	// 0 - qGAUSSLEGENDRE; 	// 1 - qGAUSSHERMITE; 	// 2 - qTRAPEZE
	int NbStep = atoi((char*)NumMethod[0]);
	bool FastStrip = true;

	double MethodNum = 1.;
	if (NumMethod.size() > 1)
	{
		NumMethod[1].toUpper();
		if ((NumMethod[1] == "l") || (NumMethod[1] == "L"))
			MethodNum = 0.;
		else if ((NumMethod[1] == "t") || (NumMethod[1] == "T"))
			MethodNum = 2.;

		if (NumMethod.size() == 3) 
		{
			string stripping = NumMethod[2].c_str();
			FastStrip = FastStripping(stripping);
		}
	}
// ***************************************************************************
// ***************************************************************************
//						CONSTRUCTION DES OBJETS

// Courbe de Taux
	ARM_ZeroFlat* CrbTaux = new ARM_ZeroFlat(AsOf, Taux);
// ---------------------------------------------------------------------------
// Courbes de défaut.
	ICM_DefaultCurve** TabCrbDefaut = new ICM_DefaultCurve*[NbIssuers+2];  // Le rajout de 2 est lié aux courbes d'indices
	
	
	if (FastStrip)
	{
		for (i = 0;i<NbIssuers;i++)
		{
			ARM_Vector* Spread = NameSpreads->TruncRowAsVector(i,nb_plot) ;
			
			TabCrbDefaut[i] = new ICM_DefCurvStr(Spread,
												 &PlotNumPtf,
												 CalibRecov[i],
												 (ARM_ZeroCurve*)CrbTaux,
												 Label[i]);
			// Delete
			MyDelete(Spread);
		}

		TabCrbDefaut[NbIssuers] = new ICM_DefCurvStr(&SpreadTraxxEur,
													 &PlotNumTraxx,
													 RecovTraxx,
													 (ARM_ZeroCurve*)CrbTaux,
													 "TRAXX_EUR");

		TabCrbDefaut[NbIssuers+1] = new ICM_DefCurvStr(&SpreadTraxxUS,
													   &PlotNumTraxx,
														RecovTraxx,
														(ARM_ZeroCurve*)CrbTaux,
														"TRAXX_US");
	}
	else
	{
		for (i = 0;i<NbIssuers;i++)
		{
			ARM_Vector* Spread = NameSpreads->TruncRowAsVector(i,nb_plot) ;
		
			TabCrbDefaut[i] = new ICM_Constant_Piecewise(AsOf,
														 TermPtf,
														 (ARM_Vector*)Spread->Clone(),
														 CalibRecov[i],
														 (ARM_ZeroCurve*)CrbTaux,
												K_ADJUSTED,	// intRule
												K_ADJUSTED,	// adjStartDate
														 qCredit_Adjust20,
														 ARM_DEFAULT_COUNTRY,
														 Label[i],
														 false,
														 NULL,//2 NULL,
														 K_QUARTERLY,qDEFCURVE_DICHO,"STD", 
														 ARM_Currency(ARM_DEFAULT_COUNTRY).GetCreditStartDateLag()
														 ,ICM_Parameters());
			// Delete
			MyDelete(Spread);
		}

		TabCrbDefaut[NbIssuers] = new ICM_Constant_Piecewise(AsOf,
															 TermTraxx,
															 (ARM_Vector*)SpreadTraxxEur.Clone(),
															 RecovTraxx,
															 (ARM_ZeroCurve*)CrbTaux,
												K_ADJUSTED,	// intRule
												K_ADJUSTED,	// adjStartDate
															 qCredit_Adjust20,
															 ARM_DEFAULT_COUNTRY,
															 "TRAXX_EUR",
															 false,
															 NULL,//2 NULL,
															 K_QUARTERLY,qDEFCURVE_DICHO,"STD",
															 ARM_Currency(ARM_DEFAULT_COUNTRY).GetCreditStartDateLag()
															 ,ICM_Parameters());

		
		TabCrbDefaut[NbIssuers+1] = new ICM_Constant_Piecewise(AsOf,
																 TermTraxx,
																 (ARM_Vector*)SpreadTraxxUS.Clone(),
																 RecovTraxx,
																 (ARM_ZeroCurve*)CrbTaux,
												K_ADJUSTED,	// intRule
												K_ADJUSTED,	// adjStartDate
																 qCredit_Adjust20,
																 ARM_DEFAULT_COUNTRY,
																 "TRAXX_US",
																 false,
																 NULL,//2 NULL,
																 K_QUARTERLY,qDEFCURVE_DICHO,"STD",
																 ARM_Currency(ARM_DEFAULT_COUNTRY).GetCreditStartDateLag()
																 ,ICM_Parameters());
	}
// ---------------------------------------------------------------------------
// Correlation
	ARM_Vector MaturityTraxx = ARM_Vector(MaturityBaseCorrel);
	int NbNomTraxx = 125;
	bool JustEuro = false;
	bool JustUS = false;
	ARM_Vector proportion = ARM_Vector(2,1);
	proportion.Elt(0) = ProportionTraxxEuro;		// Proportion Traxx EUR
	proportion.Elt(1) = 1-ProportionTraxxEuro;	// Proportion Traxx US
	if (fabs(ProportionTraxxEuro - 1)<1e-6)
		JustEuro = true;
	else if (fabs(ProportionTraxxEuro)<1e-6)
		JustUS = true;

	std::vector<std::string> LabelTraxx (2); 
	LabelTraxx[0] = "TRAXX_EUR";
	LabelTraxx[1] = "TRAXX_US";

	std::vector<std::string> LabelTraxxEur (NbNomTraxx); 
	// char** LabelTraxxEur = new char*[NbNomTraxx];
	for (i=0;i<NbNomTraxx;i++)
		LabelTraxxEur[i] = "TRAXX_EUR";
		
	
	std::vector<std::string> LabelTraxxUS (NbNomTraxx); 
	// char** LabelTraxxUS = new char*[NbNomTraxx];
	for (i=0;i<NbNomTraxx;i++)
		LabelTraxxUS[i] = "TRAXX_US";
	
// Vecteur de strike
	ARM_Vector StrikEUR = ARM_Vector(nbColCorrel,0.);
	ARM_Vector StrikUS  = ARM_Vector(nbColCorrel,0.);
	
	for (i=0;i<nbColCorrel;i++)
	{
		StrikEUR.Elt(i) = 100 * Correl[i];
		StrikUS.Elt(i) = 100 * Correl[(sizeMatBaseCorrel + 1)*nbColCorrel+ i];
	}

// On recupere les correls
	ARM_Vector CorrelEUR = ARM_Vector(sizeMatBaseCorrel*nbColCorrel,0.);
	ARM_Vector CorrelUS  = ARM_Vector(sizeMatBaseCorrel*nbColCorrel,0.);
	for (i=0;i<sizeMatBaseCorrel*nbColCorrel;i++)
	{
		CorrelEUR.Elt(i) = 100 * Correl[i + nbColCorrel];
		CorrelUS.Elt(i) = 100 * Correl[i + (2+sizeMatBaseCorrel)*nbColCorrel];
	}

	ARM_Matrix* CorrelationEur = new ARM_Matrix(MaturityTraxx.GetSize(),nbColCorrel);
	ARM_Matrix* CorrelationUS = new ARM_Matrix(MaturityTraxx.GetSize(),nbColCorrel);
	for (j = 0; j<sizeMatBaseCorrel;j++) 
	{
		for (il=0;il<nbColCorrel;il++)
		{
			CorrelationEur->Elt(j,il) = CorrelEUR.Elt(il+j*nbColCorrel);
			CorrelationUS->Elt(j,il) = CorrelUS.Elt(il+j*nbColCorrel);
		}
	}
	
// Objet Vol Curve : on doit en faire un pour chaque Traxx
	vector<const ARM_VolCurve*> Vector_VolCurve(2);
	Vector_VolCurve[0] = new ARM_VolLInterpol(AsOf,
											  (ARM_Vector*)MaturityTraxx.Clone(),
											  (ARM_Vector*)StrikEUR.Clone(),				// Strike en bp
											  (ARM_Matrix*)CorrelationEur->Clone(),			// Valeur de la correl en %
											  1,
	 										  K_SMILE_VOL);
	
	Vector_VolCurve[1] = new ARM_VolLInterpol(AsOf,
											  (ARM_Vector*)MaturityTraxx.Clone(),
											  (ARM_Vector*)StrikUS.Clone(),				// Strike en bp
											  (ARM_Matrix*)CorrelationUS->Clone(),		// Valeur de la correl en %
											  1,
	 										  K_SMILE_VOL);

// On cree un index pour chaque traxx
	vector<const ICM_Credit_Index*> VIndex(2);
	ARM_Vector YT(1); YT.Elt(0) = 5.0;
	ARM_Vector spread(1); spread.Elt(0) = 0.;
					VIndex[0] = new ICM_Credit_Index(KACTUAL_360, 
								  					  4, 
													  4, 
													  YT,
													  ARM_DEFAULT_COUNTRY,
													  LabelTraxxEur[0],
													  LabelTraxxEur,
													  qAVERAGE,spread,NULL,0,0,K_ARREARS,0,K_ARREARS,0,qCredit_Adjust20,5,2);

	// MyDeleteTab(LabelTraxxEur);


					VIndex[1]  = new ICM_Credit_Index(KACTUAL_360, 
													  4, 
													  4, 
													  YT,
													  ARM_DEFAULT_COUNTRY,
													  LabelTraxxUS[0],
													  LabelTraxxUS,
													  qAVERAGE,spread,NULL,0,0,K_ARREARS,0,K_ARREARS,0,qCredit_Adjust20,5,2);

	//MyDeleteTab(LabelTraxxUS);

// Correlation
ICM_Smile_Correlation* Corr = NULL;
ARM_Vector Vsmilestrikelow(2,DataMez[1]);			// On initialise avec les strike down
ARM_Vector Vsmilestrikehigh(2,DataMez[2]);			// On initialise avec les strike Up

if (JustEuro)
{
	Vsmilestrikelow.Elt(1) = 0.;
	Vsmilestrikehigh.Elt(1) = 0.;
}
else if (JustUS)
{
	Vsmilestrikelow.Elt(0) = 0.;
	Vsmilestrikehigh.Elt(0) = 0.;
}

	if (Rescal)
// FIXMEFRED: mig.vc8 (28/05/2007 14:21:52):cast
		Corr = new ICM_Smile_Correlation(AsOf,
										"CORRSTR",
										&(*Vector_VolCurve.begin()),
										LabelTraxx,
										&proportion,
										&(*VIndex.begin()));
	else 
		Corr = new ICM_Smile_Correlation(AsOf,
										 "CORRSTR",
										 &(*Vector_VolCurve.begin()), 
										 LabelTraxx,
										 &proportion,
										 &Vsmilestrikelow,
										 &Vsmilestrikehigh);

	// MyDelete(Vector_VolCurve,2);
	delete Vector_VolCurve[0]; delete Vector_VolCurve [1];
// ---------------------------------------------------------------------------
// Objet Vol
ARM_VolFlat* VolIndex = new ARM_VolFlat(AsOf,VolAjustConvexity*100.);	// Multiplication pour être ARM_Compliant. On passe des objets en bp.
// ---------------------------------------------------------------------------

// ---------------------------------------------------------------------------
// Model
ICM_ModelMultiCurves* Model = new ICM_ModelMultiCurves(NbIssuers+2,
													   TabCrbDefaut,
													   (ARM_ZeroCurve*) CrbTaux->Clone(), 
													   LossRecov,
													   Corr,VolIndex);
// Objet Index
ICM_Credit_Index* Index = NULL;
ICM_DefaultCurve* ForcedCurve = NULL;
ICM_DefCurvStr* CourbeTemp = NULL;
ARM_Vector SpreadFwd = ARM_Vector(SpreadCrbFwd);
double yearterm = 5.;	// Maturité du Forward
if (MezType == 4)		// Cas CMTranche
{
	if (SpreadCrbFwd[0] != CREDIT_DEFAULT_VALUE)
	{
		CourbeTemp = new ICM_DefCurvStr(&SpreadFwd,&PlotNumPtf,CalibRecov[0],CrbTaux,"CrbFwd");
		ForcedCurve = (ICM_DefaultCurve*) (CourbeTemp->ConvertStrToARM())->Clone();
		MyDelete(CourbeTemp);
	}
	else 
	{	
		if (FastStrip)
			ForcedCurve = (ICM_DefaultCurve*) AverageCurveStr(TabCrbDefaut,NbIssuers);
		else
			ForcedCurve = (ICM_DefaultCurve*) AverageCurvePWC(TabCrbDefaut,NbIssuers);
	}
	  
	YT.Elt(0) = yearterm;
	Index = new ICM_Credit_Index(KACTUAL_360, 
								 Frequency, 
								 Frequency, 
								 YT,
								 ARM_DEFAULT_COUNTRY,
								 Label[0],
								 Label,
								 qAVERAGE,
								 spread,
								 // NULL,
								 ForcedCurve,
								 //VolIndex,
								 0,0,K_ARREARS,0,K_ARREARS,0,qCredit_Adjust20,5,2);
}
// ---------------------------------------------------------------------------
// Objet Mezz
ICM_Mez* mez = NULL;
double FloatingPayerAmount = MezzAmount;
int stubrule = 1;
int FrequencyDL = Frequency;

if (MezType != 4)
{
	if (MezType == 2) // Cas Infine
			Frequency = K_ZEROCOUPON;

	mez = new ICM_Mez(AsOf,
					  EndDate,
					  (ARM_Date*)0,
					  (ARM_Date*)0,
					  FixedRate,
					  K_ADJUSTED,	// intRule
					  K_ADJUSTED,	// adjStartDate
					  SubAmount,
					  MezzAmount,
					  Label,
					  pIssuersNotionals,
					  //NbIssuers,
					  Frequency,
					  KACTUAL_360,
					  MezzAmount,
					  qACCRUED_SETTLED,
					  ARM_DEFAULT_COUNTRY,
					  FloatingPayerAmount,
					  stubrule,
					  0,
					  FrequencyDL,
				CREDIT_DEFAULT_VALUE, // const double& Binary  ,
				std::string(), // const std::string&  payCalName  
				qRunning_Leg, // const qCredit_Leg_Type& TypeFeeLeg 
				qStandart_Recovery_Leg, // const qCredit_Leg_Type& TypeDefLeg 
				INCLUDE_MATURITY // const bool& IncludeMaturity  
					  );

	mez->SetTradedCoef(Notional[1]/MezzAmount);
}
else // On est dans le cas CMTranche
{
	mez = new ICM_Mez(AsOf,
					  EndDate,
					  (ARM_Date*)0,
					  (ARM_Date*)0,
					  K_ADJUSTED,	// intRule
					  K_ADJUSTED,	// adjStartDate
					  SubAmount,
					  MezzAmount,
					  Label,
					  pIssuersNotionals,
					  //NbIssuers,
				      Index, 
					  FixedRate,
					  Frequency,
					  KACTUAL_360,
					  MezzAmount,
					  qACCRUED_SETTLED,
					  ARM_DEFAULT_COUNTRY,
					  FloatingPayerAmount,
					  stubrule,
					  0,
					  FrequencyDL,
					CREDIT_DEFAULT_VALUE, // const double& Binary   ,
					std::string(), // const std::string&  
					0, // const ARM_Date* FwdFixedDate  
					INCLUDE_MATURITY // const bool& IncludeMaturity  
					  );

	mez->SetTradedCoef(Notional[1]/MezzAmount);
}
// ---------------------------------------------------------------------------
// Construction de l'objet ICM_Matrix qui va nous permettre de passer les infos concernant la méthode numérique.
ICM_Parameters* Matrix = new ICM_Parameters();

const char* tmp0 = "COPULA";
	ARM_Vector V0 = ARM_Vector(1,1);
	Matrix->Push((ARM_Vector*)V0.Clone(), (char*)tmp0);
	
	const char* tmp1 = "INTEGRATION_STEP_1";
	ARM_Vector V1 = ARM_Vector(1,NbStep);
	Matrix->Push((ARM_Vector*)V1.Clone(), (char*)tmp1);
	
	const char* tmp2 = "FREEDOM_DEGREE";
	ARM_Vector V2 = ARM_Vector(1,1);
	Matrix->Push((ARM_Vector*)V2.Clone(), (char*)tmp2);

	const char* tmp3 = "INTEGRATION_METHOD";
	// 0 -> qGAUSSLEGENDRE; // 1 -> qGAUSSHERMITE;	// 2 -> qTRAPEZE
	ARM_Vector V3 = ARM_Vector(1,MethodNum );
	Matrix->Push((ARM_Vector*)V3.Clone(), (char*)tmp3);
// ---------------------------------------------------------------------------
// Pricer
ICM_Pricer_Distrib_Smile Pricer_Mez ; Pricer_Mez.Set((ICM_Mez*)mez->Clone(),
																	Model,
																	*Matrix,AsOf);
// ---------------------------------------------------------------------------
// Resultat
FILE *stream1 = fopen("c:\\temp\\Object.txt", "w+");
mez->View("",stream1);
fclose(stream1);

if (C_Function == "NPV")
	{
		Resultat.resize(1);
		Resultat[0] = - Pricer_Mez.Price(qCMPPRICE); 
	}
	else if (C_Function == "DL")
	{
		Resultat.resize(1);
		Resultat[0] = Pricer_Mez.Price(qCMPDEFLEGPV);
	}
	else if (C_Function == "FL")
	{
		Resultat.resize(1);
		Resultat[0] = Pricer_Mez.Price(qCMPFEELEGPV);
	}
	else if (C_Function == "BE")
	{
		Resultat.resize(1);
		if (MezType != 4)
			Resultat[0] = Pricer_Mez.ComputeSpread(0)/10000.;
		else
			Resultat[0] = Pricer_Mez.ComputeSpread(0);
	}
	else if (C_Function == "DUR")
	{
		Resultat.resize(1);
		Resultat[0] =Pricer_Mez.Price(qCMPDURATION);
	}
	else if (C_Function == "SENSI") //  Sensi au spread
	{
		Resultat.resize(1);
		Resultat[0] = - Pricer_Mez.Hedge(ICMSPREAD_TYPE,DateToShiftString,CCSTringToSTLString(LabelToShift),ValShift); 
	}
 	else if (C_Function == "ALL")
	{
		Resultat.resize(13);

		if (MezType != 4)
			Resultat[3] =Pricer_Mez.ComputeSpread(0)/10000.;	// BE
		else
			Resultat[3] =Pricer_Mez.ComputeSpread(0);			// Part Rate

		FILE *stream2 = fopen("c:\\temp\\ObjectBE.txt", "w+");
		mez->View("",stream2);
		fclose(stream2);

		Resultat[0] = Pricer_Mez.Price(qCMPDEFLEGPV);						// DL

		FILE *stream3 = fopen("c:\\temp\\ObjectDL.txt", "w+");
		mez->View("",stream3);
		fclose(stream3);

		Resultat[1] = Pricer_Mez.Price(qCMPFEELEGPV);						// FL

		FILE *stream4 = fopen("c:\\temp\\ObjectFL.txt", "w+");
		mez->View("",stream4);
		fclose(stream4);

		Resultat[2] = - Pricer_Mez.Price(qCMPPRICE);			// NPV

		FILE *stream5 = fopen("c:\\temp\\ObjectNPV.txt", "w+");
		mez->View("",stream5);
		fclose(stream5);
		
		Resultat[4] =Pricer_Mez.Price(qCMPDURATION);					// DUR

		FILE *stream6 = fopen("c:\\temp\\ObjectDUR.txt", "w+");
		mez->View("",stream6);
		fclose(stream6);

		ICM_Smile_Correlation* CorrelSmile = (ICM_Smile_Correlation*) Model->GetCorrelation();

		int nbIndex = CorrelSmile->GetSlices().size();
		int NbMatUsed = CorrelSmile->GetSlices()[0].GetSmileStrikeHigh().size();
		
		// Il faut remplir les coordonnées 5 - 13
		for (int ind = 0; ind<nbIndex;ind++)
		{
			for (int ind2 = 0; ind2<NbMatUsed;ind2++)
			{
				Resultat[5 + 2*ind + 4*ind2] = CorrelSmile->GetSlices()[ind].GetSmileStrikeHigh()[ind2];
				Resultat[6 + 2*ind + 4*ind2] = CorrelSmile->GetSlices()[ind].GetSmileStrikeLow()[ind2];
			}
		}

	}
	
	result.setDouble(0.);
// ---------------------------------------------------------------------------		
	// DELETE
	// MyDelete(Pricer_Mez);
	MyDelete(Matrix);
	MyDelete(mez);
	MyDelete(Index);
	MyDelete(ForcedCurve);
	MyDelete(VolIndex);
	MyDelete(Model);
	MyDelete(Corr);
	delete VIndex[0]; delete VIndex[1];
	MyDelete(CorrelationUS);
	MyDelete(CorrelationEur);
	//MyDeleteTab(LabelTraxx);
	MyDelete(NameSpreads);
	MyDelete(TabCrbDefaut,NbIssuers + 2);
	MyDelete(CrbTaux);
	// if (Label)
	// 	FreePointerTabChar(Label,NbIssuers);
	//MyDeleteTab(pIssuersNotionals);
	MyDelete(NameSpreads);
	MyDeleteTab(CalibRecov);
//	MyDeleteTab(LossRecov);
	MyDeleteTab(pEndDate);

	return ARM_OK;
	}
 

	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}

	catch (...)
	{
		result.setMsg ("ARM_ERR: Price : unrecognized failure");
		return ARM_KO;
	}
};

**/ 
// -----------------------------------------------------------------------------------
// CDO Smile: Interface Proto / Code ARM Damien
extern long ICMLOCAL_STR_ExpectedLossSMILE( vector<double> DataMez,
								 		    vector<double> CalibRecovery,
										    vector<double> LossRecovery,
											vector<double> YF_plot,
											vector<double> MatIssuersSpread,								
											vector<double> SpreadTraxx,
											vector<double> MatBaseCorrel,
											vector<CCString>& CorrelParStrike,
									   	    int nbRowCorrel,
										    int nbColCorrel,
											vector<double> DataTraxx,
								    		vector<CCString> NumMethod,
								    		vector<double>& Resultat,
											ARM_result& result)
{
	ICMMSG(WARN,"Using ICMLOCAL_STR_ExpectedLossSMILE" ); 
	CCString msg ("");

	int i,j = 0;
	double price = 0.;
	

	try
	{
	double Basis = 365.;
	// --------------------------------------------------------------------------------------
	// Message d'erreur
	if (DataMez.size() != 6)
	{
		result.setMsg ("ARM_ERR: Le vecteur de Data Mez est mal renseigné");
		return ARM_KO;
	}
	// --------------------------------------------------------------------------------------

	double Taux  = DataMez[4]*100.;
	double YF_Maturity  = DataMez[0];
	double FixedRate = 1.;
	ARM_Date AsOf; 
	AsOf.Today();
	double AsOfExcelDate = AsOf.GetJulian() - 2415019;
	double Maturity = YF_Maturity;
	if (YF_Maturity>100)
		Maturity = ( YF_Maturity - AsOfExcelDate)/Basis;

	ARM_Date EndDate;
	char* pEndDate =  new char[11];
	Local_XLDATE2ARMDATE(AsOf.GetJulian() - 2415019 + Basis*Maturity,pEndDate);
	EndDate = (ARM_Date) pEndDate;
	
	MyDeleteTab(pEndDate);
// ---------------------------------------------------------------------------
// Issuers Data
	int NbIssuers = (int)DataMez[5];
// ---------------------------------------------------------------------------	
// Recovery	
	double* CalibRecov =  new double[NbIssuers+2];
	memset(CalibRecov,'\0',sizeof(double)*(NbIssuers+2));

	if (CalibRecovery.size()<NbIssuers)
		for (i=0;i<NbIssuers;i++)
			CalibRecov[i] = CalibRecovery[0];
	else
		for (i=0;i<NbIssuers;i++)
			CalibRecov[i] = CalibRecovery[i];


	ARM_Vector LossRecov (NbIssuers+2);
	// memset(LossRecov,'\0',sizeof(double)*NbIssuers+2);

	if (LossRecovery.size()<NbIssuers)
		for (i=0;i<NbIssuers;i++)
			LossRecov[i] = LossRecovery[0];
	else
		for (i=0;i<NbIssuers;i++)
			LossRecov[i] = LossRecovery[i];
// ---------------------------------------------------------------------------
// Matrice de Spread
	int nb_plot = YF_plot.size();
	int NbRow = NbIssuers;
	int NbCol = nb_plot; 
	ICM_QMatrix<double>* NameSpreads = new ICM_QMatrix<double>(NbRow,NbCol,CREDIT_DEFAULT_VALUE);
	for (i=0;i<NbRow;i++)
	{
		// Cas Vraie Matrice
		if (MatIssuersSpread.size()==NbCol*NbRow)
		{
			for (j=0;j<NbCol;j++)
				NameSpreads->SetValue(i,j,MatIssuersSpread[i*NbCol+j]);
		}
		// Cas d'un input = une ligne
		else if (MatIssuersSpread.size() == NbCol)
		{
			for (j=0;j<NbCol;j++)
				NameSpreads->SetValue(i,j,MatIssuersSpread[j]);
		}
		// Cas d'une colonne
		else if (MatIssuersSpread.size() == NbIssuers)
		{
			for (j=0;j<NbCol;j++)
				NameSpreads->SetValue(i,j,MatIssuersSpread[i]);
		}
		// Cas une seule valeur
		else if (MatIssuersSpread.size() == 1)
		{
			for (j=0;j<NbCol;j++)
				NameSpreads->SetValue(i,j,MatIssuersSpread[0]);
		}
		else
		{
			result.setMsg ("ARM_ERR: Issuers Spread : Size Error");
			return ARM_KO;
		}	
	}
// ---------------------------------------------------------------------------
// Vecteur de plot
	ARM_Vector PlotNum = ARM_Vector(YF_plot);
	const int sizePlot = YF_plot.size();
	if (YF_plot[0]>100)
	{
		for (i=0;i<sizePlot;i++)
			PlotNum.Elt(i) = (YF_plot[i] - AsOfExcelDate)/Basis;
	}

	// Term

	vector<string> Term;
	for (i=0;i<sizePlot;i++)
	{
		char * t = new char[5];
		ITOA(PlotNum.Elt(i),t,10);
		string s = string(t) + string("Y");
		Term.push_back(s);
		if (t) delete t;
	}
// ---------------------------------------------------------------------------
// Frequency
	int Frequency = FindFrequency(DataMez[3],DataMez[0]);
// ---------------------------------------------------------------------------
// Vecteur d'Issuers Notionals
	// double* pIssuersNotionals = new double[NbIssuers];
	std::vector<double> IssuersNotionals  (NbIssuers,10000000); 
	// for (int il=0; il<NbIssuers; il++) 
	// 	pIssuersNotionals[il]=10000000.;
// ---------------------------------------------------------------------------
// Points d'attachement 
	double CDOAmount  = NbIssuers * 10000000.;
	double MezzAmount = (DataMez[2] - DataMez[1])*CDOAmount;
	double SubAmount  = DataMez[1] * CDOAmount;
// ---------------------------------------------------------------------------
// Traxx Data
	double MaturiteTraxx		= 5.;
	double ProportionTraxxEuro	= DataTraxx[0];
	double RecovTraxx			= DataTraxx[1];
	LossRecov[NbIssuers] = RecovTraxx;
	LossRecov[NbIssuers+1] = RecovTraxx;

	CalibRecov[NbIssuers] = RecovTraxx;
	CalibRecov[NbIssuers+1] = RecovTraxx;
// --------------------------------------------------------------------------
   // Maturity Base Correl
	vector<double> MaturityBaseCorrel;
	MaturityBaseCorrel.resize(MatBaseCorrel.size());

	// On distingue le cas Year Fraction du cas Excel Date
	if (MatBaseCorrel[0]<100)
		for (i = 0; i<MatBaseCorrel.size();i++)
			MaturityBaseCorrel[i] = MatBaseCorrel[i];
	else    
		for (i = 0; i<MatBaseCorrel.size();i++)
			MaturityBaseCorrel[i] = (MatBaseCorrel[i]  - AsOfExcelDate)/Basis;
// ---------------------------------------------------------------------------
// Correlation
	int sizeMatBaseCorrel = MaturityBaseCorrel.size();

	if (nbRowCorrel == sizeMatBaseCorrel + 1)	// On duplique les informations de l'indice Euro
	{
		for (i=0;i<(sizeMatBaseCorrel + 1)*nbColCorrel;i++)
			CorrelParStrike.push_back(CorrelParStrike[i]);
		
		ProportionTraxxEuro = 1.;				// On force la projection sur le traxx euro
	}

	vector<double> Correl;
	Correl.resize(CorrelParStrike.size());
	for (i=0;i<CorrelParStrike.size();i++)
		Correl[i] = atof((char*)CorrelParStrike[i]);
// ---------------------------------------------------------------------------
// spread
	ARM_Vector SpreadTraxxEur = ARM_Vector(sizePlot,1.);
	ARM_Vector SpreadTraxxUS =  ARM_Vector(sizePlot,1.);


	bool Rescal = true;
	if (SpreadTraxx[0] < 0.)
	{	
		Rescal = false;
		SpreadTraxx[0] = 0.01/100.;
	}


	if(SpreadTraxx.size() == 1)
	{
		for (i=0;i<sizePlot;i++)
		{
			SpreadTraxxEur.Elt(i) = SpreadTraxx[0];
			SpreadTraxxUS.Elt(i) = SpreadTraxx[0];
		}

	}
	else if(SpreadTraxx.size() == 2)
	{
		
		for (i=0;i<sizePlot;i++)
			SpreadTraxxEur.Elt(i) = SpreadTraxx[0];
		
		for (i=0;i<sizePlot;i++)
			SpreadTraxxUS.Elt(i) = SpreadTraxx[1];
	}
	else if(SpreadTraxx.size() == sizePlot)
	{
		for (i=0;i<sizePlot;i++)
		{
			SpreadTraxxEur.Elt(i) = SpreadTraxx[i];
			SpreadTraxxUS.Elt(i) = SpreadTraxx[i];	 
		}
	}
	else if(SpreadTraxx.size() == 2*sizePlot) 
	{	
		for (i=0;i<sizePlot;i++)
		{
			SpreadTraxxEur.Elt(i) = SpreadTraxx[i];
			SpreadTraxxUS.Elt(i) = SpreadTraxx[sizePlot+i];	 
		}
	}
	else
	{
		result.setMsg ("ARM_ERR: Vector Spread Traxx : Size Error");
		return ARM_KO;
	}
// ---------------------------------------------------------------------------
// Labels
	std::vector<std::string> Label (NbIssuers);
	for (int i=0;i<NbIssuers;i++)
	{
		std::stringstream sstr; sstr<<"ISSUER_"<<i; 
		sstr>>Label[i] ;
	}
// ---------------------------------------------------------------------------
// Methode NUMERIQUE
	// 0 - qGAUSSLEGENDRE; 	// 1 - qGAUSSHERMITE; 	// 2 - qTRAPEZE
	int NbStep = atoi((char*)NumMethod[0]);
	bool FastStrip = true;

	double MethodNum = 1.;
	if (NumMethod.size() > 1)
	{
		if ((NumMethod[1] == "l") || (NumMethod[1] == "L"))
			MethodNum = 0.;
		else if ((NumMethod[1] == "t") || (NumMethod[1] == "T"))
			MethodNum = 2.;

		if (NumMethod.size() == 3)
		{
			string stripping = NumMethod[2].c_str();;
			FastStrip = FastStripping(stripping);
		}
	}
// ***************************************************************************
//						CONSTRUCTION DES OBJETS

// Courbe de Taux
	ARM_ZeroFlat* CrbTaux = new ARM_ZeroFlat(AsOf, Taux);
// ---------------------------------------------------------------------------
// Courbes de défaut.
	std::vector<const ICM_DefaultCurve*> TabCrbDefaut (NbIssuers+2);  // Le rajout de 2 est lié aux courbes de traxx

	if (FastStrip)
	{
		for (i = 0;i<NbIssuers;i++)
		{
			ARM_Vector* Spread = NameSpreads->TruncRowAsVector(i,nb_plot) ; 

			TabCrbDefaut[i] = new ICM_DefCurvStr(Spread,
												 &PlotNum,
												 CalibRecov[i],
												 CrbTaux,
												 Label[i]);
			// Delete
			MyDelete(Spread);
		}

		// On cree les 2 courbes de Traxx
		TabCrbDefaut[NbIssuers] = new ICM_DefCurvStr(&SpreadTraxxEur,
													 &PlotNum,
													 RecovTraxx,
													 CrbTaux,
													 "TRAXX_EUR");

		TabCrbDefaut[NbIssuers+1] = new ICM_DefCurvStr(&SpreadTraxxUS,
													   &PlotNum,
													   RecovTraxx,
													   CrbTaux,
													   "TRAXX_US");
	}
	else
	{
		for (i = 0;i<NbIssuers;i++)
		{
			ARM_Vector* Spread = NameSpreads->TruncRowAsVector(i,nb_plot) ; 

			TabCrbDefaut[i] = new ICM_Constant_Piecewise(AsOf,
														 Term,
														(ARM_Vector*)Spread->Clone(),
														 CalibRecov[i],
														 CrbTaux,
												K_ADJUSTED,	// intRule
												K_ADJUSTED,	// adjStartDate
														 qCredit_Adjust20,
														 ARM_DEFAULT_COUNTRY,
														 Label[i],
														 false,
														 NULL,//2 NULL,
														 K_QUARTERLY,qDEFCURVE_DICHO,"STD", 
														 ARM_Currency(ARM_DEFAULT_COUNTRY).GetCreditStartDateLag()
														 ,ICM_Parameters());
			// Delete
			MyDelete(Spread);
		}

		// On cree les 2 courbes de Traxx
		TabCrbDefaut[NbIssuers] = new ICM_Constant_Piecewise(AsOf,
														 Term,
														 &SpreadTraxxEur,
														 RecovTraxx,
														 CrbTaux,
												K_ADJUSTED,	// intRule
												K_ADJUSTED,	// adjStartDate
														 qCredit_Adjust20,
														 ARM_DEFAULT_COUNTRY,
														 "TRAXX_EUR",
														 false,
														 NULL,//2 NULL,
														 K_QUARTERLY,qDEFCURVE_DICHO,"STD", 
														 ARM_Currency(ARM_DEFAULT_COUNTRY).GetCreditStartDateLag()
														 ,ICM_Parameters());

		TabCrbDefaut[NbIssuers+1] = new ICM_Constant_Piecewise(AsOf,
														 Term,
														 &SpreadTraxxUS,
														 RecovTraxx,
														 CrbTaux,
												K_ADJUSTED,	// intRule
												K_ADJUSTED,	// adjStartDate
														 qCredit_Adjust20,
														 ARM_DEFAULT_COUNTRY,
														 "TRAXX_US",
														 false,
														 NULL,//2 NULL,
														 K_QUARTERLY,qDEFCURVE_DICHO,"STD", 
														 ARM_Currency(ARM_DEFAULT_COUNTRY).GetCreditStartDateLag()
														 ,ICM_Parameters());
	}
// ---------------------------------------------------------------------------
// Correlation
	ARM_Vector MaturityTraxx = ARM_Vector(MaturityBaseCorrel);
	int NbNomTraxx = 125;
	bool JustEuro = false;
	bool JustUS = false;
	ARM_Vector proportion = ARM_Vector(2,1);
	proportion.Elt(0) = ProportionTraxxEuro;		// Proportion Traxx EUR
	proportion.Elt(1) = 1-ProportionTraxxEuro;	// Proportion Traxx US
	if (fabs(ProportionTraxxEuro - 1)<1e-6)
		JustEuro = true;
	else if (fabs(ProportionTraxxEuro)<1e-6)
		JustUS = true;

	std::vector<std::string> LabelTraxx (2); 
	// char** LabelTraxx = new char*[2];
	LabelTraxx[0] = "TRAXX_EUR";
	LabelTraxx[1] = "TRAXX_US";

	std::vector<std::string> LabelTraxxEur (NbNomTraxx); 
	// char** LabelTraxxEur = new char*[NbNomTraxx];
	for (i=0;i<NbNomTraxx;i++)
		LabelTraxxEur[i] = "TRAXX_EUR";
		
	std::vector<std::string> LabelTraxxUS (NbNomTraxx); 
	// char** LabelTraxxUS = new char*[NbNomTraxx];
	for (i=0;i<NbNomTraxx;i++)
		LabelTraxxUS[i] = "TRAXX_US";
	
// Vecteur de strike
	ARM_Vector StrikEUR = ARM_Vector(nbColCorrel,0.);
	ARM_Vector StrikUS  = ARM_Vector(nbColCorrel,0.);
	
	for (i=0;i<nbColCorrel;i++)
	{
		StrikEUR.Elt(i) = 100 * Correl[i];
		StrikUS.Elt(i) = 100 * Correl[(sizeMatBaseCorrel + 1)*nbColCorrel+ i];
	}

// On recupere les correls
	ARM_Vector CorrelEUR = ARM_Vector(sizeMatBaseCorrel*nbColCorrel,0.);
	ARM_Vector CorrelUS  = ARM_Vector(sizeMatBaseCorrel*nbColCorrel,0.);
	for (i=0;i<sizeMatBaseCorrel*nbColCorrel;i++)
	{
		CorrelEUR.Elt(i) = 100 * Correl[i + nbColCorrel];
		CorrelUS.Elt(i) = 100 * Correl[i + (2+sizeMatBaseCorrel)*nbColCorrel];
	}

	ARM_Matrix* CorrelationEur = new ARM_Matrix(MaturityTraxx.GetSize(),nbColCorrel);
	ARM_Matrix* CorrelationUS = new ARM_Matrix(MaturityTraxx.GetSize(),nbColCorrel);
	for (j = 0; j<sizeMatBaseCorrel;j++) 
	{
		for (int il=0;il<nbColCorrel;il++)
		{
			CorrelationEur->Elt(j,il) = CorrelEUR.Elt(il+j*nbColCorrel);
			CorrelationUS->Elt(j,il) = CorrelUS.Elt(il+j*nbColCorrel);
		}
	}
	
// Objet Vol Curve : on doit en faire un pour chaque Traxx
	vector<const ARM_VolCurve*> Vector_VolCurve(2);
	Vector_VolCurve[0] = new ARM_VolLInterpol(AsOf,
											  (ARM_Vector*)MaturityTraxx.Clone(),
											  (ARM_Vector*)StrikEUR.Clone(),				// Strike en bp
											  (ARM_Matrix*)CorrelationEur->Clone(),			// Valeur de la correl en %
											  1,
	 										  K_SMILE_VOL);
	
	Vector_VolCurve[1] = new ARM_VolLInterpol(AsOf,
											  (ARM_Vector*)MaturityTraxx.Clone(),
											  (ARM_Vector*)StrikUS.Clone(),				// Strike en bp
											  (ARM_Matrix*)CorrelationUS->Clone(),		// Valeur de la correl en %
											  1,
	 										  K_SMILE_VOL);

// On cree un index pour chaque traxx
	vector<const ICM_Credit_Index*> VIndex(2);
	ARM_Vector YT(1); YT.Elt(0) = 5.0;
	ARM_Vector spread(1); spread.Elt(0) = 0.;
					VIndex[0] = new ICM_Credit_Index(KACTUAL_360, 
								  					  4, 
													  4, 
													  YT,
													  ARM_DEFAULT_COUNTRY,
													  LabelTraxxEur[0],
													  LabelTraxxEur,
													  qAVERAGE,spread,NULL,0,0,K_ARREARS,0,K_ARREARS,0,qCredit_Adjust20,5,2);

	// MyDeleteTab(LabelTraxxEur);


					VIndex[1]  = new ICM_Credit_Index(KACTUAL_360, 
													  4, 
													  4, 
													  YT,
													  ARM_DEFAULT_COUNTRY,
													  LabelTraxxUS[0],
													  LabelTraxxUS,
													  qAVERAGE,spread,NULL,0,0,K_ARREARS,0,K_ARREARS,0,qCredit_Adjust20,5,2);

	//MyDeleteTab(LabelTraxxUS);

// Correlation
ICM_Smile_Correlation* Corr = NULL;
ARM_Vector Vsmilestrikelow(2,DataMez[1]);			// On initialise avec les strike down
ARM_Vector Vsmilestrikehigh(2,DataMez[2]);			// On initialise avec les strike Up

if (JustEuro)
{
	Vsmilestrikelow.Elt(1) = 0.;
	Vsmilestrikehigh.Elt(1) = 0.;
}
else if (JustUS)
{
	Vsmilestrikelow.Elt(0) = 0.;
	Vsmilestrikehigh.Elt(0) = 0.;
}

// FIXMEFRED: mig.vc8 (28/05/2007 14:22:55):cast
if (Rescal)
		Corr = new ICM_Smile_Correlation(AsOf,
										"CORRSTR",
										&(*Vector_VolCurve.begin()),
										LabelTraxx,
										&proportion,
										&(*VIndex.begin()));
	else 
		Corr = new ICM_Smile_Correlation(AsOf,
										 "CORRSTR",
										 &(*Vector_VolCurve.begin()), 
										 LabelTraxx,
										 &proportion,
										 &Vsmilestrikelow,
										 &Vsmilestrikehigh);

	// MyDelete(Vector_VolCurve,2);
	delete Vector_VolCurve[0]; delete Vector_VolCurve [1]; 
// ---------------------------------------------------------------------------
// Model
ICM_ModelMultiCurves* Model = new ICM_ModelMultiCurves(// NbIssuers+2,
													   TabCrbDefaut,
													   (ARM_ZeroCurve*) CrbTaux->Clone(), 
													   LossRecov,
													   Corr);
// ---------------------------------------------------------------------------
// Objet Mezz
ICM_Mez* mez = NULL;
double FloatingPayerAmount = MezzAmount;
int stubrule = 1;
int FrequencyDL = Frequency;

	mez = new ICM_Mez(AsOf,
					  EndDate,
					  (ARM_Date*)0,
					  (ARM_Date*)0,
					  FixedRate,
					  K_ADJUSTED,	// intRule
					  K_ADJUSTED,	// adjStartDate
					  SubAmount,
					  MezzAmount,
					  Label,
					  IssuersNotionals,
					  // NbIssuers,
					  Frequency,
					  KACTUAL_360,
					  MezzAmount,
					  qACCRUED_SETTLED,
					  ARM_DEFAULT_COUNTRY,
					  FloatingPayerAmount,
					  stubrule,
					  0,
					  FrequencyDL,
					CREDIT_DEFAULT_VALUE, // const double& Binary /* = CREDIT_DEFAULT_VALUE*/,
					std::string(), // const std::string&  payCalName /* = NULL*/,
					qRunning_Leg, // const qCredit_Leg_Type& TypeFeeLeg/* =qRunning_Leg*/,
					qStandart_Recovery_Leg, // const qCredit_Leg_Type& TypeDefLeg/* =qStandart_Recovery_Leg*/,
					INCLUDE_MATURITY // const bool& IncludeMaturity /* = INCLUDE_MATURITY*/);
					  );
// ---------------------------------------------------------------------------
// Construction de l'objet ICM_Matrix qui va nous permettre de passer les infos concernant la méthode numérique.
ICM_Parameters* Matrix = new ICM_Parameters();

const char* tmp0 = "COPULA";
	ARM_Vector V0 = ARM_Vector(1,1);
	Matrix->Push((ARM_Vector*)V0.Clone(), (char*)tmp0);
	
	const char* tmp1 = "INTEGRATION_STEP_1";
	ARM_Vector V1 = ARM_Vector(1,NbStep);
	Matrix->Push((ARM_Vector*)V1.Clone(), (char*)tmp1);
	
	const char* tmp2 = "FREEDOM_DEGREE";
	ARM_Vector V2 = ARM_Vector(1,1);
	Matrix->Push((ARM_Vector*)V2.Clone(), (char*)tmp2);

	const char* tmp3 = "INTEGRATION_METHOD";
	// 0 -> qGAUSSLEGENDRE; // 1 -> qGAUSSHERMITE; 	// 2 -> qTRAPEZE
	ARM_Vector V3 = ARM_Vector(1,MethodNum );
	Matrix->Push((ARM_Vector*)V3.Clone(), (char*)tmp3);
// ---------------------------------------------------------------------------
// Pricer
ICM_Pricer_Distrib_Smile Pricer_Mez ; Pricer_Mez.Set(mez,Model,*Matrix,AsOf);
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// Resultat
	vector<double> losses;

	if (Frequency > 0)
	{
		int nbFlow = Nb_Flux(Maturity,Frequency);
		Resultat.resize(nbFlow);
		for (i=1;i<nbFlow+1;i++)
		{
			double time = ((double)i)/Frequency;
			Resultat[i-1] = Pricer_Mez.ExpectedLossTranche(time,losses);
		}
	}
	else
	{
		Resultat.resize(1);
		Resultat[0] = Pricer_Mez.ExpectedLossTranche(Maturity,losses);
	}
	result.setDouble(0.);
// ---------------------------------------------------------------------------		
	// DELETE
	// MyDelete(Pricer_Mez);
	MyDelete(Matrix);
	MyDelete(mez);
	MyDelete(Corr);
	//MyDelete(VIndex,2);
	delete VIndex[0]; delete VIndex[1];
	//MyDelete(Vector_VolCurve,2);
	delete Vector_VolCurve[0]; delete Vector_VolCurve[1];
	MyDelete(CorrelationUS);
	MyDelete(CorrelationEur);
	// if (Label)
	// 	FreePointerTabChar(Label,NbIssuers);

	MyDeleteTab(TabCrbDefaut);
	MyDelete(CrbTaux);
	// MyDeleteTab(pIssuersNotionals);
	MyDeleteTab(CalibRecov);
	//MyDeleteTab(LossRecov);

	return ARM_OK;
	}
 

	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}

	catch (...)
	{
		result.setMsg ("ARM_ERR: Price : unrecognized failure");
		return ARM_KO;
	}
};

// ---------------------------------------------------------------------------
// ------------------------  Interpolation de la correl ----------------------
/** 
extern long ICMLOCAL_STR_LinearInterpol(vector<double> XRef,
										vector<double> YRef,
										vector<double> XTarget,
										vector<double>& YTarget,
										ARM_result& result)
{
	ICMMSG(WARN,"Using ICMLOCAL_STR_LinearInterpol" ); 
	CCString msg ("");
try
{
	InterpolLineaire(XRef,YRef,XTarget,YTarget);
	
	return ARM_OK;
}
 
	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}

	catch (...)
	{
		result.setMsg ("ARM_ERR: Linear Interpolation : unrecognized failure");
		return ARM_KO;
	}

};

**/ 
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
extern long ICMLOCAL_Str_CorrelInterTranche(long NbSimul,
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
											ARM_Vector& ProbaDefautTranche,
											ARM_Vector& VarianceTranche,
											vector<double*>& MatriceEvtDefaut,
											vector<double*>& MatriceTpsDefaut,
											ARM_result& result)
{
	ICMMSG(WARN,"Using ICMLOCAL_Str_CorrelInterTranche" ); 
	CCString msg ("");
	int i = 0;
	int j = 0;
	double Taux = taux*100.;
try
{
	double Basis = 365.;
	// ---------------------------------------------------------------------------
	// Recovery
	if (CalibRecovery.size() <NbName)
	{	
		double tempCalib = CalibRecovery[0];
		CalibRecovery.clear();
		CalibRecovery.resize(NbName);

		for (i=1;i<NbName;i++)
			CalibRecovery[i] = tempCalib;
	}
	// ---------------------------------------------------------------------------
	// Maturity
	double Maturity  = Maturite;
	ARM_Date AsOf; 
	AsOf.Today();
	double AsOfExcelDate = AsOf.GetJulian() - 2415019;
	if (Maturite>100)
		Maturity = (Maturite - AsOfExcelDate)/Basis;
	// ---------------------------------------------------------------------------
	// Matrice de Spread
	int nb_plot = YF_plot.size();
	int NbRow = NbName;
	int NbCol = nb_plot; 
	ICM_QMatrix<double>* NameSpreads = new ICM_QMatrix<double>(NbRow,NbCol,CREDIT_DEFAULT_VALUE);
	for (i=0;i<NbRow;i++)
	{
		// Cas Vraie Matrice
		if (MatIssuersSpread.size()==NbCol*NbRow)
		{
			for (j=0;j<NbCol;j++)
				NameSpreads->SetValue(i,j,MatIssuersSpread[i*NbCol+j]);
		}
		// Cas d'un input = une ligne
		else if (MatIssuersSpread.size() == NbCol)
		{
			for (j=0;j<NbCol;j++)
				NameSpreads->SetValue(i,j,MatIssuersSpread[j]);
		}
		// Cas d'une colonne
		else if (MatIssuersSpread.size() == NbName)
		{
			for (j=0;j<NbCol;j++)
				NameSpreads->SetValue(i,j,MatIssuersSpread[i]);
		}
		// Cas une seule valeur
		else if (MatIssuersSpread.size() == 1)
		{
			for (j=0;j<NbCol;j++)
				NameSpreads->SetValue(i,j,MatIssuersSpread[0]);
		}
		else
		{
			result.setMsg ("ARM_ERR: Issuers Spread : Size Error");
			return ARM_KO;
		}	
	}
	// ---------------------------------------------------------------------------
	// Vecteur de plot
	ARM_Vector PlotNum = ARM_Vector(YF_plot);
	const int sizePlot = YF_plot.size();
	if (YF_plot[0]>100)
	{
		for (i=0;i<sizePlot;i++)
			PlotNum.Elt(i) = (YF_plot[i] - AsOfExcelDate)/Basis;
	}

	// ---------------------------------------------------------------------------
	// Labels
	std::vector<std::string> Label(NbTranche); 
	// char** Label = new char*[NbTranche];
	for (i=0;i<NbTranche;i++)
	{
		std::stringstream sstr ;sstr<<"TRANCHE_"<<i;
		sstr>>Label[i] ; 
		// sprintf(Label[i],"TRANCHE_%i",i);
	}
	
	// ---------------------------------------------------------------------------
	// Matrice des strikes : Les tranches down puis les tranches up
	ICM_QMatrix<double>* TdownTup = new ICM_QMatrix<double>(2,NbTranche);
		// On remplit la matrice
		for (i=0;i<NbTranche;i++)
		{
			(*TdownTup)(0,i) = atof((char*)StrikeDownUp[i]);						// On recupere le strike Down
			(*TdownTup)(1,i) = atof((char*)StrikeDownUp[i+nbColStrikeDownUp]);	// On recupere le strike Up
		}


	// ---------------------------------------------------------------------------
	// Matrice d'appartenance
	ICM_QMatrix<int>* NameInCDO = new ICM_QMatrix<int>(NbName,NbTranche);
		// On remplit la matrice
		for (i=0;i<NbName;i++)
		{
			for (j=0;j<NbTranche;j++)
				(*NameInCDO)(i,j) = atoi((char*)Appartenance[j + nbColIssuersInCDO*i]);		// On recupere le nom j de la tranche i
		}

	// ---------------------------------------------------------------------------
	// Construction des objets qui seront retournés par la fonction
		ICM_CorrMatrix* CorrelEvtDefaut = NULL;	
		ICM_CorrMatrix* CorrelTpsDefaut = NULL;
	// ---------------------------------------------------------------------------
	// Courbe de Taux
		ARM_ZeroFlat* CrbTaux = new ARM_ZeroFlat(AsOf, Taux);
	// ---------------------------------------------------------------------------
	// Courbes de défaut.
		std::vector<const ICM_DefaultCurve*> TabCrbDefaut (NbName);  
		
		char* Label_Issuers = new char[60];
		for (i = 0;i<NbName;i++)
		{
			ARM_Vector* Spread = NameSpreads->TruncRowAsVector(i,nb_plot) ;
			sprintf(Label_Issuers,"ISSUER_%i",i);

			TabCrbDefaut[i] = new ICM_DefCurvStr(Spread,
												 &PlotNum,
												 CalibRecovery[i],
												 (ARM_ZeroCurve*)CrbTaux,
												 Label_Issuers);
			// Delete
			MyDelete(Spread);
		}	
	// ---------------------------------------------------------------------------
	// Vecteur des probas de défaut à maturité
		ARM_Vector PdefIssuers(NbName);
		ARM_Vector VLossRecovery(NbName);
		ARM_Vector VNotionals(NbName);

		
		 
		for (i=0;i<NbName;i++){
			PdefIssuers[i] = 1 - TabCrbDefaut[i]->SurvivalProba(Maturity);
			VLossRecovery[i] = LossRecovery;
			VNotionals[i]=1.;
		}
	// ---------------------------------------------------------------------------
	// Appel de la fonction
		MatriceCorrelInterTranche(NbSimul,
								  NbTranche,
								  NbName,
								  VNotionals,
								  VLossRecovery,
								  BetaIssuers,
								  TdownTup,
								  *NameInCDO,
								  PdefIssuers,
								  Label,
								  CorrelEvtDefaut,
								  CorrelTpsDefaut,
								  ProbaDefautTranche,
								  VarianceTranche,
								  true);
	
	MatriceEvtDefaut.clear();
	MatriceEvtDefaut.resize(NbTranche);

	MatriceTpsDefaut.clear();
	MatriceTpsDefaut.resize(NbTranche);

	for (i=0;i<NbTranche;i++)
	{
		MatriceEvtDefaut[i] = new double[NbTranche];
		memset(MatriceEvtDefaut[i],'\0',sizeof(double)*NbTranche);
		
		MatriceTpsDefaut[i] = new double[NbTranche];
		memset(MatriceTpsDefaut[i],'\0',sizeof(double)*NbTranche);
	}


	for (i=0;i<NbTranche;i++)
	{
		
		for (j=i+1;j<NbTranche;j++)
		{
			MatriceEvtDefaut[i][j] = CorrelEvtDefaut->GetCorrelation(Label[i],Label[j]);
			MatriceEvtDefaut[j][i] = MatriceEvtDefaut[i][j];
			MatriceTpsDefaut[i][j] = CorrelTpsDefaut->GetCorrelation(Label[i],Label[j]);
			MatriceTpsDefaut[j][i] = MatriceTpsDefaut[i][j];
		}

		MatriceEvtDefaut[i][i] = 1.;
		MatriceTpsDefaut[i][i] = 1.;
		
	}	

	// ----------------------------------------------------
	// Delete
	MyDeleteTab(TabCrbDefaut);
	MyDeleteTab(Label_Issuers);

	// if (Label)
	// 	FreePointerTabChar(Label,NbTranche);

	MyDelete(CorrelEvtDefaut);
	MyDelete(CorrelTpsDefaut);
	MyDelete(TdownTup);
	MyDelete(NameInCDO);

	return ARM_OK;
}
 
	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}

	catch (...)
	{
		result.setMsg ("ARM_ERR: Matrice Correl Inter Tranche  : unrecognized failure");
		return ARM_KO;
	}
}



/*
long ICMLOCAL_Str_PriceAndSensiForCDSOptions(CCString PricingType,
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
									  ARM_result& result)

{
	ICMMSG(WARN,"Using ICMLOCAL_Str_PriceAndSensiForCDSOptions" ); 
	try

	{
	
	int i=0;
	int j=0;
	int size = Maturities.size (); 
	double* pdRates = new double[size];
	char* pTradeDate =  new char[11];
    ARM_Vector* vRates = NULL;

	ARM_Date AsOfDate;
	ARM_Date OptMatu;
	ARM_Date CDSMatu;

	Local_XLDATE2ARMDATE(TradeDate,pTradeDate);
	AsOfDate = (ARM_Date) pTradeDate;

	char* pOptMatu = new char[11];
	char* pCDSMatu = new char[11];

	Local_XLDATE2ARMDATE(OptionMaturity,pOptMatu);
	Local_XLDATE2ARMDATE(CDSMaturity,pCDSMatu);

	OptMatu = (ARM_Date)pOptMatu;
	CDSMatu = (ARM_Date)pCDSMatu;


// ***************************************************************************
// ***************************************************************************
//						CONSTRUCTION DES OBJETS

// Courbe de Taux
	ARM_ZeroFlat* CrbTaux = new ARM_ZeroFlat(AsOfDate, Rate);


//Courbe de Défaut 
	ICM_DefaultCurve** TabCrbDefaut = new ICM_DefaultCurve*[1];  

	
	vector<string> psMatu;


	for(j = 0; j < size; j++)
		psMatu.push_back(string( (const char*)Maturities[j] ));

	for (i=0;i<size;i++)
		pdRates[i] = Spreads[i]/10000.; 
	
	vRates = new ARM_Vector(size,pdRates); 

	TabCrbDefaut[0] = new ICM_Constant_Piecewise(AsOfDate,
										psMatu,
										vRates,
										Recovery,
										CrbTaux,
												K_ADJUSTED,	// intRule
												K_ADJUSTED,	// adjStartDate
										qCredit_Adjust20,
										ARM_DEFAULT_COUNTRY,
										"CrbDefaut",
										true,
										NULL,//2 NULL,
										K_QUARTERLY,qDEFCURVE_DICHO,"STD", ARM_Currency(ARM_DEFAULT_COUNTRY).GetCreditStartDateLag());

// Objet Vol
	ARM_VolFlat* VolIndex = new ARM_VolFlat(AsOfDate,Vol*100.);	// Multiplication pour être ARM_Compliant. On passe des objets en bp.


// Security CDS 

	ICM_Cds* NewCDS = NULL;
	ICM_Credit_Index* Index = NULL; 
	
	std::vector<std::string> LabelIndex (1); 
	LabelIndex[0] = "CrbDefaut";


	// CDS 
	if (UnderlyingTp == "CDS") 
	NewCDS = new ICM_Cds(OptMatu,
			CDSMatu ,
			(ARM_Date*)0,
			(ARM_Date*)0,
			OptMatu,
			CDSMatu,
			(Spread/10000.),
			Notional,Notional,// 0,0,	// notionals 
			Frequency, 
			KACTUAL_360 , 			
			qACCRUED_SETTLED,
			ARM_DEFAULT_COUNTRY, 
			K_SHORTSTART,
			30.,
			DEFAULT_FRQ_DEFLEG, // const int& FrequencyDefLeg = DEFAULT_FRQ_DEFLEG,
			K_ADJUSTED , // const int& intRule = K_ADJUSTED ,
			INCLUDE_MATURITY, // const bool& includematurity = INCLUDE_MATURITY,
			K_ADJUSTED, // const int& adjStartDate = K_ADJUSTED,
			std::string(), // const std::string& char* payCalName = NULL*,
			qRunning_Leg,
			qStandart_Recovery_Leg,
			ISSUER_UNDEFINE, // const string& name = ISSUER_UNDEFINE ,
			CREDIT_DEFAULT_VALUE // const double& Binary = CREDIT_DEFAULT_VALUE);
			); 


	//Index
	if (UnderlyingTp == "Index")
	{
		ARM_Vector YT(1); YT.Elt(0) = 5.0;
		ARM_Vector spread(1); spread.Elt(0) = Spread5Y;
		Index = new ICM_Credit_Index(KACTUAL_360 , 
						K_QUARTERLY,
						K_QUARTERLY, 
						YT,
						ARM_DEFAULT_COUNTRY,
						LabelIndex[0],
						LabelIndex,
						qHOMOTHETIE,
						spread,
						NULL,
						K_MOD_FOLLOWING,
						K_ADJUSTED,
						K_ADVANCE,
						0,
						K_ARREARS,
						0,
						qCredit_Adjust20,5,2);

		NewCDS = new ICM_Cds(OptMatu,
			CDSMatu ,
			(ARM_Date*)0,
			(ARM_Date*)0,
			OptMatu,
			CDSMatu,
			(Spread/10000.),
			Notional,Notional,0,0,
			Index, 
			//"NULL", 
			Frequency, 
			KACTUAL_360 , 
			// Notional,
			qACCRUED_SETTLED,
			ARM_DEFAULT_COUNTRY, 
			// Notional,
			K_SHORTSTART,
			30.,
			DEFAULT_FRQ_DEFLEG, // const int& FrequencyDefLeg = DEFAULT_FRQ_DEFLEG,
			K_ADJUSTED , // const int& intRule = K_ADJUSTED ,
			INCLUDE_MATURITY, // const bool& includematurity = INCLUDE_MATURITY,
			K_ADJUSTED, // const int& adjStartDate = K_ADJUSTED,
			std::string(), // const std::string& char* payCalName = NULL,
			qCDS_INDEX, // const qSecurity_TYPE& cdstype  = qCDS_INDEX,
			ISSUER_UNDEFINE, // const string& name = ISSUER_UNDEFINE ,
			CREDIT_DEFAULT_VALUE // const double& Binary = CREDIT_DEFAULT_VALUE );
			); 

	}

// Objet Correl Matrix 

	ICM_CorrMatrix* Matrix = NULL;		
	ICM_QMatrix<double> values(1,1); values.Fill(1); 
	Matrix = new ICM_CorrMatrix(AsOfDate, "CrbDefaut", LabelIndex,values);


//MODEL 
	ICM_DefaultCurveModel* Model = NULL;
	
	if (UnderlyingTp == "CDS")
		Model = new ICM_DefaultCurveModel(TabCrbDefaut[0],CrbTaux,VolIndex);
	
	ARM_Vector Recoveries (1); 
	
	Recoveries[0] = Recovery; 


	if (UnderlyingTp == "Index")
	{
	Model = new ICM_ModelMultiCurves(1,
								   TabCrbDefaut,
								   CrbTaux,
								   Recoveries,
								   Matrix,
								   VolIndex);
	}



//PRICER CDS-Index 
	ICM_Pricer_Cds* PricerCDS = NULL;
	ICM_Pricer_Advisor advisorCDS = ICM_Pricer_Advisor();
	

	if (UnderlyingTp == "CDS")
		PricerCDS = dynamic_cast<ICM_Pricer_Cds*> (advisorCDS.GeneratePricer(NewCDS, Model, (ARM_CLASS_NAME)ICM_PRICER_CDS,
			CREDIT_DEFAULT_VALUE,(ICM_Parameters*)0,AsOfDate) );


	if (UnderlyingTp == "Index")
		PricerCDS = dynamic_cast<ICM_Pricer_CDSIndex*>  (advisorCDS.GeneratePricer(NewCDS, Model, (ARM_CLASS_NAME)ICM_PRICER_CDS_INDEX,
			CREDIT_DEFAULT_VALUE,(ICM_Parameters*)0,AsOfDate) );


//Security Option 
	ICM_Option* NewSpreadOption = NULL; 

	NewSpreadOption = new ICM_Option(((ICM_Option*)PricerCDS->GetSecurity())->GetUnderMaturityDate(),
										OptMatu,
										 Strike,
										 OptionType,
										 (qDEF_MAT)KOType) ;


//PRICER Option 
	ICM_Pricer* PricerOption = NULL;
	ICM_Pricer_Advisor advisorOpt = ICM_Pricer_Advisor();
	

	if (UnderlyingTp == "CDS")
		PricerOption = advisorOpt.GeneratePricer(NewSpreadOption, Model, (ARM_CLASS_NAME)ICM_PRICER_CDSOPTION,
			CREDIT_DEFAULT_VALUE,(ICM_Parameters*)0,AsOfDate);


	if (UnderlyingTp == "Index")
		PricerOption = advisorOpt.GeneratePricer(NewSpreadOption, Model, (ARM_CLASS_NAME)ICM_PRICER_INDEXOPTION,
			CREDIT_DEFAULT_VALUE,(ICM_Parameters*)0,AsOfDate);

	ICM_Pricer_Option* Pricer_Option = NULL;
	Pricer_Option = (ICM_Pricer_Option *) PricerOption;

//--------------------------------------------------------------------------------
// RESULTS 
	double dur = 0.;
	if (PricingType == "ALL") 
	{
	Res.resize(10);
	Res[0] = PricerCDS->ComputeSpread(0.);
	
	Res[1] = PricerCDS->Compute_Fwd_Spread(OptMatu, CDSMatu,dur);
	Res[2] = PricerCDS->Price(qCMPPRICE)*NbIssuers;
	Res[3] = PricerCDS->Price(qCMPDURATION);

	Res[4] = PricerOption->Price(qCMPPRICE);
	Res[5] = PricerOption->ComputeImpliedVol(MktPrice);
	Res[6] = Pricer_Option->ComputeBSGreeks(ICM_GREEK_DELTA_TYPE);
	Res[7] = Pricer_Option->ComputeBSGreeks(ICM_GREEK_GAMMA_TYPE);
	Res[8] = Pricer_Option->ComputeBSGreeks(ICM_GREEK_VEGA_TYPE);
	Res[9] = Pricer_Option->ComputeBSGreeks(ICM_GREEK_THETA_TYPE);
	//Res[10] = Pricer_Option->ComputeBSGreeks(ICM_GREEK_RHO_TYPE);

	}

	if (PricingType == "NOSENSI") 
	{
	Res.resize(6);

	Res[0] = PricerCDS->ComputeSpread(0.);
	Res[1] = PricerCDS->Compute_Fwd_Spread(OptMatu, CDSMatu, dur);
	Res[2] = PricerCDS->Price(qCMPPRICE)*NbIssuers;
	Res[3] = PricerCDS->Price(qCMPDURATION);

	Res[4] = PricerOption->Price(qCMPPRICE);
	Res[5] = PricerOption->ComputeImpliedVol(MktPrice);
	}

// -------------------------------------------------------------------------------		
// DELETE

	MyDeleteTab(pTradeDate);
	MyDeleteTab(pCDSMatu);
	MyDeleteTab(pOptMatu);
	MyDeleteTab(pdRates);
	//MyDeleteTab(Recoveries);

	MyDelete(VolIndex);
	MyDeleteTab(TabCrbDefaut);
	MyDelete(CrbTaux);
	MyDelete(Matrix);
	MyDelete(vRates);

	if (UnderlyingTp == "Index")
	{
	MyDelete(Index);
	}

	MyDelete(NewCDS);
	MyDelete(PricerCDS);
	MyDelete(Model);
	MyDelete(NewSpreadOption);
	MyDelete(Pricer_Option);

	
	return ARM_OK; 

	}

	catch (...)
	{
		result.setMsg ("ARM_ERR: ICMLOCAL_Str_PriceAndSensiForCDSOptions : unrecognized failure");
		return ARM_KO;
	}

}
*/