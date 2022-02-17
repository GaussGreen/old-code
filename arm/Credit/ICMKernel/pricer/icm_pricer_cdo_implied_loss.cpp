#include "ARMKernel\glob\firsttoinc.h"
#include "ICMKernel\pricer\icm_pricer_cdo_implied_loss.h"
#include "ICMKernel\pricer\icm_pricer_adviser.h"
#include "ICMKernel\pricer\icm_pricer_homogeneous_smile.h"
#include "ICMKernel\inst\icm_collateral.h"
#include "ICMKernel\mod\modelmulticurves.h"
#include "ICMKernel\inst\icm_mez.h"
#include "ICMKernel/glob/icm_mktdatamng.h"
#include "ICMKernel\crv\icm_implied_loss_tree.h"


ICM_Pricer_ImpliedLoss::~ICM_Pricer_ImpliedLoss()
	{
		if (itsPricer) delete itsPricer;
		itsPricer = NULL;
	}

void ICM_Pricer_ImpliedLoss::Set(ARM_Security *option, 
			 ARM_Object *mod,
			 const ICM_Parameters& parameters ,const ARM_Date&asof)
	{
		ICM_Pricer::Set(option,mod,parameters,&asof);
		itsTree = (ICM_ImpLossTree*)((ICM_MktDataMng*)mod)->GetLossTree("USERDEF",asof);
	}

//---------------------------------------------------------------------------
// Adjust transitions probabilities and calibrate Tree
//---------------------------------------------------------------------------
void ICM_Pricer_ImpliedLoss::CptFwdLoss(int idxstart,int idxend,int dw_loss ,int up_loss,double notional_dw,double notional_up)
{

	vector<double> EL;
	EL.resize(idxend+1);
	itsDistribLoss.clear();

	for (int k=0;k<idxend-idxstart+1;k++) {EL[k]=0.;}

/*
	itsDistribLoss[0.]=0.;
	for (int noindex=0;noindex<itsTree->GetTreeVector().size();noindex++)
	for (int i=idxstart;i<=idxend;i++)	
	{itsDistribLoss[itsYearTerms[i]]=itsTree->CptCumEL(i,loss,noindex,false)/(loss*itsTree->GetLossUnit());}
*/

	for (int noindex=0;noindex<itsTree->GetTreeVector().size();noindex++)
	{
		//pour chacun des indices
		//for (int noloss=dw_loss;noloss<=up_loss;noloss++)
		for (int noloss=0;noloss<TRX_EUR_NBNAMES;noloss++)
		{
			double value=0.;

			if (noloss<=idxstart)
				value = itsTree->GetTreeVector(noindex).data(idxstart,noloss).state;
			else
				continue;

			//pour la date de départ de la fwd start cdo
			//on met la valeur des variables d'état de l'arbre partout à 0 sauf au Kème élément qui vaut 1
			itsTree->ResetValuesForFwdStates(noindex,idxstart,noloss,value);
			itsTree->DiffuseStateLosses(idxstart,idxend,true);

			for (int slice=idxstart;slice<=idxend;slice++)	
			{
				//if (noloss<=slice)
				{EL[slice] += itsTree->CptCumEL(slice,dw_loss,up_loss,noindex,true,notional_dw,notional_up);}
			}
		}

	}

/*
	for (int noindex=0;noindex<itsTree->GetTreeVector().size();noindex++)
	{
		double value=0.;

		//pour chacun des indices
		for (int slice=idxstart;slice<=idxend;slice++)
		for (int noloss=0;noloss<=slice;noloss++)
		{
		//pour la date de départ de la fwd start cdo
		//on met la valeur des variables d'état de l'arbre partout à 0 sauf au Kème élément qui vaut 1
		itsTree->ResetValuesForFwdStates(noindex,idxstart,noloss,1.);
		itsTree->DiffuseStateLosses(idxstart,idxend,true);

		EL[slice] += itsTree->GetTreeVector(noindex).data(idxstart,noloss).state*
									itsTree->CptCumEL(slice,noloss,noindex,true);

		}

	}
*/

	itsDistribLoss[0.]=0.;

	for (int i=idxstart;i<=idxend;i++)	
	{itsDistribLoss[itsYearTerms[i]]=EL[i];}

}


//---------------------------------------------------------------------------
// Price computation
//---------------------------------------------------------------------------
double ICM_Pricer_ImpliedLoss::ComputePrice(qCMPMETH measure) 
{
	double SPLINE = 0.;

	// if (GetParameters() != NULL)
	// {	
	if (GetParameters().GetColVect("SPLINE"))
		{SPLINE = (long) GetParameters().GetColVect("SPLINE")->Elt(0);}
	// }

	ARM_Date AsOf = GetAsOfDate();
	ICM_Mez* cdo = (ICM_Mez*) GetSecurity();

	ICM_Security* feeleg = cdo->GetFeeLeg()->GetCreditInfos();

	feeleg->ComputeYF(AsOf);
	// 14514 int nbflows = feeleg->GetNumFlows();
	int nbflows = feeleg->GetAccStartDates().size();

	itsYearTerms.clear();

	double yfstart = feeleg->GetYFAccStartDates().Elt(0);
	double yfend = feeleg->GetYFAccEndDates().Elt(nbflows-1);

	int idx_start=-1;
	int idx_matu=-1;

	int noindex=0;
	double value = 0.;

	ARM_Date Maturity = feeleg->GetAccEndDates().Elt(nbflows-1);

	if (!itsTree->IsCalibrated()) 
	{
		if (SPLINE==0.)
			itsTree->calibrate(Maturity.GetJulian(),true);
		else
			itsTree->calibrate(Maturity.GetJulian(),false);
	}

	string defcurvename = cdo->GetCollateral()->GetIssuersLabels(0);
	ICM_DefaultCurve* defcurve = (ICM_DefaultCurve*)((ICM_MktDataMng*)GetModel())->GetDefaultCurve(defcurvename,AsOf);

	//création du modele
	std::vector<const ICM_DefaultCurve*> DefaultCurves(1); 
	// ICM_DefaultCurve* DefaultCurves[1]; 
	DefaultCurves[0] = defcurve;

	ICM_Correlation* correl = itsTree->GetCorrelation();
	string correlname = correl->GetStructName();

	ARM_ZeroCurve* ircurve = (ARM_ZeroCurve*)((ICM_MktDataMng*)GetModel())->find(IRCURVE);

	ICM_ModelMultiCurves mmc(DefaultCurves,ircurve,NULL,correl);

	int near_inf = 0,near_sup = 0;

	if (SPLINE==0.){

	near_sup = itsTree->GetTreeVector(noindex).depth()-1;

	for (int i=0;i<itsTree->GetTreeVector(noindex).depth();i++)
	{
		itsYearTerms.push_back(itsTree->GetTreeVector(noindex).Time(i));

		value = itsTree->GetTreeVector(noindex).Time(i);

		if ((i+1)<itsTree->GetTreeVector(noindex).depth()){
		if ((itsTree->GetTreeVector(noindex).Time(i)<=yfstart) &&
			(itsTree->GetTreeVector(noindex).Time(i+1)>yfstart))
		{near_inf = i;}}
		
		if (CHECK_EQUAL(itsTree->GetTreeVector(noindex).Time(i),yfstart))
		{idx_start = i;}

		if ((i+1)<itsTree->GetTreeVector(noindex).depth()){
		if ((itsTree->GetTreeVector(noindex).Time(i)<=yfend) &&
			(itsTree->GetTreeVector(noindex).Time(i+1)>yfend))
		{near_sup = i;}
		else
		if ((itsTree->GetTreeVector(noindex).Time(i)<yfend) &&
			(itsTree->GetTreeVector(noindex).Time(i+1)>=yfend))
		{near_sup = i+1;}}

		if ((CHECK_EQUAL(itsTree->GetTreeVector(noindex).Time(i),yfend)) ||
			((yfend<itsTree->GetTreeVector(noindex).Time(i))&&(idx_matu<0)))
		{idx_matu = i;break;}
	}

	if (idx_start<0) 
		{idx_start=near_inf;}

	if (idx_matu<0) 
		{idx_matu=near_sup;}

	}
	else
	{
		for (int i=0;i<itsTree->GetTimeStep().Getnbcols();i++)
			{itsYearTerms.push_back(itsTree->GetTimeStep()(0,i));}

/*		ARM_Vector* FinalSched = NULL;
		MergeDates(&FinalSched,&feeleg->GetYFAccStartDates(),&feeleg->GetYFAccEndDates());

		for (int i=0;i<FinalSched->GetSize();i++)
			{itsYearTerms.push_back(FinalSched->Elt(i));}
		
		if (FinalSched)
			delete (FinalSched);
*/

	}

	// double notional = cdo->GetCollateral()->GetIssuersNotional(0);
	
	double SumNot = cdo->GetCollateral()->SumNotionals(cdo->GetStartDateNA()); 

	double notional_dw = 0.;
	double notional_up = 0.;

	notional_dw = cdo->GetPercentLow(cdo->GetFeeLeg()->GetStartDate())*SumNot;
	notional_up = (cdo->GetPercentLow(cdo->GetFeeLeg()->GetStartDate()))+cdo->GetPercentHight(cdo->GetFeeLeg()->GetStartDate())*SumNot;

	double lossunit = itsTree->GetLossUnit();

	double dw_loss = floor(notional_dw/lossunit);
	double up_loss = floor(notional_up/lossunit);

	if (SPLINE==0.)
	{CptFwdLoss(idx_start,idx_matu,(int)dw_loss,(int)up_loss,notional_dw,notional_up);}
	else
	{CptFwdLoss(cdo->GetPercentLow(cdo->GetFeeLeg()->GetStartDate()),cdo->GetPercentLow(cdo->GetFeeLeg()->GetStartDate())+cdo->GetPercentHight(cdo->GetFeeLeg()->GetStartDate()));}

	ICM_Pricer_Advisor Advisor;
	if (itsPricer) delete itsPricer;
	itsPricer = NULL;
	itsPricer = (ICM_Pricer_Distrib_Smile*) Advisor.GeneratePricer(GetSecurity(),
																&mmc,
																ICM_PRICER_HOMOGENEOUS_SMILE,
																CREDIT_DEFAULT_VALUE,
																&GetParameters(),
																GetAsOfDate());

	itsPricer->Price(qCMPPRICE);
	itsPricer->ResetRootPricer();

	itsPricer->setDistribLoss(itsDistribLoss);
	double result = itsPricer->ComputePrice(measure);

	return (result);
}	


void 
ICM_Pricer_ImpliedLoss::View(char* id, FILE* ficOut)
{
	FILE* fOut;
	char  fOutName[200];

	if ( ficOut == NULL )
	{
		ARM_GetViewFile(ARM_VIEW_FILE, id, fOutName);
		(void) unlink(fOutName);
		fOut = fopen(fOutName, "w"); 
	}
	else	fOut = ficOut;

	int size =0;

	fprintf(fOut, "\t\t\t ----------------- Implied Loss Tree Pricer ----------------- \n\n");
	itsPricer->View(id,fOut); 

	if ( ficOut == NULL )fclose(fOut);
}

double ICM_Pricer_ImpliedLoss::Accrued(void) 
	{
		ICMTHROW(ERR_INVALID_ARGUMENT,"Not Implemented") ; 
	}
double ICM_Pricer_ImpliedLoss::FeeLegPV(void) 
	{
		ICMTHROW(ERR_INVALID_ARGUMENT,"Not Implemented") ; 
	}
double ICM_Pricer_ImpliedLoss::DefLegPV(void) 
	{
		ICMTHROW(ERR_INVALID_ARGUMENT,"Not Implemented") ; 
	}
double ICM_Pricer_ImpliedLoss::Compute_Fwd_Spread(const  ARM_Date &,const  ARM_Date &, double& dur)
	{
		ICMTHROW(ERR_INVALID_ARGUMENT,"Not Implemented") 
	}
double ICM_Pricer_ImpliedLoss::ComputeDuration(void) 
	{
		ICMTHROW(ERR_INVALID_ARGUMENT,"Not Implemented") ; 
	}
double ICM_Pricer_ImpliedLoss::ComputeSpread(const double &) 
	{
		ICMTHROW(ERR_INVALID_ARGUMENT,"Not Implemented") ; 
	}
double ICM_Pricer_ImpliedLoss::ComputeImpliedVol(const double &) 
	{
		ICMTHROW(ERR_INVALID_ARGUMENT,"Not Implemented") ; 
	}


//---------------------------------------------------------------------------
// Adjust transitions probabilities and calibrate Tree
//---------------------------------------------------------------------------
void ICM_Pricer_ImpliedLoss::CptFwdLoss(double dw_loss ,double up_loss)
{

	vector<double> EL;
	itsDistribLoss.clear();

	for (int noindex=0;noindex<itsYearTerms.size();noindex++)
	{EL.push_back(itsTree->CptCumELSpline(dw_loss,up_loss,itsYearTerms[noindex]));}

	itsDistribLoss[0.]=0.;

	for (int i=0;i<=itsYearTerms.size();i++)	
	{itsDistribLoss[itsYearTerms[i]]=EL[i];}

}




// 	virtual 
double ICM_Pricer_ImpliedLoss::ComputeSensitivity(const qSENSITIVITY_TYPE& typesensi, 
		const std::string& plot,
		const std::string& label,
		double epsilon , double epsilonGamma) 
{
	return ICM_Pricer::ComputeSensitivity(typesensi,plot,label,epsilon,epsilonGamma ); 
}
