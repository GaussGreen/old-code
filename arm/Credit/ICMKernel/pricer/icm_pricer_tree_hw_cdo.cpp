#include "ICMKernel\pricer\icm_pricer_tree_hw_cdo.h"
#include "ICMKernel\mod\modelmulticurves.h"
#include "ICMKernel\inst\icm_collateral.h"
#include "ICMKernel\mod\icm_lossunits.h"
#include "ICMKernel\util\icm_rootfinder1D.h"
#include "ICMKernel\random\icm_randomRan2.h"
// clude "ICMKernel\crv\icm_interpoldefcrv.h"
#include "ICMKernel\crv\icm_distriblosscalculator.h"

#include <set>

static FILE *stream_hw_tree = NULL;

ICM_QMatrix<double> ICM_Pricer_tree_hw_cdo::itsCombinations=ICM_QMatrix<double>();

void ICM_Pricer_tree_hw_cdo::ComputeFwdMatrix()
{
	int i=0,j=0;
	double nodevalue1 = 0,nodevalue2 = 0.,proba=0.;

	if (itsFwdMatrix)
		{return;}

	ARM_Vector& sched = GenCDOSchedule(false);
	int size = sched.size()-1;
	double deltat = 0.;

	itsFwdMatrix = new ICM_QMatrix<double>(size,size,0.);

	//back propagation
	for (i=size-1;i>=0;i--)
	{
//		if (i!=size)
//		{deltat =sched.Elt(i+1)-sched.Elt(i);}
		
		for (j=i;j>=0;j--)
		{

			if (i==size)
			{
				nodevalue1 = itsTree->GetTree().data(i,j).jumpsize + itsTree->GetTree().data(i,j).theta;
				itsFwdMatrix->SetValue(i,j,nodevalue1);
			}
			else
			{
				nodevalue1 = itsTree->GetTree().data(i+1,j+1).jumpsize + 
										itsTree->GetTree().data(i+1,j+1).theta;
				nodevalue2 = itsTree->GetTree().data(i+1,j).jumpsize + 
										itsTree->GetTree().data(i+1,j).theta;

				proba = deltat * itsTree->GetTree().data(i,j).lambda;
				itsFwdMatrix->SetValue(i,j,proba * nodevalue1 + (1.-proba) * nodevalue2);
			}
		}
	}


}



// ---------------------------------------------------------------------------------
// Pdef2Spreads
// ---------------------------------------------------------------------------------
double ICM_Pricer_tree_hw_cdo::Pdef2Spread(double& date1,double& date2,vector<double>&  pdef, bool spreadok)
{
	if (date1==date2)
	{return 0.;}

	int i=0,beg=0;
	ARM_Vector& sched = GenCDOSchedule(false);
	ICM_ModelMultiCurves* model = (ICM_ModelMultiCurves*) GetModel(); 
	ICM_Ftd& ftd=dynamic_cast<ICM_Ftd&>(*GetSecurity()); 
	ICM_Collateral*collat =ftd.GetCollateral(); 
	double recovery = model->GetDefaultCurve(collat->GetIssuersLabels(0))->GetRecovery();

//	ARM_Vector dates(sched.size(),sched.begin());
//	ARM_Vector defaultproba(pdef.size(),pdef.begin());


	ARM_Vector dates_(sched.size(),&(*sched.begin()));
	dates_ -= date1;
	ARM_Vector defaultproba_(pdef.size(),&(*pdef.begin()));

	for (i=0;i<sched.size();i++)
	{if CHECK_NULL(dates_[i])
		{beg=i;}}

	ARM_Vector dates(sched.size()-beg,0.);
	ARM_Vector defaultproba(sched.size()-beg,0.);
	for (i=beg;i<sched.size();i++)
	{
		dates[i-beg]=dates_[i];
		defaultproba[i-beg]=defaultproba_[i];
	}

	defaultproba -= 1;
	defaultproba /= defaultproba[0];
	defaultproba -= 1.;
	defaultproba *=-1.;


	double spr = 0.;

	if (spreadok)
		spr = QuickFwdSpread(dates,defaultproba,0.,date2-date1,model->GetZeroCurve(),recovery);
	else
		spr = QuickCdsDuration(dates,defaultproba,0.,date2-date1,model->GetZeroCurve(),recovery);

	return (spr);
}

void ICM_Pricer_tree_hw_cdo::Pdef2Spreads(vector<double>&  pdef,vector<double>&  spreads)
{

	int i=0;
	spreads.clear();
	ARM_Vector& sched = GenCDOSchedule(false);
	ICM_ModelMultiCurves* model = (ICM_ModelMultiCurves*) GetModel(); 
	ICM_Ftd& ftd=dynamic_cast<ICM_Ftd&>(*GetSecurity()); 
	ICM_Collateral*collat =ftd.GetCollateral(); 
	double recovery = model->GetDefaultCurve(collat->GetIssuersLabels(0))->GetRecovery();

	ARM_Vector dates(sched.size(),sched.begin());
	ARM_Vector dates_(sched.size(),sched.begin());

// FIXMEFRED: mig.vc8 (28/05/2007 15:29:41):cast
	ARM_Vector defaultproba(pdef.size(),&(*pdef.begin()));

	for (i=0;i<sched.size()-1;i++)
	{
		double spr = QuickFwdSpread(dates_,defaultproba,dates_[i],dates_[sched.size()-1],model->GetZeroCurve(),recovery);
		spreads.push_back(spr);
	}

}

// ---------------------------------------------------------------------------------
// Generation de l'échéancier pour l'arbre following cdo payments dates
// ---------------------------------------------------------------------------------
void ICM_Pricer_tree_hw_cdo::View(char* id, FILE* ficOut)
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

	fprintf(fOut, "\t\t\t ----------------- Hull White Credit Tree Pricer ----------------- \n\n");
	ICM_Pricer_Distrib::View(id,fOut); 
	itsTree->View(id,fOut); 

	fprintf(fOut, "\t\t\t ----------------- Statistics ----------------- \n\n");
	fprintf(fOut,"Trigger(%) : %10.4lf\n",itsPercentTrigger) ;
	fprintf(fOut,"Loss(%) : %10.4lf\n",itsPercentLoss) ;
	fprintf(fOut,"Loss & Trigger(%) : %10.4lf\n",itsPercentLosstrigger) ;

	if ( ficOut == NULL )fclose(fOut);
}

// ---------------------------------------------------------------------------------
// Generation de l'échéancier pour l'arbre following cdo payments dates
// ---------------------------------------------------------------------------------
ARM_Vector& ICM_Pricer_tree_hw_cdo::GenCDOSchedule(bool cdo)
{
	double yf = 0.;
	double yfmatu = 0.;

	if (itsMergedSchedule.size()==0)
	{
		ARM_Date AsOf = GetAsOfDate();
		int step = (int) (itsTimeStep * 365.);
		ICM_Cds* basiscds = (ICM_Cds*) GetSecurity();

		//basiscds->ResetSchedules();
		basiscds->GenerateSchedule(AsOf,step);
		itsSchedule = *(basiscds->GetSchedule());

		itsResetSchedule.clear();
		itsMergedSchedule.clear();

		ARM_Vector* vectmerged;
		yfmatu = itsSchedule.Elt(itsSchedule.size()-1);
		
		yf = 0.;
		for (int k=0;yf<=yfmatu;k++)
		{yf += itsTimeStep;
		if (yf<=yfmatu){
			itsResetSchedule.push_back(yf);}}

		MergeDates(&vectmerged,itsSchedule,itsResetSchedule);
		itsMergedSchedule = *vectmerged;

		if (vectmerged) delete vectmerged;
	}

	if (cdo)
		return itsSchedule;
	else
		return itsMergedSchedule;
}

// ---------------------------------------------------------------------------------
// Generation de l'arbre
// ---------------------------------------------------------------------------------
void ICM_Pricer_tree_hw_cdo::BuildTree()
{

	int i,j;

	if (!itsTree)
	{
		ARM_Vector sch = GenCDOSchedule(false);
		itsTree= new ICM_HWTree_intensity(sch);
	}

	double FlatIntensity=0.;

	if (itsFlatIntensity.size()>0)
	{FlatIntensity = itsFlatIntensity.Elt(0);}

	ARM_ReferenceValue RefIntensity((ARM_Vector*)itsParamsSchedule.Clone(),(ARM_Vector*)itsIntensity.Clone(),K_STEPUP_LEFT,0);
	RefIntensity.SetCalcMethod(K_STEPUP_LEFT) ;
	
	ARM_ReferenceValue RefJumpSize((ARM_Vector*)itsParamsSchedule.Clone(),(ARM_Vector*)itsJumpSize.Clone(),K_STEPUP_LEFT,0);
	RefJumpSize.SetCalcMethod(K_STEPUP_LEFT) ;

	itsTree->GetTree().data(0,0).jumpsize = 0.;
	itsTree->GetTree().data(0,0).lambda=RefIntensity.CptReferenceValue(0.);

	double T = itsMergedSchedule.Elt(itsMergedSchedule.size()-1);
	double T0 = itsMergedSchedule.Elt(0);

	int functype = 0;
	if (itsFuncType.size()>0)
		{ functype =(int) itsFuncType.Elt(0);}

	for (i=0;i<itsMergedSchedule.size();i++)
	{
		double yf = itsMergedSchedule.Elt(i);

		for (j=0;j<=i;j++)
		{
			if ((i>=1)&&(j>=1))
			{
				double jsize = 0.;

				switch (functype)
				{
				default :
				case 0 : 
				jsize = itsTree->GetTree().data(i,j-1).jumpsize +
									itsAlpha[0]*(pow(j,itsBeta[0])-pow(j-1,itsBeta[0]));
				break;
				case 1 : 
				jsize = itsTree->GetTree().data(i,j-1).jumpsize + itsAlpha[0];

				break;
				case 2 : 
				jsize = itsTree->GetTree().data(i,j-1).jumpsize +
								itsAlpha[0]*(exp(itsBeta[0]*j)-exp(itsBeta[0]*(j-1)));

				break;
				case 3 : 
				jsize = itsTree->GetTree().data(i,j-1).jumpsize +
								itsAlpha[0]*
								(exp( itsBeta[0]*pow(j,itsBeta[1]) ) - exp( itsBeta[0]*pow(j-1,itsBeta[1]) )) /exp(itsBeta[0]-1.);
				break;
				case 4 : 
				jsize = itsTree->GetTree().data(i,j-1).jumpsize +
									itsAlpha[0]*(pow(j,itsBeta[0])-pow(j-1,itsBeta[0]))*pow(T-T0,itsBeta[1]);
				break;
				}

				if (jsize>=1.)
					{jsize = 0.99;}
				itsTree->GetTree().data(i,j).jumpsize = jsize;
			}

			double Intensity = RefIntensity.CptReferenceValue(yf);

			if (FlatIntensity)
				{Intensity = FlatIntensity;}

			if (j>=1)
				itsTree->GetTree().data(i,j).lambda=Intensity * pow(j,itsBeta[1]);
			else
				itsTree->GetTree().data(i,j).lambda=Intensity;
		}
	}
}

// ---------------------------------------------------------------------------------
// Calibration des thetas
// ---------------------------------------------------------------------------------
void ICM_Pricer_tree_hw_cdo::CalibrateTheta(int& period)
{

	if (period==0)
	{
		its_call_lastTheta=0.;
		return;
	}

	ICM_ModelMultiCurves* model = (ICM_ModelMultiCurves*) GetModel(); 
	ICM_Ftd& ftd=dynamic_cast<ICM_Ftd&>(*GetSecurity()); 
	ICM_Collateral*collat =ftd.GetCollateral(); 
	ARM_Vector& sched = GenCDOSchedule(false);

	double yf = sched.Elt(period);
	its_call_proba= model->GetDefaultCurve(collat->GetIssuersLabels(0))->DefaultProba(yf);
	its_call_noperiod = period;

	/*
	double _inf = MAX(its_call_lastTheta,-itsTree->GetTree().data(period,period).jumpsize);
	double _sup = 1.-itsTree->GetTree().data(period,period).jumpsize;
	*/

	//double out = RootFinder1D(ff1::mem_call(&ICM_Pricer_tree_hw_cdo::Evaluate,(*this))).ZeroBracket(_inf,_sup);	
	//double result = RootFinder1D(ff1::mem_call(&ICM_Pricer_tree_hw_cdo::Evaluate,(*this))).Dichotomy(_inf,_sup,100,1.E-4,1.E-5,1.E-4);
	//double result = RootFinder1D(ff1::mem_call(&ICM_Pricer_tree_hw_cdo::Evaluate,(*this))).NewtonRaphsonWithBisection(_inf,_sup,100,1.E-4,1.E-5,1.E-4);

	double result = MAX(its_call_proba - itsTree->SumJumps(period),its_call_lastTheta+1.e-5);

/*
	if ((result + itsTree->GetTree().data(period,period).jumpsize)>1.)
		result = 1. - itsTree->GetTree().data(period,period).jumpsize;
*/
	if (result<0.) result = 0.;

	itsTree->SetTheta(period,result);
	its_call_lastTheta = result; 

	#ifdef _DEBUG
	double value = result + itsTree->GetTree().data(period,period).jumpsize;
	#endif
}


double ICM_Pricer_tree_hw_cdo::Evaluate(const double& x)
{
	double value = 0.;
	itsTree->SetTheta(its_call_noperiod,x);

	value = its_call_proba - itsTree->CorpEsp(its_call_noperiod,true);
	value = fabs(value);
			
	return (value);
}

// ---------------------------------------------------------------------------------
// Calcul de l'EL cas homogene pour la période i
// ---------------------------------------------------------------------------------
void ICM_Pricer_tree_hw_cdo::HomogEL(int& period,int nbnames,const double& LossUnit,const double& tranche_down,const double& tranche_up)
{

	//ARM_Vector& sched = GenCDOSchedule();

	int lup=MIN(floor(tranche_up/LossUnit),nbnames);
	int ldown=MIN(floor(tranche_down/LossUnit),nbnames);

	if ((itsCombinations.Getnbrows()==0) || (itsCombinations.Getnbrows()<MAX(nbnames+1,lup+1)))
	{
		itsCombinations.Resize(MAX(nbnames+1,lup+1),MAX(nbnames+1,lup+1));
		itsCombinations.GenCombinations();
	}

	vector<double> tmpPrbLoss;tmpPrbLoss.resize(lup+1);
	double CumtmpPrbLoss=0.,Loss = 0.,Cnp = 0;
	int i=0;
	double value1=0.,value2=0.;

//	for (int noperiod=0;noperiod<period;noperiod++)
	{
		for (int jump=0;jump<=period;jump++)
		{

			double prob = (value1 = itsTree->GetTree().data(period,jump).theta) +
							(value2 = itsTree->GetTree().data(period,jump).jumpsize);
			CumtmpPrbLoss = Loss =0.;

			for (i=0;i<=lup;i++)
			{
			Cnp = itsCombinations(nbnames,i);
			tmpPrbLoss[i]=pow(prob,i)*pow(1.-prob,nbnames-i);
			tmpPrbLoss[i]*=Cnp;
			CumtmpPrbLoss+=tmpPrbLoss[i];
			}

			for (i=ldown+1;i<=lup;i++)
			{	Loss+=tmpPrbLoss[i]*(i*LossUnit-tranche_down); }
			Loss += (1.-CumtmpPrbLoss)*(tranche_up-tranche_down);

			Loss /= (tranche_up-tranche_down);

			itsTree->GetTree().data(period,jump).EP_i_j = Loss;
		}
	}	
}

// ---------------------------------------------------------------------------------
// Calcul de la distribution d'EL
// ---------------------------------------------------------------------------------
void ICM_Pricer_tree_hw_cdo::CptExpectedLossTranche()
{

	BuildTree();
	itsTree->DiffuseTree(true);

	if (!getDistribLossObj().empty()) return; 

	ICM_ModelMultiCurves* model = (ICM_ModelMultiCurves*) GetModel(); 
	ICM_Ftd& ftd=dynamic_cast<ICM_Ftd&>(*GetSecurity()); 
	ICM_Collateral*collat =ftd.GetCollateral(); 
	const ICM_LossUnits& lossUnits= getLossUnits(); 

	ARM_Vector& sched = GenCDOSchedule(false);
	ICM_DistribLoss  DistribLoss;
	double ELT = 0.;

	for (int i=0;i<sched.size();i++)
	{
		double p =itsTree->SumProba(i);
		if (CHECK_EQUAL(p,1)==false)
		{double val=0.;}


		CalibrateTheta(i);
		double yf = sched.Elt(i);

		ARM_Date date = (ARM_Date)(model->GetStartDate().GetJulian() + K_YEAR_LEN*yf);
		double tranche_down = GetTranche_Down(date);
		double tranche_up = GetTranche_Up(date);

		tranche_down = ABS(round(tranche_down));
		tranche_up = ABS(round(tranche_up));

		HomogEL(i,collat->GetNbIssuers(),lossUnits.getLossUnit(date),tranche_down,tranche_up);
		ELT = itsTree->CorpEsp(i,false);

		if (ELT<0.) ELT = 0.;

		DistribLoss[yf]= ELT;
	}

	setDistribLoss(DistribLoss);
}

// ---------------------------------------------------------------------------------
// Calcul de la distribution d'EL
// ---------------------------------------------------------------------------------
void ICM_Pricer_tree_hw_cdo::GeneratePaths(vector<double>&	pdefaults,
											vector<double>&	spreads,
											vector<int>& namesindefault,
											vector<double>&  cumulativeloss)
{

	ICM_ModelMultiCurves* model = (ICM_ModelMultiCurves*) GetModel(); 
	ICM_Ftd& ftd=dynamic_cast<ICM_Ftd&>(*GetSecurity()); 
	ICM_Collateral*collat =ftd.GetCollateral(); 

	int jump=0,nbissuers=collat->GetNbIssuers();
	double lambda = 0.,rand = 0.,deltaT=0.,value1=0.,value2=0.;
	set<ORDER_> mOrder;

	ARM_Vector& sched = GenCDOSchedule(false);

	pdefaults.clear();pdefaults.resize(sched.size());
	spreads.clear();spreads.resize(sched.size());
	namesindefault.clear();
	cumulativeloss.clear();cumulativeloss.resize(sched.size());

	if (itsTree->GetTree().data(0,0).FwdSpread==CREDIT_DEFAULT_VALUE)
	{itsTree->EC_Node_until_maturity(0,0,pdefaults);
	spreads[0]=Pdef2Spread(sched.Elt(0),sched.Elt(sched.size()-1),pdefaults);
	itsTree->GetTree().data(0,0).FwdSpread = spreads[0];}
	else
	{spreads[0]=itsTree->GetTree().data(0,0).FwdSpread;}

	for (int period=0;period<sched.size()-1;period++)
	{
		lambda = itsTree->GetTree().data(period,jump).lambda;
		rand = its_mc_RandGen.GenerateOneRandom();
		deltaT= sched.Elt(period+1) - sched.Elt(period);

		if ( (lambda*deltaT)>rand ) 
		{jump += 1;}

		if (itsTree->GetTree().data(period+1,jump).FwdSpread==CREDIT_DEFAULT_VALUE)
		{itsTree->EC_Node_until_maturity(period+1,jump,pdefaults);
		spreads[period+1]=Pdef2Spread(sched.Elt(period+1),sched.Elt(sched.size()-1),pdefaults);
		itsTree->GetTree().data(period+1,jump).FwdSpread = spreads[period+1];
		}else
		{spreads[period+1] = itsTree->GetTree().data(period+1,jump).FwdSpread;}

		if (period)
		{cumulativeloss[period+1] = cumulativeloss[period];}

		if CHECK_NULL(pdefaults[period+1])
			continue;

		for (int nissuer=0; nissuer<nbissuers;nissuer++)
		{
			ORDER_ order(nissuer);
			if (mOrder.find(order) != mOrder.end())
				continue;

			rand = its_mc_RandGen.GenerateOneRandom();
			if ((pdefaults[period+1]-pdefaults[period])>rand)
			{
				mOrder.insert(order);
				namesindefault.push_back(nissuer);
				cumulativeloss[period+1] += (1.- model->GetRecoveryRate(collat->GetIssuersLabels(nissuer))) * 
											collat->GetIssuersNotional(nissuer);
			}

		}
	}

	#ifdef _DEBUG
	fprintf(stream_hw_tree,"\n");
	for (int y=0;y<spreads.size();y++)
	{fprintf(stream_hw_tree,"%lf,",spreads[y]);}
	#endif

}

// ---------------------------------------------------------------------------------
// OC Loss
// ---------------------------------------------------------------------------------
void ICM_Pricer_tree_hw_cdo::OCLoss(vector<double>&  cumulativeloss,
									 vector<double>&  cdoloss,
									 vector<double>&  spreads,
									 vector<double>&  pdefaults,
									 bool& trigger,
									 int& triggerBegin,
									 int& triggerend)
{

	cdoloss.clear();
	spreads.clear();
	ARM_Vector& sched = GenCDOSchedule(false);
	ICM_ModelMultiCurves* model = (ICM_ModelMultiCurves*) GetModel(); 
	double ELT = 0.,K1 = 0.,K2 = 0.;
	double ELT_inital = 0.,ELT_R = 0.,K1_R = 0.,K2_R = 0.;
	double OC=0.;
	trigger = false;

	triggerBegin = triggerend = sched.size()-1;

	ARM_Date PeriodDate;
	double yfmatu = sched[sched.size()-1]; 

	for (int i=0;i<sched.size() /* && (!trigger) */;i++)
	{
		double yf = sched.Elt(i);
		
		if ((yf>itsMaturitySub)&&(trigger==false))
			{break;}

		PeriodDate = (ARM_Date)(model->GetStartDate().GetJulian() + K_YEAR_LEN*yf);

		//standard cdo losses
		K1 = GetTranche_Down(PeriodDate);
		K2 = GetTranche_Up(PeriodDate);

		K1 = ABS(round(K1));
		K2 = ABS(round(K2));

		ELT_inital = ELT = MAX(cumulativeloss[i]-K1,0.) - MAX(cumulativeloss[i]-K2,0.);

		// -----------------------------------------------------------------
		//reset sub cdo losses
		//special for OC test
		// -----------------------------------------------------------------
		if (itsResetSub)
		{
			K1_R = K1 + itsResetSub;
			K2_R = K2 + itsResetSub;

			ELT_R = MAX(cumulativeloss[i]-K1_R,0.) - MAX(cumulativeloss[i]-K2_R,0.);
			ELT -= ELT_R; 
			//ELT = ELT_R; 

			//trigger = true;
			//triggerBegin = 0;

			if ((trigger==false)&&(cumulativeloss[i]< K1)&&(yf<itsMaturitySub))
			{
				//on calcul les spreads uniquement si l'ecart de loss des tranches est non nulle
				//if (/*ELT && */ (spreads.size()==0))
				//	{Pdef2Spreads(pdefaults,spreads);}

				if ((i<(sched.size()-1)) /*&& (ELT)*/)
				{ 
				OC = cumulativeloss[i]  + itsPtfLosse*spreads[i]* (yfmatu-yf) -K1;
				if (OC<0.) 
					{ELT = 0.;}	
				  else
					{
					 trigger = true;
					 triggerBegin = i;
					}
				}
			}

			if (trigger==false)
			{ELT = 0.;}
			else
			{
				ELT /= (K2 - K1); 
				if CHECK_EQUAL(K1,K2) {ELT=0.;}
				if (ELT<0.) ELT = 0.;
			}

			cdoloss.push_back(ELT);
		}
	}

//	if (triggerperiod<0)
//	{triggerperiod = cdoloss.size()-1;}
}


// ---------------------------------------------------------------------------------
// CDO Loss
// ---------------------------------------------------------------------------------
void ICM_Pricer_tree_hw_cdo::CDOLoss(vector<double>&  cumulativeloss,
									 vector<double>&  cdoloss,
									 vector<double>&  spreads,
									 vector<double>&  pdefaults,
									 bool& trigger,
									 int& triggerperiod)
{

	cdoloss.clear();
	spreads.clear();
	ARM_Vector& sched = GenCDOSchedule(false);
	ICM_ModelMultiCurves* model = (ICM_ModelMultiCurves*) GetModel(); 
	double ELT = 0.,K1 = 0.,K2 = 0.;
	double ELT_inital = 0.,ELT_R = 0.,K1_R = 0.,K2_R = 0.;
	double OC=0.;
	trigger = false;
	triggerperiod = -1;

	ARM_Date PeriodDate;
	double yfmatu = sched[sched.size()-1]; 

	for (int i=0;i<sched.size() /* && (!trigger) */;i++)
	{
		double yf = sched.Elt(i);
		
		PeriodDate = (ARM_Date)(model->GetStartDate().GetJulian() + K_YEAR_LEN*yf);

		//standard cdo losses
		K1 = GetTranche_Down(PeriodDate);
		K2 = GetTranche_Up(PeriodDate);

		K1 = ABS(round(K1));
		K2 = ABS(round(K2));

		ELT = MAX(cumulativeloss[i]-K1,0.) - MAX(cumulativeloss[i]-K2,0.);

		ELT /= (K2 - K1); 
		if CHECK_EQUAL(K1,K2) {ELT=0.;}
		if (ELT<0.) ELT = 0.;

		cdoloss.push_back(ELT);
	}

	if (triggerperiod<0)
	{triggerperiod = cdoloss.size()-1;}
}

// ---------------------------------------------------------------------------------
// CDO PayOff
// ---------------------------------------------------------------------------------
void ICM_Pricer_tree_hw_cdo::PayOffCDO(const vector<double>&  cdoloss,qCMPMETH mode,double& feepv,double& defpv,int end)
{
	ARM_Vector& sched = GenCDOSchedule(false);
	ICM_DistribLoss  DistribLoss;
	double ELT = 0.;
	ICM_ModelMultiCurves* model = (ICM_ModelMultiCurves*) GetModel(); 
	ICM_Ftd& ftd=dynamic_cast<ICM_Ftd&>(*GetSecurity());
	ICM_Leg* Feeleg = (ICM_Leg*) ftd.GetFeeLeg();
	ICM_Leg* Defleg = (ICM_Leg*) ftd.GetDefLeg();

	double S0 = Feeleg->GetCreditSpread()/100.;
	double notfeeleg = Feeleg->GetCreditInfos()->GetNotionals().Elt(0);
	double notdefleg = Defleg->GetCreditInfos()->GetNotionals().Elt(0);

/*	
	for (int i=0;i<sched.size();i++)
	{	double yf = sched.Elt(i);
		DistribLoss[yf]= cdoloss[i]; }

	setDistribLoss(DistribLoss);

	ICM_Pricer_Security::ResetRootPricer();
	if ((mode==qCMPACCRUED)||(mode==qCMPPRICE)||(mode==qCMPSPREAD)||(mode==qCMPFEELEGPV)||(mode==qCMPDEFLEGPV))
	{
		feepv = ICM_Pricer_Security::ComputePrice(qCMPFEELEGPV);
		defpv = ICM_Pricer_Security::ComputePrice(qCMPDEFLEGPV);
	}
*/
	
	if ((mode==qCMPACCRUED)||(mode==qCMPPRICE)||(mode==qCMPSPREAD)||(mode==qCMPFEELEGPV)||(mode==qCMPDEFLEGPV))
	{
		double NPV = FastCDOPricing(sched,0,end,S0,cdoloss,notfeeleg,notdefleg, 
					  model->GetZeroCurve(), feepv,defpv);

	}
}

void ICM_Pricer_tree_hw_cdo::PayOffCDOFwd(const vector<double>&  cdoloss,qCMPMETH mode,double& feepv,double& defpv,int begin,int end)
{
	ARM_Vector& sched = GenCDOSchedule(false);
	ICM_DistribLoss  DistribLoss;
	double ELT = 0.;
	ICM_ModelMultiCurves* model = (ICM_ModelMultiCurves*) GetModel(); 
	ICM_Ftd& ftd=dynamic_cast<ICM_Ftd&>(*GetSecurity());
	ICM_Leg* Feeleg = (ICM_Leg*) ftd.GetFeeLeg();
	ICM_Leg* Defleg = (ICM_Leg*) ftd.GetDefLeg();

	double S0 = Feeleg->GetCreditSpread()/100.;
	double notfeeleg = Feeleg->GetCreditInfos()->GetNotionals().Elt(0);
	double notdefleg = Defleg->GetCreditInfos()->GetNotionals().Elt(0);

	if ((mode==qCMPACCRUED)||(mode==qCMPPRICE)||(mode==qCMPSPREAD)||(mode==qCMPFEELEGPV)||(mode==qCMPDEFLEGPV))
	{
		double NPV = FastCDOPricing(sched, begin,end,S0,cdoloss,notfeeleg,notdefleg, 
					  model->GetZeroCurve(), feepv,defpv);

	}
}


// ---------------------------------------------------------------------------------
// Monte-Carlo Valuation
// ---------------------------------------------------------------------------------
double ICM_Pricer_tree_hw_cdo::MCPricing_OC(qCMPMETH mode)
{
	vector<double>	pdefaults;
	vector<double>	spreads;
	vector<int> namesindefault;
	vector<double>  cumulativeloss;
	vector<double>  cdoloss;
	bool trigger = false;
	int triggerBegin = 0.;
	int triggerEnd = 0.;

	itsPercentTrigger=0.;
	itsPercentLoss=0.;
	itsPercentLosstrigger=0.;
	double sumel=0.;

	double valuation = 0.;
	double feepv=0.,defpv=0.;
	double sum_feepv=0.,sum_defpv=0.;
	ICM_Ftd& ftd=dynamic_cast<ICM_Ftd&>(*GetSecurity());
	ICM_Leg* Feeleg = (ICM_Leg*) ftd.GetFeeLeg();
	double S0 = Feeleg->GetCreditSpread();

	for (int i=0;i<its_mc_nbsimuls;i++)
	{
		sumel=0.;
		trigger = false;
		GeneratePaths(pdefaults,spreads,namesindefault,cumulativeloss);
		OCLoss(cumulativeloss,cdoloss,spreads,pdefaults,trigger,triggerBegin,triggerEnd);
		PayOffCDOFwd(cdoloss,mode,feepv,defpv,triggerBegin,triggerEnd);		
		sum_feepv += feepv;
		sum_defpv += defpv;

		for (int j=0;j<cdoloss.size();j++)
		{sumel+=cdoloss[j];}

		if ((sumel) && (trigger))
		{itsPercentLosstrigger++;}

		if (sumel)
		{itsPercentLoss++;}

		if (trigger)
		{itsPercentTrigger++;}

	}

	sum_feepv /= (double)its_mc_nbsimuls;
	sum_defpv /= (double)its_mc_nbsimuls;

	itsPercentLoss/=(double)its_mc_nbsimuls;
	itsPercentLosstrigger/=(double)its_mc_nbsimuls;
	itsPercentTrigger/=(double)its_mc_nbsimuls;

	setValue(qCMPFEELEGPV,sum_feepv);
	setValue(qCMPDEFLEGPV,sum_defpv);
	setValue(qCMPPRICE,sum_feepv-sum_defpv);
	setValue(qCMPSPREAD,fabs(sum_defpv/sum_feepv*S0));
	setValue(qCMPACCRUED,0.);

	valuation = getValue(mode);
	return valuation;
}

// ---------------------------------------------------------------------------------
// Monte-Carlo Fast Valuation
// ---------------------------------------------------------------------------------
double ICM_Pricer_tree_hw_cdo::MCPricing_CDO(qCMPMETH mode)
{
	vector<double>	pdefaults;
	vector<double>	spreads;
	vector<int> namesindefault;
	vector<double>  cumulativeloss;
	vector<double>  _cumloss;
	vector<double>  cdoloss;
	bool trigger = false;
	itsPercentTrigger=0.;
	itsPercentLoss=0.;
	itsPercentLosstrigger=0.;
	double sumel=0.;
	int triggerperiod = 0.;

	int j=0;
	double valuation = 0.;
	double feepv=0.,defpv=0.;
	double sum_feepv=0.,sum_defpv=0.;
	ICM_Ftd& ftd=dynamic_cast<ICM_Ftd&>(*GetSecurity());
	ICM_Leg* Feeleg = (ICM_Leg*) ftd.GetFeeLeg();
	double S0 = Feeleg->GetCreditSpread()/100.;

	for (int i=0;i<its_mc_nbsimuls;i++)
	{
		sumel=0.;
		trigger = false;
		GeneratePaths(pdefaults,spreads,namesindefault,cumulativeloss);
		CDOLoss(cumulativeloss,cdoloss,spreads,pdefaults,trigger,triggerperiod);

		for (int j=0;j<cdoloss.size();j++)
		{sumel+=cdoloss[j];}

		if ((sumel) && (trigger))
		{itsPercentLosstrigger++;}

		if (sumel)
		{itsPercentLoss++;}

		if (trigger)
		{itsPercentTrigger++;}

		if (_cumloss.size()==0)
		{_cumloss = cdoloss;}
		else
		{	for (j=0;j<_cumloss.size();j++)
			{_cumloss[j]+=cdoloss[j];}}
	}

	for (j=0;j<_cumloss.size();j++)
	{_cumloss[j]/=(double)its_mc_nbsimuls;}

	itsPercentLoss/=(double)its_mc_nbsimuls;
	itsPercentLosstrigger/=(double)its_mc_nbsimuls;
	itsPercentTrigger/=(double)its_mc_nbsimuls;

	PayOffCDO(_cumloss,mode,feepv,defpv,triggerperiod);		

	setValue(qCMPFEELEGPV,feepv);
	setValue(qCMPDEFLEGPV,defpv);
	setValue(qCMPPRICE,feepv-defpv);
	setValue(qCMPSPREAD,fabs(defpv/feepv*S0));
	setValue(qCMPACCRUED,0.);

	valuation = getValue(mode);
	return valuation;
}


// ---------------------------------------------------------------------------------
// Pricing de la call option
// ---------------------------------------------------------------------------------
double ICM_Pricer_tree_hw_cdo::CallOptionPricing(qCMPMETH mode)
{
	ARM_Vector& sched = GenCDOSchedule(false);
	vector<double> pdefaults;pdefaults.resize(sched.size());
	double FwdSpread =0., esp=0.,V=0.,Duration=0.,Duration0=0.;
	int i,j;

	int idx = FindIndexInVector(sched,itsExerciceDate);

	for (j=0;j<=idx;j++)
	{
		itsTree->EC_Node_until_maturity(idx,j,pdefaults);
		FwdSpread=Pdef2Spread(sched.Elt(idx),sched.Elt(sched.size()-1),pdefaults);
		Duration=Pdef2Spread(sched.Elt(idx),sched.Elt(sched.size()-1),pdefaults,false);
		itsTree->GetTree().data(idx,j).FwdSpread = FwdSpread;
		itsTree->GetTree().data(idx,j).duration = Duration;
	}

	itsTree->EC_Node_until_maturity(0,0,pdefaults);
	Duration0=Pdef2Spread(sched.Elt(idx),sched.Elt(sched.size()-1),pdefaults,false);

	for (j=0;j<=idx;j++)
	{
		FwdSpread=itsTree->GetTree().data(idx,j).FwdSpread;
		Duration=itsTree->GetTree().data(idx,j).duration;
		esp += Duration * MAX(FwdSpread-itsStrike,0.) * itsTree->GetTree().data(idx,j).proba;
	}

//	esp /= Duration0;

	setValue(qCMPFEELEGPV,0.);
	setValue(qCMPDEFLEGPV,0.);
	setValue(qCMPPRICE,esp);
	setValue(qCMPSPREAD,0.);
	setValue(qCMPACCRUED,0.);

	double valuation = getValue(mode);
	return valuation;

}

// ---------------------------------------------------------------------------------
// Generic Valuation
// ---------------------------------------------------------------------------------
double ICM_Pricer_tree_hw_cdo::ComputePrice(qCMPMETH mode)
{
	if (getFlg(mode))
	{return getValue(mode);}

	#ifdef _DEBUG
		stream_hw_tree = fopen("c:\\temp\\stream_hw_tree.txt", "w+");
	#endif

	ICM_Ftd& ftd=dynamic_cast<ICM_Ftd&>(*GetSecurity()); 
	ICM_Collateral*collat =ftd.GetCollateral(); 
	ICM_ModelMultiCurves* model = (ICM_ModelMultiCurves*) GetModel(); 
	const ICM_DefaultCurve* DefCurve = model->GetDefaultCurve(collat->GetIssuersLabels(0));

	itsPtfLosse = /*(1. - DefCurve->GetRecovery()) * */ collat->SumNotionals(ftd.GetStartDateNA());
	itsResetSub = itsResetSubPercent * collat->SumNotionals(ftd.GetStartDateNA());

	double value = ICM_Pricer::ComputePrice(mode);

	if (its_mc_pricing_activate==1.)
	{	value = MCPricing_OC(mode); }
	else if (its_mc_pricing_activate==2.)
	{	value = MCPricing_CDO(mode);}
	else if (its_mc_pricing_activate==3.)
	{	value = CallOptionPricing(mode);}


	#ifdef _DEBUG
		fclose(stream_hw_tree);
	#endif

	return value;
}

