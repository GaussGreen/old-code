#include "ARMKernel/glob/firsttoinc.h" 
#include "ICMKernel\crv\icm_implied_loss_tree.h"
#include "ICMKernel\glob\icm_index_correlation.h"
#include "ICMKernel\inst\icm_mez.h"
#include "ICMKernel\pricer\icm_pricer_adviser.h"
//#include "ICMKernel\util\icm_pgcd.h"
//#include "ICMKernel\util\icm_RootFinder1D.h"
#include "ICMKernel/pricer/icm_pricer_homogeneous_smile.h"
#include "ICMKernel\mod\icm_lossunits.h"
#include "ICMKernel\mod\modelmulticurves.h"
#include "ICMKernel\crv\icm_volinterpol.h"
//#include "ARM_Kernel\util\interpol.h"

#include <nag.h>
#include <nage01.h>
#include <nage02.h>
#include <nag_stdlib.h>

class NagSplineHolder 
{ 
public:
	NagSplineHolder(); 
	~NagSplineHolder(); 
	void resize(unsigned int size); 
	Nag_Spline*  get(unsigned int i)   ;
private:
	std::vector<Nag_Spline> itsSplines; 
private:
	NagSplineHolder(const NagSplineHolder&ref); 
	NagSplineHolder & operator=(const NagSplineHolder &); 
} ;
//	-------------------------------------------------------------
NagSplineHolder::NagSplineHolder()
{}
//	-------------------------------------------------------------
NagSplineHolder::~NagSplineHolder()
{
	for (int i=0;i<itsSplines.size();i++)
	{
		NAG_FREE(itsSplines[i].lamda);
		NAG_FREE(itsSplines[i].c);
	}
}
//	-------------------------------------------------------------
void
NagSplineHolder::resize(unsigned int size)
{
	for (int i=0;i<itsSplines.size();i++)
	{
		NAG_FREE(itsSplines[i].lamda);
		NAG_FREE(itsSplines[i].c);
	}
	itsSplines.resize(size); 
}
//	-------------------------------------------------------------
Nag_Spline*  
NagSplineHolder::get(unsigned int i)  
{ 
	if (i>=itsSplines.size()) 
		ICMTHROW(ERR_INVALID_ARGUMENT,"NagSplineHolder: can't access "<<i); 
	return &itsSplines[i]; 
}
//	-------------------------------------------------------------

void ICM_ImpLossTree::Init()
{
	SetName(ICM_IMPLIED_LOSS_TREE);
	itsIdxDefaultCurve = NULL;
	itsCorrelation = NULL;
	itsTreeVector.clear();
	itsCorreltype = qCAL_BASE_CORRELATION;
	itsELoss.clear();
	itsPricesBefore.clear();
	itsPricesAfter.clear();
	itsELossAfter.clear();
	itsLossUnit=0.;

	its_cal_noindex=-1;
	its_cal_slice=-1;
	its_cal_noloss=-1;
	its_IsCalibrate = false;
	its_IsAdjusted = true;
	itsStep=30;

	
	// for(int i=0;i<itsSpleens.size();i++) 
	// 	delete itsSpleens[i]; 
	if (!itsSpleens) delete itsSpleens; 
	itsSpleens=new NagSplineHolder; 
	itsStrikes.clear();
	itsSchedCDO.clear();

};


ICM_ImpLossTree::~ICM_ImpLossTree()
{
/**	for (int i=0;i<itsSpleens.size();i++)
	{
		NAG_FREE(itsSpleens[i]->lamda);
		NAG_FREE(itsSpleens[i]->c);
		delete itsSpleens[i] ;
	} **/ 
	if (itsSpleens) delete itsSpleens ;
}

//---------------------------------------------------------------------------
// Compute cumulative loss
//---------------------------------------------------------------------------
double ICM_ImpLossTree::CptCumEL(int IndMatu,
						 int dw_IndLoss,
						 int up_IndLoss,
						 int noindex,
						 bool fwdstatus,
						 double notional_dw,
						 double notional_up)
{	
	double EL = 0.;
	//double value = 0.;
	for (int k=0;k<=IndMatu;k++)
	{

		if (itsTreeVector[noindex].data(IndMatu,k).state)
		{

		if ((notional_dw) || (notional_up))
		{
		if (!fwdstatus)
		{
			if ((k*itsLossUnit>=notional_dw)&&(k*itsLossUnit<=notional_up))
				EL+=((double)(k*itsLossUnit-notional_dw))*itsTreeVector[noindex].data(IndMatu,k).state;
			else if (k*itsLossUnit>notional_up)
				EL+=((double)(notional_up-notional_dw))*itsTreeVector[noindex].data(IndMatu,k).state;
		}
		else
		{
			if ((k*itsLossUnit>=notional_dw)&&(k*itsLossUnit<=notional_up))
				EL+=((double)(k*itsLossUnit-notional_dw))*itsTreeVector[noindex].data(IndMatu,k).fwd_state;
			else if (k*itsLossUnit>notional_up)
				EL+=((double)(notional_up-notional_dw))*itsTreeVector[noindex].data(IndMatu,k).fwd_state;
		}

		}		
		else
		{
		if (!fwdstatus)
		{
			if ((dw_IndLoss<=k) && (k<=up_IndLoss))
				EL+=((double)(k-dw_IndLoss))*itsTreeVector[noindex].data(IndMatu,k).state;
			else if (k>up_IndLoss)
				EL+=((double)(up_IndLoss-dw_IndLoss))*itsTreeVector[noindex].data(IndMatu,k).state;
		}
		else
		{
			if ((dw_IndLoss<=k) && (k<=up_IndLoss))
				EL+=((double)(k-dw_IndLoss))*itsTreeVector[noindex].data(IndMatu,k).fwd_state;
			else if (k>up_IndLoss)
				EL+=((double)(up_IndLoss-dw_IndLoss))*itsTreeVector[noindex].data(IndMatu,k).fwd_state;
		}

		}

		}
	}


	if ((notional_dw)||(notional_up))
	{
		EL /= (notional_up-notional_dw);
	}	
	else
	{
		EL /= (up_IndLoss-dw_IndLoss);
	}

	return (EL);
}

//---------------------------------------------------------------------------
// Compute initial implicit losses
//---------------------------------------------------------------------------
void ICM_ImpLossTree::CptImplicitEL(double Matu_,bool tree)
{
	itsELoss.clear();
	itsPricesBefore.clear();

	int nbslices = itsCorrelation->GetSlices().size();
	ARM_Date AsOf = itsIdxDefaultCurve->GetAsOfDate();
	ICM_VolInterpol* basecorrel=NULL;

	ICM_Parameters Params;
	ICM_Parameters ParamsTS;

	ARM_CLASS_NAME cname;

	//definition des parametres du pricer
	ARM_Vector PIntegrationStep(1,60);
	ARM_Vector PCopula(1,1.);
	ARM_Vector PIntegrationStep2(1,0.);
	ARM_Vector PFreedomDeg(1,0.);
	ARM_Vector TERMresc(1,1.);
	Params.Push(&PIntegrationStep,"INTEGRATION_STEP_1");
	Params.Push(&PCopula,"COPULA");
	Params.Push(&PIntegrationStep2,"INTEGRATION_STEP_2");
	Params.Push(&PFreedomDeg,"FREEDOM_DEGREE");

	ParamsTS =  Params ;
	ParamsTS.Push(&TERMresc,"TERMS_RESCALING");

	cname = ICM_PRICER_HOMOGENEOUS_SMILE;

	//creéation du modele
	std::vector<const ICM_DefaultCurve*> DefaultCurves(1); 
	DefaultCurves[0] = itsIdxDefaultCurve;

	ICM_ModelMultiCurves mmcTS(DefaultCurves,itsIdxDefaultCurve->GetZeroCurve(),NULL,itsCorrelation);

	//pour chacun des indices
	for (int i=0;i<nbslices;i++)
	{
		basecorrel= (ICM_VolInterpol*) itsCorrelation->GetSlices()[i].GetVolCurve();

		//on récupére les maturités standards
		ARM_Vector YFMaturities = *(*basecorrel).GetExpiryTerms();
		ARM_Vector Maturities = *(*basecorrel).GetExpiryTerms();
		Maturities *= 365.;
		Maturities += AsOf.GetJulian();

		//on récupére les strikes standards
		ARM_Vector Strikes = *(*basecorrel).GetStrikes();
		Strikes /= 100.;

		ARM_Date Start = AsOf.AddDays(1);							//départ à Asof+1j
		ARM_Date Maturity = Maturities.Elt(Maturities.size()-1);	//on choisit la maturité finale

		if ((itsCorreltype==qCAL_BASE_CORRELATION)&&(Matu_)) 
			{Maturity=Matu_;}

		// char* Collateral[TRX_EUR_NBNAMES];
		// for (int k=0; k<TRX_EUR_NBNAMES; k++)
		// {Collateral[k]=(char*)itsIdxDefaultCurve->GetLabel().c_str();}
		std::vector<std::string> Collateral(TRX_EUR_NBNAMES); 
		for (int k=0; k<TRX_EUR_NBNAMES; k++)
			Collateral[k]=itsIdxDefaultCurve->GetLabel(); 


		double strike_dw=0.;
		double strike_up=0.03;
		double amount = 1.e7;
		double amountfixed = 1.e7;

		ICM_Mez* indexCDO = CdoIndexDefinition(Start,Maturity, Collateral,
										// TRX_EUR_NBNAMES,
										0.,strike_up,
										0.01,amount,INCLUDE_MATURITY,false,
										DEFAULT_CREDIT_LAG_INDX,itsIdxDefaultCurve->GetCurrency() );

		int numflows=indexCDO->GetFeeLeg()->GetNumFlows();

		ICM_Pricer_Advisor Advisor;
		ICM_Pricer_Distrib_Smile* pricer=NULL;

		if (itsCorreltype==qCAL_BASE_CORRELATION_TS)
			pricer = (ICM_Pricer_Distrib_Smile*) Advisor.GeneratePricer(indexCDO,&mmcTS,cname,CREDIT_DEFAULT_VALUE,&ParamsTS,AsOf );
		else if (itsCorreltype==qCAL_BASE_CORRELATION)
			pricer = (ICM_Pricer_Distrib_Smile*) Advisor.GeneratePricer(indexCDO,&mmcTS,cname,CREDIT_DEFAULT_VALUE,&Params,AsOf );

		double price=pricer->Price(qCMPPRICE);

		itsLossUnit = pricer->getLossUnits().getLossUnit(AsOf);
		//nombre de losses du portefeuille
		//int sizelossunit = pgcd::u_round(((1-TRX_EUR_RECOVERY)*TRX_EUR_NBNAMES*amountfixed)/itsLossUnit);
		int sizelossunit = TRX_EUR_NBNAMES;

		//definition de la matrice d'expected loss
		double step = ((double)itsStep)/365.;
		int step_days = 30;
		double YFMaturity = (indexCDO->GetFeeLeg()->GetMaturity()-AsOf)/365.;
		int sizesched = (int)(YFMaturity/step)+1;

		const ARM_Vector& Vyfstart = indexCDO->GetFeeLeg()->GetCreditInfos()->GetYFAccStartDates();
		double yfend = indexCDO->GetFeeLeg()->GetCreditInfos()->GetYFAccEndDates().Elt(numflows-1);

		ARM_Vector* FinalSched1 = NULL;
		ARM_Vector* FinalSched = NULL;
		ARM_Vector SchedSynth(sizesched);
		for (int l=0;l<sizesched;l++) {SchedSynth.Elt(l) = (double)l*step;}
		MergeDates(&FinalSched1,&SchedSynth,&YFMaturities);
		MergeDates(&FinalSched,FinalSched1,&Vyfstart);

		sizesched=FinalSched->GetSize();
		itsTimeStep.Resize(1,(FinalSched)->GetSize());

		for (l=1;l<(FinalSched)->GetSize();l++)
		{itsTimeStep(0,l) = (FinalSched)->Elt(l-1);}

		itsTimeStep(0,0)=0.;

		if (FinalSched) delete FinalSched;
		FinalSched=NULL;

		int nbrows = sizelossunit+1;

		ICM_QMatrix<double> MatEl(nbrows,sizesched);
		ICM_QMatrix<double> PricesBefore(1.,nbrows);

		if (indexCDO) delete indexCDO;indexCDO=NULL;
		if (pricer) delete pricer;pricer=NULL;

		ICM_Smile_Correlation* correl = NULL;

		for (int j=1;j<nbrows;j++)
		{

			strike_dw = 0.;
			strike_up = ((double)(j))*itsLossUnit/(TRX_EUR_NBNAMES*amountfixed);

			if (itsCorreltype==qCAL_BASE_CORRELATION)
			{if (Matu_) {Maturity=Matu_;}}

			ICM_Mez* indexCDO = CdoIndexDefinition(Start,Maturity, Collateral,
										// TRX_EUR_NBNAMES,
										0.,strike_up,
										0.01,amount,INCLUDE_MATURITY,false,
										DEFAULT_CREDIT_LAG_INDX, itsIdxDefaultCurve->GetCurrency());

			
			correl = (ICM_Smile_Correlation*) (mmcTS.GetCorrelation());
			vector<double> strikedw;strikedw.push_back(0.);
			correl->GetSlices()[i].SetSmileStrikeLow(strikedw);
			vector<double> strikeup;strikeup.push_back(strike_up);
			correl->GetSlices()[i].SetSmileStrikeHigh(strikeup);

			if (itsCorreltype==qCAL_BASE_CORRELATION_TS)
				pricer = (ICM_Pricer_Distrib_Smile*) Advisor.GeneratePricer(indexCDO,&mmcTS,cname,CREDIT_DEFAULT_VALUE,&ParamsTS,AsOf);
			else if (itsCorreltype==qCAL_BASE_CORRELATION)
				pricer = (ICM_Pricer_Distrib_Smile*) Advisor.GeneratePricer(indexCDO,&mmcTS,cname,CREDIT_DEFAULT_VALUE,&Params,AsOf);

			double price=pricer->Price(qCMPPRICE);
			PricesBefore(0,j)=price;

			double test = 0.;
			for (int k=0;k<sizesched;k++)
			{
				MatEl(j,k)= pricer->getDistribLoss().InterpolEL(itsTimeStep(0,k));
				
				if (tree)
					{MatEl(j,k) *= itsLossUnit*((double)j);}

			}

			if (indexCDO) delete indexCDO;indexCDO=NULL;
			if (pricer) delete pricer;pricer=NULL;
		}

		itsELoss.push_back(MatEl);
		itsPricesBefore.push_back(PricesBefore);
	}
}

//---------------------------------------------------------------------------
// Build Trees
//---------------------------------------------------------------------------
void ICM_ImpLossTree::CreateTrees()
{
	itsTreeVector.clear();
	int nbslices = itsCorrelation->GetSlices().size();
	int i=0,j=0,k=0;

	ICM_QMatrix<double> disctimes(itsTimeStep.Getnbcols(),1);

	for (i=0;i<itsTimeStep.Getnbcols();i++)
	{disctimes(i,0)=itsTimeStep(0,i);}

	//pour chacun des indices
	for (k=0;k<nbslices;k++)
	{
		BinomialTree<ImpliedNode> Tree;
		Tree.setTimes(disctimes);
		itsTreeVector.push_back(Tree);
	}
}

//---------------------------------------------------------------------------
// compute state Losses
//---------------------------------------------------------------------------
void ICM_ImpLossTree::CptStateLosses(bool adjust)
{
	int i=0,j=0,k=0;

	//pour chacun des indices
	for (k=0;k<itsTreeVector.size();k++)
	{
		for (int i=0;i<itsELoss[k].Getnbcols();i++)
			for (int j=0;j<=MIN(itsELoss[k].Getnbrows()-1,i);j++)
			{
				if (j==itsELoss[k].Getnbrows()-1)
					itsTreeVector[k].data(i,j).state=
						((itsELoss[k])(j,i)-(itsELoss[k])(j-1,i))/itsLossUnit;
				else if (j)
					itsTreeVector[k].data(i,j).state=
						(2.*(itsELoss[k])(j,i)-(itsELoss[k])(j+1,i)-(itsELoss[k])(j-1,i))/itsLossUnit;
				else
					itsTreeVector[k].data(i,j).state= 1. - (itsELoss[k])(1,i)/itsLossUnit;


				//ajustement pour les valeurs invalides
				if ((j>0)&&(its_IsAdjusted))
				{
				if ((itsTreeVector[k].data(i,j).state>=1.) ||
					(itsTreeVector[k].data(i,j).state<=0.))
					//(itsTreeVector[k].data(i,j).state>itsTreeVector[k].data(i,j-1).state))
				{itsTreeVector[k].data(i,j).state = itsTreeVector[k].data(i,j-1).state/10.;}
				}
			}

	}
}

//---------------------------------------------------------------------------
// compute transitions probabilities for node
//---------------------------------------------------------------------------
void ICM_ImpLossTree::CptTransitionProbaForNode(const int noindex,
												const int slice,
												const int noloss,
												bool adjust)
{
	if (noloss>slice) 
		return;

	double value=0.;
	if (noloss)
	{
		value = itsTreeVector[noindex].data(slice,noloss).p_nodef=
				(itsTreeVector[noindex].data(slice+1,noloss).state -
				itsTreeVector[noindex].data(slice,noloss-1).state * itsTreeVector[noindex].data(slice,noloss-1).p_def) /
				itsTreeVector[noindex].data(slice,noloss).state;

		if (adjust) 
		{
			if ((itsTreeVector[noindex].data(slice,noloss).p_nodef<0.) ||
				(itsTreeVector[noindex].data(slice,noloss).p_nodef>1.))
				itsTreeVector[noindex].data(slice,noloss).p_nodef = itsTreeVector[noindex].data(slice,noloss-1).p_nodef;
		}

		itsTreeVector[noindex].data(slice,noloss).p_def=1. - itsTreeVector[noindex].data(slice,noloss).p_nodef;

	}
	else
	{
		itsTreeVector[noindex].data(slice,0).p_nodef = 
			itsTreeVector[noindex].data(slice+1,0).state /
			itsTreeVector[noindex].data(slice,0).state;


		if (adjust) 
		{
			if ((itsTreeVector[noindex].data(slice,0).p_nodef<0.) ||
				(itsTreeVector[noindex].data(slice,0).p_nodef>1.))
				itsTreeVector[noindex].data(slice,0).p_nodef = 1.e-5;
		}

		itsTreeVector[noindex].data(slice,0).p_def=1. - itsTreeVector[noindex].data(slice,0).p_nodef;
	}

}


//---------------------------------------------------------------------------
// deduce state losses from transitions probabilities for node
//---------------------------------------------------------------------------
void ICM_ImpLossTree::CptStateLossesFromProbaForNode(const int noindex,
													 const int slice,
													 const int noloss,
													 const double proba)
{
	if (noloss>slice) 
		return;

	double value=0.;

	itsTreeVector[noindex].data(slice,noloss).p_nodef=proba;
	itsTreeVector[noindex].data(slice,noloss).p_def=1. - itsTreeVector[noindex].data(slice,noloss).p_nodef;

	value=itsTreeVector[noindex].data(slice,noloss).state=
		(itsTreeVector[noindex].data(slice+1,noloss).state -
		itsTreeVector[noindex].data(slice,noloss-1).state * itsTreeVector[noindex].data(slice,noloss-1).p_def) /
		itsTreeVector[noindex].data(slice,noloss).p_nodef;

}

//----------------------------------------------------------------------------------------------
// deduce state losses from transitions probabilities for node with given transition probability
//----------------------------------------------------------------------------------------------
void ICM_ImpLossTree::CptNextStateLossesFromProbaForNode(const int noindex,
														const int slice,
														const int noloss,
														const double proba)
{
	if (noloss>slice) 
		return;

	double value=0.;

	itsTreeVector[noindex].data(slice,noloss).p_nodef=proba;
	itsTreeVector[noindex].data(slice,noloss).p_def=1. - itsTreeVector[noindex].data(slice,noloss).p_nodef;

	value = itsTreeVector[noindex].data(slice+1,noloss).state=
			itsTreeVector[noindex].data(slice,noloss).p_nodef * 
			itsTreeVector[noindex].data(slice,noloss).state +
			itsTreeVector[noindex].data(slice,noloss-1).state *
			itsTreeVector[noindex].data(slice,noloss-1).p_def;
}


//-------------------------------------------------------------
// deduce state losses from transitions probabilities for node 
//-------------------------------------------------------------
void ICM_ImpLossTree::CptNextStateLossesFromProbaForNode(const int noindex,
														const int slice,
														const int noloss,
														bool fwd)
{
	if (noloss>slice) 
		return;

	if (fwd)
	{	if (noloss<1)
			itsTreeVector[noindex].data(slice+1,noloss).fwd_state=
				itsTreeVector[noindex].data(slice,noloss).p_nodef * 
				itsTreeVector[noindex].data(slice,noloss).fwd_state;
		else
		{
			itsTreeVector[noindex].data(slice+1,noloss).fwd_state=
				itsTreeVector[noindex].data(slice,noloss).p_nodef * 
				itsTreeVector[noindex].data(slice,noloss).fwd_state +
				itsTreeVector[noindex].data(slice,noloss-1).fwd_state *
				itsTreeVector[noindex].data(slice,noloss-1).p_def;

			if (its_IsAdjusted)
			{
				if ((itsTreeVector[noindex].data(slice+1,noloss).fwd_state>=1.) ||
					(itsTreeVector[noindex].data(slice+1,noloss).fwd_state<=0.))
					//(itsTreeVector[k].data(i,j).state>itsTreeVector[k].data(i,j-1).state))
				{itsTreeVector[noindex].data(slice+1,noloss).fwd_state = itsTreeVector[noindex].data(slice+1,noloss-1).fwd_state/10.;}
			}
		}
	}
	else
	{	if (noloss<1)
			itsTreeVector[noindex].data(slice+1,noloss).state=
				itsTreeVector[noindex].data(slice,noloss).p_nodef * 
				itsTreeVector[noindex].data(slice,noloss).state;
		else
		{	
			itsTreeVector[noindex].data(slice+1,noloss).state=
				itsTreeVector[noindex].data(slice,noloss).p_nodef * 
				itsTreeVector[noindex].data(slice,noloss).state +
				itsTreeVector[noindex].data(slice,noloss-1).state *
				itsTreeVector[noindex].data(slice,noloss-1).p_def;

			if (its_IsAdjusted)
			{
				if ((itsTreeVector[noindex].data(slice+1,noloss).state>=1.) ||
					(itsTreeVector[noindex].data(slice+1,noloss).state<=0.))
					//(itsTreeVector[k].data(i,j).state>itsTreeVector[k].data(i,j-1).state))
				{itsTreeVector[noindex].data(slice+1,noloss).state = itsTreeVector[noindex].data(slice+1,noloss-1).state/10.;}
			}
		}
	}
}

//---------------------------------------------------------------------------
// deduce state losses from transitions probabilities for node
//---------------------------------------------------------------------------
double ICM_ImpLossTree::CptNextStateLossesFromProbaForNodeRoot(double proba)
{
	int noindex=its_cal_noindex;
	int slice=its_cal_slice;
	int noloss=its_cal_noloss;

	if (noloss>slice) 
		return CREDIT_DEFAULT_VALUE;

	double value=0.;

	itsTreeVector[noindex].data(slice,noloss).p_nodef=proba;
	itsTreeVector[noindex].data(slice,noloss).p_def=1. - itsTreeVector[noindex].data(slice,noloss).p_nodef;

	value = itsTreeVector[noindex].data(slice+1,noloss).state=
			itsTreeVector[noindex].data(slice,noloss).p_nodef * 
			itsTreeVector[noindex].data(slice,noloss).state +
			itsTreeVector[noindex].data(slice,noloss-1).state *
			itsTreeVector[noindex].data(slice,noloss-1).p_def;

	return value;
}


double ICM_ImpLossTree::CptNextStateLossesFromProbaForNodeRootDF(double proba)
{
	int noindex=its_cal_noindex;
	int slice=its_cal_slice;
	int noloss=its_cal_noloss;

	return itsTreeVector[noindex].data(slice,noloss).state;
}


//---------------------------------------------------------------------------
// compute transitions probabilities for slice
//---------------------------------------------------------------------------
void ICM_ImpLossTree::CptTransitionProbaForSlice(const int noindex,
												 const int slice,
												 bool adjust)
{
	double value=0.;
	for (int noloss=0;noloss<=slice;noloss++)
		{CptTransitionProbaForNode(noindex,slice,noloss,adjust);}
}

//---------------------------------------------------------------------------
// compute transitions probabilities for trees
//---------------------------------------------------------------------------
void ICM_ImpLossTree::CptTransitionProba(bool adjust)
{
	//pour chacun des indices
	for (int k=0;k<itsTreeVector.size();k++)
	{
		for (int i=0;i<itsTreeVector[k].depth()-1;i++)
			{CptTransitionProbaForSlice(k,i,adjust);}

	}
}

//---------------------------------------------------------------------------
// Adjust transitions probabilities and calibrate Tree
//---------------------------------------------------------------------------
void ICM_ImpLossTree::CalibrateTree()
{

	//pour chacun des indices
	for (int noindex=0;noindex<itsTreeVector.size();noindex++)
	{
		double value=0.;

		for (int slice=0;slice<itsTreeVector[noindex].depth()-1;slice++)
			for (int noloss=0;noloss<=slice;noloss++)
			{
				//CptTransitionProbaForNode(noindex,slice,noloss);
				its_cal_noindex=noindex;
				its_cal_slice=slice;
				its_cal_noloss=noloss;

				if ((itsTreeVector[noindex].data(slice,noloss).p_nodef<0.)||
					(itsTreeVector[noindex].data(slice,noloss).p_nodef>1.))
				{
					double _inf = 0.;
					double _sup = 1.;
					int nb_iter=100;
					double proba=0.;
					double accuracy=1.E-5;

					/*
					//on recalcule les probabilités de transitions
					proba = RootFinder1D(
						ff1::mem_call(&ICM_ImpLossTree::CptNextStateLossesFromProbaForNodeRoot,(*this)),
						ff1::mem_call(&ICM_ImpLossTree::CptNextStateLossesFromProbaForNodeRootDF,(*this))
						).NewtonRaphson(_inf,_sup,nb_iter,accuracy); */

					//proba = MIN(1.,MAX(0,1./(2.*itsTreeVector[noindex].data(slice,noloss).state)));
					proba = itsTreeVector[noindex].data(slice,noloss-1).p_nodef;

					CptNextStateLossesFromProbaForNode(noindex,slice,noloss,proba);

					if ((itsTreeVector[noindex].data(slice+1,noloss).state >1.)||
						(itsTreeVector[noindex].data(slice+1,noloss).state<0.))
					{
						proba = 0.5;
						break;
					}
				}

			}

	}
}


//---------------------------------------------------------------------------
// compute transitions probabilities for slice
//---------------------------------------------------------------------------
void ICM_ImpLossTree::CptNextStateLossesForSlice(const int noindex,
												 const int slice,
												 bool fwd)
{
	double value=0.;
	for (int noloss=0;noloss<=slice;noloss++)
		{CptNextStateLossesFromProbaForNode(noindex,slice,noloss,fwd);}
}

//---------------------------------------------------------------------------
// compute transitions probabilities for trees
//---------------------------------------------------------------------------
void ICM_ImpLossTree::DiffuseStateLosses(int begin,int end,bool fwd)
{
	//pour chacun des indices
	for (int k=0;k<itsTreeVector.size();k++)
	{
		for (int i=begin;i<end;i++)
			{CptNextStateLossesForSlice(k,i,fwd);}

	}
}


//---------------------------------------------------------------------------
// Build new expected losses
//---------------------------------------------------------------------------
void ICM_ImpLossTree::CptPrices()
{
	itsPricesAfter.clear();
	itsELossAfter.clear();

	int nbslices = itsCorrelation->GetSlices().size();
	ARM_Date AsOf = itsIdxDefaultCurve->GetAsOfDate();
	ICM_VolInterpol* basecorrel=NULL;

	ICM_Parameters Params;
	ICM_Parameters ParamsTS;

	ARM_CLASS_NAME cname;

	//definition des parametres du pricer
	ARM_Vector PIntegrationStep(1,60);
	ARM_Vector PCopula(1,1.);
	ARM_Vector PIntegrationStep2(1,0.);
	ARM_Vector PFreedomDeg(1,0.);
	ARM_Vector TERMresc(1,1.);
	Params.Push(&PIntegrationStep,"INTEGRATION_STEP_1");
	Params.Push(&PCopula,"COPULA");
	Params.Push(&PIntegrationStep2,"INTEGRATION_STEP_2");
	Params.Push(&PFreedomDeg,"FREEDOM_DEGREE");

	// ParamsTS =  (ICM_Parameters*) Params.Clone();
	ParamsTS =Params; 
	ParamsTS.Push(&TERMresc,"TERMS_RESCALING");

	cname = ICM_PRICER_HOMOGENEOUS_SMILE;

	//creéation du modele
	std::vector<const ICM_DefaultCurve*> DefaultCurves(1); 
	DefaultCurves[0] = itsIdxDefaultCurve;

	ICM_ModelMultiCurves mmcTS(DefaultCurves,itsIdxDefaultCurve->GetZeroCurve(),NULL,itsCorrelation);

	//pour chacun des indices
	for (int i=0;i<nbslices;i++)
	{
		basecorrel= (ICM_VolInterpol*) itsCorrelation->GetSlices()[i].GetVolCurve();

		//on récupére les maturités standards
		ARM_Vector YFMaturities = *(*basecorrel).GetExpiryTerms();
		ARM_Vector Maturities = *(*basecorrel).GetExpiryTerms();
		Maturities *= 365.;
		Maturities += AsOf.GetJulian();

		//on récupére les strikes standards
		ARM_Vector Strikes = *(*basecorrel).GetStrikes();
		Strikes /= 100.;

		ARM_Date Start = AsOf.AddDays(1);							//départ à Asof+1j
		ARM_Date Maturity = Maturities.Elt(Maturities.size()-1);	//on choisit la maturité finale

		// char* Collateral[TRX_EUR_NBNAMES];
		// for (int k=0; k<TRX_EUR_NBNAMES; k++)
		// {Collateral[k]=(char*)itsIdxDefaultCurve->GetLabel().c_str();}
		std::vector<std::string> Collateral(TRX_EUR_NBNAMES); 
		for (int k=0; k<TRX_EUR_NBNAMES; k++)
			Collateral[k]=itsIdxDefaultCurve->GetLabel(); 

		double strike_dw=0.;
		double strike_up=0.03;
		double amount = 1.e7;
		double amountfixed = 1.e7;

		ICM_Mez* indexCDO = CdoIndexDefinition(Start,Maturity, Collateral,
										// TRX_EUR_NBNAMES,
										0.,strike_up,
										0.01,amount,INCLUDE_MATURITY,false,
										DEFAULT_CREDIT_LAG_INDX, itsIdxDefaultCurve->GetCurrency());

		int numflows=indexCDO->GetFeeLeg()->GetNumFlows();

		ICM_Pricer_Advisor Advisor;
		ICM_Pricer_Distrib_Smile*	pricer = NULL;

		if (itsCorreltype==qCAL_BASE_CORRELATION_TS)
			pricer = (ICM_Pricer_Distrib_Smile*) Advisor.GeneratePricer(indexCDO,&mmcTS,cname,CREDIT_DEFAULT_VALUE,&ParamsTS,AsOf);
		else if (itsCorreltype==qCAL_BASE_CORRELATION)
			pricer = (ICM_Pricer_Distrib_Smile*) Advisor.GeneratePricer(indexCDO,&mmcTS,cname,CREDIT_DEFAULT_VALUE,&Params,AsOf);

		double price=pricer->Price(qCMPPRICE);

		// itsLossUnit = pricer->GetLossUnit();
		itsLossUnit = pricer->getLossUnits().getLossUnit(AsOf);
		//nombre de losses du portefeuille
		int sizelossunit = round(((1-TRX_EUR_RECOVERY)*TRX_EUR_NBNAMES*amountfixed)/itsLossUnit);

		//definition de la matrice d'expected loss
		double step = 0.25;
		double YFMaturity = (indexCDO->GetFeeLeg()->GetMaturity()-AsOf)/365.;
		int nbrows = sizelossunit+1;

		ICM_QMatrix<double> PricesAfter(1.,nbrows);

		if (indexCDO) delete indexCDO;indexCDO=NULL;
		if (pricer) delete pricer;pricer=NULL;

		ICM_DistribLoss DistribLoss;

		ICM_QMatrix<double> MatEl(nbrows,itsTimeStep.Getnbcols());

		for (int j=1;j<nbrows;j++)
		{
			//strike_dw = (j-1)*lossunit/((1-TRX_EUR_RECOVERY)*TRX_EUR_NBNAMES*amountfixed);
			strike_dw = 0.;
			strike_up = ((double)j)*itsLossUnit/((1.-TRX_EUR_RECOVERY)*TRX_EUR_NBNAMES*amountfixed);

			ICM_Mez* indexCDO = CdoIndexDefinition(Start,Maturity, Collateral,
										// TRX_EUR_NBNAMES,
										0.,strike_up,
										0.01,amount,INCLUDE_MATURITY,false,
										DEFAULT_CREDIT_LAG_INDX, itsIdxDefaultCurve->GetCurrency());

		ICM_Pricer_Distrib_Smile* pricer=NULL;

		ICM_Smile_Correlation* correl = (ICM_Smile_Correlation*) (mmcTS.GetCorrelation());
		vector<double> strikedw;strikedw.push_back(0.);
		correl->GetSlices()[i].SetSmileStrikeLow(strikedw);
		vector<double> strikeup;strikeup.push_back(strike_up);
		correl->GetSlices()[i].SetSmileStrikeHigh(strikeup);

		if (itsCorreltype==qCAL_BASE_CORRELATION_TS)
			pricer = (ICM_Pricer_Distrib_Smile*) Advisor.GeneratePricer(indexCDO,&mmcTS,cname,CREDIT_DEFAULT_VALUE,&ParamsTS,AsOf);
		else if (itsCorreltype==qCAL_BASE_CORRELATION)
			pricer = (ICM_Pricer_Distrib_Smile*) Advisor.GeneratePricer(indexCDO,&mmcTS,cname,CREDIT_DEFAULT_VALUE,&Params,AsOf);


		double price=pricer->Price(qCMPPRICE);

		for (int k=0;k<itsTimeStep.Getnbcols();k++)
		{
			//DistribLoss[(itsTimeStep(0,k))]=j*itsTreeVector[i].data(k,j).state;
			MatEl(j,k)=CptCumEL(k,0,j,i);
			DistribLoss[(itsTimeStep(0,k))]=MatEl(j,k);
		}

		pricer->setDistribLoss(DistribLoss);
		pricer->ResetRootPricer();
		PricesAfter(0,j)=pricer->Price(qCMPPRICE);

		if (indexCDO) delete indexCDO;indexCDO=NULL;
		if (pricer) delete pricer;pricer=NULL;

		}

		itsPricesAfter.push_back(PricesAfter);
		itsELossAfter.push_back(MatEl);
	}

}

//---------------------------------------------------------------------------
// View method
//---------------------------------------------------------------------------
void ICM_ImpLossTree::View(char* id, FILE* ficOut)
{
	unsigned int i=0,j=0;
	FILE* fOut;
	char  fOutName[200];

	if ( ficOut == NULL )
	{
		ARM_GetViewFile(ARM_VIEW_FILE, id, fOutName);
		(void) unlink(fOutName);
		fOut = fopen(fOutName, "w"); 
	}
	else	fOut = ficOut;
	
	for (int k=0; k<itsTreeVector.size();k++)
	{

		//itsELoss[k].View(id,fOut);

		fprintf(fOut,"\n");
		fprintf(fOut,"\n");

		std::stringstream sstr ;
		sstr<<"\t\t\t ----------------- Binomial Tree ----------------- \n\n"<<std::endl ;
		
		sstr<<"-- Dumping Tree : (state,p_def,p_nodef) "<<std::endl ;
		for(i=0;i<itsTreeVector[k].depth();i++)
		{
			sstr<<((double)(int)(1.e3*itsTimeStep(0,i)))/1.e3<<"\t"; 
		for(j=0;j<=i;j++) 
			sstr<<"(" /*<<((double)(int)(1.e3*itsTreeVector[k].data(i,j).cumEL)/1.e3)*/
				<<((double)(int)(1.e3*itsTreeVector[k].data(i,j).state)/1.e3)
				<<","<<((double)(int)(1.e3*itsTreeVector[k].data(i,j).p_def)/1.e3)
				<<","<<((double)(int)(1.e3*itsTreeVector[k].data(i,j).p_nodef)/1.e3)<<")"<<"\t"; 

		sstr<<std::endl; 
		}
		sstr<<"-- End of Dumping Tree\n"<<std::endl ;
		sstr<<std::endl; 

		sstr<<"\t\t\t ----------------- Tree Proba default ----------------- \n\n"<<std::endl ;
		
		sstr<<"-- Dumping Tree : (p_def) "<<std::endl ;
		for(i=0;i<itsTreeVector[k].depth();i++)
		{
			sstr<<((double)(int)(1.e3*itsTimeStep(0,i)))/1.e3<<"\t"; 
		for(j=0;j<=i;j++) 
			sstr<<((double)(int)(1.e3*itsTreeVector[k].data(i,j).p_def)/1.e3)<<"\t"; 

		sstr<<std::endl; 
		}
		sstr<<"-- End of Dumping Tree\n"<<std::endl ;
		sstr<<std::endl; 

		sstr<<"\t\t\t ----------------- Tree Fwd States default ----------------- \n\n"<<std::endl ;
		
		sstr<<"-- Dumping Tree : (fwd_state) "<<std::endl ;
		for(i=0;i<itsTreeVector[k].depth();i++)
		{
			sstr<<((double)(int)(1.e3*itsTimeStep(0,i)))/1.e3<<"\t"; 
		for(j=0;j<=i;j++) 
			sstr<<((double)(int)(1.e3*itsTreeVector[k].data(i,j).fwd_state)/1.e3)<<"\t"; 

		sstr<<std::endl; 
		}
		sstr<<"-- End of Dumping Tree\n"<<std::endl ;
		sstr<<std::endl; 


		sstr<<"\t\t\t ----------------- Tree States ----------------- \n\n"<<std::endl ;
		
		sstr<<"-- Dumping Tree : (state) "<<std::endl ;
		for(i=0;i<itsTreeVector[k].depth();i++)
		{
			sstr<<((double)(int)(1.e3*itsTimeStep(0,i)))/1.e3<<"\t"; 
		for(j=0;j<=i;j++) 
			sstr<<((double)(int)(1.e3*itsTreeVector[k].data(i,j).state)/1.e3)<<"\t"; 

		sstr<<std::endl; 
		}
		sstr<<"-- End of Dumping Tree\n"<<std::endl ;
		sstr<<std::endl; 

		fprintf(fOut,"%s\n",sstr.str().c_str()); 
		fprintf(fOut, "\n");

		if (itsPricesBefore.size()>0)
		{fprintf(fOut,"Price Before N°%i\n",k);
		itsPricesBefore[k].View(id,fOut);}
		
		if (itsPricesAfter.size()>0)
		{fprintf(fOut,"Price After N°%i\n\n\n",k);
		itsPricesAfter[k].View(id,fOut);}

	}

	if ( ficOut == NULL )fclose(fOut);
}

//spleens 
void ICM_ImpLossTree::BuildSplines()
{
	vector <double> times;
	vector <double> strikes;
	vector <double> el;
	
	// itsSpleens.clear();

	int i=0,j=0,k=0,l=0;
	/** for(i=0;i<itsSpleens.size();i++) 
		delete itsSpleens[i]; 
	itsSpleens.clear();  **/ 

	int nbtimes=itsTimeStep.Getnbcols(); //itsELoss[0].Getnbcols();

	itsSpleens->resize(nbtimes);
	
	for (i=0;i<nbtimes;i++)
	{times.push_back(itsTimeStep(0,i));}

	//strikes.push_back(0.);
	for (j=0;j<itsStrikes.size();j++)
	{strikes.push_back(itsStrikes[j]);}

	for (l=0;l<nbtimes;l++)
	{
		el.clear();
		for (k=0;k<itsELoss[0].Getnbrows();k++)
		{el.push_back(itsELoss[0].Getvalue(k,l));}

// FIXMEFRED: mig.vc8 (28/05/2007 15:14:08):cast
		e01bac(el.size(),&(*strikes.begin()),&(*el.begin()),itsSpleens->get(l),NAGERR_DEFAULT);
	}

	its_IsCalibrate=true;
}

//interpolation des EL
double ICM_ImpLossTree::InterpolEL(double strike, double maturity)
{

	if (maturity<=1.e-20)
		return 0.;

	int exact = -1,_inf=-1,_sup=-1;
	int nbtimes=itsTimeStep.Getnbcols(); //itsELoss[0].Getnbcols();

	for (int i=0;i<nbtimes;i++)
	{
		if (CHECK_EQUAL(maturity,itsTimeStep(0,i)))
		{exact=i;break;}
		
		if (i<nbtimes-1){
		if ((itsTimeStep(0,i)<maturity) && (maturity<itsTimeStep(0,i+1)))
		{_inf=i;_sup=i+1;}}
	}

	if (_sup==-1)
	{_sup=nbtimes-1;}

	double fit=0.,fit1=0.,fit2=0.;

	if (exact>=0)
	{e02bbc(strike,&fit,itsSpleens->get(exact),NAGERR_DEFAULT);}
	else
	{
		e02bbc(strike,&fit1,itsSpleens->get(_inf),NAGERR_DEFAULT);
		e02bbc(strike,&fit2,itsSpleens->get(_sup),NAGERR_DEFAULT);
		double deltat = (itsTimeStep(0,_sup)-itsTimeStep(0,_inf));
		fit = (maturity-itsTimeStep(0,_inf))/deltat*fit1 + 
			(itsTimeStep(0,_sup)-maturity)/deltat*fit2;
	}

	return fit;
}

double ICM_ImpLossTree::CptCumELSpline(double strikedw,
									   double strikeup,
									   double maturity)
{
	double EL_dw = InterpolEL(strikedw, maturity);
	double EL_up = InterpolEL(strikeup, maturity);

	double EL = (EL_up*strikeup-EL_dw*strikedw)/(strikeup-strikedw); 
	return EL;
}	


//---------------------------------------------------------------------------
// Compute initial implicit losses
//---------------------------------------------------------------------------
void ICM_ImpLossTree::CptImplicitELforsplines(double Matu_)
{
	int l=0;
	itsELoss.clear();
	itsPricesBefore.clear();
	itsStrikes.clear();

	int nbslices = itsCorrelation->GetSlices().size();
	ARM_Date AsOf = itsIdxDefaultCurve->GetAsOfDate();
	ICM_VolInterpol* basecorrel=NULL;

	ICM_Parameters Params;
	ICM_Parameters ParamsTS;

	ARM_CLASS_NAME cname;

	//definition des parametres du pricer
	ARM_Vector PIntegrationStep(1,60);
	ARM_Vector PCopula(1,1.);
	ARM_Vector PIntegrationStep2(1,0.);
	ARM_Vector PFreedomDeg(1,0.);
	ARM_Vector TERMresc(1,1.);
	Params.Push(&PIntegrationStep,"INTEGRATION_STEP_1");
	Params.Push(&PCopula,"COPULA");
	Params.Push(&PIntegrationStep2,"INTEGRATION_STEP_2");
	Params.Push(&PFreedomDeg,"FREEDOM_DEGREE");

	ParamsTS = Params ;
	ParamsTS.Push(&TERMresc,"TERMS_RESCALING");

	cname = ICM_PRICER_HOMOGENEOUS_SMILE;

	//creéation du modele
	std::vector<const ICM_DefaultCurve*> DefaultCurves(1); 
	DefaultCurves[0] = itsIdxDefaultCurve;

	ICM_ModelMultiCurves mmcTS(DefaultCurves,itsIdxDefaultCurve->GetZeroCurve(),NULL,itsCorrelation);

	//pour chacun des indices
	for (int i=0;i<nbslices;i++)
	{
		basecorrel= (ICM_VolInterpol*) itsCorrelation->GetSlices()[i].GetVolCurve();

		//on récupére les maturités standards
		ARM_Vector YFMaturities = *(*basecorrel).GetExpiryTerms();
		ARM_Vector Maturities = *(*basecorrel).GetExpiryTerms();
		Maturities *= 365.;
		Maturities += AsOf.GetJulian();

		//on récupére les strikes standards
		ARM_Vector Strikes = *(*basecorrel).GetStrikes();
		Strikes /= 100.;

		itsStrikes.clear();
		for (int k=0; k<Strikes.GetSize(); k++)
		{itsStrikes.push_back(Strikes.Elt(k));}

		ARM_Date Start = AsOf.AddDays(1);							//départ à Asof+1j
		ARM_Date Maturity = Maturities.Elt(Maturities.size()-1);	//on choisit la maturité finale

		if ((itsCorreltype==qCAL_BASE_CORRELATION)&&(Matu_)) 
			{Maturity=Matu_;}

		// char* Collateral[TRX_EUR_NBNAMES];
		// for (k=0; k<TRX_EUR_NBNAMES; k++)
		// {Collateral[k]=(char*)itsIdxDefaultCurve->GetLabel().c_str();}
		std::vector<std::string> Collateral(TRX_EUR_NBNAMES); 
		for (k=0; k<TRX_EUR_NBNAMES; k++)
			Collateral[k]=itsIdxDefaultCurve->GetLabel(); 


		double strike_dw=0.;
		double strike_up=0.03;
		double amount = 1.e7;
		double amountfixed = 1.e7;

		ICM_Mez* indexCDO = CdoIndexDefinition(Start,Maturity, Collateral,
										// TRX_EUR_NBNAMES,
										0.,strike_up,
										0.01,amount,INCLUDE_MATURITY,false,
										DEFAULT_CREDIT_LAG_INDX,itsIdxDefaultCurve->GetCurrency() );

		int numflows=indexCDO->GetFeeLeg()->GetNumFlows();

		ICM_Pricer_Advisor Advisor;
		ICM_Pricer_Distrib_Smile* pricer=NULL;

		if (itsCorreltype==qCAL_BASE_CORRELATION_TS)
			pricer = (ICM_Pricer_Distrib_Smile*) Advisor.GeneratePricer(indexCDO,&mmcTS,cname,CREDIT_DEFAULT_VALUE,&ParamsTS,AsOf);
		else if (itsCorreltype==qCAL_BASE_CORRELATION)
			pricer = (ICM_Pricer_Distrib_Smile*) Advisor.GeneratePricer(indexCDO,&mmcTS,cname,CREDIT_DEFAULT_VALUE,&Params,AsOf);

		double price=pricer->Price(qCMPPRICE);

		// itsLossUnit = pricer->GetLossUnit();
		itsLossUnit = pricer->getLossUnits().getLossUnit(AsOf);
		//nombre de losses du portefeuille
		int sizelossunit = TRX_EUR_NBNAMES;

		//definition de la matrice d'expected loss
		double step = ((double)itsStep)/365.;
		int step_days = 30;
		double YFMaturity = (indexCDO->GetFeeLeg()->GetMaturity()-AsOf)/365.;
		int sizesched = (int)(YFMaturity/step)+1;

		const ARM_Vector& Vyfstart = indexCDO->GetFeeLeg()->GetCreditInfos()->GetYFAccStartDates();
		double yfend = indexCDO->GetFeeLeg()->GetCreditInfos()->GetYFAccEndDates().Elt(numflows-1);
		
		
		ARM_Vector* FinalSched = NULL;
		ARM_Vector* FinalSched1 = NULL;
		ARM_Vector SchedSynth(sizesched);
		for (int l=0;l<sizesched;l++) {SchedSynth.Elt(l) = (double)l*step;}
		MergeDates(&FinalSched1,&SchedSynth,&YFMaturities);
		MergeDates(&FinalSched,FinalSched1,&Vyfstart);
		
/*		ARM_Vector* FinalSched = NULL;
		FinalSched = (ARM_Vector*) YFMaturities.Clone();
*/
		itsTimeStep.Resize(1,(FinalSched)->GetSize()+1);

		for (l=0;l<FinalSched->GetSize();l++)
		{itsTimeStep(0,l+1) = (FinalSched)->Elt(l);}
		itsTimeStep(0,0)=0.;

		sizesched=FinalSched->GetSize()+1;

		if (FinalSched) delete FinalSched;
		FinalSched=NULL;

		int nbrows = Strikes.GetSize();

		ICM_QMatrix<double> MatEl(nbrows,sizesched);
		ICM_QMatrix<double> PricesBefore(1.,nbrows);

		if (indexCDO) delete indexCDO;indexCDO=NULL;
		if (pricer) delete pricer;pricer=NULL;

		ICM_Smile_Correlation* correl = NULL;

		for (int j=0;j<nbrows;j++)
		{

			strike_dw = 0.;
			strike_up = Strikes.Elt(j);

			if (itsCorreltype==qCAL_BASE_CORRELATION)
			{if (Matu_) {Maturity=Matu_;}}
			
			ICM_Mez* indexCDO = CdoIndexDefinition(Start,Maturity, Collateral,
										// TRX_EUR_NBNAMES,
										0.,strike_up,
										0.01,amount,INCLUDE_MATURITY,false,
										DEFAULT_CREDIT_LAG_INDX, itsIdxDefaultCurve->GetCurrency());

			
			correl = (ICM_Smile_Correlation*) (mmcTS.GetCorrelation());
			vector<double> strikedw;strikedw.push_back(0.);
			correl->GetSlices()[i].SetSmileStrikeLow(strikedw);
			vector<double> strikeup;strikeup.push_back(strike_up);
			correl->GetSlices()[i].SetSmileStrikeHigh(strikeup);

			if (itsCorreltype==qCAL_BASE_CORRELATION_TS)
				pricer = (ICM_Pricer_Distrib_Smile*) Advisor.GeneratePricer(indexCDO,&mmcTS,cname,CREDIT_DEFAULT_VALUE,&ParamsTS,AsOf);
			else if (itsCorreltype==qCAL_BASE_CORRELATION)
				pricer = (ICM_Pricer_Distrib_Smile*) Advisor.GeneratePricer(indexCDO,&mmcTS,cname,CREDIT_DEFAULT_VALUE,&Params,AsOf);

			double price=pricer->Price(qCMPPRICE);
			PricesBefore(0,j)=price;

			double test = 0.;
			for (int k=0;k<sizesched;k++)
			{
				if (strike_up)
				{MatEl(j,k)= pricer->getDistribLoss().InterpolEL(itsTimeStep(0,k));}
				else
				{MatEl(j,k)= 1.;}
			}

			if (indexCDO) delete indexCDO;indexCDO=NULL;
			if (pricer) delete pricer;pricer=NULL;
		}

		itsELoss.push_back(MatEl);
		itsPricesBefore.push_back(PricesBefore);
	}
}
