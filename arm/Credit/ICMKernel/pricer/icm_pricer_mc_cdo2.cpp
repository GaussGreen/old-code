#include "firsttoinc.h"

#include "ICMKernel\pricer\icm_pricer_mc_cdo2.h"
#include "ICMKernel\crv\icm_defaultcurve.h"
#include "ICMKernel\inst\icm_cdo2.h"
#include "ICMKernel\inst\icm_nthtd.h"
#include "ICMKernel\inst\icm_collateral.h"
#include "ICMKernel\glob\icm_correlation.h"
#include "ICMKernel/util/icm_gengauss1F.h"
#include "ICMKernel/mod/modelmulticurves.h"

#include <numeric>

void ICM_Pricer_MC_Cdo2::Init()
{
	SetName(ICM_MC_PRICER_CDO2);

	itsIndCDO.clear() ;
	itsStrikesDown.clear() ;
	itsStrikesUp.clear() ;
// 17783 	itsBarriers=NULL;
	itsPortfolioLosses.clear() ;
// 17783 	itsDefaultTimes=NULL;
	itsTrancheLosses.clear() ;
	itsZcPay.clear() ;
	itsSpreads.clear() ;
	itsAccrualDates.clear() ;
	itsCouponsDates.clear();
	itsYFPayDates.clear() ;
	itsNbNames    = 0 ;
	itsNbTranches = 0 ;
	itsNbPeriods  = 0 ;
	itsStrikeDown = 0. ;
	itsStrikeUp   = 0. ;
	itsCdo2PortfolioLoss = 0. ;
	itsCdo2TrancheLoss   = 0. ;
	itsField      = NULL ;
	itsCollatLoss = NULL ;
	//itsUnionNames = NULL ;
	itsInstrumentType.clear();
	its_NTD_NumDef.clear();
	its_NTD_ActualDef.clear();
	its_sz_dt_i = 0;
	its_sz_dt_j = 0;
	itsDefaultTimes_Standart.clear();
	itsBarriers_Standart.clear();
// 17783 	itsPerturbDefaultCurves = NULL;
// 17783 	itsCorrespondanceUnionNames_Sensis.clear();
}

ICM_Pricer_MC_Cdo2::~ICM_Pricer_MC_Cdo2()
{
	int i=0,j=0;

	//FreePointerTabChar(itsUnionNames,itsNbNames);	

	if(itsField)
		delete itsField ;
	itsField = NULL;

	if(itsCollatLoss)
		delete itsCollatLoss ;
	itsCollatLoss = NULL ;

/** // 17783 
	if (GetSensiManager())
	if (GetSensiManager()->GetSensiMatrix(ICM_FAST_SPREAD_TYPE))
	{
	its_sz_dt_i = (GetSensiManager()->GetSensiMatrix(ICM_FAST_SPREAD_TYPE))->GetMktDataSize();
	its_sz_dt_j = (GetSensiManager()->GetSensiMatrix(ICM_FAST_SPREAD_TYPE))->GetBumpProfileSize();

	if (itsDefaultTimes)
	{
		for (i =0; i<its_sz_dt_i; i++)
			for (j =0; j<its_sz_dt_j; j++)
			{
				if (itsDefaultTimes->Getvalue(i,j))
				{
					(itsDefaultTimes->Getvalue(i,j))->clear();
					delete (itsDefaultTimes->Getvalue(i,j));
					itsDefaultTimes->SetValue(i,j,NULL);
				}
			}

		delete itsDefaultTimes;
		itsDefaultTimes = NULL;
	}

	if (itsBarriers)
	{
		for (i =0; i<its_sz_dt_i; i++)
			for (j =0; j<its_sz_dt_j; j++)
			{
				if (itsBarriers->Getvalue(i,j))
				{
				(itsBarriers->Getvalue(i,j))->clear();
				delete (itsBarriers->Getvalue(i,j));
				itsBarriers->SetValue(i,j,NULL);
				}
			}

		delete itsBarriers;
		itsBarriers = NULL;
	}

	}
**/ 
/**
17783 
	if (itsPerturbDefaultCurves)
	{
		for (i =0; i<its_sz_dt_i; i++)
			for (j =0; j<its_sz_dt_j; j++)
			{
				if (itsPerturbDefaultCurves->Getvalue(i,j))
				{
				delete (itsPerturbDefaultCurves->Getvalue(i,j));
				itsPerturbDefaultCurves->SetValue(i,j,NULL);
				}
			}

		delete itsPerturbDefaultCurves;
		itsPerturbDefaultCurves = NULL;
	}
**/ 
	itsDefaultTimes_Standart.clear();
	itsBarriers_Standart.clear();
}


void ICM_Pricer_MC_Cdo2::ResetPricer(void)  
{
	int i =0,j=0;

	ICM_Pricer_MC::ResetPricer();

	itsIndCDO.clear();
	itsZcPay.clear();

	if(itsField)
		delete itsField ;
	itsField = NULL;

	if(itsCollatLoss)
		delete itsCollatLoss ;
	itsCollatLoss = NULL ;

	/**if (itsUnionNames)
	{
		for (int i =0; i<itsNbNames; i++)
		{
			if (itsUnionNames[i])
				delete[] itsUnionNames[i];
			itsUnionNames[i] = NULL;
		}
	}**/ 

	itsSpreads.clear() ;
	itsAccrualDates.clear() ;
	itsCouponsDates.clear() ;

	itsYFPayDates.clear() ;
	itsStrikesDown.clear() ;
	itsStrikesUp.clear() ;
	itsPortfolioLosses.clear() ;
	itsTrancheLosses.clear() ;

	itsNbTranches = 0;
	itsNbPeriods = 0;
	itsStrikeDown = 0. ;
	itsStrikeUp   = 0. ;
	itsCdo2PortfolioLoss = 0. ;
	itsCdo2TrancheLoss   = 0. ;
	itsInstrumentType.clear();
	its_NTD_NumDef.clear();
	its_NTD_ActualDef.clear();

	itsDefaultTimes_Standart.clear();
	itsBarriers_Standart.clear();

	its_sz_dt_i = 0;
	its_sz_dt_j = 0;
}

/**

  17783 	
  void ICM_Pricer_MC_Cdo2::ResetSensiManager()
{
	int i=0,j=0;

	if (GetSensiManager())
	if (GetSensiManager()->GetSensiMatrix(ICM_FAST_SPREAD_TYPE))
	{
	its_sz_dt_i = (GetSensiManager()->GetSensiMatrix(ICM_FAST_SPREAD_TYPE))->GetMktDataSize();
	its_sz_dt_j = (GetSensiManager()->GetSensiMatrix(ICM_FAST_SPREAD_TYPE))->GetBumpProfileSize();

	if (itsDefaultTimes)
	{
		for (i =0; i<its_sz_dt_i; i++)
			for (j =0; j<its_sz_dt_j; j++)
			{
				if (itsDefaultTimes->Getvalue(i,j))
				{
					(itsDefaultTimes->Getvalue(i,j))->clear();
					delete (itsDefaultTimes->Getvalue(i,j));
					itsDefaultTimes->SetValue(i,j,NULL);
				}
			}

		delete itsDefaultTimes;
		itsDefaultTimes = NULL;
	}

	itsDefaultTimes = new ICM_QMatrix<set<Y>*>(its_sz_dt_i,its_sz_dt_j,NULL);
		
	
	if (itsBarriers)
	{
		for (i =0; i<its_sz_dt_i; i++)
			for (j =0; j<its_sz_dt_j; j++)
			{
				if (itsBarriers->Getvalue(i,j))
				{
				(itsBarriers->Getvalue(i,j))->clear();
				delete (itsBarriers->Getvalue(i,j));
				itsBarriers->SetValue(i,j,NULL);
				}
			}

		delete itsBarriers;
		itsBarriers = NULL;
	}

	itsBarriers = new ICM_QMatrix<vector<double>*>(its_sz_dt_i,its_sz_dt_j,NULL);

	}
}
17783 	
**/ 

void ICM_Pricer_MC_Cdo2::CptIntermediate(qSENSITIVITY_TYPE typesensi)	
{
	int i=0,j=0;

	ICM_Cdo2   * cdo2 = (ICM_Cdo2*) GetSecurity() ;
	ICM_ModelMultiCurves* mod = (ICM_ModelMultiCurves*) GetModel();

	//On set le nom des instrument pour gérer le cas d'un CDO² de Ntd
	itsInstrumentType.clear();
	for (i=0; i< cdo2->GetPortfolio()->GetNbSec(); i++)
		itsInstrumentType.push_back(cdo2->GetPortfolio()->GetSecurities()[i]->GetName());

	ComputeUnionNames();

	itsNbNames	  = itsUnionNames.size();
	itsNbTranches = cdo2->GetPortfolio()->GetNbSec() ;
	itsNbPeriods  = cdo2->GetFeeLeg()->GetPaymentDates()->size() ;
	itsStrikeDown = cdo2->GetSubAmount(cdo2->GetFeeLeg()->GetStartDate()) ;
	itsStrikeUp   = itsStrikeDown + cdo2->GetMezzAmount(cdo2->GetFeeLeg()->GetStartDate()) ;

	vector<double> betas ;
	betas.clear();

	ARM_Vector VBetas = mod->GetCorrelation()->GetBetaVector(itsUnionNames,itsNbNames);

	for (i=0; i<itsNbNames; i++) betas.push_back(VBetas.Elt(i));

	SetGenerator((ICM_Root_Generator *) new ICM_Gengauss1F(0,itsNbNames,betas));

// 17783 		if (typesensi == ICM_FAST_SPREAD_TYPE) ResetSensiManager();

	SetTrancheLosses(0.) ;
	SetPortfolioLosses(0.) ;
	SetCdo2PortfolioLoss(0.) ;
	SetCdo2TrancheLoss(0.) ;
	ComputeFieldAndCollatLoss() ;
	ComputeBarriers_Complete() ;
	ComputeIndCDO() ;
	ComputeZc() ;
	ComputeSpreads() ;
	ComputeAccrualDates() ;
	ComputeYFPayDates() ;
}
void ICM_Pricer_MC_Cdo2::ComputeFieldAndCollatLoss()
{
	ICM_Cdo2 * cdo2 = (ICM_Cdo2*) GetSecurity() ;

	if (itsField)
		delete itsField;
	itsField = new ICM_QMatrix<int>(itsNbNames, itsNbTranches) ;

	if (itsCollatLoss)
		delete itsCollatLoss;
	itsCollatLoss = new ICM_QMatrix<double>(itsNbNames, itsNbTranches) ;

	double binary = cdo2->GetBinary(); 
	if (binary == CREDIT_DEFAULT_VALUE) binary = 0.;

	ICM_ModelMultiCurves * mod = (ICM_ModelMultiCurves*) GetModel() ;
	
	itsStrikesDown.resize(itsNbTranches) ;
	itsStrikesUp.resize(itsNbTranches) ;

	its_NTD_NumDef.resize(itsNbTranches) ;
	its_NTD_ActualDef.resize(itsNbTranches) ;

	int nbissuers_sec = 0;
	const std::vector<std::string>& labels = cdo2->GetCollateral()->GetIssuersLabels() ;
	
	for (int i = 0; i<itsNbNames ; i++)
	{
		for(int j = 0; j< itsNbTranches; j++)
		{
			ARM_Security * Sec     = cdo2->GetPortfolio()->GetSecurity(j) ;

			ICM_Mez      * Sec_Aux = (ICM_Mez*) Sec;
			
			const std::vector<std::string>& labels_Sec      = Sec_Aux->GetCollateral()->GetIssuersLabels();
			nbissuers_sec = Sec_Aux->GetCollateral()->GetNbIssuers() ;

			ARM_Vector * Not  = new ARM_Vector(nbissuers_sec,0.) ;
			ARM_Vector * Rec = new ARM_Vector(nbissuers_sec,0.) ;

			//fabs : afin de prendre en compte les tranches filles short dans le CDO2 (nominal < 0 -> tradedcoef < 0)
			for (int il=0; il<nbissuers_sec; il++)
			{
				// Not->InitElt(il,Sec_Aux->GetCollateral()->GetIssuersNotional(labels_Sec[il])*fabs(Sec_Aux->GetTradedCoef())) ;
				Not->InitElt(il,Sec_Aux->GetCollateral()->GetIssuersNotional(il)*fabs(Sec_Aux->GetTradedCoef())) ;

			if (binary)
				Rec->InitElt(il,binary); 
			else
				Rec->InitElt(il,mod->GetRecoveryRate(labels_Sec[il])); 
			}
			
			for(int k = 0; k< nbissuers_sec; k++)
			{
				if(!(strcmp(labels[i].c_str(),labels_Sec[k].c_str())))
				{	
					itsField->SetValue(i,j,1) ;
					itsCollatLoss->SetValue(i,j,Not->Elt(k)*(1.- Rec->Elt(k))) ;
				}
			}

			if (itsInstrumentType[j] == ICM_MEZ)
			{	
			itsStrikesDown[j] = Sec_Aux->GetSubAmount(Sec_Aux->GetFeeLeg()->GetStartDate())*fabs(Sec_Aux->GetTradedCoef()) ;
			itsStrikesUp[j]   = itsStrikesDown[j] + Sec_Aux->GetMezzAmount(Sec_Aux->GetFeeLeg()->GetStartDate()) *fabs(Sec_Aux->GetTradedCoef());
			}
			else //cas Nth to default
			{	
			its_NTD_NumDef[j] = (double) ((ICM_Nthtd*) Sec)->GetFirstNumDefault() ;
			its_NTD_ActualDef[j]   = 0.;
			}

			delete Not ; delete Rec ;
		}
	}
	
}


// Retourne le nombre de payment dates avant le moment du defaut
int ICM_Pricer_MC_Cdo2::Eta(const double& tau)
{
	int cmpt = 0 ;
	int i = 0 ;
	int size = itsYFPayDates.size(); 

	if (tau >= itsYFPayDates[size-1])
		return size;

	while(tau > itsYFPayDates[i])
		{cmpt++ ; i++ ;}
	return cmpt ;
}



void ICM_Pricer_MC_Cdo2::ComputeIndCDO()
{
	itsIndCDO.resize(itsNbNames) ;
	for(unsigned int j = 0; j < itsNbNames ; j++)
	{
		for(unsigned int k = 0; k < itsNbTranches ; k++)
		{
			if( itsField->Getvalue(j,k) == 1 )
				itsIndCDO[j].push_back(k) ;
		}
	}

}

// ----------------------------------------------------------------------------
// Compute Barriers in complete case
// ----------------------------------------------------------------------------
void ICM_Pricer_MC_Cdo2::ComputeBarriers_Complete()
{
	ICM_Cdo2* cdo2 = (ICM_Cdo2*) GetSecurity() ;
	ICM_ModelMultiCurves * model = (ICM_ModelMultiCurves*)GetModel() ;
	ICM_Security* security = cdo2->GetFeeLeg()->GetCreditInfos();

	vector<double>* Barriers = NULL;
	int i = 0,j=0,k=0;

	double enddate = (security->GetEndDateNA().GetJulian() - security->GetStartDateNA().GetJulian())/365.;
	double PD = 0.;


/**
17783 	
	if (GetSensiManager())
	if (GetSensiManager()->GetSensiMatrix(ICM_FAST_SPREAD_TYPE))
	{

	its_sz_dt_i = (GetSensiManager()->GetSensiMatrix(ICM_FAST_SPREAD_TYPE))->GetMktDataSize();
	its_sz_dt_j = (GetSensiManager()->GetSensiMatrix(ICM_FAST_SPREAD_TYPE))->GetBumpProfileSize();

	CptCorrespondanceUnionNames_Sensis();

	if (itsBarriers)
	{
		for (i =0; i<its_sz_dt_i; i++)
			for (j =0; j<its_sz_dt_j; j++)
			{
			   if (itsBarriers->Getvalue(i,j))
			   {		 
			   (itsBarriers->Getvalue(i,j))->clear();
			   delete (itsBarriers->Getvalue(i,j));
			   itsBarriers->SetValue(i,j,NULL);
			   }
			}
		delete itsBarriers;
		itsBarriers = NULL;
	}
	itsBarriers = new ICM_QMatrix<vector<double>*>(its_sz_dt_i,its_sz_dt_j,NULL);
	}
17783 	**/ 

	itsBarriers_Standart.resize(itsNbNames) ;

	for(i = 0; i < itsNbNames ; i++)
	{
		PD = model->GetDefaultCurve(itsUnionNames[i])->DefaultProba(enddate);
		if (PD>1.)
			ICMTHROW(ERR_INVALID_MODEL,"Error: Default Probability >1. for "<<itsUnionNames[i]); 
		if (PD<0.0)
			ICMTHROW(ERR_INVALID_MODEL,"Error: Default Probability <0. for "<<itsUnionNames[i]); 
		if (fabs(PD - 1.)<DB_TOL) PD = 0.999999999;

		(itsBarriers_Standart)[i] = NAG_deviates_normal( PD  ) ;
	}
/**
17783 	
	ICM_Sensi2D* sensi2D = GetSensiManager()->GetSensiMatrix(ICM_FAST_SPREAD_TYPE);

	for(i = 0; i < its_sz_dt_i ; i++)
		for(j = 0; j < its_sz_dt_j ; j++)
		{
		Barriers = new vector<double>;
		Barriers->resize(itsNbNames) ;

		for(k = 0; k<itsNbNames;k++) (*Barriers)[k] = itsBarriers_Standart[k];

		PD = ((ICM_QMatrix<ICM_DefaultCurve*>*)sensi2D->GetObjectMatrix())->Getvalue(i,j)->DefaultProba(enddate);

		for(k = 0; k<itsNbNames;k++) { if (strcmp(itsUnionNames[k].c_str(),sensi2D->GetMktData(i).c_str()) == NULL) break;}

		(*Barriers)[k] = NAG_deviates_normal( PD );

		itsBarriers->SetValue(i,j,Barriers);
		}
17783 	**/ 
		
}


// ----------------------------------------------------------------------------
// Generation of default times with correlation matrix
// ----------------------------------------------------------------------------
void ICM_Pricer_MC_Cdo2::GenerateDefaultTimes_troughCorrelMatrix_Standart()
{
	itsDefaultTimes_Standart.clear() ;
	ICM_Root_Generator * Gen = GetGenerator() ;
	ICM_ModelMultiCurves * model = (ICM_ModelMultiCurves*)GetModel() ;

	std::vector<double> its_indep_normal_vector; 
	its_indep_normal_vector.resize(itsNbNames);
	std::vector<double> its_corr_normal_vector; 
	its_corr_normal_vector.resize(itsNbNames);

	double aux = 0. ;
	// génération de normales indépendantes

	for(unsigned int i = 0; i < itsNbNames; i++)
	{	
		its_corr_normal_vector[i]=0.;
		its_indep_normal_vector[i]=NAG_random_normal(0.,1.);
	}

	for( i = 0; i < itsNbNames; i++)
	{	
		for(unsigned int j = 0; j < itsNbNames; j++)
		{
			its_corr_normal_vector[i]+=its_indep_normal_vector[j]*itsCholeskyMatrix.Getvalue(i,j);
		}
	}

	for(  i = 0; i < itsNbNames; i++)
	{	
		if(its_corr_normal_vector[i] < itsBarriers_Standart[i])
		{
			aux = NAG_cumul_normal(its_corr_normal_vector[i]);
			double  t = model->GetDefaultCurve(itsUnionNames[i])->DefProbInverse(aux);
			itsDefaultTimes_Standart.insert(Y(i,t)) ;
		}
	}
}


// ----------------------------------------------------------------------------
// Compute Price 
// ----------------------------------------------------------------------------
double ICM_Pricer_MC_Cdo2::ComputePrice(qCMPMETH measure)
{
	double price = 0.;

	if (GetPriceFlg()) return CheckForPrice(measure);

	CptIntermediate(ICMSPREAD_TYPE); 

	ICM_Cdo2* Cdo2= (ICM_Cdo2*)(GetSecurity());

	if (Cdo2->IsCrossSubordinate())
		price = ComputePrice_Basket_CrossSub(measure);
	else
		price = ComputePrice_Basket(measure);

	return (price);
}

// ----------------------------------------------------------------------------
// Compute Price for single price
// ----------------------------------------------------------------------------
double ICM_Pricer_MC_Cdo2::ComputePrice_Basket(qCMPMETH measure)
{

	ICM_Cdo2 * cdo2 = (ICM_Cdo2*) GetSecurity() ;
	ICM_ModelMultiCurves * model = (ICM_ModelMultiCurves*) GetModel() ;
	ARM_Date AsOf = model->GetStartDate() ;
	ICM_Security* security = cdo2->GetFeeLeg()->GetCreditInfos();
	security->ComputeYF(AsOf);

	const ARM_Vector & AccStartDates   = security->GetAccStartDates();
	const ARM_Vector & AccEndDates     = security->GetAccEndDates();

	ICM_Correlation* Correlation = model->GetCorrelation();
	ARM_Vector Betas = Correlation->GetBetaVector(itsUnionNames,itsNbNames);

	int nPaths = GetNbPaths() ;

    int cmpt = 0 ;
	while ((AsOf.GetJulian()>AccStartDates.Elt(cmpt)) &&  (AsOf.GetJulian()>=AccEndDates.Elt(cmpt)))
			cmpt++ ;

	double Tranche  = itsStrikeUp - itsStrikeDown ;
	double Notional = Tranche ;
	double NotionalPrev = Tranche ;

	double premiumLeg = 0. ;
	double premiumLeg_Unity = 0. ;
	double period_val = 0.;

	double defaultLeg = 0. ;
	double price =0.;
	int k1=0,k2=0;

     
    {

	if(	Betas.empty())
	{
		itsCholeskyMatrix.Resize(itsNbNames,itsNbNames);
		itsCholeskyMatrix=Correlation->ComputeCholeskyMatrix();
	}

	for(unsigned int i = 0; i < nPaths; i++)
	{
		if(Betas.empty())
			GenerateDefaultTimes_troughCorrelMatrix_Standart() ;
		else
			GenerateDefaultTimes_Betas_Standart();

		double t ;
		int iName ;
		int deb = cmpt ;

		//Reinit des défauts actuels
		for (iName=0;iName<its_NTD_ActualDef.size();iName++ ) its_NTD_ActualDef[iName] = 0;

		for(set<Y>::iterator iterateur = itsDefaultTimes_Standart.begin(); iterateur != itsDefaultTimes_Standart.end(); ++ iterateur)
		{
			double em = 0 ;
			t = iterateur->tau ;
			iName = iterateur->id ;
			int e = Eta(t) ;
			for(unsigned int j = deb; j < e ; j++)
			{ 
				// Computes premimum leg
				switch (cdo2->GetFeeLeg()->GetAccruedOnDefault())
				{
				case qCONTINUE_TO_MATURITY:
					premiumLeg +=  itsZcPay[j] * Tranche * (itsCouponsDates[j] + itsAccrualDates[j]) * itsSpreads[j] ;
					premiumLeg_Unity +=  itsZcPay[j] * Tranche * (itsCouponsDates[j] + itsAccrualDates[j]);
					break;
				case qACCRUED_NOT_SETTLED :	
					premiumLeg +=  itsZcPay[j] * Notional * (itsCouponsDates[j] + itsAccrualDates[j]) * itsSpreads[j] ;
					premiumLeg_Unity +=  itsZcPay[j] * Notional * (itsCouponsDates[j] + itsAccrualDates[j]);
					break;
				case qCOMPLETE_CURRENT_PERIOD :
				case qACCRUED_SETTLED : 
				default : 
					premiumLeg +=  itsZcPay[j] * (Notional *itsCouponsDates[j] + Tranche * itsAccrualDates[j]) * itsSpreads[j] ;
					premiumLeg_Unity +=  itsZcPay[j] *(Notional *itsCouponsDates[j] + Tranche*itsAccrualDates[j]);
				}
			}

			for(unsigned int k = 0 ; k < itsIndCDO[iName].size() ; k++)
			{
				int ind = itsIndCDO[iName][k] ;

				//Récupération du traded coef : prise en compte de tranches short
				ARM_Security * Sec     = cdo2->GetPortfolio()->GetSecurity(ind) ;
				ICM_Mez      * Sec_Aux = (ICM_Mez*) Sec;
				double tranche_sign = 0.;
				if (Sec_Aux->GetTradedCoef() >= 0.)
					tranche_sign = 1.;
				else
					tranche_sign = -1.;

				//PayOff
				if (itsInstrumentType[ind] == ICM_MEZ)
				{
				itsPortfolioLosses[ind]  += itsCollatLoss->Getvalue(iName,ind) ;
				itsTrancheLosses[ind]    = tranche_sign * MAX( MIN(itsPortfolioLosses[ind],itsStrikesUp[ind]) - 
										 itsStrikesDown[ind], 0. ) ;
				}
				else
				{
				(its_NTD_ActualDef[ind])++;

				if (its_NTD_ActualDef[ind] == its_NTD_NumDef[ind])
					itsTrancheLosses[ind] = tranche_sign * itsCollatLoss->Getvalue(iName,ind);
				}
			}
			itsCdo2PortfolioLoss = accumulate(itsTrancheLosses.begin(),itsTrancheLosses.end(),0) ;
			itsCdo2TrancheLoss   = MAX( MIN(itsCdo2PortfolioLoss,itsStrikeUp) - itsStrikeDown, 0. ) ;
			if(itsCdo2TrancheLoss>0.1)
			{
				NotionalPrev = Notional ;
				Notional = Tranche - itsCdo2TrancheLoss ;
				defaultLeg			 +=  itsZcPay[e] * (NotionalPrev - Notional) ;

				if(e>0)          em	 = itsYFPayDates[e-1] ;

				//OnDefault Accrued case (end date notional is calculated on the next previous loop)
				switch (cdo2->GetFeeLeg()->GetAccruedOnDefault())
				{
				case qCONTINUE_TO_MATURITY:
				case qACCRUED_NOT_SETTLED :	
					break;
				case qCOMPLETE_CURRENT_PERIOD :
					premiumLeg			 +=  itsZcPay[e] * (NotionalPrev - Notional) * (itsYFPayDates[e]-em) * itsSpreads[e] ; 
					premiumLeg_Unity     +=  itsZcPay[e] * (NotionalPrev - Notional) * (itsYFPayDates[e]-em) ; 	
					break;
				case qACCRUED_SETTLED : 
				default : 
					premiumLeg			 +=  itsZcPay[e] * (NotionalPrev - Notional) * (t-em) * itsSpreads[e] ; 
					premiumLeg_Unity     +=  itsZcPay[e] * (NotionalPrev - Notional) * (t-em) ; 	
				}
			}
			deb = e ;
		}
		for(unsigned int l = deb; l < itsNbPeriods ; l++)
		{
			// Computes premimum leg
			switch (cdo2->GetFeeLeg()->GetAccruedOnDefault())
			{
			case qCONTINUE_TO_MATURITY:
				premiumLeg +=  itsZcPay[l] * Tranche * (itsCouponsDates[l] + itsAccrualDates[l]) * itsSpreads[l] ;
				premiumLeg_Unity +=  itsZcPay[l] * Tranche * (itsCouponsDates[l] + itsAccrualDates[l]);
				break;
			case qACCRUED_NOT_SETTLED :	
				premiumLeg +=  itsZcPay[l] * Notional * (itsCouponsDates[l] + itsAccrualDates[l]) * itsSpreads[l] ;
				premiumLeg_Unity +=  itsZcPay[l] * Notional * (itsCouponsDates[l] + itsAccrualDates[l]);
				break;
			case qCOMPLETE_CURRENT_PERIOD :
			case qACCRUED_SETTLED : 
			default : 
				premiumLeg +=  itsZcPay[l] * (Notional *itsCouponsDates[l] + Tranche * itsAccrualDates[l]) * itsSpreads[l] ;
				premiumLeg_Unity +=  itsZcPay[l] *(Notional *itsCouponsDates[l] + Tranche*itsAccrualDates[l]);
			}
		}
	
		SetPortfolioLosses(0.) ;	
		SetTrancheLosses(0.) ;
		itsCdo2PortfolioLoss = 0. ;
		itsCdo2TrancheLoss = 0. ;
		Notional = Tranche ;
		NotionalPrev = Tranche ;
	}

	defaultLeg			 =  defaultLeg / nPaths ;
	premiumLeg			 =  premiumLeg / nPaths ;
	premiumLeg_Unity     =  premiumLeg_Unity / nPaths ;
 	
	price = premiumLeg - defaultLeg ;
	
	double Flag = 1.;
	// 14514 if (cdo2->GetInitialNotional()<0) Flag = -1.;
	if (cdo2->GetFeeLeg()->GetCreditInfosRef().GetNotionals().Elt(0)<0) Flag = -1.;

	SetPrice(Flag*price*cdo2->GetTradedCoef());
	SetDefLegPrice(Flag*defaultLeg*cdo2->GetTradedCoef());
	SetFeeLegPrice(Flag*premiumLeg*cdo2->GetTradedCoef());
	SetFeeLegPrice_Unity(Flag*premiumLeg_Unity*cdo2->GetTradedCoef());

	price = CheckForPrice(measure);

    }	
    

	return (price) ;
}


// ----------------------------------------------------------------------------
// Compute Price for single price with cross subordination
// ----------------------------------------------------------------------------
double ICM_Pricer_MC_Cdo2::ComputePrice_Basket_CrossSub(qCMPMETH measure)
{

	ICM_Cdo2 * cdo2 = (ICM_Cdo2*) GetSecurity() ;
	ICM_ModelMultiCurves * model = (ICM_ModelMultiCurves*) GetModel() ;
	ARM_Date AsOf = model->GetStartDate() ;
	ICM_Security* security = cdo2->GetFeeLeg()->GetCreditInfos();
	security->ComputeYF(AsOf);

	const ARM_Vector & AccStartDates   = security->GetAccStartDates();
	const ARM_Vector & AccEndDates     = security->GetAccEndDates();

	ICM_Correlation* Correlation = model->GetCorrelation();
	ARM_Vector  Betas = Correlation->GetBetaVector(itsUnionNames,itsNbNames);

	int nPaths = GetNbPaths() ;

    int cmpt = 0 ;
	while ((AsOf.GetJulian()>AccStartDates.Elt(cmpt)) &&  (AsOf.GetJulian()>=AccEndDates.Elt(cmpt)))
			cmpt++ ;

	double Tranche  = itsStrikeUp - itsStrikeDown ;
	double Notional = Tranche ;
	double NotionalPrev = Tranche ;

	double premiumLeg = 0. ;
	double premiumLeg_Unity = 0. ;
	double period_val = 0.;

	double defaultLeg = 0. ;
	double price =0.;
	int k1=0,k2=0;

     
    {

	if(	Betas.empty())
	{
		itsCholeskyMatrix.Resize(itsNbNames,itsNbNames);
		itsCholeskyMatrix=Correlation->ComputeCholeskyMatrix();
	}

	for(unsigned int i = 0; i < nPaths; i++)
	{
		if(Betas.empty())
			GenerateDefaultTimes_troughCorrelMatrix_Standart() ;
		else
			GenerateDefaultTimes_Betas_Standart();

		double t ;
		int iName ;
		int deb = cmpt ;

		//Reinit des défauts actuels
		for (iName=0;iName<its_NTD_ActualDef.size();iName++ ) its_NTD_ActualDef[iName] = 0;

		for(set<Y>::iterator iterateur = itsDefaultTimes_Standart.begin(); iterateur != itsDefaultTimes_Standart.end(); ++ iterateur)
		{
			double em = 0 ;
			t = iterateur->tau ;
			iName = iterateur->id ;
			int e = Eta(t) ;
			for(unsigned int j = deb; j < e ; j++)
			{ 
				// Computes premimum leg
				switch (cdo2->GetFeeLeg()->GetAccruedOnDefault())
				{
				case qCONTINUE_TO_MATURITY:
					premiumLeg +=  itsZcPay[j] * Tranche * (itsCouponsDates[j] + itsAccrualDates[j]) * itsSpreads[j] ;
					premiumLeg_Unity +=  itsZcPay[j] * Tranche * (itsCouponsDates[j] + itsAccrualDates[j]);
					break;
				case qACCRUED_NOT_SETTLED :	
					premiumLeg +=  itsZcPay[j] * Notional * (itsCouponsDates[j] + itsAccrualDates[j]) * itsSpreads[j] ;
					premiumLeg_Unity +=  itsZcPay[j] * Notional * (itsCouponsDates[j] + itsAccrualDates[j]);
					break;
				case qCOMPLETE_CURRENT_PERIOD :
				case qACCRUED_SETTLED : 
				default : 
					premiumLeg +=  itsZcPay[j] * (Notional *itsCouponsDates[j] + Tranche * itsAccrualDates[j]) * itsSpreads[j] ;
					premiumLeg_Unity +=  itsZcPay[j] *(Notional *itsCouponsDates[j] + Tranche*itsAccrualDates[j]);
				}
			}

			for(unsigned int k = 0 ; k < itsIndCDO[iName].size() ; k++)
			{
				int ind = itsIndCDO[iName][k] ;

				if (itsInstrumentType[ind] == ICM_MEZ)
				{
				AffectLossForCrossSub(iName,ind,itsCollatLoss,itsStrikesDown,itsStrikesUp,itsPortfolioLosses,ind);
				//itsPortfolioLosses[ind]  += itsCollatLoss->Getvalue(iName,ind) ;
				itsTrancheLosses[ind]    = MAX( MIN(itsPortfolioLosses[ind],itsStrikesUp[ind]) - 
										 itsStrikesDown[ind], 0. ) ;
				}
				else
				{
				(its_NTD_ActualDef[ind])++;

				if (its_NTD_ActualDef[ind] == its_NTD_NumDef[ind])
					AffectLossForCrossSub(iName,ind,itsCollatLoss,itsStrikesDown,itsStrikesUp,itsPortfolioLosses,ind);
					//itsTrancheLosses[ind] = itsCollatLoss->Getvalue(iName,ind);
				}
			}
			itsCdo2PortfolioLoss = accumulate(itsTrancheLosses.begin(),itsTrancheLosses.end(),0) ;
			itsCdo2TrancheLoss   = MAX( MIN(itsCdo2PortfolioLoss,itsStrikeUp) - itsStrikeDown, 0. ) ;
			if(itsCdo2TrancheLoss>0.1)
			{
				NotionalPrev = Notional ;
				Notional = Tranche - itsCdo2TrancheLoss ;
				defaultLeg			 +=  itsZcPay[e] * (NotionalPrev - Notional) ;

				if(e>0)          em	 = itsYFPayDates[e-1] ;

				//OnDefault Accrued case (end date notional is calculated on the next previous loop)
				switch (cdo2->GetFeeLeg()->GetAccruedOnDefault())
				{
				case qCONTINUE_TO_MATURITY:
				case qACCRUED_NOT_SETTLED :	
					break;
				case qCOMPLETE_CURRENT_PERIOD :
					premiumLeg			 +=  itsZcPay[e] * (NotionalPrev - Notional) * (itsYFPayDates[e]-em) * itsSpreads[e] ; 
					premiumLeg_Unity     +=  itsZcPay[e] * (NotionalPrev - Notional) * (itsYFPayDates[e]-em) ; 	
					break;
				case qACCRUED_SETTLED : 
				default : 
					premiumLeg			 +=  itsZcPay[e] * (NotionalPrev - Notional) * (t-em) * itsSpreads[e] ; 
					premiumLeg_Unity     +=  itsZcPay[e] * (NotionalPrev - Notional) * (t-em) ; 	
				}
			}
			deb = e ;
		}
		for(unsigned int l = deb; l < itsNbPeriods ; l++)
		{
			// Computes premimum leg
			switch (cdo2->GetFeeLeg()->GetAccruedOnDefault())
			{
			case qCONTINUE_TO_MATURITY:
				premiumLeg +=  itsZcPay[l] * Tranche * (itsCouponsDates[l] + itsAccrualDates[l]) * itsSpreads[l] ;
				premiumLeg_Unity +=  itsZcPay[l] * Tranche * (itsCouponsDates[l] + itsAccrualDates[l]);
				break;
			case qACCRUED_NOT_SETTLED :	
				premiumLeg +=  itsZcPay[l] * Notional *(itsCouponsDates[l] + itsAccrualDates[l]) * itsSpreads[l] ;
				premiumLeg_Unity +=  itsZcPay[l] * Notional * (itsCouponsDates[l] + itsAccrualDates[l]);
				break;
			case qCOMPLETE_CURRENT_PERIOD :
			case qACCRUED_SETTLED : 
			default : 
				premiumLeg +=  itsZcPay[l] * (Notional *itsCouponsDates[l] + Tranche * itsAccrualDates[l]) * itsSpreads[l] ;
				premiumLeg_Unity +=  itsZcPay[l] *(Notional *itsCouponsDates[l] + Tranche*itsAccrualDates[l]);
			}
		}
	
		SetPortfolioLosses(0.) ;	
		SetTrancheLosses(0.) ;
		itsCdo2PortfolioLoss = 0. ;
		itsCdo2TrancheLoss = 0. ;
		Notional = Tranche ;
		NotionalPrev = Tranche ;
	}

	defaultLeg			 =  defaultLeg / nPaths ;
	premiumLeg			 =  premiumLeg / nPaths ;
	premiumLeg_Unity     =  premiumLeg_Unity / nPaths ;
 	
	price = premiumLeg - defaultLeg ;
	
	double Flag = 1.;
	if (cdo2->GetFeeLeg()->GetCreditInfosRef().GetNotionals().Elt(0)<0) Flag = -1.;
	// 14514 if (cdo2->GetInitialNotional()<0) Flag = -1.;

	SetPrice(Flag*price*cdo2->GetTradedCoef());
	SetDefLegPrice(Flag*defaultLeg*cdo2->GetTradedCoef());
	SetFeeLegPrice(Flag*premiumLeg*cdo2->GetTradedCoef());
	SetFeeLegPrice_Unity(Flag*premiumLeg_Unity*cdo2->GetTradedCoef());

	price = CheckForPrice(measure );

    }	
     

	return (price) ;
}


// ----------------------------------------------------------------------------
//Compute Price + Hedges for Spreads
// ----------------------------------------------------------------------------
/**JLA 
double ICM_Pricer_MC_Cdo2::ComputePrice_Basket_Complete(const double& initialprice)
{
	CptIntermediate(ICM_FAST_SPREAD_TYPE); 

	ICM_Cdo2 * cdo2 = (ICM_Cdo2*) GetSecurity() ;
	ICM_ModelMultiCurves * model = (ICM_ModelMultiCurves*) GetModel() ;
	ARM_Date AsOf = model->GetStartDate() ;
	ICM_Security* security = cdo2->GetFeeLeg()->GetCreditInfos();
	security->ComputeYF(AsOf);

	const ARM_Vector &AccStartDates   = security->GetAccStartDates();
	const ARM_Vector & AccEndDates     = security->GetAccEndDates();

	ICM_Correlation* Correlation = model->GetCorrelation();
	ARM_Vector  Betas = Correlation->GetBetaVector(itsUnionNames,itsNbNames);

	int nPaths = GetNbPaths() ;

    int cmpt = 0 ;
	while ((AsOf.GetJulian()>AccStartDates.Elt(cmpt)) &&  (AsOf.GetJulian()>=AccEndDates.Elt(cmpt)))
			cmpt++ ;

	double Tranche  = itsStrikeUp - itsStrikeDown ;
	double Notional = Tranche ;
	double NotionalPrev = Tranche ;
	double premiumLeg = 0. ;
	double defaultLeg = 0. ;
	double price =0.;
	int k1=0,k2=0;

	ICM_QCubix<double>* Cube = new ICM_QCubix<double>(its_sz_dt_i,its_sz_dt_j,nPaths,0.);

     
    {

	for(unsigned int i = 0; i < nPaths; i++)
	{
		GenerateDefaultTimes_Betas_Complete();

		for(k1 = 0; k1 < its_sz_dt_i ; k1++)
		for(k2 = 0; k2 < its_sz_dt_j ; k2++)
		{

		premiumLeg = 0. ;
		defaultLeg = 0. ;

		double t = 0. ;
		int iName = 0 ;
		int deb = cmpt ;

		//Reinit des défauts actuels
		for (iName=0;iName<its_NTD_ActualDef.size();iName++ ) its_NTD_ActualDef[iName] = 0;

		for(set<Y>::iterator iterateur = (*(itsDefaultTimes->Getvalue(k1,k2))).begin(); iterateur != (*(itsDefaultTimes->Getvalue(k1,k2))).end(); ++ iterateur)
		{
			double em = 0 ;
			t = iterateur->tau ;
			iName = iterateur->id ;
			int e = Eta(t) ;
			for(unsigned int j = deb; j < e ; j++)
				premiumLeg +=  itsZcPay[j] * Notional * itsCouponsDates[j] * itsSpreads[j] ;

			for(unsigned int k = 0 ; k < itsIndCDO[iName].size() ; k++)
			{
				int ind = itsIndCDO[iName][k] ;

				if (itsInstrumentType[ind] == ICM_MEZ)
				{
				itsPortfolioLosses[ind]  += itsCollatLoss->Getvalue(iName,ind) ;
				itsTrancheLosses[ind]    = MAX( MIN(itsPortfolioLosses[ind],itsStrikesUp[ind]) - 
										 itsStrikesDown[ind], 0. ) ;
				}
				else
				{
				(its_NTD_ActualDef[ind])++;

				if (its_NTD_ActualDef[ind] == its_NTD_NumDef[ind])
					itsTrancheLosses[ind] = itsCollatLoss->Getvalue(iName,ind);
				}
			}
			itsCdo2PortfolioLoss = accumulate(itsTrancheLosses.begin(),itsTrancheLosses.end(),0) ;
			itsCdo2TrancheLoss   = MAX( MIN(itsCdo2PortfolioLoss,itsStrikeUp) - itsStrikeDown, 0. ) ;
			if(itsCdo2TrancheLoss>0.1)
			{
				NotionalPrev = Notional ;
				Notional = Tranche - itsCdo2TrancheLoss ;
				defaultLeg			 +=  itsZcPay[e] * (NotionalPrev - Notional) ;
				if(e>0)
					em				  = itsYFPayDates[e-1] ;
				premiumLeg			 +=  itsZcPay[e] * (NotionalPrev - Notional) * (t-em) * itsSpreads[e] ; 

			}
			deb = e ;
		}
		for(unsigned int l = deb; l < itsNbPeriods ; l++)
			premiumLeg			 +=  itsZcPay[l] * Notional * itsCouponsDates[l] * itsSpreads[l] ; 
	
		SetPortfolioLosses(0.) ;	
		SetTrancheLosses(0.) ;
		itsCdo2PortfolioLoss = 0. ;
		itsCdo2TrancheLoss = 0. ;
		Notional = Tranche ;
		NotionalPrev = Tranche ;

		Cube->SetElt(k1,k2,i,premiumLeg-defaultLeg);
	}
	}

	for(k1 = 0; k1 < its_sz_dt_i ; k1++)
		for(k2 = 0; k2 < its_sz_dt_j ; k2++)
		{
			price = 0.;

			for(int l = 0; l < nPaths; l++)
			{
				price += Cube->Elt(k1,k2,l);
			}

			price /= nPaths;

		GetSensiManager()->GetSensiMatrix(ICM_FAST_SPREAD_TYPE)->GetMatrix()->SetValue(k1,k2,price-initialprice);
		}


	if (Cube) delete Cube;

    }	
   

	return (price) ;
}
**/ 

void ICM_Pricer_MC_Cdo2::View(char* id, FILE* ficOut)
{
		FILE* fOut;
		char  fOutName[200];
		int i = 0;
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
			
	
		if (itsField)
		{
			fprintf(fOut, "\n ======> Field Matrix :\n\n");
			int size = itsField->Getnbrows();
			int k =0;

			for (i = 0; i<size; i++)
			{
			fprintf(fOut, "%s\t", itsUnionNames[i].c_str()); 
				for (k = 0; k<itsNbTranches; k++)
			{
				fprintf(fOut, "\t %d ", itsField->Getvalue(i,k)); 
			}
				fprintf(fOut, "\n");
			}
		}

		if (itsCollatLoss)
		{
			fprintf(fOut, "\n ======> Collat Loss Matrix :\n\n");

			int size = itsCollatLoss->Getnbrows();
			int k =0;

			for (i = 0; i<size; i++)
			{
			fprintf(fOut, "%s\t", itsUnionNames[i].c_str()); 
				for (k = 0; k<itsNbTranches; k++)
			{
				fprintf(fOut, "\t %f ", itsCollatLoss->Getvalue(i,k)); 
			}
				fprintf(fOut, "\n");
			}
		}

		fprintf(fOut, "\n");
		fprintf(fOut, "\t\t\t -------------------------------------------------------------------- \n");
		fprintf(fOut, "\t\t\t ------------------------- average spread on CDO at Maturity -------- \n");
		fprintf(fOut, "\t\t\t -------------------------------------------------------------------- \n\n");
		ICM_Cdo2* cdo2 = dynamic_cast<ICM_Cdo2*>(GetSecurity());
		if (cdo2){
			for ( int j=0; j< cdo2->GetPortfolio()->GetNbSec(); j++){
				ICM_Mez* mez = dynamic_cast<ICM_Mez*>(cdo2->GetPortfolio()->GetSecurities()[j]);
				int NbIssuers = mez->GetCollateral()->GetNbIssuers();	
				double averageSpread = 0.0;
				for (int k = 0; k< NbIssuers; k++) {
					string name = string(mez->GetCollateral()->GetIssuersLabels(k));
					
					ICM_ModelMultiCurves* MMC = (ICM_ModelMultiCurves*) GetModel();
					double spread = MMC->GetDefaultCurve(name)->ImpliedSpreadInterpol(mez->GetEndDateNA());
					averageSpread += spread;	
				}
				averageSpread /= NbIssuers;
				fprintf(fOut, "\t CDO %d, average spread at maturity : %lf\n",j, averageSpread);
			}
		}
		fprintf(fOut, "\n");

		if ( ficOut == NULL )
		{
			fclose(fOut);
		}
}


/*----------------------------------------------------------------------------*
  Compute Credit Spread
*----------------------------------------------------------------------------*/ 
double ICM_Pricer_MC_Cdo2::ComputeSpread(const double& MtM )
{
	if (GetSpreadFlg()) return GetSpread();

	double Breakeven =0.;

	double price=  Price(qCMPPRICE);
	double FeeNPV_Unity = GetFeeLegPrice_Unity();

	//double SpdCour	= ((ICM_Cds*)GetSecurity())->GetFeeLeg()->GetCreditSpread();

	//FeeNPV /= (SpdCour/100) ;

	Breakeven = 10000.* GetDefLegPrice()/FeeNPV_Unity;
	SetSpread(Breakeven);

	return (Breakeven);
}



// -------------------------------------------------------------
// Compute DTR Method For CDO²
// -------------------------------------------------------------
double ICM_Pricer_MC_Cdo2::ComputeDTR(int NbDefaults, char* S_or_L, char** labels, double* Recoveries)
{
    double sensitivity =0.;

	ICM_ModelMultiCurves* DefModel = (ICM_ModelMultiCurves*) GetModel();

	ARM_Date ExecutionDate = GetModel()->GetStartDate();

	int i = 0, j = 0;
	double result = 0.0;
	double initialprice = 0.0;
	double modifiedprice = 0.0;

	 
    {
		if (GetInitialPriceFlg())
			initialprice = GetInitialPrice();
		else
			initialprice = Price(qCMPPRICE);

		ResetPricer();

		// if(!strcmp(S_or_L,"S")) initialprice = -initialprice ;
		if (S_or_L=="S") initialprice = -initialprice ;

		ICM_Cdo2* Security = (ICM_Cdo2*) GetSecurity();
		ICM_Portfolio* portfolio = (ICM_Portfolio*) Security->GetPortfolio();

		int NbSec = portfolio->GetNbSec();

		ICM_Cdo2* InitCDO2 = (ICM_Cdo2*) Security;
		ICM_Portfolio* InitPortfolio = (ICM_Portfolio*) InitCDO2->GetPortfolio();

		ICM_Cdo2* NewCDO2 = (ICM_Cdo2*) InitCDO2->Clone() ;
		
		double LowThreesholds  = 0. ;
		double HighThreesholds = 0. ;

		double* Losses = new double [NbDefaults] ;
		
		for (i = 0; i < NbDefaults ; i++)
		{
			Losses[i] = InitCDO2->GetCollateral()->GetIssuersNotional(
				InitCDO2->GetCollateral()->getIssuerPosition(labels[i])) * (1.0 - Recoveries[i]);
		}

		ARM_Vector* Actual_Loss = new ARM_Vector (NbSec, 0.);
		ARM_Security** securities= new ARM_Security* [NbSec];

		double MasterLoss = 0.; // Loss of the Master CDO due to the Default of the name
		bool Test = false;

		for(i = 0 ; i < NbSec ; i++)
		{
			ARM_Security* Sec = portfolio->GetSecurity(i);
			ICM_Mez* InitMez = (ICM_Mez*) Sec;
			ICM_Mez* MEZ = (ICM_Mez*) InitMez->Clone() ;
			
			Test = false ;
			
			for (int j = 0 ; j < NbDefaults ; j++)
			{
				Test = false;

				for (int k = 0; k < MEZ->GetCollateral()->GetNbIssuers(); k++)
				{
					if (!strcmp(MEZ->GetCollateral()->GetIssuersLabels(k).c_str(),labels[j]))
					{
						Test = true ;
					}
				}
				
				LowThreesholds  = MEZ->GetSubAmount(MEZ->GetFeeLeg()->GetStartDate());
				HighThreesholds = MEZ->GetSubAmount(MEZ->GetFeeLeg()->GetStartDate()) + MEZ->GetMezzAmount(MEZ->GetFeeLeg()->GetStartDate()) ;

				if(Test == true)
				{
					double aux = 0.;
					MEZ->ExcludeIssuer(labels[j]);
					if(Losses[j] > HighThreesholds)
					{
						aux = Actual_Loss->Elt(i) + HighThreesholds - LowThreesholds ;
						Actual_Loss->InitElt(i, aux );
						MEZ->SetSubAmount(0.) ;
						MEZ->SetMezzAmount(0.);
					}
					
					else
					{
						if (LowThreesholds < Losses[j]  && Losses[j]  <= HighThreesholds )
						{
							MEZ->SetSubAmount(0.);
							MEZ->SetMezzAmount(HighThreesholds - Losses[j]);
							aux = Actual_Loss->Elt(i) + Losses[j] - LowThreesholds ;
							Actual_Loss->InitElt(i, aux) ;
						}
						else
						{
							MEZ->SetSubAmount(LowThreesholds - Losses[j] );
							MEZ->SetMezzAmount(HighThreesholds - LowThreesholds );
						}
					}
					
				}
			}

			securities[i] = MEZ;
			MasterLoss += Actual_Loss->Elt(i) ;
		}


		ICM_Portfolio* NewPortfolio = new ICM_Portfolio ( securities, NbSec);
		NewCDO2->SetPortfolio(NewPortfolio);
		
		double A = InitCDO2->GetSubAmount(InitCDO2->GetFeeLeg()->GetStartDate());
		double B = InitCDO2->GetMezzAmount(InitCDO2->GetFeeLeg()->GetStartDate()) + A;
		double Effective_Loss = 0. ;
		
		if(MasterLoss > B)
		{
			Effective_Loss = B-A;
			modifiedprice = 0.;
		}
		
		else
		{
			if (A<MasterLoss && MasterLoss<=B)
			{
				NewCDO2->SetSubAmount(0.);
				NewCDO2->SetMezzAmount(B - MasterLoss);
				Effective_Loss = MasterLoss - A;
			}
			else
			{
				NewCDO2->SetSubAmount(A - MasterLoss);
				NewCDO2->SetMezzAmount(B - A );
			}
			SetSecurity(NewCDO2);
			Set(NewCDO2, DefModel, GetParameters(),GetAsOfDate(),GetNbPaths()) ;
			
			modifiedprice = ComputePrice(qCMPPRICE);
			
			// if(!strcmp(S_or_L,"S")) modifiedprice = -modifiedprice;
			if(S_or_L=="S") modifiedprice = -modifiedprice;
			
			//if (NewCDO2) delete NewCDO2; 
		}
		// memory delete
		for (int i=0; i<NbSec; i++){
			if( securities[i]!= NULL) {
				//!\ no * = NULL 


					delete securities[i];
			}
		}
		if (NewCDO2) delete NewCDO2; 

		SetSecurity(InitCDO2) ;
		InitCDO2->SetPortfolio(InitPortfolio);
		Set(InitCDO2, DefModel, GetParameters(),GetAsOfDate(),GetNbPaths()) ;				
		
		if (!result)
		{
			//Ajouter le CF Effective_Loss si la tranche est touchée
			// if(!strcmp(S_or_L,"S"))
			if (S_or_L=="S")
				result = modifiedprice - initialprice + Effective_Loss;
			else result=modifiedprice - initialprice - Effective_Loss; //On est Long par défaut
		}
		
		else result = -99999999.0;

		ResetPricer();
		
		Set(GetSecurity(), DefModel, GetParameters(),GetAsOfDate(),GetNbPaths());

		if (!result)
			result=modifiedprice - initialprice;

		if (Losses)
			delete [] Losses ;
		Losses = NULL;

		if (Actual_Loss)
			delete Actual_Loss ;
		Actual_Loss = NULL;
	}

    
	return (result);
}


// -------------------------------------------------------------
// Generate ¨Perturbates Defaults Probabilities 
// -------------------------------------------------------------

/** 
17783 
void ICM_Pricer_MC_Cdo2::PerturbDefaultCurves()
{
	int nocurve = -1;
	int i =0,j=0;

	if (GetSensiManager()==NULL) return;

	const std::vector<std::string>& labels = (GetSensiManager()->GetSensiMatrix(ICM_FAST_SPREAD_TYPE))->GetMktData();
	char** tenors = (GetSensiManager()->GetSensiMatrix(ICM_FAST_SPREAD_TYPE))->GetBumpProfile();

	const ICM_DefaultCurve** VectorDefCurves = NULL;	
	ICM_DefaultCurve** ModifVectorDefCurves = NULL;	
	ICM_ModelMultiCurves * model = (ICM_ModelMultiCurves*)GetModel() ;

	its_sz_dt_i = (GetSensiManager()->GetSensiMatrix(ICM_FAST_SPREAD_TYPE))->GetMktDataSize();
	its_sz_dt_j = (GetSensiManager()->GetSensiMatrix(ICM_FAST_SPREAD_TYPE))->GetBumpProfileSize();

	VectorDefCurves = new const ICM_DefaultCurve*[its_sz_dt_i];

	for (i=0; i<its_sz_dt_i; i++)
		VectorDefCurves[i] = model->GetDefaultCurve(labels[i]);


	ICM_QMatrix<ICM_DefaultCurve*>* ShiftedCurves = new ICM_QMatrix<ICM_DefaultCurve*>(its_sz_dt_i,its_sz_dt_j,NULL);

	

	for (i=0; i<its_sz_dt_j; i++)
	{
		ModifVectorDefCurves = BumpPortfolio_(VectorDefCurves, 
											 its_sz_dt_i,
											 ICMSPREAD_TYPE,
											 tenors[i], 
											 "NONE",
											 nocurve);								

		for (int j=0; j<its_sz_dt_i; j++) ShiftedCurves->SetValue(j,i,ModifVectorDefCurves[j]);

		delete[] ModifVectorDefCurves;
	}

	if (VectorDefCurves)
		delete[] VectorDefCurves;
	VectorDefCurves = NULL;

	(GetSensiManager()->GetSensiMatrix(ICM_FAST_SPREAD_TYPE))->SetObjectMatrix((void*)ShiftedCurves);
}
17783 
**/ 

// -----------------------------------------------------------------
// Compute Sensitivity Method For Baskets
// -----------------------------------------------------------------
double ICM_Pricer_MC_Cdo2::ComputeSensitivity(const qSENSITIVITY_TYPE& typesensi, 
											  const std::string&  plot, const std::string&  label, 
											  double  epsvalue, double  epsilonGamma // useless
											  )
{
    double sensitivity =0.;

	ICM_ModelMultiCurves* DefModel = (ICM_ModelMultiCurves*) GetModel();
	int nbcurves = DefModel->GetNbDefCurves();
	int nocurve = -1;

	ARM_Date ExecutionDate = GetModel()->GetStartDate();

	int i = 0;
	double result = 0.,initialprice = 0.,modifiedprice = 0.;

	//Perte subie effectivement par la tranche
	double Actual_Loss = 0. ;

	if (typesensi == ICMRECOVERY_BAR_TYPE)
	{
		ICMTHROW(ERR_UNIMP_METHOD_CALL,"Unimplemented ICMRECOVERY_BAR_TYPE method");
	}

     
    {

		if (typesensi != ICMBETA_WITH_SPREAD_SHIFT)
		{
		if (GetInitialPriceFlg())
			initialprice = GetInitialPrice();
		else
			initialprice = Price(qCMPPRICE);

		ResetPricer();
		}

		switch (typesensi)
		{
			case ICMRECOVERY_TYPE :
			case ICM_SAMECORRELATION_TYPE :
			case ICMCORRELATION_TYPE :
			case ICM_SAMEBETA_TYPE :
			case ICMBETA_TYPE :
			case ICMIRCURVE_TYPE :
			case ICMIRCURVE_WITHOUT_DEFCURVE_TYPE :
			case ICMSPREAD_TYPE :
			case ICMSPRELSHIFT_TYPE :
			case ICM_DTR_TYPE :
 			{
	
				ICM_ModelMultiCurves* ModelDef2 = DefModel->GenerateShiftModel(typesensi,
																			   plot, 
																			   label,
																			   nocurve,
																			   epsvalue);
				
				SetModel(ModelDef2);	

				BeforePrice(ModelDef2);
				modifiedprice = ComputePrice(qCMPPRICE);

				SetModel(DefModel); //On reset le model initial

				if (ModelDef2)
					delete ModelDef2;
				ModelDef2 = NULL;
			}
			break;
			// -------------------------------------------------------------
			case ICMBETA_WITH_SPREAD_SHIFT : //Beta parallel shift with spread shift
			// -------------------------------------------------------------
			{
				vector<double> betas ;
			}
			break;
			// -------------------------------------------------------------
			case ICM_ISSUER_DEFAULT :  //Default to recovery
			// -------------------------------------------------------------
			{
				// if(!strcmp(plot,"S")) initialprice = -initialprice ;
				if(plot=="S") initialprice = -initialprice ;

				ICM_Cdo2* Security = (ICM_Cdo2*) GetSecurity();
				ICM_Portfolio* portfolio = (ICM_Portfolio*) Security->GetPortfolio();
				int NbSec = portfolio->GetNbSec();

				ICM_Cdo2* InitCDO2 = (ICM_Cdo2*) Security;
				ICM_Portfolio* InitPortfolio = (ICM_Portfolio*) InitCDO2->GetPortfolio();

				ICM_Cdo2* NewCDO2 = (ICM_Cdo2*) InitCDO2->Clone() ;
				
				double LowThreesholds  = 0. ;
				double HighThreesholds = 0. ;

				int issuerRank = InitCDO2->GetCollateral()->getIssuerPosition(label); 
				double NOT = InitCDO2->GetCollateral()->GetIssuersNotional(issuerRank);
				double LR = 1. - epsvalue; //LR to be paid in case of default
				double alpha = NOT*LR;

				ARM_Vector* Actual_Loss = new ARM_Vector (NbSec, 0.);
				ARM_Security** securities= new ARM_Security* [NbSec];

				int Compt = 0;

				double Loss = 0.; // Loss of the Master CDO due to the Default of the name
				bool Test = false;
				bool ToInclude = true;

				for(int i = 0 ; i < NbSec ; i++)
				{
					ARM_Security* Sec = portfolio->GetSecurity(i);
					ToInclude = true ;

					if (Sec->GetName() == ICM_MEZ)
					{
						ICM_Mez* InitMez = (ICM_Mez*) Sec;
						LowThreesholds  = InitMez->GetSubAmount(InitMez->GetFeeLeg()->GetStartDate());
						HighThreesholds = InitMez->GetSubAmount(InitMez->GetFeeLeg()->GetStartDate()) + InitMez->GetMezzAmount(InitMez->GetFeeLeg()->GetStartDate()) ;
						ICM_Mez* MEZ = (ICM_Mez*) InitMez->Clone() ;
						
						Test = false ;
						for (int k = 0; k < MEZ->GetCollateral()->GetNbIssuers(); k++)
						{
							if (!strcmp(MEZ->GetCollateral()->GetIssuersLabels(k).c_str(),label.c_str()))
							{
								Test = true ;
							}
						}
						
						if(Test == true)
						{
							MEZ->ExcludeIssuer(label.c_str());
							if(alpha >= HighThreesholds)
							{
								Actual_Loss->InitElt(i, HighThreesholds - LowThreesholds);
								MEZ->SetSubAmount(0.) ;
								MEZ->SetMezzAmount(0.);
								ToInclude = false ;
								delete MEZ;
							}
							else
							{
								if (LowThreesholds < alpha && alpha <= HighThreesholds )
								{
									MEZ->SetSubAmount(0.);
									MEZ->SetMezzAmount(HighThreesholds - alpha);
									Actual_Loss->InitElt(i, alpha - LowThreesholds ) ;
								}
								else
								{
									MEZ->SetSubAmount(LowThreesholds - alpha);
									MEZ->SetMezzAmount(HighThreesholds - LowThreesholds );
								}
							}
							Loss += Actual_Loss->Elt(i) ;
						}
						if (ToInclude)
						{
							securities[Compt] = MEZ;
							Compt++;
						}
					}
					else if (Sec->GetName() == ICM_NTD)
					{
						ICM_Nthtd* InitNTD = (ICM_Nthtd*)Sec; 
						double NOT = 0.;
						ICM_Nthtd* NewNTD = (ICM_Nthtd*) InitNTD->Clone();

						Test = false ;
						for (int k = 0; k < InitNTD->GetCollateral()->GetNbIssuers(); k++)
						{
							if (!strcmp(InitNTD->GetCollateral()->GetIssuersLabels(k).c_str(),label.c_str()))
							{
								Test = true ;
							}
						}

						if(Test == true)
						{
							NOT = InitNTD->GetCollateral()->GetIssuersNotional(InitNTD->GetCollateral()->getIssuerPosition(label));

							NewNTD->ExcludeIssuer(label);
							if (InitNTD->GetFirstNumDefault() == 1)
							{
								Actual_Loss->InitElt(i, NOT * LR);
								NewNTD->SetFirstNumDefault(0);
								NewNTD->SetLastNumDefault(0) ;
								ToInclude = false ;
								delete NewNTD;
							}
							else
							{			
								NewNTD->SetFirstNumDefault(MAX(NewNTD->GetFirstNumDefault()-1,0.));
								NewNTD->SetLastNumDefault(MAX(NewNTD->GetLastNumDefault()-1,0.));
							}
							Loss += Actual_Loss->Elt(i) ;
						}

						if (ToInclude)
						{
							securities[Compt] = NewNTD;
							Compt++;
						}
					}
				}

				double Effective_Loss = 0. ;
				double A = InitCDO2->GetSubAmount(InitCDO2->GetFeeLeg()->GetStartDate());
				double B = InitCDO2->GetMezzAmount(InitCDO2->GetFeeLeg()->GetStartDate()) + A;
				
				if(Compt == 0)
				{
					modifiedprice = 0. ;
					if (Loss > B)
						Effective_Loss = B -A;
					else
					{
						if (Loss>=A)
							Effective_Loss = Loss - A;
					}						
				}
				else
				{
					ICM_Portfolio* NewPortfolio = new ICM_Portfolio (securities, Compt);
					
					NewCDO2->SetPortfolio(NewPortfolio);
					
					if(Loss > B)
					{
						Effective_Loss = B-A;
						modifiedprice = 0.;
					}
					else
					{
						if (A<Loss && Loss<=B)
						{
							NewCDO2->SetSubAmount(0.);
							NewCDO2->SetMezzAmount(B - Loss);
							Effective_Loss = Loss - A;

						}
						else
						{
							NewCDO2->SetSubAmount(A - Loss);
							NewCDO2->SetMezzAmount(B - A );
						}
						SetSecurity(NewCDO2);
						Set(NewCDO2, DefModel, GetParameters(),GetAsOfDate(),GetNbPaths()) ;
						
						modifiedprice = ComputePrice(qCMPPRICE);
						
						// if(!strcmp(plot,"S")) modifiedprice = -modifiedprice;
						if(plot=="S") modifiedprice = -modifiedprice;
						
					}

					SetSecurity(InitCDO2) ;
					InitCDO2->SetPortfolio(InitPortfolio);
					Set(InitCDO2, DefModel, GetParameters(),GetAsOfDate(),GetNbPaths()) ;				

					for (int il=0; il<Compt; il++) delete securities[il];
				}
				
				if (NewCDO2) delete NewCDO2; 
				if (Actual_Loss) delete Actual_Loss;

				if (!result)
				{
					//Ajouter le CF Effective_Loss si la tranche est touchée
					// if(!strcmp(plot,"S"))
					if(plot=="S")
						result = modifiedprice - initialprice + Effective_Loss;
					else result=modifiedprice - initialprice - Effective_Loss; //On est Long par défaut
				}
				
				else result = -99999999.0;
			}
			break;
			default :
			result = -99999999.0;
		}

	ResetPricer();

	if (!result) result=modifiedprice - initialprice;

	}
  

	return (result);
}

void ICM_Pricer_MC_Cdo2::ComputeZc()
{
	ICM_ModelMultiCurves * model = (ICM_ModelMultiCurves*)GetModel() ;
	ICM_Cdo2 * cdo2 = (ICM_Cdo2*) GetSecurity() ;
	ICM_Security* security = cdo2->GetFeeLeg()->GetCreditInfos();
	itsZcPay.resize(itsNbPeriods) ;
	for(int i = 0; i < itsNbPeriods ; i++)
		itsZcPay[i] = model->GetZeroCurve()->DiscountPrice((ARM_Date)security->GetPayDates().Elt(i)) ;
}
void ICM_Pricer_MC_Cdo2::ComputeSpreads()
{
	ICM_Cdo2 * cdo2 = (ICM_Cdo2*) GetSecurity() ;
	ICM_Security* security = cdo2->GetFeeLeg()->GetCreditInfos();
	const ARM_Vector & Spreads = security->GetCouponRates() ;
	itsSpreads.resize(itsNbPeriods) ;
	for(unsigned int p = 0; p < Spreads.size(); p++)
		itsSpreads[p] = Spreads.Elt(p)/100. ;
}
void ICM_Pricer_MC_Cdo2::ComputeAccrualDates()
{
	ICM_Cdo2 * cdo2 = (ICM_Cdo2*) GetSecurity() ;
	ICM_ModelMultiCurves * model = (ICM_ModelMultiCurves*) GetModel() ;
	ARM_Date AsOf = model->GetStartDate() ;
	ICM_Security* security = cdo2->GetFeeLeg()->GetCreditInfos();
	security->ComputeYF(AsOf);

	const ARM_Vector & InterestDays  = security->GetYFInterestDays() ;
	const ARM_Vector & AccStartDates  = security->GetAccStartDates();
	const ARM_Vector & AccEndDates  = security->GetAccEndDates();

	itsAccrualDates.resize(itsNbPeriods) ;
	itsCouponsDates.resize(itsNbPeriods) ;

	for(int i = 0; i < InterestDays.size(); i++)
	{
		itsCouponsDates[i] = InterestDays.Elt(i) ;
		itsAccrualDates[i] = 0.;

		if ((AsOf.GetJulian()>AccStartDates.Elt(i)) &&  (AsOf.GetJulian()>=AccEndDates.Elt(i)))
		{
			itsCouponsDates[i] = 0.;
			itsAccrualDates[i] = 0.;
		}
		else if ((AsOf.GetJulian()>AccStartDates.Elt(i)) &&  (AsOf.GetJulian()<AccEndDates.Elt(i))) 
		{
			itsCouponsDates[i] = CountYears(security->GetAccrualBasis(),AsOf.GetJulian(),AccEndDates.Elt(i)); //broken period
			itsAccrualDates[i] = InterestDays.Elt(i) - itsCouponsDates[i];
		}
	}
}
void ICM_Pricer_MC_Cdo2::ComputeYFPayDates()
{
	ICM_Cdo2 * cdo2 = (ICM_Cdo2*) GetSecurity() ;
	ICM_Security* security = cdo2->GetFeeLeg()->GetCreditInfos();

	ICM_ModelMultiCurves * model = (ICM_ModelMultiCurves*) GetModel() ;
	ARM_Date AsOf = model->GetStartDate() ;
	
	security->ComputeYF(AsOf);
	const ARM_Vector & YFPayDates = security->GetYFPayDates() ;	
	itsYFPayDates.resize(itsNbPeriods) ;
	for(unsigned int p = 0; p < YFPayDates.size(); p++)
		itsYFPayDates[p] = YFPayDates.Elt(p) ;
}

void ICM_Pricer_MC_Cdo2::ComputeUnionNames()
{
	ARM_Vector  notional ;
	ICM_Cdo2 * cdo2 = (ICM_Cdo2*) GetSecurity() ;

	// FreePointerTabChar(itsUnionNames,itsNbNames);	

	cdo2->GetPortfolio()->GetIssuersDatas(itsUnionNames,notional) ;
// 	delete notional;
}
/**
17783 	
void ICM_Pricer_MC_Cdo2::CptCorrespondanceUnionNames_Sensis()
{
	int i =0,k=0;
	bool status = true;
	itsCorrespondanceUnionNames_Sensis.resize(itsNbNames);
	ICM_Sensi2D* sensi2D = GetSensiManager()->GetSensiMatrix(ICM_FAST_SPREAD_TYPE);

	for(i = 0; i<sensi2D->GetMktDataSize();i++) 
	{
		status = true;

		for(k = 0; k<itsNbNames;k++) 
		{	
			if (strcmp(itsUnionNames[k].c_str(),sensi2D->GetMktData(i).c_str()) == NULL) 
			{
			status = false;
			break;
			}
		}

		if (status)
			itsCorrespondanceUnionNames_Sensis[i] = -1;
		else
			itsCorrespondanceUnionNames_Sensis[i] = k;
	}
}	
17783 	**/ 

void ICM_Pricer_MC_Cdo2::GenerateDefaultTimes_Betas_Standart()
{
	itsDefaultTimes_Standart.clear() ;
	ICM_Root_Generator* Gen = GetGenerator() ;
	ICM_ModelMultiCurves* model = (ICM_ModelMultiCurves*)GetModel() ;
	Gen->CommonFactors() ;

	double aux = 0. ;

	for(unsigned int i = 0; i < itsNbNames; i++)
	{	
		double ri = Gen->generateRandom(i) ;
		if(ri < itsBarriers_Standart[i])
		{
			aux = NAG_cumul_normal(ri);
			double  t = model->GetDefaultCurve(itsUnionNames[i])->DefProbInverse(aux);
			itsDefaultTimes_Standart.insert(Y(i,t)) ;
		}
	}
}	

// ----------------------------------------------------------------------------
// Generation of default times with Betas in complete case
// ----------------------------------------------------------------------------
void ICM_Pricer_MC_Cdo2::GenerateDefaultTimes_Betas_Complete()
{
	int i = 0,j=0,k=0,l=0,k_=0;
	set<Y>* DefaultTime = NULL;
	itsDefaultTimes_Standart.clear();
	vector<double> _t; _t.clear();_t.resize(itsNbNames);
	vector<int> _i ; _i.clear();_i.resize(itsNbNames);
	vector<double> ri_;ri_.clear();ri_.resize(itsNbNames);
	double aux = 0. ,t=0;
	// // 17783  ICM_Sensi2D* sensi2D = NULL;

	ICM_Root_Generator * Gen = GetGenerator() ;
	ICM_ModelMultiCurves * model = (ICM_ModelMultiCurves*)GetModel() ;
	Gen->CommonFactors() ;

	

	//Reinitialisation de la matrice de vecteur des temps de défauts
/** 
17783 	 
	sensi2D = GetSensiManager()->GetSensiMatrix(ICM_FAST_SPREAD_TYPE);
	if (GetSensiManager())
	if (sensi2D)
	{
	its_sz_dt_i = sensi2D->GetMktDataSize();
	its_sz_dt_j = sensi2D->GetBumpProfileSize();

		if (itsDefaultTimes)
		{
			for (i =0; i<its_sz_dt_i; i++)
				for (j =0; j<its_sz_dt_j; j++)
				{
				 if (itsDefaultTimes->Getvalue(i,j))
				   (itsDefaultTimes->Getvalue(i,j))->clear();
				 else
					{
					DefaultTime = new set<Y>;
					itsDefaultTimes->SetValue(i,j,DefaultTime);
					}	
					
				}
		}
	}
17783 	**/ 
	for(i = 0; i < itsNbNames; i++)
	{	
		ri_[i] = Gen->generateRandom(i) ;
		_t[i] = -1.;_i[i] = i;

		if(ri_[i] < itsBarriers_Standart[i])
		{
			aux = NAG_cumul_normal(ri_[i]);
			t = model->GetDefaultCurve(itsUnionNames[i])->DefProbInverse(aux);
			_t[i] = model->GetDefaultCurve(itsUnionNames[i])->DefProbInverse(aux);
			itsDefaultTimes_Standart.insert(Y(i,_t[i])) ;
		}
	}

/**
17783 	
	for(i = 0; i < its_sz_dt_i ; i++)
	for(j = 0; j < its_sz_dt_j ; j++)
	{

		DefaultTime = itsDefaultTimes->Getvalue(i,j);
		k_=itsCorrespondanceUnionNames_Sensis[i];

		for(k = 0; k<itsNbNames;k++) 
		{
			if ((k != k_) && (_t[k]>=0.))
				DefaultTime->insert(Y(_i[k],_t[k]));
		}

		if(ri_[k_] < (*(itsBarriers->Getvalue(i,j)))[k_])
		{
			aux = NAG_cumul_normal(ri_[k_]);
			t = ((ICM_QMatrix<ICM_DefaultCurve*>*)sensi2D->GetObjectMatrix())->Getvalue(i,j)->DefProbInverse(aux);
			DefaultTime->insert(Y(_i[k_],t));
		}

	}

	17783 	
	**/ 
}	


void ICM_Pricer_MC_Cdo2::AffectLossForCrossSub(const int& NoName,
													const int& NoCdo,
													ICM_QMatrix<double>* CollatLoss,
													const vector<double>& Strike_Down,
													const vector<double>& Strike_Up,
													vector<double>& Loss,
													const int& FirstNoCdo,
													const double& residu)
{
	if (!residu) return;
	double Residu = 0.;
	double Vloss = (residu == CREDIT_DEFAULT_VALUE ? CollatLoss->Getvalue(NoName,NoCdo) : residu);
	bool AllFull = true;

	for (int i=0;i<itsNbTranches;i++) {if (Loss[i]<Strike_Down[i]) AllFull=false;}

	if (AllFull) 
	{
		//cas cross subordination capé 
		if ((Loss[FirstNoCdo]+Vloss)<=Strike_Up[NoCdo])
			{Loss[FirstNoCdo]+=Vloss;}
		else
			{Loss[FirstNoCdo]=Strike_Up[NoCdo];}
		return;
	}		

	if ((Loss[NoCdo]+Vloss)<=Strike_Down[NoCdo])
	{
		Loss[NoCdo]+=Vloss;
		return;
	}
	else if ((Loss[NoCdo]<Strike_Down[NoCdo]) && ((Loss[NoCdo]+Vloss)>Strike_Down[NoCdo]))
	{
		Residu = Loss[NoCdo]+Vloss-Strike_Down[NoCdo];
		Loss[NoCdo]=Strike_Down[NoCdo];
	}
	else Residu = Vloss;

	int IncNoCdo = NoCdo;
	IncNoCdo++; //on passe au cdo suivant
	if (IncNoCdo==itsNbTranches) {IncNoCdo=0;}

	AffectLossForCrossSub(NoName, IncNoCdo, CollatLoss, Strike_Down, Strike_Up, Loss,FirstNoCdo,Residu);
}

//  
void 
ICM_Pricer_MC_Cdo2::BeforePrice(ICM_ModelMultiCurves* Model,qSENSITIVITY_TYPE type ) 
{ 
	Set(GetSecurity(), Model, GetParameters(),GetAsOfDate(),GetNbPaths()); 
}
