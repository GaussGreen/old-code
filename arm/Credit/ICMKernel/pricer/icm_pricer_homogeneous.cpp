#include "firsttoinc.h"
#include "ICMKernel\pricer\ICM_Pricer_homogeneous.h"
#include "ICMKernel\inst\icm_nthtd.h"
#include "ICMKernel\inst\icm_mez.h"
#include "ICMKernel\glob\icm_smile_correlation.h"
#include "ICMKernel/crv/icm_distriblosscalculator.h"
#include "ICMKernel\inst\icm_collateral.h"
#include "ICMKernel\mod\icm_lossunits.h"
#include "ICMKernel\mod\modelmulticurves.h"
#include "ARMKernel\crv\volint.h"
#include "ICMKernel\util\icm_rootfinder1d.h"


void ICM_Pricer_Distrib::Set(ARM_Security *sec, ARM_Model *mod,const ICM_Parameters&params,const ARM_Date&asof)
{
	if (mod->GetName()!=ICM_MODELMULTICURVES)
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "ERROR : pricer not compatible with this security or this model ");

	ICM_Pricer_Basket::Set(sec, mod,params,asof);

	const ICM_Parameters& FlowsMatrix = GetParameters();
	ARM_Vector* MINLOSSUNIT = FlowsMatrix.GetColVect("MINLOSS");
	if (MINLOSSUNIT) itsMinLossUnit = MINLOSSUNIT->Elt(0);
}

//---------------------------------------------------------------
// Détermination du % en risque
//---------------------------------------------------------------
double ICM_Pricer_Distrib::GetTranche(const ARM_Date& date)
{
	ICM_ModelMultiCurves* model = (ICM_ModelMultiCurves*) GetModel();
	ICM_Ftd* ftd = (ICM_Ftd*) GetSecurity();

	double tranche_down = 0.;
	double tranche_up = 0.;
	double Recovery = 0.;

	int nbnames = ftd->GetCollateral()->GetNbIssuers();

	if (ftd->GetBinaryFlg())
		Recovery = ftd->GetBinary();
	else
	if ((ftd->GetName() == ICM_NTD) || (ftd->GetName() == ICM_FTD))
	{
		// Gestion du cas Recovery Non Homogene pour les Ntd	
		double sum=0.;
		double* Pdefault = new double[nbnames];
		double Recov_weight = 0.;
		double yearterm = (ftd->GetEndDateNA() - model->GetStartDate())/365.;
		double RecovCoef = 1.;

		//Test Recov coef
		if (ftd->GetCollateral()->GetRecovCoefFlg())
			RecovCoef = ftd->GetCollateral()->GetRecovCoef();

		for (int i=0; i<nbnames; i++)
		{
			Pdefault[i] = model->GetDefaultCurve(ftd->GetCollateral()->GetIssuersLabels(i))->DefaultProba(yearterm);
			sum+=Pdefault[i];
			Recov_weight+=Pdefault[i]*MIN(model->GetRecoveryRate(ftd->GetCollateral()->GetIssuersLabels(i))*RecovCoef,1);
		}
		Recovery=Recov_weight/sum;

		if (Pdefault) delete[] Pdefault;
		//Recovery = AvgRecovery(ftd->GetIssuersLabels(),nbnames);
	}

	switch (ftd->GetName())
	{
		case ICM_NTD:
		{
		ICM_Nthtd* ntd = (ICM_Nthtd*) ftd;
		tranche_down = ((double)(ntd->GetFirstNumDefault()-1)) * (1. - Recovery);
		tranche_up = ((double)(ntd->GetLastNumDefault())) * (1. - Recovery);
		break;
		}
		case ICM_FTD:
		{
		tranche_up = (1. - Recovery);
		break;
		}
		case ICM_MEZ : 
		case ICM_CDO2 : 
		{
		ICM_Mez* mez = (ICM_Mez*)ftd;
		tranche_down = mez->GetSubAmount((ARM_Date)date);
		tranche_up = 1. + mez->GetSubAmount((ARM_Date)date);
		break;
		}
	}

	return (tranche_up-tranche_down);
}	

//---------------------------------------------------------------
// Détermination % Tranche Up
//---------------------------------------------------------------
double ICM_Pricer_Distrib::GetTranche_Up(const ARM_Date& date)
{
	ICM_ModelMultiCurves* model = (ICM_ModelMultiCurves*) GetModel();
	ICM_Ftd* ftd = (ICM_Ftd*) GetSecurity();

	double tranche_up = 0.,Recovery = 0.;
	ARM_Date  RefDate   = date;
	int nbnames = ftd->GetCollateral()->GetNbIssuers();
	ARM_Date AccEndDate = ftd->GetFeeLeg()->GetCreditInfos()->GetAccEndDates().Elt( ftd->GetFeeLeg()->GetCreditInfos()->GetAccEndDates().size()-1 );
 
	if (RefDate>AccEndDate) RefDate = AccEndDate;

	double Notional = ftd->GetFeeLeg()->GetCreditInfos()->Notional(RefDate);
	Notional = fabs(Notional);

	if (ftd->GetBinaryFlg())
		Recovery = ftd->GetBinary();
	else
 	if ((ftd->GetName() == ICM_NTD) || (ftd->GetName() == ICM_FTD))
	{
		// Gestion du cas Recovery Non Homogene pour les Ntd	
		double sum=0.;
		// double* Pdefault = new double[nbnames];
		double Recov_weight = 0.;
		double yearterm = (ftd->GetEndDateNA() - model->GetStartDate())/365.;
		double RecovCoef = 1.;

		//Test Recov coef
		if (ftd->GetCollateral()->GetRecovCoefFlg())
			RecovCoef = ftd->GetCollateral()->GetRecovCoef();
		for (int i=0; i<nbnames; i++)
		{
			double Pdefault  = model->GetDefaultCurve(ftd->GetCollateral()->GetIssuersLabels(i))->DefaultProba(yearterm);
			sum+=Pdefault  ;
			Recov_weight+=Pdefault  *MIN(model->GetRecoveryRate(ftd->GetCollateral()->GetIssuersLabels(i))*RecovCoef,1);
		}
		Recovery=Recov_weight/sum;

		// if (Pdefault) delete[] Pdefault;
		
		//Recovery = AvgRecovery(ftd->GetIssuersLabels(),nbnames);
	}

	switch (ftd->GetName())
	{
		case ICM_NTD:
		{
		ICM_Nthtd* ntd = (ICM_Nthtd*) ftd;
		tranche_up = ((double)(ntd->GetLastNumDefault())) * (1. - Recovery) * Notional;
		break;
		}
		case ICM_FTD:
		{
		tranche_up = (1. - Recovery) * Notional;
		break;
		}
		case ICM_MEZ : 
		case ICM_CDO2 :
		{
		ICM_Mez* mez = (ICM_Mez*)ftd;
		tranche_up = Notional + mez->GetSubAmount(RefDate);
		break;
		}
	}

	return (tranche_up);
}	


//---------------------------------------------------------------
// Détermination % Tranche Down
//---------------------------------------------------------------
double ICM_Pricer_Distrib::GetTranche_Down(const ARM_Date& date)
{
	ICM_ModelMultiCurves* model = (ICM_ModelMultiCurves*) GetModel();
	ICM_Ftd* ftd = (ICM_Ftd*) GetSecurity();

	double tranche_down = 0.,Recovery = 0.;
	ARM_Date RefDate = date;
	int nbnames = ftd->GetCollateral()->GetNbIssuers();
// 	ARM_Date AccEndDate = ftd->GetFeeLeg()->GetCreditInfos()->GetAccEndDate();
	ARM_Date AccEndDate = ftd->GetFeeLeg()->GetCreditInfos()->GetAccEndDates().Elt( ftd->GetFeeLeg()->GetCreditInfos()->GetAccEndDates().size()-1 );
	

	if (RefDate>AccEndDate) RefDate = AccEndDate;
	double Notional = ftd->GetDefLeg()->GetCreditInfos()->Notional(RefDate);

	if (ftd->GetBinaryFlg())
		Recovery = ftd->GetBinary();
	else
	if (ftd->GetName() == ICM_NTD)
	{
		// Gestion du cas Recovery Non Homogene pour les Ntd	
		double sum=0.;
		// double* Pdefault = new double[nbnames];
		double Recov_weight = 0.;
		double yearterm = (ftd->GetEndDateNA() - model->GetStartDate())/365.;

		double RecovCoef = 1.;

		//Test Recov coef
		if (ftd->GetCollateral()->GetRecovCoefFlg())
			RecovCoef = ftd->GetCollateral()->GetRecovCoef();

		for (int i=0; i<nbnames; i++)
		{
			double Pdefault = model->GetDefaultCurve(ftd->GetCollateral()->GetIssuersLabels(i))->DefaultProba(yearterm);
			sum+=Pdefault;
			Recov_weight+=Pdefault *MIN(model->GetRecoveryRate(ftd->GetCollateral()->GetIssuersLabels(i))*RecovCoef,1);
		}
		Recovery=Recov_weight/sum;

		// if (Pdefault) delete[] Pdefault;
		
		//Recovery = AvgRecovery(ftd->GetIssuersLabels(),nbnames);
	}

	switch (ftd->GetName())
	{
		case ICM_NTD:
		{
		ICM_Nthtd* ntd = (ICM_Nthtd*) ftd;
		tranche_down = ((double)(ntd->GetFirstNumDefault()-1)) * (1. - Recovery)* Notional;
		break;
		}
		case ICM_MEZ : 
		case ICM_CDO2 :
		{
		ICM_Mez* mez = (ICM_Mez*)ftd;
		tranche_down = mez->GetSubAmount(RefDate);
		break;
		}
	}

	return (tranche_down);
}	



ICM_LossUnits* ICM_Pricer_Distrib::CptLossUnits() 
{
	
	// 
	ICM_ModelMultiCurves* model = (ICM_ModelMultiCurves*) GetModel();
	ICM_Ftd* ftd = (ICM_Ftd*) GetSecurity();
	int k=0,nbnames = ftd->GetCollateral()->GetNbIssuers();


	ARM_Date date = model->GetStartDate() ;
	

	

	ICM_Vector<ARM_Date> notionalDates; 
	ftd->GetCollateral()->GetNotionalDates(notionalDates) ; 
	
	if (notionalDates.size()==0) notionalDates.push_back(date); 
	ICM_LossUnits* item = new ICM_LossUnits; 
	double initial_LU=0.;
	
	for(unsigned int iDates=0;iDates<notionalDates.size();iDates++) 
	{
		ARM_Date notionalDate = notionalDates[iDates]; 
		ICM_LossUnits::lossunit_data_t lossdata; 
		lossdata.LossRates.Resize(nbnames);
		lossdata.LSLossRates.Resize(nbnames);
		lossdata.CollatRank.resize(nbnames);
		lossdata.IsNull=false;

		//double loss_unit_min = ftd->GetCollateral()->SumNotionals(notionalDate)*itsMinLossUnit;
		double loss_unit_min = ftd->GetCollateral()->SumNotionals((ARM_Date)notionalDates[0])*itsMinLossUnit;
		double sum=0.;
		double sum_internal=0.;

		double RecovCoef = 1.;
		//Test Recov coef
		if (ftd->GetCollateral()->GetRecovCoefFlg())
			RecovCoef = ftd->GetCollateral()->GetRecovCoef();
		for (k=0; k<nbnames; k++)
		{	
			if (ftd->GetBinaryFlg())
				lossdata.LSLossRates[k] =ftd->GetCollateral()->GetIssuersNotional(k,notionalDate) * ( 1. - ftd->GetBinary());
			else if (ftd->GetCollateral()->GetRecovCoefFlg())
				lossdata.LSLossRates[k] =ftd->GetCollateral()->GetIssuersNotional(k,notionalDate) * ( 1. - MIN(model->GetRecoveryRate(ftd->GetCollateral()->GetIssuersLabels(k))*RecovCoef,1));
			else
				lossdata.LSLossRates[k] =ftd->GetCollateral()->GetIssuersNotional(k,notionalDate) * ( 1. - model->GetRecoveryRate(ftd->GetCollateral()->GetIssuersLabels(k)));

			if (!loss_unit_min) 
			{lossdata.LSLossRates[k] = round(lossdata.LSLossRates[k]/LOSS_UNIT_MIN)*LOSS_UNIT_MIN;}
			
			// case of notional too small.
			if((lossdata.LSLossRates[k]==0)&&(!loss_unit_min)){
				if(item) 
					delete item;
				if (!(ftd->GetCollateral()->GetIssuersNotional(k,notionalDate) > 100000)) {
					ICMTHROW(ERR_INVALID_DATA,"ERROR : Notional too small, set a notional >= 1000000"); 
				} else {
					ICMTHROW(ERR_INVALID_DATA,"ERROR : Notional too small for loss unit computation");
				}
				//lossdata.LSLossRates[k] = LOSS_UNIT_MIN;
			}
			sum_internal += lossdata.LSLossRates[k];

			//cas Notionels null
			//if CHECK_NULL(sum_internal)	{lossdata.LSLossRates[k]=1.;}
		}

		if (!loss_unit_min)
		{
			if (nbnames>1)
				lossdata.LU=fabs(pgcd::solve(lossdata.LSLossRates));
			else
				lossdata.LU=fabs(lossdata.LSLossRates[0]);
		}
		else if (!CHECK_NULL(sum_internal))
		{	lossdata.LU=fabs(loss_unit_min); }

		//cas Notionels null
		if (CHECK_NULL(sum_internal) && (iDates))
		{
			//lossdata.IsNull=true;
			lossdata.LU=0.;
			//lossdata.LU=initial_LU; 
		}

		// /!\ use the first date for tranche_down; tranche_up, but my not be good for "restrikable" product.
		double tranche_down = GetTranche_Down(date);
		double tranche_up = GetTranche_Up(date);

		if ((!loss_unit_min) &&  !CHECK_NULL(sum_internal))
		{
		if ((floor(tranche_up/lossdata.LU)>MAX_SIZE_LOSS_UNIT) || (floor(tranche_down/lossdata.LU)>MAX_SIZE_LOSS_UNIT)) 
			lossdata.LU = floor(tranche_up/((double)MAX_SIZE_LOSS_UNIT)); 
		}

		if (iDates==0)
		{initial_LU=lossdata.LU;}

		if (!CHECK_NULL(sum_internal)) {
		for (k=0; k<nbnames; k++)
		{
			lossdata.LSLossRates[k] /=lossdata.LU;
			lossdata.LossRates[k] = fabs(lossdata.LSLossRates[k]);
		}}

		//si nous sommes dans le cas NTD & FTD pas d'ordonnancement pour eviter un impact sur le prix
		if ((ftd->GetName()==ICM_NTD)||(ftd->GetName()==ICM_FTD)) 
		{
			for (k=0; k<nbnames; k++) lossdata.CollatRank[k]=k; 
			item->insert(notionalDate, lossdata); 
		}
		else 
		{
			ARM_Matrix matrix(nbnames,3);
			for (k=0; k<nbnames; k++)
			{	
				matrix.Elt(k,0)=lossdata.LossRates[k];
				matrix.Elt(k,1)=lossdata.LSLossRates[k];	
				matrix.Elt(k,2)=k;  
			}
			matrix.SortLines(0);

			for (k=0; k<nbnames; k++)
			{	
				lossdata.LossRates[k]=matrix.Elt(k,0);
				lossdata.LSLossRates[k]=matrix.Elt(k,1);	
				lossdata.CollatRank[k]=matrix.Elt(k,2) ; 
			}
			// ICM_LossUnits* item = new ICM_LossUnits; 
			item->insert(notionalDate, lossdata); 
		}
	}

	return item; 

}


// *************************************************************
// Computing of Avg Recovery
// *************************************************************
double ICM_Pricer_Distrib::AvgRecovery(const std::vector<std::string>&labels)
{
	ICM_ModelMultiCurves* model = (ICM_ModelMultiCurves*) GetModel();
	ICM_Ftd* ftd = (ICM_Ftd*) GetSecurity();

	if (ftd->GetBinaryFlg()) {return (ftd->GetBinary());}
	double Recovery = 0.,sum=0.;
	int size = labels.size(); 
	vector<double> Pdefault;Pdefault.resize(size);
	double Recov_weight = 0.;
	double yearterm = (ftd->GetEndDateNA() - model->GetStartDate())/365.;

	double RecovCoef = 1.;

	//Test Recov coef
	if (ftd->GetCollateral()->GetRecovCoefFlg())
		RecovCoef = ftd->GetCollateral()->GetRecovCoef();

	for (int k=0; k<size ; k++)
	{
		Pdefault[k] = model->GetDefaultCurve(labels[k])->DefaultProba(yearterm);
		sum+=Pdefault[k];
		Recov_weight+=Pdefault[k]*MIN(model->GetRecoveryRate(labels[k])*RecovCoef,1);		
	}
	Recovery=Recov_weight/sum;
	

	return (Recovery);

}	

// *************************************************************
// Computing of Expected Loss Tranche
// *************************************************************
double ICM_Pricer_Distrib::ExpectedLossTranche(const double& yearterm,vector<double>& losses)
{
	return cpt_elt_pricer_distrib(yearterm,*this,*GetModel(),*GetSecurity(),losses); 
}

//Estimation de l'espérance de l'amortissement de la tranche
//Uniquement estimée sur les tranches [x-100%]

double ICM_Pricer_Distrib::ExpectedAmortTranche(const double& yearterm)
{
	//Récupération des modèles
	ICM_ModelMultiCurves* model = (ICM_ModelMultiCurves*) GetModel();
	ICM_Ftd* ftd = (ICM_Ftd*) GetSecurity();
	
	ARM_Date date = (ARM_Date)(model->GetStartDate().GetJulian() + K_YEAR_LEN*yearterm);
	double tranche_up = GetTranche_Up(date);
	double tranche_down = GetTranche_Down(date);
	double Notional = fabs(ftd->GetCollateral()->SumNotionals(date));

	double ERT = 0., nomin=0., pdef=0., rec=0.;

	//Estimation uniquement sur les tranches de type [X% - 100%]
	if CHECK_NULL(tranche_up/Notional-1.)
	{
		//Calcul Expect Recovery
		for (int i=0; i<ftd->GetCollateral()->GetNbIssuers(); i++)
		{
			//Estimation de la proba de defaut
			pdef = model->GetDefaultCurve(ftd->GetCollateral()->GetIssuersLabels(i))->DefaultProba(yearterm);

			//Estimation du recovery
			if (ftd->GetBinaryFlg() == true)
				rec = ftd->GetBinary();
			else if (ftd->GetCollateral()->GetRecovCoefFlg())
				rec = (MIN(model->GetDefaultCurve(ftd->GetCollateral()->GetIssuersLabels(i))->GetRecovery()*ftd->GetCollateral()->GetRecovCoef(),1));
			else
				rec = model->GetDefaultCurve(ftd->GetCollateral()->GetIssuersLabels(i))->GetRecovery();

			//Estimation Expected Recovery
			ERT+= rec*pdef*ftd->GetCollateral()->GetIssuersNotional(i);
		}
		if (! CHECK_NULL(tranche_up-tranche_down))
		{
			ERT /=(tranche_up-tranche_down);
			ERT = MIN(1.,ERT);
		}
	}
	else
		ERT=0.;

	return 1.-ERT;
}

// *************************************************************
// Computing of Expected Loss Tranche
// *************************************************************
// 17783 double ICM_Pricer_Distrib::ExpectedLossTrancheForFastSpreadHedge(const double& yearterm)
// 17783 {
// 17783 	return cpt_elt_pricer_distrib_fast(yearterm,*this,*GetModel(),*GetSecurity()); 
// 17783 
// 17783 }

// *************************************************************
// Compute Fast Sensitivity Method For Baskets
// *************************************************************
double ICM_Pricer_Distrib::ComputeSensitivity(const qSENSITIVITY_TYPE& typesensi, 
											  const std::string& plot, 
											  const std::string&  label, 
											  double epsvalue, double epsilonGamma // useless
											  )
{
    double sensitivity =0.;
	int k =0;

	// keep track of Sensitivity
	//itsSensiType	=	typesensi;

	// 17783  if ((typesensi != ICM_FAST_SPREAD_TYPE) && (typesensi != ICM_FAST_RECOVERY_TYPE)&& (typesensi != ICM_SUBORDINATION))
	if ((typesensi != xICM_FAST_SPREAD_TYPE) && (typesensi != xICM_FAST_RECOVERY_TYPE)&& (typesensi != ICM_SUBORDINATION))
		return ICM_Pricer_Basket::ComputeSensitivity(typesensi,plot,label,epsvalue);

	ICM_ModelMultiCurves* DefModel = (ICM_ModelMultiCurves*) GetModel();

	ICM_Ftd* Basket = (ICM_Ftd*) GetSecurity();
	ARM_Date ExecutionDate = GetModel()->GetStartDate();
	int nbissuer = Basket->GetCollateral()->GetNbIssuers();

	int i = 0;
	double initialprice = 0., modifiedprice = 0.,result = 0.;

	// char**	itsNames=NULL;
	// int itsNbNames = 0;
 
    
    {

		// Get the initial Price if already computed
		if (GetInitialPriceFlg())
			initialprice = GetInitialPrice();
		else
			initialprice = Price(qCMPPRICE);

		// Activate Hedge Flag Computation
		// 17783 ActivateShift(true);

		switch (typesensi)
		{
			/** 
			case ICM_FAST_SPREAD_TYPE : //Fast Spread Sensitivity
			{
				
				// If already computed
				if (GetSensiManager())
				if (GetSensiManager()->GetSensiMatrix(ICM_FAST_SPREAD_TYPE))
					return 	(GetSensiManager()->GetSensiMatrix(ICM_FAST_SPREAD_TYPE))->GetValue(label,plot);

				const std::vector<std::string> & itsNames= Basket->GetCollateral()->GetIssuersLabels(); 
				// itsNbNames	= Basket->GetNbIssuers();
				// itsNames	= Basket->GetIssuersLabels();
				// if ((itsNames == NULL) || (itsNbNames == 0))
				// {
				// 	throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
			    //         "ERROR : not able to retrieve labels");
				// }
				if (itsNames.empty()) 
					ICMTHROW(ERR_INVALID_ARGUMENT,"ERROR : not able to retrieve labels"); 

				//Initialisation de la matrice "Hegdes"
				int sizetenors = 0;
				// true: in order to allow Parallel Shift
				char** tenors = DefModel->GetUnionTenorsForCollateral(itsNames,  sizetenors, true);
				
				ICM_QMatrix<ICM_DefaultCurve*>* ObjectMatrix = new ICM_QMatrix<ICM_DefaultCurve*>(itsNames.size(),sizetenors,NULL);
				ICM_Sensi2D* sensi2D = new ICM_Sensi2D(ICM_FAST_SPREAD_TYPE,itsNames,tenors,sizetenors,(void*)ObjectMatrix);

				GetSensiManager()->push_back(sensi2D);

//t1 = time(NULL);

				PerturbDefaultCurves();
//t2 = time(NULL);

//fprintf(fOut, "\nPerturbDefaultCurves:\t%.12lf\n", difftime(t2, t1));
								
				// Computes the SensiManager Matrix
				ComputePrice_Basket_Complete(ICM_FAST_SPREAD_TYPE, initialprice, epsvalue);

				result = (GetSensiManager()->GetSensiMatrix(ICM_FAST_SPREAD_TYPE))->GetValue(label,plot);
			
				break;
			}

			case ICM_FAST_RECOVERY_TYPE: // Fast Recovery Sensitivity
			{
				// plot is either "NONE", or "S" ???
				// If already computed
				if (GetSensiManager())
				if (GetSensiManager()->GetSensiMatrix(ICM_FAST_RECOVERY_TYPE))
					return 	(GetSensiManager()->GetSensiMatrix(ICM_FAST_RECOVERY_TYPE))->GetValue(label,plot);

				
				const std::vector<std::string> & itsNames= Basket->GetCollateral()->GetIssuersLabels(); 
				// itsNbNames	= Basket->GetNbIssuers();
				// itsNames	= Basket->GetIssuersLabels();
				// if ((itsNames == NULL) || (itsNbNames == 0))
				// {
				// 	throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
			    //         "ERROR : not able to retrieve labels");
				// }
				if (itsNames.empty()) 
					ICMTHROW(ERR_INVALID_ARGUMENT,"ERROR : not able to retrieve labels"); 


				//Initialisation de la matrice "Hegdes"
				int sizetenors = 1;
				// true: in order to allow Parallel Shift
				char** tenors = new char*[sizetenors];
				tenors[0]	=	new char[5] ; 
				strcpy(tenors[0],"NONE"); 
				
				ICM_QMatrix<ICM_DefaultCurve*>* ObjectMatrix = new ICM_QMatrix<ICM_DefaultCurve*>(itsNames.size(),sizetenors,NULL);
				ICM_Sensi2D* sensi2D = new ICM_Sensi2D(ICM_FAST_RECOVERY_TYPE,itsNames,tenors,sizetenors,(void*)ObjectMatrix);

				GetSensiManager()->push_back(sensi2D);
								
				// Computes the SensiManager Matrix
				ComputePrice_Basket_Complete(ICM_FAST_RECOVERY_TYPE, initialprice, epsvalue);

				result = (GetSensiManager()->GetSensiMatrix(ICM_FAST_RECOVERY_TYPE))->GetValue(label,plot);
			
				break;
			}
			**/ 
			case ICM_SUBORDINATION:
			{
				result = 0.0 ;

				ICM_Mez* NewMEZ = (ICM_Mez*) ((ICM_Mez*)GetSecurity())->Clone();

				ICM_ModelMultiCurves* ModelDef2 = (ICM_ModelMultiCurves*) DefModel ->Clone();
				
				((ICM_Smile_Correlation*) ModelDef2 ->GetCorrelation())->SetAlready_rescal(false);
				
				ICM_Pricer_Distrib* Pricer2 =(ICM_Pricer_Distrib*) CloneOnlyParams(NewMEZ,ModelDef2);
				
				Pricer2->Set(NewMEZ,ModelDef2,GetParameters(),GetAsOfDate());

				double  Dw_Bound;
				Dw_Bound = NewMEZ ->GetSubAmount(NewMEZ->GetFeeLeg()->GetStartDate());
				//Up_Bound = NewMEZ ->GetMezzAmount(NewMEZ->GetFeeLeg()->GetStartDate());
				
				double total_notionals = NewMEZ->GetCollateral()->SumNotionals(NewMEZ->GetFeeLeg()->GetStartDate());
				
				Dw_Bound = MAX(0,Dw_Bound+(epsvalue*total_notionals));

				//Up_Bound = MIN(total_notionals,Up_Bound+Dw_Bound);

				NewMEZ->SetSubAmount(Dw_Bound);
				//NewMEZ->SetMezzAmount(Up_Bound+Dw_Bound);
				

				modifiedprice =Pricer2->Price(qCMPPRICE);

				result = modifiedprice - initialprice ;

				if (ModelDef2)
					delete ModelDef2;
				if ( NewMEZ)
					delete NewMEZ;
				if( Pricer2)
					delete Pricer2;

			}
			break;
			default :
			result = -99999999.0;
		}

		// 17783 ActivateShift(false);

		ResetPricer();

	}
    

//	fclose(fOut);

	return (result);
}


// *************************************************************
// Generate ¨Perturbates Defaults Probabilities 
// *************************************************************
/**
17783 	
ICM_QMatrix<double>* ICM_Pricer_Distrib::PerturbProbasFromSensiManager(const double& yearterm)
{
	if (!GetSensiManager())
	{
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "ERROR : Pb. with PerturbProbasFromSensiManager");
	}

	int	i,j;
	double	PD;
	
	ICM_QMatrix<double>* ShiftMatrix = NULL;

	if (GetSensiManager()->GetSensiMatrix(ICM_FAST_SPREAD_TYPE))
	{
		ICM_Sensi2D* sensi2D = GetSensiManager()->GetSensiMatrix(ICM_FAST_SPREAD_TYPE);
		int its_sz_dt_i = sensi2D->GetMktDataSize();		// Names
		int its_sz_dt_j = sensi2D->GetBumpProfileSize();	// Tenor

		ShiftMatrix = new ICM_QMatrix<double>(its_sz_dt_j,its_sz_dt_i,0.);

		for (i=0;i<its_sz_dt_i;i++)
			for (j=0;j<its_sz_dt_j;j++)
			{
				PD = ((ICM_QMatrix<ICM_DefaultCurve*>*)sensi2D->GetObjectMatrix())->Getvalue(i,j)->DefaultProba(yearterm);
				ShiftMatrix->SetValue(j,i,PD);	// check the index
			}
	}
	else if (GetSensiManager()->GetSensiMatrix(ICM_FAST_RECOVERY_TYPE))
	{
		// not really optimal.
		ICM_Sensi2D* sensi2D = GetSensiManager()->GetSensiMatrix(ICM_FAST_RECOVERY_TYPE);
		int its_sz_dt_i = sensi2D->GetMktDataSize();		// Names
		int its_sz_dt_j = sensi2D->GetBumpProfileSize();	// Tenor

		ShiftMatrix = new ICM_QMatrix<double>(its_sz_dt_j,its_sz_dt_i,0.);

		for (i=0;i<its_sz_dt_i;i++)
			for (j=0;j<its_sz_dt_j;j++)
				ShiftMatrix->SetValue(j,i,1.0);	// check the index
	}

	return (ShiftMatrix);
}
17783 	
**/ 

void ICM_Pricer_Distrib::View(char* id, FILE* ficOut)
{	
	FILE* fOut;
	char  fOutName[200];
	ICM_Collateral* Collateral = ((ICM_Ftd*)GetSecurity())->GetCollateral();
	const std::vector<std::string>& Issuers = Collateral->GetIssuersLabels();

	if ( ficOut == NULL )
	{ARM_GetViewFile(ARM_VIEW_FILE, id, fOutName);
     (void) unlink(fOutName);
     fOut = fopen(fOutName, "w"); }
	else {fOut = ficOut;} 

	int size =0;

    fprintf(fOut, "\t\t\t ----------------- Homogeneous Security Pricer ----------------- \n");
// 17783 			fprintf( fOut,"TenorShift= %d, IssuerShift= %d\n",itsTenorShift,itsIssuerShift); 
	fprintf( fOut,"TargetSpreadForFlatCorrel=%f\n",itsTargetSpreadForFlatCorrel); 
// 17783 	fprintf( fOut,"IsActivateShift=%s\n",itsIsActivateShift?"TRUE":"FALSE"); 
// 17783 	fprintf( fOut,"FirstComputationForHedgeFlg=%s\n",itsFirstComputationForHedgeFlg?"TRUE":"FALSE"); 
	if (!itsLossUnits) 
	{
		fprintf( fOut,"LossUnit not defined\n"); 
	}
	else 
	{
		itsLossUnits->View(id, fOut); 
	}
	
	

//	fprintf(fOut, "\tLoss Unit : %f\n",its_lossunit);
	fprintf(fOut, "\tIssuers\t\tLoss Rate\n");

/**	for (int i=0; i<its_LossRate.size(); i++)
	{fprintf(fOut, "\t%s\t", Issuers[i].c_str());
	fprintf(fOut, "\t%f\n", its_LossRate[i]);}
	**/ 
	fprintf(fOut, "\n");
	
	ICM_Pricer_Security::View(id, fOut);

	if ( ficOut == NULL )
	{fclose(fOut);}
}


// -------------------------------------------------------------
// Generate ¨Perturbates Defaults Probabilities (new version)
// -------------------------------------------------------------
/**  17783 
void ICM_Pricer_Distrib::PerturbDefaultCurves()
{
	int nocurve = -1;
	int i =0,j=0;

	if (GetSensiManager()==NULL) return;

	const std::vector<std::string>& labels = (GetSensiManager()->GetSensiMatrix(ICM_FAST_SPREAD_TYPE))->GetMktData();
	char** tenors = (GetSensiManager()->GetSensiMatrix(ICM_FAST_SPREAD_TYPE))->GetBumpProfile();

	const ICM_DefaultCurve** VectorDefCurves = NULL;	
	ICM_DefaultCurve** ModifVectorDefCurves = NULL;	
	ICM_ModelMultiCurves * model = (ICM_ModelMultiCurves*)GetModel() ;

	int	its_sz_dt_i = (GetSensiManager()->GetSensiMatrix(ICM_FAST_SPREAD_TYPE))->GetMktDataSize();
	int its_sz_dt_j = (GetSensiManager()->GetSensiMatrix(ICM_FAST_SPREAD_TYPE))->GetBumpProfileSize();

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

// ----------------------------------------------------------------------------
//Compute Price + Hedges for Spreads
// ----------------------------------------------------------------------------
/**
17783 	
void ICM_Pricer_Distrib::ComputePrice_Basket_Complete(const ENUM_SENSITIVITY& sensitivity_type,
													  const double& initialprice, 
													  const double& epsvalue)
{
	// For Fast Spread Sensisitivity ICM_FAST_SPREAD_TYPE
	// get its_Cond

	// --------------------------------------------------------------
	// Default Proba
	ICM_Sensi2D* sensi2D = NULL;
	sensi2D = GetSensiManager()->GetSensiMatrix(sensitivity_type);

	if (sensi2D == NULL)
	{
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "ERROR : Pb. with Sensi Manager");
	}

	int	its_sz_dt_i = sensi2D->GetMktDataSize();
	int	its_sz_dt_j = sensi2D->GetBumpProfileSize();

	int i,j;

	ICM_ModelMultiCurves* DefModel = (ICM_ModelMultiCurves*) GetModel();
	ARM_Date ExecutionDate = GetModel()->GetStartDate();
	ICM_Ftd* ftd = (ICM_Ftd*) GetSecurity();

	double	modifiedprice;

	SetFirstComputationForHedgeFlg(true);
	// Clear vector
	ClearShifts();
	ResetPricer();
	for (j=0;j<its_sz_dt_j;j++)		// A Tenor (5Y, 10Y, Parallel...), Parallel has to be added
	{
		SetTenorShift(j);	// Tenor

		for (i=0;i<its_sz_dt_i;i++)	//	The given Name
		{
			if ((i) || (j))
				SetFirstComputationForHedgeFlg(false);

			SetIssuerShift(i);	// Label


			modifiedprice =	PreCheckBeforeComputePrice(sensitivity_type, ftd->GetCollateral()->GetIssuersLabels(i), initialprice, epsvalue);


			(GetSensiManager()->GetSensiMatrix(sensitivity_type))->GetMatrix()->SetValue(i, j, modifiedprice - initialprice);

			ResetPricerForHedge();
		}		
	}
}
17783 	
**/ 

/** 17783 	
double ICM_Pricer_Distrib::PreCheckBeforeComputePrice(const ENUM_SENSITIVITY& sensitivity_type, 
													  const std::string& label, 
													  const double& initialprice, 
													  const double& epsvalue)
{
	
	double	modifiedprice;

	ICM_ModelMultiCurves* DefModel = (ICM_ModelMultiCurves*) GetModel();
	ARM_Date ExecutionDate = GetModel()->GetStartDate();

	ARM_Security* Security = GetSecurity();;

	switch (sensitivity_type)
	{
		case ICM_FAST_RECOVERY_TYPE:

		if (Security->GetName() == ICM_MEZ)
		{
			modifiedprice = ComputePrice(qCMPPRICE);
		}
		else if (Security->GetName() == ICM_NTD)
		{
			modifiedprice = ComputePrice(qCMPPRICE);
		}

		break;
	
		case ICM_FAST_SPREAD_TYPE:

			modifiedprice = ComputePrice(qCMPPRICE);

			break;
	
		default:
			break;
	}

	return modifiedprice;
}
17783 	**/ 

// --------------------------------------------
// Compute Implied Curve for Tranche
// --------------------------------------------
ICM_DefaultCurve* ICM_Pricer_Distrib::GenerateImpliedCurve()
{
	ICM_ModelMultiCurves* model = (ICM_ModelMultiCurves*) GetModel();
	ICM_Ftd* basket = (ICM_Ftd*) GetSecurity() ;

	ARM_Vector* notionals = NULL ;
	ARM_Vector Tenors;
	int nbnames = basket->GetCollateral()->GetNbIssuers();
	const std::vector<std::string>& UnionIssuers=basket->GetCollateral()->GetIssuersLabels();
	model->GetUnionYFForCollateral(UnionIssuers,Tenors);
	int sizeTenors = Tenors.size();
	delete notionals;
	// FreePointerTabChar(UnionIssuers,nbnames);
	const ICM_DefaultCurve* RefCurve = model->GetDefaultCurves(0);

	int sizeTenors_Adj = 0;
	ARM_Date AsOf = RefCurve->GetAsOfDate();
	// ARM_Date Maturity = basket->GetMaturity();
	ARM_Date Maturity = basket->GetEndDateNA();
	double Matu = (Maturity-AsOf)/365.;
	for (int k=0;k<sizeTenors;k++) if (Tenors[k]>Matu) {sizeTenors_Adj=k+1;Maturity=(ARM_Date)(Tenors[k]*365.+ AsOf.GetJulian());break;}
	sizeTenors = sizeTenors_Adj;
	
	ARM_Vector spread(sizeTenors);
	ICM_Ftd* Basket_c = (ICM_Ftd*) basket->Clone();
	Basket_c->CptCashFlowDatesCredit(Maturity);
	SetSecurity(Basket_c);

	double Spread = ComputeSpread()/10000.;
	spread[sizeTenors-1]=Spread;
	double dRecovery = 0.;
	
	//For each tenor of the default curve we compute the breakeven spread
	for (int j=0;j<sizeTenors-1;j++)
	{
		ARM_Date End = Tenors[j]*365.+ AsOf.GetJulian();
		Basket_c->CptCashFlowDatesCredit(End);
		Spread = ComputeSpread()/10000.;
		spread[j]=Spread;
	}

	if (basket->GetName()==ICM_NTD) {dRecovery=AvgRecovery(UnionIssuers);}

	// char Tmp[50];
	// sprintf(Tmp,"ImpliedCurve");
	string Name("ImpliedCurve") ;

	SetSecurity(basket);
	if (Basket_c) delete Basket_c;

	//Now we can build the default curve for the underlying Cdo N°i
	//On impose un recovery = 0 pour chaque courbe sous-jacente
	return CptDefaultCurve(Tenors,spread,model,Name,dRecovery);
}

// --------------------------------------------------------------------------------------------------------------------
// Quick Generation of Default Curve
// --------------------------------------------------------------------------------------------------------------------
ICM_DefaultCurve* ICM_Pricer_Distrib::CptDefaultCurve(const ARM_Vector& dates,
													  const ARM_Vector& spreads,
													  ICM_ModelMultiCurves* model,
													  const string& label,
													  const double& recovery,
													  const string& RefCurve)
{
	int j=0,size=0;
	const ICM_DefaultCurve* curve = NULL;

	if (RefCurve == "NONE") curve = model->GetDefaultCurves(0);  	
	else curve = model->GetDefaultCurve(RefCurve);  	

	ARM_ZeroCurve* IRCurve = model->GetZeroCurve();
	ICM_DefaultCurve* DefCurve = curve->GenDefCurve(dates,spreads,recovery,label,false,IRCurve);

	return (DefCurve);
}

double ICM_Pricer_Distrib:: FlatCorrelEvaluate(double FlatCorrel)
{

	ICM_Mez* OriginalMez = (ICM_Mez*) GetSecurity();

	ICM_ModelMultiCurves* Originalmodel = (ICM_ModelMultiCurves*) GetModel();

	ICM_Mez* Mez = (ICM_Mez*) OriginalMez->Clone();
	
	ICM_ModelMultiCurves* model = (ICM_ModelMultiCurves*) Originalmodel->Clone();
	double currentspread = 0. ;

	ICM_Smile_Correlation* Correl ; 
	
	ARM_Date AsOfDate = model->GetDefaultCurves(0)->GetAsOfDate();
	ARM_Vector* SmileStrikeUp=new ARM_Vector(1,0.);
	ARM_Vector* SmileStrikeDw= new ARM_Vector(1,0.); ;

	// check function giving percentUp/down

	SmileStrikeUp->Elt(0) = Mez->GetPercentHight(Mez->GetFeeLeg()->GetStartDate()) +Mez ->GetPercentLow(Mez->GetFeeLeg()->GetStartDate()) ;
	SmileStrikeDw ->Elt(0)= Mez ->GetPercentLow(Mez->GetFeeLeg()->GetStartDate());

	// Creating Vol Obj with CorrelUp=CorrelDown=C0

	ARM_Vector* Matu = new ARM_Vector(1,0.);
	Matu->Elt(0)=5 ;

	ARM_Vector* Strikes = new ARM_Vector(2,0.);


	
	Strikes->Elt(0) = SmileStrikeDw->Elt(0) *100  ;
	Strikes->Elt(1) = SmileStrikeUp->Elt(0) *100 ;
	double * vols = new double[2];

	// start of algorithm : first Correl test

	
	vols[0]= FlatCorrel ; // TODO : voir point départ de l'algo
	vols[1] = FlatCorrel ;

	ARM_Matrix* mVol = new ARM_Matrix(1,2,vols) ;
	

	ARM_Currency currency((const char *) "EUR");
	
	// ObjVol
	ARM_VolLInterpol* newVolCrv = NULL;
	
	newVolCrv = new ARM_VolLInterpol( AsOfDate , Matu , Strikes , mVol , 1,K_ATMF_VOL,&currency);
	
	vector<const ARM_VolCurve*> VVolcurves;
	VVolcurves.resize(1);
	VVolcurves[0] = (ARM_VolCurve*) newVolCrv;

	// Obj Correl
	vector<std::string>  labels ;
	labels.resize(1);
	labels[0]="INDEX";

	ARM_Vector* proportions = new ARM_Vector(1,0.);
	proportions->Elt(0)=1;

	Correl = new ICM_Smile_Correlation( AsOfDate , "USERDEF",&(*VVolcurves.begin()) ,labels,proportions,SmileStrikeDw,SmileStrikeUp);
	
	model->SetCorrelation(Correl);


	ResetPricer();
	SetModel(model);
	computelossunit();

	currentspread = ComputeSpread();


	if (mVol)
		delete mVol;
	mVol=NULL;
	
	if (Correl)
		delete Correl;
	Correl=NULL;


	double result ;
	result = currentspread-itsTargetSpreadForFlatCorrel;
	
	return result;
}

double ICM_Pricer_Distrib::ComputeFlatCorrel()
{
	
	ICM_ModelMultiCurves* Originalmodel = (ICM_ModelMultiCurves*) GetModel();

	itsTargetSpreadForFlatCorrel = ComputeSpread();
	
	double _Up = 90 ; 
	double _Dw = 0.0;

	double epsilon = 1.E-5 ;
	double pas = 1.E-1 ;

	double result ;

	result = RootFinder1D(ff1::mem_call(&ICM_Pricer_Distrib::FlatCorrelEvaluate,(*this))).Dichotomy(_Dw,_Up,1000,epsilon,pas);

	ResetPricer();
	SetModel(Originalmodel);
	computelossunit();

	return result/100;


}

// virtual 
void ICM_Pricer_Distrib::Reset(void)
{
	ICM_Pricer_Basket::Reset();
	if (itsLossUnits) delete itsLossUnits; 
	itsLossUnits=0; 
}

ICM_Pricer_Distrib::~ICM_Pricer_Distrib(void)
{ 
	if (itsLossUnits) delete itsLossUnits; itsLossUnits=0; 
}
// 
void ICM_Pricer_Distrib::Init() 
{
	SetName(ICM_PRICERHOMOGENEOUS);

	itsLossUnits=0; 
	itsMinLossUnit=0.;
// 17783 			itsTenorShift = -1;
// 17783 			itsIssuerShift = -1;
	itsTargetSpreadForFlatCorrel = -1 ;
	// 17783 itsIsActivateShift = false;
	// 17783 itsFirstComputationForHedgeFlg	=	false;
// 17783 	LossTShifts.clear();
}

void ICM_Pricer_Distrib::ResetLossUnit()
{ 
	if (!itsLossUnits) delete itsLossUnits;
		itsLossUnits= NULL;
}
