#include "firsttoinc.h"

#include "ICMKernel\pricer\ICM_Pricer_Basket.h"
#include "ICMKernel\inst\icm_cdo2.h"
#include "ICMKernel\inst\icm_nthtd.h"
// #include "ICMKernel\mod\icm_Gauss1FSmile.h"
// #include "ICMKernel\mod\icm_MixCopulaCalibration.h"
#include "ICMKernel\glob\icm_smile_correlation.h"
#include "ICMKernel\glob\icm_betas_correlation.h"
#include "ICMKernel\inst\icm_collateral.h"
#include "ICMKernel\inst\icm_credit_index.h"
#include "ICMKernel\mod\modelmulticurves.h"
#include "ARMKernel\crv\volflat.h"

ICM_Pricer_Basket::~ICM_Pricer_Basket()
{
	if (itsBetas)
		delete itsBetas;
	itsBetas = NULL;
}
void ICM_Pricer_Basket::Init()
{
	SetName(ICM_PRICERSECURITY);

	itsBetas = NULL;
	itsEqCorrelUp =0.;
	itsEqCorrelUpFlg = false;
	itsEqCorrelDown=0.;
	itsEqCorrelDownFlg = false;
}
void ICM_Pricer_Basket::Set(ARM_Security *sec, ARM_Model *mod,const ICM_Parameters&params,const ARM_Date&asof)
{
	if ((mod->GetName()!=ICM_MODELMULTICURVES))			
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "ERROR : pricer not compatible with this security or this model ");

	ICM_Pricer_Security::Set(sec, mod,params,asof);
}
double ICM_Pricer_Basket::GetBeta(const std::string& issuer)
{
	if (!itsBetas)
		ComputeBeta() ;
	
	ICM_Ftd* ftd = (ICM_Ftd *) GetSecurity() ;	
	int index = ftd->GetCollateral()->getIssuerPosition(issuer); 
	return itsBetas->Elt(index); 

	/** ICM_ModelMultiCurves * model = (ICM_ModelMultiCurves*) GetModel() ;
	unsigned int size    = model->GetNbDefCurves() ;
	// _DefaultCurve** def = model->GetDefaultCurves();
	int indx = -1 ;

	for(unsigned int i=0; i<size;i++)
	{
		if (!strcmp(model->GetDefaultCurves(i)->GetLabel().c_str(),issuer))
		{
			indx = i ; 
			break ;
		}	
	}

	if (indx<0)	return	999.;
	return itsBetas->Elt(indx) ; 
	**/ 
} 

// virtual
void ICM_Pricer_Basket::Reset(void)
{
	ICM_Pricer_Security::Reset();
	
	if (itsBetas)
		delete itsBetas;
	itsBetas = NULL;

	ICM_Cds* cds = (ICM_Cds*) GetSecurity(); //For CM-CDS sensitivity
	if (cds->GetFeeLeg()->GetCreditIndex()) cds->GetFeeLeg()->GetCreditIndex()->ResetDefCurve();
}
// *************************************************************
// Beta Computation 
// *************************************************************
void ICM_Pricer_Basket::ComputeBeta()
{
	ICM_Ftd* ftd = (ICM_Ftd *) GetSecurity() ;	
	ICM_ModelMultiCurves* model = (ICM_ModelMultiCurves *) GetModel() ;
	ICM_Collateral*collat=ftd->GetCollateral(); 	

	int size = collat->GetNbIssuers();
	const std::vector<std::string> & labels = ftd->GetCollateral()->GetIssuersLabels();
	ARM_Date Maturity = ftd->GetEndDateNA();

	ICM_Correlation* Correlation = model->GetCorrelation();
	ARM_Vector notionals; 
	collat->GetConstantNotionals(notionals); 
	ARM_Vector* Betas = Correlation->ComputeBetas(size,labels,notionals,Maturity,model);

	if (Betas) SetBetasPointer(Betas);

	return;
}

// *************************************************************
// Compute Sensitivity Method For Baskets
// *************************************************************
double ICM_Pricer_Basket::ComputeSensitivity(const qSENSITIVITY_TYPE& typesensi, 
											 const std::string& plot, 
											 const std::string& label, 
											 double epsvalue , double epsilonGamma // useless
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
		return 0; 

		/** if (!strcmp(label,"UP")) return CptCorrSensi(epsvalue, 1);
		else
		 return CptCorrSensi(epsvalue, 0);		
		 **/ 
	}

    
    {

		if (typesensi != ICMBETA_WITH_SPREAD_SHIFT)
		{
		if (GetInitialPriceFlg())
			initialprice = GetInitialPrice();
		else
			initialprice = Price(qCMPPRICE);

		ResetPricer();
		if(typesensi == ICMRECOVERY_TYPE) ResetLU();
		}

		switch (typesensi)
		{
			case ICMRECOVERY_TYPE :
			case ICMIRCURVE_TYPE :
			case ICMIRCURVE_WITHOUT_DEFCURVE_TYPE :
			case ICMSPREAD_TYPE :
			case ICMSPRELSHIFT_TYPE :
			case ICM_DTR_TYPE :
			case ICM_SAMEBETA_TYPE :
			case ICM_SAMECORRELATION_TYPE :
			case ICMCORRELATION_TYPE :
			case ICMBETA_TYPE :
			case ICMCORREL_STRIKE_DOWN_TYPE : 
			case ICMCORREL_STRIKE_UP_TYPE : 
			case ICM_CORREL_BET_CDO_UP_TYPE :
			case ICM_CORREL_BET_CDO_DW_TYPE :
			case ICM_INFLATION_CPN_CURVE :
			case ICM_INTEREST_CPN_CURVE :
			case ICM_IRCURVE_WITH_CPN:
			case ICM_INDX_SPREAD_RESCALING :
 			{
	
				ICM_ModelMultiCurves* ModelDef2 = DefModel->GenerateShiftModel(typesensi,
																			   plot, 
																			   label,
																			   nocurve,
																			   epsvalue);
				
				SetModel(ModelDef2);	
				SetBetasPointer(NULL);
				computelossunit();
				BeforePrice(label,typesensi);
				modifiedprice = ComputePrice(qCMPPRICE);
#ifdef _DEBUG
				FILE* pFile = NULL;
				pFile = fopen("C:\\temp\\testPricerShiftedRecov.txt","w");
				View("",pFile);
				if ( pFile) fclose(pFile);
#endif
				SetModel(DefModel); //On reset le model initial
				SetBetasPointer(NULL);
				computelossunit();

				if (ModelDef2)
					delete ModelDef2;
				ModelDef2 = NULL;
			}
			break;
			case ICM_GREEK_VEGA_TYPE : 
			{
			ICM_Cds* cds = (ICM_Cds*) GetSecurity(); //For CM-CDS sensitivity
			ICM_Credit_Index* index = dynamic_cast<ICM_Credit_Index*>(cds->GetFeeLeg()->GetCreditIndex());
			if( !index) ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_Pricer_Basket::ComputeSensitivity : case ICM_GREEK_VEGA_TYPE : index is not ICM_Credit_Index");  	
			string VolName = index->GetIndexName();
			ARM_VolCurve* volCurve = dynamic_cast<ARM_VolCurve*>(DefModel->GetVolCurve());  // only One !!
			if( !volCurve) ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_Pricer_Basket::ComputeSensitivity : case ICM_GREEK_VEGA_TYPE : No vol curve in DefModel");  	
			
			std::auto_ptr<ICM_ModelMultiCurves> DefModelShifted(dyn_clone(DefModel));
			if (volCurve->GetName() == ARM_VOL_FLAT)
			{
				ARM_VolFlat* volFlat_cloned = (ARM_VolFlat*) dyn_clone(volCurve);
				volFlat_cloned->BumpVolatility(epsvalue*100.);
				DefModelShifted->SetVolCurve(volFlat_cloned);
				if (volFlat_cloned && DefModelShifted->GetFlgClone()) delete volFlat_cloned;
			} else
			{
				ARM_VolLInterpol* volInterp_cloned = (ARM_VolLInterpol*) dyn_clone(volCurve); 
				int nthLine,nthCol;
				SearchCoordForVolatility(volInterp_cloned,label,nthLine,nthCol); 
				if ((nthLine == -1) || (nthCol == -1)) break;
				volInterp_cloned->BumpVolatility(epsvalue*100.,nthLine+1,nthCol+1);
				DefModelShifted->SetVolCurve(volInterp_cloned);
				if(volInterp_cloned && DefModelShifted->GetFlgClone()) delete volInterp_cloned;
			}
				SetModel(DefModelShifted.get());	
				SetBetasPointer(NULL);
				computelossunit();
				BeforePrice(label,typesensi);
				modifiedprice = ComputePrice(qCMPPRICE);
#ifdef _DEBUG
				FILE* pFile = NULL;
				pFile = fopen("C:\\temp\\testPricerShiftedVol.txt","w");
				View("",pFile);
				if ( pFile) fclose(pFile);
#endif
				SetModel(DefModel); //On reset le model initial
				SetBetasPointer(NULL);
				computelossunit();
			}
			break;
			case ICMBETA_WITH_SPREAD_SHIFT :
			{

				ICM_ModelMultiCurves* ModelDef2 = DefModel->GenerateShiftModel(typesensi,
																	 "NONE", 
																	 "NONE",
																	 nocurve,
																	 epsvalue);

				ICM_Beta_Correlation* correl_init = (ICM_Beta_Correlation*) DefModel->GetCorrelation();

				double EPSL_BETA = 0.01;

				

				int nbcurves_corr = correl_init->GetSize(); 	

				double* NewBetas_plus = new double[nbcurves_corr];
				double* NewBetas_moins = new double[nbcurves_corr];

				ARM_Vector* InitialBetas = NULL;

				memcpy(NewBetas_plus,correl_init->GetBetas()->GetElt(),sizeof(double)*nbcurves_corr);
				memcpy(NewBetas_moins,correl_init->GetBetas()->GetElt(),sizeof(double)*nbcurves_corr);

				// if (strcmp(label,"NONE") == NULL)
				if (label=="NONE") 
				{
					for (i=0; i<nbcurves_corr; i++)
					{
						NewBetas_plus[i] = EPSL_BETA + NewBetas_plus[i];
						NewBetas_moins[i] = -EPSL_BETA + NewBetas_moins[i];
						if (NewBetas_plus[i]>=1.) NewBetas_plus[i] = 0.999;
						if (NewBetas_moins[i]<=-1.) NewBetas_moins[i] = -0.999;
					}
				}
				else
				{
						DefModel->CptNoCurve(label.c_str(),nocurve);
						NewBetas_plus[nocurve] = EPSL_BETA + NewBetas_plus[nocurve];
						NewBetas_moins[nocurve] = -EPSL_BETA + NewBetas_moins[nocurve];
						if (NewBetas_plus[nocurve]>=1.) NewBetas_plus[nocurve] = 0.999;
						if (NewBetas_moins[nocurve]<=-1.) NewBetas_moins[nocurve] = -0.999;
				}
	
				ARM_Vector V_Betas_plus (nbcurves_corr,NewBetas_plus);
				ARM_Vector V_Betas_moins (nbcurves_corr,NewBetas_moins);

				ICM_Beta_Correlation* corr_moins = new ICM_Beta_Correlation(ExecutionDate,"SENSI_CORR", V_Betas_moins, correl_init->GetLabels(),(ARM_IRIndex*)0,(ARM_IRIndex*)0);
				ICM_Beta_Correlation* corr_plus = new ICM_Beta_Correlation(ExecutionDate,"SENSI_CORR", V_Betas_plus, correl_init->GetLabels(),(ARM_IRIndex*)0,(ARM_IRIndex*)0);

				ResetPricer();

				ModelDef2->SetCorrelation(corr_moins);
		
				SetBetasPointer(NULL);
				SetModel(ModelDef2);
				
				BeforePrice(label,typesensi);

				initialprice = ComputePrice(qCMPPRICE);
				SetBetasPointer(NULL);

				ResetPricer();

				ModelDef2->SetCorrelation(corr_plus);
		
				SetBetasPointer(NULL);
				SetModel(ModelDef2);	

				BeforePrice(label,typesensi);

				modifiedprice = ComputePrice(qCMPPRICE);
				SetBetasPointer(NULL);

				SetModel(DefModel); //On reset le model initial

				if (NewBetas_plus)
					delete [] NewBetas_plus ;
				NewBetas_plus = NULL;

				if (NewBetas_moins)
					delete [] NewBetas_moins;
				NewBetas_moins = NULL;

				if (ModelDef2)
					delete ModelDef2;
				ModelDef2 = NULL;

				// if (V_Betas_plus)
				// 	delete V_Betas_plus;
				// V_Betas_plus = NULL;

				// if (V_Betas_moins)
				// 	delete V_Betas_moins;
				// V_Betas_moins = NULL;

				initialprice /= 2.;
				modifiedprice /= 2.;
			}
			break;
			case ICM_ISSUER_DEFAULT :
			{

				// if(!strcmp(plot,"S")) initialprice = -initialprice ;
				if (plot=="S") initialprice = -initialprice ;

				ARM_Security* Security = GetSecurity();

				if (Security->GetName() == ICM_MEZ)
				{

				ICM_Mez* InitMEZ = (ICM_Mez*)Security; 

				double A = InitMEZ->GetSubAmount(InitMEZ->GetFeeLeg()->GetStartDate());
				double B = InitMEZ->GetMezzAmount(InitMEZ->GetFeeLeg()->GetStartDate()) + A;
				int issuerRank =  InitMEZ->GetCollateral()->getIssuerPosition(label);
				double NOT= InitMEZ->GetCollateral()->GetIssuersNotional(issuerRank,InitMEZ->GetFeeLeg()->GetStartDate()); 
				// double NOT = InitMEZ->GetCollateral()->GetIssuersNotional(issuerRank);
				double LR = 1. - epsvalue;
				double alpha = NOT*LR;
				double FULLNOT = InitMEZ->GetCollateral()->SumNotionals(InitMEZ->GetStartDateNA());

				ICM_Mez* NewMEZ = (ICM_Mez*) ((ICM_Mez*)GetSecurity())->Clone();
				NewMEZ->ExcludeIssuer(label.c_str());

				if(alpha > B)
				{
					Actual_Loss = B-A;
					modifiedprice = 0.;
				}
				else
				{	if (A<alpha && alpha<=B)
				
					{
						NewMEZ->SetSubAmount(0.);
						NewMEZ->SetMezzAmount(B - alpha);
						Actual_Loss = alpha - A;

					}
					else
					{
						NewMEZ->SetSubAmount(A - alpha);
						NewMEZ->SetMezzAmount(B - A );
					}
					
					SetSecurity(NewMEZ);
					SetBetasPointer(NULL);
					SetModel(DefModel);	
					computelossunit();

					modifiedprice = ComputePrice(qCMPPRICE);
					
					// if(!strcmp(plot,"S")) modifiedprice = -modifiedprice;
					if (plot=="S") modifiedprice = -modifiedprice;

				}
				
				SetSecurity(InitMEZ);
				SetModel(DefModel); //On reset le model initial
				SetBetasPointer(NULL);
				computelossunit();

				if (NewMEZ) delete NewMEZ; 

				if (!result)
				{
					//Ajouter le CF alpha si la tranche est touchée
					// if(!strcmp(plot,"S")) result = modifiedprice - initialprice + Actual_Loss;
					if (plot=="S") result = modifiedprice - initialprice + Actual_Loss;
					else result=modifiedprice - initialprice - Actual_Loss; //On est Long par défaut
				}

				}
				else if (Security->GetName() == ICM_NTD)
				{
					ICM_Nthtd* InitNTD = (ICM_Nthtd*)Security; 
					int issuerRank = InitNTD->GetCollateral()->getIssuerPosition(label) ;
					double NOT = InitNTD->GetCollateral()->GetIssuersNotional(issuerRank,InitNTD->GetFeeLeg()->GetStartDate());
					double LR = 1. - epsvalue; //LR to be paid in case of default

					if (InitNTD->GetFirstNumDefault() == 1)
					{
					modifiedprice = 0.;
					Actual_Loss = NOT * LR;

					//Ajouter le CF alpha si la tranche est touchée
					// if(!strcmp(plot,"S")) modifiedprice += Actual_Loss;
					if (plot=="S") modifiedprice += Actual_Loss;
					else modifiedprice -= Actual_Loss;; //On est Long par défaut

					break;
					}

					ICM_Nthtd* NewNTD = (ICM_Nthtd*) ((ICM_Nthtd*)GetSecurity())->Clone();

					NewNTD->ExcludeIssuer(label.c_str());
					NewNTD->SetFirstNumDefault(MAX(NewNTD->GetFirstNumDefault()-1,0.));
					NewNTD->SetLastNumDefault(MAX(NewNTD->GetLastNumDefault()-1,0.));

					SetSecurity(NewNTD);
					SetModel(DefModel);	
					SetBetasPointer(NULL);

					modifiedprice = ComputePrice(qCMPPRICE);

					// if(!strcmp(plot,"S")) modifiedprice = -modifiedprice;
					if(plot=="S") modifiedprice = -modifiedprice;

					SetSecurity(InitNTD);
					SetModel(DefModel); //On reset le model initial
					SetBetasPointer(NULL);

					if (NewNTD) delete NewNTD; 
				}
				else result = -99999999.0;
			}
			break;
			default :
			result = -99999999.0;
		}

	ResetPricer();
	if(typesensi == ICMRECOVERY_TYPE) ResetLU();


	if (!result) result=modifiedprice - initialprice;

	}
     

	return (result);
}


// *************************************************************
// Compute Correlation Smile Method
// *************************************************************
/** 
double ICM_Pricer_Basket::ComputeCorrelSmile(int SmileType, double MktData, double seed, double UpfrontPay, int DataType)
{

	Gauss1FSmile* MultiDimMin = new Gauss1FSmile(this,
												 DataType, 
												 MktData,
												 UpfrontPay,
												 seed,
												 SmileType);
	
	double Result = MultiDimMin->GetitsSmile();

	return Result;
}
**/ 


// *****************************************************************
// Fonctions utilisées dans le calcule de l'At The Money Correlation 
// *****************************************************************

/**
double ICM_Pricer_Basket::CpteDL(double StrikeUp, double BC, double SmileType)
{ return 0.;}

// Compute the Equivalent Strikes for the index given the tranche to price and its Expected Loss
void ICM_Pricer_Basket::CpteStrikesEquivalents(double StrikeDown,
											   double StrikeUp,
											   double Index_Notional,
											   double EL_Tranche,
											   vector<double>& Res)
{
	double EL_Index = (1.0/Index_Notional) * (ICM_Pricer_Basket::CpteDL(Index_Notional,0.50,0)) ;
	double X1 = (StrikeDown/EL_Tranche) - 1 ;
	double X2 = (StrikeUp/EL_Tranche) - 1 ;

	Res.resize(3) ;
	Res[0] = (1+X1)*EL_Index; //StrikeDown Eq 
	Res[1] = (1+X2)*EL_Index; //StrikeUp Eq
	Res[2] = EL_Index ; //EL_Index
}


double ICM_Pricer_Basket::CpteDefLegTrancheEquivalente(double StrikeUp1,
													   double StrikeUp2,
													   double BC1,
													   double BC2,
													   double SmileType)
{
	double DefLeg = 0.;

	if (StrikeUp1 > 0)
		DefLeg = ICM_Pricer_Basket::CpteDL(StrikeUp2,BC2,SmileType) -
				 ICM_Pricer_Basket::CpteDL(StrikeUp1,BC1,SmileType) ;
	else
		DefLeg = ICM_Pricer_Basket::CpteDL(StrikeUp2,BC2,SmileType) ;
	
	return DefLeg ;
}
**/ 
void ICM_Pricer_Basket::ComputeEqCorrelUp(){
	/*	
	double lower=0.;
	double upper=0.;
	pSec->Bounds(lower, upper); */
	double asOfDate =0.;
	double Matu=0.;
	double MatuYF =0.;
	ICM_ModelMultiCurves* pModel = NULL;
	pModel = dynamic_cast<ICM_ModelMultiCurves*>(GetModel());
	if (!pModel){
		ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_Pricer_Basket::ComputeStrikeUp : pModel is not a ICM_ModelMultiCurves");  
	}
	ICM_Smile_Correlation* pCorrelation = NULL;
	ICM_Mez *pSec = dynamic_cast<ICM_Mez*>(GetSecurity());
	if (!pSec){
		ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_Pricer_Basket::ComputeCorrelUp : Security is not a Mezzanine");  
	}
	Price(qCMPPRICE);
	pCorrelation = dynamic_cast<ICM_Smile_Correlation*>(pModel->GetCorrelation()->Clone());
	if (!pCorrelation){
		ICMTHROW( ERR_INVALID_ARGUMENT,
                   "ICM_Pricer_Basket::ComputeEqCorrelUp : Model->GetCorrelation() in not an ICM_Smile_Correlation");
	}
	
	ARM_Date endDate = pSec->GetEndDateNA();
	ARM_Date startDate = pModel->GetStartDate();
	double yf_maturity = (endDate-startDate)/365;
	itsEqCorrelUp = pCorrelation->GetCorrelStrikeUp(yf_maturity);
	itsEqCorrelUpFlg = true;

	if (pCorrelation)
	delete pCorrelation;
	pCorrelation = NULL;
}
void ICM_Pricer_Basket::ComputeEqCorrelDown(){
	double asOfDate =0.;
	double Matu=0.;
	double MatuYF =0.;
	ICM_ModelMultiCurves* pModel = NULL;
	pModel = dynamic_cast<ICM_ModelMultiCurves*>(GetModel());
	if (!pModel){
		ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_Pricer_Basket::ComputeStrikeUp : pModel is not a ICM_ModelMultiCurves");  
	}
	ICM_Smile_Correlation* pCorrelation = NULL;
	ICM_Mez *pSec = dynamic_cast<ICM_Mez*>(GetSecurity());
	if (!pSec){
		ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_Pricer_Basket::ComputeEqCorrelDown : Security is not a Mezzanine");  
	}
	Price(qCMPPRICE);
	pCorrelation = dynamic_cast<ICM_Smile_Correlation*>(pModel->GetCorrelation()->Clone());
	if (!pCorrelation){
		ICMTHROW( ERR_INVALID_ARGUMENT,
                   "ICM_Pricer_Basket::ComputeStrikeDown : Model->GetCorrelation() in not an ICM_Smile_Correlation");
	}
	
	ARM_Date endDate = pSec->GetEndDateNA();
	ARM_Date startDate = pModel->GetStartDate();
	double yf_maturity = (endDate-startDate)/365;
	itsEqCorrelDown = pCorrelation->GetCorrelStrikeDown(yf_maturity);
	itsEqCorrelDownFlg = true;

	if (pCorrelation)
	delete pCorrelation;
	pCorrelation = NULL;

}

// ****************************************************************
// Fonction utilisées pour la calibration du modèle de Mixed Copula
// ****************************************************************
/**
double ICM_Pricer_Basket::Get_MixCopula_Factor (int SmileType, 
												double BC1, 
												double BC2, 
												double SeedIndep, 
												double SeedFull,
												double Accuracy,
												int FactorType)
{
	ICM_Cdo2* cdo2 = (ICM_Cdo2*) GetSecurity();
	
	ICM_Portfolio* portfolio = (ICM_Portfolio*) cdo2->GetPortfolio();
	
	MixCopulaCalibration* MultiDimMin = new MixCopulaCalibration(this,
																 portfolio,
																 SmileType, 
																 BC1, 
																 BC2, 																				
																 SeedIndep, 
																 SeedFull);
	
	double Result = MultiDimMin->Get_MixCopula_Factor(FactorType);

	return Result;
}
**/
double ICM_Pricer_Basket::DoPrice(qCMPMETH measure){
	double res =0.;
	ARM_Date ExecutionDate = GetModel()->GetStartDate();
	switch(measure) {
	case qCMPCORRELDOWN:
		if (! itsEqCorrelDownFlg) ComputeEqCorrelDown();
		res = getEqCorrelDown(); 
		break;
	case qCMPCORRELUP:
		if (! itsEqCorrelUpFlg) ComputeEqCorrelUp();
		res = getEqCorrelUp();
		break;
	case qCMPFLATCORR:
		res = ComputeFlatCorrel();
		break;

	default :
		string pmeas;
		ICM_EnumsCnv::toString(measure, pmeas);
		ICMTHROW( ERR_INVALID_ARGUMENT,
			"Pricer Basket : unable to Compute Pricing Measure :"<< pmeas);	
	}
	return res;
}

/**
double ICM_Pricer_Basket::Get_ReducedMixCopula_Factor(int SmileType,
													  double BC1,
													  double BC2,
													  double SeedIndep,
													  double SeedBeta,
													  int FactorType)
{
	ICM_Cdo2* cdo2 = (ICM_Cdo2*) GetSecurity();
	
	ICM_Portfolio* portfolio = (ICM_Portfolio*) cdo2->GetPortfolio();
	
	MixCopulaCalibration* MultiDimMin = new MixCopulaCalibration(this,
																 portfolio,																 
																 SmileType, 
																 BC1, 
																 BC2, 	
																 SeedIndep, 
																 SeedBeta);
	
	double Result = MultiDimMin->Get_ReducedMixCopula_Factor(FactorType);

	return Result;
}
**/

//********************
// Méthode de View
//********************

void ICM_Pricer_Basket::View(char* id, FILE* ficOut)
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

	int size =0;

	if (itsBetas) size = itsBetas->GetSize();

    fprintf(fOut, "-----------------------------------------------------------\n");
    fprintf(fOut, "---                      Basket Pricer                  ---\n");
	fprintf(fOut, "-----------------------------------------------------------\n");
	if (itsEqCorrelUpFlg) fprintf(fOut, "\tCorrelUp: %.2lf\t\t\t%s\n",itsEqCorrelUp,"Ok");
	if (itsEqCorrelDownFlg) fprintf(fOut, "\tCorrelDown: %.2lf\t\t\t%s\n",itsEqCorrelDown,"Ok");

	fprintf(fOut, "\n\tBeta\t\tImpliedCorr\n\n");
	for (int i=0; i<size; i++)
	{
	fprintf(fOut, "\t%f", itsBetas->Elt(i));
//	if (itsImpliedCorrelations)	
//	{
//		if (i<itsImpliedCorrelations->GetSize())
//			fprintf(fOut, "\t\t%f", itsImpliedCorrelations->Elt(i));
//	}
	fprintf(fOut, "\n");
	}

	fprintf(fOut, "\n");

//	if (itsBaseCorrelations) itsBaseCorrelations->View(id, fOut);
	
	ICM_Pricer_Security::View(id, fOut);

	if ( ficOut == NULL )
	{
		fclose(fOut);
	}
}
