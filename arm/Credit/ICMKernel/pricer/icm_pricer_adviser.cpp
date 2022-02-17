
#include "ARMKernel\glob\firsttoinc.h"

#include "ICMKernel/pricer/icm_pricer_adviser.h"

#include "ICMKernel/pricer/icm_pricer_analytic_cdo2.h"
#include "ICMKernel/pricer/icm_pricer_analytic_cdo2_smile.h"
#include "ICMKernel/pricer/icm_pricer_mc_cdo2.h"
#include "ICMKernel/pricer/icm_pricer_indexoption.h"
#include "ICMKernel/pricer/icm_pricer_capfloor_cmcds.h"
#include "ICMKernel/pricer/icm_pricer_legs_basket.h"
#include "ICMKernel/pricer/icm_pricer_tree_binomial_cds.h"
#include "ICMKernel/pricer/icm_pricer_cppi.h"
#include "ICMKernel/pricer/icm_pricer_tree_binomial_spreadoption.h"
#include "ICMKernel/pricer/icm_pricer_corridor.h"
// #include "ICMKernel/pricer/icm_pricer_indexoption_smile.h"
// #include "ICMKernel/pricer/icm_pricer_indexoption_sabr.h"
#include "ICMKernel/pricer/icm_pricer_ir_ycmod.h"
#include "ICMKernel/pricer/icm_pricer_ir_infcurve.h"
#include "ICMKernel/pricer/icm_pricer_homogeneous_smile_rf.h"
#include "ICMKernel/pricer/icm_VanillaFwdPricer.h"
#include "ICMKernel/pricer/icm_pricer_homogeneous_smile_collat_fwd.h"
// #include "ICMKernel/pricer/icm_pricer_index_capgen_smile.h"
// #include "ICMKernel/pricer/icm_pricer_index_capgen_blackscholes_strikes.h"
#include "ICMKernel/pricer/icm_pricer_cln.h"
#include "ICMKernel/pricer/icm_pricer_cdo_implied_loss.h"
#include "ICMKernel/pricer/icm_pricer_lss_gap_option.h"
#include "ICMKernel/pricer/icm_pricer_cpdo.h"
#include "ICMKernel\cair\credit_manager.h"
#include "ICMKernel\pricer\icm_basepricer.h"
#include "ICMKernel/pricer/icm_pricer_mc_cdo.h"
#include "ICMKernel/pricer/icm_pricer_tree_hw_cdo.h"
#include "ICMKernel/pricer/icm_pricer_tranche_restrikable.h"
#include "ICMKernel/inst/icm_credit_index.h"
#include "ICMKernel/mod/modelmulticurves.h"

//	----------------------------------------------------------------------------
void 
ICM_Pricer_Advisor::Init() 
{
	ICM_PRICER_TYPE type;
	type.Set(ICM_SECURITY);itsMapAdvisor.insert(type);
	type.Set(ICM_GENCF);itsMapAdvisor.insert(type);
	type.Set(ICM_LEG);itsMapAdvisor.insert(type);
	type.Set(ICM_BOND);itsMapAdvisor.insert(type);
	type.Set(ICM_FRN);itsMapAdvisor.insert(type);
	type.Set(ICM_ZCCLN);itsMapAdvisor.insert(type);
	type.Set(ICM_PF);itsMapAdvisor.insert(type);
	type.Set(ICM_CDS,ICM_PRICER_CDS);itsMapAdvisor.insert(type);
	type.Set(ICM_CDS,ICM_BINOMIAL_TREE_PRICER_CDS);itsMapAdvisor.insert(type);
	type.Set(ICM_FTD,ICM_PRICERHOMOGENEOUS);itsMapAdvisor.insert(type);
	type.Set(ICM_FTD,ICM_PRICER_HOMOGENEOUS_SMILE);itsMapAdvisor.insert(type);
	type.Set(ICM_NTD,ICM_PRICERHOMOGENEOUS);itsMapAdvisor.insert(type);
	type.Set(ICM_NTD,ICM_PRICER_HOMOGENEOUS_SMILE);itsMapAdvisor.insert(type);
	type.Set(ICM_MEZ,ICM_PRICERHOMOGENEOUS);itsMapAdvisor.insert(type);
	type.Set(ICM_MEZ,ICM_PRICER_HOMOGENEOUS_SMILE);itsMapAdvisor.insert(type);
	type.Set(ICM_CDO2,ICM_CLDPRICERCDO2);itsMapAdvisor.insert(type);
	type.Set(ICM_CDO2,ICM_MC_PRICER_CDO2);itsMapAdvisor.insert(type);
	type.Set(ICM_CAPFLOORCMCDS,ICM_PRICER_CAPFLOORCMCDS);itsMapAdvisor.insert(type);
	type.Set(ICM_CAPFLOORCMCDS,ICM_PRICER_CAPFLOORCMCDS);itsMapAdvisor.insert(type);
	type.Set(ICM_SPREADOPTION,ICM_BINOMIAL_TREE_PRICER_SPREADOPTION);itsMapAdvisor.insert(type);
	type.Set(ICM_CPPI,ICM_PRICER_CPPI);itsMapAdvisor.insert(type);
	type.Set(ICM_CORRIDORLEG,ICM_PRICER_CORRIDOR);itsMapAdvisor.insert(type);
	type.Set(ARM_INFLEG,ICM_PRICER_IR_INFCURVE);itsMapAdvisor.insert(type);
	type.Set(ARM_SWAPLEG,ICM_PRICER_IR_YCMOD);itsMapAdvisor.insert(type);
	type.Set(ARM_CMSLEG,ICM_PRICER_IR_YCMOD);itsMapAdvisor.insert(type);
	type.Set(ICM_CUSTOMIZED_CDO,ICM_CREDIT_MANAGER);itsMapAdvisor.insert(type);
	type.Set(ICM_CLN,ICM_GENERIC_NAME1);itsMapAdvisor.insert(type);
	type.Set(ICMOPTION_RESTRIKABLE,ICM_PRICER_TRANCHE_RESTRIKABLE);itsMapAdvisor.insert(type);
}

//	----------------------------------------------------------------------------
ARM_CLASS_NAME ICM_Pricer_Advisor::GetPricerType(ARM_Security* sec)
{
	ICM_PRICER_TYPE type(sec->GetName());
	if (itsMapAdvisor.find(type) == itsMapAdvisor.end()) return ARM_OBJECT;
	
	set<ICM_PRICER_TYPE>::iterator iterateur = itsMapAdvisor.find(type);

	return iterateur->m_pricer_type;
}


//	----------------------------------------------------------------------------
ICM_Pricer* 
ICM_Pricer_Advisor::GeneratePricer(ARM_Security* sec, 
									ARM_Object * mod_, 
									ARM_CLASS_NAME InputPricerType ,
									int nbpaths ,
									const ICM_Parameters* parameters_ ,
									// ICM_MktDataMng* MktDataManager ,
									const ARM_Date& asof )
{
	
	// to change:
	ICM_Parameters parameters; 
	if (parameters_) parameters=*parameters_; 

	ARM_CLASS_NAME PricerType;
	PricerType = InputPricerType;
	vector<string> DMKT;
	ARM_Model* mod=(ARM_Model*)mod_; // ugly... 
	ICM_MktDataMng* MktMng = dynamic_cast<ICM_MktDataMng*>(mod);
	

	if (InputPricerType == ARM_OBJECT) PricerType=GetPricerType(sec);
	if (PricerType == ARM_OBJECT) return NULL;

	ICM_Pricer* pricer = NULL;

	switch (PricerType)
	{
	case ICM_PRICER_CPPI :
	{
		if (!MktMng) 
			ICMTHROW(ERR_INVALID_ARGUMENT,"MarketDataManager Required"); 

		ICM_Pricer_Cppi tmp;tmp.MarketData(sec,DMKT);
		if (!(MktMng->CheckMktData(DMKT))) return NULL;

		ICM_Pricer_Cppi* pricer_tmp = new ICM_Pricer_Cppi; pricer_tmp->Set(sec, 
														  mod, 
														  parameters,
														  asof); 
		pricer = (ICM_Pricer*) pricer_tmp;
	}
	break;
	case ICM_MC_PRICER_CDO2 :
	{
		ICM_Pricer_MC_Cdo2* pricer_tmp = new ICM_Pricer_MC_Cdo2; pricer_tmp->Set(sec, 
																mod, 
																parameters,
																asof,
																nbpaths); 
		pricer = (ICM_Pricer*) pricer_tmp;
	}
	break;
	case ICM_PRICER_ANALYTIC_CDO2 :
	{
		ICM_Pricer_Analytic_Cdo2* pricer_tmp = new ICM_Pricer_Analytic_Cdo2; pricer_tmp->Set(sec,
																  mod,
																  parameters,asof);
		pricer = (ICM_Pricer*) pricer_tmp;
	}
	break;
	case ICM_PRICER_ANALYTIC_CDO2_STRIKE :
	{
		ICM_Pricer_Analytic_Cdo2_Smile	tmp;

		tmp.MarketData(sec, DMKT);

		ICM_Pricer_Analytic_Cdo2_Smile* pricer_tmp = new ICM_Pricer_Analytic_Cdo2_Smile; pricer_tmp->Set(sec,
																			mod,
																			parameters,asof);
		pricer = (ICM_Pricer*) pricer_tmp;
	}
	break;
	case ICM_PRICERHOMOGENEOUS :
	{
		ICM_Pricer_Distrib* pricer_tmp = new ICM_Pricer_Distrib; pricer_tmp->Set(sec,
																mod,parameters,
														  asof);
		pricer = (ICM_Pricer*) pricer_tmp;
	}
	break;
	case ICM_PRICER_CDS :
	{
		ICM_Pricer_Cds* pricer_tmp = new ICM_Pricer_Cds; pricer_tmp->Set(sec,
														mod,parameters,asof);
		pricer = (ICM_Pricer*) pricer_tmp;
	}
	break;
	case ICM_GENERIC_NAME1 :
	{
		ICM_Pricer_Cln* pricer_tmp = new ICM_Pricer_Cln; pricer_tmp->Set(sec,
														mod,parameters,asof);
		pricer = (ICM_Pricer*) pricer_tmp;
	}
	break;
	case ICM_PRICER_CDSOPTION :
	{
		ICM_Pricer_CdsOption* pricer_tmp = new ICM_Pricer_CdsOption; pricer_tmp->Set(sec,
																	mod,parameters,
														  asof);
		pricer = (ICM_Pricer*) pricer_tmp;
	}
	break;
	case ICM_PRICER_CDS_INDEX :
	{
		ICM_Pricer_CDSIndex* pricer_tmp = new ICM_Pricer_CDSIndex; pricer_tmp->Set(sec,
															mod,parameters,asof);
		pricer = (ICM_Pricer*) pricer_tmp;
	}
	break;
	case ICM_PRICER_INDEXOPTION :
	{
		ICM_Pricer_IndexOption* pricer_tmp = new ICM_Pricer_IndexOption; pricer_tmp->Set(sec,
																		mod,parameters,
														  asof);
		pricer = (ICM_Pricer*) pricer_tmp;
	}
	break;
	case ICM_PRICER_VANILLA_FWD:
	case ICM_PRICER_HOMOGENEOUS_SMILE:
	{
		ICM_Pricer_Distrib_Smile* pricer_tmp = new ICM_Pricer_Distrib_Smile; pricer_tmp->Set(sec,
																			mod,
																			parameters,asof);
		pricer = (ICM_Pricer*) pricer_tmp;

		if (PricerType==ICM_PRICER_VANILLA_FWD) 
		{CheckForForwardPricer(pricer);}
	}
	break;
	case ICM_PRICER_HOMOGENEOUS_SMILE_RF:
	{
		ICM_Pricer_Distrib_RandFactors* pricer_tmp = new ICM_Pricer_Distrib_RandFactors; pricer_tmp->Set(sec,
																			mod,
																			parameters,asof);
		pricer = (ICM_Pricer*) pricer_tmp;
	}
	break;
	case ICM_PRICER_HOMOGENEOUS_SMILE_COLLAT_FWD:
	{
		ICM_Pricer_Distrib_Smile_Collat_Fwd* pricer_tmp = new ICM_Pricer_Distrib_Smile_Collat_Fwd; pricer_tmp->Set(sec,
																			mod,
																			parameters,asof);
		pricer = (ICM_Pricer*) pricer_tmp;
	}
	break;
	case ICM_PRICER_CAPFLOORCMCDS:
	{
		ICM_Pricer_CapFloorCMCDS* pricer_tmp = new ICM_Pricer_CapFloorCMCDS; pricer_tmp->Set(sec,
																			mod,parameters,asof);
		pricer = (ICM_Pricer*) pricer_tmp;
	}
	break;
	case ICM_DIFFUSIONPRICER:
	{
		ICM_Pricer_Legs_Basket* pricer_tmp = new ICM_Pricer_Legs_Basket; pricer_tmp->Set(sec,
																		mod,
																		parameters,asof);
		pricer = (ICM_Pricer*) pricer_tmp;
	}
	break;
	case ICM_IMPLIED_LOSS_TREE_PRICER:
	{
		ICM_Pricer_ImpliedLoss tmp;
		//tmp.MarketData(sec,DMKT);
		if (!MktMng) 
			ICMTHROW(ERR_INVALID_ARGUMENT,"MarketDataManager Required"); 
		if (!(MktMng->CheckMktData(DMKT))) return NULL;
// 
		// ICM_Parameters localParams; 
		// if (parameters) localParams=*parameters; 
		// ARM_Date asOfDate; 
		// if (asof) asOfDate=asof;
		// else asOfDate=mod->GetAsOfDate(); 
		
		ICM_Pricer_ImpliedLoss* pricer_tmp = new ICM_Pricer_ImpliedLoss; pricer_tmp->Set(sec,
																		mod,parameters,asof);
		
		pricer = pricer_tmp;
	}
	break;
	case ICM_BINOMIAL_TREE_PRICER_CDS:
	{
		ICM_Pricer_Tree_Binomial_Cds tmp;
		//tmp.MarketData(sec,DMKT);
		if (!MktMng) 
			ICMTHROW(ERR_INVALID_ARGUMENT,"MarketDataManager Required"); 

		if (!(MktMng->CheckMktData(DMKT))) return NULL;

		// ICM_Parameters localParams; 
		// if (parameters) localParams=*parameters; 
		// ARM_Date asOfDate; 
		// if (asof) asOfDate=asof;
		// else asOfDate=mod->GetAsOfDate(); 
		
		ICM_Pricer_Tree_Binomial_Cds* pricer_tmp = new ICM_Pricer_Tree_Binomial_Cds; pricer_tmp->Set(sec,
																		mod,parameters,asof);
		
		pricer = pricer_tmp;
	}
	break;
	case ICM_BINOMIAL_TREE_PRICER_SPREADOPTION: 
	{
		// ICM_Parameters localParams; 
		// if (parameters) localParams=*parameters; 
		// ARM_Date asOfDate; 
		// if (asof) asOfDate=asof;
		ICM_Pricer_Tree_Binomial_SpreadOption * pricer_tmp = new ICM_Pricer_Tree_Binomial_SpreadOption; pricer_tmp->Set(sec,mod,parameters,asof); 
		pricer=pricer_tmp; 
	} break; 
	case ICM_PRICER_CORRIDOR :
	{
		ICM_Pricer_Corridor* pricer_tmp = new ICM_Pricer_Corridor; pricer_tmp->Set(sec,mod,parameters,asof);
		pricer = (ICM_Pricer*) pricer_tmp;
	}
	break;
/*	case ICM_PRICER_INDEXOPTION_SMILE :
	{
		ICM_Pricer_IndexOption_Smile tmp;
		if (!MktMng) 
			ICMTHROW(ERR_INVALID_ARGUMENT,"MarketDataManager Required"); 

		// ARM_Date asOfDate; 
		// if (asof) asOfDate=asof;
		// else asOfDate=mod->GetAsOfDate(); 

		tmp.SetAsOfDate(asof);
		tmp.MarketData(sec,DMKT);
		if (!(MktMng->CheckMktData(DMKT))) return NULL;

		ICM_Pricer_IndexOption_Smile* pricer_tmp = new ICM_Pricer_IndexOption_Smile ; pricer_tmp->Set(sec,
																		mod,
																		nbpaths,	// if -1 --> SmileToCompute
																		asof,
																		parameters);
		pricer = (ICM_Pricer*) pricer_tmp;
	}
	break;
	case ICM_PRICER_INDEXOPTION_SABR:
	{
		ICM_Pricer_IndexOption_SABR tmp;
		if (!MktMng) 
			ICMTHROW(ERR_INVALID_ARGUMENT,"MarketDataManager Required"); 

		// ARM_Date asOfDate; 
		// if (asof) asOfDate=asof;
		// else asOfDate=mod->GetAsOfDate(); 

		tmp.SetAsOfDate(asof);

		tmp.MarketData(sec,DMKT);
		if (!(MktMng->CheckMktData(DMKT))) return NULL;

		ICM_Pricer_IndexOption_SABR* pricer_tmp = new ICM_Pricer_IndexOption_SABR; pricer_tmp->Set(sec,
																		mod,
																		nbpaths,	// if -1 --> SmileToCompute
																		asof,
																		parameters);
		pricer = (ICM_Pricer*) pricer_tmp;
	}
	break;*/
	case ICM_PRICER_IR_INFCURVE :
	{
		ICM_Pricer_IR_InfCrv* pricer_tmp = new ICM_Pricer_IR_InfCrv; pricer_tmp->Set(sec,mod,parameters,asof);
		pricer = (ICM_Pricer*) pricer_tmp;
	}
	break;
	case ICM_PRICER_IR_YCMOD :
	{
		ICM_Pricer_IR_YcMod* pricer_tmp = new ICM_Pricer_IR_YcMod;pricer_tmp->Set(sec,mod,parameters,asof);
		pricer = (ICM_Pricer*) pricer_tmp;
	}
	break;

	case ICM_CREDIT_MANAGER :
	{
		CreditManager* pricer_tmp = new CreditManager;pricer_tmp->Set(sec, mod, parameters,asof);
		pricer = (ICM_Pricer*) pricer_tmp;
	}
	break;
/*	case ICM_PRICER_INDEX_CAPGEN_BLACKSCHOLES :
	{
		ICM_Pricer_Index_CapGen_BlackScholes tmp;
		if (!MktMng) 
			ICMTHROW(ERR_INVALID_ARGUMENT,"MarketDataManager Required"); 

		// ARM_Date asOfDate; 
		// if (asof) asOfDate=asof;
		// else asOfDate=mod->GetAsOfDate(); 

		tmp.SetAsOfDate(asof);

		tmp.MarketData(sec,DMKT);
		if (!(MktMng->CheckMktData(DMKT))) return NULL;

		ICM_Pricer_Index_CapGen_BlackScholes* pricer_tmp = new ICM_Pricer_Index_CapGen_BlackScholes;pricer_tmp->Set(sec,
																		mod,
																		nbpaths,	// if -1 --> SmileToCompute
																		asof,
																		parameters);
		pricer = (ICM_Pricer*) pricer_tmp;
	}
	break;*/
/*	case ICM_PRICER_INDEX_CAPGEN_BLACKSCHOLES_STRIKES :
	{
		ICM_Pricer_Index_CapGen_BlackScholes_Strikes tmp;
		if (!MktMng) 
			ICMTHROW(ERR_INVALID_ARGUMENT,"MarketDataManager Required"); 

		// ARM_Date asOfDate; 
		// if (asof) asOfDate=asof;
		// else asOfDate=mod->GetAsOfDate(); 

		tmp.SetAsOfDate(asof);

		tmp.MarketData(sec,DMKT);
		if (!(MktMng->CheckMktData(DMKT))) return NULL;

		ICM_Pricer_Index_CapGen_BlackScholes_Strikes* pricer_tmp = new ICM_Pricer_Index_CapGen_BlackScholes_Strikes;pricer_tmp->Set(sec,
																		mod,
																		nbpaths,	// if -1 --> SmileToCompute
																		asof,
																		parameters);
		pricer = (ICM_Pricer*) pricer_tmp;
	}
	break;*/
/*	case ICM_PRICER_INDEX_CAPGEN_SMILE :
	{
		ICM_Pricer_Index_CapGen_Smile tmp;
		if (!MktMng) 
			ICMTHROW(ERR_INVALID_ARGUMENT,"MarketDataManager Required"); 

		// ARM_Date asOfDate; 
		// if (asof) asOfDate=asof;
		// else asOfDate=mod->GetAsOfDate(); 

		tmp.SetAsOfDate(asof);

		tmp.MarketData(sec,DMKT);
		if (!(MktMng->CheckMktData(DMKT))) return NULL;

		ICM_Pricer_Index_CapGen_Smile* pricer_tmp = new ICM_Pricer_Index_CapGen_Smile;pricer_tmp->Set(sec,
																		mod,
																		nbpaths,	// if -1 --> SmileToCompute
																		asof,
																		parameters);
		pricer = (ICM_Pricer*) pricer_tmp;
	}
	break;*/
	case ICM_PRICER_LSS_GAP_OPTION :
	{
		ICM_Pricer_LssGapOption* pricer_tmp = new ICM_Pricer_LssGapOption;pricer_tmp->Set(sec,mod,parameters,nbpaths);
		pricer = (ICM_Pricer*) pricer_tmp;
	}
	break;
	case ICM_PRICER_CPDO :
	{

		// ICM_Parameters localParams; 
		// if (parameters) localParams=*parameters; 

		// ARM_Date asOfDate; 
		// if (strcmp(asof,"NULL")) 
		// 	asOfDate=asof;
		// else 
		// 	asOfDate=mod->GetAsOfDate(); 


		ICM_Pricer_CPDO* pricer_tmp = new ICM_Pricer_CPDO;pricer_tmp->Set(sec,mod,parameters, asof);
		pricer = (ICM_Pricer*) pricer_tmp;
	} 
	break; 
	case ICM_BASEPRICER:
	{
		ICM_BasePricer * pricer_tmp = new ICM_BasePricer(); 
		pricer_tmp->Set(dynamic_cast<ICM_Credit_Index*>(sec),
					dynamic_cast<ICM_ModelMultiCurves*>(mod),parameters,asof); 
		pricer=pricer_tmp; 
	} 
	break; 
	case ICM_PRICER_MC_CDO :
	{
		ICM_Pricer_MC_CDO* pricer_tmp = new ICM_Pricer_MC_CDO(); 
		pricer_tmp->Set(sec,mod,parameters,asof,nbpaths); 
		pricer = (ICM_Pricer*) pricer_tmp;
	}
	break;
	case ICM_PRICER_TREE_HW_CDO:
	{
		ICM_Pricer_tree_hw_cdo * pricer_tmp = new ICM_Pricer_tree_hw_cdo(); 
		pricer_tmp->Set(sec,
					dynamic_cast<ICM_ModelMultiCurves*>(mod),parameters,asof); 
		pricer=pricer_tmp; 
	}
	break; 

	case ICM_PRICER_TRANCHE_RESTRIKABLE :

		ICM_Pricer_Tranche_Restrikable * pricer_tmp = new ICM_Pricer_Tranche_Restrikable();
		pricer_tmp->Set(sec,mod,parameters,asof);
		pricer = pricer_tmp;

}

return (pricer);
}


// 
void ICM_Pricer_Advisor::CheckForForwardPricer(ICM_Pricer*& initpricer)
{

	ICM_Pricer* pricer = NULL;
	
	switch (initpricer->GetName())
	{
		case ICM_PRICERHOMOGENEOUS:
		case ICM_PRICER_HOMOGENEOUS_SMILE:
		{	
			pricer = new ICM_VanillaFwdPricer(initpricer);
			if (initpricer) delete initpricer; initpricer=pricer;	
		}
		break;
	}
	
}

