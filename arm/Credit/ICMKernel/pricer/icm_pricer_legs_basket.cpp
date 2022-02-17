#include "firsttoinc.h"
#include "ICMKernel\pricer\icm_pricer_legs_basket.h"
#include "ICMKernel\pricer\icm_pricer_adviser.h"
#include "ICMKernel\util\icm_macro.h"
//#include "ARMKernel\mod\model.h"
#include "ICMKernel\mod\icm_meta_model.h"
#include "ICMKernel\mod\icm_defcurvemodel.h"
#include "ICMKernel\inst\icm_leg_basket.h"

void ICM_Pricer_Legs_Basket::Set(ARM_Security *sec, ARM_Model *mod,const ICM_Parameters& parameters,const ARM_Date&asof)
{
	ICM_Pricer::Set(sec,mod,parameters,&asof);

	ICM_Meta_Model* Model = (ICM_Meta_Model*) mod;
	ICM_Leg_Basket* Security = (ICM_Leg_Basket*) sec;

	int i;
	int Nbmod = Model->GetNbModel();
	int Nbsec = Security->GetNbSec();

	if ((Nbmod>1) && (Nbmod != Nbsec))
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "ERROR : Models number <> Securities number ");

	itsModels.resize(Nbsec);	
	itsPricerType.resize(Nbsec);	
	itsSecurities.resize(Nbsec);	

	for (i=0; i<itsPricers.size(); i++) 
	{
		if (itsPricers[i])delete itsPricers[i];
			itsPricers[i] = NULL;
	}

	itsPricers.resize(Nbsec);	

	vector<ARM_CLASS_NAME> PricerType;	

	Model->GetPricerType(PricerType);

	ICM_Pricer_Advisor Advisor;

	for (i=0;i<Nbsec;i++)
	{
		if (Nbmod == 1)
		{
		itsPricerType[i] = PricerType[0];
		itsModels[i] = (ICM_DefaultCurveModel*) Model->GetModels()[0];
		}
		else
		{
		itsPricerType[i] = PricerType[i];
		itsModels[i] = (ICM_DefaultCurveModel*) Model->GetModels()[i];
		}

		itsSecurities[i] = Security->GetSecurities()[i];
		itsPricers[i] = Advisor.GeneratePricer(itsSecurities[i],itsModels[i],itsPricerType[i],CREDIT_DEFAULT_VALUE,(ICM_Parameters*)0,asof);
	}
}


double ICM_Pricer_Legs_Basket::ComputePrice(qCMPMETH measure)
{
	int i;
	int size = itsSecurities.size();
	double price = 0;

	for (i=0;i<size;i++)
		price += itsPricers[i]->Price(qCMPPRICE);

	return (price);
}


double ICM_Pricer_Legs_Basket::ComputeSensitivity(const qSENSITIVITY_TYPE& typesensi, 
												  const std::string&  plot, 
												  const std::string& label, 
												double  epsilon,  double epsilonGamma //useless 
												)
{
	int i;
	int size = itsSecurities.size();
	double price = 0;

	for (i=0;i<size;i++)
		price += itsPricers[i]->Hedge(typesensi,plot,label,epsilon);

	return (price);
}


double ICM_Pricer_Legs_Basket::Accrued(void) 
{
	ICMTHROW(ERR_INVALID_ARGUMENT,"Not Implemented") ; 
}

double ICM_Pricer_Legs_Basket::FeeLegPV(void)
{
	ICMTHROW(ERR_INVALID_ARGUMENT,"Not Implemented") ; 
}

double ICM_Pricer_Legs_Basket::DefLegPV(void)
{
	ICMTHROW(ERR_INVALID_ARGUMENT,"Not Implemented") ; 
}

double ICM_Pricer_Legs_Basket::Compute_Fwd_Spread(const ARM_Date &,const ARM_Date &, double& dur)
{
	ICMTHROW(ERR_INVALID_ARGUMENT,"Not Implemented") ; 
}

double ICM_Pricer_Legs_Basket::ComputeDuration(void)
{
	ICMTHROW(ERR_INVALID_ARGUMENT,"Not Implemented") ; 
}

double ICM_Pricer_Legs_Basket::ComputeSpread(const double &)
{
	ICMTHROW(ERR_INVALID_ARGUMENT,"Not Implemented") ; 
}

double ICM_Pricer_Legs_Basket::ComputeImpliedVol(const double &)
{
	ICMTHROW(ERR_INVALID_ARGUMENT,"Not Implemented") ; 
}

