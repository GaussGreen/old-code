#include "firsttoinc.h"
#include "icm_pricer_ir_infcurve.h" 
//#include "ARMKernel\inst\swapleg.h" 
#include "gpcalculators\pricerfactory.h" 
#include "gpinflation/infleg.h"
#include "ARMKernel/pricer/ipricer.h"

 
 
ICM_Pricer_IR_InfCrv::ICM_Pricer_IR_InfCrv()
{ Init() ;}

void ICM_Pricer_IR_InfCrv::Set(ARM_Security *sec, ARM_Model *mod,const ICM_Parameters& params,const ARM_Date& asof)
	{
		ICM_Pricer::Set(sec, mod,params,&asof);	
	}

void ICM_Pricer_IR_InfCrv::Init()
	{SetName(ICM_PRICER_IR_INFCURVE);}


// virtual 
double ICM_Pricer_IR_InfCrv::ComputePrice(qCMPMETH measure)
    {
		double price = 0.;
		/// call factory
		CC_NS(std,pair)<bool,double> pricingResult = ARM::ARM_PricerFactory::Price( GetSecurity(), GetModel() );

		if( pricingResult.first)
			price = pricingResult.second;
		else
		{
			/// this is a simple pricer!
			ARM_IFPricer pricer(GetSecurity(),GetModel());
			price = pricer.Price();
		}

		return price;
	}

// virtual 
void ICM_Pricer_IR_InfCrv::GenerateRates(ARM_Vector& rates)
	{
		double value = Price(qCMPPRICE);
		ARM::ARM_InfLeg* leg = dynamic_cast<ARM::ARM_InfLeg*>(GetSecurity());
		rates = *leg->GetDecompRates();
	}

// 	virtual 
double ICM_Pricer_IR_InfCrv::ComputeSensitivity(const qSENSITIVITY_TYPE& typesensi, 
		const std::string& plot,
		const std::string& label,
		double epsilon,  double epsilonGamma )  
	{
		return ICM_Pricer::ComputeSensitivity(typesensi,plot,label,epsilon, epsilonGamma); 
	}
 