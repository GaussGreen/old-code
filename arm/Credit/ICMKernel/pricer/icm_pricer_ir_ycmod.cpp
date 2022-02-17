#include "firsttoinc.h"
#include "icm_pricer_ir_ycmod.h" 
#include "ARMKernel\inst\swapleg.h" 
#include "gpcalculators\pricerfactory.h" 


	ICM_Pricer_IR_YcMod::ICM_Pricer_IR_YcMod(void) 
	{ Init();}

 
	void ICM_Pricer_IR_YcMod::Set(ARM_Security *sec, ARM_Model *mod,const ICM_Parameters&params,const ARM_Date&asof)
	{
		ICM_Pricer::Set(sec, mod,params,&asof);	
	}

	ICM_Pricer_IR_YcMod::~ICM_Pricer_IR_YcMod(){}

	void ICM_Pricer_IR_YcMod::Init()
	{SetName(ICM_PRICER_IR_YCMOD);}

	// void Set(ARM_Security *sec, ARM_Model *mod){}
	// virtual 
	void ICM_Pricer_IR_YcMod::Reset(void) {}
	// virtual 
	void ICM_Pricer_IR_YcMod::View(char* id, FILE* ficOut) {}

	// virtual 
	void ICM_Pricer_IR_YcMod::GenerateRates(ARM_Vector& rates)
	{
		double value = Price(qCMPPRICE);
		ARM_SwapLeg* leg = dynamic_cast<ARM_SwapLeg*>(GetSecurity());
		rates = *leg->GetDecompRates();
	}

	// virtual 
	double ICM_Pricer_IR_YcMod::ComputeSpread(const double& MtM )
	{ 
		double feelegpv = ComputePrice(qCMPFEELEGPV);
		double duration = ComputePrice(qCMPDURATION);
		double breakeven = fabs(feelegpv)/duration;

		breakeven/=10000.;
		SetSpread(breakeven);

		return (breakeven);
	}


	// virtual 
	double ICM_Pricer_IR_YcMod::ComputeDuration(void)
    {
		double price = 0.;
		/// call factory
		CC_NS(std,pair)<bool,double> pricingResult = ARM::ARM_PricerFactory::Price( GetSecurity(), GetModel() );

		if( pricingResult.first)
			price = pricingResult.second;
		else
		{
			/// this is a simple pricer!
			ARM_SwapLeg* sl = dynamic_cast<ARM_SwapLeg*>(GetSecurity());
			ARM_SwapLeg* slfixed = (ARM_SwapLeg*) sl->Clone();
			slfixed->SetLegType(K_FIXED_LEG);
			slfixed->SetFixedRate(1.);
			slfixed->SetNotionalAmount(0.);
			slfixed->SetAmount(&ARM_ReferenceValue(1.));

			ARM_IFPricer pricer(slfixed,GetModel());
			price = 100.*pricer.Price();
			if (slfixed) delete slfixed;
		}

		SetDuration(fabs(price));
		return fabs(price); 
	}

	// virtual 
	double ICM_Pricer_IR_YcMod::FeeLegPV ()
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

		SetFeeLegPrice(price);
		return price;}

	// virtual 
	double ICM_Pricer_IR_YcMod::DefLegPV ()
	{
		SetDefLegPrice(0.);
		return 0.;
	}


	// virtual 
	double ICM_Pricer_IR_YcMod::Accrued()
	{
		double price = 0.;
		double Notional = 0.;

		/// call factory
		CC_NS(std,pair)<bool,double> pricingResult = ARM::ARM_PricerFactory::Price( GetSecurity(), GetModel() );

		if( pricingResult.first)
			price = pricingResult.second;
		else
		{
			/// this is a simple pricer!
			ARM_IFPricer pricer(GetSecurity(),GetModel());
			price = pricer.Price();
			price = GetSecurity()->Accrued(GetModel()->GetStartDate());
		}

		Notional = GetSecurity()->GetAmount()->CptReferenceValue(GetModel()->GetStartDate());
		price *= Notional/100.;

		SetAccrued(price);
		return price;	 	
	}

	// virtual 
	double ICM_Pricer_IR_YcMod::ComputeSensitivity(const qSENSITIVITY_TYPE& typesensi, 
										const std::string& plot, 
										const std::string& label , 
										double  epsvalue ,  double epsilonGamma )
	{
    double sensitivity =0.;

	ARM_YCModel* ycModel = dynamic_cast<ARM_YCModel*>(GetModel());
	ARM_Date ExecutionDate = ycModel->GetStartDate();

	int i = 0;
	double result = 0.,initialprice = 0.,modifiedprice = 0.;
	double EPSL = 0.01;		//Bump sur les taux fixé à 1bp	
	bool Parallelshift = false; 

	ARM_CRV_TERMS Term;
	memset(Term,'\0',sizeof(char)*ARM_NB_TERMS*ARM_NB_MAX_CHAR_TERMS);
	ARM_Vector epsilon(1,0.0);

	if (epsvalue != CREDIT_DEFAULT_VALUE) {EPSL = epsvalue;}

	// if (!strcmp(plot,"NONE"))  //Detection d'un parallel Shift
	if (plot=="NONE") 
		Parallelshift = true;
	else
	{
		strcpy(Term[0],plot.c_str());
		switch (typesensi)
		{
			case ICMIRCURVE_WITHOUT_DEFCURVE_TYPE :
			{epsilon.InitElt(0,EPSL);}
			break;
		}
	}

    
    {

		if (GetInitialPriceFlg())
			initialprice = GetInitialPrice();
		else
			initialprice = Price(qCMPPRICE);

		ResetPricer();

		switch (typesensi)
		{
			case ICMIRCURVE_WITHOUT_DEFCURVE_TYPE :
 			{
	
				ARM_ZeroCurve* zc = (ARM_ZeroCurve*) ycModel->GetZeroCurve();
				ARM_ZeroCurve* shiftzc = NULL;

				if (Parallelshift)
				{	
					shiftzc = (ARM_ZeroCurve*)(zc->Clone());
					shiftzc->ParallelShift(EPSL);
				}
				else
				{shiftzc = (ARM_ZeroCurve*) zc->GenerateShiftCurve(Term,&epsilon);}
				
				ARM_YCModel shiftmod(shiftzc);
				SetModel(&shiftmod);	

				modifiedprice = ComputePrice(qCMPPRICE);
				SetModel(ycModel);	

				if (shiftzc) delete shiftzc;
				shiftzc = NULL;
			}
			break;
			default :
			{
			ResetPricer();
			return 0.;
			}
		}

	ResetPricer();

	if (!result) result=modifiedprice - initialprice;

	}
   

	return (result);
	}
