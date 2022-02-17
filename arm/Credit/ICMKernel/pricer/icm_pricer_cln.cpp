#include "firsttoinc.h"

#include "ICMKernel/pricer/icm_pricer_cln.h"



	ICM_Pricer_Cln::ICM_Pricer_Cln(void) 
	{ 
		Init();
	}

	/** ICM_Pricer_Cln::ICM_Pricer_Cln(ARM_Security *sec, ARM_Model *mod,const ICM_Parameters&params,const ARM_Date& asof)
	{	
		Init();
		Set(sec, mod,params,asof);	
	}
	**/ 
	ICM_Pricer_Cln::~ICM_Pricer_Cln() 
	{}

	void ICM_Pricer_Cln::Init()	
	{
		SetName(ICM_GENERIC_NAME1);
	}

	void ICM_Pricer_Cln::Set(ARM_Security *sec, ARM_Model *mod,const ICM_Parameters&params,
		const ARM_Date& asof)
	{ 	
	 //on crée un Cds fictif ayant pour feeleg le cds
	 ICM_Cds GhostCds(dynamic_cast<ICM_Leg*>(sec),dynamic_cast<ICM_Leg*>(sec));

	 int nxchg = dynamic_cast<ICM_Leg*>(sec)->GetNxFlag();
	 if ((nxchg==K_NX_START)||
		 (nxchg==K_NX_BOTH)) {SetRefNotional(0.);}

	 //on switch cln & cds
	 ICM_Pricer_Cds::Set(&GhostCds,mod,params,asof); 

	 SetSecurity(sec);
	}

	// virtual 
	double ICM_Pricer_Cln::FeeLegPV () 
	{
		if (GetFeeLegPriceFlg()) return GetFeeLegPrice();
		ICM_Cln* security =  (ICM_Cln*) GetSecurity();
		
		//on crée un Cds fictif ayant pour feeleg le cds
		ICM_Leg Ghost;Ghost.SetCreditLegStyle(qStyle_None_Leg);
		ICM_Cds GhostCds((ICM_Leg*)security,&Ghost);

		//on switch cln & cds
		SetSecurity(&GhostCds);
		double FeePv = ICM_Pricer_Cds::FeeLegPV(); 

		((ICM_Leg*)security)->BitwiseCopy(GhostCds.GetFeeLeg());

		SetSecurity(security);
		return (FeePv);
	}

	// virtual 
	double ICM_Pricer_Cln::Accrued()
	{
		// if (GetAccruedFlg()) return GetAccrued();
		ICM_Cln* security =  (ICM_Cln*) GetSecurity();
		
		//on crée un Cds fictif ayant pour feeleg le cds
		ICM_Leg Ghost;Ghost.SetCreditLegStyle(qStyle_None_Leg);
		ICM_Cds GhostCds((ICM_Leg*)security,&Ghost);

		//on switch cln & cds
		SetSecurity(&GhostCds);
		double _AccruedPV = ICM_Pricer_Cds::Accrued(); 

		SetSecurity(security);
		return (_AccruedPV);
	}

	// virtual 
	double ICM_Pricer_Cln::DefLegPV () 
	{
		if (GetDefLegPriceFlg()) return GetDefLegPrice();

		double result = 0.;
		ICM_Cln* security =  (ICM_Cln*) GetSecurity();

		//on crée un Cds fictif ayant pour feeleg le cds
		std::auto_ptr<ICM_Leg> Ghost((ICM_Leg*)security->Clone()); 
		Ghost->SetCreditLegStyle(qStyle_None_Leg);
		ICM_Cds GhostCds(Ghost.get(),(ICM_Leg*)security);

		if (GetAsOfDate() <= security->GetStartDate())
		{	int nxchg = security->GetNxFlag();
			//in cln case Reference notional is =0	
			if ((nxchg==K_NX_START)||
			    (nxchg==K_NX_BOTH)) 
				{result += security->GetRedempNotional();}}

		//on switch cln & cds
		SetSecurity(&GhostCds);
		
		result += ICM_Pricer_Cds::DefLegPV ();
		SetDefLegPrice(result);

		SetSecurity(security);

		return (result); }

	void ICM_Pricer_Cln::Reset(void)
	{ 
		ICM_Pricer_Cds::Reset();
	}

//	virtual void ComputeDuration(void)
//	{
//	}


