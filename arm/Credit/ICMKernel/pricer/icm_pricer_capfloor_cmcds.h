
#ifndef _ICM_Pricer_CapFloorCMCDS_H
#define _ICM_Pricer_CapFloorCMCDS_H

#include "ICMKernel\inst\icm_cmcds.h"
#include "ICMKernel\pricer\icm_pricer_index.h"
#include "ARMKernel\crv\volflat.h"

/*********************************************************************************/
/*! \class  ICM_Pricer_CapFloorCMCDS icm_pricer_capfloor_cmcds.h "icm_pricer_capfloor_cmcds.h"
 *  \author Fakher Ben Atig
 *	\version 1.0
 *	\date   September 2004
 *	\brief  Pricer for CMCDS with a Cap or/and a Floor */
/***********************************************************************************/

class ICM_Pricer_CapFloorCMCDS : public ICM_Pricer
{
private:

	ARM_Vector  itsAdjFwdSpreads ; //Convexity (and Payment Lag) Adjusted Fwd Spreads for CMCDS products
	ARM_Vector  itsFwdSpreads ;	//Fwd Spreads for CMCDS products (with no adjustments)
	ARM_Vector  itsFwdPV01s ;		//Fwd PV01s for CMCDS products
	

	bool AdjFwdSpreadsFlg ;			//Adjusted Fwd Spreads Flag
	bool FwdSpreadsFlg ;			//Unadjusted Fwd Spreads Flag
	bool FwdPV01sFlg ;				//Fwd PV01s Flag

	double itsCapValue ;
	double itsFloorValue ;
	
	bool CapValueFlg ;
	bool FloorValueFlg ;


public :
	ICM_Pricer_CapFloorCMCDS(void)
	{
		Init();
	}

	/** ICM_Pricer_CapFloorCMCDS(ARM_Security *capfloorcmcds, ARM_Model *mod,const ICM_Parameters&params,const ARM_Date&asof)
	{
		Init();

		Set(capfloorcmcds, mod,params,asof);
	}
	**/ 
	virtual ~ICM_Pricer_CapFloorCMCDS()
	{
	
	}

	void Init()
	{
		SetName(ICM_PRICER_CAPFLOORCMCDS);
						
		//itsAdjFwdSpreads	= NULL ;
		//itsFwdSpreads		= NULL ;
		//itsFwdPV01s			= NULL ;
		AdjFwdSpreadsFlg	= false;
		FwdSpreadsFlg		= false;
		FwdPV01sFlg			= false;
		itsCapValue = -1.;
		itsFloorValue = -1.;
		CapValueFlg = false ;
		FloorValueFlg = false ;
		
		
		
	}

	void Set(ARM_Security *capfloorcmcds, ARM_Model *mod,const ICM_Parameters&params,const ARM_Date& asof); 

	virtual void Reset(void)
	{
		ICM_Pricer::Reset();
			
		AdjFwdSpreadsFlg	= false;
		FwdSpreadsFlg		= false;
		FwdPV01sFlg			= false;
		
		itsCapValue = -1. ;
		itsFloorValue = -1.;
		CapValueFlg = false ;
		FloorValueFlg = false ;
	}

	
	void CapFloorBSPrices();

	void SetCapValue (double value)
	{
		itsCapValue = value ;
		CapValueFlg	= true ;
	}

	double GetCapValue()
	{
		if (!CapValueFlg)
			CapFloorBSPrices();
		return itsCapValue ;
	}
	
	void SetFloorValue (double value)
	{
		itsFloorValue = value ;
		FloorValueFlg	= true ;
	}

	double GetFloorValue()
	{
		if (!FloorValueFlg)
			CapFloorBSPrices();
		return itsFloorValue ;
	}

	//virtual double ComputeImpliedVol(double Price);

	const ARM_Vector& GetAdjFwdSpreads(void) ;
	
	void SetAdjFwdSpreads(const ARM_Vector& AdjFwdSpreads) 
	{ 
		// if (itsAdjFwdSpreads)
		// 	delete itsAdjFwdSpreads;
		itsAdjFwdSpreads = AdjFwdSpreads; // (ARM_Vector*)AdjFwdSpreads.Clone();

		AdjFwdSpreadsFlg	= true;
	}

	const ARM_Vector& GetFwdSpreads(void) ;
	
	void SetFwdSpreads(const ARM_Vector& FwdSpreads) 
	{ 
		// if (itsFwdSpreads)
		//	delete itsFwdSpreads;
		// itsFwdSpreads = (ARM_Vector*)FwdSpreads.Clone();
		itsFwdSpreads =FwdSpreads ;

		AdjFwdSpreadsFlg	= true;
	}

	const ARM_Vector& GetFwdPV01s(void)  ;

	void SetFwdPV01s(const ARM_Vector& FwdPV01s) 
	{ 
		// if (itsFwdPV01s)
		//	delete itsFwdPV01s;
		// itsFwdPV01s= (ARM_Vector*)FwdPV01s.Clone();
		itsFwdPV01s = FwdPV01s ;
		
		FwdPV01sFlg			= true;
	}
	
	virtual double ComputeSpread(const double& MtM = 0.);
	private:
		virtual double FeeLegPV ();
		virtual double DefLegPV ();
	protected:
	virtual double ComputeSensitivity(const qSENSITIVITY_TYPE& typesensi, 
		const std::string& plot,
		const std::string& label,
		double epsilon, double epsilonGamma =0 ); 
	public:

	 
private:
	virtual double Accrued(); 
	virtual double Compute_Fwd_Spread(const ARM_Date&,const ARM_Date&, double& dur); 
	virtual double ComputeDuration() ; 
	virtual double ComputeImpliedVol(const double&); 

};

#endif
