
#ifndef _ICM_Pricer_CDSIndex_H
#define _ICM_Pricer_CDSIndex_H

#include "ICMKernel\pricer\icm_pricer_cds.h"

/*********************************************************************************/
/*! \class  ICM_Pricer_CDSIndex icm_pricer_index.h "icm_pricer_index.h"
 *  \author Fakher Ben Atig
 *	\version 1.0
 *	\date   September 2004
 *	\file   ICM_Pricer_CDSIndex.h
 *	\brief  Pricer for Credit Index Contracts  */
/***********************************************************************************/

class ICM_Pricer_CDSIndex : public ICM_Pricer_Cds
{
private:
	ICM_DefaultCurve* itsDefCurveIndex;		//0,1, aggregated, cached. 
	bool itsFlatCalculation;
private: 
	virtual double FeeLegPV () ;
	virtual double DefLegPV () ;
protected:
	const ICM_DefaultCurve& GetDefCurveIndex(); // Compute DefCurve of the Index
	void CptPricingMeFlat(qCMPMETH measure);
	void ComputeUpfrontPay(qCMPMETH measure); 
	
public :
	ICM_Pricer_CDSIndex() ; // { itsDefCurveIndex=NULL; Init();}
	~ICM_Pricer_CDSIndex() ; // {}
	void Init();
	ICM_Pricer_CDSIndex(ARM_Security *option, ARM_Model *mod,const ICM_Parameters&params,const ARM_Date&AsOf)
	{
		Init();
		Set(option, mod,params,AsOf);
	}
	void Set(ARM_Security *option, ARM_Model *mod,const ICM_Parameters&params,const ARM_Date&AsOf);
	virtual void Reset(void);
	virtual double ComputeSpread(const double& MtM = 0.) ;
	

	double Compute_Fwd_Spread(const ARM_Date &Mty, const ARM_Date & ExpiryDate,double& FlatRbp_Fwd);
	double GetVirtualCDSSpread(const ARM_Date &Mty) ;
	virtual double ComputeDuration(void);
	void View(char* id, FILE* ficOut);
private: 
	ICM_Pricer_CDSIndex(const ICM_Pricer_CDSIndex&ref);			// NA
	ICM_Pricer_CDSIndex& operator=(const ICM_Pricer_CDSIndex&); // NA 
};

#endif

