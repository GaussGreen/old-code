
#ifndef _ICM_Pricer_SECURITY_H
#define _ICM_Pricer_SECURITY_H

#include "ICMKernel/pricer/icm_pricer.h"
#include "ICMKernel/util/icm_utils.h"
#include "ICMKernel/crv/icm_distribloss.h"
#include "ICMKernel/inst/icm_cds.h"

/*********************************************************************************/
/*! \class  ICM_Pricer_Security icm_pricer_security.h "icm_pricer_security.h"
 *  \author Damien Pouponneau 
 *	\version 1.0
 *	\date   March 2004 
 *	\file   icm_pricer_security.h
 *	\brief  Pricer for all security using ICM_Security */
/***********************************************************************************/

#ifndef STEP_DISTRIB_LOSS 
#define STEP_DISTRIB_LOSS  90
#endif

class ICM_Pricer_Security : public ICM_Pricer
{
private:

	ICM_DistribLoss  itsDistribLoss ;		//Loss Distribution
	int				 itsStepForDistribLoss ;//Minimal step for loss distrib computation	

	qTERM_STRUCTURE	 itsIsTermStructurePricing;//Is a term structure correlation
	ARM_Vector		 itsTermStructureSched; //Schedule Used for term structure pricing 	

	vector<double>	 itsScheduleForEL;		//full schedule used to compute EL
	vector<vector<double> >	itsDetailLosses;//unit expected loss matrix

public :
	ICM_Pricer_Security(void) 
	{ Init(); }

	/** ICM_Pricer_Security(ARM_Security *sec, ARM_Model *mod,const ICM_Parameters&params,const ARM_Date&asof)
	{
		Init();
		Set(sec, mod,params,asof);
	}
	**/ 
	~ICM_Pricer_Security()
	{
		itsScheduleForEL.clear();
		itsDetailLosses.clear();
	}

	void Init()
	{
		SetName(ICM_PRICERSECURITY);
		itsDistribLoss.clear(); 
		itsStepForDistribLoss=STEP_DISTRIB_LOSS;
		itsIsTermStructurePricing = qNoTermStructure;
	}

	void Set(ARM_Security *sec, ARM_Model *mod,const ICM_Parameters&,const ARM_Date&);

	// Pricer API ... 
private:
	virtual double FeeLegPV ();
	virtual double DefLegPV ();
	virtual double Accrued();
	virtual double  ComputeDuration(void);
public:
	virtual double ComputeSpread(const double& MtM = 0.);
public:
	// Pricer Security Specifics
	virtual double ExpectedLossTranche(const double& yearterm,vector<double>& losses);
	virtual double GetTranche(const ARM_Date& date) {return 1.;};

	//Expected Amortissement : uniquement estimé sur les tranhces [x-100%].
	//On considère que tout le risque d'amotissement est associé à cette tranche uniquement.
	virtual double ExpectedAmortTranche(const double& yearterm) { return 1.;}


	virtual void CptExpectedLossTranche();

	virtual void Reset(void)
	{itsDistribLoss.clear();}

	
	
	double ComputeImpliedPartRate();

	

	void View(char* id, FILE* ficOut);
	const ICM_DistribLoss& getDistribLoss()  { CptExpectedLossTranche(); return itsDistribLoss; }
	const ICM_DistribLoss& getDistribLossObj()  { return itsDistribLoss; }
	void setDistribLoss(const ICM_DistribLoss&ref) { itsDistribLoss=ref; }
	void clearDistribLoss() { itsDistribLoss.clear() ; }

	int	GetStepForDistribLoss(){ return itsStepForDistribLoss ;}	
	void SetStepForDistribLoss(const int& step){itsStepForDistribLoss=step;}	

	/** JLA: abstract class, dead code.
	virtual ICM_Pricer* CloneOnlyParams(ARM_Security* sec,ARM_Model* mod)
	{ICM_Pricer* theClone = new ICM_Pricer_Security(sec,mod);
	 return(theClone);}
	 **/ 

	ARM_Vector* GetRiskSchedule(void)
	{	ICM_Cds* cds = (ICM_Cds*) GetSecurity();
		ARM_Date AsOf = GetModel()->GetStartDate();
		cds->GenerateRiskSchedule(AsOf,itsStepForDistribLoss);
		return cds->GetRiskSchedule();	}

	ARM_Vector* GetSchedule(void)
	{	ICM_Cds* cds = (ICM_Cds*) GetSecurity();
		ARM_Date AsOf = GetModel()->GetStartDate();

		//modification pour leg zerocoupon
		double StepForDistribLoss=itsStepForDistribLoss;
		ICM_Security* security = cds->GetFeeLeg()->GetCreditInfos();
		if (security->GetPaymentFreq()==0)
			{StepForDistribLoss=-1;}

		cds->GenerateSchedule(AsOf,StepForDistribLoss);
		return cds->GetSchedule();	}

	double RiskyNotional(const int& index);

	inline void SetTermStructurePricing(qTERM_STRUCTURE value) {itsIsTermStructurePricing = value;}
	inline qTERM_STRUCTURE GetTermStructurePricing() {return itsIsTermStructurePricing;}

	vector<double>&	 GetScheduleForEL() {return itsScheduleForEL;}
	vector<vector<double> >& GetDetailLosses() {return itsDetailLosses;}


private:
	//	Unimplemented
	ICM_Pricer_Security(const ICM_Pricer_Security&ref); 
	ICM_Pricer_Security& operator=(const ICM_Pricer_Security&); 
	//	Throw
	virtual double Compute_Fwd_Spread(const ARM_Date &Mty, const ARM_Date &CDS_ExpiryDate, double& dur) ; 
 	virtual	double ComputeImpliedVol(const double& Price) ; 
};


#endif

