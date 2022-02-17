

#ifndef _ICM_Pricer_SECURITY_HOMOGENEOUS_H
#define _ICM_Pricer_SECURITY_HOMOGENEOUS_H

#include "ICMKernel/pricer/icm_pricer_basket.h"

class ICM_LossUnits; 
class ICM_DefaultCurve; 
class ICM_ModelMultiCurves ;

/*********************************************************************************/
/*! \class  ICM_Pricer_Distrib icm_pricer_homogeneous.h "icm_pricer_homogeneous.h"
 *  \author Damien Pouponneau 
 *	\version 1.0
 *	\date   June 2003
 *	\brief  ce pricer est basé sur l'utilisation d'un cache de distribution */
/***********************************************************************************/

class ICM_Pricer_Distrib : public ICM_Pricer_Basket
{

private:
	
	ICM_LossUnits*	itsLossUnits; 
	double			itsMinLossUnit;			//Minimal loss Unit user defined

	// 17783 		int				itsTenorShift;
	// 17783 		int				itsIssuerShift;
	double			itsTargetSpreadForFlatCorrel ;
	// 17783  bool			itsIsActivateShift;

protected:

	// 17783 bool		itsFirstComputationForHedgeFlg;
	// 17783 vector<ICM_QMatrix<double>*>	LossTShifts;

	
public :
	ICM_Pricer_Distrib(void) { Init();}

	
	void Init() ; 


	~ICM_Pricer_Distrib(void) ;

	void Set(ARM_Security *sec, ARM_Model *mod,const ICM_Parameters&params,const ARM_Date&asof);

	double AvgRecovery(const std::vector<std::string>&labels);


	virtual double ExpectedLossTranche(const double& yearterm,vector<double>& losses);
// 17783 	virtual double ExpectedLossTrancheForFastSpreadHedge(const double& yearterm);

	//Expected Amortissement sur la tranche : uniquement sur les tranches [x,100]
	double ExpectedAmortTranche(const double& yearterm);

	virtual void Reset(void); 


	virtual double GetTranche(const ARM_Date& date);
	virtual double GetTranche_Up(const ARM_Date& date);
	virtual double GetTranche_Down(const ARM_Date& date);
	virtual double ComputeFlatCorrel();
	double FlatCorrelEvaluate(double FlatCorrel) ;
	protected:
	virtual double ComputeSensitivity(const qSENSITIVITY_TYPE& typesensi, 
		const std::string& plot,
		const std::string& label,
		double epsilon, double epsilonGamma =0 ); 
	public:

// 17783 	void SetTenorShift(const int& nbshift) { itsTenorShift=nbshift;}
// 17783 			int GetTenorShift(void) { return itsTenorShift;}

// 17783 		void SetIssuerShift(const int& nbshift) { itsIssuerShift=nbshift;}
// 17783 			int GetIssuerShift(void) { return itsIssuerShift;}
	
	void SetTargetSpreadForFlatCorrel( const double& spread) { itsTargetSpreadForFlatCorrel = spread;}
	double GetTargetSpreadForFlatCorrel(void) { return itsTargetSpreadForFlatCorrel;}



// 17783 	void ActivateShift(const bool& val) { itsIsActivateShift = val; }
// 17783 	bool IsActivateShift(void) { return itsIsActivateShift; }


// 17783 		ICM_QMatrix<double>* PerturbProbasFromSensiManager(const double& yearterm);


// 17783 	virtual void PerturbDefaultCurves();
// 17783 		virtual void ComputePrice_Basket_Complete(const ENUM_SENSITIVITY& sensitivity_type, 
// 17783 												  const	double& initialprice, 
// 17783 												  const double& epsvalue);

// 17783	double PreCheckBeforeComputePrice(const ENUM_SENSITIVITY& sensitivity_type,
// 17783										const std::string&	label, 
// 17783										const double& initialprice, 
// 17783										const double& epsvalue);

	void View(char* id, FILE* ficOut);


	// ------------------------------------------------------------
	// FOR FAST SPREAD COMPUTATION
	// ------------------------------------------------------------

	// 17783 void ResetPricerForHedge(void)  
	// 17783 {
	// 17783 	ICM_Pricer:: ResetPricer();
	// 17783 itsFirstComputationForHedgeFlg = false;
	// 17783 }

	// ------------------------------------------------------------
	// 17783 bool GetFirstComputationForHedgeFlg(void) { return (itsFirstComputationForHedgeFlg); }
	// 17783 void SetFirstComputationForHedgeFlg(bool value) { itsFirstComputationForHedgeFlg = value; }
	// 17783 void	ClearShifts() {LossTShifts.clear();}
	// 17783 void	AppendShifts(ICM_QMatrix<double>* data) {LossTShifts.push_back(data);}
	// 17783 ICM_QMatrix<double>*	GetShifts(const int& Num) {return LossTShifts[Num];}
	// 17783 const std::vector<ICM_QMatrix<double>*>& GetLossTShifts() const { return LossTShifts; }


	ICM_DefaultCurve* CptDefaultCurve(const ARM_Vector& dates,
									  const ARM_Vector& spreads,
									  ICM_ModelMultiCurves* model,
									  const string& label,
									  const double& recovery,
									  const string& RefCurve = "NONE");

	virtual ICM_DefaultCurve* GenerateImpliedCurve();
	const ICM_LossUnits& getLossUnits()   
	{
		if (!itsLossUnits) itsLossUnits = CptLossUnits(); 
		return *itsLossUnits; 
	}
protected:
	void ResetLossUnit(); 

public:

	virtual ICM_Pricer* CloneOnlyParams(ARM_Security* sec,ARM_Model* mod)
	{
	ICM_Pricer_Distrib* theClone = new ICM_Pricer_Distrib(); 
	theClone->Set(sec,mod,ICM_Parameters(),GetAsOfDate());
    return(theClone);
	}
private:
	ICM_LossUnits*  CptLossUnits() ; 

};



#endif

