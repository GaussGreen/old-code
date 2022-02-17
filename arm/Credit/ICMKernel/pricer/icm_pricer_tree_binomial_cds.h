/********************************************************************************/
/*! \file ICM_Pricer_Tree_Binomial_Cds.h
 *  \author D. Pouponneau
 *	\version 1.0
 *	\date   Febuary 2005 */
/*
 *********************************************************************************/

#ifndef _ICM_PRICER_TREE_BINOMIAL_CDS_H
#define _ICM_PRICER_TREE_BINOMIAL_CDS_H

// #include "ICMKernel\pricer\icm_pricer_tree_cds.h"
#include "ICMKernel\pricer\icm_pricer.h"
// #include "ICMKernel\inst\icm_cds.h"
// #include "ICMKernel\util\icm_utils.h"

#include "ICMKernel\mod\icm_treeservices.h" 


class ICM_Pricer_Tree_Binomial_Cds : public ICM_Pricer
{
	typedef ICM_Pricer parent; 
private:
	StdBinomialTree itsTree ;	
	// long itsUseOldDefLegPricing; 
public:
	static std::string PARAM_FORCEVOL ;
	static std::string PARAM_USEDEFCURVEPOINTS; 
	static std::string PARAM_OPTIONDISCSTEPSIZE ;
	static std::string PARAM_OPTIONDISCSTEPNUMBER ;
	static std::string PARAM_USECONSTANTPROBA;
	// static std::string PARAM_USEOLDDEFLEGPRICING ;
	static std::string PARAM_DISCTOLERANCE ;
	static std::string PARAM_DISCMAXITER; 
	static std::string PARAM_DISCVARIANCETOL; 

public :
	ICM_Pricer_Tree_Binomial_Cds() ;
	ICM_Pricer_Tree_Binomial_Cds(ARM_Security *sec, ARM_Object *mod,const ICM_Parameters&parameters,const ARM_Date&asof) ;
	virtual ~ICM_Pricer_Tree_Binomial_Cds() ;
	void Init() ; 

	// non virtual 
	void Set(ARM_Security *option, ARM_Object *mod,const ICM_Parameters& parameters,const ARM_Date&asof) ; 

	//	--	ARM_Object
	virtual void View(char* id, FILE* ficOut) ;
	

	//	--	ICM_Pricer	--	Pricing functions entry point 
	virtual void Reset() ; 
	virtual double ComputeSpread(const double& MtM = 0.)  ;
	virtual double Accrued() ;
	virtual double FeeLegPV () ;
	virtual double DefLegPV () ;

	virtual double ComputePrice(qCMPMETH measure) ;
	virtual void MarketData(ARM_Security* sec,vector<string>& DMKT) ;

	//	--	UNIMPLEMENTED 
	//	--- ARM_Object
	virtual void Copy(const ARM_Object* src) { ICMTHROW(ERR_INVALID_OPERATION,"Not Implemented"); }
	virtual ARM_Object* Clone() { ICMTHROW(ERR_INVALID_OPERATION,"Not Implemented"); }
    virtual ARM_Object& operator = (const ARM_Object& obj){ ICMTHROW(ERR_INVALID_OPERATION,"Not Implemented"); }
	virtual ARM_CLASS_NAME GetRootName() ; 

	//	--	UNIMPLEMENTED 
	//	---	ICM_Pricer	--	
		protected:
	virtual double ComputeSensitivity(const qSENSITIVITY_TYPE& typesensi, 
		const std::string& plot,
		const std::string& label,
		double epsilon, double epsilonGamma =0 ); 
	public:

		
private:
	virtual double ComputeDuration(void) ; 
	virtual void SetModel(ARM_Model* model); 
	virtual void computelossunit() ;
// 17783 	virtual void PerturbDefaultCurves() ;
	virtual	double ComputeImpliedVol(const double& Price) ;
	virtual double Compute_Fwd_Spread(const ARM_Date &Mty, const ARM_Date & CDS_ExpiryDate, double& dur) ;


private:
	//	--	NOT IMPLEMENTED 
	//
	ICM_Pricer_Tree_Binomial_Cds(const ICM_Pricer_Tree_Binomial_Cds&) ;			
	ICM_Pricer_Tree_Binomial_Cds& operator=(const ICM_Pricer_Tree_Binomial_Cds&); 
};

#endif

