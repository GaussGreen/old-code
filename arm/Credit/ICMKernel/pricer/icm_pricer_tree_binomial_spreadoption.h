
#ifndef _ICM_PRICER_TREE_BINOMIAL_SPDOPT_H
#define _ICM_PRICER_TREE_BINOMIAL_SPDOPT_H

// #include "ICMKernel\pricer\icm_pricer_tree.h"
#include "ICMKernel\pricer\icm_pricer.h"
#include "ICMKernel\util\icm_utils.h"
#include "ICMKernel\util\icm_registry.h"
#include "ICMKernel\mod\icm_treeservices.h" 
#include "ICMKernel\inst\icm_cds.h"
class ICM_SpreadOption; 

class ICM_Pricer_Tree_Binomial_SpreadOption : public ICM_Pricer
{
	typedef ICM_Pricer parent; 
	typedef ICM_Registry<ICM_Cds,ARM_Date> cdsreg_t; 
private:
	StdBinomialTree itsTree ;		
	cdsreg_t itsCdsReg		;		// 	underlyings CDS are stored by start Date
									//	and are aggregated.
	double				itsDelta;
	double				itsGamma;
	double				itsVega;
	double				itsTheta;
	double				itsRho;

	bool				itsDeltaFlg ;
	bool				itsGammaFlg ;
	bool				itsVegaFlg ;
	bool				itsThetaFlg ;
	bool				itsRhoFlg ;

public:
	static std::string PARAM_FORCEVOL ;
	static std::string PARAM_OPTIONDISCSTEPSIZE; 
	static std::string PARAM_OPTIONDISCSTEPNUMBER; 
	static std::string PARAM_USECONSTANTPROBA; 
	static std::string PARAM_DISCTOLERANCE; 
	static std::string PARAM_DISCMAXITER ; 
	static std::string PARAM_DISCVARIANCETOL ; 
public :
	ICM_Pricer_Tree_Binomial_SpreadOption() ;
	ICM_Pricer_Tree_Binomial_SpreadOption(ARM_Security *sec, ARM_Object *mod,const ICM_Parameters& parameters,const ARM_Date&asof) ;
	virtual ~ICM_Pricer_Tree_Binomial_SpreadOption() ;
	void Init() ; 

	// non virtual 
	void Set(ARM_Security *option, ARM_Object *mod,const ICM_Parameters& parameters,const ARM_Date&asof) ; 

	//	--	ARM_Object
	virtual void View(char* id, FILE* ficOut) ;
	

	//	--	ICM_Pricer	--	Pricing functions entry point 
	virtual void Reset() ; 
	virtual double ComputePrice(qCMPMETH measure ) ;


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
		double epsilon ,  double epsilonGamma = 0); 
	public:

		
	virtual void SetModel(ARM_Model* model){ ICMTHROW(ERR_INVALID_OPERATION,"Not Implemented"); }
	virtual void computelossunit(){ ICMTHROW(ERR_INVALID_OPERATION,"Not Implemented"); }
// 17783 	virtual void PerturbDefaultCurves(){ ICMTHROW(ERR_INVALID_OPERATION,"Not Implemented"); }

	virtual	double GenericLegPrice(const int& mode,double* Data);
	void CpteSpreadsAndPv01s(ICM_QMatrix<double>&ret, double mode) ; 
	void CpteSpreadsAndPv01s(StdBinomialTree& ret) ;
	ICM_Cds&	getStdCDS(const ARM_Date& startDate,const ARM_Date&endDate);

	void SetDelta(const double& Delta)
	{
		itsDelta = Delta;
		itsDeltaFlg = true;
	}

	void SetGamma(const double& Gamma)
	{
		itsGamma = Gamma ;
		itsGammaFlg = true ;
	}

	void SetVega (const double&  Vega)
	{
		itsVega  = Vega;
		itsVegaFlg = true ;
	}

	void SetTheta(const double& Theta)
	{
		itsTheta = Theta;
		itsThetaFlg = true ;
	}

	void SetRho(const double& Rho)
	{
		itsRho   = Rho;
		itsRhoFlg = true;	
	}
	
	double GetDelta()
	{
		if(!itsDeltaFlg)
			ComputePrice(qCMPPRICE/*GetAsOfDate(),7,0*/);
		return itsDelta ;
	}

	double GetGamma()
	{
		if(!itsGammaFlg)
			ComputePrice(qCMPPRICE/*GetAsOfDate(),7,0*/);
		return itsGamma ;
	}
	
	double GetVega()
	{
		if(!itsVegaFlg)
		{
			double sensi = ComputeSensitivity(ICM_GREEK_VEGA_TYPE,"NONE","NONE",0.01);
			SetVega (sensi);
		}
		return itsVega ;
	}

	double GetTheta()
	{
		if(!itsThetaFlg)
			ComputePrice(qCMPPRICE/*GetAsOfDate(),7,0*/);
		return itsTheta ;
	}

	double GetRho()
	{
		if(!itsRhoFlg)
		{
			double sensi = ComputeSensitivity(ICMIRCURVE_TYPE,"NONE","EUR",0.01);
			SetRho (sensi);
		}
		return itsRho ;
	}

	bool GetRhoFlg() {return itsRhoFlg ;}
	bool GetVegaFlg() {return itsVegaFlg ;}

	virtual double ComputeGreeks (const int& Type);
	// virtual void ComputeGreeksData(ICM_QMatrix<double>& ret);


public:
	ICM_DefaultCurve&	getDefCurve() ; 
	ICM_SpreadOption&	getSpreadOption() ; 
	ARM_ZeroCurve&		getDiscountCurve() ;
protected: 
	double	EuropeanPrice() ; 
	double	AmericanPrice();
	double	FwdSpreadPrice();
private:
	//	--	NOT IMPLEMENTED 
	//
	ICM_Pricer_Tree_Binomial_SpreadOption(const ICM_Pricer_Tree_Binomial_SpreadOption&) ;			
	ICM_Pricer_Tree_Binomial_SpreadOption& operator=(const ICM_Pricer_Tree_Binomial_SpreadOption&); 
private: 
	virtual double ComputeSpread(const double& MtM = 0.)  ; 
	virtual double Accrued() ; 
	virtual double FeeLegPV () ; 
	virtual double DefLegPV () ; 
	virtual double ComputeDuration(void) ; 
	virtual	double ComputeImpliedVol(const double& Price) ; 
	virtual double Compute_Fwd_Spread(const ARM_Date &Mty, const ARM_Date &CDS_ExpiryDate,double& dur) ; 

};
#endif	// _ICM_PRICER_TREE_BINOMIAL_SPDOPT_H

