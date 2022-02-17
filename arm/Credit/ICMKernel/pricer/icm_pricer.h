
#ifndef _ICM_PRICER_H
#define _ICM_PRICER_H

#include <string>
#include <algorithm>

class ARM_Model; 
class ARM_Security; 
#include "ICMKernel\util\icm_matrix.h"
// #include "ICMKernel\util\icm_SensiManager.h"




/*********************************************************************************/
/*! \class  ICM_Pricer icm_pricer.h "icm_pricer.h"
 *  \author Damien Pouponneau 
 *	\version 1.0
 *	\date   June 2003
 *	\brief  Root Pricer */
/***********************************************************************************/

class ICM_MktDataMng ;
// class ICM_SensiManager ; 
class ICM_Pricer : public ARM_Object
{

private:
	ARM_Date		itsAsOfDate;    //As Of Date for valuation

	ARM_Security*	itsSecurity;	//	association -	Default Security
	ARM_Object*		itsModel;		//	association -	Default Model 
	ICM_Parameters	itsParameters; 

	double			itsInitialPrice;//Inital price (for sensitivity computing)
	bool			itsInitialPriceFlg;//Is Initial price already computed ?


	double			itsFeeLegPrice_Unity;
	bool			itsFeeLegPrice_UnityFlg;

	bool			itsFaster;		//Flag used to accelerate pricing when no view is needed
	
	//
	//	---------------------------------------
	//		This is the pricing measure cache 
	bool			itsFlags [ qCMPLAST ] ;
	double			itsValues [qCMPLAST ] ; 
	//
	//	---------------------------------------
	//		This is the Vector measur cache 
	ARM_Vector*		itsObjectValues[qVECTLAST];	
	//
	//	---------------------------------------
	//		This is the arbitrary pricing measure cache 
	std::vector<std::string>	itsMeasureNames ; 
	std::vector<double>			itsMeasureValues; 
	//	---------------------------------------
	//		This is the arbitrary vector measure cache 
	std::vector<std::string>	itsVectorMeasureNames ; 
	std::vector<ARM_Vector*>	itsVectorMeasureValues; 
protected:
	//17783  ICM_SensiManager* itsSensiManager;
	ICM_Pricer(void) { Init();}
	void Init() ;
	void Set(ARM_Security*sec,ARM_Object*mod,const ICM_Parameters&params,const ARM_Date*asof) ;
public: 
	//	Generic API. 
	double Price(qCMPMETH measure) ; 	
	double Price(const std::string&measurename); 
	const ARM_Vector& PriceVector(qVECTMETH measure); 
	const ARM_Vector& PriceVector(const std::string& measure); 
	// to be removed: use Price instead
	virtual double ComputePrice(qCMPMETH mode) ; 
	//	Computation function: to become protected. 
protected:
	virtual double DoPrice(qCMPMETH measure) ;
	virtual double DoPrice(const std::string& measure); 
	virtual void DoPriceVector(qVECTMETH measure) ;
	virtual void DoPriceVector(const std::string& measure) ;
public:
	virtual double ComputeSpread(const double& MtM = 0.) =0; 
 	virtual	double ComputeImpliedVol(const double& Price) =0 ;
	//virtual double Compute_Fwd_Spread(const ARM_Date &Mty, const ARM_Date &CDS_ExpiryDate, double& dur) =0; 
	// 
	void SetPriceFollowMode(qCMPMETH measure,const double& value);
private:
	// Computation Functions are not accessible to clients.
	// 
	virtual double Accrued() =0; 
	virtual double FeeLegPV () =0 ;
	virtual double DefLegPV () =0; 
	// cache management. 
protected:
	double CheckForPrice(qCMPMETH measure);
	// 
protected:
	virtual double ComputeSensitivity(const qSENSITIVITY_TYPE& typesensi, 
		const std::string& plot,
		const std::string& label,
		double epsilon, double epsvalueGamma = 0) =0; 
		
	
public: // wip
	bool getFlg(qCMPMETH measure)	const { return itsFlags[measure]; }
	void setFlg(qCMPMETH measure)	{ itsFlags[measure]=true; }
	void unsetFlg(qCMPMETH measure) { itsFlags[measure]=false; }
	void unsetFlgs() ;  
protected:
	bool getFlg(const std::string& measure)	const 
	{ 
	 	std::vector<std::string>::const_iterator it = std::find(itsMeasureNames.begin(),itsMeasureNames.end(),measure) ;
	 	if (it==itsMeasureNames.end()) return false; 
	 	return true; 
	}
	void unsetFlg(const std::string& measure) 
	{ 
		std::vector<std::string>::iterator it = std::find(itsMeasureNames.begin(),itsMeasureNames.end(),measure) ;
		if (it==itsMeasureNames.end()) return ; // already not set
		std::vector<double>::iterator it2 = itsMeasureValues.begin() + (it-itsMeasureNames.begin()) ;
		itsMeasureNames.erase(it);
		itsMeasureValues.erase(it2); 
	}

	bool getFlgVector(qVECTMETH measure)  const	{ return itsObjectValues[measure]!=0; }
	void unsetFlgVector(qVECTMETH  measure)		
	{ 
		if (itsObjectValues[measure]) delete itsObjectValues[measure]; 
		itsObjectValues[measure]=0; 
	}
	bool getFlgVector(const std::string& measure)
	{
	 	std::vector<std::string>::const_iterator it = std::find(itsVectorMeasureNames.begin(),itsVectorMeasureNames.end(),measure) ;
	 	if (it==itsVectorMeasureNames.end()) return false; 
	 	return true; 
	}
	void unsetValueFlg(const std::string& measure)
	{
		std::vector<std::string>::iterator it = std::find(itsMeasureNames.begin(),itsMeasureNames.end(),measure) ;
		if (it==itsMeasureNames.end()) return ; // already not set
		std::vector<double>::iterator it2 = itsMeasureValues.begin() + (it-itsMeasureNames.begin()) ;
		itsMeasureNames.erase(it);
		itsMeasureValues.erase(it2); 
	}
	void unsetVectorFlg(const std::string& measure)
	{
		std::vector<std::string>::iterator it = std::find(itsVectorMeasureNames.begin(),itsVectorMeasureNames.end(),measure) ;
		if (it==itsVectorMeasureNames.end()) return ; // already not set
		std::vector<ARM_Vector*>::iterator it2 = itsVectorMeasureValues.begin() + (it-itsVectorMeasureNames.begin()) ;
		itsVectorMeasureNames.erase(it);
		delete *it2; 
		itsVectorMeasureValues.erase(it2); 
	}
public:
	// not yet analyzed... 
	inline void SetDuration(const double& Duration) {
		setValue(qCMPDURATION,Duration); 
	}
	inline bool GetDurationFlg(void) 
	{	
		// return (itsDurationFlg);
		return getFlg(qCMPDURATION); 
	}

public:

	// public API for computing Sensitivity
	double Hedge(qSENSITIVITY_TYPE  curvetype,
		 const std::string& plot,
		 const std::string& label,
		 double epsilon, double epsilonGamma = 0) 
    {
		return ComputeSensitivity(curvetype, plot, label, epsilon, epsilonGamma);
    } 
	

    virtual ~ICM_Pricer(void) ; 

	void SetPrice(const double& price)  ;

public:
	void ResetRootPricer() ; 
public:
	virtual void Reset(void) {} //Reset Intermediate Data for each pricer
	virtual void RefreshParameters(void) {} //Reset Intermediate Data for each pricer
public:
	void ResetPricer(void)  ;

	
	inline bool GetPriceFlg(void) { return getFlg(qCMPPRICE) ; }
	


public:
	inline bool GetAccruedFlg(void) 
	{ 
		return getFlg(qCMPACCRUED); 
		
	}
	inline void SetAccrued(const double& value) 
	{
		setValue(qCMPACCRUED,value); 
	}
public:
	inline double GetInitialPrice(void) { return (itsInitialPrice); }
	inline void SetInitialPrice(const double& value) { itsInitialPrice = value;}
	inline bool GetInitialPriceFlg(void) { return (itsInitialPriceFlg); }
	inline void SetInitialPriceFlg(const bool& value) 
	{
		itsInitialPriceFlg = value; 
		itsInitialPrice = CREDIT_DEFAULT_VALUE; 	
	}

public:
	inline void SetFeeLegPrice(const double& price)  
	{ 
		setValue(qCMPFEELEGPV,price); 
	}

	inline bool GetFeeLegPriceFlg(void) 
	{ 
		return getFlg(qCMPFEELEGPV); 
	}

	inline void SetDefLegPrice(const double& price)  
	{ 
		setValue(qCMPDEFLEGPV,price); 
	}

	inline bool GetDefLegPriceFlg(void) 
	{ 
		return getFlg(qCMPDEFLEGPV); 

	}
public:
	inline void SetFeeLegPrice_Unity(const double& price)  
	{ 
		itsFeeLegPrice_Unity = price; 
		itsFeeLegPrice_UnityFlg = true; 
	}

	inline double GetFeeLegPrice_Unity(void) { return (itsFeeLegPrice_Unity); }
	inline bool GetFeeLegPrice_UnityFlg(void) { return (itsFeeLegPrice_UnityFlg) ;}

public:
	inline void SetSpread(const double& spread)  
	{ 
		setValue(qCMPSPREAD,spread); 
	}

	inline bool GetSpreadFlg(void) 
	{
		return getFlg(qCMPSPREAD); 
	}
public:
	// not defined in ARM
	virtual void SetModel(ARM_Model* model) ;

	virtual void PropagateModel(ARM_Model *model) {};

	void SetMktDataMng(ICM_MktDataMng* model) ;

	inline ARM_Model* GetModel(void) {return (ARM_Model*)itsModel;}
	inline ICM_MktDataMng* GetMktDataMng(void) {return (ICM_MktDataMng*)itsModel;}

	inline void SetSecurity(ARM_Security* sec) {itsSecurity = sec;}
	inline ARM_Security* GetSecurity(void) {return itsSecurity;}

	void View(char* id = NULL, FILE* ficOut = NULL) ;

	virtual ARM_CLASS_NAME GetRootName(void)
	{
		return(ICM_PRICER);
	}




	virtual ICM_Pricer* CloneOnlyParams(ARM_Security* sec,ARM_Model* mod)  
	{
		ICMTHROW(ERR_INVALID_ARGUMENT,"CloneOnlyParams unimplemented"); 
	} 

	// --------------
	//	Clone Method
	// --------------

	virtual void MarketData(ARM_Security* sec,vector<string>& DMKT){}


	virtual void computelossunit() {};



	inline void SetFaster(const bool& value) { itsFaster = value;};
	inline bool GetFaster(void) { return itsFaster;};

	// 17783 inline ICM_SensiManager* GetSensiManager() {return itsSensiManager;}
	// JLA Added
	// 17783 void SetSensiManager(ICM_SensiManager*item) ; 
	// 17783 virtual void PerturbDefaultCurves()
	// 17783 {
	// 17783 	ICMTHROW( ERR_UNIMP_METHOD_CALL,
	// 17783 		"Unimplemented <PerturbDefaultCurves> method");
	// 17783 }


	
	const ICM_Parameters& GetParameters() const { return itsParameters; }

	void SetParameters(const ICM_Parameters&params) 
	{ 
		itsParameters=params; 
	}
	protected:
		inline void SetAsOfDate(const ARM_Date& date) {itsAsOfDate = date;}
	public:
	inline const ARM_Date& GetAsOfDate() const {return itsAsOfDate;}

	//	JLA - Accessors for the pricer parameters
	//		returns false on 
	bool getParam(const std::string&paramName,double&ret,bool throwOnError=true) const; 
	bool getParam(const std::string&paramName,long&ret,bool throwOnError=true) const; 
	bool getParam(const std::string&paramName,ARM_Vector&ret,bool throwOnError=true) const; 


	// LJ - Get Data
	virtual	void	GetDataFromLabel(string TheLabel, double& TheValue) {TheValue = 0.0;}
	virtual	void	GetDataFromLabel(string TheLabel, string& TheValue) {TheValue = "";}

	virtual void GenerateRates(ARM_Vector& rates) //For SwapLeg pricers
	{
		ICMTHROW( ERR_UNIMP_METHOD_CALL,
			"Unimplemented <GenerateRates> method");
	}
protected:
	void setValue(qCMPMETH measure,double value) 
	{
		itsValues[measure]=value; 
		setFlg(measure); 
	}
	double getValue(qCMPMETH measure) const
	{
		return itsValues[measure]; 
	}
	void setValue(const std::string&measure,double value)
	{
		std::vector<std::string>::iterator it = std::find(itsMeasureNames.begin(),itsMeasureNames.end(),measure) ;
		if (it==itsMeasureNames.end()) 
		{
			itsMeasureNames.push_back(measure); 
			itsMeasureValues.push_back(value); 
		}
		else itsMeasureValues[it-itsMeasureNames.begin()]=value; 
	}
	double getValue(const std::string& measure) const 
	{
		std::vector<std::string>::const_iterator it = std::find(itsMeasureNames.begin(),itsMeasureNames.end(),measure) ;
		if (it==itsMeasureNames.end()) 
			ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_Pricer::getValue: can't find "<<measure); 
		return itsMeasureValues[it-itsMeasureNames.begin()] ;
	}
	const ARM_Vector& getVectorValue(qVECTMETH measure) const 
	{
		if (itsObjectValues[measure]) return *itsObjectValues[measure]; 
		ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_Pricer::getVectorValue: can't get null Object "<<measure); 
	}
	ARM_Vector& getVectorValue(qVECTMETH measure) 
	{
		if (itsObjectValues[measure]) return *itsObjectValues[measure]; 
		ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_Pricer::getVectorValue: can't get null Object "<<measure); 
	}
	const ARM_Vector& getVectorValue(const std::string& measure) const 
	{
		std::vector<std::string>::const_iterator it = std::find(itsVectorMeasureNames.begin(),itsVectorMeasureNames.end(),measure); 
		if (it==itsVectorMeasureNames.end()) 
			ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_Pricer::getVectorValue: can't get null Object "<<measure); 
		return * itsVectorMeasureValues[  it-itsVectorMeasureNames.begin()]  ;
	}
	ARM_Vector& getVectorValue(const std::string& measure) 
	{
		std::vector<std::string>::const_iterator it = std::find(itsVectorMeasureNames.begin(),itsVectorMeasureNames.end(),measure); 
		if (it==itsVectorMeasureNames.end()) 
			ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_Pricer::getVectorValue: can't get null Object "<<measure); 
		return * itsVectorMeasureValues[  it-itsVectorMeasureNames.begin()]  ;
	}
	void adoptVectorValue(qVECTMETH measure,ARM_Vector*item) 
	{
		if (!item) 
			ICMTHROW(ERR_INVALID_ARGUMENT,"Can't adopt null pointer"); 
		if (itsObjectValues[measure]) delete itsObjectValues[measure]; 
		itsObjectValues[measure]=item ;
	}
	void adoptVectorValue(const std::string&measure,ARM_Vector*item)
	{
		if (!item) 
			ICMTHROW(ERR_INVALID_ARGUMENT,"Can't adopt null pointer"); 
		std::vector<std::string>::iterator it = std::find(itsVectorMeasureNames.begin(),itsVectorMeasureNames.end(),measure) ;
		if (it==itsVectorMeasureNames.end()) 
		{
			itsVectorMeasureNames.push_back(measure); 
			itsVectorMeasureValues.push_back(item); 
		}
		else 
		{
			// if (itsMeasureValues[it-itsMeasureNames.begin()]) 
			delete itsVectorMeasureValues[it-itsVectorMeasureNames.begin()] ; 
			itsVectorMeasureValues[it-itsVectorMeasureNames.begin()]=item; 
		}
	}
	
protected:
	// non generic API. Protected, to be removed. 
	// 
	double GetPrice() const			{ return getValue(qCMPPRICE) ;}
	double GetSpread() const		{ return getValue(qCMPSPREAD);}
	double GetDefLegPrice () const	{ return getValue(qCMPDEFLEGPV);}
	double GetFeeLegPrice() const	{ return getValue(qCMPFEELEGPV); }
	virtual double ComputeDuration() =0;	
private: 
	ARM_Object* Clone(void) ; // NA 
	ICM_Pricer& operator=(const ICM_Pricer&); // NA
	ICM_Pricer(const ICM_Pricer&);  // NA 
	virtual void Copy(const ARM_Object* src) 
	{
		ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_Pricer::Copy unimplemented"); 
	}
friend class ICM_Pricer_Advisor; 
}; 
 

#endif
