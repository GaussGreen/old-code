
#if !defined(_ICM_CREDIT_INDEX_H_)
#define _ICM_CREDIT_INDEX_H_

#include "ARMKernel\inst\irindex.h"
#include "ICMKernel\glob\icm_enums.h"

/*********************************************************************************/
/*! \class  ICM_Credit_Index icm_credit_index.h "icm_credit_index.h"
 *  \author 
 *	\version 1.0
 *	\date   June 2004
 *	\brief  Defines Credit Index */
/***********************************************************************************/

class ICM_Cds;
class ICM_ModelMultiCurves ;
class ICM_DefaultCurve ;

class ICM_Leg; 
class ICM_Credit_Index : public ARM_IRIndex  
{
	friend void CptCouponsForFeeLeg(ICM_Leg& leg, ARM_Object* object ,const ARM_Date& AsOf) ;
private: 
	std::string             itsCreditIndexName;
	std::vector<std::string>	itsLabels;			// optionnel. adopted not owner if no compostion get the CreditIndexName
	qINDEX_CMPT_METHOD		itsMethod;				//Computation Method for Implied default curve 	
	ICM_DefaultCurve*		itsDefaultCurve;		//	0..1, aggreg. Implied default curve
	ICM_DefaultCurve*		itsForcedCurve;			//	0..1, aggreg. Forced Default Curve 
	ARM_Vector				itsMaturities;			// vector of maturities(YearTerms or Dates) corresponding to coupons
	bool					itsIsYT;					
	ARM_Vector				itsRunning ;			//Market Spread Quotation for index fixed evrey 6 months
	
//	ARM_VolCurve*			itsVol;					//	0-1, aggreg. Index volatility
	qCDS_ADJ				itsAdjForTenor;			//Tenor Adjustment (5 ->AIMM 20/../..)
	bool					itsDefaultCurveflg;
	int						itsCM_resetWeekDay;		//	0(Sunday)..6(Saturday). Default=5(friday)
	int						itsCM_resetOccur;		//	Default=2. "second friday following the first day of the relevant interest period".
	bool					itsIsHomogeneous;		//for full homogeneous case


public:


	ICM_Credit_Index()	{Init();}	

	void Init() ;

    ICM_Credit_Index(const int& dayCount, 
					 const int& resetFreq, 
					 const int& payFreq, 
					 const ARM_Vector& maturity,
					 const std::string &ccy /*ARM_Currency* ccyName = ARM_DEFAULT_CURRENCY*/,
					 const std::string& IndexName ,
					 const std::vector<std::string>& labels/** = NULL**/ ,
					 const qINDEX_CMPT_METHOD& Method /** = qAVERAGE**/ ,
					 const ARM_Vector& Spread /** = 0.**/ ,
					 const ICM_DefaultCurve* ForcedCurve /** = NULL**/ ,
					 //const ARM_VolCurve* Volatility /** = NULL**/ ,
					 const int& fwdRule/** =0**/ ,
					 const int& intRule/** =0**/ ,	
					 const int& resetTiming/** =K_ARREARS**/ ,	
					 const int& resetGap/** =0**/ ,	
					 const int& payTiming/** =K_ARREARS**/ ,
					 const int& payGap/** =0**/ ,
					 const qCDS_ADJ& adj /** = qCredit_Adjust20**/ ,
					 int cm_resetWeekDay/** =5**/ ,
					 int cm_resetOccur/** =2 **/ );


    void Set(const int& dayCount, 
			 const int& resetFreq, 
			 const int& payFreq, 
			 const ARM_Vector& maturity,
			 const std::string & ccy  /*ARM_Currency* ccyName  = ARM_DEFAULT_CURRENCY*/,
			 const std::string& IndexName ,
			 const std::vector<std::string>& labels /** = NULL**/ ,
			 const qINDEX_CMPT_METHOD& Method /** = qAVERAGE**/ ,
			 const ARM_Vector& Spread/** =0.**/ ,
			 const ICM_DefaultCurve* ForcedCurve /** = NULL**/ ,
			 //const ARM_VolCurve* Volatility /** = NULL**/ ,
			 const int& fwdRule/** =0**/ ,
			 const int& intRule/** =0**/ ,	
			 const int& resetTiming/** =K_ARREARS**/ ,	
			 const int& resetGap/** =0**/ ,	
			 const int& payTiming/** =K_ARREARS**/ ,
			 const int& payGap/** =0**/ ,
			 const qCDS_ADJ& adj /** = qCredit_Adjust20**/ ,
			 int cm_resetWeekDay /** =5**/ ,
			 int cm_resetOccur/** =2**/ );

	virtual ~ICM_Credit_Index() ;

private:
	void BitwiseCopy(const ARM_Object* src) ;
public:
	virtual void Copy(const ARM_Object* src) ;
	virtual ARM_Object* Clone(void) ;

	inline const std::string& GetLabel(unsigned int i) const ;
	const std::vector<std::string>& GetLabels() const {return (itsLabels);}

	int GetNbCurves() const {return itsLabels.size();}
	const ICM_DefaultCurve* GetDefCurve() const { return itsDefaultCurve;}
	void adoptDefCurve(ICM_DefaultCurve* value) ;

	const ICM_DefaultCurve* GetForcedCurve() const { return itsForcedCurve;}
	void SetForcedCurve(const ICM_DefaultCurve* value) ;

//	const ARM_VolCurve* GetVolatility() const {return itsVol;}
//	void SetVolatility(const ARM_VolCurve* value) ;
 
// private:
	//double FwdSpread(const double& yt1, const double& yt2 = -1.) const ;
	//double AjustConvexity(const double& yt1, const double& yt2 = -1., const double& FwdSpread = -1.)const ;
	//double PaymentLagAjust(const double& yt1, const double& yt2 = -1., const double& FwdSpread = -1.)const ; // Payment Lag for CMCDS
	const ICM_DefaultCurve* CptDefaultCurve(ICM_ModelMultiCurves* model);
public:
	// JLA useless void CptDefaultCurveStr(ICM_ModelMultiCurves* model);
	const ICM_DefaultCurve* CptDefaultCurveIndex(ICM_ModelMultiCurves* model, ICM_Cds* cds);

	void SetRunning(const ARM_Vector& Spread);
	const ARM_Vector& GetRunning(void)const { return itsRunning;}
	double GetSpread(void)const { if(itsRunning.size() >0) return itsRunning[0];
									else return 0;}
		bool GetIsYT() const { return itsIsYT;}
	void SetMaturities(const ARM_Vector& Maturities);
	const ARM_Vector& GetMaturities(void)const { return itsMaturities;}


	qINDEX_CMPT_METHOD GetMethod(void) const { return itsMethod ;}

	// virtual void ExcludeIssuer(const std::string&label) ;

	void View(char* id, FILE* ficOut);

	void ResetDefCurve() ;

	qCDS_ADJ GetAdjForTenor() const {return itsAdjForTenor;}
	void SetAdjForTenor(const qCDS_ADJ& value) {itsAdjForTenor=value;}
	bool GetDefaultCurveflg() const { return itsDefaultCurveflg; }
	int GetCM_resetWeekDay() const { return itsCM_resetWeekDay; }
	int GetCM_resetOccur() const { return itsCM_resetOccur; }	
private:
	void CptFullHomog(void); 
public:
	bool IsHomogeneous(void) const {return itsIsHomogeneous;}
public:
	const std::string &GetCreditIndexName() const { return itsCreditIndexName;}
	void SetCreditIndexName(const std::string& aCreditIndexName) { itsCreditIndexName = aCreditIndexName;}
private:
	ICM_Credit_Index(const ICM_Credit_Index&); //NA
	ICM_Credit_Index& operator=(const ICM_Credit_Index&); // NA
};

//	---------------------------------------------------------------------------------
inline const std::string&
ICM_Credit_Index::GetLabel(unsigned int i) const
{
	if (i>=itsLabels.size()) 
		ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_Credit_Index::GetLabel: out of bound "<<i); 
	return itsLabels[i]; 
}
#endif 
