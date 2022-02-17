
#if !defined(_ICM_SECURITY_H_)
#define _ICM_SECURITY_H_

#include <map>
#include "ICMKernel\glob\icm_enums.h"
#include "ARMKernel\inst\security.h"
 

/*********************************************************************************/
/*! \class  ICM_Security icm_security.h "icm_security.h"
 *  \author 
 *	\version 1.0
 *	\date   June 2003
 *	\file   icm_security.h
 *	\brief  generic security for credit derivatives */
/***********************************************************************************/


class ICM_Security : public ARM_Security
{
private: 

	ARM_Date itsStartDateNA;		//start date 
	ARM_Date itsEndDateNA ;			//maturity date 

	//period for interests
	ARM_Vector itsAccStartDates;		//accrual start dates
	ARM_Vector itsAccEndDates;			//accrual end dates
	ARM_Vector itsInterestDays;			// number of day of the interest period.

	//period for interests (YF from date AsOf)
	ARM_Vector itsYFAccStartDates;		//	accrual start dates (in Year Fractions)
	ARM_Vector itsYFAccEndDates;		//	accrual end dates (in Year Fractions)	
	ARM_Vector itsYFInterestDays ;		//	year fraction of the interest days (accrual basis)

	ARM_Vector itsFwdStartDates;//  Fwd start dates (deduced from ResetDates)
    ARM_Vector itsFwdEndDates;  //	Fwd end dates (deduced from ResetDates)
    ARM_Vector itsResetDates;	 // Fixing dates (deduced from ResetDates)

	ARM_Vector itsYFFwdStartDates;// start dates and end dates (in Year Fractions)
    ARM_Vector itsYFFwdEndDates;  // are used to compute (in Year Fractions) 
    ARM_Vector itsYFResetDates;   // Fixing dates (in Year Fractions)

	ARM_Vector itsStartRiskDates; // credit start risk date
	ARM_Vector itsYFStartRiskDates; // credit year fraction start risk date
	ARM_Vector itsEndRiskDates; // credit end risk date
	ARM_Vector itsYFEndRiskDates; // credit end fraction start risk date

	ARM_Vector itsPayDates;		//Payment dates
	ARM_Vector itsYFPayDates;	//Payment dates (YF from date AsOf)

	double itsFrozenMaturity;	//Frozen Maturity as a julian date (CREDIT_DEFAULT_VALUE: no frozen maturity)

	// ARM_Currency* itsCcy;		//currency of the security
	std::string itsCcy; 
	int itsAccrualBasis;		//Accrual basis		
	int itsPaymentFreq;			//Payment frequency

	ARM_Vector itsPartRates;	//Participation rate
	ARM_Vector itsCouponRates;	//Risky variable spreads
	ARM_Vector itsFwdSpreads;	//Vector added to Rates Vector
	ARM_Vector itsNotionals;	//Notionals
	ARM_Vector itsNotionalXchange; //Notional exchange

	ARM_Vector itsDefaultProbability; //default probability at period enddate 
	ARM_Vector itsDiscountRate;		//discount yield at period paydate

	ARM_Vector itsPeriodRefDefLegPV; 
	ARM_Vector itsPeriodRecDefLegPV;
	ARM_Vector itsPeriodDefLegPV;
	ARM_Vector itsPeriodFeeLegPV;
	
	ARM_Vector itsAdjFwdSpreads;	// Fwd Spread + Convexity Adj + Payment Lag Adj
	ARM_Vector itsUnAdjFwdSpreads;	// Fwd Spread withaout any adj
	ARM_Vector itsFwdPV01s;			// Fwd PV01s between Reset Date and Paydate (Settled to 1 if no Payment Lag)

	std::map<double,double>	itsPastFixings;		// the first double is a Julian date corresponding to ResetDates.
	
	qSecurity_TYPE	itsSecurityType; // Security Type (Running, Infine, ZeroCoupon...)
	vector<qCOUPON_TYPE> itsFwdsType; // Additive, Exclusive ...

	// bool itsFlgYF; //Flag YearFraction Date AsOf (if already compute)

	qPAYMENT_PREMIUM_LEG  itsAccruedOnDefault; /*!< Reference sur le calcul du recouvrement du coupon */

	//Internal use only
	double itsLastAsOf;		//doesn't compute yearfractions twice for a same date
	bool itsIncludeMaturity;//Is Maturity included ?


public:

	void Init();

	ICM_Security() {Init();	}	

	ICM_Security(const ICM_Security& security) : ARM_Security(security)
	{ Init();
	BitwiseCopy(&security);} 	



protected:
	ICM_Security(const ARM_Date& StartDate,
				const ARM_Date& EndDate ,
				const int& nbflows,
				const ARM_ReferenceValue&premiumNotionals,
				const int& AccrualBasis /* = KACTUAL_360*/ ,
				const std::string& ccy /*= ARM_DEFAULT_COUNTRY*/ )
	{
		Init();

		Set(StartDate,
			EndDate,
		nbflows,
			premiumNotionals,
			AccrualBasis,
			ccy);
	}
protected:
	void Set(const ARM_Date& StartDate,
			const ARM_Date& EndDate,
			const int& nbflows,
			const ARM_ReferenceValue&premiumNotionals,
			const int& AccrualBasis /* = KACTUAL_360*/,
			const std::string& Ccy /* = ARM_DEFAULT_COUNTRY*/);
public:
	ICM_Security(const ARM_Date& StartDate,
				const ARM_Date& EndDate,
				ARM_Vector* AccStartDates,
				ARM_Vector* AccEndDates,
				ARM_Vector* AccrualDates,
				ARM_Vector* PayDates,
				ARM_Vector* CouponRates,
				ARM_Vector* Notionals,
				ARM_Vector* NotionalXchange,
				bool includeMaturity,
				const double& InitialNotional /*=0.*/,
				//const int& CreditGap /*= 0*/,
				const int& AccrualBasis /*= KACTUAL_360*/,
				//const int& Basis /*= KACTUAL_360*/,
				const std::string& Ccy/*= ARM_DEFAULT_COUNTRY*/)
	{
		Init();

		Set(StartDate,
			EndDate,
			AccStartDates,
			AccEndDates,
			AccrualDates,
			PayDates,
			CouponRates,
			Notionals,
			NotionalXchange,
			includeMaturity,
			InitialNotional,
			AccrualBasis,
			Ccy);
	}
private:
	void Set(const ARM_Date& StartDate,
			const ARM_Date& EndDate,
			ARM_Vector* AccStartDates,
			ARM_Vector* AccEndDates,
			ARM_Vector* AccrualDates,
			ARM_Vector* PayDates,
			ARM_Vector* CouponRates,
			ARM_Vector* Notionals,
			ARM_Vector* NotionalXchange,
			bool  includeMaturity,
			const double& InitialNotional /* =0.*/ ,
			const int& AccrualBasis /* = KACTUAL_360*/ ,
			const std::string& Ccy /* = ARM_DEFAULT_COUNTRY */ );
public:

	virtual ~ICM_Security()	{}	

	void BitwiseCopy(const ARM_Object* srcleg);

	void Copy(const ARM_Object* srcleg);

	ARM_Object* Clone(void); 
	
	// JLA 
	virtual ARM_Object* Clone() const ; 

	void View(char* id, FILE* ficOut);

	const ARM_Date& GetStartDateNA() const { return itsStartDateNA;}
	const ARM_Date& GetEndDateNA() const { return itsEndDateNA;}
protected:
	void  SetStartDateNA(const ARM_Date&item)  { itsStartDateNA=item ; }
	void SetEndDateNA(const ARM_Date&item)  {  itsEndDateNA=item ; }
public:
	// do not knwow..
	inline void SetAccrualBasis(int value)	{itsAccrualBasis = value;}
	inline int GetAccrualBasis(void)const {return itsAccrualBasis;}

	void SetAccStartDates(const ARM_Vector& vector) 
	{
		itsAccStartDates = vector;
	}
	const ARM_Vector& GetAccStartDates(void) const 
	{
		return itsAccStartDates;
	}

	void SetAccEndDates(const ARM_Vector& vector)
	{	
		itsAccEndDates = vector;
	}
	const ARM_Vector& GetAccEndDates(void) const 
	{
		return itsAccEndDates;
	}
	const ARM_Vector& GetInterestDays() const 
	{ 
		return itsInterestDays ; 
	}
	void SetYFAccStartDates(const ARM_Vector&vector)
	{
		itsYFAccStartDates = vector;
	}
	const ARM_Vector& GetYFAccStartDates(void)
	{
		return itsYFAccStartDates;
	}
	const ARM_Vector& GetYFInterestDays() const 
	{ 
		return itsYFInterestDays; 
	}
	inline void SetYFAccEndDates(const ARM_Vector& vector)
	{
		itsYFAccEndDates = vector;
	}
	const ARM_Vector& GetYFAccEndDates(void)
	{
		return itsYFAccEndDates;
	}
	inline void SetPayDates(const ARM_Vector& vector)	
	{
		itsPayDates = vector;
	}
	inline const ARM_Vector& GetPayDates(void) const 
	{
		return itsPayDates;
	}
	inline void SetYFPayDates(const ARM_Vector& vector)
	{
		itsYFPayDates = vector;
	}
	inline const ARM_Vector& GetYFPayDates(void)	
	{
		return itsYFPayDates;
	}
	inline void SetCouponRates(const ARM_Vector& vector)
	{
		itsCouponRates = vector;
	}
	inline void SetPartRates(const ARM_Vector& vector)
	{
		itsPartRates = vector;
	}
	virtual void SetCreditSpread(const double& spread)
	{	//Risky variable spreads
		if (itsCouponRates.GetSize()>0)
		{for (int i=0; i<itsCouponRates.GetSize(); i++) {itsCouponRates.Elt(i)=spread;}}}

	inline const ARM_Vector& GetCouponRates(void)
	{
		return itsCouponRates;
	}
	inline const ARM_Vector& GetPartRates(void)
	{
		return itsPartRates;
	}
	virtual void SetPartRate(const double& part)
	{	//Risky variable spreads
		if (itsPartRates.GetSize()>0)
		{for (int i=0; i<itsPartRates.GetSize(); i++) {itsPartRates.Elt(i)=part;}}}

	inline const ARM_Vector& GetFwdSpreads(void) const 
	{ 
		return itsFwdSpreads;
	}
	inline const ARM_Vector& GetAdjFwdSpreads(void) const 
	{ 
		return itsAdjFwdSpreads;
	}
	inline const ARM_Vector& GetUnAdjFwdSpreads(void) const 
	{ 
		return itsUnAdjFwdSpreads;
	}
	inline const ARM_Vector& GetFwdPV01s(void) const 
	{ 
		return itsFwdPV01s;
	}
	inline void SetFwdSpreads(const double& spread)
	{	//Risky variable spreads
		if (itsFwdSpreads.GetSize()>0)
		{
			for (int i=0; i<itsFwdSpreads.GetSize(); i++) {itsFwdSpreads.Elt(i)=spread;}
		}
		updateFwdSpreadsFromPastFixings(); 
	}

	inline void SetFwdSpreads(const ARM_Vector& vector)
	{
		itsFwdSpreads = vector;
		updateFwdSpreadsFromPastFixings(); 
	}
	inline void SetNotionals(const ARM_Vector& vector)
	{
		itsNotionals = vector;
	}
	inline const ARM_Vector& GetNotionals(void) const 
	{	
		return	itsNotionals;
	}
	inline void SetFixedNotional(const double& notional)
	{ 
		if (itsNotionals.GetSize()>0)
		{for (int i=0; i<itsNotionals.GetSize(); i++) {itsNotionals.Elt(i)=notional;}}
	}
	inline void SetNotionalXchange(const ARM_Vector& vector)
	{
		itsNotionalXchange = vector;
	}
	inline const ARM_Vector& GetNotionalXchange(void) const 
	{
		return	itsNotionalXchange;
	}
	void SetCcy(const std::string& ccy) 
	{ 
		itsCcy=ccy; 
		SetCurrencyUnit(&ARM_Currency(itsCcy.c_str()));  
	}
	const std::string& GetCcy(void)	const 
	{
		return itsCcy;	
	}
	public:
	void ComputeYF(const ARM_Date &AsOf);
	int PeriodIndex(const ARM_Date& AsofDate);
	double DailyAccrued(int& i,bool includenot = true);
	double FullCoupon(const ARM_Date& date);
	double FullCoupon(const int& index);
	double FullCouponRate(const int& index);

	double Notional(const ARM_Date& date);
	double NotionalXChange(const ARM_Date& date);

	double AccruedCoupon(const ARM_Date& Asofdate,const bool& includenot =true);

	double DiffAccrued  (const ARM_Date& t1, const ARM_Date& t2,const bool& includenot = true);
	double DiffAccrued (const ARM_Date& AsOf,const double& yf1,const double& yf2,const bool& includenot = true);

	inline void SetAccruedOnDefault(qPAYMENT_PREMIUM_LEG value) {itsAccruedOnDefault = value;}
	inline qPAYMENT_PREMIUM_LEG GetAccruedOnDefault(void) {return itsAccruedOnDefault;}

	inline const ARM_Vector& GetDefaultProbability() const {return itsDefaultProbability;}
	inline const ARM_Vector& GetDiscountRate() const {return itsDiscountRate;}

	inline void SetDefaultProbability(const ARM_Vector&vector) {itsDefaultProbability=vector;}
	inline void SetDiscountRate(const ARM_Vector& vector) {itsDiscountRate=vector;}

	inline const ARM_Vector&GetPeriodRefDefLegPV() const {return itsPeriodRefDefLegPV;}
	inline const ARM_Vector&GetPeriodRecDefLegPV() const {return itsPeriodRecDefLegPV;}

	inline const ARM_Vector&GetPeriodDefLegPV() const {return itsPeriodDefLegPV;}
	inline const ARM_Vector& GetPeriodFeeLegPV() const {return itsPeriodFeeLegPV;}

	inline void SetPeriodRefDefLegPV(const ARM_Vector& vector){itsPeriodRefDefLegPV = vector;}
	inline void SetPeriodRecDefLegPV(const ARM_Vector& vector){itsPeriodRecDefLegPV = vector;}
	inline void SetPeriodDefLegPV(const ARM_Vector& vector){itsPeriodDefLegPV = vector;}
	inline void SetPeriodFeeLegPV(const ARM_Vector& vector){itsPeriodFeeLegPV = vector;}
	inline bool IsLastPeriod(const ARM_Date& AsofDate);   

	inline bool IsIncludeMaturity(void) const { return itsIncludeMaturity;}


	inline qSecurity_TYPE GetSecurityType() { return itsSecurityType;}	
	inline void SetSecurityType(const qSecurity_TYPE& type) { itsSecurityType = type;}

	inline int GetPaymentFreq(void) const { return(itsPaymentFreq);}
    inline void SetPaymentFreq(int payFreq) { itsPaymentFreq = payFreq; }


	inline const ARM_Vector& GetFwdStartDates() const {return itsFwdStartDates;}
	inline void SetFwdStartDates(const ARM_Vector&FwdStartDates){itsFwdStartDates = FwdStartDates;};
    inline const ARM_Vector& GetFwdEndDates() const {return itsFwdEndDates;}
	inline void SetFwdEndDates(const ARM_Vector& FwdEndDates){itsFwdEndDates = FwdEndDates;}
    const ARM_Vector& GetResetDates() const 
	{
		return itsResetDates;
	}
	inline void SetResetDates(const ARM_Vector&ResetDates){itsResetDates=ResetDates;}

	inline const vector<qCOUPON_TYPE>& GetFwdsType(void) const {return itsFwdsType;}
	inline void SetFwdsType(const vector<qCOUPON_TYPE>& type) {itsFwdsType = type;}
	inline void SetFwdType(const qCOUPON_TYPE& type) 
	{ for (int i =0; i<itsFwdsType.size(); i++){itsFwdsType[i]=type;}}

	inline void SetFwdType(const ARM_ReferenceValue* ref) 
	{
		itsFwdsType.resize(itsAccStartDates.size());
		for (int i =0; i<itsFwdsType.size(); i++)
		{ itsFwdsType[i]=(ENUM_COUPON_TYPE)(long)(ref->CptReferenceValue((ARM_Date)itsAccStartDates.Elt(i)));}
	}

	inline const ARM_Vector& GetYFFwdStartDates(void) const {return itsYFFwdStartDates;}
	inline const ARM_Vector& GetYFFwdEndDates(void) const {return itsYFFwdEndDates;}
	inline const ARM_Vector& GetYFResetDates(void) const {return itsYFResetDates;}

	inline void SetYFFwdStartDates(const ARM_Vector& vector) {itsYFFwdStartDates = vector;}
	inline void SetYFFwdEndDates(const ARM_Vector& vector) {itsYFFwdEndDates = vector;}
	inline void SetYFResetDates(const ARM_Vector& vector) 
	{
		itsYFResetDates = vector;
	}
	std::string ConvertFwdType(const qCOUPON_TYPE& type) const 
	{
		string retour;

		switch (type)
		{
		case qCPN_FIXED :
			retour = "Fixed";break;
		case qCPN_FIXED_AND_FWD :
			retour = "Sum";break;
		case qCPN_FWD :
			retour = "FwdOnly";break;
		default :
			retour = "???";break;
		}
		return (retour);
	}
	//	JLA 
	//	this method will count the flows whose paymentDate > asofDate
	//	and return the position of the first of these flows. 
	//	If there are no flows, then first will be meaningless.
	unsigned long countPaymentFlows(unsigned long& first,const ARM_Date&asof) const ; 

	inline void SetYFStartRiskDates(const ARM_Vector& vector){itsYFStartRiskDates = vector;}
	inline const ARM_Vector& GetYFStartRiskDates(void)	{return itsYFStartRiskDates;}

	inline void SetStartRiskDates(const ARM_Vector& vector){itsStartRiskDates = vector;}
	inline const ARM_Vector& GetStartRiskDates(void)	{return itsStartRiskDates;}

	inline void SetYFEndRiskDates(const ARM_Vector& vector){itsYFEndRiskDates = vector;}
	inline const ARM_Vector& GetYFEndRiskDates(void)	{return itsYFEndRiskDates;}

	inline void SetEndRiskDates(const ARM_Vector& vector){itsEndRiskDates = vector;}
	inline const ARM_Vector& GetEndRiskDates(void)	{return itsEndRiskDates;}

	//	Past Fixings Management
	void setPastFixings(const std::map<double,double>&fixings) ;
	bool getPastFixing(const ARM_Date&resetDate,double&value) const;
	bool getLastFixingDate(const ARM_Date&asofDate,ARM_Date&resetDate) const ;
	void updateFwdSpreadsFromPastFixings(); 

	inline bool IsFrozenMaturity() {return (itsFrozenMaturity==CREDIT_DEFAULT_VALUE ? 0 : 1) ;}
	inline ARM_Date GetFrozenMaturity() { return (ARM_Date)itsFrozenMaturity; }
	inline void SetFrozenMaturity(const ARM_Date& date) {itsFrozenMaturity=date.GetJulian();}
	int PaymentIndex(const ARM_Date& AsofDate) const ;

	// bool	IsYFAlreadyComputed() {return itsFlgYF;};

private:
	// JLA: preparation
	//	on va différencier les payDates des endDates
	//	on ne souhaite pas utiliser ces fonctions qui sont celles de la SwapLeg. 
	inline ARM_Vector* GetPaymentDates(void) const ; 
	inline void SetPaymentDates(ARM_Vector* paymentDates) ; 
	void  GetNumFlows(void) const  ; 
	void GetStartDate(); // !! not the one defined in ARM_Security.
//
//	prevent implicit use
//
	void SetMaturity(const ARM_Date& matDate)  ; 
	ARM_Date GetMaturity(void) const  ;
//
public:
	virtual	void	GetSubSchedulesDataFromLabel(vector<double*>& OutputMatrix, vector<string>& OutputLabels, int& OutputNbRows, int& OutputNbCols) {};


};


//	----------------------------------------------------------------------------------------------------------
//	----------------------------------------------------------------------------------------------------------
#endif 
