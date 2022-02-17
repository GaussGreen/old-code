
#pragma warning (disable : 4786 )

#if !defined(_ICM_CDS_H_)
#define _ICM_CDS_H_

#include "ICMKernel\inst\icm_leg.h"

/*********************************************************************************/
/*! \class  ICM_Cds icm_cds.h "icm_cds.h"
 *  \author 
 *	\version 1.0
 *	\date   June 2003
 *	\file   icm_cds.h
 *	\brief  Pricing <B> Credit Default Swap </B>
	The definition of the members is closed to the ISDA defintion. 
	 <ul> Which can be: 
			<li> Single Name </li>
			<li >First to Default </li>
			<li> Nth to Default </li>
	</ul>
			depending of the model used to price it */
/***********************************************************************************/

class ICM_Credit_Index;

class ICM_Cds : public ICM_Security  
{
private: 

	ICM_Leg*			itsFeeLeg;			//Premium Leg		
	ICM_Leg*			itsDefLeg;			//Defaultable Leg	
	double				itsTradedCoef; 
protected:

	ARM_Vector*		itsSchedule;			//Merge Schedule for Cds
	ARM_Vector*		itsRiskSchedule;		//Merge Risk Schedule for Cds
	// ARM_ReferenceValue*	itsYearTermSearch;	//used for seraches into cds schedule

private: 
	
	void Init();

public:


	ICM_Cds() {	Init();}	
	ICM_Cds(const ICM_Cds& cds) ;
	inline void SetTradedCoef(const double& value) { itsTradedCoef=value; }
	inline double GetTradedCoef(void) const { return itsTradedCoef; }

	// --------------------------------------------------------------------------
	// Constructeur à partir de 2 legs
	// --------------------------------------------------------------------------

	ICM_Cds(ICM_Leg* Feeleg, 
			ICM_Leg* Defleg,
			const int& creditlag = DEFAULT_CREDIT_LAG,
			const string& name = "default");


	void Set(ICM_Leg* Feeleg,
			 ICM_Leg* Defleg,
			const int& creditlag = DEFAULT_CREDIT_LAG,
			const string& name = "default");

	// --------------------------------------------------------------------------
	// Constructeur avec char* FirstPeriodReferenceDate = "NULL"
	// --------------------------------------------------------------------------

	ICM_Cds(const ARM_Date& EffectiveDate,
				const ARM_Date& ScheduleTerminateDate,
				const ARM_Date* FirstPeriodReferenceDate,
				const ARM_Date* FstCpnEffDate,
				const ARM_Date& ProtectionStartDate,
				const ARM_Date& ProtectionEndDate,
				const double& FixedRate,
				const ARM_ReferenceValue& premiumNotionals,
				const ARM_ReferenceValue& defaultNotionals,
				const int& FrequencyFeeLeg /* = K_QUARTERLY */ ,
				const int& DayCountBasis /* = KACTUAL_360 */ ,
				const qPAYMENT_PREMIUM_LEG& AccruedOnDefault /* = qACCRUED_SETTLED */ ,
				const std::string& ccy/* ARM_Currency *Ccy = ARM_DEFAULT_CURRENCY */ , 
				const int& stubrule /* = K_SHORTSTART*/ , /* 1 <-> SHORTSTART */
				const double& CreditLag /*=  DEFAULT_CREDIT_LAG */,
				const int& FrequencyDefLeg /* DEFAULT_FRQ_DEFLEG */,
				const int& intRule /* K_ADJUSTED */,
				const bool& includematurity /* INCLUDE_MATURITY*/ ,
				const int& adjStartDate /* K_ADJUSTED*/ ,
				const std::string& payCalName/* char* payCalName = NULL */ ,
				const qCredit_Leg_Type& TypeFeeLeg/* =qRunning_Leg*/ ,
				const qCredit_Leg_Type& TypeDefLeg/* =qStandart_Recovery_Leg*/ ,
				const string& name /* = ISSUER_UNDEFINE */ ,
				const double& Binary /* = CREDIT_DEFAULT_VALUE*/ ,
				const long l_NotionalEch_Type =  K_NX_NONE,
				const ARM_ReferenceValue* l_NotionalEchange = NULL
				);

	void Set(const ARM_Date& EffectiveDate,
				const ARM_Date& ScheduleTerminateDate,
				const ARM_Date* FirstPeriodReferenceDate,
				const ARM_Date* FstCpnEffDate,
				const ARM_Date& ProtectionStartDate,
				const ARM_Date& ProtectionEndDate,
				const double& FixedRate,
				const ARM_ReferenceValue& premiumNotionals,
				const ARM_ReferenceValue& defaultNotionals,
				const int& FrequencyFeeLeg/*  = K_QUARTERLY */ ,
				const int& DayCountBasis /* = KACTUAL_360 */ ,
				const qPAYMENT_PREMIUM_LEG& /* AccruedOnDefault = qACCRUED_SETTLED*/ ,
				const std::string& Ccy /* = ARM_DEFAULT_CURRENCY */ , 
				const int& stubrule /* = K_SHORTSTART */ ,
				const double& CreditLag /* = DEFAULT_CREDIT_LAG */ ,
				const int& FrequencyDefLeg /* = DEFAULT_FRQ_DEFLEG */ ,
				const int& intRule /* = K_ADJUSTED */ ,
				const bool& includematurity /* = INCLUDE_MATURITY */ ,
				const int& adjStartDate /* = K_ADJUSTED*/ ,
				const std::string& payCalName/* char* payCalName = NULL*/ ,
				const qCredit_Leg_Type& TypeFeeLeg/* =qRunning_Leg */ ,
				const qCredit_Leg_Type& TypeDefLeg/* =qStandart_Recovery_Leg*/ ,
				const string& name /* = ISSUER_UNDEFINE*/ ,
				const double& Binary /* = CREDIT_DEFAULT_VALUE*/,
				const long l_NotionalEch_Type =  K_NX_NONE,
				const ARM_ReferenceValue* l_NotionalEchange = NULL);

// -------------------------------------------------------------------------
// Constructor For CDS Index Product
// -------------------------------------------------------------------------

	
	ICM_Cds(const ARM_Date& EffectiveDate,
			const ARM_Date& MaturityDate,
			const ARM_Date* FirstPeriodReferenceDate ,
			const ARM_Date* FstCpnEffDate,
			const ARM_Date& ProtectionStartDate,
			const ARM_Date& ProtectionEndDate,
			const double& FixedRate,
				double premiumNotional,
				double defaultNotional,
				const ARM_ReferenceValue* premiumNotionals,
				const ARM_ReferenceValue* defaultNotionals,
			ICM_Credit_Index* Index,
			const int& FrequencyFeeLeg /* = K_QUARTERLY*/,
			const int& DayCountBasis /* = KACTUAL_360*/,
			const qPAYMENT_PREMIUM_LEG& AccruedOnDefault /* = qACCRUED_SETTLED*/,
			const std::string& ccy/* = ARM_DEFAULT_CURRENCY*/, 
			const int& stubrule /* = K_SHORTSTART*/,
			const double& CreditLag/*  = DEFAULT_CREDIT_LAG*/,
			const int& FrequencyDefLeg/*  = DEFAULT_FRQ_DEFLEG*/,
			const int& intRule /* = K_ADJUSTED*/,
			const bool& includematurity /* = INCLUDE_MATURITY*/,
			const int& adjStartDate /* = K_ADJUSTED*/,
			const std::string& payCalName /* = NULL*/,
			const qSecurity_TYPE& cdstype /* = qCDS_INDEX*/,
			const string& name /* = ISSUER_UNDEFINE*/,
			const double& Binary /* = CREDIT_DEFAULT_VALUE*/);

	void Set(const ARM_Date& EffectiveDate,
			 const ARM_Date& MaturityDate,
			 const ARM_Date* FirstPeriodReferenceDate ,
			 const ARM_Date* FstCpnEffDate,
			 const ARM_Date& ProtectionStartDate,
			 const ARM_Date& ProtectionEndDate,
			 const double& FixedRate,
				double premiumNotional,
				double defaultNotional,
				const ARM_ReferenceValue* premiumNotionals,
				const ARM_ReferenceValue* defaultNotionals,
			 ICM_Credit_Index* Index,
			 const int& FrequencyFeeLeg /* = K_QUARTERLY*/ ,
			 const int& DayCountBasis /* = KACTUAL_360*/ ,
			 const qPAYMENT_PREMIUM_LEG& AccruedOnDefault /* = qACCRUED_SETTLED*/ ,
			 const std::string& ccy /* = ARM_DEFAULT_CURRENCY */ , 
			 const int& stubrule /* = K_SHORTSTART*/ ,
			 const double& CreditLag /* = DEFAULT_CREDIT_LAG*/ ,
			 const int& FrequencyDefLeg /* = DEFAULT_FRQ_DEFLEG*/ ,
			 const int& intRule /* = K_ADJUSTED*/ ,
			 const bool& includematurity /* = INCLUDE_MATURITY*/ ,
			 const int& adjStartDate /* = K_ADJUSTED*/ ,
			 const std::string& payCalName /* = NULL*/ ,
			 const qSecurity_TYPE& cdstype /* = qCDS_INDEX*/ ,
			 const string& name /* = ISSUER_UNDEFINE*/ ,
			 const double& Binary /* = CREDIT_DEFAULT_VALUE*/ );


	virtual ~ICM_Cds() ;


	 void BitwiseCopy(const ARM_Object* srcleg);
	 void Copy(const ARM_Object* srcleg);

	ARM_Object* Clone(void);
	virtual ARM_Object* Clone() const ;

	void View(char* id, FILE* ficOut);

	inline void SetFeeLeg(ICM_Leg* FeeLeg)
	{
		if (itsFeeLeg)
			delete itsFeeLeg;
		itsFeeLeg = (ICM_Leg*) FeeLeg->Clone();
	}

	inline ICM_Leg* GetFeeLeg(void)	{return itsFeeLeg;}
	ICM_Leg& GetFeeLeg() const 
	{
		if (!itsFeeLeg) 
			ICMTHROW(ERR_INVALID_ARGUMENT,"GetFeeLeg: leg does not exists"); 
		return *itsFeeLeg; 
	}

	inline void SetDefLeg(ICM_Leg* DefLeg)
	{
		if (itsDefLeg)
			delete itsDefLeg;
		itsDefLeg = (ICM_Leg*) DefLeg->Clone();
	}

	inline ICM_Leg* GetDefLeg(void) {return itsDefLeg;}
	ICM_Leg& GetDefLeg() const
	{
		if (!itsDefLeg) 
			ICMTHROW(ERR_INVALID_ARGUMENT,"GetDefLeg: leg does not exists"); 
		return *itsDefLeg; 
	}

	/* ****************************************************************************************************************** */
	/*! \fn GetCreditSpread()
		\brief returns the rate of the premium leg 
	/* ****************************************************************************************************************** */
	virtual double GetCreditSpread()
	{
		double spread = itsFeeLeg->GetCreditSpread();
		return spread;
	}

	virtual void SetCreditSpread(const double& spread)
	{
		ICM_Security::SetCreditSpread(spread);
		itsFeeLeg->SetCreditSpread(spread);
	}

    void DisplayScheduleDates(int datesType,
                              char* id = NULL, FILE* fOut = NULL);

    void DisplayScheduleValues(int valuesType,
                              char* id = NULL, FILE* fOut = NULL);

	inline int GetCreditLag() const { return itsDefLeg->GetCreditLag();}
	inline void SetCreditLag(const int& Lag) {itsDefLeg->SetCreditLag(Lag);}

	inline void CptCashFlowDatesCredit(const ARM_Date &start,const ARM_Date &end)
	{
		SetStartDateNA(start);
		SetEndDateNA(end); 
		itsFeeLeg->CptCashFlowDatesCredit(start,end);
		itsDefLeg->CptCashFlowDatesCredit(start,end);
	}

	inline void CptCashFlowDatesCredit(const ARM_Date &end)
	{
		SetEndDateNA(end);
		itsFeeLeg->CptCashFlowDatesCredit(end);
		itsDefLeg->CptCashFlowDatesCredit(end);
	}

	inline void CptCashFlowDatesCredit(const int& Frequency)
	{
		// SetFrequency(Frequency);
		itsFeeLeg->CptCashFlowDatesCredit(Frequency);
		itsDefLeg->CptCashFlowDatesCredit(Frequency);
	}

	// -------------------------------------------------------------------------------
	//Compute the schedule following start&end periods in year fractions starting at 0
	// -------------------------------------------------------------------------------
	void GenerateSchedule(const ARM_Date& TheDate,const int& nbdays=-1);
	inline ARM_Vector*	GetSchedule() {return itsSchedule;}

	// -------------------------------------------------------------------------------------
	//Compute the schedule following Risky start&end periods in year fractions starting at 0
	// -------------------------------------------------------------------------------------
	void GenerateRiskSchedule(const ARM_Date& TheDate,const int& nbdays=-1);
	inline ARM_Vector*	GetRiskSchedule() {return itsRiskSchedule;}

	int	GetYearTermSearchInSchedule(const double& yearterm) const ;

	qCredit_Leg_Type DeduceFeeLegType(const qSecurity_TYPE& sectype);
	qCredit_Leg_Type DeduceDefLegType(const qSecurity_TYPE& sectype);

	inline void SetSingleName(const string& name) {itsDefLeg->SetSingleName(name);}

	const std::string& GetSingleName() const { return itsDefLeg->GetSingleName();}// return itsCdsSingleName; }

	inline void SetCoupons(ARM_ReferenceValue* cpnrate = NULL,
						   ARM_ReferenceValue* cpnmethod = NULL,
						   ARM_ReferenceValue* PartRates = NULL)	
	{
		if (cpnrate) GetFeeLeg()->SetVariableSpread(cpnrate);
		if (cpnmethod) GetFeeLeg()->SetFwdCalcTypes(cpnmethod);
		if (PartRates) GetFeeLeg()->SetRefPartRate(PartRates);
	}

	bool GetBinaryFlg(void) { return itsDefLeg->GetBinaryFlg();}
	inline double GetBinary(void) { return itsDefLeg->GetBinary(); }
	inline void SetBinary(const double& binary) {itsDefLeg->SetBinary(binary); }

	// inline void SetIncludeMaturity(const bool& value) {GetFeeLeg()->SetIncludeMaturity(value);} 

	void ResetSchedules() ;

};

#endif // !defined(AFX_ICM_CDS_H__0FAFBE8F_B970_476B_AE0F_A498B958D585__INCLUDED_)
