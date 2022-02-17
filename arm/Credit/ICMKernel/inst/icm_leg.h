
#if !defined(AFX_ICM_LEG_H_)
#define AFX_ICM_LEG_H_

#include "ARMKernel\inst\swapleg.h"
#include "ICMKernel\inst\icm_security.h"
//#include "ICMKernel\inst\icm_credit_index.h"
// #include "gpinflation\infcurv.h"

// using namespace ARM;

class ICM_Leg;
class ICM_Pricer; 
class ICM_Credit_Index; 

void CptCouponsForFeeLeg(ICM_Leg& leg, ARM_Object* object ,const ARM_Date& asof);

/*********************************************************************************/
/*! \class  ICM_Leg icm_leg.h "icm_leg.h"
 *  \author D. ICM_Leg
 *	\version 1.0
 *	\date   Jully 2005
 *	\brief  Description of a Credit Leg */
/***********************************************************************************/

class ICM_Leg : public ARM_SwapLeg  
{
private:
	
	qPAYMENT_PREMIUM_LEG itsAccruedOnDefault; /*!< Reference sur le calcul du recouvrement du coupon */
	bool				 itsIncludeMaturity; /*!< Include Maturity for last coupon computation */
	ARM_ReferenceValue*	 itsNotExchange;	 /*!< Notional eXchange */


	double				itsRefSpread;		// Reference Spread for CDS Index product
	qCredit_Leg_Style	itsCreditLegStyle;	//Credit Leg Style (Defaultable or Premium)
	qCredit_Leg_Type	itsCreditLegType;	//Credit Leg Type	(Funded,Constant Maturity,Forward ...)

	mutable ICM_Security*	itsCreditInfos;		//Special Infos for credit
	ICM_Credit_Index*	itsCreditIndex;		//Credit index 

	// char*				itsFwdFixedDate;	//Forward Fixed Date : after this date coupons are fixed
	ARM_Date*			itsFwdFixedDate ;
	ARM_ReferenceValue*	itsFwdCalcTypes;	//Forward Calculation Types by period

	ARM_ReferenceValue*	itsRefRiskyType;	//Credit Calculation Types by period
	ARM_ReferenceValue*	itsRefRiskyDate;	//Credit Calculation dates by period

	ARM_ReferenceValue*	itsRefPartRate;		//Participation Rate for Coupon computation

	ICM_Pricer*			itsRatesPricer;		//Pricer used for forward rates computation
	ARM_SwapLeg*		itsSwapLeg;			//Swapleg used for forward rates computation 		

	int					itsCreditLag;		//Credit Lag for Cds - in default case payment begin on j + itsCreditLag			
	double				itsBinary;			//Binary Recovery in digital case
	string				itsSingleName;		//Name of the Single Name

	int					itsIntRule;			// required here since the K_MATUNADJUSTED  is not supported by ARM_SsapLeg
	ARM_Date*			itsFstCpnEffDate;   // Used for Index CDS (full coupon on first period) 

	std::map<double,double>	itsPastFixings;		// the first double is a Julian date corresponding to ResetDates.
public:

	void BitwiseCopy(const ARM_Object* srcleg);

	void Copy(const ARM_Object* srcleg);

	virtual ARM_Object* Clone(void);
	virtual ARM_Object* Clone() const; 


/* ***************************************************************************************************************** */
/*! \fn ICM_Leg(void) 
	\brief Default Leg constructor
/* ***************************************************************************************************************** */
	ICM_Leg() : ARM_SwapLeg() { Init();}		

	ICM_Leg(const ARM_SwapLeg& swapleg,
			const qCredit_Leg_Type& LegType = qRunning_Leg,
			ICM_Credit_Index*	CreditIndex = NULL) : ARM_SwapLeg(swapleg)
	{ 
		Init();
		Set(LegType,CreditIndex);
	}		


	void Set(const qCredit_Leg_Type& LegType = qRunning_Leg,
			 ICM_Credit_Index*	CreditIndex = NULL); 
	/*{
		SetCreditLegType(LegType);
		
		if (CreditIndex)
		{
		if (itsCreditIndex)
			delete itsCreditIndex;
		itsCreditIndex = (ICM_Credit_Index*)CreditIndex->Clone();
		}
		ResetCreditInfos();
	}*/

/* ***************************************************************************************************************** */
/*! \fn ~ICM_Leg(void) 
	\brief Default Leg destructor
/* ***************************************************************************************************************** */
	~ICM_Leg() ;
	/*{
		if (itsCreditInfos)
			delete itsCreditInfos;
		itsCreditInfos = NULL;

		if (itsNotExchange)
			delete itsNotExchange;
		itsNotExchange = NULL;

		if (itsCreditIndex)
			delete itsCreditIndex;
		itsCreditIndex = NULL;

		if (itsFwdFixedDate)
			delete itsFwdFixedDate;
		itsFwdFixedDate = NULL;

		if (itsFwdCalcTypes)
			delete itsFwdCalcTypes;
		itsFwdCalcTypes = NULL;

		if (itsRefRiskyType)
			delete itsRefRiskyType;
		itsRefRiskyType = NULL;

		if (itsRefRiskyDate)
			delete itsRefRiskyDate;
		itsRefRiskyDate = NULL;

		if (itsSwapLeg)
			delete itsSwapLeg;
		itsSwapLeg = NULL;

		if (itsRefPartRate)
			delete itsRefPartRate;
		itsRefPartRate = NULL;

		if (itsFstCpnEffDate) 
			delete itsFstCpnEffDate;
		itsFstCpnEffDate = NULL; 

	}*/

/* ***************************************************************************************************************** */
/*! \fn Init(void) 
	\brief Initializer of ICM_Leg Class
/* ***************************************************************************************************************** */
	void Init(void);


/* ***************************************************************************************************************** */
/*! \fn ICM_Leg(int LesConstructeurs_pour_Jambes_à_Coupon_Fixe) 
	\brief Constructors of Fixed leg for bonds
	\note : Damien will define the best constructors to perform our needs
	\note : At the construction of the leg, the forward rates will be equal to the istlastIndexFixing. 
	\note : The pricing method will reset the forward rates and the fixed rate of the last period. 
	\note : All methods of the <B> ARM_Swapleg </B> should throw an exeption except if there are overloaded or agreed upon for ICM use*/
/* ***************************************************************************************************************** */
	ICM_Leg(const ARM_SwapLeg& SwapLegIn) : ARM_SwapLeg(SwapLegIn)
	{
		Init();
	}

/* ***************************************************************************************************************** */
/*! \fn ICM_Leg(int LesConstructeurs_pour_Jambes_à_Coupon_Fixe) 
	\brief Constructors of Fixed leg for bonds
	\note : Damien will define the best constructors to perform our needs
	\note : At the construction of the leg, the forward rates will be equal to the istlastIndexFixing. 
	\note : The pricing method will reset the forward rates and the fixed rate of the last period. 
	\note : All methods of the <B> ARM_Swapleg </B> should throw an exeption except if there are overloaded or agreed upon for ICM use*/
/* ***************************************************************************************************************** */
	ICM_Leg(const ARM_Date& startDate, 
            const ARM_Date& endDate, 
			const ARM_Date* refDate,
			const ARM_Date* FstCpnEffDate,
			const double& fixedRate /* =0. */ ,
			const qPAYMENT_PREMIUM_LEG& AccruedOnDefault /*= qACCRUED_SETTLED */,
			const int& AccruedDayCount /*= KACTUAL_365*/,
			const double& LastIndexFixing /*=0.*/,
			const int& rcvOrPay	/*= K_RCV*/, 
			const int& freq /*= K_ANNUAL*/, 
            const int& dayCount/*= K30_360*/, 
            const int& decompFreq /*= K_COMP_PROP*/,
            const int& payTiming /*= K_ARREARS*/,
            const int& intRule /*= K_ADJUSTED*/,
            const int& stubRule /*= K_SHORTSTART*/,
            const std::string& discountCcy /* ARM_Currency* discountCcy = ARM_DEFAULT_CURRENCY*/ ,
            const std::string& payCalName /*= NULL*/,
            const int& nxChange /*= K_NX_NONE*/,
			const bool& includematurity /*= EXCLUDE_MATURITY*/,
			const int& adjStartDate /*= K_ADJUSTED*/,
			const qCredit_Leg_Type& LegType /*= qRunning_Leg*/,
			const double& Binary /*= CREDIT_DEFAULT_VALUE*/,
			const string& name /*= ISSUER_UNDEFINE*/) : ARM_SwapLeg()
	{
		Init();

		Set(startDate, endDate, refDate,FstCpnEffDate,
			fixedRate, AccruedOnDefault, 
			AccruedDayCount, LastIndexFixing, rcvOrPay , freq, dayCount, decompFreq, payTiming, intRule,
			stubRule, discountCcy, payCalName, nxChange, // refDate,
			includematurity,adjStartDate,LegType,Binary,name);	
	}


	void Set(const ARM_Date& startDate, 
			 const ARM_Date& endDate, 
			 const ARM_Date* refDate,
   			 const ARM_Date* FstCpnEffDate,
			 const double& fixedRate /* =0.*/,
			 const qPAYMENT_PREMIUM_LEG& /* AccruedOnDefault = qACCRUED_SETTLED*/ ,
			 const int& AccruedDayCount /* = KACTUAL_365*/ ,
			 const double& LastIndexFixing /* =0.*/ ,
			 const int& rcvOrPay /* = K_RCV*/ , 
			 const int& freq /* = K_ANNUAL*/ , 
			 const int& dayCount/* = K30_360*/ , 
			 const int& decompFreq /* = K_COMP_PROP*/ ,
			 const int& payTiming /* = K_ARREARS*/ ,
			 const int& intRule /* = K_ADJUSTED*/ ,
			 const int& stubRule /* = K_SHORTSTART*/ ,
			 const std::string& discountCcy /* = ARM_DEFAULT_CURRENCY*/ ,
			 const std::string& payCalName /* = NULL*/ ,
			 const int& nxChange /* = K_NX_NONE*/ ,
			 const bool& includematurity /* = EXCLUDE_MATURITY*/ ,
			 const int& adjStartDate /* = K_ADJUSTED*/ ,
			 const qCredit_Leg_Type& /* LegType = qRunning_Leg*/ ,
			 const double& Binary /* = CREDIT_DEFAULT_VALUE*/ ,
			 const string& name /* = ISSUER_UNDEFINE*/ );

/* ***************************************************************************************************************** */
/*! \fn ICM_Leg(int LesConstructeurs_pour_Jambes_à_Coupon_Fixe) 
	\brief Constructors of Fixed leg for bonds
	\note : Damien will define the best constructors to perform our needs
	\note : At the construction of the leg, the forward rates will be equal to the istlastIndexFixing. 
	\note : The pricing method will reset the forward rates and the fixed rate of the last period. 
	\note : All methods of the <B> ARM_Swapleg </B> should throw an exeption except if there are overloaded or agreed upon for ICM use*/
/* ***************************************************************************************************************** */

	ICM_Leg(const ARM_Date& startDate, 
            const ARM_Date& endDate, 
			const ARM_Date* refDate,
  			const ARM_Date* FstCpnEffDate,
			const double& fixedRate /* =0.*/ ,
			const int& AccruedDayCount /* = KACTUAL_365*/ ,
			const int& freq /* = K_ANNUAL*/ , 
            const int& dayCount/* = K30_360*/ , 
            const int& payTiming /* = K_ARREARS*/ ,
            const int& intRule /* = K_ADJUSTED*/ ,
            const int& stubRule /* = K_SHORTSTART*/ ,
            const std::string& discountCcy /* = ARM_DEFAULT_CURRENCY*/ ,
            const std::string& payCalName /* = NULL*/ ,
			const bool& includematurity /* = EXCLUDE_MATURITY*/ ,
			const int& adjStartDate /* = K_ADJUSTED*/ ,
			const qCredit_Leg_Type& LegType/* LegType = qRunning_Leg*/ ,
			const double& Binary /* = CREDIT_DEFAULT_VALUE*/ ,
			const string& name /* = ISSUER_UNDEFINE*/ ) : ARM_SwapLeg()
	{
		Init();

		Set(startDate,endDate,refDate,FstCpnEffDate,fixedRate,qACCRUED_SETTLED,AccruedDayCount,0.,	K_RCV,
			freq,dayCount,K_COMP_PROP,payTiming,intRule,stubRule,discountCcy,payCalName,K_NX_NONE,// refDate,
			includematurity,adjStartDate,LegType,Binary,name);	
	}

/* ***************************************************************************************************************** */
/*! \fn ICM_Leg (Constructor for Floating Leg with char* refDate) 
	\brief  Constructors of Floating leg for bonds (Interface)
	\note : Damien will define the best constructors to perform our needs
	\note : At the construction of the leg, the forward rates will be equal to the istlastIndexFixing. 
	\note : The pricing method will reset the forward rates and the fixed rate of the last period. 
	\note : All methods of the <B> ARM_Swapleg </B> should throw an exeption except if there are overloaded or agreed upon for ICM use
	\note : Reference date is given with an ARM_Date */
/* ***************************************************************************************************************** */
	ICM_Leg(const ARM_Date& startDate, 
            const ARM_Date& endDate, 
			const ARM_Date* refDate,
   			const ARM_Date* FstCpnEffDate,
			ARM_IRIndex* irIndex,
			const double& spread /* = 0.0*/ ,
			const qPAYMENT_PREMIUM_LEG& AccruedOnDefault /*= qACCRUED_SETTLED*/,
			const int&				    AccruedDayCount /*= KACTUAL_365*/,
			const double&				InitialRate /*=0.*/,
			const double&				LastIndexFixing /*=0.*/,
			const int& rcvOrPay /*= K_PAY*/, 
			const int& dayCount/*= K30_360*/, 
            const int& decompFreq /*= K_COMP_PROP*/,
            const int& stubRule /*= K_SHORTSTART*/,
			const int& resetgap /*= 10000*/,
			const std::string& resetCalName /*= NULL*/, 
            const std::string& discountCcy /*= ARM_DEFAULT_CURRENCY*/,
            const std::string& payCalName /*= NULL*/,
            const int& nxChange /*= K_NX_NONE*/,
			const bool& includematurity /*= EXCLUDE_MATURITY*/,
			const int& adjStartDate /*= K_ADJUSTED*/,
			const qCredit_Leg_Type& LegType /*= qRunning_Leg*/,
			const double& Binary /*= CREDIT_DEFAULT_VALUE*/,
			const string& name /*= ISSUER_UNDEFINE*/) : ARM_SwapLeg() 
	{

		Init();

		Set( startDate, endDate, refDate,FstCpnEffDate,irIndex, // refDate, 
			spread, AccruedOnDefault, AccruedDayCount, InitialRate, LastIndexFixing, rcvOrPay, 
			 dayCount, decompFreq, stubRule, resetgap, resetCalName, discountCcy, payCalName, nxChange,includematurity,adjStartDate,LegType,Binary,name);
		
	}

	void Set( const ARM_Date& startDate, 
			 const ARM_Date& endDate, 
			 const ARM_Date* refDate,
   			 const ARM_Date* FstCpnEffDate,
			 ARM_IRIndex* irIndex,
			 const double& spread /* = 0.0*/,
			 const qPAYMENT_PREMIUM_LEG& AccruedOnDefault /* = qACCRUED_SETTLED*/,
			 const int&				    AccruedDayCount /* = KACTUAL_365*/,
			 const double&				InitialRate /* = 0.*/,
			 const double&				LastIndexFixing /* =0.*/,
			 const int& rcvOrPay /* = K_PAY*/, 
			 const int& dayCount/* = K30_360*/, 
			 const int& decompFreq /* = K_COMP_PROP*/,
			 const int& stubRule /* = K_SHORTSTART*/,
			 const int& resetgap /* = 10000*/,
			 const std::string& resetCalName /* = NULL*/, 
			 const std::string& discountCcy /* = ARM_DEFAULT_CURRENCY*/,
			 const std::string& payCalName /* = NULL*/,
			 const int& nxChange /* = K_NX_NONE*/,
			 const bool& includematurity /* = EXCLUDE_MATURITY*/,
			 const int& adjStartDate /* = K_ADJUSTED*/,
			 const qCredit_Leg_Type& LegType /* = qRunning_Leg*/,
			 const double& Binary /* = CREDIT_DEFAULT_VALUE*/,
			 const string& name /* = ISSUER_UNDEFINE*/ );

	/* ***************************************************************************************************************** */
	// Constructor general
	/* ***************************************************************************************************************** */
	void Set(const ARM_Date& StartDate, 
            const ARM_Date& EndDate,
			const ARM_Date* refDate,
   			const ARM_Date* FstCpnEffDate,
			const double& FixedRate,
			const double& FixedNotional,	
		ARM_ReferenceValue* Notional /*= NULL*/,
		ARM_ReferenceValue* Rates /*= NULL*/,
		ARM_ReferenceValue* Exchange /*= NULL*/,
		const int& frequency /* = K_ANNUAL*/, 
            const int& dayCount /* = K30_360*/, 
            const int& payTiming /* = K_ARREARS*/,
            const int& intRule /* = K_ADJUSTED*/,
            const int& stubRule /* = K_SHORTSTART*/,
            const std::string& discountCcy /* = ARM_DEFAULT_CURRENCY*/,
            const std::string& payCalName /* = NULL*/,
		const qCredit_Leg_Type& LegType /* = qRunning_Leg*/,
		const bool& includematurity /* = EXCLUDE_MATURITY*/,
		const int& adjStartDate /* = K_ADJUSTED*/,
		ARM_IRIndex* irIndex /* = NULL*/,
		const double& Binary /* = CREDIT_DEFAULT_VALUE*/,
		const string& name /* = ISSUER_UNDEFINE*/,
		const int& NXchange /* = K_NX_NONE*/,
		const qPAYMENT_PREMIUM_LEG& AccruedOnDefault /* = qACCRUED_SETTLED */ );


	ICM_Leg(const ARM_Date& StartDate, 
            const ARM_Date& EndDate, 
			const ARM_Date* refDate,
   			const ARM_Date* FstCpnEffDate,
			const double& FixedRate,
			const double& FixedNotional,
			ARM_ReferenceValue* Notional /*= NULL*/,
			ARM_ReferenceValue* Rates /*= NULL*/,
			ARM_ReferenceValue* Exchange /*= NULL*/,
			const int& frequency /*= K_ANNUAL*/, 
            const int& dayCount/*= K30_360*/, 
            const int& payTiming /*= K_ARREARS*/,
            const int& intRule /*= K_ADJUSTED*/,
            const int& stubRule /*= K_SHORTSTART*/,
            const std::string& discountCcy /*= ARM_DEFAULT_CURRENCY*/,
            const std::string& payCalName /*= NULL*/,
			const qCredit_Leg_Type& LegType /*= qRunning_Leg*/,
			const bool& includematurity /*= EXCLUDE_MATURITY*/,
			const int& adjStartDate /*= K_ADJUSTED*/,
			ARM_IRIndex* irIndex /*= NULL*/,
			const double& Binary /*= CREDIT_DEFAULT_VALUE*/,
			const string& name /*= ISSUER_UNDEFINE*/,
			const int& NXchange /*= K_NX_NONE*/,
			const qPAYMENT_PREMIUM_LEG& AccruedOnDefault /* =qACCRUED_SETTLED */ ) : ARM_SwapLeg()
	{
		Init();

		Set(StartDate,EndDate,refDate,FstCpnEffDate,FixedRate,FixedNotional,
				Notional,Rates,Exchange,frequency,dayCount
				,payTiming,intRule,stubRule,discountCcy,payCalName,
				LegType,includematurity,adjStartDate, irIndex,Binary,name,NXchange,AccruedOnDefault);	
	}


	/* ***************************************************************************************************************** */
	// Libor & Euribor Risky Leg
	/* ***************************************************************************************************************** */
	void Set(const ARM_Date& StartDate, 
            const ARM_Date& EndDate, 
			const ARM_Date* refDate,
   			const ARM_Date* FstCpnEffDate,
			const ARM_INDEX_TYPE& liborType /*= EURIBOR6M*/,
			const double& spread /*= 0.0*/, 
			const int& resetFreq /*= K_DEF_FREQ*/, 
			const int& payFreq /*= K_DEF_FREQ*/, 
            const int& resetTiming /*= K_ADVANCE*/, 
			const int& payTiming /*= K_ARREARS*/,
            const std::string& discountCcy /*= ARM_DEFAULT_CURRENCY*/,
            const int& intRule /*= K_ADJUSTED*/,
            const int& resetgap /*= 10000*/,
            const std::string& resetCalName /*= NULL*/,
            const std::string& payCalName /*= NULL*/,
            const int& stubRule /*= K_SHORTSTART*/,
            const int& adjStartDate /*= K_ADJUSTED*/,
			const int& dayCount /*= -1*/,
			const bool& includematurity /*= EXCLUDE_MATURITY*/,
			const double& Binary /*= CREDIT_DEFAULT_VALUE*/,
			const string& name /*= ISSUER_UNDEFINE*/);


	ICM_Leg(const ARM_Date& StartDate, 
            const ARM_Date& EndDate, 
			const ARM_Date* refDate,
   			const ARM_Date* FstCpnEffDate,
			const ARM_INDEX_TYPE& liborType /*= EURIBOR6M*/,
			const double& spread /*= 0.0*/, 
			const int& resetFreq /*= K_DEF_FREQ*/, 
			const int& payFreq /*= K_DEF_FREQ*/, 
            const int& resetTiming /*= K_ADVANCE*/, 
			const int& payTiming /*= K_ARREARS*/,
            const std::string&  discountCcy /*= ARM_DEFAULT_CURRENCY*/,
            const int& intRule /*= K_ADJUSTED*/,
            const int& resetgap /*= 10000*/,
            const std::string&  resetCalName /*= NULL*/,
            const std::string&  payCalName /*= NULL*/,
            const int& stubRule /*= K_SHORTSTART*/,
            const int& adjStartDate /*= K_ADJUSTED*/,
			const int& dayCount /*= -1*/,
			const bool& includematurity /*= EXCLUDE_MATURITY*/,
			const double& Binary /*= CREDIT_DEFAULT_VALUE*/,
			const string& name /*= ISSUER_UNDEFINE*/) : ARM_SwapLeg()
	{
		Init();

		Set(StartDate,EndDate,refDate,FstCpnEffDate,liborType,spread,
			resetFreq,payFreq,resetTiming,payTiming,discountCcy,
            intRule,resetgap,resetCalName,payCalName,stubRule,
            // refDate,
			adjStartDate,dayCount,
			includematurity,Binary,name);	
	}


	virtual void Rebuild(void);
	
	// ----------------------------------------------------------------------
	//Summit Leg
	// ----------------------------------------------------------------------

	ICM_Leg(ARM_Vector* flowStartDates, 
		    ARM_Vector* flowEndDates,
            ARM_Vector* paymentDates, 
			ARM_Vector* resetDates, 
            ARM_Vector* intDays, 
			ARM_Vector* fwdRates,
            ARM_ReferenceValue* Notional, 
			ARM_IRIndex* irIndex,
            const int& rcvOrPay /*= K_PAY*/ , 
			const double& spread /*= 0.0*/ ,
            const double& fixRate /*= K_FLOAT_RATE*/ , 
            const std::string& discountCcy /*= ARM_DEFAULT_CURRENCY*/ ,
            const int& NxId /*= K_NX_NONE*/ ,
			const int& decompPricingFlag /*= 1*/ ,
			const qCredit_Leg_Type& LegType /*= qRunning_Leg*/ ,
			const double& Binary /*= CREDIT_DEFAULT_VALUE*/ ,
			const string& name /*= ISSUER_UNDEFINE */ )
	{
		Init();

		Set(flowStartDates,flowEndDates,paymentDates,resetDates,intDays,fwdRates,Notional, 
			irIndex,rcvOrPay,spread,fixRate,discountCcy,NxId,decompPricingFlag,LegType,Binary,name);
	}

	void Set(ARM_Vector* flowStartDates, 
		    ARM_Vector* flowEndDates,
            ARM_Vector* paymentDates, 
			ARM_Vector* resetDates, 
            ARM_Vector* intDays, 
			ARM_Vector* fwdRates,
            ARM_ReferenceValue* Notional, 
			ARM_IRIndex* irIndex,
            const int& rcvOrPay /* = K_PAY*/, 
			const double& spread /* = 0.0*/,
            const double& fixRate /* = K_FLOAT_RATE*/, 
            const std::string& discountCcy /* = ARM_DEFAULT_CURRENCY*/,
            const int& NxId /* = K_NX_NONE*/,
			const int& decompPricingFlag /* = 1*/,
			const qCredit_Leg_Type& LegType /* = qRunning_Leg*/,
			const double& Binary /* = CREDIT_DEFAULT_VALUE*/,
			const string& name /* = ISSUER_UNDEFINE*/);


/* ****************************************************************************************************************** */
/*!	\fn    PeriodIndex(ARM_Date& AsofDate); 
    \brief Retrieves the current period index such that \f[ max \{ StartDate(PeriodIndex) \le t \} \f] 
	\param AsofDate
	\return integer: The period index*/
/* ****************************************************************************************************************** */
	int PeriodIndex(const ARM_Date& AsofDate);

	//Get Methods *********************************************
	inline ARM_Vector* GetCashFlowValues(void) const
	{
		
    if (ARM_SwapLeg::GetCashFlowValues())
		return(ARM_SwapLeg::GetCashFlowValues()); 
	else
       throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "ERROR : vector itsFlowValues unitialized ");
	}


	inline qPAYMENT_PREMIUM_LEG GetAccruedOnDefault(void) { return itsAccruedOnDefault;}	
 	inline void SetAccruedOnDefault(const qPAYMENT_PREMIUM_LEG& value) {itsAccruedOnDefault = value;}

	void View(char* id = NULL, FILE* fOut = NULL);

/* ****************************************************************************************************************** */
/*! \fn GetCreditSpread()
	\brief returns fixed coupon */
/* ****************************************************************************************************************** */
	virtual double GetCreditSpread()
	{
		GetCreditInfos(); // make sure ARM_Swapleg::CptDates is called.... 

		double InitialSpread = 0.;

		if (IsFixedLeg())
			InitialSpread = GetFixedRate();
		else
			InitialSpread = GetSpread();

		return (InitialSpread);
	}

/* ****************************************************************************************************************** */
/*! \fn SetCreditSpread()
	\brief set fixed coupon to fee leg*/
/* ****************************************************************************************************************** */
	virtual void SetCreditSpread(const double& spread)
	{
		double InitialSpread = 0.;

		SetLegType(K_FIXED_LEG);
		GetCreditInfos()->SetCreditSpread(spread);

		if (IsFixedLeg())
			SetFixedRate(spread);
		else
			SetSpread(spread);
	}

	double GetNotionalAmount() const 
	{
		//GetCreditInfos(); // make sure ARM_Swapleg::CptDates is called.... 
		return GetAmount()->GetDiscreteValues()->Elt(0);
	}

    void GenScheduleDates(int datesType, int viewInitExch, int modfoll, int creditgap, char* id = NULL, FILE* fOut = NULL);

	inline void ResetCreditInfos(void) 
	{	if (itsCreditInfos)	delete itsCreditInfos;
		itsCreditInfos = NULL;}

	inline ICM_Security* GetCreditInfos(void)  const 
	{if (itsCreditInfos == NULL)	
			CptCreditInfos();
	 return (itsCreditInfos);}

	const ICM_Security& GetCreditInfosRef() const 
	{if (!GetCreditInfos()) 
		ICMTHROW(ERR_INVALID_ARGUMENT,"GetCreditInfos: no value"); 
	 return *GetCreditInfos();}

	// inline ICM_Security* GetCreditInfosPointer(void) {	return (itsCreditInfos);}
	// inline void SetCreditInfosPointer(ICM_Security* p) {itsCreditInfos = p;}

	inline void ComputeYF(const ARM_Date& AsOf) 
	{
		GetCreditInfos()->ComputeYF(AsOf);
	}

	virtual void SetVariableSpread(ARM_ReferenceValue* spreads) ;
	// JLA moved to .cpp {
	// JLA moved to .cpp 	SetLegType(K_FLOATING_LEG);
	// JLA moved to .cpp	ARM_SwapLeg::SetVariableSpread(spreads);
	// JLA moved to .cpp
	// JLA moved to .cpp	
	// JLA moved to .cpp	/** JLA this logic is done in GetCreditInfos()
	// JLA moved to .cpp	ICM_Security* security = GetCreditInfosPointer();
	// JLA moved to .cpp
	// JLA moved to .cpp	if (security==NULL) CptCreditInfos();
	// JLA moved to .cpp	else
	// JLA moved to .cpp	{
	// JLA moved to .cpp	**/ 
	// JLA moved to .cppICM_Security* security = GetCreditInfos(); 
	// JLA moved to .cppint NbFlows = GetNumFlows();
	// JLA moved to .cppARM_Vector CouponRates(NbFlows,0.);
	// JLA moved to .cppfor (int i=0;i<NbFlows;i++)
	// JLA moved to .cpp{CouponRates.Elt(i) = spreads->CptReferenceValue((ARM_Date)GetFlowStartDates()->Elt(i));}
	// JLA moved to .cppsecurity->SetCouponRates(&CouponRates);
	// JLA moved to .cpp/**}**/ 
	// JLA moved to .cpp }	

	virtual void SetNotExchange(ARM_ReferenceValue* NotExchange)
	{
		ResetCreditInfos(); 
		// if (itsCreditInfos)
		// 	delete itsCreditInfos;
		// itsCreditInfos = NULL;

		if (itsNotExchange)
			delete itsNotExchange;
		itsNotExchange = NULL;

		itsNotExchange = NotExchange;
		
		// JLA. CptCreditInfos();
	}	

	inline ARM_ReferenceValue* GetNotExchange(void) { return itsNotExchange;};
	
	void CptCashFlowDatesCredit(void) ;

	inline qCredit_Leg_Style	GetCreditLegStyle() {return itsCreditLegStyle;}
	inline qCredit_Leg_Type GetCreditLegType() {return itsCreditLegType;}

	inline void SetCreditLegStyle(const qCredit_Leg_Style& style) { itsCreditLegStyle=style;}
	void SetCreditLegType(const qCredit_Leg_Type &type);

	inline double GetRefSpread() const { return itsRefSpread ; }
	inline void SetRefSpread(const double& Spread) { itsRefSpread = Spread ;}

	inline ICM_Credit_Index* GetCreditIndex() const {return itsCreditIndex;}
	void SetCreditIndex(ICM_Credit_Index* index); 

	virtual void CptCoupons(ARM_Object* object,const ARM_Date& asof) {CptCouponsForFeeLeg(*this,object,asof);}

	void ResetCoupons(ARM_Model* model);

    inline void SetPaymentFreq(const int& payFreq)
    { 
	   ARM_SwapLeg::SetPaymentFreq(payFreq); 
	   GetCreditInfos()->SetPaymentFreq(payFreq);
    }

	inline const ARM_Date* GetFwdFixedDate() const { return itsFwdFixedDate; }
	inline void SetFwdFixedDate(const ARM_Date*date) 
	{
		if (itsFwdFixedDate) delete itsFwdFixedDate; itsFwdFixedDate=0; 
		if (date) itsFwdFixedDate=dynamic_cast<ARM_Date*>(new ARM_Date(*date)); 
		ResetCreditInfos(); 
	}

	const ARM_ReferenceValue* GetRefRiskyDate() const {return itsRefRiskyDate;}
	const ARM_ReferenceValue* GetRefRiskyType() const {return itsRefRiskyType;}
	const ARM_ReferenceValue* GetRefPartRate() const  {return itsRefPartRate;}

	inline void SetRefRiskyDate(ARM_ReferenceValue* ref1,
								ARM_ReferenceValue* ref2) 
	{
		if (itsRefRiskyDate)
			delete itsRefRiskyDate;
		itsRefRiskyDate = (ARM_ReferenceValue*) ref1->Clone();

		if (itsRefRiskyType)
			delete itsRefRiskyType;
		itsRefRiskyType = (ARM_ReferenceValue*) ref2->Clone();

		// CptCreditInfos();
		ResetCreditInfos(); 
	}

	const ARM_ReferenceValue* GetFwdCalcTypes() const {return itsFwdCalcTypes;}
	inline void SetFwdCalcTypes(ARM_ReferenceValue* ref) 
	{
		if (itsFwdCalcTypes)
			delete itsFwdCalcTypes;
		itsFwdCalcTypes = (ARM_ReferenceValue*) ref->Clone();

		// CptCreditInfos();
		ResetCreditInfos(); 
	}

	inline void SetRefPartRate(ARM_ReferenceValue* ref) 
	{
		if (itsRefPartRate)
			delete itsRefPartRate;
		itsRefPartRate = (ARM_ReferenceValue*) ref->Clone();

		// CptCreditInfos();
		ResetCreditInfos(); 
	}

	void SetRatesPricer(ICM_Pricer* pricer) {itsRatesPricer=pricer;}
	ICM_Pricer* GetRatesPricer(void) {return itsRatesPricer;}

private:
	void CptCreditInfos(void) const ;	
public:

	inline int GetCreditLag() const { return itsCreditLag;}
	inline void SetCreditLag(const int& Lag) {itsCreditLag = Lag ;}

	bool GetBinaryFlg(void) 
	{ if (itsBinary == CREDIT_DEFAULT_VALUE) return false;
		else
		return true;
	}

	inline double GetBinary(void) { return itsBinary; }
	inline void SetBinary(const double& binary) { itsBinary=binary; }

	const std::string& GetSingleName() const { return itsSingleName; }
	inline void SetSingleName(const string& name){itsSingleName=name;}

	inline void CptCashFlowDatesCredit(ARM_Date start,ARM_Date end)
	{
		ARM_SwapLeg::SetStartDate(start);	
		ARM_SwapLeg::SetStartDateNA(start);
		ARM_SwapLeg::SetEndDate(end);
		ARM_SwapLeg::SetEndDateNA(end);
		CptCashFlowDatesCredit();
	}

	inline void CptCashFlowDatesCredit(ARM_Date end)
	{
		ARM_SwapLeg::SetEndDate(end);
		ARM_SwapLeg::SetEndDateNA(end);
		CptCashFlowDatesCredit();
	}

	inline void CptCashFlowDatesCredit(const int& Frequency)
	{
		SetPaymentFreq(Frequency);
		GetIRIndex()->SetPayFrequency(Frequency);
		GetIRIndex()->SetResetFrequency(Frequency);

		CptCashFlowDatesCredit();
	}

	inline void SetSwapLegForRatesComputation(ARM_SwapLeg* leg)
	{
		if (itsSwapLeg)
			delete itsSwapLeg;
		itsSwapLeg = (ARM_SwapLeg*) leg->Clone();
	}

	inline ARM_SwapLeg* GetSwapLegForRatesComputation() {return itsSwapLeg;}


	// inline const ARM_Date* GetFstCpnEffDate() const {return itsFstCpnEffDate;}

private:
	void SetFstCpnEffDate(const ARM_Date* FstCpnEffDate) 
	{
		if (itsFstCpnEffDate)
		{delete itsFstCpnEffDate;
		itsFstCpnEffDate = NULL;
		}
		if (FstCpnEffDate) 
			itsFstCpnEffDate = new ARM_Date (*FstCpnEffDate);
	}
public:
	void setPastFixing(const ARM_Date&resetDate,double value) ;
	inline bool	 GetIncludeMaturity() {return itsIncludeMaturity;}

}; // End of ICM_Leg

#endif // !defined(AFX_ICM_LEG_H_)
