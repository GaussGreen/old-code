
#ifndef _ICM_DEFAULT_CURVE_H
#define _ICM_DEFAULT_CURVE_H

#include "ARMKernel\util\refvalue.h"
#include "ICMKernel\glob\icm_enums.h"
#include "ICMKernel\util\icm_utils.h"
#include "ICMKernel\util\icm_matrix.h"
#include "interpol.h"


// #define _credit_tol 1.E-6
#define DEFAULT_NOMINAL 10000000.

ARM_Date AddPeriod(const ARM_Date& date, char* Matu, const std::string& ccy  , bool AdjorNot /*= false*/ , qCDS_ADJ adj /*= qCredit_Default */);
static inline ARM_Date AddPeriod(const ARM_Date& date, const std::string& Matu , const std::string& ccy , bool AdjorNot, qCDS_ADJ adj)
{
	return AddPeriod(date,(char*)Matu.c_str(),ccy,AdjorNot,adj); 
}

qCDS_ADJ FromSummitAdjCDSToARM( const char* AdjCode);

class ICM_DefaultCurve;
class ARM_VolCurve; 

extern ICM_DefaultCurve** BumpPortfolio_ (const ICM_DefaultCurve ** portfolio, 
								 int nbcurves,
								 qSENSITIVITY_TYPE typesensi,
								 char* plot, 
								 char* label,
								 int& nocurve,
								 double epsilon = CREDIT_DEFAULT_VALUE);

 

class ICM_Cds; 
class ICM_Pricer;
/*********************************************************************************/
/*! \class  ICM_DefaultCurve icm_defaultcurve.h "icm_defaultcurve.h"
 *  \author 
 *	\version 1.0
 *	\date   4 may 2004
 *	\file   icm_defaultcurve.h
 *	\brief  Default Curve */
/***********************************************************************************/
class ICM_DefaultCurve : public ARM_Object
{
private:
	ARM_Date		itsAsOfDate;			//AsOfDate
	ARM_ZeroCurve*	itsZeroCurve;			// 0-1, aggregated. Discount Curve
	std::vector<std::string> itsTerms ;		// optional terms. 
	double			itsRecovery;			// Recovery Rate of the curve
	std::string		itsCurrency;			// Currency of the Curve
	std::string		itsLabel; 
	ICM_Parameters	itsParameters;			// parameters (passed to pricer for instance). 
protected:
	ARM_Vector*		itsDates;			//Vector of Dates for CDS Quotations
	ARM_Vector*		itsNbDays;			//Number of days between AsOf and itsDates
	ARM_Vector*		itsYearTerms;		//Vector of YearTerms for CDS Quotations
	ARM_Vector*		itsRates;			//Vector of spreads CDS
	ARM_Vector*		itsUpFronts;		//Vector of upfront
	ARM_Vector*		itsLambdas;			//Vector of default intensities
	ARM_Vector*		itsSurvivalProba;	//Vector of survival probability for dates itsDates
	vector<ARM_ReferenceValue*>	itsAmortRefValues;//References value

private:
	qCDS_ADJ		itsCdsAdj;			//Cds date adjustment (long CDS or not)
	bool			itsAdjBusiness;		//Cds date adjustment on business dates
	int				itsLagStartDate;	//Lag for the spot Cds quotations
	int             itsLagProtectStartDate; //Lag for the spot Cds quotations for default leg

	bool			itsIsSummitCurve;	//Summit's like calibration for intensities
	ARM_VolCurve*   itsVolCurve;		//0-1, aggreg : Volatility Curve
	bool			itsIsNameInDefault;	//Used in case of default

	// RELATED to NEW INTERPOL... 
	//		Those vectors are costructed or reset during calibration
	// 
	protected:
	// ARM_ReferenceValue	itsInterpolSearch; 
	ARM_Vector			itsInterpolYF;	
	ARM_Vector			itsInterpolLambda; 
	ARM_Vector			itsInterpolDates; 
	ARM_Vector			itsInterpolRates; 
	ARM_Vector			itsInterpolSP; 
	//
	std::string				itsCalibrationData; // temp string
	qDEFCURVE_CALIB_ALGO	itsCalibrationAlgo; 
	

	// MONTE-CARLO Hedges
	ARM_Vector*		itsDirectionalLambdasShift;			//Vector of directional shift in default intensities
 
	private :
	// -- For Internal Optimisation -----------------
	int	 its_Current_indice;
	ICM_Pricer*	its_pricer;			// association during calibration
	
 
	// -- Defines Standard for Cds Quotes -----------
	ARM_Date its_STDCDS_StartDate;
	ARM_Date its_STDCDS_ProtectionStartDate;
	int		 its_STDCDS_Frequency;
	int		 its_STDCDS_Basis;
	double	 its_STDCDS_Notional;
	std::string its_STDCDS_Currency;
	int		 its_STDCDS_Stub;
	int		 its_STDCDS_CreditLag;
	int		 its_STDCDS_Adjusted;
	qPAYMENT_PREMIUM_LEG its_STDCDS_AccOnDef;
	bool	 its_STDCDS_IncludeMaturity;
	int		 its_STDCDS_AdjustedStartDateOnBusinessDay;

 
	std::vector<ICM_Cds*>	itsCalibrationInstruments ;

    public:

	//standart constructor using tenors 1Y,3Y,...
	/** 
	ICM_DefaultCurve(const ARM_Date& asOf,
					const std::vector<std::string>& terms,
					 ARM_Vector* rates,
					 double& Recovery,
					 ARM_ZeroCurve* zc,
					 int intRule, 
					 int adjStartDate,
					 qCDS_ADJ adj ,  	
					 const std::string& ccy  ,
					 const std::string& label   ,
					 bool issummitcurve  ,
 					 const ARM_VolCurve* VolCurve  ,
					 long PayFrq ,
 					 qDEFCURVE_CALIB_ALGO	calibAlgo,
					 const std::string& calibData, 
					 int LagStartDate,
					 const ICM_Parameters& params,
					 ARM_Vector* UpFronts = NULL,
					 vector<ARM_ReferenceValue*>& AmortRefValues = vector<ARM_ReferenceValue*>());
	**/ 

	void Set(const ARM_Date& asOf,
				   const std::vector<std::string>& terms,
                   ARM_Vector* rates,
				   double Recovery,
				   ARM_ZeroCurve* zc,
					 int intRule, 
					 int adjStartDate,
				   qCDS_ADJ adj ,  	
                   const std::string& ccy  ,
				   const std::string& label /*= NULL*/,
				   bool issummitcurve /* = true*/,
				   const ARM_VolCurve* VolCurve /* = NULL*/,
					 long PayFrq /* = K_QUARTERLY*/,
 					 qDEFCURVE_CALIB_ALGO	calibAlgo, 
					 const std::string& calibData,
					 int LagStartDate,
					 const ICM_Parameters& params,
					 ARM_Vector* UpFronts = NULL,
					 vector<ARM_ReferenceValue*>& AmortRefValues = vector<ARM_ReferenceValue*>());


	void Set(const ARM_Date& asOf,
					 const ARM_Vector& Dates,
					 const ARM_Vector& Intensity,
					 double  Recovery,
					 ARM_ZeroCurve* zc,  
					 int intRule, 
					 int adjStartDate,	
					 qCDS_ADJ adj ,  
					 const std::string& ccy,
					 const std::string& label,
					 const ARM_VolCurve* VolCurve ,
 					 qDEFCURVE_CALIB_ALGO	 calibAlgo,
					 const std::string& calibData,
					 int LagStartDate,
					 const ICM_Parameters&params);

protected:
	// ICM_InterpolDefCrv 
	void Set(const ARM_Date& asOf,
					 const ARM_Vector* Dates,
					 const ARM_Vector* DefCurve,
					 double& Recovery,
					 ARM_ZeroCurve* zc,
					 int intRule, 
					 int adjStartDate,
					 qCDS_ADJ adj ,  
					 bool IsDefCurveSource,
					 const std::string& ccy  ,
					 const std::string& label /*= NULL */,
					 const ARM_VolCurve* VolCurve /*= NULL*/,
 					 qDEFCURVE_CALIB_ALGO	calibAlgo,
					 const std::string& calibMethod,
					 int LagStartDate);


protected:
	void Set(const ARM_Date& asOf,
					 double* YearTerms,
					 ARM_Vector* Intensity,
					 double& Recovery,
					 ARM_ZeroCurve* zc,  	
					 int intRule, 
					 int adjStartDate,
					 qCDS_ADJ adj ,  
					 const std::string& ccy ,
					 const std::string& label ,
					 const ARM_VolCurve* VolCurve ,
					 qDEFCURVE_CALIB_ALGO calibAlgo,
					 const std::string& calibData,
					 int LagStartDate);
private:
		/** ICM_DefaultCurve(const ARM_Date& asOf,
					 const ARM_Vector& YForDates,
					 const ARM_Vector& rates,
					 const double& Recovery,
					 ARM_ZeroCurve* zc,
					 int intRule, 
					 int adjStartDate,
					 qCDS_ADJ adj ,
					 const std::string&  ,
					 const std::string&  ,
					 bool issummitcurve  ,
					 const ARM_VolCurve* VolCurve ,
					 const bool& isdatesininput  ,
					 long PayFrq ,
 					 qDEFCURVE_CALIB_ALGO	calibAlgo,
					 const std::string& calibData,
					 int LagStartDate);
					 **/ 

public:
	void Set (const ARM_Date& asOf,
                   const ARM_Vector& YForDates,
				   const ARM_Vector& rates,
				   const double& Recovery,
				   ARM_ZeroCurve* zc,
					 int intRule, 
					 int adjStartDate,
					 qCDS_ADJ adj ,
					 const std::string& ccy  ,
					 const std::string& /* = NULL */ ,
				   bool issummitcurve /*= true*/,
				   const ARM_VolCurve* VolCurve /*= NULL*/,
				   const bool& isdatesininput /*= true*/,
				   long PayFrq /*= K_QUARTERLY*/,
 					 qDEFCURVE_CALIB_ALGO	calibAlgo,
					 const std::string& calibData,
					 int LagStartDate);

public :
    ICM_DefaultCurve(void) { Init();}

	virtual ARM_CLASS_NAME GetRootName(void)
	{
		return(ICM_DEFAULTCURVE);
	}

		

  protected:
	ICM_DefaultCurve(const ICM_DefaultCurve&ref); 
public:
 	int GetLagStartDate() const {return itsLagStartDate; }
   virtual ~ICM_DefaultCurve(void) ;

    void Init(void);

	// implémentation spécifique à la classe

	double GetPWCIntensity(const ARM_Date& date) const 
	{
		double yt = (date - itsAsOfDate)/365.;
		return GetPWCIntensity(yt); 

	}
	double GetPWCIntensity(double yf) const
	{
 		return itsInterpolLambda.Elt( position(yf) ) ;
	}
public:
	unsigned int position(const ARM_Date&date) const; 
	unsigned int position(double yf ) const; 
public:
		const ARM_Vector&  GetInterpolDates() const { return itsInterpolDates ; }
		const ARM_Vector&  GetInterpolRates() const { return itsInterpolRates ; }
		const ARM_Vector&  GetInterpolYF() const { return itsInterpolYF ; }
		const ARM_Vector&  GetInterpolSP() const { return itsInterpolSP; }
		const ARM_Vector&  GetInterpolLambda() const { return itsInterpolLambda ;}
public:
		bool GetIsSummitCurve(void) const { return itsIsSummitCurve;}
		void SetIsSummitCurve(bool value) { itsIsSummitCurve = value;}
		const ICM_Parameters& GetParameters() const { return itsParameters;}
		bool GetIsNameInDefault(void) const { return itsIsNameInDefault;}
		void SetIsNameInDefault(bool value) { itsIsNameInDefault = value;}

		double GetRecovery(void) const { return itsRecovery;}
		void SetRecovery(double value) { itsRecovery = value;}

public:
	// should not be public.
	void SetAsOfDate(const ARM_Date& AsOfDate) { itsAsOfDate = AsOfDate; };
public:

		const ARM_Date&  GetAsOfDate(void) const { return itsAsOfDate;};

		void SetZeroCurve(ARM_ZeroCurve* zc) ;
		ARM_ZeroCurve* GetZeroCurve(void) const { return itsZeroCurve;}
public:
		const std::vector<std::string>& GetTerms(void) const 
		{
			return itsTerms;
		}
		const std::string& GetTerm(unsigned int i) const 
		{
			if (i>=itsTerms.size()) 
				ICMTHROW(ERR_INVALID_ARGUMENT,"GetTerm: out of bound "<<i); 
			return itsTerms[i] ;
		}

		void SetRates(ARM_Vector* vector) ;  // adoption
 		void setRates(const ARM_Vector&rates)
		{
			*itsRates=rates; 
		}
		const ARM_Vector* GetRates(void) const  { return itsRates;}
		double GetRate(int i) const { return itsRates->Elt(i);}


		const ARM_Vector* GetLambdas(void)  const { return itsLambdas;}
		double GetLambda(int i) const { return itsLambdas->Elt(i);}

		void SetLambdas(ARM_Vector* vector) 
		{
			if (itsLambdas)
				delete itsLambdas;
			itsLambdas = vector;
		}

		unsigned long GetSize() const 
		{
			if (!itsDates) 
				ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_DefaultCurve::GetSize: no dates are defined !");
			return itsDates->GetSize(); 
		}
		const ARM_Vector* GetDates(void)  const { return itsDates;}
		double GetDate(int i) const { return itsDates->Elt(i);}
		const ARM_Vector* GetYearTerms(void)  const { return itsYearTerms;}
		double GetYearTerm(int i) const { return itsYearTerms->Elt(i);}

		void SetSurvivalProba(ARM_Vector* vector) 
		{
			if (itsSurvivalProba)
				delete itsSurvivalProba;
			itsSurvivalProba = vector;
		}

		const ARM_Vector* GetSurvivalProba(void)  const { return itsSurvivalProba;}
		const ARM_Vector* GetNbDays(void)  const { return itsNbDays;}
		double GetNbDay(const int& i) const { return itsNbDays->Elt(i);}

 
		const std::string&  GetCurrency(void) const { return itsCurrency;}

		ARM_VolCurve* GetVolCurve(void) { return itsVolCurve;}
		const ARM_VolCurve* GetVolCurve(void) const { return itsVolCurve;}
		void SetVolCurve(const ARM_VolCurve* value) ;


		double SurvivalProba(const double& yearterm) const 
		{
			return SurvivalFunction(yearterm);
		}
		virtual void CptTermsSurvivalProba(void)
		{
		    throw Exception(__LINE__, __FILE__, ERR_UNIMP_METHOD_CALL,
                       "Unimplemented <CptTermsSurvivalProba> method");

		}


		ARM_Vector* GetDirectionalLambdasShift(void)  { return itsDirectionalLambdasShift;}
		double GetDirectionalLambdasShift(int i) const { return itsDirectionalLambdasShift->Elt(i);}

		void SetDirectionalLambdasShift(ARM_Vector* vector) 
		{
			if (itsDirectionalLambdasShift)
				delete itsDirectionalLambdasShift;
			itsDirectionalLambdasShift = vector;
		}

	
		
		double SurvivalProba(const ARM_Date& date) const 
		{
			double yearterm = (date.GetJulian() - itsAsOfDate.GetJulian())/K_YEAR_LEN;

			double Sprob = SurvivalFunction(yearterm);

			return (Sprob);
		}

		double DefaultProba(const double& yearterm) const 
		{
			double DefProb = 1. - SurvivalProba(yearterm);

			return (DefProb);
		}

		double DefaultProba(const ARM_Date& date) const 
		{
			double yearterm = (date.GetJulian() - itsAsOfDate.GetJulian())/K_YEAR_LEN;

			double Dprob = DefaultProba(yearterm);
			
			return (Dprob);
		}

		const std::string& GetLabel() const { return itsLabel; }

		void SetLabel(const std::string& label)
		{
			itsLabel=label ;
			if (label.find(ISSUER_IN_DEFAULT)!=-1) {itsIsNameInDefault=true;}
		}

public:
	virtual void Calibrate()=0; 	
private:	
	void OldCalibrate()  ; 
	void Calibrate_Stress_Test_Guess_Brent() ;
public: 
 
		double Evaluate(const double& x);

		void View(char* id = NULL, FILE* ficOut = NULL );
		void View(char* id  , FILE* ficOut  ) const 
		{  const_cast<ICM_DefaultCurve*>(this)->View(id,ficOut); }

		virtual double FwdSpread(const ARM_Date& t1, const ARM_Date& t2) const ;
		double FwdSpread(const double& yt1, const double& yt2) const ;

		// LJ - behaves like an Index
		double FwdSpread_AsIndex(const ARM_Date& t1, const ARM_Date& t2, double& FlatRbp_fwd, double& Fwd_RPV01) const ;
		double FwdSpread_AsIndex(const double& yt1, const double& yt2, double& FlatRbp_fwd, double& Fwd_RPV01) const ;

		double AjustConvexity(const double& yt1, const double& yt2, const double& FwdSpread, ARM_VolCurve* vol ) const ;

		double PayLagAdjustment(const double& yt1, const double& yt2, const double& FwdSpread, ARM_VolCurve* vol )const ; //Payment Lag Adjustment

		// LJ - with a given vol value
		double AjustConvexity(const double& yt1, const double& yt2, const double& FwdSpread, double vol);
		double PayLagAdjustment(const double& yt1, const double& yt2, const double& FwdSpread, double vol ); //Payment Lag Adjustment

		double RiskyPV01(const ARM_Date& t1, const ARM_Date& t2) const ;
		double RiskyPV01(const ARM_Date& t1, const std::string& tenor) const ;
		double DefLegPV(const ARM_Date& t1, const ARM_Date& t2) const ;

		qCDS_ADJ GetCdsAdj(void) const {return itsCdsAdj;}	

	// API 
	virtual ICM_DefaultCurve* GenerateShiftCurve(const std::vector<std::string> & Terms, 
												  const ARM_Vector& epsilon ,
												  qSENSITIVITY_TYPE  ) const =0 ; 
 
	public:
	virtual double DefProbInverse(const double& PriceIn) const =0
	{
    throw Exception(__LINE__, __FILE__, ERR_UNIMP_METHOD_CALL,
                "Unimplemented <DefProbInverse> method");
	}
	double ImpliedSpreadInterpol(const std::string& plot) const ; 
	double ImpliedSpreadInterpol(const ARM_Date& date) const ; 
	double RiskyDuration(const ARM_Date& date) const ;
	double RiskyDuration(const std::string& Tenor) const ;

		

	ICM_Pricer*& Get_pricer() { return its_pricer;}

	const ARM_Date&				GetSTDCDS_StartDate()						const {return its_STDCDS_StartDate;}
	const ARM_Date&				GetSTDCDS_ProtectionStartDate()				const {return its_STDCDS_ProtectionStartDate;}
	const int&					GetSTDCDS_Frequency()						const {return its_STDCDS_Frequency;}
	const double&				GetSTDCDS_Notional()						const {return its_STDCDS_Notional;}
	const std::string&			GetSTDCDS_Currency()						const {return its_STDCDS_Currency;}
	const int&					GetSTDCDS_Stub()							const {return its_STDCDS_Stub;}
	const int&					GetSTDCDS_CreditLag()						const {return its_STDCDS_CreditLag;}
	const qPAYMENT_PREMIUM_LEG& GetSTDCDS_AccOnDef()						const {return its_STDCDS_AccOnDef;}
	const bool&					GetSTDCDS_IncludeMaturity()					const {return its_STDCDS_IncludeMaturity;}
	int							GetSTDCDS_Basis()							const {return its_STDCDS_Basis; }
	int							GetSTDCDS_Adjusted()						const {return its_STDCDS_Adjusted; }
	int							GetSTDCDS_AdjustedStartDateOnBusinessDay()	const {return  its_STDCDS_AdjustedStartDateOnBusinessDay; } ;
	ARM_Security* stdCDS(const ARM_Date& t1, const ARM_Date& t2,double spread=1.,double notional=1.,bool useimpliedstartdate=false) const ;

	virtual ICM_DefaultCurve* GenDefCurve(const ARM_Vector& dates,
										  const ARM_Vector& spreads,
										  const double& recovery,
										  const std::string& label,
										  const bool& isdates,
										  ARM_ZeroCurve* ircurve) const =0 ; // {return NULL;}

	ICM_DefaultCurve* createDefCurveFromBase(const ICM_DefaultCurve* DefCurveIndex , const ARM_Vector& vBase) const ;
	int GetPayFrequency(void) const {return its_STDCDS_Frequency;}

	double GetNotional(void) const {return its_STDCDS_Notional;}

	double NPV_Implicit_Cds(ARM_Date& Maturity,
						  const string& Tenor = "NULL",	
						  const qCMPMETH& mode = qCMPPRICE) const ;

	double RiskyPV01AsSensitivity(const std::string& Tenor) const ;

	const std::string& GetCalibrationData() const { return itsCalibrationData; }
	const qDEFCURVE_CALIB_ALGO GetCalibrationAlgo() const { return itsCalibrationAlgo; }

	double ImpliedSpreadInterpolOLD(const std::string&plot , 
								const double& slope =0.,
								const double& ExactDate = 0.,
								const double& InterDate =0);
	
	ICM_DefaultCurve* createFlatCurve(const std::string& tenor) const; 
	ICM_DefaultCurve* createFlatCurve(const ARM_Date& date) const; 
	
	void SetAmortRefValues(vector<ARM_ReferenceValue*>&	AmortRefValues);//References value
	void SetUpFronts(ARM_Vector* UpFronts);		//Vector of upfront
	
	vector<ARM_ReferenceValue*>& GetAmortRefValues() { return itsAmortRefValues;}//References value
	ARM_Vector* GetUpFronts() {return itsUpFronts;}		//Vector of upfront


protected:
	virtual void ResetLambda(const int& indice,const double& lambda);
 
	int GetCurrent_indice(void) const { return  its_Current_indice;}
	void SetCurrent_indice(int i) { its_Current_indice=i; } 
	void SetLagStartDate(int lag) {itsLagStartDate = lag; }
	//
	void SetInterpolSP(unsigned int i,double v) { itsInterpolSP.Elt(i)=v;}
	void SetInterpolLambda(unsigned int i,double v) { itsInterpolLambda.Elt(i)=v;}
	void SetTerms(const std::vector<std::string>&terms) 
		{
			itsTerms=terms ;
		}; 
	void SetYearTerms(ARM_Vector* vector) 
	{
		if (itsYearTerms)
			delete itsYearTerms;
		itsYearTerms = vector;
	}
	double GetSurvivalProba(const int& i) const 
	{ 
		return itsSurvivalProba->Elt(i);
	}

private:
	void CptSetCreditMarketData(bool charTermKn);
	void AdjustDateYearsTerms(bool dateVectorKn);
	virtual double SurvivalFunction(const double& yearterm) const =0; 
private:
	ICM_DefaultCurve& operator=(const ICM_DefaultCurve&); 
};		

// 
//	------------------------------------------------------------------------------------------------
inline unsigned int 
ICM_DefaultCurve::position(const ARM_Date&date) const
{
	double yf = (date.GetJulian() - itsAsOfDate.GetJulian())/365.;
	return position(yf); 
}
//	------------------------------------------------------------------------------------------------
//
//		i such as : 
//		yf(ti) <= yf < yf(ti+1) 
inline unsigned int 
ICM_DefaultCurve::position(double yf) const
{
	static int prevLocation=0; 
	int i=locateIndex(&unconst(itsInterpolYF),yf,prevLocation); 
	if (i==-1) return 0; 
	return i; 
}
#endif /*---- End of file ----*/
