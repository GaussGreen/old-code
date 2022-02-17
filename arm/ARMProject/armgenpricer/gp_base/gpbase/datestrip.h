/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file datestrip.h
 *	\brief Date strip is an object that can generate a strip of dates
 *		Handles only dates and compute 
 *			-the reset dates, 
 *			-the payment dates
 *			-the interest terms
 *
 *		This object is a very simple object that can be used in order
 *		to compute the date ... Use rather datestrip than swap leg
 *		if you are only interested in dates generations!
 *
 *
 *	\author  Eric Benhamou
 *	\version 1.0
 *	\date August 2003
 */

#ifndef _INGPBASE_DATESTRIP_H
#define _INGPBASE_DATESTRIP_H

#include "port.h"
#include "env.h"
#include "rootobject.h"
#include "gplinalgtypedef.h"
#include "gpvector.h"
#include "assignop.h"
#include "cloneutilityfunc.h"

/// kernel headers
#include "dates.h"

CC_BEGIN_NAMESPACE( ARM )

class ARM_DateStrip : public ARM_RootObject
{
protected:
    std::vector<double>* itsFlowStartDates;	/// Flow start dates and Flow 
										/// end dates are used to compute the
    std::vector<double>* itsFlowEndDates;		/// period of interest 

    std::vector<double>* itsFwdRateStartDates;/// start dates and end dates
    std::vector<double>* itsFwdRateEndDates;	/// are used to compute 
										/// the forward rates
	std::vector<double>* itsResetDates;		/// resetDates
	std::vector<double>* itsPaymentDates;		/// paymentDates
    std::vector<double>* itsInterestDays;		/// numbers of days between 2 periods
    std::vector<double>* itsInterestTerms;	/// interest term... conversion of InterestDays 
										/// to the correctDayCount

	/// to keep track of the various arguments
	ARM_Date itsStartDate;				/// start date
	ARM_Date itsEndDate;				/// end date
	ARM_Date itsRefDate;				/// reference date.. in most cases it is the startDate
	int itsResetFreq;					/// reset frequency of the strip
	int itsDayCount;					/// dayCount method used for computing the accrued
	char* itsResetCalendar;				/// calendar used for reset
	int itsFwdRule;						/// whether fwds are with adjusted dates
	int itsIntRule;						/// whether interest itsare K_ADJUSTED
	int itsStubRule;					/// ability to have shortstart .... 
	int itsResetGap;					/// reset itsgap
	int itsPayFreq;						/// pay frequency of the strip
	int itsPayGap;						/// pay gap
	char* itsPayCalendar;				/// calendar used for payment
	int itsResetTiming;					/// whether resets are in arrears or in advance
	int itsPayTiming;					/// whether payments are in arrears or in advance
	int itsAdjFirstdate;				/// adjust the first date to business day or not
	int itsFirstDateFwdRule;			/// rule used to adjust the first date (default : itsFwdRule)
	bool itsIsBuiltFromOutside;			/// flag to tell whether this has been built from outside!
    int itsStdSpotDays;                 /// standard gap between reset and forward start dates
    int itsIndexFreq;                   /// term of the forward
    int itsAccruedOrFull;                /// Accrud mode or full mode to generat FwdDates    

	std::vector<double>* GenerateDatesFromTimingAndGap(
		int timing,
		int gap,
		const char* calendar,
		std::vector<double>* FlowStartDates,
		std::vector<double>* FlowEndDates ) const;
	
	void Init();

	///  compute the various dates
	void CptCashFlowDates();

	/// GetDefaultCalendar
	char* GetCalendar( const char* calendar );

	/// validation
	void Validate();

	/// function to account for a reference date
	std::vector<double>* ARM_DateStrip::AdjForRefDate( std::vector<double>* FlowDates, char* calendar );

	void CheckSizeAndSorted( size_t defaultSize, std::vector<double>* vec, const char* name );
	void CheckSize( size_t defaultSize, std::vector<double>* vec, const char* name );

	std::vector<double>* CptStripOfStartDates(
		ARM_Date& StartDate, ARM_Date& EndDate, 
		int frequency, int fwdRule, int TypeStub, int intRule, 
		char* ccy, int adjFirstdate);

public:

	/// standard constructor 
	ARM_DateStrip( 
		const ARM_Date& startDate,			/// startDate of the strip
		const ARM_Date& endDate,			/// end Date
		int resetFreq,						/// reset frequency
											/// no default to force the user to give one
		int dayCount,						/// dayCount method used for computing the accrued
											/// no default to force the user to give one
		const char* resetCalendar = GETDEFAULTVALUESTR,
											/// calendar used for reset... Default is the one from the ARM_DEFAULT_CURRENCY
		int fwdRule		= K_MOD_FOLLOWING,	/// whether fwds are with adjusted dates
		int intRule		= K_MOD_FOLLOWING,	/// whether interest are K_ADJUSTED
		int stubRule	= K_SHORTSTART,		/// ability to have K_SHORTSTART etc
		int resetGap	= GETDEFAULTVALUE,	/// reset gap default is -spotDay of default currency
		int payFreq		= GETDEFAULTVALUE,	/// payment frequency
		int payGap		= GETDEFAULTVALUE,	/// pay gap default is 0
		const char* payCalendar = GETDEFAULTVALUESTR,
											/// calendar used for payment
		int resetTiming	= K_ADVANCE,		/// whether reset are in arrears or in advance
		int payTiming	= K_ARREARS,		/// whether payment are in arrears or in advance
		int adjFirstdate= true,				/// adjust the first date to business day
		const char* refDateChar= GETDEFAULTVALUESTR,
											//// reference Date
        int stdSpotDays = GETDEFAULTVALUE,  /// standard spot days default is spotDay of default currency
        int indexTerm   = GETDEFAULTVALUE,
        int accruedOrFull = K_ACCRUED,
		int firstDateFwdRule = 9999
	);

	ARM_DateStrip( 
		std::vector<double>* FlowStartDates,		/// Flow start dates and Flow 
		std::vector<double>* FlowEndDates,		/// end dates are used to compute the	
		std::vector<double>* FwdRateStartDates,	/// fwd start dates
		std::vector<double>* FwdRateEndDates,		/// fwd end dates
		std::vector<double>* ResetDates,		/// resetDates
		std::vector<double>* PaymentDates,		/// paymentDates
		std::vector<double>* InterestDays,		/// numbers of days between 2 periods
		std::vector<double>* InterestTerms );	/// interest term... conversion of InterestDays

    ARM_DateStrip(); // to be used with extreme caution !!

	ARM_DateStrip(const ARM_DateStrip& rhs );
	ASSIGN_OPERATOR(ARM_DateStrip)
	virtual ~ARM_DateStrip();

	/// Standard ARM_Object support
	virtual ARM_Object* Clone() const;
    virtual string toString(const string& indent="",const string& nextIndent="") const;

	/// DateStrip is a root object
	ARM_CLASS_NAME GetRootName() { return ARM_DATESTRIP; }
	std::vector<double>* GetMemberData( int type );

	size_t inline size() const { return itsResetDates->size();}

	/// test if there is input in it!
	inline bool isBuiltFromOutside() const { return itsIsBuiltFromOutside; }

	/// various accessors
	/// we do not clone the object
	/// therefore if you need to, DO IT
	/// I repeat DO IT!
    inline std::vector<double>* GetFlowStartDates() const {return itsFlowStartDates;}
    inline std::vector<double>* GetFlowEndDates() const {return itsFlowEndDates;}
    inline std::vector<double>* GetFwdRateStartDates() const {return itsFwdRateStartDates;}
    inline std::vector<double>* GetFwdRateEndDates() const {return itsFwdRateEndDates;}
	inline std::vector<double>* GetResetDates() const {return itsResetDates;}
	inline std::vector<double>* GetPaymentDates() const {return itsPaymentDates;}
    inline std::vector<double>* GetInterestDays() const {return itsInterestDays;} 
    inline std::vector<double>* GetInterestTerms() const {return itsInterestTerms;} 
	inline const int&	  GetDayCount() const {return itsDayCount;}
	inline ARM_Date		  GetEndDate() const { return itsEndDate;}
	inline ARM_Date		  GetStartDate() const { return itsStartDate;}

    inline char* GetPayCalendar  () const { return(itsPayCalendar);		};
	inline char* GetResetCalendar() const { return(itsResetCalendar);	};


    /// set directly the std::vector<double> pointor with cloning it
	/// it has been a design decision to clone dates .. 
	/// use the noclone version if required
	void SetFlowStartDates( std::vector<double>* FlowStartDates );
	void SetFlowEndDates( std::vector<double>* FlowEndDates );
	void SetFwdRateStartDates( std::vector<double>* FwdRateStartDates );
	void SetFwdRateEndDates( std::vector<double>* FwdRateEndDates);
	void SetResetDates( std::vector<double>* ResetDates );
	void SetPaymentDates( std::vector<double>* PaymentDates );
	void SetInterestDays( std::vector<double>* InterestDays );
	void SetInterestTerms( std::vector<double>* InterestTerms );

	/// version with NoClone
	void SetFlowStartDatesNoClone( std::vector<double>* FlowStartDates );
	void SetFlowEndDatesNoClone( std::vector<double>* FlowEndDates );
	void SetFwdRateStartDatesNoClone( std::vector<double>* FwdRateStartDates );
	void SetFwdRateEndDatesNoClone( std::vector<double>* FwdRateEndDates);
	void SetResetDatesNoClone( std::vector<double>* ResetDates );
	void SetPaymentDatesNoClone( std::vector<double>* PaymentDates );
	void SetInterestDaysNoClone( std::vector<double>* InterestDays );
	void SetInterestTermsNoClone( std::vector<double>* InterestTerms );

	int GetResetFreq()	{ return itsResetFreq;	};
	int	GetPayFreq	()	{ return itsPayFreq;	};
	int GetPayGap	()	{ return itsPayGap;		};
	int GetResetGap	()	{ return itsResetGap;	};

	void InsertDate(
		int idx,
		double flowStartDate,
		double flowEndDate,
		double fwdRateStartDate,
		double fwdRateEndDate,
		double resetDate,
		double paymentDate,
		double interestDays,
		double interestTerm);

	void fill(const ARM_DateStrip& datestrip);


	/// data to put for blank data
	static const double DateStripCombiner_BlankData;

	/// resize and built a new object

	void ResizeAndBuilt( size_t begin, size_t end);
};

CC_END_NAMESPACE()

#endif
/*----------------------------------------------------------------------------*/
/*---- End of file ----*/
