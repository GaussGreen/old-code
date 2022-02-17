/*! 
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file infleg.h
 *  \brief implements the ARM_InfLeg class, a class to deal with SwapLegs
 *		Inherit from ARM_SwapLeg
 *
 *	\author  Eric Benhamou
 *	\version 1.0
 *	\date August 2003
 */

#ifndef _INGPINFLATION_INFLEG_H
#define _INGPINFLATION_INFLEG_H

/// gpbase
#include "gpbase/env.h"
#include "gpbase/port.h"
#include "gpbase/gplinalgtypedef.h"

/// kernel
#include <inst/swapleg.h>

class ARM_Date;


/*!
 * \class ARM_InfLeg
 * \author  Eric Benhamou
 * \version 1.0
 * \date August 2003
 *
 * \brief General Inflation leg 
 * that can prices anything of the type
 * multiple * CPI( Tnum_i)/CPI(Tdenom_i) + constant
 * with 2 strip of dates Tnum and Tdenom
 *
 *
 * The design choice was to
 * have a generic inflation leg inherited from ARM_Swap_Leg
 * with each inflation swap leg
 * defined as inherited object of the ARM_InfLeg
 * and that are creating a specific form of the ARM_InfLeg
 *
 *
 * 1) a Year to Year Inflation leg is a specialised 
 * version of an inflation leg that pays
 * leverage * ( CPI(Ti+1)/CPI(Ti)-1 )+ spread
 *
 * therefore:
 *		- TdenomDates   = {T0,T1,...,Tn-1}
 *		- TnumDates		= {T1,T2,....,Tn}
 *		- multiple		= leverage
 *		- constant		= -leverage + spread
 *
 *
 * 2) a Zero Coupon Inflation leg is a specialised 
 * version of an inflation leg that pays
 * leverage * ( CPI(Tn)/CPI(T0)-1 ) + spread
 *
 * therefore:
 *		- TdenomDates   = {T0}
 *		- TnumDates		= {Tn}
 *		- multiple =  leverage,
 *		- constant = -leverage + spread
 *
 *
 * 3) an OTA type Inflation leg is a specialised 
 * version of an inflation leg that pays
 * coupon * CPI(Ti+1)/CPI(T0) + spread
 *
 * therefore:
 *		- TdenomDates   = {T0,T0,... ,T0}
 *		- TnumDates		= {T1,T2,....,Tn}
 *		- multiple =  coupon,
 *		- constant =  spread
 *
 *
 * Potential types are
 *		- K_YEARTOYEAR_LEG
 *		- K_OATTYPE_LEG
 *		- K_ZEROCOUPON_LEG
 *		- K_GENERICINF_LEG
 *
 * For the notional
 * K_NX_NONE		no notional paid at the end
 * K_NX_END			100 paid at the end
 * K_NX_INFEND		CPI(Tn)/CPI(T0) * 100
 * K_NX_ASINFEND	(CPI(Tn)/CPI(T0)-1.0) * 100
 *
 *
 */

CC_BEGIN_NAMESPACE( ARM )

/// forward declaration in namespace ARM
class ARM_DateStrip;
class ARM_InfIdx;

class ARM_AverageSwapLeg : public ARM_SwapLeg{

public:
		ARM_AverageSwapLeg(	ARM_Date&		startDate, 
							ARM_Date&		endDate, 
							ARM_IRIndex*	irIndex,
							int				rcvOrPay, 
							double			spread, 
							int				stub, 
							int				decompFreq,
							ARM_Currency*	ccyName, 
							int				dayCount,
							char*			payCalName,
							int				decompPricingFlag,
							int				nxChange,
							char*			refDate,
							int				adjStartDate,
							double			couru);

		void Set	(		ARM_Date&		startDate, 
							ARM_Date&		endDate, 
							ARM_IRIndex*	irIndex, 
							int				rcvOrPay, 
							double			spread, 
							int				stub, 
							int				decompFreq,
							ARM_Currency*	ccy, 
							int				dayCount,
							char*			payCalName,
							int				decompPricingFlag,
							int				nxChange,
							char*			refDate,
							int				adjStartDate);


        ARM_AverageSwapLeg	(const ARM_AverageSwapLeg& swapLeg):ARM_SwapLeg(swapLeg){ itsCouru=swapLeg.itsCouru;}

		virtual ARM_Object* Clone(){	return new  ARM_AverageSwapLeg(*this);	}

        ARM_AverageSwapLeg  operator = (const ARM_AverageSwapLeg& swapLeg) { return  ARM_AverageSwapLeg(swapLeg); }

		virtual ~ARM_AverageSwapLeg(){}


	inline double GetCouru() const{ return itsCouru; }
	inline void		SetCouru( const double & couru){ itsCouru = couru; }

private:
	double itsCouru;

};


class ARM_InfLeg : public ARM_SwapLeg 
{
private:
	/// whether it is a generic leg
	/// a year to year, a zero coupon or an oat-type
	int itsSwapType;
	/// whether we want this to be CPILINEAR
	int itsInterpType; 

	int itsResetNumGaporRollDate;
	
	int	itsResetDenomGaporRollDate;

	/// multiple (see in the derived class to understand how the 
	/// inflation leg is built
	double itsMultiple;
	
	/// multiple ( applied to the other index mixed with the inflation index 
	double itsCoMultiple;

	/// constant of the inflation leg
	double itsConstant;					

	/// for zeroCoupon and OAT type, the concept of reset is fundamental
	double itsFirstReset;
	
	///for basculate to livret A Schedule
	bool isLivretA;

	/// numerator and denominator reset dates
	ARM_GP_Vector* itsNumResetDates;
    ARM_GP_Vector* itsDenomResetDates;

	ARM_GP_Vector* itsNumCPIRates;
	ARM_GP_Vector* itsDenomCPIRates;

	/// Stores the num Publication Dates
	ARM_GP_Vector* itsNumPublicationDates;

	ARM_GP_Vector* itsDiscountFactors;

	/// copy mehtod that only copies the member data above
	void CopyNoCleanUp( const ARM_InfLeg& infLeg );
	
	void CleanUp();
	
	void Init(); // this function is never virtual since it belongs to the object itself

	void Set(
        const ARM_InfIdx* infIdx, 
		int rcvOrPay,
		int finalNotionalType,
		const ARM_Currency*	discountCcy,
		const ARM_DateStrip* numStripDate,
		const ARM_DateStrip* denomStripDate );

	string itsPayCurrencyName;
	string itsIndexName;


public:
	/// 1) General constructor with 2 dateStrips
	/// ARM_DateStrip gives a very condensed form of the constructor
	/// this should be often preferred
    ARM_InfLeg(
        const string& indexName,
		int rcvOrPay,
		int itsInterpType,
		double multiple,
        double constant,
		int finalNotionalType,
		const ARM_DateStrip* numDateStrip,
		const ARM_DateStrip* denomDateStrip,
		const ARM_Currency* discountCcy = NULL );

	/// 2) constructor with full precision of arguments for the date strip 
	/// and a swap Type
	ARM_InfLeg(
		const ARM_Date& startDate,
		const ARM_Date& endDate,
		const string& indexName,
		int swapType,
		int rcvOrPay,
		int interpType					= K_CPILINEAR,
		double multiple					= 1.0,
		double Comultiple				= 1.0,
        double constant					= 0.0, 
		int resetFreq					= K_DEF_FREQ,
		int dayCount					= KACTUAL_ACTUAL,
		const char* resetCalendar		= "INF",				/// calendar used for reset dates
		int fwdRule						= K_MOD_FOLLOWING,		/// whether fwds are with adjusted dates
		int intRule						= K_UNADJUSTED,			/// whether period of interest are adjusted or not
		int stubRule					= K_SHORTSTART,			/// ability to have K_SHORTSTART etc
		int resetNumGaporRollDate		= 0.0,					/// reset gap default is 0 (can be either a gap or a roll date)
		int resetDenomGaporRollDate		= 0.0,					/// reset gap default is 0 (can be either a gap or a roll date)
		int payFreq						= GETDEFAULTVALUE,		/// payment frequency
		int payGap						= GETDEFAULTVALUE,		/// pay gap default is 0
		const char* payCalendar			= GETDEFAULTVALUESTR,	/// calendar used for payment
		int adjFirstDate				= 1,					/// we adjust for first date
		int finalNotionalType			= K_NX_NONE,
		double firstReset				= GETDEFAULTVALUE,		/// ability to overwritte the first reset
		const ARM_Currency*	discountCcy = NULL );

	/// copy constructor
	ARM_InfLeg( const ARM_InfLeg& infLeg);
	
	/// assignment operator
	ARM_InfLeg& operator = ( const ARM_InfLeg& infLeg );
	
	virtual ~ARM_InfLeg();

	/// methods to overwritte for pricing
	virtual double ComputePrice(int mode = 0);		
	virtual void PropagateModel(ARM_Model *model);
    virtual void CptExpectedFwdRates();		
    virtual void CptCashFlowValues();

	/// ARM_Object standard support
	virtual void View(char* id = NULL, FILE* ficOut = NULL);
	virtual ARM_Object* Clone();

	ARM_GP_Vector* GetDiscountFactors() const;
	ARM_GP_Vector* GetNumResetDates() const;
    ARM_GP_Vector* GetDenomResetDates() const;
	ARM_GP_Vector* GetNumCPIRates() const;
	ARM_GP_Vector* GetDenomCPIRates() const;
	
	ARM_GP_Vector* GetNumPublicationDates() const { return itsNumPublicationDates; }
	

	inline string& getPayCurrencyName() { return itsPayCurrencyName; }
	inline string& getIndexName() { return itsIndexName; }
	inline void setPayCurrencyName( const string& currencyName ) { itsPayCurrencyName = currencyName; }
	inline void setIndexName( const string& indexName ) { itsIndexName = indexName; }
	
	int GetSwapType() const;
	inline double GetMultiple() const { return itsMultiple; }
	inline double GetCoMultiple() const { return itsCoMultiple; }
	inline double GetConstant() const { return itsConstant; }
	inline void SetConstant( double constant ) { itsConstant=constant; }
	inline int	GetNumResetGap(){ return itsResetNumGaporRollDate; }
	inline int	GetDemResetGap(){ return itsResetDenomGaporRollDate; }
	inline int	GetInterpType(){ return itsInterpType; }
	inline double	GetFirstReset(){ return itsFirstReset;}
	/// function to convert a date to a gap
	static int GetGapFromGapOrDate( int gapInput, const char* calendar, const ARM_Date& refDate );

};


CC_END_NAMESPACE()

#endif
/*--------------------------------------------------------------------------*/
/*---- End of file ----*/


