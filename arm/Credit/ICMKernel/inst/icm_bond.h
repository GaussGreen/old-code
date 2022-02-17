
#if !defined(_ICM_BOND_H__)
#define _ICM_BOND_H__

#include "ICMKernel\inst\icm_leg.h"

/*********************************************************************************/
/*! \class  ICM_Bond icm_bond.h "icm_bond.h"
 *  \author  
 *	\version 1.0
 *	\date   June 2003
 *	\file   icm_bond.h
 *	\brief  Definition of an <B> Fixed ICM_bond </B> */
/***********************************************************************************/

class ICM_Bond : public ICM_Leg  
{

public:


	ICM_Bond() : ICM_Leg()
	{}

	virtual ~ICM_Bond()
	{}

	void BitwiseCopy(const ARM_Object* srcleg);

	void Copy(const ARM_Object* srcleg);

	ARM_Object* Clone(void);

	ICM_Bond(double CouponRate,
			 const ARM_Date& Int_Accrual_Date,
			 const ARM_Date& Maturity, 
			 const ARM_Date* First_Period_Reference_Date,
			 int	  frequency,
			 double NotionalAmount /* = 100.*/ , 
			 qPAYMENT_PREMIUM_LEG AccruedOnDefault /*= qACCRUED_SETTLED */,
			 const std::string& discountCcy /* = ARM_DEFAULT_CURRENCY*/ ,
			 const std::string& PayCalendar /* = NULL */ ,
 			 int	DayCount, /* = K30_360, /*< Si non affecté il sera attribué par défaut*/
			 int	AccruedDayCount, /* = -1000,   Si non affecté il sera attribué par défaut*/
			 int    SettlementGap,   /* = -1000,     Si non affecté il sera attribué par défaut de la devise */ 
			 int	stubrule /* = K_SHORTSTART*/ ,
			 double RedemptionValue /* = 100. */ ) : ICM_Leg() 
		{
			Init();
			
			Set( CouponRate, Int_Accrual_Date, Maturity, First_Period_Reference_Date, frequency,
			 NotionalAmount, AccruedOnDefault , discountCcy , PayCalendar, DayCount , 
			 AccruedDayCount , SettlementGap , stubrule, RedemptionValue );

		}


	void Set(double CouponRate,
			 const ARM_Date& Int_Accrual_Date,
			 const ARM_Date& Maturity, 
			 const ARM_Date* First_Period_Reference_Date,
			 int	  frequency,
			 double NotionalAmount /* = 100. */ , 
			 qPAYMENT_PREMIUM_LEG AccruedOnDefault /* = qACCRUED_SETTLED */,
			 const std::string& discountCcy /*= ARM_DEFAULT_CURRENCY */,
			 const std::string& /* char* PayCalendar = NULL*/ ,
 			 int	DayCount ,/* = K30_360, /*< Si non affecté il sera attribué par défaut*/
			 int	AccruedDayCount, /* = -1000, /*< Si non affecté il sera attribué par défaut*/
			 int    SettlementGap,   /* = -1000, /*<  Si non affecté il sera attribué par défaut de la devise */ 
			 int	stubrule /* = K_SHORTSTART*/,
			 double redemption /* =100. */);

	ICM_Bond(double CouponRate,
			 const ARM_Date& Int_Accrual_Date,
			 const ARM_Date& Maturity,
			 const ARM_Date* First_Period_Reference_Date,
			 int	  frequency,
			 double NotionalAmount /* = 100 */ , 
			 qPAYMENT_PREMIUM_LEG AccruedOnDefault /* = qACCRUED_SETTLED*/ ,
			 const std::string& discountCcy /*= ARM_DEFAULT_CURRENCY */,
			 const std::string& payCalendar /*= NULL */,
			 int	DayCount ,/* = K30_360, /*< Si non affecté il sera attribué par défaut*/
			 int	AccruedDayCount ,/* = -1000, /*< Si non affecté il sera attribué par défaut de la devise */
			 int    SettlementGap,   /* = -1000*/
			 int	stubrule ,/* = K_SHORTSTART*/ 
			 ARM_ReferenceValue* NotionalSchedule /* = NULL*/ ) : ICM_Leg()
		{

			Init();

			Set(CouponRate,
			 Int_Accrual_Date,
			 Maturity, 
			 First_Period_Reference_Date,
			 frequency,
			 NotionalAmount , 
			 AccruedOnDefault,
			 discountCcy ,
			 payCalendar,
			 DayCount , 
			 AccruedDayCount , 
			 SettlementGap ,
			 stubrule ,
			 NotionalSchedule);
		}	
		void Set(double CouponRate,
			 const ARM_Date& Int_Accrual_Date,
			 const ARM_Date& Maturity, 
			 const ARM_Date* First_Period_Reference_Date,
			 int	  frequency,
			 double NotionalAmount /*= 100*/, 
			 qPAYMENT_PREMIUM_LEG AccruedOnDefault /*= qACCRUED_SETTLED*/,
			 const std::string& discountCcy /*= ARM_DEFAULT_CURRENCY*/,
			 const std::string& payCalendar /*= NULL*/,
			 int	DayCount, /*= K30_360, /*< Si non affecté il sera attribué par défaut*/
			 int	AccruedDayCount, /*= -1000, /*< Si non affecté il sera attribué par défaut de la devise */
			 int    SettlementGap   /*= -1000*/,
			 int	stubrule, /*= K_SHORTSTART, /*<  Si non affecté il sera attribué par défaut de la devise */ 
			 ARM_ReferenceValue* /*NotionalSchedule = NULL*/ ); 

/* ****************************************************************************************************************** */
/*! \fn PriceToYield(ARM_Date& settlement, double price) 
	\brief Computes the <B> Yield to maturity (in percentage points) </B> at settlement date for a given <B> clean price (in units : Base 1) </B>
	\param settlement settlement date
	\param price Clean Price at settlement date
	\return Yield to maturity (in percentage points) */
/* ****************************************************************************************************************** */
	double PriceToYield(ARM_Date& settlement, double price);

/* ****************************************************************************************************************** */
/*! \fn YieldToPrice(ARM_Date& settlement, double yield) 
	\brief Computes the <B> clean price (in units : Base 1) </B> at settlement date for a given <B> Yield to maturity (in percentage points) </B>
	\param settlement settlement date
	\param yield Yield to maturity (in percentage points) at settlement date
	\return Clean price (in units : Base 1) */
/* ****************************************************************************************************************** */
	double YieldToPrice(ARM_Date& settlement, double yield);

/* ****************************************************************************************************************** */
/*! \fn Accrued(ARM_Date & settlement) 
	\brief Computes the <B> Accrued coupon (in units : Base 1) </B> at settlement date
	\param settlement settlement date
	\return Accrued coupon (in units : Base 1) */
/* ****************************************************************************************************************** */
	virtual double Accrued(ARM_Date & settlement);

private :

	void Init(void);
	virtual void CptCashFlowValues(void);
	double CouponPeriod(int index);

	
};

#endif // !defined(AFX_ICM_BOND_H__705C12EA_BBC1_4D17_AA0E_664E45678280__INCLUDED_)
