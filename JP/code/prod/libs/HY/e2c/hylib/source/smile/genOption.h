// -----------------------------------------------------------------------
// This proprietary software has been developed strictly for J.P. Morgan's
// internal use.  Any use or misuse, intentional or otherwise, which
// contradicts or places this policy in jeopardy is strictly forbidden.
//
// Copyright 1999 J.P. Morgan & Co. Incorporated.  All rights reserved.
// -----------------------------------------------------------------------
//
// 1.0 10/19/99 Neil Yang
// 
//

#if ! defined(_GENOPTION_)
#define _GENOPTION_
#include "genfunc.h"
#include "gtobf.h"   // GtoOPTION_CALL
#include "macros.h"
#include "genIntegral.h"


const double  INTEGRAL_LIMIT = 10.0;

class BaseOption : public CM::Object
{
	double m_strike;
	double m_yearsToExpiry;
	double m_yearsToPayment;

public:
	BaseOption(): m_strike(0),m_yearsToExpiry(0),m_yearsToPayment(0) {}
	BaseOption(double strike, double yearsToExpiry, double yearsToPayment):
		m_strike(strike),
		m_yearsToExpiry(yearsToExpiry),
		m_yearsToPayment(yearsToPayment) {}

		void set_strike(double strike) {m_strike = strike;}
		void set_yearsToExpiry(double yearsToExpiry) {m_yearsToExpiry = yearsToExpiry;}
		void set_yearsToPayment(double yearsToPayment) {m_yearsToPayment = yearsToPayment;}
		double get_strike() const { return m_strike;}
		double get_yearsToExpiry() const {return m_yearsToExpiry;}
		double get_yearsToPayment() const {return m_yearsToPayment;}

		virtual double payoff(double y) const = 0;
};


class CallOption : public BaseOption
{

public:
	CallOption(double strike, double yearsToExpiry, double yearsToPayment)
		: BaseOption(strike, yearsToExpiry, yearsToPayment) {}

	double payoff(double y) const 
	{
		return MAX(y-get_strike(),0);
	}

};

class PutOption : public BaseOption
{

public:
	PutOption(double strike, double yearsToExpiry, double yearsToPayment)
		: BaseOption(strike, yearsToExpiry, yearsToPayment) {}

	double payoff(double y) const
	{
		return MAX(get_strike()-y,0);
	}

};


class Payoff
{
	BaseOption* m_baseOption;
	BaseFunction* m_processFunc;

public:
	Payoff(BaseOption* baseOption, BaseFunction* processFunc): m_baseOption(baseOption),
	m_processFunc(processFunc) {}

	double operator()(double x) const
	{
	
		return m_baseOption->payoff( (*m_processFunc)(x) )*GtoNormalDen(x);
	
	}

	BaseFunction* get_baseFunc() { return m_processFunc;}
};




class ForwardSolverFunc
{
	double m_fwd;
	Payoff* m_payoff;

public:
	ForwardSolverFunc( double fwd, Payoff* payoff) : m_fwd(fwd), m_payoff(payoff) {}
	double operator() (double x) const 
	{
		BaseFunction* process = m_payoff->get_baseFunc();
		(process)->set_amplitude(x);

		double temp =Integral(*m_payoff,         
						 -INTEGRAL_LIMIT,          
						 INTEGRAL_LIMIT,           
						 MAX_INTEGRAL_ERROR);
		return (temp - m_fwd);
	}
};


double GenOptionPricer(double  forward,        /* Forward price */
					   BaseOption* baseOption, /* Option infor */
					   double  volatility,     /* */
					   double  discountRate,   /* */
					   BaseFunction*  xx);


#endif
