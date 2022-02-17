

//----------------------------------------------------------------------------
//
//   Group       : QR cross asset
//
//   Filename    : SRMSwapUtil.hpp
//
//   Description : swap definition, for use with Libor model mainly
//
//   Author      : Henrik Rasmussen
//
//   Date        : 16 Dec 2005
//
//
//----------------------------------------------------------------------------

#ifndef SRM_SWAP_UTIL_HPP
#define SRM_SWAP_UTIL_HPP

#include "edginc/config.hpp"
#include "edginc/VirtualDestructorBase.hpp"
#include "edginc/YieldCurve.hpp"
#include "edginc/DayCountConvention.hpp"

DRLIB_BEGIN_NAMESPACE

class SRMSwapClass : public virtual VirtualDestructorBase 
{

public:
      
	  SRMSwapClass(void);

	  SRMSwapClass(
					DateTime             startDate,      // (I) Date instrument begins at
                    bool                 stubAtEnd,
                    DayCountConvention&  dcc,
					int                  accrInMonths,
					int                  tenorInMonths);

	  ~SRMSwapClass(void);

      SRMSwapClass& operator=(const SRMSwapClass &otherSwap);
      SRMSwapClass(const SRMSwapClass &otherSwap);

	  // class is copyable, compiler generated assignment operator
	  // and copy ctor work fine

	  double parSwapRate(const IYieldCurve& zc) const;

	  // sets weights w_{i}(0)L_{i}(0)/Sr(0) which are used in
	  // Black formula for swaption implied volatility
	  void setWeights(const IYieldCurve& zc);

      // returns weights w_{i}(0)L_{i}(0)/Sr(0) which are used in
	  // Black formula for swaption implied volatility
	  const vector<double>& getWeights(void) const;

	  const DateTime& getStartDate(void) const;
      const DateTime& getMatDate(void) const;
	  const DateTimeArray& getResetDates(void) const;
      const DateTimeArray& getPayDates(void) const;

private:

	   vector<double>  
			   m_weights;  // stores w_{i}(0)L_{i}(0)/Sr(0) for all libors 
                
		DateTime             startDate;      // (I) Date instrument begins at
        DateTime             maturityDate;   // (I) Date instrument matures at
        bool                 stubAtEnd;
		string               period;               // e.g. Y, M, W, D
		int                  tenorInMonths,
		                     accrInMonths,
		                     numPeriods;
        DayCountConventionSP dcc;
        DateTimeArraySP      resetDates;
        DateTimeArraySP      payDates; 
};

DECLARE(SRMSwapClass);

DRLIB_END_NAMESPACE

#endif // SRM_SWAP_UTIL_HPP
