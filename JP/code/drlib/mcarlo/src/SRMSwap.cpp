
//----------------------------------------------------------------------------
//
//   Group       : QR cross asset
//
//   Filename    : SRMSwapUtil.cpp
//
//   Description : swap definition, for use with Libor model mainly
//
//   Author      : Henrik Rasmussen
//
//   Date        : 16 Dec 2005
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/SRMUtil.hpp"
#include "edginc/SRMSwap.hpp"
#include "edginc/SRMConstants.hpp"
#include "edginc/SwapTool.hpp"
#include "edginc/YieldCurve.hpp"
#include "edginc/MaturityPeriod.hpp"

DRLIB_BEGIN_NAMESPACE

SRMSwapClass::SRMSwapClass(void) 
               : m_weights(0),
				 startDate( DateTime() ),
				 maturityDate ( DateTime() ),
				 stubAtEnd(true),
                 period("M"),
				 tenorInMonths(12), // swap is an annual libor
				 accrInMonths(12),
				 numPeriods(1),
			     dcc(NULL), 
                 resetDates(NULL),
                 payDates(NULL) 
{
   return;
}

// assignment
SRMSwapClass& SRMSwapClass::operator=(const SRMSwapClass &otherSwap)
{
	// do nothing
	if (this==&otherSwap) return *this;

	try
	{
		m_weights = otherSwap.m_weights;
                
		startDate = otherSwap.startDate;
        maturityDate = otherSwap.maturityDate;
        stubAtEnd = otherSwap.stubAtEnd;
		period = otherSwap.period;
		tenorInMonths = otherSwap.tenorInMonths;
		accrInMonths = otherSwap.accrInMonths;
		numPeriods = otherSwap.numPeriods;
       
		dcc = DayCountConventionSP ( (otherSwap.dcc).get() );
        resetDates = DateTimeArraySP ( (otherSwap.resetDates).get() );      
        payDates = DateTimeArraySP ( (otherSwap.payDates).get() );   
	}
	catch(...)
	{	
	}

	return *this;
}

// copy 
SRMSwapClass::SRMSwapClass(const SRMSwapClass &otherSwap)
{
	// do nothing
	if (this==&otherSwap) return;

	try
	{
		m_weights = otherSwap.m_weights;
                
		startDate = otherSwap.startDate;
        maturityDate = otherSwap.maturityDate;
        stubAtEnd = otherSwap.stubAtEnd;
		period = otherSwap.period;
		tenorInMonths = otherSwap.tenorInMonths;
		accrInMonths = otherSwap.accrInMonths;
		numPeriods = otherSwap.numPeriods;
       
		dcc = DayCountConventionSP ( (otherSwap.dcc).get() );
        resetDates = DateTimeArraySP ( (otherSwap.resetDates).get() );      
        payDates = DateTimeArraySP ( (otherSwap.payDates).get() );  
	}
	catch(...)
	{	
	}

	return;
}

SRMSwapClass::SRMSwapClass(
					DateTime              _startDate,      // (I) Date instrument begins at
                    bool                  _stubAtEnd,
                    DayCountConvention&   _dcc,
					int                   _accrInMonths,
					int                   _tenorInMonths)
					: period("M")
{
	  startDate     =  _startDate;
	  stubAtEnd     =  _stubAtEnd;
	  dcc           =  DayCountConventionSP ( &_dcc);
	  accrInMonths  =  Maths::min(_accrInMonths,_tenorInMonths);
      tenorInMonths =  _tenorInMonths;

      maturityDate = startDate.rollDateInMonths(tenorInMonths);

	  int extraDays;
	  
	  SwapTool::countDates(
                    startDate,
					maturityDate,
					accrInMonths,    // interval = count periods
					"M",             // e.g. Y, M, W, D
					&numPeriods,     // (O) Answer (Quotient) 
					&extraDays);     // (O) Days left over(remainder)

	  m_weights.resize(numPeriods);
	  
	  resetDates.reset(
		  SwapTool::dateArray(
	          	  startDate,     // start here
				  accrInMonths,  // interval = count periods
				  "M",           // e.g. Y, M, W, D
				  0,             // startIdx, 0=start @ basedate, 1=start @ baseDate + interval
				  stubAtEnd,     // arrayIncrement, Usually +1 or -1
				  numPeriods)    // how many dates
		  );

    payDates.reset( new DateTimeArray(numPeriods) );

	for (int i = 0; i < numPeriods; i++)
	{
		(*payDates)[i] = (i < numPeriods - 1) ? (*resetDates)[i+1] : maturityDate;
	}

	return;
}

SRMSwapClass::~SRMSwapClass(void)
{
	 m_weights.clear();
}

double SRMSwapClass::parSwapRate(const IYieldCurve& zc) const
{
    MaturityPeriod maturityPeriod(accrInMonths, "M");
	double parSwapRate = zc.couponRate(
                startDate,       // (I) Date instrument begins at
                maturityDate,    // (I) Date instrument matures at
                maturityPeriod,  // (I) Time between payments 
                false,           // (I) stub at front
                dcc.get());

	return parSwapRate;
}

void SRMSwapClass::setWeights(const IYieldCurve& zc)
{
	double dfBegin, dfEnd, dfStart, dfMat;

	dfStart = zc.pv(startDate);
	dfMat = zc.pv(maturityDate);

	double oneOverNormFactor = 1./(dfStart - dfMat);

	for (int i = 0; i < numPeriods; i++)
	{
	    dfBegin = zc.pv((*resetDates)[i]);
	    dfEnd   = zc.pv((*payDates)[i]);
        m_weights[i] = (dfBegin - dfEnd) * oneOverNormFactor;
	}

	return;
}

const vector<double>& SRMSwapClass::getWeights(void) const
{
	return m_weights;
}

const DateTime& SRMSwapClass::getStartDate(void) const
{
    return startDate;
}

const DateTime& SRMSwapClass::getMatDate(void) const
{
    return maturityDate;
}

const DateTimeArray& SRMSwapClass::getResetDates(void) const
{
   return (*(resetDates.get()));
}

const DateTimeArray& SRMSwapClass::getPayDates(void) const
{
   return (*(payDates.get()));
}

DRLIB_END_NAMESPACE
