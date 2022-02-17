//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : HaveEquity.hpp
//
//   Description : Interface class for getEquity() method
//
//   Author      : Jay Blumenstein
//
//   Date        : 30 Oct 2002
//
//
//----------------------------------------------------------------------------

#ifndef HAVE_EQUITY_HPP
#define HAVE_EQUITY_HPP

DRLIB_BEGIN_NAMESPACE

/** abstract interface class  */
class MARKET_DLL IHaveEquity
{
public:

	/** returns a smart pointer to the equity */
	virtual EquitySP getEquity() const = 0;
};

// typedef smartPtr<IHaveEquity> IHaveEquitySP;

DRLIB_END_NAMESPACE

#endif
