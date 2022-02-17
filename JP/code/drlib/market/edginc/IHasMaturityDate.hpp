//----------------------------------------------------------------------------
//
//   Group       : Credit Hybrids QR
//
//   Filename    : IHasMaturityDate.hpp
//
//   Description : 
//
//   Date        : Dec 2006
//
//----------------------------------------------------------------------------

#ifndef IHASMATURITYDATE_HPP
#define IHASMATURITYDATE_HPP

#include "edginc/DateTime.hpp"

#include "edginc/DECLARE.hpp"


DRLIB_BEGIN_NAMESPACE

class IHasMaturityDate	
{

public:

	virtual DateTime maturityDate() const = 0;

};

DECLARE(IHasMaturityDate);

DRLIB_END_NAMESPACE

#endif
