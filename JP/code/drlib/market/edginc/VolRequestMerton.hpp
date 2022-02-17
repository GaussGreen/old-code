//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : VolRequestMerton.hpp
//
//   Description : Merton vol request 
//
//   Author      : Oliver Brockhaus
//
//   Date        : April 2003
//
//
//----------------------------------------------------------------------------

#ifndef EDG_MERTON_VOLREQUEST_HPP
#define EDG_MERTON_VOLREQUEST_HPP

#include "edginc/DateTime.hpp"
#include "edginc/VolRequest.hpp"

DRLIB_BEGIN_NAMESPACE

/** This class represents the Merton volatility */
class MARKET_DLL VolRequestMerton: public CVolRequest {
public:
    static CClassConstSP const TYPE;
    friend class VolRequestMertonHelper;
    VolRequestMerton();

protected:
private:
    VolRequestMerton(const VolRequestMerton &rhs);
    VolRequestMerton& operator=(const VolRequestMerton& rhs);
};

typedef smartConstPtr<VolRequestMerton> VolRequestMertonConstSP;
typedef smartPtr<VolRequestMerton> VolRequestMertonSP;

// array of vol request
//typedef array<VolRequestMertonSP, VolRequestMerton> VolRequestMertonArray;

// smart pointers for CVolRequestMertonArray
//typedef smartConstPtr<VolRequestMertonArray> VolRequestMertonArrayConstSP;
//typedef smartPtr<VolRequestMertonArray> VolRequestMertonArraySP;

DRLIB_END_NAMESPACE

#endif
