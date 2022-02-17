//----------------------------------------------------------------------------
//
//   Group       : QR Equities
//
//   Filename    : IPDFBoundaryProb.hpp
//
//   Description : Implemented by market data in order to use 
//                 CProcessedVolBS::defaultStrikes.
//
//   Author      : Jon Dee
//
//   Date        : 27 October 2006
//
//
//----------------------------------------------------------------------------

#ifndef EDG_IPDF_BOUNDARY_PROB_HPP
#define EDG_IPDF_BOUNDARY_PROB_HPP

#include "edginc/DECLARE.hpp"

DRLIB_BEGIN_NAMESPACE

/** Implement this interface on classes that call CProcessedVolBS::defaultStrikes() 
  * The purpose of this interface is to give a quick method to look up such classes in the market data cache. */
class MARKET_DLL IPDFBoundaryProb : public virtual IObject
{
public:

    static CClassConstSP const TYPE;
    static void load(CClassSP& clazz);

    IPDFBoundaryProb(void) {}
    virtual ~IPDFBoundaryProb(void) {}

    virtual double getPDFBoundaryProb() const = 0;
};



DRLIB_END_NAMESPACE

#endif
