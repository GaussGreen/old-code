/**
 * @file NamedRiskQuantity.hpp
 */

#ifndef DRLIB_NamedRiskObjectQuantity_H
#define DRLIB_NamedRiskObjectQuantity_H

#include "edginc/NamedRiskQuantity.hpp"
#include "edginc/RiskQuantity.hpp"

DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE(NamedRiskObjectQuantity)
//----------------------------------------------------------------------------
//
//   Group       : QR - Credit Hybrids
//
//   Filename    : NamedRiskObjectQuantity.hpp
//
//   Description : Like NamedRiskQuantity, but for returning IObjects.
//
//   Author      : Linus Thand
//
//   Date        : 26 July 2006
//
//
//----------------------------------------------------------------------------


class RISKMGR_DLL NamedRiskObjectQuantity: public NamedRiskQuantity {

    static void load(CClassSP& clazz);
    NamedRiskObjectQuantity(const NamedRiskObjectQuantity& rhs);
    NamedRiskObjectQuantity& operator=(const NamedRiskObjectQuantity& rhs);

    IObjectSP value;

public:

    static CClassConstSP const TYPE;

 
    /**
     * Constructor.
     *
     */

    NamedRiskObjectQuantity(IObjectSP value,
                      IResultsIdentifierConstSP resultsName);

    /**
     * Constructor returning smartPtr.
     */

    static NamedRiskObjectQuantitySP SP(IObjectSP value,
                                  IResultsIdentifierConstSP resultsName);

    ~NamedRiskObjectQuantity();


    /**
     * Evaluate the RiskQuantity and store its value in the Results under the
     * IResultsIdentifier.
     */

    virtual void storeResult(const CDoubleArray& vals, const CDoubleArray& dsts,
                             CResultsSP results) const;
};

DRLIB_END_NAMESPACE


#endif
