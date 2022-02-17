//----------------------------------------------------------------------------
//
//   Group       : Credit QR
//
//   Filename    : KRiskyComponent.hpp
//
//   Description : Extension of KComponent to add a credit curve for
//                 components with credit risk.
//
//   Author      : Charles Morcom
//
//   Date        : May 15, 2006
//
//----------------------------------------------------------------------------
#ifndef QR_KRISKYCOMPONENT_HPP
#define QR_KRISKYCOMPONENT_HPP

#include "edginc/FORWARD_DECLARE.hpp"
#include "edginc/FORWARD_DECLARE_WRAPPER.hpp"
#include "edginc/KComponent.hpp"

DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE(KRiskyComponent)
FORWARD_DECLARE_WRAPPER(ICDSParSpreads)

/**KRiskyComponent is just a KComponent which has a credit curve ascribed to it.*/
class KRiskyComponent: public KComponent
{
public:
    static CClassConstSP const TYPE;
    ICDSParSpreadsWrapper           cdsParSpreads;   

    virtual const string creditCurveName() const {return cdsParSpreads.getName();};

protected:
    KRiskyComponent(CClassConstSP const &type) : KComponent(type) {}
    KRiskyComponent(const string &discount, const string& cdsParSpreads,
        const string &outputName, CClassConstSP const &type) 
        : KComponent(discount, outputName, type), cdsParSpreads(cdsParSpreads) {}

    static void load(CClassSP& clazz);
    static IObject* defaultConstructor(void) { return new KRiskyComponent(TYPE); }
};

DRLIB_END_NAMESPACE

#endif
