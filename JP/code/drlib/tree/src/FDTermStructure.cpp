//----------------------------------------------------------------------------
//
//   Group       : Convertibles Research
//
//   Filename    : FDTermStructure.cpp
//
//   Description : A term structure object used in the generic 1F FD Engine
//
//   Author      : André Segger
//
//   Date        : 10 October 2003
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/Expiry.hpp"
#include "edginc/FDTermStructure.hpp"

DRLIB_BEGIN_NAMESPACE

FDTermStructure::FDTermStructure(double flatTerm):CObject(TYPE)
{ 
    term = flatTerm; 
    dimension = 0;
}

FDTermStructure::FDTermStructure(DoubleArraySP termArray):CObject(TYPE)
{ 
    if (!termArray) {
        throw ModelException("FDTermStructure::TermStructure",
            "Supplied array must not be null");
    }

    if (termArray->size() == 0) {
        throw ModelException("FDTermStructure::TermStructure",
            "Must have at least one element in term structure array");
    }

    termStructure = termArray; 
    dimension = 1;
}

double FDTermStructure::operator()(const int r) const 
{
    if ( dimension == 0 ) {
        return term;
    } else {
        if ( r >= termStructure->size() || r < 0) {
            throw ModelException("TermStructure::operator()",
                "Tried to access element outside array bounds");
        }

        return (*termStructure)[r];
    }
}

void FDTermStructure::scale(const double& scalingFactor)
{
    if ( dimension == 0 ) {
        term *= scalingFactor;
    } else {
        for(int i=0 ; i<termStructure->size() ; ++i) {
            (*termStructure)[i] *= scalingFactor;
        }
    }
}


class FDTermStructureHelper {
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(Expiry, clazz);
        SUPERCLASS(CObject);
        // no fields
    }

};

CClassConstSP const FDTermStructure::TYPE = CClass::registerClassLoadMethod(
    "FDTermStructure", typeid(FDTermStructure), FDTermStructureHelper::load);

DEFINE_TEMPLATE_TYPE(FDTermStructureArray);

DRLIB_END_NAMESPACE
