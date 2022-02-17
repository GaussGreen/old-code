#ifndef _VOLTOOLS_HPP
#define _VOLTOOLS_HPP

#include "edginc/Object.hpp"
#include "edginc/DoubleMatrix.hpp"
#include "edginc/IRVol.hpp"
#include "esl_types.h"


DRLIB_BEGIN_NAMESPACE

class VolTools : public CObject{
public:
    static CClassConstSP const TYPE;

protected:

    IRVolWrapper         irvol;
    ExpiryArraySP		 selectedTenors;
    ExpiryArraySP		 selectedExpiries;                

public:
    VolTools(ExpiryArraySP Exp, ExpiryArraySP Tenor, IRVolSP ir, CClassConstSP const &type=TYPE);

    CDoubleMatrixSP  getSwapVol();                  // calculate model swap vols
    virtual void validatePop2Object();

protected:

    VolTools(CClassConstSP const &type=TYPE) : CObject(type) {}

private:
    static void load(CClassSP& clazz);
    static IObject* defaultConstructor(void) { return new VolTools(); }
};



DRLIB_END_NAMESPACE

#endif


