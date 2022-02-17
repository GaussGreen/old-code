//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//
//----------------------------------------------------------------------------

#ifndef SENSITIVE_STRIKES_HPP
#define SENSITIVE_STRIKES_HPP

#include "edginc/config.hpp"
#include "edginc/smartPtr.hpp"
#include "edginc/Object.hpp"
#include "edginc/Model.hpp"

using namespace std;

DRLIB_BEGIN_NAMESPACE

class RISKMGR_DLL ISensitiveStrikes: virtual public IObject{
public:
    ISensitiveStrikes();  // in VegaMatrix.cpp
    virtual ~ISensitiveStrikes();  // in VegaMatrix.cpp

    /** indicates whether VEGA_MATRIX is sensible for this instrument */
    virtual bool avoidVegaMatrix(const IModel* model) = 0;

    /** returns all strikes on the vol surface to which 
        this instrument is sensitive */
    virtual DoubleArraySP getSensitiveStrikes(OutputNameConstSP outputName,
                                              const IModel*      model) = 0;

    static CClassConstSP const TYPE; // in VegaMatrix.cpp
};

// typedef for smart pointers to ISensitiveStrikes
typedef smartConstPtr<ISensitiveStrikes> ISensitiveStrikesConstSP;
typedef smartPtr<ISensitiveStrikes> ISensitiveStrikesSP;

// type registration in VegaMatrix.cpp 

DRLIB_END_NAMESPACE

#endif // SENSITIVE_STRIKES_HPP
