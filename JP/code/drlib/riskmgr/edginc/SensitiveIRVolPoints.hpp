//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : SensitiveIRVolPoints.hpp
//
//   Description : Identifies risk sensitive points on an IR vol surface
//
//   Author      : Andrew J Swain
//
//   Date        : 22 February 2002
//
//
//----------------------------------------------------------------------------

#ifndef SENSITIVEIRVOLPOINTS_HPP
#define SENSITIVEIRVOLPOINTS_HPP

#include "edginc/Object.hpp"
#include "edginc/IRGridPointAbs.hpp"
#include "edginc/FORWARD_DECLARE.hpp"

DRLIB_BEGIN_NAMESPACE
class IModel;
FORWARD_DECLARE(OutputName);

/** Should consider retiring this class once all simple IR instruments go
    through ClosedFormIRLN */
class RISKMGR_DLL ISensitiveIRVolPoints: virtual public IObject{
public:
    static CClassConstSP const TYPE;

    /** returns all strikes on the vol surface to which 
        this instrument is sensitive */
    virtual IRGridPointAbsArraySP getSensitiveIRVolPoints(
        OutputNameConstSP outputName,
        const IModel*     model) const = 0;
};

// typedef for smart pointers to ISensitiveIRVolPoints
typedef smartConstPtr<ISensitiveIRVolPoints> ISensitiveIRVolPointsConstSP;
typedef smartPtr<ISensitiveIRVolPoints> ISensitiveIRVolPointsSP;

DRLIB_END_NAMESPACE

#endif 
