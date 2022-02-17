//----------------------------------------------------------------------------
//
//   Group       : QR&D
//
//   Author      : Mark A Robson
//
//   Date        : 25 July 2005
//
//
//----------------------------------------------------------------------------

#ifndef EDR_SCALAR_TWEAK_HPP
#define EDR_SCALAR_TWEAK_HPP


DRLIB_BEGIN_NAMESPACE


/** Interface for tweaks that have a well defined shift size which is a 
    double */
class RISKMGR_DLL IScalarTweak{
public:
    IScalarTweak(); // In SensControl.cpp

    virtual ~IScalarTweak(); // In SensControl.cpp

    /** Returns size of the shift */
    virtual double getShiftSize() const = 0;
};

DRLIB_END_NAMESPACE

#endif

