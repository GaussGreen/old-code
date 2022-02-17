//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : MonteCarloLib.cpp
//
//   Description : Ensures all classes get linked in (important if there are
//                 no other dependencies upon a class)
//
//   Date        : June 2001
//
//
//----------------------------------------------------------------------------

#ifndef EDG_MONTECARLOLIB_HPP
#define EDG_MONTECARLOLIB_HPP

DRLIB_BEGIN_NAMESPACE

/** class insures all classes in the mcarlo directory are linked in. */
class MCARLO_DLL MonteCarloLib{
public:
    static void linkInClasses();
private:
    // cannot be instantiated
    MonteCarloLib();
};
        
DRLIB_END_NAMESPACE

#endif
