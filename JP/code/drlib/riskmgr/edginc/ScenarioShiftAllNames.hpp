//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : ScenarioShift.hpp
//
//   Description : Defines a scenario shift to market data
//
//   Author      : Andrew J Swain
//
//   Date        : 24 April 2001
//
//
//   $Log$
//   Revision 1.4  2005/01/19 10:17:56  mrobson
//   sensCtrl is now of type IPerturbation. Keep field name unchanged for
//   backward compatibility
//
//   Revision 1.3  2002/05/17 15:40:52  mrobson
//   Removed semi-bogus parameters from smartPtr and smartConstPtr
//
//   Revision 1.2  2001/05/17 17:41:56  mrobson
//   Simplified interface.
//
//   Revision 1.1  2001/04/26 14:42:26  aswain
//   *** empty log message ***
//
//
//
//----------------------------------------------------------------------------


#ifndef SCENARIOSHIFT_HPP
#define SCENARIOSHIFT_HPP
#include "edginc/ScenarioShift.hpp"

DRLIB_BEGIN_NAMESPACE
/** Implementation of IScenarioShift where all objects of the right type
    regardless of name are shifted */
class RISKMGR_DLL ScenarioShiftAllNames: public CObject,
                             public virtual IScenarioShift,
                             public virtual TweakTypeID{
public:
    static CClassConstSP const TYPE;

    virtual ~ScenarioShiftAllNames();

    //// returns null => tweak regardless of 'name'
    virtual ITweakNameResolver* nameResolver();

    /** does nothing */
    virtual void reset();

    /** apply this scenario shift to the supplied object */
    virtual void shift(IObjectSP object) const;
  
private:
    static void load(CClassSP& clazz);
};


DRLIB_END_NAMESPACE
#endif
