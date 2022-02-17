//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : ElemStateVariable.hpp
//
//   Description : Interface implemented by all classes that can generate
//                 'elementary' state variables.
//
//   Date        : March 2004
//
//
//----------------------------------------------------------------------------

#ifndef EDR_ELEMSTATEVARIABLEGEN_HPP
#define EDR_ELEMSTATEVARIABLEGEN_HPP

#include "edginc/StateVariableGen.hpp"
#include "edginc/ModelException.hpp"

DRLIB_BEGIN_NAMESPACE

/** Interface implemented by all classes that can generate 'elementary'
    state variables. Elementary state variables have no dependency on any
    other state variables */

class IElemStateVariableGenVisitor; // and implementation will have to include "MCPathConfigSR<GenSV.hpp"
class RISKMGR_DLL IElemStateVariableGen: public virtual IStateVariableGen {
public:
    /* probably want more methods here to help path generator understand
       how independent different SV's are.
       We need to know what dates the simulation will need to generate values
       for ? */


    /** for 'visitor' model of sorting SVGens, 
        many thanks to Vladimir for suggesting it */
    virtual void attachSVGen(IElemStateVariableGenVisitor*) const =0;
    /* NOTE: that every derived class should have this method implemented 
       the same way, e.g. SVGenDiscFactor will have

        void SVGenDiscFactor::attachSVGen(SV* sv)
        {
            sv->processSVGen(this);
        }
     */

};

typedef vector<const IElemStateVariableGen*> IElemStateVariableGenArray;

DRLIB_END_NAMESPACE

#endif
