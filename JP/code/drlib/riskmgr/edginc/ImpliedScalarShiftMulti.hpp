//----------------------------------------------------------------------------
//
//   Group       : QR - Credit Hybrids
//
//   Filename    : ImpliedScalarShiftMulti.hpp
//
//   Description : Like ImpliedScalarShift, but implemented as a "new greek".
//
//   Author      : Linus Thand
//
//   Date        : 09 June 2006
//
//----------------------------------------------------------------------------

#ifndef EDG_IMPLIED_SCALAR_SHIFT_MULTI_H
#define EDG_IMPLIED_SCALAR_SHIFT_MULTI_H

#include "edginc/RiskPropertySensitivity.hpp"
#include "edginc/ImpliedScalarShift.hpp" // Uses static methods

DRLIB_BEGIN_NAMESPACE

/** Similar to ImpliedScalarShift but is implemented as a new style sensitivity
 * which allows it to solve across an InstrumentCollection instead of just one 
 * instrument. This is primarily used for managed trades, where a scenario is
 * requested by a client (e.g. CreditNameNotionalLevel) who also want the change
 * to be PV neutral across a set of tranches (stream group). 
 * This means that some other variable (e.g. CDOParallelStrike,
 * a parallel shift in strikes) is solved for in order to make the PV the same 
 * as before the scenario is applied.
 */

class RISKMGR_DLL ImpliedScalarShiftMulti : 
    public RiskQuantityFactorySensitivity {
        
 public:
    static CClassConstSP const TYPE;
    virtual void validatePop2Object();  
    NamedRiskQuantityArraySP nameRiskQuantities(MultiTweakGroupConstSP world,
                                                RiskMappingConstSP riskMapping) const;
 private:
    ~ImpliedScalarShiftMulti();
    static void load(CClassSP& clazz);
    static IObject* defaultConstructor();
    ImpliedScalarShiftMulti();
    ImpliedScalarShiftMulti(const ImpliedScalarShiftMulti& rhs); 
    ImpliedScalarShiftMulti& operator=(const ImpliedScalarShiftMulti& rhs);
    class ObjectiveFunc;

    /** Is this sensitivity made using a discrete shift (ie a jump) or a
        an approximately continuous one. Returns true **/
    virtual bool discreteShift() const; 

    /** Registered vars **/

    //// Sensitivity whose shift size we seek to solve for, eg VEGA_PARALLEL.
    IScalarPerNameShiftSP sensToShift;                                         

    //// Target value
    double targetValue;

    //// Control that the target value refers to (eg delta). If not provided, price
    CControlSP targetControl;                     
    
    //// needed if control is not a price
    OutputNameSP targetControlOutputName;    

    /** Realistically, the only type of root finder we can use here
       (numerical computation of derivs is not generally recommended
       within 1-dimensional root finders).
       Not optional, but a default root finder class is provided. **/
    RootFinder1D::TwoInitValNoDerivSP rootFinder;     
      
    //// name under which the result is reported back
    OutputNameSP resultOutputName;                  
      
    //// lower bound for value
    double lBound;  

    //// upper (high) bound for value used by default if specified
    double hBound;                           

    //// true if the input to shift has to be always positive
    bool isPositive;                        

    /** Try to shift all supported objects, regardless of name.
        Useful when you solve for a sensitivity that is implemented
        by classes without names, e.g. CDO or WeightedInstrumentCollection.
     */
    bool shiftAllNames;

    const static string NAME;
    
    //// Control for final calculation. E.g for finding greeks at solution point.
    CControlSP ctrlAtSolution;

    //// Identifier for the results at the solution.
    string atSolutionIdentifier;

    //// Un-registered vars 
    CControlSP locCtrl;
};

FORWARD_DECLARE(ImpliedScalarShiftMulti)

DRLIB_END_NAMESPACE
#endif

