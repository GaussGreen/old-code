//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : ImpliedScalarShift.hpp
//
//   Description : Implied Scalar Shift Sensitivity
//
//   Author      : Regis Guichard
//
//   Date        : 10 Dec 2001
//
//
//----------------------------------------------------------------------------

#ifndef EDG_IMPLIED_SCALAR_SHIFT_H
#define EDG_IMPLIED_SCALAR_SHIFT_H
#include "edginc/ScalarPerNameShift.hpp"
#include "edginc/RootFinder.hpp"
#include "edginc/Control.hpp"

DRLIB_BEGIN_NAMESPACE
class OutputRequest;
typedef smartPtr<OutputRequest> OutputRequestSP;

/** ImpliedScalarShift Sensitivity. */
class RISKMGR_DLL ImpliedScalarShift: public Sensitivity{
public:
    friend class ImpliedScalarShiftHelper;
    static CClassConstSP const TYPE;
      
    /** Is this sensitivity made using a discrete shift (ie a jump) or a
        an approximately continuous one. Returns true */
    virtual bool discreteShift() const;
      
    /** Default root finder. Accuracy needs to be provided by user. */
    class RISKMGR_DLL DefaultRootFinder: public ZBrent{
    public:
        friend class ImpliedScalarShiftDefaultRootFinderHelper;
        static CClassConstSP const TYPE;
        DefaultRootFinder(double shiftAccuracy);
    private:
        DefaultRootFinder();
        DefaultRootFinder(const DefaultRootFinder &rhs);
        DefaultRootFinder& operator=(const DefaultRootFinder& rhs);
    };
      
    void validatePop2Object();
      
    /** identifies the name used for storing associated results in the output */
    virtual const string& getSensOutputName() const;
      
    /** Combines results in packet (ie does a merge) */
    virtual void addResult(Results*           results,     // (M)
                           const Results*     resultsToAdd,
                           double             scaleFactor) const;
      
      
    /** Computes implied scalar shift and record it */
    virtual void calculate(TweakGroup*     tweakGroup,
                           Results*        results);
      
      
    ImpliedScalarShift(IScalarPerNameShiftSP              sensToShift,
                       double                             targetValue,
                       RootFinder1D::TwoInitValNoDerivSP  rootFinder,
                       OutputNameSP                       resultOutputName);
      
    ImpliedScalarShift(IScalarPerNameShiftSP              sensToShift,
                       double                             targetValue,
                       RootFinder1D::TwoInitValNoDerivSP  rootFinder,
                       OutputNameSP                       resultOutputName,
                       bool                               isPositive,
                       OutputRequestSP                    initialGuessRequest);

private:
    class ObjectiveFunc;
    ImpliedScalarShift();
    ImpliedScalarShift(const ImpliedScalarShift &rhs);
    ImpliedScalarShift& operator=(const ImpliedScalarShift& rhs);
    
    /* Registered vars */
    IScalarPerNameShiftSP  sensToShift;                       // sensitivity whose shift size we seek to solve for, eg VEGA_PARALLEL.
    double          targetValue;
    CControlSP      targetControl;                     // control that the target value refers to (eg delta). If not provided, price
    OutputNameSP    targetControlOutputName;           // needed if control is not a price
    RootFinder1D::TwoInitValNoDerivSP  rootFinder;     /* Realistically, the only type of root finder we can use here
                                                          (numerical computation of derivs is not generally recommended
                                                          within 1-dimensional root finders).
                                                          Not optional, but a default root finder class is provided. */
      
    OutputNameSP    resultOutputName;                  // name under which the result is reported back
      
    double          lBound;                            // lower bound for value
    double          hBound;                            // upper (high) bound for value used by default if specified

    bool            isPositive;                        // true if the input to shift has to be always positive
      
    OutputRequestSP initialGuessRequest;               // request to get the initial guess (eg indicative volatility) 
       
    const static string NAME;
      
    /* Un-registered vars */
    CControlSP      locCtrl;
 

public: //Shared with ImpliedScalarShiftMulti
    /** calculate a range containing the root to find 
     *  Returns the pair <lower shift of the range, upper shift of the range> */
    static pair<double,double> calculateRange(const Func1D::NoDeriv & objFunc,
                                              const double &          initialGuess,
                                              const double &          lBound,
                                              const double &          hBound,
                                              const bool &            isPositive); 
};

typedef ImpliedScalarShift CImpliedScalarShift;
typedef smartConstPtr<ImpliedScalarShift> ImpliedScalarShiftConstSP;
typedef smartPtr<ImpliedScalarShift> ImpliedScalarShiftSP;
typedef ImpliedScalarShiftConstSP CImpliedScalarShiftConstSP;
typedef ImpliedScalarShiftSP CImpliedScalarShiftSP;

DRLIB_END_NAMESPACE
#endif

