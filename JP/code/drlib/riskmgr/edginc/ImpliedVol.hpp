//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : ImpliedVol.hpp
//
//   Description : Implied Volatility Sensitivity
//
//   Author      : André Segger
//
//   Date        : 23 Jan 2002
//
//
//----------------------------------------------------------------------------

#ifndef EDG_IMPLIED_VOL_H
#define EDG_IMPLIED_VOL_H
#include "edginc/Sensitivity.hpp"
#include "edginc/RootFinder.hpp"

DRLIB_BEGIN_NAMESPACE
class Results;

/** ImpliedScalarShift Sensitivity. */
class RISKMGR_DLL ImpliedVol: public Sensitivity{
public:
    friend class ImpliedVolHelper;
    static CClassConstSP const TYPE;

    void validatePop2Object();

    /** Is this sensitivity made using a discrete shift (ie a jump) or a
        an approximately continuous one. Returns true */
    virtual bool discreteShift() const;

    /** identifies the name used for storing associated results in the
        output */
    virtual const string& getSensOutputName() const;

    /** Computes implied scalar shift and record it */
    virtual void calculate(TweakGroup*     tweakGroup,
                           Results*        results);

    ImpliedVol(double        target, 
               double        guess, 
               double        tolerance, 
               const string& resultOutputName);

private:
    ImpliedVol();
    ImpliedVol(const ImpliedVol &rhs);
    ImpliedVol& operator=(const ImpliedVol& rhs);

    /* Registered vars */
    double                             targetValue;
    double                             initialGuess;
    RootFinder1D::TwoInitValNoDerivSP  rootFinder;       // Realistically, the only type of root finder we can use here
                                                         // (numerical computation of derivs is not generally recommended
                                                         // within 1-dimensional root finders).
                                                         // Not optional, but a default root finder class is provided.
    OutputNameSP                       resultOutputName; // name under which the result is reported back

    const static string NAME;
};

typedef ImpliedVol                CImpliedVol;
typedef smartConstPtr<ImpliedVol> ImpliedVolConstSP;
typedef smartPtr<ImpliedVol>      ImpliedVolSP;
typedef ImpliedVolConstSP         CImpliedVolConstSP;
typedef ImpliedVolSP              CImpliedVolSP;

DRLIB_END_NAMESPACE
#endif

