//----------------------------------------------------------------------------
//
//   Group       : Convertibles DR
//
//   Filename    : ImpliedYTM.hpp
//
//   Description : Implied Yield to Maturity
//
//   Author      : Tycho von Rosenvinge
//
//   Date        : April 16, 2002
//
//
//----------------------------------------------------------------------------

#ifndef EDG_IMPLIED_YTM_H
#define EDG_IMPLIED_YTM_H
#include "edginc/Sensitivity.hpp"

DRLIB_BEGIN_NAMESPACE
class Results;
/** ImpliedScalarShift Sensitivity. */
class RISKMGR_DLL ImpliedYTM: public Sensitivity{
public:
    friend class ImpliedYTMHelper;
    static CClassConstSP const TYPE;

    /** Is this sensitivity made using a discrete shift (ie a jump) or a
        an approximately continuous one. Returns true */
    virtual bool discreteShift() const;

    class RISKMGR_DLL IFaceYTM {
    public:
        virtual double yieldToMaturity(double price, bool clean) const = 0;
    };

    void validatePop2Object();

    /** identifies the name used for storing associated results in the output */
    virtual const string& getSensOutputName() const;

    /** Computes YTM and record it */
    virtual void calculate(TweakGroup*     tweakGroup,
                           Results*        results);

    ImpliedYTM(double        targetPrice, 
               bool          quotedClean,
               const string& resultOutputName);

private:
    ImpliedYTM();
    ImpliedYTM(const ImpliedYTM &rhs);
    ImpliedYTM& operator=(const ImpliedYTM& rhs);

    /* Registered vars */
    double                             price;
    bool                               clean;
    OutputNameSP                       resultOutputName; // name under which the result is reported back

    const static string NAME;
};

typedef ImpliedYTM                CImpliedYTM;
typedef smartConstPtr<ImpliedYTM> ImpliedYTMConstSP;
typedef smartPtr<ImpliedYTM>      ImpliedYTMSP;
typedef ImpliedYTMConstSP         CImpliedYTMConstSP;
typedef ImpliedYTMSP              CImpliedYTMSP;

DRLIB_END_NAMESPACE
#endif

