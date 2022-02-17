//----------------------------------------------------------------------------
//
//   Group       : Convertibles DR
//
//   Filename    : ImpliedYTP.hpp
//
//   Description : Implied Yield to First Put
//
//   Author      : Tycho von Rosenvinge
//
//   Date        : April 16, 2002
//
//
//----------------------------------------------------------------------------

#ifndef EDG_IMPLIED_YTP_H
#define EDG_IMPLIED_YTP_H
#include "edginc/Sensitivity.hpp"

DRLIB_BEGIN_NAMESPACE
class Results;

/** ImpliedScalarShift Sensitivity. */
class RISKMGR_DLL ImpliedYTP: public Sensitivity{
public:
    friend class ImpliedYTPHelper;
    static CClassConstSP const TYPE;

    /** Is this sensitivity made using a discrete shift (ie a jump) or a
        an approximately continuous one. Returns true */
    virtual bool discreteShift() const;

    class RISKMGR_DLL IFaceYTP {
    public:
        virtual double yieldToFirstPut(double price, bool clean) const = 0;

        virtual bool   hasPut() const = 0;
    };

    void validatePop2Object();

    /** identifies the name used for storing associated results in the output */
    virtual const string& getSensOutputName() const;

    /** Computes YTP and record it */
    virtual void calculate(TweakGroup*     tweakGroup,
                           Results*        results);

    ImpliedYTP(double        targetPrice, 
               bool          quotedClean,
               const string& resultOutputName);

private:
    ImpliedYTP();
    ImpliedYTP(const ImpliedYTP &rhs);
    ImpliedYTP& operator=(const ImpliedYTP& rhs);

    /* Registered vars */
    double                             price;
    bool                               clean;
    OutputNameSP                       resultOutputName; // name under which the result is reported back

    const static string NAME;
};

typedef ImpliedYTP                CImpliedYTP;
typedef smartConstPtr<ImpliedYTP> ImpliedYTPConstSP;
typedef smartPtr<ImpliedYTP>      ImpliedYTPSP;
typedef ImpliedYTPConstSP         CImpliedYTPConstSP;
typedef ImpliedYTPSP              CImpliedYTPSP;

DRLIB_END_NAMESPACE
#endif

