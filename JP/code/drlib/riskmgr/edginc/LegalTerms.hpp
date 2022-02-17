//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : LegalTerms.hpp
//
//   Description : Transform object into one which matches term sheet
//                 rather than one relevant for pricing/greeks
//                 e.g. use economic barrier instead of risk barrier
//
//   Date        : May 2005
//
//
//
//----------------------------------------------------------------------------

#ifndef LEGAL_TERMS__HPP
#define LEGAL_TERMS__HPP
#include "edginc/ScenarioShift.hpp"

DRLIB_BEGIN_NAMESPACE
class OutputRequestCalculator;

/** Sens Control for Legal Terms */
class RISKMGR_DLL LegalTerms: public CObject,
                  public virtual ITweakID,
                  public virtual IScenarioShift,
                  public virtual IPerturbation{ /* IPerturbation for backward
                                                   compatibility only */
public:
    static CClassConstSP const TYPE;

    // public so the LEGAL_TERMS_FV output request has access
    LegalTerms();

    /** Apply 'legal terms' */
    virtual bool applyScenario(IObjectSP object);

    // Nothing to do before market data is retrieved.
    virtual bool preapplyScenario(IObjectSP object);

    /** IPerturbation implementation - for backwards compatibility only.
        Equivalent here to applyScenario(objectToShift) */
    virtual bool findAndShift(IObjectSP         objectToShift, 
                              OutputNameConstSP name);
  
    /** What an object must implement to support LegalTerms */
    class RISKMGR_DLL Shift{
    public:
        friend class LegalTermsHelper;
        static CClassConstSP const TYPE;
        virtual ~Shift();
        /** Shifts the object using given shift. Return true to make the
            infrastructure keep tweaking the components within the object
            which implements this interface */
        virtual bool sensShift(LegalTerms* shift) = 0;
    };
    typedef Shift IShift;

    /** returns the interface identifying what an object has to do in order
        to be support the tweak that this object represents */
    CClassConstSP shiftInterface() const;

    //// implementation of ITweakID
    virtual void reset();

    //// implementation of ITweakNameResolution - return null
    virtual ITweakNameResolver* nameResolver();

    /**
     * @param obj The object to shift. The object must implement the
     LegalTerms.Shift interface. The return value indicates whether
     or not components of this object need to be tweaked ie true:
     infrastructure should continue to recurse through components
     tweaking them; false: the infrastructure shouldn't touch any
     components within this object */
    virtual bool shift(IObjectSP obj);

    /** Creates an OutputRequestCalculator which can be used for a default
        implementation of the LegalTerms output request */
    static OutputRequestCalculator* createCalculator();
private:
    class Calculator;
    friend class LegalTermsHelper;
    OutputNameArraySP toTweak; // to be removed - legacy reasons
    /** for reflection */
    LegalTerms(const LegalTerms &rhs);
    LegalTerms& operator=(const LegalTerms& rhs);
};

typedef smartConstPtr<LegalTerms> LegalTermsConstSP;
typedef smartPtr<LegalTerms> LegalTermsSP;

DRLIB_END_NAMESPACE

#endif
