//----------------------------------------------------------------------------
//
//   Group       : QR - Credit Hybrids
//
//   Filename    : CreditMultiBetaLevel.cpp
//
//   Description : Beta substitution scenario -
//                 Supply a list of names and their corresponding new betas.
//
//   Author      : Linus Thand
//
//   Date        : 24 May 2006
//
//----------------------------------------------------------------------------

#ifndef CREDITMULTIBETALEVEL__HPP
#define CREDITMULTIBETALEVEL__HPP
#include "edginc/Perturbation.hpp"
#include "edginc/AtomicArray.hpp"

DRLIB_BEGIN_NAMESPACE

/** A scenario which replaces betas for a list of names in a CDO portfolio. This is 
 * used by managed trades where clients can request a change in 
 * a porfolio, and where another variable can be solved for for PV neutrality 
 * (e.g a parallel shift in strikes). **/

class RISKMGR_DLL CreditMultiBetaLevel : public Perturbation {
 public:    
    static CClassConstSP const TYPE;

    /** What an object must implement to support  */
    class IShift {
     public:
        static CClassConstSP const TYPE;
        /** Returns the name of the asset - used to determine
            whether to shift the object */
        virtual string sensName(CreditMultiBetaLevel* shift) const = 0;
        /** Shifts the object using given shift. Return true to make the
            infrastructure keep tweaking the components within the object
            which implements this interface */
        virtual bool sensShift(CreditMultiBetaLevel* shift) = 0;
    };

    CClassConstSP shiftInterface() const;
    virtual bool nameMatches(const OutputName& name, IObjectConstSP obj);
    virtual void appendName(OutputNameArray& namesList, IObjectConstSP obj);
    virtual bool shift(IObjectSP obj);
    virtual void validatePop2Object();
    StringArrayConstSP getNames() const; 
    DoubleArrayConstSP getBetas() const; 

 private:
    ~CreditMultiBetaLevel();
    CreditMultiBetaLevel();
    CreditMultiBetaLevel(const CreditMultiBetaLevel& rhs); 
    CreditMultiBetaLevel& operator=(const CreditMultiBetaLevel& rhs);
    static void load(CClassSP& clazz);
    static IObject* defaultConstructor();
    
     /** Fields **/
    StringArraySP names; //An array of names for which betas will be substituted.
    DoubleArraySP betas; //An array of new beta values.
};

DRLIB_END_NAMESPACE

#endif
