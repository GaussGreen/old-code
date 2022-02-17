//----------------------------------------------------------------------------
//
//   Group       : Quantitative Research
//
//   Filename    : MeanCreditReversionLevel.hpp
//
//   Description : Mean credit reversion level scenario - set CDSVol mean reversion to supplied value
//
//   Author      : Antoine Gregoire
//
//   Date        : September 2005
//
//----------------------------------------------------------------------------

#ifndef MEAN_CREDIT_REVERSION_LEVEL_HPP
#define MEAN_CREDIT_REVERSION_LEVEL_HPP
#include "edginc/ScalarPerturbation.hpp"

DRLIB_BEGIN_NAMESPACE

/** Sens Control for scenario - set CDS Par Spread to supplied value */
class RISKMGR_DLL MeanCreditReversionLevel: public ScalarPerturbation {
public:
    static CClassConstSP const TYPE;

    /** What an object must implement to support MeanCreditReversionLevel */
    class RISKMGR_DLL IShift{
    public:
        friend class MeanCreditReversionLevelHelper;
        static CClassConstSP const TYPE;
        virtual ~IShift();
        /** Returns the name of the stock/asset - used to determine
            whether to shift the object */
        virtual string sensName(MeanCreditReversionLevel* shift) const = 0;
        /** Shifts the object using given shift. Return true to make the
            infrastructure keep tweaking the components within the object
            which implements this interface */
        virtual bool sensShift(MeanCreditReversionLevel* shift) = 0;
    };

    /** constructor with explicit cs level */
    MeanCreditReversionLevel(double spread);

    /** returns the interface identifying what an object has to do in order
        to be support the tweak that this object represents */
    CClassConstSP shiftInterface() const;

    // methods to shift given objects

    /** Returns true if the supplied object matches the supplied name
        for this sensitivity. The object must implement this sensitivity's
        Shift interface */
    virtual bool nameMatches(const OutputName&         name,
                             IObjectConstSP            obj);

    /** Appends the name(s) of the supplied object with respect to
        this sensitivity to the supplied list. The object must
        implement this sensitivity's Shift interface */
    virtual void appendName(OutputNameArray&          namesList,
                            IObjectConstSP            obj);

    /**
     * @param obj The object to shift. The object must implement the
     MeanCreditReversionLevel.Shift interface. The return value indicates whether
     or not components of this object need to be tweaked ie true:
     infrastructure should continue to recurse through components
     tweaking them; false: the infrastructure shouldn't touch any
     components within this object */
    virtual bool shift(IObjectSP obj);
    
private:
    friend class MeanCreditReversionLevelHelper;
    /** for reflection */
    MeanCreditReversionLevel();
    MeanCreditReversionLevel(const MeanCreditReversionLevel &rhs);
    MeanCreditReversionLevel& operator=(const MeanCreditReversionLevel& rhs);
};


typedef smartConstPtr<MeanCreditReversionLevel> MeanCreditReversionLevelConstSP;
typedef smartPtr<MeanCreditReversionLevel> MeanCreditReversionLevelSP;

DRLIB_END_NAMESPACE

#endif
