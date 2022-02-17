//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : IRVegaPointwise.hpp
//
//   Description : Pointwise IR vega sensitivity
//
//   Author      : Andrew J Swain
//
//   Date        : 25 February 2002
//
//   WARNING  WARNING  WARNING  WARNING  WARNING  WARNING  WARNING  WARNING  
//        This class is now deprecated - Use IRVegaMatrix instead.
//        See IRVegaMatrix.hpp for more information
//   WARNING  WARNING  WARNING  WARNING  WARNING  WARNING  WARNING  WARNING  
//
//----------------------------------------------------------------------------

#ifndef IRVEGAPOINTWISE_HPP
#define IRVEGAPOINTWISE_HPP
#include "edginc/IRGridShift.hpp"
#include "edginc/Additive.hpp"

DRLIB_BEGIN_NAMESPACE

class VegaProxyMatrix;

/** Sens Control for pointwise IR vega */
class RISKMGR_DLL IRVegaPointwise: public IRGridShift,
                       public virtual Additive {
public:
    static CClassConstSP const TYPE;
    const static string NAME;
    const static double DEFAULT_SHIFT;

    /** What an object must implement to be tweakable */
    class RISKMGR_DLL IShift{
    public:
        friend class IRVegaPointwiseHelper;
        static CClassConstSP const TYPE;
        virtual ~IShift();
        /** Returns the name of the vol - used to determine whether to tweak
            the object */
        virtual string sensName(IRVegaPointwise* shift) const = 0;

        /** Return the grid points which are to be tweaked */
        //virtual IRGridPointArrayConstSP sensGridPoints(IRVegaPointwise* shift) const =0;

        /** Shifts the object using given shift. Return true to make the
            infrastructure keep tweaking the components within the object
            which implements this interface */
        virtual bool sensShift(IRVegaPointwise* shift) = 0;
    };

    /** What an object must implement to be able to perform a restorable
        tweak. This is used by {@link SensMgr SensMgr} to
        determine whether it should try and tweak a particular object. */
    class RISKMGR_DLL IRestorableShift: public virtual IShift{
    public:
        friend class IRVegaPointwiseHelper;
        static CClassConstSP const TYPE;
        virtual ~IRestorableShift();
        /** Restores the object to its original form */
        virtual void sensRestore(IRVegaPointwise* shift) = 0;
    };

    /** Class to be implemented by models which use IR swaption vols. We should
        consider just putting this method on CModel directly since it just
        makes life easier if you wrap a model inside another model */
    class RISKMGR_DLL ISensitivePoints {
    public:
        virtual ~ISensitivePoints();
        /** returns all the points on the ir vol surface to which this
            model for the supplied instrument is sensitive.
            A null return value is ok - it
            will cause an exception to be thrown (rather than crash) */
        virtual IRGridPointAbsArraySP getSensitiveIRVolPoints(
            OutputNameConstSP  outputName,
            const CInstrument* inst) const = 0;
    };

    virtual IRGridPointAbsArrayConstSP getGridPoints(
        CInstrument*      instrument, 
        OutputNameConstSP outputName);


    // should instrument bother doing this greek or not?
    virtual bool avoid(CInstrument* instrument,
                       IModel*      model);

    /** constructor with explicit shift size */
    IRVegaPointwise(double shiftSize);

    /** Once used to make a shift, this reports the appropriate divisor
        for this sensitivity */
    double divisor() const;

    /** returns the interface identifying what an object has to do in order
        to be support the tweak that this object represents */
    CClassConstSP shiftInterface() const;

    /** returns the interface identifying what an object has to do in order
        to be support the tweak that this object represents which is also
        restorable */
    CClassConstSP restorableShiftInterface() const;


    // methods to shift given objects

    /** Returns true if the supplied object matches the supplied name
        for this sensitivity. The object must implement this sensitivity's
        Shift interface */
    virtual bool nameMatches(const OutputName&         name,
                             IObjectConstSP          obj);

    /** Appends the name(s) of the supplied object with respect to
        this sensitivity to the supplied list. The object must
        implement this sensitivity's Shift interface */
    virtual void appendName(OutputNameArray&          namesList,
                            IObjectConstSP          obj);

     /**
     * @param obj The object to shift. The object must implement the
     IRVegaPointwise.Shift interface. The return value indicates whether
     or not components of this object need to be tweaked ie true:
     infrastructure should continue to recurse through components
     tweaking them; false: the infrastructure shouldn't touch any
     components within this object */
    virtual bool shift(IObjectSP obj);
    
    /**
     * @param param1 The object to shift. The
     object must implement the IRVegaPointwise.RestorableShift interface
     */
    virtual void restore(IObjectSP obj);

protected:
    IRVegaPointwise(CClassConstSP clazz, const string name);
    IRVegaPointwise(CClassConstSP clazz, const string name, const double shiftSize);

private:
    friend class IRVegaPointwiseHelper;
    /** for reflection */
    IRVegaPointwise();
    IRVegaPointwise(const IRVegaPointwise &rhs);
    IRVegaPointwise& operator=(const IRVegaPointwise& rhs);
};

typedef smartConstPtr<IRVegaPointwise> IRVegaPointwiseConstSP;
typedef smartPtr<IRVegaPointwise> IRVegaPointwiseSP;


DRLIB_END_NAMESPACE

#endif
