//----------------------------------------------------------------------------
//
//   Group       : Quantitative Research & Development
//
//   Filename    : BetaSkewPointwiseTweak.hpp
//
//   Description : Tweak of each relevant point inside a skew surface
//
//   Author      : Antoine Gregoire
//
//   Date        : June 2005
//
//   WARNING  WARNING  WARNING  WARNING  WARNING  WARNING  WARNING  WARNING  
//      This class is now deprecated - Use BetaSkewMatrixTweak instead.
//      See BetaSkewMatrixTweak.hpp for more information
//   WARNING  WARNING  WARNING  WARNING  WARNING  WARNING  WARNING  WARNING  
//
//----------------------------------------------------------------------------

#ifndef QLIB_BETA_SKEW_POINTWISE_TWEAK_H
#define QLIB_BETA_SKEW_POINTWISE_TWEAK_H

#include "edginc/ScalarShift.hpp"
#include "edginc/BetaSkewGridPoint.hpp"

DRLIB_BEGIN_NAMESPACE

#define DEFAULT_VALUE_FOR_TWEAKALL       false
#define DEFAULT_VALUE_FOR_NUM_NEIGHBOURS 1

/** Tweak of each relevant point inside a skew surface */
class RISKMGR_DLL BetaSkewPointwiseTweak: public ScalarShift {
public:    
    static CClassConstSP const TYPE;
    const static string NAME;
    const static double DEFAULT_SHIFT;

    /** Called immediately after object constructed */
    void validatePop2Object();

    /** implements a one sided vector derivative for each instance of the
        market data which is sensitive to this SensControl */
    void calculate(TweakGroup*  tweakGroup,
                   CResults*    results);

    /** Returns the grid points which are to be tweaked */
    BetaSkewGridPointArrayConstSP getGridPoints(
        CInstrument* instrument, 
        OutputNameConstSP outputName);
        
    /** Returns the current grid point which is to be tweaked */        
    BetaSkewGridPoint getCurrentPointToTweak() const;        

    /** What an object must implement to be tweakable */
    class RISKMGR_DLL IShift{
    public:
        static CClassConstSP const TYPE;
        virtual ~IShift();
        /** Returns the name of the skew surface - used to determine whether to tweak
            the object */
        virtual string sensName(BetaSkewPointwiseTweak* shift) const = 0;

        /** Shifts the object using given shift. Return true to make the
            infrastructure keep tweaking the components within the object
            which implements this interface */
        virtual bool sensShift(BetaSkewPointwiseTweak* shift) = 0;
    };

    /** What an object must implement to be able to perform a restorable
        tweak. This is used by {@link SensMgr SensMgr} to
        determine whether it should try and tweak a particular object. */
    class RISKMGR_DLL IRestorableShift: public virtual IShift{
    public:
        static CClassConstSP const TYPE;
        virtual ~IRestorableShift();
     
        /** Restores the object to its original form */
        virtual void sensRestore(BetaSkewPointwiseTweak* shift) = 0;
    };

    /** Class to be implemented by models which use skew surfaces (eg: base correlation models) */
    class RISKMGR_DLL ISensitivePoints {
    public:
        static CClassConstSP const TYPE;
        virtual ~ISensitivePoints();
     
        /**
         * Returns all the points on the skew surface to which this
         * model for the supplied instrument is sensitive.
         * A null return value is ok - it will cause an exception to be thrown
         * */
        virtual BetaSkewGridPointArrayConstSP getSensitiveBetaSkewPoints(
            OutputNameConstSP  outputName,
            const CInstrument* inst,
            bool tweakAll,
            const int numberOfNeighbours) const = 0;
    };

    /** Class to be implemented by instruments which use skew surfaces (eg: CDO) */
    class RISKMGR_DLL ISensitivePointsImt {
    public:
        static CClassConstSP const TYPE;
        virtual ~ISensitivePointsImt();
     
        /**
         * Returns all the points on the skew surface to which this
         * instrument for the supplied model is sensitive.
         * A null return value is ok - it will cause an exception to be thrown
         * */
        virtual BetaSkewGridPointArrayConstSP getSensitiveBetaSkewPoints(
            OutputNameConstSP  outputName,
            const IModel* model,
            bool tweakAll) const = 0;
    };

    /** constructor with explicit shift size */
    BetaSkewPointwiseTweak(double shiftSize);

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
                             IObjectConstSP            obj);

    /** Appends the name(s) of the supplied object with respect to
        this sensitivity to the supplied list. The object must
        implement this sensitivity's Shift interface */
    virtual void appendName(OutputNameArray&          namesList,
                            IObjectConstSP            obj);

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

private:
    // for reflection
    BetaSkewPointwiseTweak();
    static IObject* defaultBetaSkewPointwiseTweak();
    static void load(CClassSP& clazz);

    BetaSkewPointwiseTweak(const BetaSkewPointwiseTweak& rhs);
    BetaSkewPointwiseTweak& operator=(const BetaSkewPointwiseTweak& rhs);

    // FIELDS
      
    /**
     * Optional field to determine if we want to do a "clever" tweak (tweak
     * only the relevant points) or if we want to tweak ALL the points (useful
     * as a debug tool to see for example the effect of strike interpolation)
     * */
    bool tweakAll;
    
    /** Optional field to determine the number of "neighbour" strikes to tweak
        at each side of the used strikes, when tweakAll is false */
    int numberOfNeighbours;

protected:
    BetaSkewPointwiseTweak(CClassConstSP clazz, const string name);
    BetaSkewPointwiseTweak(CClassConstSP clazz, 
                           const string name, 
                           const double shiftSize);
    BetaSkewPointwiseTweak(CClassConstSP clazz, 
                           const string name, 
                           const double shiftSize,
                           const bool tweakAll,
                           const int numberOfNeighbours);

    /** The current point to tweak */
    BetaSkewGridPoint currentPointToTweak;
};

DRLIB_END_NAMESPACE

#endif //QLIB_BETA_SKEW_POINTWISE_TWEAK_H
