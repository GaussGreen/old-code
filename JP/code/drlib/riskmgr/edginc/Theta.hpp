//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : Theta.hpp
//
//   Description : Theta shift
//
//   Author      : Andrew J Swain
//
//   Date        : 16 May 2001
//
//
//----------------------------------------------------------------------------


#ifndef THETA_HPP
#define THETA_HPP
#include "edginc/FORWARD_DECLARE_WRAPPER.hpp"
#include "edginc/SensControlAllNames.hpp"
#include "edginc/Additive.hpp"
// MSVC 6 can't cope if we forward declare Holiday 
#include "edginc/Holiday.hpp"

DRLIB_BEGIN_NAMESPACE
FORWARD_DECLARE(SensControlPerName);
FORWARD_DECLARE(Untweakable);
class Theta;

typedef smartConstPtr<Theta> ThetaConstSP;
typedef smartPtr<Theta> ThetaSP;

/** Theta shift */
class RISKMGR_DLL Theta: public SensControlAllNames,
             public virtual Additive {
public:
    friend class RollingTheta;
    static CClassConstSP const TYPE;
    const static string NAME;
    const static int    DEFAULT_SHIFT;

    /** Interface for sensitivities that calculate a greek some point
        in the future. This allows a central routine to calculate them
        together at the same time as theta. This class should ideally
        extend the Sensitivity interface - except that the Sensitivity
        class RISKMGR_DLL is not an interface. */
    class RISKMGR_DLL INextDaySensitivity: public virtual IObject {
    public:
        static CClassConstSP const TYPE;
        virtual ~INextDaySensitivity();
    
        /** Returns the sensitivity to be used on the date in the
            future in order to calculate the required greek. The names()
            method must return what to tweak (if anything). If hasNames()
            returns false, writeResults will still be invoked. The
            Sensitivity::removeNames method will be invoked on the return
            object to remove names that have already been calculated */
        virtual SensitivitySP requiredSensitivity(
            TweakGroup* tweakGroup) const = 0;

        /** Returns the holiday object together with an offset
            indicating number of business days to roll */
        virtual ThetaSP offsetRequired(const CInstrument* inst) const = 0;

        /** Write the results from src to dest. The untweakable object
            will be non null if the initial pricing on the rolled date
            failed (in which case src will be null). The reqdSens is that
            returned from requiredSensitivity. src will be null if there
            was nothing to price */
        virtual void writeResults(const SensitivitySP&          reqdSens,
                                  Results*                      dest,
                                  const UntweakableSP&          untweakable,
                                  const Results*                src) const = 0;
        class RISKMGR_DLL Util{
        public:
            /** Simple implementation of above method where you can specify
                a SensControl */
            static SensitivitySP requiredSensitivity(
                const Sensitivity*   nextDaySens,
                SensControlPerNameSP sensCtrl, /* eg Delta when doing
                                                  DeltaNextDay */
                TweakGroup*          tweakGroup);
            /** Simple implementation of above method where sensCtrl
                is the 'this' passed to writeResults. The
                copyWholePacket determines whether the names within
                reqdSens is ignored and the whatever results are in
                the packet are just copied over. Otherwise only the
                names listed in reqdSens are copied over */
            static void writeResults(
                const Sensitivity*            nextDaySens,
                const SensitivitySP&          reqdSens,
                bool                          copyWholePacket,
                Results*                      dest,
                const UntweakableSP&          untweakable,
                const Results*                src);
        };
    };
    typedef smartPtr<INextDaySensitivity> INextDaySensitivitySP;

    /** Is this sensitivity made using a discrete shift (ie a jump) or a
        an approximately continuous one. Returns true */
    virtual bool discreteShift() const;

    /** What an object must implement to be tweakable for THETA */
    class RISKMGR_DLL Shift{
    public:
        friend class ThetaHelper;
        static CClassConstSP const TYPE;
        Shift();
        virtual ~Shift();
        /** Shifts the object using given shift. Return true to make the
            infrastructure keep tweaking the components within the object
            which implements this interface */
        virtual bool sensShift(Theta* shift) = 0;
    };
    typedef Shift IShift;

    /** What an object must implement to be able to perform a restorable
        tweak for THETA. This is used by {@link SensMgr SensMgr} to
        determine whether it should try and tweak a particular object. */
    class RISKMGR_DLL RestorableShift: public virtual Shift{
    public:
        friend class ThetaHelper;
        static CClassConstSP const TYPE;
        virtual ~RestorableShift();
        /** Restores the object to its original form */
        virtual void sensRestore(Theta* shift) = 0;
    };
    typedef RestorableShift IRestorableShift;

    /** constructor */
    Theta(int offset, HolidaySP hols);

    /** destructor */
    virtual ~Theta();

    /** identifies the packet in which the results are stored. Theta
        results are stored in the instrument packet */
    virtual const string& getPacketName() const;

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

    /**
     * @param obj The object to shift. The object must implement the
     Theta.Shift interface. The return value indicates whether
     or not components of this object need to be tweaked ie true:
     infrastructure should continue to recurse through components
     tweaking them; false: the infrastructure shouldn't touch any
     components within this object */
    virtual bool shift(IObjectSP obj);
    
    /**
     * @param param1 The object to shift. The
     object must implement the Theta.RestorableShift interface
     */
    virtual void restore(IObjectSP obj);

    /** calculates given sensitivity - invoked by calculateSens */
    virtual void calculate(TweakGroup*  tweakGroup,
                           CResults*    results);

    /** Potentially calculates 'nexy day' sensitivities corresponding to the
        specified thetaType (eg theta, delta next day, delta proxy next day etc)
        depending on what's been asked for in the control. Optimises things by
        doing everything together. Does not repeat calculation if result has
        already been calculated. */
    static void calculateNextDaySens(const CClassConstSP& thetaType,
                                     Control*             control,
                                     TweakGroup*          tweakGroup,
                                     CResults*            results);

    /** returns date + offset business days */
    virtual DateTime rollDate(const DateTime &date) const;

    /** Returns true if assets forward value should be used (ie theta
        forward spot) */
    virtual bool useAssetFwds() const;

    /** populate from market cache */
    virtual void getMarket(const IModel* model, const MarketData* market);

    /** scale the results in the Results Object for this sensitivity
        by supplied factor. Default implementation is (provided instance
        of Sensitivity implements the Additive interface) to scale all
        the results in the packet provided by getPacketName() if this is
        the same as getSensOutputName() otherwise the result with name
        getSensOutputName() in packet getPacketName() is scaled */
    virtual void scaleResult(Results*     results,     // (M)
                             double       scaleFactor) const;

    /** Modify results in the Results Object for this sensitivity by
        adding all results in resultsToAdd as indicated by control */
    virtual void addResult(Results*           results,     // (M)
                           const Results*     resultsToAdd,
                           double             scaleFactor) const;


    /** Returns true if theta is the same as this */
    /* todo: virtual */ bool equals(const Theta* theta) const;

    /** tiny little utility class that tells you what the original and
        new value dates are when doing a theta shift. (Originally wanted
        to put these methods on Theta but it became too dangerous trying
        to ensure that theta had the right value date in it when it was
        constructed) */
    class RISKMGR_DLL Util{
    public:
        /** Returns the value date (as reported by the instrument) before the
            theta shift. This method should only be called by shift(Theta*)
            methods */
        const DateTime& getOriginalValueDate() const;

        /** Returns the shifted value date. Equivalent to
            rollDate(getOriginalValueDate()). This method should only be
            called by shift(Theta*) methods */
        const DateTime& getNewValueDate() const;

        /** Returns true if assets forward value should be used (ie theta
            forward spot) */
        bool useAssetFwds() const;

        //// returns original theta object
        const Theta* getTheta() const;

        //// constructor
        Util(const Theta* theta, const DateTime& origValueDate);
    private:
        const Theta*  theta;
        DateTime      origValueDate;
        DateTime      newValueDate;
    };
    /** Returns a ValueDates object corresponding to this shift and
        original value date */
    Util getUtil(const DateTime& origValueDate);

protected:
    /** constructor */
    Theta(int offset, HolidaySP hols,
          const CClassConstSP& clazz,
          const string&        outputName);

    /** empty constructor */
    Theta(const CClassConstSP& clazz,
          const string&        outputName);

    // fields
    int            offset;
    HolidayWrapper hols;
    OutputRequestArraySP outRequests;

    /** for reflection */
    Theta();
    Theta(const Theta &rhs);
    Theta& operator=(const Theta& rhs);
private:
    friend class ThetaHelper;
    void writeResults(Results*             dest, 
                      const UntweakableSP& untweakable,
                      Results*             src);
};

// SP typedef's before class defn

DRLIB_END_NAMESPACE

#endif
