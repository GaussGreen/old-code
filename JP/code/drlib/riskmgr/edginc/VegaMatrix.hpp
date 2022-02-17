//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : VegaMatrix.cpp
//
//   Description : Controls calculation of VEGA_MATRIX
//
//   Author      : Andre Segger
//
//   Date        : 30 March 2001
//
//----------------------------------------------------------------------------

#ifndef EDG_VEGA_MATRIX_HPP
#define EDG_VEGA_MATRIX_HPP

#include "edginc/FORWARD_DECLARE.hpp"
#include "edginc/MatrixShift.hpp"
#include "edginc/Additive.hpp"

DRLIB_BEGIN_NAMESPACE

class VegaProxyMatrix;
FORWARD_DECLARE(VegaPointwise)

/** Sens Control for vega matrix */
class RISKMGR_DLL VegaMatrix: public MatrixShift,
                  public Additive {
public:
    friend class VegaMatrixHelper;
    static CClassConstSP const TYPE;
    const static string NAME;
    const static double DEFAULT_SHIFT;

    virtual ~VegaMatrix();

    /** What an object must implement to be tweakable for VEGA_MATRIX */
    class RISKMGR_DLL IShift{
    public:
        friend class VegaMatrixHelper;
        static CClassConstSP const TYPE;
        IShift();
        virtual ~IShift();
        /** Returns the name of the vol - used to determine whether to tweak
            the object */
        virtual string sensName(VegaMatrix* shift) const = 0;

        /** Return the array of expiries (ie maturities/benchmark dates) that
            need to be tweaked for this vol */
        virtual ExpiryArrayConstSP sensExpiries(VegaMatrix* shift) const =0;

        /** Shifts the object using given shift. Return true to make the
            infrastructure keep tweaking the components within the object
            which implements this interface */
        virtual bool sensShift(VegaMatrix* shift) = 0;
    };

    /** What an object must implement to be able to perform a restorable
        tweak for VEGA_MATRIX. This is used by {@link SensMgr SensMgr} to
        determine whether it should try and tweak a particular object. */
    class RISKMGR_DLL IRestorableShift: public virtual IShift{
    public:
        friend class VegaMatrixHelper;
        static CClassConstSP const TYPE;
        IRestorableShift();
        virtual ~IRestorableShift();
        /** Restores the object to its original form */
        virtual void sensRestore(VegaMatrix* shift) = 0;
    };

    /** constructor with explicit shift size */
    VegaMatrix(double     shiftSize);

    /** constructor with extra chilli sauce */
    VegaMatrix(double     shiftSize,
               IModel*    model,
               Control*   control);

    /** the full monty */
    VegaMatrix(double                     shiftSize,
               int                        expiryIdx,
               const ExpiryArrayConstSP&  allExpiries,
               int                        strikeIdx,
               const DoubleArraySP&       strikes);

    /** returns the name identifying the market data to be shifted. Returns
        null if not set */
    virtual OutputNameConstSP getMarketDataName() const;

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
     VegaMatrix.Shift interface. The return value indicates whether
     or not components of this object need to be tweaked ie true:
     infrastructure should continue to recurse through components
     tweaking them; false: the infrastructure shouldn't touch any
     components within this object */
    virtual bool shift(IObjectSP obj);
    
    /** The supplied object is queried for the expiries array needed
        for doing vega pointwise and this array is returned. The supplied
        object must implement the VegaMatrix.Shift interface */
    virtual IObjectConstSP qualifier(IObjectConstSP obj);

    /**
     * @param param1 The object to shift. The
     object must implement the VegaMatrix.RestorableShift interface
     */
    virtual void restore(IObjectSP obj);

    /** returns all x axis elements to be tweaked */
    virtual DoubleArraySP getXAxis(CInstrument*      instrument, 
                                   OutputNameConstSP outputName);

    /** VEGA_MATRIX algorithm */
    virtual void calculate(TweakGroup*      tweakGroup,
                           CResults*        results);
    /** returns a new zero shift control */
    virtual MatrixShift* getZeroShift() const;

    static VegaMatrix* fromProxy(VegaProxyMatrix* proxy);

private:
    // fields
    VegaPointwiseSP vegaPointwise; // transient. Used when doing vega pointwise
    /** for reflection */
    VegaMatrix();
    VegaMatrix(const VegaMatrix &rhs);
    VegaMatrix& operator=(const VegaMatrix& rhs);
};

typedef smartConstPtr<VegaMatrix> VegaMatrixConstSP;
typedef smartPtr<VegaMatrix> VegaMatrixSP;


DRLIB_END_NAMESPACE

#endif
