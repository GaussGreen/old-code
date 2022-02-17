
#ifndef EDG_MATRIX_SHIFT_H
#define EDG_MATRIX_SHIFT_H
#include "edginc/ScalarShift.hpp"
#include "edginc/Expiry.hpp"
#include "edginc/TweakQualifierID.hpp"

DRLIB_BEGIN_NAMESPACE


/** Slightly more specialised ScalarShift where the shifting
    information is qualified by an expiry (eg vega pointwise).
    Note still abstract */
class RISKMGR_DLL MatrixShift: public ScalarShift,
                   public virtual ITweakQualifierID{
public:    
    friend class MatrixShiftHelper;
    static CClassConstSP const TYPE;

    virtual ~MatrixShift();
    
    /** if strike1 > strike2 && strike1/strike2 <= BETA then strikes
        are merged */
    static const double BETA;

    /** implements a one sided vector derivative for each instance of the
        market data which is sensitive to this SensControl */
    void calculate(TweakGroup*      tweakGroup,
                   CResults*        results);

    // /** Returns the expiry which is currently being tweaked */
    // virtual ExpiryConstSP getExpiry();

    /** Returns the expiries which are to be tweaked.
        This allows the expiries to be either obtained from
        the instrument or from within the SensControl. The
        default is from the instrument. Typically the object parameter would be
        the TweakGroup */
    virtual ExpiryArrayConstSP getExpiries(const IObject* tweakGroup);

    /** Returns the expiries x axis values (aka. strikes) which are 
        to be tweaked. */
    virtual DoubleArraySP getXAxis(CInstrument*      instrument, 
                                   OutputNameConstSP outputName) = 0;
    /** returns a new zero shift control */
    virtual MatrixShift* getZeroShift()  const = 0;

    /** returns the x-axis*/
    DoubleArrayConstSP getXAxisValues()  const;
    /** returns the expiries */
    ExpiryArrayConstSP  getExpiries()    const;
    /** returns the lower index to shift */
    int                 getLowerIdx()    const;
    /** returns the upper index to shift */
    int                 getUpperIdx()    const;
    /** returns the index of the expiry to shift */
    int                 getExpiryIdx()   const;

protected:
    /** Note VectorShift is abstract. Create a vector shift of type clazz,
        which uses outputName (eg VEGA_PARALLEL) to identify results and
        with given shiftSize */
    MatrixShift(const CClassConstSP& clazz,
                const string&        sensName,
                const double&        shiftSize);


    /** for reflection */
    MatrixShift(const CClassConstSP& clazz,
                const string&        sensName);

    /* access functions*/
    ResultsSP getBasePrice(TweakGroup*       tweakGroup,
                           OutputNameConstSP name);

protected:
    int           xLowerIdx;     /* lower idx of range to tweak for strikes
                                     that are too close together */
    int           xUpperIdx;     /* upper idx of range to tweak for strikes
                                     that are too close together */
    int           expiryIdx;     /* lower idx of range to tweak for strikes
                                     that are too close together */

    ExpiryArrayConstSP  expiries;   /* the y-coordinate */

    CDoubleArraySP xAxis;           /* All the xOrds that will be
                                       shifted - we need this
                                       information to allow us to
                                       shift properly in some cases -
                                       internal use only */

    MatrixShift(const MatrixShift& rhs);
    MatrixShift& operator=(const MatrixShift& rhs);

};

typedef smartPtr<MatrixShift>      MatrixShiftSP;
typedef smartConstPtr<MatrixShift> MatrixShiftConstSP;

DRLIB_END_NAMESPACE

#endif
