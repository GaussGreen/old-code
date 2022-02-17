
#ifndef EDG_VECTOR_SHIFT_H
#define EDG_VECTOR_SHIFT_H

#include "edginc/FORWARD_DECLARE.hpp"
#include "edginc/ScalarShift.hpp"
#include "edginc/VectorTweak.hpp"
#include "edginc/TweakQualifierID.hpp"

DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE(ExpiryWindow)

/** Slightly more specialised ScalarShift where the shifting
    information is qualified by an expiry (eg vega pointwise).
    Note still abstract */
class RISKMGR_DLL VectorShift: public ScalarShift,
                   public virtual ITweakQualifierID,
                   public virtual IVectorTweak{
public:    
    friend class VectorShiftHelper;
    static CClassConstSP const TYPE;

    ~VectorShift();

    /** implements a one sided vector derivative for each instance of the
        market data which is sensitive to this SensControl */
    void calculate(TweakGroup*      tweakGroup,
                   CResults*        results);

    /** Returns the expiry which is currently being tweaked */
    virtual ExpiryConstSP getExpiry()const;

    ExpiryWindowConstSP getExpiryWindow() const;

    /** sets the expiry which is currently being tweaked to a certain value*/
    virtual void setExpiry(ExpiryConstSP newExpiry);

    /** Sets which expiry to tweak within the calculate method.
     * This is VERY different from the setExpiry method above: 
     * in setExpiry we set the point CURRENTLY being tweaked (within 
     * the calculate method), whereas here we indicate what is the 
     * only point we WANT to tweak in the calculate method */ 
    virtual void setExpiryToTweak(ExpiryConstSP newExpiry);

    /** Returns the expiries which are to be tweaked.
        This allows the expiries to be either obtained from
        the instrument or from within the SensControl. The
        default is from the instrument. Typically the object parameter would be
        the TweakGroup */
    virtual ExpiryArrayConstSP getExpiries(const IObject* tweakGroup);

    /** Returns the expiries which are to be tweaked. Only valid once
        the tweaking process has started. ie value is derived during
        the tweaking */
    virtual ExpiryArrayConstSP getExpiries()const;

protected:
    /** Note VectorShift is abstract. Create a vector shift of type clazz,
        which uses outputName (eg VEGA_PARALLEL) to identify results and
        with given shiftSize */
    VectorShift(const CClassConstSP& clazz,
                const string&        sensName,
                const double&        shiftSize);

    /** for reflection */
    VectorShift(const CClassConstSP& clazz,
                const string&        sensName);

    /** Calculate shift sizes for given name and set of expiries. Default 
        implementation returns DoubleArray(expiries.size(), getShiftSize()) */
    virtual DoubleArrayConstSP calculateTweakSizes(
        IObject*           tweakGroup, 
        OutputNameConstSP  name,
        const ExpiryArray& expiries);

    /** As above but uses alternative results set for base price */
    void calculate(TweakGroup*  tweakGroup,
                   CResults*    resultsForTweak,
                   CResults*    resultsForBasePrice);

    virtual void setExpiriesToTweak(ExpiryArrayConstSP expiries);

    bool getExpiriesAreSet() const;

    ExpiryConstSP      expiry;
    ExpiryArrayConstSP cachedExpiries;

private:
    VectorShift(const VectorShift& rhs);
    VectorShift& operator=(const VectorShift& rhs);

    /* "expiriesAreSet" is used to indicate whether the setExpiryToTweak 
     * method has been used to set the expiries to tweak.
     * NOTE: there is no way to set this attribute to false once it has 
     * been set to true - which means that if the tweak is applied to
     * multiple names, it may not work */
    bool expiriesAreSet;
};

typedef smartPtr<VectorShift> VectorShiftSP;
#ifndef QLIB_VECTORSHIFT_CPP
EXTERN_TEMPLATE(class RISKMGR_DLL_SP smartPtr<VectorShift>);
#else
INSTANTIATE_TEMPLATE(class RISKMGR_DLL smartPtr<VectorShift>);
#endif


DRLIB_END_NAMESPACE

#endif
