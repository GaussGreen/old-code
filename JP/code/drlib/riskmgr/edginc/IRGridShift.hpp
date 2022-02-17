//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : IRGridShift.hpp
//
//   Description : Pointwise type IR vega sensitivities
//
//   Author      : Andrew J Swain
//
//   Date        : 22 February 2002
//
//
//----------------------------------------------------------------------------
#ifndef IRGRIDSHIFT_HPP
#define IRGRIDSHIFT_HPP
#include "edginc/ScalarShift.hpp"
#include "edginc/IRGridPointAbs.hpp"

DRLIB_BEGIN_NAMESPACE

/** Pointwise type IR vega sensitivities */
class RISKMGR_DLL IRGridShift: public ScalarShift{
public:    
    static CClassConstSP const TYPE;

    virtual ~IRGridShift();

    /** Overridden to copy over gridPoint */
    IObject* clone() const;

    /** implements a one sided vector derivative for each instance of the
        market data which is sensitive to this SensControl */
    virtual void calculate(TweakGroup*  tweakGroup,
                           CResults*    results);

    /** what point is being tweaked ? */
    ExpiryConstSP tenor() const;
    ExpiryConstSP expiry() const;

    virtual IRGridPointAbsArrayConstSP getGridPoints(
        CInstrument*      instrument, 
        OutputNameConstSP outputName) = 0;

    // should instrument bother doing this greek or not?
    virtual bool avoid(CInstrument* instrument,
                       IModel*      model) = 0;

protected:
    /** Note VectorShift is abstract. Create a vector shift of type clazz,
        which uses outputName (eg VEGA_PARALLEL) to identify results and
        with given shiftSize */
    IRGridShift(const CClassConstSP& clazz,
                const string&        sensName,
                const double         shiftSize);

    /** for reflection */
    IRGridShift(const CClassConstSP& clazz,
                const string&        sensName);

    IRGridPointAbsConstSP gridPoint; // $unregistered

private:
    friend class IRGridShiftHelper;
    IRGridShift(const IRGridShift& rhs);
    IRGridShift& operator=(const IRGridShift& rhs);
};

DRLIB_END_NAMESPACE

#endif
