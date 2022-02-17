//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : AmerSpread.hpp
//
//   Description : American spread vanilla option
//
//   Author      : Ning shen
//
//   Date        : 06 Jan 2003
//
//
//----------------------------------------------------------------------------

#ifndef EDR_AMERSPREAD_HPP
#define EDR_AMERSPREAD_HPP

#include "edginc/Generic1Factor.hpp"
#include "edginc/Tree1f.hpp"
#include "edginc/LastSensDate.hpp"

DRLIB_BEGIN_NAMESPACE

/** Shell structure for spot/hybrid/ratio payoffs */
class PRODUCTS_DLL CAmerSpread: public Generic1Factor,
                   public FDModel::IIntoProduct,
                   public virtual LastSensDate
{
public:
    static CClassConstSP const TYPE;

    /** instrument validation */
    virtual void Validate();

    /** input data validation */
    virtual void validatePop2Object();

    /** create a fd payoff tree - new for all fd/tree state variable interface */
    virtual FDProductSP createProduct(FDModel* model) const;

    /** indicates whether VEGA_MATRIX is sensible for this instrument */
    virtual bool avoidVegaMatrix(const IModel* model);

    /** Returns rolls value date and sets initial spot for Theta,
        return true if sub objects need to be tweaked */
    virtual bool sensShift(Theta* shift);
  
    virtual DateTime getValueDate() const;
 
    virtual DateTime endDate(const Sensitivity* sensControl) const;

    /** price a dead instrument until settlement - exercised, expired, knocked out etc.
    returns true if it is dead (and priced), false if it is not dead */
    virtual bool priceDeadInstrument(CControl* control, CResults* results) const;

private:
    friend class CAmerSpreadHelper;
    friend class AmerSpreadFDProd;

    CAmerSpread();
    CAmerSpread(const CAmerSpread& rhs);
    CAmerSpread& operator=(const CAmerSpread& rhs);
    static void load(CClassSP& clazz);

    void addOutputRequests(Control* control,
                           Results* results,
                           const double& fairValue,
                           const double& indVol) const;

protected:
    CAmerSpread(CClassConstSP clazz);

    bool                    isCall;
    bool                    canExerciseEarly;
    ScheduleSP              exerciseScheduleHi;
    ScheduleSP              exerciseScheduleLo;

    double                  spotAtMaturity;

    DateTime                dateExercised;
    bool                    isExercised;
    int                     noExerciseWindow;
};

typedef smartPtr<CAmerSpread> CAmerSpreadSP;

DRLIB_END_NAMESPACE

#endif
