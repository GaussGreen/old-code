//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : ExtendableNote.hpp
//
//   Description   Step down bond instrument. 
//
//
//----------------------------------------------------------------------------

#ifndef EDG_EXTENDABLE_NOTE_HPP
#define EDG_EXTENDABLE_NOTE_HPP

#include "edginc/DblBarrier.hpp"
#include "edginc/SwapLegIntFace.hpp"

DRLIB_BEGIN_NAMESPACE

class CExtendableNote : public CDblBarrier
{
public:
    static CClassConstSP const TYPE;

    virtual void Validate();
 
    void addOutputRequests(Control* control,
                           Results* results,
                           const double& fairValue,
                           const double& indVol) const;
    
    bool sensShift(Theta* shift);

	/** create a fd payoff tree - new for all fd/tree state variable interface */
    virtual FDProductSP createProduct(FDModel* model) const; 

    /** price a dead instrument until settlement - exercised, expired, knocked out etc.
        returns true if it is dead (and priced), false if it is not dead */
    bool priceDeadInstrument(CControl* control, CResults* results) const;

private:
    friend class CExtendableNoteHelper;
    friend class CExtendableNoteFDProd;

protected:
    static void load(CClassSP& clazz);
    CExtendableNote();
    CExtendableNote(CClassConstSP clazz);

    double          Notional; // $unregistered
    double          SpotFixing;

    CDateTime               CouponLockDate;        //Before this date, all coupon never KO.

    // STB member data
    LiborLegSP liborLeg;

    FixedLegSP  fixedLeg;

    bool    SwitchPlainSwap;         //Price with PlainSwap or Not.  Default = True (with PlainSwap)
    bool    RemoveSettlementOption;  //default = false, If final put/fwd not needed, make it true
    bool    RemoveCoupon;            //default = false. If only need final put/fwd value, make it true

    bool    isEquityPayer;          //Receive equity payoff for false.
    bool    isLiborPayer;           //Recieve Fixed Leg for false.
};

/////////////////////////////////////////////////////////
//           fd product
/////////////////////////////////////////////////////////
/** ExtendableNote product payoff for a FD tree */

class CExtendableNoteFDProd:  public DblBarrierFDProd
{
public:
    CExtendableNoteFDProd(const CExtendableNote* ExN, FDModel* mdl):DblBarrierFDProd(ExN,mdl), inst(ExN)
    {       
        numPrices = 4; 
    }
	
    ~CExtendableNoteFDProd() {}

	virtual void update(int& step, FDProduct::UpdateType type);

protected:
   
    //Input Array contains dammy Date
    CDoubleArray        Coupon;                 //Fixed Coupon will be payed if no KO.
    double              SettleOptionPrice;      //Store the settlement option price
    double              KOSwapPrice;            //Store the conditinal swap price

private:

	/** initialise product specific data */
    virtual void initProd();

    /** product payoff method at maturity */
	void prod_BWD_T(const TreeSlice & spot,
								   int step, 
								   int bot, 
                                   int top, 
							       int pStart, 
							       int pEnd,
                                   const vector< TreeSliceSP > & price);

    /** product payoff method at steps earlier than maturity */
	void prod_BWD(const TreeSlice & spot,
						        int step, 
						        int bot, 
								int top, 
								int pStart, 
								int pEnd,
								const vector< TreeSliceSP > & price);

	/** vanilla, quanto or struck */
    virtual string getCcyTreatment() const
	{
        return inst->ccyTreatment;
    }

	/** ignore start date if not forward starting */
    virtual DateTime getStartDate() const
    {
        return inst->fwdStarting ? inst->startDate : inst->valueDate;
    }

    /** extra output requests */
    virtual bool Positive() { return false; }

    /** make price refinement - combine Libor leg with step down coupons */
    double scalePremium(double P0, 
						double P1,
			 			YieldCurveConstSP disc);

    void    CalcStepCoupon();   //Shift coupon from pay date to accrue start dates for avoid barrier.

    /** extra output requests */
    virtual void recordOutput(Control* control, 
							  YieldCurveConstSP disc, 
							  Results* results);

   const CExtendableNote* inst; 
};

DRLIB_END_NAMESPACE
#endif
