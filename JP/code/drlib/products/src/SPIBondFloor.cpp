//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : SPIBondFloor.cpp
//
//   Description : SPI bond floor handling
//
//   Date        : Dec 2005
//
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#include "edginc/SPIBondFloor.hpp"
#include "edginc/SPIRunTime.hpp"

DRLIB_BEGIN_NAMESPACE

ISPIBondFloor::~ISPIBondFloor() {};

class FixedSPIBondFloor : virtual public ISPIBondFloor {
    // Independent of fees, lock-in : everything. 
public:
    FixedSPIBondFloor(SPIRunTime*      dynBask,
                      const Schedule*  bondFloorSched):
        bondFloorLevels(dynBask->numSteps, 0.0) {
        
        try {
            const DateTimeArray& stepDates = dynBask->dynBask->getRebalanceDates();
            for(int i=dynBask->iStepFirstRebal; i<dynBask->numSteps;i++) {
                bondFloorLevels[i] = bondFloorSched->interpolate(stepDates[i]);
            }
        } catch (exception& e){
            throw ModelException(e, "FixedSPIBondFloor - failed to compute bond floor levels");
        }            
    }

    const double getLevel(double BL,
                            int    iStep) const {
        return bondFloorLevels[iStep];
    };

    const double getLevelToday() const {
        throw ModelException("SPIBondFloor::getLevelToday",
                "Should not report fixing level for a linear bond floor");
    };

    void refresh(int iFirstFutureStep) {}

private:
    DoubleArray bondFloorLevels;
};

class StandardSPIBondFloor : virtual public ISPIBondFloor {
private:
    class IBondFloorFeeHelper {
    public:
        virtual double adjustNotionalForFees(double notional, int iStep) const = 0;
        virtual double adjustCouponForFees(double coupon, int iStep) const = 0;
        virtual void refresh(int                iFirstFutureStep,
                             const DoubleArray& DFQuot) = 0;
        virtual ~IBondFloorFeeHelper() {};
    };
    class BondFloorFeeHelper {
    private:
        class FeeOnBasketBFHelper : virtual public IBondFloorFeeHelper {
        public:
            FeeOnBasketBFHelper(SPIRunTime* dynBask,
                                double      Basis) :
            maturityAdjustment(dynBask->numSteps, 1.0) {
                const DateTime& finalRebalDate = dynBask->dynBask->getLastRebalDate();
                for (int i = dynBask->iStepFirstRebal; i < dynBask->iStepLastRebal; i++) {
                    maturityAdjustment[i] = 
                        exp((finalRebalDate.getDate()-dynBask->dynBask->getRebalanceDates()[i].getDate())/Basis * 
                            dynBask->algo->feeAtMin(dynBask->fees));
                }
            }

            double adjustNotionalForFees(double notional, int iStep) const {
                return notional * maturityAdjustment[iStep];
            };

            double adjustCouponForFees(double coupon, int iStep) const {
                return coupon * maturityAdjustment[iStep] / maturityAdjustment[iStep+1];
            };
            void refresh(int iFirstFutureStep,
                         const DoubleArray& DFQuot) {
                ; // nothing to do
            }
        private:
            DoubleArray maturityAdjustment;    
        };

        class FeeOnNotionalBFHelper : virtual public IBondFloorFeeHelper {
        // If we have coupons with a non-zero floor, or fees on
        // notional then the bond floor has an additive accounting
        // instead of a purely proportional one.
        // For the paid element I'm including a daily amount in the 
        // bond floor, rather than accumulate and then accounting for
        // a full amount actually being paid. Easier, and ok - I think.
        public:
            FeeOnNotionalBFHelper(SPIRunTime*               dynBask,
                                  const SPIFeesPerNotional* feesPN, 
                                  const DoubleArray&        DFQuot):                
            futureFees(dynBask->numSteps, 0.0),
            dynBask(dynBask), feesPN(feesPN) {
                // note - no fees after last rebal date
                for(int iStep=dynBask->iStepLastRebal-1; iStep>=dynBask->iStepFirstRebal;iStep--) {
                    futureFees[iStep] = (feesPN->getFeeAmount(iStep+1) + futureFees[iStep+1]) 
                                                * DFQuot[iStep];
                }
            }


            double adjustNotionalForFees(double notional, int iStep) const {
                return notional + futureFees[iStep];
            };

            double adjustCouponForFees(double coupon, int iStep) const {
                // do nothing - coupons do not need to be 'adjusted' for fees
                return coupon;
            };

            void refresh(int                iFirstFutureStep,
                         const DoubleArray& DFQuot) {
                // update futureFees
                // I'm a little concerned here that since futureFees[iStep] depends
                // on futureFees[iStep+1] that I need to recalc all the past too. This
                // means each iteration the fees for the PAST change!
                //for(int iStep=dynBask->iStepLastRebal-1; iStep>=iFirstFutureStep;iStep--) {
                for(int iStep=dynBask->iStepLastRebal-1; iStep>=dynBask->iStepFirstRebal;iStep--) {
                    futureFees[iStep] = (feesPN->getFeeAmount(iStep+1) + futureFees[iStep+1]) 
                                                * DFQuot[iStep];
                }
            }

        private:
            DoubleArray futureFees;
            SPIRunTime*               dynBask;
            const SPIFeesPerNotional* feesPN;
        };

        public:
        static IBondFloorFeeHelper* make(SPIRunTime* dynBaskRT, double Basis, 
                                         const DoubleArray& DFQuot) {
            // FeesPerNotional need to be accounted differently in the bond floor
            SPIFeesPerNotional* feesPN = 
                dynamic_cast<SPIFeesPerNotional*>(dynBaskRT->dynBask->feesSPI->getFeesSPI().get());
            if (feesPN) {
                return new BondFloorFeeHelper::FeeOnNotionalBFHelper(dynBaskRT, feesPN, DFQuot);
            }
            else { // fees on basket
                return new BondFloorFeeHelper::FeeOnBasketBFHelper(dynBaskRT, Basis);
            }
        };
    };
    typedef refCountPtr<IBondFloorFeeHelper> IBondFloorFeeHelperSP;

    // calculates all the fixed parts of the bond floor
    DoubleArray preCalc(SPIRunTime*           dynBask, 
                        const ICouponsSPI*    coupons, 
                        IBondFloorFeeHelperSP bondFloorFeeHelper, 
                        const DoubleArray&    DFQuot) {
        DoubleArray fixPart(dynBask->numSteps, 0.0);

        if (dynBask->dynBask->hasBFHistory) {
            int iStep;

            //no adjustment of coupons for fees beyond last rebalance date
            // need to start right at the end as we may get coupons after the last rebalance date
            for(iStep=dynBask->numSteps-2; iStep>=dynBask->iStepLastRebal;iStep--) {
                fixPart[iStep] = (fixPart[iStep+1] + coupons->getGuaranteedCoupon(iStep+1))
                                        * DFQuot[iStep];
            }
            // now adjust for fees back from last rebalance date       
            for(; iStep>=dynBask->iStepFirstRebal;iStep--) {
                fixPart[iStep] = bondFloorFeeHelper->adjustCouponForFees(fixPart[iStep+1] + 
                                                    coupons->getGuaranteedCoupon(iStep+1), iStep)
                                        * DFQuot[iStep];
            }
        }
        return fixPart;
    }

    // update for case where we have stochastic rates
    // this means updating "fixedPart" and the "bondFloorFeeHelper"
    // only deal with future steps
    void refresh(int iFirstFutureStep) {
        dynBask->refreshDiscountFactorsQuotients(DFQuot);
        bondFloorFeeHelper->refresh(iFirstFutureStep, DFQuot);
        if (dynBask->dynBask->hasBFHistory) {
            int iStep;
            const ICouponsSPI* coupons = dynBask->dynBask->couponsSPI->getCouponsSPI().get();

            //no adjustment of coupons for fees beyond last rebalance date
            // need to start right at the end as we may get coupons after the last rebalance date
            for(iStep=dynBask->numSteps-2; iStep>=dynBask->iStepLastRebal;iStep--) {
                fixedPart[iStep] = (fixedPart[iStep+1] + coupons->getGuaranteedCoupon(iStep+1))
                                        * DFQuot[iStep];
            }
            // now adjust for fees back from last rebalance date       
            for(; iStep>=dynBask->iStepFirstRebal;iStep--) {
                fixedPart[iStep] = bondFloorFeeHelper->adjustCouponForFees(fixedPart[iStep+1] + 
                                                    coupons->getGuaranteedCoupon(iStep+1), iStep)
                                        * DFQuot[iStep];
            }
        }
    }

public:
    StandardSPIBondFloor(SPIRunTime*    dynBask, 
                         double         Basis, 
                         const DateTime today) :
    dynBask(dynBask), 
    lastHistoricDateIndex(-1),
    fixingTodayIndex(-1),
    firstCalculatedDateIndex(0),
    bondFloorToday(0),
    historicBondFloor(dynBask->numSteps, 0.0){
        const ICouponsSPI* coupons = dynBask->dynBask->couponsSPI->getCouponsSPI().get();

        if (dynBask->dynBask->hasBFHistory) {
            int iStep;
            const DateTimeArray rebalDates = dynBask->dynBask->getRebalanceDates();

            for (iStep = 0; iStep < rebalDates.size() && rebalDates[iStep] < today; iStep++) {}
            firstCalculatedDateIndex = iStep;

            lastHistoricDateIndex = dynBask->dynBask->lastHistoricBFIndex;
            const CashFlowArray *BFHistory = dynBask->dynBask->bondFloorHistory.get();
            // note we've already validated the history so we know the dates match
            // the rebalance dates
            for (iStep = 0; iStep <= lastHistoricDateIndex; iStep++) {
                historicBondFloor[iStep] = (*BFHistory)[iStep].amount;
            }

            // if we're fixing today use the curve for fixings held in the bond
            if (rebalDates[firstCalculatedDateIndex].getDate()==today.getDate()) {
                fixingTodayIndex = firstCalculatedDateIndex;
                DoubleArray DFQuotFix = dynBask->dynBask->getDiscountFactorsQuotients(true, dynBask->bond->getYC());
                bondFloorFeeHelperForFixingToday = 
                    IBondFloorFeeHelperSP(BondFloorFeeHelper::make(dynBask, Basis, DFQuotFix));
                fixedPartForFixingToday = preCalc(dynBask, coupons, bondFloorFeeHelperForFixingToday, DFQuotFix);
            }
            if (lastHistoricDateIndex + 1 < firstCalculatedDateIndex) {
                // Handle the fiddly case when we roll
                DoubleArray DFQuotRoll = dynBask->dynBask->unpickDiscountFactorQuotients();

                bondFloorFeeHelperForRoll = 
                    IBondFloorFeeHelperSP(BondFloorFeeHelper::make(dynBask, Basis, DFQuotRoll));
                fixedPartForRoll = preCalc(dynBask, coupons, bondFloorFeeHelperForRoll, DFQuotRoll);
            }
        }

        // use class member for this so it is cached and the memory can be reused for possible "refresh"
        DFQuot = dynBask->dynBask->getDiscountFactorsQuotients(true, dynBask->disc);
        bondFloorFeeHelper = IBondFloorFeeHelperSP(BondFloorFeeHelper::make(dynBask, Basis, DFQuot));

        fixedPart = preCalc(dynBask, coupons, bondFloorFeeHelper, DFQuot);
    }

    const double getLevel(double BL,
                            int    iStep) const {
        if (iStep == fixingTodayIndex) {
            bondFloorToday = bondFloorFeeHelperForFixingToday->adjustNotionalForFees(BL * dynBask->Z(iStep), iStep) 
                                + fixedPartForFixingToday[iStep];
            return bondFloorToday;
        }
        if (iStep >= firstCalculatedDateIndex) {
            return bondFloorFeeHelper->adjustNotionalForFees(BL * dynBask->Z(iStep), iStep) 
                        + fixedPart[iStep];
        }
        if (iStep > lastHistoricDateIndex) {
            // we must have rolled over this date
            return bondFloorFeeHelperForRoll->adjustNotionalForFees(BL * dynBask->Z(iStep), iStep) 
                        + fixedPartForRoll[iStep];
        }
        else {
            return historicBondFloor[iStep];
        }
    };
    
    const double getLevelToday() const {
        // note if we have non bond floor history this returns 0
        return bondFloorToday;
    };

private:
    const SPIRunTime* dynBask;
    DoubleArray fixedPart;
    int lastHistoricDateIndex;
    int fixingTodayIndex;
    int firstCalculatedDateIndex;
    mutable double bondFloorToday;
    DoubleArray historicBondFloor;
    IBondFloorFeeHelperSP bondFloorFeeHelper;
    DoubleArray           DFQuot; // cache for the 'bondFloorFeeHelper'
    IBondFloorFeeHelperSP bondFloorFeeHelperForFixingToday;
    DoubleArray fixedPartForFixingToday;
    IBondFloorFeeHelperSP bondFloorFeeHelperForRoll;
    DoubleArray fixedPartForRoll;
};

ISPIBondFloor* SPIBondFloor::make(SPIRunTime* dynBaskRT,
                            double                        Basis,
                            const DateTime                today) {
    const Schedule* bondFloorSched = dynBaskRT->bond->getLinearBondFloor();
    if (bondFloorSched) {
        // this overrides everything
        return new FixedSPIBondFloor(dynBaskRT,
                                    bondFloorSched);
    }

    return new StandardSPIBondFloor(dynBaskRT, Basis, today);
}

DRLIB_END_NAMESPACE
