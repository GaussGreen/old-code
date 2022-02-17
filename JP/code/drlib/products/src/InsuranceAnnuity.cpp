//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : InsuranceAnnuity.cpp
//
//   Author      : Francois Lu
//
//   Description   Insurance annuity payoff
//
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#include "edginc/Generic1Factor.hpp"
#include "edginc/MonteCarlo.hpp"
#include "edginc/Maths.hpp"
#include "edginc/LinearStrikeVolRequest.hpp"
#include "edginc/AssetUtil.hpp"
#include "edginc/LastSensDate.hpp"
#include "edginc/SampleList.hpp"
DRLIB_BEGIN_NAMESPACE

static long Seed=-10;
static double paymentStorage=0;


// debug
//#define DEBUG_FILE_NAME "c:\\temp\\payoff_debug.txt"

//static FILE* debug_file = 0;

/** Rolling average american */
//////// instrument class //////////
class InsuranceAnnuity: public Generic1Factor, 
                        virtual public LastSensDate,
                        virtual public IMCIntoProduct{
protected:
    /// fields ////////
    SampleListSP            spotSamples;    // all historic spots
    double                  initialPayRate; // initial percentage annuity payment 
    double                  guaranteeAmount; // total benefit sum garanteed
    int                     deathBenefitAgeLimit; // age limit for MAV death benefit calculation
    int                     initialWaitPeriod;  // num of years to wait before going to normal period if step-up is made
    DoubleArray             mortalityRate;      // mortality rate as num of years
    DoubleArray             lapseRate;          // contract lapse rate
    IntArray                ageAtStart;         // age group when contract start
    DoubleArray             groupSizeAtStart;     // percentage of people in the corresponding age group
    bool                                    priceFees;                              // true=price fees, false=price structure
    double                                  feeRate;                            // the rate of fee payable as percentage of asset value
    bool                                    alwaysPayAfterWaitingPeriod;  // always pay coupon after waiting period, regardless of ITM
    bool                                    moyFee;//the fees are computed as the moy of the sum from last period and the actuel sum times the fee_Rate
    double                                  initialAmount;//initial amount of the fund
    bool                                    dyeOption; // the percentage of people who dye is given by the mortality rate
    bool                                    lapseOption;//the percentage of people who lapse is given by the lapseRate
    bool                                    stepUpOption;//if yes, we can't take coupon till waitingPeriod, but we can stepUp
    bool                                    OptionChosen;//if yes, the assumption of the pricer for the withdrawal and the lapse rate are the assumptions of Milliman

public:
    static CClassConstSP const TYPE;
    friend class InsuranceAnnuityProd;

    virtual void Validate(){
        static const string routine("InsuranceAnnuity::Validate");
        // NB don't call parent's validate - issues with fwdStart
        if (fwdStarting){
            throw ModelException(routine, "Fwd starting flag on "
                                 "Generic1Factor is not used");
        }
        if (oneContract){
            throw ModelException(routine, "oneContract flag on "
                                 "Generic1Factor is not used");
        }

        if (ageAtStart.size() <=0) {
            throw ModelException(routine, "ageAtStart array size must be >0.");
        }

        if (ageAtStart.size() != groupSizeAtStart.size()) {
            throw ModelException(routine, "ageAtStart and groupSizeAtStart arrays must be the same size.");
        }

        if (spotSamples->getDates().size() == 0) {
            throw ModelException(routine, "there must be at least one sample date.");
        }

        if (spotSamples->getFirstDate() > valueDate) {
            throw ModelException(routine, "there must be one sample date in the past (or = valueDate).");
        }

        AssetUtil::assetCrossValidate(asset.get(),
                                      false, //fwdStarting,
                                      startDate,
                                      valueDate,
                                      discount,
                                      this);
    }

    void validatePop2Object(){
        static const string routine("InsuranceAnnuity::validatePop2Object");

#ifdef DEBUG_FILE_NAME
        debug_file = fopen(DEBUG_FILE_NAME, "w");
#endif
    }
   
    /** Implementation of MonteCarlo::IntoProduct interface - the
        implementation of this is below */
    virtual IMCProduct* createProduct(
        const MonteCarlo* model) const;

    /** when to stop tweaking - default implementation assumes product can
        be priced with a MC */
    DateTime endDate(const Sensitivity*) const{
        DateTime end = valueDate.rollDate(365*(mortalityRate.size() - ageAtStart[0]));
        return end;
    }

private:
    InsuranceAnnuity(): Generic1Factor(TYPE) {}

    // for reflection
    InsuranceAnnuity(const InsuranceAnnuity& rhs); // not implemented
    InsuranceAnnuity& operator=(const InsuranceAnnuity& rhs); // not implemented

    static IObject* defaultInsuranceAnnuity(){
        return new InsuranceAnnuity();
    }

    //// roll through time (setting historic values)
    bool sensShift(Theta* theta){
        // use valueDate before it changes
        spotSamples->roll(theta->getUtil(valueDate), 0 /* iAsset */,
                          asset.get()); // then roll our past values
        Generic1Factor::sensShift(theta); // and then call parent's method
        return true; // continue to tweak components which implement Theta
    }

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(InsuranceAnnuity, clazz);
        SUPERCLASS(Generic1Factor);
        IMPLEMENTS(IMCIntoProduct);
        IMPLEMENTS(LastSensDate);   
        IMPLEMENTS(Theta::Shift);
        EMPTY_SHELL_METHOD(defaultInsuranceAnnuity);
        FIELD(spotSamples, "historical values");
        FIELD(initialPayRate,    "initial percentage annuity payment ");
        FIELD(guaranteeAmount, "total initial benefit sum garanteed");
        FIELD(deathBenefitAgeLimit, "age limit for MAV death benefit calculation");
        FIELD(initialWaitPeriod, "num of years to wait before going to normal period if step-up is made");
        FIELD(mortalityRate, "mortality rate over 1000 people as a number between 0 and 1000");
        FIELD(lapseRate, "contract lapse rate num between 0 and 1");
        FIELD(ageAtStart, "age group when contract start");
        FIELD(groupSizeAtStart, "percentage of people in the corresponding age group");
        FIELD(priceFees, "true=price fees, false=price structure");
        FIELD(feeRate, "the rate of fee payable as percentage of asset value");
        FIELD(alwaysPayAfterWaitingPeriod, "always pay coupon after waiting period, regardless of ITM");
        FIELD(moyFee, "the fees are computed as the moy of the sum from last period and the actuel sum times the fee_Rate");
        FIELD(initialAmount, "initial amount of the fund");
        FIELD(dyeOption, "the percentage of people is given by the mortalityRate");
        FIELD(lapseOption, "the percentage of people who lapse is given by the lapseRate ");
        FIELD(stepUpOption, "if yes, we can't take coupon till waitingPeriod, but we can stepUp ");
        FIELD(OptionChosen, "if yes, the assumption of the pricer for the withdrawal and the lapse rate are the assumptions of Milliman");


        clazz->setPublic(); // make visible to EAS/spreadsheet
    }
};

/* MC product class for super rainbow */
//////// product class //////////
class InsuranceAnnuityProd : public IMCProduct, virtual public IMCProductLN{
private:
    const InsuranceAnnuity*   inst; // reference to original instrument
    DateTimeArray             simulationDates;

    int             yrIndex;    // just needed for year until initialWaitPeriod years
    DoubleArray     histSpots;
    DateTimeArray   histDates;
    DoubleArray     growthFactors; // growth factor fom payout date to last date

public:
    
    /** Use this opportunity to do any LogNormal driven initialisation
        of the instrument before the main MC loop. e.g closed form barrier adjustment */
    void initialiseLN(const  IMCPathGenerator*  pathGen)const{
        // empty
    }

    /** equivalent to InstIntoMCProduct. */
    InsuranceAnnuityProd(const InsuranceAnnuity*         inst,
                         IRefLevelSP refLevel,
                         const SimSeriesSP&             simSeries) :
        IMCProduct(inst->asset.get(),
                  inst->valueDate,
                  inst->discount.get(),
                  refLevel,
                  simSeries,
                  inst->spotSamples,
                  inst->instSettle.get(),
                  simSeries->getLastDate()),
        inst(inst){
                    
        simulationDates = simSeries->getAllDates();

        histDates = inst->spotSamples->getDates();
        for (yrIndex=1; yrIndex < histDates.size(); yrIndex++)
        {
            if (histDates[yrIndex] > inst->valueDate)
                break;
        }
        histSpots = inst->spotSamples->getValues();

        // compute df for payment dates
        growthFactors.resize(histDates.size());
        DateTime lastDate = histDates[histDates.size() -1]; // this is the last date
        for (int i=0; i<histDates.size(); i++)
        {
            if(inst->valueDate <= histDates[i])
                growthFactors[i] = 1/inst->instSettle->pv(histDates[i],  
                                                          lastDate,
                                                          inst->discount.get(),
                                                          inst->asset.get());
            else
                growthFactors[i] = 99999999; // this number should never be used, or a bug
        }
    }

    /** print extra output **/
    virtual void recordExtraOutput(Control* control, Results* results, const IMCPrices& ) const;
    /** vol interp for LN */
    CVolRequestLNArray getVolInterp(const IMCPathGenerator* pathGen, int iAsset) const;


    void payoff(const IPathGenerator*  pathGen, IMCPrices&  prices);

    double calcDeathBenefit(double AV, double MAVDB, int index);

    double calcAlreadystepUp(double Spot,int index, double mySum, bool flag, double remtopay, double nextcoupon, double MAVDB, int age,const double *path, double i_end);

    double calcNostepUp(double Spot,int index, double mySum, bool flag, double remtopay, double nextcoupon, double MAVDB, int age,const double *path,bool take_coupon, double i_end);

    double alreadyTookCoupon(double Spot,int index, double mySum, bool flag, double remtopay, double nextcoupon, double MAVDB, int age,const double *path, double i_end);

    double noCouponTaken(double Spot,int index, double mySum, bool flag, double remtopay, double nextcoupon, double MAVDB, int age,const double *path, double i_end);

    double calcAfterInitialWaitPeriod(double Spot,int index, double mySum, bool flag, double remtopay, double nextcoupon, double MAVDB, int age,const double *path, double i_end);

    double calcFirstYear(int age,const double *path, int i_end);


}; 
// end of class InsuranceAnnuityProd

double InsuranceAnnuityProd::calcDeathBenefit(double AV, double MAVDB,int index)
{
    return 0.0;
    //growthFactors[index]*Maths::max(MAVDB-AV,0.0);
    //didn't take into consideration for GMWB for AMEX
}

double InsuranceAnnuityProd::calcAlreadystepUp(double Spot,int index, double mySum, bool flag, double remtopay, double nextcoupon, double MAVDB, int age,const double *path, double i_end)
{
    if (index>=(i_end-1)) return 0.0;
    double result;
    if (index>= inst->initialWaitPeriod) 
        result= calcAfterInitialWaitPeriod(Spot,index,mySum,flag,remtopay,nextcoupon,MAVDB,age,path,i_end);
    else{
        index=index+1;
        double oldSum=mySum;
        double nbShare=mySum/Spot;
        Spot=path[index];
        double AV=nbShare*Spot;
        double newFee=(inst->moyFee)?
            (
                (oldSum+AV)/2
                )
            :
            AV;
        double RBA=remtopay;
        double ITM=(RBA/AV-1);
        double ITMMilliman=AV/RBA;
        double aMilliman=(ITMMilliman>1)?0.22:2.3;
        double dead_rate=(inst->mortalityRate[age+index-1])/1000;
        if (!inst->dyeOption) dead_rate*=Maths::max(0.0,Maths::min(1.0,2-RBA/MAVDB));
        double deathBenefit = calcDeathBenefit(AV, MAVDB,index);
        if ((dead_rate<1)&&((index<histDates.size()))){
            double lapse_rate;
            if (inst->lapseOption) lapse_rate= Maths::max(0.0,Maths::min(1.0,inst->lapseRate[index-1]));
            else lapse_rate=(inst->OptionChosen)?
                     (
                         (
                             Maths::min(Maths::max(exp(aMilliman*(ITMMilliman-1)),0.2),1.25)
                             )*
                         (inst->lapseRate[index-1])
                         )
                     :
                (
                    Maths::max(0.2,Maths::min(1.0,1.0-1.5*(ITM-0.1)))*(inst->lapseRate[index-1])
                    );
            double normalCoupon=nextcoupon;
            bool stepUp=ITM<0;
            mySum=AV;
            bool old_flag=flag;
            flag=old_flag;
            double temporary=RBA; 
            remtopay=Maths::max(mySum,temporary);
            if ((age+index)<=inst->deathBenefitAgeLimit) 
                MAVDB=Maths::max(MAVDB,mySum);
            nextcoupon=Maths::min(
                (
                    (stepUp)?
                    Maths::max(0.01*(inst->initialPayRate)*remtopay,normalCoupon)
                    :
                    normalCoupon
                    )
                ,remtopay
                );
            result= inst->priceFees ? 
                (
                    (1-dead_rate)*(1-lapse_rate)*
                    (
                        (
                            (yrIndex<=index)?
                            (growthFactors[index]*newFee*inst->feeRate)
                            :
                            0.0
                            )
                        +
                        calcAlreadystepUp(Spot,index,mySum,flag,remtopay,nextcoupon,MAVDB,age,path,i_end)
                        )
                    )
                :
                (
                    (
                        (yrIndex<=index)?
                        (dead_rate*deathBenefit)
                        :
                        0.0
                        )
                    +
                    (1-dead_rate)*(1-lapse_rate)*calcAlreadystepUp(Spot,index,mySum,flag,remtopay,nextcoupon,MAVDB,age,path,i_end)
                    );
            if(index==yrIndex-1) paymentStorage +=(inst->priceFees) ? (1-dead_rate)*(1-lapse_rate)*newFee*inst->feeRate: 0.0;///0.0 because no coupon paid
                                     }
        else result=  inst->priceFees ? 
                 0.0:
            ((yrIndex<=index)?deathBenefit:0.0);
    }
    return result;
}

double InsuranceAnnuityProd::calcNostepUp(double Spot,int index, double mySum, bool flag, double remtopay, double nextcoupon, double MAVDB, int age,const double *path,bool take_coupon, double i_end)
{       if (take_coupon)
    return alreadyTookCoupon(Spot,index,mySum,flag,remtopay,nextcoupon,MAVDB,age,path,i_end);
 else
     return noCouponTaken(Spot,index,mySum,flag,remtopay,nextcoupon,MAVDB,age,path,i_end);
}

double InsuranceAnnuityProd::alreadyTookCoupon(double Spot,int index, double mySum, bool flag, double remtopay, double nextcoupon, double MAVDB, int age,const double *path, double i_end)
{       
    if (index>=(i_end-1)) return 0.0;
    double result;
    if (index>=(inst->initialWaitPeriod)) 
        result=calcAfterInitialWaitPeriod(Spot,index,mySum,flag,remtopay,nextcoupon,MAVDB,age,path,i_end); 
    else{
        index=index+1;
        double oldSum=mySum;
        double nbShare=mySum/Spot;
        Spot=path[index];
        double AV=nbShare*Spot;
        double newFee=(inst->moyFee)? 
            (
                (oldSum+AV)/2
                )
            :
            AV;
        double RBA=remtopay;
        double ITM=(RBA/AV-1);
        double ITMMilliman=AV/RBA;
        double aMilliman=(ITMMilliman>1)?0.22:2.3;
        double dead_rate=inst->mortalityRate[age+index-1]/1000;
        if (!inst->dyeOption) dead_rate*=Maths::max(0.0,Maths::min(1.0,2-RBA/MAVDB));
        double deathBenefit = calcDeathBenefit(AV, MAVDB,index);
        if ((dead_rate<1)&&((index<histDates.size()))){
            double lapse_rate;
            if (inst->lapseOption) lapse_rate= Maths::max(0.0,Maths::min(1.0,inst->lapseRate[index-1]));
            else lapse_rate=(inst->OptionChosen)?
                     (
                         (
                             Maths::min(
                                 Maths::max(exp(aMilliman*(ITMMilliman-1)),0.2)
                                 ,1.25
                                 )
                             )*
                         (inst->lapseRate[index-1])
                         )
                     :
                (
                    Maths::max(0.2,Maths::min(1.0,1.0-1.5*(ITM-0.1)))*
                    (inst->lapseRate[index-1])
                    );
            double normalCoupon=nextcoupon;
            bool take_coupon=(inst->OptionChosen)? 
                true
                :
                (
                    ran1(&Seed)<= Maths::min(1.0,Maths::max(0.2,5*ITM))
                    );
            double coupon=take_coupon?
                (
                    normalCoupon*
                    (
                        (inst->OptionChosen)? 
                        Maths::max(2.0/7.0,Maths::min(1.0,(1.0-(ITMMilliman-0.55)/0.49)))
                        :
                        1.0
                        )
                    )
                :
                0.0;
            mySum=Maths::max(AV-coupon,0.0);
            bool old_flag=flag;
            flag=mySum>0.0;
            double Payoff=take_coupon?
                (
                    old_flag ?
                    (
                        flag ? 0.0: (coupon-AV)
                        )
                    :
                    coupon
                    )
                :
                0.0;
            double temporary=RBA-coupon;
            remtopay= temporary;
            if ((age+index)<=inst->deathBenefitAgeLimit) 
                MAVDB=Maths::max(MAVDB-coupon*MAVDB/AV,mySum);
            else
                MAVDB=MAVDB-coupon*MAVDB/AV;
            if (remtopay==0.0) {
                result=inst->priceFees ? 
                    (
                        (yrIndex<=index)?
                        (
                            (1-dead_rate)*(1-lapse_rate)*
                            growthFactors[index]*newFee*inst->feeRate
                            )
                        :
                        0.0
                        )
                    : 
                    (
                        (yrIndex<=index)?
                        (
                            dead_rate*deathBenefit+
                            (1-dead_rate)*(1-lapse_rate)*growthFactors[index]*Payoff
                            )
                        :
                        0.0
                        );
                if(index==yrIndex-1) paymentStorage += (inst->priceFees) ?
                                         (
                                             (1-dead_rate)*(1-lapse_rate)*newFee*inst->feeRate
                                             )
                                         :
                    (
                        (1-dead_rate)*(1-lapse_rate)*Payoff
                        );
            }
            else{
                nextcoupon=Maths::min(normalCoupon, remtopay);
                result= inst->priceFees ? 
                    (
                        (1-dead_rate)*(1-lapse_rate)*
                        (
                            (
                                (yrIndex<=index)?
                                (
                                    growthFactors[index]*newFee*inst->feeRate
                                    )
                                :
                                0.0
                                )+
                            alreadyTookCoupon(Spot,index,mySum,flag,remtopay,nextcoupon,MAVDB,age,path,i_end)
                            )
                        )
                    :
                    (
                        (
                            (yrIndex<=index)?
                            (dead_rate*deathBenefit):0.0
                            )+
                        (1-dead_rate)*(1-lapse_rate)*
                        (
                            (
                                (yrIndex<=index)?
                                (growthFactors[index]*Payoff):0.0
                                )+
                            alreadyTookCoupon(Spot,index,mySum,flag,remtopay,nextcoupon,MAVDB,age,path,i_end)
                            )
                        );
                if(index==yrIndex-1) paymentStorage +=(inst->priceFees) ? 
                                         (
                                             (1-dead_rate)*(1-lapse_rate)*newFee*inst->feeRate
                                             )
                                         :
                    (
                        (1-dead_rate)*(1-lapse_rate)*Payoff
                        );
            }
        }
        else result=inst->priceFees ? 
                 0.0:
            ((yrIndex<=index)?deathBenefit:0.0);
    }
    return result;
}

double InsuranceAnnuityProd::noCouponTaken(double Spot,int index, double mySum, bool flag, double remtopay, double nextcoupon, double MAVDB, int age,const double *path, double i_end)
{       
    if (index>=(i_end-1)) return 0.0;
    double result;
    if (index>=(inst->initialWaitPeriod)) 
        result=calcAfterInitialWaitPeriod(Spot,index,mySum,flag,remtopay,nextcoupon,MAVDB,age,path,i_end);
    else{
        index=index+1;
        double oldSum=mySum;
        double nbShare=mySum/Spot;
        Spot=path[index];
        double AV=nbShare*Spot;
        double newFee=(inst->moyFee)? ((oldSum+AV)/2):AV;
        double RBA=remtopay;
        double ITM=(RBA/AV-1);
        double ITMMilliman=AV/RBA;
        double aMilliman=(ITMMilliman>1)?0.22:2.3;
        double dead_rate=(inst->mortalityRate[age+index-1])/1000;
        if (!inst->dyeOption) dead_rate*=Maths::max(0.0,Maths::min(1.0,2-RBA/MAVDB));
        double deathBenefit = calcDeathBenefit(AV, MAVDB,index);
        if ((dead_rate<1)&&((index<histDates.size()))){
            double lapse_rate;
            if (inst->lapseOption) lapse_rate = Maths::max(0.0,Maths::min(1.0,inst->lapseRate[index-1]));
            else lapse_rate=(inst->OptionChosen)?
                     (
                         (
                             Maths::min(
                                 Maths::max(exp(aMilliman*(ITMMilliman-1)),0.2)
                                 ,1.25
                                 )
                             )*
                         (inst->lapseRate[index-1])
                         )
                     :
                (
                    Maths::max(0.2,Maths::min(1.0,1.0-1.5*(ITM-0.1)))
                    *(inst->lapseRate[index-1])
                    );
            double normalCoupon=nextcoupon;
            bool stepUp=ITM<0;
            if(stepUp){
                mySum=AV;
                bool old_flag=flag;
                flag=old_flag;
                remtopay=mySum;
                if ((age+index)<=inst->deathBenefitAgeLimit) 
                    MAVDB=Maths::max(MAVDB,mySum);
                nextcoupon=Maths::max(0.01*(inst->initialPayRate)*remtopay,(inst->initialPayRate));
                result= inst->priceFees ? 
                    (
                        (1-dead_rate)*(1-lapse_rate)*
                        (
                            (
                                (yrIndex<=index)?
                                (growthFactors[index]*newFee*inst->feeRate):0.0
                                )+
                            calcAlreadystepUp(Spot,index,mySum,flag,remtopay,nextcoupon,MAVDB,age,path,i_end)
                            )
                        )
                    :
                    (
                        (
                            (yrIndex<=index)?
                            (dead_rate*deathBenefit):0.0
                            )+
                        (1-dead_rate)*(1-lapse_rate)*
                        calcAlreadystepUp(Spot,index,mySum,flag,remtopay,nextcoupon,MAVDB,age,path,i_end)
                        );
                if(index==yrIndex-1) paymentStorage +=(inst->priceFees) ?
                                         (1-dead_rate)*(1-lapse_rate)*newFee*inst->feeRate: 0.0;
            }
            else{
                bool take_coupon=(inst->OptionChosen)? true:
                    (ran1(&Seed)<= Maths::min(1.0,Maths::max(0.2,5*ITM)));
                double coupon=take_coupon?
                    (
                        normalCoupon*
                        (
                            (inst->OptionChosen)? 
                            Maths::max(
                                2.0/7.0,
                                Maths::min(1.0,(1.0-(ITMMilliman-0.55)/0.49))
                                )
                            :
                            1.0
                            )
                        )
                    :
                    0.0;
                mySum=Maths::max(AV-coupon,0.0);
                bool old_flag=flag;
                flag=mySum>0.0; 
                double Payoff= take_coupon?(old_flag? (flag ? 0.0: (coupon-AV)):coupon):0.0;
                double temporary=RBA-coupon;
                remtopay= temporary;
                if ((age+index)<=inst->deathBenefitAgeLimit) 
                    MAVDB=Maths::max(MAVDB-coupon*MAVDB/AV,mySum);
                else
                    MAVDB=MAVDB-coupon*MAVDB/AV;
                if (remtopay==0.0) {
                    result=inst->priceFees? 
                        (
                            (
                                (yrIndex<=index)?
                                (
                                    (1-dead_rate)*(1-lapse_rate)*
                                    growthFactors[index]*newFee*inst->feeRate
                                    )
                                :
                                0.0
                                )
                            )
                        :
                        (
                            (yrIndex<=index)?
                            (
                                dead_rate*deathBenefit+
                                (1-dead_rate)*(1-lapse_rate)*growthFactors[index]*Payoff
                                )
                            :
                            0.0
                            );
                    if(index==yrIndex-1) paymentStorage += 
                                             (inst->priceFees) ? 
                                             (1-dead_rate)*(1-lapse_rate)*newFee*inst->feeRate
                                             :
                        (1-dead_rate)*(1-lapse_rate)*Payoff;
                }
                else{
                    nextcoupon=Maths::min(normalCoupon, remtopay);
                    if (take_coupon) {
                        result=inst->priceFees ?
                            (
                                (1-dead_rate)*(1-lapse_rate)*
                                (
                                    (
                                        (yrIndex<=index)? 
                                        (growthFactors[index]*newFee*inst->feeRate):0.0
                                        )+
                                    alreadyTookCoupon(Spot,index,mySum,flag,remtopay,nextcoupon,MAVDB,age,path,i_end)
                                    )
                                )
                            :
                            (
                                (
                                    (yrIndex<=index)?
                                    (dead_rate*deathBenefit):0.0
                                    )+
                                (1-dead_rate)*(1-lapse_rate)*
                                (
                                    (
                                        (yrIndex<=index)?
                                        (growthFactors[index]*Payoff):0.0
                                        )+
                                    alreadyTookCoupon(Spot,index,mySum,flag,remtopay,nextcoupon,MAVDB,age,path,i_end)
                                    )
                                );
                        if(index==yrIndex-1) paymentStorage +=(inst->priceFees) ? (1-dead_rate)*(1-lapse_rate)*newFee*inst->feeRate: (1-dead_rate)*(1-lapse_rate)*Payoff;
                    }
                    else {
                        result=inst->priceFees ?
                            (
                                (1-dead_rate)*(1-lapse_rate)*
                                (
                                    (
                                        (yrIndex<=index)?(growthFactors[index]*newFee*inst->feeRate):0.0
                                        )+
                                    noCouponTaken(Spot,index,mySum,flag,remtopay,nextcoupon,MAVDB,age,path,i_end)
                                    )
                                )
                            :
                            (
                                (
                                    (yrIndex<=index)?
                                    (dead_rate*deathBenefit):0.0
                                    )+
                                (1-dead_rate)*(1-lapse_rate)*
                                noCouponTaken(Spot,index,mySum,flag,remtopay,nextcoupon,MAVDB,age,path,i_end)
                                );
                        if(index==yrIndex-1) paymentStorage +=(inst->priceFees) ? 
                                                 (
                                                     (1-dead_rate)*(1-lapse_rate)*newFee*inst->feeRate)
                                                 :
                            (
                                (1-dead_rate)*(1-lapse_rate)*Payoff
                                );
                    }
                }
            }
        }
        else result=inst->priceFees ? 
                 0.0:
            ((yrIndex<=index)?deathBenefit:0.0);
    }
    return result;
}

// compute when waiting period is done
double InsuranceAnnuityProd::calcAfterInitialWaitPeriod(double Spot,int index, double mySum, bool flag, double remtopay, double nextcoupon, double MAVDB, int age,const double *path,double i_end)
{
    if (index>=(i_end-1)) return 0.0;
    double result;
    index=index+1;
    double oldSum=mySum;
    double nbShare=mySum/Spot;
    Spot=path[index];
    double AV=nbShare*Spot;
    double newFee=(inst->moyFee)? 
        (
            (oldSum+AV)/2
            )
        :
        AV;
    double RBA=remtopay;
    double ITM=AV>0?(RBA/AV-1):100000000.0;
    double ITMMilliman=AV/RBA;
    double aMilliman=(ITMMilliman>1)?0.22:2.3;
    double dead_rate=(inst->mortalityRate[age+index-1])/1000;
    if (!inst->dyeOption) dead_rate*=Maths::max(0.0,Maths::min(1.0,2-RBA/MAVDB));
    double deathBenefit = calcDeathBenefit(AV, MAVDB,index);
    if ((dead_rate<1)&&((index<histDates.size()))){
        double lapse_rate;
        if(inst->lapseOption) lapse_rate=Maths::max(0.0,Maths::min(1.0,inst->lapseRate[index-1]));
        else lapse_rate=(inst->OptionChosen)?
                 (
                     (
                         Maths::min(
                             Maths::max(exp(aMilliman*(ITMMilliman-1)),0.2)
                             ,1.25
                             )
                         )*
                     (inst->lapseRate[index-1])
                     )
                 :
            (
                Maths::max(0.2,Maths::min(1.0,1.0-1.5*(ITM-0.1)))*
                (inst->lapseRate[index-1])
                );
        double normalCoupon=nextcoupon;
        bool stepUp=ITM<0;
        bool take_coupon;
        take_coupon=inst->alwaysPayAfterWaitingPeriod ? 
            true
            :
            (
                (inst->OptionChosen)?
                true
                :
                (
                    ran1(&Seed)<= Maths::min(1.0,Maths::max(0.2,5*ITM))
                    )
                );      
        double coupon=take_coupon?
            (
                normalCoupon*
                (
                    (inst->OptionChosen)?
                    Maths::max(2.0/7.0,
                               Maths::min(1.0,(1.0-(ITMMilliman-0.55)/0.49))
                        )
                    :
                    1.0
                    )
                )
            :
            0.0;    
        mySum=Maths::max(AV-coupon,0.0);
        double old_flag=flag;   
        flag=mySum>0.0; 
        double Payoff=take_coupon?
            (
                old_flag? 
                (
                    flag ? 
                    0.0
                    :
                    (coupon-AV)
                    )
                :
                coupon
                )
            :
            0.0;    
        double temporary=RBA-coupon;
        remtopay=Maths::max(mySum,temporary);
        if ((age+index)<=inst->deathBenefitAgeLimit) 
            MAVDB=Maths::max(MAVDB-coupon*MAVDB/AV,mySum);
        else 
            MAVDB=MAVDB-coupon*MAVDB/AV;
        if (remtopay==0.0) {
            result=inst->priceFees? 
                (
                    (yrIndex<=index)?
                    (
                        (1-dead_rate)*(1-lapse_rate)*growthFactors[index]*newFee*inst->feeRate
                        )
                    :
                    0.0
                    )
                :
                (
                    (yrIndex<=index)?
                    (
                        dead_rate*deathBenefit+(1-dead_rate)*(1-lapse_rate)*growthFactors[index]*Payoff
                        )
                    :0.0
                    );
            if(index==yrIndex-1) paymentStorage +=(inst->priceFees) ? 
                                     (
                                         (1-dead_rate)*(1-lapse_rate)*newFee*inst->feeRate
                                         )
                                     :
                (
                    (1-dead_rate)*(1-lapse_rate)*Payoff
                    );
        }
        else {
            nextcoupon=Maths::min(
                (
                    (stepUp)?
                    Maths::max(
                        0.01*(inst->initialPayRate)*remtopay
                        ,
                        normalCoupon
                        )
                    :
                    normalCoupon
                    ), 
                remtopay
                );
            result=inst->priceFees ? 
                (
                    (1-dead_rate)*(1-lapse_rate)*
                    (
                        (
                            (yrIndex<=index)?
                            (growthFactors[index]*newFee*inst->feeRate):0.0
                            )
                        + calcAfterInitialWaitPeriod(Spot,index,mySum,flag,remtopay,nextcoupon,MAVDB,age,path,i_end)
                        )
                    )
                :
                (
                    (
                        (yrIndex<=index)?
                        (dead_rate*deathBenefit):0.0
                        )
                    +(1-dead_rate)*(1-lapse_rate)*
                    (
                        (
                            (yrIndex<=index)?
                            (growthFactors[index]*Payoff):0.0
                            )
                        +calcAfterInitialWaitPeriod(Spot,index,mySum,flag,remtopay,nextcoupon,MAVDB,age,path,i_end)
                        )
                    );
            if(index==yrIndex-1) paymentStorage += 
                                     (inst->priceFees) ? 
                                     (
                                         (1-dead_rate)*(1-lapse_rate)*newFee*inst->feeRate
                                         )
                                     :
                (
                    (1-dead_rate)*(1-lapse_rate)*Payoff
                    );
        }
    }
    else result=inst->priceFees ? 
             0.0:
        ((yrIndex<=index)?deathBenefit:0.0);
    return result;
}

double InsuranceAnnuityProd::calcFirstYear(int age,const double *path, int i_end)
{
    // ** start from 1st year **********
    double result;
    double Spot;
    int index;
    double nbShare;
    double AV;
    double RBA;
    double ITM;
    double ITMMilliman;
    double aMilliman;
    double dead_rate;
    double deathBenefit;
    double lapse_rate;
    double normalCoupon;
    bool take_coupon;
    double coupon;
    double mySum;
    bool flag;
    double Payoff;
    double temporary;
    double remtopay;
    double nextcoupon;
    bool stepUp;                    
    double MAVDB;
    double oldSum;
    double newFee;
        
//debugging
#ifdef DEBUG_FILE_NAME
    fprintf(debug_file, "%f \n", paymentStorage);
#endif
    index=1;
    if (index>=i_end) return 0.0;
    oldSum=inst->initialAmount;
    nbShare = inst->initialAmount/histSpots[0];
    Spot=path[index];
    AV=nbShare*Spot;
    newFee=(inst->moyFee)? ((oldSum+AV)/2):AV;
    RBA=inst->guaranteeAmount;
    ITM=(RBA/AV-1);
    ITMMilliman=AV/RBA;
    aMilliman=(ITMMilliman>1)?0.22:2.3;
    MAVDB=0;
    dead_rate=inst->mortalityRate[age]/1000.0;
    if (!inst->dyeOption) dead_rate*=0;
    deathBenefit = calcDeathBenefit(AV, MAVDB,index);
    if ((dead_rate<1)&&((index<histDates.size()))){
        if(inst->lapseOption) lapse_rate=Maths::max(0.0,Maths::min(1.0,inst->lapseRate[index-1]));
        else lapse_rate=(inst->OptionChosen)?
                 (
                     (
                         Maths::min(
                             Maths::max(exp(aMilliman*(ITMMilliman-1))
                                        ,
                                        0.2
                                 )
                             ,
                             1.25
                             )
                         )*
                     (inst->lapseRate[index-1])
                     )
                 :
            (
                Maths::max(
                    0.2
                    ,
                    Maths::min(1.0,1.0-1.5*(ITM-0.1))
                    )*
                (inst->lapseRate[index-1])
                );
        normalCoupon = inst->initialPayRate;
        stepUp=ITM<0;
        if(inst->stepUpOption){
            mySum=AV;
            flag=true;
            remtopay=Maths::max(mySum,RBA);
            if ((age+index)<=inst->deathBenefitAgeLimit) 
                MAVDB=Maths::max(MAVDB,mySum);
            nextcoupon=Maths::min(
                (
                    (stepUp)?
                    Maths::max(0.01*(inst->initialPayRate)*remtopay,normalCoupon)
                    :
                    normalCoupon
                    )
                ,
                remtopay
                );
            result=inst->priceFees?
                (
                    (1-dead_rate)*(1-lapse_rate)*
                    (
                        (
                            (yrIndex<=index)?
                            (growthFactors[index]*newFee*inst->feeRate)
                            :
                            0.0
                            )+
                        calcAlreadystepUp(Spot,index,mySum,flag,remtopay,nextcoupon,MAVDB,age,path,i_end)
                        )
                    )
                :
                (
                    (
                        (yrIndex<=index)?
                        (dead_rate*deathBenefit)
                        :
                        0.0
                        )+
                    (1-dead_rate)*(1-lapse_rate)*
                    calcAlreadystepUp(Spot,index,mySum,flag,remtopay,nextcoupon,MAVDB,age,path,i_end)
                    );
            if(index==yrIndex-1) paymentStorage +=(inst->priceFees) ? 
                                     (
                                         (1-dead_rate)*(1-lapse_rate)*newFee*inst->feeRate
                                         )
                                     :
                0.0;
        }
        else{
            if (stepUp){
                take_coupon=false;
                coupon=0.0;
                mySum=AV;
                flag=true;
                Payoff=0.0;
                temporary=RBA;
                remtopay=mySum;
                if ((age+index)<=inst->deathBenefitAgeLimit) 
                    MAVDB=Maths::max(MAVDB,mySum);
                nextcoupon=Maths::min(Maths::max(0.01*(inst->initialPayRate)*remtopay,normalCoupon),remtopay);
                result=inst->priceFees?
                    (
                        (1-dead_rate)*(1-lapse_rate)*
                        (
                            (
                                (yrIndex<=index)?
                                (
                                    growthFactors[index]*newFee*inst->feeRate
                                    )
                                :
                                0.0
                                )+
                            calcAlreadystepUp(Spot,index,mySum,flag,remtopay,nextcoupon,MAVDB,age,path,i_end)
                            )
                        )
                    :
                    (
                        (
                            (yrIndex<=index)?
                            (dead_rate*deathBenefit)
                            :
                            0.0
                            )+
                        (1-dead_rate)*(1-lapse_rate)*
                        calcAlreadystepUp(Spot,index,mySum,flag,remtopay,nextcoupon,MAVDB,age,path,i_end)
                        );
                if(index==yrIndex-1) paymentStorage +=(inst->priceFees) ? 
                                         (
                                             (1-dead_rate)*(1-lapse_rate)*newFee*inst->feeRate
                                             )
                                         :
                    0.0;
            }

            else {
                take_coupon=(inst->OptionChosen)? 
                    true
                    :
                    (
                        ran1(&Seed)<= Maths::min(1.0,Maths::max(0.2,5*ITM))
                        );
                coupon=take_coupon?
                    (
                        normalCoupon*
                        (
                            (inst->OptionChosen)? 
                            Maths::max(2.0/7.0,
                                       Maths::min(1.0,(1.0-(ITMMilliman-0.55)/0.49))
                                )
                            :
                            1.0
                            )
                        )
                    :
                    0.0;
                mySum=Maths::max(AV-coupon,0.0);        
                flag=mySum>0.0 ;
                Payoff=take_coupon?(flag ? 0.0: (coupon-AV)):0.0;
                temporary=RBA-coupon;
                remtopay=temporary ;
                if ((age+index)<=inst->deathBenefitAgeLimit) 
                    MAVDB=Maths::max(MAVDB-coupon*MAVDB/AV,mySum);
                else 
                    MAVDB=MAVDB-coupon*MAVDB/AV;
                if (remtopay==0.0) {
                    result=inst->priceFees?
                        (
                            (yrIndex<=index)?
                            (
                                (1-dead_rate)*(1-lapse_rate)*
                                growthFactors[index]*newFee*inst->feeRate
                                )
                            :
                            0.0
                            )
                        :
                        (
                            (yrIndex<=index)?
                            (
                                dead_rate*deathBenefit+
                                (1-dead_rate)*(1-lapse_rate)*growthFactors[index]*Payoff
                                )
                            :
                            0.0
                            );
                    if(index==yrIndex-1) paymentStorage +=(inst->priceFees) ? 
                                             (
                                                 (1-dead_rate)*(1-lapse_rate)*newFee*inst->feeRate
                                                 )
                                             :
                        (
                            (1-dead_rate)*(1-lapse_rate)*Payoff
                            );
                }
                else{
                    nextcoupon=Maths::min(normalCoupon, remtopay);
                    result=inst->priceFees ? 
                        (
                            (1-dead_rate)*(1-lapse_rate)*
                            (
                                (
                                    (yrIndex<=index)? 
                                    (growthFactors[index]*newFee*inst->feeRate)
                                    :
                                    0.0
                                    )+
                                calcNostepUp(Spot,index,mySum,flag,remtopay,nextcoupon,MAVDB,age,path,take_coupon,i_end)
                                )
                            )
                        :
                        (
                            (
                                (yrIndex<=index)?
                                (dead_rate*deathBenefit)
                                :
                                0.0
                                )+
                            (1-dead_rate)*(1-lapse_rate)*
                            (
                                (
                                    (yrIndex<=index)?
                                    (growthFactors[index]*Payoff)
                                    :
                                    0.0
                                    )+
                                calcNostepUp(Spot,index,mySum,flag,remtopay,nextcoupon,MAVDB,age,path,take_coupon,i_end)
                                )
                            );
                    if(index==yrIndex-1) paymentStorage += (inst->priceFees) ?
                                             (
                                                 (1-dead_rate)*(1-lapse_rate)*newFee*inst->feeRate
                                                 )
                                             :
                        (
                            (1-dead_rate)*(1-lapse_rate)*Payoff
                            );
                }
            }
        }
    }
    else result=inst->priceFees ? 
             0.0:
        ((yrIndex<=index)?deathBenefit:0.0);
    return result;
}
void InsuranceAnnuityProd::payoff(const IPathGenerator*  pathGen, IMCPrices& prices) {

    const double* path = pathGen->Path(0,0); // access path
    int age;
    double pathValue = 0.0;

    
    int i_end = pathGen->end(0);

    if (pathGen->doingPast())
    {
        // compute realised ITM, Deaths, past lapses
    }
    else
    {
        // compute payoff for each age group
        for (int ageGroup =0; ageGroup<inst->ageAtStart.size(); ageGroup++)
        {
            double value = 0.0;
            age = inst->ageAtStart[ageGroup];
            value =calcFirstYear(age,path,i_end);
            pathValue+=value*inst->groupSizeAtStart[ageGroup];
        }
    }

    if (!pathGen->doingPast())
    {
        prices.add(pathValue);
        
    }
}

// for the LogNormal path generator
CVolRequestLNArray InsuranceAnnuityProd::getVolInterp(const IMCPathGenerator* pathGen,
                                                      int                     iAsset) const {
    DateTime imntStartDate = inst->fwdStarting? 
        inst->startDate: inst->valueDate;
    CVolRequestLNArray   reqarr(1); // one interp level/path per asset here

    reqarr[0] = CVolRequestLNSP(new LinearStrikeVolRequest(
        inst->asset->getSpot(),  //ATM
        imntStartDate, 
        inst->endDate(0),
        inst->fwdStarting));
    
    return reqarr;
}

// control variate is done here
void InsuranceAnnuityProd::recordExtraOutput(Control* control, Results* results, const IMCPrices& ) const
{
    OutputNameConstSP valueout(new OutputName("payment"));
    results->storeGreek(IObjectSP(CDouble::create(paymentStorage)), Results::DEBUG_PACKET, valueout);
    paymentStorage=0;
    Seed=-10;
}

/** Implementation of MonteCarlo::IntoProduct interface */
IMCProduct* InsuranceAnnuity::createProduct(const MonteCarlo* model) const {
    // we need to create a SimSeries object which says which assets need
    // which dates to be simulated
    SimSeriesSP simSeries(new SimSeries(1)); /* create empty one */
    simSeries->addDates(spotSamples->getAllDates());

    IRefLevelSP refLevel(IRefLevel::Util::makeTrivialAverage(spotSamples->getDates()[0]));
    return new InsuranceAnnuityProd(this, refLevel, simSeries);
}

CClassConstSP const InsuranceAnnuity::TYPE = CClass::registerClassLoadMethod(
    "InsuranceAnnuity", typeid(InsuranceAnnuity), InsuranceAnnuity::load);

// force linker to include this file (avoid having header file) */
extern bool InsuranceAnnuityLoad()
{
    return true && InsuranceAnnuity::TYPE;
}

DRLIB_END_NAMESPACE
