
//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : InsuranceA.cpp
//
//   Author      : Xiaolan Zhang
//
//   Description   Insurance annuity payoff
//                 common features for GMDB/GMAB/GMWB
//   Date        : Sept 2006
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#include "edginc/GenericNFBase.hpp"
#include "edginc/MonteCarlo.hpp"
#include "edginc/LinearStrikeVolRequest.hpp"
#include "edginc/AssetUtil.hpp"
#include "edginc/LastSensDate.hpp"
#include "edginc/SampleList.hpp"
#include "edginc/Random.hpp"
#include "edginc/SVGenSpot.hpp"
#include "edginc/SVGenDiscFactor.hpp"


DRLIB_BEGIN_NAMESPACE

//////// instrument class ////////
class InsuranceA: public GenericNFBase, 
                         virtual public LastSensDate,
                         virtual public IMCIntoProduct {

public:
    static CClassConstSP const TYPE;
    friend class InsuranceAProd;
    friend class InsuranceAProdSV;

    virtual void Validate();

    /** Implementation of MonteCarlo::IntoProduct interface - the
        implementation of this is below */
    virtual IMCProduct* createProduct(
        const MonteCarlo* model) const;

    /** when to stop tweaking - default implementation assumes product can
        be priced with a MC */
    DateTime endDate(const Sensitivity*) const {
        DateTime end = valueDate.rollDate(365 * (mortalityRate.size() - ageAtStart[0]));
        return end;
    }

protected:

    virtual const DateTimeArray samplingDates() const ;
    // constructor 
    InsuranceA(CClassConstSP const type=TYPE): GenericNFBase(type), 
                                feeRateIns(0.0),
                                isDynamicLapseRate(false),
                                lapseRateMultiplierType(1),                            
                                withdrawMultiplier(1.0),
                                noWithdrawMultiplier(1.0){
        featureGpSizeAtStart.resize(1,1.0);  //withdrawal group percentage.
    }

    // for reflection
    InsuranceA(const InsuranceA& rhs); // not implemented
    InsuranceA& operator=(const InsuranceA& rhs); // not implemented

    static IObject* defaultInsuranceA() {
        return new InsuranceA();
    }

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz);

    ////////////////// fields ///////////////
    
    DoubleArray       weights;              //weight for basket
    DateTimeArray     monitorDates;         // monitoring dates
    double            feeRateJPM;           // fee rate as % of the account value, charged by JPM
    double            feeRateIns;           // fee rate as % of the account value, charged by Insurance company


    /** age groups with percentages and age limit for benefit accumulation   */
    IntArray        ageAtStart;              // age of the people in the age group at start date
    DoubleArray     ageGpSizeAtStart;        // percentage of people in the corresponding group
    //to review
    DoubleArray     featureGpSizeAtStart;       // withdrawal group percentage.?????to chg to group features


    /** mortality and lapse rate definitions*/
    bool            dieOption;              // if yes, the mortality rates table is taken 
                                            // into account
                                            // if no, the mortality rate is 0%
    DoubleArray     mortalityRate;          // mortality rate as function of the client age
    bool            lapseOption;            // if yes, the lapse rate table is taken into account
                                            // if no, the mortality rate is 0%
    DoubleArray     lapseRate;              // contract lapse rate as function of the 
                                            // time from start date
    bool            isDynamicLapseRate;     // if yes, the lapse rate is dynamic 
                                            // (i.e.determined with respect to ITM)

    //to review, only for GMWB or others also?????
    int             lapseRateMultiplierType;// default: 1 (old version, power ft); 
                                            // 2: AV/Benefits if Benefit/AV >=1
    double          withdrawMultiplier;     // additional multiplier depends on the 
                                            // current status of withdrawal.
    double          noWithdrawMultiplier;

};

///////////////////// product class //////////////////////////////

//////////////////////////////////////////////////////////////////////////////////
//state variable Version
//////////////////////////////////////////////////////////////////////////////////

/* MC product class for insurance annuity */
class InsuranceAProdSV : public MCProductClient,
                        virtual public IMCProductLN {

public:
    class InsuranceAPrices;
    typedef refCountPtr<InsuranceAPrices> InsuranceAPricesSP;

    class InsuranceAPrices: public IMCPricesSimple {
    public:

        enum PricesIndices {
            OPTION_PRICE = 0,
            FEE_PRICE,
            NB_PRICES
        };

        /** adds supplied price to this set of IMCPrices */
        void add(double price, int index) {
            simplePrices[index]->add(price);
        }

        /** Returns the averaged result, together with the standard error */
        virtual void getResult(double* result, double* resultStdErr) const {
            simplePrices[OPTION_PRICE]->getResult(result, resultStdErr);
        }

        double getOptionPrice() const {
            double optionPrice, dummy;
            simplePrices[OPTION_PRICE]->getResult(&optionPrice, &dummy);
            return optionPrice;
        }

        double getFeePrice() const {
            double feePrice, dummy;
            simplePrices[FEE_PRICE]->getResult(&feePrice, &dummy);
            return feePrice;
        }

        /** Returns true if the path, identified by pathIdx, should be
        repriced. If false is returned then there will be no add()
        method called for this path and the IMCPrices object must
        take any appropriate action */
        virtual bool repriceForGreek(int pathIdx) {
            return true;
        }

        /** Support product-level caching across tweaks */
        virtual int storagePerPath(IMCProduct* product) const {
            return 0;
        }

        virtual void configureCache(const IntArray& changedAssets) {
        }

        /** Returns a deep copy of this object */
        IMCPrices* clone() const {
            InsuranceAPricesSP copy(new InsuranceAPrices());
            copy->simplePrices.resize(simplePrices.size());
            for (unsigned int iPrice = 0 ; iPrice < copy->simplePrices.size() ; iPrice++) {
                IMCPricesSP temp(this->simplePrices[iPrice]->clone());
                copy->simplePrices[iPrice] = MCPricesSimpleSP(static_cast<MCPricesSimple*>(temp.get()));
            }
            return copy.get();
        }

        InsuranceAPrices(int NbIter, int NbSubSamples):
        simplePrices(NB_PRICES){
            for (unsigned int iPrice = 0 ; iPrice < simplePrices.size() ; iPrice++){
                simplePrices[iPrice] = MCPricesSimpleSP(new MCPricesSimple(NbIter, NbSubSamples));
            }
        }

    private:
        InsuranceAPrices() {
        }

        /** adds supplied price to this set of IMCPrices */
        virtual void add(double price) {
            throw ModelException("IMCPrices::add(double)",
                     "internal error");
        }

        /** Returns the last price stored. Undefined behaviour if no
        prices yet stored */
        virtual double lastPrice() const {
            throw ModelException("IMCPrices::lastPrice",
                     "internal error");
        }

        /** On pricing run returns MAX(x, 0.0). It should be used for applying
        the 'final point of optionality' within the product. This allows
        QuickGreeks type of IMCPrices not to apply the max when doing first
        order greeks (this may sound strange but see the docs for why) */
        virtual double maxWithZero(double x) const {
            throw ModelException("IMCPrices::maxWithZero",
                     "internal error");
        }

        /** Reset this object so that it can be used for the same operation
        again. Normally, a new IMCPrices object is created for each
        pricing run. However, for quick x gamma, it is important to use
        the same one for each pair of assets */
        virtual void reset() {
            throw ModelException("IMCPrices::reset",
                     "internal error");
        }

        /** Ease cloning */
        virtual IMCPrices* emptyConstructor() const {
            throw ModelException("IMCPrices::emptyConstructor",
                     "internal error");
        }

        vector<MCPricesSimpleSP> simplePrices;

    };


public:
    
    virtual IMCPrices* createOrigPrices(int nbIter,
                             int nbSubSamples,
                             int mode) {
        //to review??????  check with basket rebalance 
        return new InsuranceAPrices(nbIter, nbSubSamples);
    }

    /** Use this opportunity to do any LogNormal driven initialisation
    of the instrument before the main MC loop. e.g closed form barrier adjustment */
    void initialiseLN(const  IMCPathGenerator*  pathGen) const {
        // empty
    }

    /** equivalent to InstIntoMCProduct. */
    InsuranceAProdSV(const InsuranceA* inst,
                            const SimSeriesSP& simSeries) ;

    /** compute fees (supposed to be paid at the begining of the period) */
    /** return this period fees */
    virtual double computeFees(double* fees, const int& timeFromStart, const int& iWDGroup);

    /** update benefits */
    virtual void updatebenefit(DoubleArraySP withdrawals, const int& iWDGroup);

    /** update step up */
    virtual void updateStepUp(const int& timeFromStart, const int& iWDGroup);

    /** compute coupons taken by clients */
    virtual void computeWithdrawals(const int& timeFromStart,
                            DoubleArraySP withdrawals,
                            const int& iWDGroup);

    /** update account value */
    virtual void updateAccountValue(DoubleArraySP withdrawals, bool alreadySettled, const int& iWDGroup);

    /** update number of clients by taking into account the lapse rate */
    virtual void updateLapseClients(const int& timeFromStart, DoubleArraySP withdrawals, const int& iWDGroup);

    /** update number of clients by taking into account the mortality rate */
    virtual void updateDeadClients(const int& timeFromStart, const int& iWDGroup);

    virtual void payoff(const IPathGenerator*  pathGen, IMCPrices&  prices);

    /** print extra output **/
    virtual void recordExtraOutput(Control* control, Results* results, const IMCPrices&) const;
    
    /** vol interp for LN */
    CVolRequestLNArray getVolInterp(const IMCPathGenerator* pathGen, int iAsset) const;

public:
    /** override the initial function. All discounting is done in payoff */
    //to review ?????do we need it?
//    virtual double pvFromPaymentDate() const {
//        return 1.0;
//    }


        /** Appends 'true' (ie non derived) state variable generators
        required to the supplied collector.*/
    virtual void collectStateVars(IStateVariableCollectorSP svCollector) const{
        // ask for a reference level State Variable
        svCollector->append(refLevelGen.get());
        svCollector->append(spotGen.get());
        svCollector->append(matDfGen.get());
        svCollector->append(dfGen.get());
    }


protected:
    /** Override default method on IMCProduct. This method is called every time
        the path generator is changed (which is, at the moment, when the
        past path generator is created, and then when the future path
        generator is created  */
    virtual void pathGenUpdated(IStateVariableGen::IStateGen* newPathGen);


protected:

    double computeBasket(const int iStep) ;

    const InsuranceA* inst;               // reference to original instrument

    RandUniformDefaultSP       uniRand;          // generator of uniform random variables
    DoubleArray                yearFrac;         // year Frac for each period
                                                 // time in years to previous date

    double                     pathValue;        // store the path value
    double                     feeValue;         // store the present value of the fees
    vector<DoubleArray>        accountValue;     // current account value;
    vector<DoubleArray>        groupSize;        // current percentage of people in the 
                                                 // corresponding age group
    vector<DoubleArray>         benefit;     // garanteed minimum withdrawal benefit for the different age groups

    /* for preservation of the past */
    double                     pathValueSoFar;   // store the path value
    double                     feeValueSoFar;    // store the present value of the fees
    vector<DoubleArray>        accountValueSoFar;// current account value;
    vector<DoubleArray>        groupSizeSoFar;   // current percentage of people in the 
                                                 // corresponding age group
    vector<DoubleArray>         benefitSoFar; // garanteed minimum withdrawal benefit for the different age groups


    SVGenSpot::IStateVarSP    spotSV;        // asset state variable
    SVGenSpotSP               spotGen;       // generator for asset


    //SVGenSpotSP                  spotGen;      //!< Generator for spot
    IRefLevel::IStateVarGenSP    refLevelGen;  //!< Generator for ref level
    //SVGenSpot::IStateVarSP       spotSV;       //!< Spot state variable
    IRefLevel::IStateVarSP       refLevelSV;   //!< Ref level state variable
    
    SVGenDiscFactorSP            matDfGen;        //!< Generator for discount factors at maturity
    SVDiscFactorSP               matDfSV;         //!< Df state variable
    SVGenDiscFactorSP            dfGen;           //!< Generator for discount factors from sim dates
    SVDiscFactorSP               dfSV;            //!< Df state variable

};
// end of class InsuranceAProd


DRLIB_END_NAMESPACE
