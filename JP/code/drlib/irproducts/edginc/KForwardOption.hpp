//----------------------------------------------------------------------------
//
//   Group       : Global QR&D
//
//   Filename    : KForwardOption.hpp
//
//   Description : forward option komponent
//
//----------------------------------------------------------------------------
#ifndef _KFORWARDOPTION_HPP
#define _KFORWARDOPTION_HPP

//#include "edginc/Schedule.hpp"

#include "edginc/FDModel.hpp"
#include "edginc/KComponent.hpp"
#include "edginc/FlexDates.hpp"
#include "edginc/IndexSpec.hpp"


DRLIB_BEGIN_NAMESPACE

//
// experimental version of a forward volatility viewer
// not used for pricing trades but for diagnostics
//
// MCProduct and ClosedFormProduct creator left for future augmentation
//
// KForwardOption: prices forward starting options with fixed strikes (test function)
class KForwardOption : public KComponent,
                  virtual public FDModel::IIntoProduct
                  // virtual public MCModel::IIntoProdct,
                  // virtual public CClosedFormLN::IIntoProduct
                  // virtual public IProdCreator, // do we need grid of grid ??? 
{
public:
    static CClassConstSP const TYPE;
    friend class KForwardOptionTree; // necessary for the optType (so far)??
    //friend class KForwardOptionClosedForm; // necessary for the optType (so far)??

    struct OptType {
        enum Enum { CALL, PUT };
    };
    typedef BoxedEnum<OptType::Enum> OptTypeBoxedEnum;

    /****************** exported fields ************/
public: // can be modified by shell instruments

    // AK: I think this should be potentially a full matrix for exp - tenor
    IProdCreatorSP  und;            // underlying index
    virtual void setup(const IModel* model, const MarketData* market) {
        und->setup(model, market);
    }

protected:
    OptType::Enum   optType;     // call, put, digicall, digiput
    //InstrumentSettlementSP  instSettle;       // instrument settlement (physical, cash etc..)

    string          outputName;      // 
    FlexDates       obsDates;        // setting (i.e. observation) of option
    DateTime        exerciseDate;    // underlying exercise date (final date in tree)
    DoubleArray     forwardObs;      // set of points to be observed on 
    DoubleArray     strikeValues;    // fixed strike for underlying options
    double          width;            // relative step for the forward-Observations


    /****************** methods ************/
public:
    // explicit constructor
    KForwardOption(
        const string &discount, 
        const string &outputName, 
        IProdCreatorSP und,  // not sure if we would not prefer the IProdCreator here
        const OptType::Enum& optType, 
        FlexDates     &obsDates,
        DateTime      &exerciseDate,
        DoubleArray   &forwardObs,
        DoubleArray   &strikeValues,
        double         width)
    :   KComponent(discount, outputName, TYPE), 
        // CInstrument(TYPE), 
        und(und), 
        optType(optType), 
        obsDates(obsDates), 
        exerciseDate(exerciseDate),
        forwardObs(forwardObs),
        strikeValues(strikeValues),
        width(width)
    { validatePop2Object(); }

    virtual void validatePop2Object(void) ;


protected:
    /* FDModel::IIntoProduct:: */
    virtual FDProductSP createProduct(FDModel * model) const;
//    virtual MCProductSP createProduct(MCModel * model) const;
//    virtual CClosedFormLN::IProduct* createProduct(CClosedFormLN * model) const;

    /* DeadPricer:: */
    virtual bool isDead(DateTime valueDate, double *price) const;

    KForwardOption(CClassConstSP const &type) 
        : KComponent(type),width(0.1) {}

private:

    static void load(CClassSP& clazz);
    static IObject* defaultConstructor(void) { return new KForwardOption(TYPE); }
};


typedef smartPtr<KForwardOption> KForwardOptionSP;
typedef smartConstPtr<KForwardOption> KForwardOptionConstSP;


DRLIB_END_NAMESPACE

#endif


