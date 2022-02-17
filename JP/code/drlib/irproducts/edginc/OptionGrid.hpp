//----------------------------------------------------------------------------
//
//   Group       : Global QR&D
//
//   Filename    : OptionGrid.hpp
//
//   Description : option grid
//
//----------------------------------------------------------------------------
#ifndef _OPTIONGRID_HPP
#define _OPTIONGRID_HPP

#include "edginc/FDModel.hpp"
#include "edginc/Schedule.hpp"
#include "edginc/KComponent.hpp"
#include "edginc/FlexDates.hpp"
#include "edginc/OptionSched.hpp"
#include "edginc/InstrumentSettlement.hpp"
#include "edginc/IndexSpec.hpp"
#include "edginc/DoubleMatrix.hpp"


DRLIB_BEGIN_NAMESPACE

// 
class OptionGrid : public KComponent,
                    // public CInstrument,
                  virtual public FDModel::IIntoProduct,
                  virtual public Callability::IEventHandler
                  // virtual public MCModel::IIntoProdct,
               //   virtual public CClosedFormLN::IIntoProduct
                  // virtual public IProdCreator, // do we need grid of grid ??? 
{
public:
    static CClassConstSP const TYPE;
    friend class OptionGridTree; // necessary for the optType (so far)??
//    friend class OptionGridClosedForm; // necessary for the optType (so far)??

    struct Type {
        enum Enum {CALL, PUT, DIGICALL, DIGIPUT};
    };
    typedef BoxedEnum<Type::Enum> TypeBoxedEnum;

    /****************** exported fields ************/
public: // can be modified by shell instruments

    // AK: I think this should be potentially a full matrix for exp - tenor
    IProdCreatorSP  und;            // underlying index
    virtual void setup(const IModel* model, const MarketData* market) ;


    /** Returns the name of the instrument's discount currency. */
  //  virtual string discountYieldCurveName() const { return discount.getName(); } ;

   // virtual DateTime getValueDate() const { return DateTime(); };
   // virtual void Validate() { } ;

protected:
    Type::Enum       optionType;     // call, put, digicall, digiput
    OptionSchedDatesSP sched; // similar to exercises
    DoubleMatrix     strikes;     // should be a "cube" of strikes (should be a matrix!)
    InstrumentSettlementSP  instSettle;       // instrument settlement (physical, cash etc..)

    /****************** methods ************/
public:
    // explicit constructor
    OptionGrid(
        const string &discount, 
        const string &outputName, 
        IProdCreatorSP und,  // not sure if we would not prefer the IProdCreator here
        const Type::Enum& optionType, 
        InstrumentSettlementSP instSettle, 
        OptionSchedDatesSP  sched,
        DoubleMatrix  &strikes)
    :   KComponent(discount, outputName, TYPE), 
        und(und), 
        optionType(optionType), 
        sched(sched), 
        strikes(strikes),
        instSettle(instSettle)
    { validatePop2Object(); }

    virtual void validatePop2Object(void) ;

public:

    /* Callability::IEventHandler */
    void getEvents( const Callability*, IModel* model, 
        const DateTime&  eventDate, EventResults* events) const;


protected:
    /* FDModel::IIntoProduct:: */
    virtual FDProductSP createProduct(FDModel * model) const;
//    virtual MCProductSP createProduct(MCModel * model) const;
 //   virtual CClosedFormLN::IProduct* createProduct(CClosedFormLN * model) const;

    /* DeadPricer:: */
    virtual bool isDead(DateTime valueDate, double *price) const;

    OptionGrid(CClassConstSP const &type) 
        : KComponent(type) {}

private:

    static void load(CClassSP& clazz);
    static IObject* defaultConstructor(void) { return new OptionGrid(TYPE); }
};

typedef smartPtr<OptionGrid> OptionGridSP;
typedef smartConstPtr<OptionGrid> OptionGridConstSP;

DRLIB_END_NAMESPACE

#endif


