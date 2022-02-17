//----------------------------------------------------------------------------
//
//   Group       : Global QR&D
//
//   Filename    : KRib.hpp
//
//   Description : Rib component
//
//----------------------------------------------------------------------------
#ifndef _KRib_HPP
#define _KRib_HPP

#include "edginc/KComponent.hpp"
#include "edginc/FlexDates.hpp"

DRLIB_BEGIN_NAMESPACE

class IRPRODUCTS_DLL KRib : public KComponent,
             virtual public FDModel::IIntoProduct
             //virtual public LastSensDate
{
public:
    static CClassConstSP const TYPE;    

	void setAcsDates(const FlexDates & dates) {accStartDates = dates;}
	void setAceDates(const FlexDates & dates) {accEndDates   = dates;}
//    void setCouponSched(CouponSchedSP sched) { this->sched = sched; }

    /* IProdCreator:: */
    virtual double getValue(DateTime date, CashflowInfo &cfi ) const;
	virtual void setup(const IModel* model, const MarketData* market);

private:
    friend class KRibTree;
    friend class KRibFloatLegTree;  // ??? just for now - look at simply making variables public instead

    /****************** exported fields ************/
    /* Observation and coupon underlyings */
    CModel::IProdCreatorSP  obsUnd;
    CModel::IProdCreatorSP  cpnUnd;

	DoubleArraySP spreads;
 
    // accEndDates and ccStartDates could be SPs to save memory - on the to-do list
    FlexDates    obsDates;
    FlexDates    obsEffDates;
//    CouponSchedSP sched;
    FlexDates    accEndDates;
    FlexDates    accStartDates; 
    bool         resetInArrears;
    
    /****************** methods ************/
    static void load(CClassSP& clazz);
    static IObject* defaultConstructor(void) { return new KRib(); }
    KRib(void) : KComponent(TYPE) {}

    /* FDModel::IIntoProduct:: */
    virtual FDProductSP createProduct(FDModel * model) const;
};

typedef smartPtr<KRib> KRibSP;

DRLIB_END_NAMESPACE

#endif


