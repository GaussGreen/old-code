//----------------------------------------------------------------------------
//
//   Group       : Global QR&D
//
//   Filename    : IndexSpecEQ.hpp
//
//   Description : specification for equity underlying/payoff index
//
//----------------------------------------------------------------------------

#ifndef INDEX_SPEC_EQ_HPP
#define INDEX_SPEC_EQ_HPP

#include "edginc/IndexSpec.hpp"
#include "edginc/Asset.hpp"
#include "edginc/BorrowCurve.hpp"
#include "edginc/InstrumentSettlement.hpp"
#include "edginc/AssetHistoryContainer.hpp"


DRLIB_BEGIN_NAMESPACE

class MARKET_DLL IndexSpecEQ : 
    virtual public IIndexSpec,
    public Asset // asslow legacy way of supplying asset
{
public:
    static CClassConstSP const TYPE;

    /*************************** methods ***************************/
    virtual void addResetDates(const DateTimeArray &resetDates);
    /** returns index name */
    virtual string getName() const { return name; }

    virtual void setup(const IModel* model, const MarketData* market);

    virtual IMarketFactorConstSP getFactor() const { return asset.getSP(); }
    virtual IMarketFactorSP getFactor() { return asset.getSP(); }

    // only used by IndexSpecEDR for bridging legacy prod
    IndexSpecEQ( const string & name, const CAssetWrapper & asset, const string& ccyTreatment) :
                Asset(TYPE), name(name), asset(asset), ccyMode(ccyTreatment){} // remove ccyTreatment

    // *** Asset:: methods ***
    virtual string getYCName() const{
        return yc->getName();
    }

    virtual double getSpot() const {return indexSpot;}

    virtual double fwdValue(const DateTime& date) const{
        // to do!!!
        throw ModelException("IndexSpecEQ::fwdValue", "not implemented");
        return 0.0;
    }
    virtual void fwdValue(const DateTimeArray&     dateList,
                          const FwdValueAlgorithm& algo,
                          CDoubleArray&            result) const{
        Asset::fwdValue(dateList, result);
    }

    /** Returns an processed vol*/
    virtual IVolProcessed* getProcessedVol(const CVolRequest* volRequest) const;

    /** Calculate the settlement date associated with a given trade date */
    virtual DateTime settleDate(const DateTime& tradeDate) const;

    virtual PDFCalculator* pdfCalculator(const PDFRequest *) const{
        throw ModelException("IndexSpecEQ::pdfCalculator", "not implemented");
    }
    //**** end of Asset::

    virtual void getEvents(
            const FixingReqEvent*, IModel* model,
            const DateTime& eDate, EventResults* events) const;

protected:
    IndexSpecEQ(const CClassConstSP & type = TYPE);
    virtual AssetHistoryConstSP getAssetHistory(const string &source) const;

	/************************ exported fields ************************/
    string                  name;
    double                  indexSpot;            
    CAssetWrapper           asset;

    CVolBaseWrapper         vol;
    BorrowCurveSP           borrowCurve;
    InstrumentSettlementSP  settle;
    /****************** transiend fields ************/
    AssetHistoryContainerSP history;
    DateTimeArray resetDates;

    // to be initialised by child classes when needed
    string                  ccyMode;
    YieldCurveWrapper       yc;

private:
    static void load(CClassSP& clazz);
    static IObject* defaultConstructor() { return new IndexSpecEQ(); }
};

typedef smartPtr<IndexSpecEQ> IndexSpecEQSP;
typedef smartConstPtr<IndexSpecEQ> IndexSpecEQConstSP;
typedef array<IndexSpecEQSP, IndexSpecEQ> IndexSpecEQArray;

DRLIB_END_NAMESPACE

#endif
