//----------------------------------------------------------------------------
//
//   Group       : QR
//
//   Filename    : IndexBasisCalcAddin.cpp
//
//   Description : Addin to compute the index basis.
//                 The reason why this addin lives here rather than in 
//                 market/CreditIndex.cpp (where it would naturally belong) 
//                 is that this addin requires ClosedFormCDSPS which also lives 
//                 in the product folder and cannot be moved because it 
//                 requires Quanto stuff which lives in mcarlo mcarlo header 
//                 files are not included when making the market folder to 
//                 avoid dependency loops - quite messy, really 
//
//   Author      : Jose Hilera
//
//   Date        : November 2005
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/Addin.hpp"
#include "edginc/CreditIndex.hpp"
#include "edginc/CreditIndexBasis.hpp"
#include "edginc/CDSIndexParSpreads.hpp"
#include "edginc/ClosedFormCDSBasket.hpp"

DRLIB_BEGIN_NAMESPACE


/******************************************************************************
 *                       INDEX_BASIS calculation add-in                       
 ******************************************************************************/

class IndexBasisCalcAddin: public CObject, 
                           public virtual ClientRunnable 
{
public:
    static CClassConstSP const TYPE;

    // Addin parameters
    CreditIndexWrapper index;  // The CreditIndex to calculate the basis for
    MarketDataSP       market;
    IModelSP           model;

    IndexBasisCalcAddin() : CObject(TYPE)
    {}

    IObjectSP run() {
        static const string method("IndexBasisCalcAddin::run");
        try {
            IModelSP mdl = model;
            if (!mdl.get()) {
                // Create a "default" model to get data out of the market
                mdl = IModelSP(new ClosedFormCDSBasket());
            }

            // Get the index data
            index.getData(mdl.get(), market);
            // and return the index basis (const_cast, because run can not 
            // return smartConstPtr
            return CreditIndexBasisSP(
                 const_cast<CreditIndexBasis*>(index->getIndexBasis().get()));
        }
        catch (exception& e) {
            throw ModelException(e, method);
        }            
    }
    
    
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(IndexBasisCalcAddin, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(ClientRunnable);
        EMPTY_SHELL_METHOD(defaultIndexBasisCalcAddin);
        FIELD(market, "Market data cache");
        FIELD(index, "Index wrapper (name)");
        FIELD(model,"(optional) model");
        FIELD_MAKE_OPTIONAL(model);
        
        Addin::registerObjectMethod("INDEX_BASIS",
                                    Addin::RISK,
                                    "Returns the index basis",
                                    true,
                                    Addin::returnHandle,
                                    &IndexBasisCalcAddin::run);
    }
    
    static IObject* defaultIndexBasisCalcAddin(){
        return new IndexBasisCalcAddin();
    }
};

CClassConstSP const IndexBasisCalcAddin::TYPE = CClass::registerClassLoadMethod(
    "IndexBasisCalcAddin", typeid(IndexBasisCalcAddin), load);



/******************************************************************************
 *                    INTERPOLATED_INDEX_CDS_CURVE add-in                       
 ******************************************************************************/

class InterpolatedIndexCDSCurve: public CObject, 
                                 public virtual ClientRunnable 
{
public:
    static CClassConstSP const TYPE;

    // Addin parameters
    CreditIndexWrapper index;  // The CreditIndex to calculate the basis for
    MarketDataSP       market;
    IModelSP           model;

    InterpolatedIndexCDSCurve() : CObject(TYPE)
    {}

    IObjectSP run() {
        static const string method("InterpolatedIndexCDSCurve::run");
        try {
            IModelSP mdl = model;
            if (!mdl.get()) {
                // Create a "default" model to get data out of the market
                mdl = IModelSP(new ClosedFormCDSBasket());
            }

            // Get the index data
            index.getData(mdl.get(), market);
            // and return the CDSIndexParSpreads (const_cast, because run can not 
            // return smartConstPtr
            return CDSIndexParSpreadsSP(
                 const_cast<CDSIndexParSpreads*>(index->getIndexCDSCurve(true).get()));
        }
        catch (exception& e) {
            throw ModelException(e, method);
        }            
    }
    
    
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(InterpolatedIndexCDSCurve, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(ClientRunnable);
        EMPTY_SHELL_METHOD(defaultInterpolatedIndexCDSCurve);
        FIELD(market, "Market data cache");
        FIELD(index, "Index wrapper (name)");
        FIELD(model,"(optional) model");
        FIELD_MAKE_OPTIONAL(model);
        
        Addin::registerObjectMethod("INTERPOLATED_INDEX_CDS_CURVE",
                                    Addin::RISK,
                                    "Returns the CDSIndexParSpreadsSP "
                                    "corresponding to the index, interpolated "
                                    "at the names' maturities if required",
                                    true,
                                    Addin::returnHandle,
                                    &InterpolatedIndexCDSCurve::run);
    }
    
    static IObject* defaultInterpolatedIndexCDSCurve() {
        return new InterpolatedIndexCDSCurve();
    }
};


CClassConstSP const InterpolatedIndexCDSCurve::TYPE = 
    CClass::registerClassLoadMethod("InterpolatedIndexCDSCurve", 
                                    typeid(InterpolatedIndexCDSCurve), 
                                    load);


bool IndexBasisCalcAddinLoad() {
    return IndexBasisCalcAddin::TYPE != NULL;
}

DRLIB_END_NAMESPACE
