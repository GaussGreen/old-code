//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : VanillaMoment.cpp
//
//   Description : Calculates moments (ATM vols, skew, convexity) of Vanilla Options
//
//   Author      : Regis Guichard
//
//   Date        : 27 March 03
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/VanillaMoment.hpp"
#include "edginc/VanillaGrid.hpp"
#include "edginc/Enum2StringListHelper.hpp"
#include "edginc/Weightmatrix.hpp"
#include "edginc/BenchmarkDate.hpp"

//#include "edginc/XMLOutputStream.hpp"

DRLIB_BEGIN_NAMESPACE
    
typedef VanillaMoment::RequestType VMReqType;
template<> string nameForType<VMReqType>(VMReqType*){
    return "VanillaMoment::RequestType";
}

typedef Enum2StringListHelper<VMReqType> VMReqTypeHelper;
template<> string VMReqTypeHelper::names[VMReqTypeHelper::EnumList::NB_ENUMS] = {
    "ATM_VOL",
    "SKEW",
    "CONVEXITY"
};

static double calcPercStrike(int iPercStrike, int nbPercStrikes){
    if (nbPercStrikes == 1){
        return 1.0;
    }
    return (0.9 + iPercStrike * 0.1);
}

static void calcMoments(const DoubleArray& vols,
                        double&            atmVol,
                        double&            skew,
                        double&            convexity){
    int nbPercSrikes = vols.size();
    if (nbPercSrikes == 1){
        atmVol = vols[0];
        return;
    }
    atmVol = vols[1];
    skew = vols[1] - vols[0];
    if (nbPercSrikes >= 3){
        convexity = vols[2] - 2.0 * vols[1] + vols[0];
    }
}

typedef pair<int, double> IntDoublePair;
class IntDoubleSortHelper{
public:
    bool operator()(const IntDoublePair& x, const IntDoublePair& y){
        return (x.second < y.second);
    }
};

DateTime VanillaMoment::getValueDate() const{
    return valueDate;
}

void VanillaMoment::validatePop2Object(){
    static const string method("VanillaMoment::ValidatePop2Object");
    try{
        // check request types
        int nbReqs = requests.size();
        if (requests.size() == 0){
            throw ModelException(method,
                                 "At least one request must be provided");
        }
        int iReq = 0;
        for (; iReq < nbReqs; ++iReq){
            // will throw exception if request type is not recognized
            VMReqTypeHelper::getIndex(requests[iReq]);
        }
    }
    catch(exception& e){
        throw ModelException(e, method);
    }
}

void VanillaMoment::GetMarket(const IModel*        model, 
                              const CMarketDataSP  market){
    static const string method("VanillaMoment::GetMarket");
    try{
        market->GetReferenceDate(valueDate);
        const IModel* modelToUse = model;
        // special treatment if CompositeModel
        // note: if not CompositeModel, will eventually fail as only CompositeModel is 
        // supported, but there is no reason to throw an exception at this very point
        if (CompositeModel::TYPE->isInstance(model)){
            const CompositeModel* mdl = dynamic_cast<const CompositeModel*>(model);
            mdl->checkNbModels(1, "VanillaMoment");
            modelToUse = mdl->getModel(0).get();
        }
        CAsset::getAssetMarketData(modelToUse, 
                                   market.get(), 
                                   ccyTreatment, 
                                   discount, 
                                   asset);
        discount.getData(modelToUse, market);
        instSettle->getMarket(modelToUse, market.get());
        if (premiumSettle.get()){
            premiumSettle->getMarket(modelToUse, market.get());
        }
    }
    catch(exception& e){
        throw ModelException(e, method);
    }
}
    
class VanillaMomentProduct: public CompositeModel::IProduct{
public:
    virtual void price(CompositeModel* model,
                       Control*        control, 
                       CResults*       results) const{
        static const string routine = "VanillaMomentProduct::price";
        try {
            // No price is returned. Results are treated as output requests
            if (!control || !control->isPricing()){
                return;
            }

            model->checkNbModels(1, "VanillaMoment");

            // control
            OutputRequestArraySP outputArray(new OutputRequestArray(1));
            (*outputArray)[0] = OutputRequestSP(new OutputRequest("IND_VOL"));
            Control ctrl(SensitivityArrayConstSP(new SensitivityArray(0)),
                         outputArray,
                         false,
                         string(""));

            // price
            CResults res;
            ctrl.calculate(model->getModel(0).get(),
                           vanillaGrid.get(),
                           &res);

            IObjectConstSP obj(res.retrieveRequestResult("IND_VOL"));
            VanillaGrid::OutputConstSP vols(VanillaGrid::OutputConstSP::dynamicCast(obj));

            // reorder vols
            int nbSteps = inst->maturities.size();
            DoubleArrayArray percVols(nbSteps);
            int iStep = 0;
            for (; iStep < nbSteps; ++iStep){
                percVols[iStep].resize(nbPercStrikes);
            }
            if (!inst->isAtmSpot){
                int iStrike = 0;
                for (; iStrike < nbStrikes; ++iStrike){
                    int iStep = idxTable[iStrike].first / nbPercStrikes;
                    int iPercStrike = idxTable[iStrike].first % nbPercStrikes;
                    percVols[iStep][iPercStrike] = vols->getValue(iStep, iStrike);
                }
            }
            else{
                int iStep = 0;
                for (; iStep < nbSteps; ++iStep){                
                    int iPercStrike = 0;
                    for (; iPercStrike < nbPercStrikes; ++iPercStrike){
                        percVols[iStep][iPercStrike] = vols->getValue(iStep, iPercStrike);
                    }
                }
            }

            // compute moments
            DoubleArray atmVols(nbSteps);
            DoubleArray skews(nbSteps);
            DoubleArray convexities(nbSteps);
            for (iStep = 0; iStep < nbSteps; ++iStep){                
                calcMoments(percVols[iStep],
                            atmVols[iStep],
                            skews[iStep],
                            convexities[iStep]);
            }

            // return results
            results->storeGreek(IObjectSP(atmVols.clone()),
                                "Moments",
                                OutputNameSP(new OutputName("atmVols")));
            if (nbPercStrikes == 1){
                return;
            }
            results->storeGreek(IObjectSP(skews.clone()),
                                "Moments",
                                OutputNameSP(new OutputName("skews")));
            if (nbPercStrikes == 2){
                return;
            }
            results->storeGreek(IObjectSP(convexities.clone()),
                                "Moments",
                                OutputNameSP(new OutputName("convexities")));
        } 
        catch (exception& e){
            throw ModelException(e, routine);
        }
    }

    VanillaMomentProduct(const VanillaMoment*         inst,
                         const VanillaGridSP&         vanillaGrid,
                         int                          nbPercStrikes,
                         int                          nbStrikes,
                         const vector<IntDoublePair>& idxTable):
    inst(inst),
    vanillaGrid(vanillaGrid),
    nbStrikes(nbStrikes),
    nbPercStrikes(nbPercStrikes),
    idxTable(idxTable){}

private:
    const VanillaMoment*    inst;
    VanillaGridSP           vanillaGrid;
    int                     nbStrikes;
    int                     nbPercStrikes;
    vector<IntDoublePair>   idxTable;
};

CompositeModel::IProduct* VanillaMoment::createProduct(CompositeModel* model) const{
    static const string method("VanillaMoment::createProduct");
    try{
        int iReq = 0;
        int nbPercStrikes = 1;
        for (; iReq < requests.size(); ++iReq){
            // count # of strikes needed
            nbPercStrikes = Maths::max(VMReqTypeHelper::getIndex(requests[iReq]) + 1, nbPercStrikes);
        }

        int nbSteps = maturities.size();
        int nbStrikes;
        DoubleArray strikes;
        CDoubleMatrixSP weights;
        vector<IntDoublePair> idxTable;
        // if atm spot, strikes are 90%, 100% and 110%
        // otherwise, strikes are different for every maturity
        if (!isAtmSpot){
            DoubleArray refs(nbSteps);
            asset->fwdValue(maturities,
                            refs);
            if (fwdStarting){
                double fwd = asset->fwdValue(startDate);
                for (int iStep = 0; iStep < nbSteps; ++iStep){
                    refs[iStep] /= fwd;
                }
            }
            nbStrikes = nbPercStrikes * nbSteps;
            idxTable.resize(nbStrikes);
            int iStep = 0;
            int iStrike = 0;
            for (; iStep < nbSteps; ++iStep){
                int iPercStrike = 0;
                for (; iPercStrike < nbPercStrikes; ++iPercStrike, ++iStrike){
                    double strike = calcPercStrike(iPercStrike, nbPercStrikes) * refs[iStep];
                    idxTable[iStrike].first = iStrike;
                    idxTable[iStrike].second = strike;
                }
            }
            sort(idxTable.begin(), idxTable.end(), IntDoubleSortHelper());
            strikes.resize(nbStrikes);
            weights = CDoubleMatrixSP(new CDoubleMatrix(nbSteps, nbStrikes));
            for (iStep = 0; iStep < nbSteps; ++iStep){
                for (iStrike = 0; iStrike < nbStrikes; ++iStrike){
                    (*weights)[iStep][iStrike] = 0.0;
                }
            }
            for (iStrike = 0; iStrike < nbStrikes; ++iStrike){
                strikes[iStrike] = idxTable[iStrike].second;
                iStep = idxTable[iStrike].first / nbPercStrikes;
                (*weights)[iStep][iStrike] = 1.0;
            }
        }
        else{
            double ref = (fwdStarting ? 1.0 : asset->getSpot());
            nbStrikes = nbPercStrikes;
            strikes.resize(nbPercStrikes);
            weights = CDoubleMatrixSP(new CDoubleMatrix(nbSteps, nbPercStrikes));
            int iPercStrike = 0;
            for (; iPercStrike < nbPercStrikes; ++iPercStrike){
                strikes[iPercStrike] = calcPercStrike(iPercStrike, nbPercStrikes) * ref;
            }
            int iStep = 0;
            for (; iStep < nbSteps; ++iStep){                
                for (iPercStrike = 0; iPercStrike < nbPercStrikes; ++iPercStrike){
                    (*weights)[iStep][iPercStrike] = 1.0;
                }
            }
        }
        
        // create a WeightMatrix for VanillaGrid
        DateTimeArraySP maturitiesSP(copy(&maturities));
        DoubleArraySP strikesSP(copy(&strikes));
        string instType = isCall ? WeightMatrix::CALL : WeightMatrix::PUT;
        WeightMatrixSP weightMatrix(new WeightMatrix("", 
                                                     maturitiesSP,
                                                     strikesSP,
                                                     weights,
                                                     WeightMatrix::ABSOLUTE,
                                                     instType,
                                                     startDate));
        WeightMatrixWrapper weightMatrixWrapper(weightMatrix);
        
        VanillaGridSP vanillaGrid(new VanillaGrid(fwdStarting,
                                                  startDate,
                                                  oneContract,
                                                  notional,
                                                  initialSpot,
                                                  valueDate,
                                                  asset,
                                                  ccyTreatment,
                                                  discount,
                                                  weightMatrixWrapper,
                                                  instSettle,
                                                  instSettle));
        vanillaGrid->Validate();

//        XMLOutputStream xml("c:/temp/xxxdebug.xml");
//        vanillaGrid->xmlWrite("vanillaGrid", xml);

        return new VanillaMomentProduct(this,
                                        vanillaGrid,
                                        nbPercStrikes,
                                        nbStrikes,
                                        idxTable);
    }
    catch(exception& e){
        throw ModelException(e, method);
    }
}


/** Returns the name of the instrument's discount currency */
string VanillaMoment::discountYieldCurveName() const {
    return discount.getName();
}



/** for reflection */
VanillaMoment::VanillaMoment():
CInstrument(TYPE),
fwdStarting(false),
isCall(true),
oneContract(true),
notional(0.0),
initialSpot(0.0),
isAtmSpot(true){}

class VanillaMomentHelper {
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(VanillaMoment, clazz);
        SUPERCLASS(CInstrument);
        IMPLEMENTS(CompositeModel::IIntoProduct);
        EMPTY_SHELL_METHOD(defaultVanillaMoment);

        FIELD(valueDate, "");
        FIELD_MAKE_OPTIONAL(valueDate);
        FIELD(asset, "Asset");
        FIELD(discount, "Discount curve");
        FIELD(ccyTreatment, "Currency treatment");
        FIELD(instSettle, "Instrument settlement at maturity");
        FIELD(premiumSettle, "Settlement of option premium");
        FIELD_MAKE_OPTIONAL(premiumSettle);

        FIELD(isCall, "Call if true; put, otherwise");
        FIELD(fwdStarting, "Fwd starting if true; Started, otherwise");
        FIELD(startDate, "Start Date");
        FIELD_MAKE_OPTIONAL(startDate);
        FIELD(oneContract, "oneContract");
        FIELD(notional, "notional");
        FIELD_MAKE_OPTIONAL(notional);
        FIELD(initialSpot, "initialSpot");
        FIELD_MAKE_OPTIONAL(initialSpot);
        FIELD(maturities, "Array of maturities");
        FIELD(isAtmSpot, "isAtmSpot");
        FIELD(requests, "Requests (" + VMReqTypeHelper::getNameList() + ")");
    }

    static IObject* defaultVanillaMoment(){
        return new VanillaMoment();
    }
};

CClassConstSP const VanillaMoment::TYPE = CClass::registerClassLoadMethod(
    "VanillaMoment", typeid(VanillaMoment), VanillaMomentHelper::load);
bool  VanillaMomentLoad() {
    return (VanillaMoment::TYPE != 0);
   }


DRLIB_END_NAMESPACE

