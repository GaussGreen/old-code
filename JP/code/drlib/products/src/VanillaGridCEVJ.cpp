//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : VanillaGridCEVJ.cpp
//
//   Description : Grid of European Vanillas for CEV+Jump Tree
//
//   Author      : Keiji Kitazawa
//
//   Date        : 07 Janualy 2003
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/UtilFuncs.hpp"
#include "edginc/LinearStrikeVolRequest.hpp"
#include "edginc/Tree1fLN.hpp"
#include "edginc/VanillaGridCEVJ.hpp"

DRLIB_BEGIN_NAMESPACE

// helpers
class VanillaGridCEVJHelper {
public:
    static void load(CClassSP& clazz){
        clazz->setPublic();
        REGISTER(VanillaGridCEVJ, clazz);
        SUPERCLASS(Generic1Factor);
        EMPTY_SHELL_METHOD(create);
        IMPLEMENTS(CTree1fCEVJ::IIntoProduct);
        FIELD(exerciseSchedule,        "Exercise Schedule");
        FIELD(StrikeArray,      "Strike Grid Array");
    }

    static IObject* create(){
        return new VanillaGridCEVJ();
    }
};


CClassConstSP const VanillaGridCEVJ::TYPE = CClass::registerClassLoadMethod(
    "VanillaGridCEVJ", typeid(VanillaGridCEVJ), VanillaGridCEVJHelper::load);
bool  VanillaGridCEVJLoad() {
    return (VanillaGridCEVJ::TYPE != 0);
   }


// constructor
VanillaGridCEVJ::VanillaGridCEVJ(): Generic1Factor(TYPE), isCall(false)
{
    isCall = true;
}


/** initialise tree1f - allow product customisation */
void VanillaGridCEVJ1fProd::InitTree(CControl* control)
{
    int i;

    // default to DollarDivTreament = false !!! to be removed after EAS/EIS can handle it correctly !!!
    tree1f->SetDivAmountTreatment(false);

    // set base step used only for default input (StepsPerYear ==0)
    model1F->TimePts.SetTreeStepBase(300);

    // ************ below is similar to Vanilla  ***************
    /** customize tree parameters here and set up the tree */
    DateTimeArray segDates;
    segDates.resize(2); // to try using more segments
    segDates[0] = instrCEVJ->valueDate;

    DateTime mat =instrCEVJ->exerciseSchedule->lastDate();
    vector<int> density (1);
    density[0] = 1;
    segDates[1] = mat;

    double minGap = 0; // insert all points
    bool useEqualTime = false; // false = equal variance steps
    int numOfInsertNode = (tree1f->GetSmoothMethod() == CTree1f::NODE_INSERTION ?2:0);

    // all exercise dates are copied to critical dates
    DateTimeArray critDates = instrCEVJ->exerciseSchedule->getDates();
    // remove exercise date from crit date
    critDates.erase(critDates.end()-1);
    // add div event dates if needed
    EventAssetMove divEvent;
    DateTimeArraySP divCritDates;
    //if (instrCEVJ->canExerciseEarly || tree1f->DEBUG_DollarDivMethod > 0)
    if (tree1f->DEBUG_DollarDivMethod > 0)
    {// American exercise or dollar div interp treatment







        // create div events, this call is very expensive for basket with lots of div dates
        // should only need once for pricing call but need to think how to store/copy for tweaks
        int numDivs = (int)(4*segDates[0].yearFrac(segDates[1]))+1; // 4 divs per year selected as critical dates

        if (AssetUtil::getJumpEvents(instrCEVJ->asset.get(), 
                                     segDates[0], 
                                     segDates[1],
                                     numDivs,
                                     divEvent)){
            // calculate critical dates
            divCritDates = divEvent.getCritDate(numDivs, instrCEVJ->isCall);
            for (i=0; i<divCritDates->size(); i++)
                critDates.push_back((*divCritDates)[i]);
        }
    }

    // call tree step set up routine
    tree1f->Setup(instrCEVJ->valueDate, segDates, density, &critDates, 
                  minGap, useEqualTime, NumOfPriceArray, numOfInsertNode);
}


/** initialise product specific data */
void VanillaGridCEVJ1fProd::InitProd()
{
    numKGrid = instrCEVJ->StrikeArray.size();

    Prices = DoubleArraySP(new DoubleArray(numKGrid));
    //Prices.resize(numKGrid);     //to store Results
    
    // set discrete adjustment factor. Estimate trading days per year using one year from now.
    // use 364 because inclusive counting
    DateTime endDate(instrCEVJ->valueDate.getDate()+364, instrCEVJ->valueDate.getTime()); 
    VolDaysPerYear = 252.0; // model1F->GetTimeMetric()->volDays(instrCEVJ->valueDate, endDate);

    // get strikes
    Strikes.resize(numKGrid);
    for (int i=0; i<numKGrid; i++)
    {
        Strikes[i] = instrCEVJ->StrikeArray[i];
    }

}

/** calculate barriers and place barriers at inserted node if needed */
void VanillaGridCEVJ1fProd::preCalcTree(int step, int idx)
{
}



/** product payoff method at steps earlier than maturity */
void VanillaGridCEVJ1fProd::PayoffAtMat(const double * s, int step, int bot, int top, int pStart, int pEnd,
                    double * const * price)
{
    static const string method = "VanillaGridCEVJ1fProd::PayoffAtMat";
    int i, j;

    // use 2 price arrays, [0] for KIKO and [1] for KO only !!!
    ASSERT(pStart == 0 && pEnd == numKGrid-1);


/* No need settlement adjustment
    double settlementPV = instrCEVJ->instSettle->pvAdjust(instrCEVJ->exerciseSchedule->lastDate(),
                                                         instrCEVJ->discount.get(), 
                                                         instrCEVJ->asset.get());
*/
    // loop the price for 
    for (j=-bot; j<=top; j++)
    {
        for (i=0; i<numKGrid; i++)
            price[i][j] = GetIntrinsic(s[j], Strikes[i], true, true);
            //price[i][j] = settlementPV*GetIntrinsic(s[j], Strikes[i], true, true);
    }

}

/** product payoff method at steps earlier than maturity */
void VanillaGridCEVJ1fProd::PayoffBeforeMat(const double * s, int step, int bot, int top, int pStart, int pEnd,
                    double * const * price)
{
    if (step == 0) {
        //add prices
        for (int i=0; i<numKGrid; i++)
        {
            (*Prices)[i] = price[i][0];     // Store the price
//            price[i][0] -= instrCEVJ->TargetPrices[i];  //NodePrice will be replaced by price diffs.
        }
    }   
}

CAssetConstSP VanillaGridCEVJ1fProd::GetAssetRef()
{
    return CAssetConstSP::dynamicCast(
        (IObjectConstSP)instrCEVJ->asset.getSP());
}

string VanillaGridCEVJ1fProd::getCcyTreatment()
{
    return instrCEVJ->ccyTreatment;
}


/** create a tree payoff product */
CTree1fCEVJ::IProduct* VanillaGridCEVJ::createProduct(CTree1fCEVJ* model) const
{
    VanillaGridCEVJ1fProd* treeProd = new VanillaGridCEVJ1fProd(this);
    treeProd->tree1f = model;
    return treeProd;
}

void VanillaGridCEVJ::Validate()
{
    static const string method = "VanillaGridCEVJ::Validate";
    // just check the things that aren't/cannot be checked in 
    // validatePop2Object
}

/** returns a vol request for log-normal vol */
CVolRequestConstSP VanillaGridCEVJ1fProd::GetLNRequest()
{
    // get strike and maturity date from instrument
    DateTime        matDate = instrCEVJ->exerciseSchedule->lastDate();

    double volStrike  = instrCEVJ->exerciseSchedule->lastValue();

    DateTime imntStartDate = instrCEVJ->fwdStarting? 
        instrCEVJ->startDate: instrCEVJ->valueDate;

    CVolRequestConstSP volRequest(
        new LinearStrikeVolRequest(volStrike, imntStartDate, 
                                   matDate, instrCEVJ->fwdStarting));
    return volRequest;
}

// initiate GetMarket 
void VanillaGridCEVJ::GetMarket(const IModel*         model, 
                            const CMarketDataSP    market)
{
    market->GetReferenceDate(valueDate);

    if (asset.usingCache() || !Tree1fLN::TYPE->isInstance(model))
    {// should always come in - just to cope with old regression convertion
        CAsset::getAssetMarketData(model, market.get(), ccyTreatment, 
                                   discount, asset);
    }

    discount.getData(model, market);
    instSettle->getMarket(model, market.get());
    if (premiumSettle.get())
    {
        premiumSettle->getMarket(model, market.get());
    }
}

/** returns the current value date */
DateTime VanillaGridCEVJ::getValueDate() const
{
   return valueDate;
}


/** when to stop tweaking */
DateTime VanillaGridCEVJ::endDate(const Sensitivity* sensControl) const {
    DateTime matDate = exerciseSchedule->lastDate();
    DateTime instEnd  = instSettle->settles(matDate, asset.get());
    DateTime assetEnd = asset->settleDate(matDate);
    DateTime end = instEnd.isGreater(assetEnd) ? instEnd : assetEnd;
    return end;
}

void VanillaGridCEVJ::addOutputRequests(Control* control,
                                 Results* results,
                                 const double& fairValue,
                                 const double& indVol) const
{
    // take care of additional outputs
    if ( control && control->isPricing() )
    {
    }
}

void VanillaGridCEVJ1fProd::recordOutputRequests(Control* control, Results* results, double fairValue)
{
    // take care of additional outputs
    if ( control && control->isPricing() )
    {
        results->storeGreek(Prices, Results::DEBUG_PACKET, OutputNameSP(new OutputName("Prices")));

        double indVol = 0;
        instrCEVJ->addOutputRequests(control,
                                    results,
                                    fairValue,
                                    indVol);
    }
}

DRLIB_END_NAMESPACE
