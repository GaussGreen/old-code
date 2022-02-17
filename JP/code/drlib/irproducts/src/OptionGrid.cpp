//----------------------------------------------------------------------------
//
//   Group       : Global QR&D
//
//   Filename    : OptionGrid.cpp
//
//   Description : option grid pricing 
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
//#include "edginc/Format.hpp"
#include "edginc/OptionGrid.hpp"
//#include "edginc/TreeSliceOper.hpp"
#include "edginc/TreeSliceOperExt.hpp"
#include "edginc/RatesUtils.hpp"
#include "edginc/Addin.hpp"
#include "edginc/PhysicalSettlement.hpp"
#include "edginc/EventResults.hpp"
#include "edginc/Results.hpp"
#include "edginc/Black.hpp"


DRLIB_BEGIN_NAMESPACE

//static const int optionWindow = 7; // in days
class PrintSlice {
public:
    static const int sliceCount = 1;
    FILE *f;
    TreeSlice &s;
    void compute() const {
        fprintf(f, "%15.13f\n", *s.iter);
    }
    template< typename S >
    const S** listInputSlices(const S** list) const
    {
        *list++ = dynamic_cast<const S*>(&s);
        return list;
    }
    template< typename S >
    S** listOutputSlices(S** list) const { return list; }
    void printDebug(char *s) const { strcat(s, "PrintSlice"); }
    PrintSlice(const string &filename, TreeSlice &s) : s(s) {
        f = fopen(filename.c_str(), "a");
        QLIB_VERIFY(f!=0, "Could not create file "+filename);
    }
    virtual ~PrintSlice() {
        if (f!=0) fclose(f);
    }
};


class OptionGridTree : public FDProduct {
public:
    /************************ variables ************************/
    OptionGridConstSP  inst;
    FDProductSP     undProd;
    int             optType;    // copied from inst->optType for convenience

    /************************ methods ************************/
    OptionGridTree(const OptionGridConstSP &inst, FDModel* model);
    DoubleArray getStrike(const DateTime &date); // get list of strikes

    virtual void init(Control*) const;
    virtual void initProd();
    virtual void update(int& step, UpdateType type);
    virtual void recordOutput(Control* ctrl, YieldCurveConstSP, Results* results);
    virtual DateTime getStartDate(void) const { return model->getDate(0);}
    // only forward induction is implemented yet
    virtual FDDirection getDirectionPref() const { return FDProduct::FWD_INDUCTION; }

    // return the corresponding dates equivalent to notification
    virtual DateTimeArray getStatePriceDates( ) const { return inst->sched->notifDate; };

protected:

    TreeSliceSP optionPrice;     // option slice

    CDoubleMatrixSP gridPrices, impliedVols;
    CDoubleArraySP  domDiscountTree, fwdTree;

    double totalSum; // sum of all the option prices
};

/******************************* OptionGridTree ********************************/

OptionGridTree::OptionGridTree(const OptionGridConstSP &inst, FDModel* model) :
    FDProduct(model), inst(inst), 
    optType(inst->optionType)
{
    undProd = model->createProduct(inst->und); // create the product

    int nbExer = inst->sched->exerciseDate.size();
    for (int i=0; i<nbExer; ++i) 
    {
        DateTime exerDate = inst->sched->exerciseDate[i];
        undProd->addModelResetDate(exerDate);
    }
}

void OptionGridTree::recordOutput(Control* ctrl, YieldCurveConstSP, Results* results) {
    string ccy = inst->discount->getCcy();

    try {

        string ccy = inst->discount->getCcy();
        //const TreeSlice &price = getValue(0, model->getDate(0));

        // store the total sum of all the options in the results
        // This does not seem the ideal solution (how to register something as price?)
        if (inst->isTopComponent()) 
        {
            if (!inst->outputName.empty()) {
                OutputRequest* request = ctrl->requestsOutput(inst->outputName);
                if (request) {
                    results->storeRequestResult(request, totalSum, inst->outputName); 
                }
                else if (ctrl->requestsOutput(OutputRequest::DBG)) {
                    // store it as a DEBUG_PACKET
                    results->storeScalarGreek(totalSum, Results::DEBUG_PACKET, 
                        OutputNameSP(new OutputName(inst->outputName)));
                } // else ignored
            }    
        }


        // always record option price
        //recordSliceToOutputName(ctrl, results, model, false, 
        //    "OPTION_PRICE", ccy, *optionPrice );
        inst->recordExtraOutput(ctrl, results);

        // save the option grid in the DEBUG_PACKET (for convenience)
        results->storeGreek( gridPrices,
                        Results::DEBUG_PACKET, 
                        OutputNameConstSP(new OutputName("GRID_PRICES")));

        results->storeGreek( impliedVols,
                        Results::DEBUG_PACKET, 
                        OutputNameConstSP(new OutputName("IMPLIED_VOLS")));

        // additional information from the tree
        results->storeGreek( domDiscountTree,
                        Results::DEBUG_PACKET, 
                        OutputNameConstSP(new OutputName("DOMESTIC_DISCOUNT_TREE")));

        results->storeGreek( fwdTree,
                        Results::DEBUG_PACKET, 
                        OutputNameConstSP(new OutputName("FWD_TREE")));
    } 
    catch (exception& e) {
        throw ModelException(e, "OptionGridTree::recordOutput, failed for outputName " +
                             getOutputName());
    }

}

// this function allows for strike interpolation for american
// exercise, when implemented
DoubleArray OptionGridTree::getStrike(const DateTime &date) {

    int i , j ;

    for (i = 0; i < inst->sched->exerciseDate.size(); i++) {
        if (inst->sched->exerciseDate[i] == date)
            break;
    }
    DoubleArray exStrikes;
    // ATTENTION: matrix class is implemented the wrong way around ( i.e. [cols,rows] )
    for (j = 0 ; j < inst->strikes.numCols(); ++j)
        exStrikes.push_back( (inst->strikes[j])[i] );
 
    return exStrikes;
}

void OptionGridTree::init(Control*) const{

    // insert critical dates
    model->addCritDates(inst->getValueDate().getFutureDates(inst->sched->notifDate));
    model->addCritDates(inst->getValueDate().getFutureDates(inst->sched->exerciseDate));
}

/** initialising and setting product variables */
void OptionGridTree::initProd(void){
    const string curve = inst->discount->getName();

    optionPrice = model->createSlice(curve) ;
    optionPrice->name = inst->outputName + "_optionPrice" ;
    *optionPrice  = 0.;
    totalSum      = 0.;
    
    // create results matrix 
    gridPrices  = CDoubleMatrixSP::attachToRef( new DoubleMatrix( inst->strikes.numCols(), inst->strikes.numRows() ) );
    impliedVols = CDoubleMatrixSP::attachToRef( new DoubleMatrix( inst->strikes.numCols(), inst->strikes.numRows() ) );
    domDiscountTree  = CDoubleArraySP::attachToRef( new DoubleArray( inst->strikes.numRows() ) );
    fwdTree     = CDoubleArraySP::attachToRef( new DoubleArray( inst->strikes.numRows() ) );

}

/** update products slice at each step */
void OptionGridTree::update(int& step, UpdateType type) {
    try {

        // this part can be used to generalise between 
        // forward and backward induction:
        // if type is somehow backward, then create a new slice
        // at each exercise possibility, otherwise, simply price
        // option (Call/Put/Digicall/Digiput)
        
        // this is kind of the "calc" function

        DateTime stepDate = model->getDate(step); // the current date
        const DateTimeArray &exerDates = inst->sched->exerciseDate.getDates();
        const DateTimeArray &notifDates = inst->sched->notifDate.getDates();
        int i, exerId, nbExer = exerDates.size();

        // on notification date, update option price
        for (exerId=0; exerId<nbExer; ++exerId) {
            if ((inst->sched->notifDate)[exerId]!=stepDate) {
                continue;
                // ignore or include notif date depending on setting
                if (notifDates[exerId] < inst->getValueDate())
                    continue;
            }

            // get the "fx" values at the exercise point
            const TreeSlice &underlying = undProd->getValue(step, stepDate);
            DoubleArray exerStrikes = getStrike(stepDate);
            DoubleArray resVec;

            for (i = 0 ; i < exerStrikes.size(); ++i)
            {
                bool isCall = true;
                switch (optType) {
                    case OptionGrid::Type::CALL: 
                         (*optionPrice) = smax( underlying - exerStrikes[i], 0. ); 
                        break;
                    case OptionGrid::Type::PUT: 
                         (*optionPrice) = smax( exerStrikes[i] - underlying , 0. ); 
                         isCall = false;
                        break;
                    case OptionGrid::Type::DIGICALL: // not really implemented yet
                        throw ModelException("DIGICALL not implemented yet");
                         //(*optionPrice) = smax( underlying - exerStrikes[i], 0. ); 
                         break;
                    case OptionGrid::Type::DIGIPUT:
                        throw ModelException("DIGIPUT not implemented yet");
                        // (*optionPrice) = smax( exerStrikes[i] - underlying , 0. ); 
                        // (*optionPrice) = digi(*exerciseValue); // currently not implemented
                        break;
                    default: throw ModelException("optType not supported");
                }
                // multiply DEV slices
                const TreeSlice & statePriceSlice = model->getStatePriceSlice( step ); // or similar

                //if (i == 0 ) {
                //    std::ostringstream fnameSP, fnameUnd;                   
                //    fnameSP << "C:/statePrices_" << Format::toString( step ) << ".dat";
                //    fnameUnd << "C:/underlying_" << Format::toString( step ) << ".dat";
                //    statePriceSlice.printSlice( fnameSP.str() );
                //    underlying.printSlice( fnameUnd.str() );
               // }

                // store information in the results matrix 
                (*gridPrices)[i][exerId] = esum( statePriceSlice, statePriceSlice * (*optionPrice));

                // add the current output to the total sum of prices
                totalSum += (*gridPrices)[i][exerId];

                if (i==0) // for each exercise date register some additional informations
                {
                    (*domDiscountTree)[exerId]     = esum( statePriceSlice, statePriceSlice );
                    (*fwdTree)[exerId]  = esum( statePriceSlice, statePriceSlice * underlying ) 
                                     /(*domDiscountTree)[exerId];
                }

                // calculate the implied vol for those strikes
                double iVar = 0.0;

                // get the Value Date and compute time to mat.
                double tau =  model->getValueDate().yearFrac( stepDate );
                switch (optType) {
                    case OptionGrid::Type::CALL: 
                    case OptionGrid::Type::PUT: 
                        Black::impliedVariance(isCall, (*fwdTree)[exerId], exerStrikes[i], (*domDiscountTree)[exerId] , 
                                    0.1*0.1*tau, (*gridPrices)[i][exerId], 1e-6, iVar);
                    break;
                    default: // do nothing
                        break; 
                }
                (*impliedVols)[i][exerId] = sqrt(iVar/tau);
                
            } // for i
        }


    } 
    catch (exception& e) {
        string errMsg = "OptionGridTree::update at date " + model->getDate(step).toString();
        throw ModelException(e, errMsg);
    }
}



/********************************** OptionGridClosedForm ****************************/

//class OptionGridClosedForm : public CClosedFormLN::IProduct {
//public:
//    /************************ variables ************************/
//    OptionGridConstSP    inst;
//    //CClosedFormLN   undProd;
//    int             optType;    // copied from inst->optType for convenience
//    //DateTime        skipExerDate;
//
//    /************************ methods ************************/
//    OptionGridClosedForm(const OptionGridConstSP &inst);
//    DoubleArray getStrike(const DateTime &date); // get list of strikes
//
//    virtual void init(Control*) const;
//    virtual void initProd();
//    virtual void recordOutput(Control* ctrl, YieldCurveConstSP, Results* results);
//    //virtual DateTime getStartDate(void) const { return model->getDate(0);}
//    // only forward induction is implemented yet
//
//    virtual void price(CClosedFormLN*    model,
//                       Control*          control, 
//                       CResults*         results) const { };
//
//protected:
//
//    CDoubleMatrixSP gridPrices, impliedVols;
//    CDoubleArraySP  domDiscount, fwds;
//};
//
///******************************* OptionGridTree ********************************/
//
//OptionGridClosedForm::OptionGridClosedForm(const OptionGridConstSP &inst ) :
//inst(inst), optType(inst->optionType)
//{
//    // no need to do anything as all information is in the OptionGrid
//    // undProd = model->createProduct( inst->und ); // create the product
//}
//
//void OptionGridClosedForm::recordOutput(Control* ctrl, YieldCurveConstSP, Results* results) {
//    string ccy = inst->discount->getCcy();
//
//    try {
//
//        //switch (optType) {
//        //        case OptionGrid::Type::CALL:
//        //        case OptionGrid::Type::PUT:
//        //        case OptionGrid::Type::DIGICALL: 
//        //        case OptionGrid::Type::DIGIPUT: 
//        //            recordSliceToOutputName(ctrl, results, model, true, 
//        //                "", ccy, *optionPrice ); 
//        //            break;
//                    // when option is top component, getValue(0) not called before this function
//        //}
//
//        // always record option price
//        //recordSliceToOutputName(ctrl, results, model, false, 
//        //    "OPTION_PRICE", ccy, *optionPrice );
//        inst->recordExtraOutput(ctrl, results);
//
//        // save the option grid in the DEBUG_PACKET (for convenience)
//        results->storeGreek( gridPrices,
//                        Results::DEBUG_PACKET, 
//                        OutputNameConstSP(new OutputName("GRID_PRICES")));
//
//        results->storeGreek( impliedVols,
//                        Results::DEBUG_PACKET, 
//                        OutputNameConstSP(new OutputName("IMPLIED_VOLS")));
//
//    } 
//    catch (exception& e) {
//        throw ModelException(e, "OptionGridClosedForm::recordOutput failed "); 
//    }
//
//}
//
//// this function allows for strike interpolation for american
//// exercise, when implemented
//DoubleArray OptionGridClosedForm::getStrike(const DateTime &date) {
//
//    int i , j ;
//
//    for (i = 0; i < inst->sched->exerciseDate.size(); i++) {
//        if (inst->sched->exerciseDate[i] == date)
//            break;
//    }
//    DoubleArray exStrikes;
//    // ATTENTION: matrix class is implemented the wrong way around ( i.e. [cols,rows] )
//    for (j = 0 ; j < inst->strikes.numCols(); ++j)
//        exStrikes.push_back( (inst->strikes[j])[i] );
// 
//    return exStrikes;
//}
//
//void OptionGridClosedForm::init(Control*) const{
//
//    // insert critical dates
//    //model->addCritDates(inst->getValueDate().getFutureDates(inst->sched->notifDate));
//   // model->addCritDates(inst->getValueDate().getFutureDates(inst->sched->exerciseDate));
//}
//
///** initialising and setting product variables */
//void OptionGridClosedForm::initProd(void){
//    const string curve = inst->discount->getName();
//    
//    // create results matrix 
//    gridPrices  = CDoubleMatrixSP::attachToRef( new DoubleMatrix( inst->strikes.numCols(), inst->strikes.numRows() ) );
//    impliedVols = CDoubleMatrixSP::attachToRef( new DoubleMatrix( inst->strikes.numCols(), inst->strikes.numRows() ) );
//  //  domDiscountTree  = CDoubleArraySP::attachToRef( new DoubleArray( inst->strikes.numRows() ) );
//  //  fxFwdTree   = CDoubleArraySP::attachToRef( new DoubleArray( inst->strikes.numRows() ) );
//
//}




/********************************** OptionGrid **********************************/

bool OptionGrid::isDead(DateTime valueDate, double *price) const {
    // here !isNotified
    if (valueDate > sched->notifDate.back())
    {    // cannot be excercised anymore
        if (price) 
            *price = 0;
        return true;
    }
    // need a real pricing
    return false;
}

/** implement FDModel product interface */
FDProductSP OptionGrid::createProduct(FDModel * model) const {
    return FDProductSP(new OptionGridTree(OptionGridConstSP(this), model));
}

//MCProductSP OptionGrid::createProduct(MCModel * model) const {
//    return MCProductSP(new OptionGridMC(OptionGridConstSP(this), model));
//}
//
// no shared pointer available??
//CClosedFormLN::IProduct * OptionGrid::createProduct(CClosedFormLN * model) const {
//    return new OptionGridClosedForm( OptionGridConstSP(this) );
//}


void OptionGrid::validatePop2Object(void) {
    try {
        if (sched->exerciseDate.size() != sched->notifDate.size()) 
            throw ModelException("OptionGrid::validatePop2Object", 
                    "exerciseDate.size()!=notifDates.size()");

        for (int i = 0; i < sched->exerciseDate.size(); ++i)
        {
            if ( sched->exerciseDate[i] != sched->notifDate[i] ) 
                throw ModelException("OptionGrid::validatePop2Object",
                "exercise and notif date must be the same at input no " 
                + Format::toString(i));
        }

        if (strikes.numRows() != sched->notifDate.size()) 
            throw ModelException("strikes.numRows()!=notifDates.size()");

    }
    catch (exception& e) {
        throw ModelException(e, "OptionGrid::validatePop2Object "+outputName);
    }
}


void OptionGrid::setup(const IModel* model, const MarketData* market) {
    try {
        //KComponent::setup(model, market); // from FloatLeg ??
        // sched->setup(model, market); // from FloatLeg ??
        KComponent::setup(model, market);

        und->addResetDates( sched->exerciseDate );
        und->setup(model, market);
    }
    catch (exception& e) {
        throw ModelException(e, "OptionGrid::setup, " + outputName);
    }
}





// Local function to map the OptionGrid Callability
// onto the proper Callability ones
static const string& typeEnumToCallabilityStr(int x) {
    static string const composite = "!Composite"; // maybe add it to Callability
    switch (x) {
        case OptionGrid::Type::CALL:
        case OptionGrid::Type::DIGICALL:
                        return Callability::CALL;
        case OptionGrid::Type::PUT:      
        case OptionGrid::Type::DIGIPUT:
                        return Callability::PUT;
        default: return composite;
    }
}

// Callability::IEventHandler interface
void OptionGrid::getEvents(const Callability*, 
                        IModel*            model, 
                        const DateTime&    eventDate,
                        EventResults*      events) const
{
    try {
        if (strikes.numRows() != sched->notifDate.size())
            throw ModelException("OptionGrid::getEvent: (internal error)",
                        " strikes.numRows()!=notifDates.size()");

        // returns ALL callability dates
        // no comment is made on the desirability of calling
        for (int i = 0; i < sched->notifDate.size(); i++) {
                events->addEvent( new Callability(sched->notifDate[i], 
                sched->exerciseDate[i],
                typeEnumToCallabilityStr(optionType),
                strikes[0][i]) );
        }
    }
    catch (exception& e) {
        throw makeException(e, __FUNCTION__);     
    }
}



void OptionGrid::load(CClassSP& clazz)
{
    clazz->setPublic();
    clazz->setDescription("Option component");
    REGISTER(OptionGrid, clazz)
    SUPERCLASS(KComponent);
    EMPTY_SHELL_METHOD(defaultConstructor);
    IMPLEMENTS(FDModel::IIntoProduct); // not
    FIELD(und, "underlying index for option payoff")
    FIELD_MAKE_OPTIONAL(und);  // why should this be optional??
    FIELD(optionType, "Type of option");
    FIELD(instSettle, "instrument settlement");
    FIELD_MAKE_OPTIONAL(instSettle);
    FIELD(sched, "Object containing option notification and exercise dates");
    FIELD(strikes, "Matrix of strikes values (set of strikes per date entry in options)");
    Addin::registerConstructor(Addin::UTILITIES, OptionGrid::TYPE);
}

CClassConstSP const OptionGrid::TYPE = CClass::registerClassLoadMethod(
    "OptionGrid", typeid(OptionGrid), OptionGrid::load);

START_PUBLIC_ENUM_DEFINITION(OptionGrid::Type::Enum, "");
ENUM_VALUE_AND_NAME(OptionGrid::Type::CALL, "CALL", "");
ENUM_VALUE_AND_NAME(OptionGrid::Type::PUT, "PUT", "");
ENUM_VALUE_AND_NAME(OptionGrid::Type::DIGICALL, "DIGICALL", "");
ENUM_VALUE_AND_NAME(OptionGrid::Type::DIGIPUT, "DIGIPUT", "");
END_ENUM_DEFINITION(OptionGrid::Type::Enum);


bool OptionGridLoad() {return OptionGrid::TYPE !=0;}

DRLIB_END_NAMESPACE
