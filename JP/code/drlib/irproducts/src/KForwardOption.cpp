//----------------------------------------------------------------------------
//
//   Group       : Global QR&D
//
//   Filename    : KForwardOption.cpp
//
//   Description : option component
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/KForwardOption.hpp"
#include "edginc/TreeSliceOperExt.hpp"
#include "edginc/RatesUtils.hpp"
#include "edginc/Addin.hpp"
#include "edginc/Results.hpp"
#include "edginc/Black.hpp"

DRLIB_BEGIN_NAMESPACE


class KForwardOptionTree : public FDProduct {

public:

    /************************ variables ************************/
    const KForwardOptionConstSP  inst;
    FDProductSP     undProd;
    int             optType;    // copied from inst->optType for convenience

    /************************ methods ************************/
    KForwardOptionTree(const KForwardOption* inst, FDModel* model);

    virtual void init(Control*) const;
    virtual void initProd();
    virtual void setCalcOff(void);
    virtual void update(int& step, UpdateType type);
    virtual void recordOutput(Control* ctrl, YieldCurveConstSP, Results* results);
    virtual DateTime getStartDate(void) const {return model->getValueDate();}
    // only forward induction is implemented yet

protected:
    vector< TreeSliceSP > slices; // containing the info on options
    TreeSliceSP domZero, fwdTree; // needed for implied vol calculation

    TreeSliceSP volSlice, indicator;  // helper slice for impliedVol calculation

    vector< CDoubleMatrixSP > impliedVols; // obsDates x (assetVals x strikes )

};

/******************************** KForwardOptionTree ********************************/

void KForwardOptionTree::recordOutput(Control* ctrl, YieldCurveConstSP, Results* results) {
    string ccy = inst->discount->getCcy();

    try {
        int nObs     = inst->obsDates.size() ;
        int nStrikes = inst->strikeValues.size() ;

        inst->recordExtraOutput(ctrl, results);

        // save the option grid in the DEBUG_PACKET (for convenience)
        for (int i = 0 ; i < nObs; ++i)        
        {
            string obsName = "VOLS_" + Format::toString(i);

            // store all the implied vols for each observation date
            results->storeGreek( impliedVols[i],
                Results::DEBUG_PACKET, 
                OutputNameConstSP(new OutputName(obsName)));

        }
    } 
    catch (exception& e) {
        throw ModelException(e, "KForwardOptionTree::recordOutput, failed for outputName " +
                             getOutputName());
    }
}


KForwardOptionTree::KForwardOptionTree(const KForwardOption* inst, FDModel* model) :
    FDProduct(model), inst(inst), optType(inst->optType)
{
    undProd = model->createProduct(inst->und); // create the product
}

void KForwardOptionTree::init(Control*) const{
    // transform the single date into array to be able to call addCritDates
    DateTimeArray exDates; 
    exDates.push_back(inst->exerciseDate);
    model->addCritDates( exDates );

    // register the different observation dates
    model->addCritDates( inst->obsDates );
}




void KForwardOptionTree::setCalcOff(void) { 
        FDProduct::setCalcOff();
        undProd->setCalcOff(); 
}

/** initialising and setting product variables */
void KForwardOptionTree::initProd(void){
    // one slice for each exercise. Used to carry "underlying - strike" from exercise to notification
    int i , nbObs  = inst->obsDates.size() ;
    int nbFwds    = inst->forwardObs.size() ;
    int nbStrikes = inst->strikeValues.size() ;
    const string curve = inst->discount->getName();

    slices.resize(nbStrikes+1);
    for (i = 0; i< nbStrikes; ++i) 
    {
        slices[i] = model->createSlice( curve ); // default all slices to no DEV
        slices[i]->name = inst->outputName + "_strikePrice" + Format::toString( inst->strikeValues[i] );
        *(slices[i])  = 0.;
    }

    // create slice for the zero Bond and the underlying 
    fwdTree = model->createSlice( curve ); 
    domZero = model->createSlice( curve ); 
    volSlice  = model->createSlice( curve ); 
    indicator = model->createSlice( curve ); 

    // create results matrix 
    for (i = 0 ; i < nbObs; ++i)
    {
        impliedVols.push_back( CDoubleMatrixSP::attachToRef( 
            new CDoubleMatrix( nbStrikes, nbFwds ) ) );
    }
}


struct bsCallVol : SliceOper< double >
{
    static Type apply( double ind, double fwd, double price, double strike, double tau) {

        double iVol = 0.0;
        if (ind > 0.0)
        {
            double iVar = 0.01*tau; // initial guess
            Black::impliedVariance( true, fwd, strike , 1.0 , 
                  iVar , price, 1e-6, iVar); 
            iVol = sqrt(iVar/tau);
        }
        return (iVol);
    }
    static const char* symbol(){return "bsCallVol";}
};

template< typename A, typename B, typename C >
inline
Slice5Expr<bsCallVol, A, B, C, double, double> impliedVol( const SliceMarker< A > & ind, 
              const SliceMarker< B > & fwd, const SliceMarker< C > & price, 
              double strike, double tau ) 
{
    return Slice5Expr<bsCallVol, A, B, C, double, double>(static_cast< const A & >(ind),
        static_cast< const B & >(fwd),static_cast< const C & >(price), strike, tau);
}


/** update products slice at each step */
//
// -- main calculation part of the product --
//
void KForwardOptionTree::update(int& step, UpdateType type) {
    try {
        int i, j, nbStrikes = inst->strikeValues.size();
        int obsId, nbObs = inst->obsDates.size();
        DateTime stepDate = model->getDate(step); // the current date

         if (inst->exerciseDate == stepDate) // update the option settlements
         {
             const TreeSlice &underlying = undProd->getValue(step, stepDate);

             for (i = 0 ; i < nbStrikes; ++i)
             {
                switch (optType) {
                    case KForwardOption::OptType::CALL: 
                        *(slices[i]) = smax(underlying - inst->strikeValues[i],0.); 
                        break;
                    case KForwardOption::OptType::PUT:  
                        // not really supported (bsvol not working yet for put - laziness)
                        *(slices[i]) = smax(inst->strikeValues[i] - underlying,0.); 
                        break;
                    default: 
                        ModelException("optType not supported");
                }
                startDEV(slices[i]);
             } // for i

             // initialise the zero Coupon and the underlying at exercise
             *fwdTree = underlying;
             *domZero = 1.0;

             startDEV(fwdTree);
             startDEV(domZero);
         } // end of setting the option slices

         // at observation date: calculate the implied vols 
         // surface at each of the prespecified vol points
         for (obsId = 0; obsId<nbObs; ++obsId) {
             if ((inst->obsDates)[obsId] != stepDate) {
                 continue; // no observation date -> test next one
             }

             // compute the time to exercise
             const double tExer = inst->obsDates[obsId].yearFrac( inst->exerciseDate  );
             const TreeSlice &underlying = undProd->getValue(step, stepDate);

             //underlying.PrintSlice("c:/fxslice", step);
             //for (i = 0; i < nbStrikes; ++i)
             //    slices[i]->PrintSlice("c:/slices_"+Format::toString(i), step);

             switch (optType) {
                 case KForwardOption::OptType::CALL: 
                     for (i = 0 ; i < nbStrikes; ++i)
                     {
                         // for each observation windoe around a forward value
                         // average the implied vol observations
                         for (j=0; j < inst->forwardObs.size() ; ++j)
                         {
                             koTernary( *indicator, underlying, inst->forwardObs[j]*(1 - inst->width),
                                 inst->forwardObs[j]*(1+inst->width), 0.0, 1.0, 0.0, false /*!!!*/, false);

                             // calculate the implied vol for each node in the tree 
                             *volSlice = impliedVol( *indicator, (*fwdTree)/(*domZero), *(slices[i])/(*domZero), 
                                 inst->strikeValues[i], tExer   );

                             const double av = esum( (*indicator), (*indicator) );
                             const double iVol = esum((*indicator), (*indicator) * (*volSlice) );

                             (*(impliedVols[obsId]))[i][j] = (av > 0) ? iVol/av : 0.0;

                         }
                     }                
                 case KForwardOption::OptType::PUT: 
                      //  (*optionPrice) = smax((*exerciseValue), (*optionPrice)); 
                        break;

                    default: throw ModelException("optType not supported");
             }
         }

    } 
    catch (exception& e) {
        throw ModelException(e, "KForwardOptionTree::update");
    }
}



/********************************** KForwardOption **********************************/

void KForwardOption::validatePop2Object(void)
{
    try {
        int i;
        int nbStrikes = strikeValues.size();
        int nbObs     = obsDates.size();
        int nbFwds    = forwardObs.size();

        if (nbStrikes == 0) throw ModelException("Strike Values are empty.");

        for (i = 1 ; i < nbObs ; ++i)
        {
            if (obsDates[i] <= obsDates[i-1])
                throw ModelException("observation dates need to be strictly increasing");
        }

        //const DateTimeArray& exerDates = exerciseSchedule->getDateArray();
        if (obsDates.back() > exerciseDate )
            throw ModelException( "last observation date is > exerciseDate");

        // check strikes are positive
        for (i = 0; i < nbStrikes; ++i) {
            if (strikeValues[i] < 0.0 )
                throw ModelException(
                    "strikeValues[" + Format::toString(i ) + "] (" 
                    + Format::toString( strikeValues[i] ) + ") is < 0.0 " );
        }

        // check forward observations are positive
        for (i = 0; i < nbFwds; ++i) {
            if (forwardObs[i] < 0.0 )
                throw ModelException(
                    "forwardObs[" + Format::toString(i ) + "] (" 
                    + Format::toString( forwardObs[i]) + ") is < 0.0 " );
        }
    }
    catch (exception& e) {
        throw ModelException(e, "KForwardOption::validatePop2Object");
    }
}

bool KForwardOption::isDead(DateTime valueDate, double *price) const {
    // need a real pricing
    return false;
}


/** implement FDModel product interface */
FDProductSP KForwardOption::createProduct(FDModel * model) const
{
    return FDProductSP(new KForwardOptionTree(this, model));
}

void KForwardOption::load(CClassSP& clazz)
{
    clazz->setPublic();
    clazz->setDescription("Forward Option component");
    REGISTER(KForwardOption, clazz)
    SUPERCLASS(KComponent);
    EMPTY_SHELL_METHOD(defaultConstructor);
    IMPLEMENTS(FDModel::IIntoProduct);
    FIELD(und, "underlying index for option payoff")
    FIELD(optType, "CALL or PUT")
    FIELD( obsDates ,"setting date for the forward");
    FIELD( exerciseDate, "exercise date of the underlying option" );    
    FIELD( forwardObs ,"forward observation points");
    FIELD( strikeValues ,"strike Values for underlying option");
    FIELD( width ,"relative width for observation window");
    FIELD_MAKE_OPTIONAL(width);
}

CClassConstSP const KForwardOption::TYPE = CClass::registerClassLoadMethod(
    "KForwardOption", typeid(KForwardOption), KForwardOption::load);


/**********************************/
bool KForwardOptionLoad()
{
    return KForwardOption::TYPE !=0;
}

START_PUBLIC_ENUM_DEFINITION(KForwardOption::OptType::Enum, "");
ENUM_VALUE_AND_NAME(KForwardOption::OptType::CALL, "CALL", "");
ENUM_VALUE_AND_NAME(KForwardOption::OptType::PUT, "PUT", "");
END_ENUM_DEFINITION(KForwardOption::OptType::Enum);


DRLIB_END_NAMESPACE

