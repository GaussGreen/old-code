//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : MatrixShift.cpp
//
//   Description : Vega Matrix tweaking
//
//   Author      : Andre Segger
//
//   Date        : 29 March 2001
//
//
//----------------------------------------------------------------------------


#include "edginc/config.hpp"
#include "edginc/MatrixResult.hpp"
#include "edginc/TweakGroup.hpp"
#include "edginc/LastSensDate.hpp"
#include "edginc/VegaMatrix.hpp"
#include "edginc/SensMgr.hpp"
#include "edginc/Results.hpp"
#include "edginc/Untweakable.hpp"
#include <algorithm>

DRLIB_BEGIN_NAMESPACE

const double MatrixShift::BETA = 1.002;

MatrixShift::~MatrixShift(){}

void MatrixShift::calculate(TweakGroup*      tweakGroup,
                            CResults*        results) 
{

    int         numTweaks  = 0;
    int         idx;
    double      firstDeriv;
    bool        emptyResult = true;

    try {
        // get list of names to calculate result for
        OutputNameArrayConstSP names(this->names(tweakGroup, results));

        if (names->empty()) {
            results->storeNotApplicable(this);            
        }

        // see if the instrument has a last sens date method
        LastSensDate* lsd = dynamic_cast<LastSensDate*>(tweakGroup->
                                                        getInstrument());
        DateTime       endDate;
        DateTime       valueDate;

        if (lsd) {
            valueDate = tweakGroup->getInstrument()->getValueDate();
            endDate   = lsd->endDate(this);
        }


        for (idx = 0; idx < names->size(); idx++) {
            // store the name of what we want to shift
            setMarketDataName((*names)[idx]);
            // skip over where result has been calculated already 
            try {
                if (!results->exists(this) )
                {
                    CDoubleArraySP revisedXAxis(new DoubleArray(0));
                    
                    // now get expiries for which we should tweak
                    expiries = getExpiries(tweakGroup);

                    // get the x axis to shift:
                    // dummy try/catch to stop nt.vc6.opt wetting its pants

					try {
                        xAxis = getXAxis(tweakGroup->getInstrument(), 
                                         (*names)[idx]);
					}
					catch (...) { throw; }
                    
                    if ( expiries->size() > 0 && xAxis->size() > 0 ) {
                        /* create matrix to hold the results */
                        CDoubleMatrixSP matrixOut(
                            new CDoubleMatrix(xAxis->size(),expiries->size()));
                        xLowerIdx = 0;
                        xUpperIdx = 0;
                        numTweaks = 0;
                        /* If we have strikes extremely close together
                           these strikes will be tweaked at the same time
                           and reported as a single sensitvity to the
                           lowest strike in the series. */

                        // sort the strikes into ascending order
                        sort(xAxis->begin(), xAxis->end());

                        while (xUpperIdx < xAxis->size() )
                        {
                            bool  pastMaturity = false; /* default */
                            numTweaks++;
                
                            xLowerIdx = xUpperIdx;

                            /* Find the upper x idx to tweak. This will do
                               nothing if adjacent strikes are
                               sufficiently far apart */
                            while( (xUpperIdx+1 < xAxis->size() ) &&
                                   (((*xAxis)[xUpperIdx+1]/(*xAxis)[xUpperIdx])
                                    <= BETA)) {
                                xUpperIdx++;
                                matrixOut->removeLastCol(); // trim matrix
                            }

                            revisedXAxis->push_back((*xAxis)[xLowerIdx]);

                            // calculate base price
                            ResultsSP baseResults = getBasePrice(tweakGroup,
                                                                 (*names)[idx]);
                        
                            /* iterate through the points on the yAxis
                               to tweak */
                            for (expiryIdx = 0; expiryIdx < expiries->size();
                                 expiryIdx++) {
                                /* only calculate if instrument is
                                   still alive */
                                if (!pastMaturity) {               
                                    // calculate sens
                                    firstDeriv = 
                                        calcOneSidedFirstDeriv(tweakGroup,
                                                               baseResults.get());
                                } else  {
                                    firstDeriv = 0.0;
                                }

                                /* store result */
                                (*matrixOut)[numTweaks-1][expiryIdx] = 
                                    firstDeriv;
                       
                                // do we need to tweak anymore ?
                                if (lsd) {
                                    pastMaturity = (*expiries)[expiryIdx]->
                                        toDate(valueDate).isGreater(endDate);
                                }                                            
                            }
                            /* This will also have the effect of
                               incrementing xLowerIdx at the begginning of
                               the iteration of the outer loop */
                            xUpperIdx++;
                        }


                        MatrixResultSP matrixResult(
                            new MatrixResult(revisedXAxis, 
                                             expiries,
                                             matrixOut));
                        // and store it
                        results->storeGreek(matrixResult, this);

                        // we have at least one result in the result
                        // set, so no need to store NotApplicable any
                        // more
                        emptyResult = false;
                    }
                }
            }
            catch (exception& e){
                results->storeGreek(IObjectSP(new Untweakable(e)), this);
                emptyResult = false;
            }
        }

        // if non of the underlyings had any sensitive strikes, we
        // report NoApplicable
        if (emptyResult && !names->empty()) {
            results->storeNotApplicable(this);            
        }

    } catch (exception& e){
        results->storeGreek(IObjectSP(new Untweakable(e)), getSensOutputName(),
                            OutputNameSP(new OutputName("")));
    }

}

/** Returns the expiries which are to be tweaked */
ExpiryArrayConstSP MatrixShift::getExpiries(const IObject* tweakGroup){
    IObjectConstSP expiriesObj = SensMgrConst(tweakGroup).qualifier(this);
    ExpiryArrayConstSP expiries = ExpiryArrayConstSP::dynamicCast(expiriesObj);
    return expiries;
}

DoubleArrayConstSP MatrixShift::getXAxisValues() const{
    if (!xAxis){
       throw ModelException("xAxis is Null",  "MatrixShift::getXAxisValues");
    }
    return xAxis;
}

/** Returns the expiries which are to be tweaked */
ExpiryArrayConstSP MatrixShift::getExpiries() const{
    if ( !expiries ) {
       throw ModelException("expiry is Null",  "MatrixShift::getExpiries");
    }
    return expiries;
}

int MatrixShift::getLowerIdx() const{
    return xLowerIdx;
}

int MatrixShift::getUpperIdx() const{
    return xUpperIdx;
}

int MatrixShift::getExpiryIdx() const{
    return expiryIdx;
}

ResultsSP MatrixShift::getBasePrice(TweakGroup*       tweakGroup,
                                    OutputNameConstSP name)
{
    // calculate base price
    MatrixShiftSP zeroControl(getZeroShift());

    zeroControl->setMarketDataName(name);
    ResultsSP baseResults = ResultsSP( new Results());
    zeroControl->xLowerIdx = 0;
    zeroControl->xUpperIdx = 0;
    zeroControl->expiryIdx = 0;
    zeroControl->expiries  = expiries;
    zeroControl->xAxis     = xAxis;
    double basePrice = zeroControl->shiftAndPrice(tweakGroup, 0.0, false);
    baseResults->storePrice(basePrice, "DUNNO");

    return baseResults;
}


MatrixShift::MatrixShift(const CClassConstSP& clazz,
                         const string&        outputName,
                         const double&        shiftSize): 
    ScalarShift(clazz, outputName, shiftSize), 
    xLowerIdx(0), xUpperIdx(0), expiryIdx(0){}

/** for reflection */
MatrixShift::MatrixShift(const CClassConstSP& clazz,
                         const string&        sensName):
    ScalarShift(clazz, sensName),
    xLowerIdx(0), xUpperIdx(0), expiryIdx(0){}


class MatrixShiftHelper{
public:
    /** Invoked when class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(MatrixShift, clazz);
        SUPERCLASS(ScalarShift);
        FIELD(xLowerIdx, "lower index to tweak");
        FIELD_MAKE_TRANSIENT(xLowerIdx);
        FIELD(xUpperIdx, "upper index to tweak");
        FIELD_MAKE_TRANSIENT(xUpperIdx);
        FIELD(expiryIdx, "expiry index to tweak");
        FIELD_MAKE_TRANSIENT(expiryIdx);
        FIELD(expiries, "the y-coordinates");
        FIELD_MAKE_TRANSIENT(expiries);
        FIELD(xAxis, "the x-coordinates");
        FIELD_MAKE_TRANSIENT(xAxis);
    }
};

CClassConstSP const MatrixShift::TYPE = CClass::registerClassLoadMethod(
    "MatrixShift", typeid(MatrixShift), MatrixShiftHelper::load);

DRLIB_END_NAMESPACE


