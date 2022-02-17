//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : FastQuote.cpp
//
//   Description : model tuned for real time vanilla options pricing.
//
//   Author      : Ning Shen
//
//   Date        : 24 June 2002
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/FastQuote.hpp"

DRLIB_BEGIN_NAMESPACE

///////////////////////////////////////////////////////
//************* FastQuote class ***********************/
///////////////////////////////////////////////////////
// EdrAction version of addin
IObjectSP FastQuote::run()
{
    return (IObjectSP) ComputeQuote();
}

/** for regression run */
IObjectSP FastQuote::runTest() const
{
    return (IObjectSP) ComputeQuote();
}

/** calculate options */
CResultsArraySP FastQuote::ComputeQuote() const
{
    static const string routine = "FastQuote::ComputeQuote";
    try {
        clock_t startTime = clock();

        int numOpt = input->size();

        // initialisation
        quoteEnv->init(spot);

        // results
        CResultsArraySP output(new CResultsArray(numOpt));

        // loop for all options
        for (int i=0; i<numOpt; i++)
        {
            (*output)[i] = CResultsSP(new CResults);
            try{
                quoteEnv->Price(*(*input)[i], *(*output)[i]);
            }
            catch(ModelException& e){
                e.errorLog(); // write log and continue
            }
        }

        clock_t endTime = clock();
        double calcTime = (double)(endTime-startTime)/CLOCKS_PER_SEC;
        // store total time in the last result's debug section
        (*output)[numOpt-1]->storeScalarGreek(calcTime, Results::DEBUG_PACKET,
                                  OutputNameSP(new OutputName("TOTAL_TIME")));

        return output;
    } 
    catch (exception& e){
        throw ModelException(e, routine);
    }
}

IObjectSP FastQuote::ComputeQuoteAddin(FastQuote* obj)
{
    return (IObjectSP) obj->ComputeQuote();
}

/** Invoked when Class is 'loaded' */
void FastQuote::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(FastQuote, clazz);
    SUPERCLASS(CObject);
    IMPLEMENTS(ClientRunnable);
    EMPTY_SHELL_METHOD(defaultFastQuote);
    FIELD(spot, "spot price");
    FIELD(input, "input data structure");
    FIELD(quoteEnv, "pricing environment");

    Addin::registerClassObjectMethod("FAST_QUOTE",
                                     Addin::RISK,
                                     "compute array of options",
                                     TYPE,
                                     false,
                                     Addin::expandSimple,
                                     (Addin::ObjMethod*)ComputeQuoteAddin);

}

CClassConstSP const FastQuote::TYPE= CClass::registerClassLoadMethod(
    "FastQuote", typeid(FastQuote), load);

bool FastQuoteLoad(){
    bool canLoad = !!FastQuote::TYPE &&
                   !!FastQuoteEnv::TYPE &&
                   !!FastQuoteInput::TYPE &&
                   !!FastQuoteInputArray::TYPE;
    return canLoad;
}

DRLIB_END_NAMESPACE

