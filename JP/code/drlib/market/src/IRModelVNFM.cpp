//----------------------------------------------------------------------------
//
//   Group       : QR Rates
//
//   Filename    : IRModelVNFM.cpp
//
//   Description : VNFM model (Vladimir's n-Factor Model).
//
//   Author      : Anwar E Sidat
//
//   Date        : 15-Aug-2006
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#define QLIB_IRModelVNFM_CPP
#include "edginc/IRModelVNFM.hpp"
#include "edginc/Expiry.hpp"


DRLIB_BEGIN_NAMESPACE

// Make sure the class links
bool IRModelVNFMLoad() { return (IRModelVNFM::TYPE != 0); }

IRModelVNFM::IRModelVNFM()
    : IRExoticParam(TYPE), Backbone(0.), i_numDates(0){}

IRModelVNFM::~IRModelVNFM() {}

IRModelVNFM::IRModelVNFM(const CClassConstSP& clazz): 
    IRExoticParam(clazz) {}

void IRModelVNFM::validatePop2Object()
{
    static const string method = "IRModelVNFM::validatePop2Object";
    try
    {
        // Check Number of Factors
        if (nbFactors < 1 || nbFactors > 3)
             throw ModelException(method, getName() + " - only 1, 2 or 3 factor model allowed!");

        // Get Number of Dates
        i_numDates = Dates.size();

        // Check matrix sizes
        if ((MeanReversion.numCols() != nbFactors)
            || (i_numDates == 0 && MeanReversion.numRows() > 1)
            || (i_numDates && MeanReversion.numRows() != i_numDates))
        {
            throw ModelException(method, getName() + " - MeanReversion matrix " +
                                 "[" + Format::toString(MeanReversion.numRows()) + "," + Format::toString(MeanReversion.numCols()) + "]" +
                                 " should have been of size " +
                                 "[" + Format::toString(i_numDates ? i_numDates : 1) + "," + Format::toString(nbFactors) + "]");
        }
        if ((Weight.numCols() != nbFactors)
            || (i_numDates == 0 && Weight.numRows() > 1)
            || (i_numDates && Weight.numRows() != i_numDates))
        {
            throw ModelException(method, getName() + " - Weight matrix " +
                                 "[" + Format::toString(Weight.numRows()) + "," + Format::toString(Weight.numCols()) + "]" +
                                 " should have been of size " +
                                 "[" + Format::toString(i_numDates ? i_numDates : 1) + "," + Format::toString(nbFactors) + "]");
        }
        int expectedCorrCols = nbFactors * (nbFactors - 1) / 2;     // number of expected correlation columns
        if ((Correlation.numCols() != expectedCorrCols)
            || (i_numDates == 0 && Correlation.numRows() > 1)
            || (nbFactors > 1 && i_numDates && Correlation.numRows() != i_numDates))
        {
            throw ModelException(method, getName() + " - Correlation matrix " +
                                 "[" + Format::toString(Correlation.numRows()) + "," + Format::toString(Correlation.numCols()) + "]" +
                                 " should have been of size " +
                                 "[" + Format::toString((nbFactors - 1) ? (i_numDates ? i_numDates : 1) : 0) + "," + Format::toString(expectedCorrCols) + "]");
        }
    }
    catch (exception& e)
    {
        throw ModelException(e, method);
    }
}

void IRModelVNFM::load(CClassSP& clazz)
{
    clazz->setPublic();
    clazz->setDescription("VNFM model (with term structure model input parameters).");
    REGISTER(IRModelVNFM, clazz);
    SUPERCLASS(IRExoticParam);
    EMPTY_SHELL_METHOD(defaultConstructor);
    FIELD(Key,           "Identifier for this object in market data cache.");
    FIELD(nbFactors,     "Enter number of factors (currently 1, 2 or 3).");
    FIELD(Dates,         "(opt) Array of Dates to represent term structure for input matrices below. Omit for non term structure case.");
    FIELD_MAKE_OPTIONAL(Dates);
    FIELD(MeanReversion, "Matrix of Mean Reversion parameters, columns as per desired number of factors (1..3) and add extra rows for term structure corresponding to Dates array.");
    FIELD(Weight,        "Matrix of factor Weights parameters, columns as per desired number of factors (1..3) and add extra rows for term structure corresponding to Dates array.");
    FIELD(Correlation,   "(opt) Matrix of Correlation values, only needed for 2 and 3 factor cases, order columns by 1&2, 1&3 and 2&3 and add extra rows for term structure corresponding to Dates array.");
    FIELD_MAKE_OPTIONAL(Correlation);
    FIELD(Backbone,      "Backbone value (defaults to 0).");
    FIELD_MAKE_OPTIONAL(Backbone);

    Addin::registerConstructor("IRModelVNFM",
                                Addin::MARKET,
                                "Creates a handle to a multi factor tree model, up to 3 factors but without term structure for input parameters.",
                                IRModelVNFM::TYPE);
}

CClassConstSP const IRModelVNFM::TYPE = CClass::registerClassLoadMethod(
    "IRModelVNFM", typeid(IRModelVNFM), IRModelVNFM::load);

/** Returns name of model */
string IRModelVNFM::getName() const
{
    return Key;
}
// definition of TYPE for MarketWrapper template class
DEFINE_TEMPLATE_TYPE(IRModelVNFMWrapper);

DRLIB_END_NAMESPACE
