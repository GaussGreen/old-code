//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : PDFLogNormal.cpp
//
//   Description : generate model pdf for log-normal
//
//   Author      : Ning Shen
//
//   Date        : 15 June 2004
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/PDFLogNormal.hpp"
#include "edginc/PDFCalculatorMaker.hpp"

DRLIB_BEGIN_NAMESPACE

// pdfCreator
static PDFCalculator* pdfCreator(const DateTime& valueDate, 
                                   const IModel* model,
                                   const Asset* asset, 
                                   const CVolProcessed* vol)
{
    if (!asset || !vol)
        throw ModelException("pdfCreator", "failed to create pdf: no asset or vol inputs");

    return new PDFLogNormal(valueDate,
                             asset,
                             vol);
}

CClassConstSP const PDFLogNormal::TYPE = 
    CClass::registerClassLoadMethod("PDFLogNormal", typeid(PDFLogNormal), load);

	/** Calculate the cumulative probability at each strike.
        Each prob is estimated as a "small" call spread (i.e., by tweaking the strike)
        Returns false if if the difference between 2 consecutive probs is negative.
    */
void PDFLogNormal::probabilities(const DoubleArray& strikes,
                               const DateTime&    maturity,
                               DoubleArray&       probs) const
{
	int nbStrikes = strikes.size();

	double fwd = asset->fwdValue(maturity);
    double v_t = sqrt(volBS->CalcVar(valueDate, maturity));

    double d2;

	int iStrike;
	for (iStrike=0; iStrike < nbStrikes; ++iStrike){
        d2 = log(fwd/strikes[iStrike])/v_t - 0.5*v_t;
		probs[iStrike]=N1(d2);
	}
}

/* external symbol to allow class to be forced to be linked in */
bool PDFLogNormalLinkIn(){
    PDFCalculatorMaker::setCreator(pdfCreator, "PDFLogNormal");
    return (PDFLogNormal::TYPE != 0);
}


DRLIB_END_NAMESPACE
