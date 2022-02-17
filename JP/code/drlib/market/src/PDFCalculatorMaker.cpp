//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : PDFCalculatorMaker.cpp
//
//   Description : generate model pdf for various processes
//
//   Author      : Ning Shen
//
//   Date        : 15 June 2004
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/PDFCalculatorMaker.hpp"

/** moved outside space name and outside class to reduce 
    the size of the name for Solaris Opt compilation 
    (bug in the compiler) */
USING_DRLIB_NAMESPACE
/** PDFCalculatorMaker::pdfCreator* has been replaced by void* to reduce the size
    of the template
    it is necessary because of a bug in the compiler for SolarisOpt */
static map<string, void*> pdfMaker;


DRLIB_BEGIN_NAMESPACE

/** Used by engine to register a method which can build a PDFCalculator */
void PDFCalculatorMaker::setCreator(
    pdfCreator* method,
    const string& type) // decides which pdfCreator to set
{
    pdfMaker[type] = (void*)method;
}

/** Used by a particular vol to build a PDFCalculator */
PDFCalculator* PDFCalculatorMaker::makePDFCalculator(
    const string& type, // decides which pdfCreator to use
    const DateTime& valueDate,
    const IModel* model,
    const Asset* asset, 
    const CVolProcessed* vol)
{
    pdfCreator* pdf = (pdfCreator*)pdfMaker[type];
    if (!pdf)
        throw ModelException("PDFCalculatorMaker::makePDFCalculator",
                             "Don't know how to build PDF Calc for " + type);
    return pdf(valueDate, model, asset, vol);
}

DRLIB_END_NAMESPACE
