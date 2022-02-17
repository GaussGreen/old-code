//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : PDFCalculatorMaker.hpp
//
//   Description : create model pdf for various processes
//
//   Author      : Ning Shen
//
//   Date        : 15 June 2004
//
//
//----------------------------------------------------------------------------

#ifndef PDFCalculatorMaker_HPP
#define PDFCalculatorMaker_HPP

#include "edginc/DateTime.hpp"
#include "edginc/Asset.hpp"
#include "edginc/VolProcessed.hpp"
#include "edginc/PDFCalculator.hpp"

DRLIB_BEGIN_NAMESPACE

/** create model pdf  */
class MARKET_DLL PDFCalculatorMaker{
public:
    typedef PDFCalculator* (pdfCreator)(const DateTime& valueDate, 
                                       const IModel* model,
                                       const Asset* asset, 
                                       const CVolProcessed* vol);

    /** Used by engine to register a method which can build a PDFCalculator */
    static void setCreator(pdfCreator* method,
                           const string& type);

    /** Used by a particular vol to build a PDFCalculator */
    static PDFCalculator* makePDFCalculator(const string& type, // decides which pdfCreator to use
                                            const DateTime& valueDate,
                                            const IModel* model,
                                            const Asset* asset, 
                                            const CVolProcessed* vol);
private:
};

DRLIB_END_NAMESPACE
#endif
