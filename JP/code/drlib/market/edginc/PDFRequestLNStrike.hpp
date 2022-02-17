//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : PDFRequestLNStrike.hpp
//
//   Description : 
//
//   Author      : Andrew J Swain
//
//   Date        : 15 March 2002
//
//
//----------------------------------------------------------------------------

#ifndef PDFREQUESTLN_HPP
#define PDFREQUESTLN_HPP

#include "edginc/PDFRequest.hpp"
#include "edginc/VolRequestLNStrike.hpp"

DRLIB_BEGIN_NAMESPACE

class MARKET_DLL PDFRequestLNStrike : public PDFRequest {
public:
    static CClassConstSP const TYPE;

    /** Returns the vol request that this PDFRequest was constructed with */
    virtual const VolRequestLNStrike* getVolRequest() const;

    virtual ~PDFRequestLNStrike();

    PDFRequestLNStrike(
        VolRequestLNStrike* request,
        DateTime            startDateTime,
        double              callSpreadWidth = default_callSpreadWidth,
        double              accuracy        = default_accuracy);

    PDFRequestLNStrike(
        VolRequestLNStrike* request,
        double              callSpreadWidth = default_callSpreadWidth,
        double              accuracy        = default_accuracy);
    

    virtual DateTime getStartDateTime()   const;
    virtual double   getCallSpreadWidth() const;
    virtual double   getAccuracy()        const;

    static const double default_callSpreadWidth;
    static const double default_accuracy;


protected:
    PDFRequestLNStrike(const CClassConstSP& clazz);
private:
    friend class PDFRequestLNStrikeHelper;
    PDFRequestLNStrike();
    PDFRequestLNStrike(const PDFRequestLNStrike &rhs);
    PDFRequestLNStrike& operator=(const PDFRequestLNStrike &rhs);

    VolRequestLNStrikeSP request;
    DateTime startDateTime;         // will be used to deduce whether CS is fwdStarting
    double   callSpreadWidth;       
    double   accuracy;              // specified as a multiple of callSpreadWidth i.e. 
                                    // accuracy = 0.001 means 0.001 * epsilon
};

typedef smartConstPtr<PDFRequestLNStrike> PDFRequestLNStrikeConstSP;
typedef smartPtr<PDFRequestLNStrike> PDFRequestLNStrikeSP;

DRLIB_END_NAMESPACE
#endif
