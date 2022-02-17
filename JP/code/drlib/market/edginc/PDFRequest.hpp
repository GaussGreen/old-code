//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : PDFRequest.hpp
//
//   Description : 
//
//   Author      : Andrew J Swain
//
//   Date        : 15 March 2002
//
//
//----------------------------------------------------------------------------

#ifndef PDFREQUEST_HPP
#define PDFREQUEST_HPP

#include "edginc/Object.hpp"
#include "edginc/smartPtr.hpp"
#include "edginc/DateTime.hpp"

DRLIB_BEGIN_NAMESPACE

class CVolRequest;

class MARKET_DLL PDFRequest : public CObject {
public:
    static CClassConstSP const TYPE;
    friend class PDFRequestHelper;
    virtual DateTime getStartDateTime() const = 0;

    virtual ~PDFRequest();

protected:
    PDFRequest(const CClassConstSP& clazz);

private:
    PDFRequest(const PDFRequest &rhs);
    PDFRequest& operator=(const PDFRequest &rhs);
};

typedef smartConstPtr<PDFRequest> PDFRequestConstSP;
typedef smartPtr<PDFRequest> PDFRequestSP;

DRLIB_END_NAMESPACE
#endif
