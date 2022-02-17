//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : PDFRequestLNStrike.cpp
//
//   Description : 
//
//   Author      : Andrew J Swain
//
//   Date        : 15 March 2002
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/PDFRequestLNStrike.hpp"

DRLIB_BEGIN_NAMESPACE

class PDFRequestLNStrikeHelper{
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic();
        REGISTER(PDFRequestLNStrike, clazz);
        SUPERCLASS(PDFRequest);
        EMPTY_SHELL_METHOD(defaultPDFRequestLNStrike);
        FIELD(request, "request");
        FIELD(startDateTime, "startDateTime");
        FIELD(callSpreadWidth, "callSpreadWidth");
        FIELD_MAKE_OPTIONAL(callSpreadWidth);
        FIELD(accuracy, "accuracy");
        FIELD_MAKE_OPTIONAL(accuracy);
    }

    static IObject* defaultPDFRequestLNStrike(){
        return new PDFRequestLNStrike();
    }
};

CClassConstSP const PDFRequestLNStrike::TYPE = CClass::registerClassLoadMethod(
    "PDFRequestLNStrike", typeid(PDFRequestLNStrike), 
    PDFRequestLNStrikeHelper::load);

const double PDFRequestLNStrike::default_callSpreadWidth = 0.0001;
const double PDFRequestLNStrike::default_accuracy        = 1.0;

const VolRequestLNStrike* PDFRequestLNStrike::getVolRequest() const {
    return request.get();
}


PDFRequestLNStrike::~PDFRequestLNStrike(){}

PDFRequestLNStrike::PDFRequestLNStrike(VolRequestLNStrike* request,
                                       DateTime            startDateTime,
                                       double              epsilon,
                                       double              accuracy):
    PDFRequest(TYPE), request(copy(request)), startDateTime(startDateTime),
    callSpreadWidth(epsilon), accuracy(accuracy) {}

PDFRequestLNStrike::PDFRequestLNStrike(VolRequestLNStrike* request,
                                       double              epsilon,
                                       double              accuracy):
    PDFRequest(TYPE), request(copy(request)), startDateTime(request->getStartDate()),
    callSpreadWidth(epsilon), accuracy(accuracy) {}

PDFRequestLNStrike::PDFRequestLNStrike(): 
callSpreadWidth(default_callSpreadWidth), accuracy(default_accuracy), PDFRequest(TYPE){}

PDFRequestLNStrike::PDFRequestLNStrike(const CClassConstSP& clazz): 
callSpreadWidth(default_callSpreadWidth), accuracy(default_accuracy), PDFRequest(clazz){}

DateTime PDFRequestLNStrike::getStartDateTime() const {
    return startDateTime;   
}

double PDFRequestLNStrike::getCallSpreadWidth() const {
    return callSpreadWidth;
}

double PDFRequestLNStrike::getAccuracy() const {
    return accuracy;
}

DRLIB_END_NAMESPACE
