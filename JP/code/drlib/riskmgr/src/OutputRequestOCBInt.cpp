//
//

#include "edginc/config.hpp"
#include "edginc/OutputRequestOCBInt.hpp"

DRLIB_BEGIN_NAMESPACE

class OutputRequestOCBIntHelper{
public:

    /** Invoked when class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic();
        REGISTER(OutputRequestOCBInt, clazz);
        SUPERCLASS(OutputRequest);
        EMPTY_SHELL_METHOD(defaultOutputRequestOCBInt);
        FIELD(mtmPrice, "Mark to market price");
        FIELD(quoteFace, "Face value mtmPrice is quoted for. Same face used for output.");
        FIELD(isClean, "Is the mtmPrice quoted clean?");
    }

    static IObject* defaultOutputRequestOCBInt(){
        return new OutputRequestOCBInt();
    }
};

//// creates deep copy - need to override default to copy data pointer over
IObject* OutputRequestOCBInt::clone() const{
    OutputRequestOCBInt* copy = new OutputRequestOCBInt();
    copy->mtmPrice = mtmPrice;
    copy->quoteFace = quoteFace;
    copy->isClean = isClean;
    return copy;
}

OutputRequestOCBInt::OutputRequestOCBInt():OutputRequest(OutputRequestOCBInt::TYPE, OutputRequest::OPTION_ON_CONVERTIBLE_INTRINSIC_VALUE),
    hasFinished(false)  {} 

CClassConstSP const OutputRequestOCBInt::TYPE = CClass::registerClassLoadMethod(
    "OutputRequestOCBInt", typeid(OutputRequestOCBInt), OutputRequestOCBIntHelper::load);



DRLIB_END_NAMESPACE
