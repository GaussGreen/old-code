#include "edginc/config.hpp"
#include "edginc/IConvolutor.hpp"
#include "edginc/ClientRunnable.hpp"

DRLIB_BEGIN_NAMESPACE



// addin for unit test of recursive convolutor
class CONVOLUTION_DLL ConvolutorUnitTest:
    public CObject,
    public ClientRunnable {
public:
    static CClassConstSP const TYPE;

    // addin parameters
    IDistribution1DArraySP distributions;
    IConvolutorSP convolutor;

    // EdrAction version of addin
    IObjectSP run() {
        return IObjectSP(convolutor->convolute(distributions)->clone());
    }

    ConvolutorUnitTest(): CObject(TYPE) {}
    
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(ConvolutorUnitTest, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(ClientRunnable);
        EMPTY_SHELL_METHOD(defaultConstructor);
        FIELD(distributions, "Array of discrete distributions");
        FIELD(convolutor, "Convolutor engine");
    }
    
    static IObject* defaultConstructor(){
        return new ConvolutorUnitTest();
    }
};

CClassConstSP const ConvolutorUnitTest::TYPE =
    CClass::registerClassLoadMethod(
        "ConvolutorUnitTest",
        typeid(ConvolutorUnitTest),
        load);

/* external symbol to allow class to be forced to be linked in */
bool ConvolutorUnitTestLoad(){
    return (ConvolutorUnitTest::TYPE != 0);
}

DRLIB_END_NAMESPACE
