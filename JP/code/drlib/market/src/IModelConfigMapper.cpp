#include "edginc/config.hpp"
#include "edginc/IModelConfigMapper.hpp"
#include "edginc/ClientRunnable.hpp"

DRLIB_BEGIN_NAMESPACE

static void IModelConfigMapperLoad(CClassSP& clazz) {
    REGISTER_INTERFACE(IModelConfigMapper, clazz);
    EXTENDS(IObject);
}

CClassConstSP const IModelConfigMapper::TYPE = 
    CClass::registerInterfaceLoadMethod(
        "IModelConfigMapper", 
        typeid(IModelConfigMapper), 
        IModelConfigMapperLoad);

DEFINE_TEMPLATE_TYPE(IModelConfigMapperArray);

/**
 * Addin for unit testing of IModelConfigMapper.
 * Just takes a IModelConfigMapper and a ICreditLossConfig as inputs and
 * return the mapped ICreditLossModelConfig.
 * */
class MARKET_DLL ModelConfigMapperUnitTest:
    public CObject,
    public ClientRunnable {
public:
    static CClassConstSP const TYPE;

    // Addin parameters
    IModelConfigMapperSP mapper;
    ICreditLossConfigSP lossConfig;

    // EdrAction version of addin
    IObjectSP run() {
        return IObjectSP(mapper->innerModel(lossConfig)->clone());
    }

    ModelConfigMapperUnitTest(): CObject(TYPE){}
    
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(ModelConfigMapperUnitTest, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(ClientRunnable);
        EMPTY_SHELL_METHOD(defaultConstructor);
        FIELD(mapper, "The mapper to test");
        FIELD(lossConfig, "Mapper input");
    }
    
    static IObject* defaultConstructor(){
        return new ModelConfigMapperUnitTest();
    }
};

CClassConstSP const ModelConfigMapperUnitTest::TYPE =
    CClass::registerClassLoadMethod(
        "ModelConfigMapperUnitTest",
        typeid(ModelConfigMapperUnitTest),
        load);

/* external symbol to allow class to be forced to be linked in */
bool IModelConfigMapperLoad(){
    return (ModelConfigMapperUnitTest::TYPE != 0 && IModelConfigMapper::TYPE != 0);
}

DRLIB_END_NAMESPACE
