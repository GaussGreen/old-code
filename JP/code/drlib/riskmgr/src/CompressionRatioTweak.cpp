//----------------------------------------------------------------------------
//
//   Group       : Quantitative Research & Development
//
//   Filename    : CompressionRatioTweak.cpp
//
//   Description : Tweak for compression ratio
//
//   Author      : Antoine Gregoire
//
//   Date        : June 2005
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/GenericAllNameScalarShift.hpp"
#include "edginc/CompressionRatioTweak.hpp"

DRLIB_BEGIN_NAMESPACE

class CompressionRatioTweak:
    public GenericAllNameScalarShift<CompressionRatioTwk, true> {

public:

    CompressionRatioTweak(double shiftSize = 0.1):
        GenericAllNameScalarShift<CompressionRatioTwk, true>(
            TYPE, NAME, shiftSize) {}

    /**
     * Reflection type of this class
     */

    static CClassConstSP const TYPE;
    
    /**
     * Output name for first derivative
     */

    const static string NAME;

    /**
     * Shift size to use if none provided
     */

    const static double DEFAULT_SHIFT;
    
    
    static IObject* defaultConstructor() {
        return new CompressionRatioTweak(DEFAULT_SHIFT);
    }

private:

    /**
     * Invoked when Class is 'loaded'
     */

    // typedef here just to avoid having "," in SUPERCLASS macro
    typedef GenericAllNameScalarShift<CompressionRatioTwk, true> __GenericCompressionRatioTwk;
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(CompressionRatioTweak, clazz);
        SUPERCLASS(__GenericCompressionRatioTwk);
        EMPTY_SHELL_METHOD(defaultConstructor);
    }
};

CClassConstSP const CompressionRatioTweak::TYPE =
    CClass::registerClassLoadMethod(
        "CompressionRatioTweak",
        typeid(CompressionRatioTweak),
        CompressionRatioTweak::load);
        
const string CompressionRatioTweak::NAME = "COMPRESSION_RATIO";

const double CompressionRatioTweak::DEFAULT_SHIFT = 0.1;


template<> CClassConstSP const GenericAllNameScalarShift<CompressionRatioTwk, true>::TYPE =
    CClass::registerClassLoadMethod(
        "GenericAllNameScalarShift<CompressionRatioTwk, true>",
        typeid(GenericAllNameScalarShift<CompressionRatioTwk, true>),
        GenericAllNameScalarShift<CompressionRatioTwk, true>::load);

template<> CClassConstSP const TweakableWith<CompressionRatioTwk>::TYPE =
    CClass::registerInterfaceLoadMethod(
        "TweakableWith<CompressionRatioTwk>",
        typeid(TweakableWith<CompressionRatioTwk>),
        TweakableWith<CompressionRatioTwk>::load);

template<> CClassConstSP const RestorableWith<CompressionRatioTwk>::TYPE =
    CClass::registerInterfaceLoadMethod(
        "RestorableWith<CompressionRatioTwk>",
        typeid(RestorableWith<CompressionRatioTwk>),
        RestorableWith<CompressionRatioTwk>::load);

/**
 * Included in RiskMgrLib::linkInClasses() to force CompressionRatioTweak
 * to get linked into the exe.
 */
bool CompressionRatioTweakLinkIn() {
    return CompressionRatioTweak::TYPE != NULL;
}

DRLIB_END_NAMESPACE

