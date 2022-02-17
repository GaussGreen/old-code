//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Description : Helper class for models to get vols for 'CAssets'
//                 (ie 'spot' type assets) out of market cache
//
//   Author      : Mark A Robson
//
//----------------------------------------------------------------------------

#ifndef EDR_MARKETDATACONVERT_HPP
#define EDR_MARKETDATACONVERT_HPP

#include "edginc/MDFAssetVol.hpp"
#include <map>


DRLIB_BEGIN_NAMESPACE

/** Interface that should be implemented by objects that "can possibly" convert to other objects */
class MARKET_DLL MarketDataConvert {
public:
    typedef MarketObjectSP (*ConvertMethod)(IObjectSP);
    
    /** Registration: CClass describing iConvert interface and pointer to convert method */
    struct Reg {
        Reg();
        Reg(CClassConstSP iConvert, ConvertMethod method);
        
        CClassConstSP   iConvert;
        ConvertMethod   method;
    };

    /** Conversion interface */
    template <class T> 
    class IConvert: virtual public IObject {
    public:
        static CClassConstSP const TYPE;
        
        /** Converts to desired type. Returns null pointer if this is not possible */
        virtual smartPtr<T> convert(T* p = 0) const = 0;
    
    private:

        /** Templated convert method that casts to desired type and invokes convert */
        static MarketObjectSP _convert(IObjectSP marketObj) {
            ASSERT(!!marketObj);
            smartPtr<IConvert<T> > convIface = smartPtr<IConvert<T> >::dynamicCast(marketObj);
            smartPtr<T> res = convIface->convert();
            return res;
        }

        /** Invoked when Class is 'loaded' */
        static void load(CClassSP& clazz) {
            clazz->setPublic();
            REGISTER_INTERFACE(IConvert, clazz);
            EXTENDS(IObject);
            MarketDataConvert::registerConversionIFace(typeid(T).name(), Reg(clazz, &_convert));
        }
    };

    /** Registers (objType, Reg>) */
    static void registerConversionIFace(const string& objType, const Reg& reg);

    /** Returns Reg info corresponding to objType */
    static const Reg& getConversionInfo(CClassConstSP objType);

private:
    /** Provides access to the conversionTable. 
        Wrapped in static function to ensure initialization of static field */
    static map<string, Reg>& getTable();

    static map<string, Reg> conversionTable;    //!< Conversion registry
};


DRLIB_END_NAMESPACE

#endif
