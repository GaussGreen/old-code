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

#include "edginc/config.hpp"
#include "edginc/MarketData.hpp"
#include "edginc/MarketDataConvert.hpp"
#include "edginc/VolSurface.hpp"
#include "edginc/FlatFXVol.hpp"

DRLIB_BEGIN_NAMESPACE

MarketDataConvert::Reg::Reg():
iConvert(0), method(0) {}

MarketDataConvert::Reg::Reg(CClassConstSP iConvert, ConvertMethod method):
iConvert(iConvert), method(method) {}


void MarketDataConvert::registerConversionIFace(const string& objType, const Reg& reg) {
    map<string, Reg>& table = getTable();
    table[objType] = reg;
}


const MarketDataConvert::Reg& MarketDataConvert::getConversionInfo(CClassConstSP objType) {
    map<string, Reg>& table = getTable();
    return table[objType->getTypeInfo().name()];
}


map<string, MarketDataConvert::Reg>& MarketDataConvert::getTable() {
    return conversionTable;
}


map<string, MarketDataConvert::Reg> MarketDataConvert::conversionTable = 
map<string, MarketDataConvert::Reg>();


DRLIB_END_NAMESPACE
