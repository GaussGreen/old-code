//----------------------------------------------------------------------------
//
//   Filename    : ObservationSource.cpp
//
//   Description : Classes for observation sources - for use with AssetHistory
//                 NB - the source can be a single source or a collection (e.g for struck)
//
//   Author      : Ian Stares   
//
//   Date        : April 5 2006
//

//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#include "edginc/ObservationSource.hpp"

DRLIB_BEGIN_NAMESPACE

// base class

ObservationSource::ObservationSource(CClassConstSP clazz) : CObject(clazz) {}

ObservationSource::ObservationSource(CClassConstSP clazz, const string& source) 
        : CObject(clazz), primarySource(source) {}

ObservationSource::ObservationSource() : CObject(TYPE) {}

ObservationSource::ObservationSource(const string& source) : CObject(TYPE),
                primarySource(source) {}

// note not public so not available in EAS/Excel
// will probably build these on fly inside QLib from lists of strings??
void ObservationSource::load(CClassSP& clazz) {
    REGISTER(ObservationSource, clazz);
    SUPERCLASS(CObject);    
    FIELD(primarySource, "The primary source for the observable");
    EMPTY_SHELL_METHOD(defaultObsSource);
}

IObject* ObservationSource::defaultObsSource() {
    return new ObservationSource();
}

string ObservationSource::getPrimarySource() const {
    return primarySource;
}   

CClassConstSP const ObservationSource::TYPE = CClass::registerClassLoadMethod(
    "ObservationSource", typeid(ObservationSource), ObservationSource::load);

DEFINE_TEMPLATE_TYPE(ObservationSourceArray);

/*****************************************************************************/
// FX Observation Source -  a pair of sources for the asset/FX
StruckObservationSource::StruckObservationSource() : ObservationSource(TYPE) {}

// note not public so not available in EAS/Excel
// will probably build these on fly inside QLib from lists of strings??
void StruckObservationSource::load(CClassSP& clazz) {
    REGISTER(StruckObservationSource, clazz);
    SUPERCLASS(ObservationSource);    
    FIELD(fxSource, "The source for the fx rate look up");
    EMPTY_SHELL_METHOD(defaultStruckObsSource);
}

IObject* StruckObservationSource::defaultStruckObsSource() {
    return new StruckObservationSource();
}

StruckObservationSource::StruckObservationSource(const string& source,
                                                 const string& fx) 
                : ObservationSource(TYPE, source) {
    fxSource = ObservationSourceSP(new ObservationSource(fx));
}

ObservationSource* StruckObservationSource::getFXSource() const {
    return fxSource.get();
}   

CClassConstSP const StruckObservationSource::TYPE = 
    CClass::registerClassLoadMethod("StruckObservationSource", 
        typeid(StruckObservationSource), StruckObservationSource::load);

DRLIB_END_NAMESPACE




