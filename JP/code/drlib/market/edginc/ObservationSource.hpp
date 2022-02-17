//----------------------------------------------------------------------------
//
//   Filename    : ObservationSource.hpp
//
//   Description : Classes for observation source - for use with AssetHistory
//                 NB - the source can be a single source or a collection (e.g for struck)
//
//   Author      : Ian Stares   
//
//   Date        : April 5 2006
//

//
//----------------------------------------------------------------------------

#ifndef OBSSOURCE_HPP
#define OBSSOURCE_HPP

#include "edginc/config.hpp"
#include "edginc/Object.hpp"
#include "edginc/Class.hpp"
#include "edginc/Array.hpp"
#include "edginc/AtomicArray.hpp"

DRLIB_BEGIN_NAMESPACE

// base class for observation source
// basically a single primary source for observations
class MARKET_DLL ObservationSource :  public CObject {
public:
    static CClassConstSP const TYPE;

    string getPrimarySource() const;

    ObservationSource(const string& source);

protected:
    // for inheritance
    ObservationSource(CClassConstSP clazz);
    ObservationSource(CClassConstSP clazz, const string& source);

private:
    ObservationSource(); // for reflection
    static void load(CClassSP& clazz);
    static IObject* defaultObsSource();

    string primarySource;
};

typedef smartPtr<ObservationSource> ObservationSourceSP;
typedef array<ObservationSourceSP, ObservationSource> ObservationSourceArray;
typedef smartPtr<ObservationSourceArray> ObservationSourceArraySP;

// struck observation source a primary source for the asset 
// and an FX source for the FX rate
class MARKET_DLL StruckObservationSource :  public ObservationSource {
public:
    static CClassConstSP const TYPE;

    ObservationSource* getFXSource() const;

    StruckObservationSource(const string& source,
                            const string& fxSource);

private:
    StruckObservationSource(); // for reflection
    static void load(CClassSP& clazz);
    static IObject* defaultStruckObsSource();

    ObservationSourceSP fxSource;
};
DRLIB_END_NAMESPACE

#endif

