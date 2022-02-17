//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : AssetNameCollector.hpp
//
//   Description : Asset name collector class
//
//   Author      : Andre Segger
//
//   Date        : 30 Apr 2001
//
//
//----------------------------------------------------------------------------

#ifndef ASSETNAME_COLLECT_HPP
#define ASSETNAME_COLLECT_HPP
#include "edginc/Collector.hpp"
#include "edginc/DateTime.hpp"

DRLIB_BEGIN_NAMESPACE

/** Validates the uniqueness of asset names. */
class MARKET_DLL AssetNameCollector: public CObject, public virtual ICollector {
public:
    static CClassConstSP const TYPE;
    /** add the supplied name to the set */
    void assetNameValidate(string   name);

    AssetNameCollector();
    ~AssetNameCollector();
private:
    AssetNameCollector(const AssetNameCollector& rhs);
    AssetNameCollector& operator=(const AssetNameCollector& rhs);
    // hide implementation in separate class
    class Imp;
    auto_ptr<Imp> my; // $unregistered

    static void load(CClassSP& clazz);
};

typedef smartPtr<AssetNameCollector> AssetNameCollectorSP;
typedef smartPtr<const AssetNameCollector> AssetNameCollectorConstSP;

DRLIB_END_NAMESPACE

#endif
