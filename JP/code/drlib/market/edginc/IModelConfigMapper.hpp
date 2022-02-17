//----------------------------------------------------------------------------
//
//   Group       : CH Quantitative Research
//
//   Date        : 01-Aug-2006
//
//----------------------------------------------------------------------------

#ifndef QLIB_IMODELCONFIGMAPPER_HPP
#define QLIB_IMODELCONFIGMAPPER_HPP

#include "edginc/ICreditLossModelConfig.hpp"
#include "edginc/ICreditLossConfig.hpp"

DRLIB_BEGIN_NAMESPACE

/**
 * Interface defining a way to map a "ICreditLossModelConfig" to a "ICreditLossConfig".
 * Typically used in ICreditLossModelConfig.
 * */
class MARKET_DLL IModelConfigMapper: public virtual IObject {
public:
    /** TYPE (for reflection) */
    static CClassConstSP const TYPE;
    
    /** Maps a "ICreditLossModelConfig" to a "ICreditLossConfig" */
    virtual ICreditLossModelConfigConstSP innerModel(
        ICreditLossConfigConstSP lossConfig) const = 0;
};

DECLARE(IModelConfigMapper);

DRLIB_END_NAMESPACE

#endif /*QLIB_IMODELCONFIGMAPPER_HPP*/
