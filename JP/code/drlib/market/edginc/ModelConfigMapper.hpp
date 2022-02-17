//----------------------------------------------------------------------------
//
//   Group       : CH Quantitative Research
//
//   Date        : 01-Aug-2006
//
//----------------------------------------------------------------------------

#ifndef QLIB_MODELCONFIGMAPPER_HPP
#define QLIB_MODELCONFIGMAPPER_HPP

#include "edginc/IModelConfigMapper.hpp"
#include "edginc/Hashtable.hpp" // for Hashtable::StringHash
#include ext_hash_map

DRLIB_BEGIN_NAMESPACE

/**
 * Implementation of IModelConfigMapper: maps
 * a "ICreditLossModelConfig" to a "ICreditLossConfig" according
 * to its name or type.
 * 
 * NB: the mapping according to the "ICreditLossConfig" type is strict i.e. if L1 is a
 *     ICreditLossConfig, L2 derives from L1 and the user defines a rule to map L1 to
 *     a ICreditLossModelConfig M1 then there will be no "type rule" to map L2.
 * */
class MARKET_DLL ModelConfigMapper:
    public CObject,
    public virtual IModelConfigMapper
{
public:
	/** TYPE (for reflection) */
    static CClassConstSP const TYPE;
    
    /** Virtual destructor */
    virtual ~ModelConfigMapper();

    /**
     * Maps a "ICreditLossModelConfig" to a "ICreditLossConfig"
     * [Implements IModelConfigMapper]
     * */
    virtual ICreditLossModelConfigConstSP innerModel(
        ICreditLossConfigConstSP lossConfig) const;

    /** Called immediately after object constructed */
    virtual void validatePop2Object();
    
    /** Overrides clone() method to copy unregistered fields */
    virtual IObject* clone() const;

    /** Explicit constructor */
    ModelConfigMapper(
        StringArraySP types,
        ICreditLossModelConfigArraySP modelsForType);
    
private:
    ModelConfigMapper();
    static void load(CClassSP& clazz);
    static IObject* defaultConstructor();
    
    /** Builds internal maps nameToModelMap and typeToModelMap */
    void buildMaps();

    /* Name - model mapping (first mapping to be checked) */
    StringArraySP names; // $optional
    ICreditLossModelConfigArraySP modelsForName; // $optional

    /* Type - model mapping (second mapping to be checked) */
    StringArraySP types; // $optional
    ICreditLossModelConfigArraySP modelsForType; // $optional

    /* Model used by default (last "mapping" to be checked) */
    ICreditLossModelConfigSP defaultModel; // $optional

    typedef hash_map<string, ICreditLossModelConfigConstSP, Hashtable::StringHash> StringToModelHashTable;
    typedef pair<string, ICreditLossModelConfigConstSP> StringModelPair;

    StringToModelHashTable nameToModelMap; // $unregistered
    StringToModelHashTable typeToModelMap; // $unregistered
};

DRLIB_END_NAMESPACE

#endif /*MODELCONFIGMAPPER*/
