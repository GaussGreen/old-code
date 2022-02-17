#include "edginc/config.hpp"
#include "edginc/ModelConfigMapper.hpp"

DRLIB_BEGIN_NAMESPACE

/** TYPE (for reflection) */        
CClassConstSP const ModelConfigMapper::TYPE =
CClass::registerClassLoadMethod(
    "ModelConfigMapper",
    typeid(ModelConfigMapper),
    ModelConfigMapper::load);

/** Virtual destructor */
ModelConfigMapper::~ModelConfigMapper(){}

/** Explicit constructor */
ModelConfigMapper::ModelConfigMapper(
    StringArraySP types,
    ICreditLossModelConfigArraySP modelsForType):
        CObject(TYPE),
        names(0),
        modelsForName(0),
        types(types),
        modelsForType(modelsForType),
        defaultModel(0)
{
    validatePop2Object();
}

/** Constructor */
ModelConfigMapper::ModelConfigMapper(): CObject(TYPE),
    names(0),
    modelsForName(0),
    types(0),
    modelsForType(0),
    defaultModel(0) {}

IObject* ModelConfigMapper::defaultConstructor()
{
    return new ModelConfigMapper();
}

void ModelConfigMapper::load(CClassSP& clazz)
{
    clazz->setPublic();
    REGISTER(ModelConfigMapper, clazz);
    SUPERCLASS(CObject);
    IMPLEMENTS(IModelConfigMapper);
    EMPTY_SHELL_METHOD(defaultConstructor);
    FIELD(names, "Names in name - model mapping (first mapping to be checked)");
    FIELD_MAKE_OPTIONAL(names);
    FIELD(modelsForName, "Models in name - model mapping (first mapping to be checked)");
    FIELD_MAKE_OPTIONAL(modelsForName);
    FIELD(types, "Types in type - model mapping (second mapping to be checked)");
    FIELD_MAKE_OPTIONAL(types);
    FIELD(modelsForType, "Models in type - model mapping (second mapping to be checked)");
    FIELD_MAKE_OPTIONAL(modelsForType);
    FIELD(defaultModel, "Model used by default (last mapping to be checked)");
    FIELD_MAKE_OPTIONAL(defaultModel);
}

/**
 * Maps a "ICreditLossModelConfig" to a "ICreditLossConfig"
 * [Implements IModelConfigMapper]
 * */
ICreditLossModelConfigConstSP ModelConfigMapper::innerModel(
    ICreditLossConfigConstSP lossConfig) const
{
    string name = lossConfig->getName();
    string type = lossConfig->getClass()->getName();
    
    StringToModelHashTable::const_iterator nameIter = nameToModelMap.find(name);
    if (nameIter != nameToModelMap.end())
    {
        // found a model corresponding to that name - use it !
        return nameIter->second;
    }
    else
    {
        // not found, try in typeToModelMap
        StringToModelHashTable::const_iterator typeIter = typeToModelMap.find(type);
        if (typeIter != typeToModelMap.end())
        {
            // found a model corresponding to that type - use it !
            return typeIter->second;
        }
    }

    if (defaultModel.get() != 0)
    {
        // use default loss model config if no match found
        return defaultModel;
    }
    else
    {
        throw ModelException(
            "ModelConfigMapper::innerModel",
            "No loss model config found for " + name + " with type " + type + ".");
    }
}

/** Called immediately after object constructed */
void ModelConfigMapper::validatePop2Object()
{
    static const string method("ModelConfigMapper::validatePop2Object");        
	try
	{
        // Checks array size
        
        // names and modelsForName
        if (names.get() != 0)
        {
            if (modelsForName.get() != 0)
            {
                if (names->size() != modelsForName->size())
                {
                    throw ModelException(
                        method,
                        "Arrays 'names' and 'modelsForName' don't have the same size.");
                }
            }
            else
            {
                throw ModelException(
                    method,
                    "Array 'modelsForName' is empty whereas array 'names' is not.");
            }    
        }
        else
        {
            if (modelsForName.get() != 0)
            {
                throw ModelException(
                    method,
                    "Array 'names' is empty whereas array 'modelsForName' is not.");
            }
        }

        // types and modelsForType
        if (types.get() != 0)
        {
            if (modelsForType.get() != 0)
            {
                if (types->size() != modelsForType->size())
                {
                    throw ModelException(
                        method,
                        "Arrays 'types' and 'modelsForType' don't have the same size.");
                }
            }
            else
            {
                throw ModelException(
                    method,
                    "Array 'modelsForType' is empty whereas array 'types' is not.");
            }    
        }
        else
        {
            if (modelsForType.get() != 0)
            {
                throw ModelException(
                    method,
                    "Array 'types' is empty whereas array 'modelsForType' is not.");
            }
        }
        
        // Checks at least one optional field is populated
        if (names.get() == 0 && types.get() == 0 && defaultModel.get() == 0)
        {
            throw ModelException(
                method,
                "At least one field should be populated "
                "('names' and 'modelsForName' or "
                "'types' and 'modelsForType' or "
                "'defaultModel')");
        }
                
        // Checks "types" are valid i.e. they derive from ICreditLossConfig
        if (types.get() != 0)
        {
            for (int i = 0; i < types->size(); ++i) {
                if (!ICreditLossConfig::TYPE->isAssignableFrom(CClass::forName((*types)[i])))
                {
                    throw ModelException(
                        method,
                        "Invalid type: " + (*types)[i] +
                        " does not derive from " +
                        ICreditLossConfig::TYPE->getName());
                }
			}
        }
        
        // Builds maps
        buildMaps();
		
	} catch (exception& e){
	    throw ModelException(e, method);
	}
}

/** Overrides clone() method to copy unregistered fields */
IObject* ModelConfigMapper::clone() const
{
    IObject* myCopy = CObject::clone();
    ModelConfigMapper& modelConfigMapper = 
        dynamic_cast<ModelConfigMapper&>(*myCopy);
    modelConfigMapper.nameToModelMap = nameToModelMap;
    modelConfigMapper.typeToModelMap = typeToModelMap;
    return myCopy;
}

/** Builds internal maps nameToModelMap and typeToModelMap */
void ModelConfigMapper::buildMaps()
{
    int i;
    static const string method("ModelConfigMapper::buildMaps");        
    
    // builds nameToModelMap
    if (names.get() != 0)
    {
        for (i = 0; i < names->size(); ++i) {
    		// insert element and test if already exists in the same line
            if (!nameToModelMap.insert(
                StringModelPair((*names)[i], (*modelsForName)[i])).second)
            {
                throw ModelException(
                    method,
                    "Name " + (*names)[i] + " appears at least twice in 'names', "
                    "making mapping ambiguous.");
            }
    	}
    }

    // builds typeToModelMap
    if (types.get() != 0)
    {
        for (i = 0; i < types->size(); ++i) {
            // insert element and test if already exists in the same line
            if (!typeToModelMap.insert(
                StringModelPair((*types)[i], (*modelsForType)[i])).second)
            {
                throw ModelException(
                    method,
                    "Type " + (*types)[i] + " appears at least twice in 'types', "
                    "making mapping ambiguous.");
            }
        }
    }
}

/* external symbol to allow class to be forced to be linked in */
bool ModelConfigMapperLoad(){
    return (ModelConfigMapper::TYPE != 0);
}

DRLIB_END_NAMESPACE
