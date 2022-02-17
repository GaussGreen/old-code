//----------------------------------------------------------------------------
//
//   Group       : QR Rates
//
//   Filename    : MarketTable.hpp
//
//   Description : MarketTable - maps key names to market object handles.
//
//   Author      : Anwar E Sidat
//
//   Date        : 07-Sep-2006
//
//----------------------------------------------------------------------------

#ifndef QLIB_MarketTable_HPP
#define QLIB_MarketTable_HPP

#include <map>
#include "edginc/Addin.hpp"
#include "edginc/MarketObject.hpp"

DRLIB_BEGIN_NAMESPACE

/** Template for family of any IObject type objects. */
template <typename X>
class MarketTable : public MarketObject
{
public:
    static CClassConstSP const TYPE;

    MarketTable();
    virtual ~MarketTable();

    /** Get object given key string. */
    virtual smartPtr<X> getMarketObject(const string& strKey) const;

    /** Returns name of model. */
    virtual string getName() const;

    /** overrides default */
    virtual void validatePop2Object();

    /** overrides clone */
    IObject* clone() const;

    typedef smartConstPtr<X> XConstSP;
    typedef smartPtr<X>      XSP;
    typedef array<XSP, X>    X_Array;

    MarketTable(const CClassConstSP& clazz);
    MarketTable(const MarketTable& irv);
    MarketTable& operator=(const MarketTable& rhs);

protected:

    // Fields
    string            Name;         // Handle name
    StringArray       KeyArray;     // array of key names
    X_Array           ObjArray;     // array of objects

private:
    static void load(CClassSP& clazz);
    static IObject* defaultConstructor(void) { return new MarketTable<X>(); }
    
    DLL_FIX_FOR_TEMPLATE_TYPE; // work around for VC71 and dlls
    typedef std::map<std::string, XSP> ObjectMap;
    ObjectMap  objMap;
};

/*
template <typename X> TEMPLATE_INLINE bool MarketTable<X>::MarketTableLoad()
{
    return (MarketTable<X>::TYPE != 0);
}
*/

template <typename X> TEMPLATE_INLINE MarketTable<X>::MarketTable() : MarketObject(MarketTable<X>::TYPE) {}

template <typename X> TEMPLATE_INLINE MarketTable<X>::~MarketTable() {}

template <typename X> TEMPLATE_INLINE MarketTable<X>::MarketTable(const CClassConstSP& clazz)
    : MarketObject(clazz)
{
    validatePop2Object();
}

template <typename X> TEMPLATE_INLINE IObject* MarketTable<X>::clone() const
{
    IObject* pIObject = MarketObject::clone();
    pIObject->validatePop2Object();
    return pIObject;
}


template <typename X> MarketTable<X>& MarketTable<X>::operator=(const MarketTable<X>& rhs)
{
    return *this;
}

template <typename X> TEMPLATE_INLINE smartPtr<X> MarketTable<X>::getMarketObject(const string& strKey) const
{
    static const string method = "MarketTable::getMarketObject";
    typename ObjectMap::const_iterator iter;
    iter = objMap.find(strKey);
    if (iter == objMap.end())
        throw ModelException(method, getName() + " - Key " + strKey + " not found in MarketTable " + getName());
    return iter->second;
}

template <typename X> TEMPLATE_INLINE void MarketTable<X>::validatePop2Object()
{
    static const string method = "MarketTable::validatePop2Object";
    try
    {
        // Check Number of Factors
        if (KeyArray.size() != ObjArray.size())
        {
             throw ModelException(method, getName() + " - KeyArray size [" + Format::toString(KeyArray.size()) +
                "] does not match size on ObjArray [" + Format::toString(ObjArray.size()) + "]");
        }

        // Set up map
        typename ObjectMap::iterator iter;
        objMap.clear();
        for (int i = 0; i < KeyArray.size(); ++i)
        {
            // No duplicates allowed.
            iter = objMap.find(KeyArray[i]);
            if (iter != objMap.end())
                throw ModelException(method, getName() + " - Duplicate key " + KeyArray[i] + " found in MarketTable " + getName());
            objMap.insert(typename ObjectMap::value_type(KeyArray[i], ObjArray[i]));
        }
    }
    catch (exception& e)
    {
        throw ModelException(e, method);
    }
}

/** Returns name of model */
template <typename X> TEMPLATE_INLINE string MarketTable<X>::getName() const
{
    return Name;
}

template <typename X> TEMPLATE_INLINE void MarketTable<X>::load(CClassSP& clazz)
{
    clazz->setPublic();
    clazz->setDescription("Object Table.");
    REGISTER(MarketTable<X>, clazz);
    SUPERCLASS(MarketObject);
    EMPTY_SHELL_METHOD(defaultConstructor);
    FIELD(Name,     "Identifier for this MarketTable object.");
    FIELD(KeyArray, "Array of key names for each corresponding object.");
    FIELD(ObjArray, "Array of objects.");

    // Addin name is MarketTable_X (angled brackets not allowed).
    string strFnName("MarketTable_");
    string strType(typeid(X).name());
    size_t pos = strType.find("::");
    strFnName += string( strType.substr(pos + 2, strType.size() ));

    Addin::registerConstructor(strFnName.c_str(),
                               Addin::MARKET,
                               "Creates a handle to an object table (maps key names to object handles).",
                               MarketTable<X>::TYPE);
}


// This definition should get moved into coreConfig.hpp once finalized.
#define DEFINE_CONTAINER_TEMPLATE_TYPE(T, X) \
template <typename X> CClassConstSP const T < X > ::TYPE = \
CClass::registerClassLoadMethod(#T "<" #X ">", typeid( T < X >), T < X > ::load);


DRLIB_END_NAMESPACE

#endif
