/**
 *  Validator.hpp
 *  Implements a mechanism by which client can check for validity of an object
 *  via DRI
 *
 *  A   Q U I C K   G U I D E
 *  For any class, say Foo, that is to be validated by a client via DRI, provide a method
 *      bool Foo::validate()
 *  and in the implementation cpp file, add this line near the end
 *      DEFINE_VALIDATOR_TYPE(Foo);
 *  so that the template class is instantiated for the class Foo.
 */

#ifndef VALIDATOR_HPP
#define VALIDATOR_HPP

#include "edginc/config.hpp"
#include "edginc/Object.hpp"
#include "edginc/Class.hpp"
#include "edginc/ClientRunnable.hpp"
#include "edginc/Atomic.hpp"

DRLIB_BEGIN_NAMESPACE

// allows client to check for validity through DRI
// requires bool T::validate() for class T
template<typename T>
class MARKET_DLL Validator : public CObject, virtual public ClientRunnable
{
public:
    static CClassConstSP const TYPE;
    virtual ~Validator() {};
    IObjectSP run()
    {
        bool result = o->validate();
        return IObjectSP(CBool::create(result));
    }

protected:
    Validator(CClassConstSP clazz) : CObject(clazz), o(0) {}
    T* o;

private:
    Validator() : CObject(TYPE), o(0) {}
    static IObject* defaultValidator() { return new Validator(); }
    static void load(CClassSP& clazz)
    {
        REGISTER(Validator<T>, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(ClientRunnable);
        EMPTY_SHELL_METHOD(defaultValidator);
        FIELD(o, "Object to be validated");
        clazz->setPublic();
    }
};

DRLIB_END_NAMESPACE

#define DEFINE_VALIDATOR_TYPE(name)             \
    template <> \
    CClassConstSP const Validator<name>::TYPE = \
        CClass::registerClassLoadMethod(#name "::Validator", typeid(Validator<name>), Validator<name>::load); \

#endif  // VALIDATOR_HPP
