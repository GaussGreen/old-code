/**
 * @file InternalInterface.cpp
 */

#include "edginc/config.hpp"
#include "edginc/Addin.hpp"
#include <sstream>

DRLIB_BEGIN_NAMESPACE

class InternalInterface: public CObject {
public:
    static CClassConstSP const TYPE;       
    static IObjectSP run(InternalInterface* params){
        try {
            ostringstream s;

            const CClassVec& allClasses = CClass::allClasses();

            for (size_t c = 0; c < allClasses.size(); ++c) {
                const CClass& clazz = *allClasses[c];
                s << "«" << clazz.getName() << "» «" <<
                     clazz.getTypeInfo().name() << "» «";
                if (clazz.getSuperClass())
                    s << clazz.getSuperClass()->getName();
                s << "» «";
                CClassVec interfaces = clazz.getInterfaces();
                for (size_t i = 0; i < interfaces.size(); ++i) {
                    if (i) s << "·";
                    s << interfaces[i]->getName();
                }
                s << "» ¡\n";
                const CFieldArray& fields = clazz.getDeclaredFields();
                for (size_t f = 0; f < fields.size(); ++f) {
                    const CField& field = *fields[f];
                    s << "  «" << field.getName() << "» «" <<
                         field.getType()->getName() << "» " <<
                         field.getModifiers() <<
                         " «" << field.getDescription() << "»\n";
                }
                s << "¡\n";
            }
            
            return CStringSP(CString::create(s.str()));
        }
        catch (exception& e) {
            throw ModelException(e, "InternalInterface");
        } 
    }

private:
    InternalInterface(): CObject(TYPE) {}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(InternalInterface, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultInternalInterface);

        Addin::registerClassObjectMethod("IFACE_INTERNAL_INTERFACE_TESTER",
                                         Addin::UTILITIES,
                                         "tests EDR InternalInterface methods",
                                         InternalInterface::TYPE,
                                         true,
                                         Addin::returnHandle,
                                         (Addin::ObjMethod*)InternalInterface::run);
    }

    static IObject* defaultInternalInterface(){
        return new InternalInterface();
    }
};

CClassConstSP const InternalInterface::TYPE = 
CClass::registerClassLoadMethod(
    "InternalInterface", typeid(InternalInterface), load);

bool InternalInterfaceLoad() {
    return InternalInterface::TYPE != NULL;
}

DRLIB_END_NAMESPACE
