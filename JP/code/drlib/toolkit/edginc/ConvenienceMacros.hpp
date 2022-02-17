#ifndef __CONVENIENCEMACROS_HPP__
#define __CONVENIENCEMACROS_HPP__

#define QUICK_CLASS_REFLECTION(className) \
	static IObject* defaultConstructor(void) { return new className(); } \
	static void load(CClassSP& clazz); /* declare what is being exposed */ \
	static CClassConstSP const TYPE;   /* type information for this class */

#define QUICK_IFACE_REFLECTION(className) \
	static void load(CClassSP& clazz); /* declare what is being exposed */ \
	static CClassConstSP const TYPE;   /* type information for this class */

#define QUICK_CLASS_LOAD(className) \
	clazz->setPublic(); /* can be seen from outside the library */ \
	REGISTER(className, clazz); /* always the same, mandatory */ \
	EMPTY_SHELL_METHOD(defaultConstructor); /* always the same, mandatory */

#define QUICK_IFACE_LOAD(className) \
	clazz->setPublic(); /* can be seen from outside the library */ \
	REGISTER_INTERFACE(className, clazz); /* always the same, mandatory */ \

#define REGISTER_MYINTERFACE(className) \
CClassConstSP const className::TYPE = CClass::registerInterfaceLoadMethod(#className, typeid(className), className::load)

#define REGISTER_MYCLASS(className) \
CClassConstSP const className::TYPE = CClass::registerClassLoadMethod(#className, typeid(className), className::load)

#define STRINGIFY(x) #x
#define DBi ErrorHandler::writeMsg("debug " STRINGIFY(__LINE__) );

#define ModEx ModelException
#define ModExS(x) ModelException(St2St(x))

#endif
