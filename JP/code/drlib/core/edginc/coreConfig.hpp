#ifndef CORE_CONFIG_HH

#    define CORE_CONFIG_HH

#if 0 && defined(LINUX64)
#include <common.hpp>
#endif

#    if defined (_MSC_VER)
#        ifdef WIN32
#            pragma warning( disable : 4786 )
					// truncation in debug info
#            pragma warning( disable : 4503 )
					// truncation in decorated name
#            pragma warning( disable : 4305 )
					// truncation from type A to type B
#            pragma warning( disable : 4244 )
					// truncation from type A to type B
// Please keep this warning - it highlights bugs/poor code
//#pragma warning( disable : 4101 ) // unreferenced local variable

// "resolved overload was found by argument-dependent lookup"
// http://msdn.microsoft.com/library/default.asp?url=/library/en-us/vclang/html/vclrfargumentdependentkoeniglookupnowsupported.asp
// http://en.wikipedia.org/wiki/Koenig_Lookup
// --- was causing lots of confusing warnings basically boiling down to
// "I used stl::exp(complex) which overloads plain old ::exp(double) because
// the argument was complex, and they're in different namespaces"
#            pragma warning( disable : 4675 )
// make this warning an error. Basically means delete operator has been
// instantiated incorrectly (because class declared but not defined)
#            pragma warning( error : 4150)
#        endif
#    endif

//// using safe mode STL for debug build should be seriously considered
//// for future testing. unfortunately the code is not fully compliant as of now
//#if defined (_MSC_VER)
//# ifdef _NDEBUG
//# undef _STLP_DEBUG
//# else
//# define _STLP_DEBUG 1
//# endif
//#endif

// Namespace supported by all compilers we use
// Core library lives in "core" namespace

#    ifndef CORE_NAMESPACE
#        define CORE_NAMESPACE core
#        define CORE_BEGIN_NAMESPACE namespace CORE_NAMESPACE{
#        define CORE_END_NAMESPACE }
#        define USING_CORE_NAMESPACE using namespace CORE_NAMESPACE
#    endif

// default namespace to drlib if NAMESPACE not defined
#    ifndef NAMEDRLIB
#        define NAMEDRLIB drlib
#    endif

#    ifndef NAMESPACE
#        define NAMESPACE ::NAMEDRLIB
#    endif


// The structure below is to identify automatically all and any cases when 
// namespaces were accidentally nested. e.g. DRLIB_BEGIN_NAMESPACE preceded 
// the #include statement
//
// If you see a message similar to:
// src\SRMBasisSpreadHJMUtil.cpp(20) : error C2556: 
//      'int drlib::drlib::namespaceChecker__(void)' : overloaded function differs only by return type
// then you know that you have drlib namespace opened within drlib namespace.
namespace NAMEDRLIB
{
  int namespaceChecker__ ();
  namespace NAMEDRLIB
  {
	char *namespaceChecker__ ();
  }
}

#    define DRLIB_BEGIN_NAMESPACE namespace NAMEDRLIB{ int namespaceChecker__();
#    define DRLIB_END_NAMESPACE }
#    define USING_DRLIB_NAMESPACE using namespace NAMESPACE;
#    define DRXML_BEGIN_NAMESPACE namespace drxml{
#    define DRXML_END_NAMESPACE }
#    define USING_DRXML_NAMESPACE using namespace drxml;



// standard workaround for the nonstandard scope of "for(int i; ... ; ... )" in VC6
// see e.g. "http://boost-consulting.com/boost/more/microsoft_vcpp.html"
#    if defined(_MSC_VER)
#        if (_MSC_VER <1300)
#            define for if(0) ; else for

// QLIB_BUILD_DLL means we're building individual dlls on windows
// (ie one per directory).
// this is for stlport-4.5.3 with VC6 (ie native io streams)
#            if defined(QLIB_BUILD_DLL)
#                define _STLP_USE_DECLSPEC
// these appear to be bogus (not completely sure though). Currently only
// switched off for VC6 and using stlport with native iostreams
#                pragma warning( disable : 4273 )
#                pragma warning( disable : 4661 )
#            endif
#        else
// force MSVC7 to do conform to C++ standard wrt for loops
#            pragma conform(forScope, on)
#        endif
#    endif

//// Legacy - this is to be removed once a couple of files are updated
#    define EDG_DLL

#    if defined (_MSC_VER) && defined(QLIB_BUILD_DLL)
#        pragma warning( disable : 4251 )
					/* warning if parent class is not exported
					   but class is - affects lots of templates */
#        define QLIB_DLL_EXPORT __declspec(dllexport)
#        define QLIB_DLL_IMPORT __declspec(dllimport)
#        if (_MSC_VER <1300)
// MSVC6 at least has a bug where templated functions within say an ordinary
// class are deemed to have been exported (although clearly impossible). Same
// is also true for methods in a templated class which use another templated
// type (eg template<class _Tp1> smartPtr(const smartPtr<_Tp1>& __r) in
// smartPtr.hpp). 
// It also has another bug which means that if you do something like
// extern template __declspec(dllexport) smartPtr<OutputNameArray>)
// and then OutputNameArraySP is used before OutputName is defined it falls over
// So when externing smartPtr we never use __declspec(dllexport) or 
// __declspec(dllimport).
#            define QLIB_DLL_IMPORT_SP
#            define QLIB_DLL_EXPORT_SP
#        else					//VC71
#            pragma warning( disable : 4661 )
					/* warning if no suitable definition found
					   for explicit template instantiation */
// MSVC71 essentially has a similar bug to the 2nd one above, but here it
// tries to instantiate the parents methods (and falls over if a class has
// been declared but not defined)
#            define QLIB_DLL_IMPORT_SP QLIB_DLL_IMPORT
#            define QLIB_DLL_EXPORT_SP

#        endif					// end of if VC6 or 7
#    else						// if not VC
#        if defined(__GNUC__) && (__GNUC__ >= 4)
#            define QLIB_DLL_IMPORT
#            define QLIB_DLL_IMPORT_SP
#            define QLIB_DLL_EXPORT __attribute__ ((visibility("default")))
#            define QLIB_DLL_EXPORT_SP __attribute__ ((visibility("default")))
#        else
#            define QLIB_DLL_IMPORT
#            define QLIB_DLL_IMPORT_SP
#            define QLIB_DLL_EXPORT
#            define QLIB_DLL_EXPORT_SP
#        endif
#    endif

#    if defined(QLIB_BUILD_DLL) && defined(_MSC_VER) && (_MSC_VER >=1300)
// ie building dlls on windows with vc71.
// VC71 also has a highly annoying feature that if there is a specialisation
// for a template of a certain type (think TYPE field) as well as a 'general'
// case, then the linker seems to choose either the wrong one or a random one
// if the specialisation is in a source file.
// Our fix for this is to have a special DLL_TYPE field on affected templates
// and initialise that.
#        define DLL_FIX_FOR_TEMPLATE_TYPE static CClassConstSP const DLL_TYPE;

#        define DEFINE_TEMPLATE_TYPE(type) \
template <> CClassConstSP const type ::DLL_TYPE = \
CClass::registerClassLoadMethod(#type, typeid(type), load);

#        define DEFINE_TEMPLATE_TYPE_WITH_NAME(name, type) \
template <> CClassConstSP const type ::DLL_TYPE = \
CClass::registerClassLoadMethod(name, typeid(type), load);

#        define DEFINE_INTERFACE_TEMPLATE_TYPE(type) \
template <> CClassConstSP const type ::DLL_TYPE = \
CClass::registerInterfaceLoadMethod(#type, typeid(type), load);

#        define DEFINE_INTERFACE_TEMPLATE_TYPE_WITH_NAME(name, type) \
template <> CClassConstSP const type ::DLL_TYPE = \
CClass::registerInterfaceLoadMethod(name, typeid(type), load);

#    else

#        define DLL_FIX_FOR_TEMPLATE_TYPE

#        define DEFINE_TEMPLATE_TYPE(type) \
template <> CClassConstSP const type ::TYPE = \
CClass::registerClassLoadMethod(#type, typeid(type), load);

#        define DEFINE_TEMPLATE_TYPE_WITH_NAME(name, type) \
template <> CClassConstSP const type ::TYPE = \
CClass::registerClassLoadMethod(name, typeid(type), load);

#        define DEFINE_INTERFACE_TEMPLATE_TYPE(type) \
template <> CClassConstSP const type ::TYPE = \
CClass::registerInterfaceLoadMethod(#type, typeid(type), load);

#        define DEFINE_INTERFACE_TEMPLATE_TYPE_WITH_NAME(name, type) \
template <> CClassConstSP const type ::TYPE = \
CClass::registerInterfaceLoadMethod(name, typeid(type), load);

#    endif

//// by default we set up the dll linkage to be importing (except for the one
//// we are actually compiling). The _SP stands for smart Ptr and is for working
//// around bug described above in definition of QLIB_DLL_IMPORT_SP
#    define CORE_DLL QLIB_DLL_IMPORT
#    define CORE_DLL_SP QLIB_DLL_IMPORT_SP
#    define RNG_DLL QLIB_DLL_IMPORT
#    define RNG_DLL_SP QLIB_DLL_IMPORT_SP
#    define TOOLKIT_DLL QLIB_DLL_IMPORT
#    define TOOLKIT_DLL_SP QLIB_DLL_IMPORT_SP
#    define UTIL_DLL QLIB_DLL_IMPORT
#    define UTIL_DLL_SP QLIB_DLL_IMPORT_SP
#    define RISKMGR_DLL QLIB_DLL_IMPORT
#    define RISKMGR_DLL_SP QLIB_DLL_IMPORT_SP
#    define MARKET_DLL QLIB_DLL_IMPORT
#    define MARKET_DLL_SP QLIB_DLL_IMPORT_SP
#    define CREDIT_DLL QLIB_DLL_IMPORT
#    define CREDIT_DLL_SP QLIB_DLL_IMPORT_SP
#    define FOURIER_DLL QLIB_DLL_IMPORT
#    define FOURIER_DLL_SP QLIB_DLL_IMPORT_SP
#    define RATESLIB_DLL QLIB_DLL_IMPORT
#    define RATESLIB_DLL_SP QLIB_DLL_IMPORT_SP
#    define TREE_DLL QLIB_DLL_IMPORT
#    define TREE_DLL_SP QLIB_DLL_IMPORT_SP
#    define MCARLO_DLL QLIB_DLL_IMPORT
#    define MCARLO_DLL_SP QLIB_DLL_IMPORT_SP
#    define CONVOLUTION_DLL QLIB_DLL_IMPORT
#    define CONVOLUTION_DLL_SP QLIB_DLL_IMPORT_SP
#    define PRODUCTS_DLL QLIB_DLL_IMPORT
#    define PRODUCTS_DLL_SP QLIB_DLL_IMPORT_SP
#    define IRPRODUCTS_DLL QLIB_DLL_IMPORT
#    define IRPRODUCTS_DLL_SP QLIB_DLL_IMPORT_SP
#    define PYPRODUCTS_DLL QLIB_DLL_IMPORT
#    define PYPRODUCTS_DLL_SP QLIB_DLL_IMPORT_SP
#    define ADDINS_DLL QLIB_DLL_IMPORT
#    define ADDINS_DLL_SP QLIB_DLL_IMPORT_SP
#    define REGTEST_DLL QLIB_DLL_IMPORT
#    define REGTEST_DLL_SP QLIB_DLL_IMPORT_SP
#    define RADAR_DLL  QLIB_DLL_IMPORT
#    define RADAR_DLL_SP QLIB_DLL_IMPORT_SP

//// then manually define each one to export where needed
#    ifdef QLIB_TOOLKIT
#        undef TOOLKIT_DLL
#        undef TOOLKIT_DLL_SP
#        define TOOLKIT_DLL QLIB_DLL_EXPORT
#        define TOOLKIT_DLL_SP QLIB_DLL_EXPORT_SP
#        if defined(_MSC_VER) && (_MSC_VER <1300) && defined(QLIB_BUILD_DLL)
// this is for stlport-4.5.3 with VC6 (ie native io streams)
#            define _STLP_DESIGNATED_DLL
#        endif
#    endif

#    ifdef QLIB_CORE
#        undef CORE_DLL
#        undef CORE_DLL_SP
#        define CORE_DLL QLIB_DLL_EXPORT
#        define CORE_DLL_SP QLIB_DLL_EXPORT_SP
#    endif

#    ifdef QLIB_RNG
#        undef RNG_DLL
#        undef RNG_DLL_SP
#        define RNG_DLL QLIB_DLL_EXPORT
#        define RNG_DLL_SP QLIB_DLL_EXPORT_SP
#    endif

#    ifdef QLIB_UTIL
#        undef UTIL_DLL
#        undef UTIL_DLL_SP
#        define UTIL_DLL QLIB_DLL_EXPORT
#        define UTIL_DLL_SP QLIB_DLL_EXPORT_SP
#    endif
#    ifdef QLIB_RISKMGR
#        undef RISKMGR_DLL
#        undef RISKMGR_DLL_SP
#        define RISKMGR_DLL QLIB_DLL_EXPORT
#        define RISKMGR_DLL_SP QLIB_DLL_EXPORT_SP
#    endif
#    ifdef QLIB_MARKET
#        undef MARKET_DLL
#        undef MARKET_DLL_SP
#        define MARKET_DLL QLIB_DLL_EXPORT
#        define MARKET_DLL_SP QLIB_DLL_EXPORT_SP
#    endif
#    ifdef QLIB_CREDIT
#        undef CREDIT_DLL
#        undef CREDIT_DLL_SP
#        define CREDIT_DLL QLIB_DLL_EXPORT
#        define CREDIT_DLL_SP QLIB_DLL_EXPORT_SP
#    endif
#    ifdef QLIB_FOURIER
#        undef FOURIER_DLL
#        undef FOURIER_DLL_SP
#        define FOURIER_DLL QLIB_DLL_EXPORT
#        define FOURIER_DLL_SP QLIB_DLL_EXPORT_SP
#    endif
#    ifdef QLIB_RATESLIB
#        undef RATESLIB_DLL
#        undef RATESLIB_DLL_SP
#        define RATESLIB_DLL QLIB_DLL_EXPORT
#        define RATESLIB_DLL_SP QLIB_DLL_EXPORT_SP
#    endif
#    ifdef QLIB_TREE
#        undef TREE_DLL
#        undef TREE_DLL_SP
#        define TREE_DLL QLIB_DLL_EXPORT
#        define TREE_DLL_SP QLIB_DLL_EXPORT_SP
#    endif
#    ifdef QLIB_MCARLO
#        undef MCARLO_DLL
#        undef MCARLO_DLL_SP
#        define MCARLO_DLL QLIB_DLL_EXPORT
#        define MCARLO_DLL_SP QLIB_DLL_EXPORT_SP
#    endif
#    ifdef QLIB_CONVOLUTION
#        undef CONVOLUTION_DLL
#        undef CONVOLUTION_DLL_SP
#        define CONVOLUTION_DLL QLIB_DLL_EXPORT
#        define CONVOLUTION_DLL_SP QLIB_DLL_EXPORT_SP
#    endif
#    ifdef QLIB_PRODUCTS
#        undef PRODUCTS_DLL
#        undef PRODUCTS_DLL_SP
#        define PRODUCTS_DLL QLIB_DLL_EXPORT
#        define PRODUCTS_DLL_SP QLIB_DLL_EXPORT_SP
#    endif
#    ifdef QLIB_IRPRODUCTS
#        undef IRPRODUCTS_DLL
#        undef IRPRODUCTS_DLL_SP
#        define IRPRODUCTS_DLL QLIB_DLL_EXPORT
#        define IRPRODUCTS_DLL_SP QLIB_DLL_EXPORT_SP
#    endif
#    ifdef QLIB_PYPRODUCTS
#        undef PYPRODUCTS_DLL
#        undef PYPRODUCTS_DLL_SP
#        define PYPRODUCTS_DLL QLIB_DLL_EXPORT
#        define PYPRODUCTS_DLL_SP QLIB_DLL_EXPORT_SP
#    endif
#    ifdef QLIB_ADDINS
#        undef ADDINS_DLL
#        undef ADDINS_DLL_SP
#        define ADDINS_DLL QLIB_DLL_EXPORT
#        define ADDINS_DLL_SP QLIB_DLL_EXPORT_SP
#    endif
#    ifdef QLIB_REGTEST
#        undef REGTEST_DLL
#        undef REGTEST_DLL_SP
#        define REGTEST_DLL QLIB_DLL_EXPORT
#        define REGTEST_DLL_SP QLIB_DLL_EXPORT_SP
#    endif
#    ifdef QLIB_RADAR
#        undef RADAR_DLL
#        undef RADAR_DLL_SP
#        define RADAR_DLL QLIB_DLL_EXPORT
#        define RADAR_DLL_SP QLIB_DLL_EXPORT_SP
#    endif


// leave it here until found a better place
#    ifdef ASSERT_OFF
#        define ASSERT(condition) (0)
#    else

#        define ASSERT(e) ((e) ? (void)0 : throw ModelException::fromAssert(__FILE__, __LINE__, #e))
#    endif

#    ifdef __GNUC__
#        if __GNUC__ < 3
// extensions are in standard place and in std space
#            define ext_hash_map <hash_map>
#            define ext_hash_set <hash_set>
#        else
#            define ext_hash_map <ext/hash_map>
#            define ext_hash_set <ext/hash_set>
#            if __GNUC__ == 3 && __GNUC_MINOR__ == 0
// extensions are in ext directory and in std space
#            else
// extensions are in ext directory and in __gnu_cxx space
// include at least one file to be able to use __gnu_cxx" namespace;
// FIXME: this part smells bad, we probably should delete next 2 lines
#			 include <ext/hash_set>
using namespace __gnu_cxx;
#            endif
#        endif
#    else
// STL port - in  standard place and in std space
#        define ext_hash_map <hash_map>
#        define ext_hash_set <hash_set>
#    endif

/** Define how to 'extern' templates. This is not ansi C++ but is supported by
    the compilers we currently use. We wrap it into a macro so we can switch
    it off (eg for non supporting compiler or more importantly for optimised
    build since, at least with MSVC, extern template overrides any inline 
    specification. */
#    ifdef DEBUG
#        if defined (_MSC_VER)
// switch off warning about non standard extension
#            pragma warning( disable : 4231 )
#            pragma warning( disable : 4275 )
					/* exported classes deriving from non-exported
					   classes - affects SPs */
// stlport already does sorts out strings for you when building dlls
#            if !defined(QLIB_BUILD_DLL)
// temporarily switch off warning about template already instantiated
#                pragma warning( disable : 4660 )
//extern template basic_string < char, char_traits < char >, allocator < char > >;
#                pragma warning( default : 4660 )
#            endif

#        endif
	   // _MSC_VER
#        define EXTERN_TEMPLATE(t) extern template t
#        define INSTANTIATE_TEMPLATE(t) template t
// don't really want anything inline for debug
// this is for templates
#        define TEMPLATE_INLINE
// this is for regular classes/methods etc
#        define OPT_INLINE
#    elif defined(QLIB_BUILD_DLL) && defined (__GNUC__) && (__GNUC__ >= 4)

#        define EXTERN_TEMPLATE(t) extern template t
#        define INSTANTIATE_TEMPLATE(t)  template QLIB_DLL_EXPORT t
#        define TEMPLATE_INLINE
#        define OPT_INLINE

#    else
// expand both to nothing
#        define EXTERN_TEMPLATE(t) 
#        define INSTANTIATE_TEMPLATE(t)
// this is for templates
#        define TEMPLATE_INLINE inline
// this is for regular classes/methods etc
#        define OPT_INLINE inline
#    endif
// handy for use in EXTERN_TEMPLATE with templated classes
#    define _COMMA_ ,

#    if defined _MSC_VER
#        define NORETURN __declspec(noreturn)
#    elif defined __GNUC__
#        define NORETURN __attribute__((noreturn))
#    endif
// introduce namespace so we can start using it
CORE_BEGIN_NAMESPACE CORE_END_NAMESPACE

/* These are here since they allow the use of precompiled headers on MSVC.
Unfortunately all compilers must see them. Note that they must come after
the various #defines above for _STLP_XXXX */
#    include <math.h>
#    include <cstddef>
#    include <cstdarg>
#    include <cstdio>
#    include <typeinfo>
#    include <cstdlib>
#    include <iostream>
#    include <vector>
#    include <exception>
#    include <memory>
#    include <string>
#include <boost/shared_ptr.hpp>
#include <boost/type_traits/is_enum.hpp>

using namespace std;
using boost::shared_ptr;


#endif
