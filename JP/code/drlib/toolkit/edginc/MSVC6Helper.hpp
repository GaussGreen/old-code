//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : MSVC6Helper.hpp
//
//   Description : Microsoft Visual C++ Version 6.0 or lower helper class
//
//   Date        : Oct 2004
//
//----------------------------------------------------------------------------

#ifndef MSVC6HELPER__HPP
#define MSVC6HELPER__HPP

#include "edginc/config.hpp"

DRLIB_BEGIN_NAMESPACE

#if defined (_MSC_VER) && !(_MSC_VER >= 1300)

/** List of helper methods and structures to work around some of msvc flaws.
    Copied from Loki VC 6.0 Port by Hume@c-plusplus.de */
struct MSVC6Helper{
//private:
////////////////////////////////////////////////////////////////////////////////
// class template AlwaysFalse
// Invocation: AlwaysFalse<T>::value
// value will always by 0 (false)
////////////////////////////////////////////////////////////////////////////////
	template< typename T >
	struct AlwaysFalse
	{
		enum { value = false };
	};

////////////////////////////////////////////////////////////////////////////////
// class template ApplyImpl1
// Invocation: ApplyImpl1<T>::template Result<T1>
// T must be a nontemplate type with a nested class template named In.
// The class template is a helper for the Apply1-Template
////////////////////////////////////////////////////////////////////////////////
	template <class TypeWithNestedTemplate>
	struct ApplyImpl1
	{
		template<bool flag>
		struct VC_WORKAROUND : public TypeWithNestedTemplate {};

        struct VC_WORKAROUND<true>
		{	template<class> struct In; };

		template< typename T1 > struct Result : public
		VC_WORKAROUND< AlwaysFalse<TypeWithNestedTemplate>::value >::template In<T1>
		{
			typedef VC_WORKAROUND< AlwaysFalse<TypeWithNestedTemplate>::value >::template In<T1> Base;
		};
	};
////////////////////////////////////////////////////////////////////////////////
// class template ApplyImpl2
// Invocation: ApplyImpl2<T>::template Result<T1, T2>
// T must be a nontemplate type with a nested class template named In.
// The class template is a helper for the Apply2-Template
////////////////////////////////////////////////////////////////////////////////
	template <class TypeWithNestedTemplate>
	struct ApplyImpl2
	{
		template<bool flag>
		struct VC_WORKAROUND : public TypeWithNestedTemplate {};

		struct VC_WORKAROUND<true>
		{template<class T, class U> struct In; };

		template< typename T1, typename T2 > struct Result : public
		VC_WORKAROUND< AlwaysFalse<TypeWithNestedTemplate>::value >::template In<T1, T2>
		{
			typedef VC_WORKAROUND< AlwaysFalse<TypeWithNestedTemplate>::value >::template In<T1, T2> Base;
		};
	};
////////////////////////////////////////////////////////////////////////////////
// class template ApplyImpl3
// Invocation: ApplyImpl3<T>::template Result<T1, T2, T3>
// T must be a nontemplate type with a nested class template named In.
// The class template is a helper for the Apply3-Template
////////////////////////////////////////////////////////////////////////////////
	template <class TypeWithNestedTemplate>
	struct ApplyImpl3 
	{
		template<bool flag>
		struct VC_WORKAROUND : public TypeWithNestedTemplate {};

		struct VC_WORKAROUND<true>
		{template<class T, class U, class V> struct In; };

		template< typename T1, typename T2, typename T3 > struct Result : public
		VC_WORKAROUND< AlwaysFalse<TypeWithNestedTemplate>::value >::template In<T1, T2, T3>
		{
			typedef VC_WORKAROUND< AlwaysFalse<TypeWithNestedTemplate>::value >::template In<T1, T2, T3> Base;
		};
	};
////////////////////////////////////////////////////////////////////////////////
// class template ApplyImpl4
// Invocation: ApplyImpl4<T>::template Result<T1, T2, T3, T4>
// T must be a nontemplate type with a nested class template named In.
// The class template is a helper for the Apply4-Template
////////////////////////////////////////////////////////////////////////////////
	template <class TypeWithNestedTemplate>
	struct ApplyImpl4 
	{
		template<bool flag>
		struct VC_WORKAROUND : public TypeWithNestedTemplate {};

		struct VC_WORKAROUND<true>
		{template<class T, class U, class V, class W> struct In; };

		template< typename T1, typename T2, typename T3, typename T4 > struct Result : public
		VC_WORKAROUND< AlwaysFalse<TypeWithNestedTemplate>::value >::template In<T1, T2, T3, T4>
		{
			typedef VC_WORKAROUND< AlwaysFalse<TypeWithNestedTemplate>::value >::template In<T1, T2, T3, T4> Base;
		};
	};

public:

////////////////////////////////////////////////////////////////////////////////
// class template Apply1
// Invocation: Apply1<T, U>
// Applies the type U to the inner template In of type T
// The class template Apply1 helps to emulate template template parameters
// i first saw this technique in boost's mpl library.
////////////////////////////////////////////////////////////////////////////////
	template<typename F, typename T1>
	struct Apply1 : ApplyImpl1<F>::template Result<T1>
	{
		typedef typename ApplyImpl1<F>::template Result<T1>::Base Base;
	};
////////////////////////////////////////////////////////////////////////////////
// class template Apply2
// Invocation: Apply2<T, U, V>
// Applies the types U and V to the inner template In of type T
// The class template Apply2 helps to emulate template template parameters
// i first saw this technique in boost's mpl library.
////////////////////////////////////////////////////////////////////////////////
	template<typename F, typename T1, typename T2>
	struct Apply2 : ApplyImpl2<F>::template Result<T1, T2> 
	{
		typedef typename ApplyImpl2<F>::template Result<T1, T2>::Base Base;
	};
////////////////////////////////////////////////////////////////////////////////
// class template Apply3
// Invocation: Apply3<T, U, V, W>
// Applies the types U, V and W to the inner template In of type T
// The class template Apply3 helps to emulate template template parameters
// i first saw this technique in boost's mpl library.
////////////////////////////////////////////////////////////////////////////////
	template<typename F, typename T1, typename T2, typename T3>
	struct Apply3 : ApplyImpl3<F>::template Result<T1, T2, T3>
	{
		typedef typename ApplyImpl3<F>::template Result<T1, T2, T3>::Base Base;
	};
////////////////////////////////////////////////////////////////////////////////
// class template Apply4
// Invocation: Apply4<T, U, V, W, X>
// Applies the types U, V, W and X to the inner template In of type T
// The class template Apply4 helps to emulate template template parameters
// i first saw this technique in boost's mpl library.
////////////////////////////////////////////////////////////////////////////////
	template<typename F, typename T1, typename T2, typename T3, typename T4>
	struct Apply4 : ApplyImpl4<F>::template Result<T1, T2, T3, T4>
	{
		typedef typename ApplyImpl4<F>::template Result<T1, T2, T3, T4>::Base Base;
	};
};

#else   // looks like gcc requires some more work around VC's workarounds (arghhhh!!!)

/** List of helper methods and structures to work around some of msvc flaws.
    Copied from Loki VC 6.0 Port by Hume@c-plusplus.de */
struct MSVC6Helper{
//private:
////////////////////////////////////////////////////////////////////////////////
// class template AlwaysFalse
// Invocation: AlwaysFalse<T>::value
// value will always by 0 (false)
////////////////////////////////////////////////////////////////////////////////
	template< typename T >
	struct AlwaysFalse
	{
		enum { value = false };
	};

////////////////////////////////////////////////////////////////////////////////
// class template ApplyImpl1
// Invocation: ApplyImpl1<T>::template Result<T1>
// T must be a nontemplate type with a nested class template named In.
// The class template is a helper for the Apply1-Template
////////////////////////////////////////////////////////////////////////////////
	// looks like gcc doesnt support specialization of nested templates
	template<class TypeWithNestedTemplate, bool flag>
	struct ApplyImpl1_VC_WORKAROUND : public TypeWithNestedTemplate {};

    template<class TypeWithNestedTemplate>
    struct ApplyImpl1_VC_WORKAROUND<TypeWithNestedTemplate, true>
	{	template<class> struct In; };

    template <class TypeWithNestedTemplate>
	struct ApplyImpl1
	{
        template< typename T1 > struct Result : public
		ApplyImpl1_VC_WORKAROUND< TypeWithNestedTemplate, AlwaysFalse<TypeWithNestedTemplate>::value >::template In<T1>
		{
            typedef typename ApplyImpl1_VC_WORKAROUND< TypeWithNestedTemplate, AlwaysFalse<TypeWithNestedTemplate>::value >::template In<T1> Base;
		};
    };
////////////////////////////////////////////////////////////////////////////////
// class template ApplyImpl2
// Invocation: ApplyImpl2<T>::template Result<T1, T2>
// T must be a nontemplate type with a nested class template named In.
// The class template is a helper for the Apply2-Template
////////////////////////////////////////////////////////////////////////////////
	// looks like gcc doesnt support specialization of nested templates
	template <class TypeWithNestedTemplate, bool flag>
	struct ApplyImpl2_VC_WORKAROUND : public TypeWithNestedTemplate {};

	template <class TypeWithNestedTemplate>
    struct ApplyImpl2_VC_WORKAROUND<TypeWithNestedTemplate, true>
    {   template<class, class> struct In; };

	template <class TypeWithNestedTemplate>
	struct ApplyImpl2
	{
		template< typename T1, typename T2 > struct Result : public
		ApplyImpl2_VC_WORKAROUND< TypeWithNestedTemplate, AlwaysFalse<TypeWithNestedTemplate>::value >::template In<T1, T2>
		{
			typedef typename ApplyImpl2_VC_WORKAROUND< TypeWithNestedTemplate, AlwaysFalse<TypeWithNestedTemplate>::value >::template In<T1, T2> Base;
		};
	};
////////////////////////////////////////////////////////////////////////////////
// class template ApplyImpl3
// Invocation: ApplyImpl3<T>::template Result<T1, T2, T3>
// T must be a nontemplate type with a nested class template named In.
// The class template is a helper for the Apply3-Template
////////////////////////////////////////////////////////////////////////////////
	template<class TypeWithNestedTemplate, bool flag>
	struct ApplyImpl3_VC_WORKAROUND : public TypeWithNestedTemplate {};

	template<class TypeWithNestedTemplate>
	struct ApplyImpl3_VC_WORKAROUND<TypeWithNestedTemplate, true>
	{template<class T, class U, class V> struct In; };

	template <class TypeWithNestedTemplate>
	struct ApplyImpl3 
	{
		template< typename T1, typename T2, typename T3 > struct Result : public
		ApplyImpl3_VC_WORKAROUND< TypeWithNestedTemplate, AlwaysFalse<TypeWithNestedTemplate>::value >::template In<T1, T2, T3>
		{
			typedef typename ApplyImpl3_VC_WORKAROUND< TypeWithNestedTemplate, AlwaysFalse<TypeWithNestedTemplate>::value >::template In<T1, T2, T3> Base;
		};
	};
////////////////////////////////////////////////////////////////////////////////
// class template ApplyImpl4
// Invocation: ApplyImpl4<T>::template Result<T1, T2, T3, T4>
// T must be a nontemplate type with a nested class template named In.
// The class template is a helper for the Apply4-Template
////////////////////////////////////////////////////////////////////////////////
	template<class TypeWithNestedTemplate, bool flag>
	struct ApplyImpl4_VC_WORKAROUND : public TypeWithNestedTemplate {};

	template<class TypeWithNestedTemplate>
	struct ApplyImpl4_VC_WORKAROUND<TypeWithNestedTemplate, true>
	{template<class T, class U, class V, class W> struct In; };

	template <class TypeWithNestedTemplate>
	struct ApplyImpl4 
	{
		template< typename T1, typename T2, typename T3, typename T4 > struct Result : public
		ApplyImpl4_VC_WORKAROUND< TypeWithNestedTemplate, AlwaysFalse<TypeWithNestedTemplate>::value >::template In<T1, T2, T3, T4>
		{
			typedef typename ApplyImpl4_VC_WORKAROUND< TypeWithNestedTemplate, AlwaysFalse<TypeWithNestedTemplate>::value >::template In<T1, T2, T3, T4> Base;
		};
	};
public:

////////////////////////////////////////////////////////////////////////////////
// class template Apply1
// Invocation: Apply1<T, U>
// Applies the type U to the inner template In of type T
// The class template Apply1 helps to emulate template template parameters
// i first saw this technique in boost's mpl library.
////////////////////////////////////////////////////////////////////////////////
	template<typename F, typename T1>
	struct Apply1 : ApplyImpl1<F>::template Result<T1>
	{
		typedef typename ApplyImpl1<F>::template Result<T1>::Base Base;
	};
////////////////////////////////////////////////////////////////////////////////
// class template Apply2
// Invocation: Apply2<T, U, V>
// Applies the types U and V to the inner template In of type T
// The class template Apply2 helps to emulate template template parameters
// i first saw this technique in boost's mpl library.
////////////////////////////////////////////////////////////////////////////////
	template<typename F, typename T1, typename T2>
	struct Apply2 : ApplyImpl2<F>::template Result<T1, T2> 
	{
		typedef typename ApplyImpl2<F>::template Result<T1, T2>::Base Base;
	};
////////////////////////////////////////////////////////////////////////////////
// class template Apply3
// Invocation: Apply3<T, U, V, W>
// Applies the types U, V and W to the inner template In of type T
// The class template Apply3 helps to emulate template template parameters
// i first saw this technique in boost's mpl library.
////////////////////////////////////////////////////////////////////////////////
	template<typename F, typename T1, typename T2, typename T3>
	struct Apply3 : ApplyImpl3<F>::template Result<T1, T2, T3>
	{
		typedef typename ApplyImpl3<F>::template Result<T1, T2, T3>::Base Base;
	};
////////////////////////////////////////////////////////////////////////////////
// class template Apply4
// Invocation: Apply4<T, U, V, W, X>
// Applies the types U, V, W and X to the inner template In of type T
// The class template Apply4 helps to emulate template template parameters
// i first saw this technique in boost's mpl library.
////////////////////////////////////////////////////////////////////////////////
	template<typename F, typename T1, typename T2, typename T3, typename T4>
	struct Apply4 : ApplyImpl4<F>::template Result<T1, T2, T3, T4>
	{
		typedef typename ApplyImpl4<F>::template Result<T1, T2, T3, T4>::Base Base;
	};
};

#endif // (_MSC_VER) && !(_MSC_VER >= 1300)

DRLIB_END_NAMESPACE

#endif
