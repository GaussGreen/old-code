/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file port.h
 *
 *  \brief portable issues
 *	\author  Eric Benhamou
 *	\version 1.0
 *	\date September 2003
 */

/*----------------------------------------------------------------------------*/

#ifndef _IN_GPABASE_PORT_H
#define _IN_GPABASE_PORT_H


#ifdef __cplusplus

#if !defined(__SUNPRO_CC) && !defined(_MSC_VER)
	#define __STANDARD_PORTABLE_MACROS
#endif


//////////////////////////////////////////////////
/// for old sun pro C++ 4.x, we use no namespace
//////////////////////////////////////////////////
#if ( defined(__SUNPRO_CC ) && __SUNPRO_CC >= 0x001 && __SUNPRO_CC < 0x500 )
	/// do not remove space! these are fundamental
	#define CC_BEGIN_NAMESPACE(name)     
	#define CC_END_NAMESPACE() 
	#define CC_NS(nameSpace,name) name   
	#define CC_USING(something)         
	#define CC_USING_NS(nameSpace,name)      
	#define CC_USING_NAMESPACE(nameSpace)  
	#define CC_STL_VECTOR( arg )  
	#define CC_USING_T(something,TEMPLATE)         
	#define CC_USING_NS_T(nameSpace,name,TEMPLATE)
	#define CC_USING_NS_BI_T( nameSpace, name, TEMPLATE1, TEMPLATE2 )
	#define CC_USING_BI_T( name, TEMPLATE1, TEMPLATE2 )
	#define CC_USING_NS_TRI_T( nameSpace, name, TEMPLATE1, TEMPLATE2, TEMPLATE3 )
	#define CC_USING_TRI_T( name, TEMPLATE1, TEMPLATE2, TEMPLATE3 )
	#define CC_TYPENAME 
	#define CC_TEMPLATE_SPECIALIZATION
	#define CC_DISTANCE distance
	#define CC_MUTABLE( type, x )  const_cast<type*>(this)->x
	#define CC_IS_MUTABLE 
#endif



//////////////////////////////////////////////////
/// for the new sun pro C++ 6.x
//////////////////////////////////////////////////
#if ( defined(__SUNPRO_CC ) && __SUNPRO_CC >= 0x500 )
	#define CC_BEGIN_NAMESPACE( name ) namespace name {
	#define CC_END_NAMESPACE() 	}
	#define CC_NS( nameSpace, name ) nameSpace::name
	#define CC_USING( something )	using something;
	#define CC_USING_NS( nameSpace, name ) using nameSpace::name;
	#define CC_USING_NAMESPACE( nameSpace )  using namespace nameSpace;
	#define CC_STL_VECTOR( arg )  CC_NS( std, vector )< arg >
	#define CC_USING_T( something, TEMPLATE )	CC_USING( something)	
	#define CC_USING_NS_T( nameSpace, name, TEMPLATE ) CC_USING_NS( nameSpace, name ) 
	#define CC_USING_NS_BI_T( nameSpace, name, TEMPLATE1, TEMPLATE2  ) CC_USING_NS( nameSpace, name ) 
	#define CC_USING_BI_T( name, TEMPLATE1, TEMPLATE2 )
	#define CC_USING_NS_TRI_T( nameSpace, name, TEMPLATE1, TEMPLATE2, TEMPLATE3 )  CC_USING_NS( nameSpace, name )
	#define CC_USING_TRI_T( name, TEMPLATE1, TEMPLATE2, TEMPLATE3 )
	#define CC_TYPENAME typename
	#define CC_TEMPLATE_SPECIALIZATION template<>

/// small distance implementation
namespace std {
	template <class Iterator>
		inline long myDistance( Iterator pos1, Iterator pos2 )
	{
		long d=0;
		distance(pos1,pos2,d);
		return d;
	}
}
	#define CC_DISTANCE myDistance
	#define CC_MUTABLE( type, x )  x
	#define CC_IS_MUTABLE mutable
#endif


////////////////////////////////////////
/// Microsoft studio 6.0 platform
////////////////////////////////////////

#if	(defined(_MSC_VER) && _MSC_VER <= 1200 )
	#define CC_BEGIN_NAMESPACE(name) namespace name {
	#define CC_END_NAMESPACE() 	}
	#define CC_NS(nameSpace,name) nameSpace::name
	#define CC_USING( something )	using something;
	#define CC_USING_NS( nameSpace, name ) using nameSpace::name;
	#define CC_USING_NAMESPACE( nameSpace )  using namespace nameSpace;
	#define CC_STL_VECTOR( arg )  CC_NS( std, vector )< arg >
	#define CC_USING_T( something, TEMPLATE )	using something< TEMPLATE >;
	#define CC_USING_NS_T( nameSpace, name, TEMPLATE ) using nameSpace::name< TEMPLATE >;
	#define CC_USING_NS_BI_T( nameSpace, name, TEMPLATE1, TEMPLATE2  ) using nameSpace::name< TEMPLATE1, TEMPLATE2 >;
	#define CC_USING_BI_T( name, TEMPLATE1, TEMPLATE2 ) name< TEMPLATE1, TEMPLATE2 >
	#define CC_USING_NS_TRI_T( nameSpace, name, TEMPLATE1, TEMPLATE2, TEMPLATE3 ) using nameSpace::name< TEMPLATE1, TEMPLATE2, TEMPLATE3 >;
	#define CC_USING_TRI_T( name, TEMPLATE1, TEMPLATE2, TEMPLATE3 ) name< TEMPLATE1, TEMPLATE2, TEMPLATE3 >
	#define CC_TYPENAME 
	#define CC_TEMPLATE_SPECIALIZATION template<>
	#define CC_DISTANCE distance
	#define CC_MUTABLE( type, x )  x
	#define CC_IS_MUTABLE mutable
#endif

////////////////////////////////////////
/// Microsoft studio .Net platform
////////////////////////////////////////
#if	( defined(_MSC_VER) && _MSC_VER >= 1300 )
	// the standard
	#define CC_BEGIN_NAMESPACE( name ) namespace name {
	#define CC_END_NAMESPACE() 	}
	#define CC_NS( nameSpace, name ) nameSpace::name
	#define CC_USING( something )	using something;
	#define CC_USING_NS( nameSpace, name ) using nameSpace::name;
	#define CC_USING_NAMESPACE( nameSpace )  using namespace nameSpace;
	#define CC_STL_VECTOR( arg )  CC_NS( std, vector )< arg >
	#define CC_USING_T( something, TEMPLATE )	CC_USING( something)	
	#define CC_USING_NS_T( nameSpace, name, TEMPLATE ) CC_USING_NS( nameSpace, name ) 
	#define CC_USING_NS_BI_T( nameSpace, name, TEMPLATE1, TEMPLATE2  ) CC_USING_NS( nameSpace, name ) 
	#define CC_USING_BI_T( name, TEMPLATE1, TEMPLATE2 )
	#define CC_USING_NS_TRI_T( nameSpace, name, TEMPLATE1, TEMPLATE2, TEMPLATE3 )  CC_USING_NS( nameSpace, name )
	#define CC_USING_TRI_T( name, TEMPLATE1, TEMPLATE2, TEMPLATE3 )
	#define CC_TYPENAME typename
	#define CC_TEMPLATE_SPECIALIZATION template<>
	#define CC_DISTANCE distance
	#define CC_MUTABLE( type, x )  x
	#define CC_IS_MUTABLE mutable
#endif


/// this is the default!
#if defined(__STANDARD_PORTABLE_MACROS)
	// the standard
	#define CC_BEGIN_NAMESPACE( name ) namespace name {
	#define CC_END_NAMESPACE()	}
	#define CC_NS( nameSpace, name ) nameSpace::name
	#define CC_USING( something )	using something;
	#define CC_USING_NS( nameSpace, name ) using nameSpace::name;
	#define CC_USING_NAMESPACE( nameSpace )  using namespace nameSpace;
	#define CC_STL_VECTOR( arg )  CC_NS( std, vector )< arg >
	#define CC_USING_T( something, TEMPLATE )	CC_USING( something)	
	#define CC_USING_NS_T( nameSpace, name, TEMPLATE ) CC_USING_NS( nameSpace, name ) 
	#define CC_USING_NS_BI_T( nameSpace, name, TEMPLATE1, TEMPLATE2  ) CC_USING_NS( nameSpace, name ) 
	#define CC_USING_BI_T( name, TEMPLATE1, TEMPLATE2 )
	#define CC_USING_NS_TRI_T( nameSpace, name, TEMPLATE1, TEMPLATE2, TEMPLATE3 )  CC_USING_NS( nameSpace, name )
	#define CC_USING_TRI_T( name, TEMPLATE1, TEMPLATE2, TEMPLATE3 )
	#define CC_TYPENAME typename
	#define CC_TEMPLATE_SPECIALIZATION template<>
	#define CC_DISTANCE distance
	#define CC_MUTABLE( type, x )  x
	#define CC_IS_MUTABLE mutable
#endif

/// for view method in unix world!
#ifdef unix
	#include <sys/types.h>
	#include <unistd.h>
#endif


#endif /*---- __cplusplus ----*/


#endif 

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/
