//	atlbase_noasset.h: We redefine ATLASSERT and include <atlbase.h>.
//					   This means that the ATLASSERT macro will write
//					   out trace output rather than display the message.
//					   This is useful in the case of the ATLASSERTS in
//					   CComPtr<> etc. where we do actually want to take
//					   the address of a non-null pointer.
//
//	author:		       David Cuin
//
//////////////////////////////////////////////////////////////////////

#ifdef ATLASSERT
	#undef ATLASSERT
#endif
#define ATLASSERT(expr)					do { if (#expr != "p==NULL") _ASSERTE(expr); } while (0)
#include <atlbase.h>
#ifdef ATLASSERT
	#undef ATLASSERT
#endif
#define ATLASSERT(expr) _ASSERTE(expr)

