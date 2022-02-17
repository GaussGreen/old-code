

#ifndef	ICM_F1_FUNCTORSDSYS_H_
#define ICM_F1_FUNCTORSDSYS_H_

#undef	TEMPLATE_F1 
#define TEMPLATE_F1

#include <functional>

namespace ff1 
{
	//	Here are a collection of classes and functions for generating functors.
	//
	//	I did this for I encountered a bug in the MSVC6 "mem_fun" collection. 
	//
	//	In their simpler forms, the binders are :
	//
	//		mem_call		: embodies the call of a member function.
	//		const_mem_call	: embodies the call of a const member function
	//		fun_call		: embodies the call of a function. 
	//	
	//		mem_call (  [&] Class::Member , ClassInstance )
	//
	//	The classes generated are functors, that is that thay have an operator() defined
	//	that takes a single argument and return the result (result type)
	//
	//	The binders store their argunemts by reference (and then allows working with
	//	polymorphic arguments). On the other hand it might be perilous to store such 
	//	a result value. 
	//
	//	We also provide a "virtual form" of the functions : 
	//
	//		mem_call_virtual 
	//		const_mem_call_virtual 
	//		fun_call_virtual 
	//	
	//	The resulting classes are all derived from the abstract class fun_virtual_t
	//
	

	//		_R		result type
	//		_Ty		class type
	//		_A		argument type

//-----------------------------------------------------------------------------------------------------------------------------------------------------
template <class _R, class _A>
class fun_virtual_t : public std::unary_function<_A,_R> 
{ 
// FIXMEFRED: mig.vc8 (28/05/2007 09:43:08):typename
#if _MSC_VER >= 1310	// Visual C++ 2005 & .net 2003
	typedef typename fun_virtual_t<_R,_A>::result_type result_type ;
#else
	typedef fun_virtual_t<_R,_A>::result_type result_type ;
#endif
	public: virtual result_type operator() (_A) = 0;
	public: virtual ~fun_virtual_t() {}
}; 
//-----------------------------------------------------------------------------------------------------------------------------------------------------
//	mem_call_t
template <class _R, class _Ty, class _A> 
class mem_call_t : public std::unary_function<_A,_R>
{
public:
	virtual ~mem_call_t(){}

	mem_call_t(_R (_Ty::* memberPtr)(_A), _Ty& objectRef):
		m_MemberPtr(memberPtr),
		m_ObjectRef(objectRef) {}

	mem_call_t(void (_Ty::* memberPtr)(_A, _R),_Ty& objectRef):
		m_MemberPtrByRef(memberPtr),
		m_ObjectRef(objectRef) {}

	result_type operator() (_A argument)			{ return (m_ObjectRef.*m_MemberPtr)(argument) ; }
	void		operator() (_A argument, _R retour) { (m_ObjectRef.*m_MemberPtrByRef)(argument, retour) ; }

private:
	_R		(_Ty::* m_MemberPtr)		(_A) ;
	void	(_Ty::* m_MemberPtrByRef)	(_A, _R) ;
	_Ty&			m_ObjectRef ; 
} ;
//-----------------------------------------------------------------------------------------------------------------------------------------------------
//	mem_call
template <class _R,class _Ty, class _A>
inline mem_call_t<_R,_Ty,_A>
mem_call(_R (_Ty::* memberPtr)(_A),_Ty& objectRef) 
{
	return mem_call_t<_R,_Ty,_A>(memberPtr,objectRef) ;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------
//	mem_call result as argument
template <class _R,class _Ty, class _A>
inline mem_call_t<_R,_Ty,_A>
mem_call(void (_Ty::* memberPtr)(_A, _R),_Ty& objectRef) 
{
	return mem_call_t<_R,_Ty,_A>(memberPtr,objectRef) ;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------
//	mem_call_virtual_t
template <class _R, class _Ty, class _A> 
class mem_call_virtual_t : public fun_virtual_t<_R,_A> 
{
public:
	virtual ~mem_call_virtual_t(){}

	mem_call_virtual_t(_R (_Ty::* memberPtr)(_A),_Ty& objectRef):
		m_MemberPtr(memberPtr),
		m_ObjectRef(objectRef) {}

	mem_call_virtual_t(void (_Ty::* memberPtr)(_A, _R),_Ty& objectRef):
		m_MemberPtrByRef(memberPtr),
		m_ObjectRef(objectRef) {}

	virtual result_type operator() (_A argument)			{ return (m_ObjectRef.*m_MemberPtr)(argument) ; }
	virtual void		operator() (_A argument, _R retour) { (m_ObjectRef.*m_MemberPtrByRef)(argument, retour) ; }

private:
	_R		(_Ty::* m_MemberPtr)		(_A) ;
	void	(_Ty::* m_MemberPtrByRef)	(_A, _R) ;
	_Ty&			m_ObjectRef ; 
}; 
//-----------------------------------------------------------------------------------------------------------------------------------------------------
//	mem_call_virtual 
template <class _R,class _Ty, class _A>
inline mem_call_virtual_t<_R,_Ty,_A> mem_call_virtual(_R (_Ty::* memberPtr)(_A),_Ty& objectRef) 
{
	return mem_call_virtual_t<_R,_Ty,_A>(memberPtr,objectRef) ;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------
//	mem_call_virtual result as argument 
template <class _R,class _Ty, class _A>
inline mem_call_virtual_t<_R,_Ty,_A> mem_call_virtual(void (_Ty::* memberPtr)(_A, _R),_Ty& objectRef) 
{
	return mem_call_virtual_t<_R,_Ty,_A>(memberPtr,objectRef) ;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------
//	const_mem_call_t
template <class _R, class _Ty, class _A> 
class const_mem_call_t : public std::unary_function<_A,_R>
{
public:
	virtual ~const_mem_call_t(){}

	const_mem_call_t(_R (_Ty::* constMemberPtr)(_A) const, const _Ty& constObjectRef):
		m_ConstMemberPtr(constMemberPtr),
		m_ConstObjectRef(constObjectRef) {}

	const_mem_call_t(void (_Ty::* constMemberPtr)(_A, _R) const, const _Ty& constObjectRef):
		m_ConstMemberPtrByRef(constMemberPtr),
		m_ConstObjectRef(constObjectRef) {}

	result_type operator() (_A argument)			{ return (m_ConstObjectRef.*m_ConstMemberPtr)(argument) ; }
	void		operator() (_A argument, _R retour) { (m_ObjectRef.*m_MemberPtrByRef)(argument, retour) ; }

private:
	_R		(_Ty::* m_ConstMemberPtr)		(_A)		const ;
	void	(_Ty::* m_ConstMemberPtrByRef)	(_A, _R)	const;
	const _Ty&		m_ConstObjectRef; 
} ;
//-----------------------------------------------------------------------------------------------------------------------------------------------------
//	const_mem_call
template <class _R,class _Ty, class _A> 
inline const_mem_call_t<_R,_Ty,_A>
const_mem_call(_R (_Ty::* constMemberPtr)(_A) const, const _Ty& constObjectRef) 
{
	return const_mem_call_t<_R,_Ty,_A>(constMemberPtr, constObjectRef) ;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------
//	const_mem_call result as argument 
template <class _R,class _Ty, class _A> 
inline const_mem_call_t<_R,_Ty,_A>
const_mem_call(void (_Ty::* constMemberPtr)(_A, _R) const, const _Ty& constObjectRef) 
{
	return const_mem_call_t<_R,_Ty,_A>(constMemberPtr, constObjectRef) ;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------
//	const_mem_call_virtual_t
template <class _R, class _Ty, class _A> 
class const_mem_call_virtual_t : public fun_virtual_t<_A,_R>
{
public:
	virtual ~const_mem_call_virtual_t(){}

	const_mem_call_virtual_t(_R (_Ty::* constMemberPtr)(_A) const, const _Ty& constObjectRef):
		m_ConstMemberPtr(constMemberPtr),
		m_ConstObjectRef(constObjectRef) {}

	const_mem_call_virtual_t(void (_Ty::* constMemberPtr)(_A, _R) const, const _Ty& constObjectRef):
		m_ConstMemberPtrByRef(constMemberPtr),
		m_ConstObjectRef(constObjectRef) {}

	virtual result_type operator() (_A argument) { return (m_ConstObjectRef.*m_ConstMemberPtr)(argument) ; }
	virtual void		operator() (_A argument, _R retour) { (m_ConstObjectRef.*m_ConstMemberPtrByRef)(argument, retour) ; }
	
private:
	_R		(_Ty::* m_ConstMemberPtr)		(_A)		const ; 
	void	(_Ty::* m_ConstMemberPtrByRef)	(_A, _R)	const;
	const _Ty& m_ConstObjectRef; 
} ;
//-----------------------------------------------------------------------------------------------------------------------------------------------------
//	const_mem_call_virtual
template <class _R,class _Ty, class _A>
inline const_mem_call_virtual_t<_R,_Ty,_A> const_mem_call_virtual(_R (_Ty::* constMemberPtr)(_A)const, const _Ty& constObjectRef) 
{
	return const_mem_call_virtual_t<_R,_Ty,_A>(constMemberPtr, constObjectRef) ;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------
//	const_mem_call_virtual result by argument
template <class _R,class _Ty, class _A>
inline const_mem_call_virtual_t<_R,_Ty,_A> const_mem_call_virtual(void (_Ty::* constMemberPtr)(_A, _R) const, const _Ty& constObjectRef) 
{
	return const_mem_call_virtual_t<_R,_Ty,_A>(constMemberPtr, constObjectRef) ;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------
//	fun_call_t
template <class _R,class _A> 
class fun_call_t : public std::unary_function<_A,_R> 
{
public:
	virtual ~fun_call_t(){}

	fun_call_t(_R (* functionPtr)(_A) ):
	  m_FunctionPtr(functionPtr) {}

	fun_call_t(void (* functionPtr)(_A, _R) ):
	  m_FunctionPtrByref(functionPtr) {}

	result_type operator() (_A argument)			{ return (*m_FunctionPtr)(argument) ;}
	void		operator() (_A argument, _R retour) { (*m_FunctionPtrByref)(argument, retour) ;}


private:
	_R		(* m_FunctionPtr)		(_A) ; 
	void	(* m_FunctionPtrByref)	(_A, _R) ; 
}; 
//-----------------------------------------------------------------------------------------------------------------------------------------------------
//	fun_call
template <class _R,class _A>
inline fun_call_t<_R,_A> fun_call(_R (* functionPtr)(_A)) 
{
	return fun_call_t<_R,_A>(functionPtr) ; 
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------
//	fun_call result by argument
template <class _R,class _A>
inline fun_call_t<_R,_A> fun_call(void (* functionPtr)(_A, _R)) 
{
	return fun_call_t<_R,_A>(functionPtr) ; 
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------
//	fun_call_virtual_t
template <class _R,class _A> 
class fun_call_virtual_t : public fun_virtual_t<_R,_A>
{
public:
	virtual ~fun_call_virtual_t(){}

	fun_call_virtual_t(_R (* functionPtr)(_A) ):
		m_FunctionPtr(functionPtr) {}

	fun_call_virtual_t(void (* functionPtr)(_A, _R) ):
		m_FunctionPtrByRef(functionPtr) {}

	virtual result_type operator() (_A argument)			{ return (*m_FunctionPtr)(argument) ;}
	virtual void		operator() (_A argument, _R retour) { return (*m_FunctionPtrByRef)(argument, retour) ;}

private:
	_R		(* m_FunctionPtr)		(_A) ; 
	void	(* m_FunctionPtrByRef)	(_A, _R) ;
};
//-----------------------------------------------------------------------------------------------------------------------------------------------------
//	fun_call_virtual
template <class _R,class _A>
inline fun_call_virtual_t<_R,_A> fun_call_virtual(_R (* functionPtr)(_A)) 
{
	return fun_call_virtual_t<_R,_A>(functionPtr) ; 
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------
//	fun_call_virtual result by argument
template <class _R,class _A>
inline fun_call_virtual_t<_R,_A> fun_call_virtual(void (* functionPtr)(_A, _R)) 
{
	return fun_call_virtual_t<_R,_A>(functionPtr) ; 
}

} // namespace f1_0 

#endif //	_F1_FUNCTORS_H_