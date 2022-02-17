#ifndef _ICM_PTR_H_ 
#define _ICM_PTR_H_ 

#include "ICMKernel/util/icm_macro.h"

//
//		aggreg_ptr<T> 	
//
//		Acts as a T* (see std::auto_ptr)
//
//		The construction is an adoption. 
//
//		On copy, the instance will be cloned (polymorphic). 
// 
//		release() will return the pointer and setting the
//		current instance to empty state. 
// 
template<class T> class aggreg_ptr 
{
private:
	T* itsPointer ;		// 0.1 , aggregated. 
public:
	explicit aggreg_ptr(T * item=0) : itsPointer(item) {}
	aggreg_ptr(const aggreg_ptr<T>&ref) : itsPointer(ref.dup()) {} 
	aggreg_ptr<T>& operator=(const aggreg_ptr<T>& ref) ; 
	~aggreg_ptr() { if (itsPointer) delete itsPointer; itsPointer=0; }
	//
	void adopt(T* item=0) { if (itsPointer) delete itsPointer; itsPointer=item; }
	T&	operator*() const 	
	{	
		if (itsPointer==0) ICMTHROW(ERR_INVALID_ARGUMENT,"aggreg_ptr<>::*"); 
		return *itsPointer; 
	}
	T*	operator->() const  
	{	
		if (itsPointer==0) ICMTHROW(ERR_INVALID_ARGUMENT,"aggreg_ptr<>::->"); 
		return itsPointer; 
	}
	T*	release()  ; 
	T*  dup() const  ;
};

//	-------------------------------------------------------------------------------------------
template <class T> 
inline 
aggreg_ptr<T>& 
aggreg_ptr<T>::operator=(const aggreg_ptr<T>&ref)
{
	if (this!=&ref)
	{
		if (itsPointer) delete itsPointer; itsPointer=0; 
		itsPointer=ref.dup(); 
	}
	return *this; 
}
//	-------------------------------------------------------------------------------------------
template<class T> 
inline 
T* 
aggreg_ptr<T>::release() 
{
	T* temp(itsPointer); 
	itsPointer=0; 
	return temp; 
} 
//	-------------------------------------------------------------------------------------------
template<class T> 
inline 
T* 
aggreg_ptr<T>::dup() const 
{
	if (itsPointer) return (T*)(itsPointer->Clone()); 
	return 0 ;
} 
//	-------------------------------------------------------------------------------------------
//	-------------------------------------------------------------------------------------------

#endif // _ICM_PTR_H_ 