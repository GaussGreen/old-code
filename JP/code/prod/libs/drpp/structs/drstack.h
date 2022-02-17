// drstack.h: interface for the drstack class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_DRSTACK_H__7FC9455C_27B6_11D2_97B5_00C04FD8EB9A__INCLUDED_)
#define AFX_DRSTACK_H__7FC9455C_27B6_11D2_97B5_00C04FD8EB9A__INCLUDED_

#if _MSC_VER >= 1000
#pragma once
#endif // _MSC_VER >= 1000

#include "drplatdep.h"
#include "drexception.h"
#include <vector>

//k  Stack class.  Exactly matches the stl spec 
//k  (ok, maybe size() should return size_t instead of int).  
//k  I wrote my own version, cause it seemed
//k  easier to write it from scratch than to get the SGI stl version to work.
//k  (Solaris's support of templates is pathetic -- no default template 
//k  arguments and templates only work for classes)

template <class T>
class drstack  
{
public:
	drstack(); //
	bool empty() const; //
	void push (const T& v); //
	void pop (); //
	const T& top () const; //
	T& top (); //
	int size() const;  //
private:
	vector < T, MYALLOC (T) > s;
};

template <class T> inline
drstack<T>::drstack() {}

template <class T> inline
bool drstack<T>::empty() const {return (s.size() == 0);}

template <class T> inline
void drstack<T>::push (const T& v) {s.push_back (v);}

template <class T> inline
void drstack<T>::pop () 
{
	if (!empty()) s.pop_back(); 
	else throw DRException ("Attempt to pop from empty stack");
}

template <class T> inline
const T& drstack<T>::top () const 
{
	if (!empty()) return s[s.size() - 1]; 
	throw DRException ("Attempt to read from empty stack"); 
	return s[0];
}

template <class T> inline
T& drstack<T>::top () 
{
	if (!empty()) return s[s.size() - 1]; 
	throw DRException ("Attempt to read from empty stack"); 
	return s[0];
}

template <class T> inline
int drstack<T>::size() const {return s.size();}


#endif // !defined(AFX_DRSTACK_H__7FC9455C_27B6_11D2_97B5_00C04FD8EB9A__INCLUDED_)
