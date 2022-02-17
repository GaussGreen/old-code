/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/win32/mfc.h
// Purpose:     various helpers for programming with MFC
// Created:     2006-02-01
// RCS-ID:      $Id: mfc.h,v 1.8 2006/05/29 19:59:31 zhang Exp $
// Copyright:   (c) 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/// @file ito33/win32/mfc.h

#ifndef _ITO33_WIN32_MFC_H_
#define _ITO33_WIN32_MFC_H_

#define TEMPLATE_1(t1)                 t1
#define TEMPLATE_2(t1, t2)             t1, t2
#define TEMPLATE_3(t1, t2, t3)         t1, t2, t3

#define TCLASS_1(theClass, t1)         theClass<t1>
#define TCLASS_2(theClass, t1, t2)     theClass<t1, t2>
#define TCLASS_3(theClass, t1, t2, t3) theClass<t1, t2, t3>

/**
    Begins the definition of message map for a template control.

    Examples
    @code
    BEGIN_TEMPLATE_MESSAGE_MAP( typename T, CMyEdit<T>, CEdit )

    BEGIN_TEMPLATE_MESSAGE_MAP( TEMPLATE_2(T, U), TCLASS_2(CMyCombo, T, U), \
                                CComboBox )
    @endcode
    
    When the template has a non-type parameter, the first parameter
    to BEGIN_TEMPLATE_MESSAGE_MAP should be the declaration of this
    non-type parameter as it appears in the class template declaration.
    
    The same applies to the arguments of TEMPLATE_X.
 */
#if _MSC_VER < 1400
  #ifdef _AFXDLL
  #define ITO_BEGIN_TEMPLATE_MESSAGE_MAP(theTemplate, theClass, baseClass) \
	  template <theTemplate> const AFX_MSGMAP* PASCAL theClass::GetThisMessageMap() \
		  { return &theClass::messageMap; } \
	  template <theTemplate> const AFX_MSGMAP* theClass::GetMessageMap() const \
		  { return &theClass::messageMap; } \
	  template <theTemplate> AFX_COMDAT const AFX_MSGMAP theClass::messageMap = \
	  { &baseClass::GetThisMessageMap, &theClass::_messageEntries[0] }; \
	  template <theTemplate> AFX_COMDAT const AFX_MSGMAP_ENTRY theClass::_messageEntries[] = \
	  { \

  #else
  #define ITO_BEGIN_TEMPLATE_MESSAGE_MAP(theTemplate, theClass, baseClass) \
	  template <theTemplate> const AFX_MSGMAP* theClass::GetMessageMap() const \
		  { return &theClass::messageMap; } \
	  template <theTemplate> AFX_COMDAT const AFX_MSGMAP theClass::messageMap = \
	  { &baseClass::messageMap, &theClass::_messageEntries[0] }; \
	  template <theTemplate> AFX_COMDAT const AFX_MSGMAP_ENTRY theClass::_messageEntries[] = \
	  { \

  #endif 
#else
  #define ITO_BEGIN_TEMPLATE_MESSAGE_MAP(theTemplate, theClass, baseClass) \
	  PTM_WARNING_DISABLE \
	  template <theTemplate> const AFX_MSGMAP* theClass::GetMessageMap() const \
		  { return GetThisMessageMap(); } \
	  template <theTemplate> const AFX_MSGMAP* PASCAL theClass::GetThisMessageMap() \
	  { \
		  typedef theClass ThisClass;						   \
		  typedef baseClass TheBaseClass;					   \
		  static const AFX_MSGMAP_ENTRY _messageEntries[] =  \
		  {
#endif // _MSC_VER < 1400

/**
    Begins the definition of message map for a fully-specialized template control

    Examples
    @code
      BEGIN_FULL_SPEC_TEMPLATE_MESSAGE_MAP(CMyEdit<T>, CEdit)

      BEGIN_FULL_SPEC_TEMPLATE_MESSAGE_MAP(CMyEdit<T, Y>, CEdit)
    
    @endcode

    Keep in mind that the second parameter cannot be a specialization.
    If you need this declare a typedef.
    For example, if you need:
      BEGIN_FULL_SPEC_TEMPLATE_MESSAGE_MAP(CMyEdit<T>, CEdit<T>)
    that statement won't compile.
    You need to do:
      typedef CEdit<T> CEditT;
    and
      BEGIN_FULL_SPEC_TEMPLATE_MESSAGE_MAP(CMyEdit<T>, CEditT)    
 */
#ifdef _AFXDLL
#define BEGIN_FULL_SPEC_TEMPLATE_MESSAGE_MAP(theClass, baseClass) \
	const AFX_MSGMAP* PASCAL theClass::GetThisMessageMap() \
		{ return &theClass::messageMap; } \
	const AFX_MSGMAP* theClass::GetMessageMap() const \
		{ return &theClass::messageMap; } \
	AFX_COMDAT const AFX_MSGMAP theClass::messageMap = \
	{ &baseClass::GetThisMessageMap, &theClass::_messageEntries[0] }; \
	AFX_COMDAT const AFX_MSGMAP_ENTRY theClass::_messageEntries[] = \
	{ \

#else
#define BEGIN_FULL_SPEC_TEMPLATE_MESSAGE_MAP(theClass, baseClass) \
	const AFX_MSGMAP* theClass::GetMessageMap() const \
		{ return &theClass::messageMap; } \
	AFX_COMDAT const AFX_MSGMAP theClass::messageMap = \
	{ &baseClass::messageMap, &theClass::_messageEntries[0] }; \
	AFX_COMDAT const AFX_MSGMAP_ENTRY theClass::_messageEntries[] = \
	{ \

#endif 


#endif // #ifndef _ITO33_WIN32_MFC_H_

