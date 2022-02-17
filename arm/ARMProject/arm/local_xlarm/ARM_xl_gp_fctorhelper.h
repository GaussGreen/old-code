/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * $Log: ARM_xl_gp_fctor_helper.h,v $
 * Revision 1.1  2004/01/13 15:08:43  ebenhamou
 * Initial version
 *
 */

/*! \file ARM_xl_gp_fctor_helper.h
 *
 *  \brief file for fctor construction
 *
 *	\author  E Benhamou
 *	\version 1.0
 *	\date September 2003
 */


#ifndef ARM_XL_GP_FCTOR_HELPER_LOCAL_H
#define ARM_XL_GP_FCTOR_HELPER_LOCAL_H

#include "ARM_xl_wrapper_local.h"

///////////////////////////////////////
//// functor with 0 argument
///////////////////////////////////////
class exportFunc0Arg : public ARMResultLong2LongFunc
{
public: 	
	typedef long (*exportFunc0)( ARM_result&, long );
	exportFunc0Arg( exportFunc0 func) : itsFunc(func)  {};
	long operator()( ARM_result& result, long objId ){	return (*itsFunc)( result, objId);}
private:
	exportFunc0 itsFunc;
};

///////////////////////////////////////
//// functor with 1 argument
///////////////////////////////////////
template <typename T1> class exportFunc1Arg : public ARMResultLong2LongFunc
{
public: 	
	typedef long (*exportFunc01)( const T1&, ARM_result&, long );
	exportFunc1Arg( const T1& arg1, exportFunc01 func): C_arg1(arg1), itsFunc(func)  {};
	long operator()( ARM_result& result, long objId ){	return (*itsFunc)( C_arg1, result, objId);}
private:
	const T1& C_arg1;
	exportFunc01 itsFunc;
};

///////////////////////////////////////
//// functor with 2 arguments
///////////////////////////////////////
template <typename T1, typename T2> class exportFunc2Args : public ARMResultLong2LongFunc
{
public: 	
	typedef long (*exportFunc02)( const T1&, const T2&,ARM_result&, long );
	exportFunc2Args ( const T1& arg1, const T2& arg2, exportFunc02 func) : C_arg1(arg1), C_arg2(arg2), itsFunc(func) {};
	long operator()( ARM_result& result, long objId ){ return (*itsFunc)( C_arg1, C_arg2, result, objId);}
private:
	const T1& C_arg1;
	const T2& C_arg2;
	exportFunc02 itsFunc;
};



///////////////////////////////////////
//// functor with 3 arguments
///////////////////////////////////////
template <typename T1, typename T2, typename T3 > class exportFunc3Args : public ARMResultLong2LongFunc
{
public: 	
	typedef long (*exportFunc03)( const T1&, const T2&,	const T3&, ARM_result&, long );
	exportFunc3Args ( const T1& arg1, const T2& arg2, const T3& arg3, exportFunc03 func) : C_arg1(arg1), C_arg2(arg2), C_arg3(arg3), itsFunc(func) {};
	long operator()( ARM_result& result, long objId ){ return (*itsFunc)( C_arg1, C_arg2, C_arg3, result, objId);}
private:
	const T1& C_arg1;
	const T2& C_arg2;
	const T3& C_arg3;
	exportFunc03 itsFunc;
};




///////////////////////////////////////
//// functor with 4 arguments
///////////////////////////////////////
template <typename T1, typename T2, typename T3,  typename T4 > class exportFunc4Args : public ARMResultLong2LongFunc
{
public: 	
	typedef long (*exportFunc04)( const T1&, const T2&,	const T3&, const T4&, ARM_result&, long );
	exportFunc4Args ( const T1& arg1, const T2& arg2, const T3& arg3, const T4& arg4, exportFunc04 func) : C_arg1(arg1), C_arg2(arg2), C_arg3(arg3), C_arg4(arg4), itsFunc(func) {};
	long operator()( ARM_result& result, long objId ){ return (*itsFunc)( C_arg1, C_arg2, C_arg3, C_arg4, result, objId);}
private:
	const T1& C_arg1;
	const T2& C_arg2;
	const T3& C_arg3;
	const T4& C_arg4;
	exportFunc04 itsFunc;
};


///////////////////////////////////////
//// functor with 5 arguments
///////////////////////////////////////
template <typename T1, typename T2, typename T3,  typename T4, typename T5 > class exportFunc5Args : public ARMResultLong2LongFunc
{
public: 	
	typedef long (*exportFunc05)( const T1&, const T2&,	const T3&, const T4&, const T5&, ARM_result&, long );
	exportFunc5Args ( const T1& arg1, const T2& arg2, const T3& arg3, const T4& arg4, const T5& arg5, exportFunc05 func) : C_arg1(arg1), C_arg2(arg2), C_arg3(arg3), C_arg4(arg4), C_arg5(arg5), itsFunc(func) {};
	long operator()( ARM_result& result, long objId ){ return (*itsFunc)( C_arg1, C_arg2, C_arg3, C_arg4, C_arg5, result, objId);}
private:
	const T1& C_arg1;
	const T2& C_arg2;
	const T3& C_arg3;
	const T4& C_arg4;
	const T5& C_arg5;
	exportFunc05 itsFunc;
};



///////////////////////////////////////
//// functor with 6 arguments
///////////////////////////////////////
template <typename T1, typename T2, typename T3,  typename T4, typename T5, typename T6 > class exportFunc6Args : public ARMResultLong2LongFunc
{
public: 	
	typedef long (*exportFunc06)( const T1&, const T2&,	const T3&, const T4&, const T5&, const T6&, ARM_result&, long );
	exportFunc6Args ( const T1& arg1, const T2& arg2, const T3& arg3, const T4& arg4, const T5& arg5, const T6& arg6, exportFunc06 func) : C_arg1(arg1), C_arg2(arg2), C_arg3(arg3), C_arg4(arg4), C_arg5(arg5), C_arg6(arg6), itsFunc(func) {};
	long operator()( ARM_result& result, long objId ){ return (*itsFunc)( C_arg1, C_arg2, C_arg3, C_arg4, C_arg5, C_arg6, result, objId);}
private:
	const T1& C_arg1;
	const T2& C_arg2;
	const T3& C_arg3;
	const T4& C_arg4;
	const T5& C_arg5;
	const T6& C_arg6;
	exportFunc06 itsFunc;
};



///////////////////////////////////////
//// functor with 7 arguments
///////////////////////////////////////
template <typename T1, typename T2, typename T3,  typename T4, typename T5, typename T6, typename T7 > class exportFunc7Args : public ARMResultLong2LongFunc
{
public: 	
	typedef long (*exportFunc07)( const T1&, const T2&,	const T3&, const T4&, const T5&, const T6&, const T7&, ARM_result&, long );
	exportFunc7Args ( const T1& arg1, const T2& arg2, const T3& arg3, const T4& arg4, const T5& arg5, const T6& arg6, const T7& arg7, exportFunc07 func) : C_arg1(arg1), C_arg2(arg2), C_arg3(arg3), C_arg4(arg4), C_arg5(arg5), C_arg6(arg6), C_arg7(arg7), itsFunc(func) {};
	long operator()( ARM_result& result, long objId ){ return (*itsFunc)( C_arg1, C_arg2, C_arg3, C_arg4, C_arg5, C_arg6, C_arg7, result, objId);}
private:
	const T1& C_arg1;
	const T2& C_arg2;
	const T3& C_arg3;
	const T4& C_arg4;
	const T5& C_arg5;
	const T6& C_arg6;
	const T7& C_arg7;
	exportFunc07 itsFunc;
};


///////////////////////////////////////
//// functor with 8 arguments
///////////////////////////////////////
template <typename T1, typename T2, typename T3,  typename T4, typename T5, typename T6, typename T7, typename T8 > class exportFunc8Args : public ARMResultLong2LongFunc
{
public: 	
	typedef long (*exportFunc08)( const T1&, const T2&,	const T3&, const T4&, const T5&, const T6&, const T7&, const T8&, ARM_result&, long );
	exportFunc8Args ( const T1& arg1, const T2& arg2, const T3& arg3, const T4& arg4, const T5& arg5, const T6& arg6, const T7& arg7, const T8& arg8, exportFunc08 func) : C_arg1(arg1), C_arg2(arg2), C_arg3(arg3), C_arg4(arg4), C_arg5(arg5), C_arg6(arg6), C_arg7(arg7), C_arg8(arg8), itsFunc(func) {};
	long operator()( ARM_result& result, long objId ){ return (*itsFunc)( C_arg1, C_arg2, C_arg3, C_arg4, C_arg5, C_arg6, C_arg7, C_arg8, result, objId);}
private:
	const T1& C_arg1;
	const T2& C_arg2;
	const T3& C_arg3;
	const T4& C_arg4;
	const T5& C_arg5;
	const T6& C_arg6;
	const T7& C_arg7;
	const T8& C_arg8;
	exportFunc08 itsFunc;
};


///////////////////////////////////////
//// functor with 9 arguments
///////////////////////////////////////
template <typename T1, typename T2, typename T3,  typename T4, typename T5, typename T6, typename T7, typename T8, typename T9 > class exportFunc9Args : public ARMResultLong2LongFunc
{
public: 	
	typedef long (*exportFunc09)( const T1&, const T2&,	const T3&, const T4&, const T5&, const T6&, const T7&, const T8& , const T9&, ARM_result&, long );
	exportFunc9Args ( const T1& arg1, const T2& arg2, const T3& arg3, const T4& arg4, const T5& arg5, const T6& arg6, const T7& arg7, const T8& arg8, const T9& arg9, exportFunc09 func) : C_arg1(arg1), C_arg2(arg2), C_arg3(arg3), C_arg4(arg4), C_arg5(arg5), C_arg6(arg6), C_arg7(arg7), C_arg8(arg8), C_arg9(arg9), itsFunc(func) {};
	long operator()( ARM_result& result, long objId )
	{ 
		return (*itsFunc)( C_arg1, C_arg2, C_arg3, C_arg4, C_arg5, C_arg6, C_arg7, C_arg8, C_arg9, result, objId);
	}
private:
	const T1& C_arg1;
	const T2& C_arg2;
	const T3& C_arg3;
	const T4& C_arg4;
	const T5& C_arg5;
	const T6& C_arg6;
	const T7& C_arg7;
	const T8& C_arg8;
	const T9& C_arg9;
	exportFunc09 itsFunc;
};

///////////////////////////////////////
//// functor with 10 arguments
///////////////////////////////////////
template <typename T1, typename T2, typename T3,  typename T4, typename T5, typename T6, typename T7, typename T8, typename T9, typename T10  > class exportFunc10Args : public ARMResultLong2LongFunc
{
public: 	
	typedef long (*exportFunc10)( const T1&, const T2&,	const T3&, const T4&, const T5&, const T6&, const T7&, const T8& , const T9&, const T10&, ARM_result&, long );
	exportFunc10Args ( const T1& arg1, const T2& arg2, const T3& arg3, const T4& arg4, const T5& arg5, const T6& arg6, const T7& arg7, const T8& arg8, const T9& arg9, const T10& arg10, exportFunc10 func)
	: C_arg1(arg1), C_arg2(arg2), C_arg3(arg3), C_arg4(arg4), C_arg5(arg5), C_arg6(arg6), C_arg7(arg7), C_arg8(arg8), C_arg9(arg9), C_arg10(arg10), itsFunc(func) {};
	long operator()( ARM_result& result, long objId ){ return (*itsFunc)( C_arg1, C_arg2, C_arg3, C_arg4, C_arg5, C_arg6, C_arg7, C_arg8, C_arg9, C_arg10, result, objId);}
private:
	const T1& C_arg1;
	const T2& C_arg2;
	const T3& C_arg3;
	const T4& C_arg4;
	const T5& C_arg5;
	const T6& C_arg6;
	const T7& C_arg7;
	const T8& C_arg8;
	const T9& C_arg9;
	const T10& C_arg10;
	exportFunc10 itsFunc;
};


///////////////////////////////////////
//// functor with 11 arguments
///////////////////////////////////////
template <typename T1, typename T2, typename T3,  typename T4, typename T5, typename T6, typename T7, typename T8, typename T9, typename T10,
		  typename T11  > class exportFunc11Args : public ARMResultLong2LongFunc
{
public: 	
	typedef long (*exportFunc11)( const T1&, const T2&,	const T3&, const T4&, const T5&, const T6&, const T7&, const T8& , const T9&, const T10&,
		 const  T11&, ARM_result&, long );
	exportFunc11Args ( const T1& arg1, const T2& arg2, const T3& arg3, const T4& arg4, const T5& arg5, const T6& arg6, const T7& arg7, const T8& arg8, const T9& arg9, const T10& arg10, const T11& arg11, exportFunc11 func)
	: C_arg1(arg1), C_arg2(arg2), C_arg3(arg3), C_arg4(arg4), C_arg5(arg5), C_arg6(arg6), C_arg7(arg7), C_arg8(arg8), C_arg9(arg9), C_arg10(arg10),
		 C_arg11(arg11), itsFunc(func) {};
	long operator()( ARM_result& result, long objId ){ return (*itsFunc)( C_arg1, C_arg2, C_arg3, C_arg4, C_arg5, C_arg6, C_arg7, C_arg8, C_arg9, C_arg10,
		 C_arg11, result, objId);}
private:
	const T1& C_arg1;
	const T2& C_arg2;
	const T3& C_arg3;
	const T4& C_arg4;
	const T5& C_arg5;
	const T6& C_arg6;
	const T7& C_arg7;
	const T8& C_arg8;
	const T9& C_arg9;
	const T10& C_arg10;
	const T11& C_arg11;
	exportFunc11 itsFunc;
};

///////////////////////////////////////
//// functor with 12 arguments
///////////////////////////////////////
template <typename T1, typename T2, typename T3,  typename T4, typename T5, typename T6, typename T7, typename T8, typename T9, typename T10,
		  typename T11, typename T12  > class exportFunc12Args : public ARMResultLong2LongFunc
{
public: 	
	typedef long (*exportFunc12)( const T1&, const T2&,	const T3&, const T4&, const T5&, const T6&, const T7&, const T8& , const T9&, const T10&,
		 const  T11&, const T12&, ARM_result&, long );
	exportFunc12Args ( const T1& arg1, const T2& arg2, const T3& arg3, const T4& arg4, const T5& arg5, const T6& arg6, const T7& arg7, const T8& arg8, const T9& arg9, const T10& arg10, const T11& arg11, const T12& arg12, exportFunc12 func)
	: C_arg1(arg1), C_arg2(arg2), C_arg3(arg3), C_arg4(arg4), C_arg5(arg5), C_arg6(arg6), C_arg7(arg7), C_arg8(arg8), C_arg9(arg9), C_arg10(arg10),
		 C_arg11(arg11), C_arg12(arg12), itsFunc(func) {};
	long operator()( ARM_result& result, long objId ){ return (*itsFunc)( C_arg1, C_arg2, C_arg3, C_arg4, C_arg5, C_arg6, C_arg7, C_arg8, C_arg9, C_arg10,
		 C_arg11, C_arg12, result, objId);}
private:
	const T1& C_arg1;
	const T2& C_arg2;
	const T3& C_arg3;
	const T4& C_arg4;
	const T5& C_arg5;
	const T6& C_arg6;
	const T7& C_arg7;
	const T8& C_arg8;
	const T9& C_arg9;
	const T10& C_arg10;
	const T11& C_arg11;
	const T12& C_arg12;
	exportFunc12 itsFunc;
};


///////////////////////////////////////
//// functor with 13 arguments
///////////////////////////////////////
template <typename T1, typename T2, typename T3,  typename T4, typename T5, typename T6, typename T7, typename T8, typename T9, typename T10,
		  typename T11, typename T12, typename T13  > class exportFunc13Args : public ARMResultLong2LongFunc
{
public: 	
	typedef long (*exportFunc13)( const T1&, const T2&,	const T3&, const T4&, const T5&, const T6&, const T7&, const T8& , const T9&, const T10&,
		 const  T11&, const T12&, const T13&, ARM_result&, long );
	exportFunc13Args ( const T1& arg1, const T2& arg2, const T3& arg3, const T4& arg4, const T5& arg5, const T6& arg6, const T7& arg7, const T8& arg8, const T9& arg9, const T10& arg10, const T11& arg11, const T12& arg12, T13& arg13, exportFunc13 func)
	: C_arg1(arg1), C_arg2(arg2), C_arg3(arg3), C_arg4(arg4), C_arg5(arg5), C_arg6(arg6), C_arg7(arg7), C_arg8(arg8), C_arg9(arg9), C_arg10(arg10),
		 C_arg11(arg11), C_arg12(arg12), C_arg13(arg13), itsFunc(func) {};
	long operator()( ARM_result& result, long objId ){ return (*itsFunc)( C_arg1, C_arg2, C_arg3, C_arg4, C_arg5, C_arg6, C_arg7, C_arg8, C_arg9, C_arg10,
		 C_arg11, C_arg12, C_arg13, result, objId);}
private:
	const T1& C_arg1;
	const T2& C_arg2;
	const T3& C_arg3;
	const T4& C_arg4;
	const T5& C_arg5;
	const T6& C_arg6;
	const T7& C_arg7;
	const T8& C_arg8;
	const T9& C_arg9;
	const T10& C_arg10;
	const T11& C_arg11;
	const T12& C_arg12;
	const T13& C_arg13;
	exportFunc13 itsFunc;
};


///////////////////////////////////////
//// functor with 14 arguments
///////////////////////////////////////
template <typename T1, typename T2, typename T3,  typename T4, typename T5, typename T6, typename T7, typename T8, typename T9, typename T10,
		  typename T11, typename T12, typename T13, typename T14  > class exportFunc14Args : public ARMResultLong2LongFunc
{
public: 	
	typedef long (*exportFunc14)( const T1&, const T2&,	const T3&, const T4&, const T5&, const T6&, const T7&, const T8& , const T9&, const T10&,
		 const  T11&, const T12&, const T13&, const T14&,  ARM_result&, long );
	exportFunc14Args ( const T1& arg1, const T2& arg2, const T3& arg3, const T4& arg4, const T5& arg5, const T6& arg6, const T7& arg7, const T8& arg8, const T9& arg9, const T10& arg10, const T11& arg11, const T12& arg12, T13& arg13, const T14& arg14, exportFunc14 func)
	: C_arg1(arg1), C_arg2(arg2), C_arg3(arg3), C_arg4(arg4), C_arg5(arg5), C_arg6(arg6), C_arg7(arg7), C_arg8(arg8), C_arg9(arg9), C_arg10(arg10),
		 C_arg11(arg11), C_arg12(arg12), C_arg13(arg13), C_arg14(arg14), itsFunc(func) {};
	long operator()( ARM_result& result, long objId ){ return (*itsFunc)( C_arg1, C_arg2, C_arg3, C_arg4, C_arg5, C_arg6, C_arg7, C_arg8, C_arg9, C_arg10,
		 C_arg11, C_arg12, C_arg13, C_arg14, result, objId);}
private:
	const T1& C_arg1;
	const T2& C_arg2;
	const T3& C_arg3;
	const T4& C_arg4;
	const T5& C_arg5;
	const T6& C_arg6;
	const T7& C_arg7;
	const T8& C_arg8;
	const T9& C_arg9;
	const T10& C_arg10;
	const T11& C_arg11;
	const T12& C_arg12;
	const T13& C_arg13;
	const T14& C_arg14;
	exportFunc14 itsFunc;
};

///////////////////////////////////////
//// functor with 15 arguments
///////////////////////////////////////
template <typename T1, typename T2, typename T3,  typename T4, typename T5, typename T6, typename T7, typename T8, typename T9, typename T10,
		  typename T11, typename T12, typename T13, typename T14, typename T15  > class exportFunc15Args : public ARMResultLong2LongFunc
{
public: 	
	typedef long (*exportFunc15)( const T1&, const T2&,	const T3&, const T4&, const T5&, const T6&, const T7&, const T8& , const T9&, const T10&,
		 const  T11&, const T12&, const T13&, const T14&, const T15&, ARM_result&, long );
	exportFunc15Args ( const T1& arg1, const T2& arg2, const T3& arg3, const T4& arg4, const T5& arg5, const T6& arg6, const T7& arg7, const T8& arg8, const T9& arg9, const T10& arg10, const T11& arg11, const T12& arg12, const T13& arg13, const T14& arg14, const T15& arg15, exportFunc15 func)
	: C_arg1(arg1), C_arg2(arg2), C_arg3(arg3), C_arg4(arg4), C_arg5(arg5), C_arg6(arg6), C_arg7(arg7), C_arg8(arg8), C_arg9(arg9), C_arg10(arg10),
		 C_arg11(arg11), C_arg12(arg12), C_arg13(arg13), C_arg14(arg14), C_arg15(arg15), itsFunc(func) {};
	long operator()( ARM_result& result, long objId ){ return (*itsFunc)( C_arg1, C_arg2, C_arg3, C_arg4, C_arg5, C_arg6, C_arg7, C_arg8, C_arg9, C_arg10,
		 C_arg11, C_arg12, C_arg13, C_arg14, C_arg15, result, objId);}
private:
	const T1& C_arg1;
	const T2& C_arg2;
	const T3& C_arg3;
	const T4& C_arg4;
	const T5& C_arg5;
	const T6& C_arg6;
	const T7& C_arg7;
	const T8& C_arg8;
	const T9& C_arg9;
	const T10& C_arg10;
	const T11& C_arg11;
	const T12& C_arg12;
	const T13& C_arg13;
	const T14& C_arg14;
	const T15& C_arg15;
	exportFunc15 itsFunc;
};

///////////////////////////////////////
//// functor with 16 arguments
///////////////////////////////////////
template <typename T1, typename T2, typename T3,  typename T4, typename T5, typename T6, typename T7, typename T8, typename T9, typename T10,
		  typename T11, typename T12, typename T13, typename T14, typename T15, typename T16  > class exportFunc16Args : public ARMResultLong2LongFunc
{
public: 	
	typedef long (*exportFunc16)( const T1&, const T2&,	const T3&, const T4&, const T5&, const T6&, const T7&, const T8& , const T9&, const T10&,
		 const  T11&, const T12&, const T13&, const T14&, const T15&, const T16&, ARM_result&, long );
	exportFunc16Args ( const T1& arg1, const T2& arg2, const T3& arg3, const T4& arg4, const T5& arg5, const T6& arg6, const T7& arg7, const T8& arg8, const T9& arg9, const T10& arg10, const T11& arg11, const T12& arg12, const T13& arg13, const T14& arg14, const T15& arg15, const T16& arg16, exportFunc16 func)
	: C_arg1(arg1), C_arg2(arg2), C_arg3(arg3), C_arg4(arg4), C_arg5(arg5), C_arg6(arg6), C_arg7(arg7), C_arg8(arg8), C_arg9(arg9), C_arg10(arg10),
		 C_arg11(arg11), C_arg12(arg12), C_arg13(arg13), C_arg14(arg14), C_arg15(arg15), C_arg16(arg16), itsFunc(func) {};
	long operator()( ARM_result& result, long objId ){ return (*itsFunc)( C_arg1, C_arg2, C_arg3, C_arg4, C_arg5, C_arg6, C_arg7, C_arg8, C_arg9, C_arg10,
		 C_arg11, C_arg12, C_arg13, C_arg14, C_arg15, C_arg16, result, objId);}
private:
	const T1& C_arg1;
	const T2& C_arg2;
	const T3& C_arg3;
	const T4& C_arg4;
	const T5& C_arg5;
	const T6& C_arg6;
	const T7& C_arg7;
	const T8& C_arg8;
	const T9& C_arg9;
	const T10& C_arg10;
	const T11& C_arg11;
	const T12& C_arg12;
	const T13& C_arg13;
	const T14& C_arg14;
	const T15& C_arg15;
	const T16& C_arg16;
	exportFunc16 itsFunc;
};

///////////////////////////////////////
//// functor with 17 arguments
///////////////////////////////////////
template <typename T1, typename T2, typename T3,  typename T4, typename T5, typename T6, typename T7, typename T8, typename T9, typename T10,
		  typename T11, typename T12, typename T13, typename T14, typename T15, typename T16, typename T17  > class exportFunc17Args : public ARMResultLong2LongFunc
{
public: 	
	typedef long (*exportFunc17)( const T1&, const T2&,	const T3&, const T4&, const T5&, const T6&, const T7&, const T8& , const T9&, const T10&,
		 const  T11&, const T12&, const T13&, const T14&, const T15&, const T16&, const T17&, ARM_result&, long );
	exportFunc17Args ( const T1& arg1, const T2& arg2, const T3& arg3, const T4& arg4, const T5& arg5, const T6& arg6, const T7& arg7, const T8& arg8, const T9& arg9, const T10& arg10, const T11& arg11, const T12& arg12, const T13& arg13, const T14& arg14, const T15& arg15, const T16& arg16, const T17& arg17, exportFunc17 func)
	: C_arg1(arg1), C_arg2(arg2), C_arg3(arg3), C_arg4(arg4), C_arg5(arg5), C_arg6(arg6), C_arg7(arg7), C_arg8(arg8), C_arg9(arg9), C_arg10(arg10),
		 C_arg11(arg11), C_arg12(arg12), C_arg13(arg13), C_arg14(arg14), C_arg15(arg15), C_arg16(arg16), C_arg17(arg17), itsFunc(func) {};
	long operator()( ARM_result& result, long objId ){ return (*itsFunc)( C_arg1, C_arg2, C_arg3, C_arg4, C_arg5, C_arg6, C_arg7, C_arg8, C_arg9, C_arg10,
		 C_arg11, C_arg12, C_arg13, C_arg14, C_arg15, C_arg16, C_arg17, result, objId);}
private:
	const T1& C_arg1;
	const T2& C_arg2;
	const T3& C_arg3;
	const T4& C_arg4;
	const T5& C_arg5;
	const T6& C_arg6;
	const T7& C_arg7;
	const T8& C_arg8;
	const T9& C_arg9;
	const T10& C_arg10;
	const T11& C_arg11;
	const T12& C_arg12;
	const T13& C_arg13;
	const T14& C_arg14;
	const T15& C_arg15;
	const T16& C_arg16;
	const T17& C_arg17;
	exportFunc17 itsFunc;
};

///////////////////////////////////////
//// functor with 18 arguments
///////////////////////////////////////
template <typename T1, typename T2, typename T3,  typename T4, typename T5, typename T6, typename T7, typename T8, typename T9, typename T10,
		  typename T11, typename T12, typename T13, typename T14, typename T15, typename T16, typename T17, typename T18  > class exportFunc18Args : public ARMResultLong2LongFunc
{
public: 	
	typedef long (*exportFunc18)( const T1&, const T2&,	const T3&, const T4&, const T5&, const T6&, const T7&, const T8& , const T9&, const T10&,
		 const  T11&, const T12&, const T13&, const T14&, const T15&, const T16&, const T17&, const T18&, ARM_result&, long );
	exportFunc18Args ( const T1& arg1, const T2& arg2, const T3& arg3, const T4& arg4, const T5& arg5, const T6& arg6, const T7& arg7, const T8& arg8, const T9& arg9, const T10& arg10, const T11& arg11, const T12& arg12, const T13& arg13, const T14& arg14, const T15& arg15, const T16& arg16, const T17& arg17, const T18& arg18, exportFunc18 func)
	: C_arg1(arg1), C_arg2(arg2), C_arg3(arg3), C_arg4(arg4), C_arg5(arg5), C_arg6(arg6), C_arg7(arg7), C_arg8(arg8), C_arg9(arg9), C_arg10(arg10),
		 C_arg11(arg11), C_arg12(arg12), C_arg13(arg13), C_arg14(arg14), C_arg15(arg15), C_arg16(arg16), C_arg17(arg17), C_arg18(arg18), itsFunc(func) {};
	long operator()( ARM_result& result, long objId ){ return (*itsFunc)( C_arg1, C_arg2, C_arg3, C_arg4, C_arg5, C_arg6, C_arg7, C_arg8, C_arg9, C_arg10,
		 C_arg11, C_arg12, C_arg13, C_arg14, C_arg15, C_arg16, C_arg17, C_arg18, result, objId);}
private:
	const T1& C_arg1;
	const T2& C_arg2;
	const T3& C_arg3;
	const T4& C_arg4;
	const T5& C_arg5;
	const T6& C_arg6;
	const T7& C_arg7;
	const T8& C_arg8;
	const T9& C_arg9;
	const T10& C_arg10;
	const T11& C_arg11;
	const T12& C_arg12;
	const T13& C_arg13;
	const T14& C_arg14;
	const T15& C_arg15;
	const T16& C_arg16;
	const T17& C_arg17;
	const T18& C_arg18;
	exportFunc18 itsFunc;
};

///////////////////////////////////////
//// functor with 20 arguments
///////////////////////////////////////
template <typename T1, typename T2, typename T3,  typename T4, typename T5, typename T6, typename T7, typename T8, typename T9, typename T10,
		  typename T11, typename T12, typename T13, typename T14, typename T15, typename T16, typename T17, typename T18, typename T19, typename T20  > class exportFunc20Args : public ARMResultLong2LongFunc
{
public: 	
	typedef long (*exportFunc20)( const T1&, const T2&,	const T3&, const T4&, const T5&, const T6&, const T7&, const T8& , const T9&, const T10&,
		 const  T11&, const T12&, const T13&, const T14&, const T15&, const T16&, const T17&, const T18&, const T19&, const T20&, ARM_result&, long );
	exportFunc20Args ( const T1& arg1, const T2& arg2, const T3& arg3, const T4& arg4, const T5& arg5, const T6& arg6, const T7& arg7, const T8& arg8, const T9& arg9, const T10& arg10, const T11& arg11, const T12& arg12, const T13& arg13, const T14& arg14, const T15& arg15, const T16& arg16, 
		const T17& arg17, const T18& arg18, const T19& arg19, const T20& arg20, exportFunc20 func)
	: C_arg1(arg1), C_arg2(arg2), C_arg3(arg3), C_arg4(arg4), C_arg5(arg5), C_arg6(arg6), C_arg7(arg7), C_arg8(arg8), C_arg9(arg9), C_arg10(arg10),
		 C_arg11(arg11), C_arg12(arg12), C_arg13(arg13), C_arg14(arg14), C_arg15(arg15), C_arg16(arg16), 
		 C_arg17(arg17), C_arg18(arg18), C_arg19(arg19), C_arg20(arg20), itsFunc(func) {};
	long operator()( ARM_result& result, long objId ){ return (*itsFunc)( C_arg1, C_arg2, C_arg3, C_arg4, C_arg5, C_arg6, C_arg7, C_arg8, C_arg9, C_arg10,
		 C_arg11, C_arg12, C_arg13, C_arg14, C_arg15, C_arg16, C_arg17, C_arg18, C_arg19, C_arg20, result, objId);}
private:
	const T1& C_arg1;
	const T2& C_arg2;
	const T3& C_arg3;
	const T4& C_arg4;
	const T5& C_arg5;
	const T6& C_arg6;
	const T7& C_arg7;
	const T8& C_arg8;
	const T9& C_arg9;
	const T10& C_arg10;
	const T11& C_arg11;
	const T12& C_arg12;
	const T13& C_arg13;
	const T14& C_arg14;
	const T15& C_arg15;
	const T16& C_arg16;
	const T17& C_arg17;
	const T18& C_arg18;
	const T19& C_arg19;
	const T20& C_arg20;

	exportFunc20 itsFunc;
};

///////////////////////////////////////
//// functor with 21 arguments  ///////
///////////////////////////////////////

template < typename T1, typename T2, typename T3,  typename T4, typename T5, typename T6, typename T7, typename T8,
typename T9, typename T10, typename T11, typename T12, typename T13, typename T14, typename T15, typename T16, 
typename T17, typename T18, typename T19, typename T20, typename T21> class exportFunc21Args : public ARMResultLong2LongFunc
{
public: 	

	typedef long (*exportFunc21)( const T1&, const T2&,	const T3&, const T4&, const T5&, const T6&, const T7&, const T8& , const T9&, 
		const T10&, const T11&, const T12&, const T13&, const T14&, const T15&,const T16&, const T17&, const T18&, const T19&, const T20&, 
		const T21&, ARM_result&, long );

	exportFunc21Args ( const T1& arg1, const T2& arg2, const T3& arg3, const T4& arg4, const T5& arg5, const T6& arg6, const T7& arg7, const T8& arg8, const T9& arg9, const T10& arg10, 
		const T11& arg11, const T12& arg12, const T13& arg13, const T14& arg14, const T15& arg15, const T16& arg16, const T17& arg17, const T18& arg18, const T19& arg19, const T20& arg20, 
		const T21& arg21, exportFunc21 func)
	: C_arg1(arg1), C_arg2(arg2), C_arg3(arg3), C_arg4(arg4), C_arg5(arg5), C_arg6(arg6), C_arg7(arg7), C_arg8(arg8), C_arg9(arg9), C_arg10(arg10), 	
	C_arg11(arg11), C_arg12(arg12), C_arg13(arg13), C_arg14(arg14), C_arg15(arg15), C_arg16(arg16), C_arg17(arg17), C_arg18(arg18), C_arg19(arg19), C_arg20(arg20), 
	C_arg21(arg21), itsFunc(func) {};

	long operator()( ARM_result& result, long objId )
	{ 
		return (*itsFunc)( C_arg1, C_arg2, C_arg3, C_arg4, C_arg5, C_arg6, C_arg7, C_arg8, C_arg9, C_arg10, C_arg11, C_arg12, C_arg13, C_arg14, C_arg15, C_arg16, 
		C_arg17, C_arg18, C_arg19, C_arg20, C_arg21, result, objId);
	}

private:

	const T1&	C_arg1;
	const T2&	C_arg2;
	const T3&	C_arg3;
	const T4&	C_arg4;
	const T5&	C_arg5;
	const T6&	C_arg6;
	const T7&	C_arg7;
	const T8&	C_arg8;
	const T9&	C_arg9;
	const T10&	C_arg10;
	const T11&	C_arg11;
	const T12&	C_arg12;
	const T13&	C_arg13;
	const T14&	C_arg14;
	const T15&	C_arg15;
	const T16&	C_arg16;
	const T17&	C_arg17;
	const T18&	C_arg18;
	const T19&	C_arg19;
	const T20&	C_arg20;
	const T21&	C_arg21;

	exportFunc21 itsFunc;
};

///////////////////////////////////////
//// functor with 22 arguments  ///////
///////////////////////////////////////
template < typename T1, typename T2, typename T3,  typename T4, typename T5, typename T6, typename T7, typename T8,
typename T9, typename T10, typename T11, typename T12, typename T13, typename T14, typename T15, typename T16, 
typename T17, typename T18, typename T19, typename T20, typename T21, typename T22> class exportFunc22Args : public ARMResultLong2LongFunc
{
public: 	

	typedef long (*exportFunc22)( const T1&, const T2&,	const T3&, const T4&, const T5&, const T6&, const T7&, const T8& , const T9&, 
		const T10&, const T11&, const T12&, const T13&, const T14&, const T15&,const T16&, const T17&, const T18&, const T19&, const T20&, 
		const T21&, const T22&, ARM_result&, long );

	exportFunc22Args ( const T1& arg1, const T2& arg2, const T3& arg3, const T4& arg4, const T5& arg5, const T6& arg6, const T7& arg7, const T8& arg8, const T9& arg9, const T10& arg10, 
		const T11& arg11, const T12& arg12, const T13& arg13, const T14& arg14, const T15& arg15, const T16& arg16, const T17& arg17, const T18& arg18, const T19& arg19, const T20& arg20, 
		const T21& arg21, const T22& arg22, exportFunc22 func)
	: C_arg1(arg1), C_arg2(arg2), C_arg3(arg3), C_arg4(arg4), C_arg5(arg5), C_arg6(arg6), C_arg7(arg7), C_arg8(arg8), C_arg9(arg9), C_arg10(arg10), 	
	C_arg11(arg11), C_arg12(arg12), C_arg13(arg13), C_arg14(arg14), C_arg15(arg15), C_arg16(arg16), C_arg17(arg17), C_arg18(arg18), C_arg19(arg19), C_arg20(arg20), 
	C_arg21(arg21), C_arg22(arg22), itsFunc(func) {};

	long operator()( ARM_result& result, long objId )
	{ 
		return (*itsFunc)( C_arg1, C_arg2, C_arg3, C_arg4, C_arg5, C_arg6, C_arg7, C_arg8, C_arg9, C_arg10, C_arg11, C_arg12, C_arg13, C_arg14, C_arg15, C_arg16, 
		C_arg17, C_arg18, C_arg19, C_arg20, C_arg21, C_arg22, result, objId);
	}

private:

	const T1&	C_arg1;
	const T2&	C_arg2;
	const T3&	C_arg3;
	const T4&	C_arg4;
	const T5&	C_arg5;
	const T6&	C_arg6;
	const T7&	C_arg7;
	const T8&	C_arg8;
	const T9&	C_arg9;
	const T10&	C_arg10;
	const T11&	C_arg11;
	const T12&	C_arg12;
	const T13&	C_arg13;
	const T14&	C_arg14;
	const T15&	C_arg15;
	const T16&	C_arg16;
	const T17&	C_arg17;
	const T18&	C_arg18;
	const T19&	C_arg19;
	const T20&	C_arg20;
	const T21&	C_arg21;
	const T22&	C_arg22;

	exportFunc22 itsFunc;
};

///////////////////////////////////////
//// functor with 23 arguments  ///////
///////////////////////////////////////
template < typename T1, typename T2, typename T3,  typename T4, typename T5, typename T6, typename T7, typename T8,
typename T9, typename T10, typename T11, typename T12, typename T13, typename T14, typename T15, typename T16, 
typename T17, typename T18, typename T19, typename T20, typename T21, typename T22, typename T23 > class exportFunc23Args : public ARMResultLong2LongFunc
{
public: 	

	typedef long (*exportFunc23)( const T1&, const T2&,	const T3&, const T4&, const T5&, const T6&, const T7&, const T8& , const T9&, 
		const T10&, const T11&, const T12&, const T13&, const T14&, const T15&,const T16&, const T17&, const T18&, const T19&, const T20&, 
		const T21&, const T22&, const T23&, ARM_result&, long );

	exportFunc23Args ( const T1& arg1, const T2& arg2, const T3& arg3, const T4& arg4, const T5& arg5, const T6& arg6, const T7& arg7, const T8& arg8, const T9& arg9, const T10& arg10, 
		const T11& arg11, const T12& arg12, const T13& arg13, const T14& arg14, const T15& arg15, const T16& arg16, const T17& arg17, const T18& arg18, const T19& arg19, const T20& arg20, 
		const T21& arg21, const T22& arg22, const T23& arg23, exportFunc23 func)
	: C_arg1(arg1), C_arg2(arg2), C_arg3(arg3), C_arg4(arg4), C_arg5(arg5), C_arg6(arg6), C_arg7(arg7), C_arg8(arg8), C_arg9(arg9), C_arg10(arg10), 	
	C_arg11(arg11), C_arg12(arg12), C_arg13(arg13), C_arg14(arg14), C_arg15(arg15), C_arg16(arg16), C_arg17(arg17), C_arg18(arg18), C_arg19(arg19), C_arg20(arg20), 
	C_arg21(arg21), C_arg22(arg22), C_arg23(arg23), itsFunc(func) {};

	long operator()( ARM_result& result, long objId )
	{ 
		return (*itsFunc)( C_arg1, C_arg2, C_arg3, C_arg4, C_arg5, C_arg6, C_arg7, C_arg8, C_arg9, C_arg10, C_arg11, C_arg12, C_arg13, C_arg14, C_arg15, C_arg16, 
		C_arg17, C_arg18, C_arg19, C_arg20, C_arg21, C_arg22, C_arg23, result, objId);
	}

private:

	const T1&	C_arg1;
	const T2&	C_arg2;
	const T3&	C_arg3;
	const T4&	C_arg4;
	const T5&	C_arg5;
	const T6&	C_arg6;
	const T7&	C_arg7;
	const T8&	C_arg8;
	const T9&	C_arg9;
	const T10&	C_arg10;
	const T11&	C_arg11;
	const T12&	C_arg12;
	const T13&	C_arg13;
	const T14&	C_arg14;
	const T15&	C_arg15;
	const T16&	C_arg16;
	const T17&	C_arg17;
	const T18&	C_arg18;
	const T19&	C_arg19;
	const T20&	C_arg20;
	const T21&	C_arg21;
	const T22&	C_arg22;
	const T23&	C_arg23;

	exportFunc23 itsFunc;
};

///////////////////////////////////////
//// functor with 24 arguments  ///////
///////////////////////////////////////
template < typename T1, typename T2, typename T3,  typename T4, typename T5, typename T6, typename T7, typename T8,
typename T9, typename T10, typename T11, typename T12, typename T13, typename T14, typename T15, typename T16, 
typename T17, typename T18, typename T19, typename T20, typename T21, typename T22, typename T23, typename T24 > class exportFunc24Args : public ARMResultLong2LongFunc
{
public: 	

	typedef long (*exportFunc24)( const T1&, const T2&,	const T3&, const T4&, const T5&, const T6&, const T7&, const T8& , const T9&, 
		const T10&, const T11&, const T12&, const T13&, const T14&, const T15&,const T16&, const T17&, const T18&, const T19&, const T20&, 
		const T21&, const T22&, const T23&, const T24&, ARM_result&, long );

	exportFunc24Args ( const T1& arg1, const T2& arg2, const T3& arg3, const T4& arg4, const T5& arg5, const T6& arg6, const T7& arg7, const T8& arg8, const T9& arg9, const T10& arg10, 
		const T11& arg11, const T12& arg12, const T13& arg13, const T14& arg14, const T15& arg15, const T16& arg16, const T17& arg17, const T18& arg18, const T19& arg19, const T20& arg20, 
		const T21& arg21, const T22& arg22, const T23& arg23, const T24& arg24, exportFunc24 func)
	: C_arg1(arg1), C_arg2(arg2), C_arg3(arg3), C_arg4(arg4), C_arg5(arg5), C_arg6(arg6), C_arg7(arg7), C_arg8(arg8), C_arg9(arg9), C_arg10(arg10), 	
	C_arg11(arg11), C_arg12(arg12), C_arg13(arg13), C_arg14(arg14), C_arg15(arg15), C_arg16(arg16), C_arg17(arg17), C_arg18(arg18), C_arg19(arg19), C_arg20(arg20), 
	C_arg21(arg21), C_arg22(arg22), C_arg23(arg23), C_arg24(arg24), itsFunc(func) {};

	long operator()( ARM_result& result, long objId )
	{ 
		return (*itsFunc)( C_arg1, C_arg2, C_arg3, C_arg4, C_arg5, C_arg6, C_arg7, C_arg8, C_arg9, C_arg10, C_arg11, C_arg12, C_arg13, C_arg14, C_arg15, C_arg16, 
		C_arg17, C_arg18, C_arg19, C_arg20, C_arg21, C_arg22, C_arg23, C_arg24, result, objId);
	}

private:

	const T1&	C_arg1;
	const T2&	C_arg2;
	const T3&	C_arg3;
	const T4&	C_arg4;
	const T5&	C_arg5;
	const T6&	C_arg6;
	const T7&	C_arg7;
	const T8&	C_arg8;
	const T9&	C_arg9;
	const T10&	C_arg10;
	const T11&	C_arg11;
	const T12&	C_arg12;
	const T13&	C_arg13;
	const T14&	C_arg14;
	const T15&	C_arg15;
	const T16&	C_arg16;
	const T17&	C_arg17;
	const T18&	C_arg18;
	const T19&	C_arg19;
	const T20&	C_arg20;
	const T21&	C_arg21;
	const T22&	C_arg22;
	const T23&	C_arg23;
	const T24&	C_arg24;

	exportFunc24 itsFunc;
};

///////////////////////////////////////
//// functor with 25 arguments  ///////
///////////////////////////////////////
template < typename T1, typename T2, typename T3,  typename T4, typename T5, typename T6, typename T7, typename T8,
typename T9, typename T10, typename T11, typename T12, typename T13, typename T14, typename T15, typename T16, 
typename T17, typename T18, typename T19, typename T20, typename T21, typename T22, typename T23, typename T24, typename T25 > class exportFunc25Args : public ARMResultLong2LongFunc
{
public: 	

	typedef long (*exportFunc25)( const T1&, const T2&,	const T3&, const T4&, const T5&, const T6&, const T7&, const T8& , const T9&, 
		const T10&, const T11&, const T12&, const T13&, const T14&, const T15&,const T16&, const T17&, const T18&, const T19&, const T20&, 
		const T21&, const T22&, const T23&, const T24&, const T25&, ARM_result&, long );

	exportFunc25Args ( const T1& arg1, const T2& arg2, const T3& arg3, const T4& arg4, const T5& arg5, const T6& arg6, const T7& arg7, const T8& arg8, const T9& arg9, const T10& arg10, 
		const T11& arg11, const T12& arg12, const T13& arg13, const T14& arg14, const T15& arg15, const T16& arg16, const T17& arg17, const T18& arg18, const T19& arg19, const T20& arg20, 
		const T21& arg21, const T22& arg22, const T23& arg23, const T24& arg24, const T25& arg25, exportFunc25 func)
	: C_arg1(arg1), C_arg2(arg2), C_arg3(arg3), C_arg4(arg4), C_arg5(arg5), C_arg6(arg6), C_arg7(arg7), C_arg8(arg8), C_arg9(arg9), C_arg10(arg10), 	
	C_arg11(arg11), C_arg12(arg12), C_arg13(arg13), C_arg14(arg14), C_arg15(arg15), C_arg16(arg16), C_arg17(arg17), C_arg18(arg18), C_arg19(arg19), C_arg20(arg20), 
	C_arg21(arg21), C_arg22(arg22), C_arg23(arg23), C_arg24(arg24), C_arg25(arg25), itsFunc(func) {};

	long operator()( ARM_result& result, long objId )
	{ 
		return (*itsFunc)( C_arg1, C_arg2, C_arg3, C_arg4, C_arg5, C_arg6, C_arg7, C_arg8, C_arg9, C_arg10, C_arg11, C_arg12, C_arg13, C_arg14, C_arg15, C_arg16, 
		C_arg17, C_arg18, C_arg19, C_arg20, C_arg21, C_arg22, C_arg23, C_arg24, C_arg25, result, objId);
	}

private:

	const T1&	C_arg1;
	const T2&	C_arg2;
	const T3&	C_arg3;
	const T4&	C_arg4;
	const T5&	C_arg5;
	const T6&	C_arg6;
	const T7&	C_arg7;
	const T8&	C_arg8;
	const T9&	C_arg9;
	const T10&	C_arg10;
	const T11&	C_arg11;
	const T12&	C_arg12;
	const T13&	C_arg13;
	const T14&	C_arg14;
	const T15&	C_arg15;
	const T16&	C_arg16;
	const T17&	C_arg17;
	const T18&	C_arg18;
	const T19&	C_arg19;
	const T20&	C_arg20;
	const T21&	C_arg21;
	const T22&	C_arg22;
	const T23&	C_arg23;
	const T24&	C_arg24;
	const T25&	C_arg25;

	exportFunc25 itsFunc;
};

///////////////////////////////////////
//// functor with 26 arguments  ///////
///////////////////////////////////////
template < typename T1, typename T2, typename T3,  typename T4, typename T5, typename T6, typename T7, typename T8,
typename T9, typename T10, typename T11, typename T12, typename T13, typename T14, typename T15, typename T16, 
typename T17, typename T18, typename T19, typename T20, typename T21, typename T22, typename T23, typename T24, typename T25, typename T26 > class exportFunc26Args : public ARMResultLong2LongFunc
{
public: 	

	typedef long (*exportFunc26)( const T1&, const T2&,	const T3&, const T4&, const T5&, const T6&, const T7&, const T8& , const T9&, 
		const T10&, const T11&, const T12&, const T13&, const T14&, const T15&,const T16&, const T17&, const T18&, const T19&, const T20&, 
		const T21&, const T22&, const T23&, const T24&, const T25&, const T26&, ARM_result&, long );

	exportFunc26Args ( const T1& arg1, const T2& arg2, const T3& arg3, const T4& arg4, const T5& arg5, const T6& arg6, const T7& arg7, const T8& arg8, const T9& arg9, const T10& arg10, 
		const T11& arg11, const T12& arg12, const T13& arg13, const T14& arg14, const T15& arg15, const T16& arg16, const T17& arg17, const T18& arg18, const T19& arg19, const T20& arg20, 
		const T21& arg21, const T22& arg22, const T23& arg23, const T24& arg24, const T25& arg25, const T26& arg26, exportFunc26 func)
	: C_arg1(arg1), C_arg2(arg2), C_arg3(arg3), C_arg4(arg4), C_arg5(arg5), C_arg6(arg6), C_arg7(arg7), C_arg8(arg8), C_arg9(arg9), C_arg10(arg10), 	
	C_arg11(arg11), C_arg12(arg12), C_arg13(arg13), C_arg14(arg14), C_arg15(arg15), C_arg16(arg16), C_arg17(arg17), C_arg18(arg18), C_arg19(arg19), C_arg20(arg20), 
	C_arg21(arg21), C_arg22(arg22), C_arg23(arg23), C_arg24(arg24), C_arg25(arg25), C_arg26(arg26), itsFunc(func) {};

	long operator()( ARM_result& result, long objId )
	{ 
		return (*itsFunc)( C_arg1, C_arg2, C_arg3, C_arg4, C_arg5, C_arg6, C_arg7, C_arg8, C_arg9, C_arg10, C_arg11, C_arg12, C_arg13, C_arg14, C_arg15, C_arg16, 
		C_arg17, C_arg18, C_arg19, C_arg20, C_arg21, C_arg22, C_arg23, C_arg24, C_arg25, C_arg26, result, objId);
	}

private:

	const T1&	C_arg1;
	const T2&	C_arg2;
	const T3&	C_arg3;
	const T4&	C_arg4;
	const T5&	C_arg5;
	const T6&	C_arg6;
	const T7&	C_arg7;
	const T8&	C_arg8;
	const T9&	C_arg9;
	const T10&	C_arg10;
	const T11&	C_arg11;
	const T12&	C_arg12;
	const T13&	C_arg13;
	const T14&	C_arg14;
	const T15&	C_arg15;
	const T16&	C_arg16;
	const T17&	C_arg17;
	const T18&	C_arg18;
	const T19&	C_arg19;
	const T20&	C_arg20;
	const T21&	C_arg21;
	const T22&	C_arg22;
	const T23&	C_arg23;
	const T24&	C_arg24;
	const T25&	C_arg25;
	const T26&	C_arg26;

	exportFunc26 itsFunc;
};

///////////////////////////////////////
//// functor with 26 arguments  ///////
///////////////////////////////////////
template < typename T1, typename T2, typename T3,  typename T4, typename T5, typename T6, typename T7, typename T8,
typename T9, typename T10, typename T11, typename T12, typename T13, typename T14, typename T15, typename T16, 
typename T17, typename T18, typename T19, typename T20, typename T21, typename T22, typename T23, typename T24, 
typename T25, typename T26, typename T27 > class exportFunc27Args : public ARMResultLong2LongFunc
{
public: 	

	typedef long (*exportFunc27)( const T1&, const T2&,	const T3&, const T4&, const T5&, const T6&, const T7&, const T8& , const T9&, 
		const T10&, const T11&, const T12&, const T13&, const T14&, const T15&,const T16&, const T17&, const T18&, const T19&, const T20&, 
		const T21&, const T22&, const T23&, const T24&, const T25&, const T26&, const T27&, ARM_result&, long );

	exportFunc27Args ( const T1& arg1, const T2& arg2, const T3& arg3, const T4& arg4, const T5& arg5, const T6& arg6, const T7& arg7, const T8& arg8, const T9& arg9, const T10& arg10, 
		const T11& arg11, const T12& arg12, const T13& arg13, const T14& arg14, const T15& arg15, const T16& arg16, const T17& arg17, const T18& arg18, const T19& arg19, const T20& arg20, 
		const T21& arg21, const T22& arg22, const T23& arg23, const T24& arg24, const T25& arg25, const T26& arg26, const T27& arg27, exportFunc27 func)
	: C_arg1(arg1), C_arg2(arg2), C_arg3(arg3), C_arg4(arg4), C_arg5(arg5), C_arg6(arg6), C_arg7(arg7), C_arg8(arg8), C_arg9(arg9), C_arg10(arg10), 	
	C_arg11(arg11), C_arg12(arg12), C_arg13(arg13), C_arg14(arg14), C_arg15(arg15), C_arg16(arg16), C_arg17(arg17), C_arg18(arg18), C_arg19(arg19), C_arg20(arg20), 
	C_arg21(arg21), C_arg22(arg22), C_arg23(arg23), C_arg24(arg24), C_arg25(arg25), C_arg26(arg26), C_arg27(arg27), itsFunc(func) {};

	long operator()( ARM_result& result, long objId )
	{ 
		return (*itsFunc)( C_arg1, C_arg2, C_arg3, C_arg4, C_arg5, C_arg6, C_arg7, C_arg8, C_arg9, C_arg10, C_arg11, C_arg12, C_arg13, C_arg14, C_arg15, C_arg16, 
		C_arg17, C_arg18, C_arg19, C_arg20, C_arg21, C_arg22, C_arg23, C_arg24, C_arg25, C_arg26, C_arg27, result, objId);
	}

private:

	const T1&	C_arg1;
	const T2&	C_arg2;
	const T3&	C_arg3;
	const T4&	C_arg4;
	const T5&	C_arg5;
	const T6&	C_arg6;
	const T7&	C_arg7;
	const T8&	C_arg8;
	const T9&	C_arg9;
	const T10&	C_arg10;
	const T11&	C_arg11;
	const T12&	C_arg12;
	const T13&	C_arg13;
	const T14&	C_arg14;
	const T15&	C_arg15;
	const T16&	C_arg16;
	const T17&	C_arg17;
	const T18&	C_arg18;
	const T19&	C_arg19;
	const T20&	C_arg20;
	const T21&	C_arg21;
	const T22&	C_arg22;
	const T23&	C_arg23;
	const T24&	C_arg24;
	const T25&	C_arg25;
	const T26&	C_arg26;
	const T27&	C_arg27;

	exportFunc27 itsFunc;
};

///////////////////////////////////////
//// functor with 29 arguments  ///////
///////////////////////////////////////
template < typename T1, typename T2, typename T3,  typename T4, typename T5, typename T6, typename T7, typename T8,
typename T9, typename T10, typename T11, typename T12, typename T13, typename T14, typename T15, typename T16, 
typename T17, typename T18, typename T19, typename T20, typename T21, typename T22, typename T23, typename T24, 
typename T25, typename T26, typename T27, typename T28, typename T29 > class exportFunc29Args : public ARMResultLong2LongFunc
{
public: 	

	typedef long (*exportFunc29)( const T1&, const T2&,	const T3&, const T4&, const T5&, const T6&, const T7&, const T8& , const T9&, 
		const T10&, const T11&, const T12&, const T13&, const T14&, const T15&,const T16&, const T17&, const T18&, const T19&, const T20&, 
		const T21&, const T22&, const T23&, const T24&, const T25&, const T26&, const T27&, const T28&, const T29&, ARM_result&, long );

	exportFunc29Args ( const T1& arg1, const T2& arg2, const T3& arg3, const T4& arg4, const T5& arg5, const T6& arg6, const T7& arg7, const T8& arg8, const T9& arg9, const T10& arg10, 
		const T11& arg11, const T12& arg12, const T13& arg13, const T14& arg14, const T15& arg15, const T16& arg16, const T17& arg17, const T18& arg18, const T19& arg19, const T20& arg20, 
		const T21& arg21, const T22& arg22, const T23& arg23, const T24& arg24, const T25& arg25, const T26& arg26, const T27& arg27, const T28& arg28, const T29& arg29, exportFunc29 func)
	: C_arg1(arg1), C_arg2(arg2), C_arg3(arg3), C_arg4(arg4), C_arg5(arg5), C_arg6(arg6), C_arg7(arg7), C_arg8(arg8), C_arg9(arg9), C_arg10(arg10), 	
	C_arg11(arg11), C_arg12(arg12), C_arg13(arg13), C_arg14(arg14), C_arg15(arg15), C_arg16(arg16), C_arg17(arg17), C_arg18(arg18), C_arg19(arg19), C_arg20(arg20), 
	C_arg21(arg21), C_arg22(arg22), C_arg23(arg23), C_arg24(arg24), C_arg25(arg25), C_arg26(arg26), C_arg27(arg27), C_arg28(arg28), C_arg29(arg29), itsFunc(func) {};

	long operator()( ARM_result& result, long objId )
	{ 
		return (*itsFunc)( C_arg1, C_arg2, C_arg3, C_arg4, C_arg5, C_arg6, C_arg7, C_arg8, C_arg9, C_arg10, C_arg11, C_arg12, C_arg13, C_arg14, C_arg15, C_arg16, 
		C_arg17, C_arg18, C_arg19, C_arg20, C_arg21, C_arg22, C_arg23, C_arg24, C_arg25, C_arg26, C_arg27, C_arg28, C_arg29, result, objId);
	}

private:

	const T1&	C_arg1;
	const T2&	C_arg2;
	const T3&	C_arg3;
	const T4&	C_arg4;
	const T5&	C_arg5;
	const T6&	C_arg6;
	const T7&	C_arg7;
	const T8&	C_arg8;
	const T9&	C_arg9;
	const T10&	C_arg10;
	const T11&	C_arg11;
	const T12&	C_arg12;
	const T13&	C_arg13;
	const T14&	C_arg14;
	const T15&	C_arg15;
	const T16&	C_arg16;
	const T17&	C_arg17;
	const T18&	C_arg18;
	const T19&	C_arg19;
	const T20&	C_arg20;
	const T21&	C_arg21;
	const T22&	C_arg22;
	const T23&	C_arg23;
	const T24&	C_arg24;
	const T25&	C_arg25;
	const T26&	C_arg26;
	const T27&	C_arg27;
	const T28&	C_arg28;
	const T29&	C_arg29;

	exportFunc29 itsFunc;
};

///////////////////////////////////////
//// functor with 31 arguments
///////////////////////////////////////
template <typename T1, typename T2, typename T3,  typename T4, typename T5, typename T6, typename T7, typename T8, typename T9, typename T10, typename T11, typename T12, typename T13, typename T14, 
typename T15, typename T16, typename T17, typename T18, typename T19, typename T20, typename T21, typename T22, typename T23, typename T24, typename T25, typename T26, typename T27, typename T28, 
typename T29, typename T30, typename T31> class exportFunc31Args : public ARMResultLong2LongFunc
{
public: 	
	typedef long (*exportFunc31)( const T1&, const T2&,	const T3&, const T4&, const T5&, const T6&, const T7&, const T8& , const T9&, 
		const T10&, const T11&, const T12&, const T13&, const T14&, const T15&,const T16&, const T17&, const T18&, const T19&, const T20&, 
		const T21&, const T22&, const T23&, const T24&, const T25&, const T26&, const T27&, const T28&, const T29&, const T30&, 
		const T31&, ARM_result&, long );
	exportFunc31Args ( const T1& arg1, const T2& arg2, const T3& arg3, const T4& arg4, const T5& arg5, const T6& arg6, const T7& arg7, const T8& arg8, const T9& arg9, const T10& arg10, 
		const T11& arg11, const T12& arg12, const T13& arg13, const T14& arg14, const T15& arg15, const T16& arg16, const T17& arg17, const T18& arg18, const T19& arg19, const T20& arg20, 
		const T21& arg21, const T22& arg22, const T23& arg23, const T24& arg24, const T25& arg25, const T26& arg26, const T27& arg27, const T28& arg28, const T29& arg29, const T30& arg30, 
		const T31& arg31, exportFunc31 func)
	: C_arg1(arg1), C_arg2(arg2), C_arg3(arg3), C_arg4(arg4), C_arg5(arg5), C_arg6(arg6), C_arg7(arg7), C_arg8(arg8), C_arg9(arg9), C_arg10(arg10), 	
	C_arg11(arg11), C_arg12(arg12), C_arg13(arg13), C_arg14(arg14), C_arg15(arg15), C_arg16(arg16), C_arg17(arg17), C_arg18(arg18), C_arg19(arg19), C_arg20(arg20), 
	C_arg21(arg21), C_arg22(arg22), C_arg23(arg23), C_arg24(arg24), C_arg25(arg25), C_arg26(arg26), C_arg27(arg27), C_arg28(arg28), C_arg29(arg29), C_arg30(arg30), 
	C_arg31(arg31), itsFunc(func) {};
	long operator()( ARM_result& result, long objId ){ return (*itsFunc)( C_arg1, C_arg2, C_arg3, C_arg4, C_arg5, C_arg6, C_arg7, C_arg8, C_arg9, C_arg10, C_arg11, C_arg12, C_arg13, C_arg14, C_arg15, C_arg16, 
		C_arg17, C_arg18, C_arg19, C_arg20, C_arg21, C_arg22, C_arg23, C_arg24, C_arg25, C_arg26, C_arg27, C_arg28, C_arg29, C_arg30, C_arg31, result, objId);}
private:
	const T1& C_arg1;
	const T2& C_arg2;
	const T3& C_arg3;
	const T4& C_arg4;
	const T5& C_arg5;
	const T6& C_arg6;
	const T7& C_arg7;
	const T8& C_arg8;
	const T9& C_arg9;
	const T10& C_arg10;
	const T11& C_arg11;
	const T12& C_arg12;
	const T13& C_arg13;
	const T14& C_arg14;
	const T15& C_arg15;
	const T16& C_arg16;
	const T17& C_arg17;
	const T18& C_arg18;
	const T19& C_arg19;
	const T20& C_arg20;
	const T21& C_arg21;
	const T22& C_arg22;
	const T23& C_arg23;
	const T24& C_arg24;
	const T25& C_arg25;
	const T26& C_arg26;
	const T27& C_arg27;
	const T28& C_arg28;
	const T29& C_arg29;
	const T30& C_arg30;
	const T31& C_arg31;
	exportFunc31 itsFunc;
};

///////////////////////////////////////
//// functor with 32 arguments
///////////////////////////////////////
template <typename T1, typename T2, typename T3,  typename T4, typename T5, typename T6, typename T7, typename T8, typename T9, typename T10, typename T11, typename T12, typename T13, typename T14, 
typename T15, typename T16, typename T17, typename T18, typename T19, typename T20, typename T21, typename T22, typename T23, typename T24, typename T25, typename T26, typename T27, typename T28, 
typename T29, typename T30, typename T31, typename T32> class exportFunc32Args : public ARMResultLong2LongFunc
{
public: 	
	typedef long (*exportFunc32)( const T1&, const T2&,	const T3&, const T4&, const T5&, const T6&, const T7&, const T8& , const T9&, 
		const T10&, const T11&, const T12&, const T13&, const T14&, const T15&,const T16&, const T17&, const T18&, const T19&, const T20&, 
		const T21&, const T22&, const T23&, const T24&, const T25&, const T26&, const T27&, const T28&, const T29&, const T30&, 
		const T31&, const T32&, ARM_result&, long );
	exportFunc32Args ( const T1& arg1, const T2& arg2, const T3& arg3, const T4& arg4, const T5& arg5, const T6& arg6, const T7& arg7, const T8& arg8, const T9& arg9, const T10& arg10, 
		const T11& arg11, const T12& arg12, const T13& arg13, const T14& arg14, const T15& arg15, const T16& arg16, const T17& arg17, const T18& arg18, const T19& arg19, const T20& arg20, 
		const T21& arg21, const T22& arg22, const T23& arg23, const T24& arg24, const T25& arg25, const T26& arg26, const T27& arg27, const T28& arg28, const T29& arg29, const T30& arg30, 
		const T31& arg31, const T32& arg32, exportFunc32 func)
	: C_arg1(arg1), C_arg2(arg2), C_arg3(arg3), C_arg4(arg4), C_arg5(arg5), C_arg6(arg6), C_arg7(arg7), C_arg8(arg8), C_arg9(arg9), C_arg10(arg10), 	
	C_arg11(arg11), C_arg12(arg12), C_arg13(arg13), C_arg14(arg14), C_arg15(arg15), C_arg16(arg16), C_arg17(arg17), C_arg18(arg18), C_arg19(arg19), C_arg20(arg20), 
	C_arg21(arg21), C_arg22(arg22), C_arg23(arg23), C_arg24(arg24), C_arg25(arg25), C_arg26(arg26), C_arg27(arg27), C_arg28(arg28), C_arg29(arg29), C_arg30(arg30), 
	C_arg31(arg31), C_arg32(arg32), itsFunc(func) {};
	long operator()( ARM_result& result, long objId ){ return (*itsFunc)( C_arg1, C_arg2, C_arg3, C_arg4, C_arg5, C_arg6, C_arg7, C_arg8, C_arg9, C_arg10, C_arg11, C_arg12, C_arg13, C_arg14, C_arg15, C_arg16, 
		C_arg17, C_arg18, C_arg19, C_arg20, C_arg21, C_arg22, C_arg23, C_arg24, C_arg25, C_arg26, C_arg27, C_arg28, C_arg29, C_arg30, C_arg31, C_arg32, result, objId);}
private:
	const T1& C_arg1;
	const T2& C_arg2;
	const T3& C_arg3;
	const T4& C_arg4;
	const T5& C_arg5;
	const T6& C_arg6;
	const T7& C_arg7;
	const T8& C_arg8;
	const T9& C_arg9;
	const T10& C_arg10;
	const T11& C_arg11;
	const T12& C_arg12;
	const T13& C_arg13;
	const T14& C_arg14;
	const T15& C_arg15;
	const T16& C_arg16;
	const T17& C_arg17;
	const T18& C_arg18;
	const T19& C_arg19;
	const T20& C_arg20;
	const T21& C_arg21;
	const T22& C_arg22;
	const T23& C_arg23;
	const T24& C_arg24;
	const T25& C_arg25;
	const T26& C_arg26;
	const T27& C_arg27;
	const T28& C_arg28;
	const T29& C_arg29;
	const T30& C_arg30;
	const T31& C_arg31;
	const T32& C_arg32;
	exportFunc32 itsFunc;
};

///////////////////////////////////////
//// functor with 33 arguments
///////////////////////////////////////
template <typename T1, typename T2, typename T3,  typename T4, typename T5, typename T6, typename T7, typename T8, typename T9, typename T10, typename T11, typename T12, typename T13, typename T14, 
typename T15, typename T16, typename T17, typename T18, typename T19, typename T20, typename T21, typename T22, typename T23, typename T24, typename T25, typename T26, typename T27, typename T28, 
typename T29, typename T30, typename T31, typename T32, typename T33> class exportFunc33Args : public ARMResultLong2LongFunc
{
public: 	
	typedef long (*exportFunc33)( const T1&, const T2&,	const T3&, const T4&, const T5&, const T6&, const T7&, const T8& , const T9&, 
		const T10&, const T11&, const T12&, const T13&, const T14&, const T15&,const T16&, const T17&, const T18&, const T19&, const T20&, 
		const T21&, const T22&, const T23&, const T24&, const T25&, const T26&, const T27&, const T28&, const T29&, const T30&, 
		const T31&, const T32&, const T33&, ARM_result&, long );
	exportFunc33Args ( const T1& arg1, const T2& arg2, const T3& arg3, const T4& arg4, const T5& arg5, const T6& arg6, const T7& arg7, const T8& arg8, const T9& arg9, const T10& arg10, 
		const T11& arg11, const T12& arg12, const T13& arg13, const T14& arg14, const T15& arg15, const T16& arg16, const T17& arg17, const T18& arg18, const T19& arg19, const T20& arg20, 
		const T21& arg21, const T22& arg22, const T23& arg23, const T24& arg24, const T25& arg25, const T26& arg26, const T27& arg27, const T28& arg28, const T29& arg29, const T30& arg30, 
		const T31& arg31, const T32& arg32, const T33& arg33, exportFunc33 func)
	: C_arg1(arg1), C_arg2(arg2), C_arg3(arg3), C_arg4(arg4), C_arg5(arg5), C_arg6(arg6), C_arg7(arg7), C_arg8(arg8), C_arg9(arg9), C_arg10(arg10), 	
	C_arg11(arg11), C_arg12(arg12), C_arg13(arg13), C_arg14(arg14), C_arg15(arg15), C_arg16(arg16), C_arg17(arg17), C_arg18(arg18), C_arg19(arg19), C_arg20(arg20), 
	C_arg21(arg21), C_arg22(arg22), C_arg23(arg23), C_arg24(arg24), C_arg25(arg25), C_arg26(arg26), C_arg27(arg27), C_arg28(arg28), C_arg29(arg29), C_arg30(arg30), 
	C_arg31(arg31), C_arg32(arg32), C_arg33(arg33), itsFunc(func) {};
	long operator()( ARM_result& result, long objId ){ return (*itsFunc)( C_arg1, C_arg2, C_arg3, C_arg4, C_arg5, C_arg6, C_arg7, C_arg8, C_arg9, C_arg10, C_arg11, C_arg12, C_arg13, C_arg14, C_arg15, C_arg16, 
		C_arg17, C_arg18, C_arg19, C_arg20, C_arg21, C_arg22, C_arg23, C_arg24, C_arg25, C_arg26, C_arg27, C_arg28, C_arg29, C_arg30, C_arg31, C_arg32, C_arg33, result, objId);}
private:
	const T1& C_arg1;
	const T2& C_arg2;
	const T3& C_arg3;
	const T4& C_arg4;
	const T5& C_arg5;
	const T6& C_arg6;
	const T7& C_arg7;
	const T8& C_arg8;
	const T9& C_arg9;
	const T10& C_arg10;
	const T11& C_arg11;
	const T12& C_arg12;
	const T13& C_arg13;
	const T14& C_arg14;
	const T15& C_arg15;
	const T16& C_arg16;
	const T17& C_arg17;
	const T18& C_arg18;
	const T19& C_arg19;
	const T20& C_arg20;
	const T21& C_arg21;
	const T22& C_arg22;
	const T23& C_arg23;
	const T24& C_arg24;
	const T25& C_arg25;
	const T26& C_arg26;
	const T27& C_arg27;
	const T28& C_arg28;
	const T29& C_arg29;
	const T30& C_arg30;
	const T31& C_arg31;
	const T32& C_arg32;
	const T33& C_arg33;
	exportFunc33 itsFunc;
};

///////////////////////////////////////
//// functor with 34 arguments
///////////////////////////////////////
template <typename T1, typename T2, typename T3,  typename T4, typename T5, typename T6, typename T7, typename T8, typename T9, typename T10, typename T11, typename T12, typename T13, typename T14, 
typename T15, typename T16, typename T17, typename T18, typename T19, typename T20, typename T21, typename T22, typename T23, typename T24, typename T25, typename T26, typename T27, typename T28, 
typename T29, typename T30, typename T31, typename T32, typename T33, typename T34  > class exportFunc34Args : public ARMResultLong2LongFunc
{
public: 	
	typedef long (*exportFunc34)( const T1&, const T2&,	const T3&, const T4&, const T5&, const T6&, const T7&, const T8& , const T9&, 
		const T10&, const T11&, const T12&, const T13&, const T14&, const T15&,const T16&, const T17&, const T18&, const T19&, const T20&, 
		const T21&, const T22&, const T23&, const T24&, const T25&, const T26&, const T27&, const T28&, const T29&, const T30&, 
		const T31&, const T32&, const T33&, const T34&, ARM_result&, long );
	exportFunc34Args ( const T1& arg1, const T2& arg2, const T3& arg3, const T4& arg4, const T5& arg5, const T6& arg6, const T7& arg7, const T8& arg8, const T9& arg9, const T10& arg10, 
		const T11& arg11, const T12& arg12, const T13& arg13, const T14& arg14, const T15& arg15, const T16& arg16, const T17& arg17, const T18& arg18, const T19& arg19, const T20& arg20, 
		const T21& arg21, const T22& arg22, const T23& arg23, const T24& arg24, const T25& arg25, const T26& arg26, const T27& arg27, const T28& arg28, const T29& arg29, const T30& arg30, 
		const T31& arg31, const T32& arg32, const T33& arg33, const T34& arg34, exportFunc34 func)
	: C_arg1(arg1), C_arg2(arg2), C_arg3(arg3), C_arg4(arg4), C_arg5(arg5), C_arg6(arg6), C_arg7(arg7), C_arg8(arg8), C_arg9(arg9), C_arg10(arg10), 	
	C_arg11(arg11), C_arg12(arg12), C_arg13(arg13), C_arg14(arg14), C_arg15(arg15), C_arg16(arg16), C_arg17(arg17), C_arg18(arg18), C_arg19(arg19), C_arg20(arg20), 
	C_arg21(arg21), C_arg22(arg22), C_arg23(arg23), C_arg24(arg24), C_arg25(arg25), C_arg26(arg26), C_arg27(arg27), C_arg28(arg28), C_arg29(arg29), C_arg30(arg30), 
	C_arg31(arg31), C_arg32(arg32), C_arg33(arg33), C_arg34(arg34),  itsFunc(func) {};
	long operator()( ARM_result& result, long objId ){ return (*itsFunc)( C_arg1, C_arg2, C_arg3, C_arg4, C_arg5, C_arg6, C_arg7, C_arg8, C_arg9, C_arg10, C_arg11, C_arg12, C_arg13, C_arg14, C_arg15, C_arg16, 
		C_arg17, C_arg18, C_arg19, C_arg20, C_arg21, C_arg22, C_arg23, C_arg24, C_arg25, C_arg26, C_arg27, C_arg28, C_arg29, C_arg30, C_arg31, C_arg32, C_arg33, C_arg34, result, objId);}
private:
	const T1& C_arg1;
	const T2& C_arg2;
	const T3& C_arg3;
	const T4& C_arg4;
	const T5& C_arg5;
	const T6& C_arg6;
	const T7& C_arg7;
	const T8& C_arg8;
	const T9& C_arg9;
	const T10& C_arg10;
	const T11& C_arg11;
	const T12& C_arg12;
	const T13& C_arg13;
	const T14& C_arg14;
	const T15& C_arg15;
	const T16& C_arg16;
	const T17& C_arg17;
	const T18& C_arg18;
	const T19& C_arg19;
	const T20& C_arg20;
	const T21& C_arg21;
	const T22& C_arg22;
	const T23& C_arg23;
	const T24& C_arg24;
	const T25& C_arg25;
	const T26& C_arg26;
	const T27& C_arg27;
	const T28& C_arg28;
	const T29& C_arg29;
	const T30& C_arg30;
	const T31& C_arg31;
	const T32& C_arg32;
	const T33& C_arg33;
	const T34& C_arg34;
	exportFunc34 itsFunc;
};

///////////////////////////////////////
//// functor with 35 arguments
///////////////////////////////////////
template <
typename T1,  typename T2,  typename T3,  typename T4,  typename T5,  typename T6,  typename T7,  typename T8,  typename T9,  typename T10, 
typename T11, typename T12, typename T13, typename T14, typename T15, typename T16, typename T17, typename T18, typename T19, typename T20, 
typename T21, typename T22, typename T23, typename T24, typename T25, typename T26, typename T27, typename T28, typename T29, typename T30, 
typename T31, typename T32, typename T33, typename T34, typename T35 > 
class exportFunc35Args : public ARMResultLong2LongFunc
{
public: 	
	typedef long (*exportFunc35)( 
		const T1&,  const T2&,	const T3&,  const T4&,  const T5&,  const T6&,  const T7&,  const T8& , const T9&,  const T10&, 
		const T11&, const T12&, const T13&, const T14&, const T15&, const T16&, const T17&, const T18&, const T19&, const T20&, 
		const T21&, const T22&, const T23&, const T24&, const T25&, const T26&, const T27&, const T28&, const T29&, const T30&, 
		const T31&, const T32&, const T33&, const T34&, const T35&, ARM_result&, long );
	exportFunc35Args ( 
		const T1& arg1,   const T2& arg2,   const T3& arg3,   const T4& arg4,   const T5& arg5,   const T6& arg6,   const T7& arg7,   const T8& arg8,   const T9& arg9,   const T10& arg10, 
		const T11& arg11, const T12& arg12, const T13& arg13, const T14& arg14, const T15& arg15, const T16& arg16, const T17& arg17, const T18& arg18, const T19& arg19, const T20& arg20, 
		const T21& arg21, const T22& arg22, const T23& arg23, const T24& arg24, const T25& arg25, const T26& arg26, const T27& arg27, const T28& arg28, const T29& arg29, const T30& arg30, 
		const T31& arg31, const T32& arg32, const T33& arg33, const T34& arg34, const T35& arg35, exportFunc35 func)
	: C_arg1(arg1), C_arg2(arg2), C_arg3(arg3), C_arg4(arg4), C_arg5(arg5), C_arg6(arg6), C_arg7(arg7), C_arg8(arg8), C_arg9(arg9), C_arg10(arg10), 	
	C_arg11(arg11), C_arg12(arg12), C_arg13(arg13), C_arg14(arg14), C_arg15(arg15), C_arg16(arg16), C_arg17(arg17), C_arg18(arg18), C_arg19(arg19), C_arg20(arg20), 
	C_arg21(arg21), C_arg22(arg22), C_arg23(arg23), C_arg24(arg24), C_arg25(arg25), C_arg26(arg26), C_arg27(arg27), C_arg28(arg28), C_arg29(arg29), C_arg30(arg30), 
	C_arg31(arg31), C_arg32(arg32), C_arg33(arg33), C_arg34(arg34), C_arg35(arg35), itsFunc(func) {};
	long operator()( ARM_result& result, long objId ){ return (*itsFunc)( 
		C_arg1, C_arg2, C_arg3, C_arg4, C_arg5, C_arg6, C_arg7, C_arg8, C_arg9, C_arg10, 
		C_arg11, C_arg12, C_arg13, C_arg14, C_arg15, C_arg16, C_arg17, C_arg18, C_arg19, C_arg20, 
		C_arg21, C_arg22, C_arg23, C_arg24, C_arg25, C_arg26, C_arg27, C_arg28, C_arg29, C_arg30, 
		C_arg31, C_arg32, C_arg33, C_arg34, C_arg35,result, objId);}
private:
	const T1& C_arg1;
	const T2& C_arg2;
	const T3& C_arg3;
	const T4& C_arg4;
	const T5& C_arg5;
	const T6& C_arg6;
	const T7& C_arg7;
	const T8& C_arg8;
	const T9& C_arg9;
	const T10& C_arg10;
	const T11& C_arg11;
	const T12& C_arg12;
	const T13& C_arg13;
	const T14& C_arg14;
	const T15& C_arg15;
	const T16& C_arg16;
	const T17& C_arg17;
	const T18& C_arg18;
	const T19& C_arg19;
	const T20& C_arg20;
	const T21& C_arg21;
	const T22& C_arg22;
	const T23& C_arg23;
	const T24& C_arg24;
	const T25& C_arg25;
	const T26& C_arg26;
	const T27& C_arg27;
	const T28& C_arg28;
	const T29& C_arg29;
	const T30& C_arg30;
	const T31& C_arg31;
	const T32& C_arg32;
	const T33& C_arg33;
	const T34& C_arg34;
	const T35& C_arg35;
	exportFunc35 itsFunc;
};


///////////////////////////////////////
//// functor with 36 arguments
///////////////////////////////////////
template <
typename T1,  typename T2,  typename T3,  typename T4,  typename T5,  typename T6,  typename T7,  typename T8,  typename T9,  typename T10, 
typename T11, typename T12, typename T13, typename T14, typename T15, typename T16, typename T17, typename T18, typename T19, typename T20, 
typename T21, typename T22, typename T23, typename T24, typename T25, typename T26, typename T27, typename T28, typename T29, typename T30, 
typename T31, typename T32, typename T33, typename T34, typename T35, typename T36 > 
class exportFunc36Args : public ARMResultLong2LongFunc
{
public: 	
	typedef long (*exportFunc36)( 
		const T1&,  const T2&,	const T3&,  const T4&,  const T5&,  const T6&,  const T7&,  const T8& , const T9&,  const T10&, 
		const T11&, const T12&, const T13&, const T14&, const T15&, const T16&, const T17&, const T18&, const T19&, const T20&, 
		const T21&, const T22&, const T23&, const T24&, const T25&, const T26&, const T27&, const T28&, const T29&, const T30&, 
		const T31&, const T32&, const T33&, const T34&, const T35&, const T36&, ARM_result&, long );
	exportFunc36Args ( 
		const T1& arg1,   const T2& arg2,   const T3& arg3,   const T4& arg4,   const T5& arg5,   const T6& arg6,   const T7& arg7,   const T8& arg8,   const T9& arg9,   const T10& arg10, 
		const T11& arg11, const T12& arg12, const T13& arg13, const T14& arg14, const T15& arg15, const T16& arg16, const T17& arg17, const T18& arg18, const T19& arg19, const T20& arg20, 
		const T21& arg21, const T22& arg22, const T23& arg23, const T24& arg24, const T25& arg25, const T26& arg26, const T27& arg27, const T28& arg28, const T29& arg29, const T30& arg30, 
		const T31& arg31, const T32& arg32, const T33& arg33, const T34& arg34, const T35& arg35, const T36& arg36, exportFunc36 func)
	: C_arg1(arg1), C_arg2(arg2), C_arg3(arg3), C_arg4(arg4), C_arg5(arg5), C_arg6(arg6), C_arg7(arg7), C_arg8(arg8), C_arg9(arg9), C_arg10(arg10), 	
	C_arg11(arg11), C_arg12(arg12), C_arg13(arg13), C_arg14(arg14), C_arg15(arg15), C_arg16(arg16), C_arg17(arg17), C_arg18(arg18), C_arg19(arg19), C_arg20(arg20), 
	C_arg21(arg21), C_arg22(arg22), C_arg23(arg23), C_arg24(arg24), C_arg25(arg25), C_arg26(arg26), C_arg27(arg27), C_arg28(arg28), C_arg29(arg29), C_arg30(arg30), 
	C_arg31(arg31), C_arg32(arg32), C_arg33(arg33), C_arg34(arg34), C_arg35(arg35),C_arg36(arg36), itsFunc(func) {};
	long operator()( ARM_result& result, long objId ){ return (*itsFunc)( 
		C_arg1, C_arg2, C_arg3, C_arg4, C_arg5, C_arg6, C_arg7, C_arg8, C_arg9, C_arg10, 
		C_arg11, C_arg12, C_arg13, C_arg14, C_arg15, C_arg16, C_arg17, C_arg18, C_arg19, C_arg20, 
		C_arg21, C_arg22, C_arg23, C_arg24, C_arg25, C_arg26, C_arg27, C_arg28, C_arg29, C_arg30, 
		C_arg31, C_arg32, C_arg33, C_arg34, C_arg35, C_arg36,result, objId);}
private:
	const T1& C_arg1;
	const T2& C_arg2;
	const T3& C_arg3;
	const T4& C_arg4;
	const T5& C_arg5;
	const T6& C_arg6;
	const T7& C_arg7;
	const T8& C_arg8;
	const T9& C_arg9;
	const T10& C_arg10;
	const T11& C_arg11;
	const T12& C_arg12;
	const T13& C_arg13;
	const T14& C_arg14;
	const T15& C_arg15;
	const T16& C_arg16;
	const T17& C_arg17;
	const T18& C_arg18;
	const T19& C_arg19;
	const T20& C_arg20;
	const T21& C_arg21;
	const T22& C_arg22;
	const T23& C_arg23;
	const T24& C_arg24;
	const T25& C_arg25;
	const T26& C_arg26;
	const T27& C_arg27;
	const T28& C_arg28;
	const T29& C_arg29;
	const T30& C_arg30;
	const T31& C_arg31;
	const T32& C_arg32;
	const T33& C_arg33;
	const T34& C_arg34;
	const T35& C_arg35;
    const T36& C_arg36;
	exportFunc36 itsFunc;
};


///////////////////////////////////////
//// functor with 38 arguments
///////////////////////////////////////
template <
typename T1,  typename T2,  typename T3,  typename T4,  typename T5,  typename T6,  typename T7,  typename T8,  typename T9,  typename T10, 
typename T11, typename T12, typename T13, typename T14, typename T15, typename T16, typename T17, typename T18, typename T19, typename T20, 
typename T21, typename T22, typename T23, typename T24, typename T25, typename T26, typename T27, typename T28, typename T29, typename T30, 
typename T31, typename T32, typename T33, typename T34, typename T35, typename T36, typename T37> 
class exportFunc37Args : public ARMResultLong2LongFunc
{
public: 	
	typedef long (*exportFunc37)( 
		const T1&,  const T2&,	const T3&,  const T4&,  const T5&,  const T6&,  const T7&,  const T8& , const T9&,  const T10&, 
		const T11&, const T12&, const T13&, const T14&, const T15&, const T16&, const T17&, const T18&, const T19&, const T20&, 
		const T21&, const T22&, const T23&, const T24&, const T25&, const T26&, const T27&, const T28&, const T29&, const T30&, 
		const T31&, const T32&, const T33&, const T34&, const T35&, const T36&, const T37&, ARM_result&, long );
	exportFunc37Args ( 
		const T1& arg1,   const T2& arg2,   const T3& arg3,   const T4& arg4,   const T5& arg5,   const T6& arg6,   const T7& arg7,   const T8& arg8,   const T9& arg9,   const T10& arg10, 
		const T11& arg11, const T12& arg12, const T13& arg13, const T14& arg14, const T15& arg15, const T16& arg16, const T17& arg17, const T18& arg18, const T19& arg19, const T20& arg20, 
		const T21& arg21, const T22& arg22, const T23& arg23, const T24& arg24, const T25& arg25, const T26& arg26, const T27& arg27, const T28& arg28, const T29& arg29, const T30& arg30, 
		const T31& arg31, const T32& arg32, const T33& arg33, const T34& arg34, const T35& arg35, const T36& arg36, const T37& arg37, exportFunc37 func)
	: C_arg1(arg1), C_arg2(arg2), C_arg3(arg3), C_arg4(arg4), C_arg5(arg5), C_arg6(arg6), C_arg7(arg7), C_arg8(arg8), C_arg9(arg9), C_arg10(arg10), 	
	C_arg11(arg11), C_arg12(arg12), C_arg13(arg13), C_arg14(arg14), C_arg15(arg15), C_arg16(arg16), C_arg17(arg17), C_arg18(arg18), C_arg19(arg19), C_arg20(arg20), 
	C_arg21(arg21), C_arg22(arg22), C_arg23(arg23), C_arg24(arg24), C_arg25(arg25), C_arg26(arg26), C_arg27(arg27), C_arg28(arg28), C_arg29(arg29), C_arg30(arg30), 
	C_arg31(arg31), C_arg32(arg32), C_arg33(arg33), C_arg34(arg34), C_arg35(arg35),C_arg36(arg36),  C_arg37(arg37), itsFunc(func) {};
	long operator()( ARM_result& result, long objId ){ return (*itsFunc)( 
		C_arg1, C_arg2, C_arg3, C_arg4, C_arg5, C_arg6, C_arg7, C_arg8, C_arg9, C_arg10, 
		C_arg11, C_arg12, C_arg13, C_arg14, C_arg15, C_arg16, C_arg17, C_arg18, C_arg19, C_arg20, 
		C_arg21, C_arg22, C_arg23, C_arg24, C_arg25, C_arg26, C_arg27, C_arg28, C_arg29, C_arg30, 
		C_arg31, C_arg32, C_arg33, C_arg34, C_arg35, C_arg36, C_arg37, result, objId);}
private:
	const T1& C_arg1;
	const T2& C_arg2;
	const T3& C_arg3;
	const T4& C_arg4;
	const T5& C_arg5;
	const T6& C_arg6;
	const T7& C_arg7;
	const T8& C_arg8;
	const T9& C_arg9;
	const T10& C_arg10;
	const T11& C_arg11;
	const T12& C_arg12;
	const T13& C_arg13;
	const T14& C_arg14;
	const T15& C_arg15;
	const T16& C_arg16;
	const T17& C_arg17;
	const T18& C_arg18;
	const T19& C_arg19;
	const T20& C_arg20;
	const T21& C_arg21;
	const T22& C_arg22;
	const T23& C_arg23;
	const T24& C_arg24;
	const T25& C_arg25;
	const T26& C_arg26;
	const T27& C_arg27;
	const T28& C_arg28;
	const T29& C_arg29;
	const T30& C_arg30;
	const T31& C_arg31;
	const T32& C_arg32;
	const T33& C_arg33;
	const T34& C_arg34;
	const T35& C_arg35;
    const T36& C_arg36;
    const T37& C_arg37;
	exportFunc37 itsFunc;
};




///////////////////////////////////////
//// functor with 38 arguments
///////////////////////////////////////
template <
typename T1,  typename T2,  typename T3,  typename T4,  typename T5,  typename T6,  typename T7,  typename T8,  typename T9,  typename T10, 
typename T11, typename T12, typename T13, typename T14, typename T15, typename T16, typename T17, typename T18, typename T19, typename T20, 
typename T21, typename T22, typename T23, typename T24, typename T25, typename T26, typename T27, typename T28, typename T29, typename T30, 
typename T31, typename T32, typename T33, typename T34, typename T35, typename T36, typename T37, typename T38 > 
class exportFunc38Args : public ARMResultLong2LongFunc
{
public: 	
	typedef long (*exportFunc38)( 
		const T1&,  const T2&,	const T3&,  const T4&,  const T5&,  const T6&,  const T7&,  const T8& , const T9&,  const T10&, 
		const T11&, const T12&, const T13&, const T14&, const T15&, const T16&, const T17&, const T18&, const T19&, const T20&, 
		const T21&, const T22&, const T23&, const T24&, const T25&, const T26&, const T27&, const T28&, const T29&, const T30&, 
		const T31&, const T32&, const T33&, const T34&, const T35&, const T36&, const T37&, const T38&, ARM_result&, long );
	exportFunc38Args ( 
		const T1& arg1,   const T2& arg2,   const T3& arg3,   const T4& arg4,   const T5& arg5,   const T6& arg6,   const T7& arg7,   const T8& arg8,   const T9& arg9,   const T10& arg10, 
		const T11& arg11, const T12& arg12, const T13& arg13, const T14& arg14, const T15& arg15, const T16& arg16, const T17& arg17, const T18& arg18, const T19& arg19, const T20& arg20, 
		const T21& arg21, const T22& arg22, const T23& arg23, const T24& arg24, const T25& arg25, const T26& arg26, const T27& arg27, const T28& arg28, const T29& arg29, const T30& arg30, 
		const T31& arg31, const T32& arg32, const T33& arg33, const T34& arg34, const T35& arg35, const T36& arg36, const T37& arg37, const T38& arg38, exportFunc38 func)
	: C_arg1(arg1), C_arg2(arg2), C_arg3(arg3), C_arg4(arg4), C_arg5(arg5), C_arg6(arg6), C_arg7(arg7), C_arg8(arg8), C_arg9(arg9), C_arg10(arg10), 	
	C_arg11(arg11), C_arg12(arg12), C_arg13(arg13), C_arg14(arg14), C_arg15(arg15), C_arg16(arg16), C_arg17(arg17), C_arg18(arg18), C_arg19(arg19), C_arg20(arg20), 
	C_arg21(arg21), C_arg22(arg22), C_arg23(arg23), C_arg24(arg24), C_arg25(arg25), C_arg26(arg26), C_arg27(arg27), C_arg28(arg28), C_arg29(arg29), C_arg30(arg30), 
	C_arg31(arg31), C_arg32(arg32), C_arg33(arg33), C_arg34(arg34), C_arg35(arg35),C_arg36(arg36),  C_arg37(arg37),C_arg38(arg38) , itsFunc(func) {};
	long operator()( ARM_result& result, long objId ){ return (*itsFunc)( 
		C_arg1, C_arg2, C_arg3, C_arg4, C_arg5, C_arg6, C_arg7, C_arg8, C_arg9, C_arg10, 
		C_arg11, C_arg12, C_arg13, C_arg14, C_arg15, C_arg16, C_arg17, C_arg18, C_arg19, C_arg20, 
		C_arg21, C_arg22, C_arg23, C_arg24, C_arg25, C_arg26, C_arg27, C_arg28, C_arg29, C_arg30, 
		C_arg31, C_arg32, C_arg33, C_arg34, C_arg35, C_arg36, C_arg37, C_arg38, result, objId);}
private:
	const T1& C_arg1;
	const T2& C_arg2;
	const T3& C_arg3;
	const T4& C_arg4;
	const T5& C_arg5;
	const T6& C_arg6;
	const T7& C_arg7;
	const T8& C_arg8;
	const T9& C_arg9;
	const T10& C_arg10;
	const T11& C_arg11;
	const T12& C_arg12;
	const T13& C_arg13;
	const T14& C_arg14;
	const T15& C_arg15;
	const T16& C_arg16;
	const T17& C_arg17;
	const T18& C_arg18;
	const T19& C_arg19;
	const T20& C_arg20;
	const T21& C_arg21;
	const T22& C_arg22;
	const T23& C_arg23;
	const T24& C_arg24;
	const T25& C_arg25;
	const T26& C_arg26;
	const T27& C_arg27;
	const T28& C_arg28;
	const T29& C_arg29;
	const T30& C_arg30;
	const T31& C_arg31;
	const T32& C_arg32;
	const T33& C_arg33;
	const T34& C_arg34;
	const T35& C_arg35;
    const T36& C_arg36;
    const T37& C_arg37;
    const T38& C_arg38;
	exportFunc38 itsFunc;
};


///////////////////////////////////////
//// functor with 39 arguments
///////////////////////////////////////
template <
typename T1,  typename T2,  typename T3,  typename T4,  typename T5,  typename T6,  typename T7,  typename T8,  typename T9,  typename T10, 
typename T11, typename T12, typename T13, typename T14, typename T15, typename T16, typename T17, typename T18, typename T19, typename T20, 
typename T21, typename T22, typename T23, typename T24, typename T25, typename T26, typename T27, typename T28, typename T29, typename T30, 
typename T31, typename T32, typename T33, typename T34, typename T35, typename T36, typename T37, typename T38, typename T39 > 
class exportFunc39Args : public ARMResultLong2LongFunc
{
public: 	
	typedef long (*exportFunc39)( 
		const T1&,  const T2&,	const T3&,  const T4&,  const T5&,  const T6&,  const T7&,  const T8& , const T9&,  const T10&, 
		const T11&, const T12&, const T13&, const T14&, const T15&, const T16&, const T17&, const T18&, const T19&, const T20&, 
		const T21&, const T22&, const T23&, const T24&, const T25&, const T26&, const T27&, const T28&, const T29&, const T30&, 
		const T31&, const T32&, const T33&, const T34&, const T35&, const T36&, const T37&, const T38&, const T39&, ARM_result&, long );
	exportFunc39Args ( 
		const T1& arg1,   const T2& arg2,   const T3& arg3,   const T4& arg4,   const T5& arg5,   const T6& arg6,   const T7& arg7,   const T8& arg8,   const T9& arg9,   const T10& arg10, 
		const T11& arg11, const T12& arg12, const T13& arg13, const T14& arg14, const T15& arg15, const T16& arg16, const T17& arg17, const T18& arg18, const T19& arg19, const T20& arg20, 
		const T21& arg21, const T22& arg22, const T23& arg23, const T24& arg24, const T25& arg25, const T26& arg26, const T27& arg27, const T28& arg28, const T29& arg29, const T30& arg30, 
		const T31& arg31, const T32& arg32, const T33& arg33, const T34& arg34, const T35& arg35, const T36& arg36, const T37& arg37, const T38& arg38, const T39& arg39, exportFunc39 func)
	: C_arg1(arg1), C_arg2(arg2), C_arg3(arg3), C_arg4(arg4), C_arg5(arg5), C_arg6(arg6), C_arg7(arg7), C_arg8(arg8), C_arg9(arg9), C_arg10(arg10), 	
	C_arg11(arg11), C_arg12(arg12), C_arg13(arg13), C_arg14(arg14), C_arg15(arg15), C_arg16(arg16), C_arg17(arg17), C_arg18(arg18), C_arg19(arg19), C_arg20(arg20), 
	C_arg21(arg21), C_arg22(arg22), C_arg23(arg23), C_arg24(arg24), C_arg25(arg25), C_arg26(arg26), C_arg27(arg27), C_arg28(arg28), C_arg29(arg29), C_arg30(arg30), 
	C_arg31(arg31), C_arg32(arg32), C_arg33(arg33), C_arg34(arg34), C_arg35(arg35), C_arg36(arg36), C_arg37(arg37), C_arg38(arg38), C_arg39(arg39), itsFunc(func) {};
	long operator()( ARM_result& result, long objId ){ return (*itsFunc)( 
		C_arg1, C_arg2, C_arg3, C_arg4, C_arg5, C_arg6, C_arg7, C_arg8, C_arg9, C_arg10, 
		C_arg11, C_arg12, C_arg13, C_arg14, C_arg15, C_arg16, C_arg17, C_arg18, C_arg19, C_arg20, 
		C_arg21, C_arg22, C_arg23, C_arg24, C_arg25, C_arg26, C_arg27, C_arg28, C_arg29, C_arg30, 
		C_arg31, C_arg32, C_arg33, C_arg34, C_arg35, C_arg36, C_arg37, C_arg38, C_arg39, result, objId);}
private:
	const T1& C_arg1;
	const T2& C_arg2;
	const T3& C_arg3;
	const T4& C_arg4;
	const T5& C_arg5;
	const T6& C_arg6;
	const T7& C_arg7;
	const T8& C_arg8;
	const T9& C_arg9;
	const T10& C_arg10;
	const T11& C_arg11;
	const T12& C_arg12;
	const T13& C_arg13;
	const T14& C_arg14;
	const T15& C_arg15;
	const T16& C_arg16;
	const T17& C_arg17;
	const T18& C_arg18;
	const T19& C_arg19;
	const T20& C_arg20;
	const T21& C_arg21;
	const T22& C_arg22;
	const T23& C_arg23;
	const T24& C_arg24;
	const T25& C_arg25;
	const T26& C_arg26;
	const T27& C_arg27;
	const T28& C_arg28;
	const T29& C_arg29;
	const T30& C_arg30;
	const T31& C_arg31;
	const T32& C_arg32;
	const T33& C_arg33;
	const T34& C_arg34;
	const T35& C_arg35;
	const T36& C_arg36;
	const T37& C_arg37;
	const T38& C_arg38;
	const T39& C_arg39;
	exportFunc39 itsFunc;
};

///////////////////////////////////////
//// functor with 40 arguments
///////////////////////////////////////
template <
typename T1,  typename T2,  typename T3,  typename T4,  typename T5,  typename T6,  typename T7,  typename T8,  typename T9,  typename T10, 
typename T11, typename T12, typename T13, typename T14, typename T15, typename T16, typename T17, typename T18, typename T19, typename T20, 
typename T21, typename T22, typename T23, typename T24, typename T25, typename T26, typename T27, typename T28, typename T29, typename T30, 
typename T31, typename T32, typename T33, typename T34, typename T35, typename T36, typename T37, typename T38, typename T39, typename T40 > 
class exportFunc40Args : public ARMResultLong2LongFunc
{
public: 	
	typedef long (*exportFunc40)( 
		const T1&,  const T2&,	const T3&,  const T4&,  const T5&,  const T6&,  const T7&,  const T8& , const T9&,  const T10&, 
		const T11&, const T12&, const T13&, const T14&, const T15&, const T16&, const T17&, const T18&, const T19&, const T20&, 
		const T21&, const T22&, const T23&, const T24&, const T25&, const T26&, const T27&, const T28&, const T29&, const T30&, 
		const T31&, const T32&, const T33&, const T34&, const T35&, const T36&, const T37&, const T38&, const T39&, const T40&, 
		ARM_result&, long );
	exportFunc40Args ( 
		const T1& arg1,   const T2& arg2,   const T3& arg3,   const T4& arg4,   const T5& arg5,   const T6& arg6,   const T7& arg7,   const T8& arg8,   const T9& arg9,   const T10& arg10, 
		const T11& arg11, const T12& arg12, const T13& arg13, const T14& arg14, const T15& arg15, const T16& arg16, const T17& arg17, const T18& arg18, const T19& arg19, const T20& arg20, 
		const T21& arg21, const T22& arg22, const T23& arg23, const T24& arg24, const T25& arg25, const T26& arg26, const T27& arg27, const T28& arg28, const T29& arg29, const T30& arg30, 
		const T31& arg31, const T32& arg32, const T33& arg33, const T34& arg34, const T35& arg35, const T36& arg36, const T37& arg37, const T38& arg38, const T39& arg39, const T40& arg40,
		exportFunc40 func)
	: C_arg1(arg1), C_arg2(arg2), C_arg3(arg3), C_arg4(arg4), C_arg5(arg5), C_arg6(arg6), C_arg7(arg7), C_arg8(arg8), C_arg9(arg9), C_arg10(arg10), 	
	C_arg11(arg11), C_arg12(arg12), C_arg13(arg13), C_arg14(arg14), C_arg15(arg15), C_arg16(arg16), C_arg17(arg17), C_arg18(arg18), C_arg19(arg19), C_arg20(arg20), 
	C_arg21(arg21), C_arg22(arg22), C_arg23(arg23), C_arg24(arg24), C_arg25(arg25), C_arg26(arg26), C_arg27(arg27), C_arg28(arg28), C_arg29(arg29), C_arg30(arg30), 
	C_arg31(arg31), C_arg32(arg32), C_arg33(arg33), C_arg34(arg34), C_arg35(arg35),C_arg36(arg36),  C_arg37(arg37), C_arg38(arg38), C_arg39(arg39), C_arg40(arg40),
	itsFunc(func) {};
	long operator()( ARM_result& result, long objId ){ return (*itsFunc)( 
		C_arg1, C_arg2, C_arg3, C_arg4, C_arg5, C_arg6, C_arg7, C_arg8, C_arg9, C_arg10, 
		C_arg11, C_arg12, C_arg13, C_arg14, C_arg15, C_arg16, C_arg17, C_arg18, C_arg19, C_arg20, 
		C_arg21, C_arg22, C_arg23, C_arg24, C_arg25, C_arg26, C_arg27, C_arg28, C_arg29, C_arg30, 
		C_arg31, C_arg32, C_arg33, C_arg34, C_arg35, C_arg36, C_arg37, C_arg38, C_arg39, C_arg40, result, objId);}
private:
	const T1& C_arg1;
	const T2& C_arg2;
	const T3& C_arg3;
	const T4& C_arg4;
	const T5& C_arg5;
	const T6& C_arg6;
	const T7& C_arg7;
	const T8& C_arg8;
	const T9& C_arg9;
	const T10& C_arg10;
	const T11& C_arg11;
	const T12& C_arg12;
	const T13& C_arg13;
	const T14& C_arg14;
	const T15& C_arg15;
	const T16& C_arg16;
	const T17& C_arg17;
	const T18& C_arg18;
	const T19& C_arg19;
	const T20& C_arg20;
	const T21& C_arg21;
	const T22& C_arg22;
	const T23& C_arg23;
	const T24& C_arg24;
	const T25& C_arg25;
	const T26& C_arg26;
	const T27& C_arg27;
	const T28& C_arg28;
	const T29& C_arg29;
	const T30& C_arg30;
	const T31& C_arg31;
	const T32& C_arg32;
	const T33& C_arg33;
	const T34& C_arg34;
	const T35& C_arg35;
    const T36& C_arg36;
    const T37& C_arg37;
    const T38& C_arg38;
    const T39& C_arg39;
    const T40& C_arg40;
	exportFunc40 itsFunc;
};

///////////////////////////////////////
//// functor with 41 arguments
///////////////////////////////////////
template <
typename T1,  typename T2,  typename T3,  typename T4,  typename T5,  typename T6,  typename T7,  typename T8,  typename T9,  typename T10, 
typename T11, typename T12, typename T13, typename T14, typename T15, typename T16, typename T17, typename T18, typename T19, typename T20, 
typename T21, typename T22, typename T23, typename T24, typename T25, typename T26, typename T27, typename T28, typename T29, typename T30, 
typename T31, typename T32, typename T33, typename T34, typename T35, typename T36, typename T37, typename T38, typename T39, typename T40, 
typename T41 > 
class exportFunc41Args : public ARMResultLong2LongFunc
{
public: 	
	typedef long (*exportFunc41)( 
		const T1&,  const T2&,	const T3&,  const T4&,  const T5&,  const T6&,  const T7&,  const T8& , const T9&,  const T10&, 
		const T11&, const T12&, const T13&, const T14&, const T15&, const T16&, const T17&, const T18&, const T19&, const T20&, 
		const T21&, const T22&, const T23&, const T24&, const T25&, const T26&, const T27&, const T28&, const T29&, const T30&, 
		const T31&, const T32&, const T33&, const T34&, const T35&, const T36&, const T37&, const T38&, const T39&, const T40&, 
		const T41&, ARM_result&, long );
	exportFunc41Args ( 
		const T1& arg1,   const T2& arg2,   const T3& arg3,   const T4& arg4,   const T5& arg5,   const T6& arg6,   const T7& arg7,   const T8& arg8,   const T9& arg9,   const T10& arg10, 
		const T11& arg11, const T12& arg12, const T13& arg13, const T14& arg14, const T15& arg15, const T16& arg16, const T17& arg17, const T18& arg18, const T19& arg19, const T20& arg20, 
		const T21& arg21, const T22& arg22, const T23& arg23, const T24& arg24, const T25& arg25, const T26& arg26, const T27& arg27, const T28& arg28, const T29& arg29, const T30& arg30, 
		const T31& arg31, const T32& arg32, const T33& arg33, const T34& arg34, const T35& arg35, const T36& arg36, const T37& arg37, const T38& arg38, const T39& arg39, const T40& arg40, 
		const T41& arg41, exportFunc41 func)
	: C_arg1(arg1), C_arg2(arg2),	C_arg3(arg3),	C_arg4(arg4),	C_arg5(arg5),	C_arg6(arg6),	C_arg7(arg7),	C_arg8(arg8),	C_arg9(arg9),	C_arg10(arg10), 	
	C_arg11(arg11), C_arg12(arg12), C_arg13(arg13), C_arg14(arg14), C_arg15(arg15), C_arg16(arg16), C_arg17(arg17), C_arg18(arg18), C_arg19(arg19), C_arg20(arg20), 
	C_arg21(arg21), C_arg22(arg22), C_arg23(arg23), C_arg24(arg24), C_arg25(arg25), C_arg26(arg26), C_arg27(arg27), C_arg28(arg28), C_arg29(arg29), C_arg30(arg30), 
	C_arg31(arg31), C_arg32(arg32), C_arg33(arg33), C_arg34(arg34), C_arg35(arg35), C_arg36(arg36), C_arg37(arg37), C_arg38(arg38), C_arg39(arg39), C_arg40(arg40), 
	C_arg41(arg41), itsFunc(func) {};
	long operator()( ARM_result& result, long objId ){ return (*itsFunc)( 
		C_arg1, C_arg2, C_arg3, C_arg4, C_arg5, C_arg6, C_arg7, C_arg8, C_arg9, C_arg10, 
		C_arg11, C_arg12, C_arg13, C_arg14, C_arg15, C_arg16, C_arg17, C_arg18, C_arg19, C_arg20, 
		C_arg21, C_arg22, C_arg23, C_arg24, C_arg25, C_arg26, C_arg27, C_arg28, C_arg29, C_arg30, 
		C_arg31, C_arg32, C_arg33, C_arg34, C_arg35, C_arg36, C_arg37, C_arg38, C_arg39, C_arg40,
		C_arg41, result, objId);}
private:
	const T1& C_arg1;
	const T2& C_arg2;
	const T3& C_arg3;
	const T4& C_arg4;
	const T5& C_arg5;
	const T6& C_arg6;
	const T7& C_arg7;
	const T8& C_arg8;
	const T9& C_arg9;
	const T10& C_arg10;
	const T11& C_arg11;
	const T12& C_arg12;
	const T13& C_arg13;
	const T14& C_arg14;
	const T15& C_arg15;
	const T16& C_arg16;
	const T17& C_arg17;
	const T18& C_arg18;
	const T19& C_arg19;
	const T20& C_arg20;
	const T21& C_arg21;
	const T22& C_arg22;
	const T23& C_arg23;
	const T24& C_arg24;
	const T25& C_arg25;
	const T26& C_arg26;
	const T27& C_arg27;
	const T28& C_arg28;
	const T29& C_arg29;
	const T30& C_arg30;
	const T31& C_arg31;
	const T32& C_arg32;
	const T33& C_arg33;
	const T34& C_arg34;
	const T35& C_arg35;
	const T36& C_arg36;
	const T37& C_arg37;
	const T38& C_arg38;
	const T39& C_arg39;
	const T40& C_arg40;
	const T41& C_arg41;
	exportFunc41 itsFunc;
};


///////////////////////////////////////
//// functor with 42 arguments
///////////////////////////////////////
template <
typename T1,  typename T2,  typename T3,  typename T4,  typename T5,  typename T6,  typename T7,  typename T8,  typename T9,  typename T10, 
typename T11, typename T12, typename T13, typename T14, typename T15, typename T16, typename T17, typename T18, typename T19, typename T20, 
typename T21, typename T22, typename T23, typename T24, typename T25, typename T26, typename T27, typename T28, typename T29, typename T30, 
typename T31, typename T32, typename T33, typename T34, typename T35, typename T36, typename T37, typename T38, typename T39, typename T40, 
typename T41, typename T42 > 
class exportFunc42Args : public ARMResultLong2LongFunc
{
public: 	
	typedef long (*exportFunc42)( 
		const T1&,  const T2&,	const T3&,  const T4&,  const T5&,  const T6&,  const T7&,  const T8& , const T9&,  const T10&, 
		const T11&, const T12&, const T13&, const T14&, const T15&, const T16&, const T17&, const T18&, const T19&, const T20&, 
		const T21&, const T22&, const T23&, const T24&, const T25&, const T26&, const T27&, const T28&, const T29&, const T30&, 
		const T31&, const T32&, const T33&, const T34&, const T35&, const T36&, const T37&, const T38&, const T39&, const T40&, 
		const T41&, const T42&, ARM_result&, long );
	exportFunc42Args ( 
		const T1& arg1,   const T2& arg2,   const T3& arg3,   const T4& arg4,   const T5& arg5,   const T6& arg6,   const T7& arg7,   const T8& arg8,   const T9& arg9,   const T10& arg10, 
		const T11& arg11, const T12& arg12, const T13& arg13, const T14& arg14, const T15& arg15, const T16& arg16, const T17& arg17, const T18& arg18, const T19& arg19, const T20& arg20, 
		const T21& arg21, const T22& arg22, const T23& arg23, const T24& arg24, const T25& arg25, const T26& arg26, const T27& arg27, const T28& arg28, const T29& arg29, const T30& arg30, 
		const T31& arg31, const T32& arg32, const T33& arg33, const T34& arg34, const T35& arg35, const T36& arg36, const T37& arg37, const T38& arg38, const T39& arg39, const T40& arg40, 
		const T41& arg41, const T42& arg42, exportFunc42 func)
	: C_arg1(arg1), C_arg2(arg2),	C_arg3(arg3),	C_arg4(arg4),	C_arg5(arg5),	C_arg6(arg6),	C_arg7(arg7),	C_arg8(arg8),	C_arg9(arg9),	C_arg10(arg10), 	
	C_arg11(arg11), C_arg12(arg12), C_arg13(arg13), C_arg14(arg14), C_arg15(arg15), C_arg16(arg16), C_arg17(arg17), C_arg18(arg18), C_arg19(arg19), C_arg20(arg20), 
	C_arg21(arg21), C_arg22(arg22), C_arg23(arg23), C_arg24(arg24), C_arg25(arg25), C_arg26(arg26), C_arg27(arg27), C_arg28(arg28), C_arg29(arg29), C_arg30(arg30), 
	C_arg31(arg31), C_arg32(arg32), C_arg33(arg33), C_arg34(arg34), C_arg35(arg35), C_arg36(arg36), C_arg37(arg37), C_arg38(arg38), C_arg39(arg39), C_arg40(arg40), 
	C_arg41(arg41), C_arg42(arg42), itsFunc(func) {};
	long operator()( ARM_result& result, long objId ){ return (*itsFunc)( 
		C_arg1, C_arg2, C_arg3, C_arg4, C_arg5, C_arg6, C_arg7, C_arg8, C_arg9, C_arg10, 
		C_arg11, C_arg12, C_arg13, C_arg14, C_arg15, C_arg16, C_arg17, C_arg18, C_arg19, C_arg20, 
		C_arg21, C_arg22, C_arg23, C_arg24, C_arg25, C_arg26, C_arg27, C_arg28, C_arg29, C_arg30, 
		C_arg31, C_arg32, C_arg33, C_arg34, C_arg35, C_arg36, C_arg37, C_arg38, C_arg39, C_arg40,
		C_arg41, C_arg42, result, objId);}
private:
	const T1& C_arg1;
	const T2& C_arg2;
	const T3& C_arg3;
	const T4& C_arg4;
	const T5& C_arg5;
	const T6& C_arg6;
	const T7& C_arg7;
	const T8& C_arg8;
	const T9& C_arg9;
	const T10& C_arg10;
	const T11& C_arg11;
	const T12& C_arg12;
	const T13& C_arg13;
	const T14& C_arg14;
	const T15& C_arg15;
	const T16& C_arg16;
	const T17& C_arg17;
	const T18& C_arg18;
	const T19& C_arg19;
	const T20& C_arg20;
	const T21& C_arg21;
	const T22& C_arg22;
	const T23& C_arg23;
	const T24& C_arg24;
	const T25& C_arg25;
	const T26& C_arg26;
	const T27& C_arg27;
	const T28& C_arg28;
	const T29& C_arg29;
	const T30& C_arg30;
	const T31& C_arg31;
	const T32& C_arg32;
	const T33& C_arg33;
	const T34& C_arg34;
	const T35& C_arg35;
	const T36& C_arg36;
	const T37& C_arg37;
	const T38& C_arg38;
	const T39& C_arg39;
	const T40& C_arg40;
	const T41& C_arg41;
	const T42& C_arg42;
	exportFunc42 itsFunc;
};

///////////////////////////////////////
//// functor with 43 arguments
///////////////////////////////////////
// to do...

///////////////////////////////////////
//// functor with 44 arguments
///////////////////////////////////////
// to do...

///////////////////////////////////////
//// functor with 45 arguments
///////////////////////////////////////
// to do...

///////////////////////////////////////
//// functor with 46 arguments
///////////////////////////////////////
// to do...

///////////////////////////////////////
//// functor with 47 arguments
///////////////////////////////////////
template <
typename T1,  typename T2,  typename T3,  typename T4,  typename T5,  typename T6,  typename T7,  typename T8,  typename T9,  typename T10, 
typename T11, typename T12, typename T13, typename T14, typename T15, typename T16, typename T17, typename T18, typename T19, typename T20, 
typename T21, typename T22, typename T23, typename T24, typename T25, typename T26, typename T27, typename T28, typename T29, typename T30, 
typename T31, typename T32, typename T33, typename T34, typename T35, typename T36, typename T37, typename T38, typename T39, typename T40, 
typename T41, typename T42 , typename T43 , typename T44 , typename T45 , typename T46, typename T47> 
class exportFunc47Args : public ARMResultLong2LongFunc
{
public: 	
	typedef long (*exportFunc47)( 
		const T1&,  const T2&,	const T3&,  const T4&,  const T5&,  const T6&,  const T7&,  const T8& , const T9&,  const T10&, 
		const T11&, const T12&, const T13&, const T14&, const T15&, const T16&, const T17&, const T18&, const T19&, const T20&, 
		const T21&, const T22&, const T23&, const T24&, const T25&, const T26&, const T27&, const T28&, const T29&, const T30&, 
		const T31&, const T32&, const T33&, const T34&, const T35&, const T36&, const T37&, const T38&, const T39&, const T40&, 
		const T41&, const T42&, const T43&, const T44&, const T45&, const T46&, const T47&, ARM_result&, long );
	exportFunc47Args ( 
		const T1& arg1,   const T2& arg2,   const T3& arg3,   const T4& arg4,   const T5& arg5,   const T6& arg6,   const T7& arg7,   const T8& arg8,   const T9& arg9,   const T10& arg10, 
		const T11& arg11, const T12& arg12, const T13& arg13, const T14& arg14, const T15& arg15, const T16& arg16, const T17& arg17, const T18& arg18, const T19& arg19, const T20& arg20, 
		const T21& arg21, const T22& arg22, const T23& arg23, const T24& arg24, const T25& arg25, const T26& arg26, const T27& arg27, const T28& arg28, const T29& arg29, const T30& arg30, 
		const T31& arg31, const T32& arg32, const T33& arg33, const T34& arg34, const T35& arg35, const T36& arg36, const T37& arg37, const T38& arg38, const T39& arg39, const T40& arg40, 
		const T41& arg41, const T42& arg42, const T43& arg43, const T44& arg44, const T45& arg45, const T46& arg46, const T47& arg47, exportFunc47 func)
	: C_arg1(arg1), C_arg2(arg2),	C_arg3(arg3),	C_arg4(arg4),	C_arg5(arg5),	C_arg6(arg6),	C_arg7(arg7),	C_arg8(arg8),	C_arg9(arg9),	C_arg10(arg10), 	
	C_arg11(arg11), C_arg12(arg12), C_arg13(arg13), C_arg14(arg14), C_arg15(arg15), C_arg16(arg16), C_arg17(arg17), C_arg18(arg18), C_arg19(arg19), C_arg20(arg20), 
	C_arg21(arg21), C_arg22(arg22), C_arg23(arg23), C_arg24(arg24), C_arg25(arg25), C_arg26(arg26), C_arg27(arg27), C_arg28(arg28), C_arg29(arg29), C_arg30(arg30), 
	C_arg31(arg31), C_arg32(arg32), C_arg33(arg33), C_arg34(arg34), C_arg35(arg35), C_arg36(arg36), C_arg37(arg37), C_arg38(arg38), C_arg39(arg39), C_arg40(arg40), 
	C_arg41(arg41), C_arg42(arg42), C_arg43(arg43), C_arg44(arg44), C_arg45(arg45), C_arg46(arg46), C_arg47(arg47), itsFunc(func) {};
	long operator()( ARM_result& result, long objId ){ return (*itsFunc)( 
		C_arg1, C_arg2, C_arg3, C_arg4, C_arg5, C_arg6, C_arg7, C_arg8, C_arg9, C_arg10, 
		C_arg11, C_arg12, C_arg13, C_arg14, C_arg15, C_arg16, C_arg17, C_arg18, C_arg19, C_arg20, 
		C_arg21, C_arg22, C_arg23, C_arg24, C_arg25, C_arg26, C_arg27, C_arg28, C_arg29, C_arg30, 
		C_arg31, C_arg32, C_arg33, C_arg34, C_arg35, C_arg36, C_arg37, C_arg38, C_arg39, C_arg40,
		C_arg41, C_arg42, C_arg43, C_arg44, C_arg45, C_arg46, C_arg47, result, objId);}
private:
	const T1& C_arg1;
	const T2& C_arg2;
	const T3& C_arg3;
	const T4& C_arg4;
	const T5& C_arg5;
	const T6& C_arg6;
	const T7& C_arg7;
	const T8& C_arg8;
	const T9& C_arg9;
	const T10& C_arg10;
	const T11& C_arg11;
	const T12& C_arg12;
	const T13& C_arg13;
	const T14& C_arg14;
	const T15& C_arg15;
	const T16& C_arg16;
	const T17& C_arg17;
	const T18& C_arg18;
	const T19& C_arg19;
	const T20& C_arg20;
	const T21& C_arg21;
	const T22& C_arg22;
	const T23& C_arg23;
	const T24& C_arg24;
	const T25& C_arg25;
	const T26& C_arg26;
	const T27& C_arg27;
	const T28& C_arg28;
	const T29& C_arg29;
	const T30& C_arg30;
	const T31& C_arg31;
	const T32& C_arg32;
	const T33& C_arg33;
	const T34& C_arg34;
	const T35& C_arg35;
	const T36& C_arg36;
	const T37& C_arg37;
	const T38& C_arg38;
	const T39& C_arg39;
	const T40& C_arg40;
	const T41& C_arg41;
	const T42& C_arg42;
	const T43& C_arg43;
	const T44& C_arg44;
	const T45& C_arg45;
	const T46& C_arg46;
	const T47& C_arg47;
	exportFunc47 itsFunc;
};

///////////////////////////////////////
//// functor with 50 arguments
///////////////////////////////////////
template <
typename T1,  typename T2,  typename T3,  typename T4,  typename T5,  typename T6,  typename T7,  typename T8,  typename T9,  typename T10, 
typename T11, typename T12, typename T13, typename T14, typename T15, typename T16, typename T17, typename T18, typename T19, typename T20, 
typename T21, typename T22, typename T23, typename T24, typename T25, typename T26, typename T27, typename T28, typename T29, typename T30, 
typename T31, typename T32, typename T33, typename T34, typename T35, typename T36, typename T37, typename T38, typename T39, typename T40, 
typename T41, typename T42 , typename T43 , typename T44 , typename T45 , typename T46, typename T47, 
typename T48, typename T49 , typename T50> 
class exportFunc50Args : public ARMResultLong2LongFunc
{
public: 	
	typedef long (*exportFunc50)( 
		const T1&,  const T2&,	const T3&,  const T4&,  const T5&,  const T6&,  const T7&,  const T8& , const T9&,  const T10&, 
		const T11&, const T12&, const T13&, const T14&, const T15&, const T16&, const T17&, const T18&, const T19&, const T20&, 
		const T21&, const T22&, const T23&, const T24&, const T25&, const T26&, const T27&, const T28&, const T29&, const T30&, 
		const T31&, const T32&, const T33&, const T34&, const T35&, const T36&, const T37&, const T38&, const T39&, const T40&, 
		const T41&, const T42&, const T43&, const T44&, const T45&, const T46&, const T47&,
		const T48&, const T49&, const T50&, ARM_result&, long );
	exportFunc50Args ( 
		const T1& arg1,   const T2& arg2,   const T3& arg3,   const T4& arg4,   const T5& arg5,   const T6& arg6,   const T7& arg7,   const T8& arg8,   const T9& arg9,   const T10& arg10, 
		const T11& arg11, const T12& arg12, const T13& arg13, const T14& arg14, const T15& arg15, const T16& arg16, const T17& arg17, const T18& arg18, const T19& arg19, const T20& arg20, 
		const T21& arg21, const T22& arg22, const T23& arg23, const T24& arg24, const T25& arg25, const T26& arg26, const T27& arg27, const T28& arg28, const T29& arg29, const T30& arg30, 
		const T31& arg31, const T32& arg32, const T33& arg33, const T34& arg34, const T35& arg35, const T36& arg36, const T37& arg37, const T38& arg38, const T39& arg39, const T40& arg40, 
		const T41& arg41, const T42& arg42, const T43& arg43, const T44& arg44, const T45& arg45, const T46& arg46, const T47& arg47,
		const T48& arg48, const T49& arg49, const T50& arg50, exportFunc50 func)
	: C_arg1(arg1), C_arg2(arg2),	C_arg3(arg3),	C_arg4(arg4),	C_arg5(arg5),	C_arg6(arg6),	C_arg7(arg7),	C_arg8(arg8),	C_arg9(arg9),	C_arg10(arg10), 	
	C_arg11(arg11), C_arg12(arg12), C_arg13(arg13), C_arg14(arg14), C_arg15(arg15), C_arg16(arg16), C_arg17(arg17), C_arg18(arg18), C_arg19(arg19), C_arg20(arg20), 
	C_arg21(arg21), C_arg22(arg22), C_arg23(arg23), C_arg24(arg24), C_arg25(arg25), C_arg26(arg26), C_arg27(arg27), C_arg28(arg28), C_arg29(arg29), C_arg30(arg30), 
	C_arg31(arg31), C_arg32(arg32), C_arg33(arg33), C_arg34(arg34), C_arg35(arg35), C_arg36(arg36), C_arg37(arg37), C_arg38(arg38), C_arg39(arg39), C_arg40(arg40), 
	C_arg41(arg41), C_arg42(arg42), C_arg43(arg43), C_arg44(arg44), C_arg45(arg45), C_arg46(arg46), C_arg47(arg47), 
	C_arg48(arg48), C_arg49(arg49), C_arg50(arg50), itsFunc(func) {};

	long operator()( ARM_result& result, long objId ){ return (*itsFunc)( 
		C_arg1, C_arg2, C_arg3, C_arg4, C_arg5, C_arg6, C_arg7, C_arg8, C_arg9, C_arg10, 
		C_arg11, C_arg12, C_arg13, C_arg14, C_arg15, C_arg16, C_arg17, C_arg18, C_arg19, C_arg20, 
		C_arg21, C_arg22, C_arg23, C_arg24, C_arg25, C_arg26, C_arg27, C_arg28, C_arg29, C_arg30, 
		C_arg31, C_arg32, C_arg33, C_arg34, C_arg35, C_arg36, C_arg37, C_arg38, C_arg39, C_arg40,
		C_arg41, C_arg42, C_arg43, C_arg44, C_arg45, C_arg46, C_arg47, 
		C_arg48, C_arg49, C_arg50, result, objId);}

private:
	const T1& C_arg1;
	const T2& C_arg2;
	const T3& C_arg3;
	const T4& C_arg4;
	const T5& C_arg5;
	const T6& C_arg6;
	const T7& C_arg7;
	const T8& C_arg8;
	const T9& C_arg9;
	const T10& C_arg10;
	const T11& C_arg11;
	const T12& C_arg12;
	const T13& C_arg13;
	const T14& C_arg14;
	const T15& C_arg15;
	const T16& C_arg16;
	const T17& C_arg17;
	const T18& C_arg18;
	const T19& C_arg19;
	const T20& C_arg20;
	const T21& C_arg21;
	const T22& C_arg22;
	const T23& C_arg23;
	const T24& C_arg24;
	const T25& C_arg25;
	const T26& C_arg26;
	const T27& C_arg27;
	const T28& C_arg28;
	const T29& C_arg29;
	const T30& C_arg30;
	const T31& C_arg31;
	const T32& C_arg32;
	const T33& C_arg33;
	const T34& C_arg34;
	const T35& C_arg35;
	const T36& C_arg36;
	const T37& C_arg37;
	const T38& C_arg38;
	const T39& C_arg39;
	const T40& C_arg40;
	const T41& C_arg41;
	const T42& C_arg42;
	const T43& C_arg43;
	const T44& C_arg44;
	const T45& C_arg45;
	const T46& C_arg46;
	const T47& C_arg47;
	const T48& C_arg48;
	const T49& C_arg49;
	const T50& C_arg50;
	exportFunc50 itsFunc;
};

///////////////////////////////////////
//// functor with 52 arguments
///////////////////////////////////////
template <
typename T1,  typename T2,  typename T3,  typename T4,  typename T5,  typename T6,  typename T7,  typename T8,  typename T9,  typename T10, 
typename T11, typename T12, typename T13, typename T14, typename T15, typename T16, typename T17, typename T18, typename T19, typename T20, 
typename T21, typename T22, typename T23, typename T24, typename T25, typename T26, typename T27, typename T28, typename T29, typename T30, 
typename T31, typename T32, typename T33, typename T34, typename T35, typename T36, typename T37, typename T38, typename T39, typename T40, 
typename T41, typename T42 , typename T43 , typename T44 , typename T45 , typename T46, typename T47, 
typename T48, typename T49 , typename T50, typename T51, typename T52> 
class exportFunc52Args : public ARMResultLong2LongFunc
{
public: 	
	typedef long (*exportFunc52)( 
		const T1&,  const T2&,	const T3&,  const T4&,  const T5&,  const T6&,  const T7&,  const T8& , const T9&,  const T10&, 
		const T11&, const T12&, const T13&, const T14&, const T15&, const T16&, const T17&, const T18&, const T19&, const T20&, 
		const T21&, const T22&, const T23&, const T24&, const T25&, const T26&, const T27&, const T28&, const T29&, const T30&, 
		const T31&, const T32&, const T33&, const T34&, const T35&, const T36&, const T37&, const T38&, const T39&, const T40&, 
		const T41&, const T42&, const T43&, const T44&, const T45&, const T46&, const T47&,
		const T48&, const T49&, const T50&,  const T51&, const T52&, ARM_result&, long );
	exportFunc52Args ( 
		const T1& arg1,   const T2& arg2,   const T3& arg3,   const T4& arg4,   const T5& arg5,   const T6& arg6,   const T7& arg7,   const T8& arg8,   const T9& arg9,   const T10& arg10, 
		const T11& arg11, const T12& arg12, const T13& arg13, const T14& arg14, const T15& arg15, const T16& arg16, const T17& arg17, const T18& arg18, const T19& arg19, const T20& arg20, 
		const T21& arg21, const T22& arg22, const T23& arg23, const T24& arg24, const T25& arg25, const T26& arg26, const T27& arg27, const T28& arg28, const T29& arg29, const T30& arg30, 
		const T31& arg31, const T32& arg32, const T33& arg33, const T34& arg34, const T35& arg35, const T36& arg36, const T37& arg37, const T38& arg38, const T39& arg39, const T40& arg40, 
		const T41& arg41, const T42& arg42, const T43& arg43, const T44& arg44, const T45& arg45, const T46& arg46, const T47& arg47,
		const T48& arg48, const T49& arg49, const T50& arg50, const T51& arg51, const T52& arg52, exportFunc52 func)
	: C_arg1(arg1), C_arg2(arg2),	C_arg3(arg3),	C_arg4(arg4),	C_arg5(arg5),	C_arg6(arg6),	C_arg7(arg7),	C_arg8(arg8),	C_arg9(arg9),	C_arg10(arg10), 	
	C_arg11(arg11), C_arg12(arg12), C_arg13(arg13), C_arg14(arg14), C_arg15(arg15), C_arg16(arg16), C_arg17(arg17), C_arg18(arg18), C_arg19(arg19), C_arg20(arg20), 
	C_arg21(arg21), C_arg22(arg22), C_arg23(arg23), C_arg24(arg24), C_arg25(arg25), C_arg26(arg26), C_arg27(arg27), C_arg28(arg28), C_arg29(arg29), C_arg30(arg30), 
	C_arg31(arg31), C_arg32(arg32), C_arg33(arg33), C_arg34(arg34), C_arg35(arg35), C_arg36(arg36), C_arg37(arg37), C_arg38(arg38), C_arg39(arg39), C_arg40(arg40), 
	C_arg41(arg41), C_arg42(arg42), C_arg43(arg43), C_arg44(arg44), C_arg45(arg45), C_arg46(arg46), C_arg47(arg47), 
	C_arg48(arg48), C_arg49(arg49), C_arg50(arg50), C_arg51(arg51), C_arg52(arg52), itsFunc(func) {};

	long operator()( ARM_result& result, long objId ){ return (*itsFunc)( 
		C_arg1, C_arg2, C_arg3, C_arg4, C_arg5, C_arg6, C_arg7, C_arg8, C_arg9, C_arg10, 
		C_arg11, C_arg12, C_arg13, C_arg14, C_arg15, C_arg16, C_arg17, C_arg18, C_arg19, C_arg20, 
		C_arg21, C_arg22, C_arg23, C_arg24, C_arg25, C_arg26, C_arg27, C_arg28, C_arg29, C_arg30, 
		C_arg31, C_arg32, C_arg33, C_arg34, C_arg35, C_arg36, C_arg37, C_arg38, C_arg39, C_arg40,
		C_arg41, C_arg42, C_arg43, C_arg44, C_arg45, C_arg46, C_arg47, 
		C_arg48, C_arg49, C_arg50,  C_arg51, C_arg52, result, objId);}

private:
	const T1& C_arg1;
	const T2& C_arg2;
	const T3& C_arg3;
	const T4& C_arg4;
	const T5& C_arg5;
	const T6& C_arg6;
	const T7& C_arg7;
	const T8& C_arg8;
	const T9& C_arg9;
	const T10& C_arg10;
	const T11& C_arg11;
	const T12& C_arg12;
	const T13& C_arg13;
	const T14& C_arg14;
	const T15& C_arg15;
	const T16& C_arg16;
	const T17& C_arg17;
	const T18& C_arg18;
	const T19& C_arg19;
	const T20& C_arg20;
	const T21& C_arg21;
	const T22& C_arg22;
	const T23& C_arg23;
	const T24& C_arg24;
	const T25& C_arg25;
	const T26& C_arg26;
	const T27& C_arg27;
	const T28& C_arg28;
	const T29& C_arg29;
	const T30& C_arg30;
	const T31& C_arg31;
	const T32& C_arg32;
	const T33& C_arg33;
	const T34& C_arg34;
	const T35& C_arg35;
	const T36& C_arg36;
	const T37& C_arg37;
	const T38& C_arg38;
	const T39& C_arg39;
	const T40& C_arg40;
	const T41& C_arg41;
	const T42& C_arg42;
	const T43& C_arg43;
	const T44& C_arg44;
	const T45& C_arg45;
	const T46& C_arg46;
	const T47& C_arg47;
	const T48& C_arg48;
	const T49& C_arg49;
	const T50& C_arg50;
	const T51& C_arg51;
	const T52& C_arg52;
	exportFunc52 itsFunc;
};

///////////////////////////////////////
//// functor with 53 arguments
///////////////////////////////////////
template <
typename T1,  typename T2,  typename T3,  typename T4,  typename T5,  typename T6,  typename T7,  typename T8,  typename T9,  typename T10, 
typename T11, typename T12, typename T13, typename T14, typename T15, typename T16, typename T17, typename T18, typename T19, typename T20, 
typename T21, typename T22, typename T23, typename T24, typename T25, typename T26, typename T27, typename T28, typename T29, typename T30, 
typename T31, typename T32, typename T33, typename T34, typename T35, typename T36, typename T37, typename T38, typename T39, typename T40, 
typename T41, typename T42 , typename T43 , typename T44 , typename T45 , typename T46, typename T47, 
typename T48, typename T49 , typename T50 , typename T51 , typename T52 , typename T53> 
class exportFunc53Args : public ARMResultLong2LongFunc
{
public: 	
	typedef long (*exportFunc53)( 
		const T1&,  const T2&,	const T3&,  const T4&,  const T5&,  const T6&,  const T7&,  const T8& , const T9&,  const T10&, 
		const T11&, const T12&, const T13&, const T14&, const T15&, const T16&, const T17&, const T18&, const T19&, const T20&, 
		const T21&, const T22&, const T23&, const T24&, const T25&, const T26&, const T27&, const T28&, const T29&, const T30&, 
		const T31&, const T32&, const T33&, const T34&, const T35&, const T36&, const T37&, const T38&, const T39&, const T40&, 
		const T41&, const T42&, const T43&, const T44&, const T45&, const T46&, const T47&,
		const T48&, const T49&, const T50&, const T51&, const T52&, const T53&, ARM_result&, long );
	exportFunc53Args ( 
		const T1& arg1,   const T2& arg2,   const T3& arg3,   const T4& arg4,   const T5& arg5,   const T6& arg6,   const T7& arg7,   const T8& arg8,   const T9& arg9,   const T10& arg10, 
		const T11& arg11, const T12& arg12, const T13& arg13, const T14& arg14, const T15& arg15, const T16& arg16, const T17& arg17, const T18& arg18, const T19& arg19, const T20& arg20, 
		const T21& arg21, const T22& arg22, const T23& arg23, const T24& arg24, const T25& arg25, const T26& arg26, const T27& arg27, const T28& arg28, const T29& arg29, const T30& arg30, 
		const T31& arg31, const T32& arg32, const T33& arg33, const T34& arg34, const T35& arg35, const T36& arg36, const T37& arg37, const T38& arg38, const T39& arg39, const T40& arg40, 
		const T41& arg41, const T42& arg42, const T43& arg43, const T44& arg44, const T45& arg45, const T46& arg46, const T47& arg47,
		const T48& arg48, const T49& arg49, const T50& arg50, const T51& arg51, const T52& arg52, const T53& arg53, exportFunc53 func)
	: C_arg1(arg1), C_arg2(arg2),	C_arg3(arg3),	C_arg4(arg4),	C_arg5(arg5),	C_arg6(arg6),	C_arg7(arg7),	C_arg8(arg8),	C_arg9(arg9),	C_arg10(arg10), 	
	C_arg11(arg11), C_arg12(arg12), C_arg13(arg13), C_arg14(arg14), C_arg15(arg15), C_arg16(arg16), C_arg17(arg17), C_arg18(arg18), C_arg19(arg19), C_arg20(arg20), 
	C_arg21(arg21), C_arg22(arg22), C_arg23(arg23), C_arg24(arg24), C_arg25(arg25), C_arg26(arg26), C_arg27(arg27), C_arg28(arg28), C_arg29(arg29), C_arg30(arg30), 
	C_arg31(arg31), C_arg32(arg32), C_arg33(arg33), C_arg34(arg34), C_arg35(arg35), C_arg36(arg36), C_arg37(arg37), C_arg38(arg38), C_arg39(arg39), C_arg40(arg40), 
	C_arg41(arg41), C_arg42(arg42), C_arg43(arg43), C_arg44(arg44), C_arg45(arg45), C_arg46(arg46), C_arg47(arg47), 
	C_arg48(arg48), C_arg49(arg49), C_arg50(arg50), C_arg51(arg51), C_arg52(arg52), C_arg53(arg53), itsFunc(func) {};

	long operator()( ARM_result& result, long objId ){ return (*itsFunc)( 
		C_arg1, C_arg2, C_arg3, C_arg4, C_arg5, C_arg6, C_arg7, C_arg8, C_arg9, C_arg10, 
		C_arg11, C_arg12, C_arg13, C_arg14, C_arg15, C_arg16, C_arg17, C_arg18, C_arg19, C_arg20, 
		C_arg21, C_arg22, C_arg23, C_arg24, C_arg25, C_arg26, C_arg27, C_arg28, C_arg29, C_arg30, 
		C_arg31, C_arg32, C_arg33, C_arg34, C_arg35, C_arg36, C_arg37, C_arg38, C_arg39, C_arg40,
		C_arg41, C_arg42, C_arg43, C_arg44, C_arg45, C_arg46, C_arg47, 
		C_arg48, C_arg49, C_arg50, C_arg51, C_arg52, C_arg53, result, objId);}

private:
	const T1& C_arg1;
	const T2& C_arg2;
	const T3& C_arg3;
	const T4& C_arg4;
	const T5& C_arg5;
	const T6& C_arg6;
	const T7& C_arg7;
	const T8& C_arg8;
	const T9& C_arg9;
	const T10& C_arg10;
	const T11& C_arg11;
	const T12& C_arg12;
	const T13& C_arg13;
	const T14& C_arg14;
	const T15& C_arg15;
	const T16& C_arg16;
	const T17& C_arg17;
	const T18& C_arg18;
	const T19& C_arg19;
	const T20& C_arg20;
	const T21& C_arg21;
	const T22& C_arg22;
	const T23& C_arg23;
	const T24& C_arg24;
	const T25& C_arg25;
	const T26& C_arg26;
	const T27& C_arg27;
	const T28& C_arg28;
	const T29& C_arg29;
	const T30& C_arg30;
	const T31& C_arg31;
	const T32& C_arg32;
	const T33& C_arg33;
	const T34& C_arg34;
	const T35& C_arg35;
	const T36& C_arg36;
	const T37& C_arg37;
	const T38& C_arg38;
	const T39& C_arg39;
	const T40& C_arg40;
	const T41& C_arg41;
	const T42& C_arg42;
	const T43& C_arg43;
	const T44& C_arg44;
	const T45& C_arg45;
	const T46& C_arg46;
	const T47& C_arg47;
	const T48& C_arg48;
	const T49& C_arg49;
	const T50& C_arg50;
	const T51& C_arg51;
	const T52& C_arg52;
	const T53& C_arg53;
	exportFunc53 itsFunc;
};

///////////////////////////////////////
//// functor with 54 arguments
///////////////////////////////////////
template <
typename T1,  typename T2,  typename T3,  typename T4,  typename T5,  typename T6,  typename T7,  typename T8,  typename T9,  typename T10, 
typename T11, typename T12, typename T13, typename T14, typename T15, typename T16, typename T17, typename T18, typename T19, typename T20, 
typename T21, typename T22, typename T23, typename T24, typename T25, typename T26, typename T27, typename T28, typename T29, typename T30, 
typename T31, typename T32, typename T33, typename T34, typename T35, typename T36, typename T37, typename T38, typename T39, typename T40, 
typename T41, typename T42 , typename T43 , typename T44 , typename T45 , typename T46, typename T47, 
typename T48, typename T49 , typename T50 , typename T51 , typename T52 , typename T53, typename T54> 
class exportFunc54Args : public ARMResultLong2LongFunc
{
public: 	
	typedef long (*exportFunc54)( 
		const T1&,  const T2&,	const T3&,  const T4&,  const T5&,  const T6&,  const T7&,  const T8& , const T9&,  const T10&, 
		const T11&, const T12&, const T13&, const T14&, const T15&, const T16&, const T17&, const T18&, const T19&, const T20&, 
		const T21&, const T22&, const T23&, const T24&, const T25&, const T26&, const T27&, const T28&, const T29&, const T30&, 
		const T31&, const T32&, const T33&, const T34&, const T35&, const T36&, const T37&, const T38&, const T39&, const T40&, 
		const T41&, const T42&, const T43&, const T44&, const T45&, const T46&, const T47&,
		const T48&, const T49&, const T50&, const T51&, const T52&, const T53&, const T54&, ARM_result&, long );
	exportFunc54Args ( 
		const T1& arg1,   const T2& arg2,   const T3& arg3,   const T4& arg4,   const T5& arg5,   const T6& arg6,   const T7& arg7,   const T8& arg8,   const T9& arg9,   const T10& arg10, 
		const T11& arg11, const T12& arg12, const T13& arg13, const T14& arg14, const T15& arg15, const T16& arg16, const T17& arg17, const T18& arg18, const T19& arg19, const T20& arg20, 
		const T21& arg21, const T22& arg22, const T23& arg23, const T24& arg24, const T25& arg25, const T26& arg26, const T27& arg27, const T28& arg28, const T29& arg29, const T30& arg30, 
		const T31& arg31, const T32& arg32, const T33& arg33, const T34& arg34, const T35& arg35, const T36& arg36, const T37& arg37, const T38& arg38, const T39& arg39, const T40& arg40, 
		const T41& arg41, const T42& arg42, const T43& arg43, const T44& arg44, const T45& arg45, const T46& arg46, const T47& arg47,
		const T48& arg48, const T49& arg49, const T50& arg50, const T51& arg51, const T52& arg52, const T53& arg53, const T54& arg54, exportFunc54 func)
	: C_arg1(arg1), C_arg2(arg2),	C_arg3(arg3),	C_arg4(arg4),	C_arg5(arg5),	C_arg6(arg6),	C_arg7(arg7),	C_arg8(arg8),	C_arg9(arg9),	C_arg10(arg10), 	
	C_arg11(arg11), C_arg12(arg12), C_arg13(arg13), C_arg14(arg14), C_arg15(arg15), C_arg16(arg16), C_arg17(arg17), C_arg18(arg18), C_arg19(arg19), C_arg20(arg20), 
	C_arg21(arg21), C_arg22(arg22), C_arg23(arg23), C_arg24(arg24), C_arg25(arg25), C_arg26(arg26), C_arg27(arg27), C_arg28(arg28), C_arg29(arg29), C_arg30(arg30), 
	C_arg31(arg31), C_arg32(arg32), C_arg33(arg33), C_arg34(arg34), C_arg35(arg35), C_arg36(arg36), C_arg37(arg37), C_arg38(arg38), C_arg39(arg39), C_arg40(arg40), 
	C_arg41(arg41), C_arg42(arg42), C_arg43(arg43), C_arg44(arg44), C_arg45(arg45), C_arg46(arg46), C_arg47(arg47), 
	C_arg48(arg48), C_arg49(arg49), C_arg50(arg50), C_arg51(arg51), C_arg52(arg52), C_arg53(arg53), C_arg54(arg54), itsFunc(func) {};

	long operator()( ARM_result& result, long objId ){ return (*itsFunc)( 
		C_arg1, C_arg2, C_arg3, C_arg4, C_arg5, C_arg6, C_arg7, C_arg8, C_arg9, C_arg10, 
		C_arg11, C_arg12, C_arg13, C_arg14, C_arg15, C_arg16, C_arg17, C_arg18, C_arg19, C_arg20, 
		C_arg21, C_arg22, C_arg23, C_arg24, C_arg25, C_arg26, C_arg27, C_arg28, C_arg29, C_arg30, 
		C_arg31, C_arg32, C_arg33, C_arg34, C_arg35, C_arg36, C_arg37, C_arg38, C_arg39, C_arg40,
		C_arg41, C_arg42, C_arg43, C_arg44, C_arg45, C_arg46, C_arg47, 
		C_arg48, C_arg49, C_arg50, C_arg51, C_arg52, C_arg53, C_arg54, result, objId);}

private:
	const T1& C_arg1;
	const T2& C_arg2;
	const T3& C_arg3;
	const T4& C_arg4;
	const T5& C_arg5;
	const T6& C_arg6;
	const T7& C_arg7;
	const T8& C_arg8;
	const T9& C_arg9;
	const T10& C_arg10;
	const T11& C_arg11;
	const T12& C_arg12;
	const T13& C_arg13;
	const T14& C_arg14;
	const T15& C_arg15;
	const T16& C_arg16;
	const T17& C_arg17;
	const T18& C_arg18;
	const T19& C_arg19;
	const T20& C_arg20;
	const T21& C_arg21;
	const T22& C_arg22;
	const T23& C_arg23;
	const T24& C_arg24;
	const T25& C_arg25;
	const T26& C_arg26;
	const T27& C_arg27;
	const T28& C_arg28;
	const T29& C_arg29;
	const T30& C_arg30;
	const T31& C_arg31;
	const T32& C_arg32;
	const T33& C_arg33;
	const T34& C_arg34;
	const T35& C_arg35;
	const T36& C_arg36;
	const T37& C_arg37;
	const T38& C_arg38;
	const T39& C_arg39;
	const T40& C_arg40;
	const T41& C_arg41;
	const T42& C_arg42;
	const T43& C_arg43;
	const T44& C_arg44;
	const T45& C_arg45;
	const T46& C_arg46;
	const T47& C_arg47;
	const T48& C_arg48;
	const T49& C_arg49;
	const T50& C_arg50;
	const T51& C_arg51;
	const T52& C_arg52;
	const T53& C_arg53;
	const T54& C_arg54;
	exportFunc54 itsFunc;
};

template <
typename T1,  typename T2,  typename T3,  typename T4,  typename T5,  typename T6,  typename T7,  typename T8,  typename T9,  typename T10, 
typename T11, typename T12, typename T13, typename T14, typename T15, typename T16, typename T17, typename T18, typename T19, typename T20, 
typename T21, typename T22, typename T23, typename T24, typename T25, typename T26, typename T27, typename T28, typename T29, typename T30, 
typename T31, typename T32, typename T33, typename T34, typename T35, typename T36, typename T37, typename T38, typename T39, typename T40, 
typename T41, typename T42 , typename T43 , typename T44 , typename T45 , typename T46, typename T47, 
typename T48, typename T49 , typename T50 , typename T51 , typename T52 , typename T53, typename T54, typename T55, typename T56> 
class exportFunc56Args : public ARMResultLong2LongFunc
{
public: 	
	typedef long (*exportFunc56)( 
		const T1&,  const T2&,	const T3&,  const T4&,  const T5&,  const T6&,  const T7&,  const T8& , const T9&,  const T10&, 
		const T11&, const T12&, const T13&, const T14&, const T15&, const T16&, const T17&, const T18&, const T19&, const T20&, 
		const T21&, const T22&, const T23&, const T24&, const T25&, const T26&, const T27&, const T28&, const T29&, const T30&, 
		const T31&, const T32&, const T33&, const T34&, const T35&, const T36&, const T37&, const T38&, const T39&, const T40&, 
		const T41&, const T42&, const T43&, const T44&, const T45&, const T46&, const T47&,	const T48&, const T49&, const T50&, 
		const T51&, const T52&, const T53&, const T54&, const T55&, const T56&,ARM_result&, long );
	exportFunc56Args ( 
		const T1& arg1,   const T2& arg2,   const T3& arg3,   const T4& arg4,   const T5& arg5,   const T6& arg6,   const T7& arg7,   const T8& arg8,   const T9& arg9,   const T10& arg10, 
		const T11& arg11, const T12& arg12, const T13& arg13, const T14& arg14, const T15& arg15, const T16& arg16, const T17& arg17, const T18& arg18, const T19& arg19, const T20& arg20, 
		const T21& arg21, const T22& arg22, const T23& arg23, const T24& arg24, const T25& arg25, const T26& arg26, const T27& arg27, const T28& arg28, const T29& arg29, const T30& arg30, 
		const T31& arg31, const T32& arg32, const T33& arg33, const T34& arg34, const T35& arg35, const T36& arg36, const T37& arg37, const T38& arg38, const T39& arg39, const T40& arg40, 
		const T41& arg41, const T42& arg42, const T43& arg43, const T44& arg44, const T45& arg45, const T46& arg46, const T47& arg47,
		const T48& arg48, const T49& arg49, const T50& arg50, const T51& arg51, const T52& arg52, const T53& arg53, const T54& arg54,const T55& arg55, const T56& arg56, exportFunc56 func)
	: C_arg1(arg1), C_arg2(arg2),	C_arg3(arg3),	C_arg4(arg4),	C_arg5(arg5),	C_arg6(arg6),	C_arg7(arg7),	C_arg8(arg8),	C_arg9(arg9),	C_arg10(arg10), 	
	C_arg11(arg11), C_arg12(arg12), C_arg13(arg13), C_arg14(arg14), C_arg15(arg15), C_arg16(arg16), C_arg17(arg17), C_arg18(arg18), C_arg19(arg19), C_arg20(arg20), 
	C_arg21(arg21), C_arg22(arg22), C_arg23(arg23), C_arg24(arg24), C_arg25(arg25), C_arg26(arg26), C_arg27(arg27), C_arg28(arg28), C_arg29(arg29), C_arg30(arg30), 
	C_arg31(arg31), C_arg32(arg32), C_arg33(arg33), C_arg34(arg34), C_arg35(arg35), C_arg36(arg36), C_arg37(arg37), C_arg38(arg38), C_arg39(arg39), C_arg40(arg40), 
	C_arg41(arg41), C_arg42(arg42), C_arg43(arg43), C_arg44(arg44), C_arg45(arg45), C_arg46(arg46), C_arg47(arg47), 
	C_arg48(arg48), C_arg49(arg49), C_arg50(arg50), C_arg51(arg51), C_arg52(arg52), C_arg53(arg53), C_arg54(arg54), C_arg55(arg55), C_arg56(arg56),itsFunc(func) {};

	long operator()( ARM_result& result, long objId ){ return (*itsFunc)( 
		C_arg1, C_arg2, C_arg3, C_arg4, C_arg5, C_arg6, C_arg7, C_arg8, C_arg9, C_arg10, 
		C_arg11, C_arg12, C_arg13, C_arg14, C_arg15, C_arg16, C_arg17, C_arg18, C_arg19, C_arg20, 
		C_arg21, C_arg22, C_arg23, C_arg24, C_arg25, C_arg26, C_arg27, C_arg28, C_arg29, C_arg30, 
		C_arg31, C_arg32, C_arg33, C_arg34, C_arg35, C_arg36, C_arg37, C_arg38, C_arg39, C_arg40,
		C_arg41, C_arg42, C_arg43, C_arg44, C_arg45, C_arg46, C_arg47, 
		C_arg48, C_arg49, C_arg50, C_arg51, C_arg52, C_arg53, C_arg54, C_arg55, C_arg56,result, objId);}

private:
	const T1& C_arg1;
	const T2& C_arg2;
	const T3& C_arg3;
	const T4& C_arg4;
	const T5& C_arg5;
	const T6& C_arg6;
	const T7& C_arg7;
	const T8& C_arg8;
	const T9& C_arg9;
	const T10& C_arg10;
	const T11& C_arg11;
	const T12& C_arg12;
	const T13& C_arg13;
	const T14& C_arg14;
	const T15& C_arg15;
	const T16& C_arg16;
	const T17& C_arg17;
	const T18& C_arg18;
	const T19& C_arg19;
	const T20& C_arg20;
	const T21& C_arg21;
	const T22& C_arg22;
	const T23& C_arg23;
	const T24& C_arg24;
	const T25& C_arg25;
	const T26& C_arg26;
	const T27& C_arg27;
	const T28& C_arg28;
	const T29& C_arg29;
	const T30& C_arg30;
	const T31& C_arg31;
	const T32& C_arg32;
	const T33& C_arg33;
	const T34& C_arg34;
	const T35& C_arg35;
	const T36& C_arg36;
	const T37& C_arg37;
	const T38& C_arg38;
	const T39& C_arg39;
	const T40& C_arg40;
	const T41& C_arg41;
	const T42& C_arg42;
	const T43& C_arg43;
	const T44& C_arg44;
	const T45& C_arg45;
	const T46& C_arg46;
	const T47& C_arg47;
	const T48& C_arg48;
	const T49& C_arg49;
	const T50& C_arg50;
	const T51& C_arg51;
	const T52& C_arg52;
	const T53& C_arg53;
	const T54& C_arg54;
	const T55& C_arg55;
	const T56& C_arg56;
	exportFunc56 itsFunc;
};

#endif /* ARM_XL_GP_FCTOR_HELPER */
