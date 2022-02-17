/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * $Log: splashwindow.h,v $
 * Revision 1.1  2004/06/08 07:52:17  ebenhamou
 * Initial revision
 *
 */

/*! \file splashwindow.h
 *
 *  \brief
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date June 2004
 */


#ifndef _INGPHELP_SPLASHWINDOW_H
#define _INGPHELP_SPLASHWINDOW_H

#include "gpbase/port.h"	/// port.h
#include <olectl.h> /// for ole
#include <crtdbg.h>	/// for assert macro

#include <string>
CC_USING_NS(std,string)

CC_BEGIN_NAMESPACE( ARM )


/// We choose to use directly the Win32 SDK instead of MFC
/// to avoid paying for the burden of MFC...
/// as this class is very simple!

class ARM_SplashWindow
{
public:
	ARM_SplashWindow( WORD ResourceHandle, HMODULE DLLHandle );
	void ShowSplashWindow();
	void CloseSplashWindow();
	void Release();

	/// static for initialisation
	static void SetResourceHandle(int ressourceHandle ) { ARM_SplashWindow::itsResourceHandle=ressourceHandle; }
	static void SetModuleHandle(HANDLE g_hInst) { ARM_SplashWindow::itsModuleHandle=(HMODULE)g_hInst; }
	static int itsResourceHandle;
	static HMODULE itsModuleHandle;

private:
	void GetPicture( WORD ResourceHandle, HMODULE DLLHandle );
	void RegisterClass();
	void UnRegisterClass();

	static const string ClassNameString;
	static const string WindowNameString;

	HWND itsWindow;
	LPPICTURE itsPicture;   /// pointor to a picture file
};


CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/

