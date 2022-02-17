/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * $Log: splashwindow.cpp,v $
 * Revision 1.1  2004/06/08 16:44:43  ebenhamou
 * Initial revision
 *
 */

/*! \file splashwindow.cpp
 *
 *  \brief
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date June 2004
 */

#include "gphelp/splashwindow.h"


CC_BEGIN_NAMESPACE( ARM )

/// static const for the class name!
const string ARM_SplashWindow::ClassNameString	= "Splash Window Class";
const string ARM_SplashWindow::WindowNameString = "Splash Window";
int ARM_SplashWindow::itsResourceHandle			= 0;
HMODULE ARM_SplashWindow::itsModuleHandle		= HMODULE();

////////////////////////////////////////////////////
///	Class  : ARM_SplashWindow
///	Routine: Constructor
///	Returns: 
///	Action : Constructor transforming the fileName into a picture
////////////////////////////////////////////////////
ARM_SplashWindow::ARM_SplashWindow( WORD ResourceHandle, HMODULE ModuleHandle )
:	itsPicture(NULL), itsWindow(NULL)
{
	GetPicture(ResourceHandle,ModuleHandle);
	RegisterClass();
}


////////////////////////////////////////////////////
///	Class  : ARM_SplashWindow
///	Routine: Constructor
///	Returns: 
///	Action : Registering the class for the window
///    This registration and its usage is only necessary if you want this code
///    to be compatible with Win32 systems prior to the 'RegisterClassEx'
///    function that was added to Windows 95. It is important to call RegisterClassEx
///    so that the application will get 'well formed' small icons associated
///    with it.
////////////////////////////////////////////////////
void ARM_SplashWindow::RegisterClass()
{
	HINSTANCE hInstance = GetModuleHandle(NULL);

	WNDCLASSEX wcex;
	wcex.cbSize = sizeof(WNDCLASSEX); 
	wcex.lpfnWndProc	= (WNDPROC) DefWindowProc;
	wcex.hInstance		= GetModuleHandle(NULL);
	wcex.lpszClassName	= ARM_SplashWindow::ClassNameString.c_str();

	wcex.style			= 0; 
	wcex.cbClsExtra		= 0;
	wcex.cbWndExtra		= 0;

	wcex.hIcon			= NULL; 
	wcex.hCursor		= NULL; 
	wcex.hbrBackground	= NULL;
	wcex.lpszMenuName	= NULL;
	wcex.hIconSm		= NULL; 
	
	RegisterClassEx(&wcex);
}


////////////////////////////////////////////////////
///	Class  : ARM_SplashWindow
///	Routine: Release
///	Returns: 
///	Action : Releases all the ressources
////////////////////////////////////////////////////
void ARM_SplashWindow::Release()
{
	/// release the picture 
	if( itsPicture)
		itsPicture->Release();
	
	/// remove the registered class
	UnregisterClass(  ARM_SplashWindow::ClassNameString.c_str(), GetModuleHandle(NULL) );

	/// free the libraries associated to the model
	FreeLibrary( GetModuleHandle(NULL) );
}


////////////////////////////////////////////////////
///	Class  : ARM_SplashWindow
///	Routine: ShowSplashWindow
///	Returns: 
///	Action : Show the window.. use a different name
///				to avoid confusing with ShowWindow from SDK
////////////////////////////////////////////////////
void ARM_SplashWindow::ShowSplashWindow()
{
	if (itsPicture)
	{
		const int HIMETRIC_PIXEL	= 26;

		/// get width and height of picture
		long hmWidth;
		long hmHeight;
		itsPicture->get_Width(&hmWidth);
		itsPicture->get_Height(&hmHeight);

		/// convert himetric to pixels
		int nHeight	= hmHeight / HIMETRIC_PIXEL;
		int nWidth	= hmWidth  / HIMETRIC_PIXEL;

		/// convert window to size
		HWND hDesktop = GetDesktopWindow();
		HDC hdcScreen = GetWindowDC(hDesktop);
		int DesktopSizeHeight = GetDeviceCaps(hdcScreen,HORZRES);
		int DesktopSizeWidth  = GetDeviceCaps(hdcScreen,VERTRES);


		int xOffset = (DesktopSizeHeight - nWidth ) / 2;
		int yOffset	= (DesktopSizeWidth  - nHeight) / 2;

		/// create the appropriate window
		itsWindow = CreateWindow( ARM_SplashWindow::ClassNameString.c_str(), 
			ARM_SplashWindow::WindowNameString.c_str(),
			WS_POPUP, xOffset,  yOffset, nWidth, nHeight, NULL, NULL, GetModuleHandle(NULL), NULL);

		/// we have to show first the window to make it active when painting
		ShowWindow(itsWindow, SW_SHOWNORMAL );

		PAINTSTRUCT ps;
		HDC hdc = BeginPaint(itsWindow, &ps);

		RECT rc;
		GetClientRect(itsWindow, &rc);
			
		/// display picture using IPicture::Render
		itsPicture->Render(hdc, 0, 0, nWidth, nHeight, 0, hmHeight, hmWidth, -hmHeight, &rc);

		EndPaint(itsWindow, &ps);
	}
}

////////////////////////////////////////////////////
///	Class  : ARM_SplashWindow
///	Routine: ShowSplashWindow
///	Returns: 
///	Action : Show the window.. use a different name
///				to avoid confusing with ShowWindow from SDK
////////////////////////////////////////////////////
void ARM_SplashWindow::CloseSplashWindow()
{
	DestroyWindow(itsWindow);
}


////////////////////////////////////////////////////
///	Class  : ARM_SplashWindow
///	Routine: GetPicture
///	Returns: 
///	Action : Load the picture knowing its fileName
///				and convert it to a bitmap object
////////////////////////////////////////////////////
void ARM_SplashWindow::GetPicture( WORD ResourceHandle, HMODULE ModuleHandle )
{
	// Open the resource
	HRSRC res = FindResource(ModuleHandle, MAKEINTRESOURCE(ResourceHandle), "BINARY");

	if (res) 
	{
		HGLOBAL mem = LoadResource(ModuleHandle, res);
		LPVOID data = LockResource(mem);
		size_t dwSize = SizeofResource(ModuleHandle, res);

		HGLOBAL hGlobal = GlobalAlloc(GMEM_MOVEABLE, dwSize);
		LPVOID pvData = GlobalLock( hGlobal );
		memcpy(pvData,data,dwSize);
		GlobalUnlock(hGlobal);
		LPSTREAM pStream = NULL;
		HRESULT hr = CreateStreamOnHGlobal( hGlobal, TRUE, &pStream );

		_ASSERTE(SUCCEEDED(hr) && pStream);

		if( itsPicture)
			itsPicture->Release();

		/// this is very strange .purify complains that there is a memory leak but when using a very
		/// large picture, we can check on the task bar that there is no memory leak,,,, NO MEMORY LEAK
		/// very very strange. It seems that purity was wrong this time... but hey, this is Windows stuff
		/// so not very transparent either
		hr = ::OleLoadPicture( pStream, dwSize, FALSE, IID_IPicture, (LPVOID *)&itsPicture);
		_ASSERTE(SUCCEEDED(hr) && itsPicture);
		if(pStream)
			pStream->Release();
	}
}

CC_END_NAMESPACE()


/*---------------------------------------------------------------------------*/
/*---- End of file ----*/