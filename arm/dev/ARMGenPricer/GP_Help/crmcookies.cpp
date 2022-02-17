/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file crmcookies.cpp
 *  \brief
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date February 2004
 */

/// this header comes firts as it includes some preprocessor constants!
#include "gpbase/removeidentifiedwarning.h"

/// gphelp
#include "gphelp/crmcookies.h"
#include "gphelp/splashwindow.h"

/// kernel
#include <glob/dates.h>

#include <fstream>
CC_USING_NS(std,fstream)
#include <iostream>
CC_USING_NS(std,ios)
#include <iomanip>
CC_USING_NS(std,setw)
CC_USING_NS(std,left)

/// uncomment this to get the splash window
///#define SHOW_GP_SPLASH_WINDOW


CC_BEGIN_NAMESPACE( ARM )


///////////////////////////////////////////////////
///	Class  : ARM_CRMCookiesImp
///	Routine: constructor
///	Returns: 
///	Action : builds the object: says that it has 
///				not been initialised
////////////////////////////////////////////////////
ARM_CRMCookiesImp::ARM_CRMCookiesImp()
:	itsStartTime(-1), itsFileName(""), itsServices()
{}


///////////////////////////////////////////////////
///	Class  : ARM_CRMCookiesImp
///	Routine: RegisterService
///	Returns: 
///	Action : register service and set up start time if not initialised
////////////////////////////////////////////////////
void ARM_CRMCookiesImp::RegisterService( const string& userName, const string& serviceName, const string& folderName )
{
	if(itsStartTime==-1)
		InitTracing(userName, folderName);
	
	/// is it a new service?
	if( itsServices.insert(serviceName ).second )
	{
		fstream file;
		/// test if the filename is not null
		if(itsFileName != "" )
		{
			file.open( itsFileName.c_str(), CC_NS(ios,out)| CC_NS(ios,app));
			file << "\t" << serviceName;
		}
	}
}

///////////////////////////////////////////////////
///	Class  : ARM_CRMCookiesImp
///	Routine: RegisterService
///	Returns: 
///	Action : close the file saving the services and the total time spent in minutes
///				save the last line nb to get fast lookup for the next connection
////////////////////////////////////////////////////
void ARM_CRMCookiesImp::EndTracing()
{
	/// static variable to guarantee that this is done only once!
	static bool done = false;

	if(itsStartTime!=-1 && !done )
	{
		fstream file;
		/// test if the filename is not null
		if(itsFileName != "" )
		{
			file.open( itsFileName.c_str(), CC_NS(ios,out)| CC_NS(ios,app));
			clock_t endTime = clock();
			file << "\t" << (endTime-itsStartTime)/CLOCKS_PER_SEC/60.0 << CC_NS(std,endl);
			file.close();
			done = true;
		}
	}
}


///////////////////////////////////////////////////
///	Class  : ARM_CRMCookiesImp
///	Routine: InitTracing
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
void ARM_CRMCookiesImp::InitTracing( const string& userName, const string& folderName )
{
	itsStartTime = clock();
	char* pPath;
	itsFileName = "\\\\cm.net\\ShareParis\\Salle\\Services\\Echanges\\CRMLog";
	pPath = getenv(itsFileName.c_str());
	if(!pPath)
		return;

	if( folderName != "" )
		itsFileName += string("\\") + folderName;
	itsFileName += string("\\") + userName + "_cookies.txt";

	fstream file;
	
	file.open( itsFileName.c_str(), CC_NS(ios,in));
	
	bool shouldInsertNewLine = false;
	
	/// does this file already exist?
	/// apparently no!
	if(file.fail())
	{
		file.clear();
		file.open(itsFileName.c_str(),CC_NS(ios,out) );

		/// are we on a computer with no access to the server...
		/// in this case, we just exist and do nothing!
		if( file.fail() )
		{
			/// no file name, hence no tracing! (typical case can be EK!)
			itsFileName = "";
			return ;
		}
		else
		{
			file  <<  left << setw(10) << "Date" << "\t" 
			  << left << setw(10)  << "Machine" << "\t" << "Xll Date & Time Stamp + Version " << "\t" << "Services" << "\t" << "TotalTime" << CC_NS(std,endl); 
		}
	}
	else
	{
		/// has to test wheter to insert some new lines!		
		file.seekg(-1,std::ios::end);
		char lastChar = file.get();

		if(lastChar != '\n' )
			shouldInsertNewLine = true;
		
		file.close();
		file.open(itsFileName.c_str(),CC_NS(ios,out)| CC_NS(ios,app) );
	
		if(shouldInsertNewLine )
			file << "\tCRASHED LAST TIME" << CC_NS(std,endl);
	}
		
	/// close everything
	file.close();

	/// reopen in write mode
	file.open(itsFileName.c_str(),CC_NS(ios,out)| CC_NS(ios,app) );
	ARM_Date now;
	file << now.toString() << "\t" <<  ARM_COMPUTERNAME << "\t" << __TIMESTAMP__;

	#if defined(_DEBUG)
		file << " Debug\t";
	#else
		file << " Release\t";
	#endif

	file.close();

	/// splash window activation
#ifdef SHOW_GP_SPLASH_WINDOW
	/// Display the splash window
	/// get the ressource
	ARM_SplashWindow window(ARM_SplashWindow::itsResourceHandle,ARM_SplashWindow::itsModuleHandle);
	window.ShowSplashWindow();
	
	/// wait for 1.5 second and close it
	Sleep(1500);
	window.CloseSplashWindow();
#endif
}


/// creation of the object ARM_CRMCookies 
ARM_SingletonHolder<ARM_CRMCookiesImp> ARM_CRMCookies;


CC_END_NAMESPACE()

//-----------------------------------------------------------------------------
/*---- End of file ----*/
