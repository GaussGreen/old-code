/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * $Log: eventviewer.h,v $
 * Revision 1.1  2004/03/01 07:52:17  ebenhamou
 * Initial revision
 *
 */


/*! \file eventviewer.h
 *
 *  \brief
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date March 2004
 */


#ifndef _INGPBASE_EVENTVIEWER_H
#define _INGPBASE_EVENTVIEWER_H

#include "port.h"

#include <string>
CC_USING_NS(std,string)

#include <cstdio>

CC_BEGIN_NAMESPACE( ARM )

/// forward declaration
template <typename T> class ARM_SingletonHolder;


/////////////////////////////////////////////////////////////////
/// \struct ARM_EventViewerImp
/// implementor class for the singleton ARM_EventViewer
/// design choice is that it should never be of size more than SIZEMAX
/// debug message gives the ability to use tracer in release mode...
/// this should be used only for debugging and removed when giving 
/// a release version to the desk.
/////////////////////////////////////////////////////////////////

struct ARM_EventViewerImp
{
	void ResetMessage() { itsMessage = ""; }
	void AddToMessage( const string& msg )
	{ 
		if( ARM_EventViewerImp::VerboseIsOn )
			AddToMessageCommon( msg, itsMessage );
	}

	void ResetDebugMessage() { itsMessage = ""; }
	void AddToDebugMessage( const string& msg )
	{ 
		if( ARM_EventViewerImp::DebugIsOn)
			AddToMessageCommon( msg, itsDebugMessage );
	}

	inline string GetMessage() const { return itsMessage; }
	inline string GetDebugMessage() const { return itsDebugMessage; }
	string toString() const;

	void WriteMessageToFile( const string& fileName, const string& txt, const char* iosMode = "a", bool addDebugReleaseToTile = true, const string& extension = "txt"  )
	{
		FILE* file;
		string finalFileName( fileName );
#if defined( _DEBUG )
		if( addDebugReleaseToTile )
			finalFileName += "_Debug."	+ extension;
#endif

#if defined( NDEBUG )
		if( addDebugReleaseToTile )
			finalFileName += "_Release."+ extension;
#endif

		file = fopen( finalFileName.c_str(), iosMode );
		fprintf( file, "%s" , txt.c_str() );
		fclose( file );
	}

	/// static variables to control the addition of text
	static bool SetVerboseMode( bool value ) { return ARM_EventViewerImp::VerboseIsOn=value; }
	static bool VerboseIsOn;
	static bool SetDebugMode(	bool value ) { return ARM_EventViewerImp::DebugIsOn=value; }
	static bool DebugIsOn;

private:
	void AddToMessageCommon( const string& msg, string& txtToAppendTo )
	{
		const size_t SIZEMAX = 10000000;
		txtToAppendTo += msg; 
		if( txtToAppendTo.size() > SIZEMAX )
			txtToAppendTo = txtToAppendTo.substr(txtToAppendTo.size()-SIZEMAX-1);
	}

	string itsMessage;
	string itsDebugMessage;


	/// to forbid client from using it except for the singleton holder
	ARM_EventViewerImp(): itsMessage(), itsDebugMessage() {};
	~ARM_EventViewerImp() {};

	/// and non implemented copy and assignment operator to avoid duplication
	ARM_EventViewerImp( const ARM_EventViewerImp& rhs );
	ARM_EventViewerImp& operator=( const ARM_EventViewerImp& rhs );

	friend class ARM_SingletonHolder<ARM_EventViewerImp>;
};

extern ARM_SingletonHolder<ARM_EventViewerImp> ARM_TheEventViewer;

CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/