/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file errviewer.cpp
 *
 *  \brief
 *
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date March 2004
 */

#include "gpbase/errviewer.h"
#include "gpbase/ostringstream.h"

/// STL
#include <fstream>

CC_BEGIN_NAMESPACE( ARM )

/////////////////////////////////////////////////////////////////
///	Class  : ARM_ErrViewer
///	Data   : filename where the truncated error are stored
/////////////////////////////////////////////////////////////////

string GetGPTmpFile(char* fileName)
{
    char fOutName[200];

    ARM_GetTmpAbsFile(fileName, fOutName);

    string absName(fOutName);

    return(absName);
}

const string ARM_ErrViewer::LogOfAllErrorsFileNameStart	= GetGPTmpFile("LogAllErrors_");
                                                                       
const string ARM_ErrViewer::LogOfAllErrorsFileNameExt	= ".txt";

/////////////////////////////////////////////////////////////////
///	Class  : ARM_ErrViewer
///	Routine: constructor,copy constructor,assignment operator and destructor
///	Returns: 
///	Action : standard action for derived object
/////////////////////////////////////////////////////////////////

ARM_ErrViewer::ARM_ErrViewer( bool reset )
:	ARM_RootObject() 
{
	if( reset)
	{
		string fileName = ARM_ErrViewer::LogOfAllErrorsFileNameStart + ARM_USERNAME + ARM_ErrViewer::LogOfAllErrorsFileNameExt;
		FILE* file = fopen( fileName.c_str(), "w+" );
		fclose( file );
	}
	SetName(ARM_ERRVIEWER);
}

ARM_ErrViewer::ARM_ErrViewer( const ARM_ErrViewer& rhs ) 
:	ARM_RootObject( rhs )
{}


ARM_ErrViewer& ARM_ErrViewer::operator=( const ARM_ErrViewer& rhs )
{
	if( this != &rhs )
		ARM_RootObject::operator=(rhs);
	return *this;
}

ARM_ErrViewer::~ARM_ErrViewer()
{}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_ErrViewer
///	Routine: Clone
///	Returns: ARM_Object*
///	Action : standard action for derived object
/////////////////////////////////////////////////////////////////
ARM_Object* ARM_ErrViewer::Clone() const
{
	return new ARM_ErrViewer(*this);
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_ErrViewer
///	Routine: toString
///	Returns: string
///	Action : stringify the object
/////////////////////////////////////////////////////////////////

string ARM_ErrViewer::toString( const string& indent, const string& nextIndent ) const
{
	string fileName = ARM_ErrViewer::LogOfAllErrorsFileNameStart + ARM_USERNAME + ARM_ErrViewer::LogOfAllErrorsFileNameExt;
	CC_NS(std,ifstream) myFile;
	myFile.open(fileName.c_str());
	CC_Ostringstream os;

	const int lineSize = 1000;
	char buffer[lineSize];
	size_t i=0;
	char firstChar;
		
	while( myFile.get(firstChar) )
	{
		myFile.getline(buffer,lineSize);
		if( buffer[0] )
			os << " line: "<< i++ << " " << firstChar << buffer << "\n\n\n";
	}
	

	myFile.close();

	return os.str();
}




CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

