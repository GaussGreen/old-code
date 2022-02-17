/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * $Log: expt.cpp,v $
 * Revision 1.11  2004/06/25 14:55:00  ranger
 * added copy constructor and handle leak in destructor
 *
 * Revision 1.10  2004/06/01 17:39:24  mab
 * Added : condition:  if (gTrace)
 * in : void Exception::DebugPrint(void)
 *
 * Revision 1.9  2004/03/15 11:14:22  ebenhamou
 * added better truncation error
 *
 * Revision 1.7  2004/03/12 09:30:45  ebenhamou
 * remove truncation as this is done at a higher level
 *
 * Revision 1.6  2004/01/07 07:39:38  jmprie
 * troncature de taille si celle du message (hors FILE & LINE) >= 170
 * et ajout des "..." seulement ds ce cas
 *
 * Revision 1.5  2003/11/18 12:08:07  jmprie
 * troncature du message a 150 char
 *
 * Revision 1.4  2003/10/14 15:06:45  ebenhamou
 * handle leak for string constructor
 *
 * Revision 1.3  2003/09/26 07:43:49  ebenhamou
 * added log
 *
 *
 *
 */



/*----------------------------------------------------------------------------*/
/*                                                                            */
/* FILE        : expt.cpp                                                     */
/*                                                                            */
/* DESCRIPTION : Exceptions classes                                           */
/*                                                                            */
/* DATE        : Tue Jun 11 1996                                              */
/*                                                                            */
/*----------------------------------------------------------------------------*/

#include "expt.h"


#include <cstdio> /// for full message
#include <stdlib.h>
#include <windows.h>	// for ::OutputDebugString



/// string support
#ifndef std
    using std::string;
#endif


extern char* ARM_STRDUP(char* inStr);




/*----------------------------------------------------------------------------*/



void InitLogFileName(char* fileName)
{
    static char* LogFile = NULL;

    if ( LogFile == NULL )
    {
       char* userName = getenv("USERNAME");

       char tmpName[100];

       char fOutName[200];

       ARM_GetTmpAbsFile("LogAllErrors_", fOutName);


       sprintf(tmpName, "%s%s%s", fOutName,
                                  userName,
                                  ".txt");

       strcpy(fileName, tmpName);

       LogFile = new char [strlen(fileName)+1];

       strcpy(LogFile, fileName);
    }
    else
    {
       strcpy(fileName, LogFile);
    }
}



void Exception::Init(void)
{
    errCode = 0;
    errLine = 0;

    errFunc = (char *) NULL;
    errFile = (char *) NULL;
    message = (char *) NULL;
    
    char theLogFile[100];

    InitLogFileName(theLogFile);
}



// Exception::Exception(void)
// {
//     Init();
// }



/// function to shorten the error Message
/// the full error message is not lost and is
/// store in the file "C:\\temp\\LogAllErrors_<UserName>.txt"




void Exception::SetGoodSizeMessage(const string& initialMessage ) 
{
    const char* initialMessageChar = initialMessage.c_str();

    if (message)
       delete message;

    message = new char[initialMessage.size()+1];
    
    strcpy(message, initialMessageChar);
}



//// copy constructor
Exception::Exception(const Exception& e)
{
    Init();

	errLine = e.errLine;

    if (e.errFunc)
       errFunc = ARM_STRDUP(e.errFunc);

    if (e.errFile)
       errFile = ARM_STRDUP(e.errFile);

    errCode = e.errCode;

    if (e.message)
       message = ARM_STRDUP(e.message);

    itsLongMessageLog = e.itsLongMessageLog;
}



//// standard constructor
Exception::Exception(long line, char* file, RET_CODE code, 
                     char* mesg, char* func)
{
    Init();

    errLine = line;

    if (func)
       errFunc = ARM_STRDUP(func); 
    else
       errFunc = ARM_STRDUP("NO_FUNC");

    if (file)
       errFile = ARM_STRDUP(file);
    else
       errFile = ARM_STRDUP("NO_FILE");

    errCode = code;

    if (mesg)
    {
       message = ARM_STRDUP(mesg);
    }
    else
    {
       message = ARM_STRDUP("NO_ERRMSG");
    }
}



/// string version of the constructor
Exception::Exception(long line, char* file, RET_CODE code, const string& mesg,
                     const string& func )
{
    Init();

    errLine = line;

    if ( func == "" )
       errFunc = ARM_STRDUP("NO_FUNC");
    else
    {
       errFunc = new char[func.size()+1];

       strcpy( errFunc, func.c_str() );
    }
 
    if (file)
       errFile = ARM_STRDUP(file);
    else
       errFile = ARM_STRDUP("NO_FILE");

    errCode = code;

    if ( mesg == "" )
       message = ARM_STRDUP("NO_ERRMSG");
    else
       SetGoodSizeMessage(mesg);
}



void Exception::DebugPrint(void)
{
    if (gTrace)
    {
       printf("\n ??===> Error : %s", message);
       printf("\n      File  : %s", errFile);

       if (errFunc)
          printf("\n      Func  : %s", errFunc); 

       printf("\n      Line  : %ld", errLine);
       printf("\n      NoErr : %ld\n", errCode);
    }
}



Exception::~Exception()
{
    if (errFunc)
       delete errFunc;
    errFunc = NULL;

    if (errFile)
       delete errFile;
    errFile = NULL;

    if (message)
       delete message;
    message = NULL;
}


//	--------------------------------------------------------------------------------
void ArmLogger::logDebugger(const std::string& file,long line,const std::string& msg)
{
	std::stringstream sstr; sstr<<file<<" ("<<line<<") :"<<msg ; 
	::OutputDebugString(sstr.str().c_str()); 
}

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/
