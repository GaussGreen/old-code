//#include "stdafx.h"

#ifdef WIN32 
	#pragma warning(disable:4786)
	#include <iostream>
	using namespace std;	
#else
	#include <iostream.h>
#endif

#include "CInitReader.h"

const char* CInitReaderExceptionMsg[CInitReaderException::MAX_EXCEPTION_NO] = {
  "File not found",
  "Cannot find section",
  "Cannot find the property",
  "Invalid Conversion",
  "Invalid syntax for a property. Should be PROPERTY=VALUE"
};



/*************************************************************************************************/
CInitReaderException::CInitReaderException(ExceptionNo aNo, const char* pTextToAppend) {
  mErrorMsg = CInitReaderExceptionMsg[aNo];
  if (pTextToAppend) {
    mErrorMsg += " ";
    mErrorMsg += pTextToAppend;
  }
}

/*************************************************************************************************/
CInitReaderException::CInitReaderException(const char* pMsg) {
  mErrorMsg = pMsg;
}


/*************************************************************************************************/
const char*  CInitReaderException::GetMessage() const {
  return (const char*)mErrorMsg.c_str();
}

/*************************************************************************************************/
 CInitReaderException:: operator const char* () const {
   return (const char*) mErrorMsg.c_str();
}


/*************************************************************************************************/
CInitReader::CInitReader(const char* pFileName, int aContinue) {
  mFileName = pFileName;
  mContinue = aContinue;
  mCurSectionPos = 0;
  mCurSectionEnd = mBuffer.length();
  mSubSection = mBuffer.substr(mCurSectionPos, mCurSectionEnd-mCurSectionPos);
}

/*************************************************************************************************/
void CInitReader::SetSection(const char* pSectionName) {
  int index;
  string tmp = "[";
  tmp += pSectionName;
  tmp += "]";
  index = mBuffer.find(pSectionName);
  
  if(index <= 0) {
    CInitReaderException* pExcep = new CInitReaderException(CInitReaderException::SECTION_NOT_FOUND, pSectionName);
    if (mContinue) {
		cout << pExcep->GetMessage();
      delete pExcep;
    }
    else
      throw pExcep;
  }
  mCurSectionPos = index + tmp.length();
  mCurSectionEnd = mBuffer.find("[", mCurSectionPos);
  if (mCurSectionEnd == -1) {
    mCurSectionEnd = mBuffer.length();
  }
  mpCurSection = ((const char*) mBuffer.c_str() ) + index;

  mSubSection = mBuffer.substr(mCurSectionPos, mCurSectionEnd-mCurSectionPos);
  
}

/*************************************************************************************************/
string CInitReader::GetString(const char* pProperty) {
  int index;
  int debut = 0;
  int lengt = 0;
  int error = 0;
  int fin = 0;
  int egualFound = 0;
  string value = "";
 
  index = mSubSection.find(pProperty);
  
  if (index < 0) {
    CInitReaderException* pExcep = new CInitReaderException(CInitReaderException::PROPERTY_NOT_FOUND, pProperty);
    if (mContinue) {
      cout << pExcep->GetMessage();
      delete pExcep;
    }
    else
      throw pExcep;
  }
  index = index + strlen(pProperty);
 
 // find the attribute value
  while(!debut && !error) {
    char ch = mSubSection[index];
    
    if (ch == '=') {
      egualFound = 1;
    }
   
    if ((ch != ' ') && (ch != '=')) {
      if (!egualFound) { 
	CInitReaderException* pExcep = new CInitReaderException(CInitReaderException::INVALID_SYNTAX, pProperty);
	if (mContinue) {
	  cout << pExcep->GetMessage();
	  delete pExcep;
	}
	else
	  throw pExcep;
      }
      else debut = index;
    }
   
    if (index == mSubSection.length() ) {
      CInitReaderException* pExcep =   new CInitReaderException(CInitReaderException::INVALID_SYNTAX, pProperty);
      if (mContinue) {
	cout << pExcep->GetMessage();
	delete pExcep;
      }
      else
	throw pExcep;
    }
    index ++;
  }

  

  while(!fin) {
    char ch = mSubSection[index];
    if ( (ch == ' ') || (ch == '\n') || (ch == 0) )
      fin = index ;

    index++;
  }
    
  value = mSubSection.substr(debut, fin-debut);
  return value;
 
}

/*************************************************************************************************/
float CInitReader::  GetFloat(const char* pProperty) {
  return atof((const char*) GetString(pProperty).c_str() );
}

/*************************************************************************************************/
int CInitReader::GetInt(const char* pProperty) {
  return atoi((const char*)GetString (pProperty).c_str());
}


/*************************************************************************************************/
void  CInitReader::Read() {
  mpFile = fopen((const char*)mFileName.c_str(), "r");
  char buffer[200];
  mBuffer = "";

  if (!mpFile) {
    CInitReaderException* pExcep= new CInitReaderException(CInitReaderException::FILE_NOT_FOUND, mFileName.c_str());
    if (mContinue) {
      cout << pExcep->GetMessage();
      delete pExcep;
    }
    else
      throw pExcep;
  }
  while(fgets(buffer,200, mpFile) != NULL) {
    if(buffer[0] != '#')
      mBuffer += buffer;
  }


  fclose(mpFile);
  mpFile = NULL;
  
  mpCurSection = mBuffer.c_str();
}

   






