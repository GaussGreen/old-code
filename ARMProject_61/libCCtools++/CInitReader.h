

#ifndef _CInitReader_H
#define _CInitReader_H



#ifdef WIN32
	#pragma warning(disable:4786)
	#include <map>
	#include<iosfwd>
    #include <string>
	using namespace std;

#else
    #include<map.h>
    #include<iostream.h>

#endif




class CInitReaderException {
 public:
	 string mErrorMsg;

  enum ExceptionNo {
    FILE_NOT_FOUND,
    SECTION_NOT_FOUND,
    PROPERTY_NOT_FOUND,
    CONVERSION_ERROR,
    INVALID_SYNTAX,
    MAX_EXCEPTION_NO
  };
  CInitReaderException(ExceptionNo aNo, const char* pTextToAppend = NULL);
  CInitReaderException(const char*);
  
  const char* GetMessage() const;
  operator const char* () const ;
 protected:
  

};


extern const char* CInitReaderExceptionMsg[CInitReaderException::MAX_EXCEPTION_NO];


class CInitReader {
 public:
  // if argument "shouldContinue" is set to 1, if the reader is not able to find a section or a property, an error;
  // message will be printed and the execution will continue. Else,by default, an Exception will be raised.
  CInitReader(const char* pFileName, int souldContinue = 0);
  virtual void Init() = NULL;
  virtual string ToString() = NULL;
  const char* GetFileName() { return (const char*)mFileName.c_str(); }

 protected:
  FILE* mpFile;
  int   mContinue;
  string mFileName;
  string mBuffer;
  const char* mpCurSection;
  int mCurSectionPos;
  int mCurSectionEnd;
  string mSubSection;

  // select a specified section into the ini file
  void SetSection(const char* pSectionName);
  string GetString(const char* pProperty);
  float GetFloat(const char* pProperty);
  int GetInt(const char* pProperty);
  void Read();

  
};





#endif



