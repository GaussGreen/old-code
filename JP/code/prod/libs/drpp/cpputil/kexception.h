// kexception.h: interface for the kexception class.
//
//////////////////////////////////////////////////////////////////////
#if !defined(KEXCEPTION_H)
/**@#-*/
#define KEXCEPTION_H
/**@#+*/

#include "kplatdep.h"
#include STDARG_H
#include STRING_H
#include STDIO_H
#include SSTREAM_H


/**	Exception class

  Stores a textual error messages where messages are appended using operator<<.

  Example:

  <code><dir>	 
  throw KException ("This ") << "is an " << "error msg.\n"; 
  </dir></code>

  Note:	EXCEPTION is a platform independent exception base class with STL interface.
  @see operator<<
*/
class KException : public EXCEPTION
{
public:
	/** Creates an empty error message */
	KException () {mMsg = new char[1]; mMsg[0] = 0;}

	/** Creates initial error message from char 
		@param msg Initial error message */
	KException (const char* msg); 


	/**@#-*/
	~KException() throw();
	
	KException (const KException& a);
	KException& operator=(const KException& a);
	
	/** Adds a text message to the current message 

		Public only due to deficiencies in SunCC-4.2 */
	KException& add(const char* s, int num);
	/**@#+*/
	
	/** Returns error message
		@returns Error message */
	const char* what() const throw();
	
protected:
	char* mMsg;
};

/**@#-*/
typedef KException DRException;
/**@#+*/

/** Appends textual string to exception message. 

	Note:  This method is not part of the exception class, since Sun-CC4.2 does not support
	templated member functions.
	@see KException
	*/
template <class T> inline
KException& operator<<(KException& e, const T& t)
{
	OSTRINGSTREAM out;
	out << t;
	
//#if defined (_MSC_VER)
	string temp = out.str();
	e.add(temp.c_str(), temp.size());
//#else
//	char* temp = out.str();
//	e.add(temp, out.pcount());
//#endif
	
//#if defined (__SUNPRO_CC)
//	delete [] temp;
//#endif
	
	return e;
}

template <class T> inline
KException& operator<<(const KException& e1, const T& t)
{
  
  KException e(e1);

	OSTRINGSTREAM out;
	out << t;

	string temp = out.str();
	e.add(temp.c_str(), temp.size());
	
	return e;
}


inline
KException::KException (const char* msg) :mMsg (NULL)
{
	mMsg = new char[1]; mMsg[0] = 0;
	add (msg, strlen(msg));
}

inline
KException::~KException() throw() {delete [] mMsg;}

inline
KException::KException (const KException& a)
{
	mMsg = new char [strlen(a.mMsg) + 1];
	strcpy (mMsg, a.mMsg);
}

inline
KException& KException::operator=(const KException& a)
{
	if (&a == this) return *this;
	
	delete [] mMsg; 
	mMsg = new char[strlen(a.mMsg)+1]; 
	strcpy (mMsg, a.mMsg); 
	return *this;
}

inline
KException& KException::add(const char* s, int num)
{
	int size1 = strlen(mMsg); 
	int size2 = num;
	char* temp = mMsg;
	mMsg = new char [size1+size2+1];
	int i;
	for (i = 0; i < size1; i++) mMsg[i] = temp[i];
	for (i = 0; i < size2; i++) mMsg[size1+i] = s[i];
	mMsg[size1+size2] ='\0';
	delete [] temp;
	
	return *this;
}


inline
const char* KException::what() const throw() {return mMsg;}


#endif // !defined(KEXCEPTION_H)



