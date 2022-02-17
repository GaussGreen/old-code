/***************************************************************
 * Module:	PenGuin
 * Submodule:	
 * File:	kexcept.h
 * Function:	Standard Definition and Include files.
 * Author:	Christian Daher
 * Revision:	$Header$
 ***************************************************************/
#ifndef	_kexcept_H
#define	_kexcept_H


/**
 * Error handling class.
 */

class	KFailure {
public:
	/** Default constructor. */
	KFailure(int level = -1); 

	/** Constructor that prints a message on the error log. */
	KFailure(const char *fmt, ...);

private:
				/** Error level. */
	int	errLevel;
};



/*
 * Useful macros.
 */



#define	IF_FAILED_THROW(statement)	\
	{if (statement != SUCCESS) throw KFailure("Condition `%s' failed.\n", \
	 #statement);}

#define	ASSERT_OR_THROW(statement)	\
	{if (!(statement)) throw KFailure("Assertion `%s' failed.\n",\
	 #statement);}

#define	ASSERT_OR_BUG(statement)	\
	{if (!(statement)) { throw KFailure(\
	"%s: `%s' failed at (%s, %d).\n",\
	 routine, #statement,__FILE__,__LINE__); }}

#define	SUCCESS_OR_THROW(statement)	\
	{if (statement != SUCCESS) {throw KFailure(\
	"%s: `%s' returned FAILURE.\n",\
	 routine, #statement);}}


#define	THROW_NA	\
	{ throw KFailure("%s, %d: not implemented.\n", __FILE__,__LINE__);}

#define	THROW_BUG	\
	{ throw KFailure("%s, %d: bug !.\n", __FILE__,__LINE__);}



#endif // _kexcept_H


