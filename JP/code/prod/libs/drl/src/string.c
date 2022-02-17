/************************************************************************
 * Module:	DRL
 * Submodule:	STR
 * File:	
 * Function:	String Utilities
 * Author:	Christian Daher
 * Revision:	$Header$
 ************************************************************************/
#include "drlstd.h"		/* platform compatibility */
#include "drlerr.h"		/* DrlErrMsg */

#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include <ctype.h>

#include "drlstr.h"

static	int	_StringCompare(const void *, const void *);

/*f---------------------------------------------------------------------
 * String routines : allocate string.
 *                                               
 * <br><br>
 * Allocates and returns a pointer to a char string
 * of length <i> n1</i>. Returns NULL if failure.
 */

char*
DrlStrAlloc(int n1)
{
	return (char*) MALLOC((size_t) n1*sizeof(char));

}

/*f---------------------------------------------------------------------
 * Frees a char string pointer <i> s</i> allocated by the previous
 * routine.
 */

void
DrlStrFree(char *s)
{
	FREE((void*) s);
}


/*f---------------------------------------------------------------------
 * Allocates and returns a matrix of char 
 * of size <i> n1</i> $x$ <i> n2</i> (i.e. vector of
 * length <i> n1</i> of char strings of length <i> n2</i>).
 * Returns NULL if failure.
 */

char**
DrlStrArrayAlloc(int n1, int n2)
{
static	char	routine[] = "DrlStrArrayAlloc";
	char	**p ;
	int	i ;

	p = (char**) MALLOC((size_t) n1 * sizeof(char*)) ;
	if (p == NULL) {
		DrlErrMsg("%s: malloc of %d char* failed \n", routine, n1);
		return(NULL) ;
	}

	for (i=0; i<=n1-1; i++) {
		p[i] = (char*) MALLOC((size_t) n2 * sizeof(char)) ;
		if (p[i] == NULL) {
			DrlErrMsg("%s: malloc of %d char (row %d)  failed\n",
				routine, n2, i);
			while(--i >= 0) FREE((void*) p[i]);
			FREE((void*) p);
			return(NULL) ;
		}
	}
	return(p) ;
}

/*f---------------------------------------------------------------------
 * Frees a matrix <i> p</i> of char of size <i> n1</i> $x$ <i> n2</i>
 * allocated by the previous routine.
 */

void
DrlStrArrayFree(char **p, int n1)
{
	int	i;

	for (i=0; i<=n1-1; i++) FREE((void*) p[i]);
	FREE((void*) p) ;
}




/*f---------------------------------------------------------------------
 * Sorts alphebetically a vector of length <i> n1</i> of char strings.
 */

int
DrlStrArraySort(char **p, int n1)
{
	qsort((void*) p, (size_t) n1, (size_t) sizeof(char*),
		_StringCompare);
	return (0);
}


static	int
_StringCompare(const void *x1, const void *x2)
{
	return strcmp(*(const char**) x1, *(const char**) x2);
}


/*f-------------------------------------------------------------
 * This routine is a similar to the standard C routine
 * <i> strtok</i> that does not use static variables.
 * <i> s0</i> and <i> s1</i> provide working space for the routine.
 * <i> s0</i> should be sufficiently large, and <i> s1</i> a pointer 
 * to char pointer.
 * It differs from <i> strtok</i> in that it identifies strings
 * delimited by quotes or double quotes as tokens.
 * %A call to <i> StrToken(s, tok, s0, s1)</i> is equivalent
 * %to  a call to the standard C <i> strtok(s, tok)</i>.
 * On first call, pass s!=NULL. Then call the routine
 * with s=NULL.
 * See C manual for details.
 * \begin{verbatim}
 *          char    sTmp0[256],
 *                  *sTmp1;
 *          ...
 *  
 *          p = StrToken(s, tok, s0, &s1);
 *          / * equivalent to 
 *          p = strtok(s, tok); * /
 *          ...
 * 
 * \end{verbatim}
 * 
 */

char*
DrlStrToken(const char *str, const char *tok, char *s0, char **s1)
{
#define	VNEW_STR

#ifdef	VNEW_STR
	/* New Version catches strings */
	char *token;
#define	olds	(*s1)

	if (str == NULL) {
		if (olds == NULL) {
			return NULL;
		} else {
			s0 = olds;
		}
	} else {
		/* first call: copy into s0 */
		strcpy(s0, str);
	}

	/* advance until delimiter */
	/* Scan leading delimiters.  */
	s0 += strspn(s0, tok);
	if (*s0 == '\0') {
		olds = NULL;
		return NULL;
	}

	/* identify strings */
	if (*s0 == '"') {
		s0++;
		token = s0;
		while ((*s0 != '"') && (*s0 != '\0')) s0++;
		if (*s0 == '"') {
			*s0++ = '\0';
		}

		if (*s0 == '\0') {
			/* This token finishes the string.  */
			olds = NULL;
		} else {
			olds = s0;
		}

	} else {
		/* Find the end of the token.  */
		token = s0;
		s0 = strpbrk(token, tok);
		if (s0 == NULL) {
			/* This token finishes the string.  */
			olds = NULL;
		} else {
			/* Terminate token and make OLDS point past it.*/
			*s0 = '\0';
			olds = s0 + 1;
		}
	}
	return(token);
#undef	olds

#elif defined(VNEW)

	/* New Version*/

	char *token;
#define	olds	(*s1)

	if (str == NULL) {
		if (olds == NULL) {
			return NULL;
		} else {
			s0 = olds;
		}
	} else {
		/* first call: copy into s0 */
		strcpy(s0, str);
	}

	/* Scan leading delimiters.  */
	s0 += strspn(s0, tok);
	if (*s0 == '\0') {
		olds = NULL;
		return NULL;
	}

	/* Find the end of the token.  */
	token = s0;
	s0 = strpbrk(token, tok);
	if (s0 == NULL) {
		/* This token finishes the string.  */
		olds = NULL;
	} else {
		/* Terminate the token and make OLDS point past it.  */
		*s0 = '\0';
		olds = s0 + 1;
	}
	return(token);
#undef	olds

#elif defined(VOLD)

	char	*p, *q;

	if (s != NULL) {
		/* first call: copy into s0 */
		strcpy(s0, s);

		/* remove first tokens */
		for(p=s0; strchr(tok, *p) != NULL; p++);
		q = s0;
		while ((*q = *p) != '\0') {p++; q++;};
		
		*s1 = s0;

	}

	if (*s1 == NULL) return NULL;
	p = *s1;
	q = p ;
	while ((strchr(tok, *q) == NULL) && (*q != '\0')) q++;

	if (*q == '\0') {
		*s1 = NULL;
		return p;
	} else {
		*q = '\0';
		q++;
		while ((strchr(tok, *q) != NULL) && (*q != '\0')) q++;
		if (*q == '\0') *s1 = NULL;
		else *s1 = q;
		return p;
	}

#endif	/* VNEW, VOLD*/
}

/*f-------------------------------------------------------------
 * This routine scans for strings
 * in a character string pointer by <i> p</i>.
 * The read strings can be enclosed in quotes (') or double quotes ("),
 * which are removed.
 * On exit, the result string is put in <i> arg</i> and returned (if <i> arg</i> 
 * is NULL, a static copy is returned).
 * The pointer <i> p</i> is advanced to the remaining unread portion 
 * the scanned string. It should not be changed between
 * successive calls to the routine.
 * For example,
 * \begin{verbatim}
 *          char    *p, *q, *string = "1.23 'hello world'";
 *          ...
 *          / * p points string to be scanned * /
 *          p = &string[0];
 *          q = DrlStrScanString(&p, NULL);	/ * q is "1.23" * /
 *          / * p now points to the unread portion of string  * /
 *          ...
 *          q = DrlStrScanString(&p, NULL);	/ * q is "hello world" * /
 *          ...
 * \end{verbatim}
 */

char*
DrlStrScanString(char **p, char *arg)
{
const	char	EOS = '\0';		/* end of string */
static	char	arg1[256];
	char	*s, *q, c;

	arg = (arg != NULL ? arg : arg1);
	q = arg;
	s = *p;

	while (((c = *s++) != EOS) && (isspace(c)));

	if (c == EOS) return(NULL);


	if (c == '"') {
	        while (((c = *s++) != EOS) && (c != '"'))
			*q++ = c;
		if (c == EOS) return(NULL);
		*q = '\0';
	} else if (c == '\'') {
	        while (((c = *s++) != EOS) && (c != '\''))
			*q++ = c;
		if (c == EOS) return(NULL);
		*q = '\0';
	} else {
		*q++ = c;
	        while (((c = *s++) != EOS) && (!isspace(c)))
			*q++ = c;
		*q = '\0';
	}
	*p = s;
	return(arg);
}


/*f-------------------------------------------------------------
 * Inserts a char <i> c</i> at position pointed by <i> p</i> in character
 * string <i> s</i>. Returns <i> s</i>.
 */

char*
DrlStrInsertChar(char *s, char *p, char c)
{
	char	*q;
	for(q=s+(int)strlen(s);q>=p;q--)
		*(q+1) = *q;
	*p = c;
	return s;
}


/*f-------------------------------------------------------------
 * Removes the char pointed by <i> p</i> in string <i> s</i>.
 * Returns <i> s</i>.
 */

char*
DrlStrDeleteChar(char *s, char *p)
{
	char	*q;
	for(q=p; *(q+1) != '\0'; q++) *q = *(q+1);
	*q = '\0';

	return s;
}

/*f-------------------------------------------------------------
 * Replaces a series of characters <i> before</i> by
 * another series of character <i> after</i> in string <i> s</i>.
 * Returns <i> s</i>.
 */

char*
DrlStrReplace(char *s, char *before, char *after)
{
register char	*p;
register int	i;

	for (p=s; *p != '\0'; p++) {
		for (i=0; before[i] != '\0'; i++) {
			if (before[i] == *p) {
				*p = after[i];
				break;
			}
		}
	}
	return(s);
}

/*f-------------------------------------------------------------
 * Replaces any occurence of string <i> before</i>
 * by the string <i> after</i>.
 * Returns <i> s</i>.
 */

char*
DrlStrReplaceStr(char *s, char *before, char *after)
{
register char	*p, *q, *r;
register int	i;
	int	j;

	j = strlen(after) - strlen(before);
	p = s;
	while ((*p != '\0') && ((p = strstr(p, before)) != NULL)) {
		r = p + strlen(before);
		if (j > 0) {
			for(q=r+(int)strlen(r);q>=r;q--)
				*(q+j) = *q;
			for(i=0; after[i] != '\0';i++)
				p[i] = after[i];
		} else {
			for(q=r; *q != '\0'; q++)
				*(q+j) = *q;
			*(q+j) = '\0';
			for(i=0; after[i] != '\0';i++)
				p[i] = after[i];
		}
		p = p + strlen(after);
	}
	return(s);
}


/*f----------------------------------------------------------------------
 * Replaces all non printable characters by ' ' in the string <i> s</i>.
 * Returns <i> s</i>.
 */

char*
DrlStrReplaceNonPrintChar(char *s)
{
	char	*q;
	for (q = s; *q != '\0'; q++) {
		if (!isprint(*q)) *q = ' ';
	}
	return (s);
}



/*f-------------------------------------------------------------
 * Remove any occurence of characters of <i> ct</i> from  <i> s</i>
 * Returns <i> s</i>.
 */

char*
DrlStrRemoveChars(char *s, char *ct)
{
	char	*p, *q;
	for (q = ct; *q != '\0'; q++) {
		while ((p = strrchr(s, *q)) != NULL)
			DrlStrDeleteChar(s, p);
	}
	return(s);
}

#ifdef	UNIX
/*--------------------------------------------------------------
 *
 */
char*
DrlStrNCpy(char *s, char *ct, int n)
{
	char	*p, *q;
	int	i;

	for (p=ct, i=0, q=s; (*p != '\0') && (i <= n-1); i++, p++, q++)
		*q = *p;
	*q = 0;


	return s;
}
#endif

/*f-------------------------------------------------------------
 * Copy the <i> n</i> rightmost character of <i> ct</i> in <i> s</i>.
 * Returns <i> s</i>.
 */

char*
DrlStrRightNCpy(char *s, char *ct, int n)
{
	char *p;

	if (strlen(ct) <= (unsigned int)n) {
		strcpy(s, ct);
	} else {
		p = ct + (strlen(ct) - n);
		strcpy(s, p);
	}
	return s;
}


/*f-------------------------------------------------------------
 * Cuts string <i> s</i> after <i> n</i> characters. If the lenngth
 * of <i> s</i> is less than <i> n</i>, does nothing.
 * Returns <i> s</i>.
 */

char*
DrlStrCut(char *s, int n)
{
	register int i;
	for (i=0; (i <= n-1) && (s[i] != '\0'); i++);
	s[i] = '\0';

	return(s);
}

/*f-------------------------------------------------------------
 * Like an <i> ssprintf</i>, but returns a pointer to a static
 * string.
 */

char*
DrlStrPrint(char *fmt, ...)
{
#undef	MAX_IDX
#define	MAX_IDX	8
	va_list	ap;
static	char	tmp[MAX_IDX][256] ;
static	int	tmpIdx=0;
	char	*s ;

	s = tmp[tmpIdx];
	tmpIdx++;
	if (tmpIdx > MAX_IDX-1) tmpIdx=0;


	va_start(ap, fmt);
	vsprintf(s, fmt, ap);
	va_end(ap);

	return(s)  ;
#undef	MAX_IDX
}




/*f---------------------------------------------------------------------
 * Copies the basename of a pathname <i> ct</i> in <i> s</i>.
 * If <i> s</i> is NULL, returns a pointer to a static string.\\
 * <i> Example:</i> for <i> /home/username/example/ex.foo</i>,
 * returns <i> ex.foo</i>.
 */

char*
DrlGetBaseName(char *s, char *ct)
{
static	char	buf[256];
	char	*p;

	s = (s != NULL ? s : buf);

	/* go back until first '/' encountered */
	p = ct + strlen(ct);
	while ((p != ct) && (*p != '/')) p--;
	if (*p == '/') p++;

	strcpy(s, p);
	return s;
}

/*f---------------------------------------------------------------------
 * Copies the suffix of a pathname <i> ct</i> in <i> s</i>.
 * Return the empty string if no suffix.
 * If <i> s</i> is NULL, returns a pointer to a static string.
 * <i> Example:</i> for <i> /home/username/example/ex.foo</i>,
 * returns <i> foo</i>.
 */

char*
DrlGetSuffix(char *s, char *ct)
{
static	char	buf[32];
	char	*p;

	s = (s != NULL ? s : buf);

	p = ct + strlen(ct);
	while ((p != ct) && (*p != '.')) p--;
	if (*p == '.') {
		p++;
		strcpy(s, p);
	} else {
		strcpy(s, "");
	}
	return s;
}



/*f---------------------------------------------------------------------
 * Copies the root of a pathname <i> ct</i> in <i> s</i>.
 * Return the basename string if no suffix (same as DrlGetBaseName).
 * If <i> s</i> is NULL, returns a pointer to a static string.
 * <i> Example:</i> for <i> /home/username/example/ex.foo</i>,
 * returns <i> ex</i>.
 */

char*
DrlGetRootName(char *s, const char *ct)
{
static	char	buf[32];
	char		*p;
	const char	*cq;
	s = (s != NULL ? s : buf);

	cq = ct;
	if ((cq = strrchr(ct, '/')) == NULL)
		cq = ct;
	else
		cq++;

	strcpy(s, cq);
	if ((p = strrchr(s, '.')) != NULL) p[0] = '\0';

	return s;
}



/*f---------------------------------------------------------------------
 * Copies the directory path of a pathname <i> ct</i> in <i> s</i>.
 * Return the empty string if no suffix.
 * If <i> s</i> is NULL, returns a pointer to a static string.
 * <i> Example:</i> for <i> /home/username/example/ex.foo</i>,
 * returns <i> /home/username/example</i>.
 */

char*
DrlGetPathDir(char *s, char *ct)
{
	char	*p;
	int	n;

	strcpy(s, ct);

	/* remove tail '/' if any */
	while ((s[n=(int)strlen(s)-1] == '/') && (n != 0)) s[n] = '\0';
	if (n == 0) return s;

	/* go back until '/' encountered */
	p = s + n;
	while ((p != s) && (*p != '/')) p--;

	if ((p != s) && (*p == '/')) 
		*p ='\0';
	else if ((p == s) && (*p == '/')) 
		strcpy(s, "");
	else
		strcpy(s, ".");
	return s;
}



/*f-------------------------------------------------------------
 * Converts a string <i> s</i> to lower case.
 * Returns <i> s</i>.
 */

char*
DrlStrToLower(char *s)
{
register char	*p;
	for (p=s; *p != '\0'; p++) 
		*p = tolower(*p);
	return(s);
}

/*f-------------------------------------------------------------
 * Converts a string <i> s</i> to upper case.
 * Returns <i> s</i>.
 */

char*
DrlStrToUpper(char *s)
{
register char	*p;
	for (p=s; *p != '\0'; p++) 
		*p = toupper(*p);
	return(s);
}

/*f-------------------------------------------------------------
 * Like a <i> strcmp</i> but case indifferent.
 */

int
DrlStrCmpCaseIndif(const char *cs, const char *ct)
{
register const char	*p, *q;
	for (p=cs, q=ct; (*p != '\0') && (*q != '\0'); p++, q++) 
	    if (toupper(*p) != toupper(*q)) return(1);
	return ((*p == '\0') && (*q == '\0') ? 0 : 1);
}


/*f-------------------------------------------------------------
 * Checks if the string <i> cs</i> starts by the string <i> ct</i>
 * (case indifferent). Ignores initial space characters in <i> cs</i>.
 * Returns 1 (TRUE) of 0 (FALSE).
 */

int
DrlStrStartsByCaseIndif(const char *cs, const char *ct)
{
register const char	*p, *q;
	for (p=cs; isspace(*p); p++);

	for (q=ct; (*p != '\0') && (*q != '\0'); p++, q++) 
	    if (toupper(*p) != toupper(*q)) return(FALSE);
	return ((*q == '\0') ? TRUE : FALSE);
}





