/*-------------------------------------------------------------------------*
 
   File: %M%
   Path: %P%
   Description: implantation de la classe CCString 
   Created: 98/06/02
   Author: Charles-Emmanuel MUSY
   Modified: %E% %U%
   Last maintained by: Charles-Emmanuel MUSY
   Revision: %I%
 *-------------------------------------------------------------------------*/
 

#include <CCcommon.h>
SCCS_ID (CCString_c_SccsId, "%W%, modified %E%");

#define CCString_c
#include "CCString.h"



/*-----------------------------------------------------------------------------*/


char* ARM_STR_DUP(const char* inStr)
{
    char* res = NULL;

    try
    {
        res = new char [strlen(inStr)+1];

        res[0] = '\0';

        strcpy(res, inStr);
    }

    catch(...)
    {
        throw "ARM_STR_DUP: Memory Problem";
    }

    return(res);
}



CCString::~CCString ()
{
	if (str)
	{
	   delete [] str;

	   str = NULL;
	}
}



void CCString::Set(const char* n_str)
{
	if (str)
	{
	   delete [] str;
		
       str = NULL;
	}

	if (n_str)
	{
	   len = strlen(n_str);
        	
       str = ARM_STR_DUP(n_str);;
	}
	else
	{
	   len = 0;
	   
       str = NULL;
	}
}



void CCString::Set(int n, const char* n_str)
{
	if (str)
	{
	   delete [] str;

	   str = NULL;
	}

	if ((n_str) && ( n > 0 ))
	{
		if ( n >= strlen(n_str) )
		{
		   Set(n_str);
		}
		else
		{
		   len = n; 
        		
           CC_ASSERT (str = new char [n+1], 
                      MSG_F_MIS_INCONSISTENT_VALUE);
			
           strncpy(str, n_str, n);
			
           str[n] = '\0';
		}
	}
	else
	{
		len = 0;
		str = NULL;
	}
}



double CCString::XL_value_convert()
{
	double result = 0.0;

	char* point_pos = strchr(str, COMMA);

    if (point_pos)
	{
		*point_pos = POINT;

		result = atof (str);
	}

	return(result);
}



void CCString::Replace (const char* n_str, char mask)
{
	char* end_str;
	char* begin_str;

	int use_n_str = 0;

	if (str)
	{
		int i = 0;
	
        while (( str[i] != mask ) && ( str[i] != '\0' )) 
		{
			i++;
		}
		
        if ( str[i] == '\0' )
		{ 
			return;
		}

        if ( i == 0 )
		{
            CC_ASSERT ((begin_str = ARM_STR_DUP (n_str)), MSG_F_MIS_INCONSISTENT_VALUE);
			
            end_str = new char [len-i];
			
            int j = 0;
			
            i++;
		
            while( str[i] != '\0' )
			{
				end_str [j] = str [i];
				
                j++; 
                i++;
			}
			
            end_str[j] = '\0';
		}
		else if ( i == len-1 )
		{
			begin_str = new char [i+1];
			
            strncpy(begin_str, str, i);
			
            begin_str[i] = '\0';
			
            CC_ASSERT((end_str = ARM_STR_DUP (n_str)), MSG_F_MIS_INCONSISTENT_VALUE);
		}
		else
		{
			begin_str = new char [i+1];

			strncpy(begin_str, str, i);
			
            begin_str[i] = '\0';
			
            end_str = new char [len-i];
			
            int j = 0;
			i++;
			
            while ( str[i] != '\0' )
			{
				end_str[j] = str [i];

				j++; 
                i++;
			}
			
            end_str[j] = '\0';
			use_n_str = 1;
		}

		Set(begin_str);
		
        if ( use_n_str == 1 )
		{
		   *this += n_str;
		}

		*this += end_str;
		
        delete [] begin_str;
		delete [] end_str;
	} 
}	



void CCString::Replace(const char* strMask, const char* n_str)
{
	char* pdest = strstr(str, strMask);

	if ( pdest == NULL )
	   return;

	int result = pdest-str;

	char* begin_str = new char [strlen(str)+strlen(n_str)];

	strncpy(begin_str, str, result);
	
    begin_str[result] = '\0';
	
	strcat(begin_str, n_str);

	char* temp;

	temp = str+result+strlen(strMask);

	strcat(begin_str, temp);


	Set(begin_str);

	delete [] begin_str;
}



CCString& CCString::operator= (const CCString& s)
{
	if ( this == &s )
	{
	   return(*this);
	}

	Set(s);

	return(*this);
}



CCString& CCString::operator= (const char* s)
{
	Set(s);

	return(*this);
}



void CCString::toUpper()
{
	for (int i = 0; i < len; i++)
	{
	    str[i] = toupper(str[i]);
	}
}



CCString& CCString::operator+= (const CCString& s)
{
	if (( len > 0 ) && ( s.len > 0 ))
	{
		len += s.len;
		
        char* p = new char [len+1];

	
        strcpy(p, str);
		
        strcat(p, s.str);

        delete [] str;

        str = p;
	}

	return(*this);
}



CCString CCString::operator+ (const CCString& s) const
{
	CCString result = *this;

	result += s;
	
    return(result);
}



int CCString::operator== (const CCString& s) const
{
	if ((!str) || (!s)) 
    {
	   return( str == (const char *) s );
	} 
    else 
    {
	   return( strcmp(str, s) == 0 );
	}
}



int CCString::operator== (const char* s) const
{
	if ((! str) || (! s)) 
    {
		return (str == s);
	} 
    else 
    {
		return( strcmp(str, s) == 0 );
	}
}



int CCString::operator== (char* s) const
{
	if ((! str) || (! s)) 
    {
	   return( str == s );
	} 
    else 
    {
	   return( strcmp(str, s) == 0 );
	}
}



int CCString::operator!= (const CCString& s) const
{
	return(! ((CCString &)str == s));
}



std::ostream& operator<< (std::ostream& os, const CCString& s)
{
	if (!s)
	{
		return os << "<empty>";
	}
	
	return os << s.str;
}



char* CCtrim_right(char* s)
{
	if (!s) return (s);

	int i = strlen (s);

	while (--i >= 0) if (isspace (s[i])) s[i] = '\0'; else return (s);

	return (s);
}



char* CCString::c_str() const
{ 
	if (str)
	{
	   char* result = new char[strlen(str)+1];
		
       strcpy(result, str);

	   return(result);
	}
	else
	   return(NULL); 
}



#ifdef STL_WIN32
void CCString::Parser (char separator, std::vector<CCString>& list)
#else
void CCString::Parser (char separator, vector<CCString>& list)
#endif	// STL_WIN32
{
	list.clear ();
	if (str)
	{
		char* p = strchr(str, separator);
		char* temp = str;
		int d = 0;

		while(p)
		{
			d = strlen(temp)-strlen(p); 
			
            if ( d > 0 )
			{
				list.push_back (CCString(d, temp));
			}

			p++;
			
            temp = p;
			
            p = strchr(temp, separator);
		}

		if ((temp) && ( strlen(temp) > 0 ))
		{
		   list.push_back (CCString (temp));
		}
	}
}


#ifdef STL_WIN32
void CCString::ParserWithEmpty (char separator, std::vector<CCString>& list)
#else
void CCString::ParserWithEmpty (char separator, vector<CCString>& list)
#endif	// STL_WIN32
{
	list.clear();

	if (str)
	{
		char* p = strchr(str, separator);

		char* temp = str;
		
        int d = 0;
		
        while(p)
		{
			d = strlen(temp)-strlen(p); 

			if ( d > 0 )
			   list.push_back(CCString(d, temp));
			else
			   list.push_back((CCString) "");

			p++;
			temp = p;

			p = strchr(temp, separator);
		}

		if (temp) 
		{
			if ( strlen (temp) > 0 )
			   list.push_back(CCString(temp));
			else
			   list.push_back((CCString) "");
		}
	}
}


//// Conversion to a standard STL string
std::string CCSTringToSTLString( const CCString& s )
{
	char* c = s.c_str();

	if (!c)
	   return(NULL);

	std::string stlString(c);

	delete [] c;
	
    return(stlString);
}



bool CCString::Contain(const CCString& chaine)
{
	CCString Tmp = *this;

	int szbef = Tmp.GetLen();
	
    Tmp.Replace(chaine, "");

	int szbef2 = Tmp.GetLen();

	if ( szbef > szbef2 ) 
       return(true);

	return(false);
}

/*------------------------------------------------------------------------------------*/
// EOF %M%