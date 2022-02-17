/*************************************************************************** 
*  Name:        ddlmaker/include/validate.cpp                                 
*  Purpose:	Implementation of functions that validate the input.                                 
*  Author:      Cosmin Cremarenco                                            
*  Created:     2005-09-05                                                   
*  RCS-ID:      $Id: validate.cpp,v 1.1 2005/09/12 08:51:36 cosmin Exp $    
*  Copyright:   (c) 2005 Trilemma LLP                                        
***************************************************************************/

#include <iostream>
#include <fstream>

#include "validate.h"

static TKeywordList ms_keywordList;

bool 
NocaseCompare(char c1, char c2)
{
	return toupper(c1) == toupper(c2);
}

void
LoadKeywordList(const std::string& fileName)
{
	std::ifstream file(fileName.c_str());

	if ( ! file )
	{
		std::cerr << "Warning: Couldn't open file keywords file: " 
			<< fileName;
	} 
	else
	{
		std::string line;
		while ( std::getline(file, line) )
		{
			ms_keywordList.push_back(line);
		}
	}
}

bool
KeywordListContainsName(const std::string& name)
{
	TKeywordList::const_iterator begin = ms_keywordList.begin();
	TKeywordList::const_iterator end = ms_keywordList.end();
	TKeywordList::const_iterator pos;
	std::string::size_type size = name.size();
	bool inSet = false;
	
	for ( pos = begin; pos != end; ++pos )
	{
		if ( size == pos->size() && 
			equal(name.begin(), name.end(),
				pos->begin(),
				NocaseCompare) )
		{
			inSet = true;
			break;
		}
	}
	
	return inSet;
}

bool
ReservedKeyword(const std::string& name)
{
	if ( KeywordListContainsName(name) )
	{
		std::cerr << "Warning: the keyword \"" << name << "\"" 
			<< " is reserved!" << std::endl;

		return true; 
	}
	else
		return false;
}

