/*************************************************************************** 
*  Name:        ddlmaker/src/column.cpp                              
*  Purpose:	Part of the ddlmaker utility                                 
*  Author:      Cosmin Cremarenco                                      
*  Created:     2005-06-07                                                   
*  RCS-ID:      $Id: column.cpp,v 1.2 2005/07/06 09:03:02 cosmin Exp $    
*  Copyright:   (c) 2005 Trilemma LLP                                        
***************************************************************************/
#include "column.h"

std::string Column::ConcatenateNames(const TColumnPtrList& list,
		const std::string& delimiter)
{
	std::string toReturn;
		
	for(TColumnPtrList::const_iterator iterator = list.begin();
			iterator != list.end(); ++iterator)
	{
		toReturn.append((*iterator)->GetName());
		if(std::distance(iterator, list.end()) > 1)
			toReturn.append(delimiter);
	}

	return toReturn;
}

void Column::CodeGenerator(Backend& backend)
{
	backend.HandleColumn(*this);
}
