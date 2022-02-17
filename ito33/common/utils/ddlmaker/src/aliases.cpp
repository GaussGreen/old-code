/*************************************************************************** 
*  Name:        ddlmaker/src/alias.cpp                              
*  Purpose:	Part of the ddlmaker utility                                 
*  Author:      Cosmin Cremarenco                                            
*  Created:     2005-06-07                                                   
*  RCS-ID:      $Id: aliases.cpp,v 1.2 2005/10/31 18:15:37 cosmin Exp $    
*  Copyright:   (c) 2005 Trilemma LLP                                        
***************************************************************************/
#include <algorithm> 
#include <iterator>
#include <iostream>

#include "knownattributes.h"
#include "fundamentaltypes.h"
#include "aliases.h"
#include "variable.h"
#include "tools.h"
#include "column.h"

void AliasesBackend::HandleColumn(Column& column)
{
	std::string columnName = column.GetName();
	std::ostringstream buffer;
	
	m_output << column.GetOwnerTable()->GetName()
		<< "." << columnName 
		<< std::endl;
}

void AliasesBackend::HandleTable(Table& table)
{
	m_output << table.GetName() << "."
		<< std::endl;

	TColumnPtrList::const_iterator tempIterator;
	for(tempIterator = table.GetColumns().begin();
			tempIterator != table.GetColumns().end();
			++tempIterator)
	{
		(*tempIterator)->CodeGenerator(*this);
	}
}

std::string AliasesBackend::GetOutput()
{
	return m_output.str();
}

