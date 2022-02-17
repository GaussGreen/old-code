/*************************************************************************** 
*  Name:        ddlmaker/src/table.cpp                              
*  Purpose:	Part of the ddlmaker utility                                 
*  Author:      Cosmin Cremarenco                                               
*  Created:     2005-06-07                                                   
*  RCS-ID:      $Id: table.cpp,v 1.5 2005/08/24 08:21:59 cosmin Exp $    
*  Copyright:   (c) 2005 Trilemma LLP                                        
***************************************************************************/
#include "table.h"

/* helper class */
class PrintTableNameToErr
{
public:
	void operator() (TTablePtr table)
	{
		std::cerr << table->GetName() << " ";
	}
};

class CallCheckCoherency
{
public:
	void operator() (TConstraintBasePtr constraint)
	{
		constraint->CheckCoherency();
	}
};

TColumnPtr Table::GetColumnByName(const std::string& name) const
{
	TColumnPtr toReturn;
	TColumnPtrList::const_iterator iterator;

	for(iterator = m_columns.begin(); iterator != m_columns.end();
			++iterator)
	{
		/* 
		 * we don't need to check that there is no other
		 * column with the same name in the table. The 
		 * later is done when the column is added to
		 * the table.
		 */
		if((*iterator)->GetName() == name)
		{
			toReturn = *iterator;
			break;
		}
	}

	return toReturn;
}

TConstraintBasePtr Table::GetConstraintByName(const std::string& name) const
{
	TConstraintBasePtr toReturn;
	TConstraintBasePtrList::const_iterator iterator;

	for(iterator = m_constraints.begin(); iterator != m_constraints.end();
			++iterator)
	{
		if((*iterator)->GetName() == name)
		{
			toReturn = *iterator;
			break;
		}
	}

	return toReturn;
}

TConstraintBasePtr Table::GetConstraintByName(const TTablePtrList& list,
	const std::string& name)
{
	TTablePtrList::const_iterator iterator;
	TConstraintBasePtr toReturn;
	
	for(iterator = list.begin(); iterator != list.end();
			++iterator)
	{
		toReturn = (*iterator)->GetConstraintByName(name);
		if(toReturn)
			break;
	}

	return toReturn;
}

TTablePtr Table::GetTableByName(const TTablePtrList& list,
	const std::string& name)
{
	TTablePtr toReturn;
	TTablePtrList::const_iterator iterator;

	for(iterator = list.begin(); iterator != list.end();
			++iterator)
	{
		if((*iterator)->GetName() == name)
		{
			toReturn = *iterator;
			break;
		}
	}

	return toReturn;
}

TTablePtr Table::GetTableByColumnName(const TTablePtrList& list,
	const std::string& name)
{
	TTablePtrList tables;
	TTablePtrList::const_iterator iterator;
	TTablePtr toReturn;

	for(iterator = list.begin(); iterator != list.end();
			++iterator)
	{
		TColumnPtr temp;

		temp = (*iterator)->GetColumnByName(name);
		
		/* remember the duplicate columns */
		if(temp)
			tables.push_back(*iterator);
	}
	
	if(tables.size() > 0)
	{
		if(tables.size() > 1)
		{
			std::cerr << "Warning: More than two tables contains same-named column "
				<< "Returning the first."
				<< std::endl;
		}
		toReturn = tables.front();
	}
	
	return toReturn;
}

TColumnPtr Table::GetColumnByName(const TTablePtrList& list,
	const std::string& name)
{
	TColumnPtrList columns;
	TTablePtrList::const_iterator iterator;
	TColumnPtr toReturn;
	
	for(iterator = list.begin(); iterator != list.end();
			++iterator)
	{
		TColumnPtr temp;

		temp = (*iterator)->GetColumnByName(name);
		
		/* remember the duplicate columns */
		if(temp)
			columns.push_back(temp);
	}

	if(columns.size() > 0)
	{
		if(columns.size() > 1)
		{
			/* 
			 * TODO: elaborate on this message,
			 * say which tables hold the duplicates
			 */
			std::cerr << "Partially qualified search of columns found "
				" two or more tables that contain columns with the same name! "
				<< "Returning the first."
				<< std::endl;
		}
		toReturn = columns.front();
	}
	
	return toReturn;
}

bool Table::AliasedTypeIsUsed(const std::string& name) const
{
	TColumnPtrList::const_iterator iterator;
	bool toReturn = false;
	
	for(iterator = m_columns.begin(); iterator != m_columns.end();
			++iterator)
	{
		if((*iterator)->GetType()->GetAlias() == name)
		{
			toReturn = true;
			break;
		}
	}

	return toReturn;
}

bool Table::AliasedTypeIsUsed(const TTablePtrList& list,
	const std::string& alias)
{
	TTablePtrList::const_iterator iterator;
	bool toReturn = false;
	
	for(iterator = list.begin(); iterator != list.end();
			++iterator)
	{
		if(toReturn = (*iterator)->AliasedTypeIsUsed(alias))
			break;
	}
	return toReturn;
}
	
bool Table::CheckCycle(TTablePtr vertex, TTablePtrList& domain,
	TTablePtrList& context)
{
	bool toReturn = true;
	TTablePtrList::iterator pos;

	pos = std::find_if(context.begin(), context.end(), 
			SameTableName(vertex->GetName())); 
	if(pos == context.end())
		context.push_back(vertex);
	else
	{
		std::cerr << "A cycle was detected: ";
		std::for_each(context.begin(), context.end(),
				PrintTableNameToErr());
		std::cerr << std::endl;
		return false;
	}
	
	TConstraintBasePtrList::iterator iterator;
	for(iterator = vertex->GetConstraints().begin();
			iterator != vertex->GetConstraints().end();
			++iterator)
	{
		/* if the GetSecondaryColumns list isn't empty that means
		 * that there is a link between two tables (a fk without 
		 * question)
		 */
		TColumnPtrList secondaryColumns = (*iterator)->GetSecondaryColumns();

		if(secondaryColumns.size() > 0)
		{
			/* 
			 * even if this fk has two son columns and two father
			 * columns, the father columns are in the same table
			 */
			TColumnPtr destColumn = secondaryColumns.front();
			
			toReturn &= CheckCycle(destColumn->GetOwnerTable(), domain, context);
		}
	}	

	/* unload context */
	context.pop_back();
	
	/* couldn't find a cycle for this table, unload it from the domain */
	std::remove_if(domain.begin(), domain.end(),
			SameTableName(vertex->GetName()));
	
	return toReturn;
}

void Table::CheckCoherency() const
{
	std::for_each(m_constraints.begin(), m_constraints.end(),
			CallCheckCoherency());
}

void Table::CodeGenerator(Backend& backend)
{
	backend.HandleTable(*this);
}
