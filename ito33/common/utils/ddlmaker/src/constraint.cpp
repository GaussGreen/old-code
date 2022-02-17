/*************************************************************************** 
*  Name:        ddlmaker/src/constraint.cpp                              
*  Purpose:	Part of the ddlmaker utility                                 
*  Author:      Cosmin Cremarenco                                              
*  Created:     2005-06-07                                                   
*  RCS-ID:      $Id: constraint.cpp,v 1.4 2005/08/24 08:21:59 cosmin Exp $    
*  Copyright:   (c) 2005 Trilemma LLP                                        
***************************************************************************/
#include "constraint.h"
#include "knownattributes.h"
#include "type.h"

void ConstraintBase::AddPrimaryColumn(TColumnPtr column)
{
	std::cerr << "You cannot add a primary column to a base constraint."
		<< "Aborting!" << std::endl;
	exit(1);
}

void ConstraintBase::AddSecondaryColumn(TColumnPtr column)
{
	std::cerr << "You cannot add a secondary column to a base constraint."
		<< "Aborting!" << std::endl;
	exit(1);
}

TColumnPtrList ConstraintBase::GetPrimaryColumns() const
{
	std::cerr << "No primary columns in a base constraint."
		<< "Aborting!" << std::endl;
	exit(1);
}

TColumnPtrList ConstraintBase::GetSecondaryColumns() const
{
	std::cerr << "No secondary columns in a base constraint."
		<< "Aborting!" << std::endl;
	exit(1);
}

void ConstraintNull::AddPrimaryColumn(TColumnPtr column)
{
	if(m_columns.size())
	{
		std::cerr << "Warning, constraint already bound to a column."
			<< "Rebounding!" << std::endl;
	}
	m_columns.push_back(column);
}

void ConstraintNull::CodeGenerator(Backend& backend)
{
	backend.HandleNotNull(*this);
}

void ConstraintNull::AddSecondaryColumn(TColumnPtr column)
{
	std::cerr << "You cannot add a secondary column to a null/not null constraint."
		<< "Aborting!" << std::endl;
	exit(1);
}

TColumnPtrList ConstraintNull::GetPrimaryColumns() const
{
	return m_columns;
}

TColumnPtrList ConstraintNull::GetSecondaryColumns() const
{
	/* 
	 * some code relies on the fact that
	 * an empty list is returned in the case when
	 * there are no secondary columns. We do this even
	 * if in the case of a NotNull constraint secondary
	 * columns don't make any sense.
	 */
	return TColumnPtrList();
}

void ConstraintNull::CheckCoherency() const
{
	TColumnPtr column = m_columns.front();

	/* we will issue a warning only when the domain of the type
	 * has the NOT NULL attribute set
	 */
	
	/* first check if a notNull attribute was set on the type */
	TVarBasePtr notNullVal = GET_ATTR(*(column->GetType()), NOT_NULL_ATTRIBUTE);
	if( notNullVal && IsTrue(*notNullVal) )
	{
		if( NotNull() )
		{
			std::cerr << "Warning: column " << column->GetName() 
				<< " has both "
				<< " a NOT NULL constraint and is of a type " 
				<< column->GetType()->GetAlias() 
				<< " which doesn't allow NULL values" << std::endl;
		}
		else
		{
			std::cerr << "Warning: column " << column->GetName() 
				<< " accepts NULL values BUT "
				<< " is of type " 
				<< column->GetType()->GetAlias() 
				<< " which doesn't allow NULL values" << std::endl;
		}
	}
}

void ConstraintUnique::AddPrimaryColumn(TColumnPtr column)
{
	m_columns.push_back(column);
}

void ConstraintUnique::AddPrimaryColumnFromList(const TColumnPtrList& list)
{
	std::copy(list.begin(), list.end(), std::back_inserter(m_columns));
}

void ConstraintUnique::AddSecondaryColumn(TColumnPtr column)
{
	std::cerr << "You cannot add a secondary column to a unique constraint."
		<< "Aborting!" << std::endl;
	exit(1);
}

TColumnPtrList ConstraintUnique::GetPrimaryColumns() const
{
	return m_columns;
}

TColumnPtrList ConstraintUnique::GetSecondaryColumns() const
{
	/* 
	 * some code relies on the fact that
	 * an empty list is returned in the case when
	 * there are no secondary columns. We do this even
	 * if in the case of a Unique constraint secondary
	 * columns don't make any sense.
	 */
	return TColumnPtrList();
}

void ConstraintUnique::CodeGenerator(Backend& backend)
{
	backend.HandleUnique(*this);
}

void ConstraintPrimaryKey::CodeGenerator(Backend& backend)
{
	backend.HandlePrimaryKey(*this);
}

void ConstraintForeignKey::AddPrimaryColumn(TColumnPtr column)
{
	m_sonColumns.push_back(column);
}

void ConstraintForeignKey::AddSecondaryColumn(TColumnPtr column)
{
	m_fatherColumns.push_back(column);
}

TColumnPtrList ConstraintForeignKey::GetPrimaryColumns() const
{
	return m_sonColumns;
}

TColumnPtrList ConstraintForeignKey::GetSecondaryColumns() const
{
	return m_fatherColumns;
}

void ConstraintForeignKey::CodeGenerator(Backend& backend)
{
	backend.HandleForeignKey(*this);
}

void ConstraintBase::SetSamePrimaryColumn(TConstraintBasePtrList& list,
	TColumnPtr column)
{
	TConstraintBasePtrList::const_iterator iterator;

	for(iterator = list.begin(); iterator != list.end(); ++iterator)
	{
		(*iterator)->AddPrimaryColumn(column);
	}
}

void ConstraintBase::CodeGenerator(Backend& backend)
{
	std::cerr << "The CodeGenerator method shouldn't be called in a ConstraintBase."
		<< "Aborting!" << std::endl;
	exit(1);
}

