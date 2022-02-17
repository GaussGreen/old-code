/*************************************************************************** 
*  Name:        ddlmaker/include/column.h                                 
*  Purpose:	Part of the ddlmaker utility                                 
*  Author:      Cosmin Cremarenco                                            
*  Created:     2005-06-07                                                   
*  $Id: column.h,v 1.3 2005/08/24 08:21:52 cosmin Exp $
*  Copyright:   (c) 2005 Trilemma LLP                                        
***************************************************************************/
#ifndef __COLUMN_H__
#define __COLUMN_H__

#include <string>
#include <boost/shared_ptr.hpp>

#include "type.h"
#include "codegeneratorinterface.h"

class Column;
typedef boost::shared_ptr<Column> TColumnPtr;
typedef std::list<TColumnPtr> TColumnPtrList;

class Table;

/*
 * We use this class to describe a table's column.
 * It inherits from Attributes because a column
 * can have its own set of attributes
 */
class Column
:public Attributes,
public CodeGeneratorInterface
{
public:
	Column(std::string name) 
	:m_name(name)
	{
	}

	Column(std::string name, TTypeBasePtr typeBase)
	:m_name(name), m_pTypeBase(typeBase)
	{
	}

	std::string GetName() const
	{
		return m_name;
	}

	TTypeBasePtr GetType() const
	{
		return m_pTypeBase;
	}

	bool operator==(const Column& column) const
	{
		return m_name == column.GetName() &&
			*m_pTypeBase == *(column.GetType());
	}

	static std::string ConcatenateNames(const TColumnPtrList& list,
			const std::string& delimiter);

	boost::shared_ptr<Table> GetOwnerTable() const
	{
		return m_pOwnerTable;
	}

	void SetOwnerTable(boost::shared_ptr<Table> table)
	{
		m_pOwnerTable = table;
	}

	void CodeGenerator(Backend& backend);

	/* default destructor is ok */
private:
	/* 
	 * this is the name of the column as it will appear
	 * in the sql dump
	 */
	std::string m_name;
	
	/*
	 * type of this column
	 */
	TTypeBasePtr m_pTypeBase;

	boost::shared_ptr<Table> m_pOwnerTable;	
};

#endif
