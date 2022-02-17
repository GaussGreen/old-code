/*************************************************************************** 
*  Name:        ddlmaker/include/constraint.h                                 
*  Purpose:	Part of the ddlmaker utility                                 
*  Author:      Cosmin Cremarenco                                            
*  Created:     2005-06-07                                                   
*  $Id: constraint.h,v 1.4 2005/09/19 16:18:43 cosmin Exp $
*  Copyright:   (c) 2005 Trilemma LLP                                        
***************************************************************************/
#ifndef __CCONSTRAINT_HPP__
#define __CCONSTRAINT_HPP__

#include <string>
#include <boost/shared_ptr.hpp>

#include "type.h"
#include "column.h"
#include "codegeneratorinterface.h"
#include "attribute.h"

class ConstraintBase;

typedef boost::shared_ptr<ConstraintBase> TConstraintBasePtr;
typedef std::list<TConstraintBasePtr> TConstraintBasePtrList;

class Table;

class ConstraintBase
:public CodeGeneratorInterface,
public Attributes
{
public:
	ConstraintBase(const std::string& name) 
	:m_name(name) 
	{
	}

	/*
	 * Every constraint is linked to at least one
	 * column. This sets it.
	 * However, for example, for foreign keys we can
	 * have more than one son column.
	 */
	virtual void AddPrimaryColumn(TColumnPtr column);

	/*
	 * Some constraints can have destination columns.
	 * This sets them.
	 * Take for example the foreign key which can
	 * point to multiple father columns.
	 */
	virtual void AddSecondaryColumn(TColumnPtr column);

	virtual TColumnPtrList GetPrimaryColumns() const;

	virtual TColumnPtrList GetSecondaryColumns() const;
		
	/* 
	 * This method will check for some coherency requirements
	 * for every constraint. For now it will only issue
	 * warnings.
	 */
	virtual void CheckCoherency() const
	{
	}
	
	static void SetSamePrimaryColumn(TConstraintBasePtrList& list,
			TColumnPtr column);
	
	virtual ~ConstraintBase()
	{
	}

	std::string GetName() const
	{
		return m_name;
	}

	void SetName(const std::string& name)
	{
		m_name = name;
	}

	boost::shared_ptr<Table> GetOwnerTable() const
	{
		return m_pOwnerTable;
	}

	void SetOwnerTable(boost::shared_ptr<Table> table)
	{
		m_pOwnerTable = table;
	}

	virtual void CodeGenerator(Backend& backend);

protected:
	std::string m_name;

	boost::shared_ptr<Table> m_pOwnerTable;
};

class ConstraintNull
:public ConstraintBase
{
public:
	ConstraintNull()
	: ConstraintBase(""),
		m_notNull(false)
	{
	}
	
	void AddPrimaryColumn(TColumnPtr column);

	void AddSecondaryColumn(TColumnPtr column);
	
	TColumnPtrList GetPrimaryColumns() const;

	TColumnPtrList GetSecondaryColumns() const;

	void CheckCoherency() const;

	void SetNotNull(bool notNull)
	{
		m_notNull = notNull;
	}

	bool NotNull() const
	{
		return m_notNull;
	}

	void CodeGenerator(Backend& backend);
	
private:
	bool m_notNull;

	/* make it a list even if it will contain a single
	 * element every time
	 */
	TColumnPtrList m_columns;
};

class ConstraintUnique
:public ConstraintBase
{
public:
	ConstraintUnique()
	: ConstraintBase("")
	{
	}

	ConstraintUnique(const std::string& name)
	: ConstraintBase(name)
	{
	}
	
	void AddPrimaryColumnFromList(const TColumnPtrList& list);
	
	virtual void AddPrimaryColumn(TColumnPtr column);

	virtual void AddSecondaryColumn(TColumnPtr column);
	
	virtual TColumnPtrList GetPrimaryColumns() const;

	virtual TColumnPtrList GetSecondaryColumns() const;
	
	virtual ~ConstraintUnique()
	{
	}

	virtual void CodeGenerator(Backend& backend);

protected:
	TColumnPtrList m_columns;	
};

typedef boost::shared_ptr<ConstraintUnique> TConstraintUniquePtr;

class ConstraintPrimaryKey
:public ConstraintUnique
{
public:
	ConstraintPrimaryKey()
	: ConstraintUnique("")
	{
	}

	ConstraintPrimaryKey(const std::string& name)
	: ConstraintUnique(name)
	{
	}

	virtual ~ConstraintPrimaryKey()
	{
	}

	void CodeGenerator(Backend& backend);
};

typedef boost::shared_ptr<ConstraintPrimaryKey> 
	TConstraintPrimaryKeyPtr;

class ConstraintForeignKey
:public ConstraintBase
{
public:
	ConstraintForeignKey()
	: ConstraintBase("")
	{
	}

	ConstraintForeignKey(const std::string& name, 
			boost::shared_ptr<Table> fatherTable)
	: ConstraintBase(name)
	{
		m_fatherTable = fatherTable;
	}

	void AddPrimaryColumn(TColumnPtr column);

	void AddSecondaryColumn(TColumnPtr column);
	
	TColumnPtrList GetPrimaryColumns() const;

	TColumnPtrList GetSecondaryColumns() const;
	
	virtual ~ConstraintForeignKey()
	{
	}

	void SetFatherTable(boost::shared_ptr<Table> fatherTable)
	{
		m_fatherTable = fatherTable;
	}

	boost::shared_ptr<Table> GetFatherTable() const
	{
		return m_fatherTable;
	}

	void CodeGenerator(Backend& backend);

protected:
	boost::shared_ptr<Table> m_fatherTable;
	
	TColumnPtrList m_sonColumns;	

	TColumnPtrList m_fatherColumns;	
};

typedef boost::shared_ptr<ConstraintForeignKey> 
	TConstraintForeignKeyPtr;

#endif
