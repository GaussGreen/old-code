#ifndef __TABLE_H__
#define __TABLE_H__

#include <string>
#include <stack>
#include <boost/shared_ptr.hpp>

#include "column.h"
#include "constraint.h"
#include "attribute.h"
#include "codegeneratorinterface.h"

class Table;

typedef boost::shared_ptr<Table> TTablePtr;
typedef boost::shared_ptr<std::list<Table> > TTableListPtr;
typedef std::list<TTablePtr> TTablePtrList;
typedef std::stack<TTablePtr> TTablePtrStack;

class Table
:public Attributes,
public CodeGeneratorInterface
{
public:
	Table(std::string name) 
	:Attributes(), m_name(name)
	{}

	void AddColumn(TColumnPtr column)
	{
		m_columns.push_back(column);
	}

	void AddConstraint(TConstraintBasePtr constraint)
	{
		m_constraints.push_back(constraint);
	}

	void AddConstraints(TConstraintBasePtrList constraints)
	{
		std::copy(constraints.begin(), constraints.end(),
				std::back_inserter(m_constraints));
	}

	TColumnPtr GetColumnByName(const std::string& name) const;
	
	TConstraintBasePtr GetConstraintByName(const std::string& name) const;

	static TConstraintBasePtr GetConstraintByName(const TTablePtrList& list,
			const std::string& name);
	
	static TColumnPtr GetColumnByName(const TTablePtrList& list,
			const std::string& name);

	static TTablePtr GetTableByColumnName(const TTablePtrList& list,
			const std::string& name);
	
	static TTablePtr GetTableByName(const TTablePtrList& list,
			const std::string& name);

	std::string GetName() const
	{
		return m_name;
	}

	void CodeGenerator(Backend& backend);
	
	bool AliasedTypeIsUsed(const std::string& alias) const;

	static bool AliasedTypeIsUsed(const TTablePtrList& list,
			const std::string& alias);
	
	/*
	 * This is just a simple function to determine if there are any
	 * cycles between tables introduced by foreign keys.
	 * This is a recursive function, at first call, "vertex"
	 * will be the root of the graph.
	 * 
	 * "domain" is just a set of tables to remove the tables
	 * that we already checked (keep in mind that
	 * a database contains multiple tables, thus can
	 * contain multiple graphs, there's no point in 
	 * checking for cycles for every table if we can
	 * eliminate some).
	 * 
	 * "context" is a depth stack to see where we are in the graph.
	 * Returns true if no cycle was detected, false if cycles
	 * were detected (outputs the cycles to std::err)
	 */
	static bool CheckCycle(TTablePtr vertex, TTablePtrList& domain,
			TTablePtrList& context);
			
	TColumnPtrList& GetColumns()
	{
		return m_columns;
	}

	TConstraintBasePtrList& GetConstraints()
	{
		return m_constraints;
	}
	
	/* This method makes some checks on the various
	 * components of the table and delivers a warning message
	 * if things do not fit (in the future might issue
	 * an exception).
	 */
	void CheckCoherency() const;
	
	class SameTableName
	{
	public:
		SameTableName(const std::string& name)
		:m_name(name)
		{
		}
		bool operator() (TTablePtr table)
		{
			return table->GetName() == m_name;
		}
	private:
		std::string m_name;
	};
	
private:
	std::string m_name;
	TColumnPtrList m_columns;
	TConstraintBasePtrList m_constraints;
};

#endif
