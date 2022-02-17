/*************************************************************************** 
*  Name:        ddlmaker/include/variable.h                                 
*  Purpose:	Part of the ddlmaker utility                                 
*  Author:      Cosmin Cremarenco                                            
*  Created:     2005-06-07                                                   
*  $Id: variable.h,v 1.3 2005/09/08 21:27:43 wang Exp $    
*  Copyright:   (c) 2005 Trilemma LLP                                        
***************************************************************************/
#ifndef __VARIABLE_H__
#define __VARIABLE_H__

#include <string>
#include <cstdlib>
#include <iostream>
#include <map>
#include <list>
#include <vector>

#include <boost/shared_ptr.hpp>

#include "tools.h"

class VarBase;

typedef boost::shared_ptr<VarBase> TVarBasePtr;
typedef std::vector<TVarBasePtr> TVectorVar;
typedef std::list<TVarBasePtr> TVarBasePtrList;

class VarBase
{
public:
	/* use the default constructor to build a pure constant */
	VarBase()
	:m_name(""), m_immutable(true)
	{}
	
	VarBase(std::string name, bool immutable)
	:m_name(name), m_immutable(immutable)
	{}	

	std::string GetName() const
	{
		return m_name;
	}

	void SetName(const std::string& name)
	{
		m_name = name;
	}

	bool IsImmutable() const
	{
		return m_immutable;
	}

	/* you cannot set immutable on a existing variable
	 * you need to set it as immutable from the start.
	 * NO SetImmutable
	 */

	virtual VarBase& operator=(const VarBase& aVarBase)
	{
		m_name = aVarBase.GetName();
		m_immutable = aVarBase.IsImmutable();

		return *this;
	}

	virtual TVarBasePtr Clone();
	virtual TVarBasePtr operator+(VarBase& value);
	virtual TVarBasePtr operator-(VarBase& value);
	virtual TVarBasePtr operator/(VarBase& value);
	virtual TVarBasePtr Mul(VarBase& value);
	virtual TVarBasePtr operator&&(VarBase& value);
	virtual TVarBasePtr operator||(VarBase& value);
	virtual TVarBasePtr operator==(VarBase& value);
	virtual TVarBasePtr operator!=(VarBase& value);
	virtual TVarBasePtr operator>(VarBase& value);
	virtual TVarBasePtr operator<(VarBase& value);
	virtual TVarBasePtr operator>=(VarBase& value);
	virtual TVarBasePtr operator<=(VarBase& value);
	virtual ~VarBase() {};
	
	virtual std::string ToString() { return ""; };
	virtual std::string ToSql() { return ""; };
protected:
	std::string m_name;
	bool m_immutable;
};

template<class T>
class Var;

template<>
class Var<int>;

template<>
class Var<bool>;

template<>
class Var<double>;

template<class T>
bool TryToCast(VarBase& value)
{
	try 
	{
		dynamic_cast<Var<T>& >(value);
		return true;
	} 
	catch(std::bad_cast&) 
	{
		return false;
	}
}

template<class T>
void FailForType(VarBase& value)
{
	try 
	{
		dynamic_cast<Var<T>& >(value);
		
		std::cerr << "FailForType: type " << typeid(T).name() 
			<< " cannot be used in this kind of expression" 
			<< std::endl; 
		exit(1);
	} 
	catch(std::bad_cast&) 
	{
		return;
	}
}

template<>
class Var<bool>
:public VarBase
{
public:
	Var()
	:VarBase("", true)
	{}	

	Var(const std::string& name)
	:VarBase(name, false)
	{}

	Var(const std::string& name, bool value)
	:VarBase(name, true)
	{
		m_value = value;
	}

	Var(bool value)
	:VarBase("", true)
	{
		m_value = value;
	}

	bool& GetValue() const
	{
		return (bool&)m_value;
	}

	void SetValue(bool value)
	{
		m_value = value;
	}

	TVarBasePtr Clone();
	VarBase& operator=(VarBase& aVarBase);
	TVarBasePtr operator&&(VarBase& value);
	TVarBasePtr operator||(VarBase& value);
	TVarBasePtr operator==(VarBase& value);
	TVarBasePtr operator!=(VarBase& value);

	std::string ToString();
	
	std::string ToSql() 
	{ 
		return m_value?"true":"false"; 
	};
private:
	bool m_value;
};

template<>
class Var<std::string>
:public VarBase
{
public:
	Var()
	:VarBase("", true)
	{}	

	Var(const std::string& name)
	:VarBase(name, false)
	{}

	Var(const std::string& name, const std::string& value)
	:VarBase(name, true)
	{
		m_value = value;
	}

	std::string& GetValue() const
	{
		return (std::string&)m_value;
	}

	void SetValue(const std::string& value)
	{
		m_value = value;
	}

	TVarBasePtr Clone();
	VarBase& operator=(VarBase& aVarBase);

	TVarBasePtr operator+(VarBase& value);

	TVarBasePtr operator==(VarBase& value);

	TVarBasePtr operator!=(VarBase& value);

	std::string ToString();

	std::string ToSql();
private:
	std::string m_value;
};

template<>
class Var<int>
:public VarBase
{
public:
	Var()
	:VarBase("", true)
	{}	

	Var(const std::string& name)
	:VarBase(name, false)
	{}

	Var(const std::string& name, int value)
	:VarBase(name, true)
	{
		m_value = value;
	}

	int& GetValue() const
	{
		return (int&)m_value;
	}

	void SetValue(int value)
	{
		m_value = value;
	}

 
	TVarBasePtr Clone();
	VarBase& operator=(VarBase& aVarBase);

	TVarBasePtr operator+(VarBase& value);

	TVarBasePtr operator-(VarBase& value);

	TVarBasePtr operator/(VarBase& value);

	TVarBasePtr Mul(VarBase& value);

	TVarBasePtr operator>(VarBase& value);

	TVarBasePtr operator<(VarBase& value);

	TVarBasePtr operator<=(VarBase& value);

	TVarBasePtr operator>=(VarBase& value);

	TVarBasePtr operator==(VarBase& value);

	TVarBasePtr operator!=(VarBase& value);

	std::string ToString();

	std::string ToSql();
private:
	int m_value;
};

template<>
class Var<double>
:public VarBase
{
public:
	Var()
	:VarBase("", true)
	{}	

	Var(const std::string& name)
	:VarBase(name, false)
	{}

	Var(const std::string& name, double value)
	:VarBase(name, true)
	{
		m_value = value;
	}

	double& GetValue() const
	{
		return (double&)m_value;
	}

	void SetValue(double value)
	{
		m_value = value;
	}

 
	TVarBasePtr Clone();
	VarBase& operator=(VarBase& aVarBase);

	TVarBasePtr operator+(VarBase& value);

	TVarBasePtr operator-(VarBase& value);

	TVarBasePtr operator/(VarBase& value);

	TVarBasePtr Mul(VarBase& value);

	TVarBasePtr operator<(VarBase& value);

	TVarBasePtr operator<=(VarBase& value);

	TVarBasePtr operator>=(VarBase& value);

	TVarBasePtr operator==(VarBase& value);

	TVarBasePtr operator!=(VarBase& value);

	std::string ToString();
	std::string ToSql();
private:
	double m_value;
};

template<>
class Var<TVectorVar>
:public VarBase
{
public:
	Var()
	:VarBase("", true)
	{}	

	Var(std::string aName)
	:VarBase(aName, false)
	{}

	Var(std::string name, TVectorVar value)
	:VarBase(name, true)
	{
		m_value = value;
	}

	TVectorVar& GetValue() const
	{
		return (TVectorVar&)m_value;
	}

	void SetValue(TVectorVar value)
	{
		m_value = value;
	}
 
	TVarBasePtr Clone();
	VarBase& operator=(VarBase& valueBase);

	TVarBasePtr GetValueCopy(int index);
	virtual std::string ToString();

	std::string ToSql();
private:
	TVectorVar m_value;
};

template<class T>
std::ostream& operator<<(std::ostream& os, const Var<T>& value)
{
	os << value.GetValue();

  return os;
}

template<class T>
bool IsTrue(T& value)
{
	bool toReturn = false;
	
	if(TryToCast<bool>(value))
		toReturn = (reinterpret_cast<Var<bool>& >(value)).GetValue();
	else
		Error("Value \"%s\" couldn't be converted to bool", value.ToString().c_str());

	return toReturn;
}

#endif
