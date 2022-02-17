/*************************************************************************** 
*  Name:        ddlmaker/src/variable.cpp                              
*  Purpose:	Part of the ddlmaker utility                                 
*  Author:      Cosmin Cremarenco                                            
*  Created:     2005-06-07                                                   
*  RCS-ID:      $Id: variable.cpp,v 1.2 2005/07/06 09:03:02 cosmin Exp $    
*  Copyright:   (c) 2005 Trilemma LLP                                        
***************************************************************************/
#include <sstream>

#include "variable.h"

TVarBasePtr VarBase::operator+( VarBase& aValue)
{
	std::cerr << "operator+ called in VarBase - Type semantic error - Aborting" << std::endl;
	exit(2);
}

TVarBasePtr VarBase::operator-( VarBase& aValue)
{
	std::cerr << "operator- called in VarBase - Type semantic error - Aborting" << std::endl;
	exit(2);
}

TVarBasePtr VarBase::operator/( VarBase& aValue)
{
	std::cerr << "operator/ called in VarBase - Type semantic error - Aborting" << std::endl;
	exit(2);
}

TVarBasePtr VarBase::Mul( VarBase& aValue)
{
	std::cerr << "operator* called in VarBase - Type semantic error - Aborting" << std::endl;
	exit(2);
}

TVarBasePtr VarBase::operator&&( VarBase& aValue)
{
	std::cerr << "operator && called in VarBase - Type semantic error - Aborting" << std::endl;
	exit(2);
}

TVarBasePtr VarBase::operator||( VarBase& aValue)
{
	std::cerr << "operator && called in VarBase - Type semantic error - Aborting" << std::endl;
	exit(2);
}

TVarBasePtr VarBase::operator==( VarBase& aValue)
{
	std::cerr << "operator == called in VarBase - Type semantic error - Aborting" << std::endl;
	exit(2);
}

TVarBasePtr VarBase::operator!=( VarBase& aValue)
{
	std::cerr << "operator != called in VarBase - Type semantic error - Aborting" << std::endl;
	exit(2);
}

TVarBasePtr VarBase::operator<( VarBase& aValue)
{
	std::cerr << "operator < called in VarBase - Type semantic error - Aborting" << std::endl;
	exit(2);
}

TVarBasePtr VarBase::operator>( VarBase& aValue)
{
	std::cerr << "operator > called in VarBase - Type semantic error - Aborting" << std::endl;
	exit(2);
}

TVarBasePtr VarBase::operator<=( VarBase& aValue)
{
	std::cerr << "operator <= called in VarBase - Type semantic error - Aborting" << std::endl;
	exit(2);
}

TVarBasePtr VarBase::operator>=( VarBase& aValue)
{
	std::cerr << "operator >= called in VarBase - Type semantic error - Aborting" << std::endl;
	exit(2);
}

TVarBasePtr VarBase::Clone()
{
	VarBase* temp = new VarBase();
	(*temp) = *this;
	
	return TVarBasePtr(temp);
}

VarBase& Var<bool>::operator=( VarBase& aVarBase)
{
	if(this->IsImmutable())
	{
		std::cerr << "cannot assign to an immutable" << std::endl;
		exit(2);
	}
	else 
	{
		Var<bool>& aTemp = reinterpret_cast<Var<bool>& >(aVarBase);
		VarBase::SetName(aTemp.GetName());
		SetValue(aTemp.GetValue());	
	}
	return *this;
}

TVarBasePtr Var<bool>::Clone()
{
	Var<bool>* temp = new Var<bool>();
	(VarBase)(*temp) = *this;
	temp->SetValue(GetValue());	
	
	return TVarBasePtr(temp);
}

TVarBasePtr Var<bool>::operator&&( VarBase& aValue)
{
	if(! TryToCast<bool>(*this) || 
		! TryToCast<bool>(aValue)) 
	{
		std::cerr << "&&: failed casting to bool" << std::endl;
		exit(2);
	}
	Var<bool>& aTemp = reinterpret_cast<Var<bool>& >(aValue);
	TVarBasePtr aReturn(new Var<bool>(std::string(""), GetValue() && aTemp.GetValue()));	
	return aReturn;
}

TVarBasePtr Var<bool>::operator||( VarBase& aValue)
{
	if(! TryToCast<bool>(*this) || 
		! TryToCast<bool>(aValue)) 
	{
		std::cerr << "||: failed casting to bool" << std::endl;
		exit(2);
	}
	Var<bool>& aTemp = reinterpret_cast<Var<bool>& >(aValue);
	TVarBasePtr aReturn(new Var<bool>(std::string(""), GetValue() || aTemp.GetValue()));	
	return aReturn;
}

TVarBasePtr Var<bool>::operator==( VarBase& aValue)
{
	if(! TryToCast<bool>(aValue)) 
	{
		std::cerr << "==: failed casting to bool" << std::endl;
		exit(2);
	}
	Var<bool>& aTemp = reinterpret_cast<Var<bool>& >(aValue);
	TVarBasePtr aReturn(new Var<bool>(std::string(""), GetValue() == aTemp.GetValue()));	
	std::cerr << "comparison result " << aReturn->ToString() << std::endl;
	return aReturn;
}

TVarBasePtr Var<bool>::operator!=( VarBase& aValue)
{
	if(! TryToCast<bool>(aValue)) 
	{
		std::cerr << "!=: failed casting to bool" << std::endl;
		exit(2);
	}
	Var<bool>& aTemp = reinterpret_cast<Var<bool>& >(aValue);
	TVarBasePtr aReturn(new Var<bool>(std::string(""), GetValue() != aTemp.GetValue()));	
	return aReturn;
}

std::string Var<bool>::ToString()
{
	std::ostringstream os;
	os << m_value;
	return os.str();
}

VarBase& Var<int>::operator=( VarBase& aVarBase)
{
	if(this->IsImmutable())
	{
		std::cerr << "cannot assign to an immutable" << std::endl;
		exit(2);
	}
	else 
	{
		Var<int>& aTemp = reinterpret_cast<Var<int>& >(aVarBase);
		VarBase::SetName(aTemp.GetName());
		SetValue(aTemp.GetValue());	
	}
	return *this;
}

TVarBasePtr Var<int>::Clone()
{
	Var<int>* temp = new Var<int>();
	(VarBase)(*temp) = *this;
	temp->SetValue(GetValue());	
	
	return TVarBasePtr(temp);
}

TVarBasePtr Var<int>::operator+( VarBase& aValue)
{
	if(TryToCast<int>(aValue)) 
	{
		Var<int>& aTemp2 = reinterpret_cast<Var<int>& >(aValue);
		return boost::shared_ptr<VarBase>(new Var<int> ("", this->GetValue() + aTemp2.GetValue()));
	} 
	else if(TryToCast<double>(aValue)) 
	{
		Var<double>& aTemp2 = reinterpret_cast<Var<double>& >(aValue);
		return boost::shared_ptr<VarBase>(new Var<double> ("", this->GetValue() + aTemp2.GetValue()));
	} 
	else 
	{
		std::cerr << "+: couldn't cast" << std::endl;
		exit(2);
	}
}

TVarBasePtr Var<int>::operator-( VarBase& aValue)
{
	if(TryToCast<int>(aValue)) 
	{
		Var<int>& aTemp2 = reinterpret_cast<Var<int>& >(aValue);
		return boost::shared_ptr<VarBase>(new Var<int> ("", this->GetValue() - aTemp2.GetValue()));
	} 
	else if(TryToCast<double>(aValue)) 
	{
		Var<double>& aTemp2 = reinterpret_cast<Var<double>& >(aValue);
		return boost::shared_ptr<VarBase>(new Var<double> ("", this->GetValue() - aTemp2.GetValue()));
	} 
	else 
	{
		std::cerr << "-: couldn't cast" << std::endl;
		exit(2);
	}
}

TVarBasePtr Var<int>::operator/( VarBase& aValue)
{
	if(TryToCast<int>(aValue)) 
	{
		Var<int>& aTemp2 = reinterpret_cast<Var<int>& >(aValue);
		return boost::shared_ptr<VarBase>(new Var<int> ("", this->GetValue() / aTemp2.GetValue()));
	} 
	else if(TryToCast<double>(aValue)) 
	{
		Var<double>& aTemp2 = reinterpret_cast<Var<double>& >(aValue);
		return boost::shared_ptr<VarBase>(new Var<double> ("", this->GetValue() / aTemp2.GetValue()));
	} 
	else 
	{
		std::cerr << "/: couldn't cast" << std::endl;
		exit(2);
	}
}

TVarBasePtr Var<int>::Mul( VarBase& aValue)
{
	if(TryToCast<int>(aValue)) 
	{
		Var<int>& aTemp2 = reinterpret_cast<Var<int>& >(aValue);
		return boost::shared_ptr<VarBase>(new Var<int> ("", this->GetValue() * aTemp2.GetValue()));
	} 
	else if(TryToCast<double>(aValue)) 
	{
		Var<double>& aTemp2 = reinterpret_cast<Var<double>& >(aValue);
		return boost::shared_ptr<VarBase>(new Var<double> ("", this->GetValue() * aTemp2.GetValue()));
	} 
	else 
	{
		std::cerr << "*: couldn't cast" << std::endl;
		exit(2);
	}
}

TVarBasePtr Var<int>::operator>( VarBase& aValue)
{
	if(TryToCast<int>(aValue)) 
	{
		Var<int>& aTemp2 = reinterpret_cast<Var<int>& >(aValue);
		return boost::shared_ptr<VarBase>(new Var<bool> ("", this->GetValue() > aTemp2.GetValue()));
	} 
	else if(TryToCast<double>(aValue)) 
	{
		Var<double>& aTemp2 = reinterpret_cast<Var<double>& >(aValue);
		return boost::shared_ptr<VarBase>(new Var<bool> ("", this->GetValue() > aTemp2.GetValue()));
	} 
	else 
	{
		std::cerr << ">: couldn't cast" << std::endl;
		exit(2);
	}
}

TVarBasePtr Var<int>::operator<( VarBase& aValue)
{
	if(TryToCast<int>(aValue)) 
	{
		Var<int>& aTemp2 = reinterpret_cast<Var<int>& >(aValue);
		return boost::shared_ptr<VarBase>(new Var<bool> ("", this->GetValue() < aTemp2.GetValue()));
	} 
	else if(TryToCast<double>(aValue)) 
	{
		Var<double>& aTemp2 = reinterpret_cast<Var<double>& >(aValue);
		return boost::shared_ptr<VarBase>(new Var<bool> ("", this->GetValue() < aTemp2.GetValue()));
	} 
	else 
	{
		std::cerr << "<: couldn't cast" << std::endl;
		exit(2);
	}
}

TVarBasePtr Var<int>::operator<=( VarBase& aValue)
{
	if(TryToCast<int>(aValue)) 
	{
		Var<int>& aTemp2 = reinterpret_cast<Var<int>& >(aValue);
		return boost::shared_ptr<VarBase>(new Var<bool> ("", this->GetValue() <= aTemp2.GetValue()));
	} 
	else if(TryToCast<double>(aValue)) 
	{
		Var<double>& aTemp2 = reinterpret_cast<Var<double>& >(aValue);
		return boost::shared_ptr<VarBase>(new Var<bool> ("", this->GetValue() <= aTemp2.GetValue()));
	} 
	else 
	{
		std::cerr << "<=: couldn't cast" << std::endl;
		exit(2);
	}
}

TVarBasePtr Var<int>::operator>=( VarBase& aValue)
{
	if(TryToCast<int>(aValue)) 
	{
		Var<int>& aTemp2 = reinterpret_cast<Var<int>& >(aValue);
		return boost::shared_ptr<VarBase>(new Var<bool> ("", this->GetValue() >= aTemp2.GetValue()));
	} 
	else if(TryToCast<double>(aValue)) 
	{
		Var<double>& aTemp2 = reinterpret_cast<Var<double>& >(aValue);
		return boost::shared_ptr<VarBase>(new Var<bool> ("", this->GetValue() >= aTemp2.GetValue()));
	} 
	else 
	{
		std::cerr << ">=: couldn't cast" << std::endl;
		exit(2);
	}
}

TVarBasePtr Var<int>::operator==( VarBase& aValue)
{
	if(TryToCast<int>(aValue)) 
	{
		Var<int>& aTemp2 = reinterpret_cast<Var<int>& >(aValue);
		return boost::shared_ptr<VarBase>(new Var<bool> ("", this->GetValue() == aTemp2.GetValue()));
	} 
	else if(TryToCast<double>(aValue)) 
	{
		Var<double>& aTemp2 = reinterpret_cast<Var<double>& >(aValue);
		return boost::shared_ptr<VarBase>(new Var<bool> ("", (this->GetValue() - aTemp2.GetValue()) < 0.00001));
	} 
	else 
	{
		std::cerr << "==: couldn't cast" << std::endl;
		exit(2);
	}
}

TVarBasePtr Var<int>::operator!=( VarBase& aValue)
{
	if(TryToCast<int>(aValue)) 
	{
		Var<int>& aTemp2 = reinterpret_cast<Var<int>& >(aValue);
		return boost::shared_ptr<VarBase>(new Var<bool> ("", this->GetValue() != aTemp2.GetValue()));
	} 
	else if(TryToCast<double>(aValue)) 
	{
		Var<double>& aTemp2 = reinterpret_cast<Var<double>& >(aValue);
		return boost::shared_ptr<VarBase>(new Var<bool> ("", (this->GetValue() - aTemp2.GetValue()) > 0.00001));
	} 
	else 
	{
		std::cerr << "!=: couldn't cast" << std::endl;
		exit(2);
	}
}

std::string Var<int>::ToString()
{
	std::ostringstream os;
	os << m_value;
	return os.str();
}

std::string Var<int>::ToSql()
{
	return ToString();
}

VarBase& Var<std::string>::operator=( VarBase& aVarBase)
{
	if(this->IsImmutable())
	{
		std::cerr << "cannot assign to an immutable" << std::endl;
		exit(2);
	}
	else 
	{
		Var<std::string>& aTemp = reinterpret_cast<Var<std::string>& >(aVarBase);
		VarBase::SetName(aTemp.GetName());
		SetValue(aTemp.GetValue());	
	}
	return *this;
}

TVarBasePtr Var<std::string>::Clone()
{
	Var<std::string>* temp = new Var<std::string>();
	(VarBase)(*temp) = *this;
	temp->SetValue(GetValue());
	
	return TVarBasePtr(temp);
}

TVarBasePtr Var<std::string>::operator+( VarBase& aValue)
{
	std::ostringstream oss;

	oss << m_value;

	if(TryToCast<bool>(aValue)) 
	{
		Var<bool>& aTemp = reinterpret_cast<Var<bool>& >(aValue);
		oss << aTemp.GetValue();
	} 
	else if(TryToCast<int>(aValue)) 
	{
		Var<int>& aTemp = reinterpret_cast<Var<int>& >(aValue);
		oss << aTemp.GetValue();
	} 
	else if(TryToCast<double>(aValue)) 
	{
		Var<double> aTemp = reinterpret_cast<Var<double>&>(aValue);
		oss << aTemp.GetValue();
	} 
	else if(TryToCast<std::string>(aValue)) 
	{
		Var<std::string>& aTemp = reinterpret_cast<Var<std::string>&>(aValue);
		oss << aTemp.GetValue();
	} 
	else 
	{
		std::cerr << "+: failed casting to string" << std::endl;
		exit(2);
	}
	
	return boost::shared_ptr<VarBase>(new Var<std::string>(std::string(""), oss.str()));
}

TVarBasePtr Var<std::string>::operator==( VarBase& aValue)
{
	if( ! TryToCast<std::string>(aValue)) 
	{
		std::cerr << "==: failed casting to string" << std::endl;
		exit(2);
	}
	Var<std::string>& aTemp = reinterpret_cast<Var<std::string>& >(aValue);
	TVarBasePtr aReturn(new Var<bool>(std::string(""), GetValue().compare(aTemp.GetValue()) == 0));
	
	return aReturn;
}

TVarBasePtr Var<std::string>::operator!=( VarBase& aValue)
{
	if(! TryToCast<std::string>(aValue)) 
	{
		std::cerr << "!=: failed casting to string" << std::endl;
		exit(2);
	}
	Var<std::string>& aTemp = reinterpret_cast<Var<std::string>& >(aValue);
	TVarBasePtr aReturn(new Var<bool>(std::string(""), GetValue().compare(aTemp.GetValue()) != 0));	
	
	return aReturn;
}

std::string Var<std::string>::ToString()
{
	return m_value;
}

std::string Var<std::string>::ToSql()
{
	std::string buffer;
	std::string::const_iterator pos;

	buffer.append("'");
	for(pos = m_value.begin(); pos != m_value.end(); ++pos)
	{
		if(*pos == '\'')
			buffer += "'";
		buffer += *pos;
	}
	buffer.append("'");

	return buffer;
}

VarBase& Var<double>::operator=( VarBase& aVarBase)
{
	if(this->IsImmutable())
	{
		std::cerr << "cannot assign to an immutable" << std::endl;
		exit(2);
	}
	else 
	{
		Var<double>& aTemp = reinterpret_cast<Var<double>& >(aVarBase);
		VarBase::SetName(aTemp.GetName());
		SetValue(aTemp.GetValue());	
	}
	return *this;
}

TVarBasePtr Var<double>::Clone()
{
	Var<double>* temp = new Var<double>();
	(VarBase)(*temp) = *this;
	temp->SetValue(GetValue());
	
	return TVarBasePtr(temp);
}
	
TVarBasePtr Var<double>::operator+( VarBase& aValue)
{
	if(TryToCast<int>(aValue)) 
	{
		Var<int>& aTemp2 = reinterpret_cast<Var<int>& >(aValue);
		return boost::shared_ptr<VarBase>(new Var<double> ("", this->GetValue() + aTemp2.GetValue()));
	} 
	else if(TryToCast<double>(aValue)) 
	{
		Var<double>& aTemp2 = reinterpret_cast<Var<double>& >(aValue);
		return boost::shared_ptr<VarBase>(new Var<double> ("", this->GetValue() + aTemp2.GetValue()));
	} 
	else 
	{
		std::cerr << "+: couldn't cast" << std::endl;
		exit(2);
	}
}

TVarBasePtr Var<double>::operator-( VarBase& aValue)
{
	if(TryToCast<int>(aValue)) 
	{
		Var<int>& aTemp2 = reinterpret_cast<Var<int>& >(aValue);
		return boost::shared_ptr<VarBase>(new Var<double> ("", this->GetValue() - aTemp2.GetValue()));
	} 
	else if(TryToCast<double>(aValue)) 
	{
		Var<double>& aTemp2 = reinterpret_cast<Var<double>& >(aValue);
		return boost::shared_ptr<VarBase>(new Var<double> ("", this->GetValue() - aTemp2.GetValue()));
	} 
	else 
	{
		std::cerr << "-: couldn't cast" << std::endl;
		exit(2);
	}
}

TVarBasePtr Var<double>::operator/( VarBase& aValue)
{
	if(TryToCast<int>(aValue)) 
	{
		Var<int>& aTemp2 = reinterpret_cast<Var<int>& >(aValue);
		return boost::shared_ptr<VarBase>(new Var<double> ("", this->GetValue() / aTemp2.GetValue()));
	} 
	else if(TryToCast<double>(aValue)) 
	{
		Var<double>& aTemp2 = reinterpret_cast<Var<double>& >(aValue);
		return boost::shared_ptr<VarBase>(new Var<double> ("", this->GetValue() / aTemp2.GetValue()));
	} 
	else 
	{
		std::cerr << "/: couldn't cast" << std::endl;
		exit(2);
	}
}

TVarBasePtr Var<double>::Mul( VarBase& aValue)
{
	if(TryToCast<int>(aValue)) 
	{
		Var<int>& aTemp2 = reinterpret_cast<Var<int>& >(aValue);
		return boost::shared_ptr<VarBase>(new Var<double> ("", this->GetValue() * aTemp2.GetValue()));
	} 
	else if(TryToCast<double>(aValue)) 
	{
		Var<double>& aTemp2 = reinterpret_cast<Var<double>& >(aValue);
		return boost::shared_ptr<VarBase>(new Var<double> ("", this->GetValue() * aTemp2.GetValue()));
	} 
	else 
	{
		std::cerr << "mul: couldn't cast" << std::endl;
		exit(2);
	}
}

TVarBasePtr Var<double>::operator<( VarBase& aValue)
{
	if(TryToCast<int>(aValue)) 
	{
		Var<int>& aTemp2 = reinterpret_cast<Var<int>& >(aValue);
		return boost::shared_ptr<VarBase>(new Var<bool> ("", (this->GetValue() - aTemp2.GetValue()) < 0.00001));
	} 
	else if(TryToCast<double>(aValue)) 
	{
		Var<double>& aTemp2 = reinterpret_cast<Var<double>& >(aValue);
		return boost::shared_ptr<VarBase>(new Var<bool> ("", (this->GetValue() - aTemp2.GetValue()) < 0.00001));
	} 
	else 
	{
		std::cerr << "<: couldn't cast" << std::endl;
		exit(2);
	}
}

TVarBasePtr Var<double>::operator<=( VarBase& aValue)
{
	if(TryToCast<int>(aValue)) 
	{
		Var<int>& aTemp2 = reinterpret_cast<Var<int>& >(aValue);
		return boost::shared_ptr<VarBase>(new Var<bool> ("", (this->GetValue() - aTemp2.GetValue()) <= 0.00001));
	} 
	else if(TryToCast<double>(aValue)) 
	{
		Var<double>& aTemp2 = reinterpret_cast<Var<double>& >(aValue);
		return boost::shared_ptr<VarBase>(new Var<bool> ("", (this->GetValue() - aTemp2.GetValue()) <= 0.00001));
	} 
	else 
	{
		std::cerr << "<=: couldn't cast" << std::endl;
		exit(2);
	}
}

TVarBasePtr Var<double>::operator>=( VarBase& aValue)
{	
	if(TryToCast<int>(aValue)) 
	{
		Var<int>& aTemp2 = reinterpret_cast<Var<int>& >(aValue);
		return boost::shared_ptr<VarBase>(new Var<bool> ("", (this->GetValue() - aTemp2.GetValue()) >= 0.00001));
	} 
	else if(TryToCast<double>(aValue)) 
	{
		Var<double>& aTemp2 = reinterpret_cast<Var<double>& >(aValue);
		return boost::shared_ptr<VarBase>(new Var<bool> ("", (this->GetValue() - aTemp2.GetValue()) >= 0.00001));
	} 
	else 
	{
		std::cerr << ">=: couldn't cast" << std::endl;
		exit(2);
	}
}

TVarBasePtr Var<double>::operator==( VarBase& aValue)
{
	if(TryToCast<int>(aValue)) 
	{
		Var<int>& aTemp2 = reinterpret_cast<Var<int>& >(aValue);
		return boost::shared_ptr<VarBase>(new Var<bool> ("", (this->GetValue() - aTemp2.GetValue()) < 0.00001));
	} 
	else if(TryToCast<double>(aValue)) 
	{
		Var<double>& aTemp2 = reinterpret_cast<Var<double>& >(aValue);
		return boost::shared_ptr<VarBase>(new Var<bool> ("", (this->GetValue() - aTemp2.GetValue()) < 0.00001));
	} 
	else 
	{
		std::cerr << "==: couldn't cast" << std::endl;
		exit(2);
	}
}

TVarBasePtr Var<double>::operator!=(VarBase& aValue)
{
	if(TryToCast<int>(aValue)) 
	{
		Var<int>& aTemp2 = reinterpret_cast<Var<int>& >(aValue);
		return boost::shared_ptr<VarBase>(new Var<bool> ("", (this->GetValue() - aTemp2.GetValue()) > 0.00001));
	} 
	else if(TryToCast<double>(aValue)) 
	{
		Var<double>& aTemp2 = reinterpret_cast<Var<double>& >(aValue);
		return boost::shared_ptr<VarBase>(new Var<bool> ("", (this->GetValue() - aTemp2.GetValue()) > 0.00001));
	} 
	else 
	{
		std::cerr << "!=: couldn't cast" << std::endl;
		exit(2);
	}
}
 
VarBase& Var<TVectorVar>::operator=(VarBase& valueBase)
{
	Var<TVectorVar>& temp = reinterpret_cast<Var<TVectorVar>& >(valueBase);
	VarBase::SetName(temp.GetName());
	SetValue(temp.GetValue());	

	return *this;
}

TVarBasePtr Var<TVectorVar>::Clone()
{
	Var<TVectorVar>* temp = new Var<TVectorVar>();
	(VarBase)(*temp) = *this;
	temp->m_value.clear();
	
	TVectorVar::const_iterator iterator;

	for(iterator = m_value.begin(); iterator != m_value.end(); ++iterator)
		temp->m_value.push_back((*iterator)->Clone());
		
	return TVarBasePtr(temp);
}

TVarBasePtr Var<TVectorVar>::GetValueCopy(int index)
{
	TVarBasePtr initVar = m_value.at(index);

	if(TryToCast<bool>(*initVar)) 
	{
		VarBase* temp = new Var<bool>("", (boost::dynamic_pointer_cast<Var<bool>, 
					VarBase>(initVar))->GetValue());
		return boost::shared_ptr<VarBase>(temp);
	} 
	else if(TryToCast<int>(*initVar)) 
	{
		VarBase* temp = new Var<int>("", (boost::dynamic_pointer_cast<Var<int>, 
					VarBase>(initVar))->GetValue());
		return boost::shared_ptr<VarBase>(temp);
	} 
	else if(TryToCast<double>(*initVar)) 
	{
		VarBase* temp = new Var<double>("", (boost::dynamic_pointer_cast<Var<double>, 
					VarBase>(initVar))->GetValue());
		return boost::shared_ptr<VarBase>(temp);
	} 
	else if(TryToCast<std::string>(*initVar)) 
	{
		VarBase* temp = new Var<std::string>("", (boost::dynamic_pointer_cast<Var<std::string>, 
					VarBase>(initVar))->GetValue());
		return boost::shared_ptr<VarBase>(temp);
	} 
	else if(TryToCast<TVectorVar>(*initVar)) 
	{
		VarBase* temp = new Var<TVectorVar>("", (boost::dynamic_pointer_cast<Var<TVectorVar>, 
					VarBase>(initVar))->GetValue());
		return boost::shared_ptr<VarBase>(temp);
	} 
	else 
	{
		std::cerr << "couldn't cast type that exists in vector" << std::endl;
		exit(1);
	}
}

std::string Var<TVectorVar>::ToString()
{
	TVectorVar::const_iterator iterator;
	std::string result;

	for(iterator = m_value.begin(); 
		iterator != m_value.end(); 
		++iterator)
	{
		result += (*iterator)->ToString();
		result += " ";
	}
	return result;
}

std::string Var<TVectorVar>::ToSql()
{
	std::string result;

	for(unsigned int i = 0; i < m_value.size(); ++i) 
	{
		result += m_value[i]->ToSql();
		
		if(i < m_value.size() - 1)
			result += ", ";
	}
	return result;
}

std::string Var<double>::ToString()
{
	std::ostringstream os;
	os << m_value;
	return os.str();
}

std::string Var<double>::ToSql()
{
	return ToString();
}
