#ifndef _ICM_REGISTRY_H_ 
#define _ICM_REGISTRY_H_ 

#include "ICMKernel/util/icm_macro.h"
#include "expt.h"

#include <string>
#include <map>
#include <vector>

///		
///		A REGISTRY is an instance holder where instances
///		are matched by an external KEY.
///
///		The REGISTRY is owner of its instances : on deletion
///		or copy, the instances are deleted or cloned. 
///
///		It can ADOPT new items (taking ownership)
///		or be completed by Cloning() instance references.
///
///			
template <class UNDERLYING,class KEY=std::string>
class ICM_Registry  
{
public:
	typedef UNDERLYING							_underlying_t ; 
	typedef	KEY									_key_t ;
	typedef	std::map<_key_t,_underlying_t*>		_map_t ;
//
private:
	template<class T> std::string _toString(const T&ref)	const { return ref.toString(); }
	std::string _toString(const std::string& ref)			const { return ref; }
public:
	//!		default constructor. 
	//!		user can define a registry name and a verbose mode (default=false)
	ICM_Registry   (const std::string& name="<noname>",bool verbose=false) : name_(name) , held_(new _map_t) , isVerbose_(verbose) 
	{
		type_= "Registry<"+std::string(typeid(UNDERLYING).name())+">"; 
		ICMLOG( "A Registry "<<name<<" of type "<<type_<<" has been created" ) ;  
		if (!held_) ICMTHROW (ERR_INVALID_ARGUMENT,type_<<"unable to create handle for "<<name_) ; 
	}
	//!		destructor is calling delete for its instances. 
	~ICM_Registry    ()
	{
		ICMLOG( type_<<": Registry "<<name_<<" is being deleted" ) 
		clear(); 
		delete held_; 
		held_=0; 
	}
	//!		empties the registry (deleting instances)
	void clear()
	{
		if (isVerbose_) ICMLOG( type_<<": Registry "<<name_<<" is being cleared" ) ; 
		_map_t::iterator it = held_->begin(); 
		while ( it!=held_->end())
		{
			delete it->second; 
			++it; 
		}
		held_->clear(); 
	}
	//!		sets the specific instance for the key
	//!		the given instance is cloned
	//!		the old instance is deleted . 
	void set(const _key_t& name,const _underlying_t& item) 
	{
		if (isVerbose_) ICMLOG( type_ <<"::add( "<<_toString(name)<<" ) into "<<name_ )   ;
		_map_t::iterator it =  held_->find(name);
		if (it!=held_->end()) delete it->second; 
		(*held_)[name]=item.Clone(); 
	}
	//!		sets the specific instance for the key
	//!		the given instance is adopted
	//!		the old instance is deleted . 
	void adopt(const _key_t& name,_underlying_t* item) 
	{
		if (isVerbose_) ICMLOG( type_ <<"::adopt( "<<_toString(name)<<" ) into "<<name_ )   ;
		if (!item) return; 
		_map_t::iterator it =  held_->find(name);
		if (it!=held_->end()) delete it->second; 
		(*held_)[name]=item; 
	}
	//!		returns the instance or NULL 
	_underlying_t* get(const _key_t& name) const 
	{
		_map_t::iterator it = held_->find(name); 
		if (it==held_->end()) return 0 ;
		return it->second; 
	}
	//!		returns the instance by ref or throws
	_underlying_t& getRef(const _key_t& name) const 
	{
		_map_t::iterator it = held_->find(name); 
		if (it==held_->end()) ICMTHROW(ERR_INVLAID_ARGUMENT,type_<<"::getByRef("<<_toString(name)<<"): nout found") ;
		return *(it->second); 
	}
	//!		true if the instance exists
	bool exists(const _key_t& name) const 
	{
		_map_t::const_iterator it = held_->find(name); 
		return (it!=held_->end()) ;
	}
	//!		remove the specific instance if found
	//!		instance is deleted
	void	remove (const _key_t& name) 
	{
		if (isVerbose_) ICMLOG( type_<<"::remove( "<<_toString(name)<<") from " <<name_  ) ;
		_underlying_t* item = detach(name); 
		if (item) delete item; 
	}
	//!		remove the specific instance if found
	//!		instance is returned to the user (null if it does not exist)
	_underlying_t* detach(const _key_t& name) 
	{
		if (isVerbose_) ICMLOG( type_<<"::detach( "<<_toString(name)<<") from " <<name_  ) ;
		_map_t::iterator it = held_->find(name); 
		if (it==held_->end()) return 0;
		_underlying_t* item = it->second; 
		held_->erase(it);
		return item; 
	}	
	//!		how many keys do we have ?
	int size() const
	{
		return held_->size(); 
	}
	//!		lists all the keys
	std::vector<_key_t> getList() const
	{
		std::vector<_key_t> ret; 
		_map_t::const_iterator it = held_->begin(); 
		while (it!=held_->end())
		{	
			ret.push_back(it->first); 
			++it; 
		}
		return ret; 
	}
	//!		standard copy constructor (deep copy)
	ICM_Registry    (const ICM_Registry&ref ) : type_(ref.type_), name(ref.name_), isVerbose(ref.isVerbose_), held_(new _map_t)
	{
		ICMLOG( "A Registry "<<name<<" of type "<<type_<<" is beeing copyConstructed" ) ;  		
		if (!held_) ICMTHROW (ERR_INVALID_ARGUMENT,type_<<"unable to create handle for "<<name_) ; 
		_map_t::const_iterator it = ref.held_->begin() ; 
		while (it!=ref.held_->end()) 
		{
			set(it->first,*it->second); 
			++it ;
		}
	}
	//!		standard operator= (deep copy)
	ICM_Registry  & operator=(const ICM_Registry&ref)
	{
		if (this!=&ref)
		{
			this->~Registry() ;
			new(this)Registry(ref) ;
		}
		return *this; 
	}
private:
	_map_t * held_;			//	m_held is aggregated, as well as helds datas
	std::string	type_;		//	automatically set at run time
	std::string	name_;		//	user defined 
	bool	isVerbose_;		//	user defined
};

//	-------------------------------------------------------------------------------------------
#endif // _ICM_REGISTRY_H_ 