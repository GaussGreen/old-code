
#ifndef _ICM_TREEBASE_H_
#define _ICM_TREEBASE_H_

#include "ICMKernel/util/ICM_QMatrix.h"
#include <algorithm>

#include <map>
//	-----------------------------------------------------------------------
//
//
//	-----------------------------------------------------------------------
template <unsigned int NCHILD,class NODEDATA>
class NodeElt 
{
public:
	typedef	NODEDATA	node_data_t; 
	typedef NodeElt		node_elt_t; 
	typedef std::pair<double,NodeElt*>	child_t; 
	typedef child_t*		child_iterator ;
	typedef const child_t *	const_child_iterator ;
	//
	typedef	node_elt_t*			iterator;
	typedef	const node_elt_t*	const_iterator;
public:
	NodeElt()
	{
		for(unsigned int i=0;i<NCHILD;i++) itsChilds[i].second=0;	// default is: no childs
	}
	static unsigned int nbChilds() { return NCHILD; }
	//
	NODEDATA& data() { return itsData; }
	const NODEDATA& data() const { return itsData; }
	//
	child_iterator	begin() { return itsChilds; }
	child_iterator	end()	{ return itsChilds+NCHILD; }
	const_child_iterator	begin() const { return itsChilds; }
	const_child_iterator	end()	const { return itsChilds+NCHILD; }
	std::ostream& dumpTo(std::ostream&o) const 
	{	
		o<<"[d="<<itsData<<"|"; 
		for(const_child_iterator it=begin();it!=end();++it)
			o<<"p="<<it->first<<":"<<it->second<<" "; 
		o<<"]" ;
		return o; 
	}
private:
	// NodeElt(const NodeElt&);				//	N/A
	// NodeElt& operator=(const NodeElt&);		//	N/A
private:
	NODEDATA		itsData;				//	 
	child_t			itsChilds[NCHILD];		//	 not owner of the childs node elements
} ;


#endif // _ICM_TREEBASE_H_