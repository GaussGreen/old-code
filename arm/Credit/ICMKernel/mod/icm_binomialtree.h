
#ifndef _ICM_BINOMIALTREE_H_
#define _ICM_BINOMIALTREE_H_

#include "ICMKernel/mod/icm_treebase.h"
#include "ICMKernel/util/ICM_QMatrix.h"
#include <algorithm>

#include <map>

//	-----------------------------------------------------------------------
template <class NODEDATA>	class BinomialTree; 
//	-----------------------------------------------------------------------
//
//	-----------------------------------------------------------------------
template <class NODEDATA> 
class	ConstBinomialSlice 
{
public:
#if _MSC_VER >= 1310	// Visual C++ 2005 & .net 2003
	typedef	typename NodeElt<2,NODEDATA>			node_elt_t; 
	typedef	typename NODEDATA						node_data_t;
	typedef	typename BinomialTree<NODEDATA>			tree_t ;
	typedef typename node_elt_t::const_iterator		const_iterator; 
#else
	typedef	NodeElt<2,NODEDATA>				node_elt_t; 
	typedef	NODEDATA						node_data_t;
	typedef	BinomialTree<NODEDATA>			tree_t ;
	typedef node_elt_t::const_iterator		const_iterator; 
#endif

	friend class BinomialTree<NODEDATA> ;
public:
	unsigned long size() const { return itsSize; }
	const node_elt_t&		node(unsigned long j) const
	{	
		if (j>=itsSize)
			ICMTHROW(ERR_INVALID_ARGUMENT,"ConstBinomialSlice::node("<<j<<"): Out of bound") ; 			;
		return *(itsFirst+j); 
	}
	const node_data_t&	data(unsigned long j)  const
	{	
		if (j>=itsSize)
			ICMTHROW(ERR_INVALID_ARGUMENT,"ConstBinomialSlice::data("<<j<<"): Out of bound") ; 			;
		return (itsFirst+j)->data(); 
	}
	const_iterator begin() const { return itsFirst; }
	const_iterator end()   const { return itsLast+1; }
protected:
#if _MSC_VER >= 1310	// Visual C++ 2005 & .net 2003
	ConstBinomialSlice(typename node_elt_t::const_iterator first, typename node_elt_t::const_iterator last) 
#else
	ConstBinomialSlice(node_elt_t::const_iterator first,node_elt_t::const_iterator last) 
#endif
	{	itsFirst=first; itsLast=last; itsSize=itsLast-itsFirst+1; }
private:
#if _MSC_VER >= 1310	// Visual C++ 2005 & .net 2003
	typename node_elt_t::const_iterator itsFirst,itsLast; 
#else
	node_elt_t::const_iterator itsFirst,itsLast; 
#endif
	unsigned long itsSize ;
}; 
//	-----------------------------------------------------------------------
//
//	-----------------------------------------------------------------------
template <class NODEDATA>
class BinomialSlice : public ConstBinomialSlice<NODEDATA>
{
public:
#if _MSC_VER >= 1310	// Visual C++ 2005 & .net 2003
	typedef	typename NodeElt<2,NODEDATA>			node_elt_t; 
	typedef	typename NODEDATA						node_data_t;
	typedef	typename BinomialTree<NODEDATA>			tree_t ;
	typedef typename node_elt_t::iterator			iterator; 
#else
	typedef	NodeElt<2,NODEDATA>				node_elt_t; 
	typedef	NODEDATA						node_data_t;
	typedef	BinomialTree<NODEDATA>			tree_t ;
	typedef node_elt_t::iterator			iterator; 
#endif
	friend class BinomialTree<NODEDATA> ;
public:
	node_elt_t&		node(unsigned long j) 
	{
		return const_cast<node_elt_t&>( ConstBinomialSlice<NODEDATA>::node(j) ); 
	}
	node_data_t&	data(unsigned long j) 
	{
		return const_cast<node_data_t&>( ConstBinomialSlice<NODEDATA>::data(j) ); 
	}
	iterator begin() { return const_cast<iterator>(ConstBinomialSlice<NODEDATA>::begin()) ; }
	iterator end()   { return const_cast<iterator>(ConstBinomialSlice<NODEDATA>::end()) ; }
protected:
#if _MSC_VER >= 1310	// Visual C++ 2005 & .net 2003
	BinomialSlice(typename node_elt_t::iterator first,typename node_elt_t::iterator last) : 
#else
	BinomialSlice(node_elt_t::iterator first,node_elt_t::iterator last) : 
#endif
	   ConstBinomialSlice<NODEDATA>(first,last) {}
}; 
//	-----------------------------------------------------------------------
//
//	-----------------------------------------------------------------------
template <class NODEDATA>
class BinomialTree
{
public:
#if _MSC_VER >= 1310	// Visual C++ 2005 & .net 2003
	typedef	typename NodeElt<2,NODEDATA>				node_elt_t; 
	typedef	typename NODEDATA						node_data_t;
	typedef typename BinomialSlice<NODEDATA>			slice_t; 
	typedef typename ConstBinomialSlice<NODEDATA>	const_slice_t; 
	typedef	typename node_elt_t::iterator			iterator; 
	typedef	typename node_elt_t::const_iterator		const_iterator; 
#else
	typedef	NodeElt<2,NODEDATA>				node_elt_t; 
	typedef	NODEDATA						node_data_t;
	typedef BinomialSlice<NODEDATA>			slice_t; 
	typedef ConstBinomialSlice<NODEDATA>	const_slice_t; 
	typedef	node_elt_t::iterator			iterator; 
	typedef	node_elt_t::const_iterator		const_iterator; 
#endif

public:
	unsigned long depth() const					{ return itsTimes.Getnbrows(); }
	double Time(unsigned long i) const			{ return itsTimes(i,0); }
	double DT(unsigned long i) const			{ return itsDT(i,0); }
	const ICM_QMatrix<double>&	Time()  const	{ return itsTimes; }
	const ICM_QMatrix<double>&	DT()	 const	{ return itsDT; }
	bool findpos(double time,unsigned long&pos) const ; 
public:
	BinomialTree();
	BinomialTree(const BinomialTree&ref);	
	BinomialTree& operator=(const BinomialTree&ref);
	void setTimes(const ICM_QMatrix<double>& times); 
public:
	bool				IsEmpty() const { return itsNodes.empty(); }	
	node_elt_t&			node(unsigned long i,unsigned long j) ; 
	node_data_t&		data(unsigned long i,unsigned long j); 
	const node_elt_t&	node(unsigned long i,unsigned long j) const ; 
	const node_data_t&	data(unsigned long i,unsigned long j) const ; 
	slice_t				slice(unsigned long i)  ; 
	const_slice_t		slice(unsigned long i) const ; 
private:
	unsigned long pos(unsigned long i,unsigned long j) const  ;
public:
	std::ostream& dumpTo(std::ostream&o) const ; 
public:
#if _MSC_VER >= 1310	// Visual C++ 2005 & .net 2003
	typename node_elt_t::iterator begin()				{ return &(*itsNodes.begin());  } 
	typename node_elt_t::iterator end()					{ return &(*itsNodes.end()) ;	}
	typename node_elt_t::const_iterator begin()	const	{ return &(*itsNodes.begin());  } 
	typename node_elt_t::const_iterator end()	const	{ return &(*itsNodes.end()) ;	}
#else
	node_elt_t::iterator begin()				{ return itsNodes.begin();  } 
	node_elt_t::iterator end()					{ return itsNodes.end() ;	}
	node_elt_t::const_iterator begin()	const	{ return itsNodes.begin();  } 
	node_elt_t::const_iterator end()	const	{ return itsNodes.end() ;	}
#endif

	unsigned long	size()				const	{ return itsNodes.size() ;  }
	void clear() ; 
private:
	void setup_relationships(); 
private:
	// node_elt_t					*itsRoot;		//	not owner
	std::vector<node_elt_t>		itsNodes; 
	ICM_QMatrix<double>			itsTimes;		//	nTime dimension
	ICM_QMatrix<double>			itsDT;			//	nTime dimension; the last value is 0.
} ;

//	-----------------------------------------------------------------------
template <class NODEDATA>
static 
inline std::ostream& operator<<(std::ostream& o,const BinomialTree<NODEDATA> &tree)
{
	tree.dumpTo(o); 
	return o; 
}
//	-----------------------------------------------------------------------
template <class NODEDATA> 
inline 
BinomialTree<NODEDATA>::BinomialTree()
{
	//	Empty Tree. Nothing to do.
}
//	-----------------------------------------------------------------------
template <class NODEDATA>
inline unsigned long 
BinomialTree<NODEDATA>::pos(unsigned long i,unsigned long j) const 
{
	return i*(i+1)/2+j; 
}
//	-----------------------------------------------------------------------
template <class NODEDATA> 
inline void
BinomialTree<NODEDATA>::setTimes(const ICM_QMatrix<double>& times)
{
	//	Special case when the Tree is empty
	if (times.IsEmpty()) 
	{
		itsNodes.resize(0); 
		itsTimes.Resize(0,0); 
		itsDT.Resize(0,0);
		return; 
	}
	//	Copy & Compute the variable time steps
	unsigned long NT = times.Getnbrows(); 
	itsTimes=times; 
	itsDT.Resize(NT,1) ;
	unsigned int i ;
	for (i=0;i<NT-1;i++) 
		itsDT(i,0)=itsTimes(i+1,0)-itsTimes(i,0); 
	itsDT(NT-1,0)=0 ; 
	//	Allocation of the node_elt_t array
	itsNodes.resize(NT*(NT+1)/2); 
	//	Setting up the node relationships
	setup_relationships(); 
}
//	-----------------------------------------------------------------------
template <class NODEDATA>
inline
void
BinomialTree<NODEDATA>::setup_relationships() 
{
	if ( itsTimes.IsEmpty()) return; 
	unsigned long NT= itsTimes.Getnbrows(); 
	unsigned long i,j ; 
	for (i=0;i<NT-1;i++)
		for (j=0;j<i+1;j++) 
		{
			node_elt_t::child_iterator it(itsNodes[pos(i,j)].begin()); 
			it->second=&itsNodes[pos(i+1,j)]; 
			++it; 
			it->second=&itsNodes[pos(i+1,j+1)]; 
		}
	for (j=0;j<NT;j++) 
	{
		node_elt_t::child_iterator it(itsNodes[pos(NT-1,j)].begin()); 
		it->second=0; 
		++it; 
		it->second=0; 
	}
}
//	-----------------------------------------------------------------------
template <class NODEDATA>
inline 
#if _MSC_VER >= 1310	// Visual C++ 2005 & .net 2003
typename BinomialTree<NODEDATA>::node_elt_t& 
#else
BinomialTree<NODEDATA>::node_elt_t& 
#endif
BinomialTree<NODEDATA>::node(unsigned long i,unsigned long j)
{
	if (i>=depth() || j>i) 
		ICMTHROW(ERR_INVALID_ARGUMENT,"BinomialTree::node("<<i<<","<<j<<"): Out of bound") ; 
	return itsNodes[pos(i,j)]; 
}
//	-----------------------------------------------------------------------
template <class NODEDATA>
inline 
#if _MSC_VER >= 1310	// Visual C++ 2005 & .net 2003
typename BinomialTree<NODEDATA>::node_data_t& 
#else
BinomialTree<NODEDATA>::node_data_t& 
#endif
BinomialTree<NODEDATA>::data(unsigned long i,unsigned long j)
{
	return node(i,j).data(); 
}
//	-----------------------------------------------------------------------
template <class NODEDATA>
inline 
#if _MSC_VER >= 1310	// Visual C++ 2005 & .net 2003
const typename BinomialTree<NODEDATA>::node_elt_t& 
#else
const BinomialTree<NODEDATA>::node_elt_t& 
#endif
BinomialTree<NODEDATA>::node(unsigned long i,unsigned long j) const
{
	if (i>=depth() || j>i) 
		ICMTHROW(ERR_INVALID_ARGUMENT,"BinomialTree::node("<<i<<","<<j<<"): Out of bound") ; 
	return itsNodes[pos(i,j)]; 
}
//	-----------------------------------------------------------------------
template <class NODEDATA>
inline 
#if _MSC_VER >= 1310	// Visual C++ 2005 & .net 2003
const typename BinomialTree<NODEDATA>::node_data_t& 
#else
const BinomialTree<NODEDATA>::node_data_t& 
#endif
BinomialTree<NODEDATA>::data(unsigned long i,unsigned long j) const 
{
	return node(i,j).data(); 
}
//	-----------------------------------------------------------------------
template <class NODEDATA>
inline 
#if _MSC_VER >= 1310	// Visual C++ 2005 & .net 2003
typename BinomialTree<NODEDATA>::slice_t 
#else
BinomialTree<NODEDATA>::slice_t
#endif
BinomialTree<NODEDATA>::slice(unsigned long i) 
{
	if (i>=depth() ) 
		ICMTHROW(ERR_INVALID_ARGUMENT,"BinomialTree::slice("<<i<<"): Out of bound") ; 
	return slice_t( &(*itsNodes.begin())+pos(i,0), &(*itsNodes.begin())+pos(i,i)) ;
}
//	-----------------------------------------------------------------------
template <class NODEDATA>
inline 
#if _MSC_VER >= 1310	// Visual C++ 2005 & .net 2003
typename BinomialTree<NODEDATA>::const_slice_t 
#else
BinomialTree<NODEDATA>::const_slice_t
#endif
BinomialTree<NODEDATA>::slice(unsigned long i) const
{
	if (i>=depth() ) 
		ICMTHROW(ERR_INVALID_ARGUMENT,"BinomialTree::slice("<<i<<"): Out of bound") ; 
// FIXMEFRED: mig.vc8 (28/05/2007 15:23:40):cast
	return const_slice_t(&(*itsNodes.begin())+pos(i,0), &(*itsNodes.begin())+pos(i,i)); 
}
//	-----------------------------------------------------------------------
template <class NODEDATA>
inline 
std::ostream& 
BinomialTree<NODEDATA>::dumpTo(std::ostream&o) const 
{
	o<<"Dumping Binomial Tree"<<std::endl; 
	o<<"depth="<<depth()<<std::endl; 
	unsigned int i,j; 
	for(i=0;i<depth();i++) 
	{
		// o<<"T("<<i<<")="<<itsTimes(i)<<std::endl ;
	}
	for(i=0;i<depth();i++) 
	{
		for (j=0;j<i+1;j++)
		{
			o<<"("<<i<<","<<j<<"):"<<&node(i,j)<<":" ; 
			node(i,j).dumpTo(o); 
			o<<std::endl; 
		}
		// o.flush(); 
	}
			
		o<<"End of Dumping Binomial Tree"<<std::endl; 
	return o; 
}
//	-----------------------------------------------------------------------
template <class NODEDATA>
inline 
BinomialTree<NODEDATA>::BinomialTree(const BinomialTree&ref) :
	itsTimes(ref.itsTimes),
	itsDT(ref.itsDT),
	itsNodes(ref.itsNodes)
{
	// we need to update the node relationships.
	setup_relationships(); 
}
//	-----------------------------------------------------------------------
template <class NODEDATA> 
inline 
BinomialTree<NODEDATA>& 
BinomialTree<NODEDATA>::operator =(const BinomialTree& ref) 
{
	if (this!=&ref) 
	{
		this->~BinomialTree(); 
		new(this)BinomialTree(ref); 
	}
	return *this; 
}
//	-----------------------------------------------------------------------
template <class NODEDATA>
inline
void
BinomialTree<NODEDATA>::clear()
{
	*this=BinomialTree<NODEDATA>(); 
}
//	-----------------------------------------------------------------------
template <class NODEDATA>
inline 
bool 
BinomialTree<NODEDATA>::findpos(double time,unsigned long&pos) const 
{
	if (depth()==0)  return false; 	
	if (lt(time,itsTimes(0,0))) return false; 
	pos=0; 
	while ( (pos<itsTimes.Getnbrows()) && neq(time,itsTimes(pos,0)) ) pos++; 
	if (pos==itsTimes.Getnbrows()) return false;
	return true; 
}
//	-----------------------------------------------------------------------
#endif // _ICM_BINOMIALTREE_H_