//----------------------------------------------------------------------------
//
//   Group       : Credit Hybrids QR
//
//   Filename    : TemplateFunctions4Containers.hpp
//
//   Description : Global template functions that wraps STL algo functions
//					and {STL containeters + QLib container Array}
//
//   Date        : March 2006
//
//
//----------------------------------------------------------------------------

#ifndef QLIB_TEMPLATEFUNCTIONS4CONTAINERS_HPP
#define QLIB_TEMPLATEFUNCTIONS4CONTAINERS_HPP

#include "edginc/Array.hpp"
#include "edginc/Format.hpp"
#include "edginc/LessEqualGreaterEps.hpp"
#include "edginc/DECLARE.hpp"
#include <algorithm> 
#include <list>

DRLIB_BEGIN_NAMESPACE
/*****************************************************************************
******************************************************************************
** Auxiliary template functions to be able to check if an instance of the   **
** container strictly decreases / decreases / contains duplicates           **
** / increases /strictly increases /is such that all its elements are more  **
** than some floor / satisfies arbitrary condition.                         **
** These are the wrappers around adjacent_find and find_if.                 **
** Efficiency =  gContainer.size(). gContainer is not obliged to be sorted. **
******************************************************************************
*****************************************************************************/	
template <class TContainer, class TElem> bool allElemsOf1AreGreaterThen2(	
						   const TContainer & gContainer, const TElem & gElem);
///////////////////////////////////////////////////////////////////////////////
template <class TContainer, class TElem> bool allElemsOf1AreLessThen2(	
						   const TContainer & gContainer, const TElem & gElem);
///////////////////////////////////////////////////////////////////////////////
template <class TContainer> bool strictlyDecreases(	
											    const TContainer & gContainer);
///////////////////////////////////////////////////////////////////////////////
template <class TContainer> bool decreases(	const TContainer & gContainer);
///////////////////////////////////////////////////////////////////////////////
template <class TContainer> bool containsDuplicates(
											const TContainer & gContainer);
///////////////////////////////////////////////////////////////////////////////
template <class TContainer> bool increases(	const TContainer & gContainer);
///////////////////////////////////////////////////////////////////////////////
template <class TContainer> bool strictlyIncreases(	
											const TContainer & gContainer);
/******************************************************************************
*******************************************************************************
** Functions merge2 merge two sorted STL containers. 
** More precisely, 1st container must increase and 
** 2nd one must strictly increase. 
** If it is not the case, an exception will be thrown.
** 
**  Union of elements of both containers might contain several subsets such 
**  that all elements there are equivalent two by two. Note that each of these 
**  subsets is
**      either one or several elements of container1
**      or one element of container2.
** 
** 	The result will increase. 
*******************************************************************************
******************************************************************************/	
template <	class TContainer1, 
			class TContainer4Result> 
                void	mergeIntoVector( const TContainer1 & g2merge,
                                         TContainer4Result *const gResult,
                                         bool gCheckThatInputsIncrease = true);
///////////////////////////////////////////////////////////////////////////////
template <	class TContainer1, 
            class TContainer4Result> 
                void	mergeIntoVector( const TContainer1 &      g2merge,
                                         TContainer4Result *const gResult,
                                         double gPrecision, 
                                         bool gCheckThatInputsIncrease = true);
///////////////////////////////////////////////////////////////////////////////
template<class TTypeInList, class TContainer2Merge>
               void mergeIntoList(  const TContainer2Merge & gContainer2Merge,
                                    list<TTypeInList> * const gResult, 
                                    bool gCheckThatInputsIncrease=true);
///////////////////////////////////////////////////////////////////////////////
template<class TTypeInList, class TContainer2Merge>
               void mergeIntoList(  const TContainer2Merge & gContainer2Merge,
                                    list<TTypeInList> * const gResult, 
                                    double gPrecision, 
                                    bool gCheckThatInputsIncrease=true);
///////////////////////////////////////////////////////////////////////////////
template <class CContainer, class CComponent> 
typename CContainer::const_iterator upperBound(	 const CContainer & gContainer,
												 const CComponent & gLimit);
///////////////////////////////////////////////////////////////////////////////
template <class CContainer, class CComponent> 
	typename CContainer::iterator 	upperBound(		   CContainer & gContainer,
												 const CComponent & gLimit);
///////////////////////////////////////////////////////////////////////////////
template <class TContainer, class TComponent> 
typename TContainer::const_iterator lowerBound(	 const TContainer & gContainer,
												 const TComponent & gLimit);
///////////////////////////////////////////////////////////////////////////////
template <class TContainer, class TComponent> 
	typename TContainer::iterator lowerBound(		TContainer & gContainer,
											 const	TComponent & gLimit);
///////////////////////////////////////////////////////////////////////////////
template <class TElemType, class TContainter> int getIndexInContainer(
									 const TElemType  & gElem, 
									 const TContainter& gInstanceOfContainter);
///////////////////////////////////////////////////////////////////////////////
template <class TContainter> 
const typename TContainter::value_type & getElemInContainerGivenIndex(
                                      const TContainter& gInstanceOfContainter,
                                      int                gIndex);
///////////////////////////////////////////////////////////////////////////////
template <class TContainter> 
typename TContainter::value_type & getElemInContainerGivenIndex(
                                            TContainter& gInstanceOfContainter,
                                            int          gIndex);
///////////////////////////////////////////////////////////////////////////////
template <class TContainter> 
refCountPtr<list<typename TContainter::value_type> > makeList(
                                                  const TContainter  & gInput);
///////////////////////////////////////////////////////////////////////////////
template <class TContainer1, class TContainer2> 
	smartPtr< array<int> >   getIndexesOf1in2(const TContainer1 & gContainer1,
									          const TContainer2 & gContainer2);
///////////////////////////////////////////////////////////////////////////////
template <class TContainer> void sortContainer(TContainer & gContainer);
///////////////////////////////////////////////////////////////////////////////
template <class TContainer, class TElem> int floorAllElemsOf1By2(
							TContainer & gContainer, const TElem & gElem);
///////////////////////////////////////////////////////////////////////////////
template<class TPointer2BaseType, class TPointer2DerivedType>
list<TPointer2DerivedType> filterDerivedTypes(
                                    list<TPointer2BaseType> & gBaseTypeList);
///////////////////////////////////////////////////////////////////////////////
template<class cType> void set1To2IfNotSet(cType ** const g1, cType * g2);
///////////////////////////////////////////////////////////////////////////////

#include "edginc/TemplateFunctions4Containers.inl"

///////////////////////////////////////////////////////////////////////////////
typedef list<string> StringList;
DECLARE_REF_COUNT(   StringList);
///////////////////////////////////////////////////////////////////////////////
DRLIB_END_NAMESPACE
#endif
