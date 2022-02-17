/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/mapiter.h
// Purpose:     helper iterators over std::map
// Author:      Vadim Zeitlin
// Created:     2006-06-28
// RCS-ID:      $Id: mapiter.h,v 1.2 2006/07/17 16:33:29 vaclav Exp $
// Copyright:   (c) 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/mapiter.h
    @brief Iterator classes adapting std::map iterators.

    Sometimes we need to iterate just over the map keys or the map values, e.g.
    to apply a standard algorithm to all values but this is impossible using
    std::map::iterator itself. So we define our own adapter iterator classes
    allowing to do this in this header.
 */

#ifndef _ITO33_MAPITER_H_
#define _ITO33_MAPITER_H_

#include <boost/iterator/iterator_facade.hpp>

namespace ito33
{

/**
    An iterator over std::map container keys.

    Notice that this is always a const iterator, it doesn't make sense to
    modify the map keys while iterating over them.

    Creating the instances of this class directly is quite verbose so usually
    the MapKeysBegin() and MapKeysEnd() functions are used.
 */
template <typename M>
class MapKeysIterator : public boost::iterator_facade
                               <
                                  MapKeysIterator<M>,
                                  const typename M::key_type,
                                  boost::bidirectional_traversal_tag
                               >
{
public:
  /// The type of the map over whose keys we iterate
  typedef M MapType;

  /// The type of the map iterator we adapt
  typedef typename MapType::const_iterator MapIterator;


  /**
      Ctor creates a new map keys iterator, possibly initialized with the given
      iterator into an existing map.
   */
  MapKeysIterator(MapIterator i = MapIterator()) : m_iter(i) { }

  // default copy ctor, assignment operator and dtor are ok

private:
  // it needs to access our dereference, equal and increment/decrement
  friend class boost::iterator_core_access;

  // we don't "inherit" value_type from base class in ISO C++, so we have to typedef
  // it in this class; this is the same as iterator_facade::value_type, but shorter:
  typedef typename MapType::key_type value_type;

  // by defining the 4 methods below we implement everything boost::iterator
  // library needs to make this class into a fully-fledged standard iterator

  const value_type& dereference() const { return m_iter->first; }

  bool equal(const MapKeysIterator<MapType>& other) const
  {
    return other.m_iter == m_iter;
  }

  void increment() { ++m_iter; }
  void decrement() { --m_iter; }

private:
  MapIterator m_iter;
};

/**
    Return an iterator to the first key of the given @a map.
 */
template <typename M>
inline MapKeysIterator<M> MapKeysBegin(const M& map)
{
  return MapKeysIterator<M>(map.begin());
}

/**
    Return an iterator to the one after last key of the given @a map.
 */
template <typename M>
inline MapKeysIterator<M> MapKeysEnd(const M& map)
{
  return MapKeysIterator<M>(map.end());
}

} // namespace ito33

#endif // _ITO33_MAPITER_H_

