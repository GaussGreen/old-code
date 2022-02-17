/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/factory.h
// Purpose:     simple generic factory class
// Author:      Vadim Zeitlin
// Created:     2004-05-11
// RCS-ID:      $Id: factory.h,v 1.10 2005/07/27 16:35:07 zeitlin Exp $
// Copyright:   (c) 2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file  ito33/factory.h
    @brief Simple generic factory class.
 */

#ifndef _ITO33_FACTORY_H_
#define _ITO33_FACTORY_H_

#include <cstdlib>

namespace ito33
{

// this class is used just as the default value for the argument template
// parameter of Factory below
namespace Private { class NoSuchClass; }

/**
    FactorySignature template defines the argument and return types of the
    factory creator functions.

    @internal
 */
template <typename A, typename R>
struct FactorySignature
{
  typedef A Argument;
  typedef R RetType;
  typedef RetType (*Creator)(const Argument *);
};

// specialization for no arguments
template <typename R>
struct FactorySignature<Private::NoSuchClass, R>
{
  typedef void Argument;
  typedef R RetType;
  typedef RetType (*Creator)();
};

/**
    Factory allows to create objects by some key instead of having to
    explicitly specify their class.

    The key may be anything at all but, condidering that the keys should be
    unique, it is usually a string (whose value may be the class name). Another
    possibility is to use a UUID as key as the factory for creating COM classes
    does. In any cas, it should be a simple class as objects of this type are
    created during static objects initialization phase and so can't rely on
    anything else working at the moment of their construction.

    To use the factory, you have to register the creators for the keys you want
    to use. This is done by using ITO33_DEFINE_FACTORY_PRODUCT() macro below. You
    should also put ITO33_IMPLEMENT_FACTORY() macro somewhere in your code (it
    should appear exactly once).

    Template parameters are:
        - K the type of the key used to index the products
        - V the type of the objects created by this factory
        - A the type of the ctor argument, none by default
        - R the type returned by this factory ("V *" by default but may also be
            a SharedPtr<V> or AutoPtr<V> for example)
 */
template
<typename K, class V, typename A = Private::NoSuchClass, typename R = V *>
class Factory
{
public:
    /// The key which identifies the classes which we create
    typedef K Key;

    /// The (base) type of classes we create
    typedef V Product;

    /// @internal
    typedef FactorySignature<A, R> Signature;

    /// The return type
    typedef typename Signature::RetType RetType;

    /// The (optional) argument passed to the constructor
    typedef typename Signature::Argument Argument;

    /// The creation function
    typedef typename Signature::Creator Creator;


    /**
        Constructor for a concrete factory implementation.

        This constructor is normally not used directly but only by
        ITO33_DEFINE_FACTORY_PRODUCT macro. You may call it directly too if the
        macro is too limited for your purposes however.
     */
    Factory(const Key& key, Creator creator)
      : // insert us in the head of the linked list (it is simpler than
        // inserting in the tail and the order doesn't matter here anyhow)
        m_next(ms_first),
        m_key(key),
        m_creator(creator)
    {
        ms_first = this;
    }


    /**
        Get the first available factory.

        Together with GetNext() these methods allow to iterate over all
        available factories. Don't use this function if you need to just create
        the object by key, use Create() instead.

        The factories are returned in an unspecified order.

        @return factory pointer (not to be deleted by caller) or @c NULL
     */
    static Factory *GetFirst() { return ms_first; }

    /**
        Return the next factory in the list.

        @sa GetFirst

        @return next factory pointer (not to be deleted by caller) or @c NULL
     */
    Factory *GetNext() const { return m_next; }


    /// Return the key of this factory
    const Key& GetKey() const { return m_key; }

    /**
        Do create an object using this factory.

        Normally this function is not used directly, use statid Create()
        instead.
     */
    RetType CreateObject() const { return (*m_creator)(); }

    /// Same as above but for factories using arguments for object creation
    RetType CreateObject(const Argument *arg) const { return (*m_creator)(arg); }


    /**
        Creates the specified product by its key.

        If no creator is found for the given key, return @c NULL. Otherwise a
        pointer to a new Product is returned and it must be deleted by the
        caller.

        This function may throw an exception if the creator function for the
        key throws it.

        @param key the key identifying the product to create
        @return pointer to new product which may be @c NULL
     */
    static RetType Create(const Key& key)
    {
        // do linear search in the linked list (ok here as the list is small)
        //
        // NB: we don't use std::find() here because this would force us to
        //     define iterators for our linked list (and we can't use
        //     std::list, of course, because it allocates memory which we don't
        //     want to do in these static objects) and it just seems as too
        //     much trouble for too little gain
        for ( const Factory *fact = GetFirst(); fact; fact = fact->GetNext() )
        {
            if ( fact->GetKey() == key )
            {
                // found the right factory
                return fact->CreateObject();
            }
        }

        return RetType();
    }

    /**
        Same as above but for factories for products which are created with an
        argument.
     */
    static RetType Create(const Key& key, const Argument *arg)
    {
      for ( const Factory *fact = GetFirst(); fact; fact = fact->GetNext() )
      {
        if ( fact->GetKey() == key )
        {
          return fact->CreateObject(arg);
        }
      }

      return RetType();
    }

private:
    // the head of the linked list of the concrete factories
    static Factory *ms_first;

    // the next object in the linked list (may be NULL)
    Factory * const m_next;

    // the key for the product we create
    const Key m_key;

    // the function which creates products of this type
    const Creator m_creator;

    // we cannot have -- and don't need --- neither assignment operator nor the
    // copy ctor
    Factory(const Factory<K, V, A, R>&);
    Factory& operator=(const Factory<K, V, A, R>&);
};

/**
    Defines a creator for the specified product.

    This macro should be used in a .cpp file after a full declaration of
    Product. You may define a shorter version of this macro, taking only key
    and Product parameters for each given specialization of Factory<K, V>
    template.

    It is only suitable for the factories without constructor arguments (i.e.
    using the default value for the last Factory template parameter).

    @param Key the type of the key
    @param key of the product to be created
    @param ProductBase the base product class
    @param Product the concrete type of product to be created for this key
 */
#define ITO33_DEFINE_FACTORY_PRODUCT(Key, key, ProductBase, Product)          \
    /* use unnamed namespace to hide the declarations inside it */            \
    namespace                                                                 \
    {                                                                         \
        ProductBase *Create()                                                 \
        {                                                                     \
          return new Product();                                               \
        }                                                                     \
                                                                              \
        static ito33::Factory< Key, ProductBase >                             \
        FactoryFor ## Product (key, Create);                                  \
    } struct Dummy /* just to force semicolon after the macro */

/**
    Implements the factory.

    This macro must appear in your code exactly once for each Factory
    specialization you use.

    @param Key the type of the key
    @param ProductBase the base product class
 */
#define ITO33_IMPLEMENT_FACTORY(Key, ProductBase)                             \
    template <>                                                               \
    ito33::Factory< Key, ProductBase > *                                      \
    ito33::Factory< Key, ProductBase >::ms_first = NULL

/**
    Implements the factory.

    This macro must appear in your code exactly once for each Factory
    specialization you use.

    @param F the factory class to be implemented
 */
#define ITO33_IMPLEMENT_THE_FACTORY(F) template<>  F * F::ms_first = NULL

} // namespace ito33

#endif // _ITO33_FACTORY_H_

