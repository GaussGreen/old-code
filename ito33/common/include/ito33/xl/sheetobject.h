/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/XL/sheetobject.h
// Purpose:     declaration of XL::SheetObject class
// Author:      Vadim Zeitlin
// Created:     2006-04-14
// RCS-ID:      $Id: sheetobject.h,v 1.7 2006/08/21 15:05:12 willy Exp $
// Copyright:   (c) 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/XL/sheetobject.h
    @brief Classes to help managing the objects created by worksheet functions.

    The problem the classes in this header help to solve is the following one:
    often, a worksheet function wants to keep around an object created during
    its execution. As a simple example, a function doing a long computation
    would want to keep the object containing the computation output around in
    order to not recompute it on each call. The problem is that we don't get
    any notification from Excel when this object is not needed any more, e.g.
    when the function is erased/overwritten by user. So the SheetObjectManager
    class remembers the cell for which the object was created and provides a
    way to check if it still contains the function which needed the object in
    the first place and destroys the object if it isn't the case. This is not
    perfect but it is the best we can with Excel.

    The major restriction of this approach is that we can have at most one
    object of a specific type in the given cell.

    @todo In fact it's not quite the best, using names instead of the cell
          coordinates would probably be better, to be investigated.
 */

#ifndef _ITO33_XL_SHEETOBJECT_H_
#define _ITO33_XL_SHEETOBJECT_H_

#include "ito33/common.h"
#include "ito33/meta/typetraits.h"

#include <map>

namespace ito33
{

namespace XL
{

/**
    Base class for SheetObjectManager template implementing type-independent
    part of it.
 */
class SheetObjectManagerBase
{
public:
  /**
      The handle used to identify the objects in Excel.

      Currently we use a string containing the numeric representation of the
      pointer to object.
   */
  typedef std::string Handle;

  /**
      Type used for passings Handle as parameters to functions.

      Using a separate typedef for it makes it simpler to update the code if
      Handle is changed to be a primitive type (as it used to be).
   */
  typedef meta::ArgType<Handle>::Value HandleArg;


  /**
      Delete all the objects which are not used by Excel any longer.

      Should be called from user code as often (or as rarely) as needed.
   */
  void PruneInactive();

protected:
  // return true if this pointer is a valid object, i.e. is known to this class
  bool DoIsValid(void *p) const
  {
    return m_formulae.find(p) != m_formulae.end();
  }

  // return an existing object for the current cell or creates a new one if
  // needed
  void *DoGetForThisCell();

  // check if we have the object with this handle and return it if we do or
  // NULL otherwise
  void *DoFromHandle(HandleArg handle) const;

  // call DoFromHandle() and throw if it returned NULL, returns its value
  // otherwise (so the return value of this method is never NULL)
  void *DoFromHandleOrThrow(HandleArg handle) const;

  // get the handle of the given object
  XLOPER *DoHandleOf(const void *p) const;

  // destroys all objects: must be called from derived class dtor
  void DoClear();


  // this map contains all the objects managed by this class indexed by their
  // location in the sheet (the key is the combination of row and column in a
  // single long)
  typedef std::map<unsigned long, void *> Objects;
  
  // this map contains the text of the formula which the given object was
  // created for, as assume it is not needed any more if the formula changed
  typedef std::map<void *, std::string> Formulae;
  Formulae m_formulae;

  // this map contains all the  map objects managed by the
  // Excel sheet they belong (the key is the sheet id as given by Excel)
  typedef std::map<unsigned long, Objects> SheetsObjects;
  SheetsObjects  m_sheetsobjects;

private:
  // create a new object
  virtual void *DoCreate() const = 0;

  // dispose of a not needed any more object
  virtual void DoDelete(void *p) const = 0;
};

/**
    The type used to pass handles from Excel to the addin functions.
 */
typedef const char *SheetObjectHandle;

/**
    Manages the lifetime of the objects created by worksheet functions.

    This is a singleton class as it doesn't make sense to have more than one
    manager for a given type, but it's not a problem to use multiple managers
    for the objects of different types.

    Template parameter @arg T must be a complete type with accessible default
    ctor and dtor.
 */
template <class T>
class SheetObjectManager : public SheetObjectManagerBase
{
public:
  /// The type of the objects managed by us.
  typedef T ObjectType;

  /**
      Destructor deletes any remaining objects.
   */
  ~SheetObjectManager() { DoClear(); }

  /**
      Check if the given object is valid, i.e. if it is still alive.
   */
  bool IsValid(ObjectType *p) const
  {
    return DoIsValid(p);
  }

  /**
      Returns the object corresponding to the given handle or @c NULL if none.

      This should be called to retrieve the object corresponding to a worksheet
      function argument if it can be a handle or another value (e.g. some
      functions take either an object or just a scalar value as argument), if
      the function argument must be an object, use operator[]() below instead.
   */
  ObjectType *FromHandle(HandleArg handle) const
  {
    return static_cast<ObjectType *>(DoFromHandle(handle));
  }

  /**
      Returns the object corresponding to the given handle or throws if the
      handle is invalid.

      This should be used when the function parameter can only be a handle.
   */
  ObjectType& operator[](HandleArg handle) const
  {
    return *static_cast<ObjectType *>(DoFromHandleOrThrow(handle));
  }

  /**
      Get the handle of the given object.

      The object pointer must be valid.
   */
  XLOPER *HandleOf(const ObjectType& ptr) const
  {
    return DoHandleOf(&ptr);
  }

  /**
      Returns the object associated with the calling cell, creating it if
      necessary.

      Usually SetForThisCell() should be used instead of this function.

      This method can only be used from the worksheet function, when the
      calling cell is well defined. It throws otherwise.
   */
  ObjectType& GetForThisCell()
  {
    return *static_cast<ObjectType *>(DoGetForThisCell());
  }

  /**
      Sets the object as being associated with the current cell and returns its
      internal handle.

      This method can only be used from the worksheet function, when the
      calling cell is well defined. It throws otherwise as well as if an error
      occurs.
   */
  XLOPER *SetForThisCell(const ObjectType& value)
  {
    ObjectType& obj = GetForThisCell();
    obj = value;
    return HandleOf(obj);
  }

private:
  virtual void *DoCreate() const { return new ObjectType; }
  virtual void DoDelete(void *p) const { delete static_cast<ObjectType *>(p); }
};

} // namespace XL

} // namespace ito33

#endif // _ITO33_XL_SHEETOBJECT_H_

