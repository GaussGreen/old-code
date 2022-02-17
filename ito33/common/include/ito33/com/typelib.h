/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/com/typelib.h
// Purpose:     TypeLibrary class wrapping ITypeLib objects
// Author:      Vadim Zeitlin
// Created:     25.03.03
// RCS-ID:      $Id: typelib.h,v 1.7 2004/11/13 10:29:14 zeitlin Exp $
// Copyright:   (c) 2002 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/com/typelib.h
    @brief Type library related stuff (mostly used inside the library only).

    This file contains TypeLibrary class which wraps ITypeLib objects and
    miscellaneous helper functions for workign with type libraries.
 */

#ifndef _ITO33_COM_TYPELIB_H_
#define _ITO33_COM_TYPELIB_H_

#include "ito33/com/ptr.h"

#include "ito33/win32/winwrap.h"

#include <ole2.h>                       // for ITypeLib

namespace ito33
{

namespace COM
{

/**
    TypeLibrary class wraps type library loading.

    This class is a usual smart pointer but it also has a ctor which loads the
    type library from file which is quite convenient.
 */
class TypeLibrary : public Ptr<ITypeLib>
{
public:
  /**
      Ctor loads type library from the main DLL itself.

      For this to succeed, DLL must have a TYPELIB section in its resources.
      This is usually achieved simply by adding a line like
          @code
              1 TYPELIB "filename.tlb"
          @endcode
      to the .rc file.

      If loading type library failed, an exception is thrown.
   */
  TypeLibrary();

  /**
      Register the type library in the system registry.

      If this function fails, an exception is thrown.
   */
  void Register();

  /**
      Unregisters the type library by removing information about it from the
      system registry.

      This function doesn't throw.

      @return true if unregisterd ok, false if an error occured
   */
  bool Unregister();

  /**
      Get the type info for the interface with the give IID.

      @return pointer to type info or @c NULL if failed
   */
  Ptr<ITypeInfo> GetTypeInfo(const GUID& guid) const
  {
    ITypeInfo *pTypeInfo = NULL;
    Get()->GetTypeInfoOfGuid(guid, &pTypeInfo);

    return Ptr<ITypeInfo>(pTypeInfo);
  }

  /**
      Get the type info for the given interface.

      @return pointer to type info or @c NULL if failed
   */
  template <class T>
  Ptr<ITypeInfo> GetTypeInfo() const
  {
    return GetTypeInfo(Traits<T>::Uuid());
  }
};

} // namespace COM

} // namespace ito33

#endif // _ITO33_COM_TYPELIB_H_

