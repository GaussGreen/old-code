/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/com/uuid.h
// Purpose:     Uuid class representing an UUID
// Author:      Vadim Zeitlin
// Created:     25.12.02
// RCS-ID:      $Id: uuid.h,v 1.7 2004/10/05 09:13:36 pedro Exp $
// Copyright:   (c) 2002 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#ifndef _ITO33_COM_UUID_H_
#define _ITO33_COM_UUID_H_

/**
    @file   ito33/com/uuid.h
    @brief  Uuid is a class wrapping the standard UUID structure.

    UUID is a universaly unique identifier and is known under various names in
    different COM areas: IID (interface id), CLSID (class id) and so on.
 */

#include "ito33/string.h"

#include "ito33/win32/winwrap.h"

#include <ole2.h>

namespace ito33
{

namespace COM
{

/**
    Uuid class is a wrapper around UUID struct adding conversion to/from string
    and more.

    @sa UuidOf
 */
class Uuid : public UUID
{
public:
  /**
      @name Constructors

      Note that default copy ctor (and assignment operator) and dtor are ok.
   */
  //@{

  /// default ctor creates an invalid (null) UUID
  Uuid();

  /// ctor from a plain UUID
  Uuid(const UUID& uuid);

  /// from string (will throw an exception if failed)
  explicit Uuid(const std::string& uuid);

  //@}

  /**
      @name Accessors

      Default comparison operator is ok.
   */
  //@{

  /// return true if the object is valid (non null), false otherwise
  bool operator!() const;

  /// returns the string representation of the UUID
  std::string AsString() const;

  /// implicit conversion to string
  operator std::string() const { return AsString(); }

  //@}
};

} // namespace COM

} // namespace ito33

#endif // _ITO33_COM_UUID_H_

