/////////////////////////////////////////////////////////////////////////////
// Name:        com/uuid.cpp
// Purpose:     implementation of Uuid class
// Author:      Vadim Zeitlin
// Created:     25.12.02
// RCS-ID:      $Id: uuid.cpp,v 1.5 2004/10/05 09:13:43 pedro Exp $
// Copyright:   (c) 2002 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

// ============================================================================
// declarations
// ============================================================================

// ----------------------------------------------------------------------------
// headers
// ----------------------------------------------------------------------------

#include "ito33/win32/exception.h"

#include "ito33/com/uuid.h"

#include <memory.h>

using namespace ito33::COM;

// ============================================================================
// implementation
// ============================================================================

// ----------------------------------------------------------------------------
// ctors
// ----------------------------------------------------------------------------

Uuid::Uuid()
{
  memset(this, 0, sizeof(UUID));
}

Uuid::Uuid(const UUID& uuid)
{
  memcpy(this, &uuid, sizeof(UUID));
}

Uuid::Uuid(const std::string& uuid)
{
  if ( ::UuidFromString(reinterpret_cast<unsigned char *>
              (const_cast<char *>(uuid.c_str())),
             this) != RPC_S_OK )
  {
    throw WIN32_EXCEPTION("UuidFromString");
  }
}

// ----------------------------------------------------------------------------
// accessors
// ----------------------------------------------------------------------------

bool Uuid::operator!() const
{
  char buf[sizeof(UUID)];
  memset(buf, 0, sizeof(UUID));

  return memcmp(this, buf, sizeof(UUID)) == 0;
}

std::string Uuid::AsString() const
{
  unsigned char *p;

  if ( ::UuidToString(const_cast<Uuid *>(this), &p) != RPC_S_OK )
  {
    throw WIN32_EXCEPTION("UuidToString");
  }

  std::string s(reinterpret_cast<char *>(p));

  ::RpcStringFree(&p);

  return s;
}

