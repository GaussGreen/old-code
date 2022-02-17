/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/net.h
// Purpose:     miscellaneous network-related stuff
// Author:      Vadim Zeitlin
// Created:     2004-11-20
// RCS-ID:      $Id: net.h,v 1.1 2004/11/20 23:54:20 zeitlin Exp $
// Copyright:   (c) 2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/net.h
    @brief Network-related classes and functions.
 */

#ifndef _ITO33_NET_H_
#define _ITO33_NET_H_

#include "ito33/common.h"

namespace ito33
{

/**
    An instance of this class must be created before using any network
    functions.

    An object of this class initializes network library in its ctor and cleans
    it up in its dtor, as usual. It is harmless to create multiple instances of
    this class, only the first one created will really do anything non trivial
    and only the last one destroyed will really shut down the network library.

    This class is only needed under Windows to initialize winsock library and
    does nothing under Unix. It is also thread-safe in the sense that its
    instances can be created in different threads without any additional
    protection.
 */
class NetworkInitializer
{
public:
  /**
      Ctor initializes the network library.

      If the library initialization fails, an exception is thrown.
   */
  NetworkInitializer();

  /**
      Dtor shuts down the network library.

      Dtor never throws.
   */
  ~NetworkInitializer();

private:
  NO_COPY_CLASS(NetworkInitializer);
};

// ----------------------------------------------------------------------------
// inline functions implementation
// ----------------------------------------------------------------------------

// NetworkInitializer does nothing under Unix
#ifndef _WIN32

inline NetworkInitializer::NetworkInitializer()
{
}

inline NetworkInitializer::~NetworkInitializer()
{
}

#endif // !_WIN32

} // namespace ito33

#endif // _ITO33_NET_H_

