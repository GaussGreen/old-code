/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/socket.h
// Purpose:     include headers needed for socket functions portably
// Author:      Vadim Zeitlin
// Created:     2005-02-16
// RCS-ID:      $Id: socket.h,v 1.8 2005/06/16 00:17:04 zeitlin Exp $
// Copyright:   (c) 2005 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/socket.h
    @brief Portable header for using sockets.

    This header simply includes the right header declaring socket functions
    depending on the platform and also defines any types/functions/constants
    missing on the current platform.
 */

#ifndef _ITO33_SOCKET_H_
#define _ITO33_SOCKET_H_

#include "ito33/string.h"

#ifdef _WIN32
  #include <winsock.h>

  typedef int socklen_t;

  namespace ito33
  {
    const int SOCKET_EAGAIN = WSAEWOULDBLOCK;
  }
#else // !_WIN32
  #include <errno.h>
  #include <unistd.h>
  #include <fcntl.h>
  #include <sys/types.h>
  #include <sys/socket.h>
  #include <sys/ioctl.h>
  #include <sys/fcntl.h>
  #include <arpa/inet.h>
  #include <netinet/in.h>

  typedef int SOCKET;
  static const SOCKET INVALID_SOCKET = -1;

  static const int SOCKET_ERROR = -1;

  #define closesocket close
  #define ioctlsocket ioctl

  namespace ito33
  {
    const int SOCKET_EAGAIN = EAGAIN;
  }
#endif // _WIN32/!_WIN32

namespace ito33
{

/**
    Socket represents a TCP/IP socket.

    This class wraps a socket descriptor ensuring that the socket is closed
    when this object is destroyed and also providing some functions abstracting
    the differences between Windows and Unix platforms.

    Note that this is a really low level class, none of its functions throw on
    error, it's up to the caller to handle this.
 */
class Socket
{
public:
  Socket(int af = AF_INET, int type = SOCK_STREAM, int protocol = PF_UNSPEC)
  {
    m_socket = socket(af, type, protocol);
  }

  ~Socket()
  {
    if ( m_socket )
      closesocket(m_socket);
  }

  /// Implicit conversion to underlying socket type.
  operator SOCKET() const { return m_socket; }

  /// Explicit test for socket validity: return true if socket is valid
  bool IsOk() const { return m_socket != INVALID_SOCKET; }

  /// Return the last error which occured on any socket.
  static int GetLastError();

  /// Switch the specified socket into non blocking mode.
  static bool MakeNonBlocking(SOCKET sock);

  /**
      Returns the address such that the specified peer can connect back to it.

      To implement this function we create a connection to the peer host on the
      specified port, so for it to work connect() should succeed.

      @param peer the peer host
      @param port the port on peer host to which we can connect
      @return this host address (in dot quad format) on the network interface
              used for connecting to the peer or empty string on error
   */
  static std::string GetAddressForPeer(const std::string& peer, int port);

private:
  // our socket descriptor
  SOCKET m_socket;
};

// ----------------------------------------------------------------------------
// implementation of inline functions
// ----------------------------------------------------------------------------

/* static */ inline
std::string
Socket::GetAddressForPeer(const std::string& peerHost, int peerPort)
{
  Socket sock(AF_INET, SOCK_STREAM, 0);

  if ( !sock )
    return "";

  struct sockaddr_in addr;
  memset(&addr, 0, sizeof(addr));
  addr.sin_family = AF_INET;
  addr.sin_addr.s_addr = inet_addr(peerHost.c_str());
  addr.sin_port = htons(static_cast<u_short>(peerPort));

  if ( connect(sock, (struct sockaddr *)&addr, sizeof(addr)) < 0 )
    return "";

  socklen_t addrlen = sizeof(addr);
  memset(&addr, 0, sizeof(addr));
  if ( getsockname(sock, (struct sockaddr*)&addr, &addrlen) < 0 )
    return "";

  return inet_ntoa(addr.sin_addr);
}

#ifdef _WIN32

/* static */ inline
int Socket::GetLastError()
{
  return WSAGetLastError();
}

/* static */
inline bool Socket::MakeNonBlocking(SOCKET sock)
{
  unsigned long flag = 1;
  return ioctlsocket(sock, FIONBIO, &flag) == 0;
}

#else // !_WIN32

/* static */ inline
int Socket::GetLastError()
{
  return errno;
}

/* static */ inline
bool Socket::MakeNonBlocking(SOCKET sock)
{
  return fcntl(sock, F_SETFL, O_NONBLOCK) == 0;
}

#endif // _WIN32/!_WIN32

} // namespace ito33

#endif // _ITO33_SOCKET_H_

