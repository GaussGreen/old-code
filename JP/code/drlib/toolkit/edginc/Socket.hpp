//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : Socket.hpp
//
//   Description : Low level socket functions
//
//   Author      : Andrew J Swain
//
//   Date        : 24 January 2003
//
//
//----------------------------------------------------------------------------

#ifndef _SOCKET_HPP
#define _SOCKET_HPP
#include <string>

using namespace std;    // string

DRLIB_BEGIN_NAMESPACE


/** Low level socket functions */
class TOOLKIT_DLL Socket {
public:
    static const string HTTP_PREFIX;

    //Sets up the connection with the server and reads in data from url
    static string urlRead(const string& url);

    // Opens a socket on the specified server to the server port and 
    // returns socket id to caller
    static int connect(const string& serverLocation, int port);

    // Sends a message to a socket
    static void sendMessage(int socketid, const string& message);

    // Get a message back from a socket
    static string receiveMessage(int socketid);

    // Close a socket
    static void close(int socket);

private:
    Socket();

    // Receive a packet of data - returns length of message received
    static int receivePacket(int socketid, char* message);
};

DRLIB_END_NAMESPACE
#endif
