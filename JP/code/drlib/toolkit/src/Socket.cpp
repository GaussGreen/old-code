//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : Socket.cpp
//
//   Description : Low level socket functions
//
//   Author      : Andrew J Swain
//
//   Date        : 24 January 2003
//
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#include "edginc/Socket.hpp"
#include "edginc/ModelException.hpp"
#include "edginc/Malloc.hpp"
#if defined(UNIX) || defined(__CYGWIN32__)
#include<sys/types.h>
#include<sys/socket.h>
#include<netinet/in.h>
#include<netdb.h>
#include <unistd.h>
#elif defined(WIN32)
#include <winsock.h>
#endif

DRLIB_BEGIN_NAMESPACE

const string Socket::HTTP_PREFIX = "http://";

static const int    HTTP_PREFIX_LEN = 7;
static const int    DEFAULT_PORT = 80;
static const string DEFAULT_FILE_NAME = "/";
static const string CMD = " HTTP/1.0\r\n\r\n";
static const int    CMD_BASE_LEN = 20;
static const int    MAX_MESG_LEN = 1024;
static const int    MSG_INCREMENT = 256000;

class URL {
public:
    int    port;   // to connect to
    string base;   // URL name
    string file;   // to get from URL

    URL(const string& urlName): port(DEFAULT_PORT) {
        static const string method("URL::URL");
        try {
            int length = urlName.length();
 
            if (length <= HTTP_PREFIX_LEN) {
                throw ModelException(method, "Malformed url name: " + urlName);
            }

            // store contents after http:// bit 
            string name = urlName.substr(HTTP_PREFIX_LEN);

            std::string::size_type colon    = name.find(":");
            std::string::size_type fwdSlash = name.find("/");

            if (colon != std::string::npos) {
                // store base name 
                base = name.substr(0,colon);
       
                // parse out file name & port number 
                if (fwdSlash != std::string::npos && fwdSlash > colon) {
                    string ptr = name.substr(colon+1, fwdSlash - colon);
                    port = atoi(ptr.c_str());

                    file = name.substr(fwdSlash);
                }
                else {
                    string ptr = name.substr(colon+1);
                    port = atoi(ptr.c_str()); 
                    file = DEFAULT_FILE_NAME;
                }
            }
            else if (fwdSlash != std::string::npos)
            {
                // no : but a /
                base = name.substr(0, fwdSlash);
                file = name.substr(fwdSlash);
            }
            else
            {
                // no : or /
                base = name.substr(0, HTTP_PREFIX_LEN);
                file = DEFAULT_FILE_NAME;
            }
        }
        catch (exception& e) {
            throw ModelException(e, method);
        }
    }
};

//Sets up the connection with the server and reads in data from url
string Socket::urlRead(const string& url) {
    static const string method("Socket::urlRead");
    try {
        URL    parsedURL(url);
        string command = "GET " + parsedURL.file + CMD;

        int socketid = connect(parsedURL.base, parsedURL.port);

        sendMessage(socketid, command);
        
        string buffer = receiveMessage(socketid);
        
        const char* start = buffer.c_str();

        while (*(start+1) != '\0' && 
               !(*start == '\n' && (*(start+1) == '\n' ||*(start+1) == 0xd)))
        {
            start++;
        }

        start++;
  
        switch (*start)  {
        case '\n':
            start+=1;
            break;
        case 0xd:
            start+=2;
            break;
        default:
            throw ModelException(method, "internal error");
        }

        string data(start);
        close(socketid);
        return data;
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

// Opens a socket on the specified server to the server port and 
// returns socket id to caller
int Socket::connect(const string& serverLocation, int port) {
    static const string method("Socket::connect");
    try {
        int    newSocket;
        struct sockaddr_in address = {0};  // Internet address of socket 
        struct hostent    *host =  NULL;   // Internet address of host to connect to

        // Get internet address of host of server 
        if (!(host=gethostbyname(serverLocation.c_str()))) {
            throw ModelException(method, "unknown host " + serverLocation);
        }

        address.sin_family=host->h_addrtype;
        address.sin_port=htons((unsigned short)port);
        memcpy((char *)&(address.sin_addr), host->h_addr, host->h_length);

        // Set up connecting socket 
        if ((newSocket = socket(AF_INET,SOCK_STREAM,0)) == -1) {
            throw ModelException(method, "unable to create connecting socket");
        }

        // Connect socket 
        if (::connect(newSocket,(struct sockaddr*)&address,sizeof(address))==-1) {
            throw ModelException(method, "unable to connect socket to address");
        }

        return newSocket;
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
};

// Sends a message to a socket
void Socket::sendMessage(int socketid, const string& message) {
    static const string method("Socket::sendMessage");
    try {
        const char *buffer = message.c_str();
        int   bytesSent;
        int   msgLen = strlen(buffer);

        while (msgLen > 0) {
            bytesSent = send(socketid, buffer, msgLen, 0);

            if (bytesSent == -1)  {
                throw ModelException(method, "Failed sending message");
            }

            buffer += bytesSent;
            msgLen -= bytesSent;
        }
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

// Get a message back from a socket
string Socket::receiveMessage(int socketid) {
    static const string method("Socket::receiveMessage");
    char* buffer = 0;
    try {
        int   bufferSize = 2 * MAX_MESG_LEN + 1;
        int   bytesReceived = 0;
        char* tmpMessage = 0;

        buffer = (char*)Malloc::allocate(bufferSize);
        int    messLength = 0;
    
        while ((bytesReceived = receivePacket(socketid,buffer + messLength)) &&
               bytesReceived > 0) {
            messLength += bytesReceived;
            if (messLength > bufferSize - MAX_MESG_LEN - 1) {
                // grow some more space 
                tmpMessage = REALLOC(buffer, 
                                     char,
                                     bufferSize + MSG_INCREMENT);
                buffer    = tmpMessage;
                bufferSize += MSG_INCREMENT;
            }
        }

        *(buffer + messLength) = '\0';
 
        string message(buffer);

        delete[] buffer;
        return message;
    }
    catch (exception& e) {
        delete[] buffer;
        throw ModelException(e, method);
    }
}

// Close a socket
void Socket::close(int socket) {
#if defined(UNIX) || defined(__CYGWIN32__) 
    ::close(socket);
#else
    closesocket(socket);
#endif
}


// Receive a packet of data - returns length of message received
// "message" is written to
int Socket::receivePacket(int socketid, char* message) {
    static const string method("Socket::receivePacket");
    try {
        int messLength = recv(socketid, message, MAX_MESG_LEN, 0);
        if (messLength == -1) {
            throw ModelException(method, "Failed to receive packet");
        }
        return messLength;
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

DRLIB_END_NAMESPACE
