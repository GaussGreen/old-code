// -----------------------------------------------------------------------
// This proprietary software has been developed strictly for J.P. Morgan's
// internal use.  Any use or misuse, intentional or otherwise, which
// contradicts or places this policy in jeopardy is strictly forbidden.
//
// Copyright 2001 J.P. Morgan & Co. Incorporated.  All rights reserved.
// -----------------------------------------------------------------------
//
// expiration.h
//
// Author: M.Huq
// License Class set up for SuperCube.
// Serves two purposes:
// 1. Provides a key to allow usage of the SuperCube library
// 2. Provides a means by which we can introduce expiration dates for usage
//    with any program really.
//
// Instance of this class must be created successfully for program to run.
// Calling the constructor is sufficient.
// -----------------------------------------------------------------------
#ifndef EXPIRATION__H
#define EXPIRATION__H

#include "edginc/coreConfig.hpp"

CORE_BEGIN_NAMESPACE

class RNG_DLL expiration
{
private:
    char publicKeyFile[128]; // Specification of license file location
    char expirationDate[12]; // MMM DD YYYY format
    void validateKey( );     // Function to read in and validate keys
    void show_expiration_message(); // Expiration message
public:
    expiration(const char *fileName) ;
    ~expiration();
    char *getExpirationDate();
};

CORE_END_NAMESPACE

#endif // EXPIRATION__H
