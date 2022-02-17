// -----------------------------------------------------------------------
// This proprietary software has been developed strictly for J.P. Morgan's
// internal use.  Any use or misuse, intentional or otherwise, which
// contradicts or places this policy in jeopardy is strictly forbidden.
//
// Copyright 2001 J.P. Morgan & Co. Incorporated.  All rights reserved.
// -----------------------------------------------------------------------
//
// expiration.cpp
//
// Author: M.Huq
// -----------------------------------------------------------------------

#include "edginc/SCException.h"

#include "edginc/expiration.h"
#include "edginc/ran2.h"

#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <cstring>

#define PAR1 16
#define PAR2 33
#define PAR3 20
#define PAR4 92

CORE_BEGIN_NAMESPACE

extern int DisassembleKeys(IUniformRNGSP uniform,
                               char *publicKey1,
                               char *publicKey2,
                               int constant1,
                               int constant2,
                               int constant3,
                               int constant4,
                               char *month,
                               char *day,
                               char *year);

//-----------------------------------
// Function to convert month from string to integer from 0 to 11.
long getNumMonth (char *month)
{
    char *months[]={"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"};
    for(int k=0; k<12; k++) {
        if(strncmp(month ,months[k],3)==0) {
            //printf("Found %s %s\n", month, months[k]);
            return k;
        }
    }
    return -1;
}//getNumMonth(char)

//-----------------------------------
// Constructor
expiration::expiration(const char *fileName)
{
    strcpy(publicKeyFile, fileName);
    validateKey();
}  // expiration::expiration(char*)

//-----------------------------------
// Destructor
expiration::~expiration()
{} // expiration::~expiration()

//-----------------------------------
// Validation method
void expiration::validateKey()
{
    // Storage for public keys
    char publicKey1[64];
    char publicKey2[64];

    // Open up the public key file
    FILE *publicFP = fopen(publicKeyFile,"r");
    if(!publicFP)
        throw SCFileIOException(__FILE__, __LINE__,
                                "License validator: Could not open public key file", publicKeyFile);
    // Read in key 1
    fscanf(publicFP,"%s",publicKey1);
    // Read in key 2
    fscanf(publicFP,"%s",publicKey2);
    // Done with file
    fclose(publicFP);
    //
    size_t keyLength = strlen(publicKey1);
    if(keyLength != strlen(publicKey2))
        throw SCException(__FILE__, __LINE__,
                          "Invalid keys found. Check validity of license file.");
    // Set up seed internally
    long seed = -39340739;
    long constant1 = PAR1;
    long constant2 = PAR2;
    long constant3 = PAR3;
    long constant4 = PAR4;

    char month[4];
    char day[3];
    char year[5];
    IUniformRNGSP uniform = SC_ran2::create(SC_ran2Gen::create(seed));
    // Get the expiration date
    DisassembleKeys(uniform,
                    publicKey1,
                    publicKey2,
                    constant1,
                    constant2,
                    constant3,
                    constant4,
                    month,
                    day,
                    year);

    // Get the current date
    time_t t;
    t = time(NULL);
    char *currentDate = asctime(localtime(&t));
    char curMonth[4];
    char curDay[3];
    char curYear[5];
    curMonth[0] = currentDate[4];
    curMonth[1] = currentDate[5];
    curMonth[2] = currentDate[6];
    curMonth[3] = '\0';
    curDay[0] = currentDate[8];
    curDay[1] = currentDate[9];
    curDay[2] = '\0';
    curYear[0] = currentDate[20];
    curYear[1] = currentDate[21];
    curYear[2] = currentDate[22];
    curYear[3] = currentDate[23];
    curYear[4] = '\0';
    // Check to see if we have exceeded the expiration date
    // Year
    long thisYear;
    sscanf(curYear,"%ld",&thisYear);
    long expYear;
    sscanf(year,"%ld",&expYear);
    if(expYear < thisYear) {
        show_expiration_message();
    }
    // Month
    long thisMonth = getNumMonth (curMonth);
    long expMonth  = getNumMonth (month);
    if( expYear == thisYear && expMonth < thisMonth)
        show_expiration_message();
    // Day
    long thisDay;
    sscanf(curDay,"%ld",&thisDay);
    long expDay;
    sscanf(day,"%ld",&expDay);
    if(expYear == thisYear && expMonth == thisMonth && expDay < thisDay) {
        show_expiration_message();
    }
}//expiration::validateKey

//--------------------------------------------
// Method to display expiration message
void expiration::show_expiration_message()
{
    cerr << "************S U P E R C U B E    L I B R A R Y************" << endl;
    cerr << " License file has expired. Please obtain new license file" << endl;
    cerr << " Contact Mijan Huq or Andrew Abrahams re getting license " << endl;
    cerr << "**********************************************************" << endl;
    throw SCException("License validation function :", 0,
                      "SuperCube license expired");
} // expiration::show_expiration_message()


CORE_END_NAMESPACE
