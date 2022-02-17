// -----------------------------------------------------------------------
// This proprietary software has been developed strictly for J.P. Morgan's
// internal use.  Any use or misuse, intentional or otherwise, which
// contradicts or places this policy in jeopardy is strictly forbidden.
//
// Copyright 2001 J.P. Morgan & Co. Incorporated.  All rights reserved.
// -----------------------------------------------------------------------
// license.cpp
//
// Author: M.Huq
//
// Functions to create and decipher licenses given particular keys
// The main objective here is introduce a license file which will have
// encoded in it an expiration date.
//
// -----------------------------------------------------------------------
//

#include "edginc/ran2.h"

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <math.h>

CORE_BEGIN_NAMESPACE


using namespace std;
//----------------------------------------------------------------------------
// ASCII character set in decimal format
//   |  0 nul |  1 soh |  2 stx |  3 etx |  4 eot |  5 enq |  6 ack |  7 bel |
//   |  8 bs  |  9 ht  | 10 nl  | 11 vt  | 12 np  | 13 cr  | 14 so  | 15 si  |
//   | 16 dle | 17 dc1 | 18 dc2 | 19 dc3 | 20 dc4 | 21 nak | 22 syn | 23 etb |
//   | 24 can | 25 em  | 26 sub | 27 esc | 28 fs  | 29 gs  | 30 rs  | 31 us  |
//   | 32 sp  | 33 !   | 34 "   | 35 #   | 36 $   | 37 %   | 38 &   | 39 '   |
//   | 40 (   | 41 )   | 42 *   | 43 +   | 44 ,   | 45 -   | 46 .   | 47 /   |
//   | 48 0   | 49 1   | 50 2   | 51 3   | 52 4   | 53 5   | 54 6   | 55 7   |
//   | 56 8   | 57 9   | 58 :   | 59 ;   | 60 <   | 61 =   | 62 >   | 63 ?   |
//   | 64 @   | 65 A   | 66 B   | 67 C   | 68 D   | 69 E   | 70 F   | 71 G   |
//   | 72 H   | 73 I   | 74 J   | 75 K   | 76 L   | 77 M   | 78 N   | 79 O   |
//   | 80 P   | 81 Q   | 82 R   | 83 S   | 84 T   | 85 U   | 86 V   | 87 W   |
//   | 88 X   | 89 Y   | 90 Z   | 91 [   | 92 \   | 93 ]   | 94 ^   | 95 _   |
//   | 96 `   | 97 a   | 98 b   | 99 c   |100 d   |101 e   |102 f   |103 g   |
//   |104 h   |105 i   |106 j   |107 k   |108 l   |109 m   |110 n   |111 o   |
//   |112 p   |113 q   |114 r   |115 s   |116 t   |117 u   |118 v   |119 w   |
//   |120 x   |121 y   |122 z   |123 {   |124 |   |125 }   |126 ~   |127 del |
//----------------------------------------------------------------------------
//
// Function to generate random character between 33 and 125 or ! and }
char getRandom(IUniformRNGSP uniform)
{
    long value = ((long)floor( 999.0 * uniform->fetch())) % 92;
    value +=33;
    return (char)value;
} //getRandom(long*)

//----------------------------------------------------------------------------
// Function to generate random stream of characters up to specified size
char * generateStream(IUniformRNGSP uniform, int streamLength)
{
    // Allocate memory for string to be returned.
    char *retString = new char [streamLength+1];
    if(!retString) {
        cerr << "generateStream: Could not allocate memory for return string"
        <<endl;
        abort();
    }
    for(int iSeq = 0; iSeq < streamLength; iSeq++) {
        retString[iSeq] = getRandom(uniform);
    } //iSeq
    // Pad with termination string. Reason for additional character at end.
    retString[streamLength] = '\0';
    // Return this stream
    return retString;
}//generateStream(long *, int)


//----------------------------------------------------------------------------
// Function to create a date string
// Expects to have date in format Jan 04 2001
char * createDateString(char *month, char *day, char *year)
{
    char *retString = new char[10];
    if(!retString) {
        cerr << "createDateString: Could not allocate memory for return string"
        <<endl;
        abort();
    }
    sprintf(retString,"%s%s%s", month,day,year);
    //printf("Date string = %s\n", retString);
    return retString;
} // createDateString

//----------------------------------------------------------------------------
// Function for assembling the final key as well as secondary keys which should
// go with the program using it.
int AssembleKeys(IUniformRNGSP uniform,
                 char **publicKey1,
                 char **publicKey2,
                 char **privateKey,
                 int padStringLength,
                 int constant2,
                 int constant3,
                 int constant4,
                 char *month,
                 char *day,
                 char *year)
{
    int k;
    int dateStringLength = 9;
    int totalLength = dateStringLength + 2*padStringLength;

    // Get the date string
    char *dateString = createDateString(month, day, year);
    // Get two sets of random streams for padding
    char *prePadStream = generateStream(uniform, padStringLength);
    char *postPadStream = generateStream(uniform, padStringLength);

    // concatencate strings
    char *string1 = new char [totalLength+1];
    if(!string1) {
        cerr << "AssembleKeys: Could not allocate memory for string1"
        <<endl;
        abort();
    }
    for( k = 0; k < padStringLength; k++) {
        string1[k] = prePadStream[k];
    }//k
    for( k = 0; k < dateStringLength; k++) {
        string1[k + padStringLength] = dateString[k];
    }//k
    for( k = 0; k < padStringLength; k++) {
        string1[k + padStringLength + dateStringLength] = postPadStream[k];
    }//k

    // Generate random stream. This will be one of the private keys
    *privateKey = generateStream(uniform, totalLength);

    // Allocate memory for public keys
    *publicKey1 = new char [totalLength+1];
    *publicKey2 = new char [totalLength+1];
    // Generate the public key by taking the sum of privateKey and string1
    // which contains the expiration date. From this we have an assumed min
    // constant of 33 and need to extract data to reconstruct the modulus.
    for( k = 0; k < totalLength; k++) {
        // Sum up each element of the stream in decimal format
        long xy = (long)string1[k] + (long)(*privateKey)[k];
        // Now take the modulus to define a string in the range 33 to 125
        long z = xy % constant4;
        if(z +constant2 < 33 || z+constant2 > 125)
            cerr << "ERROR key1" << endl;
        (*publicKey1)[k] = (char)(z + constant2);
        long c = (constant3+(long)(k*constant3/totalLength))*((xy - z)/constant4) + constant2 ;
        if(c> 125 || c<33 )
            cerr << "ERROR key2" << endl;
        (*publicKey2)[k] = (char)c;
    }
    // pad with null string
    (*publicKey1)[totalLength] = '\0';
    (*publicKey2)[totalLength] = '\0';


    delete string1;
    delete prePadStream;
    delete postPadStream;
    return 0;
}

//----------------------------------------------------------------------------
// Function for diassembling the key given the secondary keys
int DisassembleKeys(IUniformRNGSP uniform,
                    char *publicKey1,
                    char *publicKey2,
                    int padStringLength,
                    int constant2,
                    int constant3,
                    int constant4,
                    char *month,
                    char *day,
                    char *year)
{
    int k;
    int dateStringLength = 9;
    int totalLength = dateStringLength + 2*padStringLength;

    // Get two sets of random streams for padding
    char *prePadStream = generateStream(uniform, padStringLength);
    char *postPadStream = generateStream(uniform, padStringLength);

    // Generate random stream. This will be one of the private keys
    char *privateKey = generateStream(uniform, totalLength);
    //printf("Private key = %s\n", privateKey);

    // string1 to reconstruct
    char *string1 = new char [totalLength+1];
    if(!string1) {
        cerr << "AssembleKeys: Could not allocate memory for string1"
        <<endl;
        abort();
    }

    for( k = 0; k < totalLength; k++) {
        // Get back z
        long z = (long)publicKey1[k] - constant2;
        // Get back c
        long c = (long)publicKey2[k];
        c = c - constant2 ;
        c = c/(constant3+(long)(k*constant3/totalLength));
        long xy = c * constant4 + z;

        string1[k] = (char)(xy - (long)privateKey[k]);
    }//k
    // pad with null string
    string1[totalLength] = '\0';
    // Extract date string
    char *dateString = new char[dateStringLength+1];
    for( k = 0 ; k < dateStringLength; k++) {
        dateString[k] = string1[k + padStringLength];
    }// k
    month[0] = dateString[0];
    month[1] = dateString[1];
    month[2] = dateString[2];
    month[3] = '\0';
    day[0] = dateString[3];
    day[1] = dateString[4];
    day[2] = '\0';
    year[0] = dateString[5];
    year[1] = dateString[6];
    year[2] = dateString[7];
    year[3] = dateString[8];
    year[4] = '\0';
    //printf("Date string = %s/%s/%s\n",month,day,year);

    delete prePadStream;
    delete postPadStream;
    delete dateString;
    return 0;
}
CORE_END_NAMESPACE
