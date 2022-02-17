// -----------------------------------------------------------------------
// This proprietary software has been developed strictly for J.P. Morgan's
// internal use.  Any use or misuse, intentional or otherwise, which
// contradicts or places this policy in jeopardy is strictly forbidden.
//
// Copyright 2001 J.P. Morgan & Co. Incorporated.  All rights reserved.
// -----------------------------------------------------------------------
// Sequence.cpp: implementation of the Sequence class.
// Author: M.Huq, Credit DR
// -----------------------------------------------------------------------

#include "edginc/Sequence.h"

// Exception handling
#include "edginc/SCException.h"

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>

CORE_BEGIN_NAMESPACE

class RNG_DLL LicenseManager
{
public:
    LicenseManager(); // constructor checks license
    virtual ~LicenseManager()
    {}
    ; // to prevent linker from dropping the symbol
};

#ifdef USELICENSE

#include "edginc/expiration.h"
LicenseManager::LicenseManager()
{ // constructor checks license
    char fileName[80];
    sprintf(fileName, "supercube.key");
    expiration A(fileName); // FIXME put it in a static variable so checked only once.
}
#else
LicenseManager::LicenseManager()
{ // constructor checks license
    // No checks
}
#endif

static LicenseManager mgr; // checks license on the startup

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

Sequence::Sequence() :
        nDimensions(0) //,
        //    initType(DEFAULTTYPE),
        //    filePointer(NULL)
{
    // seed = 0;
    // vector = NULL;
    // inputFileName = NULL;
    // fileIndex = 0;
    // nSets = 0;


}

Sequence::Sequence(size_t nDim) :
        nDimensions(nDim),
        vector(nDim)
{
}

// Copy constructor
// FIXME: kill it
Sequence::Sequence(const Sequence & t)
{
    nDimensions = t.nDimensions;
    /*  vector = new double [nDimensions];
      for(int i =0 ; i < nDimensions; i++)
        vector[i] = t.vector[i];*/
    this->vector = t.vector;
    //  seed = t.seed;
    //  strcpy(inputFileName, t.inputFileName);
    //  initType = t.initType;
    //  fileIndex = t.fileIndex;
    //  filePointer = t.filePointer; // FIXME
    //  nSets = t.nSets;
}// Sequence::Sequence(const Sequence & t)...

// Copy assignment operator
Sequence & Sequence::operator=(const Sequence &t)
{
    // as long as we are not equating to ourselves
    if(this != &t) {
        //    if(vector){
        //   delete [] vector;
        //   vector = 0;
        //    }
        //
        //     vector = new double [nDimensions];
        //     for(int i =0 ; i < nDimensions; i++)
        //       vector[i] = t.vector[i];
        this->vector = t.vector;
        nDimensions = t.nDimensions;//    seed = t.seed;
        //    strcpy(inputFileName, t.inputFileName);
        //     initType = t.initType;
        //    fileIndex = t.fileIndex;
        //     filePointer = t.filePointer; //FIXME: this will result in double-close in dtor
        //    nSets = t.nSets;
    }
    return *this;
}// Sequence & Sequence::operator=(const Sequence &t){

Sequence::~Sequence()
{
    /* if(filePointer){
      fclose(filePointer);
     }*/
    //  if(vector){
    //    delete [] vector;
    //    // nullify
    //    vector = 0;
    //  }
}

// void Sequence::dumpSequence(const char *tag){
//  if(vector){
//   cout << "=======Sequence:dumpSequence for " << tag << " =======" << endl;
//   for(int i = 0 ; i < nDimensions; i++){
//    cout <<"vector[" << i << "] = " << vector[i] << endl;
//   }
//   cout << "======================================================" << endl;
//  }
// }


void Sequence::populateVector(int iPath)
{
    // In the base class just call the method without any path reference
    populateVector();
}

double * Sequence::getVector()
{
    //  if(vector){
    //   return vector;  // Return the address
    //  }else{
    //   return NULL;
    //  }
    // Keep vector/nDimensions consistent
    if (vector.size() < getDimension())
        vector.resize(getDimension());
    return (vector.empty() ? NULL : & vector[0]);
}



size_t Sequence::getDimension() const
{
    return nDimensions;
}

// int Sequence::getNumberSets(){
//         return nSets;
// }

/*void Sequence::dumpToFile(const char *FileName){
 ofstream outputFile(FileName, ios::out);
 if(!outputFile){
  throw SCFileIOException(__FILE__, __LINE__,
     "Unable to open file ",
     FileName);
 }
 for(int i = 0; i < nDimensions; i++)
  outputFile << vector[i] << " ";
 outputFile << endl;
 outputFile.close();
}
*/
/*void Sequence::appendToFile(const char *FileName){
 ofstream outputFile(FileName, ios::app);
 if(!outputFile){
  throw SCFileIOException(__FILE__, __LINE__,
     "Unable to open file ",
     FileName);
 }
 for(int i = 0; i < nDimensions; i++)
  outputFile << vector[i] << " ";
 outputFile << endl;
 outputFile.close();
}*/
// void Sequence::plotToFile(const char *FileName, const int index1, const int index2){
//  ofstream outputFile(FileName, ios::app);
//  if(!outputFile){
//   throw SCFileIOException(__FILE__, __LINE__,
//      "Unable to open file ",
//      FileName);
//  }
//  outputFile << vector[index1] << " " << vector[index2] << endl;
//  outputFile.close();
// }


#if 0
double * Sequence::getDeviatesbyPathDate(int iPath, int iDate)
{
    throw SCException(__FILE__, __LINE__,
                      "Sequence::getDeviatesbyPathDate - Method not implemented for base class.");
} //getDeviatebyPathDate

#endif

#if 0
void Sequence::getDeviatesbyFactorDate(int iFactor,
                                       int iDate,
                                       double *&returnVector)
{
    throw SCException(__FILE__, __LINE__,
                      "Sequence::getDeviatesbyFactorDate - Method not implemented for base class.");
} //getDeviatebyFactorDate
#endif

//======================================================
//  Method for advancing along sequence. That is, given
//  method f() for generating a sequence {x}, run f()
//  to numberElements number of elements of {x}. Discard
//  the values.

void Sequence::advanceSequence(const int numberElements)
{
    if(nDimensions == 0) {
        throw SCException (__FILE__, __LINE__,
                           "nDimensions == 0 : Sequence not initialized");
    }
    double *tmpx;
    int iSeq; // Index for an element along the sequence.
    for(iSeq = 0 ; iSeq < numberElements; iSeq++) {
        populateVector(iSeq);
        tmpx = getVector();
    }//iSeq
}//void Sequence::advanceSequence(const int numberElements)
//=====================================================
// Method for reinitializing. Dummy method unless overloaded.
void Sequence::reinitialize(const long _seed)
{
    throw SCException(__FILE__, __LINE__,
                      "reinitialize not implemented for class.");
}// void Sequence::reinitialize(const long _seed)


// Constructor that takes a file.
// Read numDimensions at a time.
FileSequence::FileSequence(const char *fileName) : Sequence()
        //    : initType(READFROMFILE)
{
    //  seed = 0;


    //  fileIndex = 0;
    filePointer = fopen(fileName,"r");
    if(!filePointer) {
        throw SCFileIOException(__FILE__, __LINE__,
                                "Unable to open file for reading",
                                fileName);
    }
    int nSets;
    if(fscanf(filePointer, "%d %d", &nDimensions, &nSets) != 2) {
        throw SCFileIOException(__FILE__, __LINE__,
                                "Unable to read from file",
                                fileName);
    }

    /*  vector = new double [nDimensions];
        if(!vector){
        throw SCException(__FILE__, __LINE__,
        "Unable to allocate memory for vector");
    }*/
    vector.assign(nDimensions, 0);
}

FileSequence::~FileSequence()
{
    if (filePointer != 0)
        std::fclose(filePointer);
}

void FileSequence::populateVector()
{

    //     if (vector == NULL)
    //         return;
    //
    //     for(int i = 0; i < nDimensions; i++)
    //         vector[i] = 0.0; // Initialize to zero
    vector.assign(nDimensions, 0.0);
    //    if (initType == READFROMFILE)
    {
        int vectorIndex = 0;
        while(vectorIndex < nDimensions  && !feof(filePointer))
        {
            if(fscanf(filePointer,"%lf", &(vector[vectorIndex]))!=1) {
                break;
            }
            vectorIndex++;
            //            fileIndex++;
        }
    }
}
CORE_END_NAMESPACE
