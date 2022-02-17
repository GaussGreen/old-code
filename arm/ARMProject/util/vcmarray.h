/*
 * $Log: vcmarray.h,v $
 * Revision 1.12  2004/05/10 16:42:28  jpriaudel
 * modif pour compil
 *
 * Revision 1.11  2004/05/10 16:09:51  jpriaudel
 * modif pour compil NT/Unix
 *
 * Revision 1.10  2003/11/19 20:03:40  ebenhamou
 * version for multi-platform compatibility
 *
 * Revision 1.9  2003/06/23 12:08:58  ebenhamou
 * dos2unix
 *
 *
 */



#ifndef _MARRAY_H
#define _MARRAY_H

#include "firsttoinc.h"
#include "vector"
#include "list"
#include "iterator"
#include "iostream"


#ifndef __ARM_NO_NAMESPACE
using namespace std;
#endif


/// Sun compiler
#if ( defined(__SUNPRO_CC ) && __SUNPRO_CC >= 0x500 )
	#define CC_TYPENAME typename
#else
	/// .Net
	#if	( defined(_MSC_VER) && _MSC_VER >= 1300 )
		#define CC_TYPENAME typename
	#else
		#define CC_TYPENAME 
	#endif
#endif

#define CC_TEMPLATE_SPECIALIZATION template<>


#define forall(__cont, __it) \
    for ((__it) = (__cont).begin(); (__it) < ((__cont).end()); __it++)

#define rforall(__cont, __it, __it_rend) \
    for (__it = (__cont).rbegin(); __it < __it_rend; __it++)



/*********************************************************/
/*                                                       */
/*                                                       */
/*                 MAGIC ARRAY CLASS                     */
/*             S. Mysona CDC Paris 2000                  */
/*                                                       */
/*                                                       */
/*********************************************************/



/*********************************************************/
/*                                                       */
/*   SUN COMPILER MANAGING TO BE EVEN WORSE THAN VC++    */
/*   THE STRUCTURE HAS BEEN EXPLODED                     */
/*   THIS VERSION IS A DIRTY WORKAROUND,                 */
/*   DO NOT USE FOR DEVELOPMENT, DO NOT IMPROVE          */
/*                                                       */
/*********************************************************/




/*********************************************************/
/*                                                       */
/* Implementation of a generic N-dim object.             */
/*                                                       */
/* The instantiation creates an 0-length N-dim object    */
/* (lenghts are 0 in every dimension)                    */
/*                                                       */
/* The lengths are set by using the method SetLength     */
/* which takes a list of integer as input and returns    */
/* true if everything is OK, false otherwise.            */
/*                                                       */
/* The riterator (for raw-iterator) gives a direct       */
/* handle to the individual elements (class T) of the    */
/* MArray without looping on the dimensions.             */
/* The most obvious way to use it is in foreach or       */
/* forall loops.                                         */
/*                                                       */
/*********************************************************/




/****************************************************************************************/
/*                                                                                      */
/*  Performance issues                                                                  */
/*                                                                                      */
/*  aggressive inlining of operators and iterator accessors (begin ,end, rbegin, rend)  */
/*                                                                                      */
/****************************************************************************************/



/****************************************************************************************/
/*                                                                                      */
/*  Robustness  issues                                                                  */
/*                                                                                      */
/* Checkings could be improved :                                                        */
/*     initialization : size checking, length validity ...                              */
/*                                                                                      */
/****************************************************************************************/


/****************************************************************************************/
/*                                                                                      */
/*  Implementation issues                                                               */
/*                                                                                      */
/*  This version of the code is quite obfuscated since VC++ 6 does not support partial  */
/*  specialization (and V7 is not likely to be conform to ISO/ANSI standard neither).   */
/*  But this form should compile on every compiler providing "basic" template support.  */
/*                                                                                      */
/*  The basic idea is to define a nested-recursive template                             */
/*                                                                                      */
/*       template <class T, int N> class MArray                                         */
/*                                                                                      */
/*  which contains a vector of <T,N-1> MArray's, and to specialize the seed             */                                                          
/*                                                                                      */
/*       template <class T, int 1> class MArray                                         */
/*                                                                                      */
/*  To circumvent the lack of partial specialization support, the trick is to define    */
/*  an intermediary "int N"-only template class                                         */
/*            template <int N> class N_MArray                                           */
/*  which itself contains a "class T" template                                          */
/*            template <class T> class T_MArray                                         */
/*  With this nested setting, we can provide a non partial N = 1 specialization.        */
/*  The final step is to derive the wished class :                                      */
/*                                                                                      */
/*  template <class T, int N> class MArray : public N_MArray<N>::T_MArray<T> {}         */
/****************************************************************************************/


/****************************************************************************************/
/*                                                                                      */
/*  Improvement needed                                                                  */
/*                                                                                      */
/* . valarray or Matlab like operators                                                  */
/* . OR designing of PETE support to generate efficient template expressions to avoid   */
/*   the performance problems generated by overloaded operators (numerous unecessary    */
/*   for loops, multiplication of temporary instanciation of classes, slowdown of       */
/*   function calls, ...                                                                */
/* . useful methods like STL vector one's (push_back, pop_back, resize, insert, ...)    */
/*                                                                                      */
/****************************************************************************************/




/****************************************************************************************/
/*                                                                                      */
/* About lengths :                                                                      */
/*                                                                                      */
/* The SetLength , GetLength methods handle rectangular hyper parrallelepiped.          */                                                                                    
/* For non regilar sizes, use Get,Set current length.                                   */
/*                                                                                      */
/****************************************************************************************/




/*********************************************************/
/*        Structure for slide movements                  */
/*********************************************************/

class SlideMove
{
public :    
    int nbDimsRemaining;
    int nbStep;

    SlideMove(int dims, int steps) : nbDimsRemaining(dims), nbStep(steps) {}

    SlideMove(void) : nbDimsRemaining(0), nbStep(0) {}

    void Set(int dims, int steps)
    {
        nbDimsRemaining = dims;
        nbStep = steps;
    }

};






struct TCoord2
{
  int x[2];

  inline int & operator[] (int k) {return x[k];}
  inline int & Elt        (int k) {return x[k];}
};





/*********************************************************/
/*        Generic N dimensional MArray                   */
/*********************************************************/

// Without template class <- member class
//template <class T, int N> class MArray;
template <class T> class M2Array;
template <class T> class M1Array;





template<class T>
        class M2riterator  {

        typedef M2Array < T > TcurrMArray;
        typedef M1Array < T > TnextMArray;
        typedef vector <T>         vectorT;



        typedef  TCoord2     TCoord;
	typedef  M2riterator         riterator;

        public :
            CC_TYPENAME TnextMArray::riterator itsNextRaw;
            CC_TYPENAME vector < TnextMArray >::iterator itsDirectPointed2;
            CC_TYPENAME vector < TnextMArray >::iterator itsFatherEnd;
            CC_TYPENAME vector < TnextMArray >::iterator itsFatherBegin;


            


            M2riterator(void) {};
            ~M2riterator(void) {};


            
            
            inline riterator & operator = (const riterator &itB)
            {
                itsNextRaw        = itB.itsNextRaw;
                itsDirectPointed2 = itB.itsDirectPointed2;
                itsFatherEnd      = itB.itsFatherEnd;
                itsFatherBegin    = itB.itsFatherBegin;

                //itsCoord          = itB.itsCoord;

                return *this;
            }


            inline riterator & operator = (const T* pointer)
            {
                itsNextRaw = pointer;
                return *this;
            }

            inline T & operator ->(void)
            {
                return itsNextRaw.operator->();
            }

    

            bool operator < (riterator &itB)
            {
                if (itsDirectPointed2 < itB.itsDirectPointed2)
                    return true;
                else
                {
                    if (itsDirectPointed2 == itB.itsDirectPointed2)
                        return (itsNextRaw < itB.itsNextRaw);
                    else
                        return false;
                }
            }

            inline bool operator > (riterator &itB);

            bool operator == (riterator &itB)
            {
                if (itsDirectPointed2 == itB.itsDirectPointed2)
                {
                    if (itsNextRaw == itB.itsNextRaw)
                        return true;
                    else
                        return false;
                }
                else
                    return false;

            }

            
            

            inline  T & operator *()
            {
                return *itsNextRaw;
            }


            //the dummy int variable enables to specify the postfix form of the operator
            inline riterator & operator ++(int i)
            {
                itsNextRaw++;
                if  (itsNextRaw.itsDirectPointed2 == itsDirectPointed2->end())
                {
                    itsDirectPointed2++;                    
                    //itsCoord++;                    
                    if (itsDirectPointed2 < itsFatherEnd)
                    {
                        itsNextRaw = itsDirectPointed2->rbegin();
                    }
                }

                return *this;
            }

            inline riterator & operator --(int i)
            {
                itsNextRaw--;
                if  (itsNextRaw.itsDirectPointed2 < itsDirectPointed2->begin())
                {
                    //itsCoord--;
                    itsDirectPointed2--;

                    if (itsDirectPointed2 >= itsFatherBegin)
                    {
                        itsNextRaw = itsDirectPointed2->rbegin();
                    }
                }

                return *this;
            }



            inline riterator operator +(int i)
            {
                riterator sum;


                sum = *this;

                while (i>0)
                {
                    sum++;
                    i--;
                }

                return sum;
            }


            inline riterator operator -(int i)
            {
                riterator sum;

                sum = *this;

                while (i>0)
                {
                    sum--;
                    i--;
                }

                return sum;
            }


            inline riterator & operator +=(int i)
            {
                while (i>0)
                {
                    (*this)++;
                    i--;
                }

                return *this;
            }


            inline riterator & operator -=(int i)
            {
                while (i>0)
                {
                    (*this)--;
                    i--;
                }

                return *this;
            }





            // This function slides the current directPointed2 of nbStep positions if nbDimsRemaining == 0
            // (ie, we are on the right dimension to ride on)
            // Otherwise, it calls the next raw slide, with nbDimsRemaining-1

            // Problem with int N=2 specialization, this is only for N=2 ...
            void slide2(int nbStep, int nbDimsRemaining = 0);


            riterator & operator << (SlideMove &move)
            {
                int d = itsNextRaw.itsDirectPointed2-itsDirectPointed2->GetVector()->begin();

                if (move.nbDimsRemaining == 0)
                {
                    itsDirectPointed2 += move.nbStep;
                    //itsCoord += move.nbStep;

                    itsNextRaw = itsDirectPointed2->rbegin() + d;

                }
                else // nbDimsRemaining == 1 or ....
                {
                    move.nbDimsRemaining--;
                    itsNextRaw << move;
                    move.nbDimsRemaining++;
                }

                return *this;
            }


            riterator & operator << (list<SlideMove> *moveList)
            {
                CC_TYPENAME list<SlideMove>::iterator it;

                
                for (it = moveList->begin(); it != moveList->end(); it++)
                    *this << *it;

                return *this;
            }

            inline int GetStep(void)
            {
                return itsDirectPointed2-itsFatherBegin;
            }

            TCoord GetCoord(void)
            {
                TCoord coord;

                coord[0] = itsDirectPointed2-itsFatherBegin;
                itsNextRaw.GetCoordInt(coord.x+1);

                return coord;
            }

            void GetCoordInt(int *coord2fill = NULL)
            {
                coord2fill[0] = itsDirectPointed2-itsFatherBegin;
                itsNextRaw.GetCoord(coord2fill+1);
            }





        };






template <class T> class M2Array
    {
        // Don't know why direct vectors crash
        // OK with this alias
    public :

        //typedef MArray < T , 2   > TcurrMArray;
        typedef M2Array < T > TcurrMArray;
        typedef M1Array < T > TnextMArray;
        typedef vector <T>         vectorT;

        typedef TCoord2   TCoord;


        



    protected :

        int itsLength;
        int itsOrder;

        vector < TnextMArray> itsVector;



    public :

        typedef M2riterator<T> riterator;


        M2Array()
        {
            itsOrder = 2;
            itsLength = 0;
        }

        ~M2Array() {}


   
    
        inline riterator rbegin()
        {
            riterator theBegin;

            theBegin.itsDirectPointed2 = itsVector.begin();
            theBegin.itsNextRaw = theBegin.itsDirectPointed2->rbegin();
            theBegin.itsFatherEnd = itsVector.end();
            theBegin.itsFatherBegin = theBegin.itsDirectPointed2;



            return theBegin;
        }

        inline riterator rend()
        {
            riterator theEnd;

            theEnd.itsDirectPointed2 = itsVector.end()-1;
            theEnd.itsNextRaw = theEnd.itsDirectPointed2->rend();
            theEnd.itsFatherEnd = itsVector.end();
            theEnd.itsFatherBegin = itsVector.begin();


            
            return theEnd;
        }




        inline CC_TYPENAME vector < TnextMArray >::iterator begin()
        {
            return itsVector.begin();
        }

        inline CC_TYPENAME vector < TnextMArray >::iterator end()
        {
            return itsVector.end();
        }




        inline TnextMArray & operator[] (int i)
        {
            return itsVector[i];
        }

        inline T & operator[] (riterator &it)
        {
            return *it;
        }



        inline TnextMArray & at (int i)
        {
            return itsVector[i];
        }

        inline T & rat (int *coord)
        {
            return itsVector[coord[0]].rat(coord+1);
        }

        riterator GetRiterator(TCoord &coord)
        {
            riterator outIt;

            outIt.itsDirectPointed2 = itsVector.begin() + coord.x[0];

            
            outIt.itsNextRaw = outIt.itsDirectPointed2->GetRiteratorInt(&(coord.x[0])+1);
            outIt.itsFatherEnd = itsVector.end();
            outIt.itsFatherBegin = itsVector.begin();

            return outIt;
        }


        riterator GetRiteratorInt(int *coord)
        {
            riterator outIt;

            outIt.itsDirectPointed2 = itsVector.begin() + coord[0];


            outIt.itsNextRaw = outIt.itsDirectPointed2->GetRiteratorInt(coord+1);
            outIt.itsFatherEnd = itsVector.end();
            outIt.itsFatherEnd = itsVector.begin();

            return outIt;
        }

        

        // inMarray should be const
        //N_MArray<N>::T_MArray<T> & operator = (TcurrMArray& inMArray)
        M2Array<T> & operator = (TcurrMArray& inMArray)
        {
            CC_TYPENAME TcurrMArray::riterator rit, ritIn, finish;

            SetLength(inMArray.GetLength());
            finish = this->rend();
            for (rit = rbegin(), ritIn = inMArray.rbegin(); rit < finish; rit++)
                *rit = *ritIn;

            return *this;
        }


    


        void reset()
        {
            CC_TYPENAME vector < TnextMArray  >::iterator next;
        
            forall(itsVector, next)
                next->reset();

            itsVector.clear();
        }



        inline void SetCurrentLength(int l)
        {
            itsLength = l;
            itsVector.resize(itsLength);

        }

        inline int GetCurrentLength(void)
        {
            return itsLength;
        }

        // Cool but unsafe (no size checking)


        bool SetLength(int *ilengths)
        {
            list <int> lLength(ilengths, ilengths+itsOrder);

            return SetLength(&lLength);
        }

        // bool SetLength(const TcurrMArray &inMArray)



        bool SetLength(list<int> *length)
        {
            CC_TYPENAME vector < TnextMArray >::iterator next;
            bool OK = true;
            int tmp;

            if (length->size() != 2)
                return false;

            itsLength = length->front();
            itsVector.resize(itsLength);
            tmp = *(length->begin());
            length->pop_front();

            forall(itsVector, next)
                OK = OK && next->SetLength(length);
            length->push_front(tmp);

            /* if there is a problem, reset the whole structure */
            if (!OK)
                reset();

            return OK;
        }

        


        void RetrieveLength(list<int> *length)
        {
            itsVector[0].RetrieveLength(length);
            length->push_front(itsLength);
        }

        list<int> *GetLength(void)
        {
            list<int> *length = new list <int>;
            RetrieveLength(length);

            return length;
        }




        void print(int skip = 0)
        {
            int i=0,j;
            CC_TYPENAME vector <TnextMArray>::iterator it;

            if (itsLength >0)
            {
                it = itsVector.begin();
                cout << "[" << i << "]";
                    it->print(skip+1);
                i++;
            }

            for(it = itsVector.begin()+1; it <itsVector.end(); it++)
            {
                for(j=0; j<skip; j++)
                    cout << "   ";
                cout << "[" << i << "]";
                it->print(skip+1);
                i++;
            }
        }


        inline vector < TnextMArray> * GetVector(void)
        {
            return &itsVector;
        }


        //template <class P>
        //void Alloc(P *dummy)
        void Alloc()
        {

            CC_TYPENAME vectorT::iterator it;

            forall(itsVector, it)
                *it = new T;
        }

        //template <class P>
        //void Alloc(P *dummy, int size)
        void Alloc(int size)
        {

            CC_TYPENAME vectorT::iterator it;

            forall(itsVector, it)
                *it = new T(size);
        }




        void Cleanit(void)
        {
            riterator theEnd = rend();
            riterator it;

            //rforall(*this, it, theEnd)
            for (it = this->rbegin(); it < theEnd; it++)
                cleanit(*it);
        }



    };


/*********************************************************/
/*                   Final leaf                          */
/*      partial specialization 1 dimension MArray        */
/*********************************************************/


        struct TCoord1
        {
            int x[1];

            inline int & operator[] (int k) {return x[k];}
            inline int & Elt        (int k) {return x[k];}
        };



template <class T>
        class M1riterator  {

        typedef vector  <T>         vectorT;
        typedef M1Array <T>   TcurrMArray;

        typedef M1riterator<T> riterator;
        typedef TCoord1     TCoord;

       

        public :
        
            CC_TYPENAME vector <T>::iterator itsDirectPointed2;
            CC_TYPENAME vector <T>::iterator itsFatherEnd;
            CC_TYPENAME vector <T>::iterator itsFatherBegin;

            //int itsCoord;

            M1riterator() {};
            ~M1riterator() {};


            riterator & operator = (const riterator &itB)
            {
                itsDirectPointed2 = itB.itsDirectPointed2;
                itsFatherEnd      = itB.itsFatherEnd;
                itsFatherBegin      = itB.itsFatherBegin;

                //itsCoord          = itB.itsCoord;

                return *this;
            }

            inline riterator & operator = (const T* pointer)
            {
                itsDirectPointed2 = const_cast<T*> (pointer);
                
                return *this;
            }

            inline bool operator < (riterator &itB)
            {
                return itsDirectPointed2 < itB.itsDirectPointed2;
            }

            inline bool operator >  (riterator &itB)
            {
                return itsDirectPointed2 > itB.itsDirectPointed2;
            }

            inline bool operator == (riterator &itB)
            {
                return (itsDirectPointed2 == itB.itsDirectPointed2);
            }


            inline T & operator ->(void)
            {
                return *itsDirectPointed2;
            }


            inline T & operator *()
            {
                return *itsDirectPointed2;
            }

            //the dummy int variable enables to specify the postfix form of the operator
            inline riterator & operator ++(int i)
            {
                itsDirectPointed2++;

                return *this;
            }

            inline riterator & operator --(int i)
            {
                itsDirectPointed2--;

                return *this;
            }

            inline riterator operator +(int i)
            {
                riterator sum;

                sum = *this;
                sum.itsDirectPointed2 += i;

                return sum;
            }


            inline riterator operator -(int i)
            {
                riterator sum;

                sum = *this;
                sum.itsDirectPointed2 -= i;

                return sum;
            }


            inline riterator & operator +=(int i)
            {
                itsDirectPointed2 += i;

                return *this;
            }


            inline riterator & operator -=(int i)
            {
                itsDirectPointed2 -= i;

                return *this;
            }


            riterator & operator << (SlideMove &move)
            {

                if (move.nbDimsRemaining == 0)
                {
                    itsDirectPointed2 += move.nbStep;
                }
                else 
                {
            		//throw Exception(__LINE__, __FILE__, ERR_CONDITION_NOT_MEET,
                    //                    "You cannot go that far !");
                }

                return *this;
            }

            riterator & operator << (list<SlideMove> *moveList)
            {

                CC_TYPENAME list<SlideMove>::iterator it;

                for (it = moveList->begin(); it != moveList->end(); it++)
                    *this = *this << *it;

				return	*this;
            }


            inline int GetStep(void)
            {
                return itsDirectPointed2-itsFatherBegin;
            }

            TCoord GetCoord(void)
            {
                TCoord coord;

                coord[0] = itsDirectPointed2-itsFatherBegin;

                return coord;
            }

            void GetCoordInt(int *coord2fill = NULL)
            {
                coord2fill[0] = itsDirectPointed2-itsFatherBegin;
            }








        };





/*template <int N> class MArray<1>
{


public :
*/
    template <class T> class M1Array 
    {
    public :

        typedef vector  <T>         vectorT;
        typedef M1Array <T>   TcurrMArray;

        typedef M1riterator<T> riterator ;


        typedef TCoord1 TCoord;


        


    protected :

        int itsLength;
        int itsOrder;


        vector < T > itsVector;






    public :






         M1Array()
        {
            //cout << "Creating M_Array of order 1" << endl;
            itsOrder = 1;
            itsLength = 0;
        }

         ~M1Array() {};


         /*           raw iterator                */




        inline riterator rbegin()
        {
            riterator theBegin;

            theBegin.itsDirectPointed2 = itsVector.begin();
            theBegin.itsFatherEnd = itsVector.end();
            theBegin.itsFatherBegin = theBegin.itsDirectPointed2;
  
            return theBegin;
        }

        inline riterator rend()
        {
            riterator ultimateEnd;

            ultimateEnd.itsDirectPointed2 = itsVector.end();
            ultimateEnd.itsFatherEnd = itsVector.end();
            ultimateEnd.itsFatherBegin = itsVector.begin();

            return ultimateEnd;
        }


        inline CC_TYPENAME vector < T >::iterator begin()
        {
            return itsVector.begin();
        }


        inline CC_TYPENAME vector < T >::iterator end()
        {
            return itsVector.end();
        }


        


        inline T & operator[] (int i)
        {
            return itsVector[i];
        }

        inline T & operator[] (riterator it)
        {
            return *it;
        }


        inline T & at (int i)
        {
            return itsVector[i];
        }

        inline T & rat (int *coord)
        {
            return itsVector[coord[0]];
        }


        inline riterator GetRiterator(TCoord &coord)
        {
            riterator outIt;

            outIt.itsDirectPointed2 = itsVector.begin() + coord.x[0];
            outIt.itsFatherEnd = itsVector.end();
            outIt.itsFatherBegin = itsVector.begin();

            return outIt;

        }

        riterator GetRiteratorInt(int *coord)
        {
            riterator outIt;

            outIt.itsDirectPointed2 = itsVector.begin() + coord[0];
            outIt.itsFatherEnd = itsVector.end();
            outIt.itsFatherBegin = itsVector.begin();

            return outIt;
        }





        M1Array<T> & operator = (const TcurrMArray& inMArray)
        {
            register int i;

            SetLength(&inMArray.itsLength);

            for (i = 0; i < itsLength; i++)
                itsVector[i] = inMArray.itsVector[i];

            return *this;
        }


        



        void reset()
        {
            itsVector.clear();
        }



        inline void SetCurrentLength(int l)
        {
            itsLength = l;
            itsVector.resize(itsLength);

        }

        inline int GetCurrentLength(void)
        {
            return itsLength;
        }





        bool SetLength(const int *ilength)
        {
            itsVector.resize(*ilength);

            return true;
        }


        bool SetLength(list<int> *length)
        {
            if (length->size() != 1)
                return false;

            itsLength = length->front();
            itsVector.resize(itsLength);

            return true;
        }


        void RetrieveLength(list<int> *length)
        {
            length->push_front(itsLength);
        }


        list<int> *GetLength(void)
        {
            list<int> *length = new list <int>;
            RetrieveLength(length);

            return length;
        }





        void print(int skip = 0)
        {
            int i=0,j;
            register CC_TYPENAME vector < T >::iterator it;

            if (itsLength >0)
            {
                it = itsVector.begin();
                cout << "[" << i << "]";
                cout << "   " << *it << endl;
                i++;
            }

            for(it = itsVector.begin()+1; it <itsVector.end(); it++)
            {
                for(j=0; j<skip; j++)
                    cout << "   ";

                cout <<"[" << i << "]" <<"   ";
                cout << *it << endl;

                i++;

            }
        }


        inline vector < T > * GetVector(void)
        {
            return &itsVector;
        }

        
        /*
        template <class P>
        void Alloc(P *dummy)
        {

            vectorT::iterator it;

            forall(itsVector, it)
                *it = new (T);
        }
        */
        /*
        template <class P>
        void Alloc(P *dummy, int size)        
        {

            vectorT::iterator it;

            forall(itsVector, it)
                *it = new T(size);
        }
        */


        void Cleanit(void)
        {
            riterator theEnd = rend();
            riterator it;

            for (it = this->rbegin(); it < theEnd; it++)
                cleanit(*it);
        }

        

    };

        template <class T>
        void Alloc(const M1Array<T*> &marray)
        {

            M1riterator it, itRend;
            itRend = marray.rend();


            forall(marray, it)
                *it = new (T);
        }



/*    
};
*/



        CC_TEMPLATE_SPECIALIZATION class M1riterator<double>  {

        typedef vector  <double>         vectorT;
        typedef M1Array <double>   TcurrMArray;

        typedef M1riterator<double> riterator;
        typedef TCoord1     TCoord;

       

        public :
        
            vector <double>::iterator itsDirectPointed2;
            vector <double>::iterator itsFatherEnd;
            vector <double>::iterator itsFatherBegin;

            //int itsCoord;

            M1riterator() {};
            ~M1riterator() {};


            riterator & operator = (const riterator &itB)
            {
                itsDirectPointed2 = itB.itsDirectPointed2;
                itsFatherEnd      = itB.itsFatherEnd;
                itsFatherBegin      = itB.itsFatherBegin;

                //itsCoord          = itB.itsCoord;

                return *this;
            }

// FIXMEFRED: mig.vc8 (21/05/2007 10:50:52): no use
			//        inline riterator & operator = (const double* pointer)
    //        {
				//itsDirectPointed2 = static_cast<std::vector<double>::iterator>(const_cast<double*>(pointer));
    //            
    //            return *this;
    //        }


                        
/*
            riterator(vector <double>::iterator it)
            {
                itsDirectPointed2 = it;
            }
*/


            inline bool operator < (riterator &itB)
            {
                return itsDirectPointed2 < itB.itsDirectPointed2;
            }

            inline bool operator >  (riterator &itB)
            {
                return itsDirectPointed2 > itB.itsDirectPointed2;
            }

            inline bool operator == (riterator &itB)
            {
                return (itsDirectPointed2 == itB.itsDirectPointed2);
            }

/*
            inline double & operator ->(void)
            {
                return *itsDirectPointed2;
            }
*/

            inline double & operator *()
            {
                return *itsDirectPointed2;
            }

            //the dummy int variable enables to specify the postfix form of the operator
            inline riterator & operator ++(int i)
            {
                itsDirectPointed2++;

                return *this;
            }

            inline riterator & operator --(int i)
            {
                itsDirectPointed2--;

                return *this;
            }

            inline riterator operator +(int i)
            {
                riterator sum;

                sum = *this;
                sum.itsDirectPointed2 += i;

                return sum;
            }


            inline riterator operator -(int i)
            {
                riterator sum;

                sum = *this;
                sum.itsDirectPointed2 -= i;

                return sum;
            }


            inline riterator & operator +=(int i)
            {
                itsDirectPointed2 += i;

                return *this;
            }


            inline riterator & operator -=(int i)
            {
                itsDirectPointed2 -= i;

                return *this;
            }


            riterator & operator << (SlideMove &move)
            {

                if (move.nbDimsRemaining == 0)
                {
                    itsDirectPointed2 += move.nbStep;
                }
                else 
                {
            		//throw Exception(__LINE__, __FILE__, ERR_CONDITION_NOT_MEET,
                    //                    "You cannot go that far !");
                }

                return *this;
            }

            riterator & operator << (list<SlideMove> *moveList)
            {

                list<SlideMove>::iterator it;

                for (it = moveList->begin(); it != moveList->end(); it++)
                    *this = *this << *it;

                return *this;
            }


            inline int GetStep(void)
            {
                return itsDirectPointed2-itsFatherBegin;
            }

            TCoord GetCoord(void)
            {
                TCoord coord;

                coord[0] = itsDirectPointed2-itsFatherBegin;

                return coord;
            }

            void GetCoordInt(int *coord2fill = NULL)
            {
                coord2fill[0] = itsDirectPointed2-itsFatherBegin;
            }








        };







/****************************************************************************************/
/*                                                                                      */
/*  Operateurs                                                                          */
/*                                                                                      */
/****************************************************************************************/


/*********************************************************************/
/*                                                                   */
/*                  UGLY MACRO                                       */
/*  Unnecessary with VC++, but needed with Sun compiler              */
/*  (no template class <- member class )                             */
/*                                                                   */
/*********************************************************************/

#define MNArray(_n,_nMinusOne)                          \
                                                        \
struct TCoord##_n                                       \
{                                                       \
  int x[_n];                                            \
                                                        \
  inline int & operator[] (int k) {return x[k];}        \
  inline int & Elt        (int k) {return x[k];}        \
};                                                      \
                                                        \
                                                        \
template <class T> class M##_n##Array;                  \
                                                        \
template <class T>                                      \
    class M##_n##riterator  {                           \
                                                        \
        typedef M##_n##Array         < T > TcurrMArray; \
        typedef M##_nMinusOne##Array < T > TnextMArray; \
        typedef vector <T>         vectorT;             \
                                                        \
                                                        \
                                                        \
        typedef  TCoord##_n             TCoord;         \
        typedef  M##_n##riterator         riterator;    \
                                                        \
        public :                                        \
            CC_TYPENAME TnextMArray::riterator itsNextRaw;          \
            CC_TYPENAME vector < TnextMArray >::iterator itsDirectPointed2;     \
            CC_TYPENAME vector < TnextMArray >::iterator itsFatherEnd;          \
            CC_TYPENAME vector < TnextMArray >::iterator itsFatherBegin;        \
                                                                    \
                                                                    \
                                                                    \
                                                                    \
                                                                    \
            M##_n##riterator(void) {};                              \
            ~M##_n##riterator(void) {};                             \
                                                                    \
                                                                    \
                                                                    \
                                                                    \
            inline riterator & operator = (const riterator &itB)    \
            {                                                       \
                itsNextRaw        = itB.itsNextRaw;                 \
                itsDirectPointed2 = itB.itsDirectPointed2;          \
                itsFatherEnd      = itB.itsFatherEnd;               \
                itsFatherBegin    = itB.itsFatherBegin;             \
                                                                    \
                return *this;                                       \
            }                                                       \
                                                                    \
                                                                    \
            inline riterator & operator = (const T* pointer)        \
            {                                                       \
                itsNextRaw = pointer;                               \
                return *this;                                       \
            }                                                       \
                                                                    \
            inline T & operator ->(void)                            \
            {                                                       \
                return itsNextRaw.operator->();                     \
            }                                                       \
                                                                    \
                                                                    \
                                                                    \
            bool operator < (riterator &itB)                        \
            {                                                       \
                if (itsDirectPointed2 < itB.itsDirectPointed2)      \
                    return true;                                    \
                else                                                \
                {                                                   \
                    if (itsDirectPointed2 == itB.itsDirectPointed2) \
                        return (itsNextRaw < itB.itsNextRaw);       \
                    else                                            \
                        return false;                               \
                }                                                   \
            }                                                       \
                                                                    \
            inline bool operator > (riterator &itB);                \
                                                                    \
            bool operator == (riterator &itB)                       \
            {                                                       \
                if (itsDirectPointed2 == itB.itsDirectPointed2)     \
                {                                                   \
                    if (itsNextRaw == itB.itsNextRaw)               \
                        return true;                                \
                    else                                            \
                        return false;                               \
                }                                                   \
                else                                                \
                    return false;                                   \
                                                                    \
            }                                                       \
                                                                    \
                                                                    \
                                                                    \
                                                                    \
            inline  T & operator *()                                \
            {                                                       \
                return *itsNextRaw;                                 \
            }                                                       \
                                                                    \
                                                                    \
            inline riterator & operator ++(int i)                   \
            {                                                       \
                itsNextRaw++;                                       \
                if  (itsNextRaw.itsDirectPointed2 == itsDirectPointed2->end())  \
                {                                                   \
                    itsDirectPointed2++;                            \
                    if (itsDirectPointed2 < itsFatherEnd)           \
                    {                                               \
                        itsNextRaw = itsDirectPointed2->rbegin();   \
                    }                                               \
                }                                                   \
                                                                    \
                return *this;                                       \
            }                                                       \
                                                                    \
            inline riterator & operator --(int i)                   \
            {                                                       \
                itsNextRaw--;                                       \
                if  (itsNextRaw.itsDirectPointed2 < itsDirectPointed2->begin()) \
                {                                                   \
                    itsDirectPointed2--;                            \
                                                                    \
                    if (itsDirectPointed2 >= itsFatherBegin)        \
                    {                                               \
                        itsNextRaw = itsDirectPointed2->rbegin();   \
                    }                                               \
                }                                                   \
                                                                    \
                return *this;                                       \
            }                                                       \
                                                                    \
                                                                    \
                                                                    \
            inline riterator operator +(int i)                      \
            {                                                       \
                riterator sum;                                      \
                                                                    \
                                                                    \
                sum = *this;\
                                                                    \
                while (i>0)\
                {\
                    sum++;\
                    i--;\
                }\
                                                                    \
                return sum;\
            }\
                                                                    \
                                                                    \
            inline riterator operator -(int i)\
            {\
                riterator sum;\
                                                                    \
                sum = *this;\
                                                                    \
                while (i>0)\
                {\
                    sum--;\
                    i--;\
                }\
                                                                    \
                return sum;\
            }\
                                                                    \
                                                                    \
            inline riterator & operator +=(int i)\
            {\
                while (i>0)\
                {\
                    (*this)++;\
                    i--;\
                }\
                                                                    \
                return *this;\
            } \
                                                                    \
                                                                    \
            inline riterator & operator -=(int i)\
            {\
                while (i>0)\
                {\
                    (*this)--;\
                    i--;\
                }\
                                                                    \
                return *this;\
            }\
                                                                    \
                                                                    \
                                                                    \
                                                                    \
                                                                    \
                                                                    \
            void slide2(int nbStep, int nbDimsRemaining = 0); \
                                                                    \
                                                                    \
            riterator & operator << (SlideMove &move) \
            { \
                int d = itsNextRaw.itsDirectPointed2-itsDirectPointed2->GetVector()->begin(); \
                                                                    \
                if (move.nbDimsRemaining == 0) \
                { \
                    itsDirectPointed2 += move.nbStep; \
                    itsNextRaw = itsDirectPointed2->rbegin() + d; \
                                                                    \
                }\
                else \
                {\
                    move.nbDimsRemaining--;\
                    itsNextRaw << move;\
                    move.nbDimsRemaining++;\
                }\
                                                                    \
                return *this;\
            }\
                                                                    \
                                                                    \
            riterator & operator << (list<SlideMove> *moveList)\
            {\
                CC_TYPENAME list<SlideMove>::iterator it;\
                                                                    \
                for (it = moveList->begin(); it != moveList->end(); it++)\
                    *this << *it; \
                                                                    \
                return *this;\
            }\
                                                                    \
            inline int GetStep(void)\
            {\
                return itsDirectPointed2-itsFatherBegin;\
            }\
                                                                    \
            TCoord GetCoord(void)\
            {\
                TCoord coord;\
                                                                    \
                coord[0] = itsDirectPointed2-itsFatherBegin;\
                itsNextRaw.GetCoordInt(coord.x+1);\
                                                                    \
                return coord;\
            }\
                                                                    \
            void GetCoordInt(int *coord2fill = NULL)\
            {\
                coord2fill[0] = itsDirectPointed2-itsFatherBegin;\
                itsNextRaw.GetCoord(coord2fill+1);\
            }\
                                                                    \
                                                                    \
                                                                    \
                                                                    \
                                                                    \
        };\
                                                                    \
    template <class T> class M##_n##Array\
    {\
    public :\
        typedef M##_n##Array < T > TcurrMArray;\
        typedef M##_nMinusOne##Array < T > TnextMArray;\
        typedef vector <T>         vectorT;\
                                                                    \
        typedef TCoord##_n   TCoord;\
                                                                    \
    protected :\
                                                                    \
        int itsLength;\
        int itsOrder;\
                                                                    \
        vector < TnextMArray> itsVector;\
                                                                    \
                                                                    \
                                                                    \
    public :\
                                                                    \
    typedef M##_n##riterator<T> riterator;\
                                                                    \
                                                                    \
        M##_n##Array()\
        {\
            itsOrder = _n;\
            itsLength = 0;\
        }\
                                                                    \
    ~M##_n##Array() {}\
                                                                    \
        inline riterator rbegin()\
        {\
            riterator theBegin;\
                                                                    \
            theBegin.itsDirectPointed2 = itsVector.begin();\
            theBegin.itsNextRaw = theBegin.itsDirectPointed2->rbegin();\
            theBegin.itsFatherEnd = itsVector.end();\
            theBegin.itsFatherBegin = theBegin.itsDirectPointed2;\
                                                                    \
                                                                    \
                                                                    \
            return theBegin;\
        }\
                                                                    \
        inline riterator rend()\
        {\
            riterator theEnd;\
                                                                    \
            theEnd.itsDirectPointed2 = itsVector.end()-1;\
            theEnd.itsNextRaw = theEnd.itsDirectPointed2->rend();\
            theEnd.itsFatherEnd = itsVector.end();\
            theEnd.itsFatherBegin = itsVector.begin();\
                                                                    \
                                                                    \
            return theEnd;\
        }\
                                                                    \
                                                                    \
                                                                    \
                                                                    \
        inline CC_TYPENAME vector < TnextMArray >::iterator begin()\
        {\
            return itsVector.begin();\
        }\
                                                                    \
        inline CC_TYPENAME vector < TnextMArray >::iterator end()\
        {\
            return itsVector.end();\
        }\
                                                                    \
                                                                    \
                                                                    \
        inline TnextMArray & operator[] (int i)\
        {\
            return itsVector[i];\
        }\
                                                                    \
        inline T & operator[] (riterator &it)\
        {\
            return *it;\
        }\
                                                                    \
                                                                    \
                                                                    \
        inline TnextMArray & at (int i)\
        {\
            return itsVector[i];\
        }\
                                                                    \
        inline T & rat (int *coord)\
        {\
            return itsVector[coord[0]].rat(coord+1);\
        }\
                                                                    \
        riterator GetRiterator(TCoord &coord)\
        {\
            riterator outIt;\
                                                                    \
            outIt.itsDirectPointed2 = itsVector.begin() + coord.x[0];\
                                                                    \
                                                                    \
            outIt.itsNextRaw = outIt.itsDirectPointed2->GetRiteratorInt(&(coord.x[0])+1);\
            outIt.itsFatherEnd = itsVector.end();\
            outIt.itsFatherBegin = itsVector.begin();\
                                                                    \
            return outIt;\
        }\
                                                                    \
                                                                    \
        riterator GetRiteratorInt(int *coord)\
        {\
            riterator outIt;\
                                                                    \
            outIt.itsDirectPointed2 = itsVector.begin() + coord[0];\
                                                                    \
                                                                    \
            outIt.itsNextRaw = outIt.itsDirectPointed2->GetRiteratorInt(coord+1);\
            outIt.itsFatherEnd = itsVector.end();\
            outIt.itsFatherEnd = itsVector.begin();\
                                                                    \
            return outIt;\
        }\
                                                                    \
                                                                    \
        M##_n##Array<T> & operator = (TcurrMArray& inMArray)\
        {\
            CC_TYPENAME TcurrMArray::riterator rit, ritIn, finish;\
                                                                    \
            SetLength(inMArray.GetLength());\
            finish = this->rend();\
            for (rit = rbegin(), ritIn = inMArray.rbegin(); rit < finish; rit++)\
                *rit = *ritIn;\
                                                                    \
            return *this;\
        }\
                                                                    \
                                                                    \
                                                                    \
                                                                    \
        void reset()\
        {\
            CC_TYPENAME vector < TnextMArray  >::iterator next;\
                                                                    \
            forall(itsVector, next)\
                next->reset();\
                                                                    \
            itsVector.clear();\
        }\
                                                                    \
                                                                    \
        inline void SetCurrentLength(int l)\
        {\
            itsLength = l;\
            itsVector.resize(itsLength);\
        }\
                                                                    \
        inline int GetCurrentLength(void)\
        {\
            return itsLength;\
        }\
                                                                    \
                                                                    \
                                                                    \
        bool SetLength(int *ilengths)\
        {\
            list <int> lLength(ilengths, ilengths+itsOrder);\
                                                                    \
            return SetLength(&lLength);\
        }\
                                                                    \
                                                                    \
                                                                    \
                                                                    \
        bool SetLength(list<int> *length)\
        {\
            CC_TYPENAME vector < TnextMArray >::iterator next;\
            bool OK = true;\
            int tmp;\
                                                                    \
            if (length->size() != _n)\
                return false;\
                                                                    \
            itsLength = length->front();\
            itsVector.resize(itsLength);\
            tmp = *(length->begin());\
            length->pop_front();\
                                                                    \
            forall(itsVector, next)\
                OK = OK && next->SetLength(length);\
            length->push_front(tmp);\
                                                                    \
            if (!OK)\
                reset();\
                                                                    \
            return OK;\
        }\
                                                                    \
        \
                                                                    \
                                                                    \
        void RetrieveLength(list<int> *length)\
        {\
            itsVector[0].RetrieveLength(length);\
            length->push_front(itsLength);\
        }\
                                                                    \
        list<int> *GetLength(void)\
        {\
            list<int> *length = new list <int>;\
            RetrieveLength(length);\
                                                                    \
            return length;\
        }\
                                                                    \
                                                                    \
                                                                    \
                                                                    \
        void print(int skip = 0)\
        {\
            int i=0,j;\
            CC_TYPENAME vector <TnextMArray>::iterator it;\
                                                                    \
            if (itsLength >0)\
            {\
                it = itsVector.begin();\
                cout << "[" << i << "]";\
                    it->print(skip+1);\
                i++;\
            }\
                                                                    \
            for(it = itsVector.begin()+1; it <itsVector.end(); it++)\
            {\
                for(j=0; j<skip; j++)\
                    cout << "   ";\
                cout << "[" << i << "]";\
                it->print(skip+1);\
                i++;\
            }\
        }\
                                                                    \
                                                                    \
        inline vector < TnextMArray> * GetVector(void)\
        {\
            return &itsVector;\
        }\
                                                                    \
                                                                    \
        void Alloc()\
        {\
                                                                    \
            CC_TYPENAME vectorT::iterator it;\
                                                                    \
            forall(itsVector, it)\
                *it = new T;\
        }\
                                                                    \
        void Alloc(int size)\
        {\
                                                                    \
            CC_TYPENAME vectorT::iterator it;\
                                                                    \
            forall(itsVector, it)\
                *it = new T(size);\
        }\
                                                                    \
                                                                    \
                                                                    \
                                                                    \
        void Cleanit(void)\
        {\
            riterator theEnd = rend();\
            riterator it;\
                                                                    \
            for (it = this->rbegin(); it < theEnd; it++)\
                cleanit(*it);\
        }\
                                                                    \
    };



MNArray(3,2);
MNArray(4,3);


#endif
