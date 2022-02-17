/*
 * $Log: containr.cpp,v $
 * Revision 1.2  2003/07/31 16:37:00  mab
 * RCS comment Added
 *
 */

/*----------------------------------------------------------------------------*/
/*                                                                            */
/* FILE        : containr.cpp                                                 */
/*                                                                            */
/* DESCRIPTION : Container class Methods File                                 */
/*                                                                            */
/* DATE        : Wed Dec  4 11:46:52 MET 1996                                 */
/*                                                                            */
/*----------------------------------------------------------------------------*/



#include "containr.h"





ARM_Container::ARM_Container(const ARM_Container& container) :
                          ARM_Object(container)
{
    Init();

    BitwiseCopy(&container);
}



ARM_Container& ARM_Container::operator = (const ARM_Container& container)
{
    (*this).ARM_Object::operator = (container);

    BitwiseCopy(&container);

    return(*this);
}



void ARM_Container::BitwiseCopy(const ARM_Object* srcCont)
{
    ARM_Container* container = (ARM_Container *) srcCont;


    if ( container == NULL )
       return;

    FreeLinks();

    int sz = container->GetSize();


    ARM_Object* obj;
 
    obj = container->Start();

    while (obj) 
    {
        Append(obj);

        obj = container->Next();
    }
}



ARM_Container* ARM_Container::DuplicateWithContent(void)
{
    ARM_Container* duplicata = new ARM_Container();  
 
 
    ARM_Object* obj;
    ARM_Object* objClone;
 
 
    obj = this->Start();
 
    while (obj) 
    {
        objClone = obj->Clone();
 
        duplicata->Append(objClone);
 
        obj = this->Next();
    }

    return(duplicata);
}



// Insert at last 

void ARM_Container::Append(ARM_Object* obj)
{
    ARM_Link* newLink;


    newLink = new ARM_Link(obj);

    if ( last == (ARM_Link *) NULL )
    {
       last = newLink;

       last->_next = last;
    }
    else
    {
       newLink->_next = last->_next;       

       last->_next = newLink;

       last = newLink;  
    } 

    size++;
}



void ARM_Container::SetItem(ARM_Object* obj, int i)
{
    ARM_Link* nthItem = (ARM_Link *) NULL;


    try
    {
       nthItem = GetNthLink(i);
    }

    catch(Exception& x)
    {
        x.DebugPrint();
 
        throw x;
    }

    nthItem->item = obj;
}










/*----------------------------------------------------------------------------*/
/*---- End of file ----*/
