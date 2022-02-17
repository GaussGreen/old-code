/**********************************************************************
*
*	Header %name:	CalypsoOptionSchedule.h %
*	Instance:		1
*	Description:	
*	%created_by:	fdegraeve %
*	%date_created:	Wed May 30 17:26:20 2007 %
*
**********************************************************************/
#ifndef _CALYPSOOPTIONSCHEDULE_H
#define _CALYPSOOPTIONSCHEDULE_H


#include <refvalue.h>
#include <exercise.h>
#include <util\fromto.h>
#include <ARM\libarm_frometk\arm_local_parsexml_util.h> 
#include "XMLTools.h"

class ARM_Date;
class ARM_ReferenceValue;

ARM_ReferenceValue* convertCalypsoOptionSchedule(MSXML2::IXMLDOMNodePtr xmlNode);
void convertCalypsoOptionSchedule(MSXML2::IXMLDOMNodePtr xmlNode, ARM_ReferenceValue*& fees, ARM_ExerciseStyle*& exerStyle);


#endif
