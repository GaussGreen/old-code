/****************************************************************
 * Module:	VirtualProduct
 * Submodule:	
 * File:	vtlfleg3Idx.h
 * Function:	
 * Author:	Robert Mattingly
 *****************************************************************/
#ifndef	_vtlfleg3Idx_H
#define	_vtlfleg3Idx_H

#include "vtlbase.h"
#include "vtlfleg.h"
#include "vpfleg3Idx.h"


//--------------------------------------------------------------
/**
 * Class for floating leg with 3 rate indices.
 */

class KVPToolFloatLeg3Idx : public KVPToolFloatLeg {
 public:
  /**
   * Constructor.
   */
  KVPToolFloatLeg3Idx(
		      SharedPointer<KVPFloatLeg3Idx> ins, // (I) inetrument to value
		      KVTree &vt);			    // (I) virtual tree

  /**
   * Destructor.
   */
  virtual	~KVPToolFloatLeg3Idx();


  /**
   * Type name.
   */
  virtual	const char*	TypeName() const {return("KVPToolFloatLeg3Idx");};

  /**
   * Returns the name for the geometry of the slice
   * (see description of method in KVPToolAtom).
   */
  virtual	const String&	GetCurveName() { return mCurveName;}


  /**
   * Returns the KVPAtom to which the KVPToolAtom is associated.
   */
  virtual	SharedPointer<KVPAtom>	Atom()
    {
      SharedPointer<KVPAtom> vp;
      SharedPointerConvertTo(mFloatLeg3Idx, vp);
      return (vp);
    }

  /**
   * Update tool.
   */
  virtual	void	Update();


  /**
   * Returns the current value of the stream.
   */
  virtual	KTSlice&	GetValue();



 private:

  /** Pointer to instrument. */
  SharedPointer<KVPFloatLeg3Idx> mFloatLeg3Idx;
 
  KVector(TDate) mModelResetDates2;
  KVector(TDate) mModelResetDates3;

};

#endif  /* _vtlfleg3Idx_H */



