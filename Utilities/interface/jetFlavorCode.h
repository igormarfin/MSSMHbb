#ifndef jetFlavorCode_h
#define jetFlavorCode_h

#include "Analysis/Utilities/interface/HbbNtuple.h"

int jetFlavorCode(const int iJet) {
  int theFlavCode = -1;
  int theFlavor = partonFlavorJet[iJet];
  if (hflContentJet[iJet] != 0) theFlavor = hflContentJet[iJet];
  switch (abs(theFlavor)) {
  case 0:
  case 1:
  case 2:
  case 3:
    theFlavCode = 0;
    break;
  case 4:
    theFlavCode = 1;
    break;
  case 5:
    theFlavCode = 2;
    break;
  case 21:
    theFlavCode = 0;
    break;
  default:
    std::cout << "bad flavor " << theFlavor << std::endl;
  }
  return theFlavCode;
}




int jetFlavorCodeNewCat(const int iJet) {
 
  int iflav = -1;
  
  const int parton_matching_mode = 3;
 
  if( (int)BContentJet3[iJet] > 0 )
    {
      if(parton_matching_mode == 3)
        {
          if((int)BContentJet3[iJet] == 1)
            iflav = 2;
          else 
            iflav = 4;  
        }
      else
        iflav = 2;
    }
  else if( (int)CContentJet3[iJet] > 0 )
    {
      if(parton_matching_mode == 3)
        {
          if( (int)CContentJet3[iJet] == 1 )
            iflav = 1;
          else
            iflav = 3; 
        }
      else
        iflav = 1;
    }
  else
    {
      //everything else is assumed to be udsg ... 
      switch (abs(partonFlavorJet[iJet])) {
      case 0: //correct?? not used in readmistag?
      case 1:
      case 2:
      case 3:
      case 4:
      case 5:
      case 21:
        iflav = 0;
        break;
      }
   /*    if(skip_undefined_jets && abs(partonFlavorJet[iJet]) == 0) */
/*         iflav = -1; */
    }
 
	if(iflav < -1 || iflav > 4) {
	std::cout << "error: wrong flav : " << iflav << std::endl;
	throw std::exception();
	}
  return iflav;
}







//this methods provides similar flavour definitions as used in the 2011 analysis
int jetFlavorCodeNewCat2011flavdef(const int iJet) {
 
  int iflav = -1;
  
  const int parton_matching_mode = 2;
 
  if( (int)BContentJet3[iJet] > 0 )
    {
      if(parton_matching_mode == 3)
        {
          if((int)BContentJet3[iJet] == 1)
            iflav = 2;
          else 
            iflav = 4;  
        }
      else
        iflav = 2;
    }
  else if( (int)CContentJet3[iJet] > 0 )
    {
      if(parton_matching_mode == 3)
        {
          if( (int)CContentJet3[iJet] == 1 )
            iflav = 1;
          else
            iflav = 3; 
        }
      else
        iflav = 1;
    }
  else
    {
      //everything else is assumed to be udsg ... 
      switch (abs(partonFlavorJet[iJet])) {
      case 0: //correct?? not used in readmistag?
      case 1:
      case 2:
      case 3:
      case 4:
      case 5:
      case 21:
        iflav = 0;
        break;
      }
   /*    if(skip_undefined_jets && abs(partonFlavorJet[iJet]) == 0) */
/*         iflav = -1; */
    }
 

  return iflav;
}




//matching code which uses a cone of size 0.5 instead of 0.3 which is the default
int jetFlavorCodeNewCat5(const int iJet) {
 
  int iflav = -1;
  
  const int parton_matching_mode = 3;
 
  if( (int)BContentJet5[iJet] > 0 )
    {
      if(parton_matching_mode == 3)
        {
          if((int)BContentJet5[iJet] == 1)
            iflav = 2;
          else 
            iflav = 4;  
        }
      else
        iflav = 2;
    }
  else if( (int)CContentJet5[iJet] > 0 )
    {
      if(parton_matching_mode == 3)
        {
          if( (int)CContentJet5[iJet] == 1 )
            iflav = 1;
          else
            iflav = 3; 
        }
      else
        iflav = 1;
    }
  else
    {
      //everything else is assumed to be udsg ... 
      switch (abs(partonFlavorJet[iJet])) {
      case 0: //correct?? not used in readmistag?
      case 1:
      case 2:
      case 3:
      case 4:
      case 5:
      case 21:
        iflav = 0;
        break;
      }
   /*    if(skip_undefined_jets && abs(partonFlavorJet[iJet]) == 0) */
/*         iflav = -1; */
    }
 

  return iflav;
}






int diJetFlavorCode(const int iJet0, const int iJet1) {
  // determine dijet flavor from MC info
  // we set a flavor pair code that does not distinguish ordering
  int theFlav[2];
//  theFlav[0] = jetFlavorCode( iJet0 );
//  theFlav[1] = jetFlavorCode( iJet1 );
  theFlav[0] = jetFlavorCodeNewCat( iJet0 );
  theFlav[1] = jetFlavorCodeNewCat( iJet1 );
    if(theFlav[0] < -1 ||  theFlav[0]>4 || theFlav[1] < -1 ||  theFlav[1]>4 ) {
    std::cout << "error diJetFlavorCode: wrong flav : " << theFlav[0]<<"   "<<theFlav[1] << std::endl;
    throw std::exception();
    }

  // now set the dijet flavor
  int theFcDijet = -1;
  int flavPatternDijet = 10*theFlav[0]+theFlav[1];
  switch (flavPatternDijet) {
  case 0: ///qq			(1)
    theFcDijet = 5;		
    break;
  case 11: ///c1c1		(2)
    theFcDijet = 3;
    break;
  case 22:///b1b1		(3)
    theFcDijet = 0;
    break;
  case 1://qc1 			(4)
  case 10:///c1q
    theFcDijet = 4;
    break;
  case 2:///qb1
  case 20:///b1q		(5)
    theFcDijet = 2;
    break;
  case 12:///c1b1		(6)	
  case 21:///b1c1
    theFcDijet = 1;
    break;

  case 32: ///c2b1		(7)
  case 23: ///b1c2
		theFcDijet = 6;
		break;

  case 31: ///c2c1		(8)
  case 13: ///c1c2
		theFcDijet = 7;
		break;

  case 30: ///c2q		(9)
  case 3: ///qc2
		theFcDijet = 8;
		break;

  case 33: ///c2c2		(10)
		theFcDijet = 9;
		break;


  case 42: ///b2b1		(11)
  case 24: ///b1b2
		theFcDijet = 10;
		break;

  case 41: ///b2c1		(12)
  case 14: ///c1b2
		theFcDijet = 11;
		break;

  case 40: ///b2q		(13)
  case 4: ///qb2
		theFcDijet = 12;
		break;

  case 44: ///b2b2		(14)
		theFcDijet = 13;
		break;


  case 43: ///b2c2  (15)
  case 34: ///c2b2
        theFcDijet = 14;
        break;


  default:
	 if(theFlav[0] < -1 ||  theFlav[0]>4 || theFlav[1] < -1 ||  theFlav[1]>4 ) {
    std::cout << "diJetFlavorCode: Bad flavor code " << theFlav[0] << " " << theFlav[1] << std::endl;
    throw std::exception();
	}
  }
  return theFcDijet;
}

#endif  // ifdef
