#ifndef HbbSelect_h
#define HbbSelect_h

   // variables in selection tree
   int selRun;    // run number for this event
   int selEvent;  // event number for this event
   int selPass;   // 1 if event passed MVA, 0 otherwise 

   double mjj; /// to validate bias in M12 from MVA selection
#endif
