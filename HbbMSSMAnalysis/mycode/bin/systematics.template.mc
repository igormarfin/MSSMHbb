##The file of systematics 
##Format:
##Type_of_systematics:nSigma
##Example for JES +2sigma
##0:2

## lines with '##' are comments
## lines with '#@' are used to create folders

## lines with '#%' are used as a line to fill the file 'systematic'

##Scale Factor of Btag
#@SF
#%1:0.0

##SF b/c + 2sigma
#@SFbc_SysUp
#%2:2.0

##SF b/c- 2sigma
#@SFbc_SysDown
#%2:-2.0

##SF udsg + 2sigma
#@SFudsg_SysUp
#%3:2.0

##SF udsg- 2sigma
#@SFudsg_SysDown
#%3:-2.0


##JES + 2sigma
#@JES_SysUp
#%0:2.0

##JES - 2sigma
#@JES_SysDown
#%0:-2.0

