## Resubmission
This is a resubmission. 
The previous problem with calling unexported from mice has been resolved
- after talking to the maintainer Stef van Buuren - by working around 
using the padModel-function on the hand and copying the other two
needed functions into an internal script on the other hand.
Software names in the DESCRIPTION file are now in single quotes.

## Test environments
* local Windows 10, R 3.3.3 and R 3.4.3, both 32/64bit

## R CMD check results
There were no ERRORs, WARNINGs or NOTEs.
  
## Downstream dependencies
There are no downstream dependencies for this package, as this is a new release.