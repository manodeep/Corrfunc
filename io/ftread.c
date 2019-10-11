/* File: ftread.c */
/*
  This file is a part of the Corrfunc package
  Copyright (C) 2015-- Manodeep Sinha (manodeep@gmail.com)
  License: MIT LICENSE. See LICENSE file under the top-level
  directory at https://github.com/manodeep/Corrfunc/
*/

/*  ftread reads unformatted data that has been written by fortran
    or in fortran convention -- i.e. an integer with the number of
    bytes in the record, the record data, and another integer with
    the number of bytes in the record.  ftread does various error
    checks, writing a message to stderr and returning a negative
    value if an error or warning occurs.  The call is identical to
    the standard i/o library routine fread.
*/

#include "ftread.h"


int my_ftread(void *ptr,size_t size,size_t nitems,FILE *stream)
{
    int err = ftread(ptr,size,nitems,stream);
    if(err != 0){
        fprintf(stderr,"ERROR: In %s> Read error in ftread while trying to read nitems = %zu "
                "of %zu bytes each\n", __FUNCTION__, nitems, size);
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}


int ftread(void *ptr,size_t size,size_t nitems,FILE *stream)
{
    unsigned int nbyte1, nbyte2;
    size_t nitem1 ;
    int errno ;
    errno = 0 ;
    if ( fread(&nbyte1,sizeof(nbyte1),1,stream) != 1 )  {
        errno = -10 ;
        fprintf(stderr,"read error, file empty ? \n") ;
    }
    nitem1 = fread(ptr,size,nitems,stream) ;
    if ( nitem1 != nitems ) {
        errno = -20 ;
        fprintf(stderr,"read error, %zu items expected, %zu items read. \n",
                nitems,nitem1) ;
    }
    if ( fread(&nbyte2,sizeof(nbyte2),1,stream) != 1 )  {
        errno = -30 ;
        fprintf(stderr,"read error, file too short ? \n") ;
    }
    if ( nbyte1 != nbyte2 ) {
        errno = errno - 1 ;
        fprintf(stderr,"read warning, byte # do not match, nbyte1 = %u, nbyte2 = %u \n",
                nbyte1,nbyte2) ;
    }
    if ( nbyte1 != size*nitems) {
        errno = errno - 2 ;
        fprintf(stderr, "read warning, byte # does not match item #, nbyte1 = %u, nitems = %zu \n",
                nbyte1,nitems) ;
    }

    return(errno) ;

}
