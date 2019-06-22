/* File: progressbar.c */
/*
  This file is a part of the Corrfunc package
  Copyright (C) 2015-- Manodeep Sinha (manodeep@gmail.com)
  License: MIT LICENSE. See LICENSE file under the top-level
  directory at https://github.com/manodeep/Corrfunc/
*/

#include "progressbar.h"
#include "utils.h"
#include <inttypes.h>

#ifndef MAXLEN
#define MAXLEN 1000
#endif


static double SMALLPRINTSTEP;
static char PROGRESSBARSTRING[MAXLEN];
static int beg_of_string_index;
static double percent=0;
static int END_INDEX_FOR_PERCENT_DONE[101];
struct timeval tstart;

void init_my_progressbar(const int64_t N,int *interrupted)
{
    if(N <= 0) {
        fprintf(stderr,"WARNING: N=%"PRId64" is not positive. Progress bar will not be printed\n",N);
        SMALLPRINTSTEP = 0.0;
    } else {
        //set the increment
        SMALLPRINTSTEP = 0.01 * N;

        //pre-fill the progress bar string
        //first the 0%
        int index=0;
        my_snprintf(&(PROGRESSBARSTRING[index]),MAXLEN-index,"%s","0%");
        index+=2;
        END_INDEX_FOR_PERCENT_DONE[0] = index;
        for(int i=1;i<100;i++) {
            if(i%10 == 0) {
                my_snprintf(&(PROGRESSBARSTRING[index]),MAXLEN-index,"%02d%%",i);
                index += 3;
            } else {
                PROGRESSBARSTRING[index++] = '.';
            }
            END_INDEX_FOR_PERCENT_DONE[i] = index;
        }

        //end with 100%
        my_snprintf(&(PROGRESSBARSTRING[index]),MAXLEN-index,"%s","100%");
        index+=4;
        END_INDEX_FOR_PERCENT_DONE[100] = index;
        PROGRESSBARSTRING[index+1] = '\0';

    }
    *interrupted=0;
    beg_of_string_index=0;
    gettimeofday(&tstart, NULL);
}

void my_progressbar(const int64_t curr_index,int *interrupted)
{
    if(SMALLPRINTSTEP > 0.0 )   {
      if(*interrupted == 1) {
        fprintf(stderr,"\n");
        *interrupted = 0;
        beg_of_string_index = 0;
      }

      percent = (curr_index+1)/SMALLPRINTSTEP;//division is in double -> C has 0-based indexing -- the +1 accounts for that.
      int integer_percent = (int) percent;

      if (integer_percent >= 0 && integer_percent <= 100)  {
        /*        for(int i=beg_of_string_index;i<END_INDEX_FOR_PERCENT_DONE[integer_percent];i++) */
        /*            fprintf(stderr,"%c",PROGRESSBARSTRING[i]); */

        fprintf(stderr,"%.*s",END_INDEX_FOR_PERCENT_DONE[integer_percent]-beg_of_string_index,&(PROGRESSBARSTRING[beg_of_string_index]));
        beg_of_string_index = END_INDEX_FOR_PERCENT_DONE[integer_percent];
      }
    }
}

void finish_myprogressbar(int *interrupted)
{
  struct timeval t1;
  gettimeofday(&t1, NULL);
  char * time_string = get_time_string(tstart, t1);
  if(SMALLPRINTSTEP > 0.0) {
    if(*interrupted == 0)   {
      if(percent < 100.0) {
        fprintf(stderr,"100%% done.");
      } else {
        fprintf(stderr," done.");
      }
    } else {
      fprintf(stderr,"\n%s done.",PROGRESSBARSTRING);
      *interrupted=0;
    }
  } else {
    fprintf(stderr," done.");
  }
  fprintf(stderr," Time taken = %s\n", time_string);
  free(time_string);

  beg_of_string_index = 0;
  my_snprintf(PROGRESSBARSTRING,MAXLEN,"%s"," ");
}
