/*
 *   segemehl - a read aligner
 *   Copyright (C) 2008-2017  Steve Hoffmann and Christian Otto
 *
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */




/*
 *  info.c
 *  nfo messages
 *
 *  @author Steve Hoffmann
 *  @email steve@bioinf.uni-leipzig.de
 *  @date 08/26/2007 06:49:02 PM CEST
 *  
 *  SVN
 *  Revision of last commit: $Rev: 19 $
 *  Author: $Author: steve $
 *  Date: $Date: 2008-05-14 15:43:29 +0200 (Wed, 14 May 2008) $
 *
 *  Id: $Id: info.c 19 2008-05-14 13:43:29Z steve $
 *  Url: $URL: file:///homes/bierdepot/steve/svn/segemehl/trunk/libs/info.c $
 *  
 */
 
 #include <stdarg.h>
 #include <stdio.h>
 #include <stdlib.h>
 #include <string.h>
 #include <time.h>
 #include "info.h"
 #include "debug.h"

 FILE *nfodevice = NULL;


 char *timestr_r(const struct tm *timeptr) {
    static const char wday_name[7][3] = {
        "Sun", "Mon", "Tue", "Wed",
        "Thu", "Fri", "Sat"
    };

    static const char mon_name[12][3] = {
        "Jan", "Feb", "Mar", "Apr", "May", "Jun",
        "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"
    };

    static char result[26];

    sprintf(result, "%.3s %.3s%3d %.2d:%.2d:%.2d %d",
        wday_name[timeptr->tm_wday], mon_name[timeptr->tm_mon],
        timeptr->tm_mday, timeptr->tm_hour, timeptr->tm_min,
        timeptr->tm_sec, 1900 + timeptr->tm_year);

    return result;
 }

 int
 infomsg( char *file, 
          int line, 
          const char *fmt, ...) {

   int ret;
   va_list ap;
   time_t rawtime;
   struct tm *timeinfo;
  
   if (mute) return 0;

    time(&rawtime);
    timeinfo = localtime (&rawtime);

   if (nfodevice == NULL) {
     nfodevice = NFODEVICE;
   }
   
   va_start(ap, fmt);
#ifdef PROGNFO   
   fprintf(nfodevice, "[%s] %s: ", "HYDI", timestr_r(timeinfo));
#endif
   ret = vfprintf(nfodevice, fmt, ap);
   va_end(ap);

   return ret; 
 }


void 
setnfodevice(char *filename) {
  FILE *fp;

  fp = fopen(filename, "w");
  if (fp == NULL) {
    DBG("Couldn't open file '%s'. Exit forced.\n", filename);
    exit(-1);
  }

  nfodevice = fp;
}

int
nfolevel( char *file,
    int line,
    int level,
    const char *fmt, ... ) {

  int ret=0;
  va_list ap;
  time_t rawtime;
  struct tm *timeinfo;
   
  if (mute) return 0;
  
  time(&rawtime);
  timeinfo = localtime (&rawtime);

   if (nfodevice == NULL) {
     nfodevice = NFODEVICE;
   }
  
   if (NFOLEVEL >= level) {

    va_start(ap, fmt);
#ifdef PROGNFO
    fprintf(nfodevice, "[%s] %s: ", "SEGEMEHL", timestr_r(timeinfo));
#endif
    ret = vfprintf(nfodevice, fmt, ap);
    va_end(ap);
  }

  return ret;
}


