#define _DEFAULT_SOURCE

#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <string.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <pthread.h>
#include "filebuffer.h"
#include "stringutils.h"
#include "basic-types.h"
#include "memory.h"
#include "mathematics.h"


unsigned char mute=0;
char *ntcode;

int 
main(int argc, char** argv) {
  unsigned int i;
  char c;
  circbuffer_t mybuffer;
  pthread_mutex_t mtx;
  FILE *fp;
  uint32_t len;
  
 
  fprintf(stderr, "opening file.\n");
  fp = fopen("test/5N_R1_1Mg7aqj2", "r");
  fprintf(stderr, "initializing filebuffer.\n");
  bl_circBufferInit(&mybuffer, 10000000, fp, NULL);
  fprintf(stderr, "getting first line.\n");
  

  
//    fprintf(stderr, "init.\n");
    // int n = fread(mybuffer.buffer, sizeof(char), mybuffer.size, mybuffer.dev);  
    //last char
    //mybuffer.end += n-1;
    char *str = NULL;

//  char *str = bl_circBufferReadLine(&mybuffer);
  //fprintf(stderr, "'%s' (len:%d)\n", str, strlen(str));
  //FREEMEMORY(NULL, str);

  
  //fprintf(stderr, "retrieving next 10 lines.\n");
  //for(i=0; i < 100; i++) {
  
  while ((str = bl_circBufferReadLine(&mybuffer, &len)) != NULL) {
    //fprintf(stderr, "*****************************************************************************\n");

    //str = bl_circBufferReadLine(&mybuffer);
    fprintf(stdout, "%s (%u;%lu)\n", str, len, strlen(str));
    assert(strlen(str) == len);

    //fprintf(stderr, "#################### \n '%s'\n ################### \n", str);
    FREEMEMORY(NULL, str);
  }
  
  fclose(fp);

  bl_circBufferDestruct(&mybuffer);
  exit(-1);
  
  fp = fopen("test/filebuffertest.txt", "wb");
  char mystring[ ] = "This is a test string with a tab\tand a line break\n";



  pthread_mutex_init(&mtx, NULL); 
  bl_circBufferInit(&mybuffer, 198, fp, &mtx);

  fprintf(mybuffer.dev, "now adding byte by byte\n");
  for(i=0; i < 50; i++) {
    bl_circBufferSaveAddByte(&mybuffer, mystring[(i%strlen(mystring))]);
  }
 
  fprintf(mybuffer.dev, "last char added '%c'\n",  mystring[((i-1)%strlen(mystring))]);

  for(i=0; i < 49; i++) {
    c= bl_circBufferRead(&mybuffer);
    if(c) fprintf(mybuffer.dev, "[%d,%d]='%c'\n", mybuffer.beg, mybuffer.end, c);
  }
    
  fprintf(mybuffer.dev, "now adding byte by byte\n");
  for(i=0; i < 25; i++) {
    bl_circBufferSaveAddByte(&mybuffer, mystring[(i%strlen(mystring))]);
  }

  bl_circBufferSaveAddByte(&mybuffer, '\n');
  
  fprintf(mybuffer.dev, "now adding data normally starting at [%d,%d]\n", mybuffer.beg, mybuffer.end);
  for(i=0; i < 3333; i++) {
    bl_circBufferAddSave(&mybuffer, mystring, strlen(mystring));
  }

  fprintf(mybuffer.dev, "empty buffer\n");
  bl_circBufferEmpty(&mybuffer);
  fclose(fp);

  bl_circBufferDestruct(&mybuffer);

/*
      if(info.bufferedwrite) {
        info.sambuffer = bl_circBufferInitArray(info.threadno, 1000000, 
            info.dev, info.mtx3);
        info.snglbuffer = bl_circBufferInitArray(info.threadno, 1000000, 
            info.singlesplitdev, info.mtx6);
        info.multbuffer = bl_circBufferInitArray(info.threadno, 1000000, 
            info.multisplitdev, info.mtx7);
        info.trnsbuffer = bl_circBufferInitArray(info.threadno, 1000000, 
            info.transsplitdev, info.mtx8);
      }



      if(info.bufferedwrite) { 
        bl_circBufferEmptyArray(info.sambuffer, info.threadno); 
        bl_circBufferEmptyArray(info.snglbuffer, info.threadno); 
        bl_circBufferEmptyArray(info.multbuffer, info.threadno); 
        bl_circBufferEmptyArray(info.trnsbuffer, info.threadno); 

        bl_circBufferDestructArray(info.sambuffer, info.threadno); 
        bl_circBufferDestructArray(info.snglbuffer, info.threadno); 
        bl_circBufferDestructArray(info.multbuffer, info.threadno); 
        bl_circBufferDestructArray(info.trnsbuffer, info.threadno); 

        FREEMEMORY(space, info.sambuffer);
        FREEMEMORY(space, info.snglbuffer);
        FREEMEMORY(space, info.multbuffer);
        FREEMEMORY(space, info.trnsbuffer);
      }
*/
}
/*
  void
sam_printSamlist (samlist_t *l, mapping_t* m, segemehl_t *nfo)

{
  char lf = '\n', alignlf = '\n';
  unsigned int i;
  FILE *dev = stdout;
  bl_fileBin_t *fx = NULL;
  MultiCharSeqAlignment *al;


  if((nfo->order || nfo->bisulfitemerging) && nfo->align){
    lf = 7;
    alignlf = 7;
  }

  for(i=0; i < l->noofrecs; i++) {
    if(!nfo->bufferedwrite || nfo->bins) { 
      dev = sam_getDevice(&l->recs[i], &fx, nfo);
      sam_printSamrec(dev, &l->recs[i], lf);

      if(nfo->align) {
        al = bl_getMapFragMCSA(&m->f[i]);  
        showAlignLF(al->al, dev, alignlf);
        fprintf(dev, "\n");
      }

      sam_closeDevice(fx, nfo);
      fx = NULL;
    } else { 
      char *line = sam_printSamrec2Buffer(&l->recs[i], lf);      
      bl_circBufferAddSave(&nfo->sambuffer[nfo->threadid], line, strlen(line));
      FREEMEMORY(NULL, line);
      if(nfo->align) {
        al = bl_getMapFragMCSA(&m->f[i]);  
        line = getAlignString(al->al, alignlf);
        bl_circBufferAddSave(&nfo->sambuffer[nfo->threadid], line, strlen(line));
        FREEMEMORY(NULL, line);
      }
    }
  }

  return ;
}
  if(!nfo->bufferedwrite || nfo->bins) {  
    dev = sam_getDevice(NULL, &fx, nfo);
    sam_printSamrec(dev, &samlist->recs[0], lf);
    sam_closeDevice(fx, nfo);
  } else {
    char *line = sam_printSamrec2Buffer(&samlist->recs[0], lf);
    bl_circBufferAddSave(&nfo->sambuffer[nfo->threadid], line, strlen(line));
    FREEMEMORY(NULL, line);
  }

      } else {
        char *tmp = bl_printMultiLocusSingle(nvrt, mseq, item);
        if(tmp) bl_circBufferAddSave(&nfo->snglbuffer[nfo->threadid], tmp, strlen(tmp));
        FREEMEMORY(NULL, tmp);

        tmp = bl_printMultiLocusJoint(nvrt, mseq, item);
        if(tmp) bl_circBufferAddSave(&nfo->multbuffer[nfo->threadid], tmp, strlen(tmp));
        FREEMEMORY(NULL, tmp);
      }

*/

