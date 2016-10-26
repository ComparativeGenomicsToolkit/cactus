/*  Last edited: Feb  1 15:35 2002 (klh) */
/**********************************************************************
 ** FILE: align.c
 ** NOTES:
 **   Functions for the manipulation of multiple sequence
 **   alignments
**********************************************************************/

#include "align.h"


static char *comment = "#";
static char *whitespace = " \t\r\n";
static char *terminator = "//";

/**********************************************************************
 FUNCTION: free_Alignment
 DESCRIPTION: 
   Frees the memory used by the given alignment
 ARGS:
   An Alignment
 RETURNS: A null pointer
 NOTES:
**********************************************************************/
void *free_Alignment( struct Alignment *al ) {
  unsigned int i;
  
  if (al != NULL) {
    if (al->seqs != NULL) {
      for( i=0; i < al->numseqs; i++) {
	al->seqs[i] = free_Sequence( al->seqs[i] );
      }
      al->seqs = free_util( al->seqs );
    }
    al = free_util( al );
  }
  return al;
}



/**********************************************************************
 FUNCTION: read_MUL_Alignment
 DESCRIPTION: 
   Reads in a muliple alignment from the given file handle and returns
   it
 ARGS: 
   A file handle
 RETURNS: A number denoting the status:
   A pointer to the struct Alignment object created, or NULL of there
   was an error parsing the file
 NOTES: It is assumed that the aligment file is in Pfam (MUL) format.
   Garbage results should be expected if the input file is not in this 
   format

   The function allocates all the memory necessary for the 
   alignment. The caller should call free_Alignment (align.h) to 
   free this memory when the alignment is no longer needed
 **********************************************************************/

struct Alignment *read_MUL_Alignment( FILE *stream ) {
  struct Alignment *aln;
  unsigned int index = 0;
  char tempname[MAX_NAME_LENGTH];
  int c;

  aln = (struct Alignment *) malloc_util( sizeof( struct Alignment ));
  aln->numseqs = 0;
  aln->seqs = NULL;
  aln->length = 0;
    

  while (( c = fgetc(stream)) != EOF) {
    if (aln->numseqs % SEQ_BLOCKSIZE == 0) {
      /* we need to allocate some more memory. But is it a malloc or realloc ? */

      if (aln->numseqs == 0) {
      	aln->seqs = (struct Sequence **) malloc_util( SEQ_BLOCKSIZE * sizeof( struct Sequence *));
      }
      else {
	aln->seqs = (struct Sequence **) 
	  realloc_util( aln->seqs, (aln->numseqs + SEQ_BLOCKSIZE) * sizeof( struct Sequence *));
      }
    }

    index = 0;
    while (! isspace(c) ) {
      tempname[index++] = c;
      c = fgetc(stream);
    }
    tempname[index] = '\0';
    while ( isspace(c = fgetc(stream)) );


    aln->seqs[aln->numseqs] = empty_Sequence();
    aln->seqs[aln->numseqs]->name = (char *) malloc_util( MAX_NAME_LENGTH * sizeof(char));
    strcpy(aln->seqs[aln->numseqs]->name, tempname);

    index = 0;
    do {
      if (index % RES_BLOCKSIZE == 0) {
      /* we need to allocate some more memory. But is it a malloc or realloc ? */
	if (index == 0) {
	  aln->seqs[aln->numseqs]->seq = (char *) malloc_util( RES_BLOCKSIZE * sizeof(char));
	}
	else {
	  aln->seqs[aln->numseqs]->seq = (char *) 
	    realloc_util( aln->seqs[aln->numseqs]->seq, (index + RES_BLOCKSIZE) * sizeof(char));
	}	
      }
      if (! isspace(c) )
	aln->seqs[aln->numseqs]->seq[index++] = c;
      c = fgetc(stream);

    } while ( c != '\n' && c != EOF);

    /* Since we now know the length of the sequence, we can do a bit of resizing */

    if (aln->numseqs == 0 || index < aln->length) {
      /* The first sequence is the trend setter */
      aln->length = index;
    }
    /* First, resize the current sequence */
    aln->seqs[aln->numseqs]->seq = (char *)
	   realloc_util( aln->seqs[aln->numseqs]->seq, aln->length * sizeof(char));

    /* if we reduced the aligment length, The earlier sequences are too big;
       Since we will only read up to the length of the smallest sequence, this
       will only prove a problem if we run out of memory; however, if the user 
       given non-flused alignments, then they are asking for everything they get...
    */

    aln->numseqs++;

  }

  /* now we can resize the seq array, and all the actual sequences themselves; this resizing
     will only save significant memory if a non-flush alignment has been given. Is it worth
     it? Well, we have to set the length of each sequence, so may as well do it while we are
     here*/

  aln->seqs = (struct Sequence **) realloc_util(  aln->seqs, aln->numseqs * sizeof( struct Sequence *));
  for (index=0; index < aln->numseqs; index++) {
    aln->seqs[index]->length = aln->length;
    aln->seqs[index]->seq = (char *) 
      realloc_util(  aln->seqs[index]->seq, aln->length * sizeof(char));
  }

  return aln;
}



/**********************************************************************
 FUNCTION: write_MUL_Alignment
 DESCRIPTION: 
   Prints a rep. of the alignment to the given handle (in MUL format)
 ARGS:
   FILE *
   struct Alignment (align.h)
 RETURNS: struct Alignment (align.h)
 NOTES:
 **********************************************************************/

void write_MUL_Alignment( FILE *handle, struct Alignment *al ) {
  unsigned int i,j;

  for( i=0; i < al->numseqs; i++ ) {
    fprintf( handle, "%-24s ", al->seqs[i]->name);
    for( j=0; j < al->length; j++) {
      fprintf( handle, "%c",  al->seqs[i]->seq[j]);
    }
    fprintf( handle, "\n");
  }
  fflush( handle );

}



/**********************************************************************
 FUNCTION: read_Stockholm_Alignment
 DESCRIPTION: 
   This function fills a simple alignment structure from an
   file assumed to be in Stockholm format. At this stage, I am ignoring
   all mark-up information (all lines beginning with '#' are ignored).
   The function also allows for wrapped alignments. Note that Pfam
   alignments in MUL format will be handled correctly by this function
 ARGS:
   FILE *
 RETURNS: struct Alignment (align.h)
 NOTES:
 **********************************************************************/

struct Alignment *read_Stockholm_Alignment( FILE *handle ) {
  struct Alignment *aln;

  char *line = (char *) malloc_util( MAX_LINE_LEN * sizeof(char) );
  unsigned int lineAllocSize = MAX_LINE_LEN; 
  char *name_ptr = NULL;
  char *seq_ptr = NULL;

  unsigned int i, j;
  unsigned int numblocks = 0;
  unsigned int thisseq = 0;
  unsigned int last_idx = 0;
  unsigned int got_line = 0;
  unsigned int saw_blank_line = 0;

  aln = (struct Alignment *) malloc_util( sizeof( struct Alignment ));
  aln->numseqs = 0;
  aln->seqs = NULL;
  aln->length = 0;

  numblocks = 0;
  thisseq = 0;

  /* skip to the start of the first block */
  got_line = 0;
  
  do {
    if (fgets(line, lineAllocSize, handle) == NULL)
      break;
    else {      
      if (strchr(line, '\n') == NULL) { /* did not read whole line! */
        do {
          lineAllocSize += MAX_LINE_LEN;
          line = (char *) realloc_util( line, lineAllocSize * sizeof(char) );
        } while (fgets(line + lineAllocSize - MAX_LINE_LEN - 1, MAX_LINE_LEN+1,
                       handle) != NULL  &&
                 strchr(line, '\n') == NULL);
      }
      
      if (strchr( comment, *line ) != NULL)
	continue;
      else if (strncmp(line, terminator, 2) == 0) {
	break;
      }
      else {
        if ( (name_ptr = strtok(line, whitespace)) != NULL &&
             (seq_ptr = strtok( NULL, whitespace )) != NULL)
          got_line = 1;        
      }
    }
  } while (! got_line);

  if (! got_line)
    fatal_util( "The alignment file appears to have no valid lines in Stockholm format\n");

  while (got_line) {

    if (numblocks == 0) {
      /* if this is the first block, set up the memory for the sequences themselves */
      if (thisseq == 0) {
	aln->seqs = (struct Sequence **) malloc_util ( sizeof( struct Sequence *));
	aln->numseqs = 1;
      }
      else { 
	aln->numseqs++;
	aln->seqs = (struct Sequence **) realloc_util( aln->seqs, aln->numseqs * sizeof( struct Sequence *));
      }
      
      aln->seqs[thisseq] = empty_Sequence();
      aln->seqs[thisseq]->length = 0;
      aln->seqs[thisseq]->name = (char *) malloc_util( (strlen( name_ptr ) + 1) * sizeof(char));
      strcpy( aln->seqs[thisseq]->name, name_ptr );
      
      /* get rid of whitespace at end of line */
      for ( last_idx = strlen(seq_ptr) - 1; strchr(whitespace, seq_ptr[last_idx]) != NULL; last_idx--);	  
      aln->seqs[thisseq]->seq = (char *) malloc_util( (last_idx + 1) * sizeof (char ) );
    }
    else {
      if (strcmp(aln->seqs[thisseq]->name, name_ptr) != 0) {
	warning_util("Your seq names are inconsistent across blocks. Using the names in the first block");
      }
      
      /* The following accounts for the fact that there may be whitepasce at the end of the line */
      for ( last_idx = strlen(seq_ptr) - 1; strchr(whitespace, seq_ptr[last_idx]) != NULL; last_idx--);

      aln->seqs[thisseq]->seq = (char *) realloc_util(aln->seqs[thisseq]->seq,	      
						      (aln->seqs[thisseq]->length + last_idx + 1) * sizeof(char)); 
      
    }
    
    /* and finally, copy the string across */
    for( j=0; j <= last_idx; j++, aln->seqs[thisseq]->length++) 
      aln->seqs[thisseq]->seq[aln->seqs[thisseq]->length] = seq_ptr[j];
    
    /* read the next line. If we come across one or more blank lines while
       we do so, then assume that we are starting a new block */

    got_line = saw_blank_line = 0;
    do {
      if (fgets(line, lineAllocSize, handle) != NULL) {
        if (strchr(line, '\n') == NULL) { /* did not read whole line! */
          do {
            lineAllocSize += MAX_LINE_LEN;
            line = (char *) realloc_util( line, lineAllocSize * sizeof(char) );
          } while (fgets(line + lineAllocSize - MAX_LINE_LEN-1, MAX_LINE_LEN+1,
                         handle) != NULL  &&
                   strchr(line, '\n') == NULL);
        }
        
	if (strchr( comment, *line ) != NULL)
	  continue;
	else if (strncmp(line, terminator, 2) == 0) {
	  break;
	}
	else {
	  if ( (name_ptr = strtok(line, whitespace)) == NULL) {
	    saw_blank_line = 1;
	  }
	  else if ( (seq_ptr = strtok(NULL, whitespace)) != NULL)
	    got_line = 1;
	}
      }
      else 
	break;

    } while (! got_line );


    if (got_line) {
      if (saw_blank_line) {
	numblocks++;
	thisseq = 0;
      }
      else {
	thisseq++;
      }
    }
  }

  free_util(line);

  /* just before we return, check the alignment */
  for(i=0; i < aln->numseqs; i++) {
    if (i==0)
      aln->length = aln->seqs[i]->length;
    else if (aln->length != aln->seqs[i]->length) {
      fatal_util( "Your alignment segments were of different sizes!\n");
    }
  }

  return aln;
}

