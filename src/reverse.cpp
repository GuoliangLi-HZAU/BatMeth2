/**********************************************************************
 * modified     program to reverse a file byte by byte                *
 *              H. J. Bernstein, yaya@dowling.edu                     *
 *              16 April 2002                                         *
 **********************************************************************/

#include <sys/types.h>
#include <sys/stat.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <stdlib.h>
#include <getopt.h>

void Parse_Command_line(int argc, char* argv[]);
void Seek_Method();

//{---------------------------- Command Line  -------------------------------------------------
option Long_Options[]=
{
{"help",0,NULL,'h'},
{"outputfile",1,NULL,'o'},
{"buffersize",1,NULL,'b'},
{"inputfile",0,NULL,'i'},
{0,0,0,0}
};
//}---------------------------- Command Line -------------------------------------------------
  FILE *in, *out;            /* The input and out file handles      */
  unsigned errflg = 0;            /* Count of argument list errors       */
  unsigned blocksize = 0;      /* The size of blocks to reverse       */
  unsigned bufsize = 0;        /* Actual buffer allocation >=512      */
  char* instr;               /* The path and name of the input      */
  char* outstr;              /* The path and name of the output     */
  char **endptr;             /* Pointer to pointer to invalid digit */
  char *buf;                 /* Pointer to buffer                   */
  char c;                    /* Option character                    */
  char tmpbuf[L_tmpnam+1];   /* Temporary file name                 */
  unsigned inlen, ipos;         /* Input file length, output position  */
  //long long inlen, ipos;         /* Input file length, output position  */
  unsigned nbytes;             /* Bytes in current transfer           */
  static struct stat64 fstatus;/* Input file stat data                */


	
int main(int argc, char * argv[]) {
       
  /* Extract options */

  instr = (char *)NULL;
  outstr = (char *)NULL;


  opterr = 0;                 /* handle bad options here            */

  Parse_Command_line(argc,argv);
  bufsize = blocksize;
  if (blocksize < 512)  bufsize = 512;
  if (blocksize < 1) blocksize = 1;

//Seek_Method();
//exit(0);
  if(! (buf=(char *)malloc(bufsize)) ) {
    fprintf(stderr,"reverse:  Insufficient memory to allocate buffer\n\n");
    exit(1);
  } 

  /* Check if the input filename string is empty or "-".  If so,
     copy to a temporary file.  If we do the copy, we have the file
     length by adding up the sizes of the transfers.  If we don't
     do the copy, we get the length with stat.                           */

  inlen = 0;
  if (!instr || strcmp(instr?instr:"-","-") == 0) {
    instr = strdup(tmpnam(&tmpbuf[0]));
    if ( (in = fopen64(instr, "w+")) == NULL) {
       fprintf(stderr,"reverse:  Can't open temporary file %s.\n", instr);
       exit(1);
    }
    while (nbytes = fread(buf, 1, bufsize, stdin)) {
      if(nbytes != fwrite(buf, 1, nbytes, in)) {
        fprintf(stderr,"reverse:  Failed to write %s.\n", instr);
        exit(1);
       }
      inlen +=nbytes;
    }
    fclose(in);
  } else {
    if ( stat64(instr, &fstatus) != 0 ) {
      fprintf(stderr,"reverse:  Unable to get length of %s.\n", instr);
    }
    inlen = fstatus.st_size;
  }

  if ( (in = fopen64(instr, "r")) == NULL ) {
    fprintf(stderr,"reverse:  Can't open input file %s.\n", instr);
    exit(1);
  }

  /* Check if the output filename string is empty or "-".  If so
     write to stdout.                                                 */  
  if (!outstr || strcmp(outstr?outstr:"-","-") == 0) {
    out = stdout;
  } else {
    if ( (out = fopen64(outstr, "a")) == NULL ) {
      fprintf(stderr,"reverse:  Can't open output file %s.\n", outstr);
      exit(1);
    }
  }

  /* Now copy from the back of the input file to the start
     of the output file, in blocks of size blocksize,
     reversing each block as it is read                               */

  for (ipos = 0; ipos < inlen; ipos += blocksize ) {
    unsigned postn;
    unsigned kbytes=blocksize;
    unsigned ii;
 
    if (ipos+blocksize > inlen) {
      kbytes += inlen-blocksize-ipos;
    }
    if (fseeko64(in, inlen-kbytes-ipos, SEEK_SET)) {printf("Seek error \n");exit(1);};
    if (kbytes != (nbytes = fread(buf,1,kbytes,in))) {
      fprintf(stderr,
        "reverse:  Read failed at position %ld for %ld bytes, got %ld.\n",
        (unsigned)(inlen-kbytes-ipos), (unsigned)kbytes, (unsigned)nbytes);
    }
    for(ii = 0; ii < nbytes/2; ii++) {
      char t;
      t = buf[nbytes-1-ii];
      buf[nbytes-1-ii] = buf[ii];
      buf[ii] = t;
    }
    fwrite(buf,1,nbytes,out);
  }

  /* close the input file, but not the output if it is stdout
     but we will do an fflush.  Note the the final exit would
     have done the same closes                                        */
  fclose(in);
  if (out != stdout) {
    fclose(out);      
  } else {
    fflush(out);      /* usually not necessary, but we'll be cautious */
  }

  exit(0);
}

//{-----------------------------  Seek Method  -------------------------------------------------
void Seek_Method()
{
  unsigned loc;
  FILE *in, *out;
  char ch;

  if((in = fopen64(instr, "rb")) == NULL) {
    printf("Cannot open input file.\n");
    exit(1);
  }
  if((out = fopen64(outstr, "wb"))==NULL) {
    //printf("Cannot open output file.\n");
    //exit(1);
    out=stdout;
  }
  if(setvbuf(in,NULL,_IOFBF,blocksize*sizeof(long))||setvbuf(out,NULL,_IOFBF,blocksize*sizeof(long)))
	{
		printf("Allocating disk buffers failed...\n");
		exit(1);
	}

  fseeko64(in, 0L, SEEK_END);
  loc = ftello64(in);
  if (out != stdout) printf("File Size %u\n",loc);
  /* copy file in reverse order */
  loc = loc-1; /* back up past end-of-file mark */
  
   for(;;){
    fseeko64(in, loc, SEEK_SET);
    ch = fgetc(in);
    fputc(ch, out);
    loc--;
    if (loc==0) break;
  }
   fseeko64(in, loc, SEEK_SET);//write last bit...
   ch = fgetc(in);
   fputc(ch, out);
   loc--;
  fclose(in);
  fclose(out);
}
//}-----------------------------  Seek Method  -------------------------------------------------

//{-----------------------------  Parse Command Line  -------------------------------------------------
void Parse_Command_line(int argc, char* argv[])
{
	int Current_Option=0;
	char* Short_Options ="ho:b:i:";//allowed options....
	char* This_Program = argv[0];//Current program name....
	char* Help_String=
"Parameters:\n"
" --help | -h\t\t\t\t Print help\n"
" --outputfile | -o <filename>\t\t Name of output file\n"
" --inputfile | -i <filename>\t\t Name of input file\n"
" --buffersize | -b <integer> \t\t Size of disk buffers\n"
;
	for(;;)	
	{
		Current_Option=getopt_long(argc, argv, Short_Options, Long_Options, NULL);
		if (Current_Option == -1 ) break;
		switch(Current_Option)
		{
			case 'h':
				printf("%s \n",Help_String);exit(0);
			case 'o':
				outstr=optarg;
				break;
			case 'i':
				instr=optarg;
				break;
			case 'b':
				blocksize=atol(optarg);
				break;
			default:
				printf("%s \n",Help_String);
				exit(0);
		}
	}	
}
