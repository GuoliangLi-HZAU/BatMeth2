#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define BUFSIZE 1000000

// g++ ./calCcontext.cpp -o calCcontext
// calCcontext genome.fa
// calculate every C or G loci, and give context, CG/CHG/CHH

int main(int argc, char **argv)
{
	FILE* iFPtr;
	FILE* oFPtr;
	char chrLenFile[500];
	
	if(argc != 2) {
		printf("Not enough args! Argc = %d.\n", argc);
		exit(1);
	}
	char* refSeqFile = argv[1];
	strcpy(chrLenFile, refSeqFile);
	strcat(chrLenFile, ".Csite.txt");
	
	fprintf(stderr, "refSeqFile = %s.\n", refSeqFile);
	//fprintf(stderr, "chrLenFile = %s.\n", chrLenFile);
	
	char* readBuffer = (char*)malloc(sizeof(char) * BUFSIZE);
	int len = 0;
	int lines = 0;
	
	char* token;
	char seps[] = " \t\n\r";
	char chrName[100];
    char context[10];
	
	if(!(iFPtr = fopen(refSeqFile, "r"))) {
		printf("File %s open error.\n", refSeqFile);
		exit(1);
	}
	
	if(!(oFPtr = fopen(chrLenFile, "w"))) {
		printf("File %s open error.\n", chrLenFile);
		exit(1);
	}
	int lineLen = 0;
	while(fgets(readBuffer, BUFSIZE, iFPtr)) {
		if(strlen(readBuffer) >= BUFSIZE - 1) {
            fprintf(stderr, "Too many characters in one row! Try to split the long row into several short rows (fewer than %d characters per row).\n", BUFSIZE);
            exit(1);
        }
		
		if(readBuffer[0] == '>') {
			/*if(lines > 0) {
				fprintf(oFPtr, "%s\t%d\n", chrName, len);
			}
            */
			// Save name
			token = strtok(readBuffer + 1, seps);	
			strcpy(chrName, token);
			len = 0;	
		}
		else {
			// Substract \n
            lineLen = strlen(readBuffer);
            for(int i=0; i<lineLen; i++){
                if(readBuffer[i]=='C' || readBuffer[i]=='c'){
                    if(i+1<lineLen){
                        if(readBuffer[i+1]=='G' || readBuffer[i+1]=='g'){
                            strcpy(context, "CG");
                        }else if(i+2<lineLen){
                            if(readBuffer[i+2]=='G' || readBuffer[i+2]=='g'){
                                strcpy(context, "CHG");
                            }else{
                                strcpy(context, "CHH");
                            }
                        }else{
                            strcpy(context, "CHH");
                        }
                    }else{
                        strcpy(context, "CHH");
                    }
                    fprintf(oFPtr, "%s\t%d\t%s\tC\n", chrName, len+i+1, context);
                }else if(readBuffer[i]=='G' || readBuffer[i]=='g'){
                    if(i>0){
                        if(readBuffer[i-1]=='C' || readBuffer[i-1]=='c'){
                            strcpy(context, "CG");
                        }else if(i>1){
                            if(readBuffer[i-2]=='C' || readBuffer[i-2]=='c'){
                                strcpy(context, "CHG");
                            }else{
                                strcpy(context, "CHH");
                            }
                        }else{
                            strcpy(context, "CHH");
                        }
                    }else{
                        strcpy(context, "CHH");
                    }
                    fprintf(oFPtr, "%s\t%d\t%s\tG\n", chrName, len+i+1, context);
                }
            }
            len += strlen(readBuffer) - 1;
		}
		
		lines++;
	}
	
	/*
    if(lines > 0) {
		fprintf(oFPtr, "%s\t%d\n", chrName, len);
	}
	*/
	fclose(iFPtr);
	fclose(oFPtr);
	
	return 0;
}
