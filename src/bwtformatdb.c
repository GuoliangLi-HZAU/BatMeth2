#include "config.h"
#ifdef MMX
/*

   BWTFormatdb.c		Build index for FASTA database

   This program builds index for FASTA database for use of BWTBlastn.

   Copyright (C) 2006, Wong Chi Kwong.

   This program is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public License
   as published by the Free Software Foundation; either version 2
   of the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

*/

#include <stdio.h>
#include <stdlib.h>
#include "TypeNLimit.h"
#include "BWTConstruct.h"
#include "MiscUtilities.h"
#include "DNACount.h"
#include "TextConverter.h"
#include "MemManager.h"
#include "iniparser.h"
#include "HSP.h"
#include "Timing.h"

// Database and ini
dictionary *ParseInput(int argc, char** argv);
void ParseIniFile(char *iniFileName);
void ProcessIni();
void ValidateIni();
void PrintIni();
void PrintShortDesc();
void PrintHelp();

void ProcessFileName(char *outputFileName, const char *inputFileName, const char *databaseName);


	// Parameters
	char IniFileName[MAX_FILENAME_LEN+1];
	int Confirmation;
	int Pacfileonly;
	
	// BuildTasks parameters
	int ParseFASTA = TRUE;
	int BuildBWT = TRUE;
	int BuildSaValue = TRUE;
	int BuildCachedSaIndex = TRUE;
	int BuildSaIndex = TRUE;

	// Memory parameters
	unsigned int PoolSize = 2097152;				// 2M  - fixed; not configurable through ini

	// Display parameters
	int ShowProgress = FALSE;

	// Database parameters
	char FASTAFileName[MAX_FILENAME_LEN+1] = "";
	char DatabaseName[MAX_FILENAME_LEN+1] = "";
	char AnnotationFileName[MAX_FILENAME_LEN+1] = "*.ann";
	char AmbiguityFileName[MAX_FILENAME_LEN+1] = "*.amb";
	char PackedDNAFileName[MAX_FILENAME_LEN+1] = "*.pac";
	char BWTCodeFileName[MAX_FILENAME_LEN+1] = "*.bwt";
	char BWTOccValueFileName[MAX_FILENAME_LEN+1] = "*.fmv";
	char SaValueFileName[MAX_FILENAME_LEN+1] = "*.sa";
	char CachedSaIndexFileName[MAX_FILENAME_LEN+1] = "*.sai";

	// Parse FASTA parameters
	unsigned int FASTARandomSeed = 0;
	int MaskLowerCase = FALSE;

	// Build BWT parameters
	unsigned int OccValueFreq = 256;
	float TargetNBit = 2.5;
	unsigned int InitialMaxBuildSize = 10000000;
	unsigned int IncMaxBuildSize = 10000000;

	// Build SA value parameters
	unsigned int SaValueFreq = 8;

	// Build SA index parameters
	unsigned int CachedSaIndexNumOfChar = 12;



int main(int argc, char** argv) {

	char c;
	MMPool *mmPool;
	dictionary *programInput;
	double startTime;
	double elapsedTime = 0, totalElapsedTime = 0;

	char filename[MAX_FILENAME_LEN+1];
	BWT *bwt = NULL;
	unsigned int textLength = 0;
	unsigned int numSeq;

	BWTInc *bwtInc = NULL;

	// Program input
	programInput = ParseInput(argc, argv);
	if (Pacfileonly)
	{
		ParseFASTA = TRUE;
		BuildBWT = FALSE;
		BuildSaValue = FALSE;
		BuildSaIndex = FALSE;
		BuildCachedSaIndex = FALSE;
	}
	PrintShortDesc();

	// Ini
	if (strcmp(argv[0] + strlen(argv[0]) - 4, ".exe") == 0) {
		*(argv[0] + strlen(argv[0]) - 4) = '\0';
	}
	sprintf(filename, "%s.ini", argv[0]);
	ParseIniFile(filename);
	printf("\n");
	ProcessIni();
	ValidateIni();
	PrintIni();

	if (Confirmation == TRUE) {
		printf("Press Y to go or N to cancel. ");
        c = (char)getchar();
		while (c != 'y' && c != 'Y' && c != 'n' && c!= 'N') {
			c = (char)getchar();
		}
		if (c == 'n' || c == 'N') {
			exit(0);
		}
	}

	startTime = setStartTime();

	MMMasterInitialize(1, 0, FALSE, NULL);
	mmPool = MMPoolCreate(PoolSize);

	// Parse FASTA file to produce packed DNA and annotation file
	if (ParseFASTA == TRUE) {

		printf("Parsing FASTA file..\n");

		numSeq = HSPParseFASTAToPacked(FASTAFileName, AnnotationFileName, PackedDNAFileName, AmbiguityFileName, FASTARandomSeed, MaskLowerCase);

		printf("Finished. Parsed %u sequences.\n", numSeq);

		elapsedTime = getElapsedTime(startTime) - totalElapsedTime;
		printf("Elapsed time = ");
		printElapsedTime(stdout, FALSE, FALSE, TRUE, 2, elapsedTime);
		totalElapsedTime += elapsedTime;
		printf("\n");

	}

	// Construct BWTInc from text
	if (BuildBWT == TRUE) {

		printf("Building BWT..\n");

		bwtInc = BWTIncConstructFromPacked(mmPool, PackedDNAFileName, ShowProgress, 
										   TargetNBit, InitialMaxBuildSize, IncMaxBuildSize);

		printf("Finished constructing BWT in %u iterations.  ", bwtInc->numberOfIterationDone);
		elapsedTime = getElapsedTime(startTime) - totalElapsedTime;
		printf("Elapsed time = ");
		printElapsedTime(stdout, FALSE, FALSE, TRUE, 2, elapsedTime);
		totalElapsedTime += elapsedTime;
		printf("\n");

		printf("Saving BWT..\n");
		BWTSaveBwtCodeAndOcc(bwtInc->bwt, BWTCodeFileName, BWTOccValueFileName);
		printf("Finished saving BWT.  ");
		elapsedTime = getElapsedTime(startTime) - totalElapsedTime;
		printf("Elapsed time = ");
		printElapsedTime(stdout, FALSE, FALSE, TRUE, 2, elapsedTime);
		totalElapsedTime += elapsedTime;
		printf("\n");

		textLength = bwtInc->bwt->textLength;

		BWTIncFree(mmPool, bwtInc);

	}

	// Load BWT
	if (BuildSaValue || BuildCachedSaIndex) {

		printf("Loading BWT...\n");

		bwt = BWTLoad(mmPool, BWTCodeFileName, BWTOccValueFileName, NULL, NULL, NULL, NULL);

		printf("Finished loading BWT.  ");

		elapsedTime = getElapsedTime(startTime) - totalElapsedTime;
		printf("Elapsed time = ");
		printElapsedTime(stdout, FALSE, FALSE, TRUE, 2, elapsedTime);
		totalElapsedTime += elapsedTime;
		printf("\n");

		textLength = bwt->textLength;

	}

	if (BuildSaValue) {

		printf("Building SA value...\n");
		
		if (ShowProgress) {
			BWTGenerateSaValue(bwt, SaValueFreq, bwt->textLength / SaValueFreq / 10);
		} else {
			BWTGenerateSaValue(bwt, SaValueFreq, 0);
		}
		BWTSaveSaValue(bwt, SaValueFileName);

		printf("Finished building SA value.  ");

		elapsedTime = getElapsedTime(startTime) - totalElapsedTime;
		printf("Elapsed time = ");
		printElapsedTime(stdout, FALSE, FALSE, TRUE, 2, elapsedTime);
		totalElapsedTime += elapsedTime;
		printf("\n");

	}

	if (BuildCachedSaIndex) {

		printf("Building cached SA index...\n");
		
		BWTGenerateCachedSaIndex(bwt, CachedSaIndexNumOfChar, CachedSaIndexFileName);

		printf("Finished building cached SA index.  ");

		elapsedTime = getElapsedTime(startTime) - totalElapsedTime;
		printf("Elapsed time = ");
		printElapsedTime(stdout, FALSE, FALSE, TRUE, 2, elapsedTime);
		totalElapsedTime += elapsedTime;
		printf("\n");

	}

	// Free BWT
	if (BuildSaValue || BuildCachedSaIndex) {
		BWTFree(mmPool, bwt);
	}

	// Finished all construction tasks
	printf("Finished all tasks.  ");
	totalElapsedTime = getElapsedTime(startTime);
	printf("Total elapsed time = ");
	printElapsedTime(stdout, FALSE, FALSE, TRUE, 2, totalElapsedTime);
	printf("\n");

	MMMasterPrintReport(stdout, FALSE, FALSE, FALSE);
	if (BuildSaValue) {
		fprintf(stdout, "Number of char   :  %u\n", textLength);
		fprintf(stdout, "Bit per char     :  %.2f\n", (float)MMMasterMaxTotalByteDispatched() * BITS_IN_BYTE / textLength);
		printf("\n");
	}

	MMPoolFree(mmPool);

	iniparser_freedict(programInput);

	return 0;

}

dictionary *ParseInput(int argc, char** argv) {

	dictionary *programInput;
	char t1[3] = "-c";	// specify that this is a boolean type parameter
	char t2[3] = "-U";	// specify that this is a boolean type parameter
	char t3[3] = "-p";	// specify that this is a boolean type parameter
	char *d[3];

	d[0] = t1;
	d[1] = t2;
	d[2] = t3;
	
	programInput = paraparser_load(argc, argv, 3, d);	// 2 boolean type parameters

	// Get database name
	if (!iniparser_find_entry(programInput, "argument:1")) {
		PrintHelp();
		exit(1);
	}
	iniparser_copystring(programInput, "argument:1", DatabaseName, DatabaseName, MAX_FILENAME_LEN);
	if (strlen(DatabaseName) + 4 > MAX_FILENAME_LEN) {
		PrintHelp();
		exit(1);
	}

	// Get FASTA file name
	iniparser_copystring(programInput, "argument:2", FASTAFileName, DatabaseName, MAX_FILENAME_LEN);
	if (strlen(FASTAFileName) > MAX_FILENAME_LEN) {
		PrintHelp();
		exit(1);
	}


	// Whether confirmation is needed
	Confirmation = iniparser_find_entry(programInput, "parameter:-c");

	MaskLowerCase = iniparser_find_entry(programInput, "parameter:-U");

	Pacfileonly = iniparser_find_entry(programInput, "parameter:-p");

	return programInput;

}

void ParseIniFile(char *iniFileName) {

	dictionary *ini;

	printf("Loading %s ..", iniFileName);
	ini = iniparser_load(iniFileName, FALSE);
	if (ini == NULL) {
		printf("not found.\n");
		return;
	}
	printf("done.\n");

	// BuildTasks parameters
	ParseFASTA = iniparser_getboolean(ini, "BuildTasks:ParseFASTA", ParseFASTA);
	BuildBWT = iniparser_getboolean(ini, "BuildTasks:BuildBWT", BuildBWT);
	BuildSaValue = iniparser_getboolean(ini, "BuildTasks:BuildSaValue", BuildSaValue);
	BuildCachedSaIndex = iniparser_getboolean(ini, "BuildTasks:BuildCachedSaIndex", BuildCachedSaIndex);

	// Display parameters
	ShowProgress = iniparser_getboolean(ini, "Display:ShowProgress", ShowProgress);

	// Parse FASTA parameters
	FASTARandomSeed = iniparser_getint(ini, "ParseFASTA:RandomSeed", FASTARandomSeed);
	if (FASTARandomSeed == 0) {
		FASTARandomSeed = getRandomSeed();
	}

	// Build BWT parameters
	OccValueFreq = iniparser_getint(ini, "BuildBWT:OccValueFreq", OccValueFreq);
	TargetNBit = (float)iniparser_getdouble(ini, "BuildBWT:TargetNBit", TargetNBit);
	InitialMaxBuildSize = iniparser_getint(ini, "BuildBWT:InitialMaxBuildSize", InitialMaxBuildSize);
	IncMaxBuildSize = iniparser_getint(ini, "BuildBWT:IncMaxBuildSize", IncMaxBuildSize);

	// Build SA value parameters
	SaValueFreq = iniparser_getint(ini, "BuildSAValue:SaValueFreq", SaValueFreq);

	// Build SA index parameters
	CachedSaIndexNumOfChar = iniparser_getint(ini, "BuildSAIndex:CachedSaIndexNumOfChar", CachedSaIndexNumOfChar);

	// Database parameters
	iniparser_copystring(ini, "Database:AnnotationFileName", AnnotationFileName, AnnotationFileName, MAX_FILENAME_LEN);
	iniparser_copystring(ini, "Database:AmbiguityFileName", AmbiguityFileName, AmbiguityFileName, MAX_FILENAME_LEN);
	iniparser_copystring(ini, "Database:PackedDNAFileName", PackedDNAFileName, PackedDNAFileName, MAX_FILENAME_LEN);
	iniparser_copystring(ini, "Database:BWTCodeFileName", BWTCodeFileName, BWTCodeFileName, MAX_FILENAME_LEN);
	iniparser_copystring(ini, "Database:BWTOccValueFileName", BWTOccValueFileName, BWTOccValueFileName, MAX_FILENAME_LEN);
	iniparser_copystring(ini, "Database:SaValueFileName", SaValueFileName, SaValueFileName, MAX_FILENAME_LEN);
	iniparser_copystring(ini, "Database:CachedSaIndexFileName", CachedSaIndexFileName, CachedSaIndexFileName, MAX_FILENAME_LEN);

	iniparser_freedict(ini);

}

void ProcessIni() {

	ProcessFileName(AnnotationFileName, AnnotationFileName, DatabaseName);
	ProcessFileName(AmbiguityFileName, AmbiguityFileName, DatabaseName);
	ProcessFileName(PackedDNAFileName, PackedDNAFileName, DatabaseName);
	ProcessFileName(BWTCodeFileName, BWTCodeFileName, DatabaseName);
	ProcessFileName(BWTOccValueFileName, BWTOccValueFileName, DatabaseName);
	ProcessFileName(SaValueFileName, SaValueFileName, DatabaseName);
	ProcessFileName(CachedSaIndexFileName, CachedSaIndexFileName, DatabaseName);

}

void ValidateIni() {

	if (!ParseFASTA && !BuildBWT && !BuildSaValue && !BuildCachedSaIndex) {
		fprintf(stderr, "No action is specified!\n");
		exit(1);
	}
	if (ParseFASTA) {
		if (PackedDNAFileName[0] == '\0') {
			fprintf(stderr, "Packed DNA file name is not specified!\n");
			exit(1);
		}
		if (AnnotationFileName[0] == '\0') {
			fprintf(stderr, "Annotation file name is not specified!\n");
			exit(1);
		}
		if (AmbiguityFileName[0] == '\0') {
			fprintf(stderr, "Ambiguity file name is not specified!\n");
			exit(1);
		}
	}
	if (BuildBWT) {
		if (PackedDNAFileName[0] == '\0') {
			fprintf(stderr, "Packed DNA file is not specified!\n");
			exit(1);
		}
		if (BWTCodeFileName[0] == '\0') {
			fprintf(stderr, "BWT code file name is not specified!\n");
			exit(1);
		}
		if (BWTOccValueFileName[0] == '\0') {
			fprintf(stderr, "BWT Occ value file name is not specified!\n");
			exit(1);
		}
		if (TargetNBit < 2.5) {
			fprintf(stderr, "Target NBit should be at least 2.5!\n");
			exit(1);
		}
	}
	if (BuildSaValue) {
		if (BWTCodeFileName[0] == '\0') {
			fprintf(stderr, "BWT code file is not specified!\n");
			exit(1);
		}
		if (BWTOccValueFileName[0] == '\0') {
			fprintf(stderr, "BWT Occ value file is not specified!\n");
			exit(1);
		}
		if (SaValueFileName[0] == '\0') {
			fprintf(stderr, "SA value file name is not specified!\n");
			exit(1);
		}
		if (SaValueFreq <= 0) {
			fprintf(stderr, "SA value frequency must > 0!\n");
			exit(1);
		}
	}

	if (BuildCachedSaIndex) {
		if (BWTCodeFileName[0] == '\0') {
			fprintf(stderr, "BWT code file is not specified!\n");
			exit(1);
		}
		if (BWTOccValueFileName[0] == '\0') {
			fprintf(stderr, "BWT Occ value file is not specified!\n");
			exit(1);
		}
		if (CachedSaIndexFileName[0] == '\0') {
			fprintf(stderr, "SA index file name is not specified!\n");
			exit(1);
		}
		if (CachedSaIndexNumOfChar <= 0) {
			fprintf(stderr, "SA index number of character must > 0!\n");
			exit(1);
		}
		if (CachedSaIndexNumOfChar > 13) {
			fprintf(stderr, "SA index number of character must <= 13!\n");
			exit(1);
		}
	}

}


void PrintIni() {

	char boolean[2];

	boolean[0] = 'N';
	boolean[1] = 'Y';

	printf("Parse FASTA file    : %c\n", boolean[ParseFASTA]);
	printf("Build BWT           : %c\n", boolean[BuildBWT]);
	printf("Build SA value      : %c\n", boolean[BuildSaValue]);
	printf("Build SA index      : %c\n", boolean[BuildCachedSaIndex]);
	printf("\n");

	printf("Show progress       : %c\n", boolean[ShowProgress]);
	printf("\n");

	if (ParseFASTA) {
		printf("Parse FASTA :\n");
		printf("Mask lower case         : %c\n", boolean[MaskLowerCase]);
		printf("Random seed             : %u\n", FASTARandomSeed);
		printf("\n");
	}

	if (BuildBWT) {
		printf("Build BWT :\n");
		printf("Target N Bits           : %.2f\n", TargetNBit);
		printf("Occ value frequency     : %u\n", OccValueFreq);
		printf("Initial Max Build Size  : %u    Inc Max Build Size : %u\n", 
				InitialMaxBuildSize, IncMaxBuildSize);
		printf("\n");
	}

	if (BuildSaValue) {
		printf("Build SA value :\n");
		printf("SA value frequency      : %u\n", SaValueFreq);
		printf("\n");
	}

	if (BuildCachedSaIndex) {
		printf("Build SA index :\n");
		printf("SA index no. of char    : %u\n", CachedSaIndexNumOfChar);
		printf("\n");
	}

	printf("Annotation file          : %s\n", AnnotationFileName);
	printf("Ambigurity file          : %s\n", AmbiguityFileName);
	printf("Packed DNA file          : %s\n", PackedDNAFileName);
	printf("BWT Code file            : %s\n", BWTCodeFileName);
	printf("BWT Occ value file       : %s\n", BWTOccValueFileName);
	printf("SA value file            : %s\n", SaValueFileName);
	printf("Cached SA index file     : %s\n", CachedSaIndexFileName);
	printf("\n");

}

void PrintShortDesc() {

	printf("BWTFormatdb v1.0, Copyright (C) 2006, Wong Chi Kwong.\n");
	printf("BWTFormatdb comes with ABSOLUTELY NO WARRENTY.\n");
	printf("BWTFormatdb is free software, and you are welcome to\n");
	printf("redistribute it under certain conditions.\n");
	printf("For details type BWTFormatdb.\n");
	printf("\n");

}

void PrintHelp() {

	printf("BWTFormatdb v1.0, Copyright (C) 2006, Wong Chi Kwong.\n");
	printf("\n");

	printf("This program is free software; you can redistribute it and/or\n");
	printf("modify it under the terms of the GNU General Public License\n");
	printf("as published by the Free Software Foundation; either version 2\n");
	printf("of the License, or (at your option) any later version.\n");
	printf("\n");

	printf("This program is distributed in the hope that it will be useful,\n");
	printf("but WITHOUT ANY WARRANTY; without even the implied warranty of\n");
	printf("MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n");
	printf("GNU General Public License for more details.\n");
	printf("\n");

	printf("You should have received a copy of the GNU General Public License\n");
	printf("along with this program; if not, write to the Free Software\n");
	printf("Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.\n");
	printf("\n");

	printf("Syntax: BWTFormatdb <database> [FASTA file; default = <database>]\n");
	printf("                               [-c Confirm]\n");
	printf("                               [-U Mask lower case]\n");

}

void ProcessFileName(char *outputFileName, const char *inputFileName, const char *databaseName) {

	char tempChar[MAX_FILENAME_LEN];
	unsigned int i;

	if (inputFileName == NULL) {
		if (outputFileName != inputFileName) {
			outputFileName[0] = '\0';
		}
		return;
	}

	if (strlen(databaseName) + strlen(inputFileName) > MAX_FILENAME_LEN) {
		fprintf(stderr, "File length is too long!\n");
		exit(1);
	}

	strncpy(tempChar, inputFileName, MAX_FILENAME_LEN);

	// locate the *
	for (i=0; i<MAX_FILENAME_LEN; i++) {
		if (tempChar[i] == '*') {
			break;
		}
	}
	if (i<MAX_FILENAME_LEN) {
		tempChar[i] = '\0';
		sprintf(outputFileName, "%s%s%s", tempChar, databaseName, tempChar + i + 1);
	} else {
		sprintf(outputFileName, "%s", tempChar);
	}

}
#else

/*

   BWTFormatdb.c		Build index for FASTA database

   This program builds index for FASTA database for use of BWTBlastn.

   Copyright (C) 2006, Wong Chi Kwong.

   This program is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public License
   as published by the Free Software Foundation; either version 2
   of the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

*/

#include <stdio.h>
#include <stdlib.h>
#include "TypeNLimit.h"
#include "BWTConstruct.h"
#include "MiscUtilities.h"
#include "DNACount.h"
#include "TextConverter.h"
#include "MemManager.h"
#include "iniparser.h"
#include "HSP.h"
#include "Timing.h"

// Database and ini
dictionary *ParseInput(int argc, char** argv);
void ParseIniFile(char *iniFileName);
void ProcessIni();
void ValidateIni();
void PrintIni();
void PrintShortDesc();
void PrintHelp();

void ProcessFileName(char *outputFileName, const char *inputFileName, const char *databaseName);


	// Parameters
	char IniFileName[MAX_FILENAME_LEN+1];
	int Confirmation;
	int Pacfileonly;
	
	// BuildTasks parameters
	int ParseFASTA = TRUE;
	int BuildBWT = TRUE;
	int BuildSaValue = TRUE;
	int BuildSaIndex = TRUE;
	int BuildCachedSaIndex = TRUE;

	// Memory parameters
	unsigned int PoolSize = 2097152;				// 2M  - fixed; not configurable through ini

	// Display parameters
	int ShowProgress = FALSE;

	// Database parameters
	char FASTAFileName[MAX_FILENAME_LEN+1] = "";
	char DatabaseName[MAX_FILENAME_LEN+1] = "";
	char AnnotationFileName[MAX_FILENAME_LEN+1] = "*.ann";
	char AmbiguityFileName[MAX_FILENAME_LEN+1] = "*.amb";
	char PackedDNAFileName[MAX_FILENAME_LEN+1] = "*.pac";
	char BWTCodeFileName[MAX_FILENAME_LEN+1] = "*.bwt";
	char BWTOccValueFileName[MAX_FILENAME_LEN+1] = "*.fmv";
	char SaValueFileName[MAX_FILENAME_LEN+1] = "*.sa";
	char SaIndexFileName[MAX_FILENAME_LEN+1] = "*.sai";

	// Parse FASTA parameters
	unsigned int FASTARandomSeed = 0;
	int MaskLowerCase = FALSE;

	// Build BWT parameters
	unsigned int OccValueFreq = 256;
	float TargetNBit = 2.5;
	unsigned int InitialMaxBuildSize = 10000000;
	unsigned int IncMaxBuildSize = 10000000;

	// Build SA value parameters
	unsigned int SaValueFreq = 8;

	// Build SA index parameters
	unsigned int SaIndexNumOfChar = 12;



int main(int argc, char** argv) {

	char c;
	MMPool *mmPool;
	dictionary *programInput;
	double startTime;
	double elapsedTime = 0, totalElapsedTime = 0;

	char filename[MAX_FILENAME_LEN+1];
	BWT *bwt = NULL;
	unsigned int textLength = 0;
	unsigned int numSeq;

	BWTInc *bwtInc = NULL;

	// Program input
	programInput = ParseInput(argc, argv);
	if (Pacfileonly)
	{
		ParseFASTA = TRUE;
		BuildBWT = FALSE;
		BuildSaValue = FALSE;
		BuildSaIndex = FALSE;
		BuildCachedSaIndex = FALSE;
	}
	PrintShortDesc();

	// Ini
	if (strcmp(argv[0] + strlen(argv[0]) - 4, ".exe") == 0) {
		*(argv[0] + strlen(argv[0]) - 4) = '\0';
	}
	sprintf(filename, "%s.ini", argv[0]);
	ParseIniFile(filename);
	printf("\n");
	ProcessIni();
	ValidateIni();
	PrintIni();

	if (Confirmation == TRUE) {
		printf("Press Y to go or N to cancel. ");
        c = (char)getchar();
		while (c != 'y' && c != 'Y' && c != 'n' && c!= 'N') {
			c = (char)getchar();
		}
		if (c == 'n' || c == 'N') {
			exit(0);
		}
	}

	startTime = setStartTime();

	MMMasterInitialize(1, 0, FALSE, NULL);
	mmPool = MMPoolCreate(PoolSize);

	// Parse FASTA file to produce packed DNA and annotation file
	if (ParseFASTA == TRUE) {

		printf("Parsing FASTA file..\n");

		numSeq = HSPParseFASTAToPacked(FASTAFileName, AnnotationFileName, PackedDNAFileName, AmbiguityFileName, FASTARandomSeed, MaskLowerCase);

		printf("Finished. Parsed %u sequences.\n", numSeq);

		elapsedTime = getElapsedTime(startTime) - totalElapsedTime;
		printf("Elapsed time = ");
		printElapsedTime(stdout, FALSE, FALSE, TRUE, 2, elapsedTime);
		totalElapsedTime += elapsedTime;
		printf("\n");

	}

	// Construct BWTInc from text
	if (BuildBWT == TRUE) {

		printf("Building BWT..\n");

		bwtInc = BWTIncConstructFromPacked(mmPool, PackedDNAFileName, ShowProgress, 
										   TargetNBit, InitialMaxBuildSize, IncMaxBuildSize);

		printf("Finished constructing BWT in %u iterations.  ", bwtInc->numberOfIterationDone);
		elapsedTime = getElapsedTime(startTime) - totalElapsedTime;
		printf("Elapsed time = ");
		printElapsedTime(stdout, FALSE, FALSE, TRUE, 2, elapsedTime);
		totalElapsedTime += elapsedTime;
		printf("\n");

		printf("Saving BWT..\n");
		BWTSaveBwtCodeAndOcc(bwtInc->bwt, BWTCodeFileName, BWTOccValueFileName);
		printf("Finished saving BWT.  ");
		elapsedTime = getElapsedTime(startTime) - totalElapsedTime;
		printf("Elapsed time = ");
		printElapsedTime(stdout, FALSE, FALSE, TRUE, 2, elapsedTime);
		totalElapsedTime += elapsedTime;
		printf("\n");

		textLength = bwtInc->bwt->textLength;

		BWTIncFree(mmPool, bwtInc);

	}

	// Load BWT
	if (BuildSaValue || BuildSaIndex) {

		printf("Loading BWT...\n");

		bwt = BWTLoad(mmPool, BWTCodeFileName, BWTOccValueFileName, NULL, NULL, NULL, NULL);

		printf("Finished loading BWT.  ");

		elapsedTime = getElapsedTime(startTime) - totalElapsedTime;
		printf("Elapsed time = ");
		printElapsedTime(stdout, FALSE, FALSE, TRUE, 2, elapsedTime);
		totalElapsedTime += elapsedTime;
		printf("\n");

		textLength = bwt->textLength;

	}

	if (BuildSaValue) {

		printf("Building SA value...\n");
		
		if (ShowProgress) {
			BWTGenerateSaValue(mmPool, bwt, SaValueFreq, bwt->textLength / SaValueFreq / 10);
		} else {
			BWTGenerateSaValue(mmPool, bwt, SaValueFreq, 0);
		}
		BWTSaveSaValue(bwt, SaValueFileName);

		printf("Finished building SA value.  ");

		elapsedTime = getElapsedTime(startTime) - totalElapsedTime;
		printf("Elapsed time = ");
		printElapsedTime(stdout, FALSE, FALSE, TRUE, 2, elapsedTime);
		totalElapsedTime += elapsedTime;
		printf("\n");

	}

	if (BuildSaIndex) {

		printf("Building SA index...\n");
		
		BWTGenerateSaRangeTable(bwt, SaIndexNumOfChar, SaIndexFileName);

		printf("Finished building SA index.  ");

		elapsedTime = getElapsedTime(startTime) - totalElapsedTime;
		printf("Elapsed time = ");
		printElapsedTime(stdout, FALSE, FALSE, TRUE, 2, elapsedTime);
		totalElapsedTime += elapsedTime;
		printf("\n");

	}

	// Free BWT
	if (BuildSaValue || BuildSaIndex) {
		BWTFree(mmPool, bwt);
	}

	// Finished all construction tasks
	printf("Finished all tasks.  ");
	totalElapsedTime = getElapsedTime(startTime);
	printf("Total elapsed time = ");
	printElapsedTime(stdout, FALSE, FALSE, TRUE, 2, totalElapsedTime);
	printf("\n");

	MMMasterPrintReport(stdout, FALSE, FALSE, FALSE);
	if (BuildSaValue) {
		fprintf(stdout, "Number of char   :  %u\n", textLength);
		fprintf(stdout, "Bit per char     :  %.2f\n", (float)MMMasterMaxTotalByteDispatched() * BITS_IN_BYTE / textLength);
		printf("\n");
	}

	MMPoolFree(mmPool);

	iniparser_freedict(programInput);

	return 0;

}

dictionary *ParseInput(int argc, char** argv) {

	dictionary *programInput;
	char t1[3] = "-c";	// specify that this is a boolean type parameter
	char t2[3] = "-U";	// specify that this is a boolean type parameter
	char t3[3] = "-p";	// specify that this is a boolean type parameter
	char *d[3];

	d[0] = t1;
	d[1] = t2;
	d[2] = t3;
	
	programInput = paraparser_load(argc, argv, 3, d);	// 2 boolean type parameters

	// Get database name
	if (!iniparser_find_entry(programInput, "argument:1")) {
		PrintHelp();
		exit(1);
	}
	iniparser_copystring(programInput, "argument:1", DatabaseName, DatabaseName, MAX_FILENAME_LEN);
	if (strlen(DatabaseName) + 4 > MAX_FILENAME_LEN) {
		PrintHelp();
		exit(1);
	}

	// Get FASTA file name
	iniparser_copystring(programInput, "argument:2", FASTAFileName, DatabaseName, MAX_FILENAME_LEN);
	if (strlen(FASTAFileName) > MAX_FILENAME_LEN) {
		PrintHelp();
		exit(1);
	}


	// Whether confirmation is needed
	Confirmation = iniparser_find_entry(programInput, "parameter:-c");

	MaskLowerCase = iniparser_find_entry(programInput, "parameter:-U");

	Pacfileonly = iniparser_find_entry(programInput, "parameter:-p");

	return programInput;

}

void ParseIniFile(char *iniFileName) {

	dictionary *ini;

	printf("Loading %s ..", iniFileName);
	ini = iniparser_load(iniFileName, FALSE);
	if (ini == NULL) {
		printf("not found.\n");
		return;
	}
	printf("done.\n");

	// BuildTasks parameters
	ParseFASTA = iniparser_getboolean(ini, "BuildTasks:ParseFASTA", ParseFASTA);
	BuildBWT = iniparser_getboolean(ini, "BuildTasks:BuildBWT", BuildBWT);
	BuildSaValue = iniparser_getboolean(ini, "BuildTasks:BuildSaValue", BuildSaValue);
	BuildSaIndex = iniparser_getboolean(ini, "BuildTasks:BuildSaIndex", BuildSaIndex);

	// Display parameters
	ShowProgress = iniparser_getboolean(ini, "Display:ShowProgress", ShowProgress);

	// Parse FASTA parameters
	FASTARandomSeed = iniparser_getint(ini, "ParseFASTA:RandomSeed", FASTARandomSeed);
	if (FASTARandomSeed == 0) {
		FASTARandomSeed = getRandomSeed();
	}

	// Build BWT parameters
	OccValueFreq = iniparser_getint(ini, "BuildBWT:OccValueFreq", OccValueFreq);
	TargetNBit = (float)iniparser_getdouble(ini, "BuildBWT:TargetNBit", TargetNBit);
	InitialMaxBuildSize = iniparser_getint(ini, "BuildBWT:InitialMaxBuildSize", InitialMaxBuildSize);
	IncMaxBuildSize = iniparser_getint(ini, "BuildBWT:IncMaxBuildSize", IncMaxBuildSize);

	// Build SA value parameters
	SaValueFreq = iniparser_getint(ini, "BuildSAValue:SaValueFreq", SaValueFreq);

	// Build SA index parameters
	SaIndexNumOfChar = iniparser_getint(ini, "BuildSAIndex:SaIndexNumOfChar", SaIndexNumOfChar);

	// Database parameters
	iniparser_copystring(ini, "Database:AnnotationFileName", AnnotationFileName, AnnotationFileName, MAX_FILENAME_LEN);
	iniparser_copystring(ini, "Database:AmbiguityFileName", AmbiguityFileName, AmbiguityFileName, MAX_FILENAME_LEN);
	iniparser_copystring(ini, "Database:PackedDNAFileName", PackedDNAFileName, PackedDNAFileName, MAX_FILENAME_LEN);
	iniparser_copystring(ini, "Database:BWTCodeFileName", BWTCodeFileName, BWTCodeFileName, MAX_FILENAME_LEN);
	iniparser_copystring(ini, "Database:BWTOccValueFileName", BWTOccValueFileName, BWTOccValueFileName, MAX_FILENAME_LEN);
	iniparser_copystring(ini, "Database:SaValueFileName", SaValueFileName, SaValueFileName, MAX_FILENAME_LEN);
	iniparser_copystring(ini, "Database:SaIndexFileName", SaIndexFileName, SaIndexFileName, MAX_FILENAME_LEN);

	iniparser_freedict(ini);

}

void ProcessIni() {

	ProcessFileName(AnnotationFileName, AnnotationFileName, DatabaseName);
	ProcessFileName(AmbiguityFileName, AmbiguityFileName, DatabaseName);
	ProcessFileName(PackedDNAFileName, PackedDNAFileName, DatabaseName);
	ProcessFileName(BWTCodeFileName, BWTCodeFileName, DatabaseName);
	ProcessFileName(BWTOccValueFileName, BWTOccValueFileName, DatabaseName);
	ProcessFileName(SaValueFileName, SaValueFileName, DatabaseName);
	ProcessFileName(SaIndexFileName, SaIndexFileName, DatabaseName);

}

void ValidateIni() {

	if (!ParseFASTA && !BuildBWT && !BuildSaValue && !BuildSaIndex) {
		fprintf(stderr, "No action is specified!\n");
		exit(1);
	}
	if (ParseFASTA) {
		if (PackedDNAFileName[0] == '\0') {
			fprintf(stderr, "Packed DNA file name is not specified!\n");
			exit(1);
		}
		if (AnnotationFileName[0] == '\0') {
			fprintf(stderr, "Annotation file name is not specified!\n");
			exit(1);
		}
		if (AmbiguityFileName[0] == '\0') {
			fprintf(stderr, "Ambiguity file name is not specified!\n");
			exit(1);
		}
	}
	if (BuildBWT) {
		if (PackedDNAFileName[0] == '\0') {
			fprintf(stderr, "Packed DNA file is not specified!\n");
			exit(1);
		}
		if (BWTCodeFileName[0] == '\0') {
			fprintf(stderr, "BWT code file name is not specified!\n");
			exit(1);
		}
		if (BWTOccValueFileName[0] == '\0') {
			fprintf(stderr, "BWT Occ value file name is not specified!\n");
			exit(1);
		}
		if (TargetNBit < 2.5) {
			fprintf(stderr, "Target NBit should be at least 2.5!\n");
			exit(1);
		}
	}
	if (BuildSaValue) {
		if (BWTCodeFileName[0] == '\0') {
			fprintf(stderr, "BWT code file is not specified!\n");
			exit(1);
		}
		if (BWTOccValueFileName[0] == '\0') {
			fprintf(stderr, "BWT Occ value file is not specified!\n");
			exit(1);
		}
		if (SaValueFileName[0] == '\0') {
			fprintf(stderr, "SA value file name is not specified!\n");
			exit(1);
		}
		if (SaValueFreq <= 0) {
			fprintf(stderr, "SA value frequency must > 0!\n");
			exit(1);
		}
	}

	if (BuildSaIndex) {
		if (BWTCodeFileName[0] == '\0') {
			fprintf(stderr, "BWT code file is not specified!\n");
			exit(1);
		}
		if (BWTOccValueFileName[0] == '\0') {
			fprintf(stderr, "BWT Occ value file is not specified!\n");
			exit(1);
		}
		if (SaIndexFileName[0] == '\0') {
			fprintf(stderr, "SA index file name is not specified!\n");
			exit(1);
		}
		if (SaIndexNumOfChar <= 0) {
			fprintf(stderr, "SA index number of character must > 0!\n");
			exit(1);
		}
		if (SaIndexNumOfChar > 13) {
			fprintf(stderr, "SA index number of character must <= 13!\n");
			exit(1);
		}
	}

}


void PrintIni() {

	char boolean[2];

	boolean[0] = 'N';
	boolean[1] = 'Y';

	printf("Parse FASTA file    : %c\n", boolean[ParseFASTA]);
	printf("Build BWT           : %c\n", boolean[BuildBWT]);
	printf("Build SA value      : %c\n", boolean[BuildSaValue]);
	printf("Build SA index      : %c\n", boolean[BuildSaIndex]);
	printf("\n");

	printf("Show progress       : %c\n", boolean[ShowProgress]);
	printf("\n");

	if (ParseFASTA) {
		printf("Parse FASTA :\n");
		printf("Mask lower case         : %c\n", boolean[MaskLowerCase]);
		printf("Random seed             : %u\n", FASTARandomSeed);
		printf("\n");
	}

	if (BuildBWT) {
		printf("Build BWT :\n");
		printf("Target N Bits           : %.2f\n", TargetNBit);
		printf("Occ value frequency     : %u\n", OccValueFreq);
		printf("Initial Max Build Size  : %u    Inc Max Build Size : %u\n", 
				InitialMaxBuildSize, IncMaxBuildSize);
		printf("\n");
	}

	if (BuildSaValue) {
		printf("Build SA value :\n");
		printf("SA value frequency      : %u\n", SaValueFreq);
		printf("\n");
	}

	if (BuildSaIndex) {
		printf("Build SA index :\n");
		printf("SA index no. of char    : %u\n", SaIndexNumOfChar);
		printf("\n");
	}

	printf("Annotation file          : %s\n", AnnotationFileName);
	printf("Ambigurity file          : %s\n", AmbiguityFileName);
	printf("Packed DNA file          : %s\n", PackedDNAFileName);
	printf("BWT Code file            : %s\n", BWTCodeFileName);
	printf("BWT Occ value file       : %s\n", BWTOccValueFileName);
	printf("SA value file            : %s\n", SaValueFileName);
	printf("SA index file            : %s\n", SaIndexFileName);
	printf("\n");

}

void PrintShortDesc() {

	printf("BWTFormatdb v1.0, Copyright (C) 2006, Wong Chi Kwong.\n");
	printf("BWTFormatdb comes with ABSOLUTELY NO WARRENTY.\n");
	printf("BWTFormatdb is free software, and you are welcome to\n");
	printf("redistribute it under certain conditions.\n");
	printf("For details type BWTFormatdb.\n");
	printf("\n");

}

void PrintHelp() {

	printf("BWTFormatdb v1.0, Copyright (C) 2006, Wong Chi Kwong.\n");
	printf("\n");

	printf("This program is free software; you can redistribute it and/or\n");
	printf("modify it under the terms of the GNU General Public License\n");
	printf("as published by the Free Software Foundation; either version 2\n");
	printf("of the License, or (at your option) any later version.\n");
	printf("\n");

	printf("This program is distributed in the hope that it will be useful,\n");
	printf("but WITHOUT ANY WARRANTY; without even the implied warranty of\n");
	printf("MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n");
	printf("GNU General Public License for more details.\n");
	printf("\n");

	printf("You should have received a copy of the GNU General Public License\n");
	printf("along with this program; if not, write to the Free Software\n");
	printf("Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.\n");
	printf("\n");

	printf("Syntax: BWTFormatdb <database> [FASTA file; default = <database>]\n");
	printf("                               [-c Confirm]\n");
	printf("                               [-U Mask lower case]\n");

}

void ProcessFileName(char *outputFileName, const char *inputFileName, const char *databaseName) {

	char tempChar[MAX_FILENAME_LEN];
	unsigned int i;

	if (inputFileName == NULL) {
		if (outputFileName != inputFileName) {
			outputFileName[0] = '\0';
		}
		return;
	}

	if (strlen(databaseName) + strlen(inputFileName) > MAX_FILENAME_LEN) {
		fprintf(stderr, "File length is too long!\n");
		exit(1);
	}

	strncpy(tempChar, inputFileName, MAX_FILENAME_LEN);

	// locate the *
	for (i=0; i<MAX_FILENAME_LEN; i++) {
		if (tempChar[i] == '*') {
			break;
		}
	}
	if (i<MAX_FILENAME_LEN) {
		tempChar[i] = '\0';
		sprintf(outputFileName, "%s%s%s", tempChar, databaseName, tempChar + i + 1);
	} else {
		sprintf(outputFileName, "%s", tempChar);
	}

}
#endif
