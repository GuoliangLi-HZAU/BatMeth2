#include "config.h"
#ifdef MMX
/*

   HSP.c		BWTBlastn functions

   This module contains miscellaneous BWTBlastn functions.

   Copyright (C) 2004, Wong Chi Kwong.

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
#include <string.h>
#include <math.h>
#include <stdint.h>
#include <emmintrin.h>
#include <mmintrin.h>
#include "TextConverter.h"
#include "MiscUtilities.h"
#include "Socket.h"
#include "r250.h"
#include "HSP.h"
#include "HSPstatistic.h"


extern  double stat_expectationValue;

void HSPFillCharMap(unsigned char charMap[255]) {

	int i;

	for (i=0; i<255; i++) {
		charMap[i] = nonMatchDnaCharIndex;
	}
	for (i=0; i<16; i++) {
		charMap[dnaChar[i]] = (unsigned char)i;
		charMap[dnaChar[i] - 'A' + 'a'] = (unsigned char)i;
	}

}

void HSPFillComplementMap(unsigned char complementMap[255]) {

	int i;

	for (i=0; i<255; i++) {
		complementMap[i] = dnaChar[nonMatchDnaCharIndex];
	}
	for (i=0; i<16; i++) {
		complementMap[dnaComplement[i]] = dnaChar[i];
		complementMap[dnaComplement[i] - 'A' + 'a'] = dnaChar[i] - 'A' + 'a';
	}

}

void HSPFillScoringMatrix(int scoringMatrix[16][16], const int matchScore, const int mismatchScore, const int leftShift) {

	int i, j, k, l;

	for (i=0; i<16; i++) {
		for (j=0; j<16; j++) {
			if (ambiguityCount[i] > 0 && ambiguityCount[j] > 0) {
				scoringMatrix[i][j] = 0;
				for (k=0; k<ambiguityCount[i]; k++) {
					for (l=0; l<ambiguityCount[j]; l++) {
						if (ambiguityMatch[i][k] == ambiguityMatch[j][l]) {
							// match
							scoringMatrix[i][j] += matchScore;
						} else {
							// mismatch
							scoringMatrix[i][j] += mismatchScore;
						}
					}
				}
				if (scoringMatrix[i][j] > 0) {
					scoringMatrix[i][j] = (scoringMatrix[i][j] + ambiguityCount[i] * ambiguityCount[j]  - 1) / (ambiguityCount[i] * ambiguityCount[j]);
				} else if (scoringMatrix[i][j] < 0) {
					scoringMatrix[i][j] = -((-scoringMatrix[i][j] + ambiguityCount[i] * ambiguityCount[j] - 1) / (ambiguityCount[i] * ambiguityCount[j]));
				} else {
					scoringMatrix[i][j] = 0;
				}
			} else {
				// Character meant to match nothing
				scoringMatrix[i][j] = mismatchScore;
			}
			scoringMatrix[i][j] = scoringMatrix[i][j] * (1 << leftShift);
		}
	}


}

HSP *HSPLoad(MMPool *mmPool, const char *PackedDNAFileName, const char *AnnotationFileName, const char *AmbiguityFileName, const unsigned int trailerBufferInWord) {

	HSP *hsp;
	unsigned int dnaLength;
	unsigned int randomSeed;
	int i;
	char c;
	unsigned char charMap[255];

	FILE *annotationFile = NULL, *ambiguityFile = NULL;

	hsp = MMPoolDispatch(mmPool, sizeof(HSP));

	// Load packed DNA
	if (PackedDNAFileName != NULL && PackedDNAFileName[0] != '\0' && PackedDNAFileName[0] != '-') {
		hsp->packedDNA = DNALoadPacked(PackedDNAFileName, &hsp->dnaLength, TRUE, trailerBufferInWord);
	} else {
		hsp->packedDNA = NULL;
		hsp->dnaLength = 0;
	}

	// Load annotation
	if (AnnotationFileName != NULL && AnnotationFileName[0] != '\0' && AnnotationFileName[0] != '-') {

		annotationFile = (FILE*)fopen64(AnnotationFileName, "r");
		if (annotationFile == NULL) {
			fprintf(stderr, "Cannot open annotation file!\n");
			exit(1);
		}

		fscanf(annotationFile, "%u %d %u\n", &dnaLength, &hsp->numOfSeq, &randomSeed);
		if (hsp->dnaLength != 0 && dnaLength != hsp->dnaLength) {
			fprintf(stderr, "Annotation database length not match!\n");
			exit(1);
		}
		hsp->dnaLength = dnaLength;
		if (hsp->numOfSeq == 0) {
			fprintf(stderr, "Annotation number of sequence = 0!\n");
			exit(1);
		}
		hsp->annotation = MMUnitAllocate((hsp->numOfSeq+1) * sizeof(Annotation));
		hsp->seqOffset = MMUnitAllocate((hsp->numOfSeq+1) * sizeof(SeqOffset));

		i = 0;
		hsp->minSeqLength = UINT32_MAX;
		while (!feof(annotationFile) && i < hsp->numOfSeq) {
			fscanf(annotationFile, "%u ", &hsp->annotation[i].gi);
			fgets(hsp->annotation[i].text, MAX_SEQ_NAME_LENGTH, annotationFile);
			fscanf(annotationFile, "%u %u %d\n", &hsp->seqOffset[i].startPos, &hsp->seqOffset[i].endPos, &hsp->seqOffset[i+1].firstAmbiguityIndex);
			hsp->seqOffset[i].lastAmbiguityIndex = hsp->seqOffset[i+1].firstAmbiguityIndex;
			if (hsp->seqOffset[i].endPos < hsp->minSeqLength) {
				hsp->minSeqLength = hsp->seqOffset[i].endPos;
			}
			hsp->seqOffset[i].endPos = hsp->seqOffset[i].startPos + hsp->seqOffset[i].endPos - 1;	// length of sequence is stored
			i++;
		}
		if (i < hsp->numOfSeq) {
			fprintf(stderr, "Annotation missing entries!\n");
			exit(1);
		}
		fclose(annotationFile);

		hsp->annotation[i].gi = 0;
		hsp->annotation[i].text[0] = '\0';

		hsp->seqOffset[i].startPos = UINT32_MAX;
		hsp->seqOffset[i].endPos = UINT32_MAX;
		hsp->seqOffset[0].firstAmbiguityIndex = 1;	// ambiguity[0] and ambiguity[numOfAmbiguity+1] are dummy
		for (i=1; i<=hsp->numOfSeq; i++) {
			hsp->seqOffset[i].firstAmbiguityIndex += hsp->seqOffset[i-1].firstAmbiguityIndex;	// number of ambiguity is stored
		}
		// hsp->seqOffset[hsp->numOfSeq].firstAmbiguityIndex = total number of ambiguity + 1 now
		hsp->seqOffset[hsp->numOfSeq].lastAmbiguityIndex = hsp->seqOffset[hsp->numOfSeq].firstAmbiguityIndex - 1;
		// hsp->seqOffset[hsp->numOfSeq].lastAmbiguityIndex = total number of ambiguity now
		for (i=hsp->numOfSeq; i>1; i--) {
			hsp->seqOffset[i-1].lastAmbiguityIndex = hsp->seqOffset[i].lastAmbiguityIndex - hsp->seqOffset[i-1].lastAmbiguityIndex;	// number of ambiguity is stored
		}
		for (i=0; i<hsp->numOfSeq; i++) {
			hsp->seqOffset[i].lastAmbiguityIndex = hsp->seqOffset[i+1].lastAmbiguityIndex;
		}

	} else {

		// Make up a dummy annotation and a dummy sequence offset
		hsp->annotation = MMUnitAllocate((1+1) * sizeof(Annotation));
		hsp->seqOffset = MMUnitAllocate((1+1) * sizeof(SeqOffset));

		hsp->numOfSeq = 1;
		hsp->numOfAmbiguity = 0;
		
		hsp->annotation[0].gi = 0;
		hsp->annotation[0].text[0] = '\0';
		hsp->annotation[1].gi = 0;
		hsp->annotation[1].text[0] = '\0';
		hsp->seqOffset[0].startPos = 0;
		hsp->seqOffset[0].endPos = hsp->dnaLength - 1;
		hsp->seqOffset[0].firstAmbiguityIndex = 1;
		hsp->seqOffset[0].lastAmbiguityIndex = 0;
		hsp->seqOffset[1].startPos = UINT32_MAX;
		hsp->seqOffset[1].endPos = UINT32_MAX;
		hsp->seqOffset[1].firstAmbiguityIndex = UINT32_MAX;	// should not be referred
		hsp->seqOffset[1].lastAmbiguityIndex = UINT32_MAX;	// should not be referred

	}

	// Load ambigity
	if (AmbiguityFileName != NULL && AmbiguityFileName[0] != '\0' && AmbiguityFileName[0] != '-') {

		// Load ambigity
		ambiguityFile = (FILE*)fopen64(AmbiguityFileName, "r");
		if (ambiguityFile == NULL) {
			fprintf(stderr, "Cannot open ambiguity file!\n");
			exit(1);
		}

		fscanf(ambiguityFile, "%u %d %d\n", &dnaLength, &i, &hsp->numOfAmbiguity);
		if (hsp->dnaLength != 0 && dnaLength != hsp->dnaLength) {
			fprintf(stderr, "Ambiguity database length not match!\n");
			exit(1);
		}
		hsp->dnaLength = dnaLength;
		if (i != hsp->numOfSeq) {
			fprintf(stderr, "Ambiguity database number of sequence not match!\n");
			exit(1);
		}
		if (hsp->numOfAmbiguity != hsp->seqOffset[hsp->numOfSeq].firstAmbiguityIndex - 1) {
			fprintf(stderr, "Ambiguity number of ambiguity not match!\n");
			exit(1);
		}

		HSPFillCharMap(charMap);

		hsp->ambiguity = MMUnitAllocate((hsp->numOfAmbiguity + 2) * sizeof(Ambiguity));
		hsp->ambiguity[0].startPos = 0;
		hsp->ambiguity[0].rightOfEndPos = 0;
		hsp->ambiguity[0].symbol = 0;
		i = 1;
		while (!feof(ambiguityFile) && i-1 < hsp->numOfAmbiguity) {
			fscanf(ambiguityFile, "%u %u %c\n", &hsp->ambiguity[i].startPos, &hsp->ambiguity[i].rightOfEndPos, &c);
			hsp->ambiguity[i].rightOfEndPos = hsp->ambiguity[i].startPos + hsp->ambiguity[i].rightOfEndPos;	// number of character is stored
			hsp->ambiguity[i].symbol = charMap[c];
			i++;
		}
		hsp->ambiguity[i].startPos = UINT32_MAX;
		hsp->ambiguity[i].rightOfEndPos = UINT32_MAX;
		hsp->ambiguity[i].symbol = 0;

		if (i-1 < hsp->numOfAmbiguity) {
			fprintf(stderr, "Ambiguity missing entries!\n");
			exit(1);
		}
		fclose(ambiguityFile);

	} else {

		if (hsp->numOfAmbiguity > 0) {
			fprintf(stderr, "Ambiguity file missing!\n");
			exit(1);
		}

		hsp->ambiguity = MMUnitAllocate((0 + 2) * sizeof(Ambiguity));
		hsp->ambiguity[0].startPos = 0;
		hsp->ambiguity[0].rightOfEndPos = 0;
		hsp->ambiguity[0].symbol = 0;
		hsp->ambiguity[1].startPos = UINT32_MAX;
		hsp->ambiguity[1].rightOfEndPos = UINT32_MAX;
		hsp->ambiguity[1].symbol = 0;

	}

	return hsp;

}

HSP *HSPConvertFromText(MMPool *mmPool, const unsigned char *text, const unsigned int textLength, 
						const unsigned int FASTARandomSeed, const int maskLowerCase,
						const int gi, const char *seqName) {

	HSP *hsp;
	Ambiguity *tempAmbiguity;
	int ambiguityAllocated =  256;
	
	char c;
	unsigned int i;
	int numAmbiguity;
	unsigned int lastAmbiguityPos;
	unsigned int numCharInBuffer, numOfConversion;
	unsigned int sequenceRandomSeed;
	unsigned char charMap[255];

	unsigned char buffer[PACKED_BUFFER_SIZE];

	HSPFillCharMap(charMap);

	hsp = MMPoolDispatch(mmPool, sizeof(HSP));

	hsp->dnaLength = textLength;
	hsp->minSeqLength = textLength;
	hsp->numOfSeq = 1;

	hsp->annotation = MMPoolDispatch(mmPool, (1+1) * sizeof(Annotation));
	hsp->annotation[0].gi = gi;
	strncpy(hsp->annotation[0].text, seqName, MAX_SEQ_NAME_LENGTH);
	hsp->annotation[1].gi = 0;
	hsp->annotation[1].text[0] = '\0';

	hsp->seqOffset = MMPoolDispatch(mmPool, (1+1) * sizeof(SeqOffset));
	hsp->seqOffset[0].startPos = 0;
	hsp->seqOffset[0].endPos = textLength - 1;
	hsp->seqOffset[0].firstAmbiguityIndex = 1;
	hsp->seqOffset[0].lastAmbiguityIndex = 0;
	hsp->seqOffset[1].startPos = UINT32_MAX;
	hsp->seqOffset[1].endPos = UINT32_MAX;
	hsp->seqOffset[1].firstAmbiguityIndex = UINT32_MAX;	// should not be referred
	hsp->seqOffset[1].lastAmbiguityIndex = UINT32_MAX;	// should not be referred

	hsp->packedDNA = MMPoolDispatch(mmPool, (textLength / CHAR_PER_WORD + 1) * sizeof(int));
	
	tempAmbiguity = MMUnitAllocate(ambiguityAllocated * sizeof(Ambiguity));

	HSPFillCharMap(charMap);

	// Set random seed for the sequence
	sequenceRandomSeed = FASTARandomSeed;
	if (gi > 0) {
		sequenceRandomSeed += gi;
	} else {
		for (i=0; i<(int)strlen(seqName); i++) {
			c = seqName[i];
			sequenceRandomSeed += c;
		}
	}
	r250_init(sequenceRandomSeed);

	numAmbiguity = -1;
	numCharInBuffer = 0;
	numOfConversion = 0;
	lastAmbiguityPos = (unsigned int)-2;

	for (i=0; i<textLength; i++) {

		c = text[i];

		if (maskLowerCase && c >= 'a' && c <= 'z') {
			c = dnaChar[lowercaseDnaCharIndex];
		}
		if (ambiguityCount[charMap[c]] == 1) {
			buffer[numCharInBuffer] = c;
		} else {
			c = charMap[c];
			if (ambiguityCount[c] > 0) {
				buffer[numCharInBuffer] = dnaChar[ambiguityMatch[c][r250() % ambiguityCount[c]]];
			} else {
				buffer[numCharInBuffer] = dnaChar[ambiguityMatch[c][r250() % ALPHABET_SIZE]];
			}
			if (i == lastAmbiguityPos + 1 && c == (char)tempAmbiguity[numAmbiguity].symbol) {
				tempAmbiguity[numAmbiguity].rightOfEndPos++;
			} else {
				numAmbiguity++;
				if (numAmbiguity >= ambiguityAllocated) {
					tempAmbiguity = MMUnitReallocate(tempAmbiguity, sizeof(Ambiguity) * ambiguityAllocated * 2, sizeof(Ambiguity) * ambiguityAllocated);
					ambiguityAllocated *= 2;
				}
				tempAmbiguity[numAmbiguity].startPos = i;
				tempAmbiguity[numAmbiguity].rightOfEndPos = i+1;
				tempAmbiguity[numAmbiguity].symbol = c;
			}
			lastAmbiguityPos = i;
		}
		numCharInBuffer++;
		if (numCharInBuffer >= PACKED_BUFFER_SIZE) {
			ConvertTextToWordPacked(buffer, hsp->packedDNA + numOfConversion * PACKED_BUFFER_SIZE / CHAR_PER_WORD, charMap, 4, PACKED_BUFFER_SIZE);
			numCharInBuffer = 0;
			numOfConversion++;
		}
	}
	if (numCharInBuffer > 0) {
		ConvertTextToWordPacked(buffer, hsp->packedDNA + numOfConversion * PACKED_BUFFER_SIZE / CHAR_PER_WORD, charMap, 4, numCharInBuffer);
	}

	numAmbiguity++;
	hsp->numOfAmbiguity = numAmbiguity;
	hsp->ambiguity = MMPoolDispatch(mmPool, (hsp->numOfAmbiguity + 2) * sizeof(Ambiguity));
	if (hsp->numOfAmbiguity > 0) {
		memcpy(hsp->ambiguity + 1, tempAmbiguity, hsp->numOfAmbiguity * sizeof(Ambiguity));
	}
	hsp->ambiguity[0].startPos = 0;
	hsp->ambiguity[0].rightOfEndPos = 0;
	hsp->ambiguity[0].symbol = 0;
	hsp->ambiguity[hsp->numOfAmbiguity+1].startPos = UINT32_MAX;
	hsp->ambiguity[hsp->numOfAmbiguity+1].rightOfEndPos = UINT32_MAX;
	hsp->ambiguity[hsp->numOfAmbiguity+1].symbol = 0;

	hsp->seqOffset[0].firstAmbiguityIndex = 1;
	hsp->seqOffset[0].lastAmbiguityIndex = numAmbiguity;

	MMUnitFree(tempAmbiguity, sizeof(Ambiguity) * ambiguityAllocated);

	return hsp;

}

void HSPFree(MMPool *mmPool, HSP *hsp, const unsigned int trailerBufferInWord) {

	if (hsp->packedDNA != NULL) {
		DNAFreePacked(hsp->packedDNA, hsp->dnaLength, trailerBufferInWord);
	}
	MMUnitFree(hsp->seqOffset, (hsp->numOfSeq+1) * sizeof(SeqOffset));
	MMUnitFree(hsp->annotation, (hsp->numOfSeq+1) * sizeof(Annotation));
	MMUnitFree(hsp->ambiguity, (hsp->numOfAmbiguity+2) * sizeof(Ambiguity));

	MMPoolReturn(mmPool, hsp, sizeof(hsp));

}

unsigned int HSPParseFASTAToPacked(const char* FASTAFileName, const char* annotationFileName, const char* packedDNAFileName, const char* ambiguityFileName,
					  const unsigned int FASTARandomSeed, const int maskLowerCase) {

	FILE *FASTAFile, *annotationFile, *packedDNAFile, *ambiguityFile;
	Annotation *annotation;
	SeqOffset *seqOffset;
	Ambiguity *ambiguity;
	int annotationAllocated = 256;
	int ambiguityAllocated =  256;
	
	char c;
	int i;
	int numAmbiguity;
	unsigned int lastAmbiguityPos;
	unsigned int numChar, numCharInBuffer, totalNumChar;
	unsigned int sequenceRandomSeed;
	int numSeq;
	unsigned char charMap[255];

	unsigned char buffer[PACKED_BUFFER_SIZE];
	unsigned char packedBuffer[PACKED_BUFFER_SIZE / 4];

	FASTAFile = (FILE*)fopen64(FASTAFileName, "r");
	if (FASTAFile == NULL) {
		fprintf(stderr, "ParseFASTToPacked() : Cannot open FASTAFileName!\n");
		exit(1);
	}

	annotationFile = (FILE*)fopen64(annotationFileName, "w");
	if (annotationFile == NULL) {
		fprintf(stderr, "ParseFASTToPacked() : Cannot open annotationFileName!\n");
		exit(1);
	}

	packedDNAFile = (FILE*)fopen64(packedDNAFileName, "wb");
	if (packedDNAFile == NULL) {
		fprintf(stderr, "ParseFASTToPacked() : Cannot open packedDNAFileName!\n");
		exit(1);
	}

	ambiguityFile = (FILE*)fopen64(ambiguityFileName, "w");
	if (ambiguityFile == NULL) {
		fprintf(stderr, "ParseFASTToPacked() : Cannot open ambiguityFileName!\n");
		exit(1);
	}

	HSPFillCharMap(charMap);

	c = (char)getc(FASTAFile);
	if (c != '>') {
		fprintf(stderr, "ParseFASTToPacked() : FASTA file does not begin with '>'!\n");
		exit(1);
	}

	totalNumChar = 0;
	numSeq = 0;
	numAmbiguity = -1;
	numCharInBuffer = 0;

	annotation = MMUnitAllocate(sizeof(Annotation) * annotationAllocated);
	seqOffset = MMUnitAllocate(sizeof(SeqOffset) * annotationAllocated);
	ambiguity = MMUnitAllocate(sizeof(Ambiguity) * ambiguityAllocated);

	while (!feof(FASTAFile)) {

		numChar = 0;
		if (numSeq >= annotationAllocated) {
			annotation = MMUnitReallocate(annotation, sizeof(Annotation) * annotationAllocated * 2, sizeof(Annotation) * annotationAllocated);
			seqOffset = MMUnitReallocate(seqOffset, sizeof(SeqOffset) * annotationAllocated * 2, sizeof(SeqOffset) * annotationAllocated);
			annotationAllocated *= 2;
		}

		annotation[numSeq].gi = 0;

		c = (char)getc(FASTAFile);
		while (!feof(FASTAFile) && c != '\n') {
			if (numChar < MAX_SEQ_NAME_LENGTH) {
				annotation[numSeq].text[numChar] = c;
				numChar++;
			}
			if (numChar == 3 && annotation[numSeq].text[0] == 'g' 
							 && annotation[numSeq].text[1] == 'i' 
							 && annotation[numSeq].text[2] == '|') {
				fscanf(FASTAFile, "%u|", &(annotation[numSeq].gi));
				numChar = 0;
			}
			c = (char)getc(FASTAFile);
		}
		annotation[numSeq].text[numChar] = '\0';

		// Set random seed for the sequence
		sequenceRandomSeed = FASTARandomSeed;
		if (annotation[numSeq].gi > 0) {
			sequenceRandomSeed += annotation[numSeq].gi;
		} else {
			for (i=0; i<(int)numChar; i++) {
				c = annotation[numSeq].text[i];
				sequenceRandomSeed += c;
			}
		}
		r250_init(sequenceRandomSeed);

		seqOffset[numSeq].startPos = totalNumChar;
		numChar = 0;
		lastAmbiguityPos = (unsigned int)-2;

		c = (char)getc(FASTAFile);
		while (!feof(FASTAFile) && c != '>') {
			// Get sequence
			if (c != '\n' && c != '\t') {
				if (maskLowerCase && c >= 'a' && c <= 'z') {
					c = dnaChar[lowercaseDnaCharIndex];
				}
				if (ambiguityCount[charMap[c]] == 1) {
					buffer[numCharInBuffer] = c;
				} else {
					c = charMap[c];
					if (ambiguityCount[c] > 0) {
						buffer[numCharInBuffer] = dnaChar[ambiguityMatch[c][r250() % ambiguityCount[c]]];
					} else {
						buffer[numCharInBuffer] = dnaChar[r250() % ALPHABET_SIZE];
					}
					if (totalNumChar + numChar == lastAmbiguityPos + 1 && c == (char)ambiguity[numAmbiguity].symbol) {
						ambiguity[numAmbiguity].rightOfEndPos++;
					} else {
						numAmbiguity++;
						if (numAmbiguity >= ambiguityAllocated) {
							ambiguity = MMUnitReallocate(ambiguity, sizeof(Ambiguity) * ambiguityAllocated * 2, sizeof(Ambiguity) * ambiguityAllocated);
							ambiguityAllocated *= 2;
						}
						ambiguity[numAmbiguity].startPos = totalNumChar + numChar;
						ambiguity[numAmbiguity].rightOfEndPos = totalNumChar + numChar;
						ambiguity[numAmbiguity].symbol = c;
					}
					lastAmbiguityPos = totalNumChar + numChar;
				}
				numCharInBuffer++;
				if (numCharInBuffer >= PACKED_BUFFER_SIZE) {
					ConvertTextToBytePacked(buffer, packedBuffer, charMap, 4, PACKED_BUFFER_SIZE);
					fwrite(packedBuffer, 1, PACKED_BUFFER_SIZE / 4, packedDNAFile);
					numCharInBuffer = 0;
				}
				numChar++;
			}
			c = (char)getc(FASTAFile);
		}

		seqOffset[numSeq].endPos = totalNumChar + numChar - 1;
		seqOffset[numSeq].firstAmbiguityIndex = numAmbiguity + 1;
		totalNumChar += numChar;
		numSeq++;

	}

	// Finish reading FASTA file
	fclose(FASTAFile);


	// Finalize packed DNA file

	numAmbiguity++;

	if (numCharInBuffer > 0) {
		ConvertTextToBytePacked(buffer, packedBuffer, charMap, 4, numCharInBuffer);
		fwrite(packedBuffer, 1, (numCharInBuffer + 3) / 4, packedDNAFile);
		numCharInBuffer = 0;
	}
	if (totalNumChar % 4 == 0) {
		c = 0;
		fwrite(&c, 1, 1, packedDNAFile);
	}
	c = (char)(totalNumChar % 4);
	fwrite(&c, 1, 1, packedDNAFile);

	fclose(packedDNAFile);

	// Output annotation file
	fprintf(annotationFile, "%u %u %u\n", totalNumChar, numSeq, FASTARandomSeed);
	for (i=0; i<numSeq; i++) {
		fprintf(annotationFile, "%u %s\n", annotation[i].gi, annotation[i].text);
		fprintf(annotationFile, "%u %u ", seqOffset[i].startPos, seqOffset[i].endPos - seqOffset[i].startPos + 1);
		// output number of ambiguity
		if (i > 0) {
			fprintf(annotationFile, "%u\n", seqOffset[i].firstAmbiguityIndex - seqOffset[i-1].firstAmbiguityIndex);
		} else {
			fprintf(annotationFile, "%u\n", seqOffset[i].firstAmbiguityIndex);
		}
	}
	fclose(annotationFile);

	MMUnitFree(annotation, sizeof(Annotation) * annotationAllocated);
	MMUnitFree(seqOffset, sizeof(SeqOffset) * annotationAllocated);

	// Output ambiguity file
	fprintf(ambiguityFile, "%u %u %u\n", totalNumChar, numSeq, numAmbiguity);
	for (i=0; i<numAmbiguity; i++) {
		// The ambiguity for the dummy length is visible
		fprintf(ambiguityFile, "%u %u %c\n", ambiguity[i].startPos, ambiguity[i].rightOfEndPos - ambiguity[i].startPos + 1, dnaChar[ambiguity[i].symbol]);
	}
	fclose(ambiguityFile);

	MMUnitFree(ambiguity, ambiguityAllocated * sizeof(Ambiguity));

	return numSeq;

}

unsigned int HSPPackedToFASTA(const char* FASTAFileName, const char* annotationFileName, const char* packedDNAFileName, const char* ambiguityFileName) {

	HSP *hsp;
	FILE *FASTAFile;

	hsp = HSPLoad(NULL, packedDNAFileName, annotationFileName, ambiguityFileName, 1);

	// Generate FASTA from packed
	FASTAFile = (FILE*)fopen64(FASTAFileName, "w");
	if (FASTAFile == NULL) {
		fprintf(stderr, "Cannot open FASTA file!\n");
		exit(1);
	}

	// to be done...
	fprintf(stderr, "HSPPackedToFASTA(): Function not complete!\n");

	fclose(FASTAFile);

	HSPFree(NULL, hsp, 1);

	return 0;

}

// The hit must be an mismatch hit anchored on the left most of the pattern
/*
unsigned int HSPShortPattern1GapRightExt(const unsigned LONG *packedDNA, const unsigned char *pattern,
						const unsigned int patternLength, const unsigned int textLength,
						HitList* __restrict hitList, const unsigned int numHit, 
						const unsigned int maxSpaceInGap, const unsigned int *maxMismatch) {

	unsigned int anchorLength;
	unsigned int maxSpaceInGap;
	unsigned int maxMismatch[MAX_SP_SPACE_IN_GAP];

	unsigned LONG anchorText, anchorPattern;
	unsigned LONG nonAnchorText, nonAnchorPattern;
	unsigned LONG nonAnchorTextReversed, nonAnchorPatternReversed;
	unsigned int pos, shift, index;
	unsigned maxRightShift;
	unsigned int anchorNumError;
	unsigned int nonAnchorError[MAX_SP_MISMATCH];
	unsigned int leftShiftError[MAX_SP_MISMATCH];
	unsigned int rightShiftError[MAX_SP_MISMATCH];

	__m128i t;
	__m128i m0, m1, m3, m0f, m00ff, m0000ffff;
	__m128i ts, tn, p, mb, mb1, mp, mneg, mb31, md;

	unsigned LONG ALIGN_16 temp[2];
	unsigned int numExtHit;
	const unsigned char *thisPattern;


	const static int bitPos[32] = {  0, 16, 14,  1, 30,  7, 12, 17, 15, 11, 10, 23, 28, 24,  2,  4, 
									31, 29, 22, 27, 26, 25,  8, 19, 13,  6,  9,  3, 21, 18,  5, 20  };

	unsigned int i, h;

	if (patternLength >= MAX_SP_NON_ANCHOR_LENGTH) {
		anchorLength = patternLength - MAX_SP_NON_ANCHOR_LENGTH;
	} else {
		anchorLength = 0;
	}

	// Determine maximum no. of spaces and maximum mismatch corresponding to the no. of spaces
	if (maxPenalty >= gapOpenScore) {
		maxSpaceInGap = (maxPenalty - gapOpenScore) / gapExtendScore;
	} else {
		maxSpaceInGap = 0;
	}

	for (i=0; i<MAX_SP_SPACE_IN_GAP; i++) {
		maxMismatch[i] = 0;
	}
	for (i=0; i<maxSpaceInGap; i++) {
		maxMismatch[i] = (maxPenalty - gapOpenScore - gapOpenScore * i) / (matchScore - mismatchScore);
	}

	numExtHit = 0;

	for (h=0; h<numHit; h++) {

		thisPattern = pattern + patternLength * hitList[h].info;
		if (hitList[h].posText + patternLength <= textLength) {
			maxRightShift = textLength - hitList[h].posText - patternLength;
			if (maxRightShift > MAX_SP_SPACE_IN_GAP) {
				maxRightShift = MAX_SP_SPACE_IN_GAP;
			}
		} else {
			continue;
		} 

		anchorNumError = 0;

		anchorPattern = 0;
		for (i=0; i<anchorLength; i++) {
			anchorPattern <<= 2;
			anchorPattern |= thisPattern[i];
		}

		nonAnchorPattern = 0;
		for (; i<patternLength; i++) {
			nonAnchorPattern <<= 2;
			nonAnchorPattern |= thisPattern[i];
		}

		nonAnchorPatternReversed = 0;	// reverse pattern
		for (i=patternLength; --i;) {	// from patternLength - 1 to 0
			nonAnchorPatternReversed <<= 2;
			nonAnchorPatternReversed |= thisPattern[i];
		}

		shift = hitList[h].posText % CHAR_PER_64;
		pos = hitList[h].posText / CHAR_PER_64;
		if (shift) {
			anchorText = (packedDNA[pos] << shift) + (packedDNA[pos + 1] >> (64 - shift));
		} else {
			anchorText = packedDNA[pos / CHAR_PER_64];
		}

		shift = (hitList[h].posText + anchorLength - MAX_SP_SPACE_IN_GAP) % CHAR_PER_64;
		pos = (hitList[h].posText + anchorLength - MAX_SP_SPACE_IN_GAP) / CHAR_PER_64;
		if (shift) {
			nonAnchorText = (packedDNA[pos] << shift) + (packedDNA[pos + 1] >> (64 - shift));
		} else {
			nonAnchorText = packedDNA[pos / CHAR_PER_64];
		}

		// Reverse nonAnchorText and nonAnchorPattern
		temp[0] = nonAnchorText;
		temp[1] = nonAnchorPattern;
		t = _mm_load_si128((__m128i*)temp);
		m3 = _mm_set1_epi32(0x33333333);
		m0f = _mm_set1_epi32(0x0f0f0f0f);
		m00ff = _mm_set1_epi32(0x00ff00ff);
		m0000ffff = _mm_set1_epi32(0x0000ffff);

		ts = _mm_srli_epi64(t, 2);
		ts = _mm_and_si128(ts, m3);
		tn = _mm_and_si128(t, m3);
		tn = _mm_srli_epi64(tn, 2);
		t = _mm_or_si128(ts, tn);

		ts = _mm_srli_epi64(t, 4);
		ts = _mm_and_si128(ts, m0f);
		tn = _mm_and_si128(t, m0f);
		tn = _mm_srli_epi64(tn, 4);
		t = _mm_or_si128(ts, tn);

		ts = _mm_srli_epi64(t, 8);
		ts = _mm_and_si128(ts, m00ff);
		tn = _mm_and_si128(t, m00ff);
		tn = _mm_srli_epi64(tn, 8);
		t = _mm_or_si128(ts, tn);

		ts = _mm_srli_epi64(t, 16);
		ts = _mm_and_si128(ts, m0000ffff);
		tn = _mm_and_si128(t, m0000ffff);
		tn = _mm_srli_epi64(tn, 16);
		t = _mm_or_si128(ts, tn);

		ts = _mm_srli_epi64(t, 32);
		tn = _mm_srli_epi64(tn, 32);
		t = _mm_or_si128(ts, tn);

		_mm_store_si128((__m128i*)temp, t);
		nonAnchorTextReversed = temp[0];
		nonAnchorPatternReversed = temp[1];

		m1 = _mm_set1_epi32(0xFFFFFFFF);
		m0 = _mm_setzero_si128();

		// Find the mismatches for anchor and non-anchor no-shift
		temp[0] = anchorText;
		temp[1] = nonAnchorTextReversed;
		t = _mm_load_si128((__m128i*)temp);
		temp[0] = anchorPattern;
		temp[1] = nonAnchorPatternReversed;
		p = _mm_load_si128((__m128i*)temp);
		mb = _mm_and_si128(t, p);
		mb1 = _mm_srli_epi64(mb, 1);
		mb = _mm_and_si128(mb, mb1);
		mb = _mm_andnot_si128(mb, m1);
		
		mneg = _mm_sub_epi64(m0, mb);
		mb = _mm_and_si128(mb, mneg);

		mb31 = _mm_srli_epi64(mb, 31);
		mb = _mm_and_si128(mb, mb31);

		md = _mm_set1_epi32(0x077CB531);
		mp = _mm_mul_epu32(mb, md);
		mp = _mm_srli_epi32(mp, 27);
		
		index = _mm_extract_epi16(mp, 0);
		anchorNumError += (index > 0);

		index = _mm_extract_epi16(mp, 4);
		nonAnchorError[0] = bitPos[index];

		// Continue with error 2 and 3

		
		// Find the mismatches for non-anchor 1-shift

		// Find the mismatches for non-anchor 2-shift

		// Find the mismatches for non-anchor 3-shift

	}


	// Find mismatches in non-anchor non-shift comparison

	// Find mismatches in non-anchor left-shift comparison

	// Find mismatches in non-anchor right-shift comparison

	// Check spaces and mismatches


	
	// reverse the text and the pattern

	// v = (v64 >> 31) | v64


	// bitPos[((v & -v) * 0x077CB531UL) >> 27];



}
*/

unsigned int HSPUngappedExtension(const unsigned int *packedDNA, const unsigned int *packedKey, const unsigned int *packedMask, const unsigned int queryPatternLength, 
						 const unsigned int subPattenLength, const unsigned int textLength,
						 HitListWithPosQuery* __restrict hitList, const unsigned int numberOfHit,
						 const HSPUngappedExtLookupTable *ungappedLookupTable,
						 const int cutoffScore, const int dropoffScore) {

	static const unsigned int oddBitMask = 0x55555555;
	static const unsigned int evenBitMask = 0xAAAAAAAA;

	unsigned int hitProcessed, numberOfUngappedHit;
	unsigned int textShift, queryShift;
	unsigned int posText, posQuery;
	unsigned int queryWordPos, textWordPos;
	unsigned int matchVector32, matchVector16;
	int numberOfValidChar;
	unsigned int currentPos;
	int currentScore, forwardMaxScore, backwardMaxScore;
	unsigned int forwardMaxScorePos, backwardMaxScorePos;
	
	hitProcessed = 0;
	numberOfUngappedHit = 0;

	while (hitProcessed < numberOfHit) {

		// Move position to mid-point of hit
		posText = hitList[hitProcessed].posText + subPattenLength / 2;
		posQuery = hitList[hitProcessed].posQuery + subPattenLength / 2;

		textShift = (posText % CHAR_PER_WORD) * BIT_PER_CHAR;
		queryShift = (posQuery % CHAR_PER_WORD) * BIT_PER_CHAR;

		// Forward extend

		textWordPos = posText / CHAR_PER_WORD;
		queryWordPos = posQuery / CHAR_PER_WORD;
		numberOfValidChar = minX(queryPatternLength - posQuery, textLength - posText);

		forwardMaxScore = 0;
		forwardMaxScorePos = posQuery;
		currentPos = posQuery;
		currentScore = 0;

		while (forwardMaxScore <= currentScore + dropoffScore && numberOfValidChar > 0) {

			matchVector32 = ((packedKey[queryWordPos] << queryShift) | ((packedKey[queryWordPos+1] >> (BITS_IN_WORD - queryShift)) * (queryShift > 0)))
							^ ((packedDNA[textWordPos] << textShift) | ((packedDNA[textWordPos+1] >> (BITS_IN_WORD - textShift)) * (textShift > 0)));

			matchVector32 |= (packedMask[queryWordPos] << queryShift) | ((packedMask[queryWordPos+1] >> (BITS_IN_WORD - queryShift)) * (queryShift > 0));

			if (numberOfValidChar < CHAR_PER_WORD) {
				// mask the invalid character as mismatch
				matchVector32 |= ALL_ONE_MASK >> (numberOfValidChar * BIT_PER_CHAR);
			}

			// Process vector to even bits for forward extend
			matchVector32 = (matchVector32 | (matchVector32 << 1)) & evenBitMask;
			matchVector16 = (ungappedLookupTable[matchVector32 >> 16].matchMismatchBitVector << 8) |
							 ungappedLookupTable[matchVector32 & 0xFFFF].matchMismatchBitVector;

			if (currentScore + ungappedLookupTable[matchVector16].maxScore > forwardMaxScore) {
				forwardMaxScore = currentScore + ungappedLookupTable[matchVector16].maxScore;
				forwardMaxScorePos = currentPos + ungappedLookupTable[matchVector16].maxScorePos;
			}
			currentScore += ungappedLookupTable[matchVector16].finalScore;
			currentPos += CHAR_PER_WORD;

			textWordPos++;
			queryWordPos++;
			numberOfValidChar -= CHAR_PER_WORD;

		}


		// Backward extend by word
		textWordPos = posText / CHAR_PER_WORD;
		queryWordPos = posQuery / CHAR_PER_WORD;
		numberOfValidChar = minX(posQuery, posText);

		backwardMaxScore = 0;
		backwardMaxScorePos = posQuery;
		currentPos = posQuery;
		currentScore = 0;

		while (backwardMaxScore <= currentScore + dropoffScore && numberOfValidChar > 0) {

			if (numberOfValidChar >= CHAR_PER_WORD) {
				matchVector32 = ((packedKey[queryWordPos-1] << queryShift) | ((packedKey[queryWordPos] >> (BITS_IN_WORD - queryShift)) * (queryShift > 0)))
								^ ((packedDNA[textWordPos-1] << textShift) | ((packedDNA[textWordPos] >> (BITS_IN_WORD - textShift)) * (textShift > 0)));
				matchVector32 |= (packedMask[queryWordPos-1] << queryShift) | ((packedMask[queryWordPos] >> (BITS_IN_WORD - queryShift)) * (queryShift > 0));
			} else {
				if (textWordPos > 0) {
					matchVector32 = ((packedKey[queryWordPos] >> (BITS_IN_WORD - queryShift)) * (queryShift > 0))
									^ ((packedDNA[textWordPos-1] << textShift) | ((packedDNA[textWordPos] >> (BITS_IN_WORD - textShift)) * (textShift > 0)));
					matchVector32 |= (packedMask[queryWordPos] >> (BITS_IN_WORD - queryShift)) * (queryShift > 0);
				} else if (queryWordPos > 0) {
					matchVector32 = ((packedKey[queryWordPos-1] << queryShift) | ((packedKey[queryWordPos] >> (BITS_IN_WORD - queryShift)) * (queryShift > 0)))
									^ ((packedDNA[textWordPos] >> (BITS_IN_WORD - textShift)) * (textShift > 0));
					matchVector32 |= (packedMask[queryWordPos-1] << queryShift) | ((packedMask[queryWordPos] >> (BITS_IN_WORD - queryShift)) * (queryShift > 0));
				} else {
					matchVector32 = ((packedKey[queryWordPos] >> (BITS_IN_WORD - queryShift)) * (queryShift > 0))
									^ ((packedDNA[textWordPos] >> (BITS_IN_WORD - textShift)) * (textShift > 0));
					matchVector32 |= (packedMask[queryWordPos] >> (BITS_IN_WORD - queryShift)) * (queryShift > 0);
				}
				matchVector32 |= ALL_ONE_MASK << (numberOfValidChar * BIT_PER_CHAR);
			}

			// Process vector to odd bits for backward extend
			matchVector32 = (matchVector32 | (matchVector32 >> 1)) & oddBitMask;
			matchVector16 = ungappedLookupTable[matchVector32 >> 16].matchMismatchBitVector |
							(ungappedLookupTable[matchVector32 & 0xFFFF].matchMismatchBitVector << 8);	// matchVecter16 is reverse of packed odd bits

			if (currentScore + ungappedLookupTable[matchVector16].maxScore > backwardMaxScore) {
				backwardMaxScore = currentScore + ungappedLookupTable[matchVector16].maxScore;
				backwardMaxScorePos = currentPos - ungappedLookupTable[matchVector16].maxScorePos;
			}
			currentScore += ungappedLookupTable[matchVector16].finalScore;
			currentPos -= CHAR_PER_WORD;

			textWordPos--;
			queryWordPos--;
			numberOfValidChar -= CHAR_PER_WORD;

		}

		// Check cut off score
		if (backwardMaxScore + forwardMaxScore >= cutoffScore) {
			hitList[numberOfUngappedHit].posQuery = (forwardMaxScorePos + backwardMaxScorePos) / 2;
			hitList[numberOfUngappedHit].posText = posText - posQuery + hitList[numberOfUngappedHit].posQuery;
			numberOfUngappedHit++;
		}

		hitProcessed++;

	}

	return numberOfUngappedHit;

}

int HSPGappedExtension(const unsigned int *packedDNA, const unsigned int textLength, const unsigned char *convertedKey, const int queryPatternLength, 
					   const HitListWithPosQuery* hitList, const int numberOfHit, 
					   GappedHitList* __restrict gappedHitList,
					   MMPool *mmPool,
					   const int matchScore, const int mismatchScore,
					   const int gapOpenScore, const int gapExtendScore,
					   const int cutoffScore, const int dropoffScore) {

	int* __restrict B;
	int M, BLast;
	int* __restrict I;	// insert and delete wrt query
	int D;				// insert and delete wrt query
	int scoringMatrix[16][16];

	int hitProcessed, numberOfGappedHit;

	int i;
	unsigned int c = 0;
	char textChar;
	int queryOffset;
	unsigned int textOffset;
	unsigned int textDecodePos;
	int charToProcess;

	int lastStartPos, startPos, nextStartPos;
	int endPos, nextEndPos;
	unsigned int textPos;

	int maxLengthOfGap;
	int currentPos;
	int forwardMaxScore, backwardMaxScore;
	unsigned int forwardMaxScorePosQuery, forwardMaxScorePosText;
	unsigned int backwardMaxScorePosQuery, backwardMaxScorePosText;

	int dpCellAllocated;

	dpCellAllocated = minX(queryPatternLength + 1, MAX_ALIGNMENT_LENGTH);
	
	// allocate working memory
	B = MMPoolDispatch(mmPool, dpCellAllocated * sizeof(int));
	I = MMPoolDispatch(mmPool, dpCellAllocated * sizeof(int));

	HSPFillScoringMatrix(scoringMatrix, matchScore, mismatchScore, 0);

	hitProcessed = 0;
	numberOfGappedHit = 0;

	maxLengthOfGap = (dropoffScore + gapOpenScore) / (-gapExtendScore);

	while (hitProcessed < numberOfHit) {

		queryOffset = hitList[hitProcessed].posQuery;
		textOffset = hitList[hitProcessed].posText;

		// Forward extend

		textDecodePos = textOffset;
		if (textDecodePos % CHAR_PER_WORD > 0) {
			// decode text to the next word boundary
			charToProcess = CHAR_PER_WORD - (textDecodePos % CHAR_PER_WORD);
			c = packedDNA[textDecodePos / CHAR_PER_WORD] << (BITS_IN_WORD - charToProcess * BIT_PER_CHAR);
			textDecodePos += charToProcess;
		}

		// Fill initial scores
		B[0] = 0;
		I[0] = gapExtendScore + gapOpenScore;
		for (i=1; i<=maxLengthOfGap; i++) {
			B[i] = i * gapExtendScore + gapOpenScore;
			I[i] = B[i] + gapExtendScore + gapOpenScore;
		}
		I[i] = DP_NEG_INFINITY;

		lastStartPos = startPos = 0;
		endPos = minX(maxLengthOfGap + 1, queryPatternLength - queryOffset);
		textPos = textOffset;
		 
		forwardMaxScore = 0;
		forwardMaxScorePosQuery = queryOffset;
		forwardMaxScorePosText = textOffset;

		while (startPos <= endPos && textPos < textLength) {

			if (endPos >= dpCellAllocated) {
				fprintf(stderr, "HSPGappedExtension(): Not enough DP cells allocated!\n");
				exit(1);
			}

			if (textPos >= textDecodePos) {
				c = packedDNA[textDecodePos / CHAR_PER_WORD];
				textDecodePos += CHAR_PER_WORD;
			}
			textChar = (char)(c >> (BITS_IN_WORD - BIT_PER_CHAR));
			c <<= BIT_PER_CHAR;

			nextEndPos = 0;
			nextStartPos = INT_MAX;

			if (startPos > lastStartPos) {
				BLast = B[startPos - 1];
				D = DP_NEG_INFINITY;
				currentPos = startPos;
			} else {
				// all scores in currentPos - 1 is DP_NEG_INFINITY
				BLast = B[startPos];
				B[startPos] = I[startPos];
				I[startPos] = I[startPos] + gapExtendScore;
				D = B[startPos] + gapExtendScore + gapOpenScore;
				if (B[startPos] + dropoffScore >= forwardMaxScore) {
					nextStartPos = startPos;
					nextEndPos = startPos;
				}
				currentPos = startPos + 1;
			}

			while (currentPos <= endPos) {

				M = BLast + scoringMatrix[textChar][convertedKey[currentPos + queryOffset - 1]];
				BLast = B[currentPos];

				if (M >= I[currentPos] && M >= D) {
					// matchScore is maximum
					B[currentPos] = M;
					I[currentPos] = maxX(M + gapOpenScore, I[currentPos]) + gapExtendScore;
					D = maxX(M + gapOpenScore, D) + gapExtendScore;
				} else {
					B[currentPos] = maxX(I[currentPos], D);
					// insert cannot follow delete and delete cannot follow insert
					I[currentPos] = I[currentPos] + gapExtendScore;
					D = D + gapExtendScore;
				}
				if (B[currentPos] + dropoffScore >= forwardMaxScore) {
					if (nextStartPos > currentPos) {
						nextStartPos = currentPos;
					}
					nextEndPos = currentPos;
				}
				if (B[currentPos] > forwardMaxScore) {
					forwardMaxScore = B[currentPos];
					forwardMaxScorePosQuery = currentPos - 1 + queryOffset;
					forwardMaxScorePosText = textPos;
				}

				currentPos++;

			}

			nextEndPos++;
			if (nextEndPos == currentPos) {
				while (nextEndPos <= queryPatternLength - queryOffset && D + dropoffScore >= forwardMaxScore)  {
					if (nextEndPos >= dpCellAllocated) {
						fprintf(stderr, "HSPGappedExtension(): Not enough DP cells allocated!\n");
						exit(1);
					}
					B[nextEndPos] = D;
					D = D + gapExtendScore;
					I[nextEndPos] = DP_NEG_INFINITY;
					nextEndPos++;
				}
				I[nextEndPos] = DP_NEG_INFINITY;
			}

			if (nextEndPos > queryPatternLength - queryOffset) {
				nextEndPos = queryPatternLength - queryOffset;
			}

			lastStartPos = startPos;
			startPos = nextStartPos;
			endPos = nextEndPos;

			textPos++;

		}

		// Backward extend

		textDecodePos = textOffset;
		if (textDecodePos % CHAR_PER_WORD > 0) {
			charToProcess = textDecodePos % CHAR_PER_WORD;
			c = packedDNA[textDecodePos / CHAR_PER_WORD] >> (BITS_IN_WORD - charToProcess * BIT_PER_CHAR);
			textDecodePos -= charToProcess;
		}

		// Fill initial scores
		B[0] = 0;
		I[0] = gapExtendScore + gapOpenScore;
		for (i=1; i<=maxLengthOfGap; i++) {
			B[i] = i * gapExtendScore + gapOpenScore;
			I[i] = B[i] + gapExtendScore + gapOpenScore;
		}
		I[i] = DP_NEG_INFINITY;

		lastStartPos = startPos = 0;
		endPos = minX(maxLengthOfGap + 1, queryOffset);
		textPos = textOffset;
		 
		backwardMaxScore = 0;
		backwardMaxScorePosQuery = queryOffset;
		backwardMaxScorePosText = textOffset;

		while (startPos <= endPos && textPos > 0) {

			if (endPos >= dpCellAllocated) {
				fprintf(stderr, "HSPGappedExtension(): Not enough DP cells allocated!\n");
				exit(1);
			}

			if (textPos <= textDecodePos) {
				c = packedDNA[textDecodePos / CHAR_PER_WORD - 1];
				textDecodePos -= CHAR_PER_WORD;
			}
			textChar = (char)(c & 0x3);
			c >>= BIT_PER_CHAR;

			nextEndPos = 0;
			nextStartPos = INT_MAX;

			if (startPos > lastStartPos) {
				BLast = B[startPos - 1];
				D = DP_NEG_INFINITY;
				currentPos = startPos;
			} else {
				// all scores in currentPos - 1 is DP_NEG_INFINITY
				BLast = B[startPos];
				B[startPos] = I[startPos];
				I[startPos] = I[startPos] + gapExtendScore;
				D = B[startPos] + gapExtendScore + gapOpenScore;
				if (B[startPos] + dropoffScore >= backwardMaxScore) {
					nextStartPos = startPos;
					nextEndPos = startPos;
				}
				currentPos = startPos + 1;
			}

			while (currentPos <= endPos) {

				M = BLast + scoringMatrix[textChar][convertedKey[queryOffset - currentPos]];
				BLast = B[currentPos];

				if (M >= I[currentPos] && M >= D) {
					// matchScore is maximum
					B[currentPos] = M;
					I[currentPos] = maxX(M + gapOpenScore, I[currentPos]) + gapExtendScore;
					D = maxX(M + gapOpenScore, D) + gapExtendScore;
				} else {
					B[currentPos] = maxX(I[currentPos], D);
					// insert cannot follow delete and delete cannot follow insert
					I[currentPos] = I[currentPos] + gapExtendScore;
					D = D + gapExtendScore;
				}
				if (B[currentPos] + dropoffScore >= backwardMaxScore) {
					if (nextStartPos > currentPos) {
						nextStartPos = currentPos;
					}
					nextEndPos = currentPos;
				}
				if (B[currentPos] > backwardMaxScore) {
					backwardMaxScore = B[currentPos];
					backwardMaxScorePosQuery =  queryOffset - currentPos;
					backwardMaxScorePosText = textPos - 1;
				}

				currentPos++;

			}

			nextEndPos++;
			if (nextEndPos == currentPos) {
				while (nextEndPos <= queryOffset && D + dropoffScore >= backwardMaxScore)  {
					if (nextEndPos >= dpCellAllocated) {
						fprintf(stderr, "HSPGappedExtension(): Not enough DP cells allocated!\n");
						exit(1);
					}
					B[nextEndPos] = D;
					D = D + gapExtendScore;
					I[nextEndPos] = DP_NEG_INFINITY;
					nextEndPos++;
				}
				I[nextEndPos] = DP_NEG_INFINITY;
			}

			if (nextEndPos > queryOffset) {
				nextEndPos = queryOffset;
			}

			lastStartPos = startPos;
			startPos = nextStartPos;
			endPos = nextEndPos;

			textPos--;

		}

		// check cutoff score
		if (forwardMaxScore + backwardMaxScore >= cutoffScore) {
			gappedHitList[numberOfGappedHit].posText = backwardMaxScorePosText;
			gappedHitList[numberOfGappedHit].posQuery = backwardMaxScorePosQuery;
			gappedHitList[numberOfGappedHit].lengthQuery = forwardMaxScorePosQuery + 1 - backwardMaxScorePosQuery;
			gappedHitList[numberOfGappedHit].score = forwardMaxScore + backwardMaxScore;
			gappedHitList[numberOfGappedHit].lengthText = forwardMaxScorePosText + 1 - backwardMaxScorePosText;
			gappedHitList[numberOfGappedHit].ungappedPosText = textOffset;
			gappedHitList[numberOfGappedHit].ungappedPosQuery = queryOffset;
			numberOfGappedHit++;
		}

		hitProcessed++;

	}

	// free working memory
	MMPoolReturn(mmPool, B, dpCellAllocated * sizeof(int));
	MMPoolReturn(mmPool, I, dpCellAllocated * sizeof(int));


	return numberOfGappedHit;

}

// mmPool can be NULL but alignmentPool must be allocated for storing alignments as they are not freed explicitly and must be freed through MMPool
// gappedHitList should be sorted in increasin posText order
int HSPGappedExtensionWithTraceback(const unsigned int *packedDNA, const unsigned char *convertedKey, const int queryPatternLength, 
					   GappedHitList* __restrict gappedHitList, const int numberOfHit, 
					   MMPool *mmPool, MMPool *alignmentPool,
					   const SeqOffset *seqOffset, const Ambiguity *ambiguity,
					   const int matchScore, const int mismatchScore,
					   const int gapOpenScore, const int gapExtendScore,
					   const double maxEvalue, const int dropoffScore) {

	char* __restrict textBuffer;
	int M, BLast;
	int* __restrict B;
	int* __restrict I;	// insert and delete wrt query
	int D;				// insert and delete wrt query
	char* __restrict IType;	// Gapped opening or gap extension
	char DType;				// Gapped opening or gap extension
	int scoringMatrix[16][16];

	char* __restrict traceback;
	int* __restrict tracebackIndex;
	char* __restrict tempForwardAlignment;
	char* __restrict tempBackwardAlignment;
	char* __restrict tempForwardAuxiliaryText;
	char* __restrict tempBackwardAuxiliaryText;
	unsigned int* __restrict alignment;
	unsigned int* __restrict auxiliaryText;

	int hitProcessed, numberOfGappedHit;

	int i, j;
	unsigned int c;
	int queryOffset;
	unsigned int textOffset;
	unsigned int sequenceStart, sequenceEnd;
	int charToProcess;
	unsigned int textDecodePos;
	double evalue;

	int lastStartPos, startPos, nextStartPos;
	int endPos, nextEndPos;
	unsigned int textPos;
	int ambiguityIndex;

	int maxLengthOfGap;
	int currentPos;
	int currentTracebackIndex;
	int forwardMaxScore, backwardMaxScore;
	unsigned int forwardMaxScorePosQuery, forwardMaxScorePosText;
	unsigned int backwardMaxScorePosQuery, backwardMaxScorePosText;
	int forwardNumOfAlignment, forwardNumOfAuxiliaryText;
	int backwardNumOfAlignment, backwardNumOfAuxiliaryText;
	int numOfAlignmentWord, numOfAuxiliaryTextWord;
	int tracebackAllocationUnit, tracebackAllocated;
	char dpType;
	int dpScore;

	int dpCellAllocated;

	dpCellAllocated = minX(queryPatternLength + 1, MAX_ALIGNMENT_LENGTH);

	if (alignmentPool == NULL) {
		fprintf(stderr, "HSPGappedExtensionWithTraceback(): alignmentPool is not allocated!\n");
		exit(1);
	}

	maxLengthOfGap = (dropoffScore + gapOpenScore) / (-gapExtendScore);

	// allocate working memory
	textBuffer = MMPoolDispatch(mmPool, dpCellAllocated * 2);
	B = MMPoolDispatch(mmPool, dpCellAllocated * sizeof(int));
	IType = MMPoolDispatch(mmPool, dpCellAllocated);
	I = MMPoolDispatch(mmPool, dpCellAllocated * sizeof(int));

	tracebackIndex = MMPoolDispatch(mmPool, dpCellAllocated * 2 * 2 * sizeof(int));

	tracebackAllocated = tracebackAllocationUnit = dpCellAllocated * maxLengthOfGap * maxLengthOfGap;
	traceback = MMUnitAllocate(tracebackAllocated);

	tempForwardAlignment = MMPoolDispatch(mmPool, dpCellAllocated * 2);
	tempBackwardAlignment = MMPoolDispatch(mmPool, dpCellAllocated * 2);
	tempForwardAuxiliaryText = MMPoolDispatch(mmPool, dpCellAllocated);
	tempBackwardAuxiliaryText = MMPoolDispatch(mmPool, dpCellAllocated);

	HSPFillScoringMatrix(scoringMatrix, matchScore, mismatchScore, 0);

	hitProcessed = 0;
	numberOfGappedHit = 0;
	ambiguityIndex = 0;

	while (hitProcessed < numberOfHit) {

		sequenceStart = seqOffset[gappedHitList[hitProcessed].dbSeqIndex].startPos;
		sequenceEnd = seqOffset[gappedHitList[hitProcessed].dbSeqIndex].endPos;
		if (seqOffset[gappedHitList[hitProcessed].dbSeqIndex].firstAmbiguityIndex > ambiguityIndex) {
			ambiguityIndex = seqOffset[gappedHitList[hitProcessed].dbSeqIndex].firstAmbiguityIndex;
		}

		textOffset = gappedHitList[hitProcessed].ungappedPosText;
		queryOffset = gappedHitList[hitProcessed].ungappedPosQuery;


		// Forward extend

		textDecodePos = textOffset;
		if (textDecodePos % CHAR_PER_WORD > 0) {
			// decode text to the next word boundary
			charToProcess = CHAR_PER_WORD - (textDecodePos % CHAR_PER_WORD);
			c = packedDNA[textDecodePos / CHAR_PER_WORD] << (BITS_IN_WORD - charToProcess * BIT_PER_CHAR);
			for (i=0; i<charToProcess; i++) {
				textBuffer[textDecodePos + i - textOffset + 1] = (unsigned char)(c >> (BITS_IN_WORD - BIT_PER_CHAR));
				c <<= BIT_PER_CHAR;
			}

			// Apply ambiguity
			while (ambiguity[ambiguityIndex].rightOfEndPos <= textOffset) {
				ambiguityIndex++;
			}
			if (ambiguity[ambiguityIndex].startPos < textOffset + charToProcess) {
				for (i=0; i<charToProcess; i++) {
					while (ambiguity[ambiguityIndex].rightOfEndPos <= textOffset + i) {
						ambiguity++;
					}
					if (ambiguity[ambiguityIndex].startPos <= textOffset + i) {
						textBuffer[textDecodePos + i - textOffset + 1] = (char)ambiguity[ambiguityIndex].symbol;
					}
				}
			}
				
			textDecodePos += charToProcess;
		}

		// Fill initial scores
		B[0] = 0;
		traceback[0] = DP_MATCH_MISMATCH;
		I[0] = gapExtendScore + gapOpenScore;
		IType[0] = DP_INSERT_OPEN;

		B[1] = gapExtendScore + gapOpenScore;
		traceback[1] = DP_DELETE | DP_DELETE_OPEN;
		I[1] = B[1] + gapExtendScore + gapOpenScore;
		IType[1] = DP_INSERT_OPEN;

		for (i=2; i<=maxLengthOfGap; i++) {
			B[i] = B[i-1] + gapExtendScore;
			traceback[i] = DP_DELETE;
			I[i] = B[i] + gapExtendScore + gapOpenScore;
			IType[i] = DP_INSERT_OPEN;
		}
		I[i] = DP_NEG_INFINITY;

		currentTracebackIndex = maxLengthOfGap + 1;

		lastStartPos = startPos = 0;
		endPos = minX((maxLengthOfGap + 1), queryPatternLength - queryOffset);
		textPos = textOffset;
		 
		forwardMaxScore = 0;
		forwardMaxScorePosQuery = queryOffset;
		forwardMaxScorePosText = textOffset;

		tracebackIndex[0] = 0;	// The first trackback of the row
		tracebackIndex[1] = 0;	// The first query position of the row

		while (startPos <= endPos && textPos <= sequenceEnd) {

			if (endPos >= dpCellAllocated) {
				fprintf(stderr, "HSPGappedExtensionWithTraceback(): Not enough DP cells allocated!\n");
				exit(1);
			}

			if (textPos < textDecodePos) {
			} else {
				// decode a word of text
				c = packedDNA[textDecodePos / CHAR_PER_WORD];
				for (j=0; j<CHAR_PER_WORD; j++) {
					textBuffer[textDecodePos + j - textOffset + 1] = (unsigned char)(c >> (BITS_IN_WORD - BIT_PER_CHAR));
					c <<= BIT_PER_CHAR;
				}

				// Apply ambiguity
				while (ambiguity[ambiguityIndex].rightOfEndPos <= textPos) {
					ambiguityIndex++;
				}
				if (ambiguity[ambiguityIndex].startPos < textPos + CHAR_PER_WORD) {
					for (j=0; j<CHAR_PER_WORD; j++) {
						while (ambiguity[ambiguityIndex].rightOfEndPos <= textPos + j) {
							ambiguity++;
						}
						if (ambiguity[ambiguityIndex].startPos <= textPos + j) {
							textBuffer[textDecodePos + j - textOffset + 1] = (char)ambiguity[ambiguityIndex].symbol;
						}
					}
				}

				textDecodePos += CHAR_PER_WORD;
			}

			// traceback
			tracebackIndex[(textPos - textOffset + 1) * 2] = currentTracebackIndex;	// The first trackback of the row
			tracebackIndex[(textPos - textOffset + 1) * 2 + 1] = startPos;			// The first query position of the row

			nextEndPos = 0;
			nextStartPos = INT_MAX;

			if (currentTracebackIndex + queryPatternLength >= tracebackAllocated) {
				traceback = MMUnitReallocate(traceback, tracebackAllocated + tracebackAllocationUnit, tracebackAllocated);
				tracebackAllocated += tracebackAllocationUnit;
			}

			if (startPos > lastStartPos) {
				BLast = B[startPos - 1];
				D = DP_NEG_INFINITY;
				DType = 0;
				currentPos = startPos;
			} else {
				// all scores in currentPos - 1 is DP_NEG_INFINITY
				BLast = B[startPos];
				B[startPos] = I[startPos];
				traceback[currentTracebackIndex] = DP_INSERT | IType[startPos];
				I[startPos] = I[startPos] + gapExtendScore;
				IType[startPos] = DP_INSERT_EXTEND;
				D = B[startPos] + gapExtendScore + gapOpenScore;
				DType = DP_DELETE_OPEN;
				if (B[startPos] + dropoffScore >= forwardMaxScore) {
					nextStartPos = startPos;
					nextEndPos = startPos;
				}
				currentPos = startPos + 1;
				currentTracebackIndex++;
			}

			while (currentPos <= endPos) {

				M = BLast + scoringMatrix[textBuffer[textPos - textOffset + 1]][convertedKey[currentPos + queryOffset - 1]];
				BLast = B[currentPos];

				if (M >= I[currentPos] && M >= D) {
					// matchScore is maximum
					B[currentPos] = M;
					traceback[currentTracebackIndex] = DP_MATCH_MISMATCH | IType[currentPos] | DType;
					if (M + gapOpenScore >= I[currentPos]) {
						I[currentPos] = M + gapOpenScore + gapExtendScore;
						IType[currentPos] = DP_INSERT_OPEN;
					} else {
						I[currentPos] = I[currentPos] + gapExtendScore;
						IType[currentPos] = DP_INSERT_EXTEND;
					}
					if (M + gapOpenScore >= D) {
						D = M + gapOpenScore + gapExtendScore;
						DType = DP_DELETE_OPEN;
					} else {
						D = D + gapExtendScore;
						DType = DP_DELETE_EXTEND;
					}
					if (B[currentPos] > forwardMaxScore) {
						forwardMaxScore = B[currentPos];
						forwardMaxScorePosQuery = currentPos + queryOffset;
						forwardMaxScorePosText = textPos + 1;
					}
				} else {
					if (I[currentPos] >= D) {
						B[currentPos] = I[currentPos];
						traceback[currentTracebackIndex] = DP_INSERT | IType[currentPos] | DType;
					} else {
						B[currentPos] = D;
						traceback[currentTracebackIndex] = DP_DELETE | IType[currentPos] | DType;
					}
					// insert cannot follow delete and delete cannot follow insert
					I[currentPos] = I[currentPos] + gapExtendScore;
					IType[currentPos] = DP_INSERT_EXTEND;
					D = D + gapExtendScore;
					DType = DP_DELETE_EXTEND;
				}
				if (B[currentPos] + dropoffScore >= forwardMaxScore) {
					if (nextStartPos > currentPos) {
						nextStartPos = currentPos;
					}
					nextEndPos = currentPos;
				}

				currentPos++;
				currentTracebackIndex++;

			}

			nextEndPos++;
			if (nextEndPos == currentPos) {
				while (nextEndPos <= queryPatternLength - queryOffset && D + dropoffScore >= forwardMaxScore)  {
					if (nextEndPos >= dpCellAllocated) {
						fprintf(stderr, "HSPGappedExtensionWithTraceback(): Not enough DP cells allocated!\n");
						exit(1);
					}
					B[nextEndPos] = D;
					traceback[currentTracebackIndex] = DP_DELETE | DType;
					D = D + gapExtendScore;
					DType = DP_DELETE_EXTEND;
					I[nextEndPos] = DP_NEG_INFINITY;
					IType[nextEndPos] = 0;
					nextEndPos++;
					currentTracebackIndex++;
				}
				I[nextEndPos] = DP_NEG_INFINITY;
				IType[nextEndPos] = 0;
			}

			if (nextEndPos > queryPatternLength - queryOffset) {
				nextEndPos = queryPatternLength - queryOffset;
			}

			lastStartPos = startPos;
			startPos = nextStartPos;
			endPos = nextEndPos;

			textPos++;

		}

		// traceback
		forwardNumOfAlignment = 0;
		forwardNumOfAuxiliaryText = 0;
		textPos = forwardMaxScorePosText - 1;
		currentPos = forwardMaxScorePosQuery - queryOffset;
		dpScore = 0;	// for verifying traceback

		currentTracebackIndex = tracebackIndex[(textPos - textOffset + 1) * 2] 
								+ currentPos - tracebackIndex[(textPos - textOffset + 1) * 2 + 1];

		while (currentTracebackIndex > 0) {

			dpType = traceback[currentTracebackIndex] & DP_MASK;

			if (dpType == DP_MATCH_MISMATCH) {
				dpScore += scoringMatrix[textBuffer[textPos - textOffset + 1]][convertedKey[currentPos + queryOffset - 1]];
				if (textBuffer[textPos - textOffset + 1] == convertedKey[currentPos + queryOffset - 1] &&
					textBuffer[textPos - textOffset + 1] < 4) {	// match and not ambiguity
					tempForwardAlignment[forwardNumOfAlignment] = ALIGN_MATCH;
				} else {
					tempForwardAlignment[forwardNumOfAlignment] = ALIGN_MISMATCH_AMBIGUITY;
					tempForwardAuxiliaryText[forwardNumOfAuxiliaryText] = textBuffer[textPos - textOffset + 1];
					forwardNumOfAuxiliaryText++;
				}
				textPos--;
				currentPos--;
				currentTracebackIndex = tracebackIndex[(textPos - textOffset + 1) * 2] 
										+ currentPos - tracebackIndex[(textPos - textOffset + 1) * 2 + 1];
			} else if (dpType ==  DP_INSERT) {
				while (!(traceback[currentTracebackIndex] & DP_INSERT_OPEN)) {
					dpScore += gapExtendScore;
					tempForwardAlignment[forwardNumOfAlignment] = ALIGN_INSERT;
					forwardNumOfAlignment++;
					tempForwardAuxiliaryText[forwardNumOfAuxiliaryText] = textBuffer[textPos - textOffset + 1];
					forwardNumOfAuxiliaryText++;
					textPos--;
					currentTracebackIndex = tracebackIndex[(textPos - textOffset + 1) * 2] 
											+ currentPos - tracebackIndex[(textPos - textOffset + 1) * 2 + 1];
				}
				dpScore += gapOpenScore + gapExtendScore;
				tempForwardAlignment[forwardNumOfAlignment] = ALIGN_INSERT;
				tempForwardAuxiliaryText[forwardNumOfAuxiliaryText] = textBuffer[textPos - textOffset + 1];
				forwardNumOfAuxiliaryText++;
				textPos--;
				currentTracebackIndex = tracebackIndex[(textPos - textOffset + 1) * 2] 
										+ currentPos - tracebackIndex[(textPos - textOffset + 1) * 2 + 1];
			} else {
				while (!(traceback[currentTracebackIndex] & DP_DELETE_OPEN)) {
					dpScore += gapExtendScore;
					tempForwardAlignment[forwardNumOfAlignment] = ALIGN_DELETE;
					forwardNumOfAlignment++;
					currentPos--;
					currentTracebackIndex--;
				}
				dpScore += gapOpenScore + gapExtendScore;
				tempForwardAlignment[forwardNumOfAlignment] = ALIGN_DELETE;
				currentPos--;
				currentTracebackIndex--;
			}

			forwardNumOfAlignment++;

		}

		if (dpScore != forwardMaxScore) {
			fprintf(stderr, "Forward gapped extension traceback error!\n");
			exit(1);
		}

		// Backward extend

		textDecodePos = textOffset;
		if (textDecodePos % CHAR_PER_WORD > 0) {
			// decode text to the next word boundary
			charToProcess = textDecodePos % CHAR_PER_WORD;
			c = packedDNA[textDecodePos / CHAR_PER_WORD] >> (BITS_IN_WORD - charToProcess * BIT_PER_CHAR);
			for (i=0; i<charToProcess; i++) {
				textBuffer[textOffset - textDecodePos + i + 1] = (unsigned char)(c & 0x3);
				c >>= BIT_PER_CHAR;
			}

			// Apply ambiguity
			while (ambiguity[ambiguityIndex].startPos >= textOffset) {
				ambiguityIndex--;
			}
			if (ambiguity[ambiguityIndex].rightOfEndPos > (textOffset - charToProcess)) {
				for (i=0; i<charToProcess; i++) {
					while (ambiguity[ambiguityIndex].startPos >= (textOffset - i)) {
						ambiguityIndex--;
					}
					if (ambiguity[ambiguityIndex].rightOfEndPos > (textOffset - i - 1)) {
						textBuffer[textOffset - textDecodePos + i + 1] = (char)ambiguity[ambiguityIndex].symbol;
					}
				}
			}

			textDecodePos -= charToProcess;
		}

		// Fill initial scores
		B[0] = 0;
		traceback[0] = DP_MATCH_MISMATCH;
		I[0] = gapExtendScore + gapOpenScore;
		IType[0] = DP_INSERT_OPEN;

		B[1] = gapExtendScore + gapOpenScore;
		traceback[1] = DP_DELETE | DP_DELETE_OPEN;
		I[1] = B[1] + gapExtendScore + gapOpenScore;
		IType[1] = DP_INSERT_OPEN;

		for (i=2; i<=maxLengthOfGap; i++) {
			B[i] = B[i-1] + gapExtendScore;
			traceback[i] = DP_DELETE;
			I[i] = B[i] + gapExtendScore + gapOpenScore;
			IType[i] = DP_INSERT_OPEN;
		}
		I[i] = DP_NEG_INFINITY;

		currentTracebackIndex = maxLengthOfGap + 1;

		lastStartPos = startPos = 0;
		endPos = minX((maxLengthOfGap + 1), queryOffset);
		textPos = textOffset;
		 
		backwardMaxScore = 0;
		backwardMaxScorePosQuery = queryOffset;
		backwardMaxScorePosText = textOffset;

		tracebackIndex[0] = 0;
		tracebackIndex[1] = 0;

		while (startPos <= endPos && textPos + 1 > sequenceStart) {

			if (endPos >= dpCellAllocated) {
				fprintf(stderr, "HSPGappedExtensionWithTraceback(): Not enough DP cells allocated!\n");
				exit(1);
			}

			if (textPos > textDecodePos) {
			} else {
				// decode a word of text
				c = packedDNA[textDecodePos / CHAR_PER_WORD - 1];
				for (j=0; j<CHAR_PER_WORD; j++) {
					textBuffer[textOffset - textDecodePos + j + 1] = (unsigned char)(c & 0x3);
					c >>= BIT_PER_CHAR;
				}

				// Apply ambiguity
				while (ambiguity[ambiguityIndex].startPos >= textPos) {
					ambiguityIndex--;
				}
				if (ambiguity[ambiguityIndex].rightOfEndPos > (textPos - CHAR_PER_WORD)) {
					for (j=0; j<CHAR_PER_WORD; j++) {
						while (ambiguity[ambiguityIndex].startPos >= (textPos - j)) {
							ambiguityIndex--;
						}
						if (ambiguity[ambiguityIndex].rightOfEndPos > (textPos - j - 1)) {
							textBuffer[textOffset - textDecodePos + j + 1] = (char)ambiguity[ambiguityIndex].symbol;
						}
					}
				}

				textDecodePos -= CHAR_PER_WORD;
			}

			// traceback
			tracebackIndex[(textOffset - textPos + 1) * 2] = currentTracebackIndex;
			tracebackIndex[(textOffset - textPos + 1) * 2 + 1] = startPos;

			nextEndPos = 0;
			nextStartPos = INT_MAX;

			if (currentTracebackIndex + queryPatternLength >= tracebackAllocated) {
				traceback = MMUnitReallocate(traceback, tracebackAllocated + tracebackAllocationUnit, tracebackAllocated);
				tracebackAllocated += tracebackAllocationUnit;
			}

			if (startPos > lastStartPos) {
				BLast = B[startPos - 1];
				D = DP_NEG_INFINITY;
				DType = 0;
				currentPos = startPos;
			} else {
				// all scores in currentPos - 1 is DP_NEG_INFINITY
				BLast = B[startPos];
				B[startPos] = I[startPos];
				traceback[currentTracebackIndex] = DP_INSERT | IType[startPos];
				I[startPos] = I[startPos] + gapExtendScore;
				IType[startPos] = DP_INSERT_EXTEND;
				D = B[startPos] + gapExtendScore + gapOpenScore;
				DType = DP_DELETE_OPEN;
				if (B[startPos] + dropoffScore >= forwardMaxScore) {
					nextStartPos = startPos;
					nextEndPos = startPos;
				}
				currentPos = startPos + 1;
				currentTracebackIndex++;
			}


			while (currentPos <= endPos) {

				M = BLast + scoringMatrix[textBuffer[textOffset - textPos + 1]][convertedKey[queryOffset - currentPos]];
				BLast = B[currentPos];

				if (M >= I[currentPos] && M >= D) {
					// matchScore is maximum
					B[currentPos] = M;
					traceback[currentTracebackIndex] = DP_MATCH_MISMATCH | IType[currentPos] | DType;
					if (M + gapOpenScore >= I[currentPos]) {
						I[currentPos] = M + gapOpenScore + gapExtendScore;
						IType[currentPos] = DP_INSERT_OPEN;
					} else {
						I[currentPos] = I[currentPos] + gapExtendScore;
						IType[currentPos] = DP_INSERT_EXTEND;
					}
					if (M + gapOpenScore >= D) {
						D = M + gapOpenScore + gapExtendScore;
						DType = DP_DELETE_OPEN;
					} else {
						D = D + gapExtendScore;
						DType = DP_DELETE_EXTEND;
					}
					if (B[currentPos] > backwardMaxScore) {
						backwardMaxScore = B[currentPos];
						backwardMaxScorePosQuery =  queryOffset - currentPos;
						backwardMaxScorePosText = textPos - 1;
					}
				} else {
					if (I[currentPos] >= D) {
						B[currentPos] = I[currentPos];
						traceback[currentTracebackIndex] = DP_INSERT | IType[currentPos] | DType;
					} else {
						B[currentPos] = D;
						traceback[currentTracebackIndex] = DP_DELETE | IType[currentPos] | DType;
					}
					// insert cannot follow delete and delete cannot follow insert
					I[currentPos] = I[currentPos] + gapExtendScore;
					IType[currentPos] = DP_INSERT_EXTEND;
					D = D + gapExtendScore;
					DType = DP_DELETE_EXTEND;
				}

				if (B[currentPos] + dropoffScore >= backwardMaxScore) {
					if (nextStartPos > currentPos) {
						nextStartPos = currentPos;
					}
					nextEndPos = currentPos;
				}

				currentPos++;
				currentTracebackIndex++;

			}

			nextEndPos++;
			if (nextEndPos == currentPos) {
				while (nextEndPos <= queryOffset && D + dropoffScore >= backwardMaxScore)  {
					if (nextEndPos >= dpCellAllocated) {
						fprintf(stderr, "HSPGappedExtensionWithTraceback(): Not enough DP cells allocated!\n");
						exit(1);
					}
					B[nextEndPos] = D;
					traceback[currentTracebackIndex] = DP_DELETE | DType;
					D = D + gapExtendScore;
					DType = DP_DELETE_EXTEND;
					I[nextEndPos] = DP_NEG_INFINITY;
					IType[nextEndPos] = 0;
					nextEndPos++;
					currentTracebackIndex++;
				}
				I[nextEndPos] = DP_NEG_INFINITY;
				IType[nextEndPos] = 0;
			}

			if (nextEndPos > queryOffset) {
				nextEndPos = queryOffset;
			}

			lastStartPos = startPos;
			startPos = nextStartPos;
			endPos = nextEndPos;

			textPos--;

		}

		// traceback
		backwardNumOfAlignment = 0;
		backwardNumOfAuxiliaryText = 0;
		textPos = backwardMaxScorePosText + 1;
		currentPos = queryOffset - backwardMaxScorePosQuery;
		dpScore = 0;	// for verifying traceback

		currentTracebackIndex = tracebackIndex[(textOffset - textPos + 1) * 2] 
								+ currentPos - tracebackIndex[(textOffset - textPos + 1) * 2 + 1];

		while (currentTracebackIndex > 0) {

			dpType = traceback[currentTracebackIndex] & DP_MASK;

			if (dpType == DP_MATCH_MISMATCH) {
				dpScore += scoringMatrix[textBuffer[textOffset - textPos + 1]][convertedKey[queryOffset - currentPos]];
				if (textBuffer[textOffset - textPos + 1] == convertedKey[queryOffset - currentPos] &&
					textBuffer[textOffset - textPos + 1] < 4) {	// match and not ambiguity
					tempBackwardAlignment[backwardNumOfAlignment] = ALIGN_MATCH;
				} else {
					tempBackwardAlignment[backwardNumOfAlignment] = ALIGN_MISMATCH_AMBIGUITY;
					tempBackwardAuxiliaryText[backwardNumOfAuxiliaryText] = textBuffer[textOffset - textPos + 1];
					backwardNumOfAuxiliaryText++;
				}
				textPos++;
				currentPos--;
				currentTracebackIndex = tracebackIndex[(textOffset - textPos + 1) * 2] 
										+ currentPos - tracebackIndex[(textOffset - textPos + 1) * 2 + 1];
			} else if (dpType == DP_INSERT) {
				while (!(traceback[currentTracebackIndex] & DP_INSERT_OPEN)) {
					dpScore += gapExtendScore;
					tempBackwardAlignment[backwardNumOfAlignment] = ALIGN_INSERT;
					backwardNumOfAlignment++;
					tempBackwardAuxiliaryText[backwardNumOfAuxiliaryText] = textBuffer[textOffset - textPos + 1];
					backwardNumOfAuxiliaryText++;
					textPos++;
					currentTracebackIndex = tracebackIndex[(textOffset - textPos + 1) * 2] 
											+ currentPos - tracebackIndex[(textOffset - textPos + 1) * 2 + 1];
				}
				dpScore += gapOpenScore + gapExtendScore;
				tempBackwardAlignment[backwardNumOfAlignment] = ALIGN_INSERT;
				tempBackwardAuxiliaryText[backwardNumOfAuxiliaryText] = textBuffer[textOffset - textPos + 1];
				backwardNumOfAuxiliaryText++;
				textPos++;
				currentTracebackIndex = tracebackIndex[(textOffset - textPos + 1) * 2] 
										+ currentPos - tracebackIndex[(textOffset - textPos + 1) * 2 + 1];
			} else {
				while (!(traceback[currentTracebackIndex] & DP_DELETE_OPEN)) {
					dpScore += gapExtendScore;
					tempBackwardAlignment[backwardNumOfAlignment] = ALIGN_DELETE;
					backwardNumOfAlignment++;
					currentPos--;
					currentTracebackIndex--;
				}
				dpScore += gapOpenScore + gapExtendScore;
				tempBackwardAlignment[backwardNumOfAlignment] = ALIGN_DELETE;
				currentPos--;
				currentTracebackIndex--;
			}

			backwardNumOfAlignment++;

		}

		if (dpScore != backwardMaxScore) {
			fprintf(stderr, "Backward gapped extension traceback error!\n");
			exit(1);
		}

		// check cutoff score
		evalue = stat_gapCalcEvalue(stat_gapNominal2normalized(forwardMaxScore + backwardMaxScore));
		if (evalue < maxEvalue) {
			gappedHitList[numberOfGappedHit].posText = backwardMaxScorePosText;
			gappedHitList[numberOfGappedHit].posQuery = backwardMaxScorePosQuery;
			gappedHitList[numberOfGappedHit].lengthQuery = forwardMaxScorePosQuery - backwardMaxScorePosQuery;
			gappedHitList[numberOfGappedHit].score = forwardMaxScore + backwardMaxScore;
			gappedHitList[numberOfGappedHit].lengthText = forwardMaxScorePosText - backwardMaxScorePosText;
			gappedHitList[numberOfGappedHit].dbSeqIndex = gappedHitList[hitProcessed].dbSeqIndex;

			// Store alignment and auxiliary text
			numOfAlignmentWord = (forwardNumOfAlignment + backwardNumOfAlignment + ALIGN_PER_WORD - 1) / ALIGN_PER_WORD;
			((GappedHitListWithAlignment*)(gappedHitList + numberOfGappedHit))->alignmentOffset = MMPoolDispatchOffset(alignmentPool, numOfAlignmentWord * sizeof(unsigned int));
			alignment = (unsigned int*)((char*)alignmentPool + ((GappedHitListWithAlignment*)(gappedHitList + numberOfGappedHit))->alignmentOffset);
			for (i=0; i<numOfAlignmentWord; i++) {
				alignment[i] = 0;
			}
			if (forwardNumOfAuxiliaryText + backwardNumOfAuxiliaryText > 0) {
				numOfAuxiliaryTextWord = (forwardNumOfAuxiliaryText + backwardNumOfAuxiliaryText + AUX_TEXT_PER_WORD - 1) / AUX_TEXT_PER_WORD;
				((GappedHitListWithAlignment*)(gappedHitList + numberOfGappedHit))->auxiliaryTextOffset = MMPoolDispatchOffset(alignmentPool, numOfAuxiliaryTextWord * sizeof(unsigned int));
				auxiliaryText = (unsigned int*)((char*)alignmentPool + ((GappedHitListWithAlignment*)(gappedHitList + numberOfGappedHit))->auxiliaryTextOffset);
				for (i=0; i<numOfAuxiliaryTextWord; i++) {
					auxiliaryText[i] = 0;
				}
			} else {
				((GappedHitListWithAlignment*)(gappedHitList + numberOfGappedHit))->auxiliaryTextOffset = 0;
			}

			// backward alignment
			for (i=0; i<backwardNumOfAlignment; i++) {
				alignment[i/ALIGN_PER_WORD] |= (tempBackwardAlignment[i] << (BITS_IN_WORD - (i % ALIGN_PER_WORD + 1) * ALIGN_BIT));
			}
			for (i=0; i<backwardNumOfAuxiliaryText; i++) {
				auxiliaryText[i/AUX_TEXT_PER_WORD] |= (tempBackwardAuxiliaryText[i] << (BITS_IN_WORD - (i % AUX_TEXT_PER_WORD + 1) * AUX_TEXT_BIT));
			}

			// forward alignment
			for (i=0; i<forwardNumOfAlignment; i++) {
				alignment[(i+backwardNumOfAlignment)/ALIGN_PER_WORD] |= (tempForwardAlignment[forwardNumOfAlignment - i - 1] 
																		<< (BITS_IN_WORD - ((i+backwardNumOfAlignment) % ALIGN_PER_WORD + 1) * ALIGN_BIT));
			}
			for (i=0; i<forwardNumOfAuxiliaryText; i++) {
				auxiliaryText[(i+backwardNumOfAuxiliaryText)/AUX_TEXT_PER_WORD] |= (tempForwardAuxiliaryText[forwardNumOfAuxiliaryText - i - 1] 
																				   << (BITS_IN_WORD - ((i+backwardNumOfAuxiliaryText) % AUX_TEXT_PER_WORD + 1) * AUX_TEXT_BIT));
			}
			
			numberOfGappedHit++;

		}

		hitProcessed++;

	}

	// free working memory
	MMPoolReturn(mmPool, textBuffer, dpCellAllocated * 2);
	MMPoolReturn(mmPool, B, dpCellAllocated * sizeof(int));
	MMPoolReturn(mmPool, IType, dpCellAllocated);
	MMPoolReturn(mmPool, I, dpCellAllocated * sizeof(int));

	MMPoolReturn(mmPool, tracebackIndex, dpCellAllocated * 2 * 2 * sizeof(int));

	MMUnitFree(traceback, tracebackAllocated);

	MMPoolReturn(mmPool, tempForwardAlignment, dpCellAllocated * 2);
	MMPoolReturn(mmPool, tempBackwardAlignment, dpCellAllocated * 2);
	MMPoolReturn(mmPool, tempForwardAuxiliaryText, dpCellAllocated);
	MMPoolReturn(mmPool, tempBackwardAuxiliaryText, dpCellAllocated);

	return numberOfGappedHit;

}

// textList must be sorted in descending posText order;
// queryPosListIndex contains the starting positions of queryPos in queryPosList within a queryPosGroup
// a dummy entry must be present at the end of queryPosListIndex

int HSPDPDBvsQuery(const unsigned int *packedText, const HitList *textList, const int numOfTextHit,
				   const SeqOffset *textSeqOffset, const Ambiguity *textAmbiguity, const int numOfTextSeq,
				   const unsigned char *convertedKey, const int queryPatternLength, 
				   const int *queryPosListIndex, const unsigned int *queryPosList,
				   GappedHitListWithAlignment* __restrict gappedHitList, const int maxNumOfGappedHit,
				   MMPool *mmPool, MMPool *alignmentPool,
				   const int matchScore, const int mismatchScore,
				   const int gapOpenScore, const int gapExtendScore,
				   const int cutoffScore, const double maxEvalue) {

	#define NUM_SOURCE_BIT		12
	#define SOURCE_BIT_MASK		0xFFF

	#define DROPOFF_MAX_ENTRY		256

	#define EARLY_DROPOFF_MIN		11
	#define EARLY_DROPOFF_FACTOR	4	// These two are heuristic values to determine the dropoff to trigger addition traceback
										// When this heuistic is not applied, about 0.01% alignments reported by Blast are
										// filtered out during traceback on experiments with query size at chromosome lengths
										// The filtered alignments are of low evalues and are close to some strong alignments

	#define TEXT_BUFFER_SIZE	16

	// DP table
	DPCell** __restrict dpCell;
	int M;
	int D;				// insert and delete wrt query

	// The following variables are all indexed by source bit

	// Max score
	DPMaxScore* __restrict maxScore;
	DPMaxScore* __restrict dropoffMaxScore;

	// Keeping track whether dropoff should be valid alignments
	int* __restrict currentMaxScore;
	int* __restrict lastDropoffMaxScoreIndex;

	// End variables indexed by source bit

	int textBuffer[TEXT_BUFFER_SIZE];

	char* __restrict tempAlignment;
	char* __restrict tempAuxiliaryText;

	int scoringMatrix[16][16];
	int gapOpenScoreShifted, gapExtendScoreShifted;
	int cutoffScoreShifted;
	int dropoffScoreShifted;
	int minPositiveScore;

	int sourceBit;

	double evalue;

	int firstTextHitBeingProcessed, numOfTextHitBeingProcessed;

	unsigned int textOffset;
	int textChar;
	int charToProcess;

	int textDbSeqIndex, textAmbiguityIndex;
	unsigned int textSeqStart, textSeqEnd;

	int i, j;
	unsigned int c;

	int numOfGappedHit;

	int sourceBitMinScore, maxSourceBit, nextSourceBit, usedSourceBit;
	int numOfDropoffMaxScore;

	int	tempAlignmentIndex;
	int	tempAuxiliaryTextIndex;
	int numOfAlignmentWord, numOfAuxiliaryTextWord;

	Traceback* __restrict traceback;
	int tracebackAllocated;
	int tracebackIndex;
	int tracebackIndexAdjustment;

	unsigned int queryListInfo;

	int lastDpCellUsed, lastDpCellIndex;
	int dpCellUsed, dpCellIndex;
	int dpCellUsedAdjustedForEdge;

	int q;
	unsigned int qPos, lastQPos, maxQPos;
	int qChar;

	int bestScore, insertScore, deleteScore;

	int minAlignmentLength, minNegAlignmentLength, maxDPGroupingDist;

	unsigned int* __restrict alignment;
	unsigned int* __restrict auxiliaryText;

	// allocate working memory
	dpCell = MMPoolDispatch(mmPool, 2 * sizeof(DPCell*));
	dpCell[0] = MMPoolDispatch(mmPool, MAX_ALIGNMENT_LENGTH * sizeof(DPCell));
	dpCell[1] = MMPoolDispatch(mmPool, MAX_ALIGNMENT_LENGTH * sizeof(DPCell));

	maxScore = MMPoolDispatch(mmPool, (1 << NUM_SOURCE_BIT) * sizeof(DPMaxScore));
	dropoffMaxScore = MMPoolDispatch(mmPool, DROPOFF_MAX_ENTRY * sizeof(DPMaxScore));

	currentMaxScore = MMPoolDispatch(mmPool, (1 << NUM_SOURCE_BIT) * sizeof(int));
	lastDropoffMaxScoreIndex = MMPoolDispatch(mmPool, (1 << NUM_SOURCE_BIT) * sizeof(int));

	tempAlignment = MMPoolDispatch(mmPool, (MAX_ALIGNMENT_LENGTH * 2) * sizeof(char));
	tempAuxiliaryText = MMPoolDispatch(mmPool, (MAX_ALIGNMENT_LENGTH * 2) * sizeof(char));

	tracebackAllocated = MMPoolByteAvailable(mmPool) / sizeof(Traceback);
	if (tracebackAllocated < 0) {
		fprintf(stderr, "HSPDPDBvsDB(): Not enough memory allocated!\n");
		exit(1);
	}
	traceback =  MMPoolDispatch(mmPool, tracebackAllocated * sizeof(Traceback));

	HSPFillScoringMatrix(scoringMatrix, matchScore, mismatchScore, NUM_SOURCE_BIT);

	gapOpenScoreShifted = gapOpenScore * (1 << NUM_SOURCE_BIT);
	gapExtendScoreShifted = gapExtendScore * (1 << NUM_SOURCE_BIT);
	cutoffScoreShifted = cutoffScore * (1 << NUM_SOURCE_BIT);
	minPositiveScore = 1 << NUM_SOURCE_BIT;
	if (cutoffScore <= EARLY_DROPOFF_MIN) {
		dropoffScoreShifted = cutoffScoreShifted;
	} else {
		dropoffScoreShifted = (cutoffScore - (cutoffScore - EARLY_DROPOFF_MIN + EARLY_DROPOFF_FACTOR - 1) / EARLY_DROPOFF_FACTOR) * (1 << NUM_SOURCE_BIT) ;
	}

	minAlignmentLength = (cutoffScore + matchScore - 1) / matchScore;
	minNegAlignmentLength = ((gapOpenScore + gapExtendScore + matchScore + 1 - minAlignmentLength * matchScore + 1) / (gapOpenScore + gapExtendScore + matchScore)) * 2 - 1;
	maxDPGroupingDist = minAlignmentLength + minNegAlignmentLength;	// minimum length for an alignment to reached cutoff and then drop back to zero 

	// maxSourceBit are initially assigned to cells
	maxSourceBit = (1 << NUM_SOURCE_BIT) - 1;
	// source bits are assigned when score >= sourceBitMinScore
	sourceBitMinScore = (1 - gapOpenScore - gapExtendScore) * (1 << NUM_SOURCE_BIT) + maxSourceBit;

	numOfGappedHit = 0;		// Number of alignments 

	firstTextHitBeingProcessed = 0;

	// Clear last dropoff
	for (j=0; j<maxSourceBit; j++) {
		lastDropoffMaxScoreIndex[j] = 0;
	}
	numOfDropoffMaxScore = 0;
	usedSourceBit = 0;

	textOffset = 0;

	while (firstTextHitBeingProcessed < numOfTextHit) {

		if (textList[firstTextHitBeingProcessed].posText > textOffset) {
			// textOffset will be larger than those processed
			textDbSeqIndex = numOfTextSeq - 1;
			textSeqStart = textSeqOffset[textDbSeqIndex].startPos;
			textSeqEnd = textSeqOffset[textDbSeqIndex].endPos;
			textAmbiguityIndex = textSeqOffset[textDbSeqIndex].lastAmbiguityIndex;
		}

		textOffset = textList[firstTextHitBeingProcessed].posText;	// posText points to the character right after the hit
		queryListInfo = textList[firstTextHitBeingProcessed].info;

		// Find corresponding DB sequence for the text
		while (textSeqOffset[textDbSeqIndex].startPos >= textOffset) {
			textDbSeqIndex--;
			textSeqStart = textSeqOffset[textDbSeqIndex].startPos;
			textSeqEnd = textSeqOffset[textDbSeqIndex].endPos;
			textAmbiguityIndex = textSeqOffset[textDbSeqIndex].lastAmbiguityIndex;
		}

		if (textOffset % CHAR_PER_WORD > 0) {

			// Decode text to the next word boundary
			charToProcess = textOffset % CHAR_PER_WORD;
			c = packedText[textOffset / CHAR_PER_WORD] >> (BITS_IN_WORD - charToProcess * BIT_PER_CHAR);
			for (i=charToProcess; i--;) {	// from charToProcess - 1 to 0
				textBuffer[i] = (c & 0x3);
				c >>= BIT_PER_CHAR;
			}

			// Apply ambiguity
			while (textAmbiguity[textAmbiguityIndex].startPos >= textOffset) {
				textAmbiguityIndex--;
			}
			if (textAmbiguity[textAmbiguityIndex].rightOfEndPos > (textOffset - charToProcess)) {
				for (i=0; i<charToProcess; i++) {
					while (textAmbiguity[textAmbiguityIndex].startPos >= (textOffset - i)) {
						textAmbiguityIndex--;
					}
					if (textAmbiguity[textAmbiguityIndex].rightOfEndPos > (textOffset - i - 1)) {
						textBuffer[(textOffset - i - 1) % CHAR_PER_WORD] = textAmbiguity[textAmbiguityIndex].symbol;
					}
				}
			}

		}

		lastDpCellUsed = 0;
		lastDpCellIndex = 0;
		dpCellIndex = 1;
		maxQPos = 0;
		tracebackIndex = 0;
		tracebackIndexAdjustment = 0;

		// Clear last dropoff
		for (j=0; j<usedSourceBit; j++) {
			lastDropoffMaxScoreIndex[j] = 0;
		}
		numOfDropoffMaxScore = 0;

		numOfTextHitBeingProcessed = 0;
		nextSourceBit = 0;

		while (textOffset > textSeqOffset[textDbSeqIndex].startPos && 
			   (lastDpCellUsed > 0 || textList[firstTextHitBeingProcessed + numOfTextHitBeingProcessed].posText == textOffset)) {

			if (textOffset % CHAR_PER_WORD == 0) {

				// Decode text to the next word boundary
				c = packedText[textOffset / CHAR_PER_WORD - 1];
				for (i=CHAR_PER_WORD; i--;) {	// from CHAR_PER_WORD - 1 to 0
					textBuffer[i] = (c & 0x3);
					c >>= BIT_PER_CHAR;
				}

				// Apply ambiguity
				while (textAmbiguity[textAmbiguityIndex].startPos >= textOffset) {
					textAmbiguityIndex--;
				}
				if (textAmbiguity[textAmbiguityIndex].rightOfEndPos > (textOffset - CHAR_PER_WORD)) {
					for (i=0; i<CHAR_PER_WORD; i++) {
						while (textAmbiguity[textAmbiguityIndex].startPos >= (textOffset - i)) {
							textAmbiguityIndex--;
						}
						if (textAmbiguity[textAmbiguityIndex].rightOfEndPos > (textOffset - i - 1)) {
							textBuffer[(textOffset - i - 1) % CHAR_PER_WORD] = textAmbiguity[textAmbiguityIndex].symbol;
						}
					}
				}

			}

			// Set text char
			textChar = textBuffer[(textOffset - 1) % CHAR_PER_WORD];

			// Check if it is a valid starting position of alignment
			if (firstTextHitBeingProcessed + numOfTextHitBeingProcessed > numOfTextHit ||
				textList[firstTextHitBeingProcessed + numOfTextHitBeingProcessed].posText != textOffset ||
				textOffset + maxDPGroupingDist < textList[firstTextHitBeingProcessed].posText) {
//				(numOfTextHitBeingProcessed > 0 && textList[firstTextHitBeingProcessed + numOfTextHitBeingProcessed].posText + maxDPSkip <
//												   textList[firstTextHitBeingProcessed + numOfTextHitBeingProcessed - 1].posText)) {
			} else {
				// Merge starting positions in query into the dpCell
				queryListInfo = textList[firstTextHitBeingProcessed + numOfTextHitBeingProcessed].info;	// the group of query hits corresponding to the text hit
				q = queryPosListIndex[queryListInfo];												// the startingindex of the group of query hits

				dpCellUsed = 0;
				i = 0;
				while ((q < queryPosListIndex[queryListInfo + 1] ||
					    i < lastDpCellUsed) && dpCellUsed < MAX_ALIGNMENT_LENGTH) {
					while (i < lastDpCellUsed && (q >= queryPosListIndex[queryListInfo + 1] || dpCell[lastDpCellIndex][i].P >= queryPosList[q] 
											  && dpCellUsed < MAX_ALIGNMENT_LENGTH)) {
						dpCell[dpCellIndex][dpCellUsed] = dpCell[lastDpCellIndex][i];
						dpCellUsed++;
						i++;
					}
					if (q < queryPosListIndex[queryListInfo + 1] && i > 0 && dpCell[lastDpCellIndex][i-1].P == queryPosList[q]) {
						// The cell is positive; ignore starting position in query
						q++;
					}
					while (q < queryPosListIndex[queryListInfo + 1] && dpCellUsed < MAX_ALIGNMENT_LENGTH &&
						   (i >= lastDpCellUsed || dpCell[lastDpCellIndex][i].P < queryPosList[q])) {
						// insert starting position in query
						if (queryPosList[q] > (unsigned int)queryPatternLength) {
							fprintf(stderr, "HSPDPDBvsQuery(): Query position is larger than query length!\n");
							exit(1);
						}
						dpCell[dpCellIndex][dpCellUsed].P = queryPosList[q];		// posText points to the character right after the hit
						dpCell[dpCellIndex][dpCellUsed].B = maxSourceBit;			// initial score for a blank cell
						dpCell[dpCellIndex][dpCellUsed].I = DP_NEG_INFINITY;		// insertion not allowed here
						dpCellUsed++;
						tracebackIndexAdjustment--;
						q++;
					}
				}
				if (q < queryPosListIndex[queryListInfo + 1] || i < lastDpCellUsed) {
					// Not enough memory
				} else {
					numOfTextHitBeingProcessed++;
					dpCellIndex ^= 1;		// Swap dpCells
					lastDpCellIndex ^= 1;	// Swap dpCells
					lastDpCellUsed = dpCellUsed;
				}
			}

			D = DP_NEG_INFINITY;

			// Clear current max score
			usedSourceBit = nextSourceBit;
			for (j=0; j<usedSourceBit; j++) {
				currentMaxScore[j] = 0;
			}

			i = 0;
			dpCellUsed = 0;
			dpCellUsedAdjustedForEdge = FALSE;
			lastQPos = ALL_ONE_MASK;
			
			// check traceback buffer
			if (tracebackIndex + MAX_ALIGNMENT_LENGTH >= tracebackAllocated) {
				fprintf(stderr, "HSPDPDBvsQuery(): Not enough traceback buffer!\n");
				exit(1);
			}

			while (i < lastDpCellUsed) {

				// DP

				// check buffer
				if (dpCellUsed + 3 >= MAX_ALIGNMENT_LENGTH) {
					// 1 cell at most give 3 positive cells
					fprintf(stderr, "HSPDPDBvsQuery(): Not enough query buffer(1)!\n");
					exit(1);
				}

				if (dpCell[lastDpCellIndex][i].B == maxSourceBit) {
					// the cell is inserted through input query positions
					tracebackIndexAdjustment++;
				}

				qPos = dpCell[lastDpCellIndex][i].P - 1;

				if (qPos != lastQPos - 1) {
					// The cell on the left is non-positive; handle the insertion + deletion only
					if (dpCell[lastDpCellIndex][i].I >= minPositiveScore || D >= minPositiveScore) {

						traceback[tracebackIndex].indexDiff = dpCellUsed - i + lastDpCellUsed + tracebackIndexAdjustment;	// indexDiff always point to the cell directly above logically
						traceback[tracebackIndex].textChar = textChar;

						if (dpCell[lastDpCellIndex][i].I >= minPositiveScore && dpCell[lastDpCellIndex][i].I == dpCell[lastDpCellIndex][i].B + gapOpenScoreShifted + gapExtendScoreShifted) {
							traceback[tracebackIndex].IOpen = 1;		// Gap opening
						} else {
							traceback[tracebackIndex].IOpen = 0;		// Not gap opening
						}
						if (D >= minPositiveScore && D == dpCell[dpCellIndex][dpCellUsed-1].B + gapOpenScoreShifted + gapExtendScoreShifted) {
							traceback[tracebackIndex].DOpen = 1;		// Gap opening
						} else {
							traceback[tracebackIndex].DOpen = 0;		// Not gap opening
						}

						if (dpCell[lastDpCellIndex][i].I >= D) {
							dpCell[dpCellIndex][dpCellUsed].B = dpCell[lastDpCellIndex][i].I;
							traceback[tracebackIndex].alignment = DP_INSERT;
						} else {
							dpCell[dpCellIndex][dpCellUsed].B = D;
							traceback[tracebackIndex].alignment = DP_DELETE;
						}
						tracebackIndex++;

						dpCell[dpCellIndex][dpCellUsed].I = dpCell[lastDpCellIndex][i].I + gapExtendScoreShifted;	// insert cannot follow delete
						D = D + gapExtendScoreShifted;	// delete cannot follow insert 
						dpCell[dpCellIndex][dpCellUsed].P = qPos + 1;
						dpCellUsed++;
					}
				}
				
				qChar = convertedKey[qPos];

				// DP - Match/mismatch
				M = dpCell[lastDpCellIndex][i].B + scoringMatrix[textChar][qChar];

				if (M >= sourceBitMinScore &&
					(M & SOURCE_BIT_MASK) == maxSourceBit) {
					// Assign source bits
					if (nextSourceBit < maxSourceBit) {
						M += nextSourceBit - maxSourceBit;
						maxScore[nextSourceBit].score = 0;
						nextSourceBit++;
					} else {
						fprintf(stderr, "HSPDPDBvsQuery(): Not enough source bit space!\n");
					}
				}

				// Store the insert and delete score before updating them
				if (i+1 < lastDpCellUsed && dpCell[lastDpCellIndex][i+1].P == qPos) {
					// insert score available
					insertScore = dpCell[lastDpCellIndex][i+1].I;
				} else {
					insertScore = DP_NEG_INFINITY;
				}
				deleteScore = D;

				// Determine the best score
				if (M >= insertScore && M >= D) {
					// match score is maximum
					dpCell[dpCellIndex][dpCellUsed].B = M;
					if (M + gapOpenScoreShifted >= insertScore) {
						dpCell[dpCellIndex][dpCellUsed].I = M + gapOpenScoreShifted + gapExtendScoreShifted;
					} else {
						dpCell[dpCellIndex][dpCellUsed].I = insertScore + gapExtendScoreShifted;
					}
					if (M + gapOpenScoreShifted >= D) {
						D = M + gapOpenScoreShifted + gapExtendScoreShifted;
					} else {
						D = D + gapExtendScoreShifted;
					}
					traceback[tracebackIndex].alignment = DP_MATCH_MISMATCH;
				} else {
					if (insertScore >= D) {
						dpCell[dpCellIndex][dpCellUsed].B = insertScore;
						traceback[tracebackIndex].alignment = DP_INSERT;
					} else {
						dpCell[dpCellIndex][dpCellUsed].B = D;
						traceback[tracebackIndex].alignment = DP_DELETE;
					}
					dpCell[dpCellIndex][dpCellUsed].I = insertScore + gapExtendScoreShifted;
					D = D + gapExtendScoreShifted;
				}

				// check if the best score is positive
				if (dpCell[dpCellIndex][dpCellUsed].B >= minPositiveScore) {

					traceback[tracebackIndex].indexDiff = dpCellUsed - i -1 + lastDpCellUsed + tracebackIndexAdjustment;	// indexDiff always point to the cell directly above logically
					traceback[tracebackIndex].textChar = textChar;
					traceback[tracebackIndex].queryChar = qChar;
					if (insertScore >= minPositiveScore && insertScore == dpCell[lastDpCellIndex][i+1].B + gapOpenScoreShifted + gapExtendScoreShifted) {
						traceback[tracebackIndex].IOpen = 1;		// Gap opening
					} else {
						traceback[tracebackIndex].IOpen = 0;		// Not gap opening
					}
					if (deleteScore >= minPositiveScore && deleteScore == dpCell[dpCellIndex][dpCellUsed-1].B + gapOpenScoreShifted + gapExtendScoreShifted) {
						traceback[tracebackIndex].DOpen = 1;		// Gap opening
					} else {
						traceback[tracebackIndex].DOpen = 0;		// Not gap opening
					}

					dpCell[dpCellIndex][dpCellUsed].P = qPos;

					sourceBit = dpCell[dpCellIndex][dpCellUsed].B & SOURCE_BIT_MASK;
					if (dpCell[dpCellIndex][dpCellUsed].B > maxScore[sourceBit].score) {
						maxScore[sourceBit].score = dpCell[dpCellIndex][dpCellUsed].B;
						maxScore[sourceBit].tracebackIndex = tracebackIndex;
						maxScore[sourceBit].posText = textOffset - 1;
						maxScore[sourceBit].posQuery = qPos;
					}
					if (sourceBit < usedSourceBit && dpCell[dpCellIndex][dpCellUsed].B > currentMaxScore[sourceBit]) {
						currentMaxScore[sourceBit] = dpCell[dpCellIndex][dpCellUsed].B;
					}

					dpCellUsed++;
					tracebackIndex++;

				}

				lastQPos = qPos;

				if ((i+1 >= lastDpCellUsed || dpCell[lastDpCellIndex][i+1].P < (qPos-1)) && qPos > 0 && D >= minPositiveScore) {

					// delete score is still positive;
					qPos--;
					dpCell[dpCellIndex][dpCellUsed].B = D;
					dpCell[dpCellIndex][dpCellUsed].I = DP_NEG_INFINITY;
					dpCell[dpCellIndex][dpCellUsed].P = qPos;
					D = D + gapExtendScoreShifted;

					traceback[tracebackIndex].indexDiff = traceback[tracebackIndex-1].indexDiff;
					traceback[tracebackIndex].textChar = textChar;
					traceback[tracebackIndex].IOpen = 0;		// Not gap opening
					if (dpCell[dpCellIndex][dpCellUsed].B == dpCell[dpCellIndex][dpCellUsed-1].B + gapOpenScoreShifted + gapExtendScoreShifted) {
						traceback[tracebackIndex].DOpen = 1;		// Gap opening
					} else {
						traceback[tracebackIndex].DOpen = 0;		// Not gap opening
					}
					traceback[tracebackIndex].alignment = DP_DELETE;

					dpCellUsed++;
					tracebackIndex++;
				}
				if (qPos == 0) {
					// Reached the start of query
					D = DP_NEG_INFINITY;
					if (dpCellUsed > 0 && dpCell[dpCellIndex][dpCellUsed - 1].P == 0) {
						// No need to process the cell on the start of query
						dpCellUsed--;
						dpCellUsedAdjustedForEdge = TRUE;
					}
					break;
				}

				i++;

			}

			if (dpCellUsedAdjustedForEdge) {
				tracebackIndexAdjustment = 1;
			} else {
				tracebackIndexAdjustment = 0;
			}

			// Add to dropoff list if maxScore > currentMaxScore + cutoff
			for (j=0; j<usedSourceBit; j++) {
				if (maxScore[j].score >= currentMaxScore[j] + dropoffScoreShifted && lastDropoffMaxScoreIndex[j] != maxScore[j].tracebackIndex) {
					if (numOfDropoffMaxScore < DROPOFF_MAX_ENTRY) {
						dropoffMaxScore[numOfDropoffMaxScore] = maxScore[j];
						lastDropoffMaxScoreIndex[j] = maxScore[j].tracebackIndex;
						numOfDropoffMaxScore++;
					} else {
						fprintf(stderr, "HSPDPDBvsQuery(): Not enough dropoff space!\n");
					}
				}
			}

			dpCellIndex ^= 1;		// Swap dpCells
			lastDpCellIndex ^= 1;	// Swap dpCells
			lastDpCellUsed = dpCellUsed;

			textOffset--;

		}

		// Move scores reached cutoff to front
		j = 0;
		for (i=0; i<nextSourceBit; i++) {
			if (maxScore[i].score >= cutoffScoreShifted) {
				maxScore[j] = maxScore[i];
				j++;
			}
		}
		nextSourceBit = j;

		// Append dropoff max score
		for (i=0; i<numOfDropoffMaxScore; i++) {
			if (nextSourceBit < (1<<NUM_SOURCE_BIT)) {
				for (j=0; j<nextSourceBit; j++) {
					if (maxScore[j].tracebackIndex == dropoffMaxScore[i].tracebackIndex) {
						break;
					}
				}
				if (j >= nextSourceBit) {
					maxScore[nextSourceBit] = dropoffMaxScore[i];
					nextSourceBit++;
				}
			} else {
				fprintf(stderr, "HSPDPDBvsQuery(): Not enough source bit for dropoff!\n");
			}
		}

		// Traceback alignments
		for (sourceBit=0; sourceBit<nextSourceBit; sourceBit++) {
			
			bestScore = maxScore[sourceBit].score;
			evalue = stat_gapCalcEvalue(stat_gapNominal2normalized(bestScore >> NUM_SOURCE_BIT));
			if (evalue <= maxEvalue) {

				// Traceback
				tempAlignmentIndex = 0;
				tempAuxiliaryTextIndex = 0;

				textOffset = maxScore[sourceBit].posText;
				qPos = maxScore[sourceBit].posQuery;

				tracebackIndex = maxScore[sourceBit].tracebackIndex;

				while (bestScore >= minPositiveScore) {

					if (traceback[tracebackIndex].alignment == DP_MATCH_MISMATCH) {
						if (tempAlignmentIndex >= MAX_ALIGNMENT_LENGTH) {
							fprintf(stderr, "HSPDPDBvsQuery(): Not enough temp alignment buffer!\n");
							exit(1);
						}
						textChar = traceback[tracebackIndex].textChar;
						qChar = traceback[tracebackIndex].queryChar;
						bestScore -= scoringMatrix[textChar][qChar];
						if (textChar == qChar && textChar < 4) {	// match and not ambiguity
							tempAlignment[tempAlignmentIndex] = ALIGN_MATCH;
						} else {
							tempAlignment[tempAlignmentIndex] = ALIGN_MISMATCH_AMBIGUITY;
							tempAuxiliaryText[tempAuxiliaryTextIndex] = (char)traceback[tracebackIndex].textChar;	// store text character if mismatch or ambiguity
							tempAuxiliaryTextIndex++;
						}
						textOffset++;
						qPos++;
						tracebackIndex -= traceback[tracebackIndex].indexDiff + 1;	// indexDiff always refer to the cell directly above logically
					} else if (traceback[tracebackIndex].alignment == DP_INSERT) {
						while (!traceback[tracebackIndex].IOpen) {
							if (tempAlignmentIndex >= MAX_ALIGNMENT_LENGTH) {
								fprintf(stderr, "HSPDPDBvsQuery(): Not enough temp alignment buffer!\n");
								exit(1);
							}
							bestScore -= gapExtendScoreShifted;
							tempAlignment[tempAlignmentIndex] = ALIGN_INSERT;
							tempAlignmentIndex++;
							tempAuxiliaryText[tempAuxiliaryTextIndex] = (char)traceback[tracebackIndex].textChar;
							tempAuxiliaryTextIndex++;
							tracebackIndex -= traceback[tracebackIndex].indexDiff;	// indexDiff always refer to the cell directly above logically
							textOffset++;
						}
						if (tempAlignmentIndex >= MAX_ALIGNMENT_LENGTH) {
							fprintf(stderr, "HSPDPDBvsQuery(): Not enough temp alignment buffer!\n");
							exit(1);
						}
						bestScore -= gapOpenScoreShifted + gapExtendScoreShifted;
						tempAlignment[tempAlignmentIndex] = ALIGN_INSERT;
						tempAuxiliaryText[tempAuxiliaryTextIndex] = (char)traceback[tracebackIndex].textChar;
						tempAuxiliaryTextIndex++;
						tracebackIndex -= traceback[tracebackIndex].indexDiff;	// indexDiff always refer to the cell directly above logically
						textOffset++;
					} else {
						while (!traceback[tracebackIndex].DOpen) {
							if (tempAlignmentIndex >= MAX_ALIGNMENT_LENGTH) {
								fprintf(stderr, "HSPDPDBvsQuery(): Not enough temp alignment buffer!\n");
								exit(1);
							}
							bestScore -= gapExtendScoreShifted;
							tempAlignment[tempAlignmentIndex] = ALIGN_DELETE;
							tempAlignmentIndex++;
							tracebackIndex--;
							qPos++;
						}
						if (tempAlignmentIndex >= MAX_ALIGNMENT_LENGTH) {
							fprintf(stderr, "HSPDPDBvsQuery(): Not enough temp alignment buffer!\n");
							exit(1);
						}
						bestScore -= gapOpenScoreShifted + gapExtendScoreShifted;
						tempAlignment[tempAlignmentIndex] = ALIGN_DELETE;
						tracebackIndex--;
						qPos++;
					}

					tempAlignmentIndex++;

				}

				if (numOfGappedHit >= maxNumOfGappedHit) {
					fprintf(stderr, "HSPDPDBvsQuery(): Not enough gapped hit buffer!\n");
					exit(1);
				}
				gappedHitList[numOfGappedHit].posText = maxScore[sourceBit].posText;
				gappedHitList[numOfGappedHit].posQuery = maxScore[sourceBit].posQuery;
				gappedHitList[numOfGappedHit].lengthText = textOffset - maxScore[sourceBit].posText;
				gappedHitList[numOfGappedHit].lengthQuery = qPos - maxScore[sourceBit].posQuery;
				gappedHitList[numOfGappedHit].score = maxScore[sourceBit].score >> NUM_SOURCE_BIT;
				gappedHitList[numOfGappedHit].dbSeqIndex = textDbSeqIndex;

				numOfAlignmentWord = (tempAlignmentIndex + ALIGN_PER_WORD - 1) / ALIGN_PER_WORD;
				gappedHitList[numOfGappedHit].alignmentOffset = MMPoolDispatchOffset(alignmentPool, numOfAlignmentWord * sizeof(unsigned int));
				alignment = (unsigned int*)((char*)alignmentPool + gappedHitList[numOfGappedHit].alignmentOffset);
				for (i=0; i<numOfAlignmentWord; i++) {
					alignment[i] = 0;
				}

				for (i=0; i<tempAlignmentIndex; i++) {
					alignment[i/ALIGN_PER_WORD] |= (tempAlignment[i] << (BITS_IN_WORD - (i % ALIGN_PER_WORD + 1) * ALIGN_BIT));
				}

				if (tempAuxiliaryTextIndex > 0) {
					numOfAuxiliaryTextWord = (tempAuxiliaryTextIndex + AUX_TEXT_PER_WORD - 1) / AUX_TEXT_PER_WORD;
					gappedHitList[numOfGappedHit].auxiliaryTextOffset = MMPoolDispatchOffset(alignmentPool, numOfAuxiliaryTextWord * sizeof(unsigned int));
					auxiliaryText = (unsigned int*)((char*)alignmentPool + gappedHitList[numOfGappedHit].auxiliaryTextOffset);
					for (i=0; i<numOfAuxiliaryTextWord; i++) {
						auxiliaryText[i] = 0;
					}
					for (i=0; i<tempAuxiliaryTextIndex; i++) {
						auxiliaryText[i/AUX_TEXT_PER_WORD] |= (tempAuxiliaryText[i] << (BITS_IN_WORD - (i % AUX_TEXT_PER_WORD + 1) * AUX_TEXT_BIT));
					}
				} else {
					gappedHitList[numOfGappedHit].auxiliaryTextOffset = 0;
				}

				numOfGappedHit++;

			}
		}

		firstTextHitBeingProcessed += numOfTextHitBeingProcessed;

	}

	// free working memory
	MMPoolReturn(mmPool, traceback, tracebackAllocated * sizeof(int));

	MMPoolReturn(mmPool, tempAuxiliaryText, (MAX_ALIGNMENT_LENGTH * 2) * sizeof(char));
	MMPoolReturn(mmPool, tempAlignment, (MAX_ALIGNMENT_LENGTH * 2) * sizeof(char));

	MMPoolReturn(mmPool, lastDropoffMaxScoreIndex, (1 << NUM_SOURCE_BIT) * sizeof(int));
	MMPoolReturn(mmPool, currentMaxScore, (1 << NUM_SOURCE_BIT) * sizeof(int));

	MMPoolReturn(mmPool, dropoffMaxScore, DROPOFF_MAX_ENTRY * sizeof(DPMaxScore));
	MMPoolReturn(mmPool, maxScore, (1 << NUM_SOURCE_BIT) * sizeof(DPMaxScore));

	MMPoolReturn(mmPool, dpCell[0], MAX_ALIGNMENT_LENGTH * sizeof(DPCell));
	MMPoolReturn(mmPool, dpCell[1], MAX_ALIGNMENT_LENGTH * sizeof(DPCell));
	MMPoolReturn(mmPool, dpCell, 2 * sizeof(DPCell*));

	return numOfGappedHit;


}


unsigned int HSPRemoveDuplicateUngappedHit(HitListWithPosQuery* __restrict ungappedHitList, const unsigned int numberOfHit) {

	unsigned int i, j;

	QSort(ungappedHitList, numberOfHit, sizeof(HitListWithPosQuery), HitListPosTextQueryOrder);

	i = j = 0;
	while (i<numberOfHit) {
		ungappedHitList[j].posText = ungappedHitList[i].posText;
		ungappedHitList[j].posQuery = ungappedHitList[i].posQuery;
		i++;
		while (i<numberOfHit && ungappedHitList[i].posText == ungappedHitList[j].posText 
							 && ungappedHitList[i].posQuery == ungappedHitList[j].posQuery) {
			i++;
		}
		j++;
	}

	return j;

}

unsigned int HSPSplitCrossBoundaryGappedHit(GappedHitList* __restrict gappedHitList, const unsigned int numberOfHit, const unsigned int queryPatternLength,
								   const SeqOffset *seqOffset,
								   const int matchScore, const int cutoffScore) {

	unsigned int i;
	unsigned int numberOfSplitHit;
	unsigned int dbSeqIndex;
	unsigned int splitDbSeqIndex;
	unsigned int minHitLength;
	unsigned int oldHitLength;
	unsigned int newHitLength;
	unsigned int hitEnd;
	LONG diagonalPos;


	minHitLength = (cutoffScore + matchScore - 1) / matchScore;

	QSort(gappedHitList, numberOfHit, sizeof(GappedHitList), GappedHitListPosTextQueryLengthScoreOrder);

	i = 0;
	numberOfSplitHit = 0;
	dbSeqIndex = 0;

	while (i<numberOfHit) {

		// Find corresponding DB sequence
		while (seqOffset[dbSeqIndex].endPos < gappedHitList[i].posText) {
			dbSeqIndex++;
		}
		gappedHitList[i].dbSeqIndex = dbSeqIndex;
		diagonalPos = (LONG)gappedHitList[i].ungappedPosText - (LONG)gappedHitList[i].ungappedPosQuery;
		oldHitLength = gappedHitList[i].lengthText;

		splitDbSeqIndex = dbSeqIndex + 1;
		while (seqOffset[splitDbSeqIndex].startPos <= gappedHitList[i].posText + gappedHitList[i].lengthText - 1 - minHitLength) {

			// the hit span across 2 db sequences; and the portion on the right > minimum hit length
			
			gappedHitList[numberOfHit + numberOfSplitHit].posText = seqOffset[splitDbSeqIndex].startPos;
			if (gappedHitList[numberOfHit + numberOfSplitHit].posText - diagonalPos < queryPatternLength) {
				gappedHitList[numberOfHit + numberOfSplitHit].posQuery = (unsigned int)((LONG)gappedHitList[numberOfHit + numberOfSplitHit].posText - diagonalPos);
			} else {
				gappedHitList[numberOfHit + numberOfSplitHit].posQuery = queryPatternLength - 1;
			}
			hitEnd = gappedHitList[i].posText + gappedHitList[i].lengthText - 1;
			if (hitEnd > seqOffset[splitDbSeqIndex].endPos) {
				hitEnd = seqOffset[splitDbSeqIndex].endPos;
			}
			newHitLength = hitEnd - gappedHitList[numberOfHit + numberOfSplitHit].posText + 1;
			gappedHitList[numberOfHit + numberOfSplitHit].lengthText = newHitLength;
			if (gappedHitList[numberOfHit + numberOfSplitHit].lengthQuery >= oldHitLength - newHitLength) {
				gappedHitList[numberOfHit + numberOfSplitHit].lengthQuery += newHitLength - oldHitLength;
			} else {
				gappedHitList[numberOfHit + numberOfSplitHit].lengthQuery = 0;
			}

			if (gappedHitList[i].ungappedPosText < seqOffset[splitDbSeqIndex].startPos) {
				gappedHitList[numberOfHit + numberOfSplitHit].ungappedPosText = seqOffset[splitDbSeqIndex].startPos;
			} else if (gappedHitList[i].ungappedPosText > hitEnd) {
				gappedHitList[numberOfHit + numberOfSplitHit].ungappedPosText = hitEnd;
			} else {
				gappedHitList[numberOfHit + numberOfSplitHit].ungappedPosText = gappedHitList[i].ungappedPosText;
			}
			if (gappedHitList[numberOfHit + numberOfSplitHit].ungappedPosText < diagonalPos) {
				gappedHitList[numberOfHit + numberOfSplitHit].ungappedPosQuery = 0;
			} else if (gappedHitList[numberOfHit + numberOfSplitHit].ungappedPosText - diagonalPos >= queryPatternLength) {
				gappedHitList[numberOfHit + numberOfSplitHit].ungappedPosQuery = queryPatternLength - 1;
			} else {
				gappedHitList[numberOfHit + numberOfSplitHit].ungappedPosQuery = (unsigned int)((LONG)gappedHitList[numberOfHit + numberOfSplitHit].ungappedPosText - diagonalPos);
			}

			if (gappedHitList[i].ungappedPosText >= seqOffset[splitDbSeqIndex].startPos &&
				gappedHitList[i].ungappedPosText <= hitEnd) {
				gappedHitList[numberOfHit + numberOfSplitHit].score = gappedHitList[i].score;
			} else {
				gappedHitList[numberOfHit + numberOfSplitHit].score = 0;
			}
			gappedHitList[numberOfHit + numberOfSplitHit].dbSeqIndex = splitDbSeqIndex;

			numberOfSplitHit++;
			splitDbSeqIndex++;

		}

		hitEnd = gappedHitList[i].posText + gappedHitList[i].lengthText - 1;
		if (hitEnd > seqOffset[splitDbSeqIndex].endPos) {
			hitEnd = seqOffset[splitDbSeqIndex].endPos;
		}
		if (gappedHitList[i].posText < seqOffset[dbSeqIndex].startPos) {
			gappedHitList[i].posText = seqOffset[dbSeqIndex].startPos;
			if (gappedHitList[i].posText - diagonalPos < queryPatternLength) {
				gappedHitList[i].posQuery = (int)((LONG)gappedHitList[i].posText - diagonalPos);
			} else {
				gappedHitList[i].posQuery = queryPatternLength - 1;
			}
		}
		newHitLength = hitEnd - gappedHitList[i].posText + 1;
		gappedHitList[i].lengthText = newHitLength;
		if (gappedHitList[i].lengthQuery >= oldHitLength - newHitLength) {
			gappedHitList[i].lengthQuery += newHitLength - oldHitLength;
		} else {
			gappedHitList[i].lengthQuery = 0;
		}

		if (gappedHitList[i].ungappedPosText < seqOffset[dbSeqIndex].startPos) {
			gappedHitList[i].ungappedPosText = seqOffset[dbSeqIndex].startPos;
		} else if (gappedHitList[i].ungappedPosText > hitEnd) {
			gappedHitList[i].ungappedPosText = hitEnd;
		} else {
			gappedHitList[i].ungappedPosText = gappedHitList[i].ungappedPosText;
		}
		if (gappedHitList[i].ungappedPosText < diagonalPos) {
			gappedHitList[i].ungappedPosQuery = 0;
		} else if (gappedHitList[i].ungappedPosText - diagonalPos >= queryPatternLength) {
			gappedHitList[i].ungappedPosQuery = queryPatternLength - 1;
		} else {
			gappedHitList[i].ungappedPosQuery = (int)((LONG)gappedHitList[i].ungappedPosText - diagonalPos);
		}

		if (gappedHitList[i].ungappedPosText >= seqOffset[dbSeqIndex].startPos &&
			gappedHitList[i].ungappedPosText <= hitEnd) {
			gappedHitList[i].score = gappedHitList[i].score;
		} else {
			gappedHitList[i].score = 0;
		}
		gappedHitList[i].dbSeqIndex = dbSeqIndex;

		if (newHitLength >= minHitLength) {
			i++;
		}

	}

	return numberOfHit + numberOfSplitHit;

}
						 
unsigned int HSPRemoveDuplicateGappedHit(GappedHitList* __restrict gappedHitList, const unsigned int numberOfHit) {

	unsigned int i, j;

	QSort(gappedHitList, numberOfHit, sizeof(GappedHitList), GappedHitListPosTextQueryLengthScoreOrder);

	i = j = 0;
	while (i<numberOfHit) {
		gappedHitList[j].posText = gappedHitList[i].posText;
		gappedHitList[j].posQuery = gappedHitList[i].posQuery;
		gappedHitList[j].lengthText = gappedHitList[i].lengthText;
		gappedHitList[j].lengthQuery = gappedHitList[i].lengthQuery;
		gappedHitList[j].score = gappedHitList[i].score;
		gappedHitList[j].dbSeqIndex = gappedHitList[i].dbSeqIndex;
		gappedHitList[j].ungappedPosText = gappedHitList[i].ungappedPosText;
		gappedHitList[j].ungappedPosQuery = gappedHitList[i].ungappedPosQuery;
		i++;
		while (i<numberOfHit && gappedHitList[i].posText == gappedHitList[j].posText 
							 && gappedHitList[i].posQuery == gappedHitList[j].posQuery
							 && gappedHitList[i].lengthText == gappedHitList[j].lengthText
							 && gappedHitList[i].lengthQuery == gappedHitList[j].lengthQuery) {
			i++;
		}
		j++;
	}

	return j;

}
						 
// alignment and auxiliary text must be allocated through MMPool
// otherwise the pointer will be lost without freeing the memory
unsigned int HSPFinalFilter(GappedHitListWithAlignment* __restrict gappedHitList, const unsigned int numberOfHit) {

	unsigned int i, j;
	unsigned int hitRemaining;

	hitRemaining = numberOfHit;

	// First filter among hits with common starting points the longer and lower scoring hits

	QSort(gappedHitList, hitRemaining, sizeof(GappedHitList), GappedHitListPosTextQueryScoreLengthOrder);

	j = 0;
	while (j<hitRemaining) {
		if (gappedHitList[j].alignmentOffset != 0) {	// not marked as discarded
			i = j + 1;
			while (i<hitRemaining && gappedHitList[i].posText == gappedHitList[j].posText
								  && gappedHitList[i].posQuery == gappedHitList[j].posQuery) {
				if (gappedHitList[i].alignmentOffset != 0) {
					if (gappedHitList[i].lengthText >= gappedHitList[j].lengthText &&
						gappedHitList[i].lengthQuery >= gappedHitList[j].lengthQuery &&
						gappedHitList[i].score <= gappedHitList[j].score) {
						// discard hit
						gappedHitList[i].alignmentOffset = 0;
					}
				}
				i++;
			}
		}
		j++;
	}

	// Pack retained hit to the beginning
	i = j = 0;
	while (i<hitRemaining) {
		if (gappedHitList[i].alignmentOffset != 0) {	// not marked as discarded
			gappedHitList[j].posText = gappedHitList[i].posText;
			gappedHitList[j].posQuery = gappedHitList[i].posQuery;
			gappedHitList[j].lengthText = gappedHitList[i].lengthText;
			gappedHitList[j].lengthQuery = gappedHitList[i].lengthQuery;
			gappedHitList[j].score = gappedHitList[i].score;
			gappedHitList[j].dbSeqIndex = gappedHitList[i].dbSeqIndex;
			gappedHitList[j].alignmentOffset = gappedHitList[i].alignmentOffset;
			gappedHitList[j].auxiliaryTextOffset = gappedHitList[i].auxiliaryTextOffset;
			j++;
		}
		i++;
	}

	hitRemaining = j;

	// Second filter among hits with common ending points the longer and lower scoring hits

	QSort(gappedHitList, hitRemaining, sizeof(GappedHitList), GappedHitListEndTextQueryScoreLengthOrder);

	j = 0;
	while (j<hitRemaining) {
		if (gappedHitList[j].alignmentOffset != 0) {	// not marked as discarded
			i = j + 1;
			while (i<hitRemaining && gappedHitList[i].posText + gappedHitList[i].lengthText == gappedHitList[j].posText + gappedHitList[j].lengthText 
							  && gappedHitList[i].posQuery + gappedHitList[i].lengthQuery == gappedHitList[j].posQuery + gappedHitList[j].lengthQuery) {
				if (gappedHitList[i].alignmentOffset != 0) {
					if (gappedHitList[i].lengthText >= gappedHitList[j].lengthText &&
						gappedHitList[i].lengthQuery >= gappedHitList[j].lengthQuery &&
						gappedHitList[i].score <= gappedHitList[j].score) {
						// discard hit
						gappedHitList[i].alignmentOffset = 0;
					}
				}
				i++;
			}
		}
		j++;
	}

	// Pack retained hit to the beginning
	i = j = 0;
	while (i<hitRemaining) {
		if (gappedHitList[i].alignmentOffset != 0) {	// not marked as discarded
			gappedHitList[j].posText = gappedHitList[i].posText;
			gappedHitList[j].posQuery = gappedHitList[i].posQuery;
			gappedHitList[j].lengthText = gappedHitList[i].lengthText;
			gappedHitList[j].lengthQuery = gappedHitList[i].lengthQuery;
			gappedHitList[j].score = gappedHitList[i].score;
			gappedHitList[j].dbSeqIndex = gappedHitList[i].dbSeqIndex;
			gappedHitList[j].alignmentOffset = gappedHitList[i].alignmentOffset;
			gappedHitList[j].auxiliaryTextOffset = gappedHitList[i].auxiliaryTextOffset;
			j++;
		}
		i++;
	}

	hitRemaining = j;

	// Third filter hits with common start points

	QSort(gappedHitList, hitRemaining, sizeof(GappedHitList), GappedHitListPosTextQueryScoreLengthOrder);

	i = j = 0;
	while (i<hitRemaining) {
		gappedHitList[j].posText = gappedHitList[i].posText;
		gappedHitList[j].posQuery = gappedHitList[i].posQuery;
		gappedHitList[j].lengthText = gappedHitList[i].lengthText;
		gappedHitList[j].lengthQuery = gappedHitList[i].lengthQuery;
		gappedHitList[j].score = gappedHitList[i].score;
		gappedHitList[j].dbSeqIndex = gappedHitList[i].dbSeqIndex;
		gappedHitList[j].alignmentOffset = gappedHitList[i].alignmentOffset;
		gappedHitList[j].auxiliaryTextOffset = gappedHitList[i].auxiliaryTextOffset;
		i++;
		while (i<hitRemaining && gappedHitList[i].posText == gappedHitList[j].posText 
							  && gappedHitList[i].posQuery == gappedHitList[j].posQuery) {
			// discard hit
			gappedHitList[i].alignmentOffset = 0;
			i++;
		}
		j++;
	}

	hitRemaining = j;

	// Fourth filter hits with common end points

	QSort(gappedHitList, hitRemaining, sizeof(GappedHitListWithAlignment), GappedHitListEndTextQueryScoreLengthOrder);

	i = j = 0;
	while (i<hitRemaining) {
		gappedHitList[j].posText = gappedHitList[i].posText;
		gappedHitList[j].posQuery = gappedHitList[i].posQuery;
		gappedHitList[j].lengthText = gappedHitList[i].lengthText;
		gappedHitList[j].lengthQuery = gappedHitList[i].lengthQuery;
		gappedHitList[j].score = gappedHitList[i].score;
		gappedHitList[j].dbSeqIndex = gappedHitList[i].dbSeqIndex;
		gappedHitList[j].alignmentOffset = gappedHitList[i].alignmentOffset;
		gappedHitList[j].auxiliaryTextOffset = gappedHitList[i].auxiliaryTextOffset;
		i++;
		while (i<hitRemaining && gappedHitList[i].posText + gappedHitList[i].lengthText == gappedHitList[j].posText + gappedHitList[j].lengthText 
							  && gappedHitList[i].posQuery + gappedHitList[i].lengthQuery == gappedHitList[j].posQuery + gappedHitList[j].lengthQuery) {
			// discard hit
			gappedHitList[i].alignmentOffset = 0;
			i++;
		}
		j++;
	}

	hitRemaining = j;

	// Fifth filter enclosed hits with lower scores

	QSort(gappedHitList, hitRemaining, sizeof(GappedHitList), GappedHitListEnclosedScoreOrder);

	j = 0;
	while (j<hitRemaining) {
		if (gappedHitList[j].alignmentOffset != 0) {	// not marked as discarded
			i = j + 1;
			while (i<hitRemaining && gappedHitList[i].posText < gappedHitList[j].posText + gappedHitList[j].lengthText) {
				if (gappedHitList[i].alignmentOffset != 0) {
					if (gappedHitList[i].posText + gappedHitList[i].lengthText <= gappedHitList[j].posText + gappedHitList[j].lengthText &&
						gappedHitList[i].posQuery >= gappedHitList[j].posQuery &&
						gappedHitList[i].posQuery + gappedHitList[i].lengthQuery <= gappedHitList[j].posQuery + gappedHitList[j].lengthQuery &&
						gappedHitList[i].score <= gappedHitList[j].score) {
						// discard hit
						gappedHitList[i].alignmentOffset = 0;
					}
				}
				i++;
			}
		}
		j++;
	}

	// Pack retained hit to the beginning
	i = j = 0;
	while (i<hitRemaining) {
		if (gappedHitList[i].alignmentOffset != 0) {	// not marked as discarded
			gappedHitList[j].posText = gappedHitList[i].posText;
			gappedHitList[j].posQuery = gappedHitList[i].posQuery;
			gappedHitList[j].lengthText = gappedHitList[i].lengthText;
			gappedHitList[j].lengthQuery = gappedHitList[i].lengthQuery;
			gappedHitList[j].score = gappedHitList[i].score;
			gappedHitList[j].dbSeqIndex = gappedHitList[i].dbSeqIndex;
			gappedHitList[j].alignmentOffset = gappedHitList[i].alignmentOffset;
			gappedHitList[j].auxiliaryTextOffset = gappedHitList[i].auxiliaryTextOffset;
			j++;
		}
		i++;
	}

	hitRemaining = j;

	return hitRemaining;

}
						 
HSPUngappedExtLookupTable *HSPGenUngappedExtLookupTable(MMPool *mmPool, const int matchScore, const int mismatchScore) {

	HSPUngappedExtLookupTable *hspUngappedExtLookupTable;

	hspUngappedExtLookupTable = MMPoolDispatch(mmPool, (1 << 16) * sizeof(HSPUngappedExtLookupTable));

	HSPCalUngappedExtLookupTable(hspUngappedExtLookupTable, matchScore, mismatchScore);

	return hspUngappedExtLookupTable;
}

void HSPCalUngappedExtLookupTable(HSPUngappedExtLookupTable *hspUngappedExtLookupTable, const int matchScore, const int mismatchScore) {

	static const unsigned int oddBitMask = 0x55555555;
	static const unsigned int evenBitMask = 0xAAAAAAAA;

	unsigned int i, j;
	char maxScore, maxScorePos, finalScore;

	// In entries where even bit = match/mismatch, odd bit = 0	-> matchMismatchBitVector is even bit compacted (used for forward extension)
	// In entries where odd bit = match/mismatch, even bit = 0	-> matchMismatchBitVector is odd bit compacted and then reversed (used for backward extension)
	// When all bits are 0, both above case can be handled without problem
	hspUngappedExtLookupTable[0].matchMismatchBitVector = 0;
	for (i=1; i<(1<<16); i++) {
		hspUngappedExtLookupTable[i].matchMismatchBitVector = 0;
		if ((i & oddBitMask) == 0) {
			for (j=0; j<8; j++) {
				hspUngappedExtLookupTable[i].matchMismatchBitVector |= ((i >> (15-j*2)) & 1) << (7-j);
			}
		}
		if ((i & evenBitMask) == 0) {
			for (j=0; j<8; j++) {
				hspUngappedExtLookupTable[i].matchMismatchBitVector |= ((i >> (14-j*2)) & 1) << j;
			}
		}
		// entries where both odd and even bits contain 1 are not used
	}

	// For each 16 bit match/mismatch vector, calculate max score, position giving max score, and final score
	// The first match/mismatch is in the leftmost bit of the 16 bit vector
	// 0 is a match and 1 is a mismatch (from result of xor)
	for (i=0; i<(1<<16); i++) {

		maxScore = 0;
		maxScorePos = 0;	// 0 -> before the 1st character; 16 -> after the 16th character
		finalScore = 0;

		for (j=0; j<16; j++) {
			if (((i >> (15-j)) & 1) == 0) {
				// match
				finalScore = finalScore + (char)matchScore;
				// favor longer hit if scores are the same
				if (finalScore >= maxScore) {
					maxScore = finalScore;
					maxScorePos = (char)(j+1);
				}
			} else {
				// mismatch
				finalScore = finalScore + (char)mismatchScore;
			}
		}

		hspUngappedExtLookupTable[i].maxScore = maxScore;
		hspUngappedExtLookupTable[i].maxScorePos = maxScorePos;
		hspUngappedExtLookupTable[i].finalScore = finalScore;
	}

}

void HSPFreeUngappedExtLookupTable(MMPool *mmPool, HSPUngappedExtLookupTable* hspUngappedExtLookupTable) {

	MMPoolReturn(mmPool, hspUngappedExtLookupTable, (1 << 16) * sizeof(HSPUngappedExtLookupTable));

}

Histogram* HSPAllocateHistogram(double minEvalue, double maxEvalue) {

	Histogram *histogram;

	histogram = MMUnitAllocate(sizeof(Histogram));
	histogram->minEvalue = minEvalue;
	histogram->maxEvalue = maxEvalue;
	histogram->histogramSize = (int)(ceil(log10(maxEvalue)) - floor(log10(minEvalue))) + 1;
	histogram->count = MMUnitAllocate(histogram->histogramSize * sizeof(int));

	return histogram;

}

void HSPInitializeHistogram(Histogram* histogram) {

	int i;

	for (i=0; i<histogram->histogramSize; i++) {
		histogram->count[i] = 0;
	}

}

void HSPFreeHistogram(Histogram* histogram) {

	MMUnitFree(histogram->count, histogram->histogramSize * sizeof(int));
	MMUnitFree(histogram, histogram->histogramSize * sizeof(int));

}

void HSPCountEvalueToHistogram(Histogram* histogram, const GappedHitListWithAlignment* gappedHitList, const int numberOfHit, const int roundEvalue) {

	int i;
	int evalueIndex;
	double tempEvalue;

	for (i=0; i<numberOfHit; i++) {

		if (roundEvalue == TRUE) {
			tempEvalue = HSPCalculateAndPrintEValue(NULL, gappedHitList[i].score);
		} else {
			tempEvalue = stat_gapCalcEvalue(stat_gapNominal2normalized(gappedHitList[i].score));
		}
		if (tempEvalue == 0.0) {
			tempEvalue = histogram->minEvalue;
		}
		evalueIndex = (int)ceil(log10(tempEvalue) - floor(log10(histogram->minEvalue)));
		if (evalueIndex < 0) {
			evalueIndex = 0;
		}
		if (evalueIndex < histogram->histogramSize) {
			histogram->count[evalueIndex]++;
		}

	}

}
	
void HSPPrintHistogram(FILE * outputFile, Histogram* histogram) {

	int i;
	double tempEvalue;

	Socketfprintf(outputFile,"\nCount histogram of gapped extension by expectation value:\n");

	tempEvalue = pow(10, (int)floor(log10(histogram->minEvalue)));
	Socketfprintf(outputFile,"(     0,%6.2g] : %d\n", tempEvalue, histogram->count[0]);

	for (i=1; i<histogram->histogramSize; i++) {
		tempEvalue = pow(10, i - 1 + (int)log10(histogram->minEvalue));
		Socketfprintf(outputFile,"(%6.2g,%6.2g] : %d\n", tempEvalue, tempEvalue*10, histogram->count[i]);
	}

}

void HSPPrintHeader(FILE *outputFile, const int outputFormat, 
					const char* databaseName, const int numOfSeq, const unsigned int databaseLength) {

	switch (outputFormat) {

	case OUTPUT_PAIRWISE :

		fprintf(outputFile, "Database         : %s\n", databaseName);
		fprintf(outputFile, "No. of sequences : %d\n", numOfSeq);
		fprintf(outputFile, "No. of letters   : %u\n", databaseLength);
		fprintf(outputFile, "\n");

		break;

	case OUTPUT_TABULAR :

		break;

	case  OUTPUT_TABULAR_COMMENT :

		fprintf(outputFile, "# Database         : %s\n", databaseName);
		fprintf(outputFile, "# No. of sequences : %d\n", numOfSeq);
		fprintf(outputFile, "# No. of letters   : %u\n", databaseLength);
		fprintf(outputFile, "#\n");

		break;

	}

}

void HSPPrintTrailer(FILE *outputFile, const int outputFormat, 
					 const char* databaseName, const int numOfSeq, const unsigned int databaseLength) {

	switch (outputFormat) {

	case OUTPUT_PAIRWISE :

		fprintf(outputFile, "Database         : %s\n", databaseName);
		fprintf(outputFile, "No. of sequences : %d\n", numOfSeq);
		fprintf(outputFile, "No. of letters   : %u\n", databaseLength);
		fprintf(outputFile, "\n\n");
	
		printHSPstatistic(outputFile);
	
		break;

	case OUTPUT_TABULAR :

		break;

	case OUTPUT_TABULAR_COMMENT :

		break;

	}

}

void HSPPrintQueryHeader(FILE *outputFile, const int outputFormat,
					   const char* queryPatternName, const int queryPatternLength,
					   const int *dbOrder, const int *dbScore, const Annotation *annotation, const int numOfSeq) {

	int i;
	int charPrinted;
	char stringFormat[8] = "%00.00s";

	switch (outputFormat) {

	case OUTPUT_PAIRWISE :

		fprintf(outputFile, "Query            : %s\n", queryPatternName);
		fprintf(outputFile, "No. of letters   : %u\n", queryPatternLength);

		fprintf(outputFile, "\n\n");
		fprintf(outputFile, "                                                                 Score    E\n");
		fprintf(outputFile, "Sequences producing significant alignments:                      (bits) Value\n");
		fprintf(outputFile, "\n");

		for (i=0; i<numOfSeq && dbScore[dbOrder[i]] > 0; i++) {
			charPrinted = 0;
			if (annotation[dbOrder[i]].gi > 0) {
				charPrinted += fprintf(outputFile, "gi|%d|", annotation[dbOrder[i]].gi);
			}
			stringFormat[1] = stringFormat[4] = (char)('0' + (64 - charPrinted) / 10);
			stringFormat[2] = stringFormat[5] = (char)('0' + (64 - charPrinted) % 10);
			fprintf(outputFile, stringFormat, annotation[dbOrder[i]].text);
			fprintf(outputFile, " ");
			HSPCalculateAndPrintBitScore(outputFile, dbScore[dbOrder[i]]);
			fprintf(outputFile, "   ");
			HSPCalculateAndPrintEValue(outputFile, dbScore[dbOrder[i]]);
			fprintf(outputFile, "\n");
		}
		fprintf(outputFile, "\n\n");
	
		break;

	case OUTPUT_TABULAR :

		break;

	case OUTPUT_TABULAR_COMMENT :

		fprintf(outputFile, "# Query            : %s\n", queryPatternName);
		fprintf(outputFile, "# Query id, Subject id, %% identity, alignment length, mismatches, gap openings, q. start, q. end, s. start, s. end, e-value, bit score\n");

	}

}

void HSPPrintAlignment(MMPool *mmPool, FILE *outputFile, MMPool *alignmentPool, const GappedHitListWithAlignment* gappedHitList, const int numOfHit, 
					   int outputFormat, const ContextInfo* contextInfo, 
					   const unsigned char *charMap, const unsigned char *complementMap,
					   const char* queryPatternName, const unsigned char* convertedQueryPattern, const int queryPatternLength,
					   const int * dbOrder, const HSP *hsp) {

	int i, j;
	int lastDbSeqIndex, dbSeqIndex, contextNum;
	char* __restrict tempDbChar;
	char* __restrict tempQueryChar;

	char queryName[MAX_SEQ_NAME_LENGTH+1];
	char seqName[MAX_SEQ_NAME_LENGTH+1];

	int matchMatrix[16][16];

	unsigned int queryIndex, dbIndex;
	int totalLength, auxiliaryTextIndex;
	unsigned int c;
	int a, lastA;

	int alignmentCharPrinted;
	int textPosNumOfChar;
	char positionFormat[6] = "%-00u";

	char lowercase[256];

	int numOfMatch, numOfMismatch, numOfSpace, numOfGap;

	unsigned int* __restrict alignment;
	unsigned int* __restrict auxiliaryText;

	// Allocate memory
	tempDbChar = MMPoolDispatch(mmPool, MAX_ALIGNMENT_LENGTH + 1);
	tempQueryChar = MMPoolDispatch(mmPool, MAX_ALIGNMENT_LENGTH + 1);

	// duplicate query name
	strncpy(queryName, queryPatternName, MAX_SEQ_NAME_LENGTH+1);

	// take the first white space as end of pattern name
	for (i=0; i<MAX_SEQ_NAME_LENGTH; i++) {
		if (queryName[i] == ' ' || queryName[i] == '\0' || queryName[i] == '\t' || queryName[i] == '\n') {
			queryName[i] = '\0';
			break;
		}
	}

	HSPFillScoringMatrix(matchMatrix, 1, -1, 0);

	// initialize lowercase
	for (i=0; i<256; i++) {
		lowercase[i] = (char)i;
	}
	for (i=0; i<26; i++) {
		lowercase['A' + i] = 'a' + (char)i;
	}

	// Process all hits
	lastDbSeqIndex = -1;
	for (i=0; i<numOfHit; i++) {

		if (gappedHitList[i].lengthQuery > MAX_ALIGNMENT_LENGTH || gappedHitList[i].lengthText > MAX_ALIGNMENT_LENGTH) {
			fprintf(stderr, "HSPPrintAlignment(): Alignment length > maximum!\n");
			exit(1);
		}

		dbSeqIndex = dbOrder[gappedHitList[i].dbSeqIndex & CONTEXT_MASK];
		contextNum = gappedHitList[i].dbSeqIndex >> (BITS_IN_WORD - CONTEXT_BIT);

		if (dbSeqIndex != lastDbSeqIndex) {

			switch(outputFormat) {
			case OUTPUT_PAIRWISE :

				// Print header for sequence
				if (hsp->annotation[dbSeqIndex].gi > 0) {
					fprintf(outputFile, ">gi|%d|%s", hsp->annotation[dbSeqIndex].gi, hsp->annotation[dbSeqIndex].text);
				} else {
					fprintf(outputFile, ">%s", hsp->annotation[dbSeqIndex].text);
				}
				fprintf(outputFile, "No. of letters   : %u\n", hsp->seqOffset[dbSeqIndex].endPos - hsp->seqOffset[dbSeqIndex].startPos + 1);
				fprintf(outputFile, "\n");

				break;

			case OUTPUT_TABULAR :
			case OUTPUT_TABULAR_COMMENT :

				// duplicate sequence name
				strncpy(seqName, hsp->annotation[dbSeqIndex].text, MAX_SEQ_NAME_LENGTH+1);
				// take the first white space as end of pattern name
				for (j=0; j<MAX_SEQ_NAME_LENGTH; j++) {
					if (seqName[j] == ' ' || seqName[j] == '\0' || seqName[j] == '\t' || seqName[j] == '\n') {
						seqName[j] = '\0';
						break;
					}
				}

				break;

			}

		}
		lastDbSeqIndex = dbSeqIndex;

		// Decode alignment and DB characters

		queryIndex = 0;
		totalLength = 0;
		auxiliaryTextIndex = 0;
		numOfMatch = 0;
		numOfMismatch = 0;
		numOfSpace = 0;
		numOfGap = 0;
		lastA = -1;
		c = 0;

		alignment = (unsigned int*)((char*)alignmentPool + gappedHitList[i].alignmentOffset);
		auxiliaryText = (unsigned int*)((char*)alignmentPool + gappedHitList[i].auxiliaryTextOffset);

		while (queryIndex < gappedHitList[i].lengthQuery) {
			if (totalLength % ALIGN_PER_WORD == 0) {
				c = alignment[totalLength / ALIGN_PER_WORD];
			}
			a = (char)(c >> (BITS_IN_WORD - ALIGN_BIT));
			c <<= ALIGN_BIT;
			switch(a) {
			case ALIGN_MATCH :
				if (contextInfo[contextNum].reversed) {
					tempQueryChar[totalLength] = dnaChar[convertedQueryPattern[queryPatternLength - 1 - queryIndex - gappedHitList[i].posQuery]];
				} else {
					tempQueryChar[totalLength] = dnaChar[convertedQueryPattern[queryIndex + gappedHitList[i].posQuery]];
				}
				tempDbChar[totalLength] = tempQueryChar[totalLength];
				numOfMatch++;
				queryIndex++;
				break;
			case ALIGN_MISMATCH_AMBIGUITY :
				tempDbChar[totalLength] = dnaChar[(auxiliaryText[auxiliaryTextIndex / AUX_TEXT_PER_WORD] << ((auxiliaryTextIndex % AUX_TEXT_PER_WORD) * AUX_TEXT_BIT) >> (BITS_IN_WORD - AUX_TEXT_BIT))];
				if (contextInfo[contextNum].reversed) {
					tempDbChar[totalLength] = complementMap[tempDbChar[totalLength]];
					tempQueryChar[totalLength] = dnaChar[convertedQueryPattern[queryPatternLength - 1 - queryIndex - gappedHitList[i].posQuery]];
				} else {
					tempQueryChar[totalLength] = dnaChar[convertedQueryPattern[queryIndex + gappedHitList[i].posQuery]];
				}
				if (matchMatrix[charMap[tempDbChar[totalLength]]][charMap[tempQueryChar[totalLength]]] <= 0) {
					// Mismatch
					numOfMismatch++;
				} else {
					// Match with ambiguity
					if (tempDbChar[totalLength] == tempQueryChar[totalLength] && charMap[tempDbChar[totalLength]] < ALPHABET_SIZE) {
						fprintf(stderr, "HSPPrintAlignment(): Alignment error!\n");
						exit(1);
					}
					numOfMatch++;
				}
				queryIndex++;
				auxiliaryTextIndex++;
				break;
			case ALIGN_INSERT :
				tempDbChar[totalLength] = dnaChar[(auxiliaryText[auxiliaryTextIndex / AUX_TEXT_PER_WORD] << ((auxiliaryTextIndex % AUX_TEXT_PER_WORD) * AUX_TEXT_BIT) >> (BITS_IN_WORD - AUX_TEXT_BIT))];
				if (contextInfo[contextNum].reversed) {
					tempDbChar[totalLength] = complementMap[tempDbChar[totalLength]];
				}
				tempQueryChar[totalLength] = '-';
				numOfSpace++;
				if (lastA != ALIGN_INSERT) {
					numOfGap++;
				}
				auxiliaryTextIndex++;
				break;
			case ALIGN_DELETE :
				tempDbChar[totalLength] = '-';
				if (contextInfo[contextNum].reversed) {
					tempQueryChar[totalLength] = dnaChar[convertedQueryPattern[queryPatternLength - 1 - queryIndex - gappedHitList[i].posQuery]];
				} else {
					tempQueryChar[totalLength] = dnaChar[convertedQueryPattern[queryIndex + gappedHitList[i].posQuery]];
				}
				numOfSpace++;
				if (lastA != ALIGN_DELETE) {
					numOfGap++;
				}
				queryIndex++;
				break;
			}
			lastA = a;
			totalLength++;
		}

		// Print alignment

		switch (outputFormat) {

		case OUTPUT_PAIRWISE :

			fprintf(outputFile, "Score = ");
			HSPCalculateAndPrintBitScore(outputFile, gappedHitList[i].score);
			fprintf(outputFile, " bits (%d), Expect = ", gappedHitList[i].score);
			HSPCalculateAndPrintEValue(outputFile, gappedHitList[i].score);
			fprintf(outputFile, "\n");
			fprintf(outputFile, "Identities = %d/%d (%.2f%%), Gaps = %d/%d (%.2f%%)\n", numOfMatch, totalLength, (float)numOfMatch / totalLength * 100, numOfSpace, totalLength, (float)numOfSpace / totalLength * 100);
			if (contextInfo[contextNum].reversed) {
				fprintf(outputFile, "Strand = Plus / Minus\n");
			} else {
				fprintf(outputFile, "Strand = Plus / Plus\n");
			}
			fprintf(outputFile, "\n\n");

			// Determine position length and format
			c = gappedHitList[i].posText + gappedHitList[i].lengthText - hsp->seqOffset[dbSeqIndex].startPos;
			textPosNumOfChar = 0;
			while (c > 0) {
				c /= 10;
				textPosNumOfChar++;
			}
			if (textPosNumOfChar >= 10) {
				positionFormat[2] = '0' + (char)(textPosNumOfChar / 10);
				positionFormat[3] = '0' + (char)(textPosNumOfChar % 10);
				positionFormat[4] = 'u';
				positionFormat[5] = '\0';
			} else {
				positionFormat[2] = '0' + (char)textPosNumOfChar;
				positionFormat[3] = 'u';
				positionFormat[4] = '\0';
			}

			alignmentCharPrinted = 0;
			if (contextInfo[contextNum].reversed) {
				queryIndex = queryPatternLength + 1 - gappedHitList[i].posQuery - gappedHitList[i].lengthQuery;
				dbIndex = gappedHitList[i].posText + gappedHitList[i].lengthText - hsp->seqOffset[dbSeqIndex].startPos;
			} else {
				queryIndex = gappedHitList[i].posQuery + 1;
				dbIndex = gappedHitList[i].posText - hsp->seqOffset[dbSeqIndex].startPos + 1;
			}

			while (alignmentCharPrinted < totalLength) {

				if (contextInfo[contextNum].reversed) {

					// Print query line
					fprintf(outputFile, "Query: ");
					fprintf(outputFile, positionFormat, queryIndex);
					fprintf(outputFile, " ");

					j = 0;
					while (j < 60 && alignmentCharPrinted + j < totalLength) {
						fprintf(outputFile, "%c", lowercase[tempQueryChar[totalLength - 1 - alignmentCharPrinted - j]]);
						if (tempQueryChar[totalLength - 1 - alignmentCharPrinted - j] != '-') {
							queryIndex++;
						}
						j++;
					}

					fprintf(outputFile, " ");
					fprintf(outputFile, "%u", queryIndex - 1);
					fprintf(outputFile, "\n");

					// Print alignment line
					for (j=0; j<textPosNumOfChar + 8; j++) {
						fprintf(outputFile, " ");
					}

					j = 0;
					while (j < 60 && alignmentCharPrinted + j < totalLength) {
						if (tempQueryChar[totalLength - 1 - alignmentCharPrinted - j] == tempDbChar[totalLength - 1 - alignmentCharPrinted - j]) {
							fprintf(outputFile, "|");
						} else {
							fprintf(outputFile, " ");
						}
						j++;
					}
					fprintf(outputFile, "\n");

					// Print database line
					fprintf(outputFile, "Sbjct: ");
					fprintf(outputFile, positionFormat, dbIndex);
					fprintf(outputFile, " ");

					j = 0;
					while (j < 60 && alignmentCharPrinted + j < totalLength) {
						fprintf(outputFile, "%c", lowercase[tempDbChar[totalLength - 1 - alignmentCharPrinted - j]]);
						if (tempDbChar[totalLength - 1 - alignmentCharPrinted - j] != '-') {
							dbIndex--;
						}
						j++;
					}

					fprintf(outputFile, " ");
					fprintf(outputFile, positionFormat, dbIndex + 1);
					fprintf(outputFile, "\n");

				} else {

					// Print query line
					fprintf(outputFile, "Query: ");
					fprintf(outputFile, positionFormat, queryIndex);
					fprintf(outputFile, " ");

					j = 0;
					while (j < 60 && alignmentCharPrinted + j < totalLength) {
						fprintf(outputFile, "%c", lowercase[tempQueryChar[alignmentCharPrinted + j]]);
						if (tempQueryChar[alignmentCharPrinted + j] != '-') {
							queryIndex++;
						}
						j++;
					}

					fprintf(outputFile, " ");
					fprintf(outputFile, "%u", queryIndex - 1);
					fprintf(outputFile, "\n");

					// Print alignment line
					for (j=0; j<textPosNumOfChar + 8; j++) {
						fprintf(outputFile, " ");
					}

					j = 0;
					while (j < 60 && alignmentCharPrinted + j < totalLength) {
						if (matchMatrix[charMap[tempQueryChar[alignmentCharPrinted + j]]][charMap[tempDbChar[alignmentCharPrinted + j]]] > 0) {
							fprintf(outputFile, "|");
						} else {
							fprintf(outputFile, " ");
						}
						j++;
					}
					fprintf(outputFile, "\n");

					// Print database line
					fprintf(outputFile, "Sbjct: ");
					fprintf(outputFile, positionFormat, dbIndex);
					fprintf(outputFile, " ");

					j = 0;
					while (j < 60 && alignmentCharPrinted + j < totalLength) {
						fprintf(outputFile, "%c", lowercase[tempDbChar[alignmentCharPrinted + j]]);
						if (tempDbChar[alignmentCharPrinted + j] != '-') {
							dbIndex++;
						}
						j++;
					}

					fprintf(outputFile, " ");
					fprintf(outputFile, positionFormat, dbIndex - 1);
					fprintf(outputFile, "\n");

				}
					
				fprintf(outputFile, "\n\n");
				alignmentCharPrinted += j;

			}

			break;

		case OUTPUT_TABULAR :
		case OUTPUT_TABULAR_COMMENT :
			// print query name
			fprintf(outputFile, "%s\t", queryName);

			// print gi
			if (hsp->annotation[dbSeqIndex].gi > 0) {
				fprintf(outputFile, "gi|%d|", hsp->annotation[dbSeqIndex].gi);
			}

			// print db sequence name
			fprintf(outputFile, "%s\t", seqName);

			// print percentage identity
			fprintf(outputFile, "%.2f\t", (float)(totalLength - numOfMismatch - numOfSpace) / (float)(totalLength) * 100);

			// print align length
			fprintf(outputFile, "%d\t", totalLength);

			// print number of mismatch
			fprintf(outputFile, "%d\t", numOfMismatch);

			// print number of gap opening
			fprintf(outputFile, "%d\t", numOfGap);

			// print query start and end
			if (contextInfo[contextNum].reversed) {
				fprintf(outputFile, "%u\t%u\t", queryPatternLength - gappedHitList[i].posQuery - gappedHitList[i].lengthQuery + 1,
											   queryPatternLength - gappedHitList[i].posQuery);
			} else {
				fprintf(outputFile, "%u\t%u\t", gappedHitList[i].posQuery + 1,
											   gappedHitList[i].posQuery + gappedHitList[i].lengthQuery);
			}

			// print text start and end
			if (contextInfo[contextNum].reversed) {
				fprintf(outputFile, "%u\t%u\t", gappedHitList[i].posText + gappedHitList[i].lengthText - hsp->seqOffset[dbSeqIndex].startPos,
											   gappedHitList[i].posText + 1 - hsp->seqOffset[dbSeqIndex].startPos);
			} else {
				fprintf(outputFile, "%u\t%u\t", gappedHitList[i].posText + 1 - hsp->seqOffset[dbSeqIndex].startPos,
											   gappedHitList[i].posText + gappedHitList[i].lengthText - hsp->seqOffset[dbSeqIndex].startPos);
			}


			// print evalue
			HSPCalculateAndPrintEValue(outputFile, gappedHitList[i].score);

			// print bit score
			HSPCalculateAndPrintBitScore(outputFile, gappedHitList[i].score);

			fprintf(outputFile, "\n");

			break;

		default :

			fprintf(stderr, "HSPPrintAlignment() : Output format is invalid!\n");
			exit(1);

		}

	}

	// Free memory
	MMPoolReturn(mmPool, tempDbChar, MAX_ALIGNMENT_LENGTH + 1);
	MMPoolReturn(mmPool, tempQueryChar, MAX_ALIGNMENT_LENGTH + 1);

}

void HSPPrintNoAlignment(FILE *outputFile, int outputFormat) {

	switch (outputFormat) {

	case OUTPUT_PAIRWISE :

		break;

	case OUTPUT_TABULAR :
	case OUTPUT_TABULAR_COMMENT :

		fprintf(outputFile, "\n");

	}

}

double HSPCalculateAndPrintBitScore(FILE *outputFile, const int score) {

	double bitScore;

	// print bit score
	bitScore = stat_gapNominal2normalized(score);

	if (outputFile != NULL) {
		if (bitScore > 9999) {
			fprintf(outputFile, "%4.3le", bitScore);
		} else if (bitScore > 99.9) {
			fprintf(outputFile, "%4.0ld", (int)bitScore);
		} else {
			fprintf(outputFile, "%4.1lf", bitScore);
		}
	}

	return bitScore;

}

double HSPCalculateAndPrintEValue(FILE *outputFile, const int score) {

	double evalue;
	char evalueString[10];

	evalue = stat_gapCalcEvalue(stat_gapNominal2normalized(score));

	if (evalue < 1.0e-180) {
		sprintf(evalueString, "0.0\t");
	} else if (evalue < 1.0e-99) {
		sprintf(evalueString, "%2.0le\t", evalue);
	} else if (evalue < 0.0009) {
		sprintf(evalueString, "%3.0le\t", evalue);
	} else if (evalue < 0.1) {
		sprintf(evalueString, "%4.3lf\t", evalue);
	} else if (evalue < 1.0) { 
		sprintf(evalueString, "%3.2lf\t", evalue);
	} else if (evalue < 10.0) {
		sprintf(evalueString, "%2.1lf\t", evalue);
	} else { 
		sprintf(evalueString, "%5.0lf\t", evalue);
	}

	if (outputFile != NULL) {
		fprintf(outputFile, "%s\t", evalueString);
	}

	sscanf(evalueString, "%g", &evalue);
	return evalue;

}

int HitListPosTextQueryOrder(const void *hitList, const int index1, const int index2) {

	if (((HitListWithPosQuery*)hitList + index1)->posText != ((HitListWithPosQuery*)hitList + index2)->posText) {
		if (((HitListWithPosQuery*)hitList + index1)->posText > ((HitListWithPosQuery*)hitList + index2)->posText) {
			return 1;
		} else {
			return -1;
		}
	} else {
		return ((HitListWithPosQuery*)hitList + index1)->posQuery - ((HitListWithPosQuery*)hitList + index2)->posQuery;
	}

}

int GappedHitListPosTextQueryLengthScoreOrder(const void *gappedHitList, const int index1, const int index2) {

	if (((GappedHitList*)gappedHitList + index1)->posText != ((GappedHitList*)gappedHitList + index2)->posText) {
		if (((GappedHitList*)gappedHitList + index1)->posText > ((GappedHitList*)gappedHitList + index2)->posText) {
			return 1;
		} else {
			return -1;
		}
	} else {
		if (((GappedHitList*)gappedHitList + index1)->posQuery != ((GappedHitList*)gappedHitList + index2)->posQuery) {
			if (((GappedHitList*)gappedHitList + index1)->posQuery > ((GappedHitList*)gappedHitList + index2)->posQuery) {
				return 1;
			} else {
				return -1;
			}
		} else {
			if (((GappedHitList*)gappedHitList + index1)->lengthText != ((GappedHitList*)gappedHitList + index2)->lengthText) {
				return ((GappedHitList*)gappedHitList + index1)->lengthText - ((GappedHitList*)gappedHitList + index2)->lengthText;
			} else {
				if (((GappedHitList*)gappedHitList + index1)->lengthQuery != ((GappedHitList*)gappedHitList + index2)->lengthQuery) {
					return ((GappedHitList*)gappedHitList + index1)->lengthQuery - ((GappedHitList*)gappedHitList + index2)->lengthQuery;
				} else {
					// higher score first
					return ((GappedHitList*)gappedHitList + index2)->score - ((GappedHitList*)gappedHitList + index1)->score;
				}
			}
		}
	}

}

int GappedHitListPosTextQueryScoreLengthOrder(const void *gappedHitList, const int index1, const int index2) {

	if (((GappedHitList*)gappedHitList + index1)->posText != ((GappedHitList*)gappedHitList + index2)->posText) {
		if (((GappedHitList*)gappedHitList + index1)->posText > ((GappedHitList*)gappedHitList + index2)->posText) {
			return 1;
		} else {
			return -1;
		}
	} else {
		if (((GappedHitList*)gappedHitList + index1)->posQuery != ((GappedHitList*)gappedHitList + index2)->posQuery) {
			if (((GappedHitList*)gappedHitList + index1)->posQuery > ((GappedHitList*)gappedHitList + index2)->posQuery) {
				return 1;
			} else {
				return -1;
			}
		} else {
			if (((GappedHitList*)gappedHitList + index1)->score != ((GappedHitList*)gappedHitList + index2)->score) {
				// higher score first
				return ((GappedHitList*)gappedHitList + index2)->score - ((GappedHitList*)gappedHitList + index1)->score;
			} else {
				// short query length first
				if (((GappedHitList*)gappedHitList + index1)->lengthQuery != ((GappedHitList*)gappedHitList + index2)->lengthQuery) {
					if (((GappedHitList*)gappedHitList + index1)->lengthQuery > ((GappedHitListWithAlignment*)gappedHitList + index2)->lengthQuery) {
						return 1;
					} else {
						return -1;
					}
				} else {
					// short text length first
					if (((GappedHitList*)gappedHitList + index1)->lengthText != ((GappedHitListWithAlignment*)gappedHitList + index2)->lengthText) {
						if (((GappedHitList*)gappedHitList + index1)->lengthText > ((GappedHitListWithAlignment*)gappedHitList + index2)->lengthText) {
							return 1;
						} else {
							return -1;
						}
					} else {
						return 0;
					}
				}
			}
		}
	}

}

int GappedHitListEndTextQueryScoreLengthOrder(const void *gappedHitList, const int index1, const int index2) {

	if (((GappedHitList*)gappedHitList + index1)->posText + ((GappedHitList*)gappedHitList + index1)->lengthText !=
		((GappedHitList*)gappedHitList + index2)->posText + ((GappedHitList*)gappedHitList + index2)->lengthText) {
		if (((GappedHitList*)gappedHitList + index1)->posText + ((GappedHitList*)gappedHitList + index1)->lengthText >
			((GappedHitList*)gappedHitList + index2)->posText + ((GappedHitList*)gappedHitList + index2)->lengthText) {
			return 1;
		} else {
			return -1;
		}
	} else {
		if (((GappedHitList*)gappedHitList + index1)->posQuery + ((GappedHitList*)gappedHitList + index1)->lengthQuery !=
			((GappedHitList*)gappedHitList + index2)->posQuery + ((GappedHitList*)gappedHitList + index2)->lengthQuery) {
			if (((GappedHitList*)gappedHitList + index1)->posQuery + ((GappedHitList*)gappedHitList + index1)->lengthQuery >
				((GappedHitList*)gappedHitList + index2)->posQuery + ((GappedHitList*)gappedHitList + index2)->lengthQuery) {
				return 1;
			} else {
				return -1;
			}
		} else {
			if (((GappedHitList*)gappedHitList + index1)->score != ((GappedHitList*)gappedHitList + index2)->score) {
				// higher score first
				return ((GappedHitList*)gappedHitList + index2)->score - ((GappedHitList*)gappedHitList + index1)->score;
			} else {
				// short query length first
				if (((GappedHitList*)gappedHitList + index1)->lengthQuery != ((GappedHitListWithAlignment*)gappedHitList + index2)->lengthQuery) {
					if (((GappedHitList*)gappedHitList + index1)->lengthQuery > ((GappedHitListWithAlignment*)gappedHitList + index2)->lengthQuery) {
						return 1;
					} else {
						return -1;
					}
				} else {
					// short text length first
					if (((GappedHitList*)gappedHitList + index1)->lengthText != ((GappedHitListWithAlignment*)gappedHitList + index2)->lengthText) {
						if (((GappedHitList*)gappedHitList + index1)->lengthText > ((GappedHitListWithAlignment*)gappedHitList + index2)->lengthText) {
							return 1;
						} else {
							return -1;
						}
					} else {
						return 0;
					}
				}
			}
		}
	}

}

int GappedHitListEnclosedScoreOrder(const void *gappedHitList, const int index1, const int index2) {

	if (((GappedHitList*)gappedHitList + index1)->posText != ((GappedHitList*)gappedHitList + index2)->posText) {
		if (((GappedHitList*)gappedHitList + index1)->posText >	((GappedHitList*)gappedHitList + index2)->posText) {
			return 1;
		} else {
			return -1;
		}
	} else {
		if (((GappedHitList*)gappedHitList + index1)->lengthText != ((GappedHitList*)gappedHitList + index2)->lengthText) {
			if (((GappedHitList*)gappedHitList + index1)->lengthText < ((GappedHitList*)gappedHitList + index2)->lengthText) {
				return 1;
			} else {
				return -1;
			}
		} else {
			if (((GappedHitList*)gappedHitList + index1)->posQuery != ((GappedHitList*)gappedHitList + index2)->posQuery) {
				if (((GappedHitList*)gappedHitList + index1)->posQuery > ((GappedHitList*)gappedHitList + index2)->posQuery) {
					return 1;
				} else {
					return -1;
				}
			} else {
				if (((GappedHitList*)gappedHitList + index1)->lengthQuery != ((GappedHitList*)gappedHitList + index2)->lengthQuery) {
					if (((GappedHitList*)gappedHitList + index1)->lengthQuery < ((GappedHitList*)gappedHitList + index2)->lengthQuery) {
						return 1;
					} else {
						return -1;
					}
				} else {
					// higher score first
					return ((GappedHitList*)gappedHitList + index2)->score - ((GappedHitList*)gappedHitList + index1)->score;
				}
			}
		}
	}

}
#else

/*

   HSP.c		BWTBlastn functions

   This module contains miscellaneous BWTBlastn functions.

   Copyright (C) 2004, Wong Chi Kwong.

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
#include <string.h>
#include <math.h>
#include <stdint.h>
#include "TextConverter.h"
#include "MiscUtilities.h"
#include "Socket.h"
#include "r250.h"
#include "HSP.h"
#include "HSPstatistic.h"


extern  double stat_expectationValue;

void HSPFillCharMap(unsigned char charMap[255]) {

	int i;

	for (i=0; i<255; i++) {
		charMap[i] = nonMatchDnaCharIndex;
	}
	for (i=0; i<16; i++) {
		charMap[dnaChar[i]] = (unsigned char)i;
		charMap[dnaChar[i] - 'A' + 'a'] = (unsigned char)i;
	}

}

void HSPFillComplementMap(unsigned char complementMap[255]) {

	int i;

	for (i=0; i<255; i++) {
		complementMap[i] = dnaChar[nonMatchDnaCharIndex];
	}
	for (i=0; i<16; i++) {
		complementMap[dnaComplement[i]] = dnaChar[i];
		complementMap[dnaComplement[i] - 'A' + 'a'] = dnaChar[i] - 'A' + 'a';
	}

}

void HSPFillScoringMatrix(int scoringMatrix[16][16], const int matchScore, const int mismatchScore, const int leftShift) {

	int i, j, k, l;

	for (i=0; i<16; i++) {
		for (j=0; j<16; j++) {
			if (ambiguityCount[i] > 0 && ambiguityCount[j] > 0) {
				scoringMatrix[i][j] = 0;
				for (k=0; k<ambiguityCount[i]; k++) {
					for (l=0; l<ambiguityCount[j]; l++) {
						if (ambiguityMatch[i][k] == ambiguityMatch[j][l]) {
							// match
							scoringMatrix[i][j] += matchScore;
						} else {
							// mismatch
							scoringMatrix[i][j] += mismatchScore;
						}
					}
				}
				if (scoringMatrix[i][j] > 0) {
					scoringMatrix[i][j] = (scoringMatrix[i][j] + ambiguityCount[i] * ambiguityCount[j]  - 1) / (ambiguityCount[i] * ambiguityCount[j]);
				} else if (scoringMatrix[i][j] < 0) {
					scoringMatrix[i][j] = -((-scoringMatrix[i][j] + ambiguityCount[i] * ambiguityCount[j] - 1) / (ambiguityCount[i] * ambiguityCount[j]));
				} else {
					scoringMatrix[i][j] = 0;
				}
			} else {
				// Character meant to match nothing
				scoringMatrix[i][j] = mismatchScore;
			}
			scoringMatrix[i][j] = scoringMatrix[i][j] * (1 << leftShift);
		}
	}


}

HSP *HSPLoad(MMPool *mmPool, const char *PackedDNAFileName, const char *AnnotationFileName, const char *AmbiguityFileName) {

	HSP *hsp;
	unsigned int dnaLength;
	unsigned int randomSeed;
	int i;
	char c;
	unsigned char charMap[255];

	FILE *annotationFile = NULL, *ambiguityFile = NULL;

	hsp = MMPoolDispatch(mmPool, sizeof(HSP));

	// Load packed DNA
	if (PackedDNAFileName != NULL && PackedDNAFileName[0] != '\0' && PackedDNAFileName[0] != '-') {
		hsp->packedDNA = DNALoadPacked(PackedDNAFileName, &hsp->dnaLength, TRUE);
	} else {
		hsp->packedDNA = NULL;
		hsp->dnaLength = 0;
	}

	// Load annotation
	if (AnnotationFileName != NULL && AnnotationFileName[0] != '\0' && AnnotationFileName[0] != '-') {

		annotationFile = (FILE*)fopen64(AnnotationFileName, "r");
		if (annotationFile == NULL) {
			fprintf(stderr, "Cannot open annotation file!\n");
			exit(1);
		}

		fscanf(annotationFile, "%u %d %u\n", &dnaLength, &hsp->numOfSeq, &randomSeed);
		if (hsp->dnaLength != 0 && dnaLength != hsp->dnaLength) {
			fprintf(stderr, "Annotation database length not match!\n");
			exit(1);
		}
		hsp->dnaLength = dnaLength;
		if (hsp->numOfSeq == 0) {
			fprintf(stderr, "Annotation number of sequence = 0!\n");
			exit(1);
		}
		hsp->annotation = MMUnitAllocate((hsp->numOfSeq+1) * sizeof(Annotation));
		hsp->seqOffset = MMUnitAllocate((hsp->numOfSeq+1) * sizeof(SeqOffset));

		i = 0;
		hsp->minSeqLength = UINT32_MAX;
		while (!feof(annotationFile) && i < hsp->numOfSeq) {
			fscanf(annotationFile, "%u ", &hsp->annotation[i].gi);
			fgets(hsp->annotation[i].text, MAX_SEQ_NAME_LENGTH, annotationFile);
			fscanf(annotationFile, "%u %u %d\n", &hsp->seqOffset[i].startPos, &hsp->seqOffset[i].endPos, &hsp->seqOffset[i+1].firstAmbiguityIndex);
			hsp->seqOffset[i].lastAmbiguityIndex = hsp->seqOffset[i+1].firstAmbiguityIndex;
			if (hsp->seqOffset[i].endPos < hsp->minSeqLength) {
				hsp->minSeqLength = hsp->seqOffset[i].endPos;
			}
			hsp->seqOffset[i].endPos = hsp->seqOffset[i].startPos + hsp->seqOffset[i].endPos - 1;	// length of sequence is stored
			i++;
		}
		if (i < hsp->numOfSeq) {
			fprintf(stderr, "Annotation missing entries!\n");
			exit(1);
		}
		fclose(annotationFile);

		hsp->annotation[i].gi = 0;
		hsp->annotation[i].text[0] = '\0';

		hsp->seqOffset[i].startPos = UINT32_MAX;
		hsp->seqOffset[i].endPos = UINT32_MAX;
		hsp->seqOffset[0].firstAmbiguityIndex = 1;	// ambiguity[0] and ambiguity[numOfAmbiguity+1] are dummy
		for (i=1; i<=hsp->numOfSeq; i++) {
			hsp->seqOffset[i].firstAmbiguityIndex += hsp->seqOffset[i-1].firstAmbiguityIndex;	// number of ambiguity is stored
		}
		// hsp->seqOffset[hsp->numOfSeq].firstAmbiguityIndex = total number of ambiguity + 1 now
		hsp->seqOffset[hsp->numOfSeq].lastAmbiguityIndex = hsp->seqOffset[hsp->numOfSeq].firstAmbiguityIndex - 1;
		// hsp->seqOffset[hsp->numOfSeq].lastAmbiguityIndex = total number of ambiguity now
		for (i=hsp->numOfSeq; i>1; i--) {
			hsp->seqOffset[i-1].lastAmbiguityIndex = hsp->seqOffset[i].lastAmbiguityIndex - hsp->seqOffset[i-1].lastAmbiguityIndex;	// number of ambiguity is stored
		}
		for (i=0; i<hsp->numOfSeq; i++) {
			hsp->seqOffset[i].lastAmbiguityIndex = hsp->seqOffset[i+1].lastAmbiguityIndex;
		}

	} else {

		// Make up a dummy annotation and a dummy sequence offset
		hsp->annotation = MMUnitAllocate((1+1) * sizeof(Annotation));
		hsp->seqOffset = MMUnitAllocate((1+1) * sizeof(SeqOffset));

		hsp->numOfSeq = 1;
		hsp->numOfAmbiguity = 0;
		
		hsp->annotation[0].gi = 0;
		hsp->annotation[0].text[0] = '\0';
		hsp->annotation[1].gi = 0;
		hsp->annotation[1].text[0] = '\0';
		hsp->seqOffset[0].startPos = 0;
		hsp->seqOffset[0].endPos = hsp->dnaLength - 1;
		hsp->seqOffset[0].firstAmbiguityIndex = 1;
		hsp->seqOffset[0].lastAmbiguityIndex = 0;
		hsp->seqOffset[1].startPos = UINT32_MAX;
		hsp->seqOffset[1].endPos = UINT32_MAX;
		hsp->seqOffset[1].firstAmbiguityIndex = UINT32_MAX;	// should not be referred
		hsp->seqOffset[1].lastAmbiguityIndex = UINT32_MAX;	// should not be referred

	}

	// Load ambigity
	if (AmbiguityFileName != NULL && AmbiguityFileName[0] != '\0' && AmbiguityFileName[0] != '-') {

		// Load ambigity
		ambiguityFile = (FILE*)fopen64(AmbiguityFileName, "r");
		if (ambiguityFile == NULL) {
			fprintf(stderr, "Cannot open ambiguity file!\n");
			exit(1);
		}

		fscanf(ambiguityFile, "%u %d %d\n", &dnaLength, &i, &hsp->numOfAmbiguity);
		if (hsp->dnaLength != 0 && dnaLength != hsp->dnaLength) {
			fprintf(stderr, "Ambiguity database length not match!\n");
			exit(1);
		}
		hsp->dnaLength = dnaLength;
		if (i != hsp->numOfSeq) {
			fprintf(stderr, "Ambiguity database number of sequence not match!\n");
			exit(1);
		}
		if (hsp->numOfAmbiguity != hsp->seqOffset[hsp->numOfSeq].firstAmbiguityIndex - 1) {
			fprintf(stderr, "Ambiguity number of ambiguity not match!\n");
			exit(1);
		}

		HSPFillCharMap(charMap);

		hsp->ambiguity = MMUnitAllocate((hsp->numOfAmbiguity + 2) * sizeof(Ambiguity));
		hsp->ambiguity[0].startPos = 0;
		hsp->ambiguity[0].rightOfEndPos = 0;
		hsp->ambiguity[0].symbol = 0;
		i = 1;
		while (!feof(ambiguityFile) && i-1 < hsp->numOfAmbiguity) {
			fscanf(ambiguityFile, "%u %u %c\n", &hsp->ambiguity[i].startPos, &hsp->ambiguity[i].rightOfEndPos, &c);
			hsp->ambiguity[i].rightOfEndPos = hsp->ambiguity[i].startPos + hsp->ambiguity[i].rightOfEndPos;	// number of character is stored
			hsp->ambiguity[i].symbol = charMap[c];
			i++;
		}
		hsp->ambiguity[i].startPos = UINT32_MAX;
		hsp->ambiguity[i].rightOfEndPos = UINT32_MAX;
		hsp->ambiguity[i].symbol = 0;

		if (i-1 < hsp->numOfAmbiguity) {
			fprintf(stderr, "Ambiguity missing entries!\n");
			exit(1);
		}
		fclose(ambiguityFile);

	} else {

		if (hsp->numOfAmbiguity > 0) {
			fprintf(stderr, "Ambiguity file missing!\n");
			exit(1);
		}

		hsp->ambiguity = MMUnitAllocate((0 + 2) * sizeof(Ambiguity));
		hsp->ambiguity[0].startPos = 0;
		hsp->ambiguity[0].rightOfEndPos = 0;
		hsp->ambiguity[0].symbol = 0;
		hsp->ambiguity[1].startPos = UINT32_MAX;
		hsp->ambiguity[1].rightOfEndPos = UINT32_MAX;
		hsp->ambiguity[1].symbol = 0;

	}

	return hsp;

}

HSP *HSPConvertFromText(MMPool *mmPool, const unsigned char *text, const unsigned int textLength, 
						const unsigned int FASTARandomSeed, const int maskLowerCase,
						const int gi, const char *seqName) {

	HSP *hsp;
	Ambiguity *tempAmbiguity;
	int ambiguityAllocated =  256;
	
	char c;
	unsigned int i;
	int numAmbiguity;
	unsigned int lastAmbiguityPos;
	unsigned int numCharInBuffer, numOfConversion;
	unsigned int sequenceRandomSeed;
	unsigned char charMap[255];

	unsigned char buffer[PACKED_BUFFER_SIZE];

	HSPFillCharMap(charMap);

	hsp = MMPoolDispatch(mmPool, sizeof(HSP));

	hsp->dnaLength = textLength;
	hsp->minSeqLength = textLength;
	hsp->numOfSeq = 1;

	hsp->annotation = MMPoolDispatch(mmPool, (1+1) * sizeof(Annotation));
	hsp->annotation[0].gi = gi;
	strncpy(hsp->annotation[0].text, seqName, MAX_SEQ_NAME_LENGTH);
	hsp->annotation[1].gi = 0;
	hsp->annotation[1].text[0] = '\0';

	hsp->seqOffset = MMPoolDispatch(mmPool, (1+1) * sizeof(SeqOffset));
	hsp->seqOffset[0].startPos = 0;
	hsp->seqOffset[0].endPos = textLength - 1;
	hsp->seqOffset[0].firstAmbiguityIndex = 1;
	hsp->seqOffset[0].lastAmbiguityIndex = 0;
	hsp->seqOffset[1].startPos = UINT32_MAX;
	hsp->seqOffset[1].endPos = UINT32_MAX;
	hsp->seqOffset[1].firstAmbiguityIndex = UINT32_MAX;	// should not be referred
	hsp->seqOffset[1].lastAmbiguityIndex = UINT32_MAX;	// should not be referred

	hsp->packedDNA = MMPoolDispatch(mmPool, (textLength / CHAR_PER_WORD + 1) * sizeof(int));
	
	tempAmbiguity = MMUnitAllocate(ambiguityAllocated * sizeof(Ambiguity));

	HSPFillCharMap(charMap);

	// Set random seed for the sequence
	sequenceRandomSeed = FASTARandomSeed;
	if (gi > 0) {
		sequenceRandomSeed += gi;
	} else {
		for (i=0; i<(int)strlen(seqName); i++) {
			c = seqName[i];
			sequenceRandomSeed += c;
		}
	}
	r250_init(sequenceRandomSeed);

	numAmbiguity = -1;
	numCharInBuffer = 0;
	numOfConversion = 0;
	lastAmbiguityPos = (unsigned int)-2;

	for (i=0; i<textLength; i++) {

		c = text[i];

		if (maskLowerCase && c >= 'a' && c <= 'z') {
			c = dnaChar[lowercaseDnaCharIndex];
		}
		if (ambiguityCount[charMap[c]] == 1) {
			buffer[numCharInBuffer] = c;
		} else {
			c = charMap[c];
			if (ambiguityCount[c] > 0) {
				buffer[numCharInBuffer] = dnaChar[ambiguityMatch[c][r250() % ambiguityCount[c]]];
			} else {
				buffer[numCharInBuffer] = dnaChar[ambiguityMatch[c][r250() % ALPHABET_SIZE]];
			}
			if (i == lastAmbiguityPos + 1 && c == (char)tempAmbiguity[numAmbiguity].symbol) {
				tempAmbiguity[numAmbiguity].rightOfEndPos++;
			} else {
				numAmbiguity++;
				if (numAmbiguity >= ambiguityAllocated) {
					tempAmbiguity = MMUnitReallocate(tempAmbiguity, sizeof(Ambiguity) * ambiguityAllocated * 2, sizeof(Ambiguity) * ambiguityAllocated);
					ambiguityAllocated *= 2;
				}
				tempAmbiguity[numAmbiguity].startPos = i;
				tempAmbiguity[numAmbiguity].rightOfEndPos = i+1;
				tempAmbiguity[numAmbiguity].symbol = c;
			}
			lastAmbiguityPos = i;
		}
		numCharInBuffer++;
		if (numCharInBuffer >= PACKED_BUFFER_SIZE) {
			ConvertTextToWordPacked(buffer, hsp->packedDNA + numOfConversion * PACKED_BUFFER_SIZE / CHAR_PER_WORD, charMap, 4, PACKED_BUFFER_SIZE);
			numCharInBuffer = 0;
			numOfConversion++;
		}
	}
	if (numCharInBuffer > 0) {
		ConvertTextToWordPacked(buffer, hsp->packedDNA + numOfConversion * PACKED_BUFFER_SIZE / CHAR_PER_WORD, charMap, 4, numCharInBuffer);
	}

	numAmbiguity++;
	hsp->numOfAmbiguity = numAmbiguity;
	hsp->ambiguity = MMPoolDispatch(mmPool, (hsp->numOfAmbiguity + 2) * sizeof(Ambiguity));
	if (hsp->numOfAmbiguity > 0) {
		memcpy(hsp->ambiguity + 1, tempAmbiguity, hsp->numOfAmbiguity * sizeof(Ambiguity));
	}
	hsp->ambiguity[0].startPos = 0;
	hsp->ambiguity[0].rightOfEndPos = 0;
	hsp->ambiguity[0].symbol = 0;
	hsp->ambiguity[hsp->numOfAmbiguity+1].startPos = UINT32_MAX;
	hsp->ambiguity[hsp->numOfAmbiguity+1].rightOfEndPos = UINT32_MAX;
	hsp->ambiguity[hsp->numOfAmbiguity+1].symbol = 0;

	hsp->seqOffset[0].firstAmbiguityIndex = 1;
	hsp->seqOffset[0].lastAmbiguityIndex = numAmbiguity;

	MMUnitFree(tempAmbiguity, sizeof(Ambiguity) * ambiguityAllocated);

	return hsp;

}

void HSPFree(MMPool *mmPool, HSP *hsp) {

	if (hsp->packedDNA != NULL) {
		DNAFreePacked(hsp->packedDNA, hsp->dnaLength);
	}
	MMUnitFree(hsp->seqOffset, (hsp->numOfSeq+1) * sizeof(SeqOffset));
	MMUnitFree(hsp->annotation, (hsp->numOfSeq+1) * sizeof(Annotation));
	MMUnitFree(hsp->ambiguity, (hsp->numOfAmbiguity+2) * sizeof(Ambiguity));

	MMPoolReturn(mmPool, hsp, sizeof(hsp));

}

unsigned int HSPParseFASTAToPacked(const char* FASTAFileName, const char* annotationFileName, const char* packedDNAFileName, const char* ambiguityFileName,
					  const unsigned int FASTARandomSeed, const int maskLowerCase) {

	FILE *FASTAFile, *annotationFile, *packedDNAFile, *ambiguityFile;
	Annotation *annotation;
	SeqOffset *seqOffset;
	Ambiguity *ambiguity;
	int annotationAllocated = 256;
	int ambiguityAllocated =  256;
	
	char c;
	int i;
	int numAmbiguity;
	unsigned int lastAmbiguityPos;
	unsigned int numChar, numCharInBuffer, totalNumChar;
	unsigned int sequenceRandomSeed;
	int numSeq;
	unsigned char charMap[255];

	unsigned char buffer[PACKED_BUFFER_SIZE];
	unsigned char packedBuffer[PACKED_BUFFER_SIZE / 4];

	FASTAFile = (FILE*)fopen64(FASTAFileName, "r");
	if (FASTAFile == NULL) {
		fprintf(stderr, "ParseFASTToPacked() : Cannot open FASTAFileName!\n");
		exit(1);
	}

	annotationFile = (FILE*)fopen64(annotationFileName, "w");
	if (annotationFile == NULL) {
		fprintf(stderr, "ParseFASTToPacked() : Cannot open annotationFileName!\n");
		exit(1);
	}

	packedDNAFile = (FILE*)fopen64(packedDNAFileName, "wb");
	if (packedDNAFile == NULL) {
		fprintf(stderr, "ParseFASTToPacked() : Cannot open packedDNAFileName!\n");
		exit(1);
	}

	ambiguityFile = (FILE*)fopen64(ambiguityFileName, "w");
	if (ambiguityFile == NULL) {
		fprintf(stderr, "ParseFASTToPacked() : Cannot open ambiguityFileName!\n");
		exit(1);
	}

	HSPFillCharMap(charMap);

	c = (char)getc(FASTAFile);
	if (c != '>') {
		fprintf(stderr, "ParseFASTToPacked() : FASTA file does not begin with '>'!\n");
		exit(1);
	}

	totalNumChar = 0;
	numSeq = 0;
	numAmbiguity = -1;
	numCharInBuffer = 0;

	annotation = MMUnitAllocate(sizeof(Annotation) * annotationAllocated);
	seqOffset = MMUnitAllocate(sizeof(SeqOffset) * annotationAllocated);
	ambiguity = MMUnitAllocate(sizeof(Ambiguity) * ambiguityAllocated);

	while (!feof(FASTAFile)) {

		numChar = 0;
		if (numSeq >= annotationAllocated) {
			annotation = MMUnitReallocate(annotation, sizeof(Annotation) * annotationAllocated * 2, sizeof(Annotation) * annotationAllocated);
			seqOffset = MMUnitReallocate(seqOffset, sizeof(SeqOffset) * annotationAllocated * 2, sizeof(SeqOffset) * annotationAllocated);
			annotationAllocated *= 2;
		}

		annotation[numSeq].gi = 0;

		c = (char)getc(FASTAFile);
		while (!feof(FASTAFile) && c != '\n') {
			if (numChar < MAX_SEQ_NAME_LENGTH) {
				annotation[numSeq].text[numChar] = c;
				numChar++;
			}
			if (numChar == 3 && annotation[numSeq].text[0] == 'g' 
							 && annotation[numSeq].text[1] == 'i' 
							 && annotation[numSeq].text[2] == '|') {
				fscanf(FASTAFile, "%u|", &(annotation[numSeq].gi));
				numChar = 0;
			}
			c = (char)getc(FASTAFile);
		}
		annotation[numSeq].text[numChar] = '\0';

		// Set random seed for the sequence
		sequenceRandomSeed = FASTARandomSeed;
		if (annotation[numSeq].gi > 0) {
			sequenceRandomSeed += annotation[numSeq].gi;
		} else {
			for (i=0; i<(int)numChar; i++) {
				c = annotation[numSeq].text[i];
				sequenceRandomSeed += c;
			}
		}
		r250_init(sequenceRandomSeed);

		seqOffset[numSeq].startPos = totalNumChar;
		numChar = 0;
		lastAmbiguityPos = (unsigned int)-2;

		c = (char)getc(FASTAFile);
		while (!feof(FASTAFile) && c != '>') {
			// Get sequence
			if (c != '\n' && c != '\t') {
				if (maskLowerCase && c >= 'a' && c <= 'z') {
					c = dnaChar[lowercaseDnaCharIndex];
				}
				if (ambiguityCount[charMap[c]] == 1) {
					buffer[numCharInBuffer] = c;
				} else {
					c = charMap[c];
					if (ambiguityCount[c] > 0) {
						buffer[numCharInBuffer] = dnaChar[ambiguityMatch[c][r250() % ambiguityCount[c]]];
					} else {
						buffer[numCharInBuffer] = dnaChar[r250() % ALPHABET_SIZE];
					}
					if (totalNumChar + numChar == lastAmbiguityPos + 1 && c == (char)ambiguity[numAmbiguity].symbol) {
						ambiguity[numAmbiguity].rightOfEndPos++;
					} else {
						numAmbiguity++;
						if (numAmbiguity >= ambiguityAllocated) {
							ambiguity = MMUnitReallocate(ambiguity, sizeof(Ambiguity) * ambiguityAllocated * 2, sizeof(Ambiguity) * ambiguityAllocated);
							ambiguityAllocated *= 2;
						}
						ambiguity[numAmbiguity].startPos = totalNumChar + numChar;
						ambiguity[numAmbiguity].rightOfEndPos = totalNumChar + numChar;
						ambiguity[numAmbiguity].symbol = c;
					}
					lastAmbiguityPos = totalNumChar + numChar;
				}
				numCharInBuffer++;
				if (numCharInBuffer >= PACKED_BUFFER_SIZE) {
					ConvertTextToBytePacked(buffer, packedBuffer, charMap, 4, PACKED_BUFFER_SIZE);
					fwrite(packedBuffer, 1, PACKED_BUFFER_SIZE / 4, packedDNAFile);
					numCharInBuffer = 0;
				}
				numChar++;
			}
			c = (char)getc(FASTAFile);
		}

		seqOffset[numSeq].endPos = totalNumChar + numChar - 1;
		seqOffset[numSeq].firstAmbiguityIndex = numAmbiguity + 1;
		totalNumChar += numChar;
		numSeq++;

	}
	numAmbiguity++;

	// Finish reading FASTA file
	fclose(FASTAFile);

	// Finalize packed DNA file
	if (numCharInBuffer > 0) {
		ConvertTextToBytePacked(buffer, packedBuffer, charMap, 4, numCharInBuffer);
		fwrite(packedBuffer, 1, (numCharInBuffer + 3) / 4, packedDNAFile);
		numCharInBuffer = 0;
	}
	if (totalNumChar % 4 == 0) {
		c = 0;
		fwrite(&c, 1, 1, packedDNAFile);
	}
	c = (char)(totalNumChar % 4);
	fwrite(&c, 1, 1, packedDNAFile);

	fclose(packedDNAFile);

	// Output annotation file
	fprintf(annotationFile, "%u %u %u\n", totalNumChar, numSeq, FASTARandomSeed);
	for (i=0; i<numSeq; i++) {
		fprintf(annotationFile, "%u %s\n", annotation[i].gi, annotation[i].text);
		fprintf(annotationFile, "%u %u ", seqOffset[i].startPos, seqOffset[i].endPos - seqOffset[i].startPos + 1);
		// output number of ambiguity
		if (i > 0) {
			fprintf(annotationFile, "%u\n", seqOffset[i].firstAmbiguityIndex - seqOffset[i-1].firstAmbiguityIndex);
		} else {
			fprintf(annotationFile, "%u\n", seqOffset[i].firstAmbiguityIndex);
		}
	}
	fclose(annotationFile);

	MMUnitFree(annotation, sizeof(Annotation) * annotationAllocated);
	MMUnitFree(seqOffset, sizeof(SeqOffset) * annotationAllocated);

	// Output ambiguity file
	fprintf(ambiguityFile, "%u %u %u\n", totalNumChar, numSeq, numAmbiguity);
	for (i=0; i<numAmbiguity; i++) {
		fprintf(ambiguityFile, "%u %u %c\n", ambiguity[i].startPos, ambiguity[i].rightOfEndPos - ambiguity[i].startPos + 1, dnaChar[ambiguity[i].symbol]);
	}
	fclose(ambiguityFile);

	MMUnitFree(ambiguity, ambiguityAllocated * sizeof(Ambiguity));

	return numSeq;

}

unsigned int HSPPackedToFASTA(const char* FASTAFileName, const char* annotationFileName, const char* packedDNAFileName, const char* ambiguityFileName) {

	HSP *hsp;
	FILE *FASTAFile;

	hsp = HSPLoad(NULL, packedDNAFileName, annotationFileName, ambiguityFileName);

	// Generate FASTA from packed
	FASTAFile = (FILE*)fopen64(FASTAFileName, "w");
	if (FASTAFile == NULL) {
		fprintf(stderr, "Cannot open FASTA file!\n");
		exit(1);
	}

	// to be done...
	fprintf(stderr, "HSPPackedToFASTA(): Function not complete!\n");

	fclose(FASTAFile);

	HSPFree(NULL, hsp);

	return 0;

}


unsigned int HSPUngappedExtension(const unsigned int *packedDNA, const unsigned int *packedKey, const unsigned int *packedMask, const unsigned int queryPatternLength, 
						 const unsigned int subPattenLength, const unsigned int textLength,
						 HitListWithPosQuery* __restrict hitList, const unsigned int numberOfHit,
						 const HSPUngappedExtLookupTable *ungappedLookupTable,
						 const int cutoffScore, const int dropoffScore) {

	static const unsigned int oddBitMask = 0x55555555;
	static const unsigned int evenBitMask = 0xAAAAAAAA;

	unsigned int hitProcessed, numberOfUngappedHit;
	unsigned int textShift, queryShift;
	unsigned int posText, posQuery;
	unsigned int queryWordPos, textWordPos;
	unsigned int matchVector32, matchVector16;
	int numberOfValidChar;
	unsigned int currentPos;
	int currentScore, forwardMaxScore, backwardMaxScore;
	unsigned int forwardMaxScorePos, backwardMaxScorePos;
	
	hitProcessed = 0;
	numberOfUngappedHit = 0;

	while (hitProcessed < numberOfHit) {

		// Move position to mid-point of hit
		posText = hitList[hitProcessed].posText + subPattenLength / 2;
		posQuery = hitList[hitProcessed].posQuery + subPattenLength / 2;

		textShift = (posText % CHAR_PER_WORD) * BIT_PER_CHAR;
		queryShift = (posQuery % CHAR_PER_WORD) * BIT_PER_CHAR;

		// Forward extend

		textWordPos = posText / CHAR_PER_WORD;
		queryWordPos = posQuery / CHAR_PER_WORD;
		numberOfValidChar = minX(queryPatternLength - posQuery, textLength - posText);

		forwardMaxScore = 0;
		forwardMaxScorePos = posQuery;
		currentPos = posQuery;
		currentScore = 0;

		while (forwardMaxScore <= currentScore + dropoffScore && numberOfValidChar > 0) {

			matchVector32 = ((packedKey[queryWordPos] << queryShift) | ((packedKey[queryWordPos+1] >> (BITS_IN_WORD - queryShift)) * (queryShift > 0)))
							^ ((packedDNA[textWordPos] << textShift) | ((packedDNA[textWordPos+1] >> (BITS_IN_WORD - textShift)) * (textShift > 0)));

			matchVector32 |= (packedMask[queryWordPos] << queryShift) | ((packedMask[queryWordPos+1] >> (BITS_IN_WORD - queryShift)) * (queryShift > 0));

			if (numberOfValidChar < CHAR_PER_WORD) {
				// mask the invalid character as mismatch
				matchVector32 |= ALL_ONE_MASK >> (numberOfValidChar * BIT_PER_CHAR);
			}

			// Process vector to even bits for forward extend
			matchVector32 = (matchVector32 | (matchVector32 << 1)) & evenBitMask;
			matchVector16 = (ungappedLookupTable[matchVector32 >> 16].matchMismatchBitVector << 8) |
							 ungappedLookupTable[matchVector32 & 0xFFFF].matchMismatchBitVector;

			if (currentScore + ungappedLookupTable[matchVector16].maxScore > forwardMaxScore) {
				forwardMaxScore = currentScore + ungappedLookupTable[matchVector16].maxScore;
				forwardMaxScorePos = currentPos + ungappedLookupTable[matchVector16].maxScorePos;
			}
			currentScore += ungappedLookupTable[matchVector16].finalScore;
			currentPos += CHAR_PER_WORD;

			textWordPos++;
			queryWordPos++;
			numberOfValidChar -= CHAR_PER_WORD;

		}


		// Backward extend by word
		textWordPos = posText / CHAR_PER_WORD;
		queryWordPos = posQuery / CHAR_PER_WORD;
		numberOfValidChar = minX(posQuery, posText);

		backwardMaxScore = 0;
		backwardMaxScorePos = posQuery;
		currentPos = posQuery;
		currentScore = 0;

		while (backwardMaxScore <= currentScore + dropoffScore && numberOfValidChar > 0) {

			if (numberOfValidChar >= CHAR_PER_WORD) {
				matchVector32 = ((packedKey[queryWordPos-1] << queryShift) | ((packedKey[queryWordPos] >> (BITS_IN_WORD - queryShift)) * (queryShift > 0)))
								^ ((packedDNA[textWordPos-1] << textShift) | ((packedDNA[textWordPos] >> (BITS_IN_WORD - textShift)) * (textShift > 0)));
				matchVector32 |= (packedMask[queryWordPos-1] << queryShift) | ((packedMask[queryWordPos] >> (BITS_IN_WORD - queryShift)) * (queryShift > 0));
			} else {
				if (textWordPos > 0) {
					matchVector32 = ((packedKey[queryWordPos] >> (BITS_IN_WORD - queryShift)) * (queryShift > 0))
									^ ((packedDNA[textWordPos-1] << textShift) | ((packedDNA[textWordPos] >> (BITS_IN_WORD - textShift)) * (textShift > 0)));
					matchVector32 |= (packedMask[queryWordPos] >> (BITS_IN_WORD - queryShift)) * (queryShift > 0);
				} else if (queryWordPos > 0) {
					matchVector32 = ((packedKey[queryWordPos-1] << queryShift) | ((packedKey[queryWordPos] >> (BITS_IN_WORD - queryShift)) * (queryShift > 0)))
									^ ((packedDNA[textWordPos] >> (BITS_IN_WORD - textShift)) * (textShift > 0));
					matchVector32 |= (packedMask[queryWordPos-1] << queryShift) | ((packedMask[queryWordPos] >> (BITS_IN_WORD - queryShift)) * (queryShift > 0));
				} else {
					matchVector32 = ((packedKey[queryWordPos] >> (BITS_IN_WORD - queryShift)) * (queryShift > 0))
									^ ((packedDNA[textWordPos] >> (BITS_IN_WORD - textShift)) * (textShift > 0));
					matchVector32 |= (packedMask[queryWordPos] >> (BITS_IN_WORD - queryShift)) * (queryShift > 0);
				}
				matchVector32 |= ALL_ONE_MASK << (numberOfValidChar * BIT_PER_CHAR);
			}

			// Process vector to odd bits for backward extend
			matchVector32 = (matchVector32 | (matchVector32 >> 1)) & oddBitMask;
			matchVector16 = ungappedLookupTable[matchVector32 >> 16].matchMismatchBitVector |
							(ungappedLookupTable[matchVector32 & 0xFFFF].matchMismatchBitVector << 8);	// matchVecter16 is reverse of packed odd bits

			if (currentScore + ungappedLookupTable[matchVector16].maxScore > backwardMaxScore) {
				backwardMaxScore = currentScore + ungappedLookupTable[matchVector16].maxScore;
				backwardMaxScorePos = currentPos - ungappedLookupTable[matchVector16].maxScorePos;
			}
			currentScore += ungappedLookupTable[matchVector16].finalScore;
			currentPos -= CHAR_PER_WORD;

			textWordPos--;
			queryWordPos--;
			numberOfValidChar -= CHAR_PER_WORD;

		}

		// Check cut off score
		if (backwardMaxScore + forwardMaxScore >= cutoffScore) {
			hitList[numberOfUngappedHit].posQuery = (forwardMaxScorePos + backwardMaxScorePos) / 2;
			hitList[numberOfUngappedHit].posText = posText - posQuery + hitList[numberOfUngappedHit].posQuery;
			numberOfUngappedHit++;
		}

		hitProcessed++;

	}

	return numberOfUngappedHit;

}

int HSPGappedExtension(const unsigned int *packedDNA, const unsigned int textLength, const unsigned char *convertedKey, const int queryPatternLength, 
					   const HitListWithPosQuery* hitList, const int numberOfHit, 
					   GappedHitList* __restrict gappedHitList,
					   MMPool *mmPool,
					   const int matchScore, const int mismatchScore,
					   const int gapOpenScore, const int gapExtendScore,
					   const int cutoffScore, const int dropoffScore) {

	int* __restrict B;
	int M, BLast;
	int* __restrict I;	// insert and delete wrt query
	int D;				// insert and delete wrt query
	int scoringMatrix[16][16];

	int hitProcessed, numberOfGappedHit;

	int i;
	unsigned int c = 0;
	char textChar;
	int queryOffset;
	unsigned int textOffset;
	unsigned int textDecodePos;
	int charToProcess;

	int lastStartPos, startPos, nextStartPos;
	int endPos, nextEndPos;
	unsigned int textPos;

	int maxLengthOfGap;
	int currentPos;
	int forwardMaxScore, backwardMaxScore;
	unsigned int forwardMaxScorePosQuery, forwardMaxScorePosText;
	unsigned int backwardMaxScorePosQuery, backwardMaxScorePosText;

	int dpCellAllocated;

	dpCellAllocated = minX(queryPatternLength + 1, MAX_ALIGNMENT_LENGTH);
	
	// allocate working memory
	B = MMPoolDispatch(mmPool, dpCellAllocated * sizeof(int));
	I = MMPoolDispatch(mmPool, dpCellAllocated * sizeof(int));

	HSPFillScoringMatrix(scoringMatrix, matchScore, mismatchScore, 0);

	hitProcessed = 0;
	numberOfGappedHit = 0;

	maxLengthOfGap = (dropoffScore + gapOpenScore) / (-gapExtendScore);

	while (hitProcessed < numberOfHit) {

		queryOffset = hitList[hitProcessed].posQuery;
		textOffset = hitList[hitProcessed].posText;

		// Forward extend

		textDecodePos = textOffset;
		if (textDecodePos % CHAR_PER_WORD > 0) {
			// decode text to the next word boundary
			charToProcess = CHAR_PER_WORD - (textDecodePos % CHAR_PER_WORD);
			c = packedDNA[textDecodePos / CHAR_PER_WORD] << (BITS_IN_WORD - charToProcess * BIT_PER_CHAR);
			textDecodePos += charToProcess;
		}

		// Fill initial scores
		B[0] = 0;
		I[0] = gapExtendScore + gapOpenScore;
		for (i=1; i<=maxLengthOfGap; i++) {
			B[i] = i * gapExtendScore + gapOpenScore;
			I[i] = B[i] + gapExtendScore + gapOpenScore;
		}
		I[i] = DP_NEG_INFINITY;

		lastStartPos = startPos = 0;
		endPos = minX(maxLengthOfGap + 1, queryPatternLength - queryOffset);
		textPos = textOffset;
		 
		forwardMaxScore = 0;
		forwardMaxScorePosQuery = queryOffset;
		forwardMaxScorePosText = textOffset;

		while (startPos <= endPos && textPos < textLength) {

			if (endPos >= dpCellAllocated) {
				fprintf(stderr, "HSPGappedExtension(): Not enough DP cells allocated!\n");
				exit(1);
			}

			if (textPos >= textDecodePos) {
				c = packedDNA[textDecodePos / CHAR_PER_WORD];
				textDecodePos += CHAR_PER_WORD;
			}
			textChar = (char)(c >> (BITS_IN_WORD - BIT_PER_CHAR));
			c <<= BIT_PER_CHAR;

			nextEndPos = 0;
			nextStartPos = INT_MAX;

			if (startPos > lastStartPos) {
				BLast = B[startPos - 1];
				D = DP_NEG_INFINITY;
				currentPos = startPos;
			} else {
				// all scores in currentPos - 1 is DP_NEG_INFINITY
				BLast = B[startPos];
				B[startPos] = I[startPos];
				I[startPos] = I[startPos] + gapExtendScore;
				D = B[startPos] + gapExtendScore + gapOpenScore;
				if (B[startPos] + dropoffScore >= forwardMaxScore) {
					nextStartPos = startPos;
					nextEndPos = startPos;
				}
				currentPos = startPos + 1;
			}

			while (currentPos <= endPos) {

				M = BLast + scoringMatrix[textChar][convertedKey[currentPos + queryOffset - 1]];
				BLast = B[currentPos];

				if (M >= I[currentPos] && M >= D) {
					// matchScore is maximum
					B[currentPos] = M;
					I[currentPos] = maxX(M + gapOpenScore, I[currentPos]) + gapExtendScore;
					D = maxX(M + gapOpenScore, D) + gapExtendScore;
				} else {
					B[currentPos] = maxX(I[currentPos], D);
					// insert cannot follow delete and delete cannot follow insert
					I[currentPos] = I[currentPos] + gapExtendScore;
					D = D + gapExtendScore;
				}
				if (B[currentPos] + dropoffScore >= forwardMaxScore) {
					if (nextStartPos > currentPos) {
						nextStartPos = currentPos;
					}
					nextEndPos = currentPos;
				}
				if (B[currentPos] > forwardMaxScore) {
					forwardMaxScore = B[currentPos];
					forwardMaxScorePosQuery = currentPos - 1 + queryOffset;
					forwardMaxScorePosText = textPos;
				}

				currentPos++;

			}

			nextEndPos++;
			if (nextEndPos == currentPos) {
				while (nextEndPos <= queryPatternLength - queryOffset && D + dropoffScore >= forwardMaxScore)  {
					if (nextEndPos >= dpCellAllocated) {
						fprintf(stderr, "HSPGappedExtension(): Not enough DP cells allocated!\n");
						exit(1);
					}
					B[nextEndPos] = D;
					D = D + gapExtendScore;
					I[nextEndPos] = DP_NEG_INFINITY;
					nextEndPos++;
				}
				I[nextEndPos] = DP_NEG_INFINITY;
			}

			if (nextEndPos > queryPatternLength - queryOffset) {
				nextEndPos = queryPatternLength - queryOffset;
			}

			lastStartPos = startPos;
			startPos = nextStartPos;
			endPos = nextEndPos;

			textPos++;

		}

		// Backward extend

		textDecodePos = textOffset;
		if (textDecodePos % CHAR_PER_WORD > 0) {
			charToProcess = textDecodePos % CHAR_PER_WORD;
			c = packedDNA[textDecodePos / CHAR_PER_WORD] >> (BITS_IN_WORD - charToProcess * BIT_PER_CHAR);
			textDecodePos -= charToProcess;
		}

		// Fill initial scores
		B[0] = 0;
		I[0] = gapExtendScore + gapOpenScore;
		for (i=1; i<=maxLengthOfGap; i++) {
			B[i] = i * gapExtendScore + gapOpenScore;
			I[i] = B[i] + gapExtendScore + gapOpenScore;
		}
		I[i] = DP_NEG_INFINITY;

		lastStartPos = startPos = 0;
		endPos = minX(maxLengthOfGap + 1, queryOffset);
		textPos = textOffset;
		 
		backwardMaxScore = 0;
		backwardMaxScorePosQuery = queryOffset;
		backwardMaxScorePosText = textOffset;

		while (startPos <= endPos && textPos > 0) {

			if (endPos >= dpCellAllocated) {
				fprintf(stderr, "HSPGappedExtension(): Not enough DP cells allocated!\n");
				exit(1);
			}

			if (textPos <= textDecodePos) {
				c = packedDNA[textDecodePos / CHAR_PER_WORD - 1];
				textDecodePos -= CHAR_PER_WORD;
			}
			textChar = (char)(c & 0x3);
			c >>= BIT_PER_CHAR;

			nextEndPos = 0;
			nextStartPos = INT_MAX;

			if (startPos > lastStartPos) {
				BLast = B[startPos - 1];
				D = DP_NEG_INFINITY;
				currentPos = startPos;
			} else {
				// all scores in currentPos - 1 is DP_NEG_INFINITY
				BLast = B[startPos];
				B[startPos] = I[startPos];
				I[startPos] = I[startPos] + gapExtendScore;
				D = B[startPos] + gapExtendScore + gapOpenScore;
				if (B[startPos] + dropoffScore >= backwardMaxScore) {
					nextStartPos = startPos;
					nextEndPos = startPos;
				}
				currentPos = startPos + 1;
			}

			while (currentPos <= endPos) {

				M = BLast + scoringMatrix[textChar][convertedKey[queryOffset - currentPos]];
				BLast = B[currentPos];

				if (M >= I[currentPos] && M >= D) {
					// matchScore is maximum
					B[currentPos] = M;
					I[currentPos] = maxX(M + gapOpenScore, I[currentPos]) + gapExtendScore;
					D = maxX(M + gapOpenScore, D) + gapExtendScore;
				} else {
					B[currentPos] = maxX(I[currentPos], D);
					// insert cannot follow delete and delete cannot follow insert
					I[currentPos] = I[currentPos] + gapExtendScore;
					D = D + gapExtendScore;
				}
				if (B[currentPos] + dropoffScore >= backwardMaxScore) {
					if (nextStartPos > currentPos) {
						nextStartPos = currentPos;
					}
					nextEndPos = currentPos;
				}
				if (B[currentPos] > backwardMaxScore) {
					backwardMaxScore = B[currentPos];
					backwardMaxScorePosQuery =  queryOffset - currentPos;
					backwardMaxScorePosText = textPos - 1;
				}

				currentPos++;

			}

			nextEndPos++;
			if (nextEndPos == currentPos) {
				while (nextEndPos <= queryOffset && D + dropoffScore >= backwardMaxScore)  {
					if (nextEndPos >= dpCellAllocated) {
						fprintf(stderr, "HSPGappedExtension(): Not enough DP cells allocated!\n");
						exit(1);
					}
					B[nextEndPos] = D;
					D = D + gapExtendScore;
					I[nextEndPos] = DP_NEG_INFINITY;
					nextEndPos++;
				}
				I[nextEndPos] = DP_NEG_INFINITY;
			}

			if (nextEndPos > queryOffset) {
				nextEndPos = queryOffset;
			}

			lastStartPos = startPos;
			startPos = nextStartPos;
			endPos = nextEndPos;

			textPos--;

		}

		// check cutoff score
		if (forwardMaxScore + backwardMaxScore >= cutoffScore) {
			gappedHitList[numberOfGappedHit].posText = backwardMaxScorePosText;
			gappedHitList[numberOfGappedHit].posQuery = backwardMaxScorePosQuery;
			gappedHitList[numberOfGappedHit].lengthQuery = forwardMaxScorePosQuery + 1 - backwardMaxScorePosQuery;
			gappedHitList[numberOfGappedHit].score = forwardMaxScore + backwardMaxScore;
			gappedHitList[numberOfGappedHit].lengthText = forwardMaxScorePosText + 1 - backwardMaxScorePosText;
			gappedHitList[numberOfGappedHit].ungappedPosText = textOffset;
			gappedHitList[numberOfGappedHit].ungappedPosQuery = queryOffset;
			numberOfGappedHit++;
		}

		hitProcessed++;

	}

	// free working memory
	MMPoolReturn(mmPool, B, dpCellAllocated * sizeof(int));
	MMPoolReturn(mmPool, I, dpCellAllocated * sizeof(int));


	return numberOfGappedHit;

}

// mmPool can be NULL but alignmentPool must be allocated for storing alignments as they are not freed explicitly and must be freed through MMPool
// gappedHitList should be sorted in increasin posText order
int HSPGappedExtensionWithTraceback(const unsigned int *packedDNA, const unsigned char *convertedKey, const int queryPatternLength, 
					   GappedHitList* __restrict gappedHitList, const int numberOfHit, 
					   MMPool *mmPool, MMPool *alignmentPool,
					   const SeqOffset *seqOffset, const Ambiguity *ambiguity,
					   const int matchScore, const int mismatchScore,
					   const int gapOpenScore, const int gapExtendScore,
					   const double maxEvalue, const int dropoffScore) {

	char* __restrict textBuffer;
	int M, BLast;
	int* __restrict B;
	int* __restrict I;	// insert and delete wrt query
	int D;				// insert and delete wrt query
	char* __restrict IType;	// Gapped opening or gap extension
	char DType;				// Gapped opening or gap extension
	int scoringMatrix[16][16];

	char* __restrict traceback;
	int* __restrict tracebackIndex;
	char* __restrict tempForwardAlignment;
	char* __restrict tempBackwardAlignment;
	char* __restrict tempForwardAuxiliaryText;
	char* __restrict tempBackwardAuxiliaryText;
	unsigned int* __restrict alignment;
	unsigned int* __restrict auxiliaryText;

	int hitProcessed, numberOfGappedHit;

	int i, j;
	unsigned int c;
	int queryOffset;
	unsigned int textOffset;
	unsigned int sequenceStart, sequenceEnd;
	int charToProcess;
	unsigned int textDecodePos;
	double evalue;

	int lastStartPos, startPos, nextStartPos;
	int endPos, nextEndPos;
	unsigned int textPos;
	int ambiguityIndex;

	int maxLengthOfGap;
	int currentPos;
	int currentTracebackIndex;
	int forwardMaxScore, backwardMaxScore;
	unsigned int forwardMaxScorePosQuery, forwardMaxScorePosText;
	unsigned int backwardMaxScorePosQuery, backwardMaxScorePosText;
	int forwardNumOfAlignment, forwardNumOfAuxiliaryText;
	int backwardNumOfAlignment, backwardNumOfAuxiliaryText;
	int numOfAlignmentWord, numOfAuxiliaryTextWord;
	int tracebackAllocationUnit, tracebackAllocated;
	char dpType;
	int dpScore;

	int dpCellAllocated;

	dpCellAllocated = minX(queryPatternLength + 1, MAX_ALIGNMENT_LENGTH);

	if (alignmentPool == NULL) {
		fprintf(stderr, "HSPGappedExtensionWithTraceback(): alignmentPool is not allocated!\n");
		exit(1);
	}

	maxLengthOfGap = (dropoffScore + gapOpenScore) / (-gapExtendScore);

	// allocate working memory
	textBuffer = MMPoolDispatch(mmPool, dpCellAllocated * 2);
	B = MMPoolDispatch(mmPool, dpCellAllocated * sizeof(int));
	IType = MMPoolDispatch(mmPool, dpCellAllocated);
	I = MMPoolDispatch(mmPool, dpCellAllocated * sizeof(int));

	tracebackIndex = MMPoolDispatch(mmPool, dpCellAllocated * 2 * 2 * sizeof(int));

	tracebackAllocated = tracebackAllocationUnit = dpCellAllocated * maxLengthOfGap * maxLengthOfGap;
	traceback = MMUnitAllocate(tracebackAllocated);

	tempForwardAlignment = MMPoolDispatch(mmPool, dpCellAllocated * 2);
	tempBackwardAlignment = MMPoolDispatch(mmPool, dpCellAllocated * 2);
	tempForwardAuxiliaryText = MMPoolDispatch(mmPool, dpCellAllocated);
	tempBackwardAuxiliaryText = MMPoolDispatch(mmPool, dpCellAllocated);

	HSPFillScoringMatrix(scoringMatrix, matchScore, mismatchScore, 0);

	hitProcessed = 0;
	numberOfGappedHit = 0;
	ambiguityIndex = 0;

	while (hitProcessed < numberOfHit) {

		sequenceStart = seqOffset[gappedHitList[hitProcessed].dbSeqIndex].startPos;
		sequenceEnd = seqOffset[gappedHitList[hitProcessed].dbSeqIndex].endPos;
		if (seqOffset[gappedHitList[hitProcessed].dbSeqIndex].firstAmbiguityIndex > ambiguityIndex) {
			ambiguityIndex = seqOffset[gappedHitList[hitProcessed].dbSeqIndex].firstAmbiguityIndex;
		}

		textOffset = gappedHitList[hitProcessed].ungappedPosText;
		queryOffset = gappedHitList[hitProcessed].ungappedPosQuery;


		// Forward extend

		textDecodePos = textOffset;
		if (textDecodePos % CHAR_PER_WORD > 0) {
			// decode text to the next word boundary
			charToProcess = CHAR_PER_WORD - (textDecodePos % CHAR_PER_WORD);
			c = packedDNA[textDecodePos / CHAR_PER_WORD] << (BITS_IN_WORD - charToProcess * BIT_PER_CHAR);
			for (i=0; i<charToProcess; i++) {
				textBuffer[textDecodePos + i - textOffset + 1] = (unsigned char)(c >> (BITS_IN_WORD - BIT_PER_CHAR));
				c <<= BIT_PER_CHAR;
			}

			// Apply ambiguity
			while (ambiguity[ambiguityIndex].rightOfEndPos <= textOffset) {
				ambiguityIndex++;
			}
			if (ambiguity[ambiguityIndex].startPos < textOffset + charToProcess) {
				for (i=0; i<charToProcess; i++) {
					while (ambiguity[ambiguityIndex].rightOfEndPos <= textOffset + i) {
						ambiguity++;
					}
					if (ambiguity[ambiguityIndex].startPos <= textOffset + i) {
						textBuffer[textDecodePos + i - textOffset + 1] = (char)ambiguity[ambiguityIndex].symbol;
					}
				}
			}
				
			textDecodePos += charToProcess;
		}

		// Fill initial scores
		B[0] = 0;
		traceback[0] = DP_MATCH_MISMATCH;
		I[0] = gapExtendScore + gapOpenScore;
		IType[0] = DP_INSERT_OPEN;

		B[1] = gapExtendScore + gapOpenScore;
		traceback[1] = DP_DELETE | DP_DELETE_OPEN;
		I[1] = B[1] + gapExtendScore + gapOpenScore;
		IType[1] = DP_INSERT_OPEN;

		for (i=2; i<=maxLengthOfGap; i++) {
			B[i] = B[i-1] + gapExtendScore;
			traceback[i] = DP_DELETE;
			I[i] = B[i] + gapExtendScore + gapOpenScore;
			IType[i] = DP_INSERT_OPEN;
		}
		I[i] = DP_NEG_INFINITY;

		currentTracebackIndex = maxLengthOfGap + 1;

		lastStartPos = startPos = 0;
		endPos = minX((maxLengthOfGap + 1), queryPatternLength - queryOffset);
		textPos = textOffset;
		 
		forwardMaxScore = 0;
		forwardMaxScorePosQuery = queryOffset;
		forwardMaxScorePosText = textOffset;

		tracebackIndex[0] = 0;	// The first trackback of the row
		tracebackIndex[1] = 0;	// The first query position of the row

		while (startPos <= endPos && textPos <= sequenceEnd) {

			if (endPos >= dpCellAllocated) {
				fprintf(stderr, "HSPGappedExtensionWithTraceback(): Not enough DP cells allocated!\n");
				exit(1);
			}

			if (textPos < textDecodePos) {
			} else {
				// decode a word of text
				c = packedDNA[textDecodePos / CHAR_PER_WORD];
				for (j=0; j<CHAR_PER_WORD; j++) {
					textBuffer[textDecodePos + j - textOffset + 1] = (unsigned char)(c >> (BITS_IN_WORD - BIT_PER_CHAR));
					c <<= BIT_PER_CHAR;
				}

				// Apply ambiguity
				while (ambiguity[ambiguityIndex].rightOfEndPos <= textPos) {
					ambiguityIndex++;
				}
				if (ambiguity[ambiguityIndex].startPos < textPos + CHAR_PER_WORD) {
					for (j=0; j<CHAR_PER_WORD; j++) {
						while (ambiguity[ambiguityIndex].rightOfEndPos <= textPos + j) {
							ambiguity++;
						}
						if (ambiguity[ambiguityIndex].startPos <= textPos + j) {
							textBuffer[textDecodePos + j - textOffset + 1] = (char)ambiguity[ambiguityIndex].symbol;
						}
					}
				}

				textDecodePos += CHAR_PER_WORD;
			}

			// traceback
			tracebackIndex[(textPos - textOffset + 1) * 2] = currentTracebackIndex;	// The first trackback of the row
			tracebackIndex[(textPos - textOffset + 1) * 2 + 1] = startPos;			// The first query position of the row

			nextEndPos = 0;
			nextStartPos = INT_MAX;

			if (currentTracebackIndex + queryPatternLength >= tracebackAllocated) {
				traceback = MMUnitReallocate(traceback, tracebackAllocated + tracebackAllocationUnit, tracebackAllocated);
				tracebackAllocated += tracebackAllocationUnit;
			}

			if (startPos > lastStartPos) {
				BLast = B[startPos - 1];
				D = DP_NEG_INFINITY;
				DType = 0;
				currentPos = startPos;
			} else {
				// all scores in currentPos - 1 is DP_NEG_INFINITY
				BLast = B[startPos];
				B[startPos] = I[startPos];
				traceback[currentTracebackIndex] = DP_INSERT | IType[startPos];
				I[startPos] = I[startPos] + gapExtendScore;
				IType[startPos] = DP_INSERT_EXTEND;
				D = B[startPos] + gapExtendScore + gapOpenScore;
				DType = DP_DELETE_OPEN;
				if (B[startPos] + dropoffScore >= forwardMaxScore) {
					nextStartPos = startPos;
					nextEndPos = startPos;
				}
				currentPos = startPos + 1;
				currentTracebackIndex++;
			}

			while (currentPos <= endPos) {

				M = BLast + scoringMatrix[textBuffer[textPos - textOffset + 1]][convertedKey[currentPos + queryOffset - 1]];
				BLast = B[currentPos];

				if (M >= I[currentPos] && M >= D) {
					// matchScore is maximum
					B[currentPos] = M;
					traceback[currentTracebackIndex] = DP_MATCH_MISMATCH | IType[currentPos] | DType;
					if (M + gapOpenScore >= I[currentPos]) {
						I[currentPos] = M + gapOpenScore + gapExtendScore;
						IType[currentPos] = DP_INSERT_OPEN;
					} else {
						I[currentPos] = I[currentPos] + gapExtendScore;
						IType[currentPos] = DP_INSERT_EXTEND;
					}
					if (M + gapOpenScore >= D) {
						D = M + gapOpenScore + gapExtendScore;
						DType = DP_DELETE_OPEN;
					} else {
						D = D + gapExtendScore;
						DType = DP_DELETE_EXTEND;
					}
					if (B[currentPos] > forwardMaxScore) {
						forwardMaxScore = B[currentPos];
						forwardMaxScorePosQuery = currentPos + queryOffset;
						forwardMaxScorePosText = textPos + 1;
					}
				} else {
					if (I[currentPos] >= D) {
						B[currentPos] = I[currentPos];
						traceback[currentTracebackIndex] = DP_INSERT | IType[currentPos] | DType;
					} else {
						B[currentPos] = D;
						traceback[currentTracebackIndex] = DP_DELETE | IType[currentPos] | DType;
					}
					// insert cannot follow delete and delete cannot follow insert
					I[currentPos] = I[currentPos] + gapExtendScore;
					IType[currentPos] = DP_INSERT_EXTEND;
					D = D + gapExtendScore;
					DType = DP_DELETE_EXTEND;
				}
				if (B[currentPos] + dropoffScore >= forwardMaxScore) {
					if (nextStartPos > currentPos) {
						nextStartPos = currentPos;
					}
					nextEndPos = currentPos;
				}

				currentPos++;
				currentTracebackIndex++;

			}

			nextEndPos++;
			if (nextEndPos == currentPos) {
				while (nextEndPos <= queryPatternLength - queryOffset && D + dropoffScore >= forwardMaxScore)  {
					if (nextEndPos >= dpCellAllocated) {
						fprintf(stderr, "HSPGappedExtensionWithTraceback(): Not enough DP cells allocated!\n");
						exit(1);
					}
					B[nextEndPos] = D;
					traceback[currentTracebackIndex] = DP_DELETE | DType;
					D = D + gapExtendScore;
					DType = DP_DELETE_EXTEND;
					I[nextEndPos] = DP_NEG_INFINITY;
					IType[nextEndPos] = 0;
					nextEndPos++;
					currentTracebackIndex++;
				}
				I[nextEndPos] = DP_NEG_INFINITY;
				IType[nextEndPos] = 0;
			}

			if (nextEndPos > queryPatternLength - queryOffset) {
				nextEndPos = queryPatternLength - queryOffset;
			}

			lastStartPos = startPos;
			startPos = nextStartPos;
			endPos = nextEndPos;

			textPos++;

		}

		// traceback
		forwardNumOfAlignment = 0;
		forwardNumOfAuxiliaryText = 0;
		textPos = forwardMaxScorePosText - 1;
		currentPos = forwardMaxScorePosQuery - queryOffset;
		dpScore = 0;	// for verifying traceback

		currentTracebackIndex = tracebackIndex[(textPos - textOffset + 1) * 2] 
								+ currentPos - tracebackIndex[(textPos - textOffset + 1) * 2 + 1];

		while (currentTracebackIndex > 0) {

			dpType = traceback[currentTracebackIndex] & DP_MASK;

			if (dpType == DP_MATCH_MISMATCH) {
				dpScore += scoringMatrix[textBuffer[textPos - textOffset + 1]][convertedKey[currentPos + queryOffset - 1]];
				if (textBuffer[textPos - textOffset + 1] == convertedKey[currentPos + queryOffset - 1] &&
					textBuffer[textPos - textOffset + 1] < 4) {	// match and not ambiguity
					tempForwardAlignment[forwardNumOfAlignment] = ALIGN_MATCH;
				} else {
					tempForwardAlignment[forwardNumOfAlignment] = ALIGN_MISMATCH_AMBIGUITY;
					tempForwardAuxiliaryText[forwardNumOfAuxiliaryText] = textBuffer[textPos - textOffset + 1];
					forwardNumOfAuxiliaryText++;
				}
				textPos--;
				currentPos--;
				currentTracebackIndex = tracebackIndex[(textPos - textOffset + 1) * 2] 
										+ currentPos - tracebackIndex[(textPos - textOffset + 1) * 2 + 1];
			} else if (dpType ==  DP_INSERT) {
				while (!(traceback[currentTracebackIndex] & DP_INSERT_OPEN)) {
					dpScore += gapExtendScore;
					tempForwardAlignment[forwardNumOfAlignment] = ALIGN_INSERT;
					forwardNumOfAlignment++;
					tempForwardAuxiliaryText[forwardNumOfAuxiliaryText] = textBuffer[textPos - textOffset + 1];
					forwardNumOfAuxiliaryText++;
					textPos--;
					currentTracebackIndex = tracebackIndex[(textPos - textOffset + 1) * 2] 
											+ currentPos - tracebackIndex[(textPos - textOffset + 1) * 2 + 1];
				}
				dpScore += gapOpenScore + gapExtendScore;
				tempForwardAlignment[forwardNumOfAlignment] = ALIGN_INSERT;
				tempForwardAuxiliaryText[forwardNumOfAuxiliaryText] = textBuffer[textPos - textOffset + 1];
				forwardNumOfAuxiliaryText++;
				textPos--;
				currentTracebackIndex = tracebackIndex[(textPos - textOffset + 1) * 2] 
										+ currentPos - tracebackIndex[(textPos - textOffset + 1) * 2 + 1];
			} else {
				while (!(traceback[currentTracebackIndex] & DP_DELETE_OPEN)) {
					dpScore += gapExtendScore;
					tempForwardAlignment[forwardNumOfAlignment] = ALIGN_DELETE;
					forwardNumOfAlignment++;
					currentPos--;
					currentTracebackIndex--;
				}
				dpScore += gapOpenScore + gapExtendScore;
				tempForwardAlignment[forwardNumOfAlignment] = ALIGN_DELETE;
				currentPos--;
				currentTracebackIndex--;
			}

			forwardNumOfAlignment++;

		}

		if (dpScore != forwardMaxScore) {
			fprintf(stderr, "Forward gapped extension traceback error!\n");
			exit(1);
		}

		// Backward extend

		textDecodePos = textOffset;
		if (textDecodePos % CHAR_PER_WORD > 0) {
			// decode text to the next word boundary
			charToProcess = textDecodePos % CHAR_PER_WORD;
			c = packedDNA[textDecodePos / CHAR_PER_WORD] >> (BITS_IN_WORD - charToProcess * BIT_PER_CHAR);
			for (i=0; i<charToProcess; i++) {
				textBuffer[textOffset - textDecodePos + i + 1] = (unsigned char)(c & 0x3);
				c >>= BIT_PER_CHAR;
			}

			// Apply ambiguity
			while (ambiguity[ambiguityIndex].startPos >= textOffset) {
				ambiguityIndex--;
			}
			if (ambiguity[ambiguityIndex].rightOfEndPos > (textOffset - charToProcess)) {
				for (i=0; i<charToProcess; i++) {
					while (ambiguity[ambiguityIndex].startPos >= (textOffset - i)) {
						ambiguityIndex--;
					}
					if (ambiguity[ambiguityIndex].rightOfEndPos > (textOffset - i - 1)) {
						textBuffer[textOffset - textDecodePos + i + 1] = (char)ambiguity[ambiguityIndex].symbol;
					}
				}
			}

			textDecodePos -= charToProcess;
		}

		// Fill initial scores
		B[0] = 0;
		traceback[0] = DP_MATCH_MISMATCH;
		I[0] = gapExtendScore + gapOpenScore;
		IType[0] = DP_INSERT_OPEN;

		B[1] = gapExtendScore + gapOpenScore;
		traceback[1] = DP_DELETE | DP_DELETE_OPEN;
		I[1] = B[1] + gapExtendScore + gapOpenScore;
		IType[1] = DP_INSERT_OPEN;

		for (i=2; i<=maxLengthOfGap; i++) {
			B[i] = B[i-1] + gapExtendScore;
			traceback[i] = DP_DELETE;
			I[i] = B[i] + gapExtendScore + gapOpenScore;
			IType[i] = DP_INSERT_OPEN;
		}
		I[i] = DP_NEG_INFINITY;

		currentTracebackIndex = maxLengthOfGap + 1;

		lastStartPos = startPos = 0;
		endPos = minX((maxLengthOfGap + 1), queryOffset);
		textPos = textOffset;
		 
		backwardMaxScore = 0;
		backwardMaxScorePosQuery = queryOffset;
		backwardMaxScorePosText = textOffset;

		tracebackIndex[0] = 0;
		tracebackIndex[1] = 0;

		while (startPos <= endPos && textPos + 1 > sequenceStart) {

			if (endPos >= dpCellAllocated) {
				fprintf(stderr, "HSPGappedExtensionWithTraceback(): Not enough DP cells allocated!\n");
				exit(1);
			}

			if (textPos > textDecodePos) {
			} else {
				// decode a word of text
				c = packedDNA[textDecodePos / CHAR_PER_WORD - 1];
				for (j=0; j<CHAR_PER_WORD; j++) {
					textBuffer[textOffset - textDecodePos + j + 1] = (unsigned char)(c & 0x3);
					c >>= BIT_PER_CHAR;
				}

				// Apply ambiguity
				while (ambiguity[ambiguityIndex].startPos >= textPos) {
					ambiguityIndex--;
				}
				if (ambiguity[ambiguityIndex].rightOfEndPos > (textPos - CHAR_PER_WORD)) {
					for (j=0; j<CHAR_PER_WORD; j++) {
						while (ambiguity[ambiguityIndex].startPos >= (textPos - j)) {
							ambiguityIndex--;
						}
						if (ambiguity[ambiguityIndex].rightOfEndPos > (textPos - j - 1)) {
							textBuffer[textOffset - textDecodePos + j + 1] = (char)ambiguity[ambiguityIndex].symbol;
						}
					}
				}

				textDecodePos -= CHAR_PER_WORD;
			}

			// traceback
			tracebackIndex[(textOffset - textPos + 1) * 2] = currentTracebackIndex;
			tracebackIndex[(textOffset - textPos + 1) * 2 + 1] = startPos;

			nextEndPos = 0;
			nextStartPos = INT_MAX;

			if (currentTracebackIndex + queryPatternLength >= tracebackAllocated) {
				traceback = MMUnitReallocate(traceback, tracebackAllocated + tracebackAllocationUnit, tracebackAllocated);
				tracebackAllocated += tracebackAllocationUnit;
			}

			if (startPos > lastStartPos) {
				BLast = B[startPos - 1];
				D = DP_NEG_INFINITY;
				DType = 0;
				currentPos = startPos;
			} else {
				// all scores in currentPos - 1 is DP_NEG_INFINITY
				BLast = B[startPos];
				B[startPos] = I[startPos];
				traceback[currentTracebackIndex] = DP_INSERT | IType[startPos];
				I[startPos] = I[startPos] + gapExtendScore;
				IType[startPos] = DP_INSERT_EXTEND;
				D = B[startPos] + gapExtendScore + gapOpenScore;
				DType = DP_DELETE_OPEN;
				if (B[startPos] + dropoffScore >= forwardMaxScore) {
					nextStartPos = startPos;
					nextEndPos = startPos;
				}
				currentPos = startPos + 1;
				currentTracebackIndex++;
			}


			while (currentPos <= endPos) {

				M = BLast + scoringMatrix[textBuffer[textOffset - textPos + 1]][convertedKey[queryOffset - currentPos]];
				BLast = B[currentPos];

				if (M >= I[currentPos] && M >= D) {
					// matchScore is maximum
					B[currentPos] = M;
					traceback[currentTracebackIndex] = DP_MATCH_MISMATCH | IType[currentPos] | DType;
					if (M + gapOpenScore >= I[currentPos]) {
						I[currentPos] = M + gapOpenScore + gapExtendScore;
						IType[currentPos] = DP_INSERT_OPEN;
					} else {
						I[currentPos] = I[currentPos] + gapExtendScore;
						IType[currentPos] = DP_INSERT_EXTEND;
					}
					if (M + gapOpenScore >= D) {
						D = M + gapOpenScore + gapExtendScore;
						DType = DP_DELETE_OPEN;
					} else {
						D = D + gapExtendScore;
						DType = DP_DELETE_EXTEND;
					}
					if (B[currentPos] > backwardMaxScore) {
						backwardMaxScore = B[currentPos];
						backwardMaxScorePosQuery =  queryOffset - currentPos;
						backwardMaxScorePosText = textPos - 1;
					}
				} else {
					if (I[currentPos] >= D) {
						B[currentPos] = I[currentPos];
						traceback[currentTracebackIndex] = DP_INSERT | IType[currentPos] | DType;
					} else {
						B[currentPos] = D;
						traceback[currentTracebackIndex] = DP_DELETE | IType[currentPos] | DType;
					}
					// insert cannot follow delete and delete cannot follow insert
					I[currentPos] = I[currentPos] + gapExtendScore;
					IType[currentPos] = DP_INSERT_EXTEND;
					D = D + gapExtendScore;
					DType = DP_DELETE_EXTEND;
				}

				if (B[currentPos] + dropoffScore >= backwardMaxScore) {
					if (nextStartPos > currentPos) {
						nextStartPos = currentPos;
					}
					nextEndPos = currentPos;
				}

				currentPos++;
				currentTracebackIndex++;

			}

			nextEndPos++;
			if (nextEndPos == currentPos) {
				while (nextEndPos <= queryOffset && D + dropoffScore >= backwardMaxScore)  {
					if (nextEndPos >= dpCellAllocated) {
						fprintf(stderr, "HSPGappedExtensionWithTraceback(): Not enough DP cells allocated!\n");
						exit(1);
					}
					B[nextEndPos] = D;
					traceback[currentTracebackIndex] = DP_DELETE | DType;
					D = D + gapExtendScore;
					DType = DP_DELETE_EXTEND;
					I[nextEndPos] = DP_NEG_INFINITY;
					IType[nextEndPos] = 0;
					nextEndPos++;
					currentTracebackIndex++;
				}
				I[nextEndPos] = DP_NEG_INFINITY;
				IType[nextEndPos] = 0;
			}

			if (nextEndPos > queryOffset) {
				nextEndPos = queryOffset;
			}

			lastStartPos = startPos;
			startPos = nextStartPos;
			endPos = nextEndPos;

			textPos--;

		}

		// traceback
		backwardNumOfAlignment = 0;
		backwardNumOfAuxiliaryText = 0;
		textPos = backwardMaxScorePosText + 1;
		currentPos = queryOffset - backwardMaxScorePosQuery;
		dpScore = 0;	// for verifying traceback

		currentTracebackIndex = tracebackIndex[(textOffset - textPos + 1) * 2] 
								+ currentPos - tracebackIndex[(textOffset - textPos + 1) * 2 + 1];

		while (currentTracebackIndex > 0) {

			dpType = traceback[currentTracebackIndex] & DP_MASK;

			if (dpType == DP_MATCH_MISMATCH) {
				dpScore += scoringMatrix[textBuffer[textOffset - textPos + 1]][convertedKey[queryOffset - currentPos]];
				if (textBuffer[textOffset - textPos + 1] == convertedKey[queryOffset - currentPos] &&
					textBuffer[textOffset - textPos + 1] < 4) {	// match and not ambiguity
					tempBackwardAlignment[backwardNumOfAlignment] = ALIGN_MATCH;
				} else {
					tempBackwardAlignment[backwardNumOfAlignment] = ALIGN_MISMATCH_AMBIGUITY;
					tempBackwardAuxiliaryText[backwardNumOfAuxiliaryText] = textBuffer[textOffset - textPos + 1];
					backwardNumOfAuxiliaryText++;
				}
				textPos++;
				currentPos--;
				currentTracebackIndex = tracebackIndex[(textOffset - textPos + 1) * 2] 
										+ currentPos - tracebackIndex[(textOffset - textPos + 1) * 2 + 1];
			} else if (dpType == DP_INSERT) {
				while (!(traceback[currentTracebackIndex] & DP_INSERT_OPEN)) {
					dpScore += gapExtendScore;
					tempBackwardAlignment[backwardNumOfAlignment] = ALIGN_INSERT;
					backwardNumOfAlignment++;
					tempBackwardAuxiliaryText[backwardNumOfAuxiliaryText] = textBuffer[textOffset - textPos + 1];
					backwardNumOfAuxiliaryText++;
					textPos++;
					currentTracebackIndex = tracebackIndex[(textOffset - textPos + 1) * 2] 
											+ currentPos - tracebackIndex[(textOffset - textPos + 1) * 2 + 1];
				}
				dpScore += gapOpenScore + gapExtendScore;
				tempBackwardAlignment[backwardNumOfAlignment] = ALIGN_INSERT;
				tempBackwardAuxiliaryText[backwardNumOfAuxiliaryText] = textBuffer[textOffset - textPos + 1];
				backwardNumOfAuxiliaryText++;
				textPos++;
				currentTracebackIndex = tracebackIndex[(textOffset - textPos + 1) * 2] 
										+ currentPos - tracebackIndex[(textOffset - textPos + 1) * 2 + 1];
			} else {
				while (!(traceback[currentTracebackIndex] & DP_DELETE_OPEN)) {
					dpScore += gapExtendScore;
					tempBackwardAlignment[backwardNumOfAlignment] = ALIGN_DELETE;
					backwardNumOfAlignment++;
					currentPos--;
					currentTracebackIndex--;
				}
				dpScore += gapOpenScore + gapExtendScore;
				tempBackwardAlignment[backwardNumOfAlignment] = ALIGN_DELETE;
				currentPos--;
				currentTracebackIndex--;
			}

			backwardNumOfAlignment++;

		}

		if (dpScore != backwardMaxScore) {
			fprintf(stderr, "Backward gapped extension traceback error!\n");
			exit(1);
		}

		// check cutoff score
		evalue = stat_gapCalcEvalue(stat_gapNominal2normalized(forwardMaxScore + backwardMaxScore));
		if (evalue < maxEvalue) {
			gappedHitList[numberOfGappedHit].posText = backwardMaxScorePosText;
			gappedHitList[numberOfGappedHit].posQuery = backwardMaxScorePosQuery;
			gappedHitList[numberOfGappedHit].lengthQuery = forwardMaxScorePosQuery - backwardMaxScorePosQuery;
			gappedHitList[numberOfGappedHit].score = forwardMaxScore + backwardMaxScore;
			gappedHitList[numberOfGappedHit].lengthText = forwardMaxScorePosText - backwardMaxScorePosText;
			gappedHitList[numberOfGappedHit].dbSeqIndex = gappedHitList[hitProcessed].dbSeqIndex;

			// Store alignment and auxiliary text
			numOfAlignmentWord = (forwardNumOfAlignment + backwardNumOfAlignment + ALIGN_PER_WORD - 1) / ALIGN_PER_WORD;
			alignment = MMPoolDispatch(alignmentPool, numOfAlignmentWord * sizeof(unsigned int));
			for (i=0; i<numOfAlignmentWord; i++) {
				alignment[i] = 0;
			}
			if (forwardNumOfAuxiliaryText + backwardNumOfAuxiliaryText > 0) {
				numOfAuxiliaryTextWord = (forwardNumOfAuxiliaryText + backwardNumOfAuxiliaryText + AUX_TEXT_PER_WORD - 1) / AUX_TEXT_PER_WORD;
				auxiliaryText = MMPoolDispatch(alignmentPool, numOfAuxiliaryTextWord * sizeof(unsigned int));
				for (i=0; i<numOfAuxiliaryTextWord; i++) {
					auxiliaryText[i] = 0;
				}
			} else {
				auxiliaryText = NULL;
			}

			// backward alignment
			for (i=0; i<backwardNumOfAlignment; i++) {
				alignment[i/ALIGN_PER_WORD] |= (tempBackwardAlignment[i] << (BITS_IN_WORD - (i % ALIGN_PER_WORD + 1) * ALIGN_BIT));
			}
			for (i=0; i<backwardNumOfAuxiliaryText; i++) {
				auxiliaryText[i/AUX_TEXT_PER_WORD] |= (tempBackwardAuxiliaryText[i] << (BITS_IN_WORD - (i % AUX_TEXT_PER_WORD + 1) * AUX_TEXT_BIT));
			}

			// forward alignment
			for (i=0; i<forwardNumOfAlignment; i++) {
				alignment[(i+backwardNumOfAlignment)/ALIGN_PER_WORD] |= (tempForwardAlignment[forwardNumOfAlignment - i - 1] 
																		<< (BITS_IN_WORD - ((i+backwardNumOfAlignment) % ALIGN_PER_WORD + 1) * ALIGN_BIT));
			}
			for (i=0; i<forwardNumOfAuxiliaryText; i++) {
				auxiliaryText[(i+backwardNumOfAuxiliaryText)/AUX_TEXT_PER_WORD] |= (tempForwardAuxiliaryText[forwardNumOfAuxiliaryText - i - 1] 
																				   << (BITS_IN_WORD - ((i+backwardNumOfAuxiliaryText) % AUX_TEXT_PER_WORD + 1) * AUX_TEXT_BIT));
			}
			
			((GappedHitListWithAlignment*)(gappedHitList + numberOfGappedHit))->alignment = alignment;
			((GappedHitListWithAlignment*)(gappedHitList + numberOfGappedHit))->auxiliaryText = auxiliaryText;

			numberOfGappedHit++;

		}

		hitProcessed++;

	}

	// free working memory
	MMPoolReturn(mmPool, textBuffer, dpCellAllocated * 2);
	MMPoolReturn(mmPool, B, dpCellAllocated * sizeof(int));
	MMPoolReturn(mmPool, IType, dpCellAllocated);
	MMPoolReturn(mmPool, I, dpCellAllocated * sizeof(int));

	MMPoolReturn(mmPool, tracebackIndex, dpCellAllocated * 2 * 2 * sizeof(int));

	MMUnitFree(traceback, tracebackAllocated);

	MMPoolReturn(mmPool, tempForwardAlignment, dpCellAllocated * 2);
	MMPoolReturn(mmPool, tempBackwardAlignment, dpCellAllocated * 2);
	MMPoolReturn(mmPool, tempForwardAuxiliaryText, dpCellAllocated);
	MMPoolReturn(mmPool, tempBackwardAuxiliaryText, dpCellAllocated);

	return numberOfGappedHit;

}

// textList must be sorted in descending posText order;
// queryPosListIndex contains the starting positions of queryPos in queryPosList within a queryPosGroup
// a dummy entry must be present at the end of queryPosListIndex

int HSPDPDBvsQuery(const unsigned int *packedText, const HitList *textList, const int numOfTextHit,
				   const SeqOffset *textSeqOffset, const Ambiguity *textAmbiguity, const int numOfTextSeq,
				   const unsigned char *convertedKey, const int queryPatternLength, 
				   const int *queryPosListIndex, const unsigned int *queryPosList,
				   GappedHitListWithAlignment* __restrict gappedHitList, const int maxNumOfGappedHit,
				   MMPool *mmPool, MMPool *alignmentPool,
				   const int matchScore, const int mismatchScore,
				   const int gapOpenScore, const int gapExtendScore,
				   const int cutoffScore, const double maxEvalue) {

	#define NUM_SOURCE_BIT		12
	#define SOURCE_BIT_MASK		0xFFF

	#define DROPOFF_MAX_ENTRY		256

	#define EARLY_DROPOFF_MIN		11
	#define EARLY_DROPOFF_FACTOR	4	// These two are heuristic values to determine the dropoff to trigger addition traceback
										// When this heuistic is not applied, about 0.01% alignments reported by Blast are
										// filtered out during traceback on experiments with query size at chromosome lengths
										// The filtered alignments are of low evalues and are close to some strong alignments

	#define TEXT_BUFFER_SIZE	16

	// DP table
	DPCell** __restrict dpCell;
	int M;
	int D;				// insert and delete wrt query

	// The following variables are all indexed by source bit

	// Max score
	DPMaxScore* __restrict maxScore;
	DPMaxScore* __restrict dropoffMaxScore;

	// Keeping track whether dropoff should be valid alignments
	int* __restrict currentMaxScore;
	int* __restrict lastDropoffMaxScoreIndex;

	// End variables indexed by source bit

	int textBuffer[TEXT_BUFFER_SIZE];

	char* __restrict tempAlignment;
	char* __restrict tempAuxiliaryText;

	int scoringMatrix[16][16];
	int gapOpenScoreShifted, gapExtendScoreShifted;
	int cutoffScoreShifted;
	int dropoffScoreShifted;
	int minPositiveScore;

	int sourceBit;

	double evalue;

	int firstTextHitBeingProcessed, numOfTextHitBeingProcessed;

	unsigned int textOffset;
	int textChar;
	int charToProcess;

	int textDbSeqIndex, textAmbiguityIndex;
	unsigned int textSeqStart, textSeqEnd;

	int i, j;
	unsigned int c;

	int numOfGappedHit;

	int sourceBitMinScore, maxSourceBit, nextSourceBit, usedSourceBit;
	int numOfDropoffMaxScore;

	int	tempAlignmentIndex;
	int	tempAuxiliaryTextIndex;
	int numOfAlignmentWord, numOfAuxiliaryTextWord;

	Traceback* __restrict traceback;
	int tracebackAllocated;
	int tracebackIndex;
	int tracebackIndexAdjustment;

	unsigned int queryListInfo;

	int lastDpCellUsed, lastDpCellIndex;
	int dpCellUsed, dpCellIndex;
	int dpCellUsedAdjustedForEdge;

	int q;
	unsigned int qPos, lastQPos, maxQPos;
	int qChar;

	int bestScore, insertScore, deleteScore;

	int minAlignmentLength, minNegAlignmentLength, maxDPGroupingDist;

	// allocate working memory
	dpCell = MMPoolDispatch(mmPool, 2 * sizeof(DPCell*));
	dpCell[0] = MMPoolDispatch(mmPool, MAX_ALIGNMENT_LENGTH * sizeof(DPCell));
	dpCell[1] = MMPoolDispatch(mmPool, MAX_ALIGNMENT_LENGTH * sizeof(DPCell));

	maxScore = MMPoolDispatch(mmPool, (1 << NUM_SOURCE_BIT) * sizeof(DPMaxScore));
	dropoffMaxScore = MMPoolDispatch(mmPool, DROPOFF_MAX_ENTRY * sizeof(DPMaxScore));

	currentMaxScore = MMPoolDispatch(mmPool, (1 << NUM_SOURCE_BIT) * sizeof(int));
	lastDropoffMaxScoreIndex = MMPoolDispatch(mmPool, (1 << NUM_SOURCE_BIT) * sizeof(int));

	tempAlignment = MMPoolDispatch(mmPool, (MAX_ALIGNMENT_LENGTH * 2) * sizeof(char));
	tempAuxiliaryText = MMPoolDispatch(mmPool, (MAX_ALIGNMENT_LENGTH * 2) * sizeof(char));

	tracebackAllocated = MMPoolByteAvailable(mmPool) / sizeof(Traceback);
	if (tracebackAllocated < 0) {
		fprintf(stderr, "HSPDPDBvsDB(): Not enough memory allocated!\n");
		exit(1);
	}
	traceback =  MMPoolDispatch(mmPool, tracebackAllocated * sizeof(Traceback));

	HSPFillScoringMatrix(scoringMatrix, matchScore, mismatchScore, NUM_SOURCE_BIT);

	gapOpenScoreShifted = gapOpenScore * (1 << NUM_SOURCE_BIT);
	gapExtendScoreShifted = gapExtendScore * (1 << NUM_SOURCE_BIT);
	cutoffScoreShifted = cutoffScore * (1 << NUM_SOURCE_BIT);
	minPositiveScore = 1 << NUM_SOURCE_BIT;
	if (cutoffScore <= EARLY_DROPOFF_MIN) {
		dropoffScoreShifted = cutoffScoreShifted;
	} else {
		dropoffScoreShifted = (cutoffScore - (cutoffScore - EARLY_DROPOFF_MIN + EARLY_DROPOFF_FACTOR - 1) / EARLY_DROPOFF_FACTOR) * (1 << NUM_SOURCE_BIT) ;
	}

	minAlignmentLength = (cutoffScore + matchScore - 1) / matchScore;
	minNegAlignmentLength = ((gapOpenScore + gapExtendScore + matchScore + 1 - minAlignmentLength * matchScore + 1) / (gapOpenScore + gapExtendScore + matchScore)) * 2 - 1;
	maxDPGroupingDist = minAlignmentLength + minNegAlignmentLength;	// minimum length for an alignment to reached cutoff and then drop back to zero 

	// maxSourceBit are initially assigned to cells
	maxSourceBit = (1 << NUM_SOURCE_BIT) - 1;
	// source bits are assigned when score >= sourceBitMinScore
	sourceBitMinScore = (1 - gapOpenScore - gapExtendScore) * (1 << NUM_SOURCE_BIT) + maxSourceBit;

	numOfGappedHit = 0;		// Number of alignments 

	firstTextHitBeingProcessed = 0;

	// Clear last dropoff
	for (j=0; j<maxSourceBit; j++) {
		lastDropoffMaxScoreIndex[j] = 0;
	}
	numOfDropoffMaxScore = 0;
	usedSourceBit = 0;

	textOffset = 0;

	while (firstTextHitBeingProcessed < numOfTextHit) {

		if (textList[firstTextHitBeingProcessed].posText > textOffset) {
			// textOffset will be larger than those processed
			textDbSeqIndex = numOfTextSeq - 1;
			textSeqStart = textSeqOffset[textDbSeqIndex].startPos;
			textSeqEnd = textSeqOffset[textDbSeqIndex].endPos;
			textAmbiguityIndex = textSeqOffset[textDbSeqIndex].lastAmbiguityIndex;
		}

		textOffset = textList[firstTextHitBeingProcessed].posText;	// posText points to the character right after the hit
		queryListInfo = textList[firstTextHitBeingProcessed].info;

		// Find corresponding DB sequence for the text
		while (textSeqOffset[textDbSeqIndex].startPos >= textOffset) {
			textDbSeqIndex--;
			textSeqStart = textSeqOffset[textDbSeqIndex].startPos;
			textSeqEnd = textSeqOffset[textDbSeqIndex].endPos;
			textAmbiguityIndex = textSeqOffset[textDbSeqIndex].lastAmbiguityIndex;
		}

		if (textOffset % CHAR_PER_WORD > 0) {

			// Decode text to the next word boundary
			charToProcess = textOffset % CHAR_PER_WORD;
			c = packedText[textOffset / CHAR_PER_WORD] >> (BITS_IN_WORD - charToProcess * BIT_PER_CHAR);
			for (i=charToProcess; i--;) {	// from charToProcess - 1 to 0
				textBuffer[i] = (c & 0x3);
				c >>= BIT_PER_CHAR;
			}

			// Apply ambiguity
			while (textAmbiguity[textAmbiguityIndex].startPos >= textOffset) {
				textAmbiguityIndex--;
			}
			if (textAmbiguity[textAmbiguityIndex].rightOfEndPos > (textOffset - charToProcess)) {
				for (i=0; i<charToProcess; i++) {
					while (textAmbiguity[textAmbiguityIndex].startPos >= (textOffset - i)) {
						textAmbiguityIndex--;
					}
					if (textAmbiguity[textAmbiguityIndex].rightOfEndPos > (textOffset - i - 1)) {
						textBuffer[(textOffset - i - 1) % CHAR_PER_WORD] = textAmbiguity[textAmbiguityIndex].symbol;
					}
				}
			}

		}

		lastDpCellUsed = 0;
		lastDpCellIndex = 0;
		dpCellIndex = 1;
		maxQPos = 0;
		tracebackIndex = 0;
		tracebackIndexAdjustment = 0;

		// Clear last dropoff
		for (j=0; j<usedSourceBit; j++) {
			lastDropoffMaxScoreIndex[j] = 0;
		}
		numOfDropoffMaxScore = 0;

		numOfTextHitBeingProcessed = 0;
		nextSourceBit = 0;

		while (textOffset > textSeqOffset[textDbSeqIndex].startPos && 
			   (lastDpCellUsed > 0 || textList[firstTextHitBeingProcessed + numOfTextHitBeingProcessed].posText == textOffset)) {

			if (textOffset % CHAR_PER_WORD == 0) {

				// Decode text to the next word boundary
				c = packedText[textOffset / CHAR_PER_WORD - 1];
				for (i=CHAR_PER_WORD; i--;) {	// from CHAR_PER_WORD - 1 to 0
					textBuffer[i] = (c & 0x3);
					c >>= BIT_PER_CHAR;
				}

				// Apply ambiguity
				while (textAmbiguity[textAmbiguityIndex].startPos >= textOffset) {
					textAmbiguityIndex--;
				}
				if (textAmbiguity[textAmbiguityIndex].rightOfEndPos > (textOffset - CHAR_PER_WORD)) {
					for (i=0; i<CHAR_PER_WORD; i++) {
						while (textAmbiguity[textAmbiguityIndex].startPos >= (textOffset - i)) {
							textAmbiguityIndex--;
						}
						if (textAmbiguity[textAmbiguityIndex].rightOfEndPos > (textOffset - i - 1)) {
							textBuffer[(textOffset - i - 1) % CHAR_PER_WORD] = textAmbiguity[textAmbiguityIndex].symbol;
						}
					}
				}

			}

			// Set text char
			textChar = textBuffer[(textOffset - 1) % CHAR_PER_WORD];

			// Check if it is a valid starting position of alignment
			if (firstTextHitBeingProcessed + numOfTextHitBeingProcessed > numOfTextHit ||
				textList[firstTextHitBeingProcessed + numOfTextHitBeingProcessed].posText != textOffset ||
				textOffset + maxDPGroupingDist < textList[firstTextHitBeingProcessed].posText) {
//				(numOfTextHitBeingProcessed > 0 && textList[firstTextHitBeingProcessed + numOfTextHitBeingProcessed].posText + maxDPSkip <
//												   textList[firstTextHitBeingProcessed + numOfTextHitBeingProcessed - 1].posText)) {
			} else {
				// Merge starting positions in query into the dpCell
				queryListInfo = textList[firstTextHitBeingProcessed + numOfTextHitBeingProcessed].info;	// the group of query hits corresponding to the text hit
				q = queryPosListIndex[queryListInfo];												// the startingindex of the group of query hits

				dpCellUsed = 0;
				i = 0;
				while ((q < queryPosListIndex[queryListInfo + 1] ||
					    i < lastDpCellUsed) && dpCellUsed < MAX_ALIGNMENT_LENGTH) {
					while (i < lastDpCellUsed && (q >= queryPosListIndex[queryListInfo + 1] || dpCell[lastDpCellIndex][i].P >= queryPosList[q] 
											  && dpCellUsed < MAX_ALIGNMENT_LENGTH)) {
						dpCell[dpCellIndex][dpCellUsed] = dpCell[lastDpCellIndex][i];
						dpCellUsed++;
						i++;
					}
					if (q < queryPosListIndex[queryListInfo + 1] && i > 0 && dpCell[lastDpCellIndex][i-1].P == queryPosList[q]) {
						// The cell is positive; ignore starting position in query
						q++;
					}
					while (q < queryPosListIndex[queryListInfo + 1] && dpCellUsed < MAX_ALIGNMENT_LENGTH &&
						   (i >= lastDpCellUsed || dpCell[lastDpCellIndex][i].P < queryPosList[q])) {
						// insert starting position in query
						if (queryPosList[q] > (unsigned int)queryPatternLength) {
							fprintf(stderr, "HSPDPDBvsQuery(): Query position is larger than query length!\n");
							exit(1);
						}
						dpCell[dpCellIndex][dpCellUsed].P = queryPosList[q];		// posText points to the character right after the hit
						dpCell[dpCellIndex][dpCellUsed].B = maxSourceBit;			// initial score for a blank cell
						dpCell[dpCellIndex][dpCellUsed].I = DP_NEG_INFINITY;		// insertion not allowed here
						dpCellUsed++;
						tracebackIndexAdjustment--;
						q++;
					}
				}
				if (q < queryPosListIndex[queryListInfo + 1] || i < lastDpCellUsed) {
					// Not enough memory
				} else {
					numOfTextHitBeingProcessed++;
					dpCellIndex ^= 1;		// Swap dpCells
					lastDpCellIndex ^= 1;	// Swap dpCells
					lastDpCellUsed = dpCellUsed;
				}
			}

			D = DP_NEG_INFINITY;

			// Clear current max score
			usedSourceBit = nextSourceBit;
			for (j=0; j<usedSourceBit; j++) {
				currentMaxScore[j] = 0;
			}

			i = 0;
			dpCellUsed = 0;
			dpCellUsedAdjustedForEdge = FALSE;
			lastQPos = ALL_ONE_MASK;
			
			// check traceback buffer
			if (tracebackIndex + MAX_ALIGNMENT_LENGTH >= tracebackAllocated) {
				fprintf(stderr, "HSPDPDBvsQuery(): Not enough traceback buffer!\n");
				exit(1);
			}

			while (i < lastDpCellUsed) {

				// DP

				// check buffer
				if (dpCellUsed + 3 >= MAX_ALIGNMENT_LENGTH) {
					// 1 cell at most give 3 positive cells
					fprintf(stderr, "HSPDPDBvsQuery(): Not enough query buffer(1)!\n");
					exit(1);
				}

				if (dpCell[lastDpCellIndex][i].B == maxSourceBit) {
					// the cell is inserted through input query positions
					tracebackIndexAdjustment++;
				}

				qPos = dpCell[lastDpCellIndex][i].P - 1;

				if (qPos != lastQPos - 1) {
					// The cell on the left is non-positive; handle the insertion + deletion only
					if (dpCell[lastDpCellIndex][i].I >= minPositiveScore || D >= minPositiveScore) {

						traceback[tracebackIndex].indexDiff = dpCellUsed - i + lastDpCellUsed + tracebackIndexAdjustment;	// indexDiff always point to the cell directly above logically
						traceback[tracebackIndex].textChar = textChar;

						if (dpCell[lastDpCellIndex][i].I >= minPositiveScore && dpCell[lastDpCellIndex][i].I == dpCell[lastDpCellIndex][i].B + gapOpenScoreShifted + gapExtendScoreShifted) {
							traceback[tracebackIndex].IOpen = 1;		// Gap opening
						} else {
							traceback[tracebackIndex].IOpen = 0;		// Not gap opening
						}
						if (D >= minPositiveScore && D == dpCell[dpCellIndex][dpCellUsed-1].B + gapOpenScoreShifted + gapExtendScoreShifted) {
							traceback[tracebackIndex].DOpen = 1;		// Gap opening
						} else {
							traceback[tracebackIndex].DOpen = 0;		// Not gap opening
						}

						if (dpCell[lastDpCellIndex][i].I >= D) {
							dpCell[dpCellIndex][dpCellUsed].B = dpCell[lastDpCellIndex][i].I;
							traceback[tracebackIndex].alignment = DP_INSERT;
						} else {
							dpCell[dpCellIndex][dpCellUsed].B = D;
							traceback[tracebackIndex].alignment = DP_DELETE;
						}
						tracebackIndex++;

						dpCell[dpCellIndex][dpCellUsed].I = dpCell[lastDpCellIndex][i].I + gapExtendScoreShifted;	// insert cannot follow delete
						D = D + gapExtendScoreShifted;	// delete cannot follow insert 
						dpCell[dpCellIndex][dpCellUsed].P = qPos + 1;
						dpCellUsed++;
					}
				}
				
				qChar = convertedKey[qPos];

				// DP - Match/mismatch
				M = dpCell[lastDpCellIndex][i].B + scoringMatrix[textChar][qChar];

				if (M >= sourceBitMinScore &&
					(M & SOURCE_BIT_MASK) == maxSourceBit) {
					// Assign source bits
					if (nextSourceBit < maxSourceBit) {
						M += nextSourceBit - maxSourceBit;
						maxScore[nextSourceBit].score = 0;
						nextSourceBit++;
					} else {
						fprintf(stderr, "HSPDPDBvsQuery(): Not enough source bit space!\n");
					}
				}

				// Store the insert and delete score before updating them
				if (i+1 < lastDpCellUsed && dpCell[lastDpCellIndex][i+1].P == qPos) {
					// insert score available
					insertScore = dpCell[lastDpCellIndex][i+1].I;
				} else {
					insertScore = DP_NEG_INFINITY;
				}
				deleteScore = D;

				// Determine the best score
				if (M >= insertScore && M >= D) {
					// match score is maximum
					dpCell[dpCellIndex][dpCellUsed].B = M;
					if (M + gapOpenScoreShifted >= insertScore) {
						dpCell[dpCellIndex][dpCellUsed].I = M + gapOpenScoreShifted + gapExtendScoreShifted;
					} else {
						dpCell[dpCellIndex][dpCellUsed].I = insertScore + gapExtendScoreShifted;
					}
					if (M + gapOpenScoreShifted >= D) {
						D = M + gapOpenScoreShifted + gapExtendScoreShifted;
					} else {
						D = D + gapExtendScoreShifted;
					}
					traceback[tracebackIndex].alignment = DP_MATCH_MISMATCH;
				} else {
					if (insertScore >= D) {
						dpCell[dpCellIndex][dpCellUsed].B = insertScore;
						traceback[tracebackIndex].alignment = DP_INSERT;
					} else {
						dpCell[dpCellIndex][dpCellUsed].B = D;
						traceback[tracebackIndex].alignment = DP_DELETE;
					}
					dpCell[dpCellIndex][dpCellUsed].I = insertScore + gapExtendScoreShifted;
					D = D + gapExtendScoreShifted;
				}

				// check if the best score is positive
				if (dpCell[dpCellIndex][dpCellUsed].B >= minPositiveScore) {

					traceback[tracebackIndex].indexDiff = dpCellUsed - i -1 + lastDpCellUsed + tracebackIndexAdjustment;	// indexDiff always point to the cell directly above logically
					traceback[tracebackIndex].textChar = textChar;
					traceback[tracebackIndex].queryChar = qChar;
					if (insertScore >= minPositiveScore && insertScore == dpCell[lastDpCellIndex][i+1].B + gapOpenScoreShifted + gapExtendScoreShifted) {
						traceback[tracebackIndex].IOpen = 1;		// Gap opening
					} else {
						traceback[tracebackIndex].IOpen = 0;		// Not gap opening
					}
					if (deleteScore >= minPositiveScore && deleteScore == dpCell[dpCellIndex][dpCellUsed-1].B + gapOpenScoreShifted + gapExtendScoreShifted) {
						traceback[tracebackIndex].DOpen = 1;		// Gap opening
					} else {
						traceback[tracebackIndex].DOpen = 0;		// Not gap opening
					}

					dpCell[dpCellIndex][dpCellUsed].P = qPos;

					sourceBit = dpCell[dpCellIndex][dpCellUsed].B & SOURCE_BIT_MASK;
					if (dpCell[dpCellIndex][dpCellUsed].B > maxScore[sourceBit].score) {
						maxScore[sourceBit].score = dpCell[dpCellIndex][dpCellUsed].B;
						maxScore[sourceBit].tracebackIndex = tracebackIndex;
						maxScore[sourceBit].posText = textOffset - 1;
						maxScore[sourceBit].posQuery = qPos;
					}
					if (sourceBit < usedSourceBit && dpCell[dpCellIndex][dpCellUsed].B > currentMaxScore[sourceBit]) {
						currentMaxScore[sourceBit] = dpCell[dpCellIndex][dpCellUsed].B;
					}

					dpCellUsed++;
					tracebackIndex++;

				}

				lastQPos = qPos;

				if ((i+1 >= lastDpCellUsed || dpCell[lastDpCellIndex][i+1].P < (qPos-1)) && qPos > 0 && D >= minPositiveScore) {

					// delete score is still positive;
					qPos--;
					dpCell[dpCellIndex][dpCellUsed].B = D;
					dpCell[dpCellIndex][dpCellUsed].I = DP_NEG_INFINITY;
					dpCell[dpCellIndex][dpCellUsed].P = qPos;
					D = D + gapExtendScoreShifted;

					traceback[tracebackIndex].indexDiff = traceback[tracebackIndex-1].indexDiff;
					traceback[tracebackIndex].textChar = textChar;
					traceback[tracebackIndex].IOpen = 0;		// Not gap opening
					if (dpCell[dpCellIndex][dpCellUsed].B == dpCell[dpCellIndex][dpCellUsed-1].B + gapOpenScoreShifted + gapExtendScoreShifted) {
						traceback[tracebackIndex].DOpen = 1;		// Gap opening
					} else {
						traceback[tracebackIndex].DOpen = 0;		// Not gap opening
					}
					traceback[tracebackIndex].alignment = DP_DELETE;

					dpCellUsed++;
					tracebackIndex++;
				}
				if (qPos == 0) {
					// Reached the start of query
					D = DP_NEG_INFINITY;
					if (dpCellUsed > 0 && dpCell[dpCellIndex][dpCellUsed - 1].P == 0) {
						// No need to process the cell on the start of query
						dpCellUsed--;
						dpCellUsedAdjustedForEdge = TRUE;
					}
					break;
				}

				i++;

			}

			if (dpCellUsedAdjustedForEdge) {
				tracebackIndexAdjustment = 1;
			} else {
				tracebackIndexAdjustment = 0;
			}

			// Add to dropoff list if maxScore > currentMaxScore + cutoff
			for (j=0; j<usedSourceBit; j++) {
				if (maxScore[j].score >= currentMaxScore[j] + dropoffScoreShifted && lastDropoffMaxScoreIndex[j] != maxScore[j].tracebackIndex) {
					if (numOfDropoffMaxScore < DROPOFF_MAX_ENTRY) {
						dropoffMaxScore[numOfDropoffMaxScore] = maxScore[j];
						lastDropoffMaxScoreIndex[j] = maxScore[j].tracebackIndex;
						numOfDropoffMaxScore++;
					} else {
						fprintf(stderr, "HSPDPDBvsQuery(): Not enough dropoff space!\n");
					}
				}
			}

			dpCellIndex ^= 1;		// Swap dpCells
			lastDpCellIndex ^= 1;	// Swap dpCells
			lastDpCellUsed = dpCellUsed;

			textOffset--;

		}

		// Move scores reached cutoff to front
		j = 0;
		for (i=0; i<nextSourceBit; i++) {
			if (maxScore[i].score >= cutoffScoreShifted) {
				maxScore[j] = maxScore[i];
				j++;
			}
		}
		nextSourceBit = j;

		// Append dropoff max score
		for (i=0; i<numOfDropoffMaxScore; i++) {
			if (nextSourceBit < (1<<NUM_SOURCE_BIT)) {
				for (j=0; j<nextSourceBit; j++) {
					if (maxScore[j].tracebackIndex == dropoffMaxScore[i].tracebackIndex) {
						break;
					}
				}
				if (j >= nextSourceBit) {
					maxScore[nextSourceBit] = dropoffMaxScore[i];
					nextSourceBit++;
				}
			} else {
				fprintf(stderr, "HSPDPDBvsQuery(): Not enough source bit for dropoff!\n");
			}
		}

		// Traceback alignments
		for (sourceBit=0; sourceBit<nextSourceBit; sourceBit++) {
			
			bestScore = maxScore[sourceBit].score;
			evalue = stat_gapCalcEvalue(stat_gapNominal2normalized(bestScore >> NUM_SOURCE_BIT));
			if (evalue <= maxEvalue) {

				// Traceback
				tempAlignmentIndex = 0;
				tempAuxiliaryTextIndex = 0;

				textOffset = maxScore[sourceBit].posText;
				qPos = maxScore[sourceBit].posQuery;

				tracebackIndex = maxScore[sourceBit].tracebackIndex;

				while (bestScore >= minPositiveScore) {

					if (traceback[tracebackIndex].alignment == DP_MATCH_MISMATCH) {
						if (tempAlignmentIndex >= MAX_ALIGNMENT_LENGTH) {
							fprintf(stderr, "HSPDPDBvsQuery(): Not enough temp alignment buffer!\n");
							exit(1);
						}
						textChar = traceback[tracebackIndex].textChar;
						qChar = traceback[tracebackIndex].queryChar;
						bestScore -= scoringMatrix[textChar][qChar];
						if (textChar == qChar && textChar < 4) {	// match and not ambiguity
							tempAlignment[tempAlignmentIndex] = ALIGN_MATCH;
						} else {
							tempAlignment[tempAlignmentIndex] = ALIGN_MISMATCH_AMBIGUITY;
							tempAuxiliaryText[tempAuxiliaryTextIndex] = (char)traceback[tracebackIndex].textChar;	// store text character if mismatch or ambiguity
							tempAuxiliaryTextIndex++;
						}
						textOffset++;
						qPos++;
						tracebackIndex -= traceback[tracebackIndex].indexDiff + 1;	// indexDiff always refer to the cell directly above logically
					} else if (traceback[tracebackIndex].alignment == DP_INSERT) {
						while (!traceback[tracebackIndex].IOpen) {
							if (tempAlignmentIndex >= MAX_ALIGNMENT_LENGTH) {
								fprintf(stderr, "HSPDPDBvsQuery(): Not enough temp alignment buffer!\n");
								exit(1);
							}
							bestScore -= gapExtendScoreShifted;
							tempAlignment[tempAlignmentIndex] = ALIGN_INSERT;
							tempAlignmentIndex++;
							tempAuxiliaryText[tempAuxiliaryTextIndex] = (char)traceback[tracebackIndex].textChar;
							tempAuxiliaryTextIndex++;
							tracebackIndex -= traceback[tracebackIndex].indexDiff;	// indexDiff always refer to the cell directly above logically
							textOffset++;
						}
						if (tempAlignmentIndex >= MAX_ALIGNMENT_LENGTH) {
							fprintf(stderr, "HSPDPDBvsQuery(): Not enough temp alignment buffer!\n");
							exit(1);
						}
						bestScore -= gapOpenScoreShifted + gapExtendScoreShifted;
						tempAlignment[tempAlignmentIndex] = ALIGN_INSERT;
						tempAuxiliaryText[tempAuxiliaryTextIndex] = (char)traceback[tracebackIndex].textChar;
						tempAuxiliaryTextIndex++;
						tracebackIndex -= traceback[tracebackIndex].indexDiff;	// indexDiff always refer to the cell directly above logically
						textOffset++;
					} else {
						while (!traceback[tracebackIndex].DOpen) {
							if (tempAlignmentIndex >= MAX_ALIGNMENT_LENGTH) {
								fprintf(stderr, "HSPDPDBvsQuery(): Not enough temp alignment buffer!\n");
								exit(1);
							}
							bestScore -= gapExtendScoreShifted;
							tempAlignment[tempAlignmentIndex] = ALIGN_DELETE;
							tempAlignmentIndex++;
							tracebackIndex--;
							qPos++;
						}
						if (tempAlignmentIndex >= MAX_ALIGNMENT_LENGTH) {
							fprintf(stderr, "HSPDPDBvsQuery(): Not enough temp alignment buffer!\n");
							exit(1);
						}
						bestScore -= gapOpenScoreShifted + gapExtendScoreShifted;
						tempAlignment[tempAlignmentIndex] = ALIGN_DELETE;
						tracebackIndex--;
						qPos++;
					}

					tempAlignmentIndex++;

				}

				if (numOfGappedHit >= maxNumOfGappedHit) {
					fprintf(stderr, "HSPDPDBvsQuery(): Not enough gapped hit buffer!\n");
					exit(1);
				}
				gappedHitList[numOfGappedHit].posText = maxScore[sourceBit].posText;
				gappedHitList[numOfGappedHit].posQuery = maxScore[sourceBit].posQuery;
				gappedHitList[numOfGappedHit].lengthText = textOffset - maxScore[sourceBit].posText;
				gappedHitList[numOfGappedHit].lengthQuery = qPos - maxScore[sourceBit].posQuery;
				gappedHitList[numOfGappedHit].score = maxScore[sourceBit].score >> NUM_SOURCE_BIT;
				gappedHitList[numOfGappedHit].dbSeqIndex = textDbSeqIndex;

				numOfAlignmentWord = (tempAlignmentIndex + ALIGN_PER_WORD - 1) / ALIGN_PER_WORD;
				gappedHitList[numOfGappedHit].alignment = MMPoolDispatch(alignmentPool, numOfAlignmentWord * sizeof(unsigned int));
				for (i=0; i<numOfAlignmentWord; i++) {
					gappedHitList[numOfGappedHit].alignment[i] = 0;
				}

				for (i=0; i<tempAlignmentIndex; i++) {
					gappedHitList[numOfGappedHit].alignment[i/ALIGN_PER_WORD] |= (tempAlignment[i] << (BITS_IN_WORD - (i % ALIGN_PER_WORD + 1) * ALIGN_BIT));
				}

				if (tempAuxiliaryTextIndex > 0) {
					numOfAuxiliaryTextWord = (tempAuxiliaryTextIndex + AUX_TEXT_PER_WORD - 1) / AUX_TEXT_PER_WORD;
					gappedHitList[numOfGappedHit].auxiliaryText = MMPoolDispatch(alignmentPool, numOfAuxiliaryTextWord * sizeof(unsigned int));
					for (i=0; i<numOfAuxiliaryTextWord; i++) {
						gappedHitList[numOfGappedHit].auxiliaryText[i] = 0;
					}
					for (i=0; i<tempAuxiliaryTextIndex; i++) {
						gappedHitList[numOfGappedHit].auxiliaryText[i/AUX_TEXT_PER_WORD] |= (tempAuxiliaryText[i] << (BITS_IN_WORD - (i % AUX_TEXT_PER_WORD + 1) * AUX_TEXT_BIT));
					}
				} else {
					gappedHitList[numOfGappedHit].auxiliaryText = NULL;
				}

				numOfGappedHit++;

			}
		}

		firstTextHitBeingProcessed += numOfTextHitBeingProcessed;

	}

	// free working memory
	MMPoolReturn(mmPool, traceback, tracebackAllocated * sizeof(int));

	MMPoolReturn(mmPool, tempAuxiliaryText, (MAX_ALIGNMENT_LENGTH * 2) * sizeof(char));
	MMPoolReturn(mmPool, tempAlignment, (MAX_ALIGNMENT_LENGTH * 2) * sizeof(char));

	MMPoolReturn(mmPool, lastDropoffMaxScoreIndex, (1 << NUM_SOURCE_BIT) * sizeof(int));
	MMPoolReturn(mmPool, currentMaxScore, (1 << NUM_SOURCE_BIT) * sizeof(int));

	MMPoolReturn(mmPool, dropoffMaxScore, DROPOFF_MAX_ENTRY * sizeof(DPMaxScore));
	MMPoolReturn(mmPool, maxScore, (1 << NUM_SOURCE_BIT) * sizeof(DPMaxScore));

	MMPoolReturn(mmPool, dpCell[0], MAX_ALIGNMENT_LENGTH * sizeof(DPCell));
	MMPoolReturn(mmPool, dpCell[1], MAX_ALIGNMENT_LENGTH * sizeof(DPCell));
	MMPoolReturn(mmPool, dpCell, 2 * sizeof(DPCell*));

	return numOfGappedHit;


}


unsigned int HSPRemoveDuplicateUngappedHit(HitListWithPosQuery* __restrict ungappedHitList, const unsigned int numberOfHit) {

	unsigned int i, j;

	QSort(ungappedHitList, numberOfHit, sizeof(HitListWithPosQuery), HitListPosTextQueryOrder);

	i = j = 0;
	while (i<numberOfHit) {
		ungappedHitList[j].posText = ungappedHitList[i].posText;
		ungappedHitList[j].posQuery = ungappedHitList[i].posQuery;
		i++;
		while (i<numberOfHit && ungappedHitList[i].posText == ungappedHitList[j].posText 
							 && ungappedHitList[i].posQuery == ungappedHitList[j].posQuery) {
			i++;
		}
		j++;
	}

	return j;

}

unsigned int HSPSplitCrossBoundaryGappedHit(GappedHitList* __restrict gappedHitList, const unsigned int numberOfHit, const unsigned int queryPatternLength,
								   const SeqOffset *seqOffset,
								   const int matchScore, const int cutoffScore) {

	unsigned int i;
	unsigned int numberOfSplitHit;
	unsigned int dbSeqIndex;
	unsigned int splitDbSeqIndex;
	unsigned int minHitLength;
	unsigned int oldHitLength;
	unsigned int newHitLength;
	unsigned int hitEnd;
	long long diagonalPos;


	minHitLength = (cutoffScore + matchScore - 1) / matchScore;

	QSort(gappedHitList, numberOfHit, sizeof(GappedHitList), GappedHitListPosTextQueryLengthScoreOrder);

	i = 0;
	numberOfSplitHit = 0;
	dbSeqIndex = 0;

	while (i<numberOfHit) {

		// Find corresponding DB sequence
		while (seqOffset[dbSeqIndex].endPos < gappedHitList[i].posText) {
			dbSeqIndex++;
		}
		gappedHitList[i].dbSeqIndex = dbSeqIndex;
		diagonalPos = (long long)gappedHitList[i].ungappedPosText - (long long)gappedHitList[i].ungappedPosQuery;
		oldHitLength = gappedHitList[i].lengthText;

		splitDbSeqIndex = dbSeqIndex + 1;
		while (seqOffset[splitDbSeqIndex].startPos <= gappedHitList[i].posText + gappedHitList[i].lengthText - 1 - minHitLength) {

			// the hit span across 2 db sequences; and the portion on the right > minimum hit length
			
			gappedHitList[numberOfHit + numberOfSplitHit].posText = seqOffset[splitDbSeqIndex].startPos;
			if (gappedHitList[numberOfHit + numberOfSplitHit].posText - diagonalPos < queryPatternLength) {
				gappedHitList[numberOfHit + numberOfSplitHit].posQuery = (unsigned int)((long long)gappedHitList[numberOfHit + numberOfSplitHit].posText - diagonalPos);
			} else {
				gappedHitList[numberOfHit + numberOfSplitHit].posQuery = queryPatternLength - 1;
			}
			hitEnd = gappedHitList[i].posText + gappedHitList[i].lengthText - 1;
			if (hitEnd > seqOffset[splitDbSeqIndex].endPos) {
				hitEnd = seqOffset[splitDbSeqIndex].endPos;
			}
			newHitLength = hitEnd - gappedHitList[numberOfHit + numberOfSplitHit].posText + 1;
			gappedHitList[numberOfHit + numberOfSplitHit].lengthText = newHitLength;
			if (gappedHitList[numberOfHit + numberOfSplitHit].lengthQuery >= oldHitLength - newHitLength) {
				gappedHitList[numberOfHit + numberOfSplitHit].lengthQuery += newHitLength - oldHitLength;
			} else {
				gappedHitList[numberOfHit + numberOfSplitHit].lengthQuery = 0;
			}

			if (gappedHitList[i].ungappedPosText < seqOffset[splitDbSeqIndex].startPos) {
				gappedHitList[numberOfHit + numberOfSplitHit].ungappedPosText = seqOffset[splitDbSeqIndex].startPos;
			} else if (gappedHitList[i].ungappedPosText > hitEnd) {
				gappedHitList[numberOfHit + numberOfSplitHit].ungappedPosText = hitEnd;
			} else {
				gappedHitList[numberOfHit + numberOfSplitHit].ungappedPosText = gappedHitList[i].ungappedPosText;
			}
			if (gappedHitList[numberOfHit + numberOfSplitHit].ungappedPosText < diagonalPos) {
				gappedHitList[numberOfHit + numberOfSplitHit].ungappedPosQuery = 0;
			} else if (gappedHitList[numberOfHit + numberOfSplitHit].ungappedPosText - diagonalPos >= queryPatternLength) {
				gappedHitList[numberOfHit + numberOfSplitHit].ungappedPosQuery = queryPatternLength - 1;
			} else {
				gappedHitList[numberOfHit + numberOfSplitHit].ungappedPosQuery = (unsigned int)((long long)gappedHitList[numberOfHit + numberOfSplitHit].ungappedPosText - diagonalPos);
			}

			if (gappedHitList[i].ungappedPosText >= seqOffset[splitDbSeqIndex].startPos &&
				gappedHitList[i].ungappedPosText <= hitEnd) {
				gappedHitList[numberOfHit + numberOfSplitHit].score = gappedHitList[i].score;
			} else {
				gappedHitList[numberOfHit + numberOfSplitHit].score = 0;
			}
			gappedHitList[numberOfHit + numberOfSplitHit].dbSeqIndex = splitDbSeqIndex;

			numberOfSplitHit++;
			splitDbSeqIndex++;

		}

		hitEnd = gappedHitList[i].posText + gappedHitList[i].lengthText - 1;
		if (hitEnd > seqOffset[splitDbSeqIndex].endPos) {
			hitEnd = seqOffset[splitDbSeqIndex].endPos;
		}
		if (gappedHitList[i].posText < seqOffset[dbSeqIndex].startPos) {
			gappedHitList[i].posText = seqOffset[dbSeqIndex].startPos;
			if (gappedHitList[i].posText - diagonalPos < queryPatternLength) {
				gappedHitList[i].posQuery = (int)((long long)gappedHitList[i].posText - diagonalPos);
			} else {
				gappedHitList[i].posQuery = queryPatternLength - 1;
			}
		}
		newHitLength = hitEnd - gappedHitList[i].posText + 1;
		gappedHitList[i].lengthText = newHitLength;
		if (gappedHitList[i].lengthQuery >= oldHitLength - newHitLength) {
			gappedHitList[i].lengthQuery += newHitLength - oldHitLength;
		} else {
			gappedHitList[i].lengthQuery = 0;
		}

		if (gappedHitList[i].ungappedPosText < seqOffset[dbSeqIndex].startPos) {
			gappedHitList[i].ungappedPosText = seqOffset[dbSeqIndex].startPos;
		} else if (gappedHitList[i].ungappedPosText > hitEnd) {
			gappedHitList[i].ungappedPosText = hitEnd;
		} else {
			gappedHitList[i].ungappedPosText = gappedHitList[i].ungappedPosText;
		}
		if (gappedHitList[i].ungappedPosText < diagonalPos) {
			gappedHitList[i].ungappedPosQuery = 0;
		} else if (gappedHitList[i].ungappedPosText - diagonalPos >= queryPatternLength) {
			gappedHitList[i].ungappedPosQuery = queryPatternLength - 1;
		} else {
			gappedHitList[i].ungappedPosQuery = (int)((long long)gappedHitList[i].ungappedPosText - diagonalPos);
		}

		if (gappedHitList[i].ungappedPosText >= seqOffset[dbSeqIndex].startPos &&
			gappedHitList[i].ungappedPosText <= hitEnd) {
			gappedHitList[i].score = gappedHitList[i].score;
		} else {
			gappedHitList[i].score = 0;
		}
		gappedHitList[i].dbSeqIndex = dbSeqIndex;

		if (newHitLength >= minHitLength) {
			i++;
		}

	}

	return numberOfHit + numberOfSplitHit;

}
						 
unsigned int HSPRemoveDuplicateGappedHit(GappedHitList* __restrict gappedHitList, const unsigned int numberOfHit) {

	unsigned int i, j;

	QSort(gappedHitList, numberOfHit, sizeof(GappedHitList), GappedHitListPosTextQueryLengthScoreOrder);

	i = j = 0;
	while (i<numberOfHit) {
		gappedHitList[j].posText = gappedHitList[i].posText;
		gappedHitList[j].posQuery = gappedHitList[i].posQuery;
		gappedHitList[j].lengthText = gappedHitList[i].lengthText;
		gappedHitList[j].lengthQuery = gappedHitList[i].lengthQuery;
		gappedHitList[j].score = gappedHitList[i].score;
		gappedHitList[j].dbSeqIndex = gappedHitList[i].dbSeqIndex;
		gappedHitList[j].ungappedPosText = gappedHitList[i].ungappedPosText;
		gappedHitList[j].ungappedPosQuery = gappedHitList[i].ungappedPosQuery;
		i++;
		while (i<numberOfHit && gappedHitList[i].posText == gappedHitList[j].posText 
							 && gappedHitList[i].posQuery == gappedHitList[j].posQuery
							 && gappedHitList[i].lengthText == gappedHitList[j].lengthText
							 && gappedHitList[i].lengthQuery == gappedHitList[j].lengthQuery) {
			i++;
		}
		j++;
	}

	return j;

}
						 
// alignment and auxiliary text must be allocated through MMPool
// otherwise the pointer will be lost without freeing the memory
unsigned int HSPFinalFilter(GappedHitListWithAlignment* __restrict gappedHitList, const unsigned int numberOfHit) {

	unsigned int i, j;
	unsigned int hitRemaining;

	hitRemaining = numberOfHit;

	// First filter among hits with common starting points the longer and lower scoring hits

	QSort(gappedHitList, hitRemaining, sizeof(GappedHitList), GappedHitListPosTextQueryScoreLengthOrder);

	j = 0;
	while (j<hitRemaining) {
		if (gappedHitList[j].alignment != NULL) {	// not marked as discarded
			i = j + 1;
			while (i<hitRemaining && gappedHitList[i].posText == gappedHitList[j].posText
								  && gappedHitList[i].posQuery == gappedHitList[j].posQuery) {
				if (gappedHitList[i].alignment != NULL) {
					if (gappedHitList[i].lengthText >= gappedHitList[j].lengthText &&
						gappedHitList[i].lengthQuery >= gappedHitList[j].lengthQuery &&
						gappedHitList[i].score <= gappedHitList[j].score) {
						// discard hit
						gappedHitList[i].alignment = NULL;
					}
				}
				i++;
			}
		}
		j++;
	}

	// Pack retained hit to the beginning
	i = j = 0;
	while (i<hitRemaining) {
		if (gappedHitList[i].alignment != NULL) {	// not marked as discarded
			gappedHitList[j].posText = gappedHitList[i].posText;
			gappedHitList[j].posQuery = gappedHitList[i].posQuery;
			gappedHitList[j].lengthText = gappedHitList[i].lengthText;
			gappedHitList[j].lengthQuery = gappedHitList[i].lengthQuery;
			gappedHitList[j].score = gappedHitList[i].score;
			gappedHitList[j].dbSeqIndex = gappedHitList[i].dbSeqIndex;
			gappedHitList[j].alignment = gappedHitList[i].alignment;
			gappedHitList[j].auxiliaryText = gappedHitList[i].auxiliaryText;
			j++;
		}
		i++;
	}

	hitRemaining = j;

	// Second filter among hits with common ending points the longer and lower scoring hits

	QSort(gappedHitList, hitRemaining, sizeof(GappedHitList), GappedHitListEndTextQueryScoreLengthOrder);

	j = 0;
	while (j<hitRemaining) {
		if (gappedHitList[j].alignment != NULL) {	// not marked as discarded
			i = j + 1;
			while (i<hitRemaining && gappedHitList[i].posText + gappedHitList[i].lengthText == gappedHitList[j].posText + gappedHitList[j].lengthText 
							  && gappedHitList[i].posQuery + gappedHitList[i].lengthQuery == gappedHitList[j].posQuery + gappedHitList[j].lengthQuery) {
				if (gappedHitList[i].alignment != NULL) {
					if (gappedHitList[i].lengthText >= gappedHitList[j].lengthText &&
						gappedHitList[i].lengthQuery >= gappedHitList[j].lengthQuery &&
						gappedHitList[i].score <= gappedHitList[j].score) {
						// discard hit
						gappedHitList[i].alignment = NULL;
					}
				}
				i++;
			}
		}
		j++;
	}

	// Pack retained hit to the beginning
	i = j = 0;
	while (i<hitRemaining) {
		if (gappedHitList[i].alignment != NULL) {	// not marked as discarded
			gappedHitList[j].posText = gappedHitList[i].posText;
			gappedHitList[j].posQuery = gappedHitList[i].posQuery;
			gappedHitList[j].lengthText = gappedHitList[i].lengthText;
			gappedHitList[j].lengthQuery = gappedHitList[i].lengthQuery;
			gappedHitList[j].score = gappedHitList[i].score;
			gappedHitList[j].dbSeqIndex = gappedHitList[i].dbSeqIndex;
			gappedHitList[j].alignment = gappedHitList[i].alignment;
			gappedHitList[j].auxiliaryText = gappedHitList[i].auxiliaryText;
			j++;
		}
		i++;
	}

	hitRemaining = j;

	// Third filter hits with common start points

	QSort(gappedHitList, hitRemaining, sizeof(GappedHitList), GappedHitListPosTextQueryScoreLengthOrder);

	i = j = 0;
	while (i<hitRemaining) {
		gappedHitList[j].posText = gappedHitList[i].posText;
		gappedHitList[j].posQuery = gappedHitList[i].posQuery;
		gappedHitList[j].lengthText = gappedHitList[i].lengthText;
		gappedHitList[j].lengthQuery = gappedHitList[i].lengthQuery;
		gappedHitList[j].score = gappedHitList[i].score;
		gappedHitList[j].dbSeqIndex = gappedHitList[i].dbSeqIndex;
		gappedHitList[j].alignment = gappedHitList[i].alignment;
		gappedHitList[j].auxiliaryText = gappedHitList[i].auxiliaryText;
		i++;
		while (i<hitRemaining && gappedHitList[i].posText == gappedHitList[j].posText 
							  && gappedHitList[i].posQuery == gappedHitList[j].posQuery) {
			// discard hit
			gappedHitList[i].alignment = NULL;
			i++;
		}
		j++;
	}

	hitRemaining = j;

	// Fourth filter hits with common end points

	QSort(gappedHitList, hitRemaining, sizeof(GappedHitListWithAlignment), GappedHitListEndTextQueryScoreLengthOrder);

	i = j = 0;
	while (i<hitRemaining) {
		gappedHitList[j].posText = gappedHitList[i].posText;
		gappedHitList[j].posQuery = gappedHitList[i].posQuery;
		gappedHitList[j].lengthText = gappedHitList[i].lengthText;
		gappedHitList[j].lengthQuery = gappedHitList[i].lengthQuery;
		gappedHitList[j].score = gappedHitList[i].score;
		gappedHitList[j].dbSeqIndex = gappedHitList[i].dbSeqIndex;
		gappedHitList[j].alignment = gappedHitList[i].alignment;
		gappedHitList[j].auxiliaryText = gappedHitList[i].auxiliaryText;
		i++;
		while (i<hitRemaining && gappedHitList[i].posText + gappedHitList[i].lengthText == gappedHitList[j].posText + gappedHitList[j].lengthText 
							  && gappedHitList[i].posQuery + gappedHitList[i].lengthQuery == gappedHitList[j].posQuery + gappedHitList[j].lengthQuery) {
			// discard hit
			gappedHitList[i].alignment = NULL;
			i++;
		}
		j++;
	}

	hitRemaining = j;

	// Fifth filter enclosed hits with lower scores

	QSort(gappedHitList, hitRemaining, sizeof(GappedHitList), GappedHitListEnclosedScoreOrder);

	j = 0;
	while (j<hitRemaining) {
		if (gappedHitList[j].alignment != NULL) {	// not marked as discarded
			i = j + 1;
			while (i<hitRemaining && gappedHitList[i].posText < gappedHitList[j].posText + gappedHitList[j].lengthText) {
				if (gappedHitList[i].alignment != NULL) {
					if (gappedHitList[i].posText + gappedHitList[i].lengthText <= gappedHitList[j].posText + gappedHitList[j].lengthText &&
						gappedHitList[i].posQuery >= gappedHitList[j].posQuery &&
						gappedHitList[i].posQuery + gappedHitList[i].lengthQuery <= gappedHitList[j].posQuery + gappedHitList[j].lengthQuery &&
						gappedHitList[i].score <= gappedHitList[j].score) {
						// discard hit
						gappedHitList[i].alignment = NULL;
					}
				}
				i++;
			}
		}
		j++;
	}

	// Pack retained hit to the beginning
	i = j = 0;
	while (i<hitRemaining) {
		if (gappedHitList[i].alignment != NULL) {	// not marked as discarded
			gappedHitList[j].posText = gappedHitList[i].posText;
			gappedHitList[j].posQuery = gappedHitList[i].posQuery;
			gappedHitList[j].lengthText = gappedHitList[i].lengthText;
			gappedHitList[j].lengthQuery = gappedHitList[i].lengthQuery;
			gappedHitList[j].score = gappedHitList[i].score;
			gappedHitList[j].dbSeqIndex = gappedHitList[i].dbSeqIndex;
			gappedHitList[j].alignment = gappedHitList[i].alignment;
			gappedHitList[j].auxiliaryText = gappedHitList[i].auxiliaryText;
			j++;
		}
		i++;
	}

	hitRemaining = j;

	return hitRemaining;

}
						 
HSPUngappedExtLookupTable *HSPGenUngappedExtLookupTable(MMPool *mmPool, const int matchScore, const int mismatchScore) {

	HSPUngappedExtLookupTable *hspUngappedExtLookupTable;

	hspUngappedExtLookupTable = MMPoolDispatch(mmPool, (1 << 16) * sizeof(HSPUngappedExtLookupTable));

	HSPCalUngappedExtLookupTable(hspUngappedExtLookupTable, matchScore, mismatchScore);

	return hspUngappedExtLookupTable;
}

void HSPCalUngappedExtLookupTable(HSPUngappedExtLookupTable *hspUngappedExtLookupTable, const int matchScore, const int mismatchScore) {

	static const unsigned int oddBitMask = 0x55555555;
	static const unsigned int evenBitMask = 0xAAAAAAAA;

	unsigned int i, j;
	char maxScore, maxScorePos, finalScore;

	// In entries where even bit = match/mismatch, odd bit = 0	-> matchMismatchBitVector is even bit compacted (used for forward extension)
	// In entries where odd bit = match/mismatch, even bit = 0	-> matchMismatchBitVector is odd bit compacted and then reversed (used for backward extension)
	// When all bits are 0, both above case can be handled without problem
	hspUngappedExtLookupTable[0].matchMismatchBitVector = 0;
	for (i=1; i<(1<<16); i++) {
		hspUngappedExtLookupTable[i].matchMismatchBitVector = 0;
		if ((i & oddBitMask) == 0) {
			for (j=0; j<8; j++) {
				hspUngappedExtLookupTable[i].matchMismatchBitVector |= ((i >> (15-j*2)) & 1) << (7-j);
			}
		}
		if ((i & evenBitMask) == 0) {
			for (j=0; j<8; j++) {
				hspUngappedExtLookupTable[i].matchMismatchBitVector |= ((i >> (14-j*2)) & 1) << j;
			}
		}
		// entries where both odd and even bits contain 1 are not used
	}

	// For each 16 bit match/mismatch vector, calculate max score, position giving max score, and final score
	// The first match/mismatch is in the leftmost bit of the 16 bit vector
	// 0 is a match and 1 is a mismatch (from result of xor)
	for (i=0; i<(1<<16); i++) {

		maxScore = 0;
		maxScorePos = 0;	// 0 -> before the 1st character; 16 -> after the 16th character
		finalScore = 0;

		for (j=0; j<16; j++) {
			if (((i >> (15-j)) & 1) == 0) {
				// match
				finalScore = finalScore + (char)matchScore;
				// favor longer hit if scores are the same
				if (finalScore >= maxScore) {
					maxScore = finalScore;
					maxScorePos = (char)(j+1);
				}
			} else {
				// mismatch
				finalScore = finalScore + (char)mismatchScore;
			}
		}

		hspUngappedExtLookupTable[i].maxScore = maxScore;
		hspUngappedExtLookupTable[i].maxScorePos = maxScorePos;
		hspUngappedExtLookupTable[i].finalScore = finalScore;
	}

}

void HSPFreeUngappedExtLookupTable(MMPool *mmPool, HSPUngappedExtLookupTable* hspUngappedExtLookupTable) {

	MMPoolReturn(mmPool, hspUngappedExtLookupTable, (1 << 16) * sizeof(HSPUngappedExtLookupTable));

}

Histogram* HSPAllocateHistogram(double minEvalue, double maxEvalue) {

	Histogram *histogram;

	histogram = MMUnitAllocate(sizeof(Histogram));
	histogram->minEvalue = minEvalue;
	histogram->maxEvalue = maxEvalue;
	histogram->histogramSize = (int)(ceil(log10(maxEvalue)) - floor(log10(minEvalue))) + 1;
	histogram->count = MMUnitAllocate(histogram->histogramSize * sizeof(int));

	return histogram;

}

void HSPInitializeHistogram(Histogram* histogram) {

	int i;

	for (i=0; i<histogram->histogramSize; i++) {
		histogram->count[i] = 0;
	}

}

void HSPFreeHistogram(Histogram* histogram) {

	MMUnitFree(histogram->count, histogram->histogramSize * sizeof(int));
	MMUnitFree(histogram, histogram->histogramSize * sizeof(int));

}

void HSPCountEvalueToHistogram(Histogram* histogram, const GappedHitListWithAlignment* gappedHitList, const int numberOfHit, const int roundEvalue) {

	int i;
	int evalueIndex;
	double tempEvalue;

	for (i=0; i<numberOfHit; i++) {

		if (roundEvalue == TRUE) {
			tempEvalue = HSPCalculateAndPrintEValue(NULL, gappedHitList[i].score);
		} else {
			tempEvalue = stat_gapCalcEvalue(stat_gapNominal2normalized(gappedHitList[i].score));
		}
		if (tempEvalue == 0.0) {
			tempEvalue = histogram->minEvalue;
		}
		evalueIndex = (int)ceil(log10(tempEvalue) - floor(log10(histogram->minEvalue)));
		if (evalueIndex < 0) {
			evalueIndex = 0;
		}
		if (evalueIndex < histogram->histogramSize) {
			histogram->count[evalueIndex]++;
		}

	}

}
	
void HSPPrintHistogram(FILE * outputFile, Histogram* histogram) {

	int i;
	double tempEvalue;

	Socketfprintf(outputFile,"\nCount histogram of gapped extension by expectation value:\n");

	tempEvalue = pow(10, (int)floor(log10(histogram->minEvalue)));
	Socketfprintf(outputFile,"(     0,%6.2g] : %d\n", tempEvalue, histogram->count[0]);

	for (i=1; i<histogram->histogramSize; i++) {
		tempEvalue = pow(10, i - 1 + (int)log10(histogram->minEvalue));
		Socketfprintf(outputFile,"(%6.2g,%6.2g] : %d\n", tempEvalue, tempEvalue*10, histogram->count[i]);
	}

}

void HSPPrintHeader(FILE *outputFile, const int outputFormat, 
					const char* databaseName, const int numOfSeq, const unsigned int databaseLength) {

	switch (outputFormat) {

	case OUTPUT_PAIRWISE :

		fprintf(outputFile, "Database         : %s\n", databaseName);
		fprintf(outputFile, "No. of sequences : %d\n", numOfSeq);
		fprintf(outputFile, "No. of letters   : %u\n", databaseLength);
		fprintf(outputFile, "\n");

		break;

	case OUTPUT_TABULAR :

		break;

	case  OUTPUT_TABULAR_COMMENT :

		fprintf(outputFile, "# Database         : %s\n", databaseName);
		fprintf(outputFile, "# No. of sequences : %d\n", numOfSeq);
		fprintf(outputFile, "# No. of letters   : %u\n", databaseLength);
		fprintf(outputFile, "#\n");

		break;

	}

}

void HSPPrintTrailer(FILE *outputFile, const int outputFormat, 
					 const char* databaseName, const int numOfSeq, const unsigned int databaseLength) {

	switch (outputFormat) {

	case OUTPUT_PAIRWISE :

		fprintf(outputFile, "Database         : %s\n", databaseName);
		fprintf(outputFile, "No. of sequences : %d\n", numOfSeq);
		fprintf(outputFile, "No. of letters   : %u\n", databaseLength);
		fprintf(outputFile, "\n\n");
	
		printHSPstatistic(outputFile);
	
		break;

	case OUTPUT_TABULAR :

		break;

	case OUTPUT_TABULAR_COMMENT :

		break;

	}

}

void HSPPrintQueryHeader(FILE *outputFile, const int outputFormat,
					   const char* queryPatternName, const int queryPatternLength,
					   const int *dbOrder, const int *dbScore, const Annotation *annotation, const int numOfSeq) {

	int i;
	int charPrinted;
	char stringFormat[8] = "%00.00s";

	switch (outputFormat) {

	case OUTPUT_PAIRWISE :

		fprintf(outputFile, "Query            : %s\n", queryPatternName);
		fprintf(outputFile, "No. of letters   : %u\n", queryPatternLength);

		fprintf(outputFile, "\n\n");
		fprintf(outputFile, "                                                                 Score    E\n");
		fprintf(outputFile, "Sequences producing significant alignments:                      (bits) Value\n");
		fprintf(outputFile, "\n");

		for (i=0; i<numOfSeq && dbScore[dbOrder[i]] > 0; i++) {
			charPrinted = 0;
			if (annotation[dbOrder[i]].gi > 0) {
				charPrinted += fprintf(outputFile, "gi|%d|", annotation[dbOrder[i]].gi);
			}
			stringFormat[1] = stringFormat[4] = (char)('0' + (64 - charPrinted) / 10);
			stringFormat[2] = stringFormat[5] = (char)('0' + (64 - charPrinted) % 10);
			fprintf(outputFile, stringFormat, annotation[dbOrder[i]].text);
			fprintf(outputFile, " ");
			HSPCalculateAndPrintBitScore(outputFile, dbScore[dbOrder[i]]);
			fprintf(outputFile, "   ");
			HSPCalculateAndPrintEValue(outputFile, dbScore[dbOrder[i]]);
			fprintf(outputFile, "\n");
		}
		fprintf(outputFile, "\n\n");
	
		break;

	case OUTPUT_TABULAR :

		break;

	case OUTPUT_TABULAR_COMMENT :

		fprintf(outputFile, "# Query            : %s\n", queryPatternName);
		fprintf(outputFile, "# Query id, Subject id, %% identity, alignment length, mismatches, gap openings, q. start, q. end, s. start, s. end, e-value, bit score\n");

	}

}

void HSPPrintAlignment(MMPool *mmPool, FILE *outputFile, const GappedHitListWithAlignment* gappedHitList, const int numOfHit, 
					   int outputFormat, const ContextInfo* contextInfo, 
					   const unsigned char *charMap, const unsigned char *complementMap,
					   const char* queryPatternName, const unsigned char* convertedQueryPattern, const int queryPatternLength,
					   const int * dbOrder, const SeqOffset *seqOffset, const Annotation *annotation) {

	int i, j;
	int lastDbSeqIndex, dbSeqIndex, contextNum;
	char* __restrict tempDbChar;
	char* __restrict tempQueryChar;

	char queryName[MAX_SEQ_NAME_LENGTH+1];
	char seqName[MAX_SEQ_NAME_LENGTH+1];

	int matchMatrix[16][16];

	unsigned int queryIndex, dbIndex;
	int totalLength, auxiliaryTextIndex;
	unsigned int c;
	int a, lastA;

	int alignmentCharPrinted;
	int textPosNumOfChar;
	char positionFormat[6] = "%-00u";

	char lowercase[256];

	int numOfMatch, numOfMismatch, numOfSpace, numOfGap;

	// Allocate memory
	tempDbChar = MMPoolDispatch(mmPool, MAX_ALIGNMENT_LENGTH + 1);
	tempQueryChar = MMPoolDispatch(mmPool, MAX_ALIGNMENT_LENGTH + 1);

	// duplicate query name
	strncpy(queryName, queryPatternName, MAX_SEQ_NAME_LENGTH+1);

	// take the first white space as end of pattern name
	for (i=0; i<MAX_SEQ_NAME_LENGTH; i++) {
		if (queryName[i] == ' ' || queryName[i] == '\0' || queryName[i] == '\t' || queryName[i] == '\n') {
			queryName[i] = '\0';
			break;
		}
	}

	HSPFillScoringMatrix(matchMatrix, 1, -1, 0);

	// initialize lowercase
	for (i=0; i<256; i++) {
		lowercase[i] = (char)i;
	}
	for (i=0; i<26; i++) {
		lowercase['A' + i] = 'a' + (char)i;
	}

	// Process all hits
	lastDbSeqIndex = -1;
	for (i=0; i<numOfHit; i++) {

		if (gappedHitList[i].lengthQuery > MAX_ALIGNMENT_LENGTH || gappedHitList[i].lengthText > MAX_ALIGNMENT_LENGTH) {
			fprintf(stderr, "HSPPrintAlignment(): Alignment length > maximum!\n");
			exit(1);
		}

		dbSeqIndex = dbOrder[gappedHitList[i].dbSeqIndex & CONTEXT_MASK];
		contextNum = gappedHitList[i].dbSeqIndex >> (BITS_IN_WORD - CONTEXT_BIT);

		if (dbSeqIndex != lastDbSeqIndex) {

			switch(outputFormat) {
			case OUTPUT_PAIRWISE :

				// Print header for sequence
				if (annotation[dbSeqIndex].gi > 0) {
					fprintf(outputFile, ">gi|%d|%s", annotation[dbSeqIndex].gi, annotation[dbSeqIndex].text);
				} else {
					fprintf(outputFile, ">%s", annotation[dbSeqIndex].text);
				}
				fprintf(outputFile, "No. of letters   : %u\n", seqOffset[dbSeqIndex].endPos - seqOffset[dbSeqIndex].startPos + 1);
				fprintf(outputFile, "\n");

				break;

			case OUTPUT_TABULAR :
			case OUTPUT_TABULAR_COMMENT :

				// duplicate sequence name
				strncpy(seqName, annotation[dbSeqIndex].text, MAX_SEQ_NAME_LENGTH+1);
				// take the first white space as end of pattern name
				for (j=0; j<MAX_SEQ_NAME_LENGTH; j++) {
					if (seqName[j] == ' ' || seqName[j] == '\0' || seqName[j] == '\t' || seqName[j] == '\n') {
						seqName[j] = '\0';
						break;
					}
				}

				break;

			}

		}
		lastDbSeqIndex = dbSeqIndex;

		// Decode alignment and DB characters

		queryIndex = 0;
		totalLength = 0;
		auxiliaryTextIndex = 0;
		numOfMatch = 0;
		numOfMismatch = 0;
		numOfSpace = 0;
		numOfGap = 0;
		lastA = -1;
		c = 0;

		while (queryIndex < gappedHitList[i].lengthQuery) {
			if (totalLength % ALIGN_PER_WORD == 0) {
				c = gappedHitList[i].alignment[totalLength / ALIGN_PER_WORD];
			}
			a = (char)(c >> (BITS_IN_WORD - ALIGN_BIT));
			c <<= ALIGN_BIT;
			switch(a) {
			case ALIGN_MATCH :
				if (contextInfo[contextNum].reversed) {
					tempQueryChar[totalLength] = dnaChar[convertedQueryPattern[queryPatternLength - 1 - queryIndex - gappedHitList[i].posQuery]];
				} else {
					tempQueryChar[totalLength] = dnaChar[convertedQueryPattern[queryIndex + gappedHitList[i].posQuery]];
				}
				tempDbChar[totalLength] = tempQueryChar[totalLength];
				numOfMatch++;
				queryIndex++;
				break;
			case ALIGN_MISMATCH_AMBIGUITY :
				tempDbChar[totalLength] = dnaChar[(gappedHitList[i].auxiliaryText[auxiliaryTextIndex / AUX_TEXT_PER_WORD] << ((auxiliaryTextIndex % AUX_TEXT_PER_WORD) * AUX_TEXT_BIT) >> (BITS_IN_WORD - AUX_TEXT_BIT))];
				if (contextInfo[contextNum].reversed) {
					tempDbChar[totalLength] = complementMap[tempDbChar[totalLength]];
					tempQueryChar[totalLength] = dnaChar[convertedQueryPattern[queryPatternLength - 1 - queryIndex - gappedHitList[i].posQuery]];
				} else {
					tempQueryChar[totalLength] = dnaChar[convertedQueryPattern[queryIndex + gappedHitList[i].posQuery]];
				}
				if (matchMatrix[charMap[tempDbChar[totalLength]]][charMap[tempQueryChar[totalLength]]] <= 0) {
					// Mismatch
					numOfMismatch++;
				} else {
					// Match with ambiguity
					if (tempDbChar[totalLength] == tempQueryChar[totalLength] && charMap[tempDbChar[totalLength]] < ALPHABET_SIZE) {
						fprintf(stderr, "HSPPrintAlignment(): Alignment error!\n");
						exit(1);
					}
					numOfMatch++;
				}
				queryIndex++;
				auxiliaryTextIndex++;
				break;
			case ALIGN_INSERT :
				tempDbChar[totalLength] = dnaChar[(gappedHitList[i].auxiliaryText[auxiliaryTextIndex / AUX_TEXT_PER_WORD] << ((auxiliaryTextIndex % AUX_TEXT_PER_WORD) * AUX_TEXT_BIT) >> (BITS_IN_WORD - AUX_TEXT_BIT))];
				if (contextInfo[contextNum].reversed) {
					tempDbChar[totalLength] = complementMap[tempDbChar[totalLength]];
				}
				tempQueryChar[totalLength] = '-';
				numOfSpace++;
				if (lastA != ALIGN_INSERT) {
					numOfGap++;
				}
				auxiliaryTextIndex++;
				break;
			case ALIGN_DELETE :
				tempDbChar[totalLength] = '-';
				if (contextInfo[contextNum].reversed) {
					tempQueryChar[totalLength] = dnaChar[convertedQueryPattern[queryPatternLength - 1 - queryIndex - gappedHitList[i].posQuery]];
				} else {
					tempQueryChar[totalLength] = dnaChar[convertedQueryPattern[queryIndex + gappedHitList[i].posQuery]];
				}
				numOfSpace++;
				if (lastA != ALIGN_DELETE) {
					numOfGap++;
				}
				queryIndex++;
				break;
			}
			lastA = a;
			totalLength++;
		}

		// Print alignment

		switch (outputFormat) {

		case OUTPUT_PAIRWISE :

			fprintf(outputFile, "Score = ");
			HSPCalculateAndPrintBitScore(outputFile, gappedHitList[i].score);
			fprintf(outputFile, " bits (%d), Expect = ", gappedHitList[i].score);
			HSPCalculateAndPrintEValue(outputFile, gappedHitList[i].score);
			fprintf(outputFile, "\n");
			fprintf(outputFile, "Identities = %d/%d (%.2f%%), Gaps = %d/%d (%.2f%%)\n", numOfMatch, totalLength, (float)numOfMatch / totalLength * 100, numOfSpace, totalLength, (float)numOfSpace / totalLength * 100);
			if (contextInfo[contextNum].reversed) {
				fprintf(outputFile, "Strand = Plus / Minus\n");
			} else {
				fprintf(outputFile, "Strand = Plus / Plus\n");
			}
			fprintf(outputFile, "\n\n");

			// Determine position length and format
			c = gappedHitList[i].posText + gappedHitList[i].lengthText - seqOffset[dbSeqIndex].startPos;
			textPosNumOfChar = 0;
			while (c > 0) {
				c /= 10;
				textPosNumOfChar++;
			}
			if (textPosNumOfChar >= 10) {
				positionFormat[2] = '0' + (char)(textPosNumOfChar / 10);
				positionFormat[3] = '0' + (char)(textPosNumOfChar % 10);
				positionFormat[4] = 'u';
				positionFormat[5] = '\0';
			} else {
				positionFormat[2] = '0' + (char)textPosNumOfChar;
				positionFormat[3] = 'u';
				positionFormat[4] = '\0';
			}

			alignmentCharPrinted = 0;
			if (contextInfo[contextNum].reversed) {
				queryIndex = queryPatternLength + 1 - gappedHitList[i].posQuery - gappedHitList[i].lengthQuery;
				dbIndex = gappedHitList[i].posText + gappedHitList[i].lengthText - seqOffset[dbSeqIndex].startPos;
			} else {
				queryIndex = gappedHitList[i].posQuery + 1;
				dbIndex = gappedHitList[i].posText - seqOffset[dbSeqIndex].startPos + 1;
			}

			while (alignmentCharPrinted < totalLength) {

				if (contextInfo[contextNum].reversed) {

					// Print query line
					fprintf(outputFile, "Query: ");
					fprintf(outputFile, positionFormat, queryIndex);
					fprintf(outputFile, " ");

					j = 0;
					while (j < 60 && alignmentCharPrinted + j < totalLength) {
						fprintf(outputFile, "%c", lowercase[tempQueryChar[totalLength - 1 - alignmentCharPrinted - j]]);
						if (tempQueryChar[totalLength - 1 - alignmentCharPrinted - j] != '-') {
							queryIndex++;
						}
						j++;
					}

					fprintf(outputFile, " ");
					fprintf(outputFile, "%u", queryIndex - 1);
					fprintf(outputFile, "\n");

					// Print alignment line
					for (j=0; j<textPosNumOfChar + 8; j++) {
						fprintf(outputFile, " ");
					}

					j = 0;
					while (j < 60 && alignmentCharPrinted + j < totalLength) {
						if (tempQueryChar[totalLength - 1 - alignmentCharPrinted - j] == tempDbChar[totalLength - 1 - alignmentCharPrinted - j]) {
							fprintf(outputFile, "|");
						} else {
							fprintf(outputFile, " ");
						}
						j++;
					}
					fprintf(outputFile, "\n");

					// Print database line
					fprintf(outputFile, "Sbjct: ");
					fprintf(outputFile, positionFormat, dbIndex);
					fprintf(outputFile, " ");

					j = 0;
					while (j < 60 && alignmentCharPrinted + j < totalLength) {
						fprintf(outputFile, "%c", lowercase[tempDbChar[totalLength - 1 - alignmentCharPrinted - j]]);
						if (tempDbChar[totalLength - 1 - alignmentCharPrinted - j] != '-') {
							dbIndex--;
						}
						j++;
					}

					fprintf(outputFile, " ");
					fprintf(outputFile, positionFormat, dbIndex + 1);
					fprintf(outputFile, "\n");

				} else {

					// Print query line
					fprintf(outputFile, "Query: ");
					fprintf(outputFile, positionFormat, queryIndex);
					fprintf(outputFile, " ");

					j = 0;
					while (j < 60 && alignmentCharPrinted + j < totalLength) {
						fprintf(outputFile, "%c", lowercase[tempQueryChar[alignmentCharPrinted + j]]);
						if (tempQueryChar[alignmentCharPrinted + j] != '-') {
							queryIndex++;
						}
						j++;
					}

					fprintf(outputFile, " ");
					fprintf(outputFile, "%u", queryIndex - 1);
					fprintf(outputFile, "\n");

					// Print alignment line
					for (j=0; j<textPosNumOfChar + 8; j++) {
						fprintf(outputFile, " ");
					}

					j = 0;
					while (j < 60 && alignmentCharPrinted + j < totalLength) {
						if (matchMatrix[charMap[tempQueryChar[alignmentCharPrinted + j]]][charMap[tempDbChar[alignmentCharPrinted + j]]] > 0) {
							fprintf(outputFile, "|");
						} else {
							fprintf(outputFile, " ");
						}
						j++;
					}
					fprintf(outputFile, "\n");

					// Print database line
					fprintf(outputFile, "Sbjct: ");
					fprintf(outputFile, positionFormat, dbIndex);
					fprintf(outputFile, " ");

					j = 0;
					while (j < 60 && alignmentCharPrinted + j < totalLength) {
						fprintf(outputFile, "%c", lowercase[tempDbChar[alignmentCharPrinted + j]]);
						if (tempDbChar[alignmentCharPrinted + j] != '-') {
							dbIndex++;
						}
						j++;
					}

					fprintf(outputFile, " ");
					fprintf(outputFile, positionFormat, dbIndex - 1);
					fprintf(outputFile, "\n");

				}
					
				fprintf(outputFile, "\n\n");
				alignmentCharPrinted += j;

			}

			break;

		case OUTPUT_TABULAR :
		case OUTPUT_TABULAR_COMMENT :
			// print query name
			fprintf(outputFile, "%s\t", queryName);

			// print gi
			if (annotation[dbSeqIndex].gi > 0) {
				fprintf(outputFile, "gi|%d|", annotation[dbSeqIndex].gi);
			}

			// print db sequence name
			fprintf(outputFile, "%s\t", seqName);

			// print percentage identity
			fprintf(outputFile, "%.2f\t", (float)(totalLength - numOfMismatch - numOfSpace) / (float)(totalLength) * 100);

			// print align length
			fprintf(outputFile, "%d\t", totalLength);

			// print number of mismatch
			fprintf(outputFile, "%d\t", numOfMismatch);

			// print number of gap opening
			fprintf(outputFile, "%d\t", numOfGap);

			// print query start and end
			if (contextInfo[contextNum].reversed) {
				fprintf(outputFile, "%u\t%u\t", queryPatternLength - gappedHitList[i].posQuery - gappedHitList[i].lengthQuery + 1,
											   queryPatternLength - gappedHitList[i].posQuery);
			} else {
				fprintf(outputFile, "%u\t%u\t", gappedHitList[i].posQuery + 1,
											   gappedHitList[i].posQuery + gappedHitList[i].lengthQuery);
			}

			// print text start and end
			if (contextInfo[contextNum].reversed) {
				fprintf(outputFile, "%u\t%u\t", gappedHitList[i].posText + gappedHitList[i].lengthText - seqOffset[dbSeqIndex].startPos,
											   gappedHitList[i].posText + 1 - seqOffset[dbSeqIndex].startPos);
			} else {
				fprintf(outputFile, "%u\t%u\t", gappedHitList[i].posText + 1 - seqOffset[dbSeqIndex].startPos,
											   gappedHitList[i].posText + gappedHitList[i].lengthText - seqOffset[dbSeqIndex].startPos);
			}


			// print evalue
			HSPCalculateAndPrintEValue(outputFile, gappedHitList[i].score);

			// print bit score
			HSPCalculateAndPrintBitScore(outputFile, gappedHitList[i].score);

			fprintf(outputFile, "\n");

			break;

		default :

			fprintf(stderr, "HSPPrintAlignment() : Output format is invalid!\n");
			exit(1);

		}

	}

	// Free memory
	MMPoolReturn(mmPool, tempDbChar, MAX_ALIGNMENT_LENGTH + 1);
	MMPoolReturn(mmPool, tempQueryChar, MAX_ALIGNMENT_LENGTH + 1);

}

void HSPPrintNoAlignment(FILE *outputFile, int outputFormat) {

	switch (outputFormat) {

	case OUTPUT_PAIRWISE :

		break;

	case OUTPUT_TABULAR :
	case OUTPUT_TABULAR_COMMENT :

		fprintf(outputFile, "\n");

	}

}

double HSPCalculateAndPrintBitScore(FILE *outputFile, const int score) {

	double bitScore;

	// print bit score
	bitScore = stat_gapNominal2normalized(score);

	if (outputFile != NULL) {
		if (bitScore > 9999) {
			fprintf(outputFile, "%4.3le", bitScore);
		} else if (bitScore > 99.9) {
			fprintf(outputFile, "%4.0ld", (int)bitScore);
		} else {
			fprintf(outputFile, "%4.1lf", bitScore);
		}
	}

	return bitScore;

}

double HSPCalculateAndPrintEValue(FILE *outputFile, const int score) {

	double evalue;
	char evalueString[10];

	evalue = stat_gapCalcEvalue(stat_gapNominal2normalized(score));

	if (evalue < 1.0e-180) {
		sprintf(evalueString, "0.0\t");
	} else if (evalue < 1.0e-99) {
		sprintf(evalueString, "%2.0le\t", evalue);
	} else if (evalue < 0.0009) {
		sprintf(evalueString, "%3.0le\t", evalue);
	} else if (evalue < 0.1) {
		sprintf(evalueString, "%4.3lf\t", evalue);
	} else if (evalue < 1.0) { 
		sprintf(evalueString, "%3.2lf\t", evalue);
	} else if (evalue < 10.0) {
		sprintf(evalueString, "%2.1lf\t", evalue);
	} else { 
		sprintf(evalueString, "%5.0lf\t", evalue);
	}

	if (outputFile != NULL) {
		fprintf(outputFile, "%s\t", evalueString);
	}

	sscanf(evalueString, "%g", &evalue);
	return evalue;

}

int HitListPosTextQueryOrder(const void *hitList, const int index1, const int index2) {

	if (((HitListWithPosQuery*)hitList + index1)->posText != ((HitListWithPosQuery*)hitList + index2)->posText) {
		if (((HitListWithPosQuery*)hitList + index1)->posText > ((HitListWithPosQuery*)hitList + index2)->posText) {
			return 1;
		} else {
			return -1;
		}
	} else {
		return ((HitListWithPosQuery*)hitList + index1)->posQuery - ((HitListWithPosQuery*)hitList + index2)->posQuery;
	}

}

int GappedHitListPosTextQueryLengthScoreOrder(const void *gappedHitList, const int index1, const int index2) {

	if (((GappedHitList*)gappedHitList + index1)->posText != ((GappedHitList*)gappedHitList + index2)->posText) {
		if (((GappedHitList*)gappedHitList + index1)->posText > ((GappedHitList*)gappedHitList + index2)->posText) {
			return 1;
		} else {
			return -1;
		}
	} else {
		if (((GappedHitList*)gappedHitList + index1)->posQuery != ((GappedHitList*)gappedHitList + index2)->posQuery) {
			if (((GappedHitList*)gappedHitList + index1)->posQuery > ((GappedHitList*)gappedHitList + index2)->posQuery) {
				return 1;
			} else {
				return -1;
			}
		} else {
			if (((GappedHitList*)gappedHitList + index1)->lengthText != ((GappedHitList*)gappedHitList + index2)->lengthText) {
				return ((GappedHitList*)gappedHitList + index1)->lengthText - ((GappedHitList*)gappedHitList + index2)->lengthText;
			} else {
				if (((GappedHitList*)gappedHitList + index1)->lengthQuery != ((GappedHitList*)gappedHitList + index2)->lengthQuery) {
					return ((GappedHitList*)gappedHitList + index1)->lengthQuery - ((GappedHitList*)gappedHitList + index2)->lengthQuery;
				} else {
					// higher score first
					return ((GappedHitList*)gappedHitList + index2)->score - ((GappedHitList*)gappedHitList + index1)->score;
				}
			}
		}
	}

}

int GappedHitListPosTextQueryScoreLengthOrder(const void *gappedHitList, const int index1, const int index2) {

	if (((GappedHitList*)gappedHitList + index1)->posText != ((GappedHitList*)gappedHitList + index2)->posText) {
		if (((GappedHitList*)gappedHitList + index1)->posText > ((GappedHitList*)gappedHitList + index2)->posText) {
			return 1;
		} else {
			return -1;
		}
	} else {
		if (((GappedHitList*)gappedHitList + index1)->posQuery != ((GappedHitList*)gappedHitList + index2)->posQuery) {
			if (((GappedHitList*)gappedHitList + index1)->posQuery > ((GappedHitList*)gappedHitList + index2)->posQuery) {
				return 1;
			} else {
				return -1;
			}
		} else {
			if (((GappedHitList*)gappedHitList + index1)->score != ((GappedHitList*)gappedHitList + index2)->score) {
				// higher score first
				return ((GappedHitList*)gappedHitList + index2)->score - ((GappedHitList*)gappedHitList + index1)->score;
			} else {
				// short query length first
				if (((GappedHitList*)gappedHitList + index1)->lengthQuery != ((GappedHitList*)gappedHitList + index2)->lengthQuery) {
					if (((GappedHitList*)gappedHitList + index1)->lengthQuery > ((GappedHitListWithAlignment*)gappedHitList + index2)->lengthQuery) {
						return 1;
					} else {
						return -1;
					}
				} else {
					// short text length first
					if (((GappedHitList*)gappedHitList + index1)->lengthText != ((GappedHitListWithAlignment*)gappedHitList + index2)->lengthText) {
						if (((GappedHitList*)gappedHitList + index1)->lengthText > ((GappedHitListWithAlignment*)gappedHitList + index2)->lengthText) {
							return 1;
						} else {
							return -1;
						}
					} else {
						return 0;
					}
				}
			}
		}
	}

}

int GappedHitListEndTextQueryScoreLengthOrder(const void *gappedHitList, const int index1, const int index2) {

	if (((GappedHitList*)gappedHitList + index1)->posText + ((GappedHitList*)gappedHitList + index1)->lengthText !=
		((GappedHitList*)gappedHitList + index2)->posText + ((GappedHitList*)gappedHitList + index2)->lengthText) {
		if (((GappedHitList*)gappedHitList + index1)->posText + ((GappedHitList*)gappedHitList + index1)->lengthText >
			((GappedHitList*)gappedHitList + index2)->posText + ((GappedHitList*)gappedHitList + index2)->lengthText) {
			return 1;
		} else {
			return -1;
		}
	} else {
		if (((GappedHitList*)gappedHitList + index1)->posQuery + ((GappedHitList*)gappedHitList + index1)->lengthQuery !=
			((GappedHitList*)gappedHitList + index2)->posQuery + ((GappedHitList*)gappedHitList + index2)->lengthQuery) {
			if (((GappedHitList*)gappedHitList + index1)->posQuery + ((GappedHitList*)gappedHitList + index1)->lengthQuery >
				((GappedHitList*)gappedHitList + index2)->posQuery + ((GappedHitList*)gappedHitList + index2)->lengthQuery) {
				return 1;
			} else {
				return -1;
			}
		} else {
			if (((GappedHitList*)gappedHitList + index1)->score != ((GappedHitList*)gappedHitList + index2)->score) {
				// higher score first
				return ((GappedHitList*)gappedHitList + index2)->score - ((GappedHitList*)gappedHitList + index1)->score;
			} else {
				// short query length first
				if (((GappedHitList*)gappedHitList + index1)->lengthQuery != ((GappedHitListWithAlignment*)gappedHitList + index2)->lengthQuery) {
					if (((GappedHitList*)gappedHitList + index1)->lengthQuery > ((GappedHitListWithAlignment*)gappedHitList + index2)->lengthQuery) {
						return 1;
					} else {
						return -1;
					}
				} else {
					// short text length first
					if (((GappedHitList*)gappedHitList + index1)->lengthText != ((GappedHitListWithAlignment*)gappedHitList + index2)->lengthText) {
						if (((GappedHitList*)gappedHitList + index1)->lengthText > ((GappedHitListWithAlignment*)gappedHitList + index2)->lengthText) {
							return 1;
						} else {
							return -1;
						}
					} else {
						return 0;
					}
				}
			}
		}
	}

}

int GappedHitListEnclosedScoreOrder(const void *gappedHitList, const int index1, const int index2) {

	if (((GappedHitList*)gappedHitList + index1)->posText != ((GappedHitList*)gappedHitList + index2)->posText) {
		if (((GappedHitList*)gappedHitList + index1)->posText >	((GappedHitList*)gappedHitList + index2)->posText) {
			return 1;
		} else {
			return -1;
		}
	} else {
		if (((GappedHitList*)gappedHitList + index1)->lengthText != ((GappedHitList*)gappedHitList + index2)->lengthText) {
			if (((GappedHitList*)gappedHitList + index1)->lengthText < ((GappedHitList*)gappedHitList + index2)->lengthText) {
				return 1;
			} else {
				return -1;
			}
		} else {
			if (((GappedHitList*)gappedHitList + index1)->posQuery != ((GappedHitList*)gappedHitList + index2)->posQuery) {
				if (((GappedHitList*)gappedHitList + index1)->posQuery > ((GappedHitList*)gappedHitList + index2)->posQuery) {
					return 1;
				} else {
					return -1;
				}
			} else {
				if (((GappedHitList*)gappedHitList + index1)->lengthQuery != ((GappedHitList*)gappedHitList + index2)->lengthQuery) {
					if (((GappedHitList*)gappedHitList + index1)->lengthQuery < ((GappedHitList*)gappedHitList + index2)->lengthQuery) {
						return 1;
					} else {
						return -1;
					}
				} else {
					// higher score first
					return ((GappedHitList*)gappedHitList + index2)->score - ((GappedHitList*)gappedHitList + index1)->score;
				}
			}
		}
	}

}

#endif
