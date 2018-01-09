#include "config.h"
#ifdef MMX
/*

   BWTConstruct.c		BWT-Index Construction

   This module constructs BWT and auxiliary data structures.

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
#include <emmintrin.h>
#include <mmintrin.h>
#include "BWTConstruct.h"
#include "MiscUtilities.h"
#include "DNACount.h"
#include "TextConverter.h"
#include "MemManager.h"
#include "QSufSort.h"
#include "r250.h"

// Static functions
static void BWTIncConstruct(BWTInc *bwtInc, const unsigned int numChar);

static void BWTIncSetBuildSizeAndTextAddr(BWTInc *bwtInc);
static void BWTIncPutPackedTextToRank(const unsigned int *packedText, unsigned int* __restrict rank, unsigned int* __restrict cumulativeCount, const unsigned int numChar);
static void BWTIncBuildPackedBwt(const unsigned int *relativeRank, unsigned int* __restrict bwt, const unsigned int numChar,
								 const unsigned int *cumulativeCount, const unsigned int *packedShift);
static unsigned int BWTIncGetAbsoluteRank(BWT *bwt, unsigned int* __restrict absoluteRank, unsigned int* __restrict seq, const unsigned int *packedText, 
								  const unsigned int numChar, const unsigned int* cumulativeCount, const unsigned int firstCharInLastIteration);
static void BWTIncSortKey(unsigned int* __restrict key, unsigned int* __restrict seq, const unsigned int numItem);
static void BWTIncBuildRelativeRank(unsigned int* __restrict sortedRank, unsigned int* __restrict seq, unsigned int* __restrict relativeRank, const unsigned int numItem, 
									unsigned int oldInverseSa0, const unsigned int *cumulativeCount);
static void BWTIncBuildBwt(unsigned int*  seq, const unsigned int *relativeRank, const unsigned int numChar, const unsigned int *cumulativeCount);	// seq is replaced with Bwt
static void BWTIncMergeBwt(const unsigned int *sortedRank, const unsigned int* oldBwt, const unsigned int *insertBwt, unsigned int* __restrict mergedBwt, 
						   const unsigned int numOldBwt, const unsigned int numInsertBwt);


BWTInc *BWTIncCreate(MMPool *mmPool, const unsigned int textLength, const float targetNBit,
					 const unsigned int initialMaxBuildSize, const unsigned int incMaxBuildSize) {

	BWTInc *bwtInc;
	unsigned int i;

	if (targetNBit == 0) {
		fprintf(stderr, "BWTIncCreate() : targetNBit = 0!\n");
		exit(1);
	}
	
	bwtInc = MMPoolDispatch(mmPool, sizeof(BWTInc));
	
	bwtInc->numberOfIterationDone = 0;

	bwtInc->bwt = BWTCreate(mmPool, textLength, NULL);

	bwtInc->initialMaxBuildSize = initialMaxBuildSize;
	bwtInc->incMaxBuildSize = incMaxBuildSize;

	bwtInc->targetNBit = targetNBit;

	bwtInc->cumulativeCountInCurrentBuild = MMPoolDispatch(mmPool, sizeof(unsigned int) * (ALPHABET_SIZE + 1));
	initializeVAL(bwtInc->cumulativeCountInCurrentBuild, ALPHABET_SIZE + 1, 0);

	// Build frequently accessed data
	bwtInc->packedShift = MMPoolDispatch(mmPool, sizeof(unsigned int) * CHAR_PER_WORD);
	for (i=0; i<CHAR_PER_WORD; i++) {
		bwtInc->packedShift[i] = BITS_IN_WORD - (i+1) * BIT_PER_CHAR;
	}

	bwtInc->targetTextLength = textLength;
	bwtInc->availableWord = (unsigned int)((textLength + OCC_INTERVAL - 1) / OCC_INTERVAL * OCC_INTERVAL / BITS_IN_WORD * bwtInc->targetNBit);
	if (bwtInc->availableWord < BWTResidentSizeInWord(textLength) + WORD_BETWEEN_OCC / 2 + BWTOccValueMinorSizeInWord(textLength)) {	// + 8 words so that the 128 bits before and after an explicit occ are in the same aligned 64 byte
		fprintf(stderr, "BWTIncCreate() : targetNBit is too low!\n");
		exit(1);
	}
	if (bwtInc->availableWord < BWTINC_MIN_MEMORY_IN_WORD) {
		bwtInc->availableWord = BWTINC_MIN_MEMORY_IN_WORD;
	}
	bwtInc->workingMemory = MMUnitAllocate(bwtInc->availableWord * BYTES_IN_WORD);

	return bwtInc;

}

void BWTIncFree(MMPool *mmPool, BWTInc *bwtInc) {

	MMUnitFree(bwtInc->workingMemory, bwtInc->availableWord * BYTES_IN_WORD);
	MMPoolReturn(mmPool, bwtInc, sizeof(BWTInc));

}

BWTInc *BWTIncConstructFromPacked(MMPool *mmPool, const char *inputFileName, const unsigned int showProgress,
								  const float targetNBit, const unsigned int initialMaxBuildSize, const unsigned int incMaxBuildSize) {

	FILE *packedFile;
	unsigned int packedFileLen;
	unsigned int totalTextLength;
	unsigned int textToLoad, textSizeInByte;
	unsigned int processedTextLength;
	unsigned char lastByteLength;

	BWTInc *bwtInc;

	packedFile = (FILE*)fopen64(inputFileName, "rb");

	if (packedFile == NULL) {
		fprintf(stderr, "BWTIncConstructFromPacked() : Cannot open inputFileName!\n");
		exit(1);
	}

	fseek(packedFile, -1, SEEK_END);
	packedFileLen = ftell(packedFile);
	if ((int)packedFileLen < 0) {
		fprintf(stderr, "BWTIncConstructFromPacked: Cannot determine file length!\n");
		exit(1);
	}
	fread(&lastByteLength, sizeof(unsigned char), 1, packedFile);
	totalTextLength = TextLengthFromBytePacked(packedFileLen, BIT_PER_CHAR, lastByteLength);

	bwtInc = BWTIncCreate(mmPool, totalTextLength, targetNBit, initialMaxBuildSize, incMaxBuildSize);

	BWTIncSetBuildSizeAndTextAddr(bwtInc);

	if (bwtInc->buildSize > totalTextLength) {
		textToLoad = totalTextLength;
	} else {
		textToLoad = totalTextLength - ((totalTextLength - bwtInc->buildSize + CHAR_PER_WORD - 1) / CHAR_PER_WORD * CHAR_PER_WORD);
	}
	textSizeInByte = textToLoad / CHAR_PER_BYTE;	// excluded the odd byte

	fseek(packedFile, -2, SEEK_CUR);
	fseek(packedFile, -((int)textSizeInByte), SEEK_CUR);
	fread(bwtInc->textBuffer, sizeof(unsigned char), textSizeInByte + 1, packedFile);
	fseek(packedFile, -((int)textSizeInByte + 1), SEEK_CUR);

	ConvertBytePackedToWordPacked(bwtInc->textBuffer, bwtInc->packedText, ALPHABET_SIZE, textToLoad);
	BWTIncConstruct(bwtInc, textToLoad);

	processedTextLength = textToLoad;

	while (processedTextLength < totalTextLength) {
		textToLoad = bwtInc->buildSize / CHAR_PER_WORD * CHAR_PER_WORD;
		if (textToLoad > totalTextLength - processedTextLength) {
			textToLoad = totalTextLength - processedTextLength;
		}
		textSizeInByte = textToLoad / CHAR_PER_BYTE;
		fseek(packedFile, -((int)textSizeInByte), SEEK_CUR);
		fread(bwtInc->textBuffer, sizeof(unsigned char), textSizeInByte, packedFile);
		fseek(packedFile, -((int)textSizeInByte), SEEK_CUR);
		ConvertBytePackedToWordPacked(bwtInc->textBuffer, bwtInc->packedText, ALPHABET_SIZE, textToLoad);
		BWTIncConstruct(bwtInc, textToLoad);
		processedTextLength += textToLoad;
		if (showProgress && bwtInc->numberOfIterationDone % 10 == 0) {
			printf("%u iterations done. %u characters processed.\n", bwtInc->numberOfIterationDone, processedTextLength);
		}
	}

	return bwtInc;

}

static void BWTIncConstruct(BWTInc *bwtInc, const unsigned int numChar) {

	unsigned int i;
	unsigned int mergedBwtSizeInWord, mergedOccSizeInWord;
	unsigned int firstCharInThisIteration;

	unsigned int *relativeRank, *seq, *sortedRank, *insertBwt, *mergedBwt;
	unsigned int newInverseSa0RelativeRank, oldInverseSa0RelativeRank, newInverseSa0;

	#ifdef DEBUG
	if (numChar > bwtInc->buildSize) {
		fprintf(stderr, "BWTIncConstruct(): numChar > buildSize!\n");
		exit(1);
	}
	#endif

	mergedBwtSizeInWord = BWTResidentSizeInWord(bwtInc->bwt->textLength + numChar);
	mergedOccSizeInWord = BWTOccValueMinorSizeInWord(bwtInc->bwt->textLength + numChar);

	initializeVAL(bwtInc->cumulativeCountInCurrentBuild, ALPHABET_SIZE + 1, 0);

	if (bwtInc->bwt->textLength == 0) {		// Initial build

		// Set address
		seq = bwtInc->workingMemory;
		relativeRank = seq + bwtInc->buildSize + 1;
		mergedBwt = insertBwt = bwtInc->workingMemory + bwtInc->availableWord - mergedBwtSizeInWord - WORD_BETWEEN_OCC / 2;	// build in place

		BWTIncPutPackedTextToRank(bwtInc->packedText, relativeRank, bwtInc->cumulativeCountInCurrentBuild, numChar);

		firstCharInThisIteration = relativeRank[0];
		relativeRank[numChar] = 0;

		// Sort suffix
		QSufSortSuffixSort((int*)relativeRank, (int*)seq, (int)numChar, (int)ALPHABET_SIZE - 1, 0, FALSE);
		newInverseSa0 = relativeRank[0];

		// Clear BWT area
		initializeVAL(insertBwt, mergedBwtSizeInWord, 0);

		// Build BWT
		BWTIncBuildPackedBwt(relativeRank, insertBwt, numChar, bwtInc->cumulativeCountInCurrentBuild, bwtInc->packedShift);

		// so that the cumulativeCount is not deducted
		bwtInc->firstCharInLastIteration = ALPHABET_SIZE;

	} else {		// Incremental build

#ifdef DEBUG
		if (numChar / CHAR_PER_WORD * CHAR_PER_WORD != numChar) {
			fprintf(stderr, "BWTIncConstruct(): numChar must be multiple of char per word!\n");
			exit(1);
		}
#endif

		// Set address
		sortedRank = bwtInc->workingMemory;
		seq = sortedRank + bwtInc->buildSize + 1;
		insertBwt = seq;
		relativeRank = seq + bwtInc->buildSize + 1;

		// Store the first character of this iteration
		firstCharInThisIteration = bwtInc->packedText[0] >> (BITS_IN_WORD - BIT_PER_CHAR);

		// Count occurrence of input text
		ForwardDNAAllOccCountNoLimit(bwtInc->packedText, numChar, bwtInc->cumulativeCountInCurrentBuild + 1, bwtInc->bwt->decodeTable);
		// Add the first character of the previous iteration to represent the inverseSa0 of the previous iteration
		bwtInc->cumulativeCountInCurrentBuild[bwtInc->firstCharInLastIteration + 1]++;
		bwtInc->cumulativeCountInCurrentBuild[2] += bwtInc->cumulativeCountInCurrentBuild[1];
		bwtInc->cumulativeCountInCurrentBuild[3] += bwtInc->cumulativeCountInCurrentBuild[2];
		bwtInc->cumulativeCountInCurrentBuild[4] += bwtInc->cumulativeCountInCurrentBuild[3];

		// Get rank of new suffix among processed suffix
		// The seq array is built into ALPHABET_SIZE + 2 groups; ALPHABET_SIZE groups + 1 group divided into 2 by inverseSa0 + inverseSa0 as 1 group
		oldInverseSa0RelativeRank = BWTIncGetAbsoluteRank(bwtInc->bwt, sortedRank, seq, bwtInc->packedText, 
														  numChar, bwtInc->cumulativeCountInCurrentBuild, bwtInc->firstCharInLastIteration);

		// Sort rank by ALPHABET_SIZE + 2 groups (or ALPHABET_SIZE + 1 groups when inverseSa0 sit on the border of a group)
		for (i=0; i<ALPHABET_SIZE; i++) {
			if (bwtInc->cumulativeCountInCurrentBuild[i] > oldInverseSa0RelativeRank ||
				bwtInc->cumulativeCountInCurrentBuild[i+1] <= oldInverseSa0RelativeRank) {
				BWTIncSortKey(sortedRank + bwtInc->cumulativeCountInCurrentBuild[i], seq + bwtInc->cumulativeCountInCurrentBuild[i], bwtInc->cumulativeCountInCurrentBuild[i+1] - bwtInc->cumulativeCountInCurrentBuild[i]);
			} else {
				if (bwtInc->cumulativeCountInCurrentBuild[i] < oldInverseSa0RelativeRank) {
					BWTIncSortKey(sortedRank + bwtInc->cumulativeCountInCurrentBuild[i], seq + bwtInc->cumulativeCountInCurrentBuild[i], oldInverseSa0RelativeRank - bwtInc->cumulativeCountInCurrentBuild[i]);
				}
				if (bwtInc->cumulativeCountInCurrentBuild[i+1] > oldInverseSa0RelativeRank + 1) {
					BWTIncSortKey(sortedRank + oldInverseSa0RelativeRank + 1, seq + oldInverseSa0RelativeRank + 1, bwtInc->cumulativeCountInCurrentBuild[i+1] - oldInverseSa0RelativeRank - 1);
				}
			}
		}

		// build relative rank; sortedRank is updated for merging to cater for the fact that $ is not encoded in bwt
		// the cumulative freq information is used to make sure that inverseSa0 and suffix beginning with different characters are kept in different unsorted groups)
		BWTIncBuildRelativeRank(sortedRank, seq, relativeRank, numChar, bwtInc->bwt->inverseSa0, bwtInc->cumulativeCountInCurrentBuild);
#ifdef DEBUG
		if (relativeRank[numChar] != oldInverseSa0RelativeRank) {
			fprintf(stderr, "BWTIncConstruct(): relativeRank[numChar] != oldInverseSa0RelativeRank!\n");
			exit(1);
		}
#endif

		// Sort suffix
		QSufSortSuffixSort((int*)relativeRank, (int*)seq, (int)numChar, (int)numChar, 1, TRUE);

		newInverseSa0RelativeRank = relativeRank[0];
		newInverseSa0 = sortedRank[newInverseSa0RelativeRank] + newInverseSa0RelativeRank;

		sortedRank[newInverseSa0RelativeRank] = 0;	// a special value so that this is skipped in the merged bwt

		// Build BWT
		BWTIncBuildBwt(seq, relativeRank, numChar, bwtInc->cumulativeCountInCurrentBuild);

		// Merge BWT
		mergedBwt = bwtInc->workingMemory + bwtInc->availableWord - mergedBwtSizeInWord - WORD_BETWEEN_OCC / 2
				    - bwtInc->numberOfIterationDone * OCC_INTERVAL / BIT_PER_CHAR;
					// minus numberOfIteration * occInterval to create a buffer for merging
		BWTIncMergeBwt(sortedRank, bwtInc->bwt->bwtCode, insertBwt, mergedBwt, bwtInc->bwt->textLength, numChar);

	}

	// Build auxiliary structure and update info and pointers in BWT
	bwtInc->bwt->textLength += numChar;
	bwtInc->bwt->bwtCode = mergedBwt;
	bwtInc->bwt->bwtSizeInWord = mergedBwtSizeInWord;
	bwtInc->bwt->occSizeInWord = mergedOccSizeInWord;
	if (mergedBwt < bwtInc->workingMemory + mergedOccSizeInWord) {
		fprintf(stderr, "BWTIncConstruct() : Not enough memory allocated!\n");
		exit(1);
	}

	bwtInc->bwt->occValue = mergedBwt - mergedOccSizeInWord;

	BWTClearTrailingBwtCode(bwtInc->bwt);
	BWTGenerateOccValueFromBwt(bwtInc->bwt->bwtCode, bwtInc->bwt->occValue, bwtInc->bwt->occValueMajor,
							   bwtInc->bwt->textLength, bwtInc->bwt->decodeTable);

	bwtInc->bwt->inverseSa0 = newInverseSa0;
	
	bwtInc->bwt->cumulativeFreq[1] += bwtInc->cumulativeCountInCurrentBuild[1] - (bwtInc->firstCharInLastIteration <= 0);
	bwtInc->bwt->cumulativeFreq[2] += bwtInc->cumulativeCountInCurrentBuild[2] - (bwtInc->firstCharInLastIteration <= 1);
	bwtInc->bwt->cumulativeFreq[3] += bwtInc->cumulativeCountInCurrentBuild[3] - (bwtInc->firstCharInLastIteration <= 2);
	bwtInc->bwt->cumulativeFreq[4] += bwtInc->cumulativeCountInCurrentBuild[4] - (bwtInc->firstCharInLastIteration <= 3);

	bwtInc->firstCharInLastIteration = firstCharInThisIteration;

	// Set build size and text address for the next build
	BWTIncSetBuildSizeAndTextAddr(bwtInc);
	bwtInc->numberOfIterationDone++;

}

static void BWTIncSetBuildSizeAndTextAddr(BWTInc *bwtInc) {

	unsigned int maxBuildSize;

	if (bwtInc->bwt->textLength == 0) {
		// initial build
		// Minus 2 because n+1 entries of seq and rank needed for n char
		maxBuildSize = (bwtInc->availableWord - 2 * CHAR_PER_WORD - OCC_INTERVAL / CHAR_PER_WORD - WORD_BETWEEN_OCC / 2)
							/ (2 * CHAR_PER_WORD + 1) * CHAR_PER_WORD;
		if (bwtInc->initialMaxBuildSize > 0) {
			bwtInc->buildSize = minX(bwtInc->initialMaxBuildSize, maxBuildSize);
		} else {
			bwtInc->buildSize = maxBuildSize;
		}
	} else {
		// Minus 3 because n+1 entries of sorted rank, seq and rank needed for n char
		// Minus numberOfIterationDone because bwt slightly shift to left in each iteration
		maxBuildSize = (bwtInc->availableWord - bwtInc->bwt->bwtSizeInWord - bwtInc->bwt->occSizeInWord - 3 * CHAR_PER_WORD
							 - bwtInc->numberOfIterationDone * OCC_INTERVAL / BIT_PER_CHAR - WORD_BETWEEN_OCC / 2) 
							 / 3;
		if (maxBuildSize < CHAR_PER_WORD) {
			fprintf(stderr, "BWTIncSetBuildSizeAndTextAddr(): Not enough space allocated to continue construction!\n");
			exit(1);
		}
		if (bwtInc->incMaxBuildSize > 0) {
            bwtInc->buildSize = minX(bwtInc->incMaxBuildSize, maxBuildSize);
		} else {
			bwtInc->buildSize = maxBuildSize;
		}
		if (bwtInc->buildSize < CHAR_PER_WORD) {
			bwtInc->buildSize = CHAR_PER_WORD;
		}
	}

	if (bwtInc->buildSize < CHAR_PER_WORD) {
		fprintf(stderr, "BWTIncSetBuildSizeAndTextAddr(): Not enough space allocated to continue construction!\n");
		exit(1);
	}

	bwtInc->buildSize = bwtInc->buildSize / CHAR_PER_WORD * CHAR_PER_WORD;

	bwtInc->packedText = bwtInc->workingMemory + 2 * (bwtInc->buildSize + CHAR_PER_WORD);
	bwtInc->textBuffer = (unsigned char*)(bwtInc->workingMemory + bwtInc->buildSize + CHAR_PER_WORD);

}

static void BWTIncPutPackedTextToRank(const unsigned int *packedText, unsigned int* __restrict rank, unsigned int* __restrict cumulativeCount, const unsigned int numChar) {

	unsigned int i, j;
	unsigned int packedMask;
	unsigned int rankIndex;
	unsigned int lastWord, numCharInLastWord;

	unsigned char ALIGN_16 temp[CHAR_PER_WORD];

	__m128i p1, p2, mask;


	lastWord = (numChar - 1) / CHAR_PER_WORD;
	numCharInLastWord = numChar - lastWord * CHAR_PER_WORD;

	packedMask = ALL_ONE_MASK >> (BITS_IN_WORD - BIT_PER_CHAR);
	rankIndex = numChar - 1;

	// Unpack word-packed text; temp[0] will be character in the least significant 2 bits
	p1 = _mm_cvtsi32_si128(packedText[lastWord]);
	p2 = _mm_srli_epi32(p1, 4);
	p1 = _mm_unpacklo_epi8(p1, p2);

	mask = _mm_set1_epi32(0x03030303);

	p2 = _mm_srli_epi32(p1, 2);
	p1 = _mm_unpacklo_epi8(p1, p2);

	p1 = _mm_and_si128(p1, mask);
	_mm_store_si128((__m128i*)temp, p1);

	for (i=CHAR_PER_WORD - numCharInLastWord; i<CHAR_PER_WORD; i++) {
		cumulativeCount[temp[i]+1]++;
		rank[rankIndex] = temp[i];
		rankIndex--;
	}

	for (i=lastWord; i--;) {	// loop from lastWord - 1 to 0

		// Unpack word-packed text; temp[0] will be character in the least significant 2 bits
		p1 = _mm_cvtsi32_si128(packedText[i]);
		p2 = _mm_srli_epi32(p1, 4);
		p1 = _mm_unpacklo_epi8(p1, p2);

		mask = _mm_set1_epi32(0x03030303);

		p2 = _mm_srli_epi32(p1, 2);
		p1 = _mm_unpacklo_epi8(p1, p2);

		p1 = _mm_and_si128(p1, mask);
		_mm_store_si128((__m128i*)temp, p1);

		for (j=0; j<CHAR_PER_WORD; j++) {
			cumulativeCount[temp[j]+1]++;
			rank[rankIndex] = temp[j];
			rankIndex--;
		}
	}

	// Convert occurrence to cumulativeCount
	cumulativeCount[2] += cumulativeCount[1];
	cumulativeCount[3] += cumulativeCount[2];
	cumulativeCount[4] += cumulativeCount[3];

}

static void BWTIncBuildPackedBwt(const unsigned int *relativeRank, unsigned int* __restrict bwt, const unsigned int numChar,
								 const unsigned int *cumulativeCount, const unsigned int *packedShift) {

	unsigned int i, c, r;
	unsigned int previousRank, currentRank;
	unsigned int wordIndex, charIndex;
	unsigned int inverseSa0;

	inverseSa0 = previousRank = relativeRank[0];

	for (i=1; i<=numChar; i++) {
		currentRank = relativeRank[i];
		// previousRank > cumulativeCount[c] because $ is one of the char
		c = (previousRank > cumulativeCount[1]) + (previousRank > cumulativeCount[2]) 
											    + (previousRank > cumulativeCount[3]);
		// set bwt for currentRank
		if (c > 0) {
			// c <> 'a'
			r = currentRank;
			if (r > inverseSa0) {
				// - 1 because $ at inverseSa0 is not encoded			
				r--;
			}
			wordIndex = r / CHAR_PER_WORD;
			charIndex = r - wordIndex * CHAR_PER_WORD;
			bwt[wordIndex] |= c << packedShift[charIndex];
		}
		previousRank = currentRank;
	}

}

static unsigned int BWTIncGetAbsoluteRank(BWT *bwt, unsigned int* __restrict absoluteRank, unsigned int* __restrict seq, const unsigned int *packedText, 
								  const unsigned int numChar, const unsigned int* cumulativeCount, const unsigned int firstCharInLastIteration) {

	unsigned int saIndex;
	unsigned int lastWord;
	unsigned int packedMask;
	unsigned int i, j;
	unsigned int c, t;
	unsigned int rankIndex;
	unsigned int shift;
	unsigned int seqIndexFromStart[ALPHABET_SIZE];
	unsigned int seqIndexFromEnd[ALPHABET_SIZE];

	for (i=0; i<ALPHABET_SIZE; i++) {
		seqIndexFromStart[i] = cumulativeCount[i];
		seqIndexFromEnd[i] = cumulativeCount[i+1] - 1;
	}

	shift = BITS_IN_WORD - BIT_PER_CHAR;
	packedMask = ALL_ONE_MASK >> shift;
	saIndex = bwt->inverseSa0;
	rankIndex = numChar - 1;

	lastWord = numChar / CHAR_PER_WORD;
	for (i=lastWord; i--;) {	// loop from lastWord - 1 to 0
		t = packedText[i];
		for (j=0; j<CHAR_PER_WORD; j++) {
			c = t & packedMask;
#ifdef DEBUG
			if (c >= ALPHABET_SIZE) {
				fprintf(stderr, "BWTIncGetAbsoluteRank() : c >= ALPHABET_SIZE!\n");
				exit(1);
			}
#endif
			saIndex = bwt->cumulativeFreq[c] + BWTOccValue(bwt, saIndex, c) + 1;
#ifdef DEBUG
			if (saIndex > bwt->textLength + 1) {
				fprintf(stderr, "BWTIncGetAbsoluteRank() : saIndex > bwt->textLength + 1!\n");
				exit(1);
			}
#endif
			// A counting sort using the first character of suffix is done here
			// If rank > inverseSa0 -> fill seq from end, otherwise fill seq from start -> to leave the right entry for inverseSa0
			if (saIndex > bwt->inverseSa0) {
				seq[seqIndexFromEnd[c]] = rankIndex;
				absoluteRank[seqIndexFromEnd[c]] = saIndex;
				seqIndexFromEnd[c]--;
			} else {
				seq[seqIndexFromStart[c]] = rankIndex;
				absoluteRank[seqIndexFromStart[c]] = saIndex;
				seqIndexFromStart[c]++;
			}
			rankIndex--;
			t >>= BIT_PER_CHAR;
		}
	}

	absoluteRank[seqIndexFromStart[firstCharInLastIteration]] = bwt->inverseSa0;	// representing the substring of all preceding characters
	seq[seqIndexFromStart[firstCharInLastIteration]] = numChar;

#ifdef DEBUG
	if (seqIndexFromStart[firstCharInLastIteration] != seqIndexFromEnd[firstCharInLastIteration]) {
		fprintf(stderr, "BWTIncGetAbsoluteRank(): seqIndexFromStart[firstCharInLastIteration] != seqIndexFromEnd[firstCharInLastIteration]!\n");
	}
#endif

	return seqIndexFromStart[firstCharInLastIteration];

}

static void BWTIncSortKey(unsigned int* __restrict key, unsigned int* __restrict seq, const unsigned int numItem) {

	#define EQUAL_KEY_THRESHOLD	4	// Partition for equal key if data array size / the number of data with equal value with pivot < EQUAL_KEY_THRESHOLD

	int lowIndex, highIndex, midIndex;
	int lowPartitionIndex, highPartitionIndex;
	int lowStack[32], highStack[32];
	int stackDepth;
	int i, j;
	unsigned int tempSeq, tempKey;
	int numberOfEqualKey;

	if (numItem < 2) {
		return;
	}

	stackDepth = 0;

    lowIndex = 0;
    highIndex = numItem - 1;

	for (;;) {

		for (;;) {

			// Sort small array of data
			if (highIndex - lowIndex < BWTINC_INSERT_SORT_NUM_ITEM) {	 // Insertion sort on smallest arrays
				for (i=lowIndex+1; i<=highIndex; i++) {
					tempSeq = seq[i];
					tempKey = key[i];
					for (j = i; j > lowIndex && key[j-1] > tempKey; j--) {
						seq[j] = seq[j-1];
						key[j] = key[j-1];
					}
					if (j != i) {
						seq[j] = tempSeq;
						key[j] = tempKey;
					}
				}
				break;
			}

			// Choose pivot as median of the lowest, middle, and highest data; sort the three data

			midIndex = average(lowIndex, highIndex);
			if (key[lowIndex] > key[midIndex]) {
				tempSeq = seq[lowIndex];
				tempKey = key[lowIndex];
				seq[lowIndex] = seq[midIndex];
				key[lowIndex] = key[midIndex];
				seq[midIndex] = tempSeq;
				key[midIndex] = tempKey;
			}
			if (key[lowIndex] > key[highIndex]) {
				tempSeq = seq[lowIndex];
				tempKey = key[lowIndex];
				seq[lowIndex] = seq[highIndex];
				key[lowIndex] = key[highIndex];
				seq[highIndex] = tempSeq;
				key[highIndex] = tempKey;
			}
			if (key[midIndex] > key[highIndex]) {
				tempSeq = seq[midIndex];
				tempKey = key[midIndex];
				seq[midIndex] = seq[highIndex];
				key[midIndex] = key[highIndex];
				seq[highIndex] = tempSeq;
				key[highIndex] = tempKey;
			}

			// Partition data

			numberOfEqualKey = 0;

			lowPartitionIndex = lowIndex + 1;
			highPartitionIndex = highIndex - 1;

			for (;;) {
				while (lowPartitionIndex <= highPartitionIndex && key[lowPartitionIndex] <= key[midIndex]) {
					numberOfEqualKey += (key[lowPartitionIndex] == key[midIndex]);
					lowPartitionIndex++;
				}
				while (lowPartitionIndex < highPartitionIndex) {
					if (key[midIndex] >= key[highPartitionIndex]) {
						numberOfEqualKey += (key[midIndex] == key[highPartitionIndex]);
						break;
					}
					highPartitionIndex--;
				}
				if (lowPartitionIndex >= highPartitionIndex) {
					break;
				}
				tempSeq = seq[lowPartitionIndex];
				tempKey = key[lowPartitionIndex];
				seq[lowPartitionIndex] = seq[highPartitionIndex];
				key[lowPartitionIndex] = key[highPartitionIndex];
				seq[highPartitionIndex] = tempSeq;
				key[highPartitionIndex] = tempKey;
				if (highPartitionIndex == midIndex) {
					// partition key has been moved
					midIndex = lowPartitionIndex;
				}
				lowPartitionIndex++;
				highPartitionIndex--;
			}

			// Adjust the partition index
			highPartitionIndex = lowPartitionIndex;
			lowPartitionIndex--;

			// move the partition key to end of low partition
			tempSeq = seq[midIndex];
			tempKey = key[midIndex];
			seq[midIndex] = seq[lowPartitionIndex];
			key[midIndex] = key[lowPartitionIndex];
			seq[lowPartitionIndex] = tempSeq;
			key[lowPartitionIndex] = tempKey;

			if (highIndex - lowIndex + BWTINC_INSERT_SORT_NUM_ITEM <= EQUAL_KEY_THRESHOLD * numberOfEqualKey) {

				// Many keys = partition key; separate the equal key data from the lower partition
		
				midIndex = lowIndex;

				for (;;) {
					while (midIndex < lowPartitionIndex && key[midIndex] < key[lowPartitionIndex]) {
						midIndex++;
					}
					while (midIndex < lowPartitionIndex && key[lowPartitionIndex] == key[lowPartitionIndex - 1]) {
						lowPartitionIndex--;
					}
					if (midIndex >= lowPartitionIndex) {
						break;
					}
					tempSeq = seq[midIndex];
					tempKey = key[midIndex];
					seq[midIndex] = seq[lowPartitionIndex - 1];
					key[midIndex] = key[lowPartitionIndex - 1];
					seq[lowPartitionIndex - 1] = tempSeq;
					key[lowPartitionIndex - 1] = tempKey;
					midIndex++;
					lowPartitionIndex--;
				}

			}

			if (lowPartitionIndex - lowIndex > highIndex - highPartitionIndex) {
				// put the larger partition to stack
				lowStack[stackDepth] = lowIndex;
				highStack[stackDepth] = lowPartitionIndex - 1;
				stackDepth++;
				// sort the smaller partition first
				lowIndex = highPartitionIndex;
			} else {
				// put the larger partition to stack
				lowStack[stackDepth] = highPartitionIndex;
				highStack[stackDepth] = highIndex;
				stackDepth++;
				// sort the smaller partition first
				if (lowPartitionIndex > lowIndex) {
					highIndex = lowPartitionIndex - 1;
				} else {
					// all keys in the partition equals to the partition key
					break;
				}
			}
			continue;

		}

		// Pop a range from stack
		if (stackDepth > 0) {
			stackDepth--;
			lowIndex = lowStack[stackDepth];
			highIndex = highStack[stackDepth];
			continue;
		} else {
			return;
		}

	}


}

static void BWTIncBuildRelativeRank(unsigned int* __restrict sortedRank, unsigned int* __restrict seq, unsigned int* __restrict relativeRank, const unsigned int numItem, 
									unsigned int oldInverseSa0, const unsigned int *cumulativeCount) {

	unsigned int i, c;
	unsigned int s, r;
	unsigned int lastRank, lastIndex;
	unsigned int oldInverseSa0RelativeRank = 0;
	unsigned int freq;

	lastIndex = numItem;
	lastRank = sortedRank[numItem];
	if (lastRank > oldInverseSa0) {
		sortedRank[numItem]--;	// to prepare for merging; $ is not encoded in bwt
	}
	s = seq[numItem];
	relativeRank[s] = numItem;
	if (lastRank == oldInverseSa0) {
		oldInverseSa0RelativeRank = numItem;
		oldInverseSa0++;	// so that this segment of code is not run again
		lastRank++;			// so that oldInverseSa0 become a sorted group with 1 item
	}

	c = ALPHABET_SIZE - 1;
	freq = cumulativeCount[c];

	for (i=numItem; i--;) {	// from numItem - 1 to 0
		r = sortedRank[i];
		if (r > oldInverseSa0) {
			sortedRank[i]--;	// to prepare for merging; $ is not encoded in bwt
		}
		s = seq[i];
		if (i < freq) {
			if (lastIndex >= freq) {
				lastRank++;	// to trigger the group across alphabet boundary to be split
			}
			c--;
			freq = cumulativeCount[c];
		}
		if (r == lastRank) {
			relativeRank[s] = lastIndex;
		} else {
			if (i == lastIndex - 1) {
				if (lastIndex < numItem && (int)seq[lastIndex + 1] < 0) {
					seq[lastIndex] = seq[lastIndex + 1] - 1;
				} else {
					seq[lastIndex] = (unsigned int)-1;
				}
			}
			lastIndex = i;
			lastRank = r;
			relativeRank[s] = i;
			if (r == oldInverseSa0) {
				oldInverseSa0RelativeRank = i;
				oldInverseSa0++;	// so that this segment of code is not run again
				lastRank++;			// so that oldInverseSa0 become a sorted group with 1 item
			}
		}
	}

}

static void BWTIncBuildBwt(unsigned int*  seq, const unsigned int *relativeRank, const unsigned int numChar, const unsigned int *cumulativeCount) {

	unsigned int i, c;
	unsigned int previousRank, currentRank;

	previousRank = relativeRank[0];

	for (i=1; i<=numChar; i++) {
		currentRank = relativeRank[i];
		c = (previousRank >= cumulativeCount[1]) + (previousRank >= cumulativeCount[2])
											  	 + (previousRank >= cumulativeCount[3]);
		seq[currentRank] = c;
		previousRank = currentRank;
	}

}

static void BWTIncMergeBwt(const unsigned int *sortedRank, const unsigned int* oldBwt, const unsigned int *insertBwt, unsigned int* __restrict mergedBwt, 
						   const unsigned int numOldBwt, const unsigned int numInsertBwt) {

	unsigned int bitsInWordMinusBitPerChar;
	unsigned int leftShift, rightShift;
	unsigned int o;
	unsigned int oIndex, iIndex, mIndex;
	unsigned int mWord, mChar, oWord, oChar;
	unsigned int numInsert;

	bitsInWordMinusBitPerChar = BITS_IN_WORD - BIT_PER_CHAR;

	oIndex = 0;
	iIndex = 0;
	mIndex = 0;

	mWord = 0;
	mChar = 0;

	mergedBwt[0] = 0;	// this can be cleared as merged Bwt slightly shift to the left in each iteration

	while (oIndex < numOldBwt) {

		// copy from insertBwt
		while (iIndex <= numInsertBwt && sortedRank[iIndex] <= oIndex) {
			if (sortedRank[iIndex] != 0) {	// special value to indicate that this is for new inverseSa0
				mergedBwt[mWord] |= insertBwt[iIndex] << (BITS_IN_WORD - (mChar + 1) * BIT_PER_CHAR);
				mIndex++;
				mChar++;
				if (mChar == CHAR_PER_WORD) {
					mChar = 0;
					mWord++;
					mergedBwt[mWord] = 0;	// no need to worry about crossing mergedBwt boundary
				}
			}
			iIndex++;
		}

		// Copy from oldBwt to mergedBwt
		if (iIndex <= numInsertBwt) {
			o = sortedRank[iIndex];
		} else {
			o = numOldBwt;
		}
		numInsert = o - oIndex;

		oWord = oIndex / CHAR_PER_WORD;
		oChar = oIndex - oWord * CHAR_PER_WORD;
		if (oChar > mChar) {
			leftShift = (oChar - mChar) * BIT_PER_CHAR;
			rightShift = (CHAR_PER_WORD + mChar - oChar) * BIT_PER_CHAR;
			mergedBwt[mWord] = mergedBwt[mWord]
								| (oldBwt[oWord] << (oChar * BIT_PER_CHAR) >> (mChar * BIT_PER_CHAR))
								| (oldBwt[oWord+1] >> rightShift);
			oIndex += minX(numInsert, CHAR_PER_WORD - mChar);
			while (o > oIndex) {
				oWord++;
				mWord++;
				mergedBwt[mWord] = (oldBwt[oWord] << leftShift) | (oldBwt[oWord+1] >> rightShift);
				oIndex += CHAR_PER_WORD;
			}
		} else if (oChar < mChar) {
			rightShift = (mChar - oChar) * BIT_PER_CHAR;
			leftShift = (CHAR_PER_WORD + oChar - mChar) * BIT_PER_CHAR;
			mergedBwt[mWord] = mergedBwt[mWord] 
								| (oldBwt[oWord] << (oChar * BIT_PER_CHAR) >> (mChar * BIT_PER_CHAR));
			oIndex += minX(numInsert, CHAR_PER_WORD - mChar);
			while (o > oIndex) {
				oWord++;
				mWord++;
				mergedBwt[mWord] = (oldBwt[oWord-1] << leftShift) | (oldBwt[oWord] >> rightShift);
				oIndex += CHAR_PER_WORD;
			}
		} else { // oChar == mChar
			mergedBwt[mWord] = mergedBwt[mWord] | truncateLeft(oldBwt[oWord], mChar * BIT_PER_CHAR);
			oIndex += minX(numInsert, CHAR_PER_WORD - mChar);
			while (o > oIndex) {
				oWord++;
				mWord++;
				mergedBwt[mWord] = oldBwt[oWord];
				oIndex += CHAR_PER_WORD;
			}
		}
		oIndex = o;
		mIndex += numInsert;

		// Clear the trailing garbage in mergedBwt
		mWord = mIndex / CHAR_PER_WORD;
		mChar = mIndex - mWord * CHAR_PER_WORD;
		if (mChar == 0) {
			mergedBwt[mWord] = 0;
		} else {
			mergedBwt[mWord] = truncateRight(mergedBwt[mWord], (BITS_IN_WORD - mChar * BIT_PER_CHAR));
		}

	}

	// copy from insertBwt
	while (iIndex <= numInsertBwt) {
		if (sortedRank[iIndex] != 0) {
			mergedBwt[mWord] |= insertBwt[iIndex] << (BITS_IN_WORD - (mChar + 1) * BIT_PER_CHAR);
			mIndex++;
			mChar++;
			if (mChar == CHAR_PER_WORD) {
				mChar = 0;
				mWord++;
				mergedBwt[mWord] = 0;	// no need to worry about crossing mergedBwt boundary
			}
		}
		iIndex++;
	}

}
unsigned int BWTGenerateOccValueToFileFromBwt(const char *bwtFileName, const char *occValueFileName, unsigned int*  decodeTable) {

	FILE *bwtFile, *occValueFile;
	unsigned int *bwt;
	unsigned int bwtFileSizeInWord, bwtResidentSizeInWord;
	unsigned int textLength;
	unsigned int inverseSa0;
	unsigned int *occValue, *occValueMajor;
	unsigned int occSizeInWord, occMajorSizeInWord;
	unsigned int i;
	unsigned int cumulativeFreq[ALPHABET_SIZE];

	bwtFile = (FILE*)fopen64(bwtFileName, "rb");
	if (bwtFile == NULL) {
		fprintf(stderr, "BWTGenerateOccValueToFileFromBwt(): Cannot open BWT file!\n");
		exit(1);
	}

	fread(&inverseSa0, sizeof(unsigned int), 1, bwtFile);

	fread(cumulativeFreq, sizeof(unsigned int), ALPHABET_SIZE, bwtFile);
	textLength = cumulativeFreq[ALPHABET_SIZE - 1];

	bwtResidentSizeInWord = BWTResidentSizeInWord(textLength);
	bwt = MMUnitAllocate(bwtResidentSizeInWord * sizeof(unsigned int));
	bwtFileSizeInWord = BWTFileSizeInWord(textLength);
	fread(bwt, sizeof(unsigned int), bwtFileSizeInWord, bwtFile);
	fclose(bwtFile);
	for (i=bwtFileSizeInWord; i<bwtResidentSizeInWord; i++) {
		bwt[i] = 0;
	}

	// occValue File

	occValueFile = (FILE*)fopen64(occValueFileName, "wb");
	if (occValueFile == NULL) {
		fprintf(stderr, "BWTGenerateOccValueToFileFromBwt(): Cannot open occ value file!\n");
		exit(1);
	}

	fwrite(&inverseSa0, sizeof(unsigned int), 1, occValueFile);
	fwrite(cumulativeFreq, sizeof(unsigned int), ALPHABET_SIZE, occValueFile);

	occSizeInWord = BWTOccValueMinorSizeInWord(textLength);
	occMajorSizeInWord = BWTOccValueMajorSizeInWord(textLength);
	occValue = MMUnitAllocate(occSizeInWord * sizeof(unsigned int));
	occValueMajor = MMUnitAllocate(occMajorSizeInWord * sizeof(unsigned int));

	if (decodeTable == NULL) {
		decodeTable = MMUnitAllocate(DNA_OCC_CNT_TABLE_SIZE_IN_WORD * sizeof(unsigned int));
		GenerateDNAOccCountTable(decodeTable);
	}

	BWTGenerateOccValueFromBwt(bwt, occValue, occValueMajor, textLength, decodeTable);

	fwrite(occValue, sizeof(unsigned int), occSizeInWord, occValueFile);
	fwrite(occValueMajor, sizeof(unsigned int), occMajorSizeInWord, occValueFile);
	fclose(occValueFile);

	MMUnitFree(occValue, occSizeInWord * sizeof(unsigned int));
	MMUnitFree(occValueMajor, occMajorSizeInWord * sizeof(unsigned int));
	MMUnitFree(bwt, bwtResidentSizeInWord * sizeof(unsigned int));
	MMUnitFree(decodeTable, DNA_OCC_CNT_TABLE_SIZE_IN_WORD * sizeof(unsigned int));

	return textLength;

}

void BWTGenerateOccValueFromBwt(const unsigned int*  bwt, unsigned int* __restrict occValue, unsigned int* __restrict occValueMajor,
								const unsigned int textLength, const unsigned int*  decodeTable) {

	unsigned int numberOfOccValueMajor, numberOfOccValue;
	unsigned int wordBetweenOccValue;
	unsigned int numberOfOccIntervalPerMajor;
	unsigned int c;
	unsigned int i, j;
	unsigned int occMajorIndex;
	unsigned int occIndex, bwtIndex;
	unsigned int sum;
	unsigned int tempOccValue0[ALPHABET_SIZE], tempOccValue1[ALPHABET_SIZE];

	wordBetweenOccValue = OCC_INTERVAL / CHAR_PER_WORD;

	// Calculate occValue

	numberOfOccValue = (textLength + OCC_INTERVAL - 1) / OCC_INTERVAL + 1;				// Value at both end for bi-directional encoding
	numberOfOccIntervalPerMajor = OCC_INTERVAL_MAJOR / OCC_INTERVAL;
	numberOfOccValueMajor = (numberOfOccValue + numberOfOccIntervalPerMajor - 1) / numberOfOccIntervalPerMajor;

	tempOccValue0[0] = 0;
	tempOccValue0[1] = 0;
	tempOccValue0[2] = 0;
	tempOccValue0[3] = 0;
	occValueMajor[0] = 0;
	occValueMajor[1] = 0;
	occValueMajor[2] = 0;
	occValueMajor[3] = 0;

	occIndex = 0;
	bwtIndex = 0;
	for (occMajorIndex=1; occMajorIndex<numberOfOccValueMajor; occMajorIndex++) {

		for (i=0; i<numberOfOccIntervalPerMajor/2; i++) {

			sum = 0;
			tempOccValue1[0] = tempOccValue0[0];
			tempOccValue1[1] = tempOccValue0[1];
			tempOccValue1[2] = tempOccValue0[2];
			tempOccValue1[3] = tempOccValue0[3];

			for (j=0; j<wordBetweenOccValue; j++) {
				c = bwt[bwtIndex];
				sum += decodeTable[c >> 16];
				sum += decodeTable[c & 0x0000FFFF];
				bwtIndex++;
			}
			if (!DNA_OCC_SUM_EXCEPTION(sum)) {
				tempOccValue1[0] += (sum & 0x000000FF);	sum >>= 8;
				tempOccValue1[1] += (sum & 0x000000FF);	sum >>= 8;
				tempOccValue1[2] += (sum & 0x000000FF);	sum >>= 8;
				tempOccValue1[3] += sum;
			} else {
				if (sum == 0x00000100) {
					tempOccValue1[0] += 256;
				} else if (sum == 0x00010000) {
					tempOccValue1[1] += 256;
				} else if (sum == 0x01000000) {
					tempOccValue1[2] += 256;
				} else {
					tempOccValue1[3] += 256;
				}
			}
			occValue[occIndex * 4 + 0] = (tempOccValue0[0] << 16) | tempOccValue1[0];
			occValue[occIndex * 4 + 1] = (tempOccValue0[1] << 16) | tempOccValue1[1];
			occValue[occIndex * 4 + 2] = (tempOccValue0[2] << 16) | tempOccValue1[2];
			occValue[occIndex * 4 + 3] = (tempOccValue0[3] << 16) | tempOccValue1[3];
			tempOccValue0[0] = tempOccValue1[0];
			tempOccValue0[1] = tempOccValue1[1];
			tempOccValue0[2] = tempOccValue1[2];
			tempOccValue0[3] = tempOccValue1[3];
			sum = 0;

			occIndex++;

			for (j=0; j<wordBetweenOccValue; j++) {
				c = bwt[bwtIndex];
				sum += decodeTable[c >> 16];
				sum += decodeTable[c & 0x0000FFFF];
				bwtIndex++;
			}
			if (!DNA_OCC_SUM_EXCEPTION(sum)) {
				tempOccValue0[0] += (sum & 0x000000FF);	sum >>= 8;
				tempOccValue0[1] += (sum & 0x000000FF);	sum >>= 8;
				tempOccValue0[2] += (sum & 0x000000FF);	sum >>= 8;
				tempOccValue0[3] += sum;
			} else {
				if (sum == 0x00000100) {
					tempOccValue0[0] += 256;
				} else if (sum == 0x00010000) {
					tempOccValue0[1] += 256;
				} else if (sum == 0x01000000) {
					tempOccValue0[2] += 256;
				} else {
					tempOccValue0[3] += 256;
				}
			}
		}

		occValueMajor[occMajorIndex * 4 + 0] = occValueMajor[(occMajorIndex - 1) * 4 + 0] + tempOccValue0[0];
		occValueMajor[occMajorIndex * 4 + 1] = occValueMajor[(occMajorIndex - 1) * 4 + 1] + tempOccValue0[1];
		occValueMajor[occMajorIndex * 4 + 2] = occValueMajor[(occMajorIndex - 1) * 4 + 2] + tempOccValue0[2];
		occValueMajor[occMajorIndex * 4 + 3] = occValueMajor[(occMajorIndex - 1) * 4 + 3] + tempOccValue0[3];
		tempOccValue0[0] = 0;
		tempOccValue0[1] = 0;
		tempOccValue0[2] = 0;
		tempOccValue0[3] = 0;

	}

	while (occIndex < (numberOfOccValue-1)/2) {
		sum = 0;
		tempOccValue1[0] = tempOccValue0[0];
		tempOccValue1[1] = tempOccValue0[1];
		tempOccValue1[2] = tempOccValue0[2];
		tempOccValue1[3] = tempOccValue0[3];
		for (j=0; j<wordBetweenOccValue; j++) {
			c = bwt[bwtIndex];
			sum += decodeTable[c >> 16];
			sum += decodeTable[c & 0x0000FFFF];
			bwtIndex++;
		}
		if (!DNA_OCC_SUM_EXCEPTION(sum)) {
			tempOccValue1[0] += (sum & 0x000000FF);	sum >>= 8;
			tempOccValue1[1] += (sum & 0x000000FF);	sum >>= 8;
			tempOccValue1[2] += (sum & 0x000000FF);	sum >>= 8;
			tempOccValue1[3] += sum;
		} else {
			if (sum == 0x00000100) {
				tempOccValue1[0] += 256;
			} else if (sum == 0x00010000) {
				tempOccValue1[1] += 256;
			} else if (sum == 0x01000000) {
				tempOccValue1[2] += 256;
			} else {
				tempOccValue1[3] += 256;
			}
		}
		occValue[occIndex * 4 + 0] = (tempOccValue0[0] << 16) | tempOccValue1[0];
		occValue[occIndex * 4 + 1] = (tempOccValue0[1] << 16) | tempOccValue1[1];
		occValue[occIndex * 4 + 2] = (tempOccValue0[2] << 16) | tempOccValue1[2];
		occValue[occIndex * 4 + 3] = (tempOccValue0[3] << 16) | tempOccValue1[3];
		tempOccValue0[0] = tempOccValue1[0];
		tempOccValue0[1] = tempOccValue1[1];
		tempOccValue0[2] = tempOccValue1[2];
		tempOccValue0[3] = tempOccValue1[3];
		sum = 0;
		occIndex++;

		for (j=0; j<wordBetweenOccValue; j++) {
			c = bwt[bwtIndex];
			sum += decodeTable[c >> 16];
			sum += decodeTable[c & 0x0000FFFF];
			bwtIndex++;
		}
		if (!DNA_OCC_SUM_EXCEPTION(sum)) {
			tempOccValue0[0] += (sum & 0x000000FF);	sum >>= 8;
			tempOccValue0[1] += (sum & 0x000000FF);	sum >>= 8;
			tempOccValue0[2] += (sum & 0x000000FF);	sum >>= 8;
			tempOccValue0[3] += sum;
		} else {
			if (sum == 0x00000100) {
				tempOccValue0[0] += 256;
			} else if (sum == 0x00010000) {
				tempOccValue0[1] += 256;
			} else if (sum == 0x01000000) {
				tempOccValue0[2] += 256;
			} else {
				tempOccValue0[3] += 256;
			}
		}
	}

	sum = 0;
	tempOccValue1[0] = tempOccValue0[0];
	tempOccValue1[1] = tempOccValue0[1];
	tempOccValue1[2] = tempOccValue0[2];
	tempOccValue1[3] = tempOccValue0[3];

	if (occIndex * 2 < numberOfOccValue - 1) {
		for (j=0; j<wordBetweenOccValue; j++) {
			c = bwt[bwtIndex];
			sum += decodeTable[c >> 16];
			sum += decodeTable[c & 0x0000FFFF];
			bwtIndex++;
		}
		if (!DNA_OCC_SUM_EXCEPTION(sum)) {
			tempOccValue1[0] += (sum & 0x000000FF);	sum >>= 8;
			tempOccValue1[1] += (sum & 0x000000FF);	sum >>= 8;
			tempOccValue1[2] += (sum & 0x000000FF);	sum >>= 8;
			tempOccValue1[3] += sum;
		} else {
			if (sum == 0x00000100) {
				tempOccValue1[0] += 256;
			} else if (sum == 0x00010000) {
				tempOccValue1[1] += 256;
			} else if (sum == 0x01000000) {
				tempOccValue1[2] += 256;
			} else {
				tempOccValue1[3] += 256;
			}
		}
	}

	occValue[occIndex * 4 + 0] = (tempOccValue0[0] << 16) | tempOccValue1[0];
	occValue[occIndex * 4 + 1] = (tempOccValue0[1] << 16) | tempOccValue1[1];
	occValue[occIndex * 4 + 2] = (tempOccValue0[2] << 16) | tempOccValue1[2];
	occValue[occIndex * 4 + 3] = (tempOccValue0[3] << 16) | tempOccValue1[3];

}

void BWTSaveBwtCodeAndOcc(const BWT *bwt, const char *bwtFileName, const char *occValueFileName) {

	FILE *bwtFile, *occValueFile;
	unsigned int bwtLength;

	bwtFile = (FILE*)fopen64(bwtFileName, "wb");
	if (bwtFile == NULL) {
		fprintf(stderr, "BWTSaveBwtCodeAndOcc(): Cannot open BWT code file!\n");
		exit(1);
	}

	fwrite(&bwt->inverseSa0, sizeof(unsigned int), 1, bwtFile);
	fwrite(bwt->cumulativeFreq + 1, sizeof(unsigned int), ALPHABET_SIZE, bwtFile);
	bwtLength = BWTFileSizeInWord(bwt->textLength);
	fwrite(bwt->bwtCode, sizeof(unsigned int), bwtLength, bwtFile);
	fclose(bwtFile);

	occValueFile = (FILE*)fopen64(occValueFileName, "wb");
	if (occValueFile == NULL) {
		fprintf(stderr, "BWTSaveBwtCodeAndOcc(): Cannot open occ value file!\n");
		exit(1);
	}

	fwrite(&bwt->inverseSa0, sizeof(unsigned int), 1, occValueFile);
	fwrite(bwt->cumulativeFreq + 1, sizeof(unsigned int), ALPHABET_SIZE, occValueFile);

	fwrite(bwt->occValue, sizeof(unsigned int), bwt->occSizeInWord, occValueFile);
	fwrite(bwt->occValueMajor, sizeof(unsigned int), bwt->occMajorSizeInWord, occValueFile);
	fclose(occValueFile);

}

void BWTGenerateSaValue(BWT *bwt, const unsigned int saValueFreq, unsigned int showProgressInterval) {

	unsigned int saValue;
	unsigned int saIndex;
	unsigned int numberOfSaValue;
	unsigned int numberOfSaValueGenerated;
	unsigned int progressInterval;
	unsigned int i;

	if (bwt->saInterval != ALL_ONE_MASK) {
		fprintf(stderr, "BWTGenerateSaValue() : saValue already exist!\n");
		exit(1);
	}

	if (saValueFreq == 1) {
		BWTGenerateFullSaValue(bwt);
		return;
	}
	saValue = bwt->textLength;
	saIndex = 0;

	numberOfSaValue = (bwt->textLength + saValueFreq) / saValueFreq;
	bwt->saValueSizeInWord = numberOfSaValue;
	bwt->saInterval = saValueFreq;

	bwt->saValue = MMUnitAllocate(bwt->saValueSizeInWord * sizeof(unsigned int));

	#ifdef DEBUG
	initializeVAL(bwt->saValue, numberOfSaValue, ALL_ONE_MASK);
	#endif

	numberOfSaValueGenerated = 0;

	if (showProgressInterval == 0 || showProgressInterval > numberOfSaValue) {
		progressInterval = numberOfSaValue;
	} else {
		progressInterval = showProgressInterval;
	}

	while (numberOfSaValueGenerated < numberOfSaValue) {
		progressInterval = minX(numberOfSaValue - numberOfSaValueGenerated, progressInterval);
		for (i=0; i<progressInterval; i++) {
			while (saIndex % saValueFreq != 0) {
				#ifdef DEBUG
				if (saValue == 0) {
					fprintf(stderr, "BWTGenerateSaValue(): saValue < 0!\n");
					exit(1);
				}
				#endif
				saValue--;
				saIndex = BWTPsiMinusValue(bwt, saIndex);
				#ifdef DEBUG
				if (saIndex > bwt->textLength) {
					fprintf(stderr, "BWTGenerateSaValue() : saIndex > textLength!\n");
					exit(1);
				}
				#endif
			}
			bwt->saValue[saIndex/saValueFreq] = saValue;
			saValue--;
			saIndex = BWTPsiMinusValue(bwt, saIndex);
		}
		numberOfSaValueGenerated += progressInterval;
		if (showProgressInterval > 0) {
			printf("SA Value generated : %u\n", numberOfSaValueGenerated);
		}
	}

	#ifdef DEBUG
	if (numberOfMatchInVAL(bwt->saValue, numberOfSaValue, ALL_ONE_MASK) > 0) {
		fprintf(stderr, "BWTGenerateSaValue() : some saValue is not filled!\n");
		exit(1);
	}
	#endif

	bwt->saValue[0] = (unsigned int)-1;	// Special handling

}

void BWTGenerateFullSaValue(BWT *bwt) {

	unsigned int i, j;
	unsigned int c, t, n;
	unsigned int freq[ALPHABET_SIZE];

	if (bwt->saInterval != ALL_ONE_MASK) {
		fprintf(stderr, "BWTGenerateSaValue() : saValue already exist!\n");
		exit(1);
	}

	bwt->saValueSizeInWord = bwt->textLength + 1;
	bwt->saInterval = 1;

	bwt->saValue = MMUnitAllocate(bwt->saValueSizeInWord * sizeof(unsigned int));

	for (i=0; i<ALPHABET_SIZE; i++) {
		freq[i] = 0;
	}

	n = 0;
	for (i=0; i<bwt->textLength / CHAR_PER_WORD; i++) {
		t = bwt->bwtCode[i];
		for (j=0; j<CHAR_PER_WORD; j++) {
			c = t >> (BITS_IN_WORD - BIT_PER_CHAR);
			t <<= BIT_PER_CHAR;
			freq[c]++;
			bwt->saValue[n + (n>=bwt->inverseSa0)] = freq[c] + bwt->cumulativeFreq[c];
			n++;
		}
	}
	t = bwt->bwtCode[i];
	for (j=0; j<bwt->textLength - i * CHAR_PER_WORD; j++) {
		c = t >> (BITS_IN_WORD - BIT_PER_CHAR);
		t <<= BIT_PER_CHAR;
		freq[c]++;
		bwt->saValue[n + (n>=bwt->inverseSa0)] = freq[c] + bwt->cumulativeFreq[c];
		n++;
	}

	t = 0;
	for (i=bwt->textLength; i>0; i--) {
		c = t;
		t = bwt->saValue[c];
		bwt->saValue[c] = i;
	}

	bwt->saValue[bwt->inverseSa0] = 0;

	bwt->saValue[0] = (unsigned int)-1;	// Special handling

}

void BWTSaveSaValue(const BWT *bwt, const char *saValueFileName) {

	FILE *saValueFile;

	saValueFile = (FILE*)fopen64(saValueFileName, "wb");
	if (saValueFile == NULL) {
		fprintf(stderr, "BWTSaveSaValue() : cannot open saValueFile!\n");
		exit(1);
	}

	fwrite(&bwt->inverseSa0, sizeof(unsigned int), 1, saValueFile);
	fwrite(bwt->cumulativeFreq + 1, sizeof(unsigned int), ALPHABET_SIZE, saValueFile);
	fwrite(&bwt->saInterval, sizeof(unsigned int), 1, saValueFile);

	fwrite(&bwt->textLength, sizeof(unsigned int), 1, saValueFile);	// Save SA values without special handling on SA[0]
	fwrite(bwt->saValue + 1, sizeof(unsigned int), bwt->saValueSizeInWord - 1, saValueFile);

	fclose(saValueFile);

}

void BWTGenerateInverseSa(BWT *bwt, const unsigned int inverseSaFreq, unsigned int showProgressInterval) {

	unsigned int saValue;
	unsigned int saIndex;
	unsigned int numberOfInverseSa;
	unsigned int numberOfInverseSaGenerated;
	unsigned int progressInterval;
	unsigned int i;
	
	if (bwt->inverseSaInterval != ALL_ONE_MASK) {
		fprintf(stderr, "BWTGenerateInverseSa() : inverseSa already exist!\n");
		exit(1);
	}

	saValue = bwt->textLength;
	saIndex = 0;

	numberOfInverseSa = (bwt->textLength + inverseSaFreq) / inverseSaFreq;
	bwt->inverseSaSizeInWord = numberOfInverseSa;
	bwt->inverseSaInterval = inverseSaFreq;

	bwt->inverseSa = MMUnitAllocate(bwt->inverseSaSizeInWord * sizeof(unsigned int));

	#ifdef DEBUG
	initializeVAL(bwt->inverseSa, numberOfInverseSa, ALL_ONE_MASK);
	#endif

	numberOfInverseSaGenerated = 0;

	if (showProgressInterval == 0 || showProgressInterval > numberOfInverseSa) {
		progressInterval = numberOfInverseSa;
	} else {
		progressInterval = showProgressInterval;
	}

	while (numberOfInverseSaGenerated < numberOfInverseSa) {
		progressInterval = minX(numberOfInverseSa - numberOfInverseSaGenerated, progressInterval);
		for (i=0; i<progressInterval; i++) {
			while (saValue % inverseSaFreq != 0) {
				saValue--;
				saIndex = BWTPsiMinusValue(bwt, saIndex);
				#ifdef DEBUG
				if (saIndex > bwt->textLength) {
					fprintf(stderr, "BWTGenerateInverseSa() : saIndex > textLength!\n");
					exit(1);
				}
				#endif
			}
			bwt->inverseSa[saValue/inverseSaFreq] = saIndex;
			saValue--;
			saIndex = BWTPsiMinusValue(bwt, saIndex);
		}
		numberOfInverseSaGenerated += progressInterval;
		if (showProgressInterval > 0) {
			printf("Inverse SA generated : %u\n", numberOfInverseSaGenerated);
		}
	}

	#ifdef DEBUG
	if (numberOfMatchInVAL(bwt->inverseSa, numberOfInverseSa, ALL_ONE_MASK) > 0) {
		fprintf(stderr, "BWTGenerateInverseSa() : some inverseSa is not filled!\n");
		exit(1);
	}
	#endif

}

void BWTSaveInverseSa(const BWT *bwt, const char *inverseSaFileName) {

	FILE *inverseSaFile;

	inverseSaFile = (FILE*)fopen64(inverseSaFileName, "wb");
	if (inverseSaFile == NULL) {
		fprintf(stderr, "BWTSaveInverseSa() : cannot open inverseSaFile!\n");
		exit(1);
	}

	fwrite(&bwt->inverseSa0, sizeof(unsigned int), 1, inverseSaFile);
	fwrite(bwt->cumulativeFreq + 1, sizeof(unsigned int), ALPHABET_SIZE, inverseSaFile);
	fwrite(&bwt->inverseSaInterval, sizeof(unsigned int), 1, inverseSaFile);
	fwrite(bwt->inverseSa, sizeof(unsigned int), bwt->inverseSaSizeInWord, inverseSaFile);

	fclose(inverseSaFile);

}

void BWTGenerateCachedSaIndex(const BWT *bwt, const unsigned int numOfChar, const char *cachedSaIndexFileName) {


	// stack
	unsigned int ALIGN_16 saIndexLeft[MAX_ARPROX_MATCH_LENGTH * ALPHABET_SIZE];
	unsigned int ALIGN_16 saIndexRight[MAX_ARPROX_MATCH_LENGTH * ALPHABET_SIZE];
	unsigned char generatedPattern[MAX_ARPROX_MATCH_LENGTH];
	unsigned int pos;
	unsigned int i;
	unsigned int numOfPattern;
	unsigned int saRangeIndex, saIndex;
	unsigned int a;
	unsigned int mask[16] = { 0x00000000, 0x00000003, 0x0000000F, 0x0000003F,
							  0x000000FF, 0x000003FF, 0x00000FFF, 0x00003FFF,
							  0x0000FFFF, 0x0003FFFF, 0x000FFFFF, 0x003FFFFF,
							  0x00FFFFFF, 0x03FFFFFF, 0x0FFFFFFF, 0x3FFFFFFF };

	unsigned int *cachedSaIndex;
	__m128i r1, r2, cf, cfp1; 

	FILE * cachedSaIndexFile;

	cachedSaIndexFile = (FILE*)fopen64(cachedSaIndexFileName, "wb");
	if (cachedSaIndexFile == NULL) {
		fprintf(stderr, "BWTGenerateSaRangeTable(): Cannot open SA range file!\n");
		exit(1);
	}

	fwrite(&bwt->inverseSa0, sizeof(unsigned int), 1, cachedSaIndexFile);
	fwrite(bwt->cumulativeFreq + 1, sizeof(unsigned int), ALPHABET_SIZE, cachedSaIndexFile);
	fwrite(&numOfChar, sizeof(unsigned int), 1, cachedSaIndexFile);

	numOfPattern = 1 << (numOfChar * 2);	// 4^numOfChar
	cachedSaIndex = MMUnitAllocate((numOfPattern+1) * sizeof(unsigned int));
	for (i=0; i<numOfPattern; i++) {
		cachedSaIndex[i] = 0;
	}
	cachedSaIndex[i] = bwt->textLength;


	// Set initial state
	a = 1;
	pos = numOfChar - 1;
	generatedPattern[pos] = (unsigned char)-1;
	saRangeIndex = (unsigned int)-1;

	cf = _mm_load_si128((__m128i*)bwt->cumulativeFreq);	// Load cumulative freq into register
	r1 = _mm_set1_epi32(1);
	cfp1 = _mm_add_epi32(cf, r1);

	saIndexLeft[pos * 4 + 0] = 1 + bwt->cumulativeFreq[0];
	saIndexLeft[pos * 4 + 1] = 1 + bwt->cumulativeFreq[1];
	saIndexLeft[pos * 4 + 2] = 1 + bwt->cumulativeFreq[2];
	saIndexLeft[pos * 4 + 3] = 1 + bwt->cumulativeFreq[3];
	saIndexRight[pos * 4 + 0] = bwt->cumulativeFreq[1];
	saIndexRight[pos * 4 + 1] = bwt->cumulativeFreq[2];
	saIndexRight[pos * 4 + 2] = bwt->cumulativeFreq[3];
	saIndexRight[pos * 4 + 3] = bwt->cumulativeFreq[4];

	while (pos < numOfChar) {
			
		generatedPattern[pos]++;
		if (generatedPattern[pos] >= ALPHABET_SIZE) {
			saRangeIndex &= mask[numOfChar - 1 - pos];
			pos++;
			a >>= 2;
			continue;
		}
		saRangeIndex += a;

		while (pos > 0 && saIndexLeft[pos * 4 + generatedPattern[pos]] <= saIndexRight[pos * 4 + generatedPattern[pos]]) {
			// Set preceding characters to 'a'
			pos--;
			a <<= 2;
			generatedPattern[pos] = 0;

			BWTAllOccValueTwoIndex(bwt, saIndexLeft[(pos+1) * 4 + generatedPattern[pos+1]], saIndexRight[(pos+1) * 4 + generatedPattern[pos+1]] + 1,
										saIndexLeft + pos * 4, saIndexRight + pos * 4);
			// Add cumulative frequency
			r1 = _mm_load_si128((__m128i*)(saIndexLeft + pos * 4));
			r2 = _mm_load_si128((__m128i*)(saIndexRight + pos * 4));
			r1 = _mm_add_epi32(r1, cfp1);
			r2 = _mm_add_epi32(r2, cf);
			_mm_store_si128((__m128i*)(saIndexLeft + pos * 4), r1);
			_mm_store_si128((__m128i*)(saIndexRight + pos * 4), r2);
		}

		if (saIndexLeft[pos * 4 + generatedPattern[pos]] <= saIndexRight[pos * 4 + generatedPattern[pos]]) {
			cachedSaIndex[saRangeIndex] = saIndexLeft[pos * 4 + generatedPattern[pos]];
		}
	}

	for (i=numOfPattern; i>0; i--) {
		if (cachedSaIndex[i-1] == 0) {
			cachedSaIndex[i-1] = cachedSaIndex[i];
		}
	}

	// Adjust the SA index ranges with consecutive 'a' as suffix so that false positive at end of text is more 'sensible'
	// False positive is in the form of X$ where |X| < numChar
	// SA index range is adjusted so that X$ appears as false positive of Xaaa.. where |Xaaa..| = numChar

	saRangeIndex = 0;
	saIndex = 0;

	for (pos=numOfChar-1; pos>0; pos--) {

		saIndex = BWTOccValueOnSpot(bwt, saIndex + 1, &a);
		saIndex += bwt->cumulativeFreq[a];
		saRangeIndex >>= BIT_PER_CHAR;
		saRangeIndex += a << ((numOfChar-1) * BIT_PER_CHAR);

		// Reuse the variables to store saRangeIndex for suffixes shorter than numOfChar
		saIndexLeft[pos] = saRangeIndex;

	}

	// Sort on saRangeIndex
	QSort(saIndexLeft + 1, numOfChar - 1, sizeof(unsigned int), QSortUnsignedIntOrder);

	for (pos=numOfChar-1; pos>0; pos--) {
		i = saIndexLeft[pos];
		a = cachedSaIndex[i];
		while (a > 1 && cachedSaIndex[i] == a) {	// No need to adjust for all 'a' suffixes
			cachedSaIndex[i]--;
			i--;
		}
	}

	fwrite(cachedSaIndex, sizeof(unsigned int), numOfPattern, cachedSaIndexFile);
	MMUnitFree(cachedSaIndex, (numOfPattern+1) * sizeof(unsigned int));

	fclose(cachedSaIndexFile);

}

void BWTGenerateSaBitmap(const BWT *bwt, const unsigned int numOfChar, const char *saBitmapFileName) {

	// stack
	unsigned int saIndexLeft[MAX_ARPROX_MATCH_LENGTH + 1];
	unsigned int saIndexRight[MAX_ARPROX_MATCH_LENGTH];
	unsigned char generatedPattern[MAX_ARPROX_MATCH_LENGTH];
	unsigned int pos;
	unsigned int i, c;
	unsigned LONG numOfPattern;
	unsigned int numOfPatternExistInText = 0;
	unsigned int bitmapSize;
	unsigned int saBitmapIndex;
	unsigned int mask[16] = { 0xFFFFFFFC, 0xFFFFFFF3, 0xFFFFFFCF, 0xFFFFFF3F,
					 0xFFFFFCFF, 0xFFFFF3FF, 0xFFFFCFFF, 0xFFFF3FFF,
					 0xFFFCFFFF, 0xFFF3FFFF, 0xFFCFFFFF, 0xFF3FFFFF,
					 0xFCFFFFFF, 0xF3FFFFFF, 0xCFFFFFFF, 0x3FFFFFFF };

	unsigned int *saBitMap;

	FILE *saBitmapFile;

	saBitmapFile = (FILE*)fopen64(saBitmapFileName, "wb");
	if (saBitmapFile == NULL) {
		fprintf(stderr, "BWTGenerateSaBitmap(): Cannot open SA bitmap file!\n");
		exit(1);
	}

	fwrite(&bwt->inverseSa0, sizeof(unsigned int), 1, saBitmapFile);
	fwrite(bwt->cumulativeFreq + 1, sizeof(unsigned int), ALPHABET_SIZE, saBitmapFile);
	fwrite(&numOfChar, sizeof(unsigned int), 1, saBitmapFile);

	bitmapSize = 1 << (numOfChar * 2 - 5);	// 4^numOfChar
	saBitMap = MMUnitAllocate(bitmapSize * sizeof(unsigned int));
	for (i=0; i<bitmapSize; i++) {
		saBitMap[i] = 0;
	}

	// Set this boundary case so that the last character of generated pattern can be located by backward search
	saIndexLeft[numOfChar] = 0;
	saIndexRight[numOfChar] = bwt->textLength;

	// Set initial state
	pos = numOfChar - 1;
	generatedPattern[pos] = (unsigned char)-1;
	saBitmapIndex = 0;

	while (pos < numOfChar) {
			
		generatedPattern[pos]++;
		if (generatedPattern[pos] >= ALPHABET_SIZE) {
			pos++;
			continue;
		}
		saBitmapIndex &= mask[numOfChar - 1 - pos];
		saBitmapIndex |= generatedPattern[pos] << ((numOfChar - 1 - pos) * BIT_PER_CHAR);

		c = generatedPattern[pos];
		saIndexLeft[pos] = bwt->cumulativeFreq[c] + BWTOccValue(bwt, saIndexLeft[pos+1], c) + 1;
		saIndexRight[pos] = bwt->cumulativeFreq[c] + BWTOccValue(bwt, saIndexRight[pos+1] + 1, c);

		while (pos > 0 && saIndexLeft[pos] <= saIndexRight[pos]) {
			// Set preceding characters to 'a'
			pos--;
			generatedPattern[pos] = 0;
			saBitmapIndex &= mask[numOfChar - 1 - pos];
			c = generatedPattern[pos];
			saIndexLeft[pos] = bwt->cumulativeFreq[c] + BWTOccValue(bwt, saIndexLeft[pos+1], c) + 1;
			saIndexRight[pos] = bwt->cumulativeFreq[c] + BWTOccValue(bwt, saIndexRight[pos+1] + 1, c);
		}

		if (saIndexLeft[pos] <= saIndexRight[pos]) {
			saBitMap[saBitmapIndex >> 5] |= FIRST_BIT_MASK >> (saBitmapIndex & 0x1F);
			numOfPatternExistInText++;
		}
	}

	fwrite(saBitMap, sizeof(unsigned int), bitmapSize, saBitmapFile);
	MMUnitFree(saBitMap, bitmapSize * sizeof(unsigned int));

	numOfPattern = (unsigned LONG)1 << (numOfChar * 2);	// 4^numOfChar

	printf("%u out of %llu length %u patterns exist in text.\n", numOfPatternExistInText, numOfPattern, numOfChar);

	fclose(saBitmapFile);

}

void BWTGenerateCompressedSaBitmap(const BWT *bwt, const unsigned int numOfChar, const char *saBitmapFileName) {

	// stack
	unsigned int saIndexLeft[MAX_ARPROX_MATCH_LENGTH + 1];
	unsigned int saIndexRight[MAX_ARPROX_MATCH_LENGTH];
	unsigned char generatedPattern[MAX_ARPROX_MATCH_LENGTH];
	unsigned int pos;
	unsigned int i, j, k, c;
	unsigned LONG numOfPattern;
	unsigned int numOfPatternExistInText = 0;
	unsigned int bitmapSize;
	unsigned LONG saBitmapIndex;
	unsigned LONG mask[32];

	unsigned int *saBitMap;

	FILE *saBitmapFile;

	unsigned LONG numOfBit = 0, numOfUnaryBit = 0, numOfGammaBit = 0;
	unsigned LONG numOfOneUnaryBit = 0, numOfOneGammaBit = 0;
	unsigned LONG numOfZeroUnaryBit = 0, numOfZeroGammaBit = 0;
	unsigned int lastBitValue, bitValue;
	unsigned int groupSize;

	unsigned int tempIntegerStream[1025];
	unsigned int numOfInteger;
	unsigned int numOfHighBit, numOfLowBit;

	unsigned int tempOneStream[1025];
	unsigned int tempZeroStream[1025];
	unsigned int numOfOne, numOfZero;


	for (i=0; i<32; i++) {
		mask[i] = ~((unsigned LONG)(0x3) << (i*2));
	}


	saBitmapFile = (FILE*)fopen64(saBitmapFileName, "wb");
	if (saBitmapFile == NULL) {
		fprintf(stderr, "BWTGenerateCompressedSaBitmap(): Cannot open SA compressed bitmap file!\n");
		exit(1);
	}

	fwrite(&bwt->inverseSa0, sizeof(unsigned int), 1, saBitmapFile);
	fwrite(bwt->cumulativeFreq + 1, sizeof(unsigned int), ALPHABET_SIZE, saBitmapFile);
	fwrite(&numOfChar, sizeof(unsigned int), 1, saBitmapFile);


	bitmapSize = 1 << (numOfChar * 2 - 5);	// 4^numOfChar
	saBitMap = MMUnitAllocate(bitmapSize * sizeof(unsigned int));
	for (i=0; i<bitmapSize; i++) {
		saBitMap[i] = 0;
	}

	// Set this boundary case so that the last character of generated pattern can be located by backward search
	saIndexLeft[numOfChar] = 0;
	saIndexRight[numOfChar] = bwt->textLength;

	// Set initial state
	pos = numOfChar - 1;
	generatedPattern[pos] = (unsigned char)-1;
	saBitmapIndex = 0;

	while (pos < numOfChar) {
			
		generatedPattern[pos]++;
		if (generatedPattern[pos] >= ALPHABET_SIZE) {
			pos++;
			continue;
		}
		saBitmapIndex &= mask[numOfChar - 1 - pos];
		saBitmapIndex |= (unsigned LONG)generatedPattern[pos] << ((numOfChar - 1 - pos) * BIT_PER_CHAR);

		c = generatedPattern[pos];
		saIndexLeft[pos] = bwt->cumulativeFreq[c] + BWTOccValue(bwt, saIndexLeft[pos+1], c) + 1;
		saIndexRight[pos] = bwt->cumulativeFreq[c] + BWTOccValue(bwt, saIndexRight[pos+1] + 1, c);

		while (pos > 0 && saIndexLeft[pos] <= saIndexRight[pos]) {
			// Set preceding characters to 'a'
			pos--;
			generatedPattern[pos] = 0;
			saBitmapIndex &= mask[numOfChar - 1 - pos];
			c = generatedPattern[pos];
			saIndexLeft[pos] = bwt->cumulativeFreq[c] + BWTOccValue(bwt, saIndexLeft[pos+1], c) + 1;
			saIndexRight[pos] = bwt->cumulativeFreq[c] + BWTOccValue(bwt, saIndexRight[pos+1] + 1, c);
		}

		if (saIndexLeft[pos] <= saIndexRight[pos]) {
			saBitMap[saBitmapIndex >> 5] |= FIRST_BIT_MASK >> ((unsigned int)saBitmapIndex & 0x1F);
			numOfPatternExistInText++;
		}
	}

	// Compress bitmap
	lastBitValue = saBitMap[0] >> (BITS_IN_WORD - 1);
	groupSize = 0;
	numOfInteger = 0;
	tempIntegerStream[0] = 0;

	tempOneStream[0] = 0;
	tempZeroStream[0] = 0;
	numOfOne = 0;
	numOfZero = 0;

	for (i=0; i<bitmapSize; i++) {
		c = saBitMap[i];
		for (j=0; j<BITS_IN_WORD; j++) {
			bitValue = c >> (BITS_IN_WORD - 1);
			c <<= 1;
			if (bitValue == lastBitValue) {
				groupSize++;
			} else {
				if (numOfInteger >= 1024) {
					for (k=1; k<=1024; k++) {
						tempIntegerStream[k] += tempIntegerStream[k-1];
					}
					numOfHighBit = 10;
					numOfLowBit = ceilLog2(tempIntegerStream[1024]) - numOfHighBit;
					for (k=1; k<1024; k++) {
						tempIntegerStream[k] >>= numOfLowBit;
						numOfUnaryBit += tempIntegerStream[k] - tempIntegerStream[k-1] + 1;
						numOfGammaBit += 2 * floorLog2(tempIntegerStream[k] - tempIntegerStream[k-1] + 1) + 1;
					}
					numOfUnaryBit += numOfLowBit * 1024;
					numOfGammaBit += numOfLowBit * 1024;
					numOfInteger = 0;
				}
				tempIntegerStream[numOfInteger+1] = groupSize;
				numOfInteger++;

				// One
				if (lastBitValue == 1) {
					if (numOfOne >= 1024) {
						for (k=1; k<=1024; k++) {
							tempOneStream[k] += tempOneStream[k-1];
						}
						numOfHighBit = 10;
						numOfLowBit = ceilLog2(tempOneStream[1024]) - numOfHighBit;
						for (k=1; k<1024; k++) {
							tempOneStream[k] >>= numOfLowBit;
							numOfOneUnaryBit += tempOneStream[k] - tempOneStream[k-1] + 1;
							numOfOneGammaBit += 2 * floorLog2(tempOneStream[k] - tempOneStream[k-1] + 1) + 1;
						}
						numOfOneUnaryBit += numOfLowBit * 1024;
						numOfOneGammaBit += numOfLowBit * 1024;
						numOfOne = 0;
					}
					tempOneStream[numOfOne+1] = groupSize;
					numOfOne++;
				} else {
					// Zero
					if (numOfZero >= 1024) {
						for (k=1; k<=1024; k++) {
							tempZeroStream[k] += tempZeroStream[k-1];
						}
						numOfHighBit = 10;
						numOfLowBit = ceilLog2(tempZeroStream[1024]) - numOfHighBit;
						for (k=1; k<1024; k++) {
							tempZeroStream[k] >>= numOfLowBit;
							numOfZeroUnaryBit += tempZeroStream[k] - tempZeroStream[k-1] + 1;
							numOfZeroGammaBit += 2 * floorLog2(tempZeroStream[k] - tempZeroStream[k-1] + 1) + 1;
						}
						numOfZeroUnaryBit += numOfLowBit * 1024;
						numOfZeroGammaBit += numOfLowBit * 1024;
						numOfZero = 0;
					}
					tempZeroStream[numOfZero+1] = groupSize;
					numOfZero++;
				}


				numOfBit += 2 * floorLog2(groupSize) + 1;
				groupSize = 1;
				lastBitValue = bitValue;
			}
		}
	}

	printf("Gamma only.\n");
	printf("Compressed bitmap size = %llu bits.\n", numOfBit);
	printf("Compressed bitmap size = %.2f bytes.\n", (double)numOfBit / BITS_IN_BYTE);

	printf("Rice + unary.\n");
	printf("Compressed bitmap size = %llu bits.\n", numOfUnaryBit);
	printf("Compressed bitmap size = %.2f bytes.\n", (double)numOfUnaryBit / BITS_IN_BYTE);

	printf("Rice + gamma.\n");
	printf("Compressed bitmap size = %llu bits.\n", numOfGammaBit);
	printf("Compressed bitmap size = %.2f bytes.\n", (double)numOfGammaBit / BITS_IN_BYTE);

	printf("Separate Ones and Zeros Rice + unary.\n");
	printf("Compressed bitmap size = %llu bits.\n", numOfOneUnaryBit + numOfZeroUnaryBit);
	printf("Compressed bitmap size = %.2f bytes.\n", (double)(numOfOneUnaryBit + numOfZeroUnaryBit) / BITS_IN_BYTE);

	printf("Separate Ones and Zeros Rice + gamma.\n");
	printf("Compressed bitmap size = %llu bits.\n", numOfOneGammaBit + numOfZeroGammaBit);
	printf("Compressed bitmap size = %.2f bytes.\n", (double)(numOfOneGammaBit + numOfZeroGammaBit) / BITS_IN_BYTE);

	fwrite(saBitMap, sizeof(unsigned int), bitmapSize, saBitmapFile);
	MMUnitFree(saBitMap, bitmapSize * sizeof(unsigned int));

	numOfPattern = (unsigned LONG)1 << (numOfChar * 2);	// 4^numOfChar

	printf("%u out of %llu length %u patterns exist in text.\n", numOfPatternExistInText, numOfPattern, numOfChar);

	fclose(saBitmapFile);

}

void BWTCountPattern(const BWT *bwt, const unsigned int numOfChar) {

	// stack
	unsigned int saIndexLeft[MAX_ARPROX_MATCH_LENGTH + 1];
	unsigned int saIndexRight[MAX_ARPROX_MATCH_LENGTH];
	unsigned char generatedPattern[MAX_ARPROX_MATCH_LENGTH];
	unsigned int pos;
	unsigned int c;
	unsigned LONG numOfPattern;
	unsigned int numOfPatternExistInText = 0;

	// Set this boundary case so that the last character of generated pattern can be located by backward search
	saIndexLeft[numOfChar] = 0;
	saIndexRight[numOfChar] = bwt->textLength;

	// Set initial state
	pos = numOfChar - 1;
	generatedPattern[pos] = (unsigned char)-1;

	while (pos < numOfChar) {
			
		generatedPattern[pos]++;
		if (generatedPattern[pos] >= ALPHABET_SIZE) {
			pos++;
			continue;
		}

		c = generatedPattern[pos];
		saIndexLeft[pos] = bwt->cumulativeFreq[c] + BWTOccValue(bwt, saIndexLeft[pos+1], c) + 1;
		saIndexRight[pos] = bwt->cumulativeFreq[c] + BWTOccValue(bwt, saIndexRight[pos+1] + 1, c);

		while (pos > 0 && saIndexLeft[pos] <= saIndexRight[pos]) {
			// Set preceding characters to 'a'
			pos--;
			generatedPattern[pos] = 0;
			c = generatedPattern[pos];
			saIndexLeft[pos] = bwt->cumulativeFreq[c] + BWTOccValue(bwt, saIndexLeft[pos+1], c) + 1;
			saIndexRight[pos] = bwt->cumulativeFreq[c] + BWTOccValue(bwt, saIndexRight[pos+1] + 1, c);
		}

		if (saIndexLeft[pos] <= saIndexRight[pos]) {
			numOfPatternExistInText++;
		}
	}

	numOfPattern = (unsigned LONG)1 << (numOfChar * 2);	// 4^numOfChar

	printf("%u out of %llu length %u patterns exist in text.\n", numOfPatternExistInText, numOfPattern, numOfChar);

}


/*
void BWTGenerateSaRangeTable(const BWT *bwt, const unsigned int numOfChar, const char *PackedDNAFileName, const char *saRangeFileName, const char *startSaIndexFileName) {

	// stack
	unsigned int saIndexLeft[MAX_ARPROX_MATCH_LENGTH + 1];
	unsigned int saIndexRight[MAX_ARPROX_MATCH_LENGTH];
	unsigned char generatedPattern[MAX_ARPROX_MATCH_LENGTH];
	unsigned int endOfTextPattern[MAX_ARPROX_MATCH_LENGTH];
	unsigned int endOfTextChar;
	unsigned int pos;
	unsigned int i, c;
	unsigned int numOfPattern;
	unsigned int saRangeIndex, lastSaRangeIndex;
	unsigned int mask[16] = { 0xFFFFFFFC, 0xFFFFFFF3, 0xFFFFFFCF, 0xFFFFFF3F,
					 0xFFFFFCFF, 0xFFFFF3FF, 0xFFFFCFFF, 0xFFFF3FFF,
					 0xFFFCFFFF, 0xFFF3FFFF, 0xFFCFFFFF, 0xFF3FFFFF,
					 0xFCFFFFFF, 0xF3FFFFFF, 0xCFFFFFFF, 0x3FFFFFFF };

	SaIndexRange *saIndexRange;
	unsigned int *startSaIndex;
	
	FILE *saRangeFile, *startSaIndexFile;
	FILE *packedDNAFile;

	unsigned int packedFileLen;
	unsigned char lastByteLength, dnaByte;
	unsigned int d;

	packedDNAFile = (FILE*)fopen64(PackedDNAFileName, "rb");
	if (packedDNAFile == NULL) {
		fprintf(stderr, "BWTGenerateSaRangeTable(): Cannot open packed DNA file!\n");
		exit(1);
	}
	fseek(packedDNAFile, -1, SEEK_END);
	packedFileLen = ftell(packedDNAFile);
	if ((int)packedFileLen < 0) {
		fprintf(stderr, "BWTGenerateSaRangeTable(): Cannot determine file length!\n");
		exit(1);
	}
	fread(&lastByteLength, sizeof(unsigned char), 1, packedDNAFile);

	fseek(packedDNAFile, -2, SEEK_CUR);
	fread(&dnaByte, sizeof(unsigned char), 1, packedDNAFile);
	d = dnaByte;
	endOfTextChar = (d >> ((CHAR_PER_BYTE - lastByteLength) * BIT_PER_CHAR));
	i = lastByteLength;

	while (i<numOfChar-1) {
		fseek(packedDNAFile, -2, SEEK_CUR);
		fread(&dnaByte, sizeof(unsigned char), 1, packedDNAFile);
		d =	dnaByte;
		endOfTextChar |= (d << (i * BIT_PER_CHAR));
		i += CHAR_PER_BYTE;
	}

	endOfTextChar &= truncateLeft(ALL_ONE_MASK, (numOfChar - 1) * BIT_PER_CHAR);

	fclose(packedDNAFile);


	// saRangeFile
	saRangeFile = (FILE*)fopen64(saRangeFileName, "wb");
	if (saRangeFile == NULL) {
		fprintf(stderr, "BWTGenerateSaRangeTable(): Cannot open SA range file!\n");
		exit(1);
	}

	fwrite(&bwt->inverseSa0, sizeof(unsigned int), 1, saRangeFile);
	fwrite(bwt->cumulativeFreq + 1, sizeof(unsigned int), ALPHABET_SIZE, saRangeFile);
	fwrite(&numOfChar, sizeof(unsigned int), 1, saRangeFile);

	// startSaIndexFile
	startSaIndexFile = (FILE*)fopen64(startSaIndexFileName, "wb");
	if (startSaIndexFile == NULL) {
		fprintf(stderr, "BWTGenerateSaRangeTable(): Cannot open Start SA index file!\n");
		exit(1);
	}

	fwrite(&bwt->inverseSa0, sizeof(unsigned int), 1, startSaIndexFile);
	fwrite(bwt->cumulativeFreq + 1, sizeof(unsigned int), ALPHABET_SIZE, startSaIndexFile);
	fwrite(&numOfChar, sizeof(unsigned int), 1, startSaIndexFile);


	numOfPattern = 1 << (numOfChar * 2);	// 4^numOfChar
	saIndexRange = MMUnitAllocate(numOfPattern * sizeof(SaIndexRange));
	for (i=0; i<numOfPattern; i++) {
		saIndexRange[i].startSaIndex = 0;
		saIndexRange[i].endSaIndex = 0;
	}

	// Set this boundary case so that the last character of generated pattern can be located by backward search
	saIndexLeft[numOfChar] = 0;
	saIndexRight[numOfChar] = bwt->textLength;

	// Set initial state
	pos = numOfChar - 1;
	generatedPattern[pos] = (unsigned char)-1;
	saRangeIndex = 0;

	while (pos < numOfChar) {
			
		generatedPattern[pos]++;
		if (generatedPattern[pos] >= ALPHABET_SIZE) {
			pos++;
			continue;
		}
		saRangeIndex &= mask[pos];
		saRangeIndex |= generatedPattern[pos] << (pos * BIT_PER_CHAR);

		c = generatedPattern[pos];
		saIndexLeft[pos] = bwt->cumulativeFreq[c] + BWTOccValue(bwt, saIndexLeft[pos+1], c) + 1;
		saIndexRight[pos] = bwt->cumulativeFreq[c] + BWTOccValue(bwt, saIndexRight[pos+1] + 1, c);

		while (pos > 0 && saIndexLeft[pos] <= saIndexRight[pos]) {
			// Set preceding characters to 'a'
			pos--;
			generatedPattern[pos] = 0;
			saRangeIndex &= mask[pos];
			c = generatedPattern[pos];
			saIndexLeft[pos] = bwt->cumulativeFreq[c] + BWTOccValue(bwt, saIndexLeft[pos+1], c) + 1;
			saIndexRight[pos] = bwt->cumulativeFreq[c] + BWTOccValue(bwt, saIndexRight[pos+1] + 1, c);
		}

		if (saIndexLeft[pos] <= saIndexRight[pos]) {
			saIndexRange[saRangeIndex].startSaIndex = saIndexLeft[pos];
			saIndexRange[saRangeIndex].endSaIndex = saIndexRight[pos];
		}

	}

	// Output saIndexRange
	fwrite(saIndexRange, sizeof(SaIndexRange), numOfPattern, saRangeFile);
	fclose(saRangeFile);

	MMUnitFree(saIndexRange, numOfPattern * sizeof(SaIndexRange));


	// Generate start SA index

	// Search again with normal order of search string

	for (i=0; i<numOfPattern; i++) {
		saIndexRange[i].startSaIndex = 0;
		saIndexRange[i].endSaIndex = 0;
	}

	// Set this boundary case so that the last character of generated pattern can be located by backward search
	saIndexLeft[numOfChar] = 0;
	saIndexRight[numOfChar] = bwt->textLength;

	// Set initial state
	pos = numOfChar - 1;
	generatedPattern[pos] = (unsigned char)-1;
	saRangeIndex = 0;

	while (pos < numOfChar) {
			
		generatedPattern[pos]++;
		if (generatedPattern[pos] >= ALPHABET_SIZE) {
			pos++;
			continue;
		}
		saRangeIndex &= mask[numberOfChar - 1 - pos];
		saRangeIndex |= generatedPattern[pos] << ((numberOfChar - 1 - pos) * BIT_PER_CHAR);

		c = generatedPattern[pos];
		saIndexLeft[pos] = bwt->cumulativeFreq[c] + BWTOccValue(bwt, saIndexLeft[pos+1], c) + 1;
		saIndexRight[pos] = bwt->cumulativeFreq[c] + BWTOccValue(bwt, saIndexRight[pos+1] + 1, c);

		while (pos > 0 && saIndexLeft[pos] <= saIndexRight[pos]) {
			// Set preceding characters to 'a'
			pos--;
			generatedPattern[pos] = 0;
			saRangeIndex &= mask[numberOfChar - 1 - pos];
			c = generatedPattern[pos];
			saIndexLeft[pos] = bwt->cumulativeFreq[c] + BWTOccValue(bwt, saIndexLeft[pos+1], c) + 1;
			saIndexRight[pos] = bwt->cumulativeFreq[c] + BWTOccValue(bwt, saIndexRight[pos+1] + 1, c);
		}

		if (saIndexLeft[pos] <= saIndexRight[pos]) {
			saIndexRange[saRangeIndex].startSaIndex = saIndexLeft[pos];
			saIndexRange[saRangeIndex].endSaIndex = saIndexRight[pos];
		}

	}

	startSaIndex = MMUnitAllocate(sizeof(unsigned int) * (numOfPattern + 1));
	for (i=0; i<numOfPattern; i++) {
		saRangeIndex = 0;
		// reverse character order
		for (pos=0; pos<numOfChar; pos++) {
			saRangeIndex |= (((i >> pos) & 0x3) << (numOfChar - pos));
		}
		startSaIndex[saRangeIndex] = saIndexRange[i].startSaIndex;
	}
	startSaIndex[numOfPattern] = bwt->textLength + 1;

	// Fill startSaIndex for pattern not found in text
	for (i=numOfPattern-1; i>0; i--) {
		if (startSaIndex[i] == 0) {
			startSaIndex[i] = startSaIndex[i+1];
		}
	}

	// Find the suffices in text with length < numOfChar

	// $
	saRangeIndex = 0;
	saIndexRange[saRangeIndex].startSaIndex--;
	endOfTextPattern[numOfChar - 1] = saRangeIndex;

	for (i=numOfChar - 1; i>0; i--) {
		for (c=0; c<ALPHABET_SIZE; c++) {
			saRangeIndex &= mask[i];
			saRangeIndex |= c << (i * BIT_PER_CHAR);
			if (saIndexRange[saRangeIndex].startSaIndex != 0 || saIndexRange[saRangeIndex].endSaIndex) {
				if (saRangeIndex != 0) {
					lastSaRangeIndex = saRangeIndex;
					pos = numOfChar - 1;
					while ((saRangeIndex & (~mask[pos])) == 0) {
						lastSaRangeIndex |= (ALPHABET_SIZE - 1) << (pos * BIT_PER_CHAR);
						pos--;
					}
					lastSaRangeIndex -= (1 << (pos * BIT_PER_CHAR));
					if (saIndexRange[saRangeIndex].startSaIndex != saIndexRange[lastSaRangeIndex].endSaIndex + 1) {
						break;
					}
				} else {
					if (saIndexRange[saRangeIndex].startSaIndex > 0) {
						break;
					}
				}
			}
		}
		saIndexRange[saRangeIndex].startSaIndex--;
		endOfTextPattern[i] = saRangeIndex;
	}

	// Output short suffices
	fwrite(endOfTextPattern, sizeof(unsigned int), numOfChar, saRangeFile);

	// Output startSaIndex
	for (i=0; i<numOfPattern; i++) {
		fwrite(&(saIndexRange[i].startSaIndex), sizeof(unsigned int), 1, saRangeFile);
	}
	MMUnitFree(saIndexRange, numOfPattern * sizeof(SaIndexRange));

	fclose(saRangeFile);

}
*/

void BWTPrintAndVerifyOccValue(const BWT *bwt, FILE *output) {

	unsigned int i, j;
	unsigned int c;
	unsigned int v;
	unsigned int occValue[ALPHABET_SIZE];

	for (i=0; i<ALPHABET_SIZE; i++) {
		occValue[i] = 0;
	}

	for (i=1; i<=bwt->textLength+1; i++) {
		if (i != bwt->inverseSa0) {
			v = BWTOccValueOnSpot(bwt, i, &c);
			if (output != NULL) {
				fprintf(output, "%u:  %u : %u\n", i, c, v);
			}
			occValue[c] = v;
			v = 0;
			for (j=0; j<ALPHABET_SIZE; j++) {
				v += occValue[j];
			}
			if ((i > bwt->inverseSa0 && v != (i-1)) || (i < bwt->inverseSa0 && v != i)) {
				fprintf(stderr, "BWTPrintAndVerifyOccValue(): occValue not summed to text length!\n");
				exit(1);
			}
		}
	}

}

void BWTPrintAndVerifySaValue(const BWT *bwt, FILE *output) {

	unsigned int i;
	unsigned int saValue;
	unsigned int *printedValue;

	printedValue = malloc(sizeof(unsigned int) * (bwt->textLength + 1));
	initializeVAL(printedValue, bwt->textLength + 1, FALSE);
    
	for (i=0; i<=bwt->textLength; i++) {
		saValue = BWTSaValue(bwt, i);
		if (output != NULL) {
			fprintf(output, "%u ", saValue);
		}
		printedValue[saValue] = TRUE;
	}

	if (output != NULL) {
		fprintf(output, "\n");
	}

	if (numberOfMatchInVAL(printedValue, bwt->textLength + 1, FALSE) > 0) {
		fprintf(stderr, "BWTPrintAndVerifySaValue() : some saValue is not printed!\n");
		exit(1);
	}

	free(printedValue);

}

void BWTPrintAndVerifyInverseSa(const BWT *bwt, FILE *output) {

	unsigned int i;
	unsigned int inverseSa;
	unsigned int *printedValue;

	printedValue = malloc(sizeof(unsigned int) * (bwt->textLength + 1));
	initializeVAL(printedValue, bwt->textLength + 1, FALSE);
    
	for (i=0; i<bwt->textLength+1; i++) {
		inverseSa = BWTInverseSa(bwt, i);
		if (output != NULL) {
			fprintf(output, "%u ", inverseSa);
		}
		printedValue[inverseSa] = TRUE;
	}

	if (output != NULL) {
		fprintf(output, "\n");
	}

	if (numberOfMatchInVAL(printedValue, bwt->textLength + 1, FALSE) > 0) {
		fprintf(stderr, "BWTPrintAndVerifyInverseSa() : some inverseSa is not printed!\n");
		exit(1);
	}

	free(printedValue);

}

#else
/*

   BWTConstruct.c		BWT-Index Construction

   This module constructs BWT and auxiliary data structures.

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
#include "BWTConstruct.h"
#include "MiscUtilities.h"
#include "DNACount.h"
#include "TextConverter.h"
#include "MemManager.h"
#include "QSufSort.h"
#include "r250.h"

// Static functions
static void BWTIncConstruct(BWTInc *bwtInc, const unsigned int numChar);

static void BWTIncSetBuildSizeAndTextAddr(BWTInc *bwtInc);
static void BWTIncPutPackedTextToRank(const unsigned int *packedText, unsigned int* __restrict rank, unsigned int* __restrict cumulativeCount, const unsigned int numChar);
static void BWTIncBuildPackedBwt(const unsigned int *relativeRank, unsigned int* __restrict bwt, const unsigned int numChar,
								 const unsigned int *cumulativeCount, const unsigned int *packedShift);
static unsigned int BWTIncGetAbsoluteRank(BWT *bwt, unsigned int* __restrict absoluteRank, unsigned int* __restrict seq, const unsigned int *packedText, 
								  const unsigned int numChar, const unsigned int* cumulativeCount, const unsigned int firstCharInLastIteration);
static void BWTIncSortKey(unsigned int* __restrict key, unsigned int* __restrict seq, const unsigned int numItem);
static void BWTIncBuildRelativeRank(unsigned int* __restrict sortedRank, unsigned int* __restrict seq, unsigned int* __restrict relativeRank, const unsigned int numItem, 
									unsigned int oldInverseSa0, const unsigned int *cumulativeCount);
static void BWTIncBuildBwt(unsigned int*  seq, const unsigned int *relativeRank, const unsigned int numChar, const unsigned int *cumulativeCount);	// seq is replaced with Bwt
static void BWTIncMergeBwt(const unsigned int *sortedRank, const unsigned int* oldBwt, const unsigned int *insertBwt, unsigned int* __restrict mergedBwt, 
						   const unsigned int numOldBwt, const unsigned int numInsertBwt);


BWTInc *BWTIncCreate(MMPool *mmPool, const unsigned int textLength, const float targetNBit,
					 const unsigned int initialMaxBuildSize, const unsigned int incMaxBuildSize) {

	BWTInc *bwtInc;
	unsigned int i;

	if (targetNBit == 0) {
		fprintf(stderr, "BWTIncCreate() : targetNBit = 0!\n");
		exit(1);
	}
	
	bwtInc = MMPoolDispatch(mmPool, sizeof(BWTInc));
	
	bwtInc->numberOfIterationDone = 0;

	bwtInc->bwt = BWTCreate(mmPool, textLength, NULL);

	bwtInc->initialMaxBuildSize = initialMaxBuildSize;
	bwtInc->incMaxBuildSize = incMaxBuildSize;

	bwtInc->targetNBit = targetNBit;

	bwtInc->cumulativeCountInCurrentBuild = MMPoolDispatch(mmPool, sizeof(unsigned int) * (ALPHABET_SIZE + 1));
	initializeVAL(bwtInc->cumulativeCountInCurrentBuild, ALPHABET_SIZE + 1, 0);

	// Build frequently accessed data
	bwtInc->packedShift = MMPoolDispatch(mmPool, sizeof(unsigned int) * CHAR_PER_WORD);
	for (i=0; i<CHAR_PER_WORD; i++) {
		bwtInc->packedShift[i] = BITS_IN_WORD - (i+1) * BIT_PER_CHAR;
	}

	bwtInc->targetTextLength = textLength;
	bwtInc->availableWord = (unsigned int)((textLength + OCC_INTERVAL - 1) / OCC_INTERVAL * OCC_INTERVAL / BITS_IN_WORD * bwtInc->targetNBit);
	if (bwtInc->availableWord < BWTResidentSizeInWord(textLength) + BWTOccValueMinorSizeInWord(textLength)) {
		fprintf(stderr, "BWTIncCreate() : targetNBit is too low!\n");
		exit(1);
	}
	bwtInc->workingMemory = MMUnitAllocate(bwtInc->availableWord * BYTES_IN_WORD);

	return bwtInc;

}

void BWTIncFree(MMPool *mmPool, BWTInc *bwtInc) {

	MMUnitFree(bwtInc->workingMemory, bwtInc->availableWord * BYTES_IN_WORD);
	MMPoolReturn(mmPool, bwtInc, sizeof(BWTInc));

}

BWTInc *BWTIncConstructFromPacked(MMPool *mmPool, const char *inputFileName, const unsigned int showProgress,
								  const float targetNBit, const unsigned int initialMaxBuildSize, const unsigned int incMaxBuildSize) {

	FILE *packedFile;
	unsigned int packedFileLen;
	unsigned int totalTextLength;
	unsigned int textToLoad, textSizeInByte;
	unsigned int processedTextLength;
	unsigned char lastByteLength;

	BWTInc *bwtInc;

	packedFile = (FILE*)fopen64(inputFileName, "rb");

	if (packedFile == NULL) {
		fprintf(stderr, "BWTIncConstructFromPacked() : Cannot open inputFileName!\n");
		exit(1);
	}

	fseek(packedFile, -1, SEEK_END);
	packedFileLen = ftell(packedFile);
	if ((int)packedFileLen < 0) {
		fprintf(stderr, "BWTIncConstructFromPacked: Cannot determine file length!\n");
		exit(1);
	}
	fread(&lastByteLength, sizeof(unsigned char), 1, packedFile);
	totalTextLength = TextLengthFromBytePacked(packedFileLen, BIT_PER_CHAR, lastByteLength);

	bwtInc = BWTIncCreate(mmPool, totalTextLength, targetNBit, initialMaxBuildSize, incMaxBuildSize);

	BWTIncSetBuildSizeAndTextAddr(bwtInc);

	if (bwtInc->buildSize > totalTextLength) {
		textToLoad = totalTextLength;
	} else {
		textToLoad = totalTextLength - ((totalTextLength - bwtInc->buildSize + CHAR_PER_WORD - 1) / CHAR_PER_WORD * CHAR_PER_WORD);
	}
	textSizeInByte = textToLoad / CHAR_PER_BYTE;	// excluded the odd byte

	fseek(packedFile, -2, SEEK_CUR);
	fseek(packedFile, -((int)textSizeInByte), SEEK_CUR);
	fread(bwtInc->textBuffer, sizeof(unsigned char), textSizeInByte + 1, packedFile);
	fseek(packedFile, -((int)textSizeInByte + 1), SEEK_CUR);

	ConvertBytePackedToWordPacked(bwtInc->textBuffer, bwtInc->packedText, ALPHABET_SIZE, textToLoad);
	BWTIncConstruct(bwtInc, textToLoad);

	processedTextLength = textToLoad;

	while (processedTextLength < totalTextLength) {
		textToLoad = bwtInc->buildSize / CHAR_PER_WORD * CHAR_PER_WORD;
		if (textToLoad > totalTextLength - processedTextLength) {
			textToLoad = totalTextLength - processedTextLength;
		}
		textSizeInByte = textToLoad / CHAR_PER_BYTE;
		fseek(packedFile, -((int)textSizeInByte), SEEK_CUR);
		fread(bwtInc->textBuffer, sizeof(unsigned char), textSizeInByte, packedFile);
		fseek(packedFile, -((int)textSizeInByte), SEEK_CUR);
		ConvertBytePackedToWordPacked(bwtInc->textBuffer, bwtInc->packedText, ALPHABET_SIZE, textToLoad);
		BWTIncConstruct(bwtInc, textToLoad);
		processedTextLength += textToLoad;
		if (showProgress && bwtInc->numberOfIterationDone % 10 == 0) {
			printf("%u iterations done. %u characters processed.\n", bwtInc->numberOfIterationDone, processedTextLength);
		}
	}

	return bwtInc;

}

static void BWTIncConstruct(BWTInc *bwtInc, const unsigned int numChar) {

	unsigned int i;
	unsigned int mergedBwtSizeInWord, mergedOccSizeInWord;
	unsigned int firstCharInThisIteration;

	unsigned int *relativeRank, *seq, *sortedRank, *insertBwt, *mergedBwt;
	unsigned int newInverseSa0RelativeRank, oldInverseSa0RelativeRank, newInverseSa0;

	#ifdef DEBUG
	if (numChar > bwtInc->buildSize) {
		fprintf(stderr, "BWTIncConstruct(): numChar > buildSize!\n");
		exit(1);
	}
	#endif

	mergedBwtSizeInWord = BWTResidentSizeInWord(bwtInc->bwt->textLength + numChar);
	mergedOccSizeInWord = BWTOccValueMinorSizeInWord(bwtInc->bwt->textLength + numChar);

	initializeVAL(bwtInc->cumulativeCountInCurrentBuild, ALPHABET_SIZE + 1, 0);

	if (bwtInc->bwt->textLength == 0) {		// Initial build

		// Set address
		seq = bwtInc->workingMemory;
		relativeRank = seq + bwtInc->buildSize + 1;
		mergedBwt = insertBwt = bwtInc->workingMemory + bwtInc->availableWord - mergedBwtSizeInWord;	// build in place

		BWTIncPutPackedTextToRank(bwtInc->packedText, relativeRank, bwtInc->cumulativeCountInCurrentBuild, numChar);

		firstCharInThisIteration = relativeRank[0];
		relativeRank[numChar] = 0;

		// Sort suffix
		QSufSortSuffixSort((int*)relativeRank, (int*)seq, (int)numChar, (int)ALPHABET_SIZE - 1, 0, FALSE);
		newInverseSa0 = relativeRank[0];

		// Clear BWT area
		initializeVAL(insertBwt, mergedBwtSizeInWord, 0);

		// Build BWT
		BWTIncBuildPackedBwt(relativeRank, insertBwt, numChar, bwtInc->cumulativeCountInCurrentBuild, bwtInc->packedShift);

		// so that the cumulativeCount is not deducted
		bwtInc->firstCharInLastIteration = ALPHABET_SIZE;

	} else {		// Incremental build

#ifdef DEBUG
		if (numChar / CHAR_PER_WORD * CHAR_PER_WORD != numChar) {
			fprintf(stderr, "BWTIncConstruct(): numChar must be multiple of char per word!\n");
			exit(1);
		}
#endif

		// Set address
		sortedRank = bwtInc->workingMemory;
		seq = sortedRank + bwtInc->buildSize + 1;
		insertBwt = seq;
		relativeRank = seq + bwtInc->buildSize + 1;

		// Store the first character of this iteration
		firstCharInThisIteration = bwtInc->packedText[0] >> (BITS_IN_WORD - BIT_PER_CHAR);

		// Count occurrence of input text
		ForwardDNAAllOccCountNoLimit(bwtInc->packedText, numChar, bwtInc->cumulativeCountInCurrentBuild + 1, bwtInc->bwt->decodeTable);
		// Add the first character of the previous iteration to represent the inverseSa0 of the previous iteration
		bwtInc->cumulativeCountInCurrentBuild[bwtInc->firstCharInLastIteration + 1]++;
		bwtInc->cumulativeCountInCurrentBuild[2] += bwtInc->cumulativeCountInCurrentBuild[1];
		bwtInc->cumulativeCountInCurrentBuild[3] += bwtInc->cumulativeCountInCurrentBuild[2];
		bwtInc->cumulativeCountInCurrentBuild[4] += bwtInc->cumulativeCountInCurrentBuild[3];

		// Get rank of new suffix among processed suffix
		// The seq array is built into ALPHABET_SIZE + 2 groups; ALPHABET_SIZE groups + 1 group divided into 2 by inverseSa0 + inverseSa0 as 1 group
		oldInverseSa0RelativeRank = BWTIncGetAbsoluteRank(bwtInc->bwt, sortedRank, seq, bwtInc->packedText, 
														  numChar, bwtInc->cumulativeCountInCurrentBuild, bwtInc->firstCharInLastIteration);

		// Sort rank by ALPHABET_SIZE + 2 groups (or ALPHABET_SIZE + 1 groups when inverseSa0 sit on the border of a group)
		for (i=0; i<ALPHABET_SIZE; i++) {
			if (bwtInc->cumulativeCountInCurrentBuild[i] > oldInverseSa0RelativeRank ||
				bwtInc->cumulativeCountInCurrentBuild[i+1] <= oldInverseSa0RelativeRank) {
				BWTIncSortKey(sortedRank + bwtInc->cumulativeCountInCurrentBuild[i], seq + bwtInc->cumulativeCountInCurrentBuild[i], bwtInc->cumulativeCountInCurrentBuild[i+1] - bwtInc->cumulativeCountInCurrentBuild[i]);
			} else {
				if (bwtInc->cumulativeCountInCurrentBuild[i] < oldInverseSa0RelativeRank) {
					BWTIncSortKey(sortedRank + bwtInc->cumulativeCountInCurrentBuild[i], seq + bwtInc->cumulativeCountInCurrentBuild[i], oldInverseSa0RelativeRank - bwtInc->cumulativeCountInCurrentBuild[i]);
				}
				if (bwtInc->cumulativeCountInCurrentBuild[i+1] > oldInverseSa0RelativeRank + 1) {
					BWTIncSortKey(sortedRank + oldInverseSa0RelativeRank + 1, seq + oldInverseSa0RelativeRank + 1, bwtInc->cumulativeCountInCurrentBuild[i+1] - oldInverseSa0RelativeRank - 1);
				}
			}
		}

		// build relative rank; sortedRank is updated for merging to cater for the fact that $ is not encoded in bwt
		// the cumulative freq information is used to make sure that inverseSa0 and suffix beginning with different characters are kept in different unsorted groups)
		BWTIncBuildRelativeRank(sortedRank, seq, relativeRank, numChar, bwtInc->bwt->inverseSa0, bwtInc->cumulativeCountInCurrentBuild);
#ifdef DEBUG
		if (relativeRank[numChar] != oldInverseSa0RelativeRank) {
			fprintf(stderr, "BWTIncConstruct(): relativeRank[numChar] != oldInverseSa0RelativeRank!\n");
			exit(1);
		}
#endif

		// Sort suffix
		QSufSortSuffixSort((int*)relativeRank, (int*)seq, (int)numChar, (int)numChar, 1, TRUE);

		newInverseSa0RelativeRank = relativeRank[0];
		newInverseSa0 = sortedRank[newInverseSa0RelativeRank] + newInverseSa0RelativeRank;

		sortedRank[newInverseSa0RelativeRank] = 0;	// a special value so that this is skipped in the merged bwt

		// Build BWT
		BWTIncBuildBwt(seq, relativeRank, numChar, bwtInc->cumulativeCountInCurrentBuild);

		// Merge BWT
		mergedBwt = bwtInc->workingMemory + bwtInc->availableWord - mergedBwtSizeInWord 
				    - bwtInc->numberOfIterationDone * OCC_INTERVAL / BIT_PER_CHAR;
					// minus numberOfIteration * occInterval to create a buffer for merging
		BWTIncMergeBwt(sortedRank, bwtInc->bwt->bwtCode, insertBwt, mergedBwt, bwtInc->bwt->textLength, numChar);

	}

	// Build auxiliary structure and update info and pointers in BWT
	bwtInc->bwt->textLength += numChar;
	bwtInc->bwt->bwtCode = mergedBwt;
	bwtInc->bwt->bwtSizeInWord = mergedBwtSizeInWord;
	bwtInc->bwt->occSizeInWord = mergedOccSizeInWord;
	if (mergedBwt < bwtInc->workingMemory + mergedOccSizeInWord) {
		fprintf(stderr, "BWTIncConstruct() : Not enough memory allocated!\n");
		exit(1);
	}

	bwtInc->bwt->occValue = mergedBwt - mergedOccSizeInWord;

	BWTClearTrailingBwtCode(bwtInc->bwt);
	BWTGenerateOccValueFromBwt(bwtInc->bwt->bwtCode, bwtInc->bwt->occValue, bwtInc->bwt->occValueMajor,
							   bwtInc->bwt->textLength, bwtInc->bwt->decodeTable);

	bwtInc->bwt->inverseSa0 = newInverseSa0;
	
	bwtInc->bwt->cumulativeFreq[1] += bwtInc->cumulativeCountInCurrentBuild[1] - (bwtInc->firstCharInLastIteration <= 0);
	bwtInc->bwt->cumulativeFreq[2] += bwtInc->cumulativeCountInCurrentBuild[2] - (bwtInc->firstCharInLastIteration <= 1);
	bwtInc->bwt->cumulativeFreq[3] += bwtInc->cumulativeCountInCurrentBuild[3] - (bwtInc->firstCharInLastIteration <= 2);
	bwtInc->bwt->cumulativeFreq[4] += bwtInc->cumulativeCountInCurrentBuild[4] - (bwtInc->firstCharInLastIteration <= 3);

	bwtInc->firstCharInLastIteration = firstCharInThisIteration;

	// Set build size and text address for the next build
	BWTIncSetBuildSizeAndTextAddr(bwtInc);
	bwtInc->numberOfIterationDone++;

}

static void BWTIncSetBuildSizeAndTextAddr(BWTInc *bwtInc) {

	unsigned int maxBuildSize;

	if (bwtInc->bwt->textLength == 0) {
		// initial build
		// Minus 2 because n+1 entries of seq and rank needed for n char
		maxBuildSize = (bwtInc->availableWord - 2 - OCC_INTERVAL / CHAR_PER_WORD)
							/ (2 * CHAR_PER_WORD + 1) * CHAR_PER_WORD;
		if (bwtInc->initialMaxBuildSize > 0) {
			bwtInc->buildSize = min(bwtInc->initialMaxBuildSize, maxBuildSize);
		} else {
			bwtInc->buildSize = maxBuildSize;
		}
	} else {
		// Minus 3 because n+1 entries of sorted rank, seq and rank needed for n char
		// Minus numberOfIterationDone because bwt slightly shift to left in each iteration
		maxBuildSize = (bwtInc->availableWord - bwtInc->bwt->bwtSizeInWord - bwtInc->bwt->occSizeInWord - 3
							 - bwtInc->numberOfIterationDone * OCC_INTERVAL / BIT_PER_CHAR) 
							 / 3;
		if (maxBuildSize < CHAR_PER_WORD) {
			fprintf(stderr, "BWTIncSetBuildSizeAndTextAddr(): Not enough space allocated to continue construction!\n");
			exit(1);
		}
		if (bwtInc->incMaxBuildSize > 0) {
            bwtInc->buildSize = min(bwtInc->incMaxBuildSize, maxBuildSize);
		} else {
			bwtInc->buildSize = maxBuildSize;
		}
		if (bwtInc->buildSize < CHAR_PER_WORD) {
			bwtInc->buildSize = CHAR_PER_WORD;
		}
	}

	if (bwtInc->buildSize < CHAR_PER_WORD) {
		fprintf(stderr, "BWTIncSetBuildSizeAndTextAddr(): Not enough space allocated to continue construction!\n");
		exit(1);
	}

	bwtInc->buildSize = bwtInc->buildSize / CHAR_PER_WORD * CHAR_PER_WORD;

	bwtInc->packedText = bwtInc->workingMemory + 2 * (bwtInc->buildSize + 1);
	bwtInc->textBuffer = (unsigned char*)(bwtInc->workingMemory + bwtInc->buildSize + 1);

}

static void BWTIncPutPackedTextToRank(const unsigned int *packedText, unsigned int* __restrict rank, unsigned int* __restrict cumulativeCount, const unsigned int numChar) {

	unsigned int i, j;
	unsigned int c, t;
	unsigned int packedMask;
	unsigned int rankIndex;
	unsigned int lastWord, numCharInLastWord;

	lastWord = (numChar - 1) / CHAR_PER_WORD;
	numCharInLastWord = numChar - lastWord * CHAR_PER_WORD;

	packedMask = ALL_ONE_MASK >> (BITS_IN_WORD - BIT_PER_CHAR);
	rankIndex = numChar - 1;

	t = packedText[lastWord] >> (BITS_IN_WORD - numCharInLastWord * BIT_PER_CHAR);
	for (i=0; i<numCharInLastWord; i++) {
		c = t & packedMask;
#ifdef DEBUG
		if (c >= ALPHABET_SIZE) {
			fprintf(stderr, "BWTIncPutPackedTextToRank() : c >= ALPHABET_SIZE!\n");
			exit(1);
		}
#endif
		cumulativeCount[c+1]++;
		rank[rankIndex] = c;
		rankIndex--;
		t >>= BIT_PER_CHAR;
	}

	for (i=lastWord; i--;) {	// loop from lastWord - 1 to 0
		t = packedText[i];
		for (j=0; j<CHAR_PER_WORD; j++) {
			c = t & packedMask;
#ifdef DEBUG
			if (c >= ALPHABET_SIZE) {
				fprintf(stderr, "BWTIncPutPackedTextToRank() : c >= ALPHABET_SIZE!\n");
				exit(1);
			}
#endif
			cumulativeCount[c+1]++;
			rank[rankIndex] = c;
			rankIndex--;
			t >>= BIT_PER_CHAR;
		}
	}

	// Convert occurrence to cumulativeCount
	cumulativeCount[2] += cumulativeCount[1];
	cumulativeCount[3] += cumulativeCount[2];
	cumulativeCount[4] += cumulativeCount[3];

}

static void BWTIncBuildPackedBwt(const unsigned int *relativeRank, unsigned int* __restrict bwt, const unsigned int numChar,
								 const unsigned int *cumulativeCount, const unsigned int *packedShift) {

	unsigned int i, c, r;
	unsigned int previousRank, currentRank;
	unsigned int wordIndex, charIndex;
	unsigned int inverseSa0;

	inverseSa0 = previousRank = relativeRank[0];

	for (i=1; i<=numChar; i++) {
		currentRank = relativeRank[i];
		// previousRank > cumulativeCount[c] because $ is one of the char
		c = (previousRank > cumulativeCount[1]) + (previousRank > cumulativeCount[2]) 
											    + (previousRank > cumulativeCount[3]);
		// set bwt for currentRank
		if (c > 0) {
			// c <> 'a'
			r = currentRank;
			if (r > inverseSa0) {
				// - 1 because $ at inverseSa0 is not encoded			
				r--;
			}
			wordIndex = r / CHAR_PER_WORD;
			charIndex = r - wordIndex * CHAR_PER_WORD;
			bwt[wordIndex] |= c << packedShift[charIndex];
		}
		previousRank = currentRank;
	}

}

static unsigned int BWTIncGetAbsoluteRank(BWT *bwt, unsigned int* __restrict absoluteRank, unsigned int* __restrict seq, const unsigned int *packedText, 
								  const unsigned int numChar, const unsigned int* cumulativeCount, const unsigned int firstCharInLastIteration) {

	unsigned int saIndex;
	unsigned int lastWord;
	unsigned int packedMask;
	unsigned int i, j;
	unsigned int c, t;
	unsigned int rankIndex;
	unsigned int shift;
	unsigned int seqIndexFromStart[ALPHABET_SIZE];
	unsigned int seqIndexFromEnd[ALPHABET_SIZE];

	for (i=0; i<ALPHABET_SIZE; i++) {
		seqIndexFromStart[i] = cumulativeCount[i];
		seqIndexFromEnd[i] = cumulativeCount[i+1] - 1;
	}

	shift = BITS_IN_WORD - BIT_PER_CHAR;
	packedMask = ALL_ONE_MASK >> shift;
	saIndex = bwt->inverseSa0;
	rankIndex = numChar - 1;

	lastWord = numChar / CHAR_PER_WORD;
	for (i=lastWord; i--;) {	// loop from lastWord - 1 to 0
		t = packedText[i];
		for (j=0; j<CHAR_PER_WORD; j++) {
			c = t & packedMask;
#ifdef DEBUG
			if (c >= ALPHABET_SIZE) {
				fprintf(stderr, "BWTIncGetAbsoluteRank() : c >= ALPHABET_SIZE!\n");
				exit(1);
			}
#endif
			saIndex = bwt->cumulativeFreq[c] + BWTOccValue(bwt, saIndex, c) + 1;
#ifdef DEBUG
			if (saIndex > bwt->textLength + 1) {
				fprintf(stderr, "BWTIncGetAbsoluteRank() : saIndex > bwt->textLength + 1!\n");
				exit(1);
			}
#endif
			// A counting sort using the first character of suffix is done here
			// If rank > inverseSa0 -> fill seq from end, otherwise fill seq from start -> to leave the right entry for inverseSa0
			if (saIndex > bwt->inverseSa0) {
				seq[seqIndexFromEnd[c]] = rankIndex;
				absoluteRank[seqIndexFromEnd[c]] = saIndex;
				seqIndexFromEnd[c]--;
			} else {
				seq[seqIndexFromStart[c]] = rankIndex;
				absoluteRank[seqIndexFromStart[c]] = saIndex;
				seqIndexFromStart[c]++;
			}
			rankIndex--;
			t >>= BIT_PER_CHAR;
		}
	}

	absoluteRank[seqIndexFromStart[firstCharInLastIteration]] = bwt->inverseSa0;	// representing the substring of all preceding characters
	seq[seqIndexFromStart[firstCharInLastIteration]] = numChar;

#ifdef DEBUG
	if (seqIndexFromStart[firstCharInLastIteration] != seqIndexFromEnd[firstCharInLastIteration]) {
		fprintf(stderr, "BWTIncGetAbsoluteRank(): seqIndexFromStart[firstCharInLastIteration] != seqIndexFromEnd[firstCharInLastIteration]!\n");
	}
#endif

	return seqIndexFromStart[firstCharInLastIteration];

}

static void BWTIncSortKey(unsigned int* __restrict key, unsigned int* __restrict seq, const unsigned int numItem) {

	#define EQUAL_KEY_THRESHOLD	4	// Partition for equal key if data array size / the number of data with equal value with pivot < EQUAL_KEY_THRESHOLD

	int lowIndex, highIndex, midIndex;
	int lowPartitionIndex, highPartitionIndex;
	int lowStack[32], highStack[32];
	int stackDepth;
	int i, j;
	unsigned int tempSeq, tempKey;
	int numberOfEqualKey;

	if (numItem < 2) {
		return;
	}

	stackDepth = 0;

    lowIndex = 0;
    highIndex = numItem - 1;

	for (;;) {

		for (;;) {

			// Sort small array of data
			if (highIndex - lowIndex < BWTINC_INSERT_SORT_NUM_ITEM) {	 // Insertion sort on smallest arrays
				for (i=lowIndex+1; i<=highIndex; i++) {
					tempSeq = seq[i];
					tempKey = key[i];
					for (j = i; j > lowIndex && key[j-1] > tempKey; j--) {
						seq[j] = seq[j-1];
						key[j] = key[j-1];
					}
					if (j != i) {
						seq[j] = tempSeq;
						key[j] = tempKey;
					}
				}
				break;
			}

			// Choose pivot as median of the lowest, middle, and highest data; sort the three data

			midIndex = average(lowIndex, highIndex);
			if (key[lowIndex] > key[midIndex]) {
				tempSeq = seq[lowIndex];
				tempKey = key[lowIndex];
				seq[lowIndex] = seq[midIndex];
				key[lowIndex] = key[midIndex];
				seq[midIndex] = tempSeq;
				key[midIndex] = tempKey;
			}
			if (key[lowIndex] > key[highIndex]) {
				tempSeq = seq[lowIndex];
				tempKey = key[lowIndex];
				seq[lowIndex] = seq[highIndex];
				key[lowIndex] = key[highIndex];
				seq[highIndex] = tempSeq;
				key[highIndex] = tempKey;
			}
			if (key[midIndex] > key[highIndex]) {
				tempSeq = seq[midIndex];
				tempKey = key[midIndex];
				seq[midIndex] = seq[highIndex];
				key[midIndex] = key[highIndex];
				seq[highIndex] = tempSeq;
				key[highIndex] = tempKey;
			}

			// Partition data

			numberOfEqualKey = 0;

			lowPartitionIndex = lowIndex + 1;
			highPartitionIndex = highIndex - 1;

			for (;;) {
				while (lowPartitionIndex <= highPartitionIndex && key[lowPartitionIndex] <= key[midIndex]) {
					numberOfEqualKey += (key[lowPartitionIndex] == key[midIndex]);
					lowPartitionIndex++;
				}
				while (lowPartitionIndex < highPartitionIndex) {
					if (key[midIndex] >= key[highPartitionIndex]) {
						numberOfEqualKey += (key[midIndex] == key[highPartitionIndex]);
						break;
					}
					highPartitionIndex--;
				}
				if (lowPartitionIndex >= highPartitionIndex) {
					break;
				}
				tempSeq = seq[lowPartitionIndex];
				tempKey = key[lowPartitionIndex];
				seq[lowPartitionIndex] = seq[highPartitionIndex];
				key[lowPartitionIndex] = key[highPartitionIndex];
				seq[highPartitionIndex] = tempSeq;
				key[highPartitionIndex] = tempKey;
				if (highPartitionIndex == midIndex) {
					// partition key has been moved
					midIndex = lowPartitionIndex;
				}
				lowPartitionIndex++;
				highPartitionIndex--;
			}

			// Adjust the partition index
			highPartitionIndex = lowPartitionIndex;
			lowPartitionIndex--;

			// move the partition key to end of low partition
			tempSeq = seq[midIndex];
			tempKey = key[midIndex];
			seq[midIndex] = seq[lowPartitionIndex];
			key[midIndex] = key[lowPartitionIndex];
			seq[lowPartitionIndex] = tempSeq;
			key[lowPartitionIndex] = tempKey;

			if (highIndex - lowIndex + BWTINC_INSERT_SORT_NUM_ITEM <= EQUAL_KEY_THRESHOLD * numberOfEqualKey) {

				// Many keys = partition key; separate the equal key data from the lower partition
		
				midIndex = lowIndex;

				for (;;) {
					while (midIndex < lowPartitionIndex && key[midIndex] < key[lowPartitionIndex]) {
						midIndex++;
					}
					while (midIndex < lowPartitionIndex && key[lowPartitionIndex] == key[lowPartitionIndex - 1]) {
						lowPartitionIndex--;
					}
					if (midIndex >= lowPartitionIndex) {
						break;
					}
					tempSeq = seq[midIndex];
					tempKey = key[midIndex];
					seq[midIndex] = seq[lowPartitionIndex - 1];
					key[midIndex] = key[lowPartitionIndex - 1];
					seq[lowPartitionIndex - 1] = tempSeq;
					key[lowPartitionIndex - 1] = tempKey;
					midIndex++;
					lowPartitionIndex--;
				}

			}

			if (lowPartitionIndex - lowIndex > highIndex - highPartitionIndex) {
				// put the larger partition to stack
				lowStack[stackDepth] = lowIndex;
				highStack[stackDepth] = lowPartitionIndex - 1;
				stackDepth++;
				// sort the smaller partition first
				lowIndex = highPartitionIndex;
			} else {
				// put the larger partition to stack
				lowStack[stackDepth] = highPartitionIndex;
				highStack[stackDepth] = highIndex;
				stackDepth++;
				// sort the smaller partition first
				if (lowPartitionIndex > lowIndex) {
					highIndex = lowPartitionIndex - 1;
				} else {
					// all keys in the partition equals to the partition key
					break;
				}
			}
			continue;

		}

		// Pop a range from stack
		if (stackDepth > 0) {
			stackDepth--;
			lowIndex = lowStack[stackDepth];
			highIndex = highStack[stackDepth];
			continue;
		} else {
			return;
		}

	}


}

static void BWTIncBuildRelativeRank(unsigned int* __restrict sortedRank, unsigned int* __restrict seq, unsigned int* __restrict relativeRank, const unsigned int numItem, 
									unsigned int oldInverseSa0, const unsigned int *cumulativeCount) {

	unsigned int i, c;
	unsigned int s, r;
	unsigned int lastRank, lastIndex;
	unsigned int oldInverseSa0RelativeRank = 0;
	unsigned int freq;

	lastIndex = numItem;
	lastRank = sortedRank[numItem];
	if (lastRank > oldInverseSa0) {
		sortedRank[numItem]--;	// to prepare for merging; $ is not encoded in bwt
	}
	s = seq[numItem];
	relativeRank[s] = numItem;
	if (lastRank == oldInverseSa0) {
		oldInverseSa0RelativeRank = numItem;
		oldInverseSa0++;	// so that this segment of code is not run again
		lastRank++;			// so that oldInverseSa0 become a sorted group with 1 item
	}

	c = ALPHABET_SIZE - 1;
	freq = cumulativeCount[c];

	for (i=numItem; i--;) {	// from numItem - 1 to 0
		r = sortedRank[i];
		if (r > oldInverseSa0) {
			sortedRank[i]--;	// to prepare for merging; $ is not encoded in bwt
		}
		s = seq[i];
		if (i < freq) {
			if (lastIndex >= freq) {
				lastRank++;	// to trigger the group across alphabet boundary to be split
			}
			c--;
			freq = cumulativeCount[c];
		}
		if (r == lastRank) {
			relativeRank[s] = lastIndex;
		} else {
			if (i == lastIndex - 1) {
				if (lastIndex < numItem && (int)seq[lastIndex + 1] < 0) {
					seq[lastIndex] = seq[lastIndex + 1] - 1;
				} else {
					seq[lastIndex] = (unsigned int)-1;
				}
			}
			lastIndex = i;
			lastRank = r;
			relativeRank[s] = i;
			if (r == oldInverseSa0) {
				oldInverseSa0RelativeRank = i;
				oldInverseSa0++;	// so that this segment of code is not run again
				lastRank++;			// so that oldInverseSa0 become a sorted group with 1 item
			}
		}
	}

}

static void BWTIncBuildBwt(unsigned int*  seq, const unsigned int *relativeRank, const unsigned int numChar, const unsigned int *cumulativeCount) {

	unsigned int i, c;
	unsigned int previousRank, currentRank;

	previousRank = relativeRank[0];

	for (i=1; i<=numChar; i++) {
		currentRank = relativeRank[i];
		c = (previousRank >= cumulativeCount[1]) + (previousRank >= cumulativeCount[2])
											  	 + (previousRank >= cumulativeCount[3]);
		seq[currentRank] = c;
		previousRank = currentRank;
	}

}

static void BWTIncMergeBwt(const unsigned int *sortedRank, const unsigned int* oldBwt, const unsigned int *insertBwt, unsigned int* __restrict mergedBwt, 
						   const unsigned int numOldBwt, const unsigned int numInsertBwt) {

	unsigned int bitsInWordMinusBitPerChar;
	unsigned int leftShift, rightShift;
	unsigned int o;
	unsigned int oIndex, iIndex, mIndex;
	unsigned int mWord, mChar, oWord, oChar;
	unsigned int numInsert;

	bitsInWordMinusBitPerChar = BITS_IN_WORD - BIT_PER_CHAR;

	oIndex = 0;
	iIndex = 0;
	mIndex = 0;

	mWord = 0;
	mChar = 0;

	mergedBwt[0] = 0;	// this can be cleared as merged Bwt slightly shift to the left in each iteration

	while (oIndex < numOldBwt) {

		// copy from insertBwt
		while (iIndex <= numInsertBwt && sortedRank[iIndex] <= oIndex) {
			if (sortedRank[iIndex] != 0) {	// special value to indicate that this is for new inverseSa0
				mergedBwt[mWord] |= insertBwt[iIndex] << (BITS_IN_WORD - (mChar + 1) * BIT_PER_CHAR);
				mIndex++;
				mChar++;
				if (mChar == CHAR_PER_WORD) {
					mChar = 0;
					mWord++;
					mergedBwt[mWord] = 0;	// no need to worry about crossing mergedBwt boundary
				}
			}
			iIndex++;
		}

		// Copy from oldBwt to mergedBwt
		if (iIndex <= numInsertBwt) {
			o = sortedRank[iIndex];
		} else {
			o = numOldBwt;
		}
		numInsert = o - oIndex;

		oWord = oIndex / CHAR_PER_WORD;
		oChar = oIndex - oWord * CHAR_PER_WORD;
		if (oChar > mChar) {
			leftShift = (oChar - mChar) * BIT_PER_CHAR;
			rightShift = (CHAR_PER_WORD + mChar - oChar) * BIT_PER_CHAR;
			mergedBwt[mWord] = mergedBwt[mWord]
								| (oldBwt[oWord] << (oChar * BIT_PER_CHAR) >> (mChar * BIT_PER_CHAR))
								| (oldBwt[oWord+1] >> rightShift);
			oIndex += minX(numInsert, CHAR_PER_WORD - mChar);
			while (o > oIndex) {
				oWord++;
				mWord++;
				mergedBwt[mWord] = (oldBwt[oWord] << leftShift) | (oldBwt[oWord+1] >> rightShift);
				oIndex += CHAR_PER_WORD;
			}
		} else if (oChar < mChar) {
			rightShift = (mChar - oChar) * BIT_PER_CHAR;
			leftShift = (CHAR_PER_WORD + oChar - mChar) * BIT_PER_CHAR;
			mergedBwt[mWord] = mergedBwt[mWord] 
								| (oldBwt[oWord] << (oChar * BIT_PER_CHAR) >> (mChar * BIT_PER_CHAR));
			oIndex += minX(numInsert, CHAR_PER_WORD - mChar);
			while (o > oIndex) {
				oWord++;
				mWord++;
				mergedBwt[mWord] = (oldBwt[oWord-1] << leftShift) | (oldBwt[oWord] >> rightShift);
				oIndex += CHAR_PER_WORD;
			}
		} else { // oChar == mChar
			mergedBwt[mWord] = mergedBwt[mWord] | truncateLeft(oldBwt[oWord], mChar * BIT_PER_CHAR);
			oIndex += minX(numInsert, CHAR_PER_WORD - mChar);
			while (o > oIndex) {
				oWord++;
				mWord++;
				mergedBwt[mWord] = oldBwt[oWord];
				oIndex += CHAR_PER_WORD;
			}
		}
		oIndex = o;
		mIndex += numInsert;

		// Clear the trailing garbage in mergedBwt
		mWord = mIndex / CHAR_PER_WORD;
		mChar = mIndex - mWord * CHAR_PER_WORD;
		if (mChar == 0) {
			mergedBwt[mWord] = 0;
		} else {
			mergedBwt[mWord] = truncateRight(mergedBwt[mWord], (BITS_IN_WORD - mChar * BIT_PER_CHAR));
		}

	}

	// copy from insertBwt
	while (iIndex <= numInsertBwt) {
		if (sortedRank[iIndex] != 0) {
			mergedBwt[mWord] |= insertBwt[iIndex] << (BITS_IN_WORD - (mChar + 1) * BIT_PER_CHAR);
			mIndex++;
			mChar++;
			if (mChar == CHAR_PER_WORD) {
				mChar = 0;
				mWord++;
				mergedBwt[mWord] = 0;	// no need to worry about crossing mergedBwt boundary
			}
		}
		iIndex++;
	}

}

unsigned int BWTGenerateOccValueToFileFromBwt(const char *bwtFileName, const char *occValueFileName, unsigned int*  decodeTable) {

	FILE *bwtFile, *occValueFile;
	unsigned int *bwt;
	unsigned int bwtFileSizeInWord, bwtResidentSizeInWord;
	unsigned int textLength;
	unsigned int inverseSa0;
	unsigned int *occValue, *occValueMajor;
	unsigned int occSizeInWord, occMajorSizeInWord;
	unsigned int i;
	unsigned int cumulativeFreq[ALPHABET_SIZE];

	bwtFile = (FILE*)fopen64(bwtFileName, "rb");
	if (bwtFile == NULL) {
		fprintf(stderr, "BWTGenerateOccValueToFileFromBwt(): Cannot open BWT file!\n");
		exit(1);
	}

	fread(&inverseSa0, sizeof(unsigned int), 1, bwtFile);

	fread(cumulativeFreq, sizeof(unsigned int), ALPHABET_SIZE, bwtFile);
	textLength = cumulativeFreq[ALPHABET_SIZE - 1];

	bwtResidentSizeInWord = BWTResidentSizeInWord(textLength);
	bwt = MMUnitAllocate(bwtResidentSizeInWord * sizeof(unsigned int));
	bwtFileSizeInWord = BWTFileSizeInWord(textLength);
	fread(bwt, sizeof(unsigned int), bwtFileSizeInWord, bwtFile);
	fclose(bwtFile);
	for (i=bwtFileSizeInWord; i<bwtResidentSizeInWord; i++) {
		bwt[i] = 0;
	}

	// occValue File

	occValueFile = (FILE*)fopen64(occValueFileName, "wb");
	if (occValueFile == NULL) {
		fprintf(stderr, "BWTGenerateOccValueToFileFromBwt(): Cannot open occ value file!\n");
		exit(1);
	}

	fwrite(&inverseSa0, sizeof(unsigned int), 1, occValueFile);
	fwrite(cumulativeFreq, sizeof(unsigned int), ALPHABET_SIZE, occValueFile);

	occSizeInWord = BWTOccValueMinorSizeInWord(textLength);
	occMajorSizeInWord = BWTOccValueMajorSizeInWord(textLength);
	occValue = MMUnitAllocate(occSizeInWord * sizeof(unsigned int));
	occValueMajor = MMUnitAllocate(occMajorSizeInWord * sizeof(unsigned int));

	if (decodeTable == NULL) {
		decodeTable = MMUnitAllocate(DNA_OCC_CNT_TABLE_SIZE_IN_WORD * sizeof(unsigned int));
		GenerateDNAOccCountTable(decodeTable);
	}

	BWTGenerateOccValueFromBwt(bwt, occValue, occValueMajor, textLength, decodeTable);

	fwrite(occValue, sizeof(unsigned int), occSizeInWord, occValueFile);
	fwrite(occValueMajor, sizeof(unsigned int), occMajorSizeInWord, occValueFile);
	fclose(occValueFile);

	MMUnitFree(occValue, occSizeInWord * sizeof(unsigned int));
	MMUnitFree(occValueMajor, occMajorSizeInWord * sizeof(unsigned int));
	MMUnitFree(bwt, bwtResidentSizeInWord * sizeof(unsigned int));
	MMUnitFree(decodeTable, DNA_OCC_CNT_TABLE_SIZE_IN_WORD * sizeof(unsigned int));

	return textLength;

}

void BWTGenerateOccValueFromBwt(const unsigned int*  bwt, unsigned int* __restrict occValue, unsigned int* __restrict occValueMajor,
								const unsigned int textLength, const unsigned int*  decodeTable) {

	unsigned int numberOfOccValueMajor, numberOfOccValue;
	unsigned int wordBetweenOccValue;
	unsigned int numberOfOccIntervalPerMajor;
	unsigned int c;
	unsigned int i, j;
	unsigned int occMajorIndex;
	unsigned int occIndex, bwtIndex;
	unsigned int sum;
	unsigned int tempOccValue0[ALPHABET_SIZE], tempOccValue1[ALPHABET_SIZE];

	wordBetweenOccValue = OCC_INTERVAL / CHAR_PER_WORD;

	// Calculate occValue

	numberOfOccValue = (textLength + OCC_INTERVAL - 1) / OCC_INTERVAL + 1;				// Value at both end for bi-directional encoding
	numberOfOccIntervalPerMajor = OCC_INTERVAL_MAJOR / OCC_INTERVAL;
	numberOfOccValueMajor = (numberOfOccValue + numberOfOccIntervalPerMajor - 1) / numberOfOccIntervalPerMajor;

	tempOccValue0[0] = 0;
	tempOccValue0[1] = 0;
	tempOccValue0[2] = 0;
	tempOccValue0[3] = 0;
	occValueMajor[0] = 0;
	occValueMajor[1] = 0;
	occValueMajor[2] = 0;
	occValueMajor[3] = 0;

	occIndex = 0;
	bwtIndex = 0;
	for (occMajorIndex=1; occMajorIndex<numberOfOccValueMajor; occMajorIndex++) {

		for (i=0; i<numberOfOccIntervalPerMajor/2; i++) {

			sum = 0;
			tempOccValue1[0] = tempOccValue0[0];
			tempOccValue1[1] = tempOccValue0[1];
			tempOccValue1[2] = tempOccValue0[2];
			tempOccValue1[3] = tempOccValue0[3];

			for (j=0; j<wordBetweenOccValue; j++) {
				c = bwt[bwtIndex];
				sum += decodeTable[c >> 16];
				sum += decodeTable[c & 0x0000FFFF];
				bwtIndex++;
			}
			if (!DNA_OCC_SUM_EXCEPTION(sum)) {
				tempOccValue1[0] += (sum & 0x000000FF);	sum >>= 8;
				tempOccValue1[1] += (sum & 0x000000FF);	sum >>= 8;
				tempOccValue1[2] += (sum & 0x000000FF);	sum >>= 8;
				tempOccValue1[3] += sum;
			} else {
				if (sum == 0x00000100) {
					tempOccValue1[0] += 256;
				} else if (sum == 0x00010000) {
					tempOccValue1[1] += 256;
				} else if (sum == 0x01000000) {
					tempOccValue1[2] += 256;
				} else {
					tempOccValue1[3] += 256;
				}
			}
			occValue[occIndex * 4 + 0] = (tempOccValue0[0] << 16) | tempOccValue1[0];
			occValue[occIndex * 4 + 1] = (tempOccValue0[1] << 16) | tempOccValue1[1];
			occValue[occIndex * 4 + 2] = (tempOccValue0[2] << 16) | tempOccValue1[2];
			occValue[occIndex * 4 + 3] = (tempOccValue0[3] << 16) | tempOccValue1[3];
			tempOccValue0[0] = tempOccValue1[0];
			tempOccValue0[1] = tempOccValue1[1];
			tempOccValue0[2] = tempOccValue1[2];
			tempOccValue0[3] = tempOccValue1[3];
			sum = 0;

			occIndex++;

			for (j=0; j<wordBetweenOccValue; j++) {
				c = bwt[bwtIndex];
				sum += decodeTable[c >> 16];
				sum += decodeTable[c & 0x0000FFFF];
				bwtIndex++;
			}
			if (!DNA_OCC_SUM_EXCEPTION(sum)) {
				tempOccValue0[0] += (sum & 0x000000FF);	sum >>= 8;
				tempOccValue0[1] += (sum & 0x000000FF);	sum >>= 8;
				tempOccValue0[2] += (sum & 0x000000FF);	sum >>= 8;
				tempOccValue0[3] += sum;
			} else {
				if (sum == 0x00000100) {
					tempOccValue0[0] += 256;
				} else if (sum == 0x00010000) {
					tempOccValue0[1] += 256;
				} else if (sum == 0x01000000) {
					tempOccValue0[2] += 256;
				} else {
					tempOccValue0[3] += 256;
				}
			}
		}

		occValueMajor[occMajorIndex * 4 + 0] = occValueMajor[(occMajorIndex - 1) * 4 + 0] + tempOccValue0[0];
		occValueMajor[occMajorIndex * 4 + 1] = occValueMajor[(occMajorIndex - 1) * 4 + 1] + tempOccValue0[1];
		occValueMajor[occMajorIndex * 4 + 2] = occValueMajor[(occMajorIndex - 1) * 4 + 2] + tempOccValue0[2];
		occValueMajor[occMajorIndex * 4 + 3] = occValueMajor[(occMajorIndex - 1) * 4 + 3] + tempOccValue0[3];
		tempOccValue0[0] = 0;
		tempOccValue0[1] = 0;
		tempOccValue0[2] = 0;
		tempOccValue0[3] = 0;

	}

	while (occIndex < (numberOfOccValue-1)/2) {
		sum = 0;
		tempOccValue1[0] = tempOccValue0[0];
		tempOccValue1[1] = tempOccValue0[1];
		tempOccValue1[2] = tempOccValue0[2];
		tempOccValue1[3] = tempOccValue0[3];
		for (j=0; j<wordBetweenOccValue; j++) {
			c = bwt[bwtIndex];
			sum += decodeTable[c >> 16];
			sum += decodeTable[c & 0x0000FFFF];
			bwtIndex++;
		}
		if (!DNA_OCC_SUM_EXCEPTION(sum)) {
			tempOccValue1[0] += (sum & 0x000000FF);	sum >>= 8;
			tempOccValue1[1] += (sum & 0x000000FF);	sum >>= 8;
			tempOccValue1[2] += (sum & 0x000000FF);	sum >>= 8;
			tempOccValue1[3] += sum;
		} else {
			if (sum == 0x00000100) {
				tempOccValue1[0] += 256;
			} else if (sum == 0x00010000) {
				tempOccValue1[1] += 256;
			} else if (sum == 0x01000000) {
				tempOccValue1[2] += 256;
			} else {
				tempOccValue1[3] += 256;
			}
		}
		occValue[occIndex * 4 + 0] = (tempOccValue0[0] << 16) | tempOccValue1[0];
		occValue[occIndex * 4 + 1] = (tempOccValue0[1] << 16) | tempOccValue1[1];
		occValue[occIndex * 4 + 2] = (tempOccValue0[2] << 16) | tempOccValue1[2];
		occValue[occIndex * 4 + 3] = (tempOccValue0[3] << 16) | tempOccValue1[3];
		tempOccValue0[0] = tempOccValue1[0];
		tempOccValue0[1] = tempOccValue1[1];
		tempOccValue0[2] = tempOccValue1[2];
		tempOccValue0[3] = tempOccValue1[3];
		sum = 0;
		occIndex++;

		for (j=0; j<wordBetweenOccValue; j++) {
			c = bwt[bwtIndex];
			sum += decodeTable[c >> 16];
			sum += decodeTable[c & 0x0000FFFF];
			bwtIndex++;
		}
		if (!DNA_OCC_SUM_EXCEPTION(sum)) {
			tempOccValue0[0] += (sum & 0x000000FF);	sum >>= 8;
			tempOccValue0[1] += (sum & 0x000000FF);	sum >>= 8;
			tempOccValue0[2] += (sum & 0x000000FF);	sum >>= 8;
			tempOccValue0[3] += sum;
		} else {
			if (sum == 0x00000100) {
				tempOccValue0[0] += 256;
			} else if (sum == 0x00010000) {
				tempOccValue0[1] += 256;
			} else if (sum == 0x01000000) {
				tempOccValue0[2] += 256;
			} else {
				tempOccValue0[3] += 256;
			}
		}
	}

	sum = 0;
	tempOccValue1[0] = tempOccValue0[0];
	tempOccValue1[1] = tempOccValue0[1];
	tempOccValue1[2] = tempOccValue0[2];
	tempOccValue1[3] = tempOccValue0[3];

	if (occIndex * 2 < numberOfOccValue - 1) {
		for (j=0; j<wordBetweenOccValue; j++) {
			c = bwt[bwtIndex];
			sum += decodeTable[c >> 16];
			sum += decodeTable[c & 0x0000FFFF];
			bwtIndex++;
		}
		if (!DNA_OCC_SUM_EXCEPTION(sum)) {
			tempOccValue1[0] += (sum & 0x000000FF);	sum >>= 8;
			tempOccValue1[1] += (sum & 0x000000FF);	sum >>= 8;
			tempOccValue1[2] += (sum & 0x000000FF);	sum >>= 8;
			tempOccValue1[3] += sum;
		} else {
			if (sum == 0x00000100) {
				tempOccValue1[0] += 256;
			} else if (sum == 0x00010000) {
				tempOccValue1[1] += 256;
			} else if (sum == 0x01000000) {
				tempOccValue1[2] += 256;
			} else {
				tempOccValue1[3] += 256;
			}
		}
	}

	occValue[occIndex * 4 + 0] = (tempOccValue0[0] << 16) | tempOccValue1[0];
	occValue[occIndex * 4 + 1] = (tempOccValue0[1] << 16) | tempOccValue1[1];
	occValue[occIndex * 4 + 2] = (tempOccValue0[2] << 16) | tempOccValue1[2];
	occValue[occIndex * 4 + 3] = (tempOccValue0[3] << 16) | tempOccValue1[3];

}

void BWTSaveBwtCodeAndOcc(const BWT *bwt, const char *bwtFileName, const char *occValueFileName) {

	FILE *bwtFile, *occValueFile;
	unsigned int bwtLength;

	bwtFile = (FILE*)fopen64(bwtFileName, "wb");
	if (bwtFile == NULL) {
		fprintf(stderr, "BWTSaveBwtCodeAndOcc(): Cannot open BWT code file!\n");
		exit(1);
	}

	fwrite(&bwt->inverseSa0, sizeof(unsigned int), 1, bwtFile);
	fwrite(bwt->cumulativeFreq + 1, sizeof(unsigned int), ALPHABET_SIZE, bwtFile);
	bwtLength = BWTFileSizeInWord(bwt->textLength);
	fwrite(bwt->bwtCode, sizeof(unsigned int), bwtLength, bwtFile);
	fclose(bwtFile);

	occValueFile = (FILE*)fopen64(occValueFileName, "wb");
	if (occValueFile == NULL) {
		fprintf(stderr, "BWTSaveBwtCodeAndOcc(): Cannot open occ value file!\n");
		exit(1);
	}

	fwrite(&bwt->inverseSa0, sizeof(unsigned int), 1, occValueFile);
	fwrite(bwt->cumulativeFreq + 1, sizeof(unsigned int), ALPHABET_SIZE, occValueFile);

	fwrite(bwt->occValue, sizeof(unsigned int), bwt->occSizeInWord, occValueFile);
	fwrite(bwt->occValueMajor, sizeof(unsigned int), bwt->occMajorSizeInWord, occValueFile);
	fclose(occValueFile);

}

void BWTGenerateSaValue(MMPool *mmPool, BWT *bwt, const unsigned int saValueFreq, unsigned int showProgressInterval) {

	unsigned int saValue;
	unsigned int saIndex;
	unsigned int numberOfSaValue;
	unsigned int numberOfSaValueGenerated;
	unsigned int progressInterval;
	unsigned int i;

	if (bwt->saInterval != ALL_ONE_MASK) {
		fprintf(stderr, "BWTGenerateSaValue() : saValue already exist!\n");
		exit(1);
	}

	saValue = bwt->textLength;
	saIndex = 0;

	numberOfSaValue = (bwt->textLength + saValueFreq) / saValueFreq;
	if (numberOfSaValue > ALL_ONE_MASK / sizeof(unsigned int)) {
		fprintf(stderr, "BWTGenerateSaValue() : numberOfSaValue > limit!\n");
		exit(1);
	}

	bwt->saValueSize = sizeof(unsigned int) * numberOfSaValue;
	bwt->saInterval = saValueFreq;

	bwt->saValue = MMUnitAllocate(bwt->saValueSize);

	#ifdef DEBUG
	initializeVAL(bwt->saValue, numberOfSaValue, ALL_ONE_MASK);
	#endif

	numberOfSaValueGenerated = 0;

	if (showProgressInterval == 0 || showProgressInterval > numberOfSaValue) {
		progressInterval = numberOfSaValue;
	} else {
		progressInterval = showProgressInterval;
	}

	while (numberOfSaValueGenerated < numberOfSaValue) {
		progressInterval = minX(numberOfSaValue - numberOfSaValueGenerated, progressInterval);
		for (i=0; i<progressInterval; i++) {
			while (saIndex % saValueFreq != 0) {
				#ifdef DEBUG
				if (saValue == 0) {
					fprintf(stderr, "BWTGenerateSaValue(): saValue < 0!\n");
					exit(1);
				}
				#endif
				saValue--;
				saIndex = BWTPsiMinusValue(bwt, saIndex);
				#ifdef DEBUG
				if (saIndex > bwt->textLength) {
					fprintf(stderr, "BWTGenerateSaValue() : saIndex > textLength!\n");
					exit(1);
				}
				#endif
			}
			bwt->saValue[saIndex/saValueFreq] = saValue;
			saValue--;
			saIndex = BWTPsiMinusValue(bwt, saIndex);
		}
		numberOfSaValueGenerated += progressInterval;
		if (showProgressInterval > 0) {
			printf("SA Value generated : %u\n", numberOfSaValueGenerated);
		}
	}

	#ifdef DEBUG
	if (numberOfMatchInVAL(bwt->saValue, numberOfSaValue, ALL_ONE_MASK) > 0) {
		fprintf(stderr, "BWTGenerateSaValue() : some saValue is not filled!\n");
		exit(1);
	}
	#endif

	bwt->saValue[0] = (unsigned int)-1;	// Special handling
	BWTGenerateSaValueOnBoundary(mmPool, bwt);

}

void BWTSaveSaValue(const BWT *bwt, const char *saValueFileName) {

	FILE *saValueFile;

	saValueFile = (FILE*)fopen64(saValueFileName, "wb");
	if (saValueFile == NULL) {
		fprintf(stderr, "BWTSaveSaValue() : cannot open saValueFile!\n");
		exit(1);
	}

	fwrite(&bwt->inverseSa0, sizeof(unsigned int), 1, saValueFile);
	fwrite(bwt->cumulativeFreq + 1, sizeof(unsigned int), ALPHABET_SIZE, saValueFile);
	fwrite(&bwt->saInterval, sizeof(unsigned int), 1, saValueFile);

	fwrite(&bwt->textLength, sizeof(unsigned int), 1, saValueFile);	// Save SA values without special handling on SA[0]
	fwrite(bwt->saValue + 1, 1, bwt->saValueSize - sizeof(unsigned int), saValueFile);

	fclose(saValueFile);

}

void BWTGenerateInverseSa(BWT *bwt, const unsigned int inverseSaFreq, unsigned int showProgressInterval) {

	unsigned int saValue;
	unsigned int saIndex;
	unsigned int numberOfInverseSa;
	unsigned int numberOfInverseSaGenerated;
	unsigned int progressInterval;
	unsigned int i;
	
	if (bwt->inverseSaInterval != ALL_ONE_MASK) {
		fprintf(stderr, "BWTGenerateInverseSa() : inverseSa already exist!\n");
		exit(1);
	}

	saValue = bwt->textLength;
	saIndex = 0;

	numberOfInverseSa = (bwt->textLength + inverseSaFreq) / inverseSaFreq;
	bwt->inverseSaSize = sizeof(unsigned int) * numberOfInverseSa;
	bwt->inverseSaInterval = inverseSaFreq;

	bwt->inverseSa = MMUnitAllocate(bwt->inverseSaSize);

	#ifdef DEBUG
	initializeVAL(bwt->inverseSa, numberOfInverseSa, ALL_ONE_MASK);
	#endif

	numberOfInverseSaGenerated = 0;

	if (showProgressInterval == 0 || showProgressInterval > numberOfInverseSa) {
		progressInterval = numberOfInverseSa;
	} else {
		progressInterval = showProgressInterval;
	}

	while (numberOfInverseSaGenerated < numberOfInverseSa) {
		progressInterval = minX(numberOfInverseSa - numberOfInverseSaGenerated, progressInterval);
		for (i=0; i<progressInterval; i++) {
			while (saValue % inverseSaFreq != 0) {
				saValue--;
				saIndex = BWTPsiMinusValue(bwt, saIndex);
				#ifdef DEBUG
				if (saIndex > bwt->textLength) {
					fprintf(stderr, "BWTGenerateInverseSa() : saIndex > textLength!\n");
					exit(1);
				}
				#endif
			}
			bwt->inverseSa[saValue/inverseSaFreq] = saIndex;
			saValue--;
			saIndex = BWTPsiMinusValue(bwt, saIndex);
		}
		numberOfInverseSaGenerated += progressInterval;
		if (showProgressInterval > 0) {
			printf("Inverse SA generated : %u\n", numberOfInverseSaGenerated);
		}
	}

	#ifdef DEBUG
	if (numberOfMatchInVAL(bwt->inverseSa, numberOfInverseSa, ALL_ONE_MASK) > 0) {
		fprintf(stderr, "BWTGenerateInverseSa() : some inverseSa is not filled!\n");
		exit(1);
	}
	#endif

}

void BWTSaveInverseSa(const BWT *bwt, const char *inverseSaFileName) {

	FILE *inverseSaFile;

	inverseSaFile = (FILE*)fopen64(inverseSaFileName, "wb");
	if (inverseSaFile == NULL) {
		fprintf(stderr, "BWTSaveInverseSa() : cannot open inverseSaFile!\n");
		exit(1);
	}

	fwrite(&bwt->inverseSa0, sizeof(unsigned int), 1, inverseSaFile);
	fwrite(bwt->cumulativeFreq + 1, sizeof(unsigned int), ALPHABET_SIZE, inverseSaFile);
	fwrite(&bwt->inverseSaInterval, sizeof(unsigned int), 1, inverseSaFile);
	fwrite(bwt->inverseSa, 1, bwt->inverseSaSize, inverseSaFile);

	fclose(inverseSaFile);

}

void BWTGenerateSaRangeTable(const BWT *bwt, const unsigned int numOfChar, const char *saRangeFileName) {

	// stack
	unsigned int saIndexLeft[MAX_ARPROX_MATCH_LENGTH + 1];
	unsigned int saIndexRight[MAX_ARPROX_MATCH_LENGTH];
	unsigned char generatedPattern[MAX_ARPROX_MATCH_LENGTH];
	unsigned int pos;
	unsigned int i, c;
	unsigned int numOfPattern;
	unsigned int saRangeIndex;
	unsigned int mask[16] = { 0xFFFFFFFC, 0xFFFFFFF3, 0xFFFFFFCF, 0xFFFFFF3F,
					 0xFFFFFCFF, 0xFFFFF3FF, 0xFFFFCFFF, 0xFFFF3FFF,
					 0xFFFCFFFF, 0xFFF3FFFF, 0xFFCFFFFF, 0xFF3FFFFF,
					 0xFCFFFFFF, 0xF3FFFFFF, 0xCFFFFFFF, 0x3FFFFFFF };

	SaIndexRange *saIndexRange;

	FILE * saRangeFile;

	saRangeFile = (FILE*)fopen64(saRangeFileName, "wb");
	if (saRangeFile == NULL) {
		fprintf(stderr, "BWTGenerateSaRangeTable(): Cannot open SA range file!\n");
		exit(1);
	}

	fwrite(&bwt->inverseSa0, sizeof(unsigned int), 1, saRangeFile);
	fwrite(bwt->cumulativeFreq + 1, sizeof(unsigned int), ALPHABET_SIZE, saRangeFile);
	fwrite(&numOfChar, sizeof(unsigned int), 1, saRangeFile);

	numOfPattern = 1 << (numOfChar * 2);	// 4^numOfChar
	saIndexRange = MMUnitAllocate(numOfPattern * sizeof(SaIndexRange));
	for (i=0; i<numOfPattern; i++) {
		saIndexRange[i].startSaIndex = 1;
		saIndexRange[i].endSaIndex = 0;
	}

	// Set this boundary case so that the last character of generated pattern can be located by backward search
	saIndexLeft[numOfChar] = 0;
	saIndexRight[numOfChar] = bwt->textLength;

	// Set initial state
	pos = numOfChar - 1;
	generatedPattern[pos] = (unsigned char)-1;
	saRangeIndex = 0;

	while (pos < numOfChar) {
			
		generatedPattern[pos]++;
		if (generatedPattern[pos] >= ALPHABET_SIZE) {
			pos++;
			continue;
		}
		saRangeIndex &= mask[pos];
		saRangeIndex |= generatedPattern[pos] << (pos * BIT_PER_CHAR);

		c = generatedPattern[pos];
		saIndexLeft[pos] = bwt->cumulativeFreq[c] + BWTOccValue(bwt, saIndexLeft[pos+1], c) + 1;
		saIndexRight[pos] = bwt->cumulativeFreq[c] + BWTOccValue(bwt, saIndexRight[pos+1] + 1, c);

		while (pos > 0 && saIndexLeft[pos] <= saIndexRight[pos]) {
			// Set preceding characters to 'a'
			pos--;
			generatedPattern[pos] = 0;
			saRangeIndex &= mask[pos];
			c = generatedPattern[pos];
			saIndexLeft[pos] = bwt->cumulativeFreq[c] + BWTOccValue(bwt, saIndexLeft[pos+1], c) + 1;
			saIndexRight[pos] = bwt->cumulativeFreq[c] + BWTOccValue(bwt, saIndexRight[pos+1] + 1, c);
		}

		if (saIndexLeft[pos] <= saIndexRight[pos]) {
			saIndexRange[saRangeIndex].startSaIndex = saIndexLeft[pos];
			saIndexRange[saRangeIndex].endSaIndex = saIndexRight[pos];
		}
	}

	fwrite(saIndexRange, sizeof(SaIndexRange), numOfPattern, saRangeFile);
	MMUnitFree(saIndexRange, numOfPattern * sizeof(SaIndexRange));


	fclose(saRangeFile);

}

void BWTGenerateSaBitmap(const BWT *bwt, const unsigned int numOfChar, const char *saBitmapFileName) {

	// stack
	unsigned int saIndexLeft[MAX_ARPROX_MATCH_LENGTH + 1];
	unsigned int saIndexRight[MAX_ARPROX_MATCH_LENGTH];
	unsigned char generatedPattern[MAX_ARPROX_MATCH_LENGTH];
	unsigned int pos;
	unsigned int i, c;
	unsigned long long numOfPattern;
	unsigned int numOfPatternExistInText = 0;
	unsigned int bitmapSize;
	unsigned int saBitmapIndex;
	unsigned int mask[16] = { 0xFFFFFFFC, 0xFFFFFFF3, 0xFFFFFFCF, 0xFFFFFF3F,
					 0xFFFFFCFF, 0xFFFFF3FF, 0xFFFFCFFF, 0xFFFF3FFF,
					 0xFFFCFFFF, 0xFFF3FFFF, 0xFFCFFFFF, 0xFF3FFFFF,
					 0xFCFFFFFF, 0xF3FFFFFF, 0xCFFFFFFF, 0x3FFFFFFF };

	unsigned int *saBitMap;

	FILE *saBitmapFile;

	saBitmapFile = (FILE*)fopen64(saBitmapFileName, "wb");
	if (saBitmapFile == NULL) {
		fprintf(stderr, "BWTGenerateSaBitmap(): Cannot open SA bitmap file!\n");
		exit(1);
	}

	fwrite(&bwt->inverseSa0, sizeof(unsigned int), 1, saBitmapFile);
	fwrite(bwt->cumulativeFreq + 1, sizeof(unsigned int), ALPHABET_SIZE, saBitmapFile);
	fwrite(&numOfChar, sizeof(unsigned int), 1, saBitmapFile);

	bitmapSize = 1 << (numOfChar * 2 - 5);	// 4^numOfChar
	saBitMap = MMUnitAllocate(bitmapSize * sizeof(unsigned int));
	for (i=0; i<bitmapSize; i++) {
		saBitMap[i] = 0;
	}

	// Set this boundary case so that the last character of generated pattern can be located by backward search
	saIndexLeft[numOfChar] = 0;
	saIndexRight[numOfChar] = bwt->textLength;

	// Set initial state
	pos = numOfChar - 1;
	generatedPattern[pos] = (unsigned char)-1;
	saBitmapIndex = 0;

	while (pos < numOfChar) {
			
		generatedPattern[pos]++;
		if (generatedPattern[pos] >= ALPHABET_SIZE) {
			pos++;
			continue;
		}
		saBitmapIndex &= mask[numOfChar - 1 - pos];
		saBitmapIndex |= generatedPattern[pos] << ((numOfChar - 1 - pos) * BIT_PER_CHAR);

		c = generatedPattern[pos];
		saIndexLeft[pos] = bwt->cumulativeFreq[c] + BWTOccValue(bwt, saIndexLeft[pos+1], c) + 1;
		saIndexRight[pos] = bwt->cumulativeFreq[c] + BWTOccValue(bwt, saIndexRight[pos+1] + 1, c);

		while (pos > 0 && saIndexLeft[pos] <= saIndexRight[pos]) {
			// Set preceding characters to 'a'
			pos--;
			generatedPattern[pos] = 0;
			saBitmapIndex &= mask[numOfChar - 1 - pos];
			c = generatedPattern[pos];
			saIndexLeft[pos] = bwt->cumulativeFreq[c] + BWTOccValue(bwt, saIndexLeft[pos+1], c) + 1;
			saIndexRight[pos] = bwt->cumulativeFreq[c] + BWTOccValue(bwt, saIndexRight[pos+1] + 1, c);
		}

		if (saIndexLeft[pos] <= saIndexRight[pos]) {
			saBitMap[saBitmapIndex >> 5] |= FIRST_BIT_MASK >> (saBitmapIndex & 0x1F);
			numOfPatternExistInText++;
		}
	}

	fwrite(saBitMap, sizeof(unsigned int), bitmapSize, saBitmapFile);
	MMUnitFree(saBitMap, bitmapSize * sizeof(unsigned int));

	numOfPattern = (unsigned long long)1 << (numOfChar * 2);	// 4^numOfChar

	printf("%u out of %llu length %u patterns exist in text.\n", numOfPatternExistInText, numOfPattern, numOfChar);

	fclose(saBitmapFile);

}

void BWTGenerateCompressedSaBitmap(const BWT *bwt, const unsigned int numOfChar, const char *saBitmapFileName) {

	// stack
	unsigned int saIndexLeft[MAX_ARPROX_MATCH_LENGTH + 1];
	unsigned int saIndexRight[MAX_ARPROX_MATCH_LENGTH];
	unsigned char generatedPattern[MAX_ARPROX_MATCH_LENGTH];
	unsigned int pos;
	unsigned int i, j, k, c;
	unsigned long long numOfPattern;
	unsigned int numOfPatternExistInText = 0;
	unsigned int bitmapSize;
	unsigned long long saBitmapIndex;
	unsigned long long mask[32];

	unsigned int *saBitMap;

	FILE *saBitmapFile;

	unsigned long long numOfBit = 0, numOfUnaryBit = 0, numOfGammaBit = 0;
	unsigned long long numOfOneUnaryBit = 0, numOfOneGammaBit = 0;
	unsigned long long numOfZeroUnaryBit = 0, numOfZeroGammaBit = 0;
	unsigned int lastBitValue, bitValue;
	unsigned int groupSize;

	unsigned int tempIntegerStream[1025];
	unsigned int numOfInteger;
	unsigned int numOfHighBit, numOfLowBit;

	unsigned int tempOneStream[1025];
	unsigned int tempZeroStream[1025];
	unsigned int numOfOne, numOfZero;


	for (i=0; i<32; i++) {
		mask[i] = ~((unsigned long long)(0x3) << (i*2));
	}


	saBitmapFile = (FILE*)fopen64(saBitmapFileName, "wb");
	if (saBitmapFile == NULL) {
		fprintf(stderr, "BWTGenerateCompressedSaBitmap(): Cannot open SA compressed bitmap file!\n");
		exit(1);
	}

	fwrite(&bwt->inverseSa0, sizeof(unsigned int), 1, saBitmapFile);
	fwrite(bwt->cumulativeFreq + 1, sizeof(unsigned int), ALPHABET_SIZE, saBitmapFile);
	fwrite(&numOfChar, sizeof(unsigned int), 1, saBitmapFile);


	bitmapSize = 1 << (numOfChar * 2 - 5);	// 4^numOfChar
	saBitMap = MMUnitAllocate(bitmapSize * sizeof(unsigned int));
	for (i=0; i<bitmapSize; i++) {
		saBitMap[i] = 0;
	}

	// Set this boundary case so that the last character of generated pattern can be located by backward search
	saIndexLeft[numOfChar] = 0;
	saIndexRight[numOfChar] = bwt->textLength;

	// Set initial state
	pos = numOfChar - 1;
	generatedPattern[pos] = (unsigned char)-1;
	saBitmapIndex = 0;

	while (pos < numOfChar) {
			
		generatedPattern[pos]++;
		if (generatedPattern[pos] >= ALPHABET_SIZE) {
			pos++;
			continue;
		}
		saBitmapIndex &= mask[numOfChar - 1 - pos];
		saBitmapIndex |= (unsigned long long)generatedPattern[pos] << ((numOfChar - 1 - pos) * BIT_PER_CHAR);

		c = generatedPattern[pos];
		saIndexLeft[pos] = bwt->cumulativeFreq[c] + BWTOccValue(bwt, saIndexLeft[pos+1], c) + 1;
		saIndexRight[pos] = bwt->cumulativeFreq[c] + BWTOccValue(bwt, saIndexRight[pos+1] + 1, c);

		while (pos > 0 && saIndexLeft[pos] <= saIndexRight[pos]) {
			// Set preceding characters to 'a'
			pos--;
			generatedPattern[pos] = 0;
			saBitmapIndex &= mask[numOfChar - 1 - pos];
			c = generatedPattern[pos];
			saIndexLeft[pos] = bwt->cumulativeFreq[c] + BWTOccValue(bwt, saIndexLeft[pos+1], c) + 1;
			saIndexRight[pos] = bwt->cumulativeFreq[c] + BWTOccValue(bwt, saIndexRight[pos+1] + 1, c);
		}

		if (saIndexLeft[pos] <= saIndexRight[pos]) {
			saBitMap[saBitmapIndex >> 5] |= FIRST_BIT_MASK >> ((unsigned int)saBitmapIndex & 0x1F);
			numOfPatternExistInText++;
		}
	}

	// Compress bitmap
	lastBitValue = saBitMap[0] >> (BITS_IN_WORD - 1);
	groupSize = 0;
	numOfInteger = 0;
	tempIntegerStream[0] = 0;

	tempOneStream[0] = 0;
	tempZeroStream[0] = 0;
	numOfOne = 0;
	numOfZero = 0;

	for (i=0; i<bitmapSize; i++) {
		c = saBitMap[i];
		for (j=0; j<BITS_IN_WORD; j++) {
			bitValue = c >> (BITS_IN_WORD - 1);
			c <<= 1;
			if (bitValue == lastBitValue) {
				groupSize++;
			} else {
				if (numOfInteger >= 1024) {
					for (k=1; k<=1024; k++) {
						tempIntegerStream[k] += tempIntegerStream[k-1];
					}
					numOfHighBit = 10;
					numOfLowBit = ceilLog2(tempIntegerStream[1024]) - numOfHighBit;
					for (k=1; k<1024; k++) {
						tempIntegerStream[k] >>= numOfLowBit;
						numOfUnaryBit += tempIntegerStream[k] - tempIntegerStream[k-1] + 1;
						numOfGammaBit += 2 * floorLog2(tempIntegerStream[k] - tempIntegerStream[k-1] + 1) + 1;
					}
					numOfUnaryBit += numOfLowBit * 1024;
					numOfGammaBit += numOfLowBit * 1024;
					numOfInteger = 0;
				}
				tempIntegerStream[numOfInteger+1] = groupSize;
				numOfInteger++;

				// One
				if (lastBitValue == 1) {
					if (numOfOne >= 1024) {
						for (k=1; k<=1024; k++) {
							tempOneStream[k] += tempOneStream[k-1];
						}
						numOfHighBit = 10;
						numOfLowBit = ceilLog2(tempOneStream[1024]) - numOfHighBit;
						for (k=1; k<1024; k++) {
							tempOneStream[k] >>= numOfLowBit;
							numOfOneUnaryBit += tempOneStream[k] - tempOneStream[k-1] + 1;
							numOfOneGammaBit += 2 * floorLog2(tempOneStream[k] - tempOneStream[k-1] + 1) + 1;
						}
						numOfOneUnaryBit += numOfLowBit * 1024;
						numOfOneGammaBit += numOfLowBit * 1024;
						numOfOne = 0;
					}
					tempOneStream[numOfOne+1] = groupSize;
					numOfOne++;
				} else {
					// Zero
					if (numOfZero >= 1024) {
						for (k=1; k<=1024; k++) {
							tempZeroStream[k] += tempZeroStream[k-1];
						}
						numOfHighBit = 10;
						numOfLowBit = ceilLog2(tempZeroStream[1024]) - numOfHighBit;
						for (k=1; k<1024; k++) {
							tempZeroStream[k] >>= numOfLowBit;
							numOfZeroUnaryBit += tempZeroStream[k] - tempZeroStream[k-1] + 1;
							numOfZeroGammaBit += 2 * floorLog2(tempZeroStream[k] - tempZeroStream[k-1] + 1) + 1;
						}
						numOfZeroUnaryBit += numOfLowBit * 1024;
						numOfZeroGammaBit += numOfLowBit * 1024;
						numOfZero = 0;
					}
					tempZeroStream[numOfZero+1] = groupSize;
					numOfZero++;
				}


				numOfBit += 2 * floorLog2(groupSize) + 1;
				groupSize = 1;
				lastBitValue = bitValue;
			}
		}
	}

	printf("Gamma only.\n");
	printf("Compressed bitmap size = %llu bits.\n", numOfBit);
	printf("Compressed bitmap size = %.2f bytes.\n", (double)numOfBit / BITS_IN_BYTE);

	printf("Rice + unary.\n");
	printf("Compressed bitmap size = %llu bits.\n", numOfUnaryBit);
	printf("Compressed bitmap size = %.2f bytes.\n", (double)numOfUnaryBit / BITS_IN_BYTE);

	printf("Rice + gamma.\n");
	printf("Compressed bitmap size = %llu bits.\n", numOfGammaBit);
	printf("Compressed bitmap size = %.2f bytes.\n", (double)numOfGammaBit / BITS_IN_BYTE);

	printf("Separate Ones and Zeros Rice + unary.\n");
	printf("Compressed bitmap size = %llu bits.\n", numOfOneUnaryBit + numOfZeroUnaryBit);
	printf("Compressed bitmap size = %.2f bytes.\n", (double)(numOfOneUnaryBit + numOfZeroUnaryBit) / BITS_IN_BYTE);

	printf("Separate Ones and Zeros Rice + gamma.\n");
	printf("Compressed bitmap size = %llu bits.\n", numOfOneGammaBit + numOfZeroGammaBit);
	printf("Compressed bitmap size = %.2f bytes.\n", (double)(numOfOneGammaBit + numOfZeroGammaBit) / BITS_IN_BYTE);

	fwrite(saBitMap, sizeof(unsigned int), bitmapSize, saBitmapFile);
	MMUnitFree(saBitMap, bitmapSize * sizeof(unsigned int));

	numOfPattern = (unsigned long long)1 << (numOfChar * 2);	// 4^numOfChar

	printf("%u out of %llu length %u patterns exist in text.\n", numOfPatternExistInText, numOfPattern, numOfChar);

	fclose(saBitmapFile);

}

void BWTCountPattern(const BWT *bwt, const unsigned int numOfChar) {

	// stack
	unsigned int saIndexLeft[MAX_ARPROX_MATCH_LENGTH + 1];
	unsigned int saIndexRight[MAX_ARPROX_MATCH_LENGTH];
	unsigned char generatedPattern[MAX_ARPROX_MATCH_LENGTH];
	unsigned int pos;
	unsigned int c;
	unsigned long long numOfPattern;
	unsigned int numOfPatternExistInText = 0;

	// Set this boundary case so that the last character of generated pattern can be located by backward search
	saIndexLeft[numOfChar] = 0;
	saIndexRight[numOfChar] = bwt->textLength;

	// Set initial state
	pos = numOfChar - 1;
	generatedPattern[pos] = (unsigned char)-1;

	while (pos < numOfChar) {
			
		generatedPattern[pos]++;
		if (generatedPattern[pos] >= ALPHABET_SIZE) {
			pos++;
			continue;
		}

		c = generatedPattern[pos];
		saIndexLeft[pos] = bwt->cumulativeFreq[c] + BWTOccValue(bwt, saIndexLeft[pos+1], c) + 1;
		saIndexRight[pos] = bwt->cumulativeFreq[c] + BWTOccValue(bwt, saIndexRight[pos+1] + 1, c);

		while (pos > 0 && saIndexLeft[pos] <= saIndexRight[pos]) {
			// Set preceding characters to 'a'
			pos--;
			generatedPattern[pos] = 0;
			c = generatedPattern[pos];
			saIndexLeft[pos] = bwt->cumulativeFreq[c] + BWTOccValue(bwt, saIndexLeft[pos+1], c) + 1;
			saIndexRight[pos] = bwt->cumulativeFreq[c] + BWTOccValue(bwt, saIndexRight[pos+1] + 1, c);
		}

		if (saIndexLeft[pos] <= saIndexRight[pos]) {
			numOfPatternExistInText++;
		}
	}

	numOfPattern = (unsigned long long)1 << (numOfChar * 2);	// 4^numOfChar

	printf("%u out of %llu length %u patterns exist in text.\n", numOfPatternExistInText, numOfPattern, numOfChar);

}


/*
void BWTGenerateSaRangeTable(const BWT *bwt, const unsigned int numOfChar, const char *PackedDNAFileName, const char *saRangeFileName, const char *startSaIndexFileName) {

	// stack
	unsigned int saIndexLeft[MAX_ARPROX_MATCH_LENGTH + 1];
	unsigned int saIndexRight[MAX_ARPROX_MATCH_LENGTH];
	unsigned char generatedPattern[MAX_ARPROX_MATCH_LENGTH];
	unsigned int endOfTextPattern[MAX_ARPROX_MATCH_LENGTH];
	unsigned int endOfTextChar;
	unsigned int pos;
	unsigned int i, c;
	unsigned int numOfPattern;
	unsigned int saRangeIndex, lastSaRangeIndex;
	unsigned int mask[16] = { 0xFFFFFFFC, 0xFFFFFFF3, 0xFFFFFFCF, 0xFFFFFF3F,
					 0xFFFFFCFF, 0xFFFFF3FF, 0xFFFFCFFF, 0xFFFF3FFF,
					 0xFFFCFFFF, 0xFFF3FFFF, 0xFFCFFFFF, 0xFF3FFFFF,
					 0xFCFFFFFF, 0xF3FFFFFF, 0xCFFFFFFF, 0x3FFFFFFF };

	SaIndexRange *saIndexRange;
	unsigned int *startSaIndex;
	
	FILE *saRangeFile, *startSaIndexFile;
	FILE *packedDNAFile;

	unsigned int packedFileLen;
	unsigned char lastByteLength, dnaByte;
	unsigned int d;

	packedDNAFile = (FILE*)fopen64(PackedDNAFileName, "rb");
	if (packedDNAFile == NULL) {
		fprintf(stderr, "BWTGenerateSaRangeTable(): Cannot open packed DNA file!\n");
		exit(1);
	}
	fseek(packedDNAFile, -1, SEEK_END);
	packedFileLen = ftell(packedDNAFile);
	if ((int)packedFileLen < 0) {
		fprintf(stderr, "BWTGenerateSaRangeTable(): Cannot determine file length!\n");
		exit(1);
	}
	fread(&lastByteLength, sizeof(unsigned char), 1, packedDNAFile);

	fseek(packedDNAFile, -2, SEEK_CUR);
	fread(&dnaByte, sizeof(unsigned char), 1, packedDNAFile);
	d = dnaByte;
	endOfTextChar = (d >> ((CHAR_PER_BYTE - lastByteLength) * BIT_PER_CHAR));
	i = lastByteLength;

	while (i<numOfChar-1) {
		fseek(packedDNAFile, -2, SEEK_CUR);
		fread(&dnaByte, sizeof(unsigned char), 1, packedDNAFile);
		d =	dnaByte;
		endOfTextChar |= (d << (i * BIT_PER_CHAR));
		i += CHAR_PER_BYTE;
	}

	endOfTextChar &= truncateLeft(ALL_ONE_MASK, (numOfChar - 1) * BIT_PER_CHAR);

	fclose(packedDNAFile);


	// saRangeFile
	saRangeFile = (FILE*)fopen64(saRangeFileName, "wb");
	if (saRangeFile == NULL) {
		fprintf(stderr, "BWTGenerateSaRangeTable(): Cannot open SA range file!\n");
		exit(1);
	}

	fwrite(&bwt->inverseSa0, sizeof(unsigned int), 1, saRangeFile);
	fwrite(bwt->cumulativeFreq + 1, sizeof(unsigned int), ALPHABET_SIZE, saRangeFile);
	fwrite(&numOfChar, sizeof(unsigned int), 1, saRangeFile);

	// startSaIndexFile
	startSaIndexFile = (FILE*)fopen64(startSaIndexFileName, "wb");
	if (startSaIndexFile == NULL) {
		fprintf(stderr, "BWTGenerateSaRangeTable(): Cannot open Start SA index file!\n");
		exit(1);
	}

	fwrite(&bwt->inverseSa0, sizeof(unsigned int), 1, startSaIndexFile);
	fwrite(bwt->cumulativeFreq + 1, sizeof(unsigned int), ALPHABET_SIZE, startSaIndexFile);
	fwrite(&numOfChar, sizeof(unsigned int), 1, startSaIndexFile);


	numOfPattern = 1 << (numOfChar * 2);	// 4^numOfChar
	saIndexRange = MMUnitAllocate(numOfPattern * sizeof(SaIndexRange));
	for (i=0; i<numOfPattern; i++) {
		saIndexRange[i].startSaIndex = 0;
		saIndexRange[i].endSaIndex = 0;
	}

	// Set this boundary case so that the last character of generated pattern can be located by backward search
	saIndexLeft[numOfChar] = 0;
	saIndexRight[numOfChar] = bwt->textLength;

	// Set initial state
	pos = numOfChar - 1;
	generatedPattern[pos] = (unsigned char)-1;
	saRangeIndex = 0;

	while (pos < numOfChar) {
			
		generatedPattern[pos]++;
		if (generatedPattern[pos] >= ALPHABET_SIZE) {
			pos++;
			continue;
		}
		saRangeIndex &= mask[pos];
		saRangeIndex |= generatedPattern[pos] << (pos * BIT_PER_CHAR);

		c = generatedPattern[pos];
		saIndexLeft[pos] = bwt->cumulativeFreq[c] + BWTOccValue(bwt, saIndexLeft[pos+1], c) + 1;
		saIndexRight[pos] = bwt->cumulativeFreq[c] + BWTOccValue(bwt, saIndexRight[pos+1] + 1, c);

		while (pos > 0 && saIndexLeft[pos] <= saIndexRight[pos]) {
			// Set preceding characters to 'a'
			pos--;
			generatedPattern[pos] = 0;
			saRangeIndex &= mask[pos];
			c = generatedPattern[pos];
			saIndexLeft[pos] = bwt->cumulativeFreq[c] + BWTOccValue(bwt, saIndexLeft[pos+1], c) + 1;
			saIndexRight[pos] = bwt->cumulativeFreq[c] + BWTOccValue(bwt, saIndexRight[pos+1] + 1, c);
		}

		if (saIndexLeft[pos] <= saIndexRight[pos]) {
			saIndexRange[saRangeIndex].startSaIndex = saIndexLeft[pos];
			saIndexRange[saRangeIndex].endSaIndex = saIndexRight[pos];
		}

	}

	// Output saIndexRange
	fwrite(saIndexRange, sizeof(SaIndexRange), numOfPattern, saRangeFile);
	fclose(saRangeFile);

	MMUnitFree(saIndexRange, numOfPattern * sizeof(SaIndexRange));


	// Generate start SA index

	// Search again with normal order of search string

	for (i=0; i<numOfPattern; i++) {
		saIndexRange[i].startSaIndex = 0;
		saIndexRange[i].endSaIndex = 0;
	}

	// Set this boundary case so that the last character of generated pattern can be located by backward search
	saIndexLeft[numOfChar] = 0;
	saIndexRight[numOfChar] = bwt->textLength;

	// Set initial state
	pos = numOfChar - 1;
	generatedPattern[pos] = (unsigned char)-1;
	saRangeIndex = 0;

	while (pos < numOfChar) {
			
		generatedPattern[pos]++;
		if (generatedPattern[pos] >= ALPHABET_SIZE) {
			pos++;
			continue;
		}
		saRangeIndex &= mask[numberOfChar - 1 - pos];
		saRangeIndex |= generatedPattern[pos] << ((numberOfChar - 1 - pos) * BIT_PER_CHAR);

		c = generatedPattern[pos];
		saIndexLeft[pos] = bwt->cumulativeFreq[c] + BWTOccValue(bwt, saIndexLeft[pos+1], c) + 1;
		saIndexRight[pos] = bwt->cumulativeFreq[c] + BWTOccValue(bwt, saIndexRight[pos+1] + 1, c);

		while (pos > 0 && saIndexLeft[pos] <= saIndexRight[pos]) {
			// Set preceding characters to 'a'
			pos--;
			generatedPattern[pos] = 0;
			saRangeIndex &= mask[numberOfChar - 1 - pos];
			c = generatedPattern[pos];
			saIndexLeft[pos] = bwt->cumulativeFreq[c] + BWTOccValue(bwt, saIndexLeft[pos+1], c) + 1;
			saIndexRight[pos] = bwt->cumulativeFreq[c] + BWTOccValue(bwt, saIndexRight[pos+1] + 1, c);
		}

		if (saIndexLeft[pos] <= saIndexRight[pos]) {
			saIndexRange[saRangeIndex].startSaIndex = saIndexLeft[pos];
			saIndexRange[saRangeIndex].endSaIndex = saIndexRight[pos];
		}

	}

	startSaIndex = MMUnitAllocate(sizeof(unsigned int) * (numOfPattern + 1));
	for (i=0; i<numOfPattern; i++) {
		saRangeIndex = 0;
		// reverse character order
		for (pos=0; pos<numOfChar; pos++) {
			saRangeIndex |= (((i >> pos) & 0x3) << (numOfChar - pos));
		}
		startSaIndex[saRangeIndex] = saIndexRange[i].startSaIndex;
	}
	startSaIndex[numOfPattern] = bwt->textLength + 1;

	// Fill startSaIndex for pattern not found in text
	for (i=numOfPattern-1; i>0; i--) {
		if (startSaIndex[i] == 0) {
			startSaIndex[i] = startSaIndex[i+1];
		}
	}

	// Find the suffices in text with length < numOfChar

	// $
	saRangeIndex = 0;
	saIndexRange[saRangeIndex].startSaIndex--;
	endOfTextPattern[numOfChar - 1] = saRangeIndex;

	for (i=numOfChar - 1; i>0; i--) {
		for (c=0; c<ALPHABET_SIZE; c++) {
			saRangeIndex &= mask[i];
			saRangeIndex |= c << (i * BIT_PER_CHAR);
			if (saIndexRange[saRangeIndex].startSaIndex != 0 || saIndexRange[saRangeIndex].endSaIndex) {
				if (saRangeIndex != 0) {
					lastSaRangeIndex = saRangeIndex;
					pos = numOfChar - 1;
					while ((saRangeIndex & (~mask[pos])) == 0) {
						lastSaRangeIndex |= (ALPHABET_SIZE - 1) << (pos * BIT_PER_CHAR);
						pos--;
					}
					lastSaRangeIndex -= (1 << (pos * BIT_PER_CHAR));
					if (saIndexRange[saRangeIndex].startSaIndex != saIndexRange[lastSaRangeIndex].endSaIndex + 1) {
						break;
					}
				} else {
					if (saIndexRange[saRangeIndex].startSaIndex > 0) {
						break;
					}
				}
			}
		}
		saIndexRange[saRangeIndex].startSaIndex--;
		endOfTextPattern[i] = saRangeIndex;
	}

	// Output short suffices
	fwrite(endOfTextPattern, sizeof(unsigned int), numOfChar, saRangeFile);

	// Output startSaIndex
	for (i=0; i<numOfPattern; i++) {
		fwrite(&(saIndexRange[i].startSaIndex), sizeof(unsigned int), 1, saRangeFile);
	}
	MMUnitFree(saIndexRange, numOfPattern * sizeof(SaIndexRange));

	fclose(saRangeFile);

}
*/

void BWTPrintAndVerifyOccValue(const BWT *bwt, FILE *output) {

	unsigned int i, j;
	unsigned int c;
	unsigned int v;
	unsigned int occValue[ALPHABET_SIZE];

	for (i=0; i<ALPHABET_SIZE; i++) {
		occValue[i] = 0;
	}

	for (i=1; i<=bwt->textLength+1; i++) {
		if (i != bwt->inverseSa0) {
			v = BWTOccValueOnSpot(bwt, i, &c);
			if (output != NULL) {
				fprintf(output, "%u:  %u : %u\n", i, c, v);
			}
			occValue[c] = v;
			v = 0;
			for (j=0; j<ALPHABET_SIZE; j++) {
				v += occValue[j];
			}
			if ((i > bwt->inverseSa0 && v != (i-1)) || (i < bwt->inverseSa0 && v != i)) {
				fprintf(stderr, "BWTPrintAndVerifyOccValue(): occValue not summed to text length!\n");
				exit(1);
			}
		}
	}

}

void BWTPrintAndVerifySaValue(const BWT *bwt, FILE *output) {

	unsigned int i;
	unsigned int saValue;
	unsigned int *printedValue;

	printedValue = malloc(sizeof(unsigned int) * (bwt->textLength + 1));
	initializeVAL(printedValue, bwt->textLength + 1, FALSE);
    
	for (i=0; i<=bwt->textLength; i++) {
		saValue = BWTSaValue(bwt, i);
		if (output != NULL) {
			fprintf(output, "%u ", saValue);
		}
		printedValue[saValue] = TRUE;
	}

	if (output != NULL) {
		fprintf(output, "\n");
	}

	if (numberOfMatchInVAL(printedValue, bwt->textLength + 1, FALSE) > 0) {
		fprintf(stderr, "BWTPrintAndVerifySaValue() : some saValue is not printed!\n");
		exit(1);
	}

	free(printedValue);

}

void BWTPrintAndVerifyInverseSa(const BWT *bwt, FILE *output) {

	unsigned int i;
	unsigned int inverseSa;
	unsigned int *printedValue;

	printedValue = malloc(sizeof(unsigned int) * (bwt->textLength + 1));
	initializeVAL(printedValue, bwt->textLength + 1, FALSE);
    
	for (i=0; i<bwt->textLength+1; i++) {
		inverseSa = BWTInverseSa(bwt, i);
		if (output != NULL) {
			fprintf(output, "%u ", inverseSa);
		}
		printedValue[inverseSa] = TRUE;
	}

	if (output != NULL) {
		fprintf(output, "\n");
	}

	if (numberOfMatchInVAL(printedValue, bwt->textLength + 1, FALSE) > 0) {
		fprintf(stderr, "BWTPrintAndVerifyInverseSa() : some inverseSa is not printed!\n");
		exit(1);
	}

	free(printedValue);

}

#endif
