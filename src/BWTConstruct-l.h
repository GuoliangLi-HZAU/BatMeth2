#include "config.h"
#ifdef MMX
/*

   BWTConstruct.h		BWT-Index Construction

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

#ifndef __BWTCONSTRUCT_H__
#define __BWTCONSTRUCT_H__

#include "TypeNLimit-l.h"
#include "BWT-l.h"

#define BWTINC_INSERT_SORT_NUM_ITEM 7
#define BWTINC_MIN_MEMORY_IN_WORD	1048576	// 4M RAM 

typedef struct BWTInc {
	BWT *bwt;
	unsigned int numberOfIterationDone;
	unsigned int *cumulativeCountInCurrentBuild;
	unsigned int availableWord;
	unsigned int targetTextLength;
	float targetNBit;
	unsigned int buildSize;
	unsigned int initialMaxBuildSize;
	unsigned int incMaxBuildSize;
	unsigned int firstCharInLastIteration;
	unsigned int *workingMemory;
	unsigned int *packedText;
	unsigned char *textBuffer;
	unsigned int *packedShift;
} BWTInc;

// This is how the workingMemory is used (initial build)

// 1. When ABS Rank is generated and sorting is being done
// |-----------Seq-----------|--------ABS Rank--------|---Text----|
//            4nbyte                   4nbyte             2nbit
// 2. When seq are regenerated and BWT is generated
// |-----------Seq-----------|--------ABS Rank--------|----BWT----|
//            4nbyte                   4nbyte             2nbit


// This is how the workingMemory is used (incremental build)

// 1. When ABS Rank is generated
// |--ABS Rank--|----Seq-----|----Text----|----Occ----|----BWT----|
//     4nbyte      4nbyte     4nbyte(2nbit)
// 2. When sorting is being done
// |sorted ARank|----Seq-----|--REL-Rank--|----Occ----|----BWT----|
//     4nbyte      4nbyte        4nbyte
// 3. When new BWT is built
// |sorted ARank|---new BWT--|--REL-Rank--|----Occ----|----BWT----|
//     4nbyte      4nbyte        4nbyte
// 4. When insertion is being done
// |sorted ARank|---new BWT--|---new Occ---|------merged BWT------|
//     4nbyte      4nbyte     (4nbyte freed for new Occ and merged BWT)

// BWT Construction
BWTInc *BWTIncCreate(MMPool *mmPool, const unsigned int textLength, const float targetNBit, 
					 const unsigned int initialMaxBuildSize, const unsigned int incMaxBuildSize);
void BWTIncFree(MMPool *mmPool, BWTInc *fmiInc);
BWTInc *BWTIncConstructFromPacked(MMPool *mmPool, const char *inputFileName, const unsigned int showProgress,
								  const float targetNBit, const unsigned int initialMaxBuildSize, const unsigned int incMaxBuildSize);
unsigned int BWTGenerateOccValueToFileFromBwt(const char *bwtFileName, const char *occValueFileName, unsigned int*  decodeTable);
void BWTGenerateOccValueFromBwt(const unsigned int*  bwt, unsigned int* __restrict occValue, unsigned int* __restrict occValueMajor,
								const unsigned int textLength, const unsigned int*  decodeTable);

// Auxliary data structures
void BWTSaveBwtCodeAndOcc(const BWT *bwt, const char *bwtFileName, const char *occValueFileName);
void BWTGenerateSaValue(BWT *bwt, const unsigned int saValueFreq, unsigned int showProgressInterval);
void BWTGenerateFullSaValue(BWT *bwt);
void BWTSaveSaValue(const BWT *bwt, const char *saValueFileName);
void BWTGenerateInverseSa(BWT *bwt, const unsigned int inverseSaFreq, unsigned int showProgressInterval);
void BWTSaveInverseSa(const BWT *bwt, const char *inverseSaFileName);
void BWTGenerateCachedSaIndex(const BWT *bwt, const unsigned int numOfChar, const char *cachedSaIndexFileName);
void BWTGenerateSaBitmap(const BWT *bwt, const unsigned int numOfChar, const char *saBitmapFileName);
void BWTGenerateCompressedSaBitmap(const BWT *bwt, const unsigned int numOfChar, const char *saBitmapFileName);
void BWTCountPattern(const BWT *bwt, const unsigned int numOfChar);

// verification functions
void BWTPrintAndVerifyOccValue(const BWT *bwt, FILE *output);
void BWTPrintAndVerifySaValue(const BWT *bwt, FILE *output);
void BWTPrintAndVerifyInverseSa(const BWT *bwt, FILE *output);





#endif
#else

/*

   BWTConstruct.h		BWT-Index Construction

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

#ifndef __BWTCONSTRUCT_H__
#define __BWTCONSTRUCT_H__

#include "TypeNLimit-l.h"
#include "BWT-l.h"

#define BWTINC_INSERT_SORT_NUM_ITEM 7


typedef struct BWTInc {
	BWT *bwt;
	unsigned int numberOfIterationDone;
	unsigned int *cumulativeCountInCurrentBuild;
	unsigned int availableWord;
	unsigned int targetTextLength;
	float targetNBit;
	unsigned int buildSize;
	unsigned int initialMaxBuildSize;
	unsigned int incMaxBuildSize;
	unsigned int firstCharInLastIteration;
	unsigned int *workingMemory;
	unsigned int *packedText;
	unsigned char *textBuffer;
	unsigned int *packedShift;
} BWTInc;

// This is how the workingMemory is used (initial build)

// 1. When ABS Rank is generated and sorting is being done
// |-----------Seq-----------|--------ABS Rank--------|---Text----|
//            4nbyte                   4nbyte             2nbit
// 2. When seq are regenerated and BWT is generated
// |-----------Seq-----------|--------ABS Rank--------|----BWT----|
//            4nbyte                   4nbyte             2nbit


// This is how the workingMemory is used (incremental build)

// 1. When ABS Rank is generated
// |--ABS Rank--|----Seq-----|----Text----|----Occ----|----BWT----|
//     4nbyte      4nbyte     4nbyte(2nbit)
// 2. When sorting is being done
// |sorted ARank|----Seq-----|--REL-Rank--|----Occ----|----BWT----|
//     4nbyte      4nbyte        4nbyte
// 3. When new BWT is built
// |sorted ARank|---new BWT--|--REL-Rank--|----Occ----|----BWT----|
//     4nbyte      4nbyte        4nbyte
// 4. When insertion is being done
// |sorted ARank|---new BWT--|---new Occ---|------merged BWT------|
//     4nbyte      4nbyte     (4nbyte freed for new Occ and merged BWT)

// BWT Construction
BWTInc *BWTIncCreate(MMPool *mmPool, const unsigned int textLength, const float targetNBit, 
					 const unsigned int initialMaxBuildSize, const unsigned int incMaxBuildSize);
void BWTIncFree(MMPool *mmPool, BWTInc *fmiInc);
BWTInc *BWTIncConstructFromPacked(MMPool *mmPool, const char *inputFileName, const unsigned int showProgress,
								  const float targetNBit, const unsigned int initialMaxBuildSize, const unsigned int incMaxBuildSize);
unsigned int BWTGenerateOccValueToFileFromBwt(const char *bwtFileName, const char *occValueFileName, unsigned int*  decodeTable);
void BWTGenerateOccValueFromBwt(const unsigned int*  bwt, unsigned int* __restrict occValue, unsigned int* __restrict occValueMajor,
								const unsigned int textLength, const unsigned int*  decodeTable);

// Auxliary data structures
void BWTSaveBwtCodeAndOcc(const BWT *bwt, const char *bwtFileName, const char *occValueFileName);
void BWTGenerateSaValue(MMPool *mmPool, BWT *bwt, const unsigned int saValueFreq, unsigned int showProgressInterval);
void BWTSaveSaValue(const BWT *bwt, const char *saValueFileName);
void BWTGenerateInverseSa(BWT *bwt, const unsigned int inverseSaFreq, unsigned int showProgressInterval);
void BWTSaveInverseSa(const BWT *bwt, const char *inverseSaFileName);
void BWTGenerateSaRangeTable(const BWT *bwt, const unsigned int numOfChar, const char *saRangeFileName);
void BWTGenerateSaBitmap(const BWT *bwt, const unsigned int numOfChar, const char *saBitmapFileName);
void BWTGenerateCompressedSaBitmap(const BWT *bwt, const unsigned int numOfChar, const char *saBitmapFileName);
void BWTCountPattern(const BWT *bwt, const unsigned int numOfChar);

// verification functions
void BWTPrintAndVerifyOccValue(const BWT *bwt, FILE *output);
void BWTPrintAndVerifySaValue(const BWT *bwt, FILE *output);
void BWTPrintAndVerifyInverseSa(const BWT *bwt, FILE *output);





#endif
#endif
