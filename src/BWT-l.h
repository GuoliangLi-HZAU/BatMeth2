#include "config.h"

#ifdef MMX
/*

   BWT.h	BWT-Index

   This module contains an implementation of BWT-index for alphabet size = 4.
   The functions provided include:
    Load functions for loading BWT to memory;
    Core functions for accessing core Inverse Psi values;
	Search functions for searching patterns from text;
	Text retrieval functions for retrieving text from BWT.

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

#ifndef __BWT_H__
#define __BWT_H__

#include "TypeNLimit-l.h"
#include "MemManager-l.h"
#include "TextConverter-l.h"
#include "HSP-l.h"

#define BITS_PER_OCC_VALUE			16
#define OCC_VALUE_PER_WORD			2
#define OCC_INTERVAL				256
#define WORD_BETWEEN_OCC			16
#define OCC_INTERVAL_MAJOR			65536

#define SORT_ALL					0
#define SORT_16_BIT					1
#define SORT_NONE					2

#define BUCKET_BIT					16
#define NUM_BUCKET					65536

#define MAX_APPROX_MATCH_ERROR	7
#define MAX_ARPROX_MATCH_LENGTH	64

#define BWTDP_MAX_SUBSTRING_LENGTH	512

#define ESTIMATED_OCC_DIFF			32	// 128 / 4
#define MAX_OCC_DIFF				128



typedef struct BWT {
	unsigned int textLength;			// length of the text
	unsigned int saInterval;			// interval between two SA values stored explicitly
	unsigned int inverseSaInterval;		// interval between two inverse SA stored explicitly
	unsigned int inverseSa0;			// SA-1[0]
	unsigned int *cumulativeFreq;		// cumulative frequency
	unsigned int *bwtCode;				// BWT code
	unsigned int *occValue;				// Occurrence values stored explicitly
	unsigned int *occValueMajor;		// Occurrence values stored explicitly
	unsigned int *saValue;				// SA values stored explicitly
	unsigned int *inverseSa;			// Inverse SA stored explicitly
	unsigned int *cachedSaIndex;		// Cached SA index
	unsigned int cachedSaIndexNumOfChar;	// Number of characters indexed in SA index range
	unsigned int *saValueOnBoundary;	// Pre-calculated frequently referred data
	unsigned int *decodeTable;			// For decoding BWT by table lookup
	unsigned int decodeTableGenerated;	// == TRUE if decode table is generated on load and will be freed
	unsigned int bwtSizeInWord;			// Temporary variable to hold the memory allocated
	unsigned int occSizeInWord;			// Temporary variable to hold the memory allocated
	unsigned int occMajorSizeInWord;	// Temporary variable to hold the memory allocated
	unsigned int saValueSizeInWord;		// Temporary variable to hold the memory allocated
	unsigned int inverseSaSizeInWord;	// Temporary variable to hold the memory allocated
	unsigned int cachedSaIndexSizeInWord;	// Temporary variable to hold the memory allocated
} BWT;

#define MAX_DIAGONAL_LEVEL 4				// Number of sub-pattern to keep for detecting diagonal hit

// Error information is stored as:
// 1. bitVector
//	  After hamming distance match
// 2. count
//    After edit distance match
// 3. score
//    After the hits are processed with scoring functions

typedef struct SaIndexGroupNew {	// SA index range and information of a particular error arrangement of a matched sub-pattern
	unsigned int startSaIndex;			// starting SA index
	unsigned int numOfMatch;				// number of match
	unsigned int posQuery;				// position in query; used for detecting diagonal hits
	unsigned int info;					// extra hit information; to be copied to hitList.info
} SaIndexGroupNew;

typedef struct SaIndexGroupTemp {	// SA index range and information of a particular error arrangement of a matched sub-pattern
	unsigned int startSaIndex1;			// starting SA index
	unsigned int numOfMatch1;			// number of match
	unsigned int startSaIndex2;			// position in query; used for detecting diagonal hits
	unsigned int numOfMatch2;			// extra hit information; to be copied to hitList.info
} SaIndexGroupTemp;

typedef struct SaIndexGroupOld {	// SA index range and information of a particular error arrangement of a matched sub-pattern
	unsigned int startSaIndex;			// starting SA index
	unsigned int numOfMatch;				// number of match
	unsigned int info;					// extra hit information; to be copied to hitList.info
} SaIndexGroupOld;

typedef struct SaIndexGroup {	// SA index range and information of a particular error arrangement of a matched sub-pattern
	unsigned int startSaIndex;			// starting SA index
	unsigned int numOfMatch;			// number of match
	unsigned int info;					// extra hit information
} SaIndexGroup;

typedef struct SaIndexGroupWithErrorBitVector {	// SA index range and information of a particular error arrangement of a matched sub-pattern
	unsigned int startSaIndex;			// starting SA index
	unsigned int numOfMatch;			// number of match
	unsigned int errorBitVector;			// error bit vector
} SaIndexGroupWithErrorBitVector;

typedef struct SaIndexGroupWithLengthError {	// SA index range and information of a particular error arrangement of a matched sub-pattern
	unsigned int startSaIndex;			// starting SA index
	unsigned int numOfMatch;			// number of match
	unsigned posQuery : 16;		// position in query
	unsigned length   : 8;		// length of hit
	unsigned error    : 8;		// error in hit
} SaIndexGroupWithLengthError;

typedef struct SaIndexGroupProcessed {	// Alternative usage of SaIndexGroup - once processed, error bit vector is replaced by index to text position
	unsigned int startSaIndex;			// starting SA index
	unsigned int numOfMatch;			// number of match
	unsigned int textPositionIndex;		// storing the pointer to text position
} SaIndexGroupProcessed;

typedef struct DupSaIndexGroup {	// Alternative usage of SaIndexGroup - the group duplicates another group
	unsigned int lastDupSaIndexGroupIndex;	// index to last duplicated group
	unsigned int saIndexGroupIndex;			// index to the first SA into group among the duplicates
	unsigned int textPositionIndex;			// storing the pointer to text position
} DupSaIndexGroup;

typedef struct SaIndexGroupHash {	// Hash table for checking duplicate SA index group
	unsigned int startSaIndex;
	unsigned int saIndexGroupIndex;
} SaIndexGroupHash;

typedef struct BWTSaRetrievalStatistics {
	unsigned int bwtSaRetrieved;
	unsigned int saDiagonalLinked;
	unsigned int saDuplicated;
	unsigned int cachedSaRetrieved;
} BWTSaRetrievalStatistics;

typedef struct BWTDPStatistics {
	int maxDepth;
	int maxDPCell;
	int maxDPMemoryInWord;
	int totalMaxDepth;
	int totalMaxDPCell;
	int totalMaxDPMemoryInWord;
	LONG acceptedPathDepth;
	LONG acceptedPath;
	LONG rejectedPathDepth;
	LONG rejectedPath;
	LONG* __restrict totalNode;
	LONG* __restrict rejectedNode;
	LONG* __restrict totalDPCell;
} BWTDPStatistics;

typedef struct SaIndexList {
	unsigned int saIndex;
	unsigned int textPositionIndex;
} SaIndexList;

typedef struct HitCombination {
	int numOfCombination;
	int maxError;
	int keyLength;
	int skipTableWidth;
	int *errorPos;
	int *skip;
	int *skipErrorIndex;
} HitCombination;

typedef struct DPText {
	int charBeingProcessed;
	int dpCellIndex;
	int numOfDpCellSegment;
	unsigned int dummy1;	// Must not be removed; so that saIndexLeft and saIndexRight are aligned to 16 byte boundary
	unsigned int saIndexLeft[ALPHABET_SIZE];
	unsigned int saIndexRight[ALPHABET_SIZE];
} DPText;

typedef struct DPScanDepth {
	unsigned P				:	31;
	unsigned withAmbiguity	:	1;
} DPScanDepth;


// Load / unload functions
BWT *BWTCreate(MMPool *mmPool, const unsigned int textLength, unsigned int *decodeTable);
BWT *BWTLoad(MMPool *mmPool, const char *bwtCodeFileName, const char *occValueFileName, 
			 const char *saValueFileName, const char *inverseSaFileName, const char *saIndexRangeFileName,
			 unsigned int *decodeTable);
void BWTFree(MMPool *mmPool, BWT *bwt);
void BWTPrintMemoryUsage(const BWT *bwt, FILE *output, const unsigned int packedDNASize);

// Precalculate frequenctly accessed data
void BWTGenerateSaValueOnBoundary(MMPool *mmPool, BWT *bwt);

// Core functions
// The following must be customized for differenet compression schemes ***
unsigned int BWTDecode(const BWT *bwt, const unsigned int index1, const unsigned int index2, const unsigned int character);
void BWTDecodeAll(const BWT *bwt, const unsigned int index1, const unsigned int index2, unsigned int* __restrict occValue);
unsigned int BWTOccValue(const BWT *bwt, unsigned int index, const unsigned int character);
void BWTOccValueTwoIndex(const BWT *bwt, unsigned int index1, unsigned int index2, const unsigned int character, unsigned int* __restrict occValue);
void BWTAllOccValue(const BWT *bwt, unsigned int index, unsigned int* __restrict occValue);
void BWTAllOccValueTwoIndex(const BWT *bwt, unsigned int index1, unsigned int index2, unsigned int* __restrict occValue1, unsigned int* __restrict occValue2);
unsigned int BWTOccValueOnSpot(const BWT *bwt, unsigned int index, unsigned int* __restrict character);
unsigned int BWTSearchOccValue(const BWT *bwt, const unsigned int character, const unsigned int searchOccValue);


// Utility functions for no compression only
unsigned int BWTResidentSizeInWord(const unsigned int numChar);
unsigned int BWTFileSizeInWord(const unsigned int numChar);
void BWTClearTrailingBwtCode(BWT *bwt);

// These are generic to different compression schemes (and generic to no compression as well)
unsigned int BWTPsiMinusValue(const BWT *bwt, const unsigned int index);
unsigned int BWTPsiPlusValue(const BWT *bwt, const unsigned int index);
unsigned int BWTSaValue(const BWT *bwt, unsigned int index);
unsigned int BWTInverseSa(const BWT *bwt, unsigned int saValue);
unsigned int BWTOccIntervalMajor(const unsigned int occInterval);
unsigned int BWTOccValueMinorSizeInWord(const unsigned int numChar);
unsigned int BWTOccValueMajorSizeInWord(const unsigned int numChar);

// Search functions
// packedText should be allocated with at least 1 Word buffer initialized to zero
int BWTForwardSearch(const unsigned int *packedKey, const unsigned int keyLength, const BWT *bwt, const unsigned int *packedText);
int BWTForwardSearchSaIndex(const unsigned int *packedKey, const unsigned int keyLength, const BWT *bwt, const unsigned int *packedText, 
					 unsigned int *resultSaIndexLeft, unsigned int *resultSaIndexRight);
int BWTSaBinarySearch(const unsigned char *convertedKey, const unsigned int keyLength, const BWT *bwt, const unsigned int *packedText, 
					  unsigned int *resultSaIndexLeft, unsigned int *resultSaIndexRight, unsigned int *tempKey);	// tempKey = buffer large enough to hold packed key
int BWTBackwardSearch(const unsigned char *convertedKey, const unsigned int keyLength, const BWT *bwt, 
					  unsigned int *resultSaIndexLeft, unsigned int *resultSaIndexRight);
int BWTBackwardSearchCheckWithText(const unsigned char *convertedKey, const unsigned int *packedKey, const unsigned int keyLength,
							  const BWT *bwt, const unsigned int *packedText, const unsigned int textCheckingCostFactor,
							  const unsigned int maxnumOfTextPosition, HitList* __restrict hitList,
							  unsigned int *resultSaIndexLeft, unsigned int *resultSaIndexRight);

// Approximate match functions - brute force deep first search by backward search is used
unsigned int BWTHammingDistMaxSaIndexGroup(const unsigned int keyLength, const unsigned int maxError);
unsigned int BWTHammingDistCountOcc(const unsigned char *convertedKey, const unsigned int keyLength, const BWT *bwt, const unsigned int maxError);
unsigned int BWTHammingDistMatch(const BWT *bwt, const unsigned char *convertedKey, const HitCombination *hitCombination,
						const unsigned int *cachedSaIndex, const unsigned int cachedSaIndexNumOfChar,
						SaIndexGroupNew* __restrict saIndexGroup, const unsigned int maxSaIndexGroup);
int BWTHammingDistMatchOld(const unsigned char *convertedKey, const unsigned int keyLength, const BWT *bwt, const unsigned int maxError,
						   SaIndexGroupNew* __restrict saIndexGroup, const unsigned int maxSaIndexGroup,
						   const unsigned int posQuery, const unsigned int info);
unsigned int BWTEditDistMaxSaIndexGroup(const unsigned int keyLength, const unsigned int maxError);
// Does not insert characters on pattern boundary
unsigned int BWTEditDistMatch(const unsigned char *convertedKey, const unsigned int keyLength, const BWT *bwt, const unsigned int maxError,
					 SaIndexGroupWithLengthError* __restrict saIndexGroup, const unsigned int maxSaIndexGroup);
unsigned int BWTEditDistMatchOld(const unsigned char *convertedKey, const unsigned int keyLength, const BWT *bwt, const unsigned int maxError,
					 SaIndexGroupWithLengthError* __restrict saIndexGroup, const unsigned int maxSaIndexGroup);


unsigned int BWTEliminateDupSaIndexGroup(SaIndexGroupWithLengthError* __restrict saIndexGroup, const unsigned int numOfSaGroup);

unsigned int BWTSubPatternHammingDistCountOcc(const unsigned char *convertedKey, const unsigned int keyLength, const unsigned int subPatternLength, const BWT *bwt, 
								     const unsigned int maxError, const unsigned int skip);
int BWTSubPatternHammingDistSaIndex(const BWT *bwt, const unsigned char *convertedKey, const int keyLength, const int skip, 
									const HitCombination *hitCombination, 
									const unsigned int *cachedSaIndex, const unsigned int cachedSaIndexNumOfChar,
									SaIndexGroupNew* __restrict saIndexGroup, const int maxnumOfSaIndexGroup,
									int* __restrict firstSaIndexGroupForSubPattern);
int BWTSubPatternHammingDistSaIndexOld(const unsigned char *convertedKey, const int keyLength, const int subPatternLength, const BWT *bwt, 
								    const int maxError, const int skip, const int lengthProcessed, int* __restrict lengthInCurrentRound,
									SaIndexGroupNew* __restrict saIndexGroup, const int maxnumOfSaIndexGroup);
int BWTDPHit(const BWT *bwt, SaIndexGroupNew* __restrict saIndexGroup, const int numOfSaIndexGroup, 
			 const int firstSaIndexGrouptoProcess, int* __restrict saIndexGroupProcessed,
			 const int discardDiagonalHit,
			 char* workingMemory, const int workingMemorySize,
			 BWTSaRetrievalStatistics* __restrict bwtSaRetrievalStatistics);
//unsigned int BWTSubPatternHammingDistMatch(const unsigned char *convertedKey, const unsigned int keyLength, const unsigned int subPatternLength, const BWT *bwt, 
//								  const unsigned int maxError, const unsigned int skip, const unsigned int matchBitVector,
//								  const SaIndexRange *saIndexRange, const unsigned int saIndexRangeNumOfChar,
//								  char *workingMemory, const unsigned int workingMemorySize, const unsigned int sortOption,
//								  BWTSaRetrievalStatistics* __restrict bwtSaRetrievalStatistics);
unsigned int BWTSubPatternEditDistMatch(const unsigned char *convertedKey, const unsigned int keyLength, const unsigned int subPatternLength, const BWT *bwt, 
								  const unsigned int maxError, const unsigned int skip, const unsigned int maxnumOfHit, 
								  HitListWithPosQueryLengthError* __restrict hitList, BWTSaRetrievalStatistics* __restrict bwtSaRetrievalStatistics, 
								  const unsigned int eliminateDuplicateStartingPos);
int BWTGappedDPDBvsQuery(BWT *bwt, const unsigned char *convertedKey, const int queryPatternLength, 
						 char* __restrict workingMemory, const int workingMemorySize, int* __restrict totalNumOfQueryPos,
						 const int matchScore, const int mismatchScore,
						 const int gapOpenScore, const int gapExtendScore,
						 const int cutoffScore,
						 BWTDPStatistics* __restrict bwtDPStatistics,
						 const int printProgressDepth);

// Text retrieval functions
// Position in text will be placed at the first word of hitListSizeInWord

// startSaIndex + resultInfo must be sorted in increasing order; there must be no overlapping groups except that one group can completely enclose another
unsigned int BWTTextPosition(const BWT *bwt, const SaIndexGroupNew *saIndexGroup, const unsigned int numOfSaIndexGroups, 
				    HitList* __restrict hitList, 
					SaIndexList* __restrict tempSaIndexList1, SaIndexList* __restrict tempSaIndexList2, 
					BWTSaRetrievalStatistics* __restrict bwtSaRetrievalStatistics, const unsigned int eliminateDuplicateStartingPos);

unsigned int BWTDecodeTextPosition(const BWT *bwt, SaIndexList* __restrict evenIterationSaIndexList, SaIndexList* __restrict oddIterationSaIndexList,
						  const unsigned int numOfSaIndex, const unsigned int numOfIterationProcessed, const unsigned int maxnumOfIteration, 
						  HitList* __restrict hitList, BWTSaRetrievalStatistics* __restrict bwtSaRetrievalStatistics);



unsigned int BWTDecompressText(const BWT *bwt, const unsigned int endSaIndex, const unsigned int length, const unsigned char *reverseCharMap,
					  unsigned char *decompressedText);
unsigned int BWTDecompressTextAsWordPacked(const BWT *bwt, const unsigned int endSaIndex, const unsigned int length, unsigned int *decompressedText);

void BWTPrintDPStatistics(FILE * outFile, const BWTDPStatistics* bwtDPStatistics);
void BWTInitializeSaRetrievalStatistics(BWTSaRetrievalStatistics *bwtSaRetrievalStatistics);
void BWTAllocateDPStatistics(BWTDPStatistics *bwtDPStatistics);
void BWTInitializeDPStatistics(BWTDPStatistics *bwtDPStatistics);
void BWTFreeDPStatistics(BWTDPStatistics *bwtDPStatistics);

// QSort comparison functions
int SaIndexGroupStartSaIndexOrder(const void *saIndexGroup, const int index1, const int index2);
int SaIndexGroupStartSaIndexLengthErrorOrder(const void *saIndexGroup, const int index1, const int index2);
int HitListPosTextErrorLengthOrder(const void *hitList, const int index1, const int index2);
int HitListPosText16BitOrder(const void *hitList, const int index1, const int index2);
int HitListPosTextOrder(const void *hitList, const int index1, const int index2);
int GappedHitListScorePosTextOrder(const void *gappedHitList, const int index1, const int index2);
int GappedHitListDbSeqIndexScorePosTextOrder(const void *gappedHitList, const int index1, const int index2);


#endif
#else
/*

   BWT.h	BWT-Index

   This module contains an implementation of BWT-index for alphabet size = 4.
   The functions provided include:
    Load functions for loading BWT to memory;
    Core functions for accessing core Inverse Psi values;
	Search functions for searching patterns from text;
	Text retrieval functions for retrieving text from BWT.

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

#ifndef __BWT_H__
#define __BWT_H__

#include "TypeNLimit-l.h"
#include "MemManager-l.h"
#include "TextConverter-l.h"
#include "HSP-l.h"

#define BITS_PER_OCC_VALUE			16
#define OCC_VALUE_PER_WORD			2
#define OCC_INTERVAL				256
#define OCC_INTERVAL_MAJOR			65536

#define SORT_ALL					0
#define SORT_16_BIT					1
#define SORT_NONE					2

#define BUCKET_BIT					16
#define NUM_BUCKET					65536

#define MAX_APPROX_MATCH_ERROR	7
#define MAX_ARPROX_MATCH_LENGTH	32

#define BWTDP_MAX_SUBSTRING_LENGTH	512

typedef struct SaIndexRange {
	unsigned int startSaIndex;
	unsigned int endSaIndex;
} SaIndexRange;


typedef struct BWT {
	unsigned int textLength;			// length of the text
	unsigned int saInterval;			// interval between two SA values stored explicitly
	unsigned int inverseSaInterval;		// interval between two inverse SA stored explicitly
	unsigned int inverseSa0;			// SA-1[0]
	unsigned int *cumulativeFreq;		// cumulative frequency
	unsigned int *bwtCode;				// BWT code
	unsigned int *occValue;				// Occurrence values stored explicitly
	unsigned int *occValueMajor;		// Occurrence values stored explicitly
	unsigned int *saValue;				// SA values stored explicitly
	unsigned int *inverseSa;			// Inverse SA stored explicitly
	SaIndexRange *saIndexRange;			// SA index range
	int saIndexRangeNumOfChar;			// Number of characters indexed in SA index range
	unsigned int *saValueOnBoundary;	// Pre-calculated frequently referred data
	unsigned int *decodeTable;			// For decoding BWT by table lookup
	unsigned int decodeTableGenerated;	// == TRUE if decode table is generated on load and will be freed
	unsigned int bwtSizeInWord;			// Temporary variable to hold the memory allocated
	unsigned int occSizeInWord;			// Temporary variable to hold the memory allocated
	unsigned int occMajorSizeInWord;	// Temporary variable to hold the memory allocated
	unsigned int saValueSize;			// Temporary variable to hold the memory allocated
	unsigned int inverseSaSize;			// Temporary variable to hold the memory allocated
	unsigned int saIndexRangeSize;		// Temporary variable to hold the memory allocated
} BWT;

#define MAX_DIAGONAL_LEVEL 4				// Number of sub-pattern to keep for detecting diagonal hit

// Error information is stored as:
// 1. bitVector
//	  After hamming distance match
// 2. count
//    After edit distance match
// 3. score
//    After the hits are processed with scoring functions

typedef struct SaIndexGroupNew {	// SA index range and information of a particular error arrangement of a matched sub-pattern
	unsigned int startSaIndex;			// starting SA index
	unsigned int numOfMatch;				// number of match
	unsigned int posQuery;				// position in query; used for detecting diagonal hits
	unsigned int info;					// extra hit information; to be copied to hitList.info
} SaIndexGroupNew;

typedef struct SaIndexGroupOld {	// SA index range and information of a particular error arrangement of a matched sub-pattern
	unsigned int startSaIndex;			// starting SA index
	unsigned int numOfMatch;				// number of match
	unsigned int info;					// extra hit information; to be copied to hitList.info
} SaIndexGroupOld;

typedef struct SaIndexGroup {	// SA index range and information of a particular error arrangement of a matched sub-pattern
	unsigned int startSaIndex;			// starting SA index
	unsigned int numOfMatch;			// number of match
	unsigned int info;					// extra hit information
} SaIndexGroup;

typedef struct SaIndexGroupWithErrorBitVector {	// SA index range and information of a particular error arrangement of a matched sub-pattern
	unsigned int startSaIndex;			// starting SA index
	unsigned int numOfMatch;			// number of match
	unsigned int errorBitVector;			// error bit vector
} SaIndexGroupWithErrorBitVector;

typedef struct SaIndexGroupWithLengthError {	// SA index range and information of a particular error arrangement of a matched sub-pattern
	unsigned int startSaIndex;			// starting SA index
	unsigned int numOfMatch;			// number of match
	unsigned posQuery : 16;		// position in query
	unsigned length   : 8;		// length of hit
	unsigned error    : 8;		// error in hit
} SaIndexGroupWithLengthError;

typedef struct SaIndexGroupProcessed {	// Alternative usage of SaIndexGroup - once processed, error bit vector is replaced by index to text position
	unsigned int startSaIndex;			// starting SA index
	unsigned int numOfMatch;			// number of match
	unsigned int textPositionIndex;		// storing the pointer to text position
} SaIndexGroupProcessed;

typedef struct DupSaIndexGroup {	// Alternative usage of SaIndexGroup - the group duplicates another group
	unsigned int lastDupSaIndexGroupIndex;	// index to last duplicated group
	unsigned int saIndexGroupIndex;			// index to the first SA into group among the duplicates
	unsigned int textPositionIndex;			// storing the pointer to text position
} DupSaIndexGroup;

typedef struct SaIndexGroupHash {	// Hash table for checking duplicate SA index group
	unsigned int startSaIndex;
	unsigned int saIndexGroupIndex;
} SaIndexGroupHash;

typedef struct BWTSaRetrievalStatistics {
	unsigned int bwtSaRetrieved;
	unsigned int saDiagonalLinked;
	unsigned int saDiagonalFiltered;
	unsigned int saDuplicated;
} BWTSaRetrievalStatistics;

typedef struct BWTDPStatistics {
	int maxDepth;
	int maxDPCell;
	int maxDPMemoryInWord;
	int totalMaxDepth;
	int totalMaxDPCell;
	int totalMaxDPMemoryInWord;
	long long acceptedPathDepth;
	long long acceptedPath;
	long long rejectedPathDepth;
	long long rejectedPath;
	long long* __restrict totalNode;
	long long* __restrict rejectedNode;
	long long* __restrict totalDPCell;
} BWTDPStatistics;

typedef struct SaIndexList {
	unsigned int saIndex;
	unsigned int textPositionIndex;
} SaIndexList;

typedef struct HitCombination {
	int numOfCombination;
	int maxError;
	int keyLength;
	int skipTableWidth;
	int *errorPos;
	int *skip;
	int *skipErrorIndex;
} HitCombination;

typedef struct DPText {
	int charBeingProcessed;
	unsigned int saIndexLeft[ALPHABET_SIZE];
	unsigned int saIndexRight[ALPHABET_SIZE];
	int dpCellIndex;
	int numOfDpCellSegment;
} DPText;

typedef struct DPScanDepth {
	unsigned P				:	31;
	unsigned withAmbiguity	:	1;
} DPScanDepth;


// Load / unload functions
BWT *BWTCreate(MMPool *mmPool, const unsigned int textLength, unsigned int *decodeTable);
BWT *BWTLoad(MMPool *mmPool, const char *bwtCodeFileName, const char *occValueFileName, 
			 const char *saValueFileName, const char *inverseSaFileName, const char *saIndexRangeFileName,
			 unsigned int *decodeTable);
void BWTFree(MMPool *mmPool, BWT *bwt);
void BWTPrintMemoryUsage(const BWT *bwt, FILE *output, const unsigned int packedDNASize);

// Precalculate frequenctly accessed data
void BWTGenerateSaValueOnBoundary(MMPool *mmPool, BWT *bwt);

// Core functions
// The following must be customized for differenet compression schemes ***
unsigned int BWTOccValue(const BWT *bwt, unsigned int index, const unsigned int character);
void BWTOccValueTwoIndex(const BWT *bwt, unsigned int index1, unsigned int index2, const unsigned int character, unsigned int* __restrict occValue);
void BWTAllOccValue(const BWT *bwt, unsigned int index, unsigned int* __restrict occValue);
void BWTAllOccValueTwoIndex(const BWT *bwt, unsigned int index1, unsigned int index2, unsigned int* __restrict occValue1, unsigned int* __restrict occValue2);
unsigned int BWTOccValueOnSpot(const BWT *bwt, unsigned int index, unsigned int* __restrict character);
unsigned int BWTSearchOccValue(const BWT *bwt, const unsigned int character, const unsigned int searchOccValue);


// Utility functions for no compression only
unsigned int BWTResidentSizeInWord(const unsigned int numChar);
unsigned int BWTFileSizeInWord(const unsigned int numChar);
void BWTClearTrailingBwtCode(BWT *bwt);

// These are generic to different compression schemes (and generic to no compression as well)
unsigned int BWTPsiMinusValue(const BWT *bwt, const unsigned int index);
unsigned int BWTPsiPlusValue(const BWT *bwt, const unsigned int index);
unsigned int BWTSaValue(const BWT *bwt, unsigned int index);
unsigned int BWTInverseSa(const BWT *bwt, unsigned int saValue);
unsigned int BWTOccIntervalMajor(const unsigned int occInterval);
unsigned int BWTOccValueMinorSizeInWord(const unsigned int numChar);
unsigned int BWTOccValueMajorSizeInWord(const unsigned int numChar);

// Search functions
// packedText should be allocated with at least 1 Word buffer initialized to zero
int BWTForwardSearch(const unsigned int *packedKey, const unsigned int keyLength, const BWT *bwt, const unsigned int *packedText);
int BWTForwardSearchSaIndex(const unsigned int *packedKey, const unsigned int keyLength, const BWT *bwt, const unsigned int *packedText, 
					 unsigned int *resultSaIndexLeft, unsigned int *resultSaIndexRight);
int BWTForwardSearchNoText(const unsigned char *convertedKey, const unsigned int keyLength, const BWT *bwt);
int BWTForwardSearchSaIndexNoText(const unsigned char *convertedKey, const unsigned int keyLength, const BWT *bwt, 
					 unsigned int *resultSaIndexLeft, unsigned int *resultSaIndexRight);
int BWTBackwardSearch(const unsigned char *convertedKey, const unsigned int keyLength, const BWT *bwt, 
					  unsigned int *resultSaIndexLeft, unsigned int *resultSaIndexRight);
int BWTBackwardSearchCheckWithText(const unsigned char *convertedKey, const unsigned int *packedKey, const unsigned int keyLength,
							  const BWT *bwt, const unsigned int *packedText, const unsigned int textCheckingCostFactor,
							  const unsigned int maxnumOfTextPosition, HitList* __restrict hitList,
							  unsigned int *resultSaIndexLeft, unsigned int *resultSaIndexRight);

// Approximate match functions - brute force deep first search by backward search is used
unsigned int BWTHammingDistMaxSaIndexGroup(const unsigned int keyLength, const unsigned int maxError);
unsigned int BWTHammingDistCountOcc(const unsigned char *convertedKey, const unsigned int keyLength, const BWT *bwt, const unsigned int maxError, const unsigned int matchBitVector);
unsigned int BWTHammingDistMatch(const BWT *bwt, const unsigned char *convertedKey, const HitCombination *hitCombination,
						const SaIndexRange *saIndexRange, const unsigned int saIndexRangeNumOfChar,
						SaIndexGroupNew* __restrict saIndexGroup, const unsigned int maxSaIndexGroup);
int BWTHammingDistMatchOld(const unsigned char *convertedKey, const int keyLength, const BWT *bwt, const int maxError,
						   SaIndexGroupNew* __restrict saIndexGroup, const int maxSaIndexGroup,
						   const unsigned int posQuery, const unsigned int info);
unsigned int BWTEditDistMaxSaIndexGroup(const unsigned int keyLength, const unsigned int maxError);
// Does not insert characters on pattern boundary
unsigned int BWTEditDistMatch(const unsigned char *convertedKey, const unsigned int keyLength, const BWT *bwt, const unsigned int maxError,
					 const SaIndexRange *saIndexRange, const unsigned int saIndexRangeNumOfChar,
					 SaIndexGroupWithLengthError* __restrict saIndexGroup, const unsigned int maxSaIndexGroup);
unsigned int BWTEditDistMatchOld(const unsigned char *convertedKey, const unsigned int keyLength, const BWT *bwt, const unsigned int maxError,
					 SaIndexGroupWithLengthError* __restrict saIndexGroup, const unsigned int maxSaIndexGroup);


unsigned int BWTEliminateDupSaIndexGroup(SaIndexGroupWithLengthError* __restrict saIndexGroup, const unsigned int numOfSaGroup);

unsigned int BWTSubPatternHammingDistCountOcc(const unsigned char *convertedKey, const unsigned int keyLength, const unsigned int subPatternLength, const BWT *bwt, 
								     const unsigned int maxError, const unsigned int skip, const unsigned int matchBitVector);
int BWTSubPatternHammingDistSaIndex(const BWT *bwt, const unsigned char *convertedKey, const int keyLength, const int skip, 
									const HitCombination *hitCombination, 
									const SaIndexRange *saIndexRange, const unsigned int saIndexRangeNumOfChar,
									SaIndexGroupNew* __restrict saIndexGroup, const int maxnumOfSaIndexGroup,
									int* __restrict firstSaIndexGroupForSubPattern);
int BWTSubPatternHammingDistSaIndexOld(const unsigned char *convertedKey, const int keyLength, const int subPatternLength, const BWT *bwt, 
								    const int maxError, const int skip, const int lengthProcessed, int* __restrict lengthInCurrentRound,
									SaIndexGroupNew* __restrict saIndexGroup, const int maxnumOfSaIndexGroup);
int BWTDPHit(const BWT *bwt, SaIndexGroupNew* __restrict saIndexGroup, const int numOfSaIndexGroup, 
			 const int firstSaIndexGrouptoProcess, int* __restrict saIndexGroupProcessed,
			 const int discardDiagonalHit,
			 char* workingMemory, const int workingMemorySize,
			 BWTSaRetrievalStatistics* __restrict bwtSaRetrievalStatistics);
//unsigned int BWTSubPatternHammingDistMatch(const unsigned char *convertedKey, const unsigned int keyLength, const unsigned int subPatternLength, const BWT *bwt, 
//								  const unsigned int maxError, const unsigned int skip, const unsigned int matchBitVector,
//								  const SaIndexRange *saIndexRange, const unsigned int saIndexRangeNumOfChar,
//								  char *workingMemory, const unsigned int workingMemorySize, const unsigned int sortOption,
//								  BWTSaRetrievalStatistics* __restrict bwtSaRetrievalStatistics);
unsigned int BWTSubPatternEditDistMatch(const unsigned char *convertedKey, const unsigned int keyLength, const unsigned int subPatternLength, const BWT *bwt, 
								  const unsigned int maxError, const unsigned int skip, const unsigned int maxnumOfHit, 
								  const SaIndexRange *saIndexRange, const unsigned int saIndexRangeNumOfChar,
								  HitListWithPosQueryLengthError* __restrict hitList, BWTSaRetrievalStatistics* __restrict bwtSaRetrievalStatistics, 
								  const unsigned int eliminateDuplicateStartingPos);
int BWTGappedDPDBvsQuery(BWT *bwt, const unsigned char *convertedKey, const int queryPatternLength, 
						 char* __restrict workingMemory, const int workingMemorySize, int* __restrict totalNumOfQueryPos,
						 const int matchScore, const int mismatchScore,
						 const int gapOpenScore, const int gapExtendScore,
						 const int cutoffScore,
						 BWTDPStatistics* __restrict bwtDPStatistics,
						 const int printProgressDepth);

// Text retrieval functions
// Position in text will be placed at the first word of hitListSizeInWord

// startSaIndex + resultInfo must be sorted in increasing order; there must be no overlapping groups except that one group can completely enclose another
unsigned int BWTTextPosition(const BWT *bwt, const SaIndexGroupNew *saIndexGroup, const unsigned int numOfSaIndexGroups, 
				    HitList* __restrict hitList, 
					SaIndexList* __restrict tempSaIndexList1, SaIndexList* __restrict tempSaIndexList2, 
					BWTSaRetrievalStatistics* __restrict bwtSaRetrievalStatistics, const unsigned int eliminateDuplicateStartingPos);

unsigned int BWTDecodeTextPosition(const BWT *bwt, SaIndexList* __restrict evenIterationSaIndexList, SaIndexList* __restrict oddIterationSaIndexList,
						  const unsigned int numOfSaIndex, const unsigned int numOfIterationProcessed, const unsigned int maxnumOfIteration, 
						  HitList* __restrict hitList, BWTSaRetrievalStatistics* __restrict bwtSaRetrievalStatistics);



unsigned int BWTDecompressText(const BWT *bwt, const unsigned int endSaIndex, const unsigned int length, const unsigned char *reverseCharMap,
					  unsigned char *decompressedText);
unsigned int BWTDecompressTextAsWordPacked(const BWT *bwt, const unsigned int endSaIndex, const unsigned int length, unsigned int *decompressedText);

void BWTPrintDPStatistics(FILE * outFile, const BWTDPStatistics* bwtDPStatistics);
void BWTInitializeSaRetrievalStatistics(BWTSaRetrievalStatistics *bwtSaRetrievalStatistics);
void BWTAllocateDPStatistics(BWTDPStatistics *bwtDPStatistics);
void BWTInitializeDPStatistics(BWTDPStatistics *bwtDPStatistics);
void BWTFreeDPStatistics(BWTDPStatistics *bwtDPStatistics);

// QSort comparison functions
int SaIndexGroupStartSaIndexOrder(const void *saIndexGroup, const int index1, const int index2);
int SaIndexGroupStartSaIndexLengthErrorOrder(const void *saIndexGroup, const int index1, const int index2);
int HitListPosTextErrorLengthOrder(const void *hitList, const int index1, const int index2);
int HitListPosText16BitOrder(const void *hitList, const int index1, const int index2);
int HitListPosTextOrder(const void *hitList, const int index1, const int index2);
int GappedHitListScorePosTextOrder(const void *gappedHitList, const int index1, const int index2);
int GappedHitListDbSeqIndexScorePosTextOrder(const void *gappedHitList, const int index1, const int index2);


#endif
#endif
