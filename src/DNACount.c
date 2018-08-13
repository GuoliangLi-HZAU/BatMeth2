/*

   DNACount.c		DNA Count

   This module contains DNA occurrence counting functions. The DNA must be
   in word-packed format.

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
#include "DNACount.h"
#include "MiscUtilities.h"


void GenerateDNAOccCountTable(unsigned int *dnaDecodeTable) {

	unsigned int i, j, c, t;

	for (i=0; i<DNA_OCC_CNT_TABLE_SIZE_IN_WORD; i++) {
		dnaDecodeTable[i] = 0;
		c = i;
		for (j=0; j<8; j++) {
			t = c & 0x00000003;
			dnaDecodeTable[i] += 1 << (t * 8);
			c >>= 2;
		}
	}

}

unsigned int ForwardDNAOccCount(const unsigned int*  dna, const unsigned int index, const unsigned int character, const unsigned int*  dnaDecodeTable) {

	static const unsigned int truncateRightMask[16] = { 0x00000000, 0xC0000000, 0xF0000000, 0xFC000000,
											   0xFF000000, 0xFFC00000, 0xFFF00000, 0xFFFC0000,
											   0xFFFF0000, 0xFFFFC000, 0xFFFFF000, 0xFFFFFC00,
											   0xFFFFFF00, 0xFFFFFFC0, 0xFFFFFFF0, 0xFFFFFFFC };

	unsigned int wordToCount, charToCount;
	unsigned int i, c;
	unsigned int sum = 0;

#ifdef DEBUG
	if (index >= 256) {
		fprintf(stderr, "ForwardDNAOccCount() : index >= 256!\n");
		exit(1);
	}
#endif

	wordToCount = index / 16;
	charToCount = index - wordToCount * 16;

	for (i=0; i<wordToCount; i++) {
		sum += dnaDecodeTable[dna[i] >> 16];
		sum += dnaDecodeTable[dna[i] & 0x0000FFFF];
	}

	if (charToCount > 0) {
		c = dna[i] & truncateRightMask[charToCount];	// increase count of 'a' by 16 - c;
		sum += dnaDecodeTable[c >> 16];
		sum += dnaDecodeTable[c & 0xFFFF];
		sum += charToCount - 16;	// decrease count of 'a' by 16 - positionToProcess
	}

	return (sum >> (character * 8)) & 0x000000FF;

}

unsigned int BackwardDNAOccCount(const unsigned int*  dna, const unsigned int index, const unsigned int character, const unsigned int*  dnaDecodeTable) {

	static const unsigned int truncateLeftMask[16] =  { 0x00000000, 0x00000003, 0x0000000F, 0x0000003F,
											   0x000000FF, 0x000003FF, 0x00000FFF, 0x00003FFF,
											   0x0000FFFF, 0x0003FFFF, 0x000FFFFF, 0x003FFFFF,
											   0x00FFFFFF, 0x03FFFFFF, 0x0FFFFFFF, 0x3FFFFFFF };

	unsigned int wordToCount, charToCount;
	unsigned int i, c;
	unsigned int sum = 0;

#ifdef DEBUG
	if (index >= 256) {
		fprintf(stderr, "ForwardDNAOccCount() : index >= 256!\n");
		exit(1);
	}
#endif

	wordToCount = index / 16;
	charToCount = index - wordToCount * 16;

	dna -= wordToCount + 1;

	if (charToCount > 0) {
		c = *dna & truncateLeftMask[charToCount];	// increase count of 'a' by 16 - c;
		sum += dnaDecodeTable[c >> 16];
		sum += dnaDecodeTable[c & 0xFFFF];
		sum += charToCount - 16;	// decrease count of 'a' by 16 - positionToProcess
	}
	
	for (i=0; i<wordToCount; i++) {
		dna++;
		sum += dnaDecodeTable[*dna >> 16];
		sum += dnaDecodeTable[*dna & 0x0000FFFF];
	}

	return (sum >> (character * 8)) & 0x000000FF;

}

void ForwardDNAAllOccCount(const unsigned int*  dna, const unsigned int index, unsigned int*  __restrict occCount, const unsigned int*  dnaDecodeTable) {

	static const unsigned int truncateRightMask[16] = { 0x00000000, 0xC0000000, 0xF0000000, 0xFC000000,
											   0xFF000000, 0xFFC00000, 0xFFF00000, 0xFFFC0000,
											   0xFFFF0000, 0xFFFFC000, 0xFFFFF000, 0xFFFFFC00,
											   0xFFFFFF00, 0xFFFFFFC0, 0xFFFFFFF0, 0xFFFFFFFC };

	unsigned int wordToCount, charToCount;
	unsigned int i, c;
	unsigned int sum = 0;

#ifdef DEBUG
	if (index >= 256) {
		fprintf(stderr, "ForwardDNAOccCount() : index >= 256!\n");
		exit(1);
	}
#endif

	wordToCount = index / 16;
	charToCount = index - wordToCount * 16;

	for (i=0; i<wordToCount; i++) {
		sum += dnaDecodeTable[dna[i] >> 16];
		sum += dnaDecodeTable[dna[i] & 0x0000FFFF];
	}

	if (charToCount > 0) {
		c = dna[i] & truncateRightMask[charToCount];	// increase count of 'a' by 16 - c;
		sum += dnaDecodeTable[c >> 16];
		sum += dnaDecodeTable[c & 0xFFFF];
		sum += charToCount - 16;	// decrease count of 'a' by 16 - positionToProcess
	}

	occCount[0] = sum & 0x000000FF;	sum >>= 8;
	occCount[1] = sum & 0x000000FF;	sum >>= 8;
	occCount[2] = sum & 0x000000FF;	sum >>= 8;
	occCount[3] = sum;

}

void BackwardDNAAllOccCount(const unsigned int*  dna, const unsigned int index, unsigned int* __restrict occCount, const unsigned int*  dnaDecodeTable) {

	static const unsigned int truncateLeftMask[16] =  { 0x00000000, 0x00000003, 0x0000000F, 0x0000003F,
											   0x000000FF, 0x000003FF, 0x00000FFF, 0x00003FFF,
											   0x0000FFFF, 0x0003FFFF, 0x000FFFFF, 0x003FFFFF,
											   0x00FFFFFF, 0x03FFFFFF, 0x0FFFFFFF, 0x3FFFFFFF };

	unsigned int wordToCount, charToCount;
	unsigned int i, c;
	unsigned int sum = 0;

#ifdef DEBUG
	if (index >= 256) {
		fprintf(stderr, "ForwardDNAOccCount() : index >= 256!\n");
		exit(1);
	}
#endif

	wordToCount = index / 16;
	charToCount = index - wordToCount * 16;

	dna -= wordToCount + 1;

	if (charToCount > 0) {
		c = *dna & truncateLeftMask[charToCount];	// increase count of 'a' by 16 - c;
		sum += dnaDecodeTable[c >> 16];
		sum += dnaDecodeTable[c & 0xFFFF];
		sum += charToCount - 16;	// decrease count of 'a' by 16 - positionToProcess
	}

	for (i=0; i<wordToCount; i++) {
		dna++;
		sum += dnaDecodeTable[*dna >> 16];
		sum += dnaDecodeTable[*dna & 0x0000FFFF];
	}
	
	occCount[0] = sum & 0x000000FF;	sum >>= 8;
	occCount[1] = sum & 0x000000FF;	sum >>= 8;
	occCount[2] = sum & 0x000000FF;	sum >>= 8;
	occCount[3] = sum;

}

unsigned int Forward1OccCount(const unsigned int*  bitVector, const unsigned int index, const unsigned int*  dnaDecodeTable) {

	static const unsigned int truncateRightMask[32] = { 0x00000000, 0x80000000, 0xC0000000, 0xE0000000, 0xF0000000, 0xF8000000, 0xFC000000, 0xFE000000,
											   0xFF000000, 0xFF800000, 0xFFC00000, 0xFFE00000, 0xFFF00000, 0xFFF80000, 0xFFFC0000, 0xFFFE0000,
											   0xFFFF0000, 0xFFFF8000, 0xFFFFC000, 0xFFFFE000, 0xFFFFF000, 0xFFFFF800, 0xFFFFFC00, 0xFFFFFE00,
											   0xFFFFFF00, 0xFFFFFF80, 0xFFFFFFC0, 0xFFFFFFE0, 0xFFFFFFF0, 0xFFFFFFF8, 0xFFFFFFFC, 0xFFFFFFFE};

	unsigned int wordToCount, bitToCount;
	unsigned int i, c;
	unsigned int sum = 0;
	unsigned int numberOf1;

#ifdef DEBUG
	if (index >= 256) {
		fprintf(stderr, "Forward1OccCount() : index >= 256!\n");
		exit(1);
	}
#endif

	wordToCount = index / 32;
	bitToCount = index - wordToCount * 32;

	for (i=0; i<wordToCount; i++) {
		sum += dnaDecodeTable[bitVector[i] >> 16];
		sum += dnaDecodeTable[bitVector[i] & 0x0000FFFF];
	}

	if (bitToCount > 0) {
		c = bitVector[i] & truncateRightMask[bitToCount];
		sum += dnaDecodeTable[c >> 16];
		sum += dnaDecodeTable[c & 0x0000FFFF];
	}

	sum >>= 8;
	numberOf1 = sum & 0x000000FF;
	sum >>= 8;
	numberOf1 = sum & 0x000000FF;
	sum >>= 8;
	numberOf1 = sum * 2;

	return numberOf1;

}

unsigned int Backward1OccCount(const unsigned int*  bitVector, const unsigned int index, const unsigned int*  dnaDecodeTable) {

	static const unsigned int truncateLeftMask[32] =  { 0x00000000, 0x00000001, 0x00000003, 0x00000007, 0x0000000F, 0x0000001F, 0x0000003F, 0x0000007F,
											   0x000000FF, 0x000001FF, 0x000003FF, 0x000007FF, 0x00000FFF, 0x00001FFF, 0x00003FFF, 0x00007FFF,
											   0x0000FFFF, 0x0001FFFF, 0x0003FFFF, 0x0007FFFF, 0x000FFFFF, 0x001FFFFF, 0x003FFFFF, 0x007FFFFF,
											   0x00FFFFFF, 0x01FFFFFF, 0x03FFFFFF, 0x07FFFFFF, 0x0FFFFFFF, 0x1FFFFFFF, 0x3FFFFFFF, 0x7FFFFFFF};

	unsigned int wordToCount, bitToCount;
	unsigned int i, c;
	unsigned int sum = 0;
	unsigned int numberOf1;

#ifdef DEBUG
	if (index >= 256) {
		fprintf(stderr, "ForwardDNAOccCount() : index >= 256!\n");
		exit(1);
	}
#endif

	wordToCount = index / 32;
	bitToCount = index - wordToCount * 32;

	bitVector -= wordToCount + 1;

	if (bitToCount > 0) {
		c = *bitVector & truncateLeftMask[bitToCount];
		sum += dnaDecodeTable[c >> 16];
		sum += dnaDecodeTable[c & 0xFFFF];
	}
	
	for (i=0; i<wordToCount; i++) {
		bitVector++;
		sum += dnaDecodeTable[*bitVector >> 16];
		sum += dnaDecodeTable[*bitVector & 0x0000FFFF];
	}

	sum >>= 8;
	numberOf1 = sum & 0x000000FF;
	sum >>= 8;
	numberOf1 = sum & 0x000000FF;
	sum >>= 8;
	numberOf1 = sum * 2;

	return numberOf1;

}

unsigned int ForwardDNAOccCountNoLimit(const unsigned int*  dna, const unsigned int index, const unsigned int character, const unsigned int*  dnaDecodeTable) {

	static const unsigned int truncateRightMask[16] = { 0x00000000, 0xC0000000, 0xF0000000, 0xFC000000,
											   0xFF000000, 0xFFC00000, 0xFFF00000, 0xFFFC0000,
											   0xFFFF0000, 0xFFFFC000, 0xFFFFF000, 0xFFFFFC00,
											   0xFFFFFF00, 0xFFFFFFC0, 0xFFFFFFF0, 0xFFFFFFFC };

	unsigned int iteration, wordToCount, charToCount;
	unsigned int i, j, c;
	unsigned int sum;
	unsigned int occCount = 0;

	iteration = index / 256;
	wordToCount = (index - iteration * 256) / 16;
	charToCount = index - iteration * 256 - wordToCount * 16;

	for (i=0; i<iteration; i++) {

		sum = 0;
		for (j=0; j<16; j++) {
			sum += dnaDecodeTable[*dna >> 16];
			sum += dnaDecodeTable[*dna & 0x0000FFFF];
			dna++;
		}
		if (!DNA_OCC_SUM_EXCEPTION(sum)) {
            occCount += (sum >> (character * 8)) & 0x000000FF;
		} else {
			// only some or all of the 3 bits are on
			// in reality, only one of the four cases are possible
			if (sum == 0x00000100) {
				if (character == 0) {
					occCount += 256;
				}
			} else if (sum == 0x00010000) {
				if (character == 1) {
					occCount += 256;
				}
			} else if (sum == 0x01000000) {
				if (character == 2) {
					occCount += 256;
				}
			} else if (sum == 0x00000000) {
				if (character == 3) {
					occCount += 256;
				}
			} else {
				fprintf(stderr, "ForwardDNAOccCountNoLimit(): DNA occ sum exception!\n");
				exit(1);
			}
		}

	}

	sum = 0;
	for (j=0; j<wordToCount; j++) {
		sum += dnaDecodeTable[*dna >> 16];
		sum += dnaDecodeTable[*dna & 0x0000FFFF];
		dna++;
	}

	if (charToCount > 0) {
		c = *dna & truncateRightMask[charToCount];	// increase count of 'a' by 16 - c;
		sum += dnaDecodeTable[c >> 16];
		sum += dnaDecodeTable[c & 0xFFFF];
		sum += charToCount - 16;	// decrease count of 'a' by 16 - positionToProcess
	}

	occCount += (sum >> (character * 8)) & 0x000000FF;

	return occCount;

}

unsigned int BackwardDNAOccCountNoLimit(const unsigned int*  dna, const unsigned int index, const unsigned int character, const unsigned int*  dnaDecodeTable) {

	static const unsigned int truncateLeftMask[16] =  { 0x00000000, 0x00000003, 0x0000000F, 0x0000003F,
											   0x000000FF, 0x000003FF, 0x00000FFF, 0x00003FFF,
											   0x0000FFFF, 0x0003FFFF, 0x000FFFFF, 0x003FFFFF,
											   0x00FFFFFF, 0x03FFFFFF, 0x0FFFFFFF, 0x3FFFFFFF };

	unsigned int iteration, wordToCount, charToCount;
	unsigned int i, j, c;
	unsigned int sum = 0;
	unsigned int occCount;

	dna -= index / 16 + 1;

	iteration = index / 256;
	wordToCount = (index - iteration * 256) / 16;
	charToCount = index - iteration * 256 - wordToCount * 16;

	if (charToCount > 0) {
		c = *dna & truncateLeftMask[charToCount];	// increase count of 'a' by 16 - c;
		sum += dnaDecodeTable[c >> 16];
		sum += dnaDecodeTable[c & 0xFFFF];
		sum += charToCount - 16;	// decrease count of 'a' by 16 - positionToProcess
	}

	for (j=0; j<wordToCount; j++) {
		dna++;
		sum += dnaDecodeTable[*dna >> 16];
		sum += dnaDecodeTable[*dna & 0x0000FFFF];
	}

	occCount = (sum >> (character * 8)) & 0x000000FF;


	for (i=0; i<iteration; i++) {

		sum = 0;
		for (j=0; j<16; j++) {
			dna++;
			sum += dnaDecodeTable[*dna >> 16];
			sum += dnaDecodeTable[*dna & 0x0000FFFF];
		}
		if (!DNA_OCC_SUM_EXCEPTION(sum)) {
            occCount += (sum >> (character * 8)) & 0x000000FF;
		} else {
			// only some or all of the 3 bits are on
			// in reality, only one of the four cases are possible
			if (sum == 0x00000100) {
				if (character == 0) {
					occCount += 256;
				}
			} else if (sum == 0x00010000) {
				if (character == 1) {
					occCount += 256;
				}
			} else if (sum == 0x01000000) {
				if (character == 2) {
					occCount += 256;
				}
			} else if (sum == 0x00000000) {
				if (character == 3) {
					occCount += 256;
				}
			} else {
				fprintf(stderr, "BackwardDNAOccCountNoLimit(): DNA occ sum exception!\n");
				exit(1);
			}
		}

	}

	return occCount;

}

void ForwardDNAAllOccCountNoLimit(const unsigned int*  dna, const unsigned int index, unsigned int* __restrict occCount, const unsigned int*  dnaDecodeTable) {

	static const unsigned int truncateRightMask[16] = { 0x00000000, 0xC0000000, 0xF0000000, 0xFC000000,
											   0xFF000000, 0xFFC00000, 0xFFF00000, 0xFFFC0000,
											   0xFFFF0000, 0xFFFFC000, 0xFFFFF000, 0xFFFFFC00,
											   0xFFFFFF00, 0xFFFFFFC0, 0xFFFFFFF0, 0xFFFFFFFC };

	unsigned int iteration, wordToCount, charToCount;
	unsigned int i, j, c;
	unsigned int sum;

	occCount[0] = 0;
	occCount[1] = 0;
	occCount[2] = 0;
	occCount[3] = 0;

	iteration = index / 256;
	wordToCount = (index - iteration * 256) / 16;
	charToCount = index - iteration * 256 - wordToCount * 16;

	for (i=0; i<iteration; i++) {

		sum = 0;
		for (j=0; j<16; j++) {
			sum += dnaDecodeTable[*dna >> 16];
			sum += dnaDecodeTable[*dna & 0x0000FFFF];
			dna++;
		}
		if (!DNA_OCC_SUM_EXCEPTION(sum)) {
			occCount[0] += sum & 0x000000FF;	sum >>= 8;
			occCount[1] += sum & 0x000000FF;	sum >>= 8;
			occCount[2] += sum & 0x000000FF;	sum >>= 8;
			occCount[3] += sum;
		} else {
			// only some or all of the 3 bits are on
			// in reality, only one of the four cases are possible
			if (sum == 0x00000100) {
				occCount[0] += 256;
			} else if (sum == 0x00010000) {
				occCount[1] += 256;
			} else if (sum == 0x01000000) {
				occCount[2] += 256;
			} else if (sum == 0x00000000) {
				occCount[3] += 256;
			} else {
				fprintf(stderr, "ForwardDNAAllOccCountNoLimit(): DNA occ sum exception!\n");
				exit(1);
			}
		}

	}

	sum = 0;
	for (j=0; j<wordToCount; j++) {
		sum += dnaDecodeTable[*dna >> 16];
		sum += dnaDecodeTable[*dna & 0x0000FFFF];
		dna++;
	}

	if (charToCount > 0) {
		c = *dna & truncateRightMask[charToCount];	// increase count of 'a' by 16 - c;
		sum += dnaDecodeTable[c >> 16];
		sum += dnaDecodeTable[c & 0xFFFF];
		sum += charToCount - 16;	// decrease count of 'a' by 16 - positionToProcess
	}

	occCount[0] += sum & 0x000000FF;	sum >>= 8;
	occCount[1] += sum & 0x000000FF;	sum >>= 8;
	occCount[2] += sum & 0x000000FF;	sum >>= 8;
	occCount[3] += sum;

}

void BackwardDNAAllOccCountNoLimit(const unsigned int*  dna, const unsigned int index, unsigned int* __restrict occCount, const unsigned int*  dnaDecodeTable) {

	static const unsigned int truncateLeftMask[16] =  { 0x00000000, 0x00000003, 0x0000000F, 0x0000003F,
											   0x000000FF, 0x000003FF, 0x00000FFF, 0x00003FFF,
											   0x0000FFFF, 0x0003FFFF, 0x000FFFFF, 0x003FFFFF,
											   0x00FFFFFF, 0x03FFFFFF, 0x0FFFFFFF, 0x3FFFFFFF };

	unsigned int iteration, wordToCount, charToCount;
	unsigned int i, j, c;
	unsigned int sum;

	dna -= index / 16 + 1;

	iteration = index / 256;
	wordToCount = (index - iteration * 256) / 16;
	charToCount = index - iteration * 256 - wordToCount * 16;

	sum = 0;
	if (charToCount > 0) {
		c = *dna & truncateLeftMask[charToCount];	// increase count of 'a' by 16 - c;
		sum += dnaDecodeTable[c >> 16];
		sum += dnaDecodeTable[c & 0xFFFF];
		sum += charToCount - 16;	// decrease count of 'a' by 16 - positionToProcess
	}

	for (j=0; j<wordToCount; j++) {
		dna++;
		sum += dnaDecodeTable[*dna >> 16];
		sum += dnaDecodeTable[*dna & 0x0000FFFF];
	}

	occCount[0] = sum & 0x000000FF;	sum >>= 8;
	occCount[1] = sum & 0x000000FF;	sum >>= 8;
	occCount[2] = sum & 0x000000FF;	sum >>= 8;
	occCount[3] = sum;

	for (i=0; i<iteration; i++) {

		sum = 0;
		for (j=0; j<16; j++) {
			dna++;
			sum += dnaDecodeTable[*dna >> 16];
			sum += dnaDecodeTable[*dna & 0x0000FFFF];
		}
		if (!DNA_OCC_SUM_EXCEPTION(sum)) {
			occCount[0] += sum & 0x000000FF;	sum >>= 8;
			occCount[1] += sum & 0x000000FF;	sum >>= 8;
			occCount[2] += sum & 0x000000FF;	sum >>= 8;
			occCount[3] += sum;
		} else {
			// only some or all of the 3 bits are on
			// in reality, only one of the four cases are possible
			if (sum == 0x00000100) {
				occCount[0] += 256;
			} else if (sum == 0x00010000) {
				occCount[1] += 256;
			} else if (sum == 0x01000000) {
				occCount[2] += 256;
			} else if (sum == 0x00000000) {
				occCount[3] += 256;
			} else {
				fprintf(stderr, "BackwardDNAAllOccCountNoLimit(): DNA occ sum exception!\n");
				exit(1);
			}
		}

	}

}


void GenerateDNA_NOccCountTable(unsigned int *dnaDecodeTable) {

	unsigned int i, j, c, t;

	for (i=0; i<DNA_N_OCC_CNT_TABLE_SIZE_IN_WORD; i++) {
		dnaDecodeTable[i] = 0;
		c = i;
		for (j=0; j<5; j++) {
			t = c & 0x00000007;
			if (t != 4) {
				// Count of 'n' is to be derived from acgt
				dnaDecodeTable[i] += 1 << (t * 8);
			}
			c >>= 3;
		}
	}

}

unsigned int ForwardDNA_NOccCount(const unsigned int*  dna, const unsigned int index, const unsigned int character, const unsigned int*  dnaDecodeTable) {

	static const unsigned int truncateRightMask[10] = { 0x00000000, 0xE0000000, 0xFC000000, 0xFF800000,
											   0xFFF00000, 0xFFFE0000, 0xFFFFC000, 0xFFFFF800,
											   0xFFFFFF00, 0xFFFFFFE0};

	unsigned int wordToCount, charToCount;
	unsigned int i, c;
	unsigned int sum = 0;
	unsigned int occCount;

#ifdef DEBUG
	if (index > 250) {
		fprintf(stderr, "ForwardDNA_NOccCount() : index > 250!\n");
		exit(1);
	}
#endif

	wordToCount = index / 10;
	charToCount = index - wordToCount * 10;

	for (i=0; i<wordToCount; i++) {
		sum += dnaDecodeTable[dna[i] >> 17];
		sum += dnaDecodeTable[(dna[i] >> 2) & 0x00007FFF];
	}

	if (charToCount > 0) {
		c = dna[i] & truncateRightMask[charToCount];	// increase count of 'a' by 10 - charToCount;
		sum += dnaDecodeTable[c >> 17];
		sum += dnaDecodeTable[(c >> 2) & 0x00007FFF];
		sum += charToCount - 10;	// decrease count of 'a' by 10 - charToCount
	}

	if (character != 4) {
		occCount = (sum >> (character * 8)) & 0x000000FF;
	} else {
		occCount = index;
		occCount -= sum & 0x000000FF;	sum >>= 8;
		occCount -= sum & 0x000000FF;	sum >>= 8;
		occCount -= sum & 0x000000FF;	sum >>= 8;
		occCount -= sum;
	}

	return occCount;

}

unsigned int BackwardDNA_NOccCount(const unsigned int*  dna, const unsigned int index, const unsigned int character, const unsigned int*  dnaDecodeTable) {

	static const unsigned int truncateLeftMask[10] =  { 0x00000000, 0x0000001C, 0x000000FC, 0x000007FC,
											   0x00003FFC, 0x0001FFFC, 0x000FFFFC, 0x007FFFFC,
											   0x03FFFFFC, 0x1FFFFFFC};

	unsigned int wordToCount, charToCount;
	unsigned int j, c;
	unsigned int sum = 0;
	unsigned int occCount;

#ifdef DEBUG
	if (index > 250) {
		fprintf(stderr, "BackwardDNA_NOccCount() : index >= 250!\n");
		exit(1);
	}
#endif

	wordToCount = index / 10;
	charToCount = index - wordToCount * 10;

	dna -= wordToCount + 1;

	if (charToCount > 0) {
		c = *dna & truncateLeftMask[charToCount];	// increase count of 'a' by 10 - charToCount;
		sum += dnaDecodeTable[c >> 17];
		sum += dnaDecodeTable[(c >> 2) & 0x00007FFF];
		sum += charToCount - 10;	// decrease count of 'a' by 10 - charToCount
	}

	for (j=0; j<wordToCount; j++) {
		dna++;
		sum += dnaDecodeTable[*dna >> 17];
		sum += dnaDecodeTable[(*dna >> 2) & 0x00007FFF];
	}

	if (character != 4) {
		occCount = (sum >> (character * 8)) & 0x000000FF;
	} else {
		occCount = index;
		occCount -= sum & 0x000000FF;	sum >>= 8;
		occCount -= sum & 0x000000FF;	sum >>= 8;
		occCount -= sum & 0x000000FF;	sum >>= 8;
		occCount -= sum;
	}

#ifdef DEBUG
	if (occCount > index + 1) {
		fprintf(stderr, "BackwardDNA_NOccCount() : occCount > index + 1!\n");
		exit(1);
	}
#endif

	return occCount;

}

void ForwardDNA_NAllOccCount(const unsigned int*  dna, const unsigned int index, unsigned int* __restrict occCount, const unsigned int*  dnaDecodeTable) {

	static const unsigned int truncateRightMask[10] = { 0x00000000, 0xE0000000, 0xFC000000, 0xFF800000,
											   0xFFF00000, 0xFFFE0000, 0xFFFFC000, 0xFFFFF800,
											   0xFFFFFF00, 0xFFFFFFE0};

	unsigned int wordToCount, charToCount;
	unsigned int i, c;
	unsigned int sum = 0;

#ifdef DEBUG
	if (index > 250) {
		fprintf(stderr, "ForwardDNA_NAllOccCount() : index >= 250!\n");
		exit(1);
	}
#endif

	wordToCount = index / 10;
	charToCount = index - wordToCount * 10;

	for (i=0; i<wordToCount; i++) {
		sum += dnaDecodeTable[dna[i] >> 17];
		sum += dnaDecodeTable[(dna[i] >> 2) & 0x00007FFF];
	}

	if (charToCount > 0) {
		c = dna[i] & truncateRightMask[charToCount];	// increase count of 'a' by 10 - charToCount;
		sum += dnaDecodeTable[c >> 17];
		sum += dnaDecodeTable[(c >> 2) & 0x00007FFF];
		sum += charToCount - 10;	// decrease count of 'a' by 10 - charToCount
	}

	occCount[0] = sum & 0x000000FF;	sum >>= 8;
	occCount[1] = sum & 0x000000FF;	sum >>= 8;
	occCount[2] = sum & 0x000000FF;	sum >>= 8;
	occCount[3] = sum;

}

void BackwardDNA_NAllOccCount(const unsigned int*  dna, const unsigned int index, unsigned int* __restrict occCount, const unsigned int*  dnaDecodeTable) {

	static const unsigned int truncateLeftMask[10] =  { 0x00000000, 0x0000001C, 0x000000FC, 0x000007FC,
											   0x00003FFC, 0x0001FFFC, 0x000FFFFC, 0x007FFFFC,
											   0x03FFFFFC, 0x1FFFFFFC};

	unsigned int wordToCount, charToCount;
	unsigned int j, c;
	unsigned int sum = 0;

#ifdef DEBUG
	if (index > 250) {
		fprintf(stderr, "BackwardDNA_NAllOccCount() : index >= 250!\n");
		exit(1);
	}
#endif

	wordToCount = index / 10;
	charToCount = index - wordToCount * 10;

	dna -= wordToCount + 1;

	if (charToCount > 0) {
		c = *dna & truncateLeftMask[charToCount];	// increase count of 'a' by 10 - charToCount;
		sum += dnaDecodeTable[c >> 17];
		sum += dnaDecodeTable[(c >> 2) & 0x00007FFF];
		sum += charToCount - 10;	// decrease count of 'a' by 16 - charToCount
	}

	for (j=0; j<wordToCount; j++) {
		dna++;
		sum += dnaDecodeTable[*dna >> 17];
		sum += dnaDecodeTable[(*dna >> 2) & 0x00007FFF];
	}

	occCount[0] = sum & 0x000000FF;	sum >>= 8;
	occCount[1] = sum & 0x000000FF;	sum >>= 8;
	occCount[2] = sum & 0x000000FF;	sum >>= 8;
	occCount[3] = sum;

}

unsigned int ForwardDNA_NOccCountNoLimit(const unsigned int*  dna, const unsigned int index, const unsigned int character, const unsigned int*  dnaDecodeTable) {

	static const unsigned int truncateRightMask[10] = { 0x00000000, 0xE0000000, 0xFC000000, 0xFF800000,
											   0xFFF00000, 0xFFFE0000, 0xFFFFC000, 0xFFFFF800,
											   0xFFFFFF00, 0xFFFFFFE0};

	unsigned int iteration, wordToCount, charToCount;
	unsigned int i, j, c;
	unsigned int sum;
	unsigned int occCount = 0;

	iteration = index / 250;
	wordToCount = (index - iteration * 250) / 10;
	charToCount = index - iteration * 250 - wordToCount * 10;

	for (i=0; i<iteration; i++) {

		sum = 0;
		for (j=0; j<25; j++) {
			sum += dnaDecodeTable[*dna >> 17];
			sum += dnaDecodeTable[(*dna >> 2) & 0x00007FFF];
			dna++;
		}
		if (character != 4) {
			occCount += (sum >> (character * 8)) & 0x000000FF;
		} else {
			occCount -= sum & 0x000000FF;	sum >>= 8;
			occCount -= sum & 0x000000FF;	sum >>= 8;
			occCount -= sum & 0x000000FF;	sum >>= 8;
			occCount -= sum;
		}
	}

	sum = 0;
	for (j=0; j<wordToCount; j++) {
		sum += dnaDecodeTable[*dna >> 17];
		sum += dnaDecodeTable[(*dna >> 2) & 0x00007FFF];
		dna++;
	}

	if (charToCount > 0) {
		c = *dna & truncateRightMask[charToCount];	// increase count of 'a' by 10 - charToCount;
		sum += dnaDecodeTable[c >> 17];
		sum += dnaDecodeTable[(c >> 2) & 0x00007FFF];
		sum += charToCount - 10;	// decrease count of 'a' by 10 - charToCount
	}

	if (character != 4) {
		occCount += (sum >> (character * 8)) & 0x000000FF;
	} else {
		occCount += index;
		occCount -= sum & 0x000000FF;	sum >>= 8;
		occCount -= sum & 0x000000FF;	sum >>= 8;
		occCount -= sum & 0x000000FF;	sum >>= 8;
		occCount -= sum;
	}

	return occCount;

}

unsigned int BackwardDNA_NOccCountNoLimit(const unsigned int*  dna, const unsigned int index, const unsigned int character, const unsigned int*  dnaDecodeTable) {

	static const unsigned int truncateLeftMask[10] =  { 0x00000000, 0x0000001C, 0x000000FC, 0x000007FC,
											   0x00003FFC, 0x0001FFFC, 0x000FFFFC, 0x007FFFFC,
											   0x03FFFFFC, 0x1FFFFFFC};

	unsigned int iteration, wordToCount, charToCount;
	unsigned int i, j, c;
	unsigned int sum = 0;
	unsigned int occCount = 0;

	dna -= index / 10 + 1;

	iteration = index / 250;
	wordToCount = (index - iteration * 250) / 10;
	charToCount = index - iteration * 250 - wordToCount * 10;

	if (charToCount > 0) {
		c = *dna & truncateLeftMask[charToCount];	// increase count of 'a' by 10 - charToCount;
		sum += dnaDecodeTable[c >> 17];
		sum += dnaDecodeTable[(c >> 2) & 0x00007FFF];
		sum += charToCount - 10;	// decrease count of 'a' by 10 - charToCount
	}

	for (j=0; j<wordToCount; j++) {
		dna++;
		sum += dnaDecodeTable[*dna >> 17];
		sum += dnaDecodeTable[(*dna >> 2) & 0x00007FFF];
	}

	if (character != 4) {
		occCount = (sum >> (character * 8)) & 0x000000FF;
	} else {
		occCount = index;
		occCount -= sum & 0x000000FF;	sum >>= 8;
		occCount -= sum & 0x000000FF;	sum >>= 8;
		occCount -= sum & 0x000000FF;	sum >>= 8;
		occCount -= sum;
	}

	for (i=0; i<iteration; i++) {

		sum = 0;
		for (j=0; j<25; j++) {
			dna++;
			sum += dnaDecodeTable[*dna >> 17];
			sum += dnaDecodeTable[(*dna >> 2) & 0x00007FFF];
		}
		if (character != 4) {
			occCount += (sum >> (character * 8)) & 0x000000FF;
		} else {
			occCount -= sum & 0x000000FF;	sum >>= 8;
			occCount -= sum & 0x000000FF;	sum >>= 8;
			occCount -= sum & 0x000000FF;	sum >>= 8;
			occCount -= sum;
		}
	}

	return occCount;

}

void ForwardDNA_NAllOccCountNoLimit(const unsigned int*  dna, const unsigned int index, unsigned int* __restrict occCount, const unsigned int*  dnaDecodeTable) {

	static const unsigned int truncateRightMask[10] = { 0x00000000, 0xE0000000, 0xFC000000, 0xFF800000,
											   0xFFF00000, 0xFFFE0000, 0xFFFFC000, 0xFFFFF800,
											   0xFFFFFF00, 0xFFFFFFE0};

	unsigned int iteration, wordToCount, charToCount;
	unsigned int i, j, c;
	unsigned int sum;

	occCount[0] = 0;
	occCount[1] = 0;
	occCount[2] = 0;
	occCount[3] = 0;

	iteration = index / 250;
	wordToCount = (index - iteration * 250) / 10;
	charToCount = index - iteration * 250 - wordToCount * 10;

	for (i=0; i<iteration; i++) {

		sum = 0;
		for (j=0; j<25; j++) {
			sum += dnaDecodeTable[*dna >> 17];
			sum += dnaDecodeTable[(*dna >> 2) & 0x00007FFF];
			dna++;
		}
		occCount[0] += sum & 0x000000FF;	sum >>= 8;
		occCount[1] += sum & 0x000000FF;	sum >>= 8;
		occCount[2] += sum & 0x000000FF;	sum >>= 8;
		occCount[3] += sum;
	}

	sum = 0;
	for (j=0; j<wordToCount; j++) {
		sum += dnaDecodeTable[*dna >> 17];
		sum += dnaDecodeTable[(*dna >> 2) & 0x00007FFF];
		dna++;
	}

	if (charToCount > 0) {
		c = *dna & truncateRightMask[charToCount];	// increase count of 'a' by 10 - charToCount;
		sum += dnaDecodeTable[c >> 17];
		sum += dnaDecodeTable[(c >> 2) & 0x00007FFF];
		sum += charToCount - 10;	// decrease count of 'a' by 10 - charToCount
	}

	occCount[0] += sum & 0x000000FF;	sum >>= 8;
	occCount[1] += sum & 0x000000FF;	sum >>= 8;
	occCount[2] += sum & 0x000000FF;	sum >>= 8;
	occCount[3] += sum;

}

void BackwardDNA_NAllOccCountNoLimit(const unsigned int*  dna, const unsigned int index, unsigned int* __restrict occCount, const unsigned int*  dnaDecodeTable) {

	static const unsigned int truncateLeftMask[10] =  { 0x00000000, 0x0000001C, 0x000000FC, 0x000007FC,
											   0x00003FFC, 0x0001FFFC, 0x000FFFFC, 0x007FFFFC,
											   0x03FFFFFC, 0x1FFFFFFC};

	unsigned int iteration, wordToCount, charToCount;
	unsigned int i, j, c;
	unsigned int sum = 0;

	dna -= index / 10 + 1;

	iteration = index / 250;
	wordToCount = (index - iteration * 250) / 10;
	charToCount = index - iteration * 250 - wordToCount * 10;

	if (charToCount > 0) {
		c = *dna & truncateLeftMask[charToCount];	// increase count of 'a' by 10 - charToCount;
		sum += dnaDecodeTable[c >> 17];
		sum += dnaDecodeTable[(c >> 2) & 0x00007FFF];
		sum += charToCount - 10;	// decrease count of 'a' by 16 - charToCount
	}

	for (j=0; j<wordToCount; j++) {
		dna++;
		sum += dnaDecodeTable[*dna >> 17];
		sum += dnaDecodeTable[(*dna >> 2) & 0x00007FFF];
	}

	occCount[0] = sum & 0x000000FF;	sum >>= 8;
	occCount[1] = sum & 0x000000FF;	sum >>= 8;
	occCount[2] = sum & 0x000000FF;	sum >>= 8;
	occCount[3] = sum;

	for (i=0; i<iteration; i++) {

		sum = 0;
		for (j=0; j<25; j++) {
			dna++;
			sum += dnaDecodeTable[*dna >> 17];
			sum += dnaDecodeTable[(*dna >> 2) & 0x00007FFF];
		}
		occCount[0] += sum & 0x000000FF;	sum >>= 8;
		occCount[1] += sum & 0x000000FF;	sum >>= 8;
		occCount[2] += sum & 0x000000FF;	sum >>= 8;
		occCount[3] += sum;
	}

}

unsigned int ForwardOccCount(const unsigned int*  packed, const unsigned int index, const unsigned int character, const unsigned int alphabetSize) {

	unsigned int wordToCount, charToCount;
	unsigned int bitPerChar, charPerWord;
	unsigned int i, j, c;
	unsigned int occCount = 0;

	bitPerChar = ceilLog2(alphabetSize);
	charPerWord = BITS_IN_WORD / bitPerChar;

	wordToCount = index / charPerWord;
	charToCount = index - wordToCount * charPerWord;

	for (i=0; i<wordToCount; i++) {
		c = packed[i];
		for (j=0; j<charPerWord; j++) {
			if (c >> (BITS_IN_WORD - bitPerChar) == character) {
				occCount++;
			}
			c <<= bitPerChar;
		}
	}
	if (charToCount > 0) {
		c = packed[i];
		for (j=0; j<charToCount; j++) {
			if (c >> (BITS_IN_WORD - bitPerChar) == character) {
				occCount++;
			}
			c <<= bitPerChar;
		}
	}

	return occCount;

}


unsigned int BackwardOccCount(const unsigned int*  packed, const unsigned int index, const unsigned int character, const unsigned int alphabetSize) {

	unsigned int wordToCount, charToCount;
	unsigned int bitPerChar, charPerWord;
	unsigned int i, j, c;
	unsigned int occCount = 0;

	bitPerChar = ceilLog2(alphabetSize);
	charPerWord = BITS_IN_WORD / bitPerChar;

	wordToCount = index / charPerWord;
	charToCount = index - wordToCount * charPerWord;

	packed -= wordToCount + 1;

	if (charToCount > 0) {
		c = *packed << (bitPerChar * (charPerWord - charToCount));
		for (j=0; j<charToCount; j++) {
			if (c >> (BITS_IN_WORD - bitPerChar) == character) {
				occCount++;
			}
			c <<= bitPerChar;
		}
	}

	for (i=1; i<=wordToCount; i++) {
		packed++;
		c = *packed;
		for (j=0; j<charPerWord; j++) {
			if (c >> (BITS_IN_WORD - bitPerChar) == character) {
				occCount++;
			}
			c <<= bitPerChar;
		}
	}

	return occCount;
}

void ForwardAllOccCount(const unsigned int*  packed, const unsigned int index, const unsigned int alphabetSize, unsigned int*  occCount) {

	unsigned int wordToCount, charToCount;
	unsigned int bitPerChar, charPerWord;
	unsigned int i, j, c;

	bitPerChar = ceilLog2(alphabetSize);
	charPerWord = BITS_IN_WORD / bitPerChar;

	wordToCount = index / charPerWord;
	charToCount = index - wordToCount * charPerWord;

	for (i=0; i<wordToCount; i++) {
		c = packed[i];
		for (j=0; j<charPerWord; j++) {
			occCount[c >> (BITS_IN_WORD - bitPerChar)]++;
			c <<= bitPerChar;
		}
	}
	if (charToCount > 0) {
		c = packed[i];
		for (j=0; j<charToCount; j++) {
			occCount[c >> (BITS_IN_WORD - bitPerChar)]++;
			c <<= bitPerChar;
		}
	}

}

void BackwardAllOccCount(const unsigned int*  packed, const unsigned int index, const unsigned int alphabetSize, unsigned int*  occCount) {

	unsigned int wordToCount, charToCount;
	unsigned int bitPerChar, charPerWord;
	unsigned int i, j, c;

	bitPerChar = ceilLog2(alphabetSize);
	charPerWord = BITS_IN_WORD / bitPerChar;

	wordToCount = index / charPerWord;
	charToCount = index - wordToCount * charPerWord;

	packed -= wordToCount + 1;

	if (charToCount > 0) {
		c = *packed << (bitPerChar * (charPerWord - charToCount));
		for (j=0; j<charToCount; j++) {
			occCount[c >> (BITS_IN_WORD - bitPerChar)]++;
			c <<= bitPerChar;
		}
	}

	for (i=1; i<=wordToCount; i++) {
		packed++;
		c = *packed;
		for (j=0; j<charPerWord; j++) {
			occCount[c >> (BITS_IN_WORD - bitPerChar)]++;
			c <<= bitPerChar;
		}
	}

}

