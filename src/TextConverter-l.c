#include "config.h"
#ifdef MMX
/*

   TextConverter.c		Text Converter

   This module contains miscellaneous text conversion functions.

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
#include "TextConverter-l.h"
#include "MiscUtilities-l.h"
#include "r250-l.h"


unsigned int GetWordPackedText(const unsigned int *packedText, const unsigned int index, const unsigned int shift, const unsigned int numberOfBit, const unsigned int vacantBit) {

	unsigned int text;
	const static unsigned int mask[32] = { 0x00000000, 0x80000000, 0xC0000000, 0xE0000000,
								  0xF0000000, 0xF8000000, 0xFC000000, 0xFE000000,
								  0xFF000000, 0xFF800000, 0xFFC00000, 0xFFE00000,
								  0xFFF00000, 0xFFF80000, 0xFFFC0000, 0xFFFE0000,
								  0xFFFF0000, 0xFFFF8000, 0xFFFFC000, 0xFFFFE000,
								  0xFFFFF000, 0xFFFFF800, 0xFFFFFC00, 0xFFFFFE00,
								  0xFFFFFF00, 0xFFFFFF80, 0xFFFFFFC0, 0xFFFFFFE0,
								  0xFFFFFFF0, 0xFFFFFFF8, 0xFFFFFFFC, 0xFFFFFFFE };

	if (shift > 0) {
		// packedText should be allocated with at least 1 Word buffer initialized to zero
#ifdef DNA_ONLY
		text = (packedText[index] << shift) | (packedText[index + 1] >> (BITS_IN_WORD - shift));
#else
		text = (packedText[index] << shift) | (packedText[index + 1] >> (BITS_IN_WORD - shift) << vacantBit);
#endif
	} else {
		text = packedText[index];
	}

	if (numberOfBit < BITS_IN_WORD) {
		// Fill unused bit with zero
		text &= mask[numberOfBit];
	}

	return text;
}


unsigned int ReadCharMap(unsigned char *charMap, const char *inputFileName, const unsigned char defaultMapping) {

	FILE *inputFile;
	char c;
	unsigned int v, alphabetSize;

	inputFile = (FILE*)fopen64(inputFileName, "r");

	if (inputFile == NULL) {
		fprintf(stderr, "ReadCharMap() : Cannot open character map!\n");
		exit(1);
	}

	for (v=0; v<CHAR_MAP_SIZE; v++) {
		charMap[v] = defaultMapping;
	}

	alphabetSize = 0;

	while (!feof(inputFile)) {
		fscanf(inputFile, " %c %u \n", &c, &v);
		if (v > CHAR_MAP_SIZE) {
			fprintf(stderr, "ReadCharMap() : Invalid charMap!\n");
			return 0;
		}
		charMap[(unsigned int)c] = (unsigned char)v;
		if (v > alphabetSize) {
			alphabetSize = v;
		}
	}

	fclose(inputFile);

	alphabetSize++;

	return alphabetSize;

}

void GenerateReverseCharMap(const unsigned char *charMap, unsigned char *reverseCharMap) {

	unsigned int i, j;

	for (i=0; i<CHAR_MAP_SIZE; i++) {
		reverseCharMap[i] = INVALID_CHAR;
		for (j=0; j<CHAR_MAP_SIZE; j++) {
			if (charMap[j] == i) {
				reverseCharMap[i] = (unsigned char)j;
				break;
			}
		}
	}

}

unsigned int BitPerWordPackedChar(const unsigned int alphabetSize) {

	#ifdef DEBUG
	if (alphabetSize < 2) {
		fprintf(stderr, "BitPerWordPackedChar() : alphabetSize < 2!\n");
		exit(1);
	}
	#endif

	return ceilLog2(alphabetSize);

}

unsigned int TextLengthFromWordPacked(unsigned int wordPackedLength, unsigned int bitPerChar, unsigned int lastWordLength) {

	return (wordPackedLength - 1) * (BITS_IN_WORD / bitPerChar) + lastWordLength;

}

unsigned int WordPackedLengthFromText(unsigned int textLength, unsigned int bitPerChar) {

	return (textLength + (BITS_IN_WORD / bitPerChar) - 1) / (BITS_IN_WORD / bitPerChar);

}

unsigned int LastWordLength(unsigned int textLength, unsigned int bitPerChar) {

	return textLength % (BITS_IN_WORD / bitPerChar);

}

unsigned int BitPerBytePackedChar(const unsigned int alphabetSize) {

	unsigned int bitPerChar;

	#ifdef DEBUG
	if (alphabetSize < 2) {
		fprintf(stderr, "BitPerBytePackedChar() : alphabetSize < 2!\n");
		exit(1);
	}
	#endif

	bitPerChar = ceilLog2(alphabetSize);

	#ifdef DEBUG
	if (bitPerChar > BITS_IN_BYTE) {
		fprintf(stderr, "BitPerBytePackedChar() : bitPerChar > BITS_IN_BYTE!\n");
		exit(1);
	}
	#endif

	// Return the largest number of bit that does not affect packing efficiency
	if (BITS_IN_BYTE / (BITS_IN_BYTE / bitPerChar) > bitPerChar) {
		bitPerChar = BITS_IN_BYTE / (BITS_IN_BYTE / bitPerChar);
	}
	return bitPerChar;
}

unsigned int TextLengthFromBytePacked(unsigned int bytePackedLength, unsigned int bitPerChar, unsigned int lastByteLength) {

	if (bytePackedLength > ALL_ONE_MASK / (BITS_IN_BYTE / bitPerChar)) {
		fprintf(stderr, "TextLengthFromBytePacked(): text length > 2^32!\n");
		exit(1);
	}
	return (bytePackedLength - 1) * (BITS_IN_BYTE / bitPerChar) + lastByteLength;

}

unsigned int BytePackedLengthFromText(unsigned int textLength, unsigned int bitPerChar) {

	return (textLength + (BITS_IN_BYTE / bitPerChar) - 1) / (BITS_IN_BYTE / bitPerChar);

}

unsigned char LastByteLength(unsigned int textLength, unsigned int bitPerChar) {

	return (unsigned char)(textLength % (BITS_IN_BYTE / bitPerChar));

}

void ConvertTextToWordPacked(const unsigned char *input, unsigned int *output, const unsigned char *charMap, const unsigned int alphabetSize, const unsigned int textLength) {

	unsigned int bitPerChar, charPerWord;
	unsigned int i, j, k;
	unsigned int c;
	unsigned int charValue;

	bitPerChar = BitPerWordPackedChar(alphabetSize);
	charPerWord = BITS_IN_WORD / bitPerChar;

	for (i=0; i<textLength/charPerWord; i++) {
		c = 0;
		j = i * charPerWord;
		for (k=0; k<charPerWord; k++) {
			charValue = charMap[input[j+k]];
			if (charValue >= alphabetSize) {
				charValue = 0;
			}
			c = c | (charValue << (BITS_IN_WORD - (k+1) * bitPerChar));
		}
		output[i] = c;
	}
	if (i * charPerWord < textLength) {
		c = 0;
		j = i * charPerWord;
		for (k=0; j+k < textLength; k++) {
			charValue = charMap[input[j+k]];
			if (charValue >= alphabetSize) {
				charValue = 0;
			}
			c = c | (charValue << (BITS_IN_WORD - (k+1) * bitPerChar));
		}
		output[i] = c;
	}

}

void ConvertTextToBytePacked(const unsigned char *input, unsigned char *output, const unsigned char *charMap, const unsigned int alphabetSize, const unsigned int textLength) {

	unsigned int bitPerChar, charPerByte;
	unsigned int i, j, k;
	unsigned char c;

	bitPerChar = BitPerBytePackedChar(alphabetSize);
	charPerByte = BITS_IN_BYTE / bitPerChar;

	for (i=0; i<textLength/charPerByte; i++) {
		c = 0;
		j = i * charPerByte;
		for (k=0; k<charPerByte; k++) {
			c = c | (unsigned char)(charMap[input[j+k]] << (BITS_IN_BYTE - (k+1) * bitPerChar));
		}
		output[i] = c;
	}
	if (i * charPerByte < textLength) {
		c = 0;
		j = i * charPerByte;
		for (k=0; j+k < textLength; k++) {
			c = c | (unsigned char)(charMap[input[j+k]] << (BITS_IN_BYTE - (k+1) * bitPerChar));
		}
		output[i] = c;
	}

}

void ConvertWordPackedToText(const unsigned int *input, unsigned char *output, const unsigned char *reverseCharMap, const unsigned int alphabetSize, const unsigned int textLength) {

	unsigned int bitPerChar, charPerWord;
	unsigned int i, j, k;
	unsigned int c;

	bitPerChar = BitPerWordPackedChar(alphabetSize);
	charPerWord = BITS_IN_WORD / bitPerChar;

	for (i=0; i<textLength/charPerWord; i++) {
		c = input[i];
		j = i * charPerWord;
		for (k=0; k<charPerWord; k++) {
			output[j+k] = reverseCharMap[c >> (BITS_IN_WORD - bitPerChar)];
			c <<= bitPerChar;
		}
	}
	if (i * charPerWord < textLength) {
		c = input[i];
		j = i * charPerWord;
		for (k=0; j+k<textLength; k++) {
			output[j+k] = reverseCharMap[c >> (BITS_IN_WORD - bitPerChar)];
			c <<= bitPerChar;
		}
	}

}

void ConvertBytePackedToText(const unsigned char *input, unsigned char *output, const unsigned char *reverseCharMap, const unsigned int alphabetSize, const unsigned int textLength) {

	unsigned int bitPerChar, charPerByte;
	unsigned int i, j, k;
	unsigned char c;

	bitPerChar = BitPerBytePackedChar(alphabetSize);
	charPerByte = BITS_IN_BYTE / bitPerChar;

	for (i=0; i<textLength/charPerByte; i++) {
		c = input[i];
		j = i * charPerByte;
		for (k=0; k<charPerByte; k++) {
			output[j+k] = reverseCharMap[c >> (BITS_IN_BYTE - bitPerChar)];
			c <<= bitPerChar;
		}
	}
	if (i * charPerByte < textLength) {
		c = input[i];
		j = i * charPerByte;
		for (k=0; j+k<textLength; k++) {
			output[j+k] = reverseCharMap[c >> (BITS_IN_BYTE - bitPerChar)];
			c <<= bitPerChar;
		}
	}

}

void ConvertBytePackedToCode(const unsigned char *input, unsigned char *output, const unsigned int alphabetSize, const unsigned int textLength) {

	unsigned int bitPerChar, charPerByte;
	unsigned int i, j, k;
	unsigned char c;

	bitPerChar = BitPerBytePackedChar(alphabetSize);
	charPerByte = BITS_IN_BYTE / bitPerChar;

	for (i=0; i<textLength/charPerByte; i++) {
		c = input[i];
		j = i * charPerByte;
		for (k=0; k<charPerByte; k++) {
			output[j+k] = c >> (unsigned char)(BITS_IN_BYTE - bitPerChar);
			c <<= bitPerChar;
		}
	}
	if (i * charPerByte < textLength) {
		c = input[i];
		j = i * charPerByte;
		for (k=0; j+k<textLength; k++) {
			output[j+k] = c >> (unsigned char)(BITS_IN_BYTE - bitPerChar);
			c <<= bitPerChar;
		}
	}

}

void ConvertWordPackedToBytePacked(const unsigned int *input, unsigned char *output, const unsigned int alphabetSize, const unsigned int textLength) {

	unsigned int i, j, k;
	unsigned int c;
	unsigned int bitPerBytePackedChar;
	unsigned int bitPerWordPackedChar;
	unsigned int charPerWord;
	unsigned int charPerByte;
	unsigned int bytePerIteration;
	unsigned int byteProcessed = 0;
	unsigned int wordProcessed = 0;
	unsigned int mask, shift;
	
	unsigned int buffer[BITS_IN_WORD];

	bitPerBytePackedChar = BitPerBytePackedChar(alphabetSize);
	bitPerWordPackedChar = BitPerWordPackedChar(alphabetSize);
	charPerWord = BITS_IN_WORD / bitPerBytePackedChar;
	charPerByte = BITS_IN_BYTE / bitPerWordPackedChar;

	bytePerIteration = charPerWord / charPerByte;
	mask = truncateRight(ALL_ONE_MASK, BITS_IN_WORD - bitPerWordPackedChar);
	shift = BITS_IN_WORD - bitPerWordPackedChar;

	while ((wordProcessed + 1) * charPerWord < textLength) {

		c = input[wordProcessed];
		for (i=0; i<charPerWord; i++) {
			buffer[i] = c >> shift;
			c <<= bitPerWordPackedChar;
		}
		wordProcessed++;

		k = 0;
		for (i=0; i<bytePerIteration; i++) {
			c = 0;
			for (j=0; j<charPerByte; j++) {
				c |= buffer[k] << (BITS_IN_BYTE - (j+1) * bitPerBytePackedChar);
				k++;
			}
			output[byteProcessed] = (unsigned char)c;
			byteProcessed++;
		}

	}

	c = input[wordProcessed];
	for (i=0; i < textLength - wordProcessed * charPerWord; i++) {
		buffer[i] = c >> shift;
		c <<= bitPerWordPackedChar;
	}

	k = 0;
	while (byteProcessed * charPerByte < textLength) {
		c = 0;
		for (j=0; j < textLength - wordProcessed * charPerWord; j++) {
			c |= buffer[k] << (BITS_IN_BYTE - (j+1) * bitPerBytePackedChar);
			k++;
		}
		output[byteProcessed] = (unsigned char)c;
		byteProcessed++;
	}

}

void ConvertBytePackedToWordPacked(const unsigned char *input, unsigned int *output, const unsigned int alphabetSize, const unsigned int textLength) {

	unsigned int i, j, k;
	unsigned int c;
	unsigned int bitPerBytePackedChar;
	unsigned int bitPerWordPackedChar;
	unsigned int charPerWord;
	unsigned int charPerByte;
	unsigned int bytePerIteration;
	unsigned int byteProcessed = 0;
	unsigned int wordProcessed = 0;
	unsigned int mask, shift;
	
	unsigned int buffer[BITS_IN_WORD];

	bitPerBytePackedChar = BitPerBytePackedChar(alphabetSize);
	bitPerWordPackedChar = BitPerWordPackedChar(alphabetSize);
	charPerByte = BITS_IN_BYTE / bitPerBytePackedChar;
	charPerWord = BITS_IN_WORD / bitPerWordPackedChar;

	bytePerIteration = charPerWord / charPerByte;
	mask = truncateRight(ALL_ONE_MASK, BITS_IN_WORD - bitPerWordPackedChar);
	shift = BITS_IN_WORD - BITS_IN_BYTE + bitPerBytePackedChar - bitPerWordPackedChar;

	while ((wordProcessed + 1) * charPerWord < textLength) {

		k = 0;
		for (i=0; i<bytePerIteration; i++) {
			c = (unsigned int)input[byteProcessed] << shift;
			for (j=0; j<charPerByte; j++) {
				buffer[k] = c & mask;
				c <<= bitPerBytePackedChar;
				k++;
			}
			byteProcessed++;
		}

		c = 0;
		for (i=0; i<charPerWord; i++) {
			c |= buffer[i] >> bitPerWordPackedChar * i;
		}
		output[wordProcessed] = c;
		wordProcessed++;

	}

	k = 0;
	for (i=0; i < (textLength - wordProcessed * charPerWord - 1) / charPerByte + 1; i++) {
		c = (unsigned int)input[byteProcessed] << shift;
		for (j=0; j<charPerByte; j++) {
			buffer[k] = c & mask;
			c <<= bitPerBytePackedChar;
			k++;
		}
		byteProcessed++;
	}

	c = 0;
	for (i=0; i<textLength - wordProcessed * charPerWord; i++) {
		c |= buffer[i] >> bitPerWordPackedChar * i;
	}
	output[wordProcessed] = c;

}

void ConvertTextToCode(const unsigned char *input, unsigned char *output, const unsigned char *charMap, const unsigned int textLength) {

	unsigned int i;

	for (i=0; i< textLength; i++) {
		output[i] = charMap[input[i]];
	}

}

void ConvertCodeToText(const unsigned char *input, unsigned char *output, const unsigned char *reverseCharMap, const unsigned int textLength) {

	unsigned int i;

	for (i=0; i< textLength; i++) {
		output[i] = reverseCharMap[input[i]];
	}

}

void PackTextWithAllShift(const unsigned char *input, unsigned int **output, const unsigned char *charMap, const unsigned int alphabetSize, const unsigned int textLength) {

	unsigned int bitPerChar;
	unsigned int numberOfShift;
	unsigned int numberOfWord;
	unsigned int shift;

	unsigned int i, j;

	bitPerChar = BitPerWordPackedChar(alphabetSize);
	numberOfShift = BITS_IN_WORD / bitPerChar;
	numberOfWord = WordPackedLengthFromText(textLength, bitPerChar);

	ConvertTextToWordPacked(input, output[0], charMap, alphabetSize, textLength);

	for (i=1; i<numberOfShift; i++) {
		shift = i * bitPerChar;
		output[i][0] = output[0][0] >> shift;
		for (j=1; j<=numberOfWord; j++) {
			output[i][j] = (output[0][j] >> shift) | (output[0][j-1] << (BITS_IN_WORD - shift));
		}
	}

}


unsigned int ReadTextAsWordPacked(const char *inputFileName, const unsigned char *charMap, const unsigned int alphabetSize, unsigned int *targetAddress, const unsigned int maxTextLength) {

	FILE *inputFile;
	unsigned char *buffer;
	unsigned int charPerWord;
	unsigned int charRead;
	unsigned int charProcessed = 0, wordProcessed = 0;
	unsigned int charPerBuffer;

	inputFile = (FILE*)fopen64(inputFileName, "rb");

	if (inputFile == NULL) {
		fprintf(stderr, "ReadTextAsWordPacked() : Cannot open inputFileName!\n");
		exit(1);
	}

	charPerWord = BITS_IN_WORD / BitPerWordPackedChar(alphabetSize);
	charPerBuffer = PACKED_BUFFER_SIZE / charPerWord * charPerWord;

	buffer = MMUnitAllocate(charPerBuffer);

	charRead = (unsigned int)fread(buffer, 1, charPerBuffer, inputFile);
	while (charRead > 0 && charProcessed + charRead < maxTextLength) {
		ConvertTextToWordPacked(buffer, targetAddress + wordProcessed, charMap, alphabetSize, charRead);
		wordProcessed += charRead / charPerWord;
		charProcessed += charRead;
		charRead = (unsigned int)fread(buffer, 1, charPerBuffer, inputFile);
	}

	if (charRead > 0 && charProcessed < maxTextLength) {
		ConvertTextToWordPacked(buffer, targetAddress + wordProcessed, charMap, alphabetSize, minX(charRead, maxTextLength - charProcessed));
		charProcessed += charRead;
	}

	MMUnitFree(buffer, charPerBuffer);

	fclose(inputFile);

	return charProcessed;

}

unsigned int ReadBytePackedAsWordPacked(const char *inputFileName, const unsigned int alphabetSize, unsigned int *targetAddress, const unsigned int maxTextLength) {

	FILE *inputFile;
	unsigned char *buffer1, *buffer2;
	unsigned int charPerByte, charPerWord;
	unsigned int charPerBuffer, wordPerBuffer;
	unsigned int charProcessed = 0, wordProcessed = 0;
	unsigned int byteRead, tempByteRead;
	unsigned int charInLastBuffer;
	unsigned int bufferSize;

	inputFile = (FILE*)fopen64(inputFileName, "rb");

	if (inputFile == NULL) {
		fprintf(stderr, "ReadBytePackedAsWordPacked() : Cannot open inputFileName!\n");
		exit(1);
	}

	charPerByte = BITS_IN_BYTE / BitPerBytePackedChar(alphabetSize);
	charPerWord = BITS_IN_WORD / BitPerWordPackedChar(alphabetSize);
	bufferSize = PACKED_BUFFER_SIZE / charPerByte / charPerWord * charPerByte * charPerWord;

	charPerBuffer = bufferSize * charPerByte;
	wordPerBuffer = charPerBuffer / charPerWord;

	buffer1 = MMUnitAllocate(bufferSize);
	buffer2 = MMUnitAllocate(bufferSize);

	byteRead = (unsigned int)fread(buffer1, 1, bufferSize, inputFile);
	tempByteRead = (unsigned int)fread(buffer2, 1, bufferSize, inputFile);

	while (tempByteRead > 1 && charProcessed + charPerBuffer < maxTextLength) {
		ConvertBytePackedToWordPacked(buffer1, targetAddress + wordProcessed, alphabetSize, charPerBuffer);
		charProcessed += charPerBuffer;
		wordProcessed += wordPerBuffer;
		memcpy(buffer1, buffer2, bufferSize);
		byteRead = tempByteRead;
		tempByteRead = (unsigned int)fread(buffer2, 1, bufferSize, inputFile);
	}

	if (tempByteRead > 1) {
		ConvertBytePackedToWordPacked(buffer1, targetAddress + wordProcessed, alphabetSize, maxTextLength - charProcessed);
		charProcessed += charPerBuffer;
	} else {
		if (tempByteRead == 1) {
			charInLastBuffer = charPerBuffer - charPerByte + buffer2[0];
		} else {
			charInLastBuffer = (byteRead - 2) * charPerByte + buffer1[byteRead - 1];
		}
		ConvertBytePackedToWordPacked(buffer1, targetAddress + wordProcessed, alphabetSize, minX(maxTextLength - charProcessed, charInLastBuffer));
		charProcessed += charInLastBuffer;
	}

	MMUnitFree(buffer1, bufferSize);
	MMUnitFree(buffer2, bufferSize);

	fclose(inputFile);

	return charProcessed;

}

// Alphabet size of DNA must be 4
void *DNALoadPacked(const char *inputFileName, unsigned int *textLength, const unsigned int convertToWordPacked, const unsigned int trailerBufferInWord) {

	FILE *inputFile;
	unsigned char tempChar[4];
	unsigned int *packedText;
	unsigned int packedFileLen;
	unsigned char lastByteLength;
	unsigned int wordToProcess;
	unsigned int i;
	unsigned int trailerBufferIn128;

	trailerBufferIn128 = (trailerBufferInWord + 3) / 4 * 4;

	inputFile = (FILE*)(FILE*)fopen64(inputFileName, "rb");

	if (inputFile == NULL) {
		fprintf(stderr, "DNALoadPacked() : Cannot open inputFileName!\n");
		exit(1);
	}

	fseek(inputFile, -1, SEEK_END);
	packedFileLen = ftell(inputFile);
	if ((int)packedFileLen < 0) {
		fprintf(stderr, "DNALoadPacked(): Cannot determine file length!\n");
		exit(1);
	}
	fread(&lastByteLength, sizeof(unsigned char), 1, inputFile);

	*textLength = (packedFileLen - 1) * 4 + lastByteLength;

	wordToProcess = (*textLength + 64 - 1) / 64 * 4 + trailerBufferIn128 * 4;		// allocate multiple of 128 bit + trailer buffer

	packedText = MMUnitAllocate(wordToProcess * sizeof(unsigned int));
	for (i=(*textLength)/16; i<wordToProcess; i++) {
		packedText[i] = 0;
	}

	fseek(inputFile, 0, SEEK_SET);
	fread(packedText, 1, packedFileLen, inputFile);
	fclose(inputFile);

	if (convertToWordPacked) {

		for (i=0; i<wordToProcess; i++) {
	
			*(unsigned int*)tempChar = packedText[i];
			packedText[i] = (tempChar[0] << 24) | (tempChar[1] << 16) | (tempChar[2] << 8) | tempChar[3];

		}

	}

	return (void*)packedText;

}

void DNAFreePacked(void* packedDNA, const unsigned int textLength, const unsigned int trailerBufferInWord) {

	unsigned int trailerBufferIn128;

	trailerBufferIn128 = (trailerBufferInWord + 3) / 4 * 4;

	MMUnitFree(packedDNA, ((textLength + 64 - 1) / 64 * 4 + trailerBufferIn128 * 4) * sizeof(unsigned int));

}

void SaveText(const char *outputFileName, const unsigned char *text, const unsigned int textLength) {

	FILE *outputFile;

	outputFile = (FILE*)fopen64(outputFileName, "wb");

	if (outputFile == NULL) {
		fprintf(stderr, "SaveText() : Cannot open output file!\n");
		exit(1);
	}

	fwrite(text, sizeof(unsigned char), textLength, outputFile);
	fclose(outputFile);

}

void SaveBytePacked(const char *outputFileName, const unsigned char *bytePacked, const unsigned int textLength, const unsigned int alphabetSize) {

	FILE *outputFile;
	unsigned int bitPerChar, charPerByte, bytePackedLen;
	unsigned char lastByteLen;
	unsigned char zero = 0;

	outputFile = (FILE*)fopen64(outputFileName, "wb");

	if (outputFile == NULL) {
		fprintf(stderr, "SaveBytePacked() : Cannot open output file!\n");
		exit(1);
	}

	bitPerChar = BitPerBytePackedChar(alphabetSize);
	charPerByte = BITS_IN_BYTE / bitPerChar;

	bytePackedLen = BytePackedLengthFromText(textLength, bitPerChar);
	lastByteLen = LastByteLength(textLength, bitPerChar);

	fwrite(bytePacked, sizeof(unsigned char), bytePackedLen, outputFile);
	if (lastByteLen == 0) {
		fwrite(&zero, sizeof(unsigned char), 1, outputFile);
	}
	fwrite(&lastByteLen, sizeof(unsigned char), 1, outputFile);
	fclose(outputFile);

}

void SaveWordPacked(const char *outputFileName, const unsigned int *wordPacked, const unsigned int textLength, const unsigned int alphabetSize) {

	FILE *outputFile;
	unsigned int bitPerChar, charPerWord, wordPackedLen;
	unsigned int lastWordLen;
	unsigned int zero = 0;

	outputFile = (FILE*)fopen64(outputFileName, "wb");

	if (outputFile == NULL) {
		fprintf(stderr, "SaveWordPacked() : Cannot open output file!\n");
		exit(1);
	}

	bitPerChar = BitPerWordPackedChar(alphabetSize);
	charPerWord = BITS_IN_WORD / bitPerChar;

	wordPackedLen = WordPackedLengthFromText(textLength, bitPerChar);
	lastWordLen = LastWordLength(textLength, bitPerChar);

	fwrite(wordPacked, sizeof(unsigned int), wordPackedLen, outputFile);
	if (lastWordLen == 0) {
		fwrite(&zero, sizeof(unsigned int), 1, outputFile);
	}
	fwrite(&lastWordLen, sizeof(unsigned int), 1, outputFile);
	fclose(outputFile);

}

FILE *InitialLoadPackedIncFromEnd(const char* inputFileName, unsigned char *packedOutput, const unsigned int alphabetSize, 
								  const unsigned int packedLengthPerLoad, unsigned int *textLength, unsigned int *textLengthForThisLoad) {

	FILE *packedFile;
	unsigned int len, packedFileLenForThisLoad, packedFileLen;
	unsigned char lastByteLength;
	unsigned int bitPerChar, charPerWord;

	packedFile = (FILE*)fopen64(inputFileName, "rb");

	if (packedFile == NULL) {
		fprintf(stderr, "InitialLoadPackedIncFromEnd() : Cannot open inputFileName!\n");
		exit(1);
	}

	bitPerChar = BitPerBytePackedChar(alphabetSize);
	charPerWord = BITS_IN_WORD / bitPerChar;

	fseek(packedFile, -1, SEEK_END);
	packedFileLen = ftell(packedFile);
	if ((int)packedFileLen < 0) {
		fprintf(stderr, "InitialLoadPackedIncFromEnd(): Cannot determine file length!\n");
		exit(1);
	}
	fread(&lastByteLength, sizeof(unsigned char), 1, packedFile);

	len = TextLengthFromBytePacked(packedFileLen, bitPerChar, lastByteLength);

	if (lastByteLength == 0 && (packedFileLen - 1) % packedLengthPerLoad == 0) {
		packedFileLenForThisLoad = 0;
		fseek(packedFile, -((int)(2+packedLengthPerLoad)), SEEK_END);
		*textLength = len;
		*textLengthForThisLoad = 0;
		return packedFile;
	}

	if (packedFileLen % packedLengthPerLoad == 0) {
		packedFileLenForThisLoad = packedLengthPerLoad;
	} else {
		packedFileLenForThisLoad = packedFileLen % packedLengthPerLoad;
	}
	fseek(packedFile, -1, SEEK_END);

	fseek(packedFile, -((int)packedFileLenForThisLoad), SEEK_CUR);
	fread(packedOutput, sizeof(unsigned char), packedFileLenForThisLoad, packedFile);
	fseek(packedFile, -((int)packedFileLenForThisLoad), SEEK_CUR);
	if (packedFileLen > packedFileLenForThisLoad) {
		fseek(packedFile, -((int)packedLengthPerLoad), SEEK_CUR);
	}

	*textLength = len;
	*textLengthForThisLoad = TextLengthFromBytePacked(packedFileLenForThisLoad, bitPerChar, lastByteLength);

	return packedFile;

}

void LoadPackedIncFromEnd(FILE *packedFile, unsigned char *packedOutput, const unsigned int packedLengthPerLoad) {
	
	fread(packedOutput, sizeof(unsigned char), packedLengthPerLoad, packedFile);
	fseek(packedFile, -(2*(int)packedLengthPerLoad), SEEK_CUR);

}


FILE *InitialLoadTextIncFromEnd(const char* inputFileName, unsigned char *textOutput, const unsigned int textLengthPerLoad, unsigned int *textLength, unsigned int *textLengthForThisLoad) {

	FILE *textFile;
	unsigned int len, textLenForThisLoad;

	textFile = (FILE*)fopen64(inputFileName, "rb");

	if (textFile == NULL) {
		fprintf(stderr, "InitialLoadTextIncFromEnd() : Cannot open inputFileName!\n");
		exit(1);
	}

	fseek(textFile, 0, SEEK_END);
	len = ftell(textFile);
	if ((int)len < 0) {
		fprintf(stderr, "InitialLoadTextIncFromEnd(): Cannot determine file length!\n");
		exit(1);
	}

	textLenForThisLoad = len % textLengthPerLoad;

	if (textLenForThisLoad > 0) {
		fseek(textFile, -((int)textLenForThisLoad), SEEK_END);
		fread(textOutput, sizeof(unsigned char), textLenForThisLoad, textFile);
		fseek(textFile, -((int)textLenForThisLoad), SEEK_END);
	}

	*textLength = len;
	*textLengthForThisLoad = textLenForThisLoad;

	return textFile;
}

void LoadTextIncFromEnd(FILE *textFile, unsigned char *textOutput, const unsigned int textLengthPerLoad) {

	if (ftell(textFile) < (int)textLengthPerLoad) {
		fprintf(stderr, "LoadTextIncFromEnd(): file pointer is not correctly placed!\n");
		exit(1);
	}

	fseek(textFile, -((int)textLengthPerLoad), SEEK_CUR);
	fread(textOutput, sizeof(unsigned char), textLengthPerLoad, textFile);
	fseek(textFile, -((int)textLengthPerLoad), SEEK_CUR);

}
#else

/*

   TextConverter.c		Text Converter

   This module contains miscellaneous text conversion functions.

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
#include "TextConverter-l.h"
#include "MiscUtilities-l.h"
#include "r250-l.h"


unsigned int GetWordPackedText(const unsigned int *packedText, const unsigned int index, const unsigned int shift, const unsigned int numberOfBit, const unsigned int vacantBit) {

	unsigned int text;
	const static unsigned int mask[32] = { 0x00000000, 0x80000000, 0xC0000000, 0xE0000000,
								  0xF0000000, 0xF8000000, 0xFC000000, 0xFE000000,
								  0xFF000000, 0xFF800000, 0xFFC00000, 0xFFE00000,
								  0xFFF00000, 0xFFF80000, 0xFFFC0000, 0xFFFE0000,
								  0xFFFF0000, 0xFFFF8000, 0xFFFFC000, 0xFFFFE000,
								  0xFFFFF000, 0xFFFFF800, 0xFFFFFC00, 0xFFFFFE00,
								  0xFFFFFF00, 0xFFFFFF80, 0xFFFFFFC0, 0xFFFFFFE0,
								  0xFFFFFFF0, 0xFFFFFFF8, 0xFFFFFFFC, 0xFFFFFFFE };

	if (shift > 0) {
		// packedText should be allocated with at least 1 Word buffer initialized to zero
#ifdef DNA_ONLY
		text = (packedText[index] << shift) | (packedText[index + 1] >> (BITS_IN_WORD - shift));
#else
		text = (packedText[index] << shift) | (packedText[index + 1] >> (BITS_IN_WORD - shift) << vacantBit);
#endif
	} else {
		text = packedText[index];
	}

	if (numberOfBit < BITS_IN_WORD) {
		// Fill unused bit with zero
		text &= mask[numberOfBit];
	}

	return text;
}


unsigned int ReadCharMap(unsigned char *charMap, const char *inputFileName, const unsigned char defaultMapping) {

	FILE *inputFile;
	char c;
	unsigned int v, alphabetSize;

	inputFile = (FILE*)fopen64(inputFileName, "r");

	if (inputFile == NULL) {
		fprintf(stderr, "ReadCharMap() : Cannot open character map!\n");
		exit(1);
	}

	for (v=0; v<CHAR_MAP_SIZE; v++) {
		charMap[v] = defaultMapping;
	}

	alphabetSize = 0;

	while (!feof(inputFile)) {
		fscanf(inputFile, " %c %u \n", &c, &v);
		if (v > CHAR_MAP_SIZE) {
			fprintf(stderr, "ReadCharMap() : Invalid charMap!\n");
			return 0;
		}
		charMap[(unsigned int)c] = (unsigned char)v;
		if (v > alphabetSize) {
			alphabetSize = v;
		}
	}

	fclose(inputFile);

	alphabetSize++;

	return alphabetSize;

}

void GenerateReverseCharMap(const unsigned char *charMap, unsigned char *reverseCharMap) {

	unsigned int i, j;

	for (i=0; i<CHAR_MAP_SIZE; i++) {
		reverseCharMap[i] = INVALID_CHAR;
		for (j=0; j<CHAR_MAP_SIZE; j++) {
			if (charMap[j] == i) {
				reverseCharMap[i] = (unsigned char)j;
				break;
			}
		}
	}

}

unsigned int BitPerWordPackedChar(const unsigned int alphabetSize) {

	#ifdef DEBUG
	if (alphabetSize < 2) {
		fprintf(stderr, "BitPerWordPackedChar() : alphabetSize < 2!\n");
		exit(1);
	}
	#endif

	return ceilLog2(alphabetSize);

}

unsigned int TextLengthFromWordPacked(unsigned int wordPackedLength, unsigned int bitPerChar, unsigned int lastWordLength) {

	return (wordPackedLength - 1) * (BITS_IN_WORD / bitPerChar) + lastWordLength;

}

unsigned int WordPackedLengthFromText(unsigned int textLength, unsigned int bitPerChar) {

	return (textLength + (BITS_IN_WORD / bitPerChar) - 1) / (BITS_IN_WORD / bitPerChar);

}

unsigned int LastWordLength(unsigned int textLength, unsigned int bitPerChar) {

	return textLength % (BITS_IN_WORD / bitPerChar);

}

unsigned int BitPerBytePackedChar(const unsigned int alphabetSize) {

	unsigned int bitPerChar;

	#ifdef DEBUG
	if (alphabetSize < 2) {
		fprintf(stderr, "BitPerBytePackedChar() : alphabetSize < 2!\n");
		exit(1);
	}
	#endif

	bitPerChar = ceilLog2(alphabetSize);

	#ifdef DEBUG
	if (bitPerChar > BITS_IN_BYTE) {
		fprintf(stderr, "BitPerBytePackedChar() : bitPerChar > BITS_IN_BYTE!\n");
		exit(1);
	}
	#endif

	// Return the largest number of bit that does not affect packing efficiency
	if (BITS_IN_BYTE / (BITS_IN_BYTE / bitPerChar) > bitPerChar) {
		bitPerChar = BITS_IN_BYTE / (BITS_IN_BYTE / bitPerChar);
	}
	return bitPerChar;
}

unsigned int TextLengthFromBytePacked(unsigned int bytePackedLength, unsigned int bitPerChar, unsigned int lastByteLength) {

	if (bytePackedLength > ALL_ONE_MASK / (BITS_IN_BYTE / bitPerChar)) {
		fprintf(stderr, "TextLengthFromBytePacked(): text length > 2^32!\n");
		exit(1);
	}
	return (bytePackedLength - 1) * (BITS_IN_BYTE / bitPerChar) + lastByteLength;

}

unsigned int BytePackedLengthFromText(unsigned int textLength, unsigned int bitPerChar) {

	return (textLength + (BITS_IN_BYTE / bitPerChar) - 1) / (BITS_IN_BYTE / bitPerChar);

}

unsigned char LastByteLength(unsigned int textLength, unsigned int bitPerChar) {

	return (unsigned char)(textLength % (BITS_IN_BYTE / bitPerChar));

}

void ConvertTextToWordPacked(const unsigned char *input, unsigned int *output, const unsigned char *charMap, const unsigned int alphabetSize, const unsigned int textLength) {

	unsigned int bitPerChar, charPerWord;
	unsigned int i, j, k;
	unsigned int c;
	unsigned int charValue;

	bitPerChar = BitPerWordPackedChar(alphabetSize);
	charPerWord = BITS_IN_WORD / bitPerChar;

	for (i=0; i<textLength/charPerWord; i++) {
		c = 0;
		j = i * charPerWord;
		for (k=0; k<charPerWord; k++) {
			charValue = charMap[input[j+k]];
			if (charValue >= alphabetSize) {
				charValue = 0;
			}
			c = c | (charValue << (BITS_IN_WORD - (k+1) * bitPerChar));
		}
		output[i] = c;
	}
	if (i * charPerWord < textLength) {
		c = 0;
		j = i * charPerWord;
		for (k=0; j+k < textLength; k++) {
			charValue = charMap[input[j+k]];
			if (charValue >= alphabetSize) {
				charValue = 0;
			}
			c = c | (charValue << (BITS_IN_WORD - (k+1) * bitPerChar));
		}
		output[i] = c;
	}

}

void ConvertTextToBytePacked(const unsigned char *input, unsigned char *output, const unsigned char *charMap, const unsigned int alphabetSize, const unsigned int textLength) {

	unsigned int bitPerChar, charPerByte;
	unsigned int i, j, k;
	unsigned char c;

	bitPerChar = BitPerBytePackedChar(alphabetSize);
	charPerByte = BITS_IN_BYTE / bitPerChar;

	for (i=0; i<textLength/charPerByte; i++) {
		c = 0;
		j = i * charPerByte;
		for (k=0; k<charPerByte; k++) {
			c = c | (unsigned char)(charMap[input[j+k]] << (BITS_IN_BYTE - (k+1) * bitPerChar));
		}
		output[i] = c;
	}
	if (i * charPerByte < textLength) {
		c = 0;
		j = i * charPerByte;
		for (k=0; j+k < textLength; k++) {
			c = c | (unsigned char)(charMap[input[j+k]] << (BITS_IN_BYTE - (k+1) * bitPerChar));
		}
		output[i] = c;
	}

}

void ConvertWordPackedToText(const unsigned int *input, unsigned char *output, const unsigned char *reverseCharMap, const unsigned int alphabetSize, const unsigned int textLength) {

	unsigned int bitPerChar, charPerWord;
	unsigned int i, j, k;
	unsigned int c;

	bitPerChar = BitPerWordPackedChar(alphabetSize);
	charPerWord = BITS_IN_WORD / bitPerChar;

	for (i=0; i<textLength/charPerWord; i++) {
		c = input[i];
		j = i * charPerWord;
		for (k=0; k<charPerWord; k++) {
			output[j+k] = reverseCharMap[c >> (BITS_IN_WORD - bitPerChar)];
			c <<= bitPerChar;
		}
	}
	if (i * charPerWord < textLength) {
		c = input[i];
		j = i * charPerWord;
		for (k=0; j+k<textLength; k++) {
			output[j+k] = reverseCharMap[c >> (BITS_IN_WORD - bitPerChar)];
			c <<= bitPerChar;
		}
	}

}

void ConvertBytePackedToText(const unsigned char *input, unsigned char *output, const unsigned char *reverseCharMap, const unsigned int alphabetSize, const unsigned int textLength) {

	unsigned int bitPerChar, charPerByte;
	unsigned int i, j, k;
	unsigned char c;

	bitPerChar = BitPerBytePackedChar(alphabetSize);
	charPerByte = BITS_IN_BYTE / bitPerChar;

	for (i=0; i<textLength/charPerByte; i++) {
		c = input[i];
		j = i * charPerByte;
		for (k=0; k<charPerByte; k++) {
			output[j+k] = reverseCharMap[c >> (BITS_IN_BYTE - bitPerChar)];
			c <<= bitPerChar;
		}
	}
	if (i * charPerByte < textLength) {
		c = input[i];
		j = i * charPerByte;
		for (k=0; j+k<textLength; k++) {
			output[j+k] = reverseCharMap[c >> (BITS_IN_BYTE - bitPerChar)];
			c <<= bitPerChar;
		}
	}

}

void ConvertBytePackedToCode(const unsigned char *input, unsigned char *output, const unsigned int alphabetSize, const unsigned int textLength) {

	unsigned int bitPerChar, charPerByte;
	unsigned int i, j, k;
	unsigned char c;

	bitPerChar = BitPerBytePackedChar(alphabetSize);
	charPerByte = BITS_IN_BYTE / bitPerChar;

	for (i=0; i<textLength/charPerByte; i++) {
		c = input[i];
		j = i * charPerByte;
		for (k=0; k<charPerByte; k++) {
			output[j+k] = c >> (unsigned char)(BITS_IN_BYTE - bitPerChar);
			c <<= bitPerChar;
		}
	}
	if (i * charPerByte < textLength) {
		c = input[i];
		j = i * charPerByte;
		for (k=0; j+k<textLength; k++) {
			output[j+k] = c >> (unsigned char)(BITS_IN_BYTE - bitPerChar);
			c <<= bitPerChar;
		}
	}

}

void ConvertWordPackedToBytePacked(const unsigned int *input, unsigned char *output, const unsigned int alphabetSize, const unsigned int textLength) {

	unsigned int i, j, k;
	unsigned int c;
	unsigned int bitPerBytePackedChar;
	unsigned int bitPerWordPackedChar;
	unsigned int charPerWord;
	unsigned int charPerByte;
	unsigned int bytePerIteration;
	unsigned int byteProcessed = 0;
	unsigned int wordProcessed = 0;
	unsigned int mask, shift;
	
	unsigned int buffer[BITS_IN_WORD];

	bitPerBytePackedChar = BitPerBytePackedChar(alphabetSize);
	bitPerWordPackedChar = BitPerWordPackedChar(alphabetSize);
	charPerWord = BITS_IN_WORD / bitPerBytePackedChar;
	charPerByte = BITS_IN_BYTE / bitPerWordPackedChar;

	bytePerIteration = charPerWord / charPerByte;
	mask = truncateRight(ALL_ONE_MASK, BITS_IN_WORD - bitPerWordPackedChar);
	shift = BITS_IN_WORD - bitPerWordPackedChar;

	while ((wordProcessed + 1) * charPerWord < textLength) {

		c = input[wordProcessed];
		for (i=0; i<charPerWord; i++) {
			buffer[i] = c >> shift;
			c <<= bitPerWordPackedChar;
		}
		wordProcessed++;

		k = 0;
		for (i=0; i<bytePerIteration; i++) {
			c = 0;
			for (j=0; j<charPerByte; j++) {
				c |= buffer[k] << (BITS_IN_BYTE - (j+1) * bitPerBytePackedChar);
				k++;
			}
			output[byteProcessed] = (unsigned char)c;
			byteProcessed++;
		}

	}

	c = input[wordProcessed];
	for (i=0; i < textLength - wordProcessed * charPerWord; i++) {
		buffer[i] = c >> shift;
		c <<= bitPerWordPackedChar;
	}

	k = 0;
	while (byteProcessed * charPerByte < textLength) {
		c = 0;
		for (j=0; j < textLength - wordProcessed * charPerWord; j++) {
			c |= buffer[k] << (BITS_IN_BYTE - (j+1) * bitPerBytePackedChar);
			k++;
		}
		output[byteProcessed] = (unsigned char)c;
		byteProcessed++;
	}

}

void ConvertBytePackedToWordPacked(const unsigned char *input, unsigned int *output, const unsigned int alphabetSize, const unsigned int textLength) {

	unsigned int i, j, k;
	unsigned int c;
	unsigned int bitPerBytePackedChar;
	unsigned int bitPerWordPackedChar;
	unsigned int charPerWord;
	unsigned int charPerByte;
	unsigned int bytePerIteration;
	unsigned int byteProcessed = 0;
	unsigned int wordProcessed = 0;
	unsigned int mask, shift;
	
	unsigned int buffer[BITS_IN_WORD];

	bitPerBytePackedChar = BitPerBytePackedChar(alphabetSize);
	bitPerWordPackedChar = BitPerWordPackedChar(alphabetSize);
	charPerByte = BITS_IN_BYTE / bitPerBytePackedChar;
	charPerWord = BITS_IN_WORD / bitPerWordPackedChar;

	bytePerIteration = charPerWord / charPerByte;
	mask = truncateRight(ALL_ONE_MASK, BITS_IN_WORD - bitPerWordPackedChar);
	shift = BITS_IN_WORD - BITS_IN_BYTE + bitPerBytePackedChar - bitPerWordPackedChar;

	while ((wordProcessed + 1) * charPerWord < textLength) {

		k = 0;
		for (i=0; i<bytePerIteration; i++) {
			c = (unsigned int)input[byteProcessed] << shift;
			for (j=0; j<charPerByte; j++) {
				buffer[k] = c & mask;
				c <<= bitPerBytePackedChar;
				k++;
			}
			byteProcessed++;
		}

		c = 0;
		for (i=0; i<charPerWord; i++) {
			c |= buffer[i] >> bitPerWordPackedChar * i;
		}
		output[wordProcessed] = c;
		wordProcessed++;

	}

	k = 0;
	for (i=0; i < (textLength - wordProcessed * charPerWord - 1) / charPerByte + 1; i++) {
		c = (unsigned int)input[byteProcessed] << shift;
		for (j=0; j<charPerByte; j++) {
			buffer[k] = c & mask;
			c <<= bitPerBytePackedChar;
			k++;
		}
		byteProcessed++;
	}

	c = 0;
	for (i=0; i<textLength - wordProcessed * charPerWord; i++) {
		c |= buffer[i] >> bitPerWordPackedChar * i;
	}
	output[wordProcessed] = c;

}

void ConvertTextToCode(const unsigned char *input, unsigned char *output, const unsigned char *charMap, const unsigned int textLength) {

	unsigned int i;

	for (i=0; i< textLength; i++) {
		output[i] = charMap[input[i]];
	}

}

void ConvertCodeToText(const unsigned char *input, unsigned char *output, const unsigned char *reverseCharMap, const unsigned int textLength) {

	unsigned int i;

	for (i=0; i< textLength; i++) {
		output[i] = reverseCharMap[input[i]];
	}

}

void PackTextWithAllShift(const unsigned char *input, unsigned int **output, const unsigned char *charMap, const unsigned int alphabetSize, const unsigned int textLength) {

	unsigned int bitPerChar;
	unsigned int numberOfShift;
	unsigned int numberOfWord;
	unsigned int shift;

	unsigned int i, j;

	bitPerChar = BitPerWordPackedChar(alphabetSize);
	numberOfShift = BITS_IN_WORD / bitPerChar;
	numberOfWord = WordPackedLengthFromText(textLength, bitPerChar);

	ConvertTextToWordPacked(input, output[0], charMap, alphabetSize, textLength);

	for (i=1; i<numberOfShift; i++) {
		shift = i * bitPerChar;
		output[i][0] = output[0][0] >> shift;
		for (j=1; j<=numberOfWord; j++) {
			output[i][j] = (output[0][j] >> shift) | (output[0][j-1] << (BITS_IN_WORD - shift));
		}
	}

}


unsigned int ReadTextAsWordPacked(const char *inputFileName, const unsigned char *charMap, const unsigned int alphabetSize, unsigned int *targetAddress, const unsigned int maxTextLength) {

	FILE *inputFile;
	unsigned char *buffer;
	unsigned int charPerWord;
	unsigned int charRead;
	unsigned int charProcessed = 0, wordProcessed = 0;
	unsigned int charPerBuffer;

	inputFile = (FILE*)fopen64(inputFileName, "rb");

	if (inputFile == NULL) {
		fprintf(stderr, "ReadTextAsWordPacked() : Cannot open inputFileName!\n");
		exit(1);
	}

	charPerWord = BITS_IN_WORD / BitPerWordPackedChar(alphabetSize);
	charPerBuffer = PACKED_BUFFER_SIZE / charPerWord * charPerWord;

	buffer = MMUnitAllocate(charPerBuffer);

	charRead = (unsigned int)fread(buffer, 1, charPerBuffer, inputFile);
	while (charRead > 0 && charProcessed + charRead < maxTextLength) {
		ConvertTextToWordPacked(buffer, targetAddress + wordProcessed, charMap, alphabetSize, charRead);
		wordProcessed += charRead / charPerWord;
		charProcessed += charRead;
		charRead = (unsigned int)fread(buffer, 1, charPerBuffer, inputFile);
	}

	if (charRead > 0 && charProcessed < maxTextLength) {
		ConvertTextToWordPacked(buffer, targetAddress + wordProcessed, charMap, alphabetSize, minX(charRead, maxTextLength - charProcessed));
		charProcessed += charRead;
	}

	MMUnitFree(buffer, charPerBuffer);

	fclose(inputFile);

	return charProcessed;

}

unsigned int ReadBytePackedAsWordPacked(const char *inputFileName, const unsigned int alphabetSize, unsigned int *targetAddress, const unsigned int maxTextLength) {

	FILE *inputFile;
	unsigned char *buffer1, *buffer2;
	unsigned int charPerByte, charPerWord;
	unsigned int charPerBuffer, wordPerBuffer;
	unsigned int charProcessed = 0, wordProcessed = 0;
	unsigned int byteRead, tempByteRead;
	unsigned int charInLastBuffer;
	unsigned int bufferSize;

	inputFile = (FILE*)fopen64(inputFileName, "rb");

	if (inputFile == NULL) {
		fprintf(stderr, "ReadBytePackedAsWordPacked() : Cannot open inputFileName!\n");
		exit(1);
	}

	charPerByte = BITS_IN_BYTE / BitPerBytePackedChar(alphabetSize);
	charPerWord = BITS_IN_WORD / BitPerWordPackedChar(alphabetSize);
	bufferSize = PACKED_BUFFER_SIZE / charPerByte / charPerWord * charPerByte * charPerWord;

	charPerBuffer = bufferSize * charPerByte;
	wordPerBuffer = charPerBuffer / charPerWord;

	buffer1 = MMUnitAllocate(bufferSize);
	buffer2 = MMUnitAllocate(bufferSize);

	byteRead = (unsigned int)fread(buffer1, 1, bufferSize, inputFile);
	tempByteRead = (unsigned int)fread(buffer2, 1, bufferSize, inputFile);

	while (tempByteRead > 1 && charProcessed + charPerBuffer < maxTextLength) {
		ConvertBytePackedToWordPacked(buffer1, targetAddress + wordProcessed, alphabetSize, charPerBuffer);
		charProcessed += charPerBuffer;
		wordProcessed += wordPerBuffer;
		memcpy(buffer1, buffer2, bufferSize);
		byteRead = tempByteRead;
		tempByteRead = (unsigned int)fread(buffer2, 1, bufferSize, inputFile);
	}

	if (tempByteRead > 1) {
		ConvertBytePackedToWordPacked(buffer1, targetAddress + wordProcessed, alphabetSize, maxTextLength - charProcessed);
		charProcessed += charPerBuffer;
	} else {
		if (tempByteRead == 1) {
			charInLastBuffer = charPerBuffer - charPerByte + buffer2[0];
		} else {
			charInLastBuffer = (byteRead - 2) * charPerByte + buffer1[byteRead - 1];
		}
		ConvertBytePackedToWordPacked(buffer1, targetAddress + wordProcessed, alphabetSize, minX(maxTextLength - charProcessed, charInLastBuffer));
		charProcessed += charInLastBuffer;
	}

	MMUnitFree(buffer1, bufferSize);
	MMUnitFree(buffer2, bufferSize);

	fclose(inputFile);

	return charProcessed;

}

// Alphabet size of DNA must be 4
void *DNALoadPacked(const char *inputFileName, unsigned int *textLength, const unsigned int convertToWordPacked) {

	FILE *inputFile;
	unsigned char tempChar[4];
	unsigned int *packedText;
	unsigned int packedFileLen;
	unsigned char lastByteLength;
	unsigned int wordToProcess;
	unsigned int i;

	inputFile = (FILE*)(FILE*)fopen64(inputFileName, "rb");

	if (inputFile == NULL) {
		fprintf(stderr, "DNALoadPacked() : Cannot open inputFileName!\n");
		exit(1);
	}

	fseek(inputFile, -1, SEEK_END);
	packedFileLen = ftell(inputFile);
	if ((int)packedFileLen < 0) {
		fprintf(stderr, "DNALoadPacked(): Cannot determine file length!\n");
		exit(1);
	}
	fread(&lastByteLength, sizeof(unsigned char), 1, inputFile);

	*textLength = (packedFileLen - 1) * 4 + lastByteLength;

	wordToProcess = (*textLength + 16 - 1) / 16;
	packedText = MMUnitAllocate((wordToProcess + 1) * sizeof(unsigned int));	// allocate 1 more word at end
	packedText[wordToProcess - 1] = 0;
	packedText[wordToProcess] = 0;

	fseek(inputFile, 0, SEEK_SET);
	fread(packedText, 1, packedFileLen, inputFile);
	fclose(inputFile);

	if (convertToWordPacked) {

		for (i=0; i<wordToProcess; i++) {
	
			*(unsigned int*)tempChar = packedText[i];
			packedText[i] = (tempChar[0] << 24) | (tempChar[1] << 16) | (tempChar[2] << 8) | tempChar[3];

		}

	}

	return (void*)packedText;

}

void DNAFreePacked(void* packedDNA, const unsigned int textLength) {

	MMUnitFree(packedDNA, ((textLength + 16 - 1) / 16 + 1) * sizeof(unsigned int));

}

void SaveText(const char *outputFileName, const unsigned char *text, const unsigned int textLength) {

	FILE *outputFile;

	outputFile = (FILE*)fopen64(outputFileName, "wb");

	if (outputFile == NULL) {
		fprintf(stderr, "SaveText() : Cannot open output file!\n");
		exit(1);
	}

	fwrite(text, sizeof(unsigned char), textLength, outputFile);
	fclose(outputFile);

}

void SaveBytePacked(const char *outputFileName, const unsigned char *bytePacked, const unsigned int textLength, const unsigned int alphabetSize) {

	FILE *outputFile;
	unsigned int bitPerChar, charPerByte, bytePackedLen;
	unsigned char lastByteLen;
	unsigned char zero = 0;

	outputFile = (FILE*)fopen64(outputFileName, "wb");

	if (outputFile == NULL) {
		fprintf(stderr, "SaveBytePacked() : Cannot open output file!\n");
		exit(1);
	}

	bitPerChar = BitPerBytePackedChar(alphabetSize);
	charPerByte = BITS_IN_BYTE / bitPerChar;

	bytePackedLen = BytePackedLengthFromText(textLength, bitPerChar);
	lastByteLen = LastByteLength(textLength, bitPerChar);

	fwrite(bytePacked, sizeof(unsigned char), bytePackedLen, outputFile);
	if (lastByteLen == 0) {
		fwrite(&zero, sizeof(unsigned char), 1, outputFile);
	}
	fwrite(&lastByteLen, sizeof(unsigned char), 1, outputFile);
	fclose(outputFile);

}

void SaveWordPacked(const char *outputFileName, const unsigned int *wordPacked, const unsigned int textLength, const unsigned int alphabetSize) {

	FILE *outputFile;
	unsigned int bitPerChar, charPerWord, wordPackedLen;
	unsigned int lastWordLen;
	unsigned int zero = 0;

	outputFile = (FILE*)fopen64(outputFileName, "wb");

	if (outputFile == NULL) {
		fprintf(stderr, "SaveWordPacked() : Cannot open output file!\n");
		exit(1);
	}

	bitPerChar = BitPerWordPackedChar(alphabetSize);
	charPerWord = BITS_IN_WORD / bitPerChar;

	wordPackedLen = WordPackedLengthFromText(textLength, bitPerChar);
	lastWordLen = LastWordLength(textLength, bitPerChar);

	fwrite(wordPacked, sizeof(unsigned int), wordPackedLen, outputFile);
	if (lastWordLen == 0) {
		fwrite(&zero, sizeof(unsigned int), 1, outputFile);
	}
	fwrite(&lastWordLen, sizeof(unsigned int), 1, outputFile);
	fclose(outputFile);

}

FILE *InitialLoadPackedIncFromEnd(const char* inputFileName, unsigned char *packedOutput, const unsigned int alphabetSize, 
								  const unsigned int packedLengthPerLoad, unsigned int *textLength, unsigned int *textLengthForThisLoad) {

	FILE *packedFile;
	unsigned int len, packedFileLenForThisLoad, packedFileLen;
	unsigned char lastByteLength;
	unsigned int bitPerChar, charPerWord;

	packedFile = (FILE*)fopen64(inputFileName, "rb");

	if (packedFile == NULL) {
		fprintf(stderr, "InitialLoadPackedIncFromEnd() : Cannot open inputFileName!\n");
		exit(1);
	}

	bitPerChar = BitPerBytePackedChar(alphabetSize);
	charPerWord = BITS_IN_WORD / bitPerChar;

	fseek(packedFile, -1, SEEK_END);
	packedFileLen = ftell(packedFile);
	if ((int)packedFileLen < 0) {
		fprintf(stderr, "InitialLoadPackedIncFromEnd(): Cannot determine file length!\n");
		exit(1);
	}
	fread(&lastByteLength, sizeof(unsigned char), 1, packedFile);

	len = TextLengthFromBytePacked(packedFileLen, bitPerChar, lastByteLength);

	if (lastByteLength == 0 && (packedFileLen - 1) % packedLengthPerLoad == 0) {
		packedFileLenForThisLoad = 0;
		fseek(packedFile, -((int)(2+packedLengthPerLoad)), SEEK_END);
		*textLength = len;
		*textLengthForThisLoad = 0;
		return packedFile;
	}

	if (packedFileLen % packedLengthPerLoad == 0) {
		packedFileLenForThisLoad = packedLengthPerLoad;
	} else {
		packedFileLenForThisLoad = packedFileLen % packedLengthPerLoad;
	}
	fseek(packedFile, -1, SEEK_END);

	fseek(packedFile, -((int)packedFileLenForThisLoad), SEEK_CUR);
	fread(packedOutput, sizeof(unsigned char), packedFileLenForThisLoad, packedFile);
	fseek(packedFile, -((int)packedFileLenForThisLoad), SEEK_CUR);
	if (packedFileLen > packedFileLenForThisLoad) {
		fseek(packedFile, -((int)packedLengthPerLoad), SEEK_CUR);
	}

	*textLength = len;
	*textLengthForThisLoad = TextLengthFromBytePacked(packedFileLenForThisLoad, bitPerChar, lastByteLength);

	return packedFile;

}

void LoadPackedIncFromEnd(FILE *packedFile, unsigned char *packedOutput, const unsigned int packedLengthPerLoad) {
	
	fread(packedOutput, sizeof(unsigned char), packedLengthPerLoad, packedFile);
	fseek(packedFile, -(2*(int)packedLengthPerLoad), SEEK_CUR);

}


FILE *InitialLoadTextIncFromEnd(const char* inputFileName, unsigned char *textOutput, const unsigned int textLengthPerLoad, unsigned int *textLength, unsigned int *textLengthForThisLoad) {

	FILE *textFile;
	unsigned int len, textLenForThisLoad;

	textFile = (FILE*)fopen64(inputFileName, "rb");

	if (textFile == NULL) {
		fprintf(stderr, "InitialLoadTextIncFromEnd() : Cannot open inputFileName!\n");
		exit(1);
	}

	fseek(textFile, 0, SEEK_END);
	len = ftell(textFile);
	if ((int)len < 0) {
		fprintf(stderr, "InitialLoadTextIncFromEnd(): Cannot determine file length!\n");
		exit(1);
	}

	textLenForThisLoad = len % textLengthPerLoad;

	if (textLenForThisLoad > 0) {
		fseek(textFile, -((int)textLenForThisLoad), SEEK_END);
		fread(textOutput, sizeof(unsigned char), textLenForThisLoad, textFile);
		fseek(textFile, -((int)textLenForThisLoad), SEEK_END);
	}

	*textLength = len;
	*textLengthForThisLoad = textLenForThisLoad;

	return textFile;
}

void LoadTextIncFromEnd(FILE *textFile, unsigned char *textOutput, const unsigned int textLengthPerLoad) {

	if (ftell(textFile) < (int)textLengthPerLoad) {
		fprintf(stderr, "LoadTextIncFromEnd(): file pointer is not correctly placed!\n");
		exit(1);
	}

	fseek(textFile, -((int)textLengthPerLoad), SEEK_CUR);
	fread(textOutput, sizeof(unsigned char), textLengthPerLoad, textFile);
	fseek(textFile, -((int)textLengthPerLoad), SEEK_CUR);

}
#endif
