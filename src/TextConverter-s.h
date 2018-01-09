#include "config.h"
#ifdef MMX
/*

   TextConverter.h		Text Converter

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

#ifndef __TEXTCONVERTOR_H__
#define __TEXTCONVERTOR_H__

#include "TypeNLimit-s.h"
#include "MemManager-s.h"

#define INVALID_CHAR 0xFF
#define CHAR_MAP_SIZE 256
#define PACKED_BUFFER_SIZE			(PACKED_BUFFER_SIZE_IN_WORD * BYTES_IN_WORD)
#define PACKED_BUFFER_SIZE_IN_WORD	65536
#define MAX_SEQ_NAME_LENGTH			256
#define RANDOM_SUBSTITUTE			'R'

// charMap is a char array of size 256. The index of the array is the input text value
// and the content of the array is the output text value. e.g. A -> 0, C -> 1
// If the value of an entry = INVALID_CHAR, the indexed text value is an invalid input

// Retrieve word packed text
unsigned int GetWordPackedText(const unsigned int *packedText, const unsigned int index, const unsigned int shift, const unsigned int numberOfBit, const unsigned int vacantBit);

// Character map functions
unsigned int ReadCharMap(unsigned char *charMap, const char *inputFileName, const unsigned char defaultMapping);
void GenerateReverseCharMap(const unsigned char *charMap, unsigned char *reverseCharMap);

// Word packed text functions
unsigned int BitPerWordPackedChar(const unsigned int alphabetSize);
unsigned int TextLengthFromWordPacked(unsigned int wordPackedLength, unsigned int bitPerChar, unsigned int lastWordLength);
unsigned int WordPackedLengthFromText(unsigned int textLength, unsigned int bitPerChar);
unsigned int LastWordLength(unsigned int textLength, unsigned int bitPerChar);

// Byte packed text functions
unsigned int BitPerBytePackedChar(const unsigned int alphabetSize);
unsigned int TextLengthFromBytePacked(unsigned int bytePackedLength, unsigned int bitPerChar, unsigned int lastByteLength);
unsigned int BytePackedLengthFromText(unsigned int textLength, unsigned int bitPerChar);
unsigned char LastByteLength(unsigned int textLength, unsigned int bitPerChar);

// Conversion functions
void ConvertTextToWordPacked(const unsigned char *input, unsigned int *output, const unsigned char *charMap, const unsigned int alphabetSize, const unsigned int textLength);
void ConvertTextToBytePacked(const unsigned char *input, unsigned char *output, const unsigned char *charMap, const unsigned int alphabetSize, const unsigned int textLength);
void ConvertWordPackedToText(const unsigned int *input, unsigned char *output, const unsigned char *reverseCharMap, const unsigned int alphabetSize, const unsigned int textLength);
void ConvertBytePackedToText(const unsigned char *input, unsigned char *output, const unsigned char *reverseCharMap, const unsigned int alphabetSize, const unsigned int textLength);
void ConvertBytePackedToCode(const unsigned char *input, unsigned char *output, const unsigned int alphabetSize, const unsigned int textLength);
void ConvertWordPackedToBytePacked(const unsigned int *input, unsigned char *output, const unsigned int alphabetSize, const unsigned int textLength);
void ConvertBytePackedToWordPacked(const unsigned char *input, unsigned int *output, const unsigned int alphabetSize, const unsigned int textLength);
void ConvertTextToCode(const unsigned char *input, unsigned char *output, const unsigned char *charMap, const unsigned int textLength);
void ConvertCodeToText(const unsigned char *input, unsigned char *output, const unsigned char *reverseCharMap, const unsigned int textLength);

// Pack text with all shift
void PackTextWithAllShift(const unsigned char *input, unsigned int **output, const unsigned char *charMap, const unsigned int alphabetSize, const unsigned int textLength);

// Full load function
unsigned int ReadTextAsWordPacked(const char *inputFileName, const unsigned char *charMap, const unsigned int alphabetSize, unsigned int *targetAddress, const unsigned int maxTextLength);
unsigned int ReadBytePackedAsWordPacked(const char *inputFileName, const unsigned int alphabetSize, unsigned int *targetAddress, const unsigned int maxTextLength);
void *DNALoadPacked(const char *inputFileName, unsigned int *textLength, const unsigned int convertToWordPacked, const unsigned int trailerBufferInWord);
void DNAFreePacked(void* packedDna, const unsigned int textLength, const unsigned int trailerBufferInWord);

// Save functions
void SaveText(const char *outputFileName, const unsigned char *text, const unsigned int textLength);
void SaveBytePacked(const char *outputFileName, const unsigned char *wordPacked, const unsigned int textLength, const unsigned int alphabetSize);
void SaveWordPacked(const char *outputFileName, const unsigned int *wordPacked, const unsigned int textLength, const unsigned int alphabetSize);

// Incremental load functions (start from end of text)
FILE *InitialLoadPackedIncFromEnd(const char* inputFileName, unsigned char *packedOutput, const unsigned int alphabetSize, const unsigned int packedLengthPerLoad, unsigned int *textLength, unsigned int *textLengthForThisLoad);
void LoadPackedIncFromEnd(FILE *packedFile, unsigned char *packedOutput, const unsigned int packedLengthPerLoad);
FILE *InitialLoadTextIncFromEnd(const char* inputFileName, unsigned char *textOutput, const unsigned int textLengthPerLoad, unsigned int *textLength, unsigned int *textLengthForThisLoad);
void LoadTextIncFromEnd(FILE *textFile, unsigned char *textOutput, const unsigned int textLengthPerLoad);


#endif
#else

/*

   TextConverter.h		Text Converter

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

#ifndef __TEXTCONVERTOR_H__
#define __TEXTCONVERTOR_H__

#include "TypeNLimit-s.h"
#include "MemManager-s.h"

#define INVALID_CHAR 0xFF
#define CHAR_MAP_SIZE 256
#define PACKED_BUFFER_SIZE			(PACKED_BUFFER_SIZE_IN_WORD * BYTES_IN_WORD)
#define PACKED_BUFFER_SIZE_IN_WORD	65536
#define MAX_SEQ_NAME_LENGTH			256
#define RANDOM_SUBSTITUTE			'R'

// charMap is a char array of size 256. The index of the array is the input text value
// and the content of the array is the output text value. e.g. A -> 0, C -> 1
// If the value of an entry = INVALID_CHAR, the indexed text value is an invalid input

// Retrieve word packed text
unsigned int GetWordPackedText(const unsigned int *packedText, const unsigned int index, const unsigned int shift, const unsigned int numberOfBit, const unsigned int vacantBit);

// Character map functions
unsigned int ReadCharMap(unsigned char *charMap, const char *inputFileName, const unsigned char defaultMapping);
void GenerateReverseCharMap(const unsigned char *charMap, unsigned char *reverseCharMap);

// Word packed text functions
unsigned int BitPerWordPackedChar(const unsigned int alphabetSize);
unsigned int TextLengthFromWordPacked(unsigned int wordPackedLength, unsigned int bitPerChar, unsigned int lastWordLength);
unsigned int WordPackedLengthFromText(unsigned int textLength, unsigned int bitPerChar);
unsigned int LastWordLength(unsigned int textLength, unsigned int bitPerChar);

// Byte packed text functions
unsigned int BitPerBytePackedChar(const unsigned int alphabetSize);
unsigned int TextLengthFromBytePacked(unsigned int bytePackedLength, unsigned int bitPerChar, unsigned int lastByteLength);
unsigned int BytePackedLengthFromText(unsigned int textLength, unsigned int bitPerChar);
unsigned char LastByteLength(unsigned int textLength, unsigned int bitPerChar);

// Conversion functions
void ConvertTextToWordPacked(const unsigned char *input, unsigned int *output, const unsigned char *charMap, const unsigned int alphabetSize, const unsigned int textLength);
void ConvertTextToBytePacked(const unsigned char *input, unsigned char *output, const unsigned char *charMap, const unsigned int alphabetSize, const unsigned int textLength);
void ConvertWordPackedToText(const unsigned int *input, unsigned char *output, const unsigned char *reverseCharMap, const unsigned int alphabetSize, const unsigned int textLength);
void ConvertBytePackedToText(const unsigned char *input, unsigned char *output, const unsigned char *reverseCharMap, const unsigned int alphabetSize, const unsigned int textLength);
void ConvertBytePackedToCode(const unsigned char *input, unsigned char *output, const unsigned int alphabetSize, const unsigned int textLength);
void ConvertWordPackedToBytePacked(const unsigned int *input, unsigned char *output, const unsigned int alphabetSize, const unsigned int textLength);
void ConvertBytePackedToWordPacked(const unsigned char *input, unsigned int *output, const unsigned int alphabetSize, const unsigned int textLength);
void ConvertTextToCode(const unsigned char *input, unsigned char *output, const unsigned char *charMap, const unsigned int textLength);
void ConvertCodeToText(const unsigned char *input, unsigned char *output, const unsigned char *reverseCharMap, const unsigned int textLength);

// Pack text with all shift
void PackTextWithAllShift(const unsigned char *input, unsigned int **output, const unsigned char *charMap, const unsigned int alphabetSize, const unsigned int textLength);

// Full load function
unsigned int ReadTextAsWordPacked(const char *inputFileName, const unsigned char *charMap, const unsigned int alphabetSize, unsigned int *targetAddress, const unsigned int maxTextLength);
unsigned int ReadBytePackedAsWordPacked(const char *inputFileName, const unsigned int alphabetSize, unsigned int *targetAddress, const unsigned int maxTextLength);
void *DNALoadPacked(const char *inputFileName, unsigned int *textLength, const unsigned int convertToWordPacked);
void DNAFreePacked(void* packedDna, const unsigned int textLength);

// Save functions
void SaveText(const char *outputFileName, const unsigned char *text, const unsigned int textLength);
void SaveBytePacked(const char *outputFileName, const unsigned char *wordPacked, const unsigned int textLength, const unsigned int alphabetSize);
void SaveWordPacked(const char *outputFileName, const unsigned int *wordPacked, const unsigned int textLength, const unsigned int alphabetSize);

// Incremental load functions (start from end of text)
FILE *InitialLoadPackedIncFromEnd(const char* inputFileName, unsigned char *packedOutput, const unsigned int alphabetSize, const unsigned int packedLengthPerLoad, unsigned int *textLength, unsigned int *textLengthForThisLoad);
void LoadPackedIncFromEnd(FILE *packedFile, unsigned char *packedOutput, const unsigned int packedLengthPerLoad);
FILE *InitialLoadTextIncFromEnd(const char* inputFileName, unsigned char *textOutput, const unsigned int textLengthPerLoad, unsigned int *textLength, unsigned int *textLengthForThisLoad);
void LoadTextIncFromEnd(FILE *textFile, unsigned char *textOutput, const unsigned int textLengthPerLoad);


#endif
#endif
