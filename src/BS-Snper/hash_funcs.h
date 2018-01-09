#ifndef _HASH_FUNCS_H
#define _HASH_FUNCS_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define HASH_TABLE_MAX_SIZE 10000

typedef struct HashNode_Struct HashNode;
struct HashNode_Struct
{
	char* sKey;
    int nValue;
    HashNode* pNext;
};

void hash_table_init(HashNode** hashTable, int* hash_table_size);
unsigned int hash_table_hash_str(const char* skey);
void hash_table_insert(HashNode** hashTable, int* hash_table_size, const char* skey, int nvalue);
void hash_table_remove(HashNode** hashTable, int* hash_table_size, const char* skey);
int hash_table_lookup(HashNode** hashTable, const char* skey);
void hash_table_print(HashNode** hashTable);
void hash_table_release(HashNode** hashTable);

#endif
