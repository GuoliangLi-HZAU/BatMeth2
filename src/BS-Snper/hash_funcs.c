#include "hash_funcs.h"

///////////////////////////////////////////////////////////////////////////////////////////////
// Basic hash functions
///////////////////////////////////////////////////////////////////////////////////////////////
// {{
// Initialize hash table
void hash_table_init(HashNode** hashTable, int* hash_table_size)
{
    *hash_table_size = 0;
    memset(hashTable, 0, sizeof(HashNode*) * HASH_TABLE_MAX_SIZE);
}

// String hash function
unsigned int hash_table_hash_str(const char* skey)
{
    const signed char *p = (const signed char*)skey;
    unsigned int h = *p;
    if(h) {
		for(p += 1; *p != '\0'; ++p)
			h = (h << 5) - h + *p;
    }
    return h;
}

// Insert key-value into hash table
void hash_table_insert(HashNode** hashTable, int* hash_table_size, const char* skey, int nvalue)
{
    if(*hash_table_size >= HASH_TABLE_MAX_SIZE) {
        printf("out of hash table memory!\n");
        exit(1);
    }

    unsigned int pos = hash_table_hash_str(skey) % HASH_TABLE_MAX_SIZE;
    HashNode* pHead =  hashTable[pos];
    while(pHead) {
        if(strcmp(pHead->sKey, skey) == 0) {
            printf("%s already exists!\n", skey);
            exit(1);
        }
        pHead = pHead->pNext;
    }

    HashNode* pNewNode = (HashNode*)malloc(sizeof(HashNode));
    memset(pNewNode, 0, sizeof(HashNode));
    pNewNode->sKey = (char*)malloc(sizeof(char) * (strlen(skey) + 1));
    strcpy(pNewNode->sKey, skey);
    pNewNode->nValue = nvalue;

    pNewNode->pNext = hashTable[pos];
    hashTable[pos] = pNewNode;

    *hash_table_size++;
}

// Remove key-value frome the hash table
void hash_table_remove(HashNode** hashTable, int* hash_table_size, const char* skey)
{
    unsigned int pos = hash_table_hash_str(skey) % HASH_TABLE_MAX_SIZE;
    if(hashTable[pos]) {
        HashNode* pHead = hashTable[pos];
        HashNode* pLast = NULL;
        HashNode* pRemove = NULL;
        while(pHead) {
            if(strcmp(skey, pHead->sKey) == 0) {
                pRemove = pHead;
                break;
            }
            pLast = pHead;
            pHead = pHead->pNext;
        }
        if(pRemove) {
            if(pLast)
                pLast->pNext = pRemove->pNext;
            else
                hashTable[pos] = NULL;

            free(pRemove->sKey);
            free(pRemove);
			
			*hash_table_size--;
        }
    }
}

// Lookup a key in the hash table
int hash_table_lookup(HashNode** hashTable, const char* skey)
{
    unsigned int pos = hash_table_hash_str(skey) % HASH_TABLE_MAX_SIZE;
    if(hashTable[pos]) {
        HashNode* pHead = hashTable[pos];
        while(pHead) {
            if(strcmp(skey, pHead->sKey) == 0)
                return pHead->nValue;
            pHead = pHead->pNext;
        }
    }
    return -1;
}

// Print the content in the hash table
void hash_table_print(HashNode** hashTable)
{
    printf("===========Content of hash table=================\n");
    int i;
    for(i = 0; i < HASH_TABLE_MAX_SIZE; ++i) {
        if(hashTable[i]) {
            HashNode* pHead = hashTable[i];
            printf("%d=>", i);
            while(pHead) {
                printf("%s:%d  ", pHead->sKey, pHead->nValue);
                pHead = pHead->pNext;
            }
            printf("\n");
        }
	}
}

// Free the memory of the hash table
void hash_table_release(HashNode** hashTable)
{
    int i;
    for(i = 0; i < HASH_TABLE_MAX_SIZE; ++i) {
        if(hashTable[i]) {
            HashNode* pHead = hashTable[i];
            while(pHead) {
                HashNode* pTemp = pHead;
                pHead = pHead->pNext;
                if(pTemp) {
                    free(pTemp->sKey);
                    free(pTemp);
                }
            }
        }
    }
}
// }}

