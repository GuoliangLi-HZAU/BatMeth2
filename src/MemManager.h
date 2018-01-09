#include "config.h"
#ifdef MMX
/*

   MemManager.h		Memory Manager

   This module provides memory management functions.

   Copyright (C) 2004, Wong Chi Kwong.

   This program is FREEALIGN software; you can redistribute it and/or
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

#ifndef __MEM_MANAGER_H__
#define __MEM_MANAGER_H__

#include "TypeNLimit.h"

#define MAX_ALIGN	64	// All memory except pool memory are aligned to MAX_ALIGN; pool memory is aligned to finer boundary for small memory size
#define MIN_ALIGN	1

#define RECORD_GRAND_TOTAL

//	Memory type:
//
//		unit memory:	allocation managed by malloc() individually;
//						to be used for large and less frequently accessed items
//						allocation can be freed individually at any time
//		pool memory:	pre-allocated memory pool for items with varying sizes
//						allocation cannot be freed by individually
//						to be used for small and frequently accessed items
//		temp memory:	temporary use granted from pool memory
//						allocation is allocated and freed like the items in a stack
//						pool memory allocation is disabled while temporary memory is in use
//		bulk memory:	pre-allocated memory pool for items with the same size
//						to be used for massively numbered items
//						memory address of dispatched items can be calculated by dispatch index


#ifdef DEBUG
#define Mem(mmBulk, index)  MMBulkAddress(mmBulk, index)
#else
#define Mem(mmBulk, index) 	(void*)&(mmBulk->directory[index >> mmBulk->itemPerAllocationInPowerOf2][(index & mmBulk->indexMask) * mmBulk->itemSize])
#endif

typedef struct MMPool {
	unsigned int poolSize;						// Size of memory pool; the beginning of the pool holds the MMPool structure
	unsigned int poolByteDispatched;			// Includes any spillover and memory skipped for align
	unsigned int poolByteSpillover;				// Exclude spillover pointers
	unsigned int currentTempByteDispatched;		// Includes any spillover
	unsigned int currentTempByteSpillover;		// Exclude spillover pointers
	unsigned int maxTotalByteDispatched;		// The max of pool memory + temp memory dispatched
	void *firstSpillOverAddress;				// if pool is freed, = address of mmPool
} MMPool;


typedef struct MMBulk {
	unsigned int itemSize;
	unsigned int itemPerAllocationInPowerOf2;
	unsigned int boundaryCushionSize;			// boundary cushion is a piece of memory allocated so that the memory around items can be safely referenced
	unsigned int indexMask;
	unsigned int currentDirectoryEntry;
	unsigned int nextUnusedItem;
	unsigned int directorySize;
	unsigned char **directory;			// if bulk is freed, = NULL
} MMBulk;

typedef struct MMMaster {
	size_t currentUnitByteAllocated;
	size_t maxUnitByteAllocated;
	unsigned int maxNumberOfPools;
	MMPool **mmPool;
	unsigned int maxNumberOfBulks;
	MMBulk **mmBulk;
	size_t maxTotalByteAllocated;
	size_t maxTotalByteDispatched;
	int traceUnitByteAllocation;
	FILE *unitByteTraceFile;
} MMMaster;

void *MMMalloc(const size_t memSize);
void MMFree(void *address);
void MMMasterInitialize(const unsigned int maxNumberOfPools, const unsigned int maxNumberOfBulks,
						const int traceUnitByteAllocation, FILE *unitByteTraceFile);
void MMMasterFreeAll();
size_t MMMasterCurrentTotalByteAllocated();
size_t MMMasterCurrentTotalByteDispatched();
size_t MMMasterMaxTotalByteAllocated();
size_t MMMasterMaxTotalByteDispatched();
void MMMasterSetMaxTotalByteAllocated();
void MMMasterSetMaxTotalByteDispatched();
void MMMasterPrintReport(FILE *output, const unsigned int withUnitDetails, const unsigned int withPoolDetails, const unsigned int withBulkDetails);

void *MMUnitAllocate(const size_t memSize);
void *MMUnitReallocate(void *address, const size_t newMemSize, const size_t oldMemSize);
void MMUnitFree(void *address, const size_t memSize);
size_t MMUnitCurrentByteAllocated();
size_t MMUnitMaxByteAllocated();
void MMUnitPrintReport(FILE *output);

MMPool *MMPoolCreate(const unsigned int poolSize);
unsigned int MMPoolIsActive(const MMPool *mmPool);
void MMPoolSetInactive(MMPool *mmPool);
unsigned int MMPoolCurrentTotalByteAllocated(const MMPool *mmPool);
unsigned int MMPoolCurrentTotalByteDispatched(const MMPool *mmPool);
unsigned int MMPoolMaxTotalByteDispatched(const MMPool *mmPool);
unsigned int MMPoolByteAvailable(const MMPool *mmPool);
MMPool *MMPoolFree(MMPool *mmPool);
void MMPoolReset(MMPool *mmPool);
void MMPoolDestory(MMPool *mmPool);
void *MMPoolDispatch(MMPool *mmPool, const unsigned int memSize);
unsigned int MMPoolDispatchOffset(MMPool *mmPool, const unsigned int memSize);
void MMPoolReturn(MMPool *mmPool, void *address, const unsigned int memSize);		// Dummy function
void MMPoolPrintReport(MMPool *mmPool, FILE *output);

void *MMTempDispatch(MMPool *mmPool, const unsigned int memsize);
void MMTempReturn(MMPool *mmPool, void *address, const unsigned int memSize);
void MMTempPrintReport(MMPool *mmPool, FILE *output);

MMBulk *MMBulkCreate(MMPool *mmPool, const unsigned int itemSize, const unsigned int itemPerAllocationInPowerOf2, 
					 unsigned int const boundaryCushionSize, unsigned int const directorySize);
unsigned int MMBulkIsActive(const MMBulk *mmBulk);
void MMBulkSetInactive(MMBulk *mmBulk);
unsigned int MMBulkByteAllocated(const MMBulk *mmBulk);
unsigned int MMBulkByteDispatched(const MMBulk *mmBulk);
unsigned int MMBulkUnitDispatched(const MMBulk *mmBulk);
void MMBulkFree(MMBulk *mmBulk);
void MMBulkDestory(MMBulk *mmBulk);
unsigned int MMBulkDispatch(MMBulk *mmBulk);
void *MMBulkAddress(const MMBulk *mmBulk, const unsigned int index);
MMPool *MMBulkFindPoolUsed(const MMBulk *mmBulk);
void MMBulkPrintReport(MMBulk *mmBulk, FILE *output);

void MMBulkSave(MMBulk *mmBulk, FILE *output);
MMBulk *MMBulkLoad(MMPool *mmPool, FILE *input);


#endif
#else

/*

   MemManager.h		Memory Manager

   This module provides memory management functions.

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

#ifndef __MEM_MANAGER_H__
#define __MEM_MANAGER_H__

#include "TypeNLimit.h"

#define MAX_ALIGN	16
#define MIN_ALIGN	1

#define RECORD_GRAND_TOTAL

//	Memory type:
//
//		unit memory:	allocation managed by malloc() individually;
//						to be used for large and less frequently accessed items
//						allocation can be freed individually at any time
//		pool memory:	pre-allocated memory pool for items with varying sizes
//						allocation cannot be freed by individually
//						to be used for small and frequently accessed items
//		temp memory:	temporary use granted from pool memory
//						allocation is allocated and freed like the items in a stack
//						pool memory allocation is disabled while temporary memory is in use
//		bulk memory:	pre-allocated memory pool for items with the same size
//						to be used for massively numbered items
//						memory address of dispatched items can be calculated by dispatch index


#ifdef DEBUG
#define Mem(mmBulk, index)  MMBulkAddress(mmBulk, index)
#else
#define Mem(mmBulk, index) 	(void*)&(mmBulk->directory[index >> mmBulk->itemPerAllocationInPowerOf2][(index & mmBulk->indexMask) * mmBulk->itemSize])
#endif

typedef struct MMPool {
	unsigned int poolSize;
	unsigned int poolByteDispatched;				// Includes any spillover and memory skipped for align
	unsigned int poolByteSkippedForAlign;
	unsigned int poolByteSpillover;				// Exclude spillover pointers
	void *firstSpillOverAddress;		// if pool is freed, = address of mmPool
	unsigned int currentTempByteDispatched;		// Includes any spillover
	unsigned int currentTempByteSpillover;		// Exclude spillover pointers
	unsigned int maxTotalByteDispatched;			// The max of pool memory + temp memory dispatched
} MMPool;


typedef struct MMBulk {
	unsigned int itemSize;
	unsigned int itemPerAllocationInPowerOf2;
	unsigned int boundaryCushionSize;			// boundary cushion is a piece of memory allocated so that the memory around items can be safely referenced
	unsigned int indexMask;
	unsigned int currentDirectoryEntry;
	unsigned int nextUnusedItem;
	unsigned int directorySize;
	unsigned char **directory;			// if bulk is freed, = NULL
} MMBulk;

typedef struct MMMaster {
	unsigned int currentUnitByteAllocated;
	unsigned int maxUnitByteAllocated;
	unsigned int maxNumberOfPools;
	MMPool **mmPool;
	unsigned int maxNumberOfBulks;
	MMBulk **mmBulk;
	unsigned int maxTotalByteAllocated;
	unsigned int maxTotalByteDispatched;
	int traceUnitByteAllocation;
	FILE *unitByteTraceFile;
} MMMaster;

void *MMMalloc(const unsigned int memSize);
void MMFree(void *address);
void MMMasterInitialize(const unsigned int maxNumberOfPools, const unsigned int maxNumberOfBulks,
						const int traceUnitByteAllocation, FILE *unitByteTraceFile);
void MMMasterFreeAll();
unsigned int MMMasterCurrentTotalByteAllocated();
unsigned int MMMasterCurrentTotalByteDispatched();
unsigned int MMMasterMaxTotalByteAllocated();
unsigned int MMMasterMaxTotalByteDispatched();
void MMMasterSetMaxTotalByteAllocated();
void MMMasterSetMaxTotalByteDispatched();
void MMMasterPrintReport(FILE *output, const unsigned int withUnitDetails, const unsigned int withPoolDetails, const unsigned int withBulkDetails);

void *MMUnitAllocate(const unsigned int memSize);
void *MMUnitReallocate(void *address, const unsigned int newMemSize, const unsigned int oldMemSize);
void MMUnitFree(void *address, const unsigned int memSize);
unsigned int MMUnitCurrentByteAllocated();
unsigned int MMUnitMaxByteAllocated();
void MMUnitPrintReport(FILE *output);

MMPool *MMPoolCreate(const unsigned int poolSize);
unsigned int MMPoolIsActive(const MMPool *mmPool);
void MMPoolSetInactive(MMPool *mmPool);
unsigned int MMPoolCurrentTotalByteAllocated(const MMPool *mmPool);
unsigned int MMPoolCurrentTotalByteDispatched(const MMPool *mmPool);
unsigned int MMPoolMaxTotalByteDispatched(const MMPool *mmPool);
unsigned int MMPoolByteAvailable(const MMPool *mmPool);
MMPool *MMPoolFree(MMPool *mmPool);
void MMPoolReset(MMPool *mmPool);
void MMPoolDestory(MMPool *mmPool);
void *MMPoolDispatch(MMPool *mmPool, const unsigned int memSize);
void MMPoolReturn(MMPool *mmPool, void *address, const unsigned int memSize);		// Dummy function
void MMPoolPrintReport(MMPool *mmPool, FILE *output);

void *MMTempDispatch(MMPool *mmPool, const unsigned int memsize);
void MMTempReturn(MMPool *mmPool, void *address, const unsigned int memSize);
void MMTempPrintReport(MMPool *mmPool, FILE *output);

MMBulk *MMBulkCreate(MMPool *mmPool, const unsigned int itemSize, const unsigned int itemPerAllocationInPowerOf2, 
					 unsigned int const boundaryCushionSize, unsigned int const directorySize);
unsigned int MMBulkIsActive(const MMBulk *mmBulk);
void MMBulkSetInactive(MMBulk *mmBulk);
unsigned int MMBulkByteAllocated(const MMBulk *mmBulk);
unsigned int MMBulkByteDispatched(const MMBulk *mmBulk);
unsigned int MMBulkUnitDispatched(const MMBulk *mmBulk);
void MMBulkFree(MMBulk *mmBulk);
void MMBulkDestory(MMBulk *mmBulk);
unsigned int MMBulkDispatch(MMBulk *mmBulk);
void *MMBulkAddress(const MMBulk *mmBulk, const unsigned int index);
MMPool *MMBulkFindPoolUsed(const MMBulk *mmBulk);
void MMBulkPrintReport(MMBulk *mmBulk, FILE *output);

void MMBulkSave(MMBulk *mmBulk, FILE *output);
MMBulk *MMBulkLoad(MMPool *mmPool, FILE *input);


#endif
#endif
