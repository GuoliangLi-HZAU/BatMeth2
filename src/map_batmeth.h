#ifndef MAP_H
#define MAP_H
#include "common.h"
unsigned Map(BWT* fwfmi,BWT* revfmi,FILE* Output_File,LEN & L, BATPARAMETERS P, MEMLOOK MLook,INFILE & F);
void *Map_t(void* TP);



//copy arguments to thread arg to be passed to threading ...
#define Thread_Arg_Copy(Tx,fwfmix,revfmix,Fx,Px,Output_Filex,Lx,MLookx) \
	{	\
		Tx.fwfmi=fwfmix;\
		Tx.revfmi=revfmix;\
		Tx.In=Fx;\
		Tx.Bat_Para=Px;\
		Tx.Output_File=Output_Filex;\
		Tx.L=Lx; \
		Tx.MLook=MLookx;\
	}

#define Init_Decode(Ox,Px,Genome_Offsetsx,Location_Arrayx,Lx,Fx,LOCATIONFILEx)\
{	\
	Offset_Record Genome_Offsetsx[80];\
	unsigned Location_Arrayx[80];\
	Ox.Length_Array[1]=Lx.STRINGLENGTH_ORG;Ox.Length_Array[2]=0;\
	Ox.SAM=FALSE;Ox.PLUSSTRAND=FALSE;Ox.MaxHits=P.MAXHITS;Ox.Offset=0;Ox.Genome_Offsets=Genome_Offsetsx;Ox.Location_Array=Location_Arrayx;Ox.FILETYPE=F.FILETYPE;\
	Ox.Genome_Count= Load_Location(LOCATIONFILEx,Genome_Offsetsx,Location_Arrayx);\
}

#endif
