#include "swroutines.h"
#include "assert.h"
#include <queue>
#include "filters.h"
#include "math.h"
#include "fastsw.h"
extern long File_size;
extern bool FASTDECODE;
extern int SW_STRING_BUFFER;
extern int Inter_MM;
extern int MODE;
extern int QUALITYCONVERSIONFACTOR;
extern int INDELGAP;//15 good for //17 for 8, 19=g00d for 9
//extern READ RawRead,revRawR,RawMead,revRawM;
extern int MAX_MISMATCHES;
extern char Char_To_R[255];
int Ref_clip = 8;

//extern int match;
//extern int mismatch;
//extern int gap_open;
//extern int gap_extension;


int Do_SW_Pair(char* Pattern_raw,READ & RawR,unsigned char* Original_Text,Pair* Pairs,char* Pattern,int Flank_Size,int StringLength,int & Err,int Shift,char Sign,std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> & Alignments,int Current_Score,int & Total_SW_Scans)
{
	s_align* Aln;
	Alignment A;
	unsigned Ref_Location;
	char Org_String[SW_STRING_BUFFER];
	char Org_String_Ori[SW_STRING_BUFFER];
	int HITS=0;
	Cigar_Info Cig_Info;
	s_profile* p = ssw_init((int8_t*)Pattern, StringLength/2, mata, n, 1);

	for (int i=0;Pairs[i].Head && MAX_SW > HITS;i++)
	{
		if(Total_SW_Scans>MAX_SW) 
		{
			Err++;break;
		}
		Total_SW_Scans++;
		Ref_Location=Pairs[i].Head+Shift;//StringLength/3;
		Get_Bases(Original_Text,Ref_Location,Flank_Size,Org_String);
		Get_Bases(Original_Text_Ori,Ref_Location,Flank_Size,Org_String_Ori);
		Aln=mengyao_ssw_core(Org_String,StringLength/2, Pattern,Flank_Size,0,(MODE==VERYFAST)?1:0, p);
		if(Aln->score1 >= ACC_SCORE)
		{
			A.SW_Score=Aln->score1;
			A.swscore = A.SW_Score;
			if(MODE>=VERYFAST) ssw_cigar_process(Sign,Aln->read_begin1,Pattern_raw+Aln->read_begin1,Org_String_Ori+Aln->ref_begin1,Aln,Cig_Info,Org_String+Aln->ref_begin1,Pattern+Aln->read_begin1,StringLength/2);//momo
			HITS++;
			Ref_Location-=Shift;
			//if(Aln->ref_begin1 || Aln->read_begin1) Ref_Location=0;
			A.Loc=Ref_Location;
			A.Realigned=NO;
			if(MODE==VERYFAST)
			{
				int Clip_H,Clip_T;Clip_H=Aln->read_begin1;Clip_T=0;
				if(Aln->read_end1!=StringLength/2-1) Clip_T=StringLength/2-1-Aln->read_end1;
				int SC_Penalty=(Clip_T+Clip_H)? gap_open+(Clip_T+Clip_H)*gap_extension:0;
				A.Score= -(-SC_Penalty+Aln->score1+Current_Score);
				A.Score= -A.Score;
			}
			else
			{
				A.Score= -Cig_Info.Score-Current_Score;
			}
			A.QScore=INT_MAX;
			A.Sign=Sign;
			Alignments.push(A);
		}
		align_destroy(Aln);
	}
	init_destroy(p); 
	return HITS;
}

int Do_SW(char* Pattern_raw,READ & RawR,RQINDEX & RQHALF,unsigned char* Original_Text,BWT* revfmi,SARange* Head_Hits,char* Pattern,int Flank_Size,int StringLength,int & Err,int Shift,char Sign,int Current_Score,std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> & Alignments,int & Total_SW_Scans,int & Filter)
{
	//filter=Filter;
	unsigned Ref_Location;
	unsigned S,L;
	char Org_String[SW_STRING_BUFFER];
	char Org_String_Ori[SW_STRING_BUFFER];
	int HITS=0;
	s_align* Aln;
	Alignment A;
	Cigar_Info Cig_Info;
	int Max_Offset_Size;
	Tag_Info Head_Info;

	if (!Head_Hits) return 0;
	//printf("O-G--------------------\n");
	s_profile* p = ssw_init((int8_t*)Pattern, StringLength/2, mata, n, 1);
	for (int i=0;Head_Hits[i].Start && MAX_SW > HITS;i++)
	{
	//printf("\tO-G--------------------\n");
		S=Head_Hits[i].Start;L=Head_Hits[i].End;
		assert(L>=S);
		bool Do_Fast=false;
		int Offset=0;
		if(L!=S)
		{
			Total_SW_Scans+=(L-S+1);
			if (L-S>MAX_SW) 
			{
				Err=L-S;continue;//Too many hits..
			}
			else if(FASTDECODE)
			{
				assert(L-S>0);
				if((L-S>MAXGAP)||(L-S<=SAGAP_CUTOFF))//too low hits, not cached..too many hits, not cached..
				{
					Head_Hits[i].Start=Head_Hits[i].End=revfmi->textLength-StringLength/2-BWTSaValue(revfmi,S);
					//printf("G-%u\n",Head_Hits[i].Start);
				}
				else
				{
					Load_Info(RQHALF,Head_Info,Head_Hits[i],Entries_Half);
					Head_Hits[i].Start=Head_Hits[i].End=Head_Info.First;//H2=Head_Info.Last + d;
					Max_Offset_Size=Head_Info.Gap-1;
					Do_Fast=true;
					//printf("G-%u\n",Head_Hits[i].Start);
				}
			}
			else 
			{
				Head_Hits[i].Start=Head_Hits[i].End=revfmi->textLength-StringLength/2-BWTSaValue(revfmi,S);
				//printf("O-%u\n",Head_Hits[i].Start);
			}
		}
		else
		{
			Total_SW_Scans++;
		}
		if(Total_SW_Scans>MAX_SW) 
		{
			Err=L-S;
			break;
		}

		do
		{
			int Begin=0;
			Ref_Location=Head_Hits[i].Start+Shift;

			Get_Bases(Original_Text,Ref_Location,Flank_Size,Org_String);
			Get_Bases(Original_Text_Ori,Ref_Location,Flank_Size,Org_String_Ori);
			bool Do_Smith_Waterman=true;
			/*if(Sign=='+')
			{
				char Org_StringX[SW_STRING_BUFFER+1];
				char PatternX[SW_STRING_BUFFER+1];
				Org_StringX[Flank_Size]=0;
				for(int i=0;i<Flank_Size;i++)
				{
					Org_StringX[i]="ACGT"[Org_String[i]];
				}
				for(int i=0;i<StringLength/2;i++)
				{
					PatternX[i]="ACGT"[Pattern[i]];
				}
				PatternX[StringLength/2]=0;
				A=Fast_SWX(Org_StringX,PatternX,0);
				A.QScore=INT_MAX;
				A.Sign='+';
				if(A.Score!=INT_MAX)
				{
					if(MODE==VERYFAST && A.SW_Score >Filter) 
					{
						Filter=A.SW_Score+2;
					}
					A.Score= A.SW_Score+Current_Score;
					A.Loc=Ref_Location-Shift;
					Do_Smith_Waterman=false;
					assert(A.Score>=0);
					A.Realigned=NO;
					A.Clip_T=A.Clip_H=0;
					Alignments.push(A);
				}
			}*/
			if(Do_Smith_Waterman)
			{
				Aln=mengyao_ssw_core(Org_String,StringLength/2, Pattern,Flank_Size,0,(MODE==VERYFAST)? 1:0, p);
				if(Aln->score1 >=Filter)//ACC_SCORE)
				{
					int Clip_H=Aln->read_begin1;int Clip_T=0;
					if(Aln->read_end1!=StringLength-1) Clip_T=StringLength-1-Aln->read_end1;
					if(Clip_T+Clip_H)
						A.QScore= -1;
					else 
						A.QScore=INT_MAX;

					A.SW_Score=Aln->score1;
					A.swscore=A.SW_Score;
					A.Realigned=NO;
					if(MODE==VERYFAST && Aln->score1 >Filter) 
					{
						Filter=Aln->score1+2;
					}
					if(MODE>=VERYFAST)
					{
						ssw_cigar_process(Sign,Aln->read_begin1,Pattern_raw+Aln->read_begin1,Org_String_Ori+Aln->ref_begin1,Aln,Cig_Info,Org_String+Aln->ref_begin1,Pattern+Aln->read_begin1,StringLength/2);
					}
					HITS++;
					if(Shift<0)
					{
						Ref_Location+=Aln->ref_begin1;
					}
					else Ref_Location-=Shift;
					A.Loc=Ref_Location;
					if(MODE==VERYFAST)
					{
						int ClipH,ClipT;ClipH=Aln->read_begin1;ClipT=0;
						if(Aln->read_end1!=StringLength/2-1) ClipT=StringLength/2-1-Aln->read_end1;
						int SC_Penalty=(ClipT+ClipH)? gap_open+(ClipT+ClipH)*gap_extension:0;
						A.Score= -(-SC_Penalty+Aln->score1)-Current_Score;
						A.Score= -A.Score;
					}
					else
						A.Score= -Cig_Info.Score-Current_Score;
					A.Sign=Sign;

					Alignments.push(A);
				}
				align_destroy(Aln);
			}
			S++;
			Offset++;
			if (S<=L) 
			{
				if(Do_Fast)
				{
					Head_Hits[i].Start=Head_Hits[i].End=Get_Location(revfmi,RQHALF,Head_Info,Offset,revfmi->textLength-StringLength/2);
				}
				else
				{
					Head_Hits[i].Start=Head_Hits[i].End=revfmi->textLength-StringLength/2-BWTSaValue(revfmi,S);
				}
			}
		}
		while (S<=L);

	}
	init_destroy(p); 
	return HITS;
}
//}-------------------------------------------------- SMITH WATERMAN ----------------------------------------------------------------------

void Get_Bases (unsigned char* Original_Text,unsigned Location,int StringLength,char* Org_String) 
{
	if ((long)(--Location)<0) Location=0;
	if(Location >= File_size-StringLength) Location=0;
	//if (unsigned(--Location)<0) Location=0;
	//if(int(Location) < 0) Location=0;
	assert (StringLength<SW_STRING_BUFFER);
	{
		for (int i=0;i<=StringLength;i++)
		{
			unsigned char L= (unsigned char)(Original_Text[(Location+i)/4]<< (((Location+i) % 4) * 2)) >>6;
			Org_String[i]=L;
		}
	}
}
//}---------------------------------------- Smith Waterman ------------------------------------------------------------------------------

void ssw_cigar_print(char Sign,int ReadSpan,char* rawRead,char* rawRef,s_align* a,char* Cigar,Cigar_Info & C,char* Ref,char* Pattern,int Clip_H,int Clip_T,int StringLength) { //print the cigar out
	int c = 0,Tot_Length=0;char* Cigar_Ptr=Cigar;
	bool Cig_Err=false;
	C.M=0;C.I=0;C.D=0,C.Mis=0;
	if(Clip_H) Cigar_Ptr+=sprintf(Cigar_Ptr,"%dS", Clip_H); 
	for (c = 0; c < a->cigarLen; ++c) {//rawRef 3=T 0=G
		int32_t letter = 0xf&*(a->cigar + c);
		int32_t length = (0xfffffff0&*(a->cigar + c))>>4;
		Cigar_Ptr+=sprintf(Cigar_Ptr,"%d", length);
		if (letter == 0) 
		{
			Tot_Length+=length;
			*Cigar_Ptr='M';
			for(int i=0;i<length;i++)
			{
				assert(*Ref<4 && *Ref>=0);
				assert(*Pattern<4 && *Pattern>=0);
				assert(*rawRef<4 && *rawRef>=0);
				if(*Ref!=*Pattern || (*Ref==*Pattern && *rawRef==3 && *rawRead==1  ) ||  (*Ref==*Pattern && *rawRef==0 && *rawRead==2 ) ) {C.Mis++;}
				Ref++;Pattern++,rawRef++,rawRead++;
			}
			Cigar_Ptr++;C.M+=length;
		}
		else if (letter == 1)
		{
			Tot_Length+=length;
			*Cigar_Ptr='I';Cigar_Ptr++;C.I+=length;
			Pattern+=length,rawRead+=length;
		} 
		else 
		{
			*Cigar_Ptr='D';Cigar_Ptr++;C.D+=length;
			Ref+=length,rawRef+=length;
		}
		if(Cigar_Ptr-Cigar>=MAX_SIGLEN-6) 
		{
			assert(StringLength-Tot_Length-Clip_H>=0);
			Cigar_Ptr+=sprintf(Cigar_Ptr,"%dS", StringLength-Tot_Length-Clip_H); 
			Cig_Err=true;
			break;
		}
	}
	if(Clip_T && !Cig_Err) Cigar_Ptr+=sprintf(Cigar_Ptr,"%dS", Clip_T); 
	*Cigar_Ptr=0;
	//if(Cig_Err) *Cigar=0;
	C.Length=Tot_Length;
}

void ssw_cigar_process(char Sign,int ReadSpan,char* rawRead,char* rawRef,s_align* a,Cigar_Info & C,char* Ref,char* Pattern,int StringLength) { //print the cigar out

	int gap_openX=40,gap_extensionX=6;
	if(MODE>=FAST)
		gap_openX=40,gap_extensionX=6;
	else
		gap_openX=gap_open,gap_extensionX=gap_extension;

	int c = 0,Tot_Length=0;
	C.M=0;C.I=0;C.D=0,C.Mis=0;C.Score=0;C.Indel_Count=0;
	for (c = 0; c < a->cigarLen; ++c) {
		int32_t letter = 0xf&*(a->cigar + c);
		int32_t length = (0xfffffff0&*(a->cigar + c))>>4;
		if (letter == 0) 
		{
			Tot_Length+=length;
			for(int i=0;i<length;i++)
			{
				assert(*Ref<4 && *Ref>=0  );
				assert(*Pattern<4 && *Pattern>=0);
				if(*Ref!=*Pattern || (*Ref==*Pattern && *rawRef==3 && *rawRead==1  ) ||  (*Ref==*Pattern && *rawRef==0 && *rawRead==2 )  )
				{
					//std::cout << " ";//*Ref<< *Pattern;
					C.Mis++;
					if(MODE==VERYFAST) C.D-=mismatch;
					//C.Score-=mismatch;
					else
						C.Score+=Mis_FA_Score;
				}
				else //same
				{
				//	if(*Ref_Ori!=*Pattern)
				//	{
						
			//	}
			//		else
			//		{
						//C.Score+=match;
						if(MODE==VERYFAST) C.D+=match;
						else	
							C.Score+=Match_FA_Score;
			//	}
				}
				Ref++;Pattern++,rawRef++,rawRead++;
			}
			C.M+=length;
		}
		else if (letter == 1)
		{
			Tot_Length+=length;
			C.Indel_Count++;
			C.I+=length;
			if(MODE>=FAST)
				C.Score+=(gap_openX+length*gap_extensionX);
			else
				C.D-=(gap_openX+(length-1)*gap_extensionX);
			Pattern+=length;rawRead+=length;
		} 
		else 
		{
			C.Indel_Count++;
			if(MODE>=FAST)
				C.Score+=(gap_openX+length*gap_extensionX);
			else
				C.D-=(gap_openX+(length-1)*gap_extensionX);
			//C.D+=length;
			Ref+=length;rawRef+=length;
		} 
	}
	C.Length=Tot_Length;
	assert(Tot_Length<=StringLength);
	if(Tot_Length<StringLength)
	{
		if(MODE>=FAST)
			C.Score+=gap_openX+(StringLength-Tot_Length)*gap_extensionX;
		else
			C.D-=(gap_openX+(StringLength-Tot_Length)*gap_extensionX);
	}
}


void ssw_cigar_processQ(char Sign,char* rawRead,char* rawRef,s_align* a,Cigar_Info & C,char* Ref,int Ref_Off,char* Pattern,int Pat_Off,int StringLength,const char* Qual,char* Cigar,int Clip_H,int Clip_T,bool Hard_Penalty) { //print the cigar out //,READ & R
	bool Cig_Err=false;
	const int gap_openN=40,gap_extensionN=6; //changed by qiangwei and temp
	const int MX=6,MN=2,BOPEN=6,BEXT=3,MATCH_BONUS=0; //0;//2;
	int c = 0,Tot_Length=0;char* Cigar_Ptr=Cigar;
	float Score=0,QScore=0,BQScore=0;
	float TScore=0,TQScore=0,TBQScore=0;
	int Tclipswscore = 0, Hclipswscore = 0;
	int swscore=0;
	C.M=0;C.I=0;C.D=0,C.Mis=0;C.Score=0;C.Indel_Count=0;
	//if(Cigar && Clip_H) Cigar_Ptr+=sprintf(Cigar_Ptr,"%dS", Clip_H); 
	Ref+=Ref_Off;Pattern+=Pat_Off;Qual+=Pat_Off;
	rawRef+=Ref_Off;rawRead+=Pat_Off; //momo
	int Tot_Ref_Length=0;
/////////////////////////////////////////////////////////////////////	
	int HClip_Mis=0;
	float HClip_QScore=0;
	float HClip_Score=0;
	float HClip_BQScore=0;
	bool HClip_Resc=false;
	if(Qual==NULL) return;//Segment fault
	int lenM=0;
//printf("\nCCC %d %s\n",Clip_H, Cigar_Ptr);
	if(Clip_H && Clip_H<=Ref_clip && Clip_H <= Ref_Off )//Ref_Off)
	{
		Ref-=Clip_H; rawRef-=Clip_H;//momo
		Pattern-=Clip_H;rawRead-=Clip_H;
		Qual-=Clip_H;
		HClip_Resc=true;
//if(!ESTIMATE)  
		for(int i=Clip_H-1;i>=0;i--)
		{
			if(*Ref>4 || *Ref<0 || *Pattern>=4 || *Pattern<0)
			{
				Cig_Err=true;
				break;
			}
			assert(*Ref<4 && *Ref>=0);
			assert(*Pattern<4 && *Pattern>=0);
			assert(*rawRef<4 && *rawRef>=0);
			assert(*rawRead<4 && *rawRead>=0);
			//if(*(Ref+i)!=*(Pattern+i) || (*(Ref+i)==*(Pattern+i) && *(rawRef+i)==3 && *(rawRead+i)==1  ) ||  (*(Ref+i)==*(Pattern+i) && *(rawRef+i)==0 && *(rawRead+i)==2 ) )
			if(*Ref!=*Pattern || (*Ref==*Pattern && *rawRef==3 && *rawRead==1  ) ||  (*Ref==*Pattern && *rawRef==0 && *rawRead==2 ) )
			{
//printf("- %d%d%d%d =",*rawRead,*Pattern,*Ref,*rawRef);
				HClip_Mis++;
				Hclipswscore -= mismatch;
				if(INPUT_FILE_TYPE==FA) 
				{
					HClip_Score+=Mis_FA_Score;
					HClip_QScore+=std::min(QLIMIT_FLOAT,Mis_FA_Score/3);
				}
				else
				{
					float Q_Value=*Qual-QUALITYCONVERSIONFACTOR;//make quality to integer..
					assert(Q_Value>=0);
					//HClip_BQScore-= MN + floor( (MX-MN)*(std::min(Q_Value, 40.0f)/40.0f) );
					float Penalty= -10*log10((1-Pr(Q_Value))/3);
					Penalty=std::min(QLIMIT_FLOAT,Penalty);
					HClip_Score+= Penalty;//Convert to probability of base being wrong =1-10^(Q_Value/10)

					Penalty= Q_Value;///3;//prob II..
					Penalty=std::min(QLIMIT_FLOAT,Penalty/3);
					HClip_QScore+=Penalty;
				}
			}
			else
			{
				Hclipswscore += match;
				if( i+lenM == Clip_H-1) lenM++;
				if(INPUT_FILE_TYPE==FA) 
					HClip_Score+=Match_FA_Score;
				else
				{
					float Q_Value=*Qual-QUALITYCONVERSIONFACTOR;//make quality to integer..
					assert(Q_Value>=0);
					float Penalty= -10*log10(Pr(Q_Value));
					Penalty=std::min(QLIMIT_FLOAT,Penalty);
					assert(Penalty<=QLIMIT); 
					HClip_Score+= Penalty;
					//HClip_BQScore+=MATCH_BONUS;
					//QScore+=Penalty;
				}
			}
			Ref++;Pattern++;Qual++,rawRead++,rawRef++;
		}
	}
	if(Cigar && Clip_H) Cigar_Ptr+=sprintf(Cigar_Ptr,"%dS", Clip_H); //-lenM); 
/////////////////////////////////////////////////////////////////////

	for (c = 0; c < a->cigarLen; ++c) {
		int32_t letter = 0xf&*(a->cigar + c);
		int32_t length = (0xfffffff0&*(a->cigar + c))>>4;
		if(Cigar && c==1) {
			Cigar_Ptr+=sprintf(Cigar_Ptr,"%d", length+lenM);
		}
		else if(Cigar) Cigar_Ptr+=sprintf(Cigar_Ptr,"%d", length);
//printf("\nlen %d %d %s %s %d",length,lenM,Cigar_Ptr,Cigar, c);
		if (letter == 0) 
		{
			Tot_Length+=length;
			Tot_Ref_Length+=length;
			if(Cigar) *Cigar_Ptr='M';
			for(int i=0;i<length;i++)
			{
				if(*Ref>4 || *Ref<0 || *Pattern>=4 || *Pattern<0)
				{
					Cig_Err=true;
					break;
				}
				assert(*Ref<4 && *Ref>=0);
				assert(*Pattern<4 && *Pattern>=0);
				assert(*rawRef<4 && *rawRef>=0);//printf("- %d%d%d%d -",*rawRead,*Pattern,*Ref,*rawRef);
				assert(*rawRead<4 && *rawRead>=0);
				if(*Ref!=*Pattern || (*Ref==*Pattern && *rawRef==3 && *rawRead==1  ) ||  (*Ref==*Pattern && *rawRef==0 && *rawRead==2 ) )
				{
//printf("- %d%d%d%d -",*rawRead,*Pattern,*Ref,*rawRef);
					C.Mis++;
					swscore -= mismatch;
					if(C.Mis>MAX_MISMATCHES && false) 
					{
						//C.Score==INT_MAX;
						return;
					}
					//C.Score-=mismatch;
					if(INPUT_FILE_TYPE==FA) 
					{
						Score+=Mis_FA_Score;
						QScore+=std::min(QLIMIT_FLOAT,Mis_FA_Score/3);
					}
					else
					{
						float Q_Value=*Qual-QUALITYCONVERSIONFACTOR;//make quality to integer..
						assert(Q_Value>=0);
						BQScore-= MN + floor( (MX-MN)*(std::min(Q_Value, 40.0f)/40.0f) );
						float Penalty= -10*log10((1-Pr(Q_Value))/3);
						Penalty=std::min(QLIMIT_FLOAT,Penalty);
						Score+= Penalty;//Convert to probability of base being wrong =1-10^(Q_Value/10)

						Penalty= Q_Value;///3;//prob II..
						Penalty=std::min(QLIMIT_FLOAT,Penalty/3);
						QScore+=Penalty;
					}
				}
				else
				{
					swscore += match;
					if(INPUT_FILE_TYPE==FA) 
						Score+=Match_FA_Score;
					else
					{
						float Q_Value=*Qual-QUALITYCONVERSIONFACTOR;//make quality to integer..
						assert(Q_Value>=0);
						float Penalty= -10*log10(Pr(Q_Value));
						Penalty=std::min(QLIMIT_FLOAT,Penalty);
						assert(Penalty<=QLIMIT); 
						Score+= Penalty;
						BQScore+=MATCH_BONUS;
						//QScore+=Penalty;
					}
				}
				Ref++;Pattern++;Qual++,rawRead++,rawRef++;
			}
			Cigar_Ptr++;C.M+=length;
		}
		else if (letter == 1)
		{
			Tot_Length+=length;
			C.Indel_Count++;
			swscore -= (gap_open+length*gap_extension);
			if(Cigar) *Cigar_Ptr='I';Cigar_Ptr++;C.I+=length;
			Score+=(gap_openN+length*gap_extensionN);
			BQScore-= BOPEN + length*BEXT;
			Pattern+=length;Qual+=length;rawRead+=length;
		} 
		else 
		{
			C.Indel_Count++;
			swscore -= (gap_open+length*gap_extension);
			Score+=(gap_openN+length*gap_extensionN);
			BQScore-= BOPEN + length*BEXT;
			if(Cigar) *Cigar_Ptr='D';Cigar_Ptr++;C.D+=length;
			Ref+=length; rawRef+=length;//momo
			Tot_Ref_Length+=length;
		}
		//printf("\nlen %d %d %s %s",length,lenM,Cigar_Ptr,Cigar);
		if(Cigar_Ptr-Cigar>=MAX_SIGLEN-6) 
		{
			Cig_Err=true;
			break;
		}
	}

	int TClip_Mis=0;
	float TClip_QScore=0;
	float TClip_Score=0;
	float TClip_BQScore=0;
	bool TClip_Resc=false;
//OOOOO 22 128 150 17 1
//printf("\nOOOOO %d %d %d %d %d\n", Clip_T, Tot_Ref_Length, StringLength, INDELGAP, Ref_Off);
	if(Clip_T && Clip_T<Ref_clip && Clip_T+Tot_Ref_Length<StringLength+INDELGAP-Ref_Off) //+INDELGAP-Ref_Off)
	{
		TClip_Resc=true;
		for(int i=0;i<Clip_T;i++)
		{
			if(*Ref>4 || *Ref<0 || *Pattern>=4 || *Pattern<0)
			{
				Cig_Err=true;
				break;
			}
			assert(*Ref<4 && *Ref>=0);
			assert(*Pattern<4 && *Pattern>=0);
			assert(*rawRef<4 && *rawRef>=0);

			if(*Ref!=*Pattern || (*Ref==*Pattern && *rawRef==3 && *rawRead==1  ) ||  (*Ref==*Pattern && *rawRef==0 && *rawRead==2 ) )
			{ 
//printf("= %d%d%d%d -",*rawRead,*Pattern,*Ref,*rawRef);
				Tclipswscore -= mismatch;
				TClip_Mis++;
				if(INPUT_FILE_TYPE==FA) 
				{
					TClip_Score+=Mis_FA_Score;
					TClip_QScore+=std::min(QLIMIT_FLOAT,Mis_FA_Score/3);
				}
				else
				{
					float Q_Value=*Qual-QUALITYCONVERSIONFACTOR;//make quality to integer..
					assert(Q_Value>=0);
					TClip_BQScore-= MN + floor( (MX-MN)*(std::min(Q_Value, 40.0f)/40.0f) );
					float Penalty= -10*log10((1-Pr(Q_Value))/3);
					Penalty=std::min(QLIMIT_FLOAT,Penalty);
					TClip_Score+= Penalty;//Convert to probability of base being wrong =1-10^(Q_Value/10)

					Penalty= Q_Value;///3;//prob II..
					Penalty=std::min(QLIMIT_FLOAT,Penalty/3);
					TClip_QScore+=Penalty;
				}
			}
			else
			{
				Tclipswscore += match;
				if(INPUT_FILE_TYPE==FA) 
					TClip_Score+=Match_FA_Score;
				else
				{
					float Q_Value=*Qual-QUALITYCONVERSIONFACTOR;//make quality to integer..
					assert(Q_Value>=0);
					float Penalty= -10*log10(Pr(Q_Value));
					Penalty=std::min(QLIMIT_FLOAT,Penalty);
					assert(Penalty<=QLIMIT); 
					TClip_Score+= Penalty;
					TClip_BQScore+=MATCH_BONUS;
					//QScore+=Penalty;
				}
			}
			Ref++;Pattern++;Qual++,rawRead++,rawRef++;
		}
	}

	if(Cigar)
	{
		if(Clip_T) Cigar_Ptr+=sprintf(Cigar_Ptr,"%dS", Clip_T); 
		*Cigar_Ptr=0;
		if(Cig_Err) *Cigar=0;
	}

	C.Length=Tot_Length;
	if(Clip_T){
		TClip_Score=Clip_T*30>150?150:Clip_T*30;
		if(!TClip_Resc){
			Score+=TClip_Score;
		}
	}
	if(TClip_Resc)// && TClip_QScore<2*QLIMIT_FLOAT)
	{
		if(Hard_Penalty)
		{
			QScore+=TClip_QScore;
			Score+=TClip_Score;
			BQScore+=TClip_BQScore;
			swscore += Tclipswscore;
		}
		Tot_Length+=Clip_T;
		C.Mis+=TClip_Mis;
	}
//printf("\nHHHs %d %d", HClip_Resc, Clip_H);
	if(Clip_H) {
		HClip_Score=Clip_H*30>150?150:Clip_H*30;
		if(!HClip_Resc) {
			Score+=HClip_Score;
		}
	}
	if(HClip_Resc)// && HClip_QScore<2*QLIMIT_FLOAT)
	{
		if(Hard_Penalty)
		{
			QScore+=HClip_QScore;
			Score+=HClip_Score;
			BQScore+=HClip_BQScore;
			swscore += Hclipswscore;
		}
		Tot_Length+=Clip_H;
		C.Mis+=HClip_Mis;
	}
	assert(Tot_Length<=StringLength);

	if(Tot_Length<StringLength)
	{
		Score+= gap_openN+(StringLength-Tot_Length)*gap_extensionN;
		QScore+= std::min(3,StringLength-Tot_Length)*Match_FA_Score/3;
	}
	C.Score=int(Score);
	C.QScore=int(QScore);
	C.BQScore=int(BQScore);
	C.swscore = int(swscore);
//if(C.BQScore==INT_MAX)
//{
  //      printf("\n22222222222 ssw\n");
//	        exit(0);
//	}
}
