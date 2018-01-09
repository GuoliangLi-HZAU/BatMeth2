#include "swroutines-l.h"
#include "assert.h"
#include <queue>
#include "filters-l.h"
#include "math.h"
#include "fastsw-l.h"
extern long File_size;
extern bool FASTDECODE;
extern int SW_STRING_BUFFER;
extern int Inter_MM;
extern int MODE;
extern int QUALITYCONVERSIONFACTOR;
extern int INDELGAP;//15 good for //17 for 8, 19=g00d for 9
int Do_SW_Pair(char* Pattern_raw,READ & RawR,unsigned char* Original_Text,Pair* Pairs,char* Pattern,int Flank_Size,int StringLength,int & Err,int Shift,char Sign,std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> & Alignments,int Current_Score,int & Total_SW_Scans)
{
	s_align* Aln;
	Alignment A;
	unsigned Ref_Location;
	char Org_String[SW_STRING_BUFFER];
	char Org_String_Ori[SW_STRING_BUFFER];
	int HITS=0;
	Cigar_Info Cig_Info;
	int8_t* mata_tmp; mata_tmp = (int8_t*)calloc(25, sizeof(int8_t));
	char state=RawR.state;if(Sign=='-') state=RawR.state=='1'?'2':'1';
	if(state=='1')
		init_SSW_BS_c2t(mata_tmp);
	else if(state=='2')
		init_SSW_BS_g2a(mata_tmp);
	else init_SSW(mata_tmp);
	
	s_profile* p = ssw_init((int8_t*)Pattern_raw, StringLength/2, mata_tmp, n, 1);

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
		
		Aln=mengyao_ssw_core(Org_String_Ori,StringLength/2, Pattern_raw,Flank_Size,0,(MODE==VERYFAST)?1:0, p);
		if(Aln->score1 >= ACC_SCORE)
		{
			A.SW_Score=Aln->score1;
			if(MODE>=VERYFAST) ssw_cigar_process(Sign,Aln->read_begin1,RawR,Org_String_Ori+Aln->ref_begin1,Aln,Cig_Info,Org_String+Aln->ref_begin1,Pattern+Aln->read_begin1,StringLength/2);
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
	int8_t* mata_tmp;mata_tmp = (int8_t*)calloc(25, sizeof(int8_t));
	char state=RawR.state;if(Sign=='-') state=RawR.state=='1'?'2':'1';
	if(state=='1')
		init_SSW_BS_c2t(mata_tmp);
	else if(state=='2')
		init_SSW_BS_g2a(mata_tmp);
	else init_SSW(mata_tmp);
	
	s_profile* p = ssw_init((int8_t*)Pattern_raw, StringLength/2, mata_tmp, n, 1);
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
				Aln=mengyao_ssw_core(Org_String_Ori,StringLength/2, Pattern_raw,Flank_Size,0,(MODE==VERYFAST)? 1:0, p);
				if(Aln->score1 >=Filter)//ACC_SCORE)
				{
					int Clip_H=Aln->read_begin1;int Clip_T=0;
					if(Aln->read_end1!=StringLength-1) Clip_T=StringLength-1-Aln->read_end1;
					if(Clip_T+Clip_H)
						A.QScore= -1;
					else 
						A.QScore=INT_MAX;

					A.SW_Score=Aln->score1;
					A.Realigned=NO;
					if(MODE==VERYFAST && Aln->score1 >Filter) 
					{
						Filter=Aln->score1+2;
					}
					if(MODE>=VERYFAST) ssw_cigar_process(Sign,Aln->read_begin1,RawR,Org_String_Ori+Aln->ref_begin1,Aln,Cig_Info,Org_String+Aln->ref_begin1,Pattern+Aln->read_begin1,StringLength/2);
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
	if(Location >= File_size || (--Location)<0) Location=0;
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

void ssw_cigar_print(char Sign,int ReadSpan,READ & RawR,char* rawRef,s_align* a,char* Cigar,Cigar_Info & C,char* Ref,char* Pattern,int Clip_H,int Clip_T,int StringLength) { //print the cigar out
	READ RawRs;
	if(Sign=='-')
	{
		int i;int Lens=RawR.Real_Len;
		for (i=0;i<=Lens-1;i++){RawRs.Tag_Copy[Lens-1-i]=Char_To_CharC[RawR.Tag_Copy[i]];};
	}
	else
	{
		RawRs=RawR;
	}
			char* rawRead =RawRs.Tag_Copy;
			rawRead+=ReadSpan;//momo
	int c = 0,Tot_Length=0;char* Cigar_Ptr=Cigar;
	bool Cig_Err=false;
	C.M=0;C.I=0;C.D=0,C.Mis=0;
	if(Clip_H) Cigar_Ptr+=sprintf(Cigar_Ptr,"%dS", Clip_H); 
	for (c = 0; c < a->cigarLen; ++c) {
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
				if(*Ref!=*Pattern )//|| (*Ref==*rawRef && *rawRef==3 && (*rawRead=='c' ||*rawRead=='C' ) ) ||  (*Ref==*rawRef && *rawRef==0 && (*rawRead=='g' ||*rawRead=='G' ) ) ) 
				{C.Mis++;}
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

void ssw_cigar_process(char Sign,int ReadSpan,READ & RawR,char* rawRef,s_align* a,Cigar_Info & C,char* Ref,char* Pattern,int StringLength) { //print the cigar out
	READ RawRs;
	if(Sign=='-')
	{
		int i;int Lens=RawR.Real_Len;
		for (i=0;i<=Lens-1;i++){RawRs.Tag_Copy[Lens-1-i]=Char_To_CharC[RawR.Tag_Copy[i]];};
	}
	else
	{
		RawRs=RawR;
	}
	int gap_openX=40,gap_extensionX=6;
	if(MODE>=FAST)
		gap_openX=40,gap_extensionX=6;
	else
		gap_openX=gap_open,gap_extensionX=gap_extension;

			char* rawRead =RawRs.Tag_Copy;
			rawRead+=ReadSpan;
			
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
				if(*Ref!=*Pattern )//|| (*Ref==*rawRef && *rawRef==3 && (*rawRead=='c' ||*rawRead=='C' ) ) ||  (*Ref==*rawRef && *rawRef==0 && (*rawRead=='g' ||*rawRead=='G' ) )  )
				{
					C.Mis++;
					if(MODE==VERYFAST) C.D-=mismatch;
					//C.Score-=mismatch;
					else
						C.Score+=Mis_FA_Score;
				}
				else
				{
					//C.Score+=match;
					if(MODE==VERYFAST) C.D+=match;
					else	
						C.Score+=Match_FA_Score;
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


void ssw_cigar_processQ(char Sign,READ & RawR,char* rawRef,s_align* a,Cigar_Info & C,char* Ref,int Ref_Off,char* Pattern,int Pat_Off,int StringLength,const char* Qual,char* Cigar,int Clip_H,int Clip_T) { //print the cigar out
	READ RawRs;
	if(Sign=='-')
	{
		int i;int Lens=RawR.Real_Len;
		for (i=0;i<=Lens-1;i++){RawRs.Tag_Copy[Lens-1-i]=Char_To_CharC[RawR.Tag_Copy[i]];};
	}
	else
	{
		RawRs=RawR;
	}
	bool Cig_Err=false;
	const int gap_open=40,gap_extension=6;
	const int MX=6,MN=2,BOPEN=6,BEXT=3,MATCH_BONUS=0;//2;
	int c = 0,Tot_Length=0;char* Cigar_Ptr=Cigar;
	float Score=0,QScore=0,BQScore=0;
	float TScore=0,TQScore=0,TBQScore=0;
	C.M=0;C.I=0;C.D=0,C.Mis=0;C.Score=0;C.Indel_Count=0;

	char* rawRead =RawRs.Tag_Copy;
	Ref+=Ref_Off;Pattern+=Pat_Off;Qual+=Pat_Off;
	rawRef+=Ref_Off;rawRead+=Pat_Off; //momo
	int Tot_Ref_Length=0;
/////////////////////////////////////////////////////////////////////	
	int HClip_Mis=0;
	float HClip_QScore=0;
	float HClip_Score=0;
	float HClip_BQScore=0;
	bool HClip_Resc=false;
	if(Qual==NULL) return;
	if(Cigar && Clip_H) Cigar_Ptr+=sprintf(Cigar_Ptr,"%dS", Clip_H); 
	if(Clip_H && Clip_H<=Ref_Off)
	{
		Ref-=Clip_H;rawRef-=Clip_H;
		Pattern-=Clip_H;rawRead-=Clip_H;
		Qual-=Clip_H;
		HClip_Resc=true;
		for(int i=0;i<Clip_H;i++)
		{
				assert(*Ref<4 && *Ref>=0);
				assert(*Pattern<4 && *Pattern>=0);
			//	assert(*rawRef<4 && *rawRef>=0);
			//	assert(*rawRead!='\0');
			if(*Ref!=*Pattern )//|| ( *Ref==*rawRef &&  *rawRef==3 && (*rawRead=='c' ||*rawRead=='C' ) )  ||  (*Ref==*rawRef && *rawRef==0 && (*rawRead=='g' ||*rawRead=='G' ) )  )
			{
				HClip_Mis++;
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

/////////////////////////////////////////////////////////////////////

	for (c = 0; c < a->cigarLen; ++c) {
		int32_t letter = 0xf&*(a->cigar + c);
		int32_t length = (0xfffffff0&*(a->cigar + c))>>4;
		if(Cigar) Cigar_Ptr+=sprintf(Cigar_Ptr,"%d", length);
		if (letter == 0) 
		{
			Tot_Length+=length;
			Tot_Ref_Length+=length;
			if(Cigar) *Cigar_Ptr='M';
			for(int i=0;i<length;i++)
			{
				assert(*Ref<4 && *Ref>=0);
				assert(*Pattern<4 && *Pattern>=0);
		//		assert(*rawRef<4 && *rawRef>=0);//printf("%c-%d%d ",*rawRead,*Pattern,*rawRef);
		//		assert(*rawRead!='\0');
				if(*Ref!=*Pattern )//|| ( *Ref==*rawRef &&  *rawRef==3 && (*rawRead=='c' ||*rawRead=='C' ) )  ||  (*Ref==*rawRef && *rawRef==0 && (*rawRead=='g' ||*rawRead=='G' ) )  )
				{
					C.Mis++;
					//C.Score-=mismatch;
					if(INPUT_FILE_TYPE==FA) 
					{
						Score+=Mis_FA_Score;
						QScore+=std::min(QLIMIT_FLOAT,Mis_FA_Score/3);
					}
					else
					{
						float Q_Value=*Qual-QUALITYCONVERSIONFACTOR;//make quality to integer..
						if(Q_Value<0) return;//momo
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
					if(INPUT_FILE_TYPE==FA) 
						Score+=Match_FA_Score;
					else
					{
						float Q_Value=*Qual-QUALITYCONVERSIONFACTOR;//make quality to integer..
						if(Q_Value<0) return;//momo
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
			if(Cigar) *Cigar_Ptr='I';Cigar_Ptr++;C.I+=length;
			Score+=(gap_open+length*gap_extension);
			BQScore-= BOPEN + length*BEXT;
			Pattern+=length;Qual+=length;rawRead+=length;
		} 
		else 
		{
			C.Indel_Count++;
			Score+=(gap_open+length*gap_extension);
			BQScore-= BOPEN + length*BEXT;
			if(Cigar) *Cigar_Ptr='D';Cigar_Ptr++;C.D+=length;
			Ref+=length; rawRef+=length;//momo
			Tot_Ref_Length+=length;
		} 
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
	if(Clip_T && Clip_T+Tot_Ref_Length<StringLength+INDELGAP-Ref_Off)
	{
		TClip_Resc=true;
		for(int i=0;i<Clip_T;i++)
		{
			assert(*Ref<4 && *Ref>=0);
			assert(*Pattern<4 && *Pattern>=0);
			if(*Ref!=*Pattern )//|| ( *Ref==*rawRef &&  *rawRef==3 && (*rawRead=='c' ||*rawRead=='C' ) )  ||  (*Ref==*rawRef && *rawRef==0 && (*rawRead=='g' ||*rawRead=='G' ) )  )
			{
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
	if(TClip_Resc)// && TClip_QScore<2*QLIMIT_FLOAT)
	{
		QScore+=TClip_QScore;
		Score+=TClip_Score;
		BQScore+=TClip_BQScore;
		Tot_Length+=Clip_T;
		C.Mis+=TClip_Mis;
	}
	if(HClip_Resc)// && HClip_QScore<2*QLIMIT_FLOAT)
	{
		QScore+=HClip_QScore;
		Score+=HClip_Score;
		BQScore+=HClip_BQScore;
		Tot_Length+=Clip_H;
		C.Mis+=HClip_Mis;
	}
	assert(Tot_Length<=StringLength);

	if(Tot_Length<StringLength)
	{
		Score+= gap_open+(StringLength-Tot_Length)*gap_extension;
		QScore+= std::min(3,StringLength-Tot_Length)*Match_FA_Score/3;
	}
	C.Score=int(Score);
	C.QScore=int(QScore);
	C.BQScore=int(BQScore);
}
