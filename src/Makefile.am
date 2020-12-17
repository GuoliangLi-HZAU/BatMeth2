EXTRA_DIST =  batmeth2
AM_CFLAGS = -w -O3 -funroll-loops -maccumulate-outgoing-args -msse2 -lz -lm -pthread -no-pie
AM_CXXFLAGS = -w -O3 -funroll-loops -maccumulate-outgoing-args -msse2 -lz -lm -pthread -no-pie
RANLIB=ranlib
#AM_CPPFLAGS = -O3 -funroll-loops -maccumulate-outgoing-args -msse2 -lm

noinst_LIBRARIES = libz.a
libz_a_SOURCES = adler32.c compress.c crc32.c gzio.c uncompr.c deflate.c trees.c \
       zutil.c inflate.c infback.c inftrees.c inffast.c\
       zlib.h zconf.h crc32.h  deflate.h zutil.h inftrees.h inflate.h inffast.h trees.h inffixed.h

bin_PROGRAMS = batmethindex penguin penguin-a bwtformatdb reverse filter methyGff build_index_rrbs
batmethindex_SOURCES= bfix.cpp index.cpp BWT.c MiscUtilities.c MemManager.c TextConverter.c r250.c QSufSort.c\
 iniparser.c inistrlib.c dictionary.c DNACount.c Timing.c Socket.c HSP.c HSPstatistic.c karlin.c\
 BWT.h TypeNLimit.h MemManager.h TextConverter.h HSP.h MiscUtilities.h DNACount.h r250.h HSPstatistic.h\
 BWTConstruct.h QSufSort.h r250.h dictionary.h inistrlib.h iniparser.h Timing.h Socket.h karlin.h\
 bfix.h 
penguin_SOURCES = print.cpp filters.cpp utils.cpp batlib.cpp rqindex.cpp penguin.cpp map.cpp ssw.h command.h global.h swroutines.h kseq.h\
 bfix.cpp BWT.c MiscUtilities.c MemManager.c TextConverter.c r250.c QSufSort.c ssw.c command.cpp swroutines.cpp\
 iniparser.c inistrlib.c dictionary.c DNACount.c Timing.c Socket.c HSP.c HSPstatistic.c karlin.c fastsw.cpp\
 print.h filters.h BWT.h TypeNLimit.h MemManager.h TextConverter.h HSP.h MiscUtilities.h DNACount.h r250.h HSPstatistic.h\
 BWTConstruct.h QSufSort.h r250.h dictionary.h inistrlib.h iniparser.h Timing.h Socket.h karlin.h\
 bfix.h rqindex.h utils.h batlib.h common.h map.h fastsw.h
penguin_a_SOURCES = print.cpp filters.cpp utils.cpp batlib.cpp rqindex.cpp penguin-a.cpp map.cpp ssw.h command.h global.h swroutines.h kseq.h\
 bfix.cpp BWT.c MiscUtilities.c MemManager.c TextConverter.c r250.c QSufSort.c ssw.c command.cpp swroutines.cpp\
 iniparser.c inistrlib.c dictionary.c DNACount.c Timing.c Socket.c HSP.c HSPstatistic.c karlin.c fastsw.cpp\
 print.h filters.h BWT.h TypeNLimit.h MemManager.h TextConverter.h HSP.h MiscUtilities.h DNACount.h r250.h HSPstatistic.h\
 BWTConstruct.h QSufSort.h r250.h dictionary.h inistrlib.h iniparser.h Timing.h Socket.h karlin.h\
 bfix.h rqindex.h utils.h batlib.h common.h map.h fastsw.h
bwtformatdb_SOURCES= bwtformatdb.c BWT.c BWTConstruct.c MiscUtilities.c MemManager.c TextConverter.c r250.c QSufSort.c\
 iniparser.c inistrlib.c dictionary.c DNACount.c Timing.c Socket.c HSP.c HSPstatistic.c karlin.c\
 BWT.h TypeNLimit.h MemManager.h TextConverter.h HSP.h MiscUtilities.h DNACount.h r250.h HSPstatistic.h\
 BWTConstruct.h QSufSort.h r250.h dictionary.h inistrlib.h iniparser.h Timing.h Socket.h karlin.h
reverse_SOURCES = reverse.cpp 
filter_SOURCES = filter.cpp 
methyGff_SOURCES = methyGff.cpp
build_index_rrbs_SOURCES = build_index_rrbs.cpp
#batman_SOURCES= batman.cpp BWT.c BWTConstruct.c MiscUtilities.c MemManager.c TextConverter.c r250.c QSufSort.c\
# iniparser.c inistrlib.c dictionary.c DNACount.c Timing.c Socket.c HSP.c HSPstatistic.c karlin.c 
#batman_LDADD = libz.a
#decode_SOURCES= decode.cpp BWT.c BWTConstruct.c MiscUtilities.c MemManager.c TextConverter.c r250.c QSufSort.c\
# iniparser.c inistrlib.c dictionary.c DNACount.c Timing.c Socket.c HSP.c HSPstatistic.c karlin.c
#decode_LDADD = libz.a
#reverse_SOURCES = reverse.cpp 
#copy:
#	cp bwtformatdb ../bin
#	cp reverse ../bin
#	cp batman ../bin
#	cp decode ../bin
