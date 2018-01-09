EXTRA_DIST =  batalign
AM_CFLAGS = -w -O3 -funroll-loops -maccumulate-outgoing-args -msse2 -lm -pthread
AM_CXXFLAGS = -w -O3 -funroll-loops -maccumulate-outgoing-args -msse2 -lm -pthread
RANLIB=ranlib
#AM_CPPFLAGS = -O3 -funroll-loops -maccumulate-outgoing-args -msse2 -lm

noinst_LIBRARIES = libz.a
libz_a_SOURCES = adler32.c compress.c crc32.c gzio.c uncompr.c deflate.c trees.c \
       zutil.c inflate.c infback.c inftrees.c inffast.c\
       zlib.h zconf.h crc32.h  deflate.h zutil.h inftrees.h inflate.h inffast.h trees.h inffixed.h

bin_PROGRAMS = batmethindex penguin penguin-a bwtformatdb reverse filter penguin-s penguin-l 
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
penguin_s_SOURCES = print-s.cpp filters-s.cpp utils-s.cpp batlib-s.cpp rqindex-s.cpp penguin-s.cpp map-s.cpp ssw-s.h command-s.h global-s.h swroutines-s.h kseq.h\
 bfix.cpp BWT-s.c MiscUtilities-s.c MemManager-s.c TextConverter-s.c r250-s.c QSufSort-s.c ssw-s.c command-s.cpp swroutines-s.cpp\
 iniparser-s.c inistrlib-s.c dictionary-s.c DNACount-s.c Timing-s.c Socket-s.c HSP-s.c HSPstatistic-s.c karlin-s.c fastsw-s.cpp\
 print-s.h filters-s.h BWT-s.h TypeNLimit-s.h MemManager-s.h TextConverter-s.h HSP-s.h MiscUtilities-s.h DNACount-s.h r250-s.h HSPstatistic-s.h\
 BWTConstruct-s.h QSufSort-s.h r250-s.h dictionary-s.h inistrlib-s.h iniparser-s.h Timing-s.h Socket-s.h karlin-s.h\
 bfix.h rqindex-s.h utils-s.h batlib-s.h common-s.h map-s.h fastsw-s.h
penguin_l_SOURCES = print-l.cpp filters-l.cpp utils-l.cpp batlib-l.cpp rqindex-l.cpp penguin-l.cpp map-l.cpp ssw-l.h command-l.h global-l.h swroutines-l.h kseq-l.h\
 bfix.cpp BWT-l.c MiscUtilities-l.c MemManager-l.c TextConverter-l.c r250-l.c QSufSort-l.c ssw-l.c command-l.cpp swroutines-l.cpp\
 iniparser-l.c inistrlib-l.c dictionary-l.c DNACount-l.c Timing-l.c Socket-l.c HSP-l.c HSPstatistic-l.c karlin-l.c fastsw-l.cpp\
 print-l.h filters-l.h BWT-l.h TypeNLimit-l.h MemManager-l.h TextConverter-l.h HSP-l.h MiscUtilities-l.h DNACount-l.h r250-l.h HSPstatistic-l.h\
 BWTConstruct-l.h QSufSort-l.h r250-l.h dictionary-l.h inistrlib-l.h iniparser-l.h Timing-l.h Socket-l.h karlin-l.h\
 bfix-l.h rqindex-l.h utils-l.h batlib-l.h common-l.h map-l.h fastsw-l.h
reverse_SOURCES = reverse.cpp 
filter_SOURCES = filter.cpp 
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
