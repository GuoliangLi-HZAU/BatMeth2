AUTOMAKE_OPTIONS = foreign
SUBDIRS = src scripts src/batDMR
script:
	g++ ./src/calmeth.cpp -o ./src/calmeth -m64 -I./src/samtools-0.1.18/ -L./src/samtools-0.1.18/ -lbam -lz
	g++ ./src/splitSam.cpp -o ./src/splitSam -m64 -I./src/samtools-0.1.18/ -L./src/samtools-0.1.18/ -lbam -lz -pthread
copy:
	if [ -d "bin" ]; then echo bin exists; else mkdir bin; fi
	cp scripts/strip.pl bin 
	cp scripts/b2c.pl bin 
	cp scripts/build_complement bin 
	cp scripts/filter.pl bin
	cp scripts/build_indexX bin 
	cp scripts/build_all bin
	cp scripts/ann2loc.pl bin 
	cp scripts/build_location.pl bin 
	cp scripts/build_revcmp bin 
	cp scripts/complement.pl bin
	cp scripts/ReverseComplteFQ bin
	cp src/bwtformatdb bin 
	cp src/reverse bin 
	cp src/penguin bin 
	cp src/penguin-a bin
	cp src/calmeth bin
	cp src/batmethindex bin
	cp src/filter bin
	cp bwtformatdb.ini bin
	cp src/splitSam bin
	cp src/methyGff bin
	cp src/methyPlot bin
	cp src/*.r bin
	cp scripts/*.r bin
	cp src/DMCannotation* bin
	cp src/GeneMethHeatmap ./bin/
	cp scripts/chrLenExtract ./bin
	cp scripts/combined.element* bin
	cp scripts/batmeth2-align bin
	cp scripts/BatMeth2 bin
	cp src/batDMR/batDMR ./bin
	cp src/genome_filter bin
	cp src/build_index_rrbs bin
