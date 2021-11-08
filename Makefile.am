AUTOMAKE_OPTIONS = foreign
SUBDIRS = src/mealign src/samtools-0.1.18 src scripts src/batDMR
script:
	g++ $(CXXFLAGS) ./scripts/report2html.cpp -o ./scripts/report2html 
myinstall:
	if [ -d "bin" ]; then echo bin exists; else mkdir bin; fi
	cp scripts/batmeth2_to_bigwig.py bin
	cp scripts/bedGraphToBigWig bin
	cp scripts/bedSort bin
	g++ -o ./scripts/BatMeth2 ./scripts/BatMeth2.cpp -lpthread
	g++ ./scripts/report2html.cpp -o ./scripts/report2html
	cp ./scripts/BatMeth2 ./bin/batmeth2
	cp scripts/strip.pl bin 
	cp scripts/report2html bin
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
	cp scripts/GeneMethHeatmap ./bin/
	cp scripts/chrLenExtract ./bin
	cp scripts/combined.element* bin
	cp scripts/BatMeth2 bin
	cp scripts/build_ann_location.pl bin
	cp scripts/preGenome bin
	cp src/batDMR/batDMR ./bin
	cp src/genome2cg bin
	cp src/genomebinLen bin
	cp src/memalign bin
	cp src/build_index_rrbs bin
	cp src/mealign/bwame bin
	cp src/build_index_rrbs bin
copy:
	if [ -d "bin" ]; then echo bin exists; else mkdir bin; fi
	g++ -o ./scripts/BatMeth2 ./scripts/BatMeth2.cpp -lpthread
	g++ ./scripts/report2html.cpp -o ./scripts/report2html
	cp ./scripts/BatMeth2 ./bin/batmeth2
	cp scripts/batmeth2_to_bigwig.py bin
	cp scripts/bedGraphToBigWig bin
	cp scripts/bedSort bin
	cp scripts/strip.pl bin 
	cp scripts/report2html bin
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
	cp scripts/GeneMethHeatmap ./bin/
	cp scripts/chrLenExtract ./bin
	cp scripts/combined.element* bin
	cp scripts/BatMeth2 bin
	cp scripts/build_ann_location.pl bin
	cp scripts/preGenome bin
	cp src/batDMR/batDMR ./bin
	cp src/genome2cg bin
	cp src/genomebinLen bin
	cp src/memalign bin
	cp src/build_index_rrbs bin
	cp src/mealign/bwame bin
	cp src/build_index_rrbs bin
