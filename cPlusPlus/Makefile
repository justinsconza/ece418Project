CC = g++
LD = g++
CCFLAGS = -c -O3 
LDFLAGS = -lm -lpng -lz -L/usr/local/lib
BASE = image.o imagefft.o helpers.o

digitalZoom:	digitalZoom.o helpers.h helpers.o $(BASE)
	$(LD) digitalZoom.o $(BASE) $(LDFLAGS) -o digitalZoom

dZCropFirst:	dZCropFirst.o helpers.h helpers.o $(BASE)
	$(LD) dZCropFirst.o $(BASE) $(LDFLAGS) -o dZCropFirst

playground: 	playground.o helpers.h helpers.o $(BASE)
	$(LD) playground.o $(BASE) $(LDFLAGS) -o playground

clean:
	rm  $(BASE) digitalZoom digitalZoom.o dZCropFirst dZCropFirst.o playground playground.o

.cc.o:
	$(CC) $(CCFLAGS) $<
