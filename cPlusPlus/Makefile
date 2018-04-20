CC = g++
LD = g++
CCFLAGS = -c -O3 
LDFLAGS = -lm -lpng -lz -L/usr/local/lib
BASE = image.o imagefft.o helpers.o

digitalZoom:	digitalZoom.o helpers.h helpers.o $(BASE)
	$(LD) digitalZoom.o $(BASE) $(LDFLAGS) -o digitalZoom

playground: 	playground.o helpers.h helpers.o $(BASE)
	$(LD) playground.o $(BASE) $(LDFLAGS) -o playground

clean:
	rm  $(BASE) digitalZoom digitalZoom.o playground playground.o

.cc.o:
	$(CC) $(CCFLAGS) $<
