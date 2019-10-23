CFLAGS = -std=c++0x

all: filterBam
	
filterBam: filterBam.cpp 
	$(CC) $(CFLAGS) -o $@ $<
  
clean:
	rm -f filterBam 

.PHONY: all clean

