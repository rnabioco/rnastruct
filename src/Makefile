CC = g++
CFLAGS = -std=c++0x
SRC_DIR = mpileup2readcounts
TARGET = $(SRC_DIR)/bin/mpileupToReadCounts
SOURCE = $(SRC_DIR)/mpileup2readcounts.cc

all: mpileupToReadCounts filterBam
	
filterBam: filterBam.cpp 
	$(CC) $(CFLAGS) -o $@ $<
  
mpileupToReadCounts: $(TARGET) 
	install $< $@

$(TARGET): $(SOURCE)
	cd $(SRC_DIR); make; cd ..

clean:
	rm -f mpileupToReadCounts

.PHONY: all clean

