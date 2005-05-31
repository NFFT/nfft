NFFT_DIR = /home/kunis/Nfft3/nfft/trunk

# ----- Compiler ---------------------------------------------------------------- 
CC = gcc
CCFLAGS = -Wall -O3 -finline


# -------------------------------------------------------------------------------
all: util.o

util.o: util.c $(NFFT_OPT)
	$(CC) $(CCFLAGS) -c util.c -o util.o -I$(NFFT_DIR)/include

clean:
	rm util.o
