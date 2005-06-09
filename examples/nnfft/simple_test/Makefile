# ----- NFFT/FFTW library and path ---------------------------------------
NFFT_DIR = /home/ana/nfft2
FFTW_DIR = /home/kunis/DevelopNfft/dagobert_fftw3
EXAMPLE_NAME = simple_test
# ----- Compiler --------------------------------------------------------------- 
CC = gcc
CCFLAGS = -Wall -O3 -finline

# ----- NFFT library and path --------------------------------------------------
NFFT_LIB = $(NFFT_DIR)/lib/libnfft.a

# ----- Link flags  ------------------------------------------------------------
LINK_FLAGS = -lnfft -lfftw3 -lm

# ----- Example programs -------------------------------------------------------
PROG_SRC = $(EXAMPLE_NAME).c
PROG_OBJ = $(EXAMPLE_NAME).o
PROG_EXE = $(EXAMPLE_NAME)

# ------------------------------------------------------------------------------
all: $(PROG_EXE)

$(EXAMPLE_NAME): $(EXAMPLE_NAME).o $(NFFT_LIB)
	$(CC) $(CCFLAGS) -o $(EXAMPLE_NAME) $(EXAMPLE_NAME).o $(LINK_FLAGS) -L$(NFFT_DIR)/lib -L$(FFTW_DIR)/lib 

$(EXAMPLE_NAME).o: $(EXAMPLE_NAME).c
	$(CC) $(CCFLAGS) -c $(EXAMPLE_NAME).c -o $(EXAMPLE_NAME).o -I$(NFFT_DIR)/include  -I$(FFTW_DIR)/include

clean:
	rm $(EXAMPLE_NAME) $(EXAMPLE_NAME).o
