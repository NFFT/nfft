# ----- NFFT/FFTW library and path ---------------------------------------
NFFT_DIR = /home/knopp/Rest/nfft2
FFTW_DIR = /home/include

OBJECTS = construct_data.o reconstruct_data_2+1d.o reconstruct_data_3d.o reconstruct_data_gridding.o

# ----- Compiler --------------------------------------------------------------- 
CC = gcc
CCFLAGS = -Wall -g -O3 -finline

# ----- NFFT library and path --------------------------------------------------
NFFT_LIB = $(NFFT_DIR)/lib/libnfft.a

# ----- Link flags  ------------------------------------------------------------
LINK_FLAGS = -lnfft -lfftw3 -lm

.PHONY:all clean
all:	construct_data reconstruct_data_2+1d reconstruct_data_3d reconstruct_data_gridding
clean: 
	-rm -f $(OBJECTS)


construct_data: construct_data.o Makefile
	$(CC) $(CCFLAGS) -o construct_data construct_data.o $(LINK_FLAGS) -L$(NFFT_DIR)/lib -L$(FFTW_DIR)/lib 

reconstruct_data_2+1d: reconstruct_data_2+1d.o Makefile
	$(CC) $(CCFLAGS) -o reconstruct_data_2+1d reconstruct_data_2+1d.o $(LINK_FLAGS) -L$(NFFT_DIR)/lib -L$(FFTW_DIR)/lib

reconstruct_data_3d: reconstruct_data_3d.o Makefile
	$(CC) $(CCFLAGS) -o reconstruct_data_3d reconstruct_data_3d.o $(LINK_FLAGS) -L$(NFFT_DIR)/lib -L$(FFTW_DIR)/lib
        
reconstruct_data_gridding: reconstruct_data_gridding.o Makefile
	$(CC) $(CCFLAGS) -o reconstruct_data_gridding reconstruct_data_gridding.o $(LINK_FLAGS) -L$(NFFT_DIR)/lib -L$(FFTW_DIR)/lib
                
%.o: %.c Makefile
	$(CC) $(CCFLAGS) -o $@ -c $< -I$(NFFT_DIR)/include  -I$(FFTW_DIR)/include

