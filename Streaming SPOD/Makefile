FC = gfortran
FFTW_DIR = /Your/Path/To/fftw3mpi
FFLAGS = -O2 -Wall -I$(FFTW_DIR)/include 
LDFLAGS = -L$(FFTW_DIR)/lib -lfftw3_mpi -lfftw3 -lm

TARGET = streamingSPOD
SRC = streamingSPODModule.o

$(TARGET): $(SRC)
	$(FLINKER) $(FFLAGS) -o $(TARGET) $(SRC) ${SLEPC_SYS_LIB} $(LDFLAGS) 
 
include ${SLEPC_DIR}/lib/slepc/conf/slepc_common
