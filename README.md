# Parallelized Streaming SPOD


A parallel implementation of Streaming Spectral Proper Orthogonal Decomposition (SPOD) algorithm presented in [Schmidt and Towne, 2019](https://doi.org/10.1016/j.cpc.2018.11.009) developed at the **Laboratory for Simulations of Flows using Numerics (Lab SoFun)**, Department of Aerospace Engineering, IIT Kanpur, India.

## About

This code implements a parallelized approach to streaming Spectral Proper Orthogonal Decomposition (SPOD), designed for analyzing large-scale fluid flow datasets. Developed by the Laboratory for Simulations of Flows using Numerics (Lab SoFun) headed by Dr. Pradeep Moise at the Department of Aerospace Engineering, IIT Kanpur, India.

SPOD is a modal decomposition technique that identifies coherent structures in fluid flows based on space-time correlations. This implementation leverages parallel computing to efficiently handle large datasets.

## Prerequisites

The following libraries are required:

- MPICH (Message Passing Interface for parallel computing)
- FFTW with MPI support (Fast Fourier Transform library)
- PETSc (Portable, Extensible Toolkit for Scientific Computation)
- SLEPc (Scalable Library for Eigenvalue Problem Computations)

## Installation

### 1. MPICH Installation

```bash
# Download and extract MPICH
tar -xvzf mpich-4.2.3.tar.gz
cd mpich-4.2.3/

# Configure and install
./configure --prefix=/path/to/save/mpich
make -j$(nproc)
make install

# Add to environment (modify path as needed)
echo 'export PATH=/path/to/save/mpich/bin:$PATH' >> ~/.bashrc
echo 'export LD_LIBRARY_PATH=/path/to/save/mpich/lib:$LD_LIBRARY_PATH' >> ~/.bashrc
echo 'export MANPATH=/path/to/save/mpich/share/man:$MANPATH' >> ~/.bashrc
source ~/.bashrc
```

### 2. FFTW with MPI Support

```bash
# Download and extract FFTW
tar -xvzf fftw-3.3.10.tar.gz
cd fftw-3.3.10/

# Configure and install with MPI support
./configure --enable-mpi --enable-shared --prefix=/path/to/save/fftwmpi
make -j$(nproc)
make install

# Add to environment
echo 'export LD_LIBRARY_PATH=/path/to/save/fftwmpi/lib:$LD_LIBRARY_PATH' >> ~/.bashrc
source ~/.bashrc
```

### 3. PETSc Installation

```bash
# Download and extract PETSc
tar -xvzf petsc-with-docs-3.22.2.tar.gz
cd petsc-3.22.2/

# Configure with MPI directory and complex scalar type
./configure --with-mpi-dir=/path/to/save/mpich --download-fblaslapack --with-scalar-type=complex

# Build and test
make PETSC_DIR=/path/to/save/petsc-3.22.2 PETSC_ARCH=arch-linux-c-debug all
make PETSC_DIR=/path/to/save/petsc-3.22.2 PETSC_ARCH=arch-linux-c-debug check

# Set environment variables
export PETSC_DIR=/path/to/save/petsc-3.22.2
export PETSC_ARCH=arch-linux-c-debug
```

### 4. SLEPc Installation

```bash
# Download and extract SLEPc
tar -xvzf slepc-with-docs-3.22.2.tar.gz
cd slepc-3.22.2/

# Configure and build
export SLEPC_DIR=/path/to/save/slepc-3.22.2
./configure
make SLEPC_DIR=/path/to/save/slepc-3.22.2 PETSC_DIR=/path/to/save/petsc-3.22.2 PETSC_ARCH=arch-linux-c-debug
make SLEPC_DIR=/path/to/save/slepc-3.22.2 PETSC_DIR=/path/to/save/petsc-3.22.2 check

# Set environment variable
export SLEPC_DIR=/path/to/save/slepc-3.22.2
```

## Environment Setup

After installation, you need to set these environment variables every time you restart the server:

```bash
export PETSC_DIR=/path/to/save/petsc-3.22.2
export PETSC_ARCH=arch-linux-c-debug
export SLEPC_DIR=/path/to/save/slepc-3.22.2
export LD_LIBRARY_PATH=/path/to/save/mpich/lib:$LD_LIBRARY_PATH
```

For convenience, you can add these lines to your `.bashrc` file or create a setup script.


## Usage

### Compilation

1. First, modify the `Makefile` to point to your FFTW installation:

```makefile
FC = gfortran
FFTW_DIR = /Your/Path/To/fftw3mpi
FFLAGS = -O2 -Wall -I$(FFTW_DIR)/include 
LDFLAGS = -L$(FFTW_DIR)/lib -lfftw3_mpi -lfftw3 -lm
TARGET = streamingSPOD
SRC = streamingSPODModule.o
$(TARGET): $(SRC)
	$(FLINKER) $(FFLAGS) -o $(TARGET) $(SRC) ${SLEPC_SYS_LIB} $(LDFLAGS) 
 
include ${SLEPC_DIR}/lib/slepc/conf/slepc_common
```

2. Compile the code:
```bash
make
```

### Configuration

Before running the code, you need to modify the parameters in the `streamingSPOD.f90` file. Open the file and locate the `userinput` subroutine. Update the following parameters according to your dataset:

```fortran
! Number of data points in each snapshot
numberOfDataInSnapShot = 374740

! Total number of snapshots
numberOfSnapShot = 1400

! Number of data points for DFT window
nDFT = 700

! K modes
k = 3

! Overlap percentage between DFT windows (0.0 to 1.0)
percentageOfOverall = 0.5
! User can change this value as needed

! Filename format for output files
! Example: '/path/to/your/data/', I0, '.bin'
filenameFormat = '("/path/to/your/data/", I0, ".bin")'
! ^^^ Change this path and pattern to match your data files

! Weight matrix filepath
! NOTE: IF PATH IS SET AS EMPTY i.e. '' THEN 1 WILL BE SET AS WEIGHT i.e. EQUAL WEIGHTS
weightMatrixFilePath = 'weightMatrix.csv'

! Hamming window or not
useHamming = .false.   ! Set to .true. to use Hamming window

! Number of frequencies to save.
! This will be used to save the nFreqToSave frequencies corresponding to the nFreqToSave largest eigenvalues.
nFreqToSave = 5
```

Key parameters to configure:
- `numberOfDataInSnapShot`: The size of each snapshot
- `numberOfSnapShot`: Total number of snapshots in your dataset
- `nDFT`: Size of the DFT window
- `k`: Number of modes to extract
- `percentageOfOverall`: Overlap percentage between windows
- `filenameFormat`: Path to your binary data files
- `weightMatrixFilePath`: Path to weight matrix (leave empty for equal weights)
- `nFreqToSave`: Number of frequencies to save

### Running the Code

After compilation and configuration, run the code using MPI:

```bash
mpirun -np 16 ./streamingSPOD
```

You can adjust the number of processes (`-np 16`) based on your system's capabilities.


### Performance Guidelines

For optimal performance, consider the following recommendations for different dataset sizes:

- **Small datasets** (e.g., laminar separation bubble, ~34,000 grid points): Use 2-4 cores
- **Medium-scale problems** (e.g., 2D transonic buffet, ~374,740 grid points): Use 8-16 cores
- **Large-scale simulations** (e.g., 3D vortex breakdown, ~20,268,100 grid points): Use 8-16 cores

While the code maintains speedup compared to MATLAB implementation across 1-64 cores, the above configurations are recommended for maximum performance.

### Input Data Format

**IMPORTANT NOTES:**
1. **Snapshot files** must be in binary format with all data stored as a single vector
2. **Weight matrix** must be in CSV format

Sample snapshot and weight matrix files are provided for reference.

You can visualize binary files in MATLAB to better understand their structure.

### Output Format

The SPOD modes extracted are saved as `.dat` files and can be opened in MATLAB using the PETSc binary reader provided.


## Citation

If you use this code in your research, please cite:

```bibtex
@misc{labsofun_spod,
  author = {{Laboratory for Simulations of Flows using Numerics}},
  title = {Parallelized Streaming SPOD},
  year = {2025},
  publisher = {GitHub},
  howpublished = {\url{https://github.com/yourusername/parallelized-streaming-spod}}
}
```

For the methodology, please also cite:

```bibtex
@article{schmidt2019efficient,
  title={An efficient streaming algorithm for spectral proper orthogonal decomposition},
  author={Schmidt, Oliver T and Towne, Aaron},
  journal={Computer Physics Communications},
  volume={237},
  pages={98--109},
  year={2019},
  publisher={Elsevier}
}
```

## Acknowledgments

- Laboratory for Simulations of Flows using Numerics (Lab SoFun)
- Dr. Pradeep Moise, Head of Lab SoFun
- Department of Aerospace Engineering, IIT Kanpur, India


