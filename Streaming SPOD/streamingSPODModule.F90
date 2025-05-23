
Module streamingSPODModule
#include <slepc/finclude/slepcsvd.h>
                            
                              use, intrinsic :: iso_c_binding
                              use slepc
                              use slepcsvd
                            
                              implicit none


                            ! ================================================================
                            ! USER-DEFINED CONSTANTS
                            ! ================================================================

                            integer :: numberOfDataInSnapShot 
                            integer :: numberOfSnapShot      
                            integer :: nDFT     
                            integer :: k                               
                            real :: percentageOfOverall 
                            logical :: useHamming
                            integer :: nFreqToSave 

                            ! User can change these values as needed

                            character(len=500) :: filenameFormat
                            character(len=500) :: weightMatrixFilePath

                            ! User can change this path and pattern as needed

                            
                            ! ================================================================
                            ! VARIABLES SET INTERNALLY
                            ! ================================================================
                            
                            integer :: nFreq                  
                            integer :: nOvlp                 ! Number of overlapping points
                            integer :: currentSnapshotNumber ! Current snapshot number
                            integer :: fileIndex             ! File index for output
                            character(len=200) :: filename   ! Output filename

                            ! ================================================================
                            ! VARIABLES FOR WINDOW
                            ! ================================================================
                            real :: winWeight
                            PetscScalar, allocatable :: window(:)


                            ! ================================================================
                            !  MPI VARIABLES AND ERROR VARIABLES
                            ! ================================================================
                            ! Declare PETSc and MPI-related variables
                            PetscMPIInt    :: rank, size   ! rank: process rank, size: total number of processes
                            PetscErrorCode :: ierr         ! ierr: error code returned by PETSc/MPI routines

                            ! ================================================================
                            !  VARIABLES USED TO DECIDE WHICH BLOCK TO UPDATE 
                            ! ================================================================
                            integer, allocatable :: t_idx(:)
                            integer :: dn

                            ! ================================================================
                            !  VARIABLES USED TO DFT
                            ! ================================================================
                            PetscScalar, allocatable :: fourierF90Array(:,:)

                            ! ================================================================
                            !  VARIABLES USED TO MEAN UPDATE
                            ! ================================================================
                            integer :: t_i

                            ! ================================================================
                            !  VARIABLES USED TO CONSTRUCT K MATRIX
                            ! ================================================================
                            integer :: current_Block

                            ! ================================================================
                            !  PETSC MATRIX 
                            ! ================================================================
                            Mat :: X_sum_BlockA
                            Mat :: X_sum_BlockB
                            Mat :: mu
                            Mat :: x_new
                            Mat :: U_hat
                            Mat :: S_hat
                            Mat :: FourierRange ! Only used for range

                            ! ================================================================
                            !  PETSC VECTORS
                            ! ================================================================
                            Vec :: sqrtW


                              contains
                              Subroutine UserInput
                                implicit none
                                ! =======================
                                ! USER-DEFINED VARIABLES
                                ! =======================

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
                                ! Example: '/home/mitul/S1.5_alph1_new/', I0, '.bin'
                                filenameFormat = '("/home/mitul/Trying SLEPC/slepc fortran example 1/parallelJuelData/juelbinary/", I0, ".bin")'
                                ! ^^^ User can change this path and pattern as needed

                                ! Weight matrix filepath
                                ! NOTE: IF PATH IS SET AS EMPTY i,e '' THEN 1 WILL BE SET AS WEIGHT i,e EQUAL WEIGHTS
                                weightMatrixFilePath = 'weightMatrix.csv'

                                ! Hamming window or not
                                useHamming = .false.   ! User can set this to .false. if desired

                                ! Number of frequencies to save.
                                ! This will be used to save the nFreqToSave frequencies corresponding to the nFreqToSave largest eigenvalues.
                                nFreqToSave = 5


                              End Subroutine UserInput

                            
                                Subroutine InitializeBlocks
                                  implicit none
                    
                                  ! =========================
                                  ! LOCAL VARIABLES OF InitializeBlocks
                                  ! =========================
                                  integer ::  blk_i, nBlks

                                  ! =========================
                                  ! LOCAL VARIABLES FOR SHAT RANGES AND SETTING IT AS 0
                                  ! =========================
                                  PetscInt :: rstart, rend, cstart, cend, row, col 

                                  



                                  ! Calculate number of overlapping points
                                  nOvlp = int(percentageOfOverall * nDFT)
                                  nFreq = nDFT / 2 + 1

                                  ! =========================
                                  ! OUTPUT FOR CONFIRMATION
                                  ! =========================
                                  if (rank == 0) then

                                      print *, "numberOfDataInSnapShot = ", numberOfDataInSnapShot
                                      print *, "numberOfSnapShot      = ", numberOfSnapShot
                                      print *, "nDFT                 = ", nDFT
                                      print *, "nFreq                 = ", nFreq
                                      print *, "percentageOfOverall  = ", percentageOfOverall
                                      print *, "nOvlp                = ", nOvlp
                                      print *, "filename             = ", filenameFormat

                                  end if

                                  ! =========================
                                  ! LOGIC TO DECIDE WHEN THE BLOCKS WILL BE UPDATES 
                                  ! BECAUSE SNAPSHOT CAN BE IN OVERLAP
                                  ! =========================
                                  current_Block = 1
                                  nBlks = 2
                                  dn = nDFT - nOvlp
                                  allocate(t_idx(nBlks))
                                  t_idx = 0
                                  do blk_i = 1, nBlks
                                    t_idx(blk_i) = t_idx(blk_i) - (blk_i - 1)*dn + 1
                                  end do


                                  ! =========================
                                  ! FIRST BLOCK
                                  ! =========================
                                  call MatCreate(PETSC_COMM_WORLD, X_sum_BlockA, ierr)
                                  call MatSetSizes(X_sum_BlockA, PETSC_DECIDE, PETSC_DECIDE, numberOfDataInSnapShot, nFreq, ierr)
                                  call MatSetType(X_sum_BlockA, MATMPIDENSE , ierr)
                                  call MatSetUp(X_sum_BlockA, ierr)

                                  ! =========================
                                  ! SECOND BLOCK
                                  ! =========================
                                  if (nOvlp > 0) then
                                    call MatCreate(PETSC_COMM_WORLD, X_sum_BlockB, ierr)
                                    call MatSetSizes(X_sum_BlockB, PETSC_DECIDE, PETSC_DECIDE, numberOfDataInSnapShot, nFreq, ierr)
                                    call MatSetType(X_sum_BlockB, MATMPIDENSE , ierr)
                                    call MatSetUp(X_sum_BlockB, ierr)
                                  end if

                                  ! =========================
                                  ! USED TO STORE THE SNAPSHOT
                                  ! =========================
                                  call MatCreate(PETSC_COMM_WORLD, x_new, ierr)
                                  call MatSetSizes(x_new, PETSC_DECIDE, PETSC_DECIDE, numberOfDataInSnapShot, 1, ierr)
                                  call MatSetType(x_new, MATMPIDENSE, ierr)
                                  call MatSetFromOptions(x_new, ierr)
                                  call MatSetUp(x_new, ierr)

                                  ! =========================
                                  ! USED TO STORE THE MEAN
                                  ! =========================
                                  call MatCreate(PETSC_COMM_WORLD, mu, ierr)
                                  call MatSetSizes(mu, PETSC_DECIDE, PETSC_DECIDE, numberOfDataInSnapShot, 1, ierr)
                                  call MatSetType(mu, MATDENSE, ierr)
                                  call MatSetFromOptions(mu, ierr)
                                  call MatSetUp(mu, ierr)
                                  call MatZeroEntries(mu, ierr)
                                  call MatAssemblyBegin(mu, MAT_FINAL_ASSEMBLY, ierr)        
                                  call MatAssemblyEnd(mu, MAT_FINAL_ASSEMBLY, ierr)
 
                                  ! =========================
                                  ! USED TO STORE THE LEFT SINGULAR VECTOR
                                  ! =========================
                                  call MatCreate(PETSC_COMM_WORLD, U_hat, ierr)
                                  call MatSetSizes(U_hat, PETSC_DECIDE, PETSC_DECIDE, numberOfDataInSnapShot, k*nFreq, ierr)
                                  call MatSetType(U_hat, MATDENSE, ierr)
                                  call MatSetUp(U_hat, ierr)

                                  ! =========================
                                  ! USED TO STORE THE EIGEN VALUES
                                  ! =========================
                                  call MatCreate(PETSC_COMM_WORLD, S_hat, ierr)
                                  call MatSetSizes(S_hat, PETSC_DECIDE, PETSC_DECIDE, k, nFreq, ierr)
                                  call MatSetType(S_hat, MATAIJ, ierr)
                                  call MatSetUp(S_hat, ierr)
                                  call MatGetOwnershipRange(S_hat, rstart, rend, ierr)
                                  call MatGetOwnershipRangeColumn(S_hat, cstart, cend, ierr)
                                  do row = rstart, rend - 1
                                    do col = 0, nFreq - 1
                                        call MatSetValue(S_hat, row, col, (0.0d0,0.0d0), INSERT_VALUES, ierr)
                                    end do
                                  end do
                                  call MatAssemblyBegin(S_hat, MAT_FINAL_ASSEMBLY, ierr)        
                                  call MatAssemblyEnd(S_hat, MAT_FINAL_ASSEMBLY, ierr)
                


                                End Subroutine InitializeBlocks
                            
                                
                                Subroutine ReadWeightMatrix
                                  implicit none
                    
                                  ! =========================
                                  ! LOCAL VARIABLES OF ReadWeightMatrix
                                  ! =========================
                                  integer :: unitno
                                  real(kind=8), allocatable :: w(:)
                                  PetscScalar, allocatable :: sqrtWArrayF90(:)
                                  integer, allocatable :: indices(:)
                                  integer :: i

                                  call VecCreate(PETSC_COMM_WORLD, sqrtW, ierr)
                                  call VecSetSizes(sqrtW, PETSC_DECIDE, numberOfDataInSnapShot, ierr)
                                  call VecSetFromOptions(sqrtW, ierr)

                                  if (trim(weightMatrixFilePath) /= '') then
                                    ! Only rank 0 reads and processes the weight matrix file
                                    if (rank == 0) then
                                        allocate(indices(numberOfDataInSnapShot))
                                        allocate(w(numberOfDataInSnapShot))
                                        allocate(sqrtWArrayF90(numberOfDataInSnapShot))

                                        ! Open and read the weight matrix file
                                        unitno = 10
                                        open(unit=unitno, file=trim(weightMatrixFilePath), status='old')
                                        read(unitno, *) w(:)
                                        close(unitno)
                                
                                        ! Compute the square root of the weights
                                        sqrtWArrayF90 = sqrt(w)
                                
                                        ! Fill indices array (Fortran is 1-based, PETSc is 0-based)
                                        do i = 1, numberOfDataInSnapShot
                                            indices(i) = i - 1
                                        end do
              
                                        ! Set values in the PETSc vector
                                        call VecSetValues(sqrtW, numberOfDataInSnapShot, indices, sqrtWArrayF90, INSERT_VALUES, ierr)
                                
                                        ! Deallocate arrays
                                        deallocate(w, sqrtWArrayF90, indices)
                                      end if
                                    else
                                        ! If no file is provided, set all entries to 1.0 (complex)
                                        call VecSet(sqrtW, cmplx(1.0, 0.0, kind=c_double), ierr)
                                    end if

                                  call VecAssemblyBegin(sqrtW, ierr)
                                  call VecAssemblyEnd(sqrtW, ierr)
                                 
                                End Subroutine ReadWeightMatrix

                              
                                subroutine WindowingFunction()
                                  implicit none

                                  integer :: i
                                  real :: sum_window

                                  ! =========================
                                  ! LOCAL VARIABLES OF WindowingFunction
                                  ! =========================
                            
                                  allocate (window(nDFT))


                                  if (rank == 0) then
                                    do i = 1, nDFT
                                      if (useHamming) then
                                        window(i) = 0.54 - 0.46 * cos(2.0 * 3.141592653589793 * (i-1) / (nDFT-1))
                                      else
                                        window(i) = 1.0
                                      end if
            
                                    end do
                                    sum_window = sum(window)
                            
                                    winWeight = 1.0 / (sum_window / real(nDFT))
            
                                    print *, 'winWeight on RANK 0', winWeight
                    
                                    call MPI_Bcast(window , nDFT * 2, MPI_DOUBLE_PRECISION, 0, PETSC_COMM_WORLD, ierr)
                                        if (ierr /= 0) then
                                            print *, 'Error broadcasting matrix size.'
                                            call PetscFinalize(ierr)
                                            call SlepcFinalize(ierr)
                                            stop
                                        end if
                    
                                    call MPI_Bcast(winWeight, 1, MPI_INT, 0, PETSC_COMM_WORLD, ierr)
                                        if (ierr /= 0) then
                                            print *, 'Error broadcasting matrix size.'
                                            call PetscFinalize(ierr)
                                            call SlepcFinalize(ierr)
                                            stop
                                        end if
            
                                  else 
                                    call MPI_Bcast(window , nDFT * 2, MPI_DOUBLE_PRECISION, 0, PETSC_COMM_WORLD, ierr)
                                    if (ierr /= 0) then
                                        print *, 'Error broadcasting matrix size.'
                                        call PetscFinalize(ierr)
                                        call SlepcFinalize(ierr)
                                        stop
                                    end if
                    
                                    call MPI_Bcast(winWeight, 1, MPI_INT, 0, PETSC_COMM_WORLD, ierr)
                                    if (ierr /= 0) then
                                        print *, 'Error broadcasting matrix size.'
                                        call PetscFinalize(ierr)
                                        call SlepcFinalize(ierr)
                                        stop
                                    end if
                    
                                  end if 
                                end subroutine WindowingFunction
                            
                            
                                subroutine PerformDFT()
                                  use, intrinsic :: iso_c_binding
                            
                                  implicit none
                                  include 'fftw3-mpi.f03'
                                 
                                  ! FFTW variables
                                  complex(c_double), allocatable :: data_in(:,:)
                                  complex(c_double), allocatable :: data_out(:,:)
                                  type(C_PTR) :: plan_ptr


                                  integer :: local_start, local_end
                                  integer :: i, j

                                  call fftw_mpi_init()

                                  call MatCreate(PETSC_COMM_WORLD, FourierRange, ierr)
                                  call MatSetSizes(FourierRange, PETSC_DECIDE, PETSC_DECIDE, nDFT, 1, ierr)
                                  call MatSetType(FourierRange, MATAIJ, ierr)  
                                  call MatSetFromOptions(FourierRange, ierr)
                                  call MatSetUp(FourierRange, ierr)       
                                  call MatGetOwnershipRange(FourierRange, local_start, local_end, ierr)
  

                                  allocate(fourierF90Array(local_end - local_start,nFreq))
        
                                  fourierF90Array = 0.0
  
                                  allocate(data_in(nDFT, local_end - local_start))
                                  allocate(data_out(nDFT, local_end - local_start))
                                  data_in = (0.0_C_DOUBLE, 0.0_C_DOUBLE)
                                  do i = 1, nDFT
                                      do j = local_start + 1, local_end
                                        if (i == j) then
                                            data_in(i, j - local_start) = 1.0_C_DOUBLE  ! Assigning diagonal elements
                                        end if
                                      end do
                                  end do
                                  do j = 1, local_end - local_start
                                    plan_ptr = fftw_plan_dft_1d(nDFT, data_in(:, j), data_out(:, j), FFTW_FORWARD, FFTW_ESTIMATE)
                                    call dfftw_execute(plan_ptr)
                                    call fftw_destroy_plan(plan_ptr)
                                  end do

                                  do i = local_start, local_end - 1
                                    do j = 1, nFreq
                                       
                                       if (j == 1) then
  
                                          fourierF90Array(i - local_start + 1,j) = cmplx(real(data_out(j,i - local_start + 1), kind=kind(0.0)), &
                                          aimag(data_out(j,i - local_start + 1)), kind=kind((0.0, 0.0)))
  
                                       else
  
                                          fourierF90Array(i - local_start + 1,j) = 2 * cmplx(real(data_out(j,i - local_start + 1), kind=kind(0.0)), &
                                          aimag(data_out(j,i - local_start + 1)), kind=kind((0.0, 0.0)))
                                          
                                       end if
                                    end do
                                  end do
                                  ! if (rank == 0) then

                                  !   print *, fourierF90Array(1,:)
                                  !   stop
                                  ! end if
                                  deallocate(data_out,data_in)

                                end subroutine PerformDFT
                            
                             
            
                                Subroutine Load_Snapshot
                                    use mpi
                                    implicit none


                                    character(len=300) :: filename
                                    integer ::  status(MPI_STATUS_SIZE), fh
                                    integer :: globalRows, globalCols
                                    integer :: myRows, remainder, start_row, end_row
                                    integer :: i
                                    double precision, allocatable :: localData(:,:)
                                    integer :: nprocs

                                    integer, dimension(2) :: sizes, subsizes, starts
                                    integer :: subarrayType
                                    integer(kind=MPI_OFFSET_KIND) :: offseta

                                    PetscScalar, pointer :: local_x_new(:,:)

                                    offseta = 0_MPI_OFFSET_KIND

                                    call MPI_Comm_rank(PETSC_COMM_WORLD, rank, ierr)
                                    call MPI_Comm_size(PETSC_COMM_WORLD, nprocs, ierr)


                                    write(filename, filenameFormat) currentSnapshotNumber
                                    call MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_RDONLY, MPI_INFO_NULL, fh, ierr)
                                    if (ierr /= 0) then
                                      print *, 'Error opening file on rank', rank, filename
                                      call MPI_Abort(MPI_COMM_WORLD, ierr, ierr)
                                    end if
              
                                    globalRows = numberOfDataInSnapShot
                                    globalCols = 1

                                    ! Partition the rows among processes.
                                    myRows = globalRows / nprocs
                                    remainder = mod(globalRows, nprocs)
                                    if (rank < remainder) then
                                      myRows = myRows + 1
                                      start_row = rank * myRows
                                    else
                                      start_row = rank * myRows + remainder
                                    end if
                                    end_row = start_row + myRows - 1
                              
                                    
                                    allocate(localData(myRows, globalCols))

                                    ! Create an MPI subarray type that represents a block of rows.
                                    sizes    = (/ globalRows, globalCols /)   ! Full matrix dimensions (93685 x 8)
                                    subsizes = (/ myRows, globalCols /)         ! Block dimensions: myRows x 8
                                    starts   = (/ start_row, 0 /)               ! Starting indices (0-based) in each dimension
                                    call MPI_Type_create_subarray(2, sizes, subsizes, starts, MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, subarrayType, ierr)
                                    call MPI_Type_commit(subarrayType, ierr)
                                  
                                    ! Set the file view. With the derived type the file “view” is the subarray of rows.
                                    call MPI_File_set_view(fh, offseta, MPI_DOUBLE_PRECISION, subarrayType, 'native', MPI_INFO_NULL, ierr)
                                  
                                    ! Read the subarray into localData.
                                    call MPI_File_read_all(fh, localData, myRows*globalCols, MPI_DOUBLE_PRECISION, status, ierr)
                                    if (ierr /= 0) then
                                      print *, 'Error reading file on rank', rank
                                      call MPI_Abort(MPI_COMM_WORLD, ierr, ierr)
                                    end if
                                  
                                    call MPI_File_close(fh, ierr)
                              
                                    call MatDenseGetArrayF90(x_new, local_x_new, ierr)
                                    do i = 1, myRows
                                    
                                      local_x_new(i, 1) = cmplx(localData(i, 1), kind=8)

                                    end do

                              

                                    call MatDenseRestoreArrayF90(x_new, local_x_new, ierr)
                                    call MatAssemblyBegin(x_new, MAT_FINAL_ASSEMBLY, ierr)
                                    call MatAssemblyEnd(x_new, MAT_FINAL_ASSEMBLY, ierr)
                                    deallocate(localData)
                            
                                    call MPI_Barrier(PETSC_COMM_WORLD, ierr)
                                    call handle_Latest_Snapshot
                                End Subroutine Load_Snapshot
              
                                
                                subroutine handle_Latest_Snapshot()
                              
                                  implicit none
                            
                                  PetscScalar, pointer :: local_x_new(:,:)
                                  PetscScalar, pointer :: local_x_sum(:,:)
                                  PetscScalar, pointer :: local_mu(:,:)

                                  PetscInt :: localRows,localCols
                                  PetscInt :: i

                                  ! To decide which process has current DFT row
                                  logical :: i_own_row
                                  integer :: row_owner, global_owner
                                  integer :: global_row_idx, local_row_idx
                                  PetscInt :: local_start, local_end
                                  integer ::  j, ld
                                  PetscInt :: rstart, rend
                                  PetscScalar, allocatable ::  global_row_data(:)

                                  PetscScalar :: alpha, beta

                                  integer ::  min_t_idx


                                  if (rank == 0) then
                                    t_i = t_i + 1
                
                                    call MPI_Bcast(t_i, 1, MPI_INT, 0, PETSC_COMM_WORLD, ierr)
                                        if (ierr /= 0) then
                                           print *, 'Error broadcasting matrix size.'
                                           call PetscFinalize(ierr)
                                           call SlepcFinalize(ierr)
                                           stop
                                        end if
                                  else 
                                    call MPI_Bcast(t_i, 1, MPI_INT, 0, PETSC_COMM_WORLD, ierr)
                                        if (ierr /= 0) then
                                          print *, 'Error broadcasting matrix size.'
                                          call PetscFinalize(ierr)
                                          call SlepcFinalize(ierr)
                                          stop
                                        end if
                                  end if
                                  
                                  call MatGetLocalSize(mu, localRows, localCols, ierr)

                                  call MatDenseGetArrayF90(mu, local_mu, ierr)
                                  call MatDenseGetArrayF90(x_new, local_x_new, ierr)
                                  alpha = t_i - 1.0d0
                                  beta = 1.0d0 / t_i

                                  do i = 1, localRows
                                    local_mu(i, 1) = beta * (alpha * local_mu(i, 1) + local_x_new(i, 1))
                                  end do

                                  call MatDenseRestoreArrayF90(mu, local_mu, ierr)
                                  call MatDenseRestoreArrayF90(x_new, local_x_new, ierr)
                



                                  if (t_idx(1) > 0) then

                                    call MatGetOwnershipRange(FourierRange, local_start, local_end, ierr)
        
                                    global_row_idx = t_idx(1) - 1

                                    i_own_row = (global_row_idx >= local_start) .and. (global_row_idx < local_end)

                                    allocate(global_row_data(nFreq))
                                    if (i_own_row) then
                                        row_owner = rank
                                    else
                                        row_owner = -1
                                    end if
                                    call MPI_Allreduce(row_owner, global_owner, 1, MPI_INTEGER, MPI_MAX, PETSC_COMM_WORLD, ierr)

                                    if (i_own_row) then
                                      local_row_idx = global_row_idx - local_start + 1
                                      global_row_data = fourierF90Array(local_row_idx, :)
                                    else
                                        global_row_data = (0.0_C_DOUBLE, 0.0_C_DOUBLE)
                                    end if

                                    call MPI_Bcast(global_row_data, 2*nFreq, MPI_DOUBLE_PRECISION, global_owner, PETSC_COMM_WORLD, ierr)
                                    call MPI_Barrier(PETSC_COMM_WORLD, ierr)

                                    call MPI_Allreduce(row_owner, global_owner, 1, MPI_INTEGER, MPI_MAX, PETSC_COMM_WORLD, ierr)



                                    call MatDenseGetArrayF90(x_new, local_x_new, ierr)
                                    call MatDenseGetArrayF90(X_sum_BlockA, local_x_sum, ierr)
                                    
                                    call MatGetOwnershipRange(X_sum_BlockA, rstart, rend, ierr)
                                    ld = rend - rstart  

                                    do j = 1, nFreq
                                      do i = 1, ld
                                        local_x_sum(i, j) = local_x_sum(i, j) + (local_x_new(i, 1) * window(t_idx(1))) * global_row_data(j)
                                      end do
                                    end do

                                    call MatDenseRestoreArrayF90(x_new, local_x_new, ierr)
                                    call MatDenseRestoreArrayF90(X_sum_BlockA, local_x_sum, ierr)
                                    deallocate(global_row_data)


                                  end if


                                  if (nOvlp > 0) then
                                    if (t_idx(2) > 0) then

                                      call MatGetOwnershipRange(FourierRange, local_start, local_end, ierr)
          
                                      global_row_idx = t_idx(2) - 1
                                      i_own_row = (global_row_idx >= local_start) .and. (global_row_idx < local_end)

                                      allocate(global_row_data(nFreq))
                                      if (i_own_row) then
                                          row_owner = rank
                                      else
                                          row_owner = -1
                                      end if
                                      
                                      call MPI_Allreduce(row_owner, global_owner, 1, MPI_INTEGER, MPI_MAX, PETSC_COMM_WORLD, ierr)
                                      
                                      if (i_own_row) then
                                          ! Calculate local index of the global row
                                          local_row_idx = global_row_idx - local_start + 1
                                          global_row_data = fourierF90Array(local_row_idx, :)
                                      
                                      else
                                          global_row_data = (0.0_C_DOUBLE, 0.0_C_DOUBLE)
                                      end if

                                      call MPI_Bcast(global_row_data, 2*nFreq, MPI_DOUBLE_PRECISION, global_owner, PETSC_COMM_WORLD, ierr)
                                      call MPI_Barrier(PETSC_COMM_WORLD, ierr)

                                      
                                      call MatDenseGetArrayF90(x_new, local_x_new, ierr)
                                      call MatDenseGetArrayF90(X_sum_BlockB, local_x_sum, ierr)
                                      
                                      call MatGetOwnershipRange(X_sum_BlockB, rstart, rend, ierr)
                                      ld = rend - rstart  ! Number of locally owned rows
          
                                      do j = 1, nFreq
                                        do i = 1, ld
                                          local_x_sum(i, j) = local_x_sum(i, j) + (local_x_new(i, 1) * window(t_idx(2))) * global_row_data(j)
                                        end do
                                      end do
                                    
          
                                      call MatDenseRestoreArrayF90(x_new, local_x_new, ierr)
                                      call MatDenseRestoreArrayF90(X_sum_BlockB, local_x_sum, ierr)
                                      deallocate(global_row_data)
                                    end if

                                  end if
         

                                  if (t_idx(1) == nDFT) then      
                                    min_t_idx = minval(t_idx)
                                    t_idx(1) = min_t_idx - dn
                                    call handle_Latest_Block(1,X_sum_BlockA)
                                    call MatZeroEntries(X_sum_BlockA, ierr)
                                    call MatAssemblyBegin(X_sum_BlockA, MAT_FINAL_ASSEMBLY, ierr)        
                                    call MatAssemblyEnd(X_sum_BlockA, MAT_FINAL_ASSEMBLY, ierr)
                                  else
                                      t_idx(1) = t_idx(1) + 1
                                  end if

                                  if (nOvlp > 0) then
                                    if (t_idx(2) == nDFT) then
                                      min_t_idx = minval(t_idx)
                                      t_idx(2) = min_t_idx - dn
                                      call handle_Latest_Block(2,X_sum_BlockB)
                                      call MatZeroEntries(X_sum_BlockB, ierr)
                                      call MatAssemblyBegin(X_sum_BlockB, MAT_FINAL_ASSEMBLY, ierr)        
                                      call MatAssemblyEnd(X_sum_BlockB, MAT_FINAL_ASSEMBLY, ierr)
                                    else
                                        t_idx(2) = t_idx(2) + 1
                                    end if
                                  end if


                                end subroutine handle_Latest_Snapshot
                              
                            
                                subroutine handle_Latest_Block(blk_j,XHat)
                            
                                  implicit none
                                  
                                  ! =========================
                                  ! blk_j and XHat WILL BE PASSED BY handle_Latest_Snapshot
                                  ! =========================
                                  integer :: blk_j
                                  Mat, intent(inout) :: XHat 

                                  ! =========================
                                  ! LOCAL VARIABLES FOR MEAN REMOVAL
                                  ! =========================
                                  integer ::   ld
                                  PetscScalar, pointer :: local_b(:,:)
                                  PetscScalar, pointer :: local_xhat(:,:)
                                  PetscScalar, allocatable :: global_row_data(:)
                                  integer :: row_idx, i,j , rstart, rend, local_start, local_end
                                  PetscScalar :: scaling_factor
                                  logical :: i_own_row
                                  integer :: row_owner, global_owner
                                  integer :: global_row_idx, local_row_idx

                                  ! =========================
                                  ! LOCAL VARIABLES FOR FREQUENCY
                                  ! =========================
                                  integer :: fi
                                  PetscInt :: col_start_fi, rows_local_fi
                                  PetscInt ::  local_size
                                  PetscInt :: row, col
                                  PetscScalar :: valueOfCmatrixArray(1)
                                  PetscScalar, pointer :: u_hat_array(:,:)
                                  PetscScalar, pointer :: u_array_copy(:,:)
                                  PetscScalar, pointer :: x_hat_col_array(:)
                                  PetscScalar, pointer :: u_p_array(:)
                                  PetscScalar, pointer :: u_new_array(:)

                                  PetscReal abs_up
                                  PetscReal :: factor
                                  PetscScalar ::  scale_factor

                                  ! =========================
                                  ! MATRIX VARIABLES FOR FREQUENCY
                                  ! =========================
                                  Mat :: U
                                  Mat :: Kmode   

                                  ! =========================
                                  ! VECTORS VARIABLES FOR FREQUENCY
                                  ! =========================
                                  Vec :: X_hat_column
                                  Vec :: x       
                                  Vec :: u_p      
                                  Vec :: u_new 
                                  Vec :: Ux      
                                  Vec :: temp 
                                  Vec :: U_col
                                  Vec :: V_col

                                  ! =========================
                                  ! LOCAL VARIABLES FOR SVD
                                  ! =========================
                                  SVD :: svd
                                  PetscReal :: tolerance
                                  PetscInt :: maxIterations
                                  PetscInt :: nconverged, local_array_i
                                  PetscScalar, pointer :: array(:)
                                  PetscReal :: sigma
                                  Mat :: Umat

                                  
                                  ! =========================
                                  ! LOCAL VARIABLES FOR RECONSTRUCTION OF SPOD MODE
                                  ! =========================
                                  PetscInt ::  n, m_local, n_local, vec_local_size
                                  Mat :: Utemp  
                                  Mat :: UtempMult 
                                  PetscInt :: cstart
                                  PetscScalar, pointer :: u_array(:,:), utemp_array(:,:)
                                  PetscScalar, pointer :: vec_array(:)
                                  PetscScalar, pointer ::  uhat_array(:,:)
                                  integer :: nrowsForUtemp


                                  ! =========================
                                  ! VIEWER
                                  ! =========================
                                  PetscViewer :: viewer


                                  current_Block = current_Block + 1
                                  tolerance = 1.0d-30
                                  maxIterations = 10000

                                  if (rank == 0) then
                                    print *, '-------------------------------------------'                                  
                                    print *, 'Processing Block'
                                    print *, '-------------------------------------------'                                  
                                  end if

                                  call MatDenseGetArrayF90(mu, local_b, ierr)

                                  do row_idx = 1, nDFT
                                    if (rank == 0) then
                                      print *, 'Processing Block Realization', row_idx
                                    end if
                                    call MatGetOwnershipRange(FourierRange, local_start, local_end, ierr)
        
                                    global_row_idx = row_idx - 1
    
                                    i_own_row = (global_row_idx >= local_start) .and. (global_row_idx < local_end)
    
                                    allocate(global_row_data(nFreq))
                                    if (i_own_row) then
                                        row_owner = rank
                                    else
                                        row_owner = -1
                                    end if
                                    
                                    call MPI_Allreduce(row_owner, global_owner, 1, MPI_INTEGER, MPI_MAX, PETSC_COMM_WORLD, ierr)
                                    
                                    if (i_own_row) then
                                        ! Calculate local index of the global row
                                        local_row_idx = global_row_idx - local_start + 1
                                        global_row_data = fourierF90Array(local_row_idx, :)
                                     
                                    else
                                        global_row_data = (0.0_C_DOUBLE, 0.0_C_DOUBLE)
                                    end if
    
                                    call MPI_Bcast(global_row_data, 2*nFreq, MPI_DOUBLE_PRECISION, global_owner, PETSC_COMM_WORLD, ierr)
                                    call MPI_Barrier(PETSC_COMM_WORLD, ierr)
                                  
                                    call MatDenseGetArrayF90(XHat, local_xhat, ierr)
                                    scaling_factor = -window(row_idx)
                                    call MatGetOwnershipRange(XHat, rstart, rend, ierr)
                                    ld = rend - rstart  ! Number of locally owned rows
            
                                    do j = 1, nFreq
                                      do i = 1, ld
                                        local_xhat(i, j) = local_xhat(i, j) + ((scaling_factor * local_b(i, 1)) * global_row_data(j))
                                      end do
                                    end do
                                    
                                    call MatDenseRestoreArrayF90(XHat, local_xhat, ierr)
          
                                    call MPI_Barrier(PETSC_COMM_WORLD, ierr)

                                    deallocate(global_row_data)

                                  end do

                                  call MatDenseRestoreArrayF90(mu, local_b, ierr)
                                  call MatDenseGetArrayF90(XHat, local_xhat, ierr)
                                  scaling_factor = cmplx(winWeight / real(nDFT), 0.0d0, kind=8)

                                  call MatGetOwnershipRange(XHat, rstart, rend, ierr)
                                  ld = rend - rstart  ! Number of locally owned rows

                                  do j = 1, nFreq
                                    do i = 1, ld
                                      local_xhat(i, j) = local_xhat(i, j) * scaling_factor
                                    end do
                                  end do

                                  call MatDenseRestoreArrayF90(XHat, local_xhat, ierr)
                                  call MatDenseGetArrayF90(XHat, local_xhat, ierr)
                                  if (rank == 0) then
                                    print *, '-------------------------------------------'                                  
                                    print *, 'DONE : Processing Block'
                                    print *, '-------------------------------------------'                                  
                                  end if


                                  ! =========================
                                  ! FREQUENCY VARIABLE INITIALISE
                                  ! =========================
                                  call MatCreate(PETSC_COMM_WORLD, U, ierr)
                                  call MatSetSizes(U, PETSC_DECIDE, PETSC_DECIDE, numberOfDataInSnapShot, k, ierr)
                                  call MatSetType(U, MATDENSE, ierr)
                                  call MatSetUp(U, ierr)
                                  
                                  call VecCreate(PETSC_COMM_WORLD, X_hat_column, ierr)
                                  call VecSetSizes(X_hat_column, PETSC_DECIDE, numberOfDataInSnapShot, ierr)
                                  call VecSetType(X_hat_column, VECSTANDARD, ierr)

                                  call VecCreate(PETSC_COMM_WORLD, x, ierr)
                                  call VecSetSizes(x, PETSC_DECIDE, numberOfDataInSnapShot, ierr)
                                  call VecSetType(x, VECSTANDARD, ierr)

                                  call VecCreate(PETSC_COMM_WORLD, u_new, ierr)
                                  call VecSetSizes(u_new, PETSC_DECIDE, numberOfDataInSnapShot, ierr)
                                  call VecSetType(u_new, VECSTANDARD, ierr)

                                  call VecCreate(PETSC_COMM_WORLD, u_p, ierr)
                                  call VecSetSizes(u_p, PETSC_DECIDE, numberOfDataInSnapShot, ierr)
                                  call VecSetType(u_p, VECSTANDARD, ierr)

                                  call VecCreate(PETSC_COMM_WORLD, Ux, ierr)
                                  call VecSetSizes(Ux, PETSC_DECIDE, k, ierr)
                                  call VecSetType(Ux, VECSTANDARD, ierr)
          
                                  call VecCreate(PETSC_COMM_WORLD, temp, ierr)
                                  call VecSetSizes(temp, PETSC_DECIDE, numberOfDataInSnapShot, ierr)
                                  call VecSetType(temp, VECSTANDARD, ierr)

                                  call MatCreate(PETSC_COMM_WORLD, Kmode, ierr)
                                  call MatSetSizes(Kmode, PETSC_DECIDE, PETSC_DECIDE, k + 1, k + 1, ierr)
                                  call MatSetType(Kmode, MATAIJ, ierr)
                                  call MatSetUp(Kmode, ierr)

                                  call MatCreate(PETSC_COMM_WORLD, Umat, ierr)
                                  call MatSetSizes(Umat, PETSC_DECIDE, PETSC_DECIDE, k + 1, k + 1, ierr)
                                  call MatSetType(Umat, MATDENSE, ierr)
                                  call MatSetUp(Umat, ierr)

                                  call MatCreate(PETSC_COMM_WORLD, Utemp, ierr)
                                  call MatSetSizes(Utemp, PETSC_DECIDE, PETSC_DECIDE, numberOfDataInSnapShot, k+1, ierr)
                                  call MatSetType(Utemp, MATDENSE, ierr)
                                  call MatSetUp(Utemp, ierr)
            
                                  call MatCreate(PETSC_COMM_WORLD, UtempMult, ierr)
                                  call MatSetSizes(UtempMult, PETSC_DECIDE, PETSC_DECIDE, numberOfDataInSnapShot, k+1, ierr)
                                  call MatSetType(UtempMult, MATDENSE, ierr)
                                  call MatSetUp(UtempMult, ierr)

                                  ! =========================
                                  ! FREQUENCY CALCULATIONS
                                  ! =========================

                                  if (rank == 0) then
                                    print *, '-------------------------------------------'                                  
                                    print *, 'Processing Each Frequency'
                                    print *, '-------------------------------------------'                                  
                                  end if
                                  call MatDenseGetArrayF90(XHat, local_xhat, ierr)


                                  do fi = 1, nFreq
                                    if (rank == 0) then
                                      print *, 'Processing Frequency', fi
                                    end if
                                  ! =========================
                                  ! Get U from U_hat and X column from Xhat realization
                                  ! =========================
                                    call MatGetOwnershipRange(U_hat, rstart, rend, ierr)
                                    rows_local_fi = rend - rstart
                                    col_start_fi = (fi-1) * k
          
                                    call MatDenseGetArrayF90(U_hat, u_hat_array, ierr)
                                    call MatDenseGetArrayF90(U, u_array_copy, ierr)
                                    
                                    do j = 1, k
                                      do i = 1, rows_local_fi
                                        u_array_copy(i, j) = u_hat_array(i, col_start_fi + j)
                                      end do
                                    end do
                                    
                                    call MatDenseRestoreArrayF90(U_hat, u_hat_array, ierr)
                                    call MatDenseRestoreArrayF90(U, u_array_copy, ierr)
                                    call MatAssemblyBegin(U, MAT_FINAL_ASSEMBLY, ierr)        
                                    call MatAssemblyEnd(U, MAT_FINAL_ASSEMBLY, ierr)

                                    call VecGetArrayF90(X_hat_column, x_hat_col_array, ierr)
                                    call MatGetOwnershipRange(XHat, rstart, rend, ierr)
                                    ld = rend - rstart 
                                    
                                    do i = 1, ld
                                      x_hat_col_array(i) = local_xhat(i, fi)
                                    end do
                                    
                                    call VecRestoreArrayF90(X_hat_column, x_hat_col_array, ierr)
        

                                    call VecPointwiseMult(x, X_hat_column, sqrtW, ierr)
                                    call MatMultHermitianTranspose(U, x, Ux, ierr)
                                    call MatMult(U, Ux, temp, ierr)
                                    call VecWAXPY(u_p, (-1.0d0, 0.0d0), temp, x, ierr) 
                                    call VecNorm(u_p, NORM_2, abs_up, ierr)

                                    call VecGetArrayF90(u_p, u_p_array, ierr)
                                    call VecGetArrayF90(u_new, u_new_array, ierr)
                                    call VecGetLocalSize(u_new, local_size, ierr)
          
                                    do i = 1, local_size
                                        u_new_array(i) = u_p_array(i)
                                    end do
                                    
                                    call VecRestoreArrayF90(u_p, u_p_array, ierr)
                                    call VecRestoreArrayF90(u_new, u_new_array, ierr)
          
                                    scale_factor = cmplx(1.0d0 / abs_up, 0.0d0, kind=8)
                                    call VecScale(u_new,scale_factor , ierr)


                                  ! ===================================================
                                  ! CONSTRUCT K MATRIX
                                  ! ===================================================
                                    factor = sqrt(real(current_Block) / real((current_Block + 1)**2)) * sqrt(real(current_Block + 1))
    
                                    call MatGetOwnershipRange(S_hat, rstart, rend, ierr)
            
                                    do row = rstart, rend - 1
                                        if (row < k) then
                                          do col = 0, k - 1
                                            if (row == col) then
                                            call MatGetValues(S_hat, 1, [row], 1, [fi - 1], valueOfCmatrixArray, ierr)
                                            
                                              call MatSetValue(Kmode, row, row, factor * valueOfCmatrixArray(1), INSERT_VALUES, ierr)
                                            else
                                              call MatSetValue(Kmode, row, col, cmplx(0.0d0, 0.0d0, kind=8), INSERT_VALUES, ierr)
                    
                                            end if
                                          end do
                                          call VecGetValues(Ux, 1, [row], valueOfCmatrixArray, ierr)
                                         
                                          scale_factor = sqrt(real(current_Block) / real((current_Block + 1)**2)) * valueOfCmatrixArray(1)
                                          call MatSetValue(Kmode, row, k,scale_factor , INSERT_VALUES, ierr)
                                        else if (row == k) then
                                          do col = 0, k - 1
                                              call MatSetValue(Kmode, row, col,cmplx(0.0d0, 0.0d0, kind=8), INSERT_VALUES, ierr)
                                          end do
                                          scale_factor = sqrt(real(current_Block) / real((current_Block + 1)**2)) * abs_up
                                          call MatSetValue(Kmode, row, k,scale_factor , INSERT_VALUES, ierr)
                                        end if
                                    end do
            
                                    if ((rend) == k) then
                                      do col = 0, k - 1
                                        call MatSetValue(Kmode, rend, col,cmplx(0.0d0, 0.0d0, kind=8), INSERT_VALUES, ierr)
                                      end do
                                      scale_factor = sqrt(real(current_Block) / real((current_Block + 1)**2)) * abs_up
                                      call MatSetValue(Kmode, rend, k,scale_factor , INSERT_VALUES, ierr)
            
                                    end if
                                    call MatAssemblyBegin(Kmode, MAT_FINAL_ASSEMBLY, ierr)
                                    call MatAssemblyEnd(Kmode, MAT_FINAL_ASSEMBLY, ierr)
            
                                  
            
                                    ! call PetscViewerASCIIGetStdout(PETSC_COMM_WORLD, viewer, ierr)
                                    ! call MatView(Kmode, viewer, ierr)
                    
                                  ! ===================================================
                                  ! PERFORM SVD USING SLEPC
                                  ! ===================================================
                                    call MatCreateVecs(Kmode, U_col, V_col, ierr)
                                    call VecAssemblyBegin(U_col, ierr)
                                    call VecAssemblyEnd(U_col, ierr)

                                    call SVDCreate(PETSC_COMM_WORLD, svd, ierr)
                                    if (ierr /= 0) then
                                      print *, "Error in SVDCreate: ", ierr
                                      stop
                                    end if
                                    call SVDSetTolerances(svd,tolerance,maxIterations,ierr)
                                    if (ierr /= 0) then
                                      print *, "Error in SVDSetTolerances: ", ierr
                                      stop
                                    end if
                                    PetscCallA(SVDSetOperators(svd,Kmode,PETSC_NULL_MAT,ierr))
                                    PetscCallA(SVDSetProblemType(svd,SVD_STANDARD,ierr))
                                    PetscCallA(SVDSetType(svd,SVDLAPACK ,ierr))
                                    call SVDSetFromOptions(svd, ierr)
                                    call SVDSolve(svd, ierr)
                                    call SVDGetConverged(svd, nconverged, ierr)

                                  ! ===================================================
                                  ! EXTRACT EIGEN VALUES AND VECTORS
                                  ! ===================================================

                                    do i = 0, nconverged - 1
                                      call PetscViewerASCIIGetStdout(PETSC_COMM_WORLD, viewer, ierr)
                                      call SVDGetSingularTriplet(svd, i, sigma, U_col, V_col, ierr)


                                  ! ===================================================
                                  ! UPDATE SHAT WITH LATEST EIGEN VALUES
                                  ! ===================================================
                                      call MatGetOwnershipRange(S_hat, rstart, rend, ierr)
                                      if (rank == 0) then
                                        print*, sigma
                                      end if
                                      do row = rstart, rend - 1
                                        if (row == i) then
                                            call MatSetValue(S_hat, row, fi - 1, cmplx(sigma, 0.0d0, kind=8), INSERT_VALUES, ierr)
                                        end if
                                      end do
                                      call MatAssemblyBegin(S_hat, MAT_FINAL_ASSEMBLY, ierr)
                                      call MatAssemblyEnd(S_hat, MAT_FINAL_ASSEMBLY, ierr)

                                  ! ===================================================
                                  ! STORE EIGEN VECTORS IN UMAT
                                  ! ===================================================
                                      call VecGetOwnershipRange(U_col, rstart, rend, ierr)
                    
                                      allocate(array(rend - rstart))
                                      call VecGetArrayF90(U_col, array, ierr)
                    
                                      local_array_i = 1
                                      do row = rstart, rend - 1
                                
                                      call MatSetValue(Umat, row, i, array(local_array_i), INSERT_VALUES, ierr) 
                                      local_array_i = local_array_i + 1
                                    
                                      end do
                                      call VecRestoreArrayF90(U_col, array, ierr)
                
                                    end do

                                    call MatAssemblyBegin(Umat, MAT_FINAL_ASSEMBLY, ierr)
                                    call MatAssemblyEnd(Umat, MAT_FINAL_ASSEMBLY, ierr)
                                    call SVDDestroy(svd, ierr)





                                  ! ===================================================
                                  ! RECONSTRUCT SPOD MODES AND ASSIGN TO GLOBAL MATRIX
                                  ! ===================================================
                                    call MatGetOwnershipRange(Utemp, rstart, rend, ierr)
                                    call MatGetLocalSize(U, m_local, n_local, ierr)  
                                    nrowsForUtemp = rend - rstart
                                    n = k

                                    call MatDenseGetArrayF90(U, u_array, ierr)
                                    call MatDenseGetArrayF90(Utemp, utemp_array, ierr)

                                    utemp_array(1:m_local, 1:n) = u_array(1:m_local, 1:n)

                                    call VecGetArrayReadF90(u_new, vec_array, ierr)
                                    call VecGetLocalSize(u_new, vec_local_size, ierr)

                                    utemp_array(1:vec_local_size, n+1) = vec_array(1:vec_local_size)

                                    call VecRestoreArrayReadF90(u_new, vec_array, ierr)
                                    call MatDenseRestoreArrayF90(U, u_array, ierr)
                                    call MatDenseRestoreArrayF90(Utemp, utemp_array, ierr)

                                    call MatAssemblyBegin(Utemp, MAT_FINAL_ASSEMBLY, ierr)
                                    call MatAssemblyEnd(Utemp, MAT_FINAL_ASSEMBLY, ierr)

                                    call MatMatMult(Utemp, Umat, MAT_REUSE_MATRIX, PETSC_DEFAULT_REAL, UtempMult, ierr)

                                    call MatGetOwnershipRange(UtempMult, rstart, rend, ierr)
                                    m_local = rend - rstart
                                    cstart = (fi - 1) * k

                                    call MatDenseGetArrayF90(UtempMult, utemp_array, ierr)
                                    call MatDenseGetArrayF90(U_hat, uhat_array, ierr)

                                    do j = 1, k
                                      uhat_array(1:m_local, cstart+j) = utemp_array(1:m_local, j)
                                    end do

                                    call MatDenseRestoreArrayF90(UtempMult, utemp_array, ierr)
                                    call MatDenseRestoreArrayF90(U_hat, uhat_array, ierr)
            
                                    call MatAssemblyBegin(U_hat, MAT_FINAL_ASSEMBLY, ierr)
                                    call MatAssemblyEnd(U_hat, MAT_FINAL_ASSEMBLY, ierr)


                                    call VecDestroy(U_col, ierr)
                                    call VecDestroy(V_col, ierr)


                                  end do


                                  call MatDenseRestoreArrayF90(XHat, local_xhat, ierr)
                                  call MPI_Barrier(MPI_COMM_WORLD, ierr)

                                  call MatDestroy(U, ierr)
                                  call MatDestroy(Kmode, ierr)
                                  call MatDestroy(Umat, ierr)
                                  call MatDestroy(Utemp, ierr)
                                  call MatDestroy(UtempMult, ierr)

                                  call VecDestroy(X_hat_column, ierr)
                                  call VecDestroy(x, ierr)
                                  call VecDestroy(u_p, ierr)
                                  call VecDestroy(u_new, ierr)
                                  call VecDestroy(Ux, ierr)
                                  call VecDestroy(temp, ierr)


                            
                                end subroutine handle_Latest_Block
                                
                                subroutine SPODMode()
                                  implicit none

                                  Vec :: invSqrtW

                                  call VecDuplicate(sqrtW, invSqrtW, ierr)  ! Create a new PETSc vector
                                  call VecCopy(sqrtW, invSqrtW, ierr)       ! Copy sqrtW to invSqrtW
                                  call VecReciprocal(invSqrtW, ierr)        ! Compute 1/sqrt(w)
            
                                  call MatDiagonalScale(U_hat, invSqrtW, PETSC_NULL_VEC, ierr)  ! Element-wise row division
                                  call MatAssemblyBegin(U_hat, MAT_FINAL_ASSEMBLY, ierr)
                                  call MatAssemblyEnd(U_hat, MAT_FINAL_ASSEMBLY, ierr)
    
                                end subroutine SPODMode

                                subroutine SaveEnergySpectraAndSPODMode()
                                  implicit none

                                end subroutine SaveEnergySpectraAndSPODMode

                                ! Subroutine to find and save the z highest SPOD modes
                                ! Subroutine to find and save the z highest SPOD modes
                                subroutine SaveHighestSPODModes()
                                  implicit none
                                  
                                  integer :: z
                                  
                                  PetscInt :: i, j, nrows, rank, size
                                  PetscInt :: cstart, cend
                                  PetscInt, allocatable :: indices(:)
                                  PetscScalar, allocatable :: values(:), local_value
                                  Vec :: XSPOD_VEC
                                  PetscViewer :: viewer
                                  character(len=300) :: filename
                                  PetscInt :: temp_idx
                                  PetscInt :: actual_col

                                  call MPI_Comm_rank(PETSC_COMM_WORLD, rank, ierr)
                                  call MPI_Comm_size(PETSC_COMM_WORLD, size, ierr)

                                  z = nFreqToSave
                                  ! Get MPI rank and size
                                  call PetscViewerASCIIOpen(PETSC_COMM_WORLD, 'EnegySpectra.txt', viewer, ierr)
                                  call MatView(S_hat, viewer, ierr)
                                  call PetscViewerDestroy(viewer, ierr)
                                  if (rank == 0) then
                                    print *, 'Saved Enegy spectra to File EnegySpectra.txt '
                                  end if

                                  
                                  
                                  call MatGetSize(S_hat, nrows, nFreq, ierr)
                                 
                                  allocate(values(nFreq))
                                  allocate(indices(nFreq))
                                  values = 0.0
                                  
                              
                                  if (rank == 0) then
                                    call MatGetOwnershipRangeColumn(S_hat, cstart, cend, ierr)
                                    do i = 1, nFreq
                                      indices(i) = i - 1  ! 0-based indices for PETSc
                                    end do
                                   
                                    do j = 1, nFreq
                                        call MatGetValues(S_hat, 1, [0], 1, [j-1], values(j), ierr)

                                    end do
                                 
                                    do i = 1, nFreq
                                      do j = 1, nFreq - i
                                        if (abs(values(j)) < abs(values(j+1))) then
                                          ! Swap values
                                          local_value = values(j)
                                          values(j) = values(j+1)
                                          values(j+1) = local_value
                                          
                                          ! Swap indices
                                          temp_idx = indices(j)
                                          indices(j) = indices(j+1)
                                          indices(j+1) = temp_idx
                                        end if
                                      end do
                                    end do



                                  end if

                                 
                                  call MPI_Barrier(PETSC_COMM_WORLD, ierr)


                                  ! Broadcast the sorted indices and values from rank 0 to all processors
                                  call MPI_Bcast(values, nFreq, MPI_DOUBLE_PRECISION, 0, PETSC_COMM_WORLD, ierr)
                                  call MPI_Bcast(indices, nFreq, MPI_INTEGER, 0, PETSC_COMM_WORLD, ierr)
                                  
                                  ! Now all processors have the same sorted list
                                  
                                  ! Create a vector to hold the mode
                                  call MatCreateVecs(U_hat, PETSC_NULL_VEC, XSPOD_VEC, ierr)
                                  if (ierr /= 0) then
                                    print *, "Error in XSPOD_VEC MatCreateVecs : ", ierr
                                    return
                                  end if
                                  
                                  ! Save the top z modes
                                  do i = 1, z
                                    ! Get the column vector from U_hat corresponding to the ith highest mode
                                    actual_col = indices(i) * k 

                                    call MatGetColumnVector(U_hat, XSPOD_VEC, actual_col, ierr)
                                    if (ierr /= 0) then
                                      print *, "Error in MatGetColumnVector: ", ierr
                                      return
                                    end if
                                    
                                    ! Create filename with mode index and corresponding sigma value
                                    write(filename, '(A,A,I0,A,E12.5,A)') trim('SPODMODE'), "_freq_", i, "_sigma_", abs(values(i)), ".dat"
                                    
                                    ! Save the vector to a binary file
                                    call PetscViewerBinaryOpen(PETSC_COMM_WORLD, trim(filename), FILE_MODE_WRITE, viewer, ierr)
                                    if (ierr /= 0) then
                                      print *, "Error in PetscViewerBinaryOpen: ", ierr
                                      return
                                    end if
                                    
                                    call VecView(XSPOD_VEC, viewer, ierr)
                                    if (ierr /= 0) then
                                      print *, "Error in VecView: ", ierr
                                      return
                                    end if
                                    
                                    call PetscViewerDestroy(viewer, ierr)
                                    
                                    ! Print information about the saved mode (only from rank 0 to avoid duplicate messages)
                                    if (rank == 0) then
                                      print *, 'Eigen values at index ',indices(i) , ' Globally in Uhat at ', actual_col
                                      print '(A,I0,A,E12.5,A,A)', "Saved SPOD mode ", i, " with sigma value ", abs(values(i)), " to file: ", trim(filename)
                                    end if
                                  end do
                                  
                                  call VecDestroy(XSPOD_VEC, ierr)
                                  deallocate(values)
                                  deallocate(indices)
                                  
                                end subroutine SaveHighestSPODModes

                                subroutine performStreamingSPOD()
                                  implicit none

                                  integer :: i

                                  do i = 1,numberOfSnapShot
                                    call MPI_Barrier(PETSC_COMM_WORLD, ierr)
                                      if (rank == 0) then
                                          print *, '-------------------------------------------'                                  
                                          print *, 'Processing Snapshot', i
                                          print *, '-------------------------------------------'                                  
                                      end if
                                    currentSnapshotNumber = i
                                    call Load_Snapshot
                                end do

                                end subroutine performStreamingSPOD
                            
                              End Module streamingSPODModule
                              
                              Program streamingSPOD
#include <slepc/finclude/slepcsvd.h>
                            
                                use slepc
                                use slepcsvd
                                use streamingSPODModule
                            
                                implicit none
                            

                                call SlepcInitialize(ierr)
                                if (ierr /= 0) then
                                  print *, 'Error initializing SLEPc.'
                                  call SlepcFinalize(ierr)
                                  stop
                                end if
                                ! Get the rank of the calling process in PETSC_COMM_WORLD
                                call MPI_Comm_rank(PETSC_COMM_WORLD, rank, ierr)

                                ! Get the total number of processes in PETSC_COMM_WORLD
                                call MPI_Comm_size(PETSC_COMM_WORLD, size, ierr)

                                ! Print the number of processes from rank 0
                                if (rank == 0) then
                                  print *, "Number of MPI processes (size) =", size
                                end if

                                call UserInput
                                call InitializeBlocks
                                call ReadWeightMatrix
                                call WindowingFunction
                                call PerformDFT
                                call performStreamingSPOD
                                call MPI_Barrier(PETSC_COMM_WORLD, ierr)
                                call SPODMode
                                call SaveHighestSPODModes
                          
                               
                                call SlepcFinalize(ierr)

                              End Program streamingSPOD
                            
                            
                    
                    
