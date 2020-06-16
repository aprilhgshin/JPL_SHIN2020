c With Python, create initial temperature and location mask and save into a binary file
c With Fortran, load the binary files
c In a header file, initialize the domain and tile dimensions and number of processes
c With MPI_SEND, send a piece of the domain to each process through a do loop
c Have each process MPI_RECV and do their job: listed in notepad
c After process does what it should to the masks, output a .bin file of each mask it receives and creates

c 1. 4 processor job.  each processor will have a 10x10 array of floats, representing temperature
c 2. the "domain" will be 40x10,
c 3. a 40x10 array of ints will define the 'mask',  all zeros except for 1s through the 20th row, columns 5:25
c 4. a 40x10 array of floats defining the 'temperature' will be read in and distributed to the processors,
c    each will get 10x10
c 5. the 40x10 mask array will be read in and distributed to the processes
c 6. each process will determine whether they have any of the mask values in them
c 7. where they do have the mask they will print out the 'temperature' at the mask points --> save all to .bin files
c 8. define temperature as T(row,col) = row*col

c Wednesday, JUNE 17:
c TODO: output bin files - construct corresponding filenames and write data
c TODO: research for other read functions that specify dtype of data read
c TODO: accomodate for processes along y-axis in loop
c TODO: organize/clean code and better documentation

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c234567890


      program create_output

      include "SIZE.h"
      include "mpif.h"
      !include "/home/mitgcm/Work/MITgcm/pkg/mdsio/MDSIO_OPTIONS.h"

      integer, dimension(Ny,Nx) :: mask
      real, dimension(Ny,Nx) :: temp
      real, dimension(Ny,Nx) :: final_mask

c     Arrays to check if mask and temp are distributed to processes correctly:
c     -1 for not being sent, 0,1,2,3 for corresponding process to which piece of array was sent
      integer, dimension(Ny, Nx) :: mask_correct
      integer, dimension(Ny, Nx) :: temp_correct
      integer, dimension(Ny, Nx) :: final_correct

c     In Parent Process 0: Used as arrays to store pieces of mask and temp to send to each process
c     All Processes 0-3: Used as buffer arrays to store pieces of mask and temp each process receives
      integer, dimension(sNy,sNx) :: sub_mask
      real, dimension(sNy, sNx) :: sub_temp

c     Store temperature values only on elements == 1 on each sub_mask
      real, dimension(sNy, sNx) :: sub_tempOnMask

c      character(len=14) :: filename
      integer row, col, p_count, p_x, p_y, ierror, proc_id, num_ones
      integer status(MPI_STATUS_SIZE), final_col, correct
      character*18 filename
      integer, parameter :: debug = 0
      correct = 0

c     Initialize all arrays
        do row=1, Ny, 1
          do col=1, Nx, 1
            mask(row,col) = 0
            temp(row,col) = 0.0
            final_mask(row,col) = 0
            mask_correct(row,col) = 0
            temp_correct(row,col) = 0
            final_correct(row,col) = 0
            if ((row .le. sNy) .and. (col .le. sNx)) then
              sub_mask(row,col) = 0
              sub_temp(row,col) = 0
              sub_tempOnMask = 0.0
            end if
          end do
        end do

c     Read mask and temp data from binary files
       open(1, FILE="input_domains/domain_mask.bin", FORM="unformatted"
     & ,ACCESS="stream")
       read(1) mask

       open(2, FILE="input_domains/domain_temp.bin", FORM="unformatted"
     & ,ACCESS="stream")
       read(2) temp

c     Initialize MPI environment
      call MPI_INIT(ierror)

c     Return (current) process id of the process that called the function
      call MPI_COMM_RANK(MPI_COMM_WORLD, proc_id, ierror)


c     Parent Process' task:
      if (proc_id .eq. 0) then

c     NOTE: Tile number = process id. Therefore, tile numbers and process id's start from 0.
c           In this do loop: when separately counting tiles(procesess) along the x axis and y axis, denoted as p_x and p_y, we start from 0
c           in order to correctly compute the row and col of mask and temp corresponding to the row and col of sub_mask and sub_temp.
c     NOTE: However, the Fortran convention of starting all array indices from 1 is held.
c           Additionally, when writing to binary files, we index sub_masks and sub_temps by treating each of them
c           as a single massive grid cell in the full domain mask and temp. These indices start from 1 as per the Fortran convention.

c      Distribute pieces of mask and temp to each process
        do p_y=0, nPy-1, 1
          do p_x=0, nPx-1, 1

            do row=1, sNy, 1
              do col=1, sNx, 1
c               Adding p_x*sNx to col and p_y*sNy to row allows us to move onto next portion
c               of mask and temp to save into sub_mask and sub_temp for each process.
                sub_mask(row, col) = mask(int(row+p_y*sNy),
     & int(col+p_x*sNx))
                sub_temp(row, col) = temp(int(row+p_y*sNy),
     & int(col+p_x*sNx))
              end do
            end do

c           p_count here is equivalent to process id assigned to each tile in increasing order to the right and downwards
            p_count = p_x + (nPx*p_y)

            !MPI_SEND params: data_to_send, send_count, send_type, destination_ID, tag, comm, ierror
            call MPI_SEND(sub_temp, sNx*sNy, MPI_REAL, p_count, 200,
     & MPI_COMM_WORLD, ierror)
            call MPI_SEND(sub_mask, sNx*sNy, MPI_INT, p_count, 100,
     & MPI_COMM_WORLD, ierror)

c         Printing sub masks sent through MPI_SENT
          print *, "Sent to process", p_count, "with sub_mask ",
     & "dimensions: (", shape(sub_mask), ")"
          print *, "For column indices of full mask: ", 1+p_count*sNx,
     & "to ", 10+p_count*sNx, "and row indices: 1 to 10"

          end do
        end do
      end if

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

c     ALL Processes' task:

c     MPI_RECV params: received_data, receive_count, receive_type, sender_ID, tag, comm, status, ierr
      call MPI_RECV(sub_mask, sNx*sNy, MPI_INT, 0, 100, MPI_COMM_WORLD,
     & status, ierror)
      call MPI_RECV(sub_temp, sNx*sNy, MPI_REAL, 0, 200, MPI_COMM_WORLD
     & ,status, ierror)

      call VERIFY_SUBMASKS(sub_mask, sub_temp, mask, temp, proc_id,
     & mask_correct, temp_correct, correct)

      call CREATE_FILENAME(proc_id, 1, 1, filename)

      print *, "Process",proc_id, "Received from parent with tag",
     & status(MPI_TAG)

      num_ones = 0
      do row=1, sNy, 1
        do col=1, sNx, 1
          if (sub_mask(row,col) .eq. 1) then
            sub_tempOnMask(row,col) = sub_temp(row,col)
            num_ones = num_ones + 1
          else
            sub_tempOnMask(row,col) = 0.0
          end if
        end do
      end do
      print*, "Process", proc_id, "has", num_ones, "number of ones"

      if (num_ones .gt. 0) then
        !MPI_SEND params: data_to_send, send_count, send_type, destination_ID, tag, comm, ierror
        call MPI_SEND(sub_tempOnMask, sNx*sNy, MPI_REAL, 0, 1,
     & MPI_COMM_WORLD, ierror)
        print *, "process", proc_id, "sends tag 1 to 0"

      else
        call MPI_SEND(sub_tempOnMask, sNx*sNy, MPI_REAL, 0, 0,
     & MPI_COMM_WORLD, ierror)
        print *, "process", proc_id, "sends tag 0 to 0"
      end if

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

      if (proc_id .eq. 0) then

        do p_count=0, 3, 1
c         MPI_RECV params: received_data, receive_count, receive_type, sender_ID, tag, comm, status, ierr
          call MPI_RECV(sub_tempOnMask, sNx*sNy, MPI_REAL, p_count,
     & MPI_ANY_TAG ,MPI_COMM_WORLD, status, ierror)
          print *, "Parent received from process", status(MPI_SOURCE),
     & "with tag", status(MPI_TAG)

          final_col = 0
          if (status(MPI_TAG) .eq. 1) then
            print *, "Tag is 1. Time to append"
            do row=1, sNy, 1
              do col=1, sNx, 1
c             Adding p_count*sNx to col allows us to move onto next portion of mask and temp
c             to save into sub_mask and sub_temp for each process
                final_col = int(col+status(MPI_SOURCE)*sNx)
                final_mask(row, final_col) = sub_tempOnMask(row, col)
              end do
            end do

          else
            print *, "Tag is 0. Nothing to append to final_arr"
          end if

        end do
        call VERIFY_FINAL(mask, final_mask, final_correct)
        print *, "Final updated temperature array over mask: "
        call PRINT_FLOAT_ARR(final_mask, Ny*Nx)

        print *, "Mask array: "
        call PRINT_INT_ARR(mask, Ny*Nx)

      end if


      call MPI_FINALIZE(ierror)

      stop
      end

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

      SUBROUTINE PRINT_INT_ARR(array, arr_length)
      integer arr_length
      integer, dimension(arr_length) :: array
      do, i=1, arr_length
        write(*,"(I5,$)") array(i)
        !write(*,"(I1,$)") array(i)
          if (i .EQ. arr_length) then
             write(*, '(A,/)') ''
          endif
      enddo
      end


      SUBROUTINE PRINT_FLOAT_ARR(array, arr_length)
c     Edit Descriptor:
c       w: the number of positions to be used
c       m: the minimum number of positions to be used
c       d: the number of digits to the right of the decimal point
c       e: the number of digits in the exponent part

c       Reading/writing REALs Decimal form Fw.d
      integer arr_length
      real, dimension(arr_length) :: array
      do, i=1, arr_length
        write(*,"(F5.0,$)") array(i)
          if (i .EQ. arr_length) then
             write(*, '(A,/)') ''
          endif
      enddo
      end

      SUBROUTINE CREATE_FILENAME(proc_id, is_temp, is_child, filename)
        include 'SIZE.h'
        integer, intent(in) :: proc_id
        integer, intent(in) :: is_temp
        integer, intent(in) :: is_child
        character*18, intent(out) :: filename

        integer p_x, p_y
        character*4 str_px, str_py

        p_y = proc_id / nPx      ! Integer division does not round.
        p_x = mod(proc_id, nPx)  ! The remainder is equivalent to tile(process) number along x_axis.

        print *, "is child: ", is_child, "is temp: ", is_temp
        if (is_child .eq. 1) then
          if (p_x .lt. 10) then
            write(str_px, "(A3, I1)") "000", p_x
          else if (p_x .ge. 10) then
            write(str_px, "(A2, I2)") "00", p_x
          else if ((p_x .ge. 100) .and. p_x .le. 999) then
            write(str_px, "(A1, I3)") "0", p_x
          else if ((p_x .ge. 1000) .and. (p_x .le. 9999)) then
            write(str_px, "(I4)") p_x
          end if

          if (p_y .lt. 10) then
            write(str_py, "(A3, I1)") "000", p_y
          else if (p_y .ge. 10) then
            write(str_py, "(A2, I2)") "00", p_y
          else if ((p_y .ge. 100) .and. p_y .le. 999) then
            write(str_py, "(A1, I3)") "0", p_y
          else if ((p_y .ge. 1000) .and. (p_y .le. 9999)) then
            write(str_py, "(I4)") p_y
          end if

          if (is_temp .eq. 1) then
            write(filename, "(A2, A4, A1, A4, A4)") "T.",
     &  str_px(1:4), ".", str_py(1:4), ".bin"
          else
            write(filename, "(A5, A4, A1, A4, A4)") "mask.",
     &  str_px(1:4), ".", str_py(1:4), ".bin"
          end if

          print *, "filename:",filename

        else
          if (is_temp .eq. 1) then
            write(filename, "(A12)") "T_global.bin"
          else
            write(filename, "(A15)") "mask_global.bin"
          end if
          
        end if
      end


      SUBROUTINE VERIFY_SUBMASKS(sub_mask, sub_temp, mask, temp,
     & proc_id, mask_correct, temp_correct, correct)

c     NOTE: Tile number = process id. Therefore, tile numbers and process id's start from 0.
c           In this subroutine: when separately counting tiles(procesess) along the x axis and y axis, we start from 0
c           in order to correctly compute the row and col of mask and temp corresponding to the row and col of sub_mask and sub_temp
c     NOTE: However, the Fortran convention of starting all array indices from 1 is held.

        include "SIZE.h"

        integer, dimension(sNy,sNx) :: sub_mask
        real, dimension(sNy,sNx) :: sub_temp
        integer, dimension(Ny,Nx) :: mask
        real, dimension(Ny,Nx) :: temp
        integer, dimension(Ny,Nx) :: mask_correct
        integer, dimension(Ny,Nx) :: temp_correct
        integer proc_id, big_col, big_row, row, col
c       p_x is the process number in x axis, p_y is the tile(process) number in y axis
        integer p_x, p_y
        integer, intent(out) :: correct

c       Determining p_x and p_y using proc_id using integer division and modulus:
        p_y = proc_id / nPx      ! Integer division does not round.
        p_x = mod(proc_id, nPx)  ! The remainder is equivalent to tile(process) number along x_axis.

        big_col = 0
        big_row = 0

        do row=1, sNy, 1
          do col=1, sNx, 1
            big_col = int(col+p_x*sNx)
            big_row = int(row+p_y*sNy)
            if (sub_mask(row, col) .eq. mask(big_row, big_col)) then
              correct = 1
              mask_correct(big_row, big_col) = 1
            else
              correct = 0
              mask_correct(big_row, big_col) = 0
            end if

            if (sub_temp(row, col) .eq. temp(row, big_col)) then
              correct = 1
              mask_correct(big_row, big_col) = 1
            else
              correct = 0
              mask_correct(big_row, big_col) = 0
            end if

          end do
        end do

        if (correct .eq. 1) then
          print *, "Process", proc_id, "receieved correct data."
        else
          print *, "INCORRECTLY DISTRIBUTED FOR PROCESS", proc_id
          print *, "mask: 1 for correct. 0 for incorrect: "
          call PRINT_INT_ARR(mask_correct, Ny*Nx)
          print *, "temp: 1 for correct. 0 for incorrect: "
          call PRINT_INT_ARR(temp_correct, Ny*Nx)
        end if
      end


      SUBROUTINE VERIFY_FINAL(mask, final_mask, final_correct)
c       This subroutine is only valid under the assumption that the temperature will not be 0.0
        include "SIZE.h"
        integer, dimension(Ny,Nx) :: mask
        real, dimension(Ny,Nx) :: final_mask
        integer, dimension(Ny,Nx) :: final_correct
        integer row, col, correct
        correct = 0

        do row=1, Ny, 1
          do col=1, Nx, 1
            if (((mask(row,col) .eq. 1) .and.
     & (final_mask(row,col) .ge. 0)) .or.
     & ((mask(row,col) .eq. 0) .and. (final_mask(row,col) .eq. 0.0)))
     & then
              final_correct(row,col)=1
              correct = 1
            else
              print *, "INCORRECT FINAL"
              correct = 0
              final_correct(row,col)=0
            end if
          end do
        end do

        if (correct .ne. 1) then
          print *, "final_mask: 1 for correct, 0 for incorrect: "
          call PRINT_INT_ARR(final_correct, Ny*Nx)
        else
          print *, "Final mask is correct."
        end if
      end
