c External files:
c With Python, create initial temperature and location mask and save into a binary file.
c In a header file, initialize the domain and tile dimensions and number of processes.

c In this program:
c Open and read the binary files.c to extract global arrays, temp and mask
c With MPI_SEND, send a piece of the global arrays to each process.
c Have each process MPI_RECV and save temperature values only on mask points.
c Output a .bin file of each sub mask.
c Each process calls MPI SEND to send the masked sub temperature values to Parent process
c Parent combines all messages with tag 1 and include in global final mask. Write to binary file.
c NOTE: All binary files will be saved in the output_masks folder

c In more detail:
c 1. 4 processor job.  each processor will have a 10x10 array of floats, representing temperature
c 2. the "domain" will be 40x10,
c 3. a 40x10 array of ints will define the 'mask',  all zeros except for 1s through the 20th row, columns 5:25
c 4. a 40x10 array of floats defining the 'temperature' will be read in and distributed to the processors,
c    each will get 10x10
c 5. the 40x10 mask array will be read in and distributed to the processes
c 6. each process will determine whether they have any of the mask values in them
c 7. where they do have the mask they will print out the 'temperature' at the mask points --> save all to .bin files
c 8. define temperature as T(row,col) = row*col


C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c234567890


      program create_output

      include "SIZE.h"
      include "mpif.h"

      integer, dimension(Ny, Nx) :: mask
      real, dimension(Ny, Nx) :: temp
      real, dimension ((sNy+sNx)*(nPx*nPy)) :: temp_on_mask
      integer, dimension((sNy+sNx)*(nPx*nPy)) :: global_indices

c     In Parent Process 0: Used as arrays to store pieces of mask and temp to send to each process
c     All Processes 0-3: Used as buffer arrays to store pieces of mask and temp each process receives
      integer, dimension(sNy,sNx) :: sub_mask
      real, dimension(sNy, sNx) :: sub_temp

c     Store temperature values only on elements == 1 on each sub_mask
      real, dimension(sNy+sNx+1) :: sub_tempOnMask
      integer, dimension(sNy+sNx+1) :: sub_global_indices

      integer row, col, p_count, p_x, p_y, ierror, proc_id, num_ones
      integer status(MPI_STATUS_SIZE), final_col
      integer s_row, e_row, s_col, e_col, sub_count, sub_num_ones
      integer, parameter :: debug = 1
      integer, parameter :: print_MPI_SEND = 0
      integer, parameter :: print_MPI_RECV = 0

      character*10 output_file

      common /sub_arr_pid/ sub_mask, sub_temp, proc_id

c     Initialize arrays used by all processes
      do row=1, Ny, 1
        do col=1, Nx, 1
          mask(row,col) = 0
          if ((row .le. sNy) .and. (col .le. sNx)) then
            sub_mask(row,col) = 0
            sub_temp(row,col) = 0
          end if
        end do
      end do

      do col=1, sNy+sNx+1, 1
        sub_tempOnMask(col) = 0.0
        sub_global_indices = 0
      end do



c     Read mask and temp data from binary files
       open(1, FILE="input_domains/domain_mask.bin", FORM="unformatted"
     & ,ACCESS="stream")
       read(1) mask
       close(1)
c      For double checking:
c       call PRINT_INT_ARR(mask, Ny*Nx)
c       row = 5
c       do col=1, Nx, 1
c         print *, mask(row,col)
c       end do


c     Initialize MPI environment
      call MPI_INIT(ierror)

c     Return (current) process id of the process that called the function
      call MPI_COMM_RANK(MPI_COMM_WORLD, proc_id, ierror)


C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

c     Parent Process' task:
      if (proc_id .eq. 0) then

c       Initialize global temp array
        do row=1, Ny, 1
          do col=1, Nx, 1
            temp(row,col) = 0
          end do
        end do

c       Read temp binary file
        open(2, FILE="input_domains/domain_temp.bin", FORM="unformatted"
     &       ,ACCESS="stream")
        read(2) temp
        close(1)

c     Creating output file to store all print statements in file instead of printing to terminal
      output_file = "output.txt"
      open(10, FILE=output_file, FORM="formatted")

c       For double checking:
c       call PRINT_FLOAT_ARR(temp, Ny*Nx)
c       print *, temp(1,12)

        call TEST_TEMP_READ(temp)

c     NOTE: Tile number = process id. Therefore, tile numbers and process id's start from 0.
c           When separately counting tiles(procesess) along the x axis and y axis, denoted as p_x and p_y, we start from 0
c           in order to correctly compute the row and col of mask and temp corresponding to the row and col of sub_mask and sub_temp.
c     NOTE: However, the Fortran convention of starting all array indices from 1 is held.
c           Additionally, when writing to binary files, we index sub_masks and sub_temps by treating each of them
c           as a single massive grid cell in the full domain mask and temp. These indices start from 1 as per the Fortran convention.
c     NOTE: We do not have to loop throught the global arrays. We can simply index to get tiles from global arrays

c      Distribute pieces of mask and temp to each process

       do p_y=0, nPy-1, 1
         do p_x=0, nPx-1, 1
           do row=1, sNy, 1
             do col=1, sNx, 1
c               Adding p_x*sNx to col and p_y*sNy to row allows us to move onto next portion
c               of mask and temp to save into sub_mask and sub_temp for each process.
                sub_mask(row, col) = mask(int(row+p_y*sNy),
     &                                    int(col+p_x*sNx))
                sub_temp(row, col) = temp(int(row+p_y*sNy),
     &                                    int(col+p_x*sNx))
              end do
            end do

c           p_count here is equivalent to process id assigned to each tile in increasing order to the right and downwards
            p_count = p_x + (nPx*p_y)
            write(10, FMT=*) "Message SEND status from Parent to ",
     &                        "Process", p_count, ": "

c           Parent sends pieces of temp and mask to each process
            !MPI_SEND params: data_to_send, send_count, send_type, destination_ID, tag, comm, ierror
            call MPI_SEND(sub_temp, sNx*sNy, MPI_REAL, p_count, 200,
     &                    MPI_COMM_WORLD, ierror)
c            print *, "Parent sent sub_temp to process", p_count,
c     &               "with tag 200"
            write(10, FMT=*) "Parent sent sub_temp to process", p_count,
     &               "with tag 200"

            call MPI_SEND(sub_mask, sNx*sNy, MPI_INT, p_count, 100,
     &                    MPI_COMM_WORLD, ierror)

            write(10, FMT=*) "Parent sent sub_mask to process", p_count,
     &               "with tag 100"

            !call MPI_BARRIER(MPI_COMM_WORLD, ierror)

c           Printing more MPI_SENT details
            if (print_MPI_SEND .eq. 0) then
              write(10, FMT=*) "Parent sent to process", p_count,
     &              "with sub_mask dimensions: (", shape(sub_mask), ")"
              write(10, FMT=*) "with column indices of full mask: ",
     & 1+p_count*sNx, "to ", 10+p_count*sNx, "and row indices: 1 to 10"
              write(10, FMT=*) "  "
            end if

          end do
        end do

        write(10, FMT=*) "Parent completed distribution to each ",
     &                   "process."
        write(10, FMT=*) " "
        write(10, FMT=*) " "
        close(10)
      end if

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

c     ALL Processes' task:
c      open(11, FILE=output_file, status="old", position="append"
c     &       ,action="write")

c     Each process receives a sub_mask and sub_temp from Parent.
c     MPI_RECV params: received_data, receive_count, receive_type, sender_ID, tag, comm, status, ierr
      call MPI_RECV(sub_mask, sNx*sNy, MPI_INT, 0, 100, MPI_COMM_WORLD,
     &              status, ierror)
      if (print_MPI_RECV .eq. 1) then
        print*,"Process",proc_id, "Received from parent "
     &                ,"with tag", status(MPI_TAG)
      end if

      call MPI_RECV(sub_temp, sNx*sNy, MPI_REAL, 0, 200, MPI_COMM_WORLD
     &             ,status, ierror)
      if (print_MPI_RECV .eq. 1) then
        print *,"Process",proc_id, "Received from parent with"
     &      ," tag" ,status(MPI_TAG)
      end if

c     Verify if each proess received the correct sub_mask and sub_temp
      call TEST_SUBMASKS(mask)
c     Write all sub_mask and sub_temp to binary file
      call WRITE_SUB_BIN()

      p_y = proc_id / nPx
      p_x = mod(proc_id, nPx)

c     Storing temperature values from sub_temp only on mask points (value 1) in sub_mask
      num_ones = 0
      sub_count = 2
      big_row = 0
      big_col = 0

      ! Looping column-major order

      do col=1, sNx, 1
        do row=1, sNy, 1
          if (sub_mask(row,col) .eq. 1) then
            big_row = int(row+p_y*sNy)
            big_col = int(col+p_x*sNx)

            sub_tempOnMask(sub_count) = sub_temp(row,col)
c           Saving corresponding global index
            sub_global_indices(sub_count) = int((Nx*(big_row-1)) !Index to be corrected
     &                                           + big_col)
            num_ones = num_ones + 1
            sub_count = sub_count + 1
          end if
        end do
      end do

c     Printing number of ones sub_mask has for each process
      print *,"Process", proc_id, "has", num_ones,
     &                "number of ones"

c     Saving number of valid elements in each sub array
      sub_tempOnMask(1) = num_ones
      sub_global_indices(1) = num_ones
c      print *, "Global indices: "
c      call PRINT_INT_ARR(sub_global_indices, sNy+sNx+1)

c     If there are mask points (value 1) in the sub__mask, send sub_tempOnMask with message tag 1. Otherwise, send with message tag 0.
      if (num_ones .gt. 0) then
        !MPI_SEND params: data_to_send, send_count, send_type, destination_ID, tag, comm, ierror
        print *,"A tag one subarray from", proc_id
        call PRINT_FLOAT_ARR(sub_tempOnMask, sNx+sNy)
        call MPI_SEND(sub_tempOnMask, sNx+sNy+1, MPI_REAL, 0, 1,
     &                MPI_COMM_WORLD, ierror)
        print *,"process", proc_id, "sends tag 1 to Parent"
      else
        call MPI_SEND(sub_tempOnMask, sNx+sNy+1, MPI_REAL, 0, 0,
     &                MPI_COMM_WORLD, ierror)
        print *,"process", proc_id, "sends tag 0 to Parent"
      end if

c     Send sub global indices array to parent process.
      call MPI_SEND(sub_global_indices, sNx+sNy, MPI_INT, 0, 2,
     &              MPI_COMM_WORLD, ierror)






C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

c     Final tasks for parent process:
c     - Receive sub arrays from all processes. Only include sub arrays with message tag 1 into global final mask array.
c     - Write global mask, temp, and final mask array to binary file.
      if (proc_id .eq. 0) then
        open(100, FILE=output_file, status="old", position="append",
     &       action="write")


c       Initialize arrays
        do col=1, ((sNy+sNx)*(nPx*nPy)), 1
          temp_on_mask(col) = 0.0
          global_indices(col) = 0
        end do

c       Initializing temp_on_mask and global_indices array counter
        final_col = 1
        write(100,FMT=*) "Parent Message Recieve Status: "
        write(100,FMT=*) " "

c       Receiving messages from all processes
        do p_count=0, nPx*nPy-1, 1
c         MPI_RECV params: received_data, receive_count, receive_type, sender_ID, tag, comm, status, ierr
c         Receive message tag 2 first.
          call MPI_RECV(sub_global_indices, sNx*sNy+1, MPI_REAL, p_count
     &                  ,2 ,MPI_COMM_WORLD, status, ierror)
          write(100,FMT=*)"Parent received from process",
     &              status(MPI_SOURCE), "with tag", 2
c         Then receive all other messages sent to parent process.
          call MPI_RECV(sub_tempOnMask, sNx*sNy+1, MPI_REAL, p_count,
     &                  MPI_ANY_TAG ,MPI_COMM_WORLD, status, ierror)
          write(100,FMT=*)"Parent received from process",
     &              status(MPI_SOURCE), "with tag", status(MPI_TAG)


c         If sub array comes with a message tag 1, include in global final mask.
          if (status(MPI_TAG) .eq. 1) then
            write(100,FMT=*)"Parent receives Tag 1 from Process",
     &                status(MPI_SOURCE)!, "with sub array:"
            !call PRINT_FLOAT_ARR(sub_tempOnMask, sNx+sNy)
            write(100,FMT=*)"Time to put Process",status(MPI_SOURCE),
     &               "'s sub array in global final mask."
            write(100,FMT=*)"There are ", sub_global_indices(1),
     &               "valid elements to save into global final mask"


c           Values used to compute the corresponding indices in final_mask to indices of sub array
            p_y = (status(MPI_SOURCE) / nPx)     ! Integer division does not round.
            p_x = (mod(status(MPI_SOURCE), nPx))  ! The remainder is equivalent to tile(process) number along x_axis.


            sub_num_ones = sub_global_indices(1)
            num_ones = num_ones + sub_num_ones
            do col=2, sub_num_ones+1, 1
              temp_on_mask(final_col) = sub_tempOnMask(col)
              global_indices(final_col) = sub_global_indices(col)
              final_col = final_col + 1
            end do

          else
            write(100,FMT=*)"Parent receives Tag 0 from Process",
     &        proc_id, ". Nothing to include in global temp on mask"
          end if

        end do

c       Verifying whether final mask is correctly generated.
        call TEST_FINAL(num_ones, temp_on_mask, global_indices)

c       Printing global final mask array and global mask array:
        print *, "Total number of valid elements in final:", num_ones
        print *, "Final updated temperature array over mask: "
        call PRINT_FLOAT_ARR(temp_on_mask, (sNy+sNx)*(nPx*nPy))
        print *, "Final global index array over mask: "
        call PRINT_INT_ARR(global_indices, (sNy+sNx)*(nPx*nPy))

c       Writing final T_on_mask global array to binary files
        call WRITE_GLOBAL_BIN(temp_on_mask)


      end if
      close(100)


c     End all MPI communication. Clean MPI environment.
      call MPI_FINALIZE(ierror)

      stop
      end

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
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


      SUBROUTINE CREATE_FILENAME(proc_id, prefix, is_child, filename)
c     Creating filenames for all binary files.

        include 'SIZE.h'
        integer, intent(in) :: proc_id
        integer, intent(in) :: is_child
        character*50, intent(out) :: filename

        integer p_x, p_y
        character*3 str_px, str_py
        character*30 prefix

c       Indices for each process with respect to global mask divided by the number of processes
        p_y = (proc_id / nPx) + 1     ! Integer division does not round.
        p_x = (mod(proc_id, nPx)) + 1  ! The remainder is equivalent to tile(process) number along x_axis.

c       If creating a filename for a child process, create a name in the format: T.001.002.bin
c       T for temp or mask for mask
c       0001 for index 1 in x axis
c       0002 for index 2 in x axis
        if (is_child .eq. 1) then
          if (p_x .lt. 10) then
            write(str_px, "(I3.3)") p_x
          else if (p_x .ge. 10) then
            write(str_px, "(I3.2)") p_x
          else if ((p_x .ge. 100) .and. p_x .le. 999) then
            write(str_px, "(I3.1)") p_x
          end if

          if (p_y .lt. 10) then
            write(str_py, "(I3.3)") p_y
          else if (p_y .ge. 10) then
            write(str_py, "(I3.2)") p_y
          else if ((p_y .ge. 100) .and. p_y .le. 999) then
            write(str_py, "(I3.1)") p_y
          end if
          write(filename, "(A, A3, A1, A3, A4)") prefix, str_px(1:3),
     &                                          ".", str_py(1:3), ".bin"

        else
c         Creating a filename for global arrays
          write(filename, "(A, A11)") prefix, "_global.bin"

        end if
      end


      SUBROUTINE WRITE_SUB_BIN()
c       Writing all sub arrays that each process receives from the Parent to binary files
        include 'SIZE.h'

        integer proc_id
        real, dimension(sNy,sNx) :: sub_temp
        integer, dimension(sNy,sNx) :: sub_mask
        real, dimension(sNy,sNx) :: read_sub_temp
        integer, dimension(sNy,sNx) :: read_sub_mask
        integer, parameter :: print_write = 1
        integer, parameter :: print_read = 0

        common /sub_arr_pid/ sub_mask, sub_temp, proc_id

        character*30 prefix
        character*50 filename

c       Writing sub_temp to binary file:
c       CREATE_FILENAME params: proc_id, prefix, is_child, filename
        prefix = "./output_masks/T."
        call CREATE_FILENAME(proc_id, prefix, 1, filename)
        open(1, FILE=filename, FORM="unformatted", ACCESS="stream")
        write(1) sub_temp
        close(1)

        if (print_write .eq. 1) then
          print *, "Process", proc_id, "'s sub_temp written to ",
     &             filename
        end if

        if (print_read .eq. 1) then
c         Reading binary files to make sure sub temp is correctly written
          open(1, FILE=filename, FORM="unformatted", ACCESS="stream")
          print *, "Reading sub_temp back from ", filename
          read(1) read_sub_temp
          call PRINT_FLOAT_ARR(read_sub_temp, sNx*sNy)
          close(1)
        end if

c       Writing sub_mask to binary file
        prefix = "./output_masks/mask."
        call CREATE_FILENAME(proc_id, prefix, 1, filename)
        open(2, FILE=filename, FORM="unformatted", ACCESS="stream")
        write(2) sub_mask
        close(2)

        if (print_write .eq. 1) then
          print *, "Process", proc_id, "'s sub_mask written to ",
     &             filename
        end if

c         Reading binary files to make sure sub mask is correctly written
        if (print_read .eq. 1) then
          open(2, FILE=filename, FORM="unformatted", ACCESS="stream")
          print *, "Reading sub_mask back from ", filename
          read(2) read_sub_mask
          call PRINT_INT_ARR(read_sub_mask, sNx*sNy)
          close(2)
        end if


      end


      SUBROUTINE WRITE_GLOBAL_BIN(temp_on_mask)
c       Writing all global arrays to binary files

        include "SIZE.h"
        real, dimension((sNy+sNx)*(nPx*nPy)) :: temp_on_mask
        real, dimension((sNy+sNx)*(nPx*nPy)) :: read_temp_on_mask
        character*50 filename
        character*30 prefix
        integer, parameter :: print_output = 0

c       Writing global final mask to binary file
        prefix = "./output_masks/T_on_mask"
        call CREATE_FILENAME(0,prefix, 0, filename)
        open(3, FILE=filename, FORM="unformatted", ACCESS="stream")
        write(3) temp_on_mask
        close(3)

        if (print_output .eq. 1) then
          print *, "Wrote final_mask to ", filename
c         Reading file to make sure global final mask is correctly written to binary file
          open(3, FILE=filename, FORM="unformatted", ACCESS="stream")
          print *, "Reading final_mask back from ", filename
          read(3) read_temp_on_mask
          call PRINT_FLOAT_ARR(read_temp_on_mask, (sNy+sNx)*(nPx*nPy))
          close(3)
        end if

      end

      SUBROUTINE TEST_TEMP_READ(temp)
        include "SIZE.h"
        real, dimension(Ny, Nx) :: temp
        integer row, col
        integer correct
        integer, parameter :: print_correctness=1

c      Array to check if mask and temp are distributed to processes correctly:
c     -1 for not being sent, 1 for being correct, 0 for incorrect.
        integer, dimension(Ny, Nx) :: temp_correct
        do row=1, Ny, 1
          do col=1, Nx, 1
            temp_correct(row,col) = 0
          end do
        end do

        do row=1, Ny, 1
          do col=1, Nx, 1
            if (temp(row,col) .eq. row*col) then
              correct = 1
              temp_correct(row,col) = 1
            else
              correct = 0
              temp_correct(row,col) = 0
              if (print_correctness .eq. 1) then
                print *, "TEST_TEMP_READ: global temp read ",
     &                   "incorrectly."
              end if
            end if
          end do
        end do
        if (print_correctness .eq. 1) then
          if (correct .eq. 1) then
            print *, "TEST_TEMP_READ: global temp is read ",
     &               "correctly."
          end if
        end if
      end


      SUBROUTINE TEST_SUBMASKS(mask)
c    Checking if all sub_masks are correctly distributed to processes
c    by comparing distributed sub arrays to manually created sub arrays from the global arrays

c     NOTE: Tile number = process id. Therefore, tile numbers and process id's start from 0.
c           In this subroutine: when separately counting tiles(procesess) along the x axis and y axis, we start from 0
c           in order to correctly compute the row and col of mask and temp corresponding to the row and col of sub_mask and sub_temp
c     NOTE: However, the Fortran convention of starting all array indices from 1 is held.

        include "SIZE.h"

        integer, dimension(sNy,sNx) :: sub_mask
        real, dimension(sNy,sNx) :: sub_temp
        integer, dimension(Ny,Nx) :: mask

c       Arrays to check if mask and temp are distributed to processes correctly:
c       1 for being correct, 0 for incorrect.
        integer, dimension(sNy,sNx) :: sub_mask_correct
        integer, dimension(sNy,sNx) :: sub_temp_correct

        integer proc_id, big_col, big_row, row, col
        common /sub_arr_pid/ sub_mask, sub_temp, proc_id
c       p_x is the process number in x axis, p_y is the tile(process) number in y axis
        integer p_x, p_y
        integer correct
        integer, parameter :: print_correctness = 0

c       Initialize mask_correct and temp_correct
        do row=1, sNy, 1
          do col=1, sNx, 1
              sub_mask_correct(row,col) = -1
            sub_temp_correct(row,col) = -1
          end do
        end do

c       Determining p_x and p_y using proc_id using integer division and modulus:
        p_y = proc_id / nPx      ! Integer division does not round.
        p_x = mod(proc_id, nPx)  ! The remainder is equivalent to tile(process) number along x_axis.

        big_col = 0
        big_row = 0

c       Looping through sub arrays and comparing each element to corresponding elements of global arrays
        do row=1, sNy, 1
          do col=1, sNx, 1
            big_col = int(col+p_x*sNx)
            big_row = int(row+p_y*sNy)
c           Comparing sub mask to corresponding piece of global mask
            if (sub_mask(row, col) .eq. mask(big_row, big_col)) then
              correct = 1
              sub_mask_correct(row, col) = 1
            else
              correct = 0
              sub_mask_correct(row, col) = 0
            end if
c           Comparing sub temp to corresponding piece of global temp which  is defined by the function global row * global col
            if (sub_temp(row, col) .eq. big_row*big_col) then
              correct = 1
              sub_temp_correct(row, col) = 1
            else
              correct = 0
              sub_temp_correct(row, col) = 0
            end if

          end do
        end do

c       correct == 1 indicates sub arrays were correctly generated.
c       correct == 0 indicates sub arrays are incorrect.
        if (correct .eq. 1) then
          if (print_correctness .eq. 1) then
            print *, "TEST_SUBMASKS: Process", proc_id, "received ",
     &               "correct sub_mask and sub_temp"
          end if
        else
          print *, "TEST_SUBMASKS: INCORRECTLY DISTRIBUTED FOR PROCESS"
     &              , proc_id
          print *, "TEST_SUBMASKS: sub_mask: 1 for correct. 0 for ",
     &             "incorrect: "
          call PRINT_INT_ARR(sub_mask_correct, Ny*Nx)
          print *, "TEST_SUBMASKS: sub_temp: 1 for correct. 0 for ",
     &             "incorrect: "
          call PRINT_INT_ARR(sub_temp_correct, Ny*Nx)
        end if
      end


      SUBROUTINE TEST_FINAL(num_elem, temp_on_mask, global_indices)
c       Checking if global temp on mask is correctly generated

        include "SIZE.h"
        integer, dimension((sNy+sNx)*(nPx*nPy)) :: global_indices
        real, dimension((sNy+sNx)*(nPx*nPy)) :: temp_on_mask

c       Array to check if mask and temp are distributed to processes correctly:
c       -1 for not being sent, 1 for being correct, 0 for incorrect.
        integer, dimension((sNy+sNx)*(nPx*nPy)) :: tempOnMask_correct
        integer i, correct, num_elem, global_index
        integer row, col
        integer, parameter :: print_correctness = 0
        correct = 0

c       initialize tempOnMask_correct
        do i=1, (sNy+sNx)*(nPx*nPy), 1
          tempOnMask_correct = -1
        end do

        do i=1, num_elem, 1
          global_index = global_indices(i)
          row = int(global_index/Nx) + 1
          col = mod(global_index, Nx)
          !print *, "TEST_FINAL row col: ", row, col
          if (temp_on_mask(i) .eq. row*col) then
            correct = 1
            tempOnMask_correct(i) = 1
            if (print_correctness .eq. 1) then
              print *, "TEST_FINAL Co: ", temp_on_mask(i), row*col
            end if
          else
            correct = 0
            tempOnMask_correct(i) = 0
            print *, "TEST_FINAL: INCORRECT", temp_on_mask(i),
     &                row*col
          end if
        end do

        if (correct .ne. 1) then
          print *, "TEST_FINAL: tempOnMask_mask: 1 for correct, ",
     & "0 for incorrect:"
          call PRINT_INT_ARR(tempOnMask_correct,(sNy+sNx)*(nPx*nPy))
        else
          print *, "TEST_FINAL: Global temp on mask is correct."
        end if
      end
