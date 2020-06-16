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
c 7. where they do have the mask they will print out the 'temperature' at the mask points
c 8. define temperature as T(row,col) = row*col

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c234567890




      program create_output

      include "SIZE.h"
      include "mpif.h"
      !include "/home/mitgcm/Work/MITgcm/pkg/mdsio/MDSIO_OPTIONS.h"
      integer, dimension(Ny,Nx) :: mask
      real, dimension(Ny,Nx) :: temp

c     Arrays to check if mask and temp are distributed to processes correctly:
c     -1 for not being sent, 0,1,2,3 for corresponding process to which piece of array was sent
      integer, dimension(Ny, Nx) :: mask_distr
      integer, dimension(Ny, Nx) :: temp_distr

c     Arrays to store pieces of mask and temp to send to each process
c
      integer, dimension(sNy,sNx) :: sub_mask
      real, dimension(sNy, sNx) :: sub_temp

c     The arrays above are used only in Process 0,
c     and, therefore, are only initialized in Process 0. (Under first if statement)


      integer row, col, p_count
      integer ierror, proc_id

c     Initialize MPI environment
      call MPI_INIT(ierror)

c     return process id (process_rank) of the process that called the function
      call MPI_COMM_RANK(MPI_COMM_WORLD, proc_id, ierror)


c     Parent Process' task:

c     Distribute pieces of mask and temp to each process
c     Looping through tiles in mask and temp
      if (proc_id .eq. 0) then

c     Initialize all arrays
        do row=1, Ny, 1
          do col=1, Nx, 1
            mask(row,col) = 0
            temp(row,col) = 0.0
            mask_distr(row,col) = -1
            temp_distr(row,col) = -1
            if ((row .le. sNy) .and. (col .le. sNx)) then
              sub_mask(row,col) = 0
              sub_temp(row,col) = 0.0
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


        do p_count=0, nPx-1, 1

          do row=1, sNy, 1
            do col=1, sNx, 1
              sub_mask(row, col) = mask(row, int(col+p_count*sNx))
              sub_temp(row, col) = temp(row, int(col+p_count*sNx))
            end do
          end do

        !MPI_SEND params: data_to_send, send_count, send_type, destination_ID, tag, comm, ierror
        call MPI_SEND(sub_mask, sNx*sNy, MPI_INT, p_count, 0,
     & MPI_COMM_WORLD, ierror)
c       Printing sub_masks sent through MPI_SENT
        print *, "Sent to process", p_count, "with sub_mask ",
     & "dimensions: (", shape(sub_mask), ")"
        print *, "For column indices of full mask: ", 1+p_count*sNx,
     & "to ", 10+p_count*sNx, "and row indices: 1 to 10"
        print *, "To get the sub_mask: "
        call PRINT_INT_ARR(sub_mask, sNx*sNy)

        end do

c        call PRINT_INT_ARR(mask(1:10, 21:30), 100)

      end if


c     ALL Processes' task:
c       MPI_RECV params: received_data, receive_count, receive_type, sender_ID,
c                        tag, comm, status, ierr
c        call MPI_RECV(message_item, 1, MPI_INT, num_procs-1, num_procs-1
c     & ,MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierror)




      call MPI_FINALIZE(ierror)
      end


C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

      SUBROUTINE PRINT_INT_ARR(array, arr_length)
      integer arr_length
      integer, dimension(arr_length) :: array
      do, i=1, arr_length
        write(*,"(I1,$)") array(i)
          if (i .EQ. arr_length) then
             write(*, '(A,/)') ''
          endif
      enddo

      end



      SUBROUTINE PRINT_FLOAT_ARR(array, arr_length)
      integer arr_length
      real, dimension(arr_length) :: array
      do, i=1, arr_length
        write(*,"(100g15.5,$)") array(i)
          if (i .EQ. arr_length) then
             write(*, '(A,/)') ''
          endif
      enddo

      end
