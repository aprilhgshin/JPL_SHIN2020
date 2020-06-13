C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c234567890
c This program sends pieces of a global array of 1's and 0's to 4 processes.
c For each process, if its piece contains a 1, it will send its piece back to the parent process with a tag 1.
c Only messages with tag 1 sent to the parent process will be appended to the final array.

      program select_tiles
      include "./mpif.h"

      integer ierror, num_procs, proc_id, count, final_count, i
      integer status(MPI_STATUS_SIZE), j
c     Note: status is an int array of size MPI_STATUS_SIZE
c           It has indices MPI_SOURCE, MPI_TAG, and more.
      integer, parameter :: GLOBAL_SIZE = 16
      integer, parameter :: SUB_SIZE = 4
      integer, dimension(GLOBAL_SIZE) :: global_arr
      integer, dimension(SUB_SIZE) :: sub_arr
      integer, dimension(GLOBAL_SIZE) :: final_arr
      integer, parameter :: seed = 1000


c     initialize MPI environment
      call MPI_INIT(ierror)
c     return number of processes
      call MPI_COMM_SIZE(MPI_COMM_WORLD, num_procs, ierror)
c     return process id (process_rank) of the process that called the function
      call MPI_COMM_RANK(MPI_COMM_WORLD, proc_id, ierror)

c     initialize global_arr with 0's and 1's
c     initialize final_arr with 0's
      call srand(seed)
      do count=1, GLOBAL_SIZE, 1
        global_arr(count) = mod(irand(),2)
        final_arr(count) = 0
      end do


c      MPI_Scatter(
c              disperse_arr,   !Array we are scattering from
c              num_elem,       !Number of items we are sending each processor
c              MPI_INT,        !MPI Datatype of elements in array to scatter.
c              scattered_data, !Pre-declared variable in which piece of array will be saved.
c              num_elem,       !Amount of data each process will receive.
c              MPI_INT,        !MPI Datatype of receiver array.
c              0,              !Process ID that will distribute the data.
c              MPI_COMM_WORLD  !MPI Communicator.
c              ierror          !Error Handling Variable
c      )
      call MPI_Scatter(global_arr, SUB_SIZE, MPI_INT, sub_arr,
     & SUB_SIZE, MPI_INT, 0, MPI_COMM_WORLD, ierror)

c     Checking if sub array contains any piece of open boundary i.e. contains a 1 anywhere
c     If it does, then send message with tag 1. Otherwise, send with tag 0.
      if (sum(sub_arr) .gt. 0) then
c       MPI_SEND params: data_to_send, send_count, send_type, destination_ID,
c                        tag, comm, ierror
        call MPI_SEND(sub_arr, SUB_SIZE, MPI_INT, 0, 1, MPI_COMM_WORLD,
     & ierror)

      else
        call MPI_SEND(sub_arr, SUB_SIZE, MPI_INT, 0, 0, MPI_COMM_WORLD,
     & ierror)
      end if

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

c     If the process ID is 0, then receive messages from all processes of any tag.
      if (proc_id .eq. 0) then
        final_count = 1
c       Receiving messages sent from processes 1, 2 and 3
        do j = 0, 3, 1
c        MPI_RECV params: received_data, receive_count, receive_type, sender_ID,
c                         tag, comm, status, ierr
          call MPI_RECV(sub_arr, SUB_SIZE, MPI_INT, MPI_ANY_SOURCE,
     & MPI_ANY_TAG ,MPI_COMM_WORLD, status, ierror)
          print *, "Received from process", status(MPI_SOURCE),
     & "with tag", status(MPI_TAG), "and sub array: "
          call PRINT_INT_ARR(sub_arr, SUB_SIZE)

c         if received tag is 1, then append sub array to final array to save to disk
          if (status(MPI_TAG) .eq. 1) then
            do i=1, SUB_SIZE, 1
              final_arr(final_count) = sub_arr(i)
              final_count = final_count + 1
            end do
            print *, "updated final_arr: "
            call PRINT_INT_ARR(final_arr, final_count)
          else
            print *, "Tag is 0. Nothing to append to final_arr"
          end if

        end do

        print *, "Final array to save to disk with empty ",
     & "elements: "
        call PRINT_INT_ARR(final_arr, GLOBAL_SIZE)

      end if


      call MPI_FINALIZE(ierror)

      end


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
