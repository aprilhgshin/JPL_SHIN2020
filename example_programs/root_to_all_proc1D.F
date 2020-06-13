c Created a program where parent (root) communicates to child
c processes and child processes communicate back to root using MPI gather and scatter.
c Two child processes communicate with each other.
c Array passed from process 1 to process 2. Properly work with pointers

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c234567890

      program root_to_all_proc_communication
      include "./mpif.h"

      integer ierror, num_procs, proc_id, count
      integer, parameter :: seed = 1000

c     Specifying maximum number of processors because arrays cannot dynamically change
c     Change number of processors accordingly
      integer, parameter :: max_num_proc = 4

c     Number of elements to pass to each processor
      integer, parameter :: num_elem = 2

c     One dimensional arrays
      integer, dimension(max_num_proc*num_elem) :: disperse_arr
      integer, dimension(max_num_proc*num_elem) :: gathered_arr
      integer, dimension(num_elem) :: scattered_data
      integer, dimension(num_elem) :: message_arr

c     Choose whether you would like modify the scattered array of each process:
      logical :: modify = .false.

c     Initializing array with random values to scatter (disperse to child processes)
      call srand(seed)
      do count=1, max_num_proc*num_elem, 1
        disperse_arr(count) = irand()/10000000 !irand: integer random
      end do

c     initialize MPI environment
      call MPI_INIT(ierror)

c     return number of processes
      call MPI_COMM_SIZE(MPI_COMM_WORLD, num_procs, ierror)

c     return process id (process_rank) of the process that called the function
      call MPI_COMM_RANK(MPI_COMM_WORLD, proc_id, ierror)

c      MPI_Scatter(
c              disperse_arr,   !Array we are scattering from
c              num_elem,       !Number of items we are sending each processor
c              MPI_INT,        !MPI Datatype of elements in array to scatter.
c              scattered_data, !Pre-declared variable at which piece of array to scatter will be saved.
c              num_elem,       !Amount of data each process will receive.
c              MPI_INT,        !MPI Datatype of receiver array.
c              0,              !Process ID that will distribute the data.
c              MPI_COMM_WORLD  !MPI Communicator.
c              ierror          !Error Handling Variable
c      )
      if (proc_id .eq. 0) then
        print *, "Array before scattering: ", disperse_arr
        print *, "Now scattered array: "
      end if

c     Dispersing a 1 by num_procs*2 array so that each process receives a 1 by 2 array
      call MPI_Scatter(disperse_arr, num_elem, MPI_INT, scattered_data,
     & num_elem, MPI_INT, 0, MPI_COMM_WORLD, ierror)
      print *, "Process ", proc_id, "received", scattered_data


C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c234567890

c     Modifying scattered array of each child process
      if (modify) then
        do count = 1, 2, 1
          scattered_data(count) = scattered_data(count) + 1000
        end do

        print *, "Modified scattered_data of process", proc_id, "to",
     & scattered_data

c       Process 1 sends its modified scattered_data to Process 2.
c       However, no data should be modified during this process.
c       Goal: Work with pointers and don't accidentally modify pointer values.
        if (proc_id .eq. 1) then
          message_arr = scattered_data
          call MPI_SEND(message_arr, num_elem, MPI_INT, num_elem, 1,
     & MPI_COMM_WORLD ,ierror)
          print *, "Process", proc_id ,"Sending the message: ",
     & message_arr,"to process", proc_id + 1
        else if (proc_id .eq. 2) then
          call MPI_RECV(message_arr, num_elem, MPI_INT, proc_id-1,
     & proc_id-1 ,MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierror)
          print *, "Process", proc_id, "Received message from ",
     & "process",proc_id - 1, "with message:",
     & message_arr
          print *, "Process", proc_id, "has its own scattered_data",
     & scattered_data
        end if

      end if


c In the gather operation, each process (root process included) sends scount(2) elements
c of type MPI_INT of scattered_data to the root process. The root process receives the messages
c and stores them in rank order in the gathered_data. For scatter, the reverse holds.
c The root process sends a buffer of N chunks of data (N = number of processes
c in the group) so that process 1 gets the first element(chunk), process 2 gets the
c second element(chunk), etc.
c Explanation source: https://cvw.cac.cornell.edu/MPIcc/gatherscatter
c and https://mpitutorial.com/tutorials/mpi-scatter-gather-and-allgather/

c     MPI_Gather(
c       scattered_data      !Array with values we will send from child process.
c       num_elem            !Number of data elements that will sent .
c       MPI_INT             !MPI Datatype of the data that is sent.
c       gathered_data       !Variable that will store the gathered data.
c       num_elem            !Number of data elements per process that will be received.
c       MPI_INT             !MPI Datatype of the data that will be received.
c       0                   !The rank of the process rank that will gather the information.
c       MPI_COMM_WORLD      !The MPI_Communicator.
c       ierror              !An error handling variable.
c     )

      call MPI_Gather(scattered_data, num_elem, MPI_INT, gathered_arr,
     & num_elem, MPI_INT, 0, MPI_COMM_WORLD, ierror)
      if (proc_id .eq. 0) then
        print *, "In Parent Process: Combined arrays from child",
     & " processes: ", gathered_arr
      end if

      call MPI_FINALIZE(ierror)
      end
