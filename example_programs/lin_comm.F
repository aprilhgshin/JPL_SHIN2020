c A program where each process communicates to the next process
c in a linear manner until the last process is reached and sends message
c to the parent process.

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c234567890


      program linear_communication
      include './mpif.h'

      integer sum, next_proc_id, message_item, parent_id
      integer ierror, num_procs, proc_id

c     parent process id set to 0
      parent_id = 0

c     initialize MPI environment
      call MPI_INIT(ierror)

c     return number of processes
      call MPI_COMM_SIZE(MPI_COMM_WORLD, num_procs, ierror)

c     return process id (process_rank) of the process that called the function
      call MPI_COMM_RANK(MPI_COMM_WORLD, proc_id, ierror)

c     edge case: if there is only one process, there is nothing to be done.
      if (num_procs .eq. 1) then
        print *, "There is only process in this computer. Therefore,",
     & " total sum of process ID's is 0."
        call MPI_FINALIZE(ierror)
        stop
      end if

c     parent process tasks:
      if(proc_id .eq. parent_id) then
        message_item = parent_id
c       MPI_SEND params: data_to_send, send_count, send_type, destination_ID,
c                        tag, comm, ierror
        call MPI_SEND(message_item, 1, MPI_INT, 1, 0, MPI_COMM_WORLD,
     & ierror)
        print *, "Process", parent_id ,"Sending the message: ",
     & message_item,"to process", parent_id + 1
c       MPI_RECV params: received_data, receive_count, receive_type, sender_ID,
c                        tag, comm, status, ierr
        call MPI_RECV(message_item, 1, MPI_INT, num_procs-1, num_procs-1
     & ,MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierror)
        print *, "Process", parent_id, "Received message from ",
     & "process",num_procs - 1, "with message containing total sum:",
     & message_item

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

      else
c       We expect the sender processor to have an id one less than current proc
        call MPI_RECV(message_item, 1, MPI_INT, proc_id - 1, proc_id - 1
     & ,MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierror)
        print *, "Process",proc_id,"Received the message: ",
     & message_item, "from process", proc_id - 1, "with tag",
     & proc_id - 1
        sum = proc_id + message_item
        next_proc_id = proc_id + 1

        if (next_proc_id < num_procs) then
          call MPI_SEND(sum, 1, MPI_INT, next_proc_id, proc_id,
     &  MPI_COMM_WORLD, ierror)
          print *, "Process", proc_id, "Sending the message: ",
     &  sum, "to process", next_proc_id

        else
          call MPI_SEND(sum, 1, MPI_INT, parent_id, proc_id,
     &  MPI_COMM_WORLD, ierror)
          print *, "Process",proc_id,"Sending the message: ",
     & sum, "to process", parent_id

        end if

      end if


c     clean up the MPI environment and end MPI communications
      call MPI_FINALIZE(ierror)

      end
