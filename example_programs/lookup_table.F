C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c234567890

      program lookup_table

      include "SIZE.h"
      include "mpif.h"

      integer, parameter :: num_openB = 2
c     Each open boundary will have a 3D array with the following info:
c     1) open boundary mask,  2) temp
      integer, parameter :: num_info_masks = 2

      real temp_on_mask(num_openB,(sNy+sNx)*(nPx*nPy))
      integer global_indices(num_openB,(sNy+sNx)*(nPx*nPy))

      integer lookup(num_openB, Ny*Nx)
c     global_ob_masks array is a 3D array. First dimension corresponds to different open boundaries i.e. mask
      integer global_ob_masks(num_openB, Ny, Nx)
      real global_temps(num_openB, Ny, Nx)

c     In Parent Process 0: Used as arrays to store pieces of mask and temp to send to each process
c     All Processes 0-3: Used as buffer arrays to store pieces of mask and temp each process receives
c     NOTE: sub_masks1 and sub_masks2 are 3D arrays. First dimension of both corresponds to different info masks:
c           1) open boundary mask
c           2) temp
      real sub_masks1(num_info_masks, sNy, sNx) !open boundary 1
      real sub_masks2(num_info_masks, sNy, sNx) !open boundary 2
      !integer global_ob_index(num_openB)
      integer global_ob_index

c     Store temperature values only on elements == 1 on each sub_mask
      real sub_tempOnMask(num_openB, sNy+sNx+1)
      integer sub_global_indices(num_openB, sNy+sNx+1)

      integer proc_id, ierror, status(MPI_STATUS_SIZE)
      integer row, col, p_x, p_y, p_count, ob, lookup_col, infomask_id
      integer global_row, global_col, counter1, counter2
      integer num_valid1, num_valid2, counter, ob_index

      integer, parameter :: PRINT_PARENT_SEND = 0

      global_row = 1
      global_col = 1
      global_ob_index = 1


c     Initialize arrays used by all processes
      do ob=1, num_openB, 1
        do row=1, Ny, 1
          do col=1, Nx, 1
            global_ob_masks(ob, row, col) = 0
          end do
        end do
      end do

      do infomask_id=1, num_info_masks, 1
        do row=1, sNy, 1
          do col=1, sNx, 1      !if ((row .le. sNy) .and. (col .le. sNx)) then
            sub_masks1(infomask_id,row,col) = 0.0
            sub_masks2(infomask_id,row,col) = 0.0
          end do
        end do
      end do

      do ob=1, num_openB, 1
        do col=1, sNy+sNx+1, 1
          sub_tempOnMask(ob, col) = 0.0
          sub_global_indices(ob, col) = 0
          !if (col .le. 2) then
          !  global_ob_index(col) = 0
          !end if
        end do
      end do




c     Read global_ob_masks and temp data from binary files
       open(1, FILE="input_domains/domain_mask1.bin", FORM="unformatted"
     & ,ACCESS="stream")
       read(1) global_ob_masks(1, 1:Ny, 1:Nx)
       close(1)
       !call PRINT_INT_ARR(global_ob_masks(1,1:,1:), Ny*Nx)


       open(2, FILE="input_domains/domain_mask2.bin", FORM="unformatted"
     & ,ACCESS="stream")
       read(2) global_ob_masks(2, 1:Ny, 1:Nx)
       close(2)
       !call PRINT_INT_ARR(global_ob_masks(2, 1:Ny, 1:Nx), Ny*Nx)



c     Initialize MPI environment
      call MPI_INIT(ierror)

c     Return (current) process id of the process that called the function
      call MPI_COMM_RANK(MPI_COMM_WORLD, proc_id, ierror)


C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

c     Parent Process' task:
      if (proc_id .eq. 0) then


       !Initialize lookup table
       do row=1, num_openB, 1
         do col=1, Ny*Nx, 1
           lookup(row,col)=0
         end do
       end do

       !Initialize global_temps array
       do infomask_id=1, num_info_masks, 1
         if (ob .eq. 1) then
           global_temps(infomask_id, row,col) = 0.0
         end if
       end do

       !Store global_ob_masks into lookup table
       do ob=1, num_openB, 1
         lookup_col = 1
         do row=1, Ny, 1
           do col=1, Nx, 1
             lookup(ob, lookup_col)
     &            = global_ob_masks(ob, row, col)
             lookup_col = lookup_col + 1
           end do
         end do
         !print *, "Global Mask 1: "
         !call PRINT_INT_ARR(lookup(ob, 1:Ny*Nx), Ny*Nx)
       end do


       !Read temps binary file - for this program, all temp masks will be exactly the same.
       open(3, FILE="input_domains/domain_temp.bin", FORM="unformatted"
     &       ,ACCESS="stream")
       read(3) global_temps(1, 1:Ny, 1:Nx)
       close(3)

       open(4, FILE="input_domains/domain_temp.bin", FORM="unformatted"
     &       ,ACCESS="stream")
       read(4) global_temps(2, 1:Ny, 1:Nx)
       close(4)

       !call PRINT_FLOAT_ARR(global_temps(1, 1:Ny, 1:Nx), Ny*Nx)
       !call PRINT_FLOAT_ARR(global_temps(2, 1:Ny, 1:Nx), Ny*Nx)
       !call TEST_TEMP_READ(global_temps(1, 1:Ny, 1:Nx), Ny*Nx)
       !call TEST_TEMP_READ(global_temps(2, 1:Ny, 1:Nx), Ny*Nx)


c      Distribute pieces of mask and temp of each open boundary to each process

        do p_y=0, 0,1!nPy-1, 1
          do p_x=0, 1,1!nPx-1, 1
            do row=1, sNy, 1
              do col=1, sNx, 1
c               Adding p_x*sNx to col and p_y*sNy to row allows us to move onto next portion
c               of mask and temp to save into sub_mask and sub_temp for each process.
c               Each process will split each open boundary on global domain.
c               In this program there are 2 open boundaries. Therefore, we define two sub_masks buffers.
                global_row = int(row+p_y*sNy)
                global_col = int(col+p_x*sNx)

                sub_masks1(1, row, col) = global_ob_masks(1,
     &                                    global_row, global_col)
                sub_masks1(2, row, col) = global_temps(1,
     &                                    global_row, global_col)

                sub_masks2(1, row, col) = global_ob_masks(2,
     &                                    global_row, global_col)
                sub_masks2(2, row, col) = global_temps(2,
     &                                    global_row, global_col)

              end do
            end do
            !global_ob_index defined at row=1, col=1 of each open boundary
            global_ob_index = int((Nx*(p_y*sNy)) + (1+p_x*sNx))
            p_count = p_x + (2*p_y)!(nPx*p_y)

            if (PRINT_PARENT_SEND .eq. 1) then
              print *, "Piece ", p_count
              print *, "global_ob_index: ", global_ob_index
              print *, "sub_mask1: ob_mask"
              call PRINT_FLOAT_ARR(sub_masks1(1, 1:sNy, 1:sNx), sNy*sNx)

              print *, "sub_mask1: temp"
              call PRINT_FLOAT_ARR(sub_masks1(2, 1:sNy, 1:sNx), sNy*sNx)

              print *, "sub_mask2: ob_mask"
              call PRINT_FLOAT_ARR(sub_masks2(1, 1:sNy, 1:sNx), sNy*sNx)

              print *, "sub_mask2: temp"
              call PRINT_FLOAT_ARR(sub_masks2(2, 1:sNy, 1:sNx), sNy*sNx)

              print *, " "
            end if


c           Parent sends pieces of masks of each open boundary to each process
            !MPI_SEND params: data_to_send, send_count, send_type, destination_ID, tag, comm, ierror
            call MPI_SEND(global_ob_index, 1, MPI_INT, p_count, 300,
     &                    MPI_COMM_WORLD, ierror)
            call MPI_SEND(sub_masks1, num_info_masks*sNx*sNy, MPI_REAL,
     &                    p_count, 100, MPI_COMM_WORLD, ierror)
            call MPI_SEND(sub_masks2, num_info_masks*sNx*sNy, MPI_REAL,
     &                    p_count, 200, MPI_COMM_WORLD, ierror)

          end do
        end do

      end if

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

c     ALL Processes' task:

c     Each process receives masks of each open boundary from Parent.
c     MPI_RECV params: received_data, receive_count, receive_type, sender_ID, tag, comm, status, ierr
      call MPI_RECV(global_ob_index, 1, MPI_INT,0,300,MPI_COMM_WORLD,
     &              status, ierror)
      call MPI_RECV(sub_masks1, num_info_masks*sNx*sNy, MPI_REAL,0,100,
     &              MPI_COMM_WORLD, status, ierror)
      call MPI_RECV(sub_masks2, num_info_masks*sNx*sNy, MPI_REAL,0,200,
     &              MPI_COMM_WORLD, status, ierror)

      print *, "Received sub_masks1 mask for piece", proc_id
      call PRINT_FLOAT_ARR(sub_masks1(1, 1:, 1:), sNx*sNy)
      print *, "Received global_ob_index mask1 for piece", proc_id
      print *, global_ob_index
      print *, "Received sub_masks1 temp for piece", proc_id
      call PRINT_FLOAT_ARR(sub_masks1(2, 1:, 1:), sNx*sNy)
      print *, "Received sub_masks2 mask for piece", proc_id
      call PRINT_FLOAT_ARR(sub_masks2(1, 1:, 1:), sNx*sNy)
      print *, "Received sub_masks2 temp for piece", proc_id
      call PRINT_FLOAT_ARR(sub_masks2(2, 1:, 1:), sNx*sNy)

      print *, " "
      counter1 = 2
      counter2 = 2
      num_valid1 = 0
      num_valid2 = 0


      do row=1, sNy, 1
        do col=1, sNx, 1
          if (sub_masks1(1,row,col) .ge. 1) then
            sub_tempOnMask(1,counter1) = sub_masks1(2,row,col)
            sub_global_indices(1,counter1) = global_ob_index
            counter1 = counter1 + 1
            num_valid1 = num_valid1 + 1
          end if
          if (sub_masks2(1,row,col) .ge. 1) then
            sub_tempOnMask(2,counter2) = sub_masks2(2,row,col)
            sub_global_indices(2,counter2) = global_ob_index
            counter2 = counter2 + 1
            num_valid2 = num_valid2 + 1
          end if

          global_ob_index = global_ob_index + 1

        end do
      end do

      sub_tempOnMask(1,1) = num_valid1
      sub_global_indices(1,1) = num_valid1
      sub_tempOnMask(2,1) = num_valid2
      sub_global_indices(2,1) = num_valid2

      print *, "sub_tempOnMask mask1 for piece", proc_id
      call PRINT_FLOAT_ARR(sub_tempOnMask(1, 1:), 1+sNx+sNy)
      print *, "sub_global_indices mask1 for piece", proc_id
      call PRINT_INT_ARR(sub_global_indices(1, 1:), 1+sNx+sNy)
      print *, "sub_tempOnMask mask2 for piece", proc_id
      call PRINT_FLOAT_ARR(sub_tempOnMask(2, 1:), 1+sNx+sNy)
      print *, "sub_global_indices mask2 for piece", proc_id
      call PRINT_INT_ARR(sub_global_indices(2, 1:), 1+sNx+sNy)


      call MPI_SEND(sub_tempOnMask, num_openB*(sNx+sNy+1), MPI_REAL,
     &              0, 1, MPI_COMM_WORLD, ierror)
      call MPI_SEND(sub_global_indices, num_openB*(sNx+sNy+1),
     &              MPI_INT, 0, 2, MPI_COMM_WORLD, ierror)



C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c     Final tasks for parent process:
c     - Receive sub arrays from all processes. Only include sub arrays with message tag 1 into global final mask array.
      if (proc_id .eq. 0) then

        !Initialize arrays
        do ob=1, num_openB, 1
          do col=1, ((sNy+sNx)*(nPx*nPy)), 1
            temp_on_mask(ob,col) = 0.0
            global_indices(ob,col) = 0
          end do
        end do

        ob_index = 1
c       Receiving messages from all processes
        do p_count=0, 1,1!nPx*nPy-1, 1
c         MPI_RECV params: received_data, receive_count, receive_type, sender_ID, tag, comm, status, ierr
          call MPI_RECV(sub_global_indices, num_openB*(sNx*sNy+1),
     &                  MPI_REAL, p_count, 2, MPI_COMM_WORLD, status,
     &                  ierror)
          call MPI_RECV(sub_tempOnMask, num_openB*(sNx*sNy+1), MPI_REAL,
     &                  p_count, 1 ,MPI_COMM_WORLD, status, ierror)

          print *, "Parent received sub_tempOnMask 1 from piece",
     &              p_count
          call PRINT_FLOAT_ARR(sub_tempOnMask(1, 1:), sNx*sNy+1)
          print *, "Parent received sub_tempOnMask 2 from piece",
     &              p_count
          call PRINT_FLOAT_ARR(sub_tempOnMask(2, 1:), sNx*sNy+1)


          do ob=1, num_openB, 1
            num_valid = sub_global_indices(ob, 1)
            !if there exists open boundary points in the mask, then append values to final array
            if (sub_global_indices(ob,1) .ge. 1) then
              do col=2,num_valid+1, 1
                ob_index = lookup(ob, sub_global_indices(ob, col))
                print *, "ob_index", ob_index
                temp_on_mask(ob,ob_index) = sub_tempOnMask(ob,col)
              end do
            end if
          end do
        end do

      print *, "temp_on_mask: 1"
      call PRINT_FLOAT_ARR(temp_on_mask(1,1:),(sNy+sNx)*(nPx*nPy))
      print *, "temp_on_mask: 2"
      call PRINT_FLOAT_ARR(temp_on_mask(2,1:),(sNy+sNx)*(nPx*nPy))
      end if



      call MPI_FINALIZE(ierror)

      end


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
