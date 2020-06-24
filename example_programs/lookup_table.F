C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c234567890


      module constants
        integer, parameter :: num_openB = 2
c     Each open boundary will have a 3D array with the following info:
c     1) open boundary mask,  2) temp
        integer, parameter :: num_info_masks = 2
      end module constants


      program lookup_table
      use constants

      include "SIZE.h"
      include "mpif.h"



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

c     global_ob_index defined at row=1, col=1 of each open boundary
      integer global_ob_index

c     Store temperature values only on elements == 1 on each sub_mask
      real sub_tempOnMask(num_openB, sNy+sNx+1)
      integer sub_global_indices(num_openB, sNy+sNx+1)

      integer proc_id, ierror, status(MPI_STATUS_SIZE)
      integer row, col, p_count, ob, lookup_col, infomask_id
      integer global_row, global_col
      integer total_valid1, total_valid2
      integer  subtile_x, subtile_y, subtile_id

      integer, parameter :: PRINT_PARENT_SEND = 0
      integer, parameter :: PRINT_SUB_RECV = 0
      integer ,parameter :: PRINT_SUB_SEND = 0
      integer, parameter :: PRINT_PARENT_RECV = 0
      integer, parameter :: debug = 1

      common /sub_arr_info/ sub_masks1, sub_masks2, global_ob_index
      common /sub_global/ sub_tempOnMask, sub_global_indices
      common /global_arr/ lookup, total_valid1, total_valid2,
     &                    temp_on_mask


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

        !Initialize arrays
        do ob=1, num_openB, 1
          do col=1, ((sNy+sNx)*(nPx*nPy)), 1
            temp_on_mask(ob,col) = 0.0
            global_indices(ob,col) = 0
          end do
        end do

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
       call TEST_TEMP_READ(global_temps(1, 1:, 1:))

       open(4, FILE="input_domains/domain_temp.bin", FORM="unformatted"
     &       ,ACCESS="stream")
       read(4) global_temps(2, 1:Ny, 1:Nx)
       close(4)
       call TEST_TEMP_READ(global_temps(2, 1:, 1:))

       !call PRINT_FLOAT_ARR(global_temps(1, 1:Ny, 1:Nx), Ny*Nx)
       !call PRINT_FLOAT_ARR(global_temps(2, 1:Ny, 1:Nx), Ny*Nx)



c      Distribute pieces(tiles) of mask and temp of each open boundary to each process
        p_count = 0
        do subtile_y=0, nPy-1, 1
          do subtile_x=0, nPx-1, 1
            do row=1, sNy, 1
              do col=1, sNx, 1
c               Adding subtile_x*sNx to col and subtile_y*sNy to row allows us to move onto next portion
c               of mask and temp to save into sub_mask and sub_temp for each process.
c               Each process will split each open boundary on global domain.
c               In this program there are 2 open boundaries. Therefore, we define two sub_masks buffers.
                global_row = int(row+subtile_y*sNy)
                global_col = int(col+subtile_x*sNx)

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
            global_ob_index = int((Nx*(subtile_y*sNy)) +
     &                        (1+subtile_x*sNx))
            subtile_id = subtile_x + (nPx*subtile_y)

            if (PRINT_PARENT_SEND .eq. 1) then
              print *, "subtile_id ", subtile_id
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
            if (p_count .gt. 0) then
c             Parent sends pieces of masks of each open boundary to each process other than itself - zero.
              !MPI_SEND params: data_to_send, send_count, send_type, destination_ID, tag, comm, ierror
              call MPI_SEND(global_ob_index, 1, MPI_INT, p_count, 300,
     &                      MPI_COMM_WORLD, ierror)
              call MPI_SEND(sub_masks1, num_info_masks*sNx*sNy, MPI_REAL
     &                      ,p_count, 100, MPI_COMM_WORLD, ierror)
              call MPI_SEND(sub_masks2, num_info_masks*sNx*sNy, MPI_REAL
     &                      ,p_count, 200, MPI_COMM_WORLD, ierror)
            else if (p_count .eq. 0) then
              print *, "Doing Process 0 tasks locally. No MPI."
              call CHILD_TASKS()
              call WRITE_SUB_BIN()
              call APPEND_OB_VALUES()

            end if

            !We are sending each subtile to each process by incrementing the destination process id at each loop.
            !No computation is dependent on the process id.
            p_count = p_count + 1

          end do
        end do

      end if

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

c     ALL Processes' task other than Process 0:
      if (proc_id .gt. 0) then

c       Each process receives masks of each open boundary from Parent.
c       MPI_RECV params: received_data, receive_count, receive_type, sender_ID, tag, comm, status, ierr
        call MPI_RECV(global_ob_index, 1, MPI_INT,0,300,MPI_COMM_WORLD,
     &                status, ierror)
        call MPI_RECV(sub_masks1, num_info_masks*sNx*sNy, MPI_REAL,0,100
     &                ,MPI_COMM_WORLD, status, ierror)
        call MPI_RECV(sub_masks2, num_info_masks*sNx*sNy, MPI_REAL,0,200
     &                ,MPI_COMM_WORLD, status, ierror)

        print *, "Testing sub_masks of all ob from process", proc_id,":"
        call TEST_SUBMASKS(global_ob_masks)


        if (PRINT_SUB_RECV .eq. 1) then
          print *, "Received sub_masks1 ob 1 with global_ob_index",
     &              global_ob_index
          call PRINT_FLOAT_ARR(sub_masks1(1, 1:, 1:), sNx*sNy)
          print *, "Received sub_masks1 temp for ob 1"
          call PRINT_FLOAT_ARR(sub_masks1(2, 1:, 1:), sNx*sNy)

          print *, "Received sub_masks2 ob 2 with global_ob_index",
     &              global_ob_index
          call PRINT_FLOAT_ARR(sub_masks2(1, 1:, 1:), sNx*sNy)
          print *, "Received sub_masks2 temp for ob 2"
          call PRINT_FLOAT_ARR(sub_masks2(2, 1:, 1:), sNx*sNy)
          print *, " "
        end if

        call CHILD_TASKS()
        call WRITE_SUB_BIN()


        if (PRINT_SUB_SEND .eq. 1) then
          print *, "Sub arrays with global_ob_index",global_ob_index,
     &             "will be sent to Parent:"
          print *, "sub_tempOnMask ob 1:"
          call PRINT_FLOAT_ARR(sub_tempOnMask(1, 1:), 1+sNx+sNy)
          print *, "sub_global_indices ob 1:"
          call PRINT_INT_ARR(sub_global_indices(1, 1:), 1+sNx+sNy)
          print *, "sub_tempOnMask ob 2:"
          call PRINT_FLOAT_ARR(sub_tempOnMask(2, 1:), 1+sNx+sNy)
          print *, "sub_global_indices ob 2:"
          call PRINT_INT_ARR(sub_global_indices(2, 1:), 1+sNx+sNy)
          print *, " "
        end if


        call MPI_SEND(sub_tempOnMask, num_openB*(sNx+sNy+1), MPI_REAL,
     &                0, 1, MPI_COMM_WORLD, ierror)
        call MPI_SEND(sub_global_indices, num_openB*(sNx+sNy+1),
     &                MPI_INT, 0, 2, MPI_COMM_WORLD, ierror)
      end if




C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c     Final tasks for parent process:
c     - Receive sub arrays from all processes. Only include sub arrays with message tag 1 into global final mask array.
      if (proc_id .eq. 0) then

        !total_valid1 = 0
        !total_valid2 = 0

c       Receiving messages from all processes other than itself
        do p_count=1,nPx*nPy-1, 1
c         MPI_RECV params: received_data, receive_count, receive_type, sender_ID, tag, comm, status, ierr
          call MPI_RECV(sub_global_indices, num_openB*(sNx+sNy+1),
     &                  MPI_REAL, p_count, 2, MPI_COMM_WORLD, status,
     &                  ierror)
          call MPI_RECV(sub_tempOnMask, num_openB*(sNx+sNy+1), MPI_REAL,
     &                  p_count, 1 ,MPI_COMM_WORLD, status, ierror)

          if (PRINT_PARENT_RECV .eq. 1) then
            print *, "Parent received sub_tempOnMask ob 1"
            call PRINT_FLOAT_ARR(sub_tempOnMask(1, 1:), sNx+sNy+1)
            print *, "Parent received sub_tempOnMask ob 2"
            call PRINT_FLOAT_ARR(sub_tempOnMask(2, 1:), sNx+sNy+1)
          end if

          call APPEND_OB_VALUES()

        end do

      print *, "FINAL ARRAY constructed in Parent process: "
      print *, "temp_on_mask: open boundary #1"
      print *, "Total valid elements: ", total_valid1
      call PRINT_FLOAT_ARR(temp_on_mask(1,1:),(sNy+sNx)*(nPx*nPy))
      print *, "temp_on_mask: open boundary #2"
      print *, "Total valid elements: ", total_valid2
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
      print *, " "
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
      print *, " "
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

      SUBROUTINE CHILD_TASKS()
        use constants
        include "SIZE.h"

        real sub_masks1(num_info_masks, sNy, sNx) !open boundary 1
        real sub_masks2(num_info_masks, sNy, sNx) !open boundary 2
        real sub_tempOnMask(num_openB, sNy+sNx+1)
        integer sub_global_indices(num_openB, sNy+sNx+1)

        integer global_ob_index, counter1, counter2
        integer num_valid1, num_valid2, row, col
        integer new_goi
        common /sub_arr_info/ sub_masks1, sub_masks2, global_ob_index
        common /sub_global/ sub_tempOnMask, sub_global_indices

        integer, parameter :: PRINT_DEBUG = 0

        counter1 = 2
        counter2 = 2
        num_valid1 = 0
        num_valid2 = 0
        new_goi = global_ob_index

        do row=1, sNy, 1
          do col=1, sNx, 1
            if (sub_masks1(1,row,col) .ge. 1) then
              sub_tempOnMask(1,counter1) = sub_masks1(2,row,col)
              sub_global_indices(1,counter1) = new_goi
              counter1 = counter1 + 1
              num_valid1 = num_valid1 + 1
            end if
            if (sub_masks2(1,row,col) .ge. 1) then
              sub_tempOnMask(2,counter2) = sub_masks2(2,row,col)
              sub_global_indices(2,counter2) = new_goi
              counter2 = counter2 + 1
              num_valid2 = num_valid2 + 1
            end if

            new_goi = new_goi + 1

          end do
        end do

        sub_tempOnMask(1,1) = num_valid1
        sub_global_indices(1,1) = num_valid1
        sub_tempOnMask(2,1) = num_valid2
        sub_global_indices(2,1) = num_valid2

        if (PRINT_DEBUG .eq. 1) then
          print *, "CHILD_TASKS: sub_tempOnMask ob1, global_ob_index:",
     &              global_ob_index
          call PRINT_FLOAT_ARR(sub_tempOnMask(1,1:),
     &                       (sNy+sNx+1))
          print *, "CHILD_TASKS: sub_tempOnMask ob2"
          call PRINT_FLOAT_ARR(sub_tempOnMask(2,1:),
     &                       (sNy+sNx+1))
        end if

      end


      SUBROUTINE APPEND_OB_VALUES()
        use constants
        include "SIZE.h"
        real sub_tempOnMask(num_openB, sNy+sNx+1)
        integer sub_global_indices(num_openB, sNy+sNx+1)

        real temp_on_mask(num_openB,(sNy+sNx)*(nPx*nPy))
        integer lookup(num_openB, Ny*Nx)

        integer ob, ob_index, total_valid1, total_valid2, col
        integer num_valid

        common /sub_global/ sub_tempOnMask, sub_global_indices
        common /global_arr/ lookup, total_valid1,total_valid2,
     &                        temp_on_mask


        ob_index = 1
        num_valid = 0
        do ob=1, num_openB, 1
          num_valid = sub_global_indices(ob, 1)
          if (ob .eq. 1) then
            total_valid1 = total_valid1 + num_valid
          else if (ob .eq. 2) then
            total_valid2 = total_valid2 + num_valid
          end if

          !if there exists open boundary points in the mask, then append values to final array
          if (num_valid .ge. 1) then
            print *, "APPEND_OB_VALUES: OB points exist in subtile ob",
     &                ob,"with starting global element",
     &               sub_global_indices(ob,2),". Append values."
            do col=2,num_valid+1, 1
              ob_index = lookup(ob, sub_global_indices(ob, col))
              temp_on_mask(ob,ob_index) = sub_tempOnMask(ob,col)
            end do
            print *, "APPEND_OB_VALUES: Updated temp_on_mask for ob:"
     &               ,ob
            call PRINT_FLOAT_ARR(temp_on_mask(1,1:),
     &                          (sNy+sNx)*(nPx*nPy))
          end if
        end do

      end


      SUBROUTINE TEST_SUBMASKS(global_ob_masks)
c    Checking if all sub_masks are correctly distributed to processes
c    by comparing distributed sub arrays to manually created sub arrays from the global arrays
        use constants
        include "SIZE.h"

        real sub_masks1(num_info_masks, sNy, sNx) !open boundary 1
        real sub_masks2(num_info_masks, sNy, sNx) !open boundary 2
        integer global_ob_masks(num_openB,Ny,Nx)

c       Arrays to check if mask and temp are distributed to processes correctly:
c       1 for being correct, 0 for incorrect.
        integer, dimension(sNy,sNx) :: sub_mask_correct
        integer, dimension(sNy,sNx) :: sub_temp_correct

        integer s_global_col, s_global_row, row, col
        integer global_ob_index, global_row, global_col
        integer correct
        integer, parameter :: print_correctness = 0

        common /sub_arr_info/ sub_masks1, sub_masks2, global_ob_index


c       Initialize mask_correct and temp_correct
        do row=1, sNy, 1
          do col=1, sNx, 1
            sub_mask_correct(row,col) = -1
            sub_temp_correct(row,col) = -1
          end do
        end do

        s_global_col = (global_ob_index/sNx) * sNx + 1
        s_global_row = (global_ob_index/Nx) + 1
        !print*, "s_global_col",s_global_col
        !print*, "s_global_row",s_global_row

        correct = 1

c       Looping through sub arrays and comparing each element to corresponding elements of global arrays
        do row=1, sNy, 1
          do col=1, sNx, 1
            global_col = int(s_global_col + col-1)
            global_row = int(s_global_row + row-1)
            !print *, global_row, global_col
c           Comparing sub mask to corresponding piece of global mask
            if (int(sub_masks1(1,row, col)) .ne.
     &          global_ob_masks(1,global_row, global_col)) then
              correct = 0
              sub_mask_correct(row, col) = 0
            end if
            if (int(sub_masks2(1,row, col)) .ne.
     &          global_ob_masks(2,global_row, global_col)) then
              correct = 0
              sub_mask_correct(row, col) = 0
            end if

c           Comparing sub temp to corresponding piece of global temp which  is defined by the function global row * global col
            if (sub_masks1(2,row, col) .ne. real(global_row*global_col))
     &      then
              correct = 0
              sub_temp_correct(row, col) = 0
            end if
            if (sub_masks2(2,row, col) .ne. real(global_row*global_col))
     &      then

              correct = 0
              sub_temp_correct(row, col) = 0
            end if

          end do
        end do

c        if (correct .eq. 0) then
c          print *, "submask1 incorrect: "
c          call PRINT_FLOAT_ARR(sub_masks1(1,1:,1:),sNy*sNx)
c          print *, "from global:"
c          call PRINT_INT_ARR(global_ob_masks(1,
c     &                         s_global_row:s_global_row+sNy-1,
c     &                         s_global_col:s_global_col+sNx-1),
c     &                         sNy*sNx)
c        end if

c       correct == 1 indicates sub arrays were correctly generated.
c       correct == 0 indicates sub arrays are incorrect.
        if (correct .eq. 1) then
          print *, "TEST_SUBMASKS: correct sub_mask and sub_temp"
        else
          print *, "TEST_SUBMASKS: INCORRECTLY DISTRIBUTED"
c          print *, "TEST_SUBMASKS: sub_mask: 1 for correct. 0 for ",
c     &             "incorrect: "
c          call PRINT_INT_ARR(sub_mask_correct, Ny*Nx)
c          print *, "TEST_SUBMASKS: sub_temp: 1 for correct. 0 for ",
c     &             "incorrect: "
c          call PRINT_INT_ARR(sub_temp_correct, Ny*Nx)
        end if
      end


      SUBROUTINE CREATE_FILENAME(s_global_index, prefix, is_child,
     &                           filename)
c       Creating filenames for all binary files.

        include 'SIZE.h'
        integer is_child, s_global_index
        character*50, intent(out) :: filename

        integer subtile_x, subtile_y
        character*3 str_px, str_py
        character*24 prefix

c       tile number along x and y axis:
        subtile_x = s_global_index/sNx + 1
        subtile_y = ((s_global_index/Nx + 1)/sNy) + 1 !row # divided by sNy then + 1

c       If creating a filename for a child process, create a name in the format: T.001.002.bin
c       T for temp or mask for mask
c       0001 for index 1 in x axis
c       0002 for index 2 in x axis
        if (is_child .eq. 1) then
          write(str_px, "(I3.3)") subtile_x
          write(str_py, "(I3.3)") subtile_y

          write(filename, "(5A)") prefix, str_px(1:3),
     &                            ".", str_py(1:3), ".bin"

        else
c         Creating a filename for global arrays
          write(filename, "(2A)") prefix, "_global.bin"

        end if
      end


      SUBROUTINE WRITE_SUB_BIN()
c       Writing all sub arrays that each process receives from the Parent to binary files
        use constants
        include 'SIZE.h'
        real sub_masks1(num_info_masks, sNy, sNx) !open boundary 1
        real sub_masks2(num_info_masks, sNy, sNx) !open boundary 2
        integer global_ob_index

        !real, dimension(sNy,sNx) :: read_sub_temp
        real, dimension(sNy,sNx) :: read_sub_mask

        integer, parameter :: print_write = 1
        integer, parameter :: print_read = 0

        common /sub_arr_info/ sub_masks1, sub_masks2, global_ob_index

        character*24 prefix
        character*50 filename

c       Writing sub mask of ob1 with global_ob_index to binary file:
c       CREATE_FILENAME params: s_global_index, prefix, is_child, filename
        prefix = "./output_masks/ob1.mask."
        call CREATE_FILENAME(global_ob_index, prefix, 1, filename)
        open(1, FILE=filename, FORM="unformatted", ACCESS="stream")
        write(1) sub_masks1(1,1:,1:)
        close(1)
        if (print_write .eq. 1) then
          print *, "sub mask for ob1 with global_ob_index",
     &              global_ob_index,"is written to ",filename
        end if

        if (print_read .eq. 1) then
c         DOUBLE CHECKING: Reading binary files to make sure sub temp is correctly written
          open(1, FILE=filename, FORM="unformatted", ACCESS="stream")
          print *, "Reading sub mask for ob1 with global_ob_index",
     &              global_ob_index,"back from ",filename
          read(1) read_sub_mask
          call PRINT_FLOAT_ARR(read_sub_mask, sNx*sNy)
          close(1)
        end if

c       Writing sub_ temp of ob1 with global_ob_index to binary file
        prefix = "./output_masks/ob1.temp."
        call CREATE_FILENAME(global_ob_index, prefix, 1, filename)
        open(2, FILE=filename, FORM="unformatted", ACCESS="stream")
        write(2) sub_masks1(2,1:,1:)
        close(2)

        if (print_write .eq. 1) then
          print *, "sub temp for ob1 with global_ob_index",
     &              global_ob_index,"is written to ",filename
        end if


        prefix = "./output_masks/ob2.mask."
        call CREATE_FILENAME(global_ob_index, prefix, 1, filename)
        open(3, FILE=filename, FORM="unformatted", ACCESS="stream")
        write(3) sub_masks2(1,1:,1:)
        close(3)
        if (print_write .eq. 1) then
          print *, "sub mask for ob2 with global_ob_index",
     &              global_ob_index,"is written to ",filename
        end if

        prefix = "./output_masks/ob2.temp."
        call CREATE_FILENAME(global_ob_index, prefix, 1, filename)
        open(4, FILE=filename, FORM="unformatted", ACCESS="stream")
        write(4) sub_masks2(2,1:,1:)
        close(4)

        if (print_write .eq. 1) then
          print *, "sub temp for ob2 with global_ob_index",
     &              global_ob_index,"is written to ",filename
        end if

      end
