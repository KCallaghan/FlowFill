!This version with an area implementation seems to be working. 
!This version of the code sorts cells by priority before moving water
!It is significantly faster than just running through the array start to finish. 


!Subroutines are used to sort the array from min to max:

subroutine Merge(A,NA,B,NB,C,NC)
 
   integer, intent(in) :: NA,NB,NC         ! Normal usage: NA+NB = NC
   integer, intent(in out) :: A(NA)        ! B overlays C(NA+1:NC)
   integer, intent(in)     :: B(NB)
   integer, intent(in out) :: C(NC)
 
   integer :: I,J,K
 
   I = 1; J = 1; K = 1;
   do while(I <= NA .and. J <= NB)
      if (A(I) <= B(J)) then
         C(K) = A(I)
         I = I+1
      else
         C(K) = B(J)
         J = J+1
      endif
      K = K + 1
   enddo
   do while (I <= NA)
      C(K) = A(I)
      I = I + 1
      K = K + 1
   enddo
   return
 
end subroutine merge
 
recursive subroutine MergeSort(A,N,T,indices)           !A is the hz array, N is the length of the array (or a shorter length), T is an empty array, indices is the indices array to also be sorted. 
 
   integer, intent(in) :: N
   integer, dimension(N), intent(in out) :: A
   real,dimension(2,N) :: indices
   integer, dimension((N+1)/2), intent (out) :: T
 
   integer :: NA,NB,V
 
   if (N < 2) return
   if (N == 2) then
      if (A(1) > A(2)) then                 !Both the hz array with topo + h info and an array recording the indices are sorted
         V = A(1)
         A(1) = A(2)
         A(2) = V

         V_row = indices(1,1)
         v_col = indices(2,1)
         indices(1,1) = indices(1,2)
         indices(1,2) = V_row

         indices(2,1) = indices(2,2)
         indices(2,1) = V_col

      endif
      return
   endif      
   NA=(N+1)/2
   NB=N-NA
 
   call MergeSort(A,NA,T,indices)
   call MergeSort(A(NA+1),NB,T,indices)
 
   if (A(NA) > A(NA+1)) then
      T(1:NA)=A(1:NA)
      call Merge(T,NA,A(NA+1),NB,A,N)
   endif
   return
 
end subroutine MergeSort



!*************************************************************************************************************


program surface_water

use netcdf
use mpi

implicit none

integer :: i,j,n2,n3,ncid,varid,error,iret,col,row,converged,counter

integer :: nmax,n,counter1,counter2

integer*8 :: ntotal

integer ierr,pid,numtasks,tasktype,status(MPI_STATUS_SIZE),rc,columntype,columntypeint

real :: watercount,upvalue,downvalue,leftvalue,rightvalue,water,lost,water1,water2

real diff_total,maxdiff,dx,dy,xn,xs

character*20 :: surfdatadir,time_start,initdatadir
character*100 :: filetopo_start,filemask

real,allocatable,dimension(:,:) :: topo,h,hz,hold,diff,topo_read,h_values,&
hz_read,h_read,hold_read,mask,mask_read

real,allocatable,dimension(:) :: hz_1D,area

integer :: Merge_size
integer,allocatable,dimension(:)  :: T

real,allocatable,dimension(:,:) :: arr

REAL,PARAMETER :: UNDEF = -1.0E+7
REAL(KIND=8) :: SEDGE
REAL(KIND=8),PARAMETER :: pi=3.141592653589793D0+0
real :: dltxy = 120 !there are 120 30 arc-second pieces in one degree

integer,allocatable :: domblock(:),domblocksmall(:),domblockint(:),nini(:),nend(:)



call MPI_INIT(ierr)       !Initialise MPI
if (ierr .ne. MPI_SUCCESS) then       !Error catching - this part should never run.
    print *,'Error starting MPI program. Terminating.'
    call MPI_ABORT(MPI_COMM_WORLD,rc,ierr)
end if

!setup data ***************************************************************************************************


n2 = 2000            !number of columns in the rotated version of madagascar. Fortran thinks these are ROWS
n3 = 1000            !number of rows in the rotated version of madagascar. Fortran thinks these are COLUMNS.
 
converged = 0          !initialising values to 0
watercount = 0.
upvalue = 0.
downvalue = 0.
leftvalue = 0.
rightvalue = 0.
counter = 0
lost = 0.

counter1=0
counter2=0

time_start = 'Mad_021000'                    !Folder names
surfdatadir = 'surfdata/'
initdatadir='initdata/'


filetopo_start = trim(surfdatadir)//trim(time_start)//'_topo_rotated.nc'         !filenames
filemask = filetopo_start!trim(initdatadir)//trim(time_start)//'_mask.nc'       !For now just using the topo as the mask

SEDGE=-60               !Used for cell area calculations. This is the southern-most latitude of your area

call MPI_COMM_RANK(MPI_COMM_WORLD,pid,ierr)
call MPI_COMM_SIZE(MPI_COMM_WORLD,numtasks,ierr)
print *,'Number of tasks=',numtasks,'My rank=',pid

allocate(nini(1:numtasks-1))
allocate(nend(1:numtasks-1))
allocate(domblock(1:numtasks-1))
allocate(domblockint(1:numtasks-1))
allocate(domblocksmall(1:numtasks-1))

if(pid.eq.0) then

!import data:

    allocate(h(n2,n3))
    allocate(hold(n2,n3))
    h(:,:) = 1.
    hold = h

    allocate(topo(n2,n3))
     
      iret = nf90_open(filetopo_start,0,ncid) !reading in the topo
      call check_err(iret)

      iret = nf90_inq_varid(ncid,'value',varid)
      call check_err(iret)

      iret = nf90_get_var(ncid,varid,topo)
      call check_err(iret)

      iret = nf90_close(ncid)
      call check_err(iret)


    where(topo .le. UNDEF) topo = 0. !change undefined cells to 0


    allocate(mask(n2,n3)) 
     
      iret = nf90_open(filemask,0,ncid) !reading in the mask
      call check_err(iret)

      iret = nf90_inq_varid(ncid,'value',varid)
      call check_err(iret)

      iret = nf90_get_var(ncid,varid,mask)
      call check_err(iret)

      iret = nf90_close(ncid)
      call check_err(iret)

    where(mask.eq.0) topo=0.
    where(mask.eq.0) h=0.

endif

call MPI_TYPE_CONTIGUOUS(numtasks-1,MPI_INTEGER,tasktype,ierr) !creates a contiguous datatype. The new datatype is called tasktype
  
call MPI_Type_commit(tasktype,ierr)

!divide the domain among tasks:

if(pid .eq. 0) then
    write(6,*) 'PID = 0'
    call dividedomain(n2,n3,numtasks,nini,filemask,ntotal) !Divides up the work equally among all of the ranks, by number of defined land cells.
    write (6,*) 'Dividedomain done'


    nend(numtasks-1) = n3+1 !define where each task must finish
    do n=2,numtasks-1
        nend(n-1) = nini(n) -1 !moving everything in the list along one space - the column to end in for each task.
    end do

      
    do n=1,numtasks-1
        call MPI_send(nini(1),1,tasktype,n,1,MPI_COMM_WORLD,ierr) !because only PID=0 has these values right now, so we have to send them out. 
        call MPI_send(nend(1),1,tasktype,n,20,MPI_COMM_WORLD,ierr) 
    end do

else
    call MPI_recv(nini(1),1,tasktype,0,1,MPI_COMM_WORLD,status,ierr) !receive what was sent above.
    call MPI_recv(nend(1),1,tasktype,0,20,MPI_COMM_WORLD,status,ierr)
endif   
write(6,*)'sent here',pid


do n=1,numtasks-1
    nmax = nend(n) - nini(n) + 4 !max number of columns we have in this group
    print *,nmax

    call MPI_TYPE_CONTIGUOUS(n2*nmax,MPI_REAL,domblock(n),ierr)
    call MPI_type_commit(domblock(n),ierr)

    call MPI_TYPE_CONTIGUOUS(n2*nmax,MPI_INTEGER1,domblockint(n),ierr)
    call MPI_type_commit(domblockint(n),ierr)

    nmax = nmax-3
    call MPI_TYPE_CONTIGUOUS(n2*nmax,MPI_REAL,domblocksmall(n),ierr)
    call MPI_type_commit(domblocksmall(n),ierr)
end do


call MPI_TYPE_CONTIGUOUS(n2,MPI_REAL,columntype,ierr)
call MPI_type_commit(columntype,ierr)


call MPI_TYPE_CONTIGUOUS(n2,MPI_INTEGER1,columntypeint,ierr)
call MPI_type_commit(columntypeint,ierr)



dy = 6370000.*pi/(180.*dltxy) !radius of the earth * pi / number of possible cells in the y-direction. This should equal the height of each cell in the N-S direction.
dx=dy
  
if(pid .gt. 0) then
   !   nmax = nend(pid) - nini(pid) +1
  
    allocate(area(n2))
    
    do j=1,n2  !changing area of cell depending on its latitude. Not totally sure what each of the different variables here represents...
   
        xs = (float(2*(j-1))/(dltxy*2.)+SEDGE)*pi/180. !Latitude on southern cell edge in radians
        xn = (float(2*(j+1))/(dltxy*2.)+SEDGE)*pi/180. !latitude on northern cell edge in radians
        area(j) = dy * 6370000.*(sin(xn)-sin(xs))/2. !final cell area for that latitude: trapezoid dy * dx
    end do

end if
  
  
if(pid.eq.0) then                    !send out topo, mask, h, hold. 
    do n=1,numtasks-2
        call MPI_send(topo(1,nini(n)-1),1,domblock(n),n,1,MPI_COMM_WORLD,ierr)
        call MPI_send(mask(1,nini(n)-1),1,domblock(n),n,2,MPI_COMM_WORLD,ierr)
        call MPI_send(h(1,nini(n)-1),1,domblock(n),n,5,MPI_COMM_WORLD,ierr)
        call MPI_send(hold(1,nini(n)-1),1,domblock(n),n,10,MPI_COMM_WORLD,ierr)
    end do

  
    call MPI_send(topo(1,nini(numtasks-1)-1),1,domblocksmall(numtasks-1),numtasks-1,1,MPI_COMM_WORLD,ierr)
    call MPI_send(mask(1,nini(numtasks-1)-1),1,domblocksmall(numtasks-1),numtasks-1,2,MPI_COMM_WORLD,ierr)
    call MPI_send(h(1,nini(numtasks-1)-1),1,domblocksmall(numtasks-1),numtasks-1,5,MPI_COMM_WORLD,ierr)
    call MPI_send(hold(1,nini(numtasks-1)-1),1,domblocksmall(numtasks-1),numtasks-1,10,MPI_COMM_WORLD,ierr) 

    call MPI_send(topo(1,1),1,columntype,numtasks-1,1,MPI_COMM_WORLD,ierr)
    call MPI_send(mask(1,1),1,columntype,numtasks-1,2,MPI_COMM_WORLD,ierr)
    call MPI_send(h(1,1),1,columntype,numtasks-1,5,MPI_COMM_WORLD,ierr)
    call MPI_send(hold(1,1),1,columntype,numtasks-1,10,MPI_COMM_WORLD,ierr)

    call MPI_send(topo(1,2),1,columntype,numtasks-1,1,MPI_COMM_WORLD,ierr)
    call MPI_send(mask(1,2),1,columntype,numtasks-1,2,MPI_COMM_WORLD,ierr)
    call MPI_send(h(1,2),1,columntype,numtasks-1,5,MPI_COMM_WORLD,ierr)
    call MPI_send(hold(1,2),1,columntype,numtasks-1,10,MPI_COMM_WORLD,ierr)


    deallocate(mask)

else
    nmax = nend(pid) - nini(pid) +4
    allocate(topo_read(n2,nmax))
    allocate(mask_read(n2,nmax))
    allocate(h_read(n2,nmax))
    allocate(hold_read(n2,nmax))

    if(pid.lt.numtasks-1)then

        call MPI_recv(topo_read(1,1),1,domblock(pid),0,1,MPI_COMM_WORLD,status,ierr)
        call MPI_recv(mask_read(1,1),1,domblock(pid),0,2,MPI_COMM_WORLD,status,ierr)
        call MPI_recv(h_read(1,1),1,domblock(pid),0,5,MPI_COMM_WORLD,status,ierr)
        call MPI_recv(hold_read(1,1),1,domblock(pid),0,10,MPI_COMM_WORLD,status,ierr)

    else
        call MPI_recv(topo_read(1,1),1,domblocksmall(pid),0,1,MPI_COMM_WORLD,status,ierr)
        call MPI_recv(mask_read(1,1),1,domblocksmall(pid),0,2,MPI_COMM_WORLD,status,ierr)
        call MPI_recv(h_read(1,1),1,domblocksmall(pid),0,5,MPI_COMM_WORLD,status,ierr)
        call MPI_recv(hold_read(1,1),1,domblocksmall(pid),0,10,MPI_COMM_WORLD,status,ierr)

        call MPI_recv(topo_read(1,nmax-1),1,columntype,0,1,MPI_COMM_WORLD,status,ierr)
        call MPI_recv(mask_read(1,nmax-1),1,columntype,0,2,MPI_COMM_WORLD,status,ierr)
        call MPI_recv(h_read(1,nmax-1),1,columntype,0,5,MPI_COMM_WORLD,status,ierr)
        call MPI_recv(hold_read(1,nmax-1),1,columntype,0,10,MPI_COMM_WORLD,status,ierr)

        call MPI_recv(topo_read(1,nmax),1,columntype,0,1,MPI_COMM_WORLD,status,ierr)
        call MPI_recv(mask_read(1,nmax),1,columntype,0,2,MPI_COMM_WORLD,status,ierr)
        call MPI_recv(h_read(1,nmax),1,columntype,0,5,MPI_COMM_WORLD,status,ierr)
        call MPI_recv(hold_read(1,nmax),1,columntype,0,10,MPI_COMM_WORLD,status,ierr)

    endif
endif


!done with data setup ************************************************************************************************


diff_total = 0

if(pid.eq.0)then
    allocate(h_values(n2,n3))
    h_values = h
else
    nmax = nend(pid) - nini(pid) +3  
 
    Merge_size = n2*nmax           !Note that this may be faster with a smaller value. Largely because the actual sort takes longer, although the number of iterations may sometimes also be higher despite that not making much sense
    !Also note if unsure - a smaller value may only partially sort the array, but a too-large value may fill with zeros! 

    allocate(T((Merge_size+1)/2))

print *,'nmax',nmax,'nini',nini(pid),'nend',nend(pid),'pid',pid,'shape',shape(topo_read)
endif



MAIN: do while(converged .eq. 0)           !Main loop for moving water
    watercount = 0.

    counter = counter + 1


    if(counter .eq.5000)then!.gt.1 .and. diff_total.lt. 1) then           !select threshold here. Consider a better way to threshold or a max number of iterations if it isn't reaching that. BUT note it'll sometimes level out for a while and then still be able to take a big jump to improvement, so it's not so simple as just looking for when it levels out! 
        print *,'success',diff_total
        converged = 1
    endif



    if (pid .eq. 0) then
 

        if (mod(counter,100) .eq. 0) then       !Doing part-way filewrites. Maybe do these less often for global.   
  
            open(23,file = 'Mad_parallel.dat',form='unformatted',access='stream')!access='direct',recl=n2*n3)!access='stream')!do the final write - create the file
            print *,'the file has been opened'
            write(23,rec=1)((h_values(i,j),i=1,n2),j=1,n3) !and write it to file
            print *,'the file has been written'
            close(23)
            print *,'written to file',pid

        endif
      
        print *,'counter',counter
        print *,'max',diff_total
  

    else   !Any PID that is not 0
        hold_read = h_read   !To compare the change for thresholding
        nmax = nend(pid) - nini(pid) +3  


        allocate(hz_read(n2,nmax+1))
        allocate(diff(n2,nmax+1))

        hz_read = topo_read+h_read


        allocate(hz_1D(n2*nmax))  !Convert the array to a 1D for the mergesort

        do i=1,n2
            hz_1D(((i-1)*nmax)+1:i*nmax) = hz_read(i,:)
        end do


        allocate(arr(2,n2*nmax))    !Create an array of the indices, for picking which cells to process first
        do row=1,n2
            arr(1,((row-1)*nmax)+1:row*nmax)=row
        end do

        do row=1,n2
            do col = 1,nmax
                arr(2,(row-1)*nmax+col)=col
            end do
        end do


        call MergeSort(hz_1D,Merge_size,T,arr)  !Sort to obtain order for processing the array

!print *,arr(2,:)
!arr(2) is until nmax, arr(1) is until n2

        COLS1: do i=1,nmax
            ROWS1: do j=1,n2
!counter1 = counter1+1
                row = arr(1,n2*nmax-(j-1)-(i-1)*n2) !get the next item in the sorted list to be processed. 
                col = arr(2,n2*nmax-(j-1)-(i-1)*n2)

                if(col.ge.nmax-1)then !Doing the end two columns separately, so skip them here
!counter1 = counter1-1
                    CYCLE
                endif

                if(col.le.2)then !Doing the end two columns separately, so skip them here
!counter1 = counter1-1
                    CYCLE
                endif

                if(row.eq.1 .or.row.eq.n2)then !Doing the end two columns separately, so skip them here
!counter1 = counter1-1
                    CYCLE
                endif


                if(mask_read(row,col) .eq. 0) then !h is 0 over the ocean
                    h_read(row,col) = 0
                endif


            if(h_read(row,col) .eq. 0) then !skip cells with no water
                CYCLE

            elseif(hz_read(row,col) .le. hz_read(row+1,col) .and. hz_read(row,col) & 
            .le. hz_read(row-1,col) .and. hz_read(row,col) .le. hz_read(row,col+1) .and. &
            hz_read(row,col) .le. hz_read(row,col-1)) then !skip cells whos neighbours are all higher
                CYCLE
            else

                upvalue = hz_read(row,col) - hz_read(row-1,col) !Get steepest slope direction
                downvalue = hz_read(row,col) - hz_read(row+1,col)
                leftvalue = hz_read(row,col) - hz_read(row,col-1)
                rightvalue = hz_read(row,col) - hz_read(row,col+1)
            
                if(max(upvalue,downvalue,leftvalue,rightvalue) .le. 0.) then  !skip cells with all negative slopes (this should never happen, should have already been skipped above)
                    CYCLE   
                elseif(max(upvalue,downvalue,leftvalue,rightvalue) .eq. upvalue) then
                    water = min(h_read(row,col)*area(row),upvalue*area(row)/2.)     !water is the minumum of the total water available, or of half the difference between 2 cells
                    water1 = water/area(row-1)   !Area calculation - something needs to be fixed
                    water2 = water/area(row)
                    h_read(row,col) = h_read(row,col) - water2
                    h_read(row-1,col) = h_read(row-1,col) + water1


                    hz_read(row,col) = hz_read(row,col) - water2
                    hz_read(row-1,col) = hz_read(row-1,col) + water1

                elseif(max(upvalue,downvalue,leftvalue,rightvalue) .eq. downvalue) then
                    water = min(h_read(row,col)*area(row),downvalue*area(row)/2.)       !water is the minumum of the total water available, or of half the difference between 2 cells
                    water1 = water/area(row+1)   !Area calculation - something needs to be fixed
                    water2 = water/area(row)

                    h_read(row,col) = h_read(row,col) - water2
                    h_read(row+1,col) = h_read(row+1,col) + water1

                    hz_read(row,col) = hz_read(row,col) - water2
                    hz_read(row+1,col) = hz_read(row+1,col) + water1

                elseif(max(upvalue,downvalue,leftvalue,rightvalue) .eq. rightvalue) then
                   water = min(h_read(row,col)*area(row),rightvalue*area(row)/2.)      !water is the minumum of the total water available, or of half the difference between 2 cells
                    water1 = water/area(row)   !Area calculation - something needs to be fixed
                    water2 = water/area(row)

                    h_read(row,col) = h_read(row,col) - water2
                    h_read(row,col+1) = h_read(row,col+1) + water1

                    hz_read(row,col) = hz_read(row,col) - water2
                    hz_read(row,col+1) = hz_read(row,col+1) + water1

                elseif(max(upvalue,downvalue,leftvalue,rightvalue) .eq. leftvalue) then
                   water = min(h_read(row,col)*area(row),leftvalue*area(row)/2.)        /area(row) !water is the minumum of the total water available, or of half the difference between 2 cells
                    water1 = water/area(row)   !Area calculation - something needs to be fixed
                    water2 = water/area(row)
                    h_read(row,col) = h_read(row,col) - water2
                    h_read(row,col-1) = h_read(row,col-1) + water1

                    hz_read(row,col) = hz_read(row,col) - water2
                    hz_read(row,col-1) = hz_read(row,col-1) + water1

                endif
            watercount = watercount + water  !I don't think watercount actually does anything right now. 
            endif

            end do ROWS1
        end do COLS1

!print *,'counter1',counter1

        if(pid.eq.1)then        !send & receive the edge columns
!print *,'try this',shape(h_read),nmax
            call MPI_recv(h_read(1,nmax+1),1,columntype,pid+1,8,MPI_COMM_WORLD,status,ierr)
            call MPI_recv(h_read(1,nmax),1,columntype,pid+1,9,MPI_COMM_WORLD,status,ierr)
            call MPI_send(h_read(1,3),1,columntype,numtasks-1,8,MPI_COMM_WORLD,ierr)
            call MPI_send(h_read(1,2),1,columntype,numtasks-1,9,MPI_COMM_WORLD,ierr)

           call MPI_recv(hz_read(1,nmax+1),1,columntype,pid+1,8,MPI_COMM_WORLD,status,ierr)
            call MPI_recv(hz_read(1,nmax),1,columntype,pid+1,9,MPI_COMM_WORLD,status,ierr)
            call MPI_send(hz_read(1,3),1,columntype,numtasks-1,8,MPI_COMM_WORLD,ierr)
            call MPI_send(hz_read(1,2),1,columntype,numtasks-1,9,MPI_COMM_WORLD,ierr)

        elseif(pid.eq.numtasks-1)then
            call MPI_send(h_read(1,3),1,columntype,pid-1,8,MPI_COMM_WORLD,ierr)   
            call MPI_send(h_read(1,2),1,columntype,pid-1,9,MPI_COMM_WORLD,ierr)
            call MPI_recv(h_read(1,nmax+1),1,columntype,1,8,MPI_COMM_WORLD,status,ierr)
            call MPI_recv(h_read(1,nmax),1,columntype,1,9,MPI_COMM_WORLD,status,ierr)

           call MPI_send(hz_read(1,3),1,columntype,pid-1,8,MPI_COMM_WORLD,ierr)   
            call MPI_send(hz_read(1,2),1,columntype,pid-1,9,MPI_COMM_WORLD,ierr)
            call MPI_recv(hz_read(1,nmax+1),1,columntype,1,8,MPI_COMM_WORLD,status,ierr)
            call MPI_recv(hz_read(1,nmax),1,columntype,1,9,MPI_COMM_WORLD,status,ierr)

        else
            if(mod(pid,2).eq.0)then
                call MPI_recv(h_read(1,nmax+1),1,columntype,pid+1,8,MPI_COMM_WORLD,status,ierr)
                call MPI_recv(h_read(1,nmax),1,columntype,pid+1,9,MPI_COMM_WORLD,status,ierr)
                call MPI_send(h_read(1,3),1,columntype,pid-1,8,MPI_COMM_WORLD,ierr)
                call MPI_send(h_read(1,2),1,columntype,pid-1,9,MPI_COMM_WORLD,ierr)

                call MPI_recv(hz_read(1,nmax+1),1,columntype,pid+1,8,MPI_COMM_WORLD,status,ierr)
                call MPI_recv(hz_read(1,nmax),1,columntype,pid+1,9,MPI_COMM_WORLD,status,ierr)
                call MPI_send(hz_read(1,3),1,columntype,pid-1,8,MPI_COMM_WORLD,ierr)
                call MPI_send(hz_read(1,2),1,columntype,pid-1,9,MPI_COMM_WORLD,ierr)

            else
                call MPI_send(h_read(1,3),1,columntype,pid-1,8,MPI_COMM_WORLD,ierr)
                call MPI_send(h_read(1,2),1,columntype,pid-1,9,MPI_COMM_WORLD,ierr)
                call MPI_recv(h_read(1,nmax+1),1,columntype,pid+1,8,MPI_COMM_WORLD,status,ierr)
                call MPI_recv(h_read(1,nmax),1,columntype,pid+1,9,MPI_COMM_WORLD,status,ierr)

                call MPI_send(hz_read(1,3),1,columntype,pid-1,8,MPI_COMM_WORLD,ierr)
                call MPI_send(hz_read(1,2),1,columntype,pid-1,9,MPI_COMM_WORLD,ierr)
                call MPI_recv(hz_read(1,nmax+1),1,columntype,pid+1,8,MPI_COMM_WORLD,status,ierr)
                call MPI_recv(hz_read(1,nmax),1,columntype,pid+1,9,MPI_COMM_WORLD,status,ierr)

            endif
        endif


        if(pid.lt.numtasks) then


            COLS3: do col=nmax-1,nmax  !process the edge columns. I'm doing these unsorted for now, but could probably place them in a new array when they are skipped above and do them sorted. 
                ROWS3:do row=2,n2
counter2=counter2+1

                if(row.eq.1 .or.row.eq.n2)then !Doing the end two columns separately, so skip them here
!counter1 = counter1-1
                    CYCLE
                endif


                if(mask_read(row,col) .eq. 0) then !h is 0 over the ocean
                    h_read(row,col) = 0
                endif


                    if(h_read(row,col) .eq. 0) then
                       CYCLE
    
                    elseif(hz_read(row,col) .le. hz_read(row+1,col) .and. hz_read(row,col) &
                    .le. hz_read(row-1,col) .and. hz_read(row,col) .le. hz_read(row,col+1) .and. &
                    hz_read(row,col) .le. hz_read(row,col-1)) then
    
                        CYCLE
                    else

                        upvalue = hz_read(row,col) - hz_read(row-1,col)
                        downvalue = hz_read(row,col) - hz_read(row+1,col)
                        leftvalue = hz_read(row,col) - hz_read(row,col-1)
                        rightvalue = hz_read(row,col) - hz_read(row,col+1)

                        if(max(upvalue,downvalue,leftvalue,rightvalue) .le. 0.) then         
                            CYCLE   
                        elseif(max(upvalue,downvalue,leftvalue,rightvalue) .eq. upvalue) then
                            water = min(h_read(row,col)*area(row),upvalue*area(row)/2.)       !water is the minumum of the total water available, or of half the difference between 2 cells
                            water1 = water/area(row-1)   !Area calculation - something needs to be fixed
                            water2 = water/area(row)

                            h_read(row,col) = h_read(row,col) - water2
                            h_read(row-1,col) = h_read(row-1,col) + water1

                            hz_read(row,col) = hz_read(row,col) - water2
                            hz_read(row-1,col) = hz_read(row-1,col) + water1


                        elseif(max(upvalue,downvalue,leftvalue,rightvalue) .eq. downvalue) then
                            water = min(h_read(row,col)*area(row),downvalue*area(row)/2.)     !water is the minumum of the total water available, or of half the difference between 2 cells
                            water1 = water/area(row+1)   !Area calculation - something needs to be fixed
                            water2 = water/area(row)

                            h_read(row,col) = h_read(row,col) - water2
                            h_read(row+1,col) = h_read(row+1,col) + water1

                            hz_read(row,col) = hz_read(row,col) - water2
                            hz_read(row+1,col) = hz_read(row+1,col) + water1

                        elseif(max(upvalue,downvalue,leftvalue,rightvalue) .eq. rightvalue) then
                            water = min(h_read(row,col)*area(row),rightvalue*area(row)/2.)   !water is the minumum of the total water available, or of half the difference between 2 cells
                            water1 = water/area(row)   !Area calculation - something needs to be fixed
                            water2 = water/area(row)
                          
                            h_read(row,col) = h_read(row,col) - water2
                            h_read(row,col+1) = h_read(row,col+1) + water1

                            hz_read(row,col) = hz_read(row,col) - water2
                            hz_read(row,col+1) = hz_read(row,col+1) + water1

                        elseif(max(upvalue,downvalue,leftvalue,rightvalue) .eq. leftvalue) then
                            water = min(h_read(row,col)*area(row),leftvalue*area(row)/2.)   !water is the minumum of the total water available, or of half the difference between 2 cells
                            water1 = water/area(row)   !Area calculation - something needs to be fixed
                            water2 = water/area(row)

                            h_read(row,col) = h_read(row,col) - water2
                            h_read(row,col-1) = h_read(row,col-1) + water1

                            hz_read(row,col) = hz_read(row,col) - water2
                            hz_read(row,col-1) = hz_read(row,col-1) + water1

                        endif
                        watercount = watercount + water
                    endif

                end do ROWS3
            end do COLS3
!print *,'counter2',counter2
        endif

        if(pid.eq.1)then  !Once again send & receive the edge columns
            call MPI_send(h_read(1,nmax+1),1,columntype,pid+1,8,MPI_COMM_WORLD,ierr)
            call MPI_send(h_read(1,nmax),1,columntype,pid+1,9,MPI_COMM_WORLD,ierr)
            call MPI_recv(h_read(1,3),1,columntype,numtasks-1,8,MPI_COMM_WORLD,status,ierr)
            call MPI_recv(h_read(1,2),1,columntype,numtasks-1,9,MPI_COMM_WORLD,status,ierr)

            call MPI_send(hz_read(1,nmax+1),1,columntype,pid+1,8,MPI_COMM_WORLD,ierr)
            call MPI_send(hz_read(1,nmax),1,columntype,pid+1,9,MPI_COMM_WORLD,ierr)
            call MPI_recv(hz_read(1,3),1,columntype,numtasks-1,8,MPI_COMM_WORLD,status,ierr)
            call MPI_recv(hz_read(1,2),1,columntype,numtasks-1,9,MPI_COMM_WORLD,status,ierr)

        elseif(pid.eq.numtasks-1)then
            call MPI_recv(h_read(1,3),1,columntype,pid-1,8,MPI_COMM_WORLD,status,ierr)   
            call MPI_recv(h_read(1,2),1,columntype,pid-1,9,MPI_COMM_WORLD,status,ierr)
            call MPI_send(h_read(1,nmax+1),1,columntype,1,8,MPI_COMM_WORLD,ierr)
            call MPI_send(h_read(1,nmax),1,columntype,1,9,MPI_COMM_WORLD,ierr)

            call MPI_recv(hz_read(1,3),1,columntype,pid-1,8,MPI_COMM_WORLD,status,ierr)   
            call MPI_recv(hz_read(1,2),1,columntype,pid-1,9,MPI_COMM_WORLD,status,ierr)
            call MPI_send(hz_read(1,nmax+1),1,columntype,1,8,MPI_COMM_WORLD,ierr)
            call MPI_send(hz_read(1,nmax),1,columntype,1,9,MPI_COMM_WORLD,ierr)

        else
            if(mod(pid,2).eq.0)then
                call MPI_recv(h_read(1,3),1,columntype,pid-1,8,MPI_COMM_WORLD,status,ierr)
                call MPI_recv(h_read(1,2),1,columntype,pid-1,9,MPI_COMM_WORLD,status,ierr)
                call MPI_send(h_read(1,nmax+1),1,columntype,pid+1,8,MPI_COMM_WORLD,ierr)
                call MPI_send(h_read(1,nmax),1,columntype,pid+1,9,MPI_COMM_WORLD,ierr)

                call MPI_recv(hz_read(1,3),1,columntype,pid-1,8,MPI_COMM_WORLD,status,ierr)
                call MPI_recv(hz_read(1,2),1,columntype,pid-1,9,MPI_COMM_WORLD,status,ierr)
                call MPI_send(hz_read(1,nmax+1),1,columntype,pid+1,8,MPI_COMM_WORLD,ierr)
                call MPI_send(hz_read(1,nmax),1,columntype,pid+1,9,MPI_COMM_WORLD,ierr)

            else
                call MPI_send(h_read(1,nmax+1),1,columntype,pid+1,8,MPI_COMM_WORLD,ierr)
                call MPI_send(h_read(1,nmax),1,columntype,pid+1,9,MPI_COMM_WORLD,ierr)
                call MPI_recv(h_read(1,3),1,columntype,pid-1,8,MPI_COMM_WORLD,status,ierr)
                call MPI_recv(h_read(1,2),1,columntype,pid-1,9,MPI_COMM_WORLD,status,ierr)

                call MPI_send(hz_read(1,nmax+1),1,columntype,pid+1,8,MPI_COMM_WORLD,ierr)
                call MPI_send(hz_read(1,nmax),1,columntype,pid+1,9,MPI_COMM_WORLD,ierr)
                call MPI_recv(hz_read(1,3),1,columntype,pid-1,8,MPI_COMM_WORLD,status,ierr)
                call MPI_recv(hz_read(1,2),1,columntype,pid-1,9,MPI_COMM_WORLD,status,ierr)

            endif
        endif

        diff = abs(hold_read-h_read)   !for thresholding purposes
        maxdiff = 0
        maxdiff = maxval(diff)
        deallocate(hz_read)
        deallocate(diff)
        deallocate(hz_1D)
        deallocate(arr)


    endif

call MPI_ALLREDUCE(maxdiff,diff_total,1,MPI_REAL,mpi_max,MPI_COMM_WORLD,ierr)  !to get overall threshold
   
end do MAIN

print *,'**************************************************************************************************'
print *,'this is the final nmax',nmax


if (pid .eq. 0) then
    print *,'done'

  
    do n=1,numtasks-1
        call MPI_recv(h_values(1,nini(n)),1,domblocksmall(n),n,4,MPI_COMM_WORLD,status,ierr) !receive the final values and save the file
    end do
  
 !   call MPI_recv(h_values(1,nini(numtasks-1)-1),1,domblocksmall(numtasks-1),numtasks-1,4,MPI_COMM_WORLD,status,ierr)
 
    
    open(23,file = 'Mad_parallel.dat',form='unformatted',access='stream')!access='direct',recl=n2*n3)!access='stream')!do the final write - create the file

    print *,'the file has been opened'

    write(23,rec=1)((h_values(i,j),i=1,n2),j=1,n3) !and write it to file

    print *,'the file has been written'

    close(23)

    print *,'written to file',pid

else

    call MPI_send(h_read(1,2),1,domblocksmall(pid),0,4,MPI_COMM_WORLD,ierr) !everyone sends out the final result to pid 0! 

endif

print *,'about to try deallocating',pid

!deallocate(h_values)


deallocate(topo,stat=error)
if (error.ne.0) then
    print *, 'topo error'
endif

deallocate(mask,stat=error)
if (error.ne.0) then
    print *,'mask error'
endif

deallocate(nini,stat=error)
if (error.ne.0) then
    print *,'nini error'
endif

deallocate(nend,stat=error)
if (error.ne.0) then
    print *,'nend error'
endif

deallocate(domblock,stat=error)
if (error.ne.0) then
    print *,'domblock error'
endif

deallocate(domblockint,stat=error)
if (error.ne.0) then
    print *,'domblockint error'
endif

deallocate(domblocksmall,stat=error)
if (error.ne.0) then
    print *,'domblocksmall error'
end if




print *,'done',pid


call MPI_FINALIZE(ierr)

print *,'finished',pid



end program surface_water


!********************************************************************************************************


subroutine dividedomain(n2,n3,numtasks,nini,filemask,ntotal)
  use netcdf
  implicit none

 
  integer :: n2,n3,numtasks
  integer :: nini(1:numtasks-1)
  real,allocatable,dimension(:,:) :: varread
  integer,allocatable,dimension(:) :: ncells
  integer :: iret,ncid,varid,ncount,n,j
  integer*8 :: ntotal
  character*100 :: filemask

  allocate(varread(n2,n3))

  write(6,*)'reading in the mask to divide the domain'

  iret = nf90_open(filemask,0,ncid)  !open the mask file
  call check_err(iret)
  write(6,*)'first call'

  iret = nf90_inq_varid(ncid,'value',varid) !get the ID of the value layer
  call check_err(iret)
  write(6,*) 'second call'

  iret = nf90_get_var(ncid,varid,varread) !read the actual values into the array called varread
  call check_err(iret)
  write(6,*) 'third call'

  iret = nf90_close(ncid) !close the mask file
  call check_err(iret)
  write(6,*)'fourth call'

  ntotal = count(varread>-100) !count the number of land cells. I changed this slightly since I am using mask here rather than topo; all cells with a value of 1 should be included. 

print *,'ntotal',ntotal,numtasks,'numtasks'

  allocate(ncells(n3))

  ncells = count(varread>-100,1) !The number of cells which are defined in the 1 dimension

  ncount=0

  nini(1) = 2

  n=2 !counter
  
  do j=1,n3
      ncount=ncount+ncells(j) !add the number of cells defined in the current column
      if (ncount .ge. ntotal/(numtasks-1)) then !>= total number of defined cells/number of threads
          nini(n) = j-1 !Telling it which row it needs to start in, to divide up the work equally, so each task handles the same number of cells.
         ncount = ncells(j) !reset the ncount for the next task
          n = n+1
      endif
      if(n .eq. numtasks) exit
  end do

print *,'down here in the subroutine',nini


deallocate(ncells)
deallocate(varread)

  return

end subroutine dividedomain


!********************************************************************************************************

subroutine check_err(statusnc)

  use netcdf
  integer statusnc

  if(statusnc.ne. nf90_noerr) then
      stop 'Stopped due to a catch in the check_err subroutine'
  endif

end subroutine check_err


!try outputting data from each individual MPI processor
