! FlowFill fills or partially fills depressions in a landscape using a prescribed runoff value in a way that conserves water mass.

! Copyright (C) 2019 Kerry Callaghan

! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.

! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.

! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.


!cleaning up this version ahead of using for surface-water focussed paper. 
!This version is for the stormwater style runs, runoff only, same value everywhere. Code for water balance style runs (E, P-ET) still needs to be added. 
!This version of the code sorts cells by priority before moving water. This is significantly faster than just running through the array start to finish. 
!It is paralellised using MPI. Note that a minimum of 3 processors should be used. 
!The only required input data is topography in a netcdf format. Mask is optional, it should be provided if your domain includes ocean cells; otherwise it should not be provided. 
!Any provided mask should have values 1 for land cells, 0 for ocean cells. 
!This is now for smaller study areas, so NO wraparound of water E-W or N-S is needed! 

!Code outline:
!1) Subroutines: We start with subroutines for sorting the array by priority; as well as one to divide the domain for parallelisation; and an error-checking subroutine. 
!We also have subroutines for finding the steepest direction for water to flow, and for the amount of water that will move in a given iteration. 
!2) variable declarations
!3) initialise changing variables - this should be the only place for edits for future runs
!4) MPI setup - divide the domain, create needed datatypes - this is all overhead and can be more or less ignored. 
!5) Open data files: your supplied topography file. If you are supplying a mask file, this should be indicated in the changing variables above;
! otherwise the entire domain is assumed to be land surface (provide a mask if your domain includes ocean cells)
!6) Data processing loop: water moves downslope until equilibrium is reached. The amount of water moving in a single iteration is min(water available, half the difference between target cell and steepest downslope cell). 
!7) The data processing loop is repeated to allow flow across boundaries between processors. 
!8) In order to avoid issues for flow accumulation, any tiny depressions that occur in regions where the water surface should actually be flat are removed.
!9) Finally, we output the new files indicating surface water thickness, and topography + water thickness. 


!SUBROUTINES *******************************************************************************************************************************************************************************************************

!used to sort the array from min to max for more efficient processing:
subroutine Merge(A,NA,B,NB,C,NC)
  integer, intent(in)     :: NA,NB,NC   ! Normal usage: NA+NB = NC
  integer, intent(in out) :: A(NA)      ! B overlays C(NA+1:NC)
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

!used to sort the array from min to max for more efficient processing:

recursive subroutine MergeSort(A,N,T,indices)           !A is the hz array, N is the length of the array (or a shorter length), T is an empty array, indices is the indices array to also be sorted. 
  integer, intent(in)                       :: N
  integer, dimension(N), intent(in out)     :: A
  real,    dimension(2,N)                   :: indices
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

!Divide the domain into x equal regions, not including masked cells, where x is the number of processors used:


subroutine dividedomain(n2,n3,numtasks,nini,filemask,ntotal)
  use netcdf
  implicit none

 
  integer                          :: n2,n3,numtasks,nini(1:numtasks-1),iret,ncid,varid,ncount,n,j
  real,allocatable,dimension(:,:)  :: varread
  integer,allocatable,dimension(:) :: ncells
  integer*8                        :: ntotal
  character*100                    :: filemask

  allocate(varread(n2,n3))

  write(6,*)'reading in the mask to divide the domain'

  iret = nf90_open(filemask,0,ncid)  !open the mask file
  call check_err(iret)

  iret = nf90_inq_varid(ncid,'value',varid) !get the ID of the value layer
  call check_err(iret)

  iret = nf90_get_var(ncid,varid,varread) !read the actual values into the array called varread
  call check_err(iret)

  iret = nf90_close(ncid) !close the mask file
  call check_err(iret)

  ntotal = count(varread>-500) !count the number of land cells. 

  allocate(ncells(n3))

  ncells = count(varread>-500,1) !The number of cells which are defined in the 1 dimension

  ncount=0

  nini(1) = 1

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

  deallocate(ncells)
  deallocate(varread)

  return

end subroutine dividedomain


!********************************************************************************************************

!error message subroutine:

subroutine check_err(statusnc)

  use netcdf
  integer statusnc

  if(statusnc.ne. nf90_noerr) then
    stop 'Stopped due to a catch in the check_err subroutine'
  endif

end subroutine check_err


!********************************************************************************************************

!get steepest slope direction subroutine:

subroutine steepest(x,y,my_val,N_val,NE_val,E_val,SE_val,S_val,SW_val,W_val,NW_val,target_cell)

  integer, intent(in) :: x,y
  integer                :: target_cell(3)
  real   , intent(in) :: my_val,N_val,S_val,E_val,W_val, NE_val, SE_val,NW_val,SW_val
  !real                :: N,S,E,W,NE,NW,SE,SW
  real,allocatable,dimension(:) :: direction_array, x_array, y_array,test_arr



 allocate(test_arr(8))
  test_arr(:) = my_val
  test_arr(1) = N_val
  test_arr(2) = NE_val
  test_arr(3) = E_val
  test_arr(4) = SE_val
  test_arr(5) = S_val
  test_arr(6) = SW_val
  test_arr(7) = W_val
  test_arr(8) = NW_val


  allocate(direction_array(8))
  direction_array(:) = my_val
  direction_array(1) = my_val - N_val
  direction_array(2) = my_val - NE_val
  direction_array(3) = my_val - E_val
  direction_array(4) = my_val - SE_val
  direction_array(5) = my_val - S_val
  direction_array(6) = my_val - SW_val
  direction_array(7) = my_val - W_val
  direction_array(8) = my_val - NW_val

  allocate(x_array(8))
  x_array(1) = x-1
  x_array(2) = x-1
  x_array(3) = x
  x_array(4) = x+1
  x_array(5) = x+1
  x_array(6) = x+1
  x_array(7) = x
  x_array(8) = x-1

  allocate(y_array(8))
  y_array(1) = y
  y_array(2) = y+1
  y_array(3) = y+1
  y_array(4) = y+1
  y_array(5) = y
  y_array(6) = y-1
  y_array(7) = y-1
  y_array(8) = y-1




  
  if(max(direction_array(1),direction_array(2),direction_array(3),direction_array(4),&
    direction_array(5),direction_array(6),direction_array(7),direction_array(8)) .le.0) then
  !  print *,'no lower cells'
    target_cell(1) = -100
    target_cell(2) = -100
    target_cell(3) = -100
    
    return
  endif

  do n= 1,8
    if(max(direction_array(1),direction_array(2),direction_array(3),&
      direction_array(4),direction_array(5),direction_array(6),direction_array(7),&
      direction_array(8)) .eq. direction_array(n)) then
   !   print *,'lower cell',my_val,' ',test_arr(n),' ',direction_array(n)
      target_cell(1) = x_array(n)
      target_cell(2) = y_array(n)
    endif
  end do 


  target_cell(3) = 1

  equal_count = 0

!  do n=1,8
  
!  end do


  !note that this may still have more than correct water moving in some direction since the cell it has moved to, might be processed again before this cell is processed again. Ideas to correct this?
 ! if(upvalue .eq. downvalue .or. upvalue.eq.leftvalue .or.upvalue.eq.rightvalue &
  !  .or.downvalue.eq.leftvalue .or.downvalue .eq. rightvalue .or. leftvalue.eq.rightvalue) then
   ! if(upvalue.eq.downvalue .and. upvalue.eq.leftvalue .and. upvalue.eq.rightvalue) then
    !  target_cell(3) = 4
   ! elseif(upvalue.eq.downvalue .and. upvalue.eq.leftvalue .or.upvalue.eq.downvalue &
    !  .and.upvalue.eq.rightvalue .or. downvalue.eq.leftvalue.and.downvalue.eq.rightvalue) then
    !  target_cell(3) = 3 
   ! else
    !  target_cell(3) = 2
  !  endif
 ! endif

  return


end subroutine steepest



!uses outputs from subroutine steepest, gets the amount of water to move from the target cell. 
subroutine water_to_move(directions,my_val,&
                        to_val,my_water,water)

  integer, intent(in) :: directions                                       !directions here refers to the number of directions with equal steepness, used to move only part of the water in the selected direction. 
  real                :: water
  real   , intent(in) :: my_val,to_val,my_water                           !my_val and to_val are the topo+water thickness of those respective cells; my_water is the total water available to move
  

  to_value = my_val - to_val                  !height difference between the two cells


  !equal heights of water being moved, 2 adjacent cells have equal height after the move. We are not taking area into account for the small regions processed with this code. 

  water = min(my_water, to_value/2.)
  if(directions.eq.2)then
    water = water/2.
  elseif(directions.eq.3)then
    water = water/3.
  elseif(directions.eq.4)then
    water=water/4.
  endif

  return

end subroutine water_to_move



!VARIABLE DECLARATIONS: *************************************************************************************************************


program surface_water

use netcdf
use mpi

implicit none

integer   :: i,j,n2,n3,ncid,varid,error,iret,col,row,converged,counter,nmax,n, Merge_size
integer   :: ierr,pid,numtasks,tasktype,status(MPI_STATUS_SIZE),rc,columntype,columntypeint
integer   :: mask_supplied,target_cell(3),depressions,iters ,my_current_diff
integer*8 :: ntotal

real :: diff_total,maxdiff,water,water_threshold,depression_threshold,starting_water, start, finish ,time_increment,y
real :: total_adjust,max_adjust

character*100 :: filetopo,filemask,output_string,outfile,waterfile

real,allocatable,dimension(:,:) :: topo,h,hold,diff,topo_read,h_values,&
hz_read,h_read,hold_read,mask,mask_read, arr,topo_import,mask_import,add,hz_values

real,allocatable,dimension(:) :: hz_1D,threshold_array

integer,allocatable,dimension(:)  :: T

REAL,PARAMETER :: UNDEF = -1.0E+7

integer,allocatable :: domblock(:),domblocksmall(:),nini(:),nend(:)

integer::narg,cptArg !#of arg & counter of arg
character(len=100)::name,my_name !Arg name

logical ::supplied_runoff=.false.,supplied_file=.false.,ready_to_exit=.false.


!INITIALISE MPI: ********************************************************************************************************************

call MPI_INIT(ierr)       !Initialise MPI
if (ierr .ne. MPI_SUCCESS) then       !Error catching - this part should never run.
  print *,'Error starting MPI program. Terminating.'
  call MPI_ABORT(MPI_COMM_WORLD,rc,ierr)
end if

!SETUP DATA - CHANGE VARIABLE VALUES: ***************************************************************************************************

!User can change these variables as needed:

!Check if any arguments are found
 narg=command_argument_count()
!Loop over the arguments
 if(narg>0)then
!loop across options
 do cptArg=1,narg



  call get_command_argument(cptArg,name)
 
     

      if(cptArg==1)then
        my_name = adjustl(name)
        Read( my_name , '(f5.0)' )  y      ! reads the value from string and assigns it to y
        starting_water=y        
     
      elseif(cptArg==2)then
        filetopo = adjustl(name)

      elseif(cptArg==3)then
        outfile = trim(adjustl(name))//'.dat'
        waterfile = trim(adjustl(name))//'_water.dat'
      
      endif

 end do
 end if


  print *, "you used a runoff value of ",starting_water

  print *,"your supplied topography is in file ",filetopo





call cpu_time(start)


n2 = 497!932!600!469            !number of columns in the topography. Fortran thinks these are ROWS
n3 = 600!1125!900!416          !number of rows in the topography. Fortran thinks these are COLUMNS.
time_increment = 0.0
output_string = '_'//trim(my_name)//'m_runoff_text_output.txt'  
open (15,file='Argentina_text_output_'//trim(output_string)) !creating a text file to store the values of iterations etc


mask_supplied = 0 !change to 1 if you are including a mask with ocean cells. 

!starting_water = 0.2

water_threshold = starting_water/10000.  


allocate(threshold_array(2000))
threshold_array(:) = 100.

depression_threshold = water_threshold



!Don't touch these variables:
 
converged  = 0
counter    = 0
my_current_diff = 1


write(15,*) 'filetopo ',filetopo 
print *, 'filetopo ',filetopo 

if(mask_supplied .eq. 1) then
  filemask = trim(my_name)//'_mask.nc'
endif


!MPI SETUP AND BEHIND-THE-SCENES: ****************************************************************************************************************

call MPI_COMM_RANK(MPI_COMM_WORLD,pid,ierr)
call MPI_COMM_SIZE(MPI_COMM_WORLD,numtasks,ierr)
write(15,*)"you used a runoff value of ",starting_water

write(15,*) 'Number of tasks=',numtasks,'My rank=',pid

print *, 'Number of tasks=',numtasks,'My rank=',pid


allocate(nini(1:numtasks-1))
allocate(nend(1:numtasks-1))
allocate(domblock(1:numtasks-1))
allocate(domblocksmall(1:numtasks-1))


call MPI_TYPE_CONTIGUOUS(numtasks-1,MPI_INTEGER,tasktype,ierr) !creates a contiguous datatype. The new datatype is called tasktype
  
call MPI_Type_commit(tasktype,ierr)

!divide the domain among tasks:


if(pid .eq. 0) then
  call dividedomain(n2,n3,numtasks,nini,filetopo,ntotal) !Divides up the work equally among all of the ranks, by number of defined land cells.
  
  nend(numtasks-1) = n3+2 !define where each task must finish
  do n=2,numtasks-1
    nend(n-1) = nini(n) -1 !moving everything in the list along one space - the column to end in for each task.
  end do
write(15,*)'nini',nini,'nend',nend
      
  do n=1,numtasks-1
    call MPI_send(nini(1),1,tasktype,n,1, MPI_COMM_WORLD,ierr) !because only PID=0 has these values right now, so we have to send them out. 
    call MPI_send(nend(1),1,tasktype,n,20,MPI_COMM_WORLD,ierr) 
  end do

else
  call MPI_recv(nini(1),1,tasktype,0,1, MPI_COMM_WORLD,status,ierr) !receive what was sent above.
  call MPI_recv(nend(1),1,tasktype,0,20,MPI_COMM_WORLD,status,ierr)
endif   

do n=1,numtasks-1
  nmax = nend(n) - nini(n) + 4 !max number of columns we have in this group
  write(15,*) nmax,n2*nmax ,nini

  call MPI_TYPE_CONTIGUOUS((n2+2)*nmax,MPI_REAL,domblock(n),ierr)
  call MPI_type_commit(domblock(n),ierr)


  nmax = nmax-3
  call MPI_TYPE_CONTIGUOUS((n2+2)*nmax,MPI_REAL,domblocksmall(n),ierr)
  call MPI_type_commit(domblocksmall(n),ierr)
end do

call MPI_TYPE_CONTIGUOUS(n2,MPI_REAL,columntype,ierr)
call MPI_type_commit(columntype,ierr)

call MPI_TYPE_CONTIGUOUS(n2,MPI_INTEGER1,columntypeint,ierr)
call MPI_type_commit(columntypeint,ierr)



!IMPORT DATA: ***********************************************************************************************************

if(pid.eq.0) then

  allocate(hold(n2,n3))
  hold(:,:) = starting_water

  allocate(topo_import(n2,n3))
     
  iret = nf90_open(filetopo,0,ncid) !reading in the topo
  call check_err(iret)

  iret = nf90_inq_varid(ncid,'value',varid)
  call check_err(iret)

  iret = nf90_get_var(ncid,varid,topo_import)
  call check_err(iret)

  iret = nf90_close(ncid)
  call check_err(iret)

 
  where(topo_import .le. UNDEF) topo_import = 0. !change undefined cells to 0

  allocate(mask_import(n2,n3)) 
  if(mask_supplied .eq. 1) then
   
     
    iret = nf90_open(filemask,0,ncid) !reading in the mask
    call check_err(iret)

    iret = nf90_inq_varid(ncid,'value',varid)
    call check_err(iret)

    iret = nf90_get_var(ncid,varid,mask_import)
    call check_err(iret)

    iret = nf90_close(ncid)
    call check_err(iret)
    

    where(mask_import.eq.0) topo_import=0.
    where(mask_import.eq.0) hold=0.
  else
  
    mask_import = 1
  endif

  allocate(topo(n2+2,n3+2))
  allocate(mask(n2+2,n3+2))
  allocate(h(n2+2,n3+2))


  do i = 1,n2+2   !add the border around the domain to allow processing to happen in inland regions (allow water to flow off the topography). 
    do j = 1,n3+2
      if(i.eq.1)then
        mask(i,j) = 0
        h(i,j) = starting_water
        if(j.eq.1)then 
          topo(i,j) = topo_import(i,j) 
        elseif(j.eq.n3+2) then
          topo(i,j) = topo_import(i,j-2)
    !    elseif(j.eq.n3+1)then
     !     topo(i,j) = topo_import(i,j-1)
        else
          topo(i,j)= topo_import(i,j-1)  !leftmost column
        endif 

      elseif(i.eq.n2+2)then 
        mask(i,j) = 0
        h(i,j) = starting_water
        if(j.eq.1)then 
          topo(i,j) = topo_import(i-2,j) 
        elseif(j.eq.n3+2) then
          topo(i,j) = topo_import(i-2,j-2)
      !  elseif(j.eq.n3+1)then
       !   topo(i,j) = topo_import(i-2,j-1)
        else
          topo(i,j) = topo_import(i-2,j-1) !rightmost column 
        endif 
      
      elseif(j.eq.1) then  
        mask(i,j) = 0
        h(i,j) = starting_water
       ! if(i.eq.n2+1)then
        !  topo(i,j) = topo_import(i-1,j+1)
     !   else
          topo(i,j)= topo_import(i-1,j) !top row
      !  endif 
      elseif(j.eq.n3+2)then 
        mask(i,j) = 0
        h(i,j) = starting_water
    !    if(i.eq.n2+1)then
     !     topo(i,j) = topo_import(i-1,j-2)
       ! else
          topo(i,j)= topo_import(i-1,j-2) !bottom row
       ! endif 
      else
        topo(i,j) = topo_import(i-1,j-1)
        mask(i,j) = mask_import(i-1,j-1)
        h(i,j) = hold(i-1,j-1)
      endif
    end do 
  end do 

endif

!SEND & RECEIVE IMPORTED DATA: ****************************************************************************************

    
if(pid.eq.0) then                    !send out topo, mask, h, hold. 
  do n=1,numtasks-1
    call MPI_send(topo(1,nini(n)),1,domblock(n),n,1, MPI_COMM_WORLD,ierr)
    call MPI_send(mask(1,nini(n)),1,domblock(n),n,2, MPI_COMM_WORLD,ierr)
    call MPI_send(h   (1,nini(n)),1,domblock(n),n,5, MPI_COMM_WORLD,ierr)
  end do

  deallocate(mask)
  

else
  nmax = nend(pid) - nini(pid) + 4
  allocate(topo_read(n2+2,nmax))
  allocate(mask_read(n2+2,nmax))
  allocate(h_read   (n2+2,nmax))
  allocate(hold_read(n2+2,nmax))
  allocate(hz_read  (n2+2,nmax))

  call MPI_recv(topo_read(1,1),1,domblock(pid),0,1, MPI_COMM_WORLD,status,ierr)
  call MPI_recv(mask_read(1,1),1,domblock(pid),0,2, MPI_COMM_WORLD,status,ierr)
  call MPI_recv(h_read   (1,1),1,domblock(pid),0,5, MPI_COMM_WORLD,status,ierr)

endif

!done with data setup ************************************************************************************************  

diff_total = 0.


if(pid.eq.0)then
  allocate(h_values(n2+2,n3+2))  !to save the final data later on
  allocate(hz_values(n2+2,n3+2))
  h_values = h
else
  nmax = nend(pid) - nini(pid) +3  
 
  Merge_size = n2*nmax           !Note that this may be faster with a smaller value. Largely because the actual sort takes longer, although the number of iterations may sometimes also be higher despite that not making much sense
  !Also note if unsure - a smaller value may only partially sort the array, but a too-large value may fill with zeros! 

  allocate(T((Merge_size+1)/2))
endif


!START THE MAIN LOOP: **********************************************************************************************************************************************

write(15,*)'nini',nini,'nend',nend,'nmax',nmax,'pid',pid,'size',size(topo_read)
write(15,*)'starting the main loop'


MAIN: do while(converged .eq. 0)           !Main loop for moving water
    
  flush(15)
  counter = counter + 1


!  if(counter .gt. 10 .and. diff_total.lt. water_threshold)then      !select threshold here. Consider a better way to threshold or a max number of iterations if it isn't reaching that. BUT note it'll sometimes level out for a while and then still be able to take a big jump to improvement, so it's not so simple as just looking for when it levels out! 
!    write(15,*)'success',diff_total
!    converged = 1
!  endif

  if (mod(counter,500) .eq. 0)then
    ready_to_exit = .true.
    do n=1,2000
      if(abs(diff_total - threshold_array(n)) .gt. 0.005)then
        ready_to_exit = .false.
      endif
    end do
    if(ready_to_exit .eqv. .true.)then
      print *,'checking if ready_to_exit',ready_to_exit

    if(diff_total < starting_water)then 
      print *,diff_total
      print *, threshold_array
            converged = 1

    endif 

    endif
  endif


  if(counter .eq. 1000000)then               !unlikely to ever reach the selected threshold, almost certainly in an endless loop by this stage. 
    write(15,*)'timed out',diff_total
    converged = 1
  endif 


  if (pid .eq. 0) then

 
    if (mod(counter,500) .eq. 0) then       !Doing part-way filewrites. 
    write(15,*)'counter',counter,'max',diff_total  

    call cpu_time(finish)
    time_increment = time_increment + finish-start
    write(15,*) '("Time = ",f6.3," seconds.")',time_increment
    call cpu_time(start)

      do n=1,numtasks-1
        call MPI_recv(h_values(1,nini(n)),1,domblocksmall(n),n,4,MPI_COMM_WORLD,status,ierr) !receive the final values and save the file
      end do


      open(23,file = outfile,form='unformatted',access='stream')!access='direct',recl=n2*n3)!access='stream')!do the final write - create the file
      write(15,*)'the file has been opened'
      write(23,rec=1)((h_values(i,j),i=1,n2+2),j=1,n3+2) !and write it to file
      write(15,*)'the file has been written'
      close(23)
      write(15,*)'written to file',pid

    endif
      
    
  

  else   !Any PID that is not 0
    if(mod(counter,500).eq.0)then
      call MPI_send(hz_read(1,2),1,domblocksmall(pid),0,4,MPI_COMM_WORLD,ierr) !everyone sends out the current result to pid 0! 

    endif


    hold_read = h_read   !To compare the change for thresholding
    nmax = nend(pid) - nini(pid) +3  

    allocate(diff(n2+2,nmax+1))

    hz_read = topo_read+h_read

    allocate(hz_1D((n2+2)*(nmax+1)))  !Convert the array to a 1D for the mergesort

    do i=1,n2+2
      hz_1D(((i-1)*(nmax+1))+1:i*(nmax+1)) = hz_read(i,:)
    end do


    allocate(arr(2,(n2+2)*(nmax+1)))    !Create an array of the indices, for picking which cells to process first
    do row=1,n2+2
      arr(1,((row-1)*(nmax+1))+1:row*(nmax+1))=row
    end do

    do row=1,n2+2
      do col = 1,nmax+1
        arr(2,(row-1)*(nmax+1)+col)=col
      end do
    end do
 
    call MergeSort(hz_1D,Merge_size,T,arr)  !Sort to obtain order for processing the array



    COLS1: do i=1,nmax+1
      ROWS1: do j=1,n2+2
        row = arr(1,(n2+2)*(nmax+1)-(j-1)-(i-1)*(n2+2)) !get the next item in the sorted list to be processed. 
        col = arr(2,(n2+2)*(nmax+1)-(j-1)-(i-1)*(n2+2))

        if(pid.ne.1 .and. col.le.2)then!.ge.nmax-1)then !Doing the end two columns separately, so skip them here
          CYCLE
        endif

        if(pid.ne.numtasks-1 .and. col .ge.nmax-1)then
          CYCLE 
        elseif(pid .eq.numtasks-1 .and. col.gt.nmax-2)then !Doing the end two columns separately, so skip them here
            hz_read(row,col) = topo_read(row,col)
            h_read(row,col) = 0
          CYCLE
        endif

        if(pid.eq.1.and.col.lt.2) then
          hz_read(row,col) = topo_read(row,col)
          h_read(row,col) = 0
          CYCLE
        endif


        if(mask_read(row,col) .eq. 0) then !h is 0 over the ocean
          hz_read(row,col) = topo_read(row,col)
          h_read(row,col) = 0
          CYCLE
        endif

        if(h_read(row,col) .eq. 0) then !skip cells with no water
          CYCLE

        else
          


   
          call steepest(row,col,hz_read(row,col),hz_read(row-1,col),&
             hz_read(row-1,col+1),hz_read(row,col+1),hz_read(row+1,col+1),&
             hz_read(row+1,col),hz_read(row+1,col-1),hz_read(row,col-1),&
             hz_read(row-1,col-1),target_cell)

          if(target_cell(1).eq.-100)then
            
            CYCLE
          endif

         ! call water_to_move(target_cell(3),hz_read(row,col),hz_read(target_cell(1),target_cell(2)),&
         !      h_read(row,col),water)
         water = min(h_read(row,col), (hz_read(row,col) - hz_read(target_cell(1),target_cell(2)))/2)
        endif

        if(h_read(row,col) - water .lt. 0) then  !this should never happen 
    !      print *,'first if'
          water = h_read(row,col)
        endif
        if(hz_read(row,col)-water .lt. topo_read(row,col))then  !this should never happen 
     !     print *,'second if'
          water = hz_read(row,col)-topo_read(row,col)
        endif 

        if(water .lt. 0)then
      !    print *,'we have negative water',water,' ',h_read(row,col),' ',hz_read(row,col),' ',hz_read(target_cell(1),target_cell(2))
        endif 
          
  
        h_read(row,col) = h_read(row,col) - water                      !adjust the target cell and neighbour cell, both h and hz 
        h_read(target_cell(1),target_cell(2)) = h_read(target_cell(1),target_cell(2)) + water
        hz_read(row,col) = hz_read(row,col) - water
        hz_read(target_cell(1),target_cell(2)) = hz_read(target_cell(1),target_cell(2)) + water

 !       if(h_read(row,col) .lt. 0 .or. h_read(target_cell(1),target_cell(2)) .lt. 0 ) then
  !        print *, 'negative water'
   !     endif

    !    if(hz_read(row,col) .lt. topo_read(row,col) .or. hz_read(target_cell(1),target_cell(2))&
     !    .lt. topo_read(target_cell(1),target_cell(2)) )then
      !    print *, 'negative topo'
       ! endif
        
       
        if(mask_read(row,col) .eq. 0) then !h is 0 over the ocean
          hz_read(row,col) = topo_read(row,col)
          h_read(row,col) = 0
          CYCLE
        endif

        if(hz_read(row,col) .lt. topo_read(row,col))then
          print *, 'hz_read has become an impossible value',water," ",hz_read(row,col)," ",topo(row,col)
        endif

        if(hz_read(row,col) .lt. 0)then
          print *, 'hz_read has become a negative value',water," ",hz_read(row,col)," ",topo(row,col)
        endif


      end do ROWS1
    end do COLS1


    if(pid.eq.1)then        !send & receive the edge columns
      call MPI_recv(h_read(1,nmax+1),1,columntype,pid+1,8,MPI_COMM_WORLD,status,ierr)
      call MPI_recv(h_read(1,nmax),1,columntype,pid+1,9,MPI_COMM_WORLD,status,ierr)

      call MPI_recv(hz_read(1,nmax+1),1,columntype,pid+1,8,MPI_COMM_WORLD,status,ierr)
      call MPI_recv(hz_read(1,nmax),1,columntype,pid+1,9,MPI_COMM_WORLD,status,ierr)

    elseif(pid.eq.numtasks-1)then
      call MPI_send(h_read(1,3),1,columntype,pid-1,8,MPI_COMM_WORLD,ierr)   
      call MPI_send(h_read(1,2),1,columntype,pid-1,9,MPI_COMM_WORLD,ierr)

      call MPI_send(hz_read(1,3),1,columntype,pid-1,8,MPI_COMM_WORLD,ierr)   
      call MPI_send(hz_read(1,2),1,columntype,pid-1,9,MPI_COMM_WORLD,ierr)

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
 

    if(pid.lt.numtasks-1) then

      COLS3: do col=nmax-1,nmax  !process the edge columns. I'm doing these unsorted for now, but could probably place them in a new array when they are skipped above and do them sorted. 
        ROWS3:do row=1,n2+2


          if(mask_read(row,col) .eq. 0) then !h is 0 over the ocean
            hz_read(row,col) = topo_read(row,col)
            h_read(row,col) = 0
            CYCLE
          endif
                   
          if(h_read(row,col) .eq. 0) then !skip cells with no water
            CYCLE
          else
           call steepest(row,col,hz_read(row,col),hz_read(row-1,col),&
             hz_read(row-1,col+1),hz_read(row,col+1),hz_read(row+1,col+1),&
             hz_read(row+1,col),hz_read(row+1,col-1),hz_read(row,col-1),&
             hz_read(row-1,col-1),target_cell)

            if(target_cell(1).eq.-100)then
              CYCLE
            endif

            call water_to_move(target_cell(3),hz_read(row,col),hz_read(target_cell(1),target_cell(2)),&
               h_read(row,col),water)
          endif 

          if(h_read(row,col) - water .lt. 0) then
            water = h_read(row,col)
          endif
        
          if(hz_read(row,col)-water .lt. topo_read(row,col))then 
            water = hz_read(row,col)-topo_read(row,col)
          endif 
          
          h_read(row,col) = h_read(row,col) - water
          h_read(target_cell(1),target_cell(2)) = h_read(target_cell(1),target_cell(2)) + water
          hz_read(row,col) = hz_read(row,col) - water
          hz_read(target_cell(1),target_cell(2)) = hz_read(target_cell(1),target_cell(2)) + water

          if(mask_read(row,col) .eq. 0) then !h is 0 over the ocean
            hz_read(row,col) = topo_read(row,col)
            h_read(row,col) = 0
            CYCLE
          endif
        
        end do ROWS3
      end do COLS3

    endif


    if(pid.eq.1)then  !Once again send & receive the edge columns
      call MPI_send(h_read(1,nmax+1),1,columntype,pid+1,8,MPI_COMM_WORLD,ierr)
      call MPI_send(h_read(1,nmax),1,columntype,pid+1,9,MPI_COMM_WORLD,ierr)
    
      call MPI_send(hz_read(1,nmax+1),1,columntype,pid+1,8,MPI_COMM_WORLD,ierr)
      call MPI_send(hz_read(1,nmax),1,columntype,pid+1,9,MPI_COMM_WORLD,ierr)
    
    elseif(pid.eq.numtasks-1)then
      call MPI_recv(h_read(1,3),1,columntype,pid-1,8,MPI_COMM_WORLD,status,ierr)   
      call MPI_recv(h_read(1,2),1,columntype,pid-1,9,MPI_COMM_WORLD,status,ierr)
    
      call MPI_recv(hz_read(1,3),1,columntype,pid-1,8,MPI_COMM_WORLD,status,ierr)   
      call MPI_recv(hz_read(1,2),1,columntype,pid-1,9,MPI_COMM_WORLD,status,ierr)
    
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
 
    deallocate(diff)
    deallocate(hz_1D)
    deallocate(arr)


  endif

  call MPI_ALLREDUCE(maxdiff,diff_total,1,MPI_REAL,mpi_max,MPI_COMM_WORLD,ierr)  !to get overall threshold

  if(mod(counter,10).eq.0)then

    threshold_array(my_current_diff) = diff_total

    my_current_diff = my_current_diff + 1
    if(my_current_diff == 2001) then
      my_current_diff = 1
    endif
  endif
   
end do MAIN


if(pid.eq.0)then   !send & receive data for the final local depression check 
  do n=1,numtasks-1
    call MPI_recv(hz_values(1,nini(n)),1,domblocksmall(n),n,4,MPI_COMM_WORLD,status,ierr) !receive the final values and save the file
  end do
  

  do n=1,numtasks-1
    call MPI_recv(h_values(1,nini(n)),1,domblocksmall(n),n,5,MPI_COMM_WORLD,status,ierr) !receive the final values and save the file
  end do

  do n=1,numtasks-1
    call MPI_recv(topo(1,nini(n)),1,domblocksmall(n),n,6,MPI_COMM_WORLD,status,ierr) !receive the final values and save the file
  end do
  
  
else

  call MPI_send(hz_read(1,2),1,domblocksmall(pid),0,4,MPI_COMM_WORLD,ierr) !everyone sends out the final result to pid 0! 

  call MPI_send(h_read(1,2),1,domblocksmall(pid),0,5,MPI_COMM_WORLD,ierr) !everyone sends out the final result to pid 0! 

call MPI_send(topo_read(1,2),1,domblocksmall(pid),0,6,MPI_COMM_WORLD,ierr) !everyone sends out the final result to pid 0! 
endif

!***************************************************************************************************************************************************8
!This section gets rid of tiny spurious depressions created as a result of the way the algorithm works. 
!Adjacent lake cells may have an extremely small difference between them rather than perfectly flat water, which causes problems for flow routing. 
!This looks for cells which are local depressions, and have neighbours containing water, and corrects this issue. 

max_adjust = 0.
total_adjust = 0.

if(pid.eq.0)then

  allocate(add(n2+2,n3+2))
  depressions = 100
  converged = 0
  iters = 0

  do while (converged .eq. 0)

    if(depressions .eq. 0 .or. iters .eq. 10000)then !no depressions are left, or we've checked many times and the depressions are supposed to be there. 
      converged = 1
    endif 

    iters = iters + 1 
    depressions = 0
    add = 0.
    do col=2,n3+1  
      do row=2,n2+1
    
        if(hz_values(row,col) .le. hz_values(row+1,col) .and. hz_values(row,col) .le. hz_values(row-1,col) .and.&
          hz_values(row,col) .le. hz_values(row,col+1) .and. hz_values(row,col) .le. hz_values(row,col-1)) then         !if this cell is a local depression
          if(h_values(row+1,col).le.0 .and.h_values(row-1,col).le.0 &                                                 !if no neighbour cells have any water, skip it (technically this shouldn't happen)
            .and.h_values(row,col+1).le.0 .and.h_values(row,col-1).le.0) then                 
            CYCLE

          else         
            if(((hz_values(row,col+1)-hz_values(row,col)).lt.depression_threshold) .and. (h_values(row,col+1) .gt. 0)) then  !if a neighbour has water, is less than some threshold above the cell, and is more than other neighbours
              add(row,col) = hz_values(row,col+1)                                                                        !then record the neighbour's height for later 
            endif
            if(((hz_values(row,col-1)-hz_values(row,col)).lt.depression_threshold) .and. (hz_values(row,col-1) .gt. add(row,col)) &  !check this for all four neighbours
              .and. (h_values(row,col-1) .gt.0)) then
              add(row,col) = hz_values(row,col-1)
            endif
            if(((hz_values(row-1,col)-hz_values(row,col)).lt.depression_threshold) .and. (hz_values(row-1,col) .gt. add(row,col)) &
              .and. (h_values(row-1,col) .gt. 0)) then
              add(row,col) = hz_values(row-1,col)
            endif
            if(((hz_values(row+1,col)-hz_values(row,col)).lt.depression_threshold) .and. (hz_values(row+1,col) .gt. add(row,col)) &
              .and. (h_values(row+1,col) .gt.0)) then
              add(row,col) = hz_values(row+1,col)
            endif
          endif
        endif
      end do
    end do    !the 'add' array should be completely populated with changes to be made 

    do col=2,n3 +1 
      do row=2,n2+1
        if(add(row,col).gt.hz_values(row,col))then    !if the add array has a higher neighbour for you 
          if(add(row,col)-hz_values(row,col) .gt. max_adjust)then
            max_adjust = add(row,col)-hz_values(row,col)
          endif
          total_adjust = total_adjust + (add(row,col)-hz_values(row,col))

          depressions = depressions + 1               !then you were a depression, get counted
          h_values(row,col) = h_values(row,col) + (add(row,col)-hz_values(row,col))  !change h_values as well for the next iteration (shouldn't matter, but just in case. And does matter if you want to look at the output for this layer)
          hz_values(row,col)= max(hz_values(row,col),add(row,col)) !we already know add is greater but checking again and only changing hz if add is greater
          add(row,col) = 0.  !resetting add for the next iteration. 
        endif 
      end do 
    end do 
write(15,*)pid,'depressions',depressions


  end do
endif

print *, "the final max_adjust was ",max_adjust
print *, "the final total_adjust was ",total_adjust



    write(15,*) '("Total Time = ",f6.3," seconds.")',finish-start
           
!WRITE AND CLEAN UP DATA: *************************************************************************************************************************



write(15,*)'**************************************************************************************************'
write(15,*)'this is the final nmax',nmax


if (pid .eq. 0) then
  write(15,*)'done'
   
  open(23,file = outfile,form='unformatted',access='stream')!access='direct',recl=n2*n3)!access='stream')!do the final write - create the file
  write(15,*)'the file has been opened'
  write(23,rec=1)((hz_values(i,j),i=1,n2+2),j=1,n3+2) !and write it to file
  write(15,*)'the file has been written'
  close(23)
  write(15,*)'written to file',pid


   
  open(25,file = waterfile,form='unformatted',access='stream')!access='direct',recl=n2*n3)!access='stream')!do the final write - create the file
  write(15,*)'the file has been opened'
  write(25,rec=1)((h_values(i,j),i=1,n2+2),j=1,n3+2) !and write it to file
  write(15,*)'the file has been written'
  close(25)
  write(15,*)'written to file',pid



endif


write(15,*)'about to try deallocating',pid

!deallocate(h_values)
!close(15)


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

deallocate(domblocksmall,stat=error)
if (error.ne.0) then
    print *,'domblocksmall error'
end if




print *,'done',pid


call MPI_FINALIZE(ierr)

print *,'finished',pid



end program surface_water


!********************************************************************************************************


