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

!This version of FlowFill allows the user to select whether they would like to use a constant runoff everywhere in the domain, or to use variable runoff.
!It is possible to allow infiltration, evaporation, etc, and these may be added at a later stage if there is a need for them. 
!This version of the code does not sort cells by priority before moving water, but simply iterates through the array left to right. With current coding and D8 
!direction selection, this was found to be faster as the time taken by a sort was longer than any gain in processing higher cells first. 
!It is paralellised using MPI. Note that a minimum of 3 processors should be used. 
!The only required input data is topography in a netcdf format. A runoff netcdf is optional.
!Other input values are required for starting runoff depth, x and y dimensions of your input DEM, number of processors to use, selection of method to use for ties, and thresholding value.
!We recommend using the supplied user_inputs and run_me files to make selection of these input values and running of the code easy. 
!This is intended for smaller study areas, so NO wraparound of water E-W or N-S is needed (as one would use for global modelling)! 

!Code outline:
!1) Subroutines: We start with subroutines to divide the domain for parallelisation; and an error-checking subroutine. 
!We also have a subroutine for finding the steepest direction for water to flow.
!2) variable declarations
!3) initialise changing variables - read in values supplied by the user at runtime.
!4) MPI setup - divide the domain, create needed datatypes - this is all overhead and can be more or less ignored. 
!5) Open data files: your supplied topography file. It is possible to supply a mask file, but this is not implemented in code right now. 
! For now, the entire domain is assumed to be land surface (future updates may allow you to provide a mask if your domain includes ocean cells)
!6) Data processing loop: water moves downslope until equilibrium is reached. The amount of water moving in a single iteration is min(water available, half the difference between target cell and steepest downslope cell). 
!7) The data processing loop is repeated to allow flow across boundaries between processors. 
!8) In order to avoid issues for flow accumulation, any tiny depressions that occur in regions where the water surface should actually be flat are removed.
!9) Finally, we output the new files indicating surface water thickness, and topography + water thickness. 


!SUBROUTINES *******************************************************************************************************************************************************************************************************


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

  write(6,*)'Reading in the raster map and dividing the domain among processors'

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


subroutine init_random_seed()
  integer :: i, n, clock
  integer, dimension(:), allocatable :: seed
         
  call RANDOM_SEED(size = n)
  allocate(seed(n))
         
  call SYSTEM_CLOCK(COUNT=clock)
          
  seed = clock + 37 * (/ (i - 1, i = 1, n) /)
  call RANDOM_SEED(PUT = seed)
          
  deallocate(seed)
end subroutine

!********************************************************************************************************

!get steepest slope direction subroutine:
!This subroutine takes the elevation of the current cell and of its 8 neighbours and returns
!the steepest downslope direction for water to flow. In the case of a tie, water moves preferentially
!NW, W, SW, S, SE, E, NE, then N. 
subroutine steepest(x,y,my_val,N_val,NE_val,E_val,SE_val,S_val,SW_val,W_val,NW_val,target_cell,treat_ties)

  integer, intent(in)           :: x,y
  integer                       :: ties,random_int
  integer                       :: target_cell(2),cell_num(8)
  real   , intent(in)           :: my_val,N_val,S_val,E_val,W_val, NE_val, SE_val,NW_val,SW_val
  real                          :: random
  real,allocatable,dimension(:) :: direction_array, x_array, y_array
  character(len=10) ,intent(in) :: treat_ties


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


  !If there is no downslope direction, water will not leave this cell, and we return -100 values to indicate this. 
  if(max(direction_array(1),direction_array(2),direction_array(3),direction_array(4),&
    direction_array(5),direction_array(6),direction_array(7),direction_array(8)) .le.0) then
    target_cell(1) = -100
    target_cell(2) = -100
    return
  endif

  ties = 0 !There are no ties yet, we will check below if there are
  !There is a downslope direction, so we assign the location of the steepest slope to target_cell.
  do n= 1,8
    if(max(direction_array(1),direction_array(2),direction_array(3),&
      direction_array(4),direction_array(5),direction_array(6),direction_array(7),&
      direction_array(8)) .eq. direction_array(n)) then
      ties = ties + 1 !Add to ties so we can see if there was a unique cell with that value, or multiples
      cell_num(ties) = n
      target_cell(1) = x_array(n)
      target_cell(2) = y_array(n)
    endif
  end do 

  !Check for ties and switch to a random direction if these exist:

  if(treat_ties == 'RAND')then
  if(ties .gt. 1) then
    CALL init_random_seed()         
    CALL RANDOM_NUMBER(random)
    random_int = FLOOR(ties*random) + 1  !e.g. if ties = 3, we will get a number 1, 2 or 3, which correspond to the first 3 positions in cell_num
    target_cell(1) = x_array(cell_num(random_int))
    target_cell(2) = y_array(cell_num(random_int))
  endif
  endif




  return

end subroutine steepest






!VARIABLE DECLARATIONS: *************************************************************************************************************


program surface_water

use netcdf
use mpi

implicit none

integer   :: i,j,n2,n3,ncid,varid,error,iret,col,row,converged,counter,nmax,n, Merge_size
integer   :: ierr,pid,numtasks,tasktype,status(MPI_STATUS_SIZE),rc,columntype,columntypeint
integer   :: target_cell(3),depressions,iters ,my_current_diff
integer*8 :: ntotal

real :: diff_total,maxdiff,water,water_threshold,depression_threshold,starting_water, start, finish ,time_increment,y
real :: total_adjust,max_adjust,diffsum,diffsum_total, thresh_mean,user_selected_threshold
real :: hz_me,hz_N,hz_S,hz_W,hz_E,hz_NE,hz_NW,hz_SW,hz_SE,h_N,h_S,h_W,h_E,h_NE,h_NW,h_SW,h_SE

character*100 :: filetopo,output_string,outfile,waterfile,textfile,bool_runoff,runoff_file

real,allocatable,dimension(:,:) :: topo,diff,topo_read,h_values,&
hz_read,h_read,hold_read,arr,topo_import,add,hz_values,runoff_import

real,allocatable,dimension(:) :: hz_1D,threshold_array

integer,allocatable,dimension(:)  :: T

REAL,PARAMETER :: UNDEF = -1.0E+7

integer,allocatable :: domblock(:),domblocksmall(:),nini(:),nend(:)

integer::narg,cptArg !#of arg & counter of arg
character(len=100)::name,my_name !Arg name
character(len=10)::treat_ties


logical ::supplied_runoff=.false.,supplied_file=.false.,ready_to_exit=.false.


!INITIALISE MPI: ********************************************************************************************************************

call MPI_INIT(ierr)       !Initialise MPI
if (ierr .ne. MPI_SUCCESS) then       !Error catching - this part should never run.
  print *,'Error starting MPI program. Terminating.'
  call MPI_ABORT(MPI_COMM_WORLD,rc,ierr)
end if

!SETUP DATA - CHANGE VARIABLE VALUES: ***************************************************************************************************

!These are the arguments supplied by the user at runtime. 
!User arguments should be *starting runoff depth* *DEM file* *number of columns* *number of rows* *threshold value* *output file name*, in that order.

!Check if any arguments are found
 narg=command_argument_count()
!Loop over the arguments
 if(narg>0)then
!loop across options
 do cptArg=1,narg

  call get_command_argument(cptArg,name)
     
      if(cptArg==1)then                    !runoff depth to add to the landscape. Always use a decimal point e.g. 1.0 not 1 for one metre of runoff depth.
        my_name = adjustl(name)
        Read( my_name , '(f5.5)' )  y      ! reads the value from string and assigns it to y
        starting_water=y        
     
      elseif(cptArg==2)then                !Input topography file
        filetopo = adjustl(name)

      elseif(cptArg==3)then                !number of columns in the DEM
        my_name = adjustl(name)
        Read( my_name , '(f5.0)' )  y      
        n2=y        

      elseif(cptArg==4)then                !number of rows in the DEM
        my_name = adjustl(name)
        Read( my_name , '(f5.0)' )  y      
        n3=y        

      elseif(cptArg==5)then                !threshold for when to exit the program. See the readme for more information.
        my_name = adjustl(name)
        Read( my_name , '(f8.8)' )  y      
        user_selected_threshold=y        


      elseif(cptArg==6)then                !output file name. 
        outfile = trim(adjustl(name))//'.dat'
        waterfile = trim(adjustl(name))//'_water.dat'
        textfile = trim(adjustl(name))//'_text_output.txt'

      elseif(cptArg==7)then
        bool_runoff = adjustl(name)

      elseif(cptArg==8)then
        runoff_file = adjustl(name)

      elseif(cptArg==9)then
        treat_ties = adjustl(name)
      
      endif

 end do
 end if


call cpu_time(start)

time_increment = 0.0
output_string = '_'//trim(my_name)//'m_runoff_text_output.txt'  
open (15,file=textfile) !creating a text file to store the values of iterations etc


water_threshold = starting_water/10000.  
depression_threshold = water_threshold


!Don't touch these variables:
 
converged  = 0
counter    = 0
my_current_diff = 1
thresh_mean = 0

write(15,*) 'filetopo ',filetopo 

!MPI SETUP AND BEHIND-THE-SCENES: ****************************************************************************************************************

call MPI_COMM_RANK(MPI_COMM_WORLD,pid,ierr)
call MPI_COMM_SIZE(MPI_COMM_WORLD,numtasks,ierr)
write(15,*)"you used a runoff value of ",starting_water

write(15,*) 'Number of tasks=',numtasks,'My rank=',pid

allocate(nini(1:numtasks-1))
allocate(nend(1:numtasks-1))
allocate(domblock(1:numtasks-1))
allocate(domblocksmall(1:numtasks-1))

call MPI_TYPE_CONTIGUOUS(numtasks-1,MPI_INTEGER,tasktype,ierr) !creates a contiguous datatype. The new datatype is called tasktype
  
call MPI_Type_commit(tasktype,ierr)

!divide the domain among tasks:
!The DEM is split between the processors to speed computation. 
!The D8 flow of water works by passing edge columns between processors each 
!iteration and processing these separately.

if (pid .eq. 0) then
  print *,''
  print *,'Beginning FlowFill.'
  print *,'Convergence threshold: ',user_selected_threshold
  print *,'Topography file: ',filetopo 
  call dividedomain(n2,n3,numtasks,nini,filetopo,ntotal) !Divides up the work equally among all of the ranks, by number of defined land cells.
  
  nend(numtasks-1) = n3+2 !define where each task must finish
  do n=2,numtasks-1
    nend(n-1) = nini(n) -1 !moving everything in the list along one space - the column to end in for each task.
  end do
      
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
  
  allocate(threshold_array(2000))
  threshold_array(:) = 100.

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


  allocate(topo(n2+2,n3+2))
  allocate(h_values(n2+2,n3+2))
  h_values(:,:) = starting_water

  if(bool_runoff .eq. "Y")then
    allocate(runoff_import(n2,n3))
     
    iret = nf90_open(runoff_file,0,ncid) !reading in the topo
    call check_err(iret)

    iret = nf90_inq_varid(ncid,'value',varid)
    call check_err(iret)

    iret = nf90_get_var(ncid,varid,runoff_import)
    call check_err(iret)

    iret = nf90_close(ncid)
    call check_err(iret)

    do i = 1,n2+2   !add the border around the domain to allow processing to happen in inland regions (allow water to flow off the topography). 
    do j = 1,n3+2
      if(i.eq.1)then
        if(j.eq.1)then 
          h_values(i,j) = runoff_import(i,j) 
        elseif(j.eq.n3+2) then
          h_values(i,j) = runoff_import(i,j-2)
        else
          h_values(i,j)= runoff_import(i,j-1)  !leftmost column
        endif 

      elseif(i.eq.n2+2)then 
        if(j.eq.1)then 
          h_values(i,j) = runoff_import(i-2,j) 
        elseif(j.eq.n3+2) then
          h_values(i,j) = runoff_import(i-2,j-2)
        else
          h_values(i,j) = runoff_import(i-2,j-1) !rightmost column 
        endif 
      
      elseif(j.eq.1) then  
        h_values(i,j)= runoff_import(i-1,j) !top row
      elseif(j.eq.n3+2)then 
        h_values(i,j)= runoff_import(i-1,j-2) !bottom row
      else
        h_values(i,j) = runoff_import(i-1,j-1)
      endif
    end do 
    end do 

  endif

!Water in cells at the edge of the domain may flow into or out of the domain. We do not want it to indiscriminately
!flow out if it may be able to flow inwards, so we add a one-cell border with the same elevation as the current
!edge cells. 

  do i = 1,n2+2   !add the border around the domain to allow processing to happen in inland regions (allow water to flow off the topography). 
    do j = 1,n3+2
      if(i.eq.1)then
        if(j.eq.1)then 
          topo(i,j) = topo_import(i,j) 
        elseif(j.eq.n3+2) then
          topo(i,j) = topo_import(i,j-2)
        else
          topo(i,j)= topo_import(i,j-1)  !leftmost column
        endif 

      elseif(i.eq.n2+2)then 
        if(j.eq.1)then 
          topo(i,j) = topo_import(i-2,j) 
        elseif(j.eq.n3+2) then
          topo(i,j) = topo_import(i-2,j-2)
        else
          topo(i,j) = topo_import(i-2,j-1) !rightmost column 
        endif 
      
      elseif(j.eq.1) then  
        topo(i,j)= topo_import(i-1,j) !top row
      elseif(j.eq.n3+2)then 
        topo(i,j)= topo_import(i-1,j-2) !bottom row
      else
        topo(i,j) = topo_import(i-1,j-1)
      endif
    end do 
  end do 

endif

!SEND & RECEIVE IMPORTED DATA: ****************************************************************************************

    
if(pid.eq.0) then                    !send out topo, h. 
  do n=1,numtasks-1
    call MPI_send(topo(1,nini(n)),1,domblock(n),n,1, MPI_COMM_WORLD,ierr)
    call MPI_send(h_values(1,nini(n)),1,domblock(n),n,5, MPI_COMM_WORLD,ierr)
  end do

else
  nmax = nend(pid) - nini(pid) + 4
  allocate(topo_read(n2+2,nmax))
  allocate(h_read   (n2+2,nmax))
  allocate(hold_read(n2+2,nmax))
  allocate(hz_read  (n2+2,nmax))

  call MPI_recv(topo_read(1,1),1,domblock(pid),0,1, MPI_COMM_WORLD,status,ierr)
  call MPI_recv(h_read   (1,1),1,domblock(pid),0,5, MPI_COMM_WORLD,status,ierr)

endif

!done with data setup ************************************************************************************************  

diff_total = 0.

if(pid.eq.0)then
  allocate(hz_values(n2+2,n3+2)) !to save the final data later on
  diffsum_total = 0
  diff_total = 0
  maxdiff = 0
  diffsum = 0
else
  nmax = nend(pid) - nini(pid) +3  
  allocate(diff(n2+2,nmax+1))
  diff(:,:) = 0.
endif


!START THE MAIN LOOP: **********************************************************************************************************************************************

write(15,*)'nini',nini,'nend',nend,'nmax',nmax,'pid',pid,'size',size(topo_read)
write(15,*)'starting the main loop'

!Now we are going to do the work! This loop iterates across the domain, moving water downslope in each iteration. 
!The amount of water moved is the minimum of all the water in a cell, or half the difference between that cell's 
!elevation and the elevation of its steepest downslope neighbour. 

MAIN: do while(converged .eq. 0)           !Main loop for moving water
  flush(15)
  counter = counter + 1

  if(counter .eq. 1000000)then               !unlikely to ever reach the selected threshold, almost certainly in an endless loop by this stage. If your DEM takes longer than this for a run to complete, comment this out or select a larger value.
    write(15,*)'timed out',diff_total
    converged = 1
  endif 

  if (pid .eq. 0) then
    if (mod(counter,500) .eq. 0) then       !Doing part-way filewrites. 
      ready_to_exit = .true.             !Every 500 iterations, we check to see whether we have reached that plateau that means we have finished the calculation.
      do n=1,2000
        if(abs(diff_total - threshold_array(n)) .gt. user_selected_threshold)then  !check to see if the max amount of water moving has not changed by more than the threshold in the past 20000 itertions.
          ready_to_exit = .false.
        endif
      end do
      if(ready_to_exit .eqv. .true.)then
        print *,'checking if ready_to_exit',ready_to_exit,pid
        print *,diff_total,starting_water
        if(diff_total < starting_water)then                                        !The max water moving must at least be less than the initial runoff supplied to be allowed to exit.
          print *,'exiting the main loop',diff_total
          converged = 1

        endif 
      endif
       call MPI_Bcast(converged,1,MPI_INT,0,MPI_COMM_WORLD,ierr)


    write(15,*)'counter',counter,'max',diff_total  ,'sum',diffsum_total

    call cpu_time(finish)
    time_increment = time_increment + finish-start
    write(15,*) "Time in seconds:",time_increment
    call cpu_time(start)

      do n=1,numtasks-1
        call MPI_recv(h_values(1,nini(n)),1,domblocksmall(n),n,4,MPI_COMM_WORLD,status,ierr) !receive the final values and save the file
      end do

      open(23,file = outfile,form='unformatted',access='stream')!access='direct',recl=n2*n3)!access='stream')!do the final write - create the file
      write(23)((h_values(i,j),i=1,n2+2),j=1,n3+2) !and write it to file
      close(23)
    endif
   diff_total = 0
   diffsum_total = 0   

  else   !Any PID that is not 0
    if(mod(counter,500).eq.0)then
      call MPI_Bcast(converged,1,MPI_INT,0,MPI_COMM_WORLD,ierr)
      call MPI_send(hz_read(1,2),1,domblocksmall(pid),0,4,MPI_COMM_WORLD,ierr) !everyone sends out the current result to pid 0! 
    endif

    hold_read = h_read   !To compare the change for thresholding
    nmax = nend(pid) - nini(pid) +3  

    hz_read = topo_read+h_read

    COLS1: do col=1,nmax!+1
      ROWS1: do row=1,n2+2
 
        if((pid.ne.1 .and. col.le.2) .or. (pid.ne.numtasks-1 .and. col .ge. nmax-1) )then!.ge.nmax-1)then !Doing the end two columns separately, so skip them here
          CYCLE
        elseif((pid.eq.1 .and. col.eq.1) .or. (pid .eq.numtasks-1 .and. col.ge.nmax-2) &
          .or. (row.eq.1 .or. row.eq.n2+2))then
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
             hz_read(row-1,col-1),target_cell,treat_ties)

          if(target_cell(1).eq.-100)then     !there is no downslope cell            
            CYCLE
          endif

         water = min(h_read(row,col), (hz_read(row,col) - hz_read(target_cell(1),target_cell(2)))/2)
        endif

        if(h_read(row,col) - water .lt. 0) then  !this should never happen, except for correcting floating point errors
          water = h_read(row,col)
        endif
        if(hz_read(row,col)-water .lt. topo_read(row,col))then  !this should never happen, except for correcting floating point errors
          water = hz_read(row,col)-topo_read(row,col)
        endif 
          
        h_read(row,col) = h_read(row,col) - water                      !adjust the target cell and neighbour cell, both h and hz 
        h_read(target_cell(1),target_cell(2)) = h_read(target_cell(1),target_cell(2)) + water
        hz_read(row,col) = hz_read(row,col) - water
        hz_read(target_cell(1),target_cell(2)) = hz_read(target_cell(1),target_cell(2)) + water
        
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

       if(row.eq.1 .or. row.eq.n2+2)then
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
             hz_read(row-1,col-1),target_cell,treat_ties)

            if(target_cell(1).eq.-100)then
              CYCLE
            endif

            water = min(h_read(row,col), (hz_read(row,col) - hz_read(target_cell(1),target_cell(2)))/2)

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
    diffsum = 0
    diffsum = sum(diff)

  endif

  call MPI_REDUCE(maxdiff,diff_total,1,MPI_REAL,mpi_max,0,MPI_COMM_WORLD,ierr)  !to get overall threshold

  call MPI_REDUCE(diffsum,diffsum_total,1,MPI_REAL,mpi_sum,0,MPI_COMM_WORLD,ierr)  !to get overall threshold

if(pid.eq.0)then
  thresh_mean = thresh_mean + diff_total
  if(mod(counter,10).eq.0)then
    threshold_array(my_current_diff) = thresh_mean/10
    thresh_mean = 0
    my_current_diff = my_current_diff + 1
    if(my_current_diff == 2001) then
      my_current_diff = 1
    endif
  endif
endif
   
end do MAIN


if(pid.eq.0)then   !send & receive data for the final local depression check 
  do n=1,numtasks-1
    call MPI_recv(hz_values(1,nini(n)),1,domblocksmall(n),n,4,MPI_COMM_WORLD,status,ierr) !receive the final values and save the file
    call MPI_recv(h_values(1,nini(n)),1,domblocksmall(n),n,5,MPI_COMM_WORLD,status,ierr) !receive the final values and save the file
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

        hz_me = hz_values(row,col)
        hz_N  = hz_values(row-1,col)
        hz_NE = hz_values(row-1,col+1)
        hz_E  = hz_values(row,col+1)
        hz_SE = hz_values(row+1,col+1)
        hz_S  = hz_values(row+1,col)
        hz_SW = hz_values(row+1,col-1)
        hz_W  = hz_values(row,col-1)
        hz_NW = hz_values(row-1,col-1)

        h_N  = h_values(row-1,col)
        h_NE = h_values(row-1,col+1)
        h_E  = h_values(row,col+1)
        h_SE = h_values(row+1,col+1)
        h_S  = h_values(row+1,col)
        h_SW = h_values(row+1,col-1)
        h_W  = h_values(row,col-1)
        h_NW = h_values(row-1,col-1)

    
        if(hz_me .le. hz_N .and. hz_me .le. hz_NE .and. hz_me .le. hz_E .and. hz_me .le. hz_SE .and. &
          hz_me .le. hz_S .and. hz_me .le. hz_SW .and. hz_me .le. hz_W .and. hz_me .le. hz_NW) then         !if this cell is a local depression

          if(h_N .le. 0 .and. h_NE .le. 0 .and. h_E .le. 0 .and. h_SE .le. 0 .and. h_S .le.0 .and. &   !if no neighbour cells have any water, skip it (technically this shouldn't happen)
            h_SW .le. 0 .and. h_W  .le. 0 .and. h_NW .le.0 ) then                 
            CYCLE

          else        

            if(((hz_E-hz_me).lt.depression_threshold) .and. (h_E .gt. 0)) then  !if a neighbour has water, is less than some threshold above the cell, and is more than other neighbours
              add(row,col) = hz_E                                                                        !then record the neighbour's height for later 
            endif
            if(((hz_W-hz_me).lt.depression_threshold) .and. (hz_W .gt. add(row,col)) &  !check this for all eight neighbours
              .and. (h_W .gt.0)) then
              add(row,col) = hz_W
            endif
            if(((hz_N-hz_me).lt.depression_threshold) .and. (hz_N .gt. add(row,col)) &
              .and. (h_N .gt. 0)) then
              add(row,col) = hz_N
            endif
            if(((hz_S-hz_me).lt.depression_threshold) .and. (hz_S .gt. add(row,col)) &
              .and. (h_S .gt.0)) then
              add(row,col) = hz_S
            endif
            if(((hz_SE-hz_me).lt.depression_threshold) .and. (hz_SE .gt. add(row,col)) &
              .and. (h_SE .gt.0)) then
              add(row,col) = hz_SE
            endif
            if(((hz_SW-hz_me).lt.depression_threshold) .and. (hz_SW .gt. add(row,col)) &
              .and. (h_SW .gt.0)) then
              add(row,col) = hz_SW
            endif
            if(((hz_NE-hz_me).lt.depression_threshold) .and. (hz_NE .gt. add(row,col)) &
              .and. (h_NE .gt.0)) then
              add(row,col) = hz_NE
            endif
            if(((hz_NW-hz_me).lt.depression_threshold) .and. (hz_NW .gt. add(row,col)) &
              .and. (h_NW .gt.0)) then
              add(row,col) = hz_NW

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

write(15,*)"the final max_adjust was ",max_adjust
write(15,*)"the final total_adjust was ",total_adjust
write(15,*) '("Total Time = ",f6.3," seconds.")',finish-start
           
!WRITE AND CLEAN UP DATA: *************************************************************************************************************************

write(15,*)'**************************************************************************************************'
write(15,*)'this is the final nmax',nmax


if (pid .eq. 0) then
  write(15,*)'FlowFill Complete.'
  print *,'FlowFill Complete.'
  print *,''
   
  open(23,file = outfile,form='unformatted',access='stream')!access='direct',recl=n2*n3)!access='stream')!do the final write - create the file
  write(23)((hz_values(i,j),i=1,n2+2),j=1,n3+2) !and write it to file
  close(23)

  open(25,file = waterfile,form='unformatted',access='stream')!access='direct',recl=n2*n3)!access='stream')!do the final write - create the file
  write(25)((h_values(i,j),i=1,n2+2),j=1,n3+2) !and write it to file
  close(25)

endif

call MPI_FINALIZE(ierr)

print *,'Process ',pid,' complete.'
write (15,*)'Process ',pid,' complete.'

if (pid .eq. 0) then
  print *,''
  print *,'Output file written. FlowFill finalized.'
  print *,''
endif

end program surface_water

!********************************************************************************************************
