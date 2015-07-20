Program driver

 ! --------------------------------------------------------------------------
 ! MPI Version of the Cartesian Acoustic Sun Simulator. 
 ! Copyright 2006, Shravan Hanasoge
 
 ! Hansen Experimental Physics Laboratory
 ! 455 Via Palou way, Stanford
 ! CA 94305, USA
 ! Email: shravan@stanford.edu 

 ! The driver routine. It utilizes the other subroutines to perform time-
 ! advancements. The time advancement is performed with the explicit Low
 ! Dissipation, Low Dispersion Runge Kutta 4th order method, described in
 ! Hu et al. (1996). Radial derivatives are computed using compact finite
 ! difference (Lele 1992). Horizontal variations are extracted by projecting
 ! the desired quantity onto Fourier harmonic space followed by an
 ! appropriate inverse fourier transform.
 

 !
 ! BASIC CHECKS:
 ! 1. If init is all right.
 ! 2. If the time-length of the simulation matches
 ! 3. If the directories are OK
 ! 4. If the forcing function is being read from the right directory
 ! 5. If the variables are initialized OK
 ! 6. If the flow profiles are allright   
 ! 7. If output and status are deleted or appropriately moved
 !
 
 ! ---------------------------------------------------------------------------
  use initialize
  use all_modules
  use mp_physics
  use mp_physics_2d
  use time_step
  use kernels

  implicit none

    integer i,j, init, ierr, t, k, index(1000,2), randinit(2), aind, nsteps_given,reclmax, loc
    real*8 start_time,end_time,tempf,t1,T00,nor1,nor2,nor3,nor4,nor5,e, mp_time
    real*8 total_time, start_mp_time, avg_time,zz, tempxy, rand1, rand2
    logical saved, iteration, init_variables 
    character*1 ci
    character*5 tempc
 20 format(8f12.4)


 ! Initializing the computation
 ! --------------------------------------------

 call Initialize_all
! do i=1,11
 ! call convert_to_string(i,contrib,2)
!  call adjoint_source_comp(211)
 !enddo
!stop
! if (COMPUTE_DATA) then
!  do k=1,nz
!   nor3 = (z(k) - 1.0)*695.9894
!   zz = exp((nor3 + 2.5)/0.9)  
!   zz = 12.*zz/(zz + 1.) + 4.
! !e =0.25 - 0.5/(1. + exp((nor3+2.5)/0.9))
!  e =0.12/(1.+exp((nor3 + 5.)/0.2)) + 0.25/(1. + exp(-(nor3+2.5)/0.9))
!  e = e/(1. + exp((nor3 + 0.2)/0.1))
!  print *,e,zz,nor3,z(k),k
  !if (nor3 > -0.25 .and. nor3 < 0.0) then
!   do i=1,nx
!    nor1 = (x(i)-0.5)*xlength*10.0**(-8.)
!    c2(i,1,k) = c2(i,1,k) * (1. + e * exp(-nor1**2./(2.*zz**2.)) + &
!	0.1*exp(-(nor3+3.)**2./(2.*1.5**2.) - (nor1-13.)**2./(2.*3.**2.)))
!   enddo
!  !endif
!  enddo
!  call writefits_3d('true_anomaly.fits',c2**0.5*dimc,nz)
!  !call writefits_3d('model_00.fits',c2**0.5*dimc,nz)
!! stop

! else
!!! stop
!  inquire(file=directory//'model.fits', exist = iteration)
!  if (iteration) then
!   call readfits(directory//'model.fits',c2,nz)
!   c2 = (c2/dimc)**2.
!  endif
! 
! endif

! v0_x = 5.0/dimc
! v0_z = 0.0

 if (CONSTRUCT_KERNELS) call PRODUCE_KERNELS 

 start_mp_time = MPI_WTIME()
 start_time = MPI_WTIME()

 !!!!!!! SOUND-SPEED PERTURBATIONS HERE !!!!!!

 if (.not. kernel_mode) then
  maxtime = floor(solartime*3600.0/timestep) + 1
  maxtime = 2*(maxtime/2) 
 endif

 if (.not. TEST_IN_2D) call mp_initialize_RHS
 if (TEST_IN_2D) call mp_initialize_RHS_2d 

 call Initialize_step
 
 end_time = MPI_WTIME()

 if (rank == 0) &
  print *, 'Initialization took: ', (end_time -start_time)

 !--------------------------------------------- 


 if ((time0 .ne.0) .and. (.not. kernel_mode)) &
   call read_in_initial_condition(directory, time0)


!c2 = (gradp0_x - boz*curlboy) !*diml*10.**(-8.) /( dimc**2. * dimrho) !
!c2 = (box**2. + boz**2.)**0.5/rho0**0.5 * dimc * 10.**(-5.)
!do k=1,nz

!c2(:,1,k) = (gradp0_z(:,1,k) + rho0(:,1,k)*g(k) + boz(:,1,k)*curlboy(:,1,k))!/(rho0(1,1,k)*g(k))
!c2(:,1,k) = g(k)*(gradp0_z(:,1,k)/c2(:,1,k) -gradrho0_z(:,1,k))/rho0(:,1,k) * (dimc/diml)**2.!+ rho0(:,1,k)*g(k) + boz(:,1,k)*curlboy(:,1,k))!/(rho0(1,1,k)*g(k))
!enddo
!c2 = (box**2. + boz**2.)**0.5/(rho0)**0.5 * 10.0
!c2(:,:,1) = 0.
!c2(:,:,nz) = 0.
!call writefits_3d('ca.fits',c2,nz)!p0*dimc**2.*dimrho,nz)
!stop
!do k=1,nz
! c2(:,1,k) = c2(:,1,k) * 1.01!(1. + 0.01*exp(-(x -0.5 + 10./300.)**2./(64.*(x(2)-x(1))**2.)))
!enddo
!c2rho0 = c2*rho0
!c_speed = c2**0.5
! c2 = c2*1.01 !* 1.01
! c2rho0 = c2rho0*1.01 !* 1.01
! c_speed = c2**0.5

!!!! OTHER INITIAL CONDITIONS!!!!!!
! a(:,:,:,1) - density
! a(:,:,:,2) - v_x
! and so on - v_y, v_z, p, b_x, b_y, b_z
! a(:,:,1,1) = orad_2d 
! call writefits_3d('orad_2d.fits',a(:,:,1,1),1)
!stop
 ! Length of simulation

  total_time = wall_time * 3600.

  steps = floor(outputcad/timestep )
  freq_filtering = floor(60./timestep)

  if (rank ==0) then
    print *,'FREQUENCY OF OUTPUT AT', steps, 'TIMESTEPS'
    print *, 'TIMESTEP =', timestep 
!    print *,'XXXXXX ------ NOT WRITING OUT XI VARIABLES -------- XXXXXXXX'
    print *,'OBSERVATION GRID POINT: ', o_rad, 'AND RADIUS:', z(o_Rad)
  endif

  ! Excitation series computation -----------------------------
  num_steps = FLOOR(DBLE(maxtime)*timestep/cadforcing) + 2
  num_steps = 2*(num_steps/2) + 2
  saved = .false.
  cadforcing_step = FLOOR(cadforcing/timestep)

!   if (kernel_mode) then  
!    maxtime = floor((final_time - timeline)/timestep) 
!    maxtime = cadforcing_step * floor(maxtime/cadforcing_step*1.) 
!    nsteps_given = floor(maxtime/cadforcing_step*1.) + 1
!   endif

  if (time0 == 0) init = 0
  if (time0 .NE. 0) init = time0+1

  if (rank ==0)  &
    print *, 'NUMBER OF TIMESTEPS =', maxtime


IF (kernel_mode) then
   call DETERMINE_STATUS(init, maxtime)
   if (rank==0) call system('rm Instruction')

  if (compute_forward) then

   indexglob = 0
   if (rank == 0) then 
     open(19,file=directory//'forward'//contrib//'/rms_hist',status='replace',form='formatted',position='append') 
     open(745, file=directory//'forward'//contrib//'/kernel_info',status='unknown', action='write')
     write(745,*) init
     write(745,*) maxtime
     close(745)
   endif

  delta_width = 2.0d0/(z(e_rad+1) - z(e_rad-1))

  if (magnetic) then
   do k=1,dim2(rank)
    do j=1,nx  
      delta_width_2d(j,k) = 2.0d0/(z(erad_2d(j,k)+1) - z(erad_2d(j,k)-1))
    enddo
   enddo
  endif

   call timestepping(init)

   if (rank==0) then 
    close(28)
    close(1244)
   endif

  call read_binary_writefits(directory//'forward'//contrib//'/vz_cc.bin',&
		indexglob, directory//'forward'//contrib//'/vz_cc.fits')

  call adjoint_source_comp(indexglob)

 endif

 if (compute_adjoint) then
   if (rank == 0) then 
     open(19,file=directory//'adjoint'//contrib//'/rms_hist',status='replace',form='formatted',position='append') 
     open(745, file=directory//'adjoint'//contrib//'/kernel_info',status='unknown', action='write')
     write(745,*) floor(maxtime/cadforcing_step*1.) + 1
   endif

   call timestepping(init)

   if (rank==0) then 
    write(745,*) floor(time/cadforcing_step*1.)*cadforcing_step
    close(745)
    close(28)
    close(1244)

   endif


  open(44,file=directory//'adjoint'//contrib//'/kernel_info',form='formatted',status='unknown')
  read(44,*) init!nt_kern
  read(44,*) maxtime
  close(44)
 endif


ENDIF  ! CORRESPONDING TO "IF KERNEL_MODE"



 call MPI_BARRIER(MPI_COMM_WORLD, ierr)

 if (.not. kernel_mode) call write_out_full_state(directory) 
! if (kernel_mode) then
!!  if (compute_adjoint) call write_out_full_state(directory//'adjoint'//contrib//'/')
!  if (compute_forward) call write_out_full_state(directory//'forward'//contrib//'/')
! endif

 call dfftw_destroy_plan(fftw_plan_fwd_x)
 call dfftw_destroy_plan(fftw_plan_inv_x)

 call dfftw_destroy_plan(fftw_plan_fwd_y)
 call dfftw_destroy_plan(fftw_plan_inv_y)

 if (rank == 0) then
   close(19)
 endif

 call MPI_BARRIER(MPI_COMM_WORLD, ierr) 
 call MPI_FINALIZE(ierr)


end program driver

! --------------------------------------------------------------------------------------------



SUBROUTINE TIMESTEPPING(init)
 
 use time_step
 use mp_physics
 use mp_physics_2d
 use initialize
 use all_modules

 implicit none

 integer i,j,k, init, nsteps_given
 real*8 start_time, end_time, T00, t1, nor1
 real*8 tempstore(nx, dim2(rank))

  do i = init,maxtime
    time = i
  ! code to be timed

    call cpu_time(t1)
    T00 = t1
    start_time = MPI_WTIME() 
    if (.not. kernel_mode) call lagrange_interp

    if ( kernel_mode .and. compute_adjoint .and. &
	(i .le. forcing_length*cadforcing_step) ) call lagrange_interp

    call step()
    end_time = MPI_WTIME() 
    call cpu_time(t1) 
   
    if (.not. displ) nor1 = norm(a(:,:,:,4)) !* dimc
    if (displ .and. test_in_2d) nor1 = sum(abs(a(:,:,:,6)**2.))**0.5
  
 
  ! Timing and 'convergence' statistics 
  ! The criterion used to see if the results are 'convergent' is the
  ! L2 norm of the radial velocity. 

   if (rank ==0) then
 
    print *, '----------------------------------------------------------'
    print *, '	Iteration #,  Cpu Time and Real Time:'
    print *,       time,       	(t1 - T00),        (time+1)*timestep
    print *
    print *,'	Wall clock time', (end_time - start_time)
    print *
    print *,nor1
    print *, '----------------------------------------------------------'
    print *,'' 

    write(19,*) nor1

   endif
 

   timeline = timeline + timestep

   if (mod(time,steps) == 0) then
    
     write(28, *) time, timeline
     flush(28)

       if (compute_adjoint) tempstore = a(:,:,o_rad,6) !* rho0(:,:,o_rad)**0.5
       if (compute_forward .and. (.not. magnetic)) tempstore = a(:,:,o_rad,6)
      
            
      if (compute_forward .and. magnetic) then
        do k=1,dim2(rank) 
         do j=1,nx
 	  tempstore(j,k) = a(j,k,orad_2d(j,k),6) 
  	 enddo
        enddo
      endif

      indexglob = indexglob +1 
      if ( compute_forward) then
        call write_binary_seq_2D(1244, tempstore, 1, indexglob)      
        flush(1244)
      endif
       if (rank == 0) print *, 'Writing output'

       if (compute_forward) then
         if (.not. LINESEARCH .and. (.NOT. COMPUTE_DATA)) call write_out_partial_state(directory//'forward'//contrib//'/')
	 !call write_out_slice(directory//'forward'//contrib//'/')
       endif
       if (compute_adjoint) call write_out_partial_state(directory//'adjoint'//contrib//'/') 

   endif


   if (mod(time+1,200) == 0)  then

     if (.not. kernel_mode) call write_out_full_state(directory) 
     if (kernel_mode) then

       if (rank==0) open(99, file=directory//'unfinished_calc_'//contrib,action='write', status='replace')

       if (compute_adjoint) then
          call write_out_full_state(directory//'adjoint'//contrib//'/') 

          if (rank==0) then 
            if (.not. generate_wavefield) write(99,"(A)") 'adjoint'
            if (generate_wavefield) write(99,"(A)") 'spectral'
            write(99,*) time
            write(99,"(A)") contrib
            write(99,*) indexglob
            close(99)
          endif
                 
       else
          call write_out_full_state(directory//'forward'//contrib//'/')
          if (rank==0) then 
            write(99,"(A)") 'forward'
            write(99,*) time
            write(99,"(A)") contrib
            write(99,*) indexglob
            close(99)
          endif
       endif

     endif


     if (rank ==0 ) print *,'Full output done, index:',time

   endif


 end do

 if (rank==0) then

  call system('rm '//directory//'unfinished_calc_'//contrib)

  if (compute_forward) then
   open(223, file=directory//'status/'//'forward'//contrib,status='unknown')
   close(223)
   !call system('rm -rf '//directory//'forward'//contrib//'/*full*')
  endif

  if (compute_adjoint) then
   open(223, file=directory//'status/'//'adjoint'//contrib,status='unknown')
   close(223)
   !call system('rm -rf '//directory//'adjoint'//contrib//'/*full*')
  endif

 endif

END SUBROUTINE TIMESTEPPING


! --------------------------------------------------------------------------------------------

SUBROUTINE DETERMINE_STATUS(init, nsteps_given)

 use time_step
 use mp_physics
 use mp_physics_2d
 use initialize
 use all_modules

 implicit none

 integer i,j,k, init, nsteps_given, indexnum
 real*8 start_time, end_time, T00, t1, nor1, reclmax, xloc
 character*80 calctype
 character*7 keyword
 character*1 ci
 integer inde
 logical lexist

 reclmax = nx * ny * 2
 init =0
  
!  ci = '1' 
 ! if (contrib == '1') ci = '2'

  if (generate_wavefield) then
    COMPUTE_ADJOINT = .FALSE.
    COMPUTE_FORWARD = .TRUE.
    generate_wavefield = .true.
    keyword = 'forward'
  endif

  if (compute_adjoint) then
    COMPUTE_ADJOINT = .TRUE.
    COMPUTE_FORWARD = .FALSE.
    generate_wavefield = .false.
    keyword = 'adjoint'
  endif
  initcond = 0

  if (compute_forward) then
   allocate(vr(nx,dim2(rank),1), fwdsource(nsteps_given+8))
   vr = 0.0 
   call forward_source(nsteps_given + 8)



   !if (rank==0) print *,'Need a file of: ', forcing_length, ' temporal slices.'

!   call read_binary_reverse(directory//'forward'//contrib//'/vz_2D_data.bin', &
!	vr(:,:,1:nsteps_given),nsteps_given)

   !forcing_length = nsteps_given
!    call readfits(directory//'forward'//contrib//'/source.fits',& !'_'//contrib//
!	vr(:,:,1:forcing_length),forcing_length)

!   if (contrib=='1') then
     read(contrib,*) indexnum

     open(356,file=directory//'master.pixels',action='read', position='rewind')
     do i=1,nmasters
      if (indexnum .ne. i) read(356,*) 
      if (indexnum .eq. i) read(356,*) xloc
     enddo
     close(356)

    dx = (x(2) - x(1)) * xlength*10.**(-8)
    do j=1,dim2(rank)
      vr(:,j,1) = exp(-((x-0.5)*xlength*10.**(-8) - xloc)**2./(2.*dx**2.))
    enddo

    dx = dx/(xlength*10.**(-8))
!    timeline = timeline + (nsteps_given -1)*cadforcing_step*floor(timestep)
!    timeline = -timeline

!    maxtime = floor((abs(timeline))*2./timestep) 
!    nsteps_given = floor(maxtime/cadforcing_step*1.) + 1

    if (rank==0) open(28, file=directory//'forward'//contrib//'/timeline',&
		status='replace', action='write',position='append')


    reclmax = nx * ny * 2
    open(1244,file=directory//'forward'//contrib//'/vz_cc.bin',&
        form='unformatted',status='replace', action ='write',&
	access='direct',recordtype='fixed',recl=reclmax)


    delta_width = 2.0d0/(z(e_rad+1) - z(e_rad-1))

    if (magnetic) then
     do k=1,dim2(rank)
      do j=1,nx  
       delta_width_2d(j,k) = 2.0/(z(erad_2d(j,k)+1) - z(erad_2d(j,k)-1))
      enddo
     enddo
    endif


   endif

   if (compute_adjoint) then

     !forcing_length = floor(maxtime/cadforcing_step*1.) + 1
!    print *,forcing_length,maxtime, cadforcing_step
    allocate(vr(nx,dim2(rank),forcing_length + 8))

    vr = 0.0 

    if (rank==0) print *,'Need a file of: ', forcing_length, ' temporal slices.'
    call readfits(directory//keyword//contrib//'/source.fits',& !'_'//contrib//
	vr(:,:,1:forcing_length),forcing_length)

 
    if (rank==0 .and. (init .gt. 1)) open(28, file=directory//keyword//contrib//'/timeline',&
		status='old', action='write', position='append')

    if (rank==0 .and. (init .eq. 1)) open(28, file=directory//keyword//contrib//'/timeline',&
		status='replace', action='write', position='rewind')

    delta_width = 2.0d0/(z(o_rad+1) - z(o_rad-1))

    if (magnetic) then
     do k=1,dim2(rank)
      do j=1,nx  
       delta_width_2d(j,k) = 2.0/(z(orad_2d(j,k)+1) - z(orad_2d(j,k)-1))
      enddo
     enddo
    endif


   endif
 indexglob = 0
 time_old = 0
 init=0
END SUBROUTINE DETERMINE_STATUS

!================================================================================

 SUBROUTINE convert_binary_to_fits(input_filename, nxs, nys, dim3, output_filename)
  use all_modules
  implicit none
  integer nxs, nys, dim3, i, j, k, reclmax
  real*8 array(nxs, nys, dim3)
  character*(*) input_filename, output_filename

  reclmax = nxs * nys * 2

  open(44,file=input_filename,recordtype='fixed', &
        form='unformatted', action='read', access='direct',status='old', &
        recl = reclmax)

   do k=1,dim3
    read(44, rec=k) ((array(i,j,k),i=1,nxs),j=1,nys)
   enddo
   close(44)

   call writefits_local(output_filename, array, nxs, nys, dim3)

 END SUBROUTINE CONVERT_BINARY_TO_FITS

!================================================================================================


SUBROUTINE ADJOINT_SOURCE_COMP(nt)
 use initialize
 use all_modules
 implicit none
 integer *8 fwdplantemp, invplantemp, invplantemp2, fwdplandata, invplandata
 logical lexist
 integer i, nt, indexnum, pord, timesmax, halftime!, jj
 integer loc, leng, lef, rig, j, nmeasurements, nmeas, ierr
 character*1 ord
 real*8 mindist, maxdist, window, distances(nx), x00, t(nt), dat(nx, dim2(rank),nt), tau
 real*8 pcoef(5,4), adj(nx,dim2(rank),nt), acc(nx,dim2(rank),nt), windows(nt), filt(nt),dt
 real*8 freqnu(nt), leftcorner, rightcorner, ccdot(nx,dim2(rank),nt), dnu, con, misfit,misfit_tot
 complex*16, dimension(nx/2+1,dim2(rank),nt) :: filtout, temp, filtdat

 filtout = 1.0
 pcoef = 0.0
 misfit = 0.0 
 read(contrib,*) indexnum

 if (rank==0) print *,'Doing Adjoint Source'

 open(356,file=directory//'master.pixels',action='read', position='rewind')
 do i=1,nmasters
  if (indexnum .ne. i) read(356,*) 
  if (indexnum .eq. i) read(356,*) x00
 enddo
 close(356)

 distances = abs((x-0.5)*xlength*10.0**(-8.) - x00)

 open(44,file=directory//'forward'//contrib//'/timeline',action='read')
 do i=1,nt
  read(44,*) loc, t(nt-i+1)
 enddo
 close(44)  
 t = -t/60.0
 dt = t(2)-t(1)

 call distmat(nt,1,freqnu)
 freqnu = freqnu/(nt*dt*60.)
 dnu = freqnu(2) - freqnu(1)
 
 call dfftw_plan_dft_r2c_3d(fwdplantemp, nx, dim2(rank), nt, acc, filtout, FFTW_ESTIMATE)
 call dfftw_plan_dft_c2r_3d(invplantemp, nx, dim2(rank), nt, filtout, acc, FFTW_ESTIMATE)
 call dfftw_plan_dft_c2r_3d(invplantemp2, nx, dim2(rank), nt, filtout, ccdot, FFTW_ESTIMATE)

 call dfftw_plan_dft_r2c_3d(fwdplandata, nx, dim2(rank), nt, dat, filtdat, FFTW_ESTIMATE)
 call dfftw_plan_dft_c2r_3d(invplandata, nx, dim2(rank), nt, filtdat, dat, FFTW_ESTIMATE)


 call readfits(directory//'forward'//contrib//'/vz_cc.fits', acc, nt)
! call readfits(directory//'data/'//contrib//'.fits', dat, nt)

 call dfftw_execute(fwdplantemp)
 call dfftw_execute(fwdplandata)

 call dfftw_destroy_plan(fwdplantemp)
 call dfftw_destroy_plan(fwdplandata)

 inquire(file=directory//'filter.params', exist=lexist)
 if (lexist) then

  open(94,file=directory//'filter.params', action='read',position='rewind')
   read(94,*) leftcorner ! mHz
   read(94,*) rightcorner  ! mHz
  close(94)

  filt = 1./(1. + exp((leftcorner - abs(freqnu))*2./dnu)) &
       - 1./(1. + exp((rightcorner - abs(freqnu))*2./dnu))

  do i=1,nt
   filtout(:,1,i) = filtout(:,1,i) * filt(i)
   filtdat(:,1,i) = filtdat(:,1,i) * filt(i)
  enddo
  temp = filtout
  call dfftw_execute(invplantemp)
  call dfftw_execute(invplandata)
  filtout = temp

  call dfftw_destroy_plan(invplantemp)
  call dfftw_destroy_plan(invplandata)

 endif
 con = 2.0*pi/(dble(nt)*dble(nx))
 do i=1,nt
  filtout(:,1,i) = filtout(:,1,i) * eye * freqnu(i) * con
 enddo

 call dfftw_execute(invplantemp2)
 call dfftw_destroy_plan(invplantemp2)


 pcoef = 0.0
 pcoef(1,1) = 11.5414
 pcoef(2,1) = 0.742286
 pcoef(3,1) = -0.00599439
 pcoef(4,1) = 2.59565e-5
 pcoef(5,1) = -2.69437e-8

 pcoef(1,2) = 16.9286
 pcoef(2,2) = 0.840476
 pcoef(3,2) = -0.00238095

 pcoef(1,3) = 25.0
 pcoef(2,3) = 0.9
 pcoef(3,3) = -0.002222

 call readfits(directory//'quiet.fits', acc, nt)
 call readfits(directory//'flow.fits', dat, nt)

 adj = 0.0
 nmeasurements = 0
 do pord=1,4
  call convert_to_string(pord, ord, 1)

  inquire(file=directory//'params.p'//ord, exist=lexist) 
  if (lexist) then
    open(97, file = directory//'params.p'//ord, action = 'read', status='old')
    read(97,*) mindist
    read(97,*) maxdist
    read(97,*) window
    close(97)
!    open(993, file=directory//'forward'//contrib//'/times.p'//ord,status='unknown',position='rewind',action='write')
!    open(45, file='temp',status='unknown',position='rewind',action='write')
    halftime = nint(window/(2.*dt))
    leng = 2*halftime+1
    do i=1,nx

     if ((distances(i) > mindist) .and. (distances(i) < maxdist)) then
!     if (i==nx/2+15) then

      timesmax = nint((pcoef(1,pord) - t(1) + pcoef(2,pord)*distances(i) + pcoef(3,pord)*distances(i)**2. &
	 + pcoef(4,pord)*distances(i)**3. + pcoef(5,pord)*distances(i)**4.)/dt) + 1

      loc = maxloc(acc(i,1,(timesmax-6):(timesmax+6)),1) + timesmax - 7
      lef = loc - halftime
      rig = loc + halftime
      call compute_tt(acc(i,1,lef:rig),dat(i,1,lef:rig),tau,dt,leng)
!        open(45,file='temp',action='write')
!       do jj=1,nt !lef,rig
!       write(45,*) acc(i,1,jj), dat(i,1,jj)
!       enddo
!       close(45)
!	print *,lef, rig,tau, i
! 	stop
!      write(993, *) (x(i)-0.5)*xlength*10.**(-8.),tau
      misfit = misfit + tau**2.
      nmeasurements = nmeasurements + 1
      tau = 1.0 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!HERERE WWWW
      con = -tau/(sum(ccdot(i,1,lef:rig)**2.)*dt)

      print *,sum(con*ccdot(i,1,lef:rig)*(dat(i,1,lef:rig) - acc(i,1,lef:rig)))*dt,'=TT DIFF'

      do j=lef,rig
       adj(i,1,nt-j+1) = ccdot(i,1,j) * con + adj(i,1,nt-j+1)
      enddo
     endif
    enddo
 !   close(993)
  endif
 enddo
! adj = adj * nx
 inquire(directory=directory//'adjoint'//contrib, exist=lexist)
 if (.not. lexist .and. (rank==0)) &
	call system('mkdir '//directory//'adjoint'//contrib)

 call writefits_3d(directory//'adjoint'//contrib//'/source.fits',adj,nt)

 call MPI_REDUCE(misfit, misfit_tot, 1, MPI_DOUBLE_PRECISION, &
				MPI_SUM, 0, MPI_COMM_WORLD, ierr)

 call MPI_REDUCE(nmeasurements, nmeas, 1, MPI_INTEGER, &
				MPI_SUM, 0, MPI_COMM_WORLD, ierr)

! misfit_tot = misfit_tot * 0.5
! print *,'Total',misfit_tot, misfit_tot/nmeasurements
 if (rank==0) then
   open(34, file=directory//'adjoint'//contrib//'/information', &
		action='write', position='rewind', status='replace')
   write(34,*) nt
   write(34,*) t(1)*60.
   write(34,*) t(nt)*60.
   if (test_in_2d) write(34,*) '2D'
   close(34)
   
   inquire(file=directory//'kernel/misfit',exist=lexist)
   if (.not. lexist) &
   open(88, file=directory//'kernel/misfit', status='unknown',&
		action='write')

   if (lexist) &
   open(88, file=directory//'kernel/misfit', status='old',&
		action='write', position='append')

   misfit_tot = misfit_tot * 0.5
   write(88,*) indexnum, nmeas, misfit_tot
   close(88)

 endif

END SUBROUTINE ADJOINT_SOURCE_COMP

!================================================================================

SUBROUTINE COMPUTE_TT(u0, u, tau, dt, nt)

 implicit none
 integer nt,l1,l2, loc, i
 real*8 u0(nt), u(nt), cc(-nt+1:nt-1), t(-nt+1:nt-1)
 real*8 times(3), tau, invmat(3,3), p1, p2, dt, mat(3,3)

 do i=(-nt+1),0
  t(i) = i*dt
  l1 = -i + 1
  l2 = nt
  cc(i) = sum(u0(l1:l2)*u(1:(l2-l1+1))) * dt
 enddo

 do i=1,nt-1
  t(i) = i*dt
  l1 = i + 1
  l2 = nt-i
  cc(i) = sum(u0(1:l2)*u(l1:nt)) * dt
 enddo

 loc = maxloc(cc,1)-nt

 times(1) = t(loc-1)
 times(2) = t(loc)
 times(3) = t(loc+1)

 mat(:,1) = 1
 mat(:,2) = times
 mat(:,3) = times**2.

 call inverse(mat, invmat, 3)

 p1 = invmat(2,1)*cc(loc-1) + invmat(2,2) * cc(loc) + invmat(2,3) * cc(loc+1)
 p2 = invmat(3,1)*cc(loc-1) + invmat(3,2) * cc(loc) + invmat(3,3) * cc(loc+1)

 tau = -p1*0.5/p2

 !print *,loc,-nt+1, nt-1, tau
END SUBROUTINE COMPUTE_TT

!================================================================================

  subroutine inverse(a,c,n)
!============================================================
! Inverse matrix
! Method: Based on Doolittle LU factorization for Ax=b
! Alex G. December 2009
!-----------------------------------------------------------
! input ...
! a(n,n) - array of coefficients for matrix A
! n      - dimension
! output ...
! c(n,n) - inverse matrix of A
! comments ...
! the original matrix a(n,n) will be destroyed 
! during the calculation
!===========================================================
implicit none 
integer n
double precision a(n,n), c(n,n)
double precision L(n,n), U(n,n), b(n), d(n), x(n)
double precision coeff
integer i, j, k

! step 0: initialization for matrices L and U and b
! Fortran 90/95 aloows such operations on matrices
L=0.0
U=0.0
b=0.0

! step 1: forward elimination
do k=1, n-1
   do i=k+1,n
      coeff=a(i,k)/a(k,k)
      L(i,k) = coeff
      do j=k+1,n
         a(i,j) = a(i,j)-coeff*a(k,j)
      end do
   end do
end do

! Step 2: prepare L and U matrices 
! L matrix is a matrix of the elimination coefficient
! + the diagonal elements are 1.0
do i=1,n
  L(i,i) = 1.0
end do
! U matrix is the upper triangular part of A
do j=1,n
  do i=1,j
    U(i,j) = a(i,j)
  end do
end do

! Step 3: compute columns of the inverse matrix C
do k=1,n
  b(k)=1.0
  d(1) = b(1)
! Step 3a: Solve Ld=b using the forward substitution
  do i=2,n
    d(i)=b(i)
    do j=1,i-1
      d(i) = d(i) - L(i,j)*d(j)
    end do
  end do
! Step 3b: Solve Ux=d using the back substitution
  x(n)=d(n)/U(n,n)
  do i = n-1,1,-1
    x(i) = d(i)
    do j=n,i+1,-1
      x(i)=x(i)-U(i,j)*x(j)
    end do
    x(i) = x(i)/u(i,i)
  end do
! Step 3c: fill the solutions x(n) into column k of C
  do i=1,n
    c(i,k) = x(i)
  end do
  b(k)=0.0
end do
end subroutine inverse

!================================================================================

SUBROUTINE FORWARD_SOURCE(nt)
 use initialize
 use all_modules
 implicit none
 integer nt, k
 real*8 ttime


 do k=1,nt
  ttime = (k- 1.)*timestep
  fwdsource(k) = cos(2.*pi*nupeak*ttime) * exp(-2.*(ttime-2./nupeak)**2.*nuwidth**2.)
  if (abs(fwdsource(k)) .lt. 0.00001) fwdsource(k) = 0.0
 enddo
 fwdsource = fwdsource * 10.**(-5.)


END SUBROUTINE FORWARD_SOURCE
!================================================================================
