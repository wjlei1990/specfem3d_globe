!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  7 . 0
!          --------------------------------------------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!                and CNRS / University of Marseille, France
!                 (there are currently many more authors!)
! (c) Princeton University and CNRS / University of Marseille, April 2014
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along
! with this program; if not, write to the Free Software Foundation, Inc.,
! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
!
!=====================================================================
!
! United States and French Government Sponsorship Acknowledged.

! we switch between vectorized and non-vectorized version by using pre-processor flag FORCE_VECTORIZATION
! and macros INDEX_IJK, DO_LOOP_IJK, ENDDO_LOOP_IJK defined in config.fh
#include "config.fh"

  subroutine iterate_time()

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore
  use specfem_par_outercore
  use specfem_par_movie

  implicit none

  ! timing
  double precision, external :: wtime
  
  ! Track maximum norm of displacement, velocity in crust/mantle, 
  ! outer core, and inner core
  ! Only works in the CPU version of the code
  real(kind=CUSTOM_REAL), dimension(NGLOB_CRUST_MANTLE) :: ndispl_max_cm,nveloc_max_cm
  real(kind=CUSTOM_REAL), dimension(NGLOB_INNER_CORE) :: ndispl_max_ic,nveloc_max_ic
  real(kind=CUSTOM_REAL), dimension(NGLOB_OUTER_CORE) :: ndispl_max_oc,nveloc_max_oc
  
  ! Track maximum norm of strain and stress in the crust/mantle at the element level
  ! Only works in the CPU version of the code
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE) :: nstrain_max_cm,nstress_max_cm
  ! Global level
  real(kind=CUSTOM_REAL), dimension(NGLOB_CRUST_MANTLE) :: nstrain_max_cm_glob,nstress_max_cm_glob
  
  ! Do we want to track the peak norm values?
  logical :: GET_PEAK_NORMS
  

  ! for EXACT_UNDOING_TO_DISK
  integer :: ispec,iglob,i,j,k
  
  
  ! Initialize variables
  ndispl_max_cm(:)=0
  nveloc_max_cm(:)=0
  ndispl_max_ic(:)=0
  nveloc_max_ic(:)=0
  ndispl_max_oc(:)=0
  nveloc_max_oc(:)=0
  nstress_max_cm(:,:,:,:)=0
  nstrain_max_cm(:,:,:,:)=0
  GET_PEAK_NORMS = .true. 
  

  !----  create a Gnuplot script to display the energy curve in log scale
  if (OUTPUT_ENERGY .and. myrank == 0) then
    open(unit=IOUT_ENERGY,file=trim(OUTPUT_FILES)//'plot_energy.gnu',status='unknown',action='write')
    write(IOUT_ENERGY,*) 'set term wxt'
    write(IOUT_ENERGY,*) '#set term postscript landscape color solid "Helvetica" 22'
    write(IOUT_ENERGY,*) '#set output "energy.ps"'
    write(IOUT_ENERGY,*) 'set logscale y'
    write(IOUT_ENERGY,*) 'set xlabel "Time step number"'
    write(IOUT_ENERGY,*) 'set ylabel "Energy (J)"'
    write(IOUT_ENERGY,'(a152)') '#plot "energy.dat" us 1:2 t "Kinetic Energy" w l lc 1, "energy.dat" us 1:3 &
                         &t "Potential Energy" w l lc 2, "energy.dat" us 1:4 t "Total Energy" w l lc 4'
    write(IOUT_ENERGY,*) '#pause -1 "Hit any key..."'
    write(IOUT_ENERGY,*) '#plot "energy.dat" us 1:2 t "Kinetic Energy" w l lc 1'
    write(IOUT_ENERGY,*) '#pause -1 "Hit any key..."'
    write(IOUT_ENERGY,*) '#plot "energy.dat" us 1:3 t "Potential Energy" w l lc 2'
    write(IOUT_ENERGY,*) '#pause -1 "Hit any key..."'
    write(IOUT_ENERGY,*) 'plot "energy.dat" us 1:4 t "Total Energy" w l lc 4'
    write(IOUT_ENERGY,*) 'pause -1 "Hit any key..."'
    close(IOUT_ENERGY)
  endif

  ! open the file in which we will store the energy curve
  if (OUTPUT_ENERGY .and. myrank == 0) &
    open(unit=IOUT_ENERGY,file=trim(OUTPUT_FILES)//'energy.dat',status='unknown',action='write')

!
!   s t a r t   t i m e   i t e r a t i o n s
!

  ! synchronize all processes to make sure everybody is ready to start time loop
  call synchronize_all()
  if (myrank == 0) then
    write(IMAIN,*) 'All processes are synchronized before time loop'
    write(IMAIN,*)
    write(IMAIN,*) 'Starting time iteration loop...'
    write(IMAIN,*)
    call flush_IMAIN()
  endif
  call synchronize_all()

  ! create an empty file to monitor the start of the simulation
  if (myrank == 0) then
    open(unit=IOUT,file=trim(OUTPUT_FILES)//'/starttimeloop.txt',status='unknown',action='write')
    write(IOUT,*) 'hello, starting time loop'
    close(IOUT)
  endif

  ! initialize variables for writing seismograms
  seismo_offset = it_begin-1
  seismo_current = 0

  ! get MPI starting time
  time_start = wtime()

  ! *********************************************************
  ! ************* MAIN LOOP OVER THE TIME STEPS *************
  ! *********************************************************

  if (EXACT_UNDOING_TO_DISK) call setup_exact_undoing_to_disk()

  ! time loop
  do it = it_begin,it_end

    ! simulation status output and stability check
    if (mod(it,NTSTEP_BETWEEN_OUTPUT_INFO) == 0 .or. it == it_begin + 4 .or. it == it_end) then
      call check_stability()
      if (I_am_running_on_a_slow_node) goto 100
    endif

    do istage = 1, NSTAGE_TIME_SCHEME ! is equal to 1 if Newmark because only one stage then

      if (USE_LDDRK) then
        ! update displacement using Runge-Kutta time scheme
        call update_displ_lddrk()
      else
        ! update displacement using Newmark time scheme
        call update_displ_Newmark()
      endif

      ! acoustic solver for outer core
      ! (needs to be done first, before elastic one)
      call compute_forces_acoustic()

      ! elastic solver for crust/mantle and inner core
      call compute_forces_viscoelastic()

    enddo ! end of very big external loop on istage for all the stages of the LDDRK time scheme (only one stage if Newmark)

    ! save the forward run to disk for the alpha kernel only
    if (EXACT_UNDOING_TO_DISK .and. SIMULATION_TYPE == 1) then
      do ispec = 1, NSPEC_CRUST_MANTLE
        do k = 1, NGLLZ
          do j = 1, NGLLY
            do i = 1, NGLLX
              iglob = ibool_crust_mantle(i,j,k,ispec)
              if (integer_mask_ibool_exact_undo(iglob) /= -1) &
                buffer_for_disk(integer_mask_ibool_exact_undo(iglob),it_of_this_exact_subset) = &
                  eps_trace_over_3_crust_mantle(i,j,k,ispec)
            enddo
          enddo
        enddo
      enddo
      if (it_of_this_exact_subset == it_exact_subset_end) then
        do it_of_this_exact_subset = 1, it_exact_subset_end
          write(IFILE_FOR_EXACT_UNDOING,rec=it_exact_subset_offset+it_of_this_exact_subset) &
            buffer_for_disk(:,it_of_this_exact_subset)
        enddo
        it_of_this_exact_subset = 1
        it_exact_subset_offset = it_exact_subset_offset + it_exact_subset_end
        it_exact_subset_end = min(NSTEP_FOR_EXACT_UNDOING, it_end - it_exact_subset_offset)
      else
        it_of_this_exact_subset = it_of_this_exact_subset + 1
      endif
    endif

    ! kernel simulations (forward and adjoint wavefields)
    if (SIMULATION_TYPE == 3) then

      if (.not. EXACT_UNDOING_TO_DISK) then
        ! note: we step back in time (using time steps - DT ), i.e. wavefields b_displ_..() are time-reversed here

        ! reconstructs forward wavefields based on last stored wavefield data

        ! note: NSTAGE_TIME_SCHEME is equal to 1 if Newmark because only one stage then
        do istage = 1, NSTAGE_TIME_SCHEME

          if (USE_LDDRK) then
            ! update displacement using Runge-Kutta time scheme
            call update_displ_lddrk_backward()
          else
            ! update displacement using Newmark time scheme
            call update_displ_Newmark_backward()
          endif

          ! acoustic solver for outer core
          ! (needs to be done first, before elastic one)
          call compute_forces_acoustic_backward()

          ! elastic solver for crust/mantle and inner core
          call compute_forces_viscoelastic_backward()

        enddo

        ! restores last time snapshot saved for backward/reconstruction of wavefields
        ! note: this is done here after the Newmark time scheme, otherwise the indexing for sources
        !          and adjoint sources will become more complicated
        !          that is, index it for adjoint sources will match index NSTEP - 1 for backward/reconstructed wavefields
        if (it == 1) then
          call read_forward_arrays()
        endif

      else ! of if (.not. EXACT_UNDOING_TO_DISK)

        ! read the forward run from disk for the alpha kernel only
        it_of_this_exact_subset = it_of_this_exact_subset + 1
        if (it_of_this_exact_subset > it_exact_subset_end) then
          it_exact_subset_offset = it_exact_subset_offset + it_exact_subset_end
          it_exact_subset_end = min(NSTEP_FOR_EXACT_UNDOING, it_end - it_exact_subset_offset)
          do it_of_this_exact_subset = 1, it_exact_subset_end
            ! here we time revert the forward run by reading time step NSTEP - it + 1
            ! but here, it == it_exact_subset_offset + it_of_this_exact_subset
            read(IFILE_FOR_EXACT_UNDOING,rec=NSTEP-it_exact_subset_offset-it_of_this_exact_subset+1) &
              buffer_for_disk(:,it_of_this_exact_subset)
          enddo
          it_of_this_exact_subset = 1
        endif

        do ispec = 1, NSPEC_CRUST_MANTLE
          do k = 1, NGLLZ
            do j = 1, NGLLY
              do i = 1, NGLLX
                iglob = ibool_crust_mantle(i,j,k,ispec)
                if (integer_mask_ibool_exact_undo(iglob) /= -1) then
                  b_eps_trace_over_3_crust_mantle(i,j,k,ispec) = &
                    buffer_for_disk(integer_mask_ibool_exact_undo(iglob),it_of_this_exact_subset)
                endif
              enddo
            enddo
          enddo
        enddo

      endif ! of if (.not. EXACT_UNDOING_TO_DISK)

    endif ! kernel simulations

    ! calculating gravity field at current timestep
    if (GRAVITY_SIMULATION) call gravity_timeseries()

    ! write the seismograms with time shift (GPU_MODE transfer included)
    call write_seismograms()

    ! adjoint simulations: kernels
    ! attention: for GPU_MODE and ANISOTROPIC_KL it is necessary to use resort_array (see lines 265-268)
    if (SIMULATION_TYPE == 3) then
      call compute_kernels()
    endif

    ! outputs movie files
    if (MOVIE_SURFACE .or. MOVIE_VOLUME) call write_movie_output()

    ! first step of noise tomography, i.e., save a surface movie at every time step
    ! modified from the subroutine 'write_movie_surface'
    if (NOISE_TOMOGRAPHY == 1) then
      call noise_save_surface_movie()
    endif

    ! updates VTK window
    if (VTK_MODE) then
      call it_update_vtkwindow()
    endif
	
    ! Call subroutine to compute the maximum norms of displacement and velocity
    ! Only for CPU version
    if (.not. GPU_MODE) then
      if (GET_PEAK_NORMS) then
        if (SIMULATION_TYPE == 1) then
           call max_norms_dispvel(ndispl_max_cm,nveloc_max_cm,ndispl_max_ic,nveloc_max_ic, &
                ndispl_max_oc,nveloc_max_oc)
		  
           ! Compute the strain (global) in the crust and mantle from displacement
           ! and store the peak norms
           if (COMPUTE_AND_STORE_STRAIN) then
              call max_norms_strain_cm(nstrain_max_cm)
           endif
           call max_norms_stress_cm(nstress_max_cm)
        endif
      endif
    endif

  !
  !---- end of time iteration loop
  !
  enddo   ! end of main time loop

 100 continue


  if (SIMULATION_TYPE == 3 .and. GPU_MODE) then
    ! attention: cijkl_kl_crust_mantle is sorted differently on GPU and CPU
    call resort_array(Mesh_pointer)
  endif

  ! close the huge file that contains a dump of all the time steps to disk
  if (EXACT_UNDOING_TO_DISK) call finish_exact_undoing_to_disk()

  ! user output of runtime
  call print_elapsed_time()

  ! Transfer fields from GPU card to host for further analysis
  if (GPU_MODE) call it_transfer_from_GPU()

!----  close energy file
  if (OUTPUT_ENERGY .and. myrank == 0) close(IOUT_ENERGY)
  
  
  ! Write the peak values for norm displacement, velocity to .bin files, which can be
  ! converted to VTK using xcombine_vol_data_vtk
  ! Only for CPU version
  if (.not. GPU_MODE) then
    if (GET_PEAK_NORMS) then
      if (SIMULATION_TYPE == 1) then
        call write_bin_ndispvel(ndispl_max_cm,nveloc_max_cm,ndispl_max_ic,nveloc_max_ic, &
          ndispl_max_oc,nveloc_max_oc)
        if (COMPUTE_AND_STORE_STRAIN) then	
          call write_bin_strain_cm(nstrain_max_cm)
        endif
        call write_bin_stress_cm(nstress_max_cm)
      endif
    endif
  endif
  

  end subroutine iterate_time

!
!-------------------------------------------------------------------------------------------------
!

  subroutine it_transfer_from_GPU()

! transfers fields on GPU back onto CPU

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore
  use specfem_par_outercore
  use specfem_par_noise

  implicit none

  ! to store forward wave fields
  if (SIMULATION_TYPE == 1) then
    if (SAVE_FORWARD .or. (NUMBER_OF_RUNS > 1 .and. NUMBER_OF_THIS_RUN < NUMBER_OF_RUNS)) then
      ! wavefield
      call transfer_fields_cm_from_device(NDIM*NGLOB_CRUST_MANTLE, &
                                          displ_crust_mantle,veloc_crust_mantle,accel_crust_mantle,Mesh_pointer)
      call transfer_fields_ic_from_device(NDIM*NGLOB_INNER_CORE, &
                                          displ_inner_core,veloc_inner_core,accel_inner_core,Mesh_pointer)
      call transfer_fields_oc_from_device(NGLOB_OUTER_CORE, &
                                          displ_outer_core,veloc_outer_core,accel_outer_core,Mesh_pointer)
      ! strain
      call transfer_strain_cm_from_device(Mesh_pointer,eps_trace_over_3_crust_mantle, &
                                          epsilondev_xx_crust_mantle,epsilondev_yy_crust_mantle, &
                                          epsilondev_xy_crust_mantle,epsilondev_xz_crust_mantle, &
                                          epsilondev_yz_crust_mantle)
      call transfer_strain_ic_from_device(Mesh_pointer,eps_trace_over_3_inner_core, &
                                          epsilondev_xx_inner_core,epsilondev_yy_inner_core, &
                                          epsilondev_xy_inner_core,epsilondev_xz_inner_core, &
                                          epsilondev_yz_inner_core)
      ! rotation
      if (ROTATION_VAL) then
        call transfer_rotation_from_device(Mesh_pointer,A_array_rotation,B_array_rotation)
      endif

      ! attenuation memory variables
      if (ATTENUATION_VAL .and. .not. PARTIAL_PHYS_DISPERSION_ONLY_VAL) then
        call transfer_rmemory_cm_from_device(Mesh_pointer,R_xx_crust_mantle,R_yy_crust_mantle,R_xy_crust_mantle, &
                                             R_xz_crust_mantle,R_yz_crust_mantle)
        call transfer_rmemory_ic_from_device(Mesh_pointer,R_xx_inner_core,R_yy_inner_core,R_xy_inner_core, &
                                             R_xz_inner_core,R_yz_inner_core)
      endif
    endif

  else if (SIMULATION_TYPE == 3) then

    ! note: for kernel simulations (SIMULATION_TYPE == 3), attenuation is
    !          only mimicking effects on phase shifts, but not on amplitudes.
    !          flag PARTIAL_PHYS_DISPERSION_ONLY will have to be set to true in this case.
    !
    ! arrays b_R_xx, ... are not used when PARTIAL_PHYS_DISPERSION_ONLY is set,
    ! therefore no need to transfer arrays from GPU to CPU
    !if (ATTENUATION) then
    !endif

    ! to store kernels
    ! crust/mantle
    call transfer_kernels_cm_to_host(Mesh_pointer, &
                                     rho_kl_crust_mantle,alpha_kl_crust_mantle,beta_kl_crust_mantle, &
                                     NSPEC_CRUST_MANTLE)

    ! full anisotropic kernel
    if (ANISOTROPIC_KL) then
      call transfer_kernels_ani_cm_to_host(Mesh_pointer,cijkl_kl_crust_mantle,NSPEC_CRUST_MANTLE)
    endif

    ! specific noise strength kernel
    if (NOISE_TOMOGRAPHY == 3) then
      call transfer_kernels_noise_to_host(Mesh_pointer,sigma_kl_crust_mantle,NSPEC_CRUST_MANTLE)
    endif

    ! approximative Hessian for preconditioning kernels
    if (APPROXIMATE_HESS_KL) then
      call transfer_kernels_hess_cm_tohost(Mesh_pointer, &
                                           hess_kl_crust_mantle, &
                                           hess_rho_kl_crust_mantle, &
                                           hess_kappa_kl_crust_mantle, &
                                           hess_mu_kl_crust_mantle, &
                                           NSPEC_CRUST_MANTLE)
    endif

    ! outer core
    if (SAVE_KERNELS_OC) then
      call transfer_kernels_oc_to_host(Mesh_pointer, &
                                       rho_kl_outer_core, &
                                       alpha_kl_outer_core,NSPEC_OUTER_CORE)
    endif

    ! inner core
    if (SAVE_KERNELS_IC) then
      call transfer_kernels_ic_to_host(Mesh_pointer, &
                                       rho_kl_inner_core, &
                                       alpha_kl_inner_core, &
                                       beta_kl_inner_core,NSPEC_INNER_CORE)
    endif
  endif

  end subroutine it_transfer_from_GPU

!
!-------------------------------------------------------------------------------------------------
!

  subroutine it_update_vtkwindow()

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_movie

  implicit none

  real :: currenttime
  integer :: iglob,inum,data_size
  real, dimension(1) :: dummy

  ! VTK rendering at frame interval
  if (mod(it,NTSTEP_BETWEEN_FRAMES) == 0) then

    ! user output
    !if (myrank == 0 ) print *,"  VTK rendering..."

    ! updates time
    currenttime = sngl((it-1)*DT-t0)

    ! transfers fields from GPU to host
    if (GPU_MODE) then
      !if (myrank == 0 ) print *,"  VTK: transferring velocity from GPU"
      call transfer_veloc_cm_from_device(NDIM*NGLOB_CRUST_MANTLE,veloc_crust_mantle,Mesh_pointer)
    endif

    ! updates wavefield
    !if (myrank == 0 ) print *,"  VTK: it = ",it," out of ",it_end," - norm of velocity field"
    inum = 0
    vtkdata(:) = 0.0
    do iglob = 1,NGLOB_CRUST_MANTLE
      if (vtkmask(iglob) .eqv. .true.) then
        inum = inum + 1
        ! stores norm of velocity vector
        vtkdata(inum) = sqrt(veloc_crust_mantle(1,iglob)**2 &
                           + veloc_crust_mantle(2,iglob)**2 &
                           + veloc_crust_mantle(3,iglob)**2)
      endif
    enddo

    ! updates for multiple MPI process
    if (NPROCTOT_VAL > 1) then
      data_size = size(vtkdata)
      if (myrank == 0) then
        ! gather data
        call gatherv_all_r(vtkdata,data_size, &
                            vtkdata_all,vtkdata_points_all,vtkdata_offset_all, &
                            vtkdata_numpoints_all,NPROCTOT_VAL)
        ! updates VTK window
        call visualize_vtkdata(it,currenttime,vtkdata_all)
      else
        ! all other process just send data
        call gatherv_all_r(vtkdata,data_size, &
                            dummy,vtkdata_points_all,vtkdata_offset_all, &
                            1,NPROCTOT_VAL)
      endif
    else
      ! serial run
      ! updates VTK window
      call visualize_vtkdata(it,currenttime,vtkdata)
    endif

  endif

  end subroutine it_update_vtkwindow

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gravity_timeseries()

  implicit none

  stop 'gravity_timeseries() not implemented in this code yet'

  end subroutine gravity_timeseries


!
!-------------------------------------------------------------------------------------------------
!

  subroutine setup_exact_undoing_to_disk()

  use specfem_par
  use specfem_par_crustmantle

  implicit none

  ! local parameters
  integer :: ispec,iglob,i,j,k
  integer :: counter,record_length
  real(kind=CUSTOM_REAL) :: radius
  character(len=MAX_STRING_LEN) :: outputname

  ! checks if anything to do
  if (.not. EXACT_UNDOING_TO_DISK) return

  ! checks flags
  if (GPU_MODE) call exit_MPI(myrank,'EXACT_UNDOING_TO_DISK not supported for GPUs')

  if (UNDO_ATTENUATION) &
    call exit_MPI(myrank,'EXACT_UNDOING_TO_DISK needs UNDO_ATTENUATION to be off because it computes the kernel directly instead')

  if (SIMULATION_TYPE == 1 .and. .not. SAVE_FORWARD) &
    call exit_MPI(myrank,'EXACT_UNDOING_TO_DISK requires SAVE_FORWARD if SIMULATION_TYPE == 1')

  if (ANISOTROPIC_KL) call exit_MPI(myrank,'EXACT_UNDOING_TO_DISK requires ANISOTROPIC_KL to be turned off')

  if (SIMULATION_TYPE /= 1 .and. SIMULATION_TYPE /= 3) &
    call exit_MPI(myrank,'EXACT_UNDOING_TO_DISK can only be used with SIMULATION_TYPE == 1 or SIMULATION_TYPE == 3')


!! DK DK determine the largest value of iglob that we need to save to disk,
!! DK DK since we save the upper part of the mesh only in the case of surface-wave kernels
  ! crust_mantle
  allocate(integer_mask_ibool_exact_undo(NGLOB_CRUST_MANTLE))
  integer_mask_ibool_exact_undo(:) = -1

  counter = 0
  do ispec = 1, NSPEC_CRUST_MANTLE
    do k = 1, NGLLZ
      do j = 1, NGLLY
        do i = 1, NGLLX
          iglob = ibool_crust_mantle(i,j,k,ispec)
          ! xstore ystore zstore have previously been converted to r theta phi in rstore
          radius = rstore_crust_mantle(1,iglob) ! radius r (normalized)
          ! save that element only if it is in the upper part of the mesh
          if (radius >= R670 / R_EARTH) then
            ! if this point has not yet been found before
            if (integer_mask_ibool_exact_undo(iglob) == -1) then
              ! create a new unique point
              counter = counter + 1
              integer_mask_ibool_exact_undo(iglob) = counter
            endif
          endif
        enddo
      enddo
    enddo
  enddo

  ! allocate the buffer used to dump a single time step
  allocate(buffer_for_disk(counter,NSTEP_FOR_EXACT_UNDOING))

  ! open the file in which we will dump all the time steps (in a single file)
  write(outputname,"('huge_dumps/proc',i6.6,'_huge_dump_of_all_time_steps.bin')") myrank
  inquire(iolength=record_length) buffer_for_disk(:,1)
  ! we write to or read from the file depending on the simulation type
  if (SIMULATION_TYPE == 1) then
    open(file=outputname, unit=IFILE_FOR_EXACT_UNDOING, action='write', status='unknown', &
                    form='unformatted', access='direct', recl=record_length)
  else if (SIMULATION_TYPE == 3) then
    open(file=outputname, unit=IFILE_FOR_EXACT_UNDOING, action='read', status='old', &
                    form='unformatted', access='direct', recl=record_length)
  endif

  if (SIMULATION_TYPE == 1) then
    it_of_this_exact_subset = 1
    it_exact_subset_offset = it_begin - 1
    it_exact_subset_end = min(NSTEP_FOR_EXACT_UNDOING, it_end - it_begin + 1)
  else if (SIMULATION_TYPE == 3) then
    ! Trigger a read at the start of the loop
    it_of_this_exact_subset = 0
    it_exact_subset_offset = it_begin - 1
    it_exact_subset_end = 0
  endif

  end subroutine setup_exact_undoing_to_disk

!
!-------------------------------------------------------------------------------------------------
!

  subroutine finish_exact_undoing_to_disk()

  use specfem_par
  use specfem_par_crustmantle

  implicit none

  ! checks if anything to do
  if (.not. EXACT_UNDOING_TO_DISK) return

  ! frees memory
  deallocate(integer_mask_ibool_exact_undo)
  deallocate(buffer_for_disk)

  ! close the huge file that contains a dump of all the time steps to disk
  close(IFILE_FOR_EXACT_UNDOING)

  end subroutine finish_exact_undoing_to_disk

!
!-------------------------------------------------------------------------------------------------
!
  
  subroutine max_norms_dispvel(maxnormdisp_cm,maxnormvel_cm,maxnormdisp_ic, &
				maxnormvel_ic,maxnormdisp_oc,maxnormvel_oc)
  
  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore
  use specfem_par_outercore
  
  real(kind=CUSTOM_REAL) :: ndispl_cm,nveloc_cm,ndispl_ic,nveloc_ic
  real(kind=CUSTOM_REAL), dimension(NGLOB_CRUST_MANTLE),intent(inout) :: maxnormdisp_cm,maxnormvel_cm
  real(kind=CUSTOM_REAL), dimension(NGLOB_INNER_CORE),intent(inout) :: maxnormdisp_ic,maxnormvel_ic
  real(kind=CUSTOM_REAL), dimension(NGLOB_OUTER_CORE),intent(inout) :: maxnormdisp_oc,maxnormvel_oc
  
  integer :: iglob
  
  ! crust and mantle
  do iglob = 1, NGLOB_CRUST_MANTLE
    ndispl_cm=scale_displ*sqrt(displ_crust_mantle(1,iglob)**2 &
		+ displ_crust_mantle(2,iglob)**2 &
		+ displ_crust_mantle(3,iglob)**2)
    maxnormdisp_cm(iglob)=max(ndispl_cm,maxnormdisp_cm(iglob))

    nveloc_cm=scale_veloc*sqrt(veloc_crust_mantle(1,iglob)**2 &
		+ veloc_crust_mantle(2,iglob)**2 &
		+ veloc_crust_mantle(3,iglob)**2)
    maxnormvel_cm(iglob)=max(nveloc_cm,maxnormvel_cm(iglob))	
  enddo
	  
  ! inner core
  do iglob = 1, NGLOB_INNER_CORE
    ndispl_ic=scale_displ*sqrt(displ_inner_core(1,iglob)**2 &
		+ displ_inner_core(2,iglob)**2 &
		+ displ_inner_core(3,iglob)**2)
    maxnormdisp_ic(iglob)=max(ndispl_ic,maxnormdisp_ic(iglob))

    nveloc_ic=scale_veloc*sqrt(veloc_inner_core(1,iglob)**2 &
		+ veloc_inner_core(2,iglob)**2 &
		+ veloc_inner_core(3,iglob)**2)
    maxnormvel_ic(iglob)=max(nveloc_ic,maxnormvel_ic(iglob))	
  enddo
	  
  ! outer core
  do iglob = 1, NGLOB_OUTER_CORE
    maxnormdisp_oc(iglob)=max(abs(displ_outer_core(iglob)),maxnormdisp_oc(iglob))
    maxnormvel_oc(iglob)=max(abs(veloc_outer_core(iglob)),maxnormvel_oc(iglob))
  enddo
  
  
  end subroutine max_norms_dispvel

!
!-------------------------------------------------------------------------------------------------
!
  
  subroutine max_norms_strain_cm(maxnormstrain_cm)
  
  use constants_solver
  use specfem_par
  use specfem_par_crustmantle
  
  ! Norm of strain
  real(kind=CUSTOM_REAL) :: nstrain_cm
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE), &
       intent(inout) :: maxnormstrain_cm
  
  ! Components of the strain tensor
  real(kind=CUSTOM_REAL) :: epsilon_xx,epsilon_yy,epsilon_zz,epsilon_xy,epsilon_xz,epsilon_yz 
  real(kind=CUSTOM_REAL), dimension(5,NGLLX,NGLLY,NGLLZ) :: epsilondev_loc_matrix
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: eps_trace_over_3_loc_matrix
  
  integer :: ispec,iglob,i,j,k
  
  ! Compute the strain from displacement in the crust/mantle
  do ispec = 1, NSPEC_CRUST_MANTLE_STR_OR_ATT				
    do k=1, NGLLZ
      do j=1, NGLLY
        do i=1, NGLLX
          epsilon_xx=eps_trace_over_3_crust_mantle(i,j,k,ispec) & 
						+ epsilondev_xx_crust_mantle(i,j,k,ispec)
          epsilon_yy=eps_trace_over_3_crust_mantle(i,j,k,ispec) & 
						+ epsilondev_yy_crust_mantle(i,j,k,ispec)
          epsilon_zz=eps_trace_over_3_crust_mantle(i,j,k,ispec) &
					   - epsilondev_xx_crust_mantle(i,j,k,ispec) &
					   - epsilondev_yy_crust_mantle(i,j,k,ispec)
          epsilon_xy = epsilondev_xy_crust_mantle(i,j,k,ispec)
          epsilon_xz = epsilondev_xz_crust_mantle(i,j,k,ispec)
          epsilon_yz = epsilondev_yz_crust_mantle(i,j,k,ispec)
		  
          ! Compute the norm of strain, and track the maximum norm across all timesteps
          nstrain_cm=sqrt(epsilon_xx**2 & 
                + epsilon_yy**2 &
                + epsilon_zz**2 &
                + epsilon_xy**2 &
                + epsilon_xz**2 &
                + epsilon_yz**2 & 
                + epsilon_xy**2 &
                + epsilon_xz**2 &
                + epsilon_yz**2)
		  
          maxnormstrain_cm(i,j,k,ispec)=max(nstrain_cm,maxnormstrain_cm(i,j,k,ispec))
        enddo
      enddo
    enddo
  enddo
  
  
  end subroutine max_norms_strain_cm
  
!
!--------------------------------------------------------------------------------------------
!
  subroutine compute_element_stress_cm_Dev(ispec,nglob,nspec, &
                                           displ,ibool, &
                                           hprime_xxl,hprime_xxTl, &
                                           deriv, &
                                           sigma_xx,sigma_yy,sigma_zz, &
                                           sigma_xy,sigma_xz,sigma_yz, &
                                           sigma_yx,sigma_zx,sigma_zy)

! computes stress for single element based on Deville routine setup (NGLLX == NGLLY == NGLLZ == 5)
! in the crust/mantle

  use constants
  use specfem_par
  use specfem_par_crustmantle
  

  implicit none

  integer,intent(in) :: ispec,nglob,nspec

  real(kind=CUSTOM_REAL), dimension(NDIM,nglob),intent(in) :: displ
  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec),intent(in) :: ibool

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLX),intent(in) :: hprime_xxl,hprime_xxTl

  real(kind=CUSTOM_REAL), dimension(9,NGLLX,NGLLY,NGLLZ,nspec),intent(in) :: deriv
  
  ! Components of the stress tensor
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ),intent(out) :: sigma_xx,sigma_yy,sigma_zz, &
       sigma_xy,sigma_xz,sigma_yz, &
       sigma_yx,sigma_zx,sigma_zy

  !  local variable
  integer :: iglob
  integer i_SLS

  !real(kind=CUSTOM_REAL) :: templ
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: &
       tempx1,tempx2,tempx3,tempy1,tempy2,tempy3,tempz1,tempz2,tempz3

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: dummyx_loc,dummyy_loc,dummyz_loc

  real(kind=CUSTOM_REAL) :: xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl
  real(kind=CUSTOM_REAL) :: duxdxl,duydyl,duzdzl,duxdyl,duydxl,duzdxl,duxdzl,duzdyl,duydzl, &
                            duxdxl_plus_duydyl,duxdxl_plus_duzdzl,duydyl_plus_duzdzl, &
                            duxdyl_plus_duydxl,duzdxl_plus_duxdzl,duzdyl_plus_duydzl
							
  ! Components of the stress tensor at each i,j,k
  real(kind=CUSTOM_REAL) sigma_xxl,sigma_yyl,sigma_zzl,sigma_xyl,sigma_xzl,sigma_yzl, &
       sigma_yxl,sigma_zxl,sigma_zyl
  ! the 21 coefficients for an anisotropic medium in reduced notation
  real(kind=CUSTOM_REAL) c11,c22,c33,c44,c55,c66,c12,c13,c23,c14,c24,c34,c15,c25,c35,c45,c16,c26,c36,c46,c56
  real(kind=CUSTOM_REAL) lambdal,mul,lambdalplus2mul
  real(kind=CUSTOM_REAL) kappal
  real(kind=CUSTOM_REAL) R_xx_val,R_yy_val
  !real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,5) :: epsilondev_loc
  real(kind=CUSTOM_REAL) sx_l,sy_l,sz_l,gxl,gyl,gzl

#ifdef FORCE_VECTORIZATION
  integer :: ijk
#else
  integer :: i,j,k
#endif

  DO_LOOP_IJK

    iglob = ibool(INDEX_IJK,ispec)
    dummyx_loc(INDEX_IJK) = displ(1,iglob)
    dummyy_loc(INDEX_IJK) = displ(2,iglob)
    dummyz_loc(INDEX_IJK) = displ(3,iglob)

  ENDDO_LOOP_IJK

  ! Deville optimizations
  ! computes 1. matrix multiplication for tempx1,..
  call mxm5_3comp_singleA(hprime_xxl,m1,dummyx_loc,dummyy_loc,dummyz_loc,tempx1,tempy1,tempz1,m2)
  ! computes 2. matrix multiplication for tempx2,..
  call mxm5_3comp_3dmat_singleB(dummyx_loc,dummyy_loc,dummyz_loc,m1,hprime_xxTl,m1,tempx2,tempy2,tempz2,NGLLX)
  ! computes 3. matrix multiplication for tempx1,..
  call mxm5_3comp_singleB(dummyx_loc,dummyy_loc,dummyz_loc,m2,hprime_xxTl,tempx3,tempy3,tempz3,m1)

  DO_LOOP_IJK

    ! get derivatives of ux, uy and uz with respect to x, y and z
    xixl = deriv(1,INDEX_IJK,ispec)
    xiyl = deriv(2,INDEX_IJK,ispec)
    xizl = deriv(3,INDEX_IJK,ispec)
    etaxl = deriv(4,INDEX_IJK,ispec)
    etayl = deriv(5,INDEX_IJK,ispec)
    etazl = deriv(6,INDEX_IJK,ispec)
    gammaxl = deriv(7,INDEX_IJK,ispec)
    gammayl = deriv(8,INDEX_IJK,ispec)
    gammazl = deriv(9,INDEX_IJK,ispec)

    duxdxl = xixl*tempx1(INDEX_IJK) + etaxl*tempx2(INDEX_IJK) + gammaxl*tempx3(INDEX_IJK)
    duxdyl = xiyl*tempx1(INDEX_IJK) + etayl*tempx2(INDEX_IJK) + gammayl*tempx3(INDEX_IJK)
    duxdzl = xizl*tempx1(INDEX_IJK) + etazl*tempx2(INDEX_IJK) + gammazl*tempx3(INDEX_IJK)

    duydxl = xixl*tempy1(INDEX_IJK) + etaxl*tempy2(INDEX_IJK) + gammaxl*tempy3(INDEX_IJK)
    duydyl = xiyl*tempy1(INDEX_IJK) + etayl*tempy2(INDEX_IJK) + gammayl*tempy3(INDEX_IJK)
    duydzl = xizl*tempy1(INDEX_IJK) + etazl*tempy2(INDEX_IJK) + gammazl*tempy3(INDEX_IJK)

    duzdxl = xixl*tempz1(INDEX_IJK) + etaxl*tempz2(INDEX_IJK) + gammaxl*tempz3(INDEX_IJK)
    duzdyl = xiyl*tempz1(INDEX_IJK) + etayl*tempz2(INDEX_IJK) + gammayl*tempz3(INDEX_IJK)
    duzdzl = xizl*tempz1(INDEX_IJK) + etazl*tempz2(INDEX_IJK) + gammazl*tempz3(INDEX_IJK)

    ! precompute some sums to save CPU time
    duxdxl_plus_duydyl = duxdxl + duydyl
    duxdxl_plus_duzdzl = duxdxl + duzdzl
    duydyl_plus_duzdzl = duydyl + duzdzl
    duxdyl_plus_duydxl = duxdyl + duydxl
    duzdxl_plus_duxdzl = duzdxl + duxdzl
    duzdyl_plus_duydzl = duzdyl + duydzl


    if (ANISOTROPIC_3D_MANTLE_VAL) then

      c11 = c11store_crust_mantle(INDEX_IJK,ispec)
      c12 = c12store_crust_mantle(INDEX_IJK,ispec)
      c13 = c13store_crust_mantle(INDEX_IJK,ispec)
      c14 = c14store_crust_mantle(INDEX_IJK,ispec)
      c15 = c15store_crust_mantle(INDEX_IJK,ispec)
      c16 = c16store_crust_mantle(INDEX_IJK,ispec)
      c22 = c22store_crust_mantle(INDEX_IJK,ispec)
      c23 = c23store_crust_mantle(INDEX_IJK,ispec)
      c24 = c24store_crust_mantle(INDEX_IJK,ispec)
      c25 = c25store_crust_mantle(INDEX_IJK,ispec)
      c26 = c26store_crust_mantle(INDEX_IJK,ispec)
      c33 = c33store_crust_mantle(INDEX_IJK,ispec)
      c34 = c34store_crust_mantle(INDEX_IJK,ispec)
      c35 = c35store_crust_mantle(INDEX_IJK,ispec)
      c36 = c36store_crust_mantle(INDEX_IJK,ispec)
      c44 = c44store_crust_mantle(INDEX_IJK,ispec)
      c45 = c45store_crust_mantle(INDEX_IJK,ispec)
      c46 = c46store_crust_mantle(INDEX_IJK,ispec)
      c55 = c55store_crust_mantle(INDEX_IJK,ispec)
      c56 = c56store_crust_mantle(INDEX_IJK,ispec)
      c66 = c66store_crust_mantle(INDEX_IJK,ispec)
      
      sigma_xxl = c11*duxdxl + c16*duxdyl_plus_duydxl + c12*duydyl + &
           c15*duzdxl_plus_duxdzl + c14*duzdyl_plus_duydzl + c13*duzdzl
      
      sigma_yyl = c12*duxdxl + c26*duxdyl_plus_duydxl + c22*duydyl + &
           c25*duzdxl_plus_duxdzl + c24*duzdyl_plus_duydzl + c23*duzdzl

      sigma_zzl = c13*duxdxl + c36*duxdyl_plus_duydxl + c23*duydyl + &
           c35*duzdxl_plus_duxdzl + c34*duzdyl_plus_duydzl + c33*duzdzl

      sigma_xyl = c16*duxdxl + c66*duxdyl_plus_duydxl + c26*duydyl + &
           c56*duzdxl_plus_duxdzl + c46*duzdyl_plus_duydzl + c36*duzdzl

      sigma_xzl = c15*duxdxl + c56*duxdyl_plus_duydxl + c25*duydyl + &
           c55*duzdxl_plus_duxdzl + c45*duzdyl_plus_duydzl + c35*duzdzl

      sigma_yzl = c14*duxdxl + c46*duxdyl_plus_duydxl + c24*duydyl + &
           c45*duzdxl_plus_duxdzl + c44*duzdyl_plus_duydzl + c34*duzdzl

   else

     ! do not use transverse isotropy except if element is between d220 and Moho
      if (.not. ispec_is_tiso_crust_mantle(ispec)) then

        ! isotropic element
         
        ! layer with no transverse isotropy, use kappav and muv
        kappal = kappavstore_crust_mantle(INDEX_IJK,ispec)
        mul = muvstore_crust_mantle(INDEX_IJK,ispec)

        ! use unrelaxed parameters if attenuation
        ! already done in prepare_timerun...
        !if (ATTENUATION_VAL) mul = mul * one_minus_sum_beta_use
        
        lambdalplus2mul = kappal + FOUR_THIRDS * mul
        lambdal = lambdalplus2mul - 2.*mul

        ! compute stress sigma
        sigma_xxl = lambdalplus2mul*duxdxl + lambdal*duydyl_plus_duzdzl
        sigma_yyl = lambdalplus2mul*duydyl + lambdal*duxdxl_plus_duzdzl
        sigma_zzl = lambdalplus2mul*duzdzl + lambdal*duxdxl_plus_duydyl
        sigma_xyl = mul*duxdyl_plus_duydxl
        sigma_xzl = mul*duzdxl_plus_duxdzl
        sigma_yzl = mul*duzdyl_plus_duydzl

      else
	! transverse isotropic element
	c11 = c11store_crust_mantle(INDEX_IJK,ispec)
	c12 = c12store_crust_mantle(INDEX_IJK,ispec)
	c13 = c13store_crust_mantle(INDEX_IJK,ispec)
	c14 = c14store_crust_mantle(INDEX_IJK,ispec)
	c15 = c15store_crust_mantle(INDEX_IJK,ispec)
	c16 = c16store_crust_mantle(INDEX_IJK,ispec)
	c22 = c22store_crust_mantle(INDEX_IJK,ispec)
	c23 = c23store_crust_mantle(INDEX_IJK,ispec)
	c24 = c24store_crust_mantle(INDEX_IJK,ispec)
	c25 = c25store_crust_mantle(INDEX_IJK,ispec)
	c26 = c26store_crust_mantle(INDEX_IJK,ispec)
	c33 = c33store_crust_mantle(INDEX_IJK,ispec)
	c34 = c34store_crust_mantle(INDEX_IJK,ispec)
	c35 = c35store_crust_mantle(INDEX_IJK,ispec)
	c36 = c36store_crust_mantle(INDEX_IJK,ispec)
	c44 = c44store_crust_mantle(INDEX_IJK,ispec)
	c45 = c45store_crust_mantle(INDEX_IJK,ispec)
	c46 = c46store_crust_mantle(INDEX_IJK,ispec)
	c55 = c55store_crust_mantle(INDEX_IJK,ispec)
	c56 = c56store_crust_mantle(INDEX_IJK,ispec)
	c66 = c66store_crust_mantle(INDEX_IJK,ispec)

	! general expression of stress tensor for full Cijkl with 21 coefficients
	sigma_xxl = c11*duxdxl + c16*duxdyl_plus_duydxl + c12*duydyl + &
			   c15*duzdxl_plus_duxdzl + c14*duzdyl_plus_duydzl + c13*duzdzl
	sigma_yyl = c12*duxdxl + c26*duxdyl_plus_duydxl + c22*duydyl + &
			   c25*duzdxl_plus_duxdzl + c24*duzdyl_plus_duydzl + c23*duzdzl
	sigma_zzl = c13*duxdxl + c36*duxdyl_plus_duydxl + c23*duydyl + &
			   c35*duzdxl_plus_duxdzl + c34*duzdyl_plus_duydzl + c33*duzdzl
	sigma_xyl = c16*duxdxl + c66*duxdyl_plus_duydxl + c26*duydyl + &
			   c56*duzdxl_plus_duxdzl + c46*duzdyl_plus_duydzl + c36*duzdzl
	sigma_xz = c15*duxdxl + c56*duxdyl_plus_duydxl + c25*duydyl + &
                           c55*duzdxl_plus_duxdzl + c45*duzdyl_plus_duydzl + c35*duzdzl
        sigma_yz = c14*duxdxl + c46*duxdyl_plus_duydxl + c24*duydyl + &
			   c45*duzdxl_plus_duxdzl + c44*duzdyl_plus_duydzl + c34*duzdzl

      endif

   endif   ! end of test whether isotropic or anisotropic element

   ! subtract memory variables if attenuation
   if (ATTENUATION_VAL .and. .not. PARTIAL_PHYS_DISPERSION_ONLY_VAL) then
     do i_SLS = 1,N_SLS
        R_xx_val = R_xx_crust_mantle(INDEX_IJK,i_SLS,ispec)
        R_yy_val = R_yy_crust_mantle(INDEX_IJK,i_SLS,ispec)
        sigma_xxl = sigma_xxl - R_xx_val
        sigma_yyl = sigma_yyl - R_yy_val
        sigma_zzl = sigma_zzl + R_xx_val + R_yy_val
        sigma_xyl = sigma_xyl - R_xy_crust_mantle(INDEX_IJK,i_SLS,ispec)
        sigma_xzl = sigma_xzl - R_xz_crust_mantle(INDEX_IJK,i_SLS,ispec)
        sigma_yzl = sigma_yzl - R_yz_crust_mantle(INDEX_IJK,i_SLS,ispec)
     enddo
   endif

   ! define symmetric components of sigma for gravity
   sigma_yxl = sigma_xyl
   sigma_zxl = sigma_xzl
   sigma_zyl = sigma_yzl

   ! compute non-symmetric terms for gravity
   if (GRAVITY_VAL) then
     ! use mesh coordinates to get theta and phi
     ! x y and z contain r theta and phi
     iglob = ibool(INDEX_IJK,ispec)

     ! Cartesian components of the gravitational acceleration
     gxl = gravity_pre_store_crust_mantle(1,iglob) ! minus_g*sin_theta*cos_phi * rho
     gyl = gravity_pre_store_crust_mantle(2,iglob) ! minus_g*sin_theta*sin_phi * rho
     gzl = gravity_pre_store_crust_mantle(3,iglob) ! minus_g*cos_theta * rho

     ! Cartesian components of gradient of gravitational acceleration
     ! get displacement and multiply by density to compute G tensor
     sx_l = displ_crust_mantle(1,iglob)
     sy_l = displ_crust_mantle(2,iglob)
     sz_l = displ_crust_mantle(3,iglob)
     
     ! compute G tensor from s . g and add to sigma (not symmetric)
     sigma_xxl = sigma_xxl + sy_l * gyl + sz_l * gzl
     sigma_yyl = sigma_yyl + sx_l * gxl + sz_l * gzl
     sigma_zzl = sigma_zzl + sx_l * gxl + sy_l * gyl
     
     sigma_xyl = sigma_xyl - sx_l * gyl
     sigma_yxl = sigma_yxl - sy_l * gxl
     
     sigma_xzl = sigma_xzl - sx_l * gzl
     sigma_zxl = sigma_zxl - sz_l * gxl

     sigma_yzl = sigma_yzl - sy_l * gzl
     sigma_zyl = sigma_zyl - sz_l * gyl
     
   endif  ! end of section with gravity terms
	
   ! Populate the stress tensor components
   sigma_xx(INDEX_IJK)=sigma_xxl
   sigma_xy(INDEX_IJK)=sigma_xyl
   sigma_xz(INDEX_IJK)=sigma_xzl
   sigma_yx(INDEX_IJK)=sigma_yxl
   sigma_yy(INDEX_IJK)=sigma_yyl
   sigma_yz(INDEX_IJK)=sigma_yzl
   sigma_zx(INDEX_IJK)=sigma_zxl
   sigma_zy(INDEX_IJK)=sigma_zyl
   sigma_zz(INDEX_IJK)=sigma_zzl

	
  ENDDO_LOOP_IJK
  
  
  contains

!--------------------------------------------------------------------------------------------
!
! matrix-matrix multiplications
!
! subroutines adapted from Deville, Fischer and Mund, High-order methods
! for incompressible fluid flow, Cambridge University Press (2002),
! pages 386 and 389 and Figure 8.3.1
!
!--------------------------------------------------------------------------------------------
!
! note: the matrix-matrix multiplications are used for very small matrices ( 5 x 5 x 5 elements);
!       thus, calling external optimized libraries for these multiplications are in general slower
!
! please leave the routines here to help compilers inlining the code

  subroutine mxm5_3comp_singleA(A,n1,B1,B2,B3,C1,C2,C3,n3)

! 3 different arrays for x/y/z-components, 2-dimensional arrays (25,5)/(5,25), same B matrix for all 3 component arrays

  use constants_solver, only: CUSTOM_REAL

  implicit none

  integer,intent(in) :: n1,n3
  real(kind=CUSTOM_REAL),dimension(n1,5) :: A
  real(kind=CUSTOM_REAL),dimension(5,n3) :: B1,B2,B3
  real(kind=CUSTOM_REAL),dimension(n1,n3) :: C1,C2,C3

  ! local parameters
  integer :: i,j

  ! matrix-matrix multiplication
  do j = 1,n3
    do i = 1,n1
      C1(i,j) =  A(i,1) * B1(1,j) &
               + A(i,2) * B1(2,j) &
               + A(i,3) * B1(3,j) &
               + A(i,4) * B1(4,j) &
               + A(i,5) * B1(5,j)

      C2(i,j) =  A(i,1) * B2(1,j) &
               + A(i,2) * B2(2,j) &
               + A(i,3) * B2(3,j) &
               + A(i,4) * B2(4,j) &
               + A(i,5) * B2(5,j)

      C3(i,j) =  A(i,1) * B3(1,j) &
               + A(i,2) * B3(2,j) &
               + A(i,3) * B3(3,j) &
               + A(i,4) * B3(4,j) &
               + A(i,5) * B3(5,j)
    enddo
  enddo

  end subroutine mxm5_3comp_singleA


!--------------------------------------------------------------------------------------------

  subroutine mxm5_3comp_singleB(A1,A2,A3,n1,B,C1,C2,C3,n3)

! 3 different arrays for x/y/z-components, 2-dimensional arrays (25,5)/(5,25), same B matrix for all 3 component arrays

  use constants_solver, only: CUSTOM_REAL

  implicit none

  integer,intent(in) :: n1,n3
  real(kind=CUSTOM_REAL),dimension(n1,5) :: A1,A2,A3
  real(kind=CUSTOM_REAL),dimension(5,n3) :: B
  real(kind=CUSTOM_REAL),dimension(n1,n3) :: C1,C2,C3

  ! local parameters
  integer :: i,j

  ! matrix-matrix multiplication
  do j = 1,n3
    do i = 1,n1
      C1(i,j) =  A1(i,1) * B(1,j) &
               + A1(i,2) * B(2,j) &
               + A1(i,3) * B(3,j) &
               + A1(i,4) * B(4,j) &
               + A1(i,5) * B(5,j)

      C2(i,j) =  A2(i,1) * B(1,j) &
               + A2(i,2) * B(2,j) &
               + A2(i,3) * B(3,j) &
               + A2(i,4) * B(4,j) &
               + A2(i,5) * B(5,j)

      C3(i,j) =  A3(i,1) * B(1,j) &
               + A3(i,2) * B(2,j) &
               + A3(i,3) * B(3,j) &
               + A3(i,4) * B(4,j) &
               + A3(i,5) * B(5,j)
    enddo
  enddo

  end subroutine mxm5_3comp_singleB


!--------------------------------------------------------------------------------------------

  subroutine mxm5_3comp_3dmat_singleB(A1,A2,A3,n1,B,n2,C1,C2,C3,n3)

! 3 different arrays for x/y/z-components, 3-dimensional arrays (5,5,5), same B matrix for all 3 component arrays

  use constants_solver, only: CUSTOM_REAL

  implicit none

  integer,intent(in) :: n1,n2,n3
  real(kind=CUSTOM_REAL),dimension(n1,5,n3) :: A1,A2,A3
  real(kind=CUSTOM_REAL),dimension(5,n2) :: B
  real(kind=CUSTOM_REAL),dimension(n1,n2,n3) :: C1,C2,C3

  ! local parameters
  integer :: i,j,k

  ! matrix-matrix multiplication
  do j = 1,n2
    do i = 1,n1
      ! for efficiency it is better to leave this loop on k inside, it leads to slightly faster code
      do k = 1,n3
        C1(i,j,k) =  A1(i,1,k) * B(1,j) &
                   + A1(i,2,k) * B(2,j) &
                   + A1(i,3,k) * B(3,j) &
                   + A1(i,4,k) * B(4,j) &
                   + A1(i,5,k) * B(5,j)

        C2(i,j,k) =  A2(i,1,k) * B(1,j) &
                   + A2(i,2,k) * B(2,j) &
                   + A2(i,3,k) * B(3,j) &
                   + A2(i,4,k) * B(4,j) &
                   + A2(i,5,k) * B(5,j)

        C3(i,j,k) =  A3(i,1,k) * B(1,j) &
                   + A3(i,2,k) * B(2,j) &
                   + A3(i,3,k) * B(3,j) &
                   + A3(i,4,k) * B(4,j) &
                   + A3(i,5,k) * B(5,j)
      enddo
    enddo
  enddo

  end subroutine mxm5_3comp_3dmat_singleB


  end subroutine compute_element_stress_cm_Dev


!
!--------------------------------------------------------------------------------------------
!

  subroutine compute_element_stress_cm_noDev(ispec,nglob,nspec, &
                                          displ, &
                                          hprime_xxl,hprime_yyl,hprime_zzl, &
                                          ibool, &
                                          xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz, &
                                          sigma_xx,sigma_yy,sigma_zz, &
                                          sigma_xy,sigma_xz,sigma_yz, &
                                          sigma_yx,sigma_zx,sigma_zy)

! computes stress for single element in the crust/mantle

  use constants
  use specfem_par
  use specfem_par_crustmantle

  implicit none

  ! element id
  integer,intent(in) :: ispec

  integer,intent(in) :: NSPEC,NGLOB

  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB),intent(in) :: displ

  ! array with derivatives of Lagrange polynomials and precalculated products
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLX),intent(in) :: hprime_xxl
  real(kind=CUSTOM_REAL), dimension(NGLLY,NGLLY),intent(in) :: hprime_yyl
  real(kind=CUSTOM_REAL), dimension(NGLLZ,NGLLZ),intent(in) :: hprime_zzl

  ! arrays with mesh parameters per slice
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC),intent(in) :: ibool
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC),intent(in) :: &
    xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz
	
  ! Components of the stress tensor
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ),intent(out) :: sigma_xx,sigma_yy,sigma_zz, &
       sigma_xy,sigma_xz,sigma_yz, &
       sigma_yx,sigma_zx,sigma_zy
  

  ! local parameters
  integer :: iglob
  integer :: i,j,k,l
  integer i_SLS

  real(kind=CUSTOM_REAL) tempx1l,tempx2l,tempx3l
  real(kind=CUSTOM_REAL) tempy1l,tempy2l,tempy3l
  real(kind=CUSTOM_REAL) tempz1l,tempz2l,tempz3l

  real(kind=CUSTOM_REAL) hp1,hp2,hp3
  real(kind=CUSTOM_REAL) xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl
  real(kind=CUSTOM_REAL) duxdxl,duxdyl,duxdzl,duydxl,duydyl,duydzl,duzdxl,duzdyl,duzdzl
  real(kind=CUSTOM_REAL) duxdxl_plus_duydyl,duxdxl_plus_duzdzl,duydyl_plus_duzdzl
  real(kind=CUSTOM_REAL) duxdyl_plus_duydxl,duzdxl_plus_duxdzl,duzdyl_plus_duydzl
  
  ! Components of the stress tensor at each i,j,k
  real(kind=CUSTOM_REAL) sigma_xxl,sigma_yyl,sigma_zzl,sigma_xyl,sigma_xzl,sigma_yzl, &
       sigma_yxl,sigma_zxl,sigma_zyl
  ! the 21 coefficients for an anisotropic medium in reduced notation
  real(kind=CUSTOM_REAL) c11,c22,c33,c44,c55,c66,c12,c13,c23,c14,c24,c34,c15,c25,c35,c45,c16,c26,c36,c46,c56
  real(kind=CUSTOM_REAL) lambdal,mul,lambdalplus2mul
  real(kind=CUSTOM_REAL) kappal
  real(kind=CUSTOM_REAL) R_xx_val,R_yy_val
  !real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,5) :: epsilondev_loc
  real(kind=CUSTOM_REAL) sx_l,sy_l,sz_l,gxl,gyl,gzl
  !real(kind=CUSTOM_REAL) templ
  

  do k = 1,NGLLZ
    do j = 1,NGLLY
      do i = 1,NGLLX
         
        tempx1l = 0._CUSTOM_REAL
        tempx2l = 0._CUSTOM_REAL
        tempx3l = 0._CUSTOM_REAL

        tempy1l = 0._CUSTOM_REAL
        tempy2l = 0._CUSTOM_REAL
        tempy3l = 0._CUSTOM_REAL

        tempz1l = 0._CUSTOM_REAL
        tempz2l = 0._CUSTOM_REAL
        tempz3l = 0._CUSTOM_REAL

        do l = 1,NGLLX
          hp1 = hprime_xxl(i,l)
          iglob = ibool(l,j,k,ispec)
          tempx1l = tempx1l + displ(1,iglob)*hp1
          tempy1l = tempy1l + displ(2,iglob)*hp1
          tempz1l = tempz1l + displ(3,iglob)*hp1
!!! can merge these loops because NGLLX = NGLLY = NGLLZ          enddo

!!! can merge these loops because NGLLX = NGLLY = NGLLZ          do l = 1,NGLLY
          hp2 = hprime_yyl(j,l)
          iglob = ibool(i,l,k,ispec)
          tempx2l = tempx2l + displ(1,iglob)*hp2
          tempy2l = tempy2l + displ(2,iglob)*hp2
          tempz2l = tempz2l + displ(3,iglob)*hp2
!!! can merge these loops because NGLLX = NGLLY = NGLLZ          enddo

!!! can merge these loops because NGLLX = NGLLY = NGLLZ          do l = 1,NGLLZ
          hp3 = hprime_zzl(k,l)
          iglob = ibool(i,j,l,ispec)
          tempx3l = tempx3l + displ(1,iglob)*hp3
          tempy3l = tempy3l + displ(2,iglob)*hp3
          tempz3l = tempz3l + displ(3,iglob)*hp3
        enddo

        ! get derivatives of ux, uy and uz with respect to x, y and z
        xixl = xix(i,j,k,ispec)
        xiyl = xiy(i,j,k,ispec)
        xizl = xiz(i,j,k,ispec)
        etaxl = etax(i,j,k,ispec)
        etayl = etay(i,j,k,ispec)
        etazl = etaz(i,j,k,ispec)
        gammaxl = gammax(i,j,k,ispec)
        gammayl = gammay(i,j,k,ispec)
        gammazl = gammaz(i,j,k,ispec)

        ! compute the Jacobian
        !jacobianl = 1._CUSTOM_REAL / (xixl*(etayl*gammazl-etazl*gammayl) &
        !              - xiyl*(etaxl*gammazl-etazl*gammaxl) &
        !              + xizl*(etaxl*gammayl-etayl*gammaxl))

        duxdxl = xixl*tempx1l + etaxl*tempx2l + gammaxl*tempx3l
        duxdyl = xiyl*tempx1l + etayl*tempx2l + gammayl*tempx3l
        duxdzl = xizl*tempx1l + etazl*tempx2l + gammazl*tempx3l

        duydxl = xixl*tempy1l + etaxl*tempy2l + gammaxl*tempy3l
        duydyl = xiyl*tempy1l + etayl*tempy2l + gammayl*tempy3l
        duydzl = xizl*tempy1l + etazl*tempy2l + gammazl*tempy3l

        duzdxl = xixl*tempz1l + etaxl*tempz2l + gammaxl*tempz3l
        duzdyl = xiyl*tempz1l + etayl*tempz2l + gammayl*tempz3l
        duzdzl = xizl*tempz1l + etazl*tempz2l + gammazl*tempz3l

        ! precompute some sums to save CPU time
        duxdxl_plus_duydyl = duxdxl + duydyl
        duxdxl_plus_duzdzl = duxdxl + duzdzl
        duydyl_plus_duzdzl = duydyl + duzdzl
        duxdyl_plus_duydxl = duxdyl + duydxl
        duzdxl_plus_duxdzl = duzdxl + duxdzl
        duzdyl_plus_duydzl = duzdyl + duydzl
		
		
        if (ANISOTROPIC_3D_MANTLE_VAL) then

          c11 = c11store_crust_mantle(i,j,k,ispec)
          c12 = c12store_crust_mantle(i,j,k,ispec)
          c13 = c13store_crust_mantle(i,j,k,ispec)
          c14 = c14store_crust_mantle(i,j,k,ispec)
          c15 = c15store_crust_mantle(i,j,k,ispec)
          c16 = c16store_crust_mantle(i,j,k,ispec)
          c22 = c22store_crust_mantle(i,j,k,ispec)
          c23 = c23store_crust_mantle(i,j,k,ispec)
          c24 = c24store_crust_mantle(i,j,k,ispec)
          c25 = c25store_crust_mantle(i,j,k,ispec)
          c26 = c26store_crust_mantle(i,j,k,ispec)
          c33 = c33store_crust_mantle(i,j,k,ispec)
          c34 = c34store_crust_mantle(i,j,k,ispec)
          c35 = c35store_crust_mantle(i,j,k,ispec)
          c36 = c36store_crust_mantle(i,j,k,ispec)
          c44 = c44store_crust_mantle(i,j,k,ispec)
          c45 = c45store_crust_mantle(i,j,k,ispec)
          c46 = c46store_crust_mantle(i,j,k,ispec)
          c55 = c55store_crust_mantle(i,j,k,ispec)
          c56 = c56store_crust_mantle(i,j,k,ispec)
          c66 = c66store_crust_mantle(i,j,k,ispec)
          
          sigma_xxl = c11*duxdxl + c16*duxdyl_plus_duydxl + c12*duydyl + &
               c15*duzdxl_plus_duxdzl + c14*duzdyl_plus_duydzl + c13*duzdzl

          sigma_yyl = c12*duxdxl + c26*duxdyl_plus_duydxl + c22*duydyl + &
               c25*duzdxl_plus_duxdzl + c24*duzdyl_plus_duydzl + c23*duzdzl

          sigma_zzl = c13*duxdxl + c36*duxdyl_plus_duydxl + c23*duydyl + &
               c35*duzdxl_plus_duxdzl + c34*duzdyl_plus_duydzl + c33*duzdzl
          
          sigma_xyl = c16*duxdxl + c66*duxdyl_plus_duydxl + c26*duydyl + &
               c56*duzdxl_plus_duxdzl + c46*duzdyl_plus_duydzl + c36*duzdzl

          sigma_xzl = c15*duxdxl + c56*duxdyl_plus_duydxl + c25*duydyl + &
               c55*duzdxl_plus_duxdzl + c45*duzdyl_plus_duydzl + c35*duzdzl

          sigma_yzl = c14*duxdxl + c46*duxdyl_plus_duydxl + c24*duydyl + &
               c45*duzdxl_plus_duxdzl + c44*duzdyl_plus_duydzl + c34*duzdzl

        else

          ! do not use transverse isotropy except if element is between d220 and Moho
          if (.not. ispec_is_tiso_crust_mantle(ispec)) then
             
            ! isotropic element
             
            ! layer with no transverse isotropy, use kappav and muv
            kappal = kappavstore_crust_mantle(i,j,k,ispec)
            mul = muvstore_crust_mantle(i,j,k,ispec)
            
            ! use unrelaxed parameters if attenuation
            ! already done in prepare_timerun...
            !if (ATTENUATION_VAL) mul = mul * one_minus_sum_beta_use
            
            lambdalplus2mul = kappal + FOUR_THIRDS * mul
            lambdal = lambdalplus2mul - 2.*mul

            ! compute stress sigma
            
            sigma_xxl = lambdalplus2mul*duxdxl + lambdal*duydyl_plus_duzdzl
            sigma_yyl = lambdalplus2mul*duydyl + lambdal*duxdxl_plus_duzdzl
            sigma_zzl = lambdalplus2mul*duzdzl + lambdal*duxdxl_plus_duydyl
            
            sigma_xyl = mul*duxdyl_plus_duydxl
            sigma_xzl = mul*duzdxl_plus_duxdzl
            sigma_yzl = mul*duzdyl_plus_duydzl

          else

            ! transverse isotropic element
            c11 = c11store_crust_mantle(i,j,k,ispec)
            c12 = c12store_crust_mantle(i,j,k,ispec)
            c13 = c13store_crust_mantle(i,j,k,ispec)
            c14 = c14store_crust_mantle(i,j,k,ispec)
            c15 = c15store_crust_mantle(i,j,k,ispec)
            c16 = c16store_crust_mantle(i,j,k,ispec)
            c22 = c22store_crust_mantle(i,j,k,ispec)
            c23 = c23store_crust_mantle(i,j,k,ispec)
            c24 = c24store_crust_mantle(i,j,k,ispec)
            c25 = c25store_crust_mantle(i,j,k,ispec)
            c26 = c26store_crust_mantle(i,j,k,ispec)
            c33 = c33store_crust_mantle(i,j,k,ispec)
            c34 = c34store_crust_mantle(i,j,k,ispec)
            c35 = c35store_crust_mantle(i,j,k,ispec)
            c36 = c36store_crust_mantle(i,j,k,ispec)
            c44 = c44store_crust_mantle(i,j,k,ispec)
            c45 = c45store_crust_mantle(i,j,k,ispec)
            c46 = c46store_crust_mantle(i,j,k,ispec)
            c55 = c55store_crust_mantle(i,j,k,ispec)
            c56 = c56store_crust_mantle(i,j,k,ispec)
            c66 = c66store_crust_mantle(i,j,k,ispec)
            
            ! general expression of stress tensor for full Cijkl with 21 coefficients
            sigma_xxl = c11*duxdxl + c16*duxdyl_plus_duydxl + c12*duydyl + &
                 c15*duzdxl_plus_duxdzl + c14*duzdyl_plus_duydzl + c13*duzdzl
            
            sigma_yyl = c12*duxdxl + c26*duxdyl_plus_duydxl + c22*duydyl + &
                 c25*duzdxl_plus_duxdzl + c24*duzdyl_plus_duydzl + c23*duzdzl

            sigma_zzl = c13*duxdxl + c36*duxdyl_plus_duydxl + c23*duydyl + &
                 c35*duzdxl_plus_duxdzl + c34*duzdyl_plus_duydzl + c33*duzdzl

            sigma_xyl = c16*duxdxl + c66*duxdyl_plus_duydxl + c26*duydyl + &
                 c56*duzdxl_plus_duxdzl + c46*duzdyl_plus_duydzl + c36*duzdzl
            
            sigma_xzl = c15*duxdxl + c56*duxdyl_plus_duydxl + c25*duydyl + &
                 c55*duzdxl_plus_duxdzl + c45*duzdyl_plus_duydzl + c35*duzdzl

            sigma_yzl = c14*duxdxl + c46*duxdyl_plus_duydxl + c24*duydyl + &
                 c45*duzdxl_plus_duxdzl + c44*duzdyl_plus_duydzl + c34*duzdzl
            
          endif

        endif   ! end of test whether isotropic or anisotropic element

        ! subtract memory variables if attenuation
        if (ATTENUATION_VAL .and. .not. PARTIAL_PHYS_DISPERSION_ONLY_VAL) then
           do i_SLS = 1,N_SLS
              R_xx_val = R_xx_crust_mantle(i,j,k,i_SLS,ispec)
              R_yy_val = R_yy_crust_mantle(i,j,k,i_SLS,ispec)
              sigma_xxl = sigma_xxl - R_xx_val
              sigma_yyl = sigma_yyl - R_yy_val
              sigma_zzl = sigma_zzl + R_xx_val + R_yy_val
              sigma_xyl = sigma_xyl - R_xy_crust_mantle(i,j,k,i_SLS,ispec)
              sigma_xzl = sigma_xzl - R_xz_crust_mantle(i,j,k,i_SLS,ispec)
              sigma_yzl = sigma_yzl - R_yz_crust_mantle(i,j,k,i_SLS,ispec)
           enddo
        endif

        ! define symmetric components of sigma for gravity
        sigma_yxl = sigma_xyl
        sigma_zxl = sigma_xzl
        sigma_zyl = sigma_yzl
        
        ! compute non-symmetric terms for gravity
        if (GRAVITY_VAL) then
           ! use mesh coordinates to get theta and phi
           ! x y and z contain r theta and phi
           iglob = ibool(i,j,k,ispec)
           
           ! Cartesian components of the gravitational acceleration
           gxl = gravity_pre_store_crust_mantle(1,iglob) ! minus_g*sin_theta*cos_phi * rho
           gyl = gravity_pre_store_crust_mantle(2,iglob) ! minus_g*sin_theta*sin_phi * rho
           gzl = gravity_pre_store_crust_mantle(3,iglob) ! minus_g*cos_theta * rho
           
           ! Cartesian components of gradient of gravitational acceleration
           ! get displacement and multiply by density to compute G tensor
           sx_l = displ_crust_mantle(1,iglob)
           sy_l = displ_crust_mantle(2,iglob)
           sz_l = displ_crust_mantle(3,iglob)
           
           ! compute G tensor from s . g and add to sigma (not symmetric)
           sigma_xxl = sigma_xxl + sy_l * gyl + sz_l * gzl
           sigma_yyl = sigma_yyl + sx_l * gxl + sz_l * gzl
           sigma_zzl = sigma_zzl + sx_l * gxl + sy_l * gyl
           
           sigma_xyl = sigma_xyl - sx_l * gyl
           sigma_yxl = sigma_yxl - sy_l * gxl
           
           sigma_xzl = sigma_xzl - sx_l * gzl
           sigma_zxl = sigma_zxl - sz_l * gxl
           
           sigma_yzl = sigma_yzl - sy_l * gzl
           sigma_zyl = sigma_zyl - sz_l * gyl
           
        endif  ! end of section with gravity terms

        ! Populate the stress tensor components
        sigma_xx(i,j,k)=sigma_xxl
        sigma_xy(i,j,k)=sigma_xyl
        sigma_xz(i,j,k)=sigma_xzl
        sigma_yx(i,j,k)=sigma_yxl
        sigma_yy(i,j,k)=sigma_yyl
        sigma_yz(i,j,k)=sigma_yzl
        sigma_zx(i,j,k)=sigma_zxl
        sigma_zy(i,j,k)=sigma_zyl
        sigma_zz(i,j,k)=sigma_zzl

      enddo ! NGLLX
    enddo ! NGLLY
  enddo ! NGLLZ
  

  end subroutine compute_element_stress_cm_noDev
  
!
!-------------------------------------------------------------------------------------------------
!

  subroutine max_norms_stress_cm(maxnormstress_cm)
  
  use specfem_par
  use specfem_par_crustmantle
  
  ! Norm of stress
  real(kind=CUSTOM_REAL) :: nstress_cm
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE), &
       intent(inout) :: maxnormstress_cm
	
  ! Components of the stress tensor
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: sigma_xx,sigma_yy,sigma_zz, &
		sigma_xy,sigma_xz,sigma_yz, &
		sigma_yx,sigma_zx,sigma_zy
							
  ! Scaling of the stress tensor to bars
  real(kind=CUSTOM_REAL) :: scaleval,scale_bar
  
  integer :: ispec,i,j,k
  
  ! Define scaling factor for stress
  scaleval = dsqrt(PI*GRAV*RHOAV)
  scale_bar = (RHOAV/1000.d0)*((R_PLANET*scaleval/1000.d0)**2)*10000.d0

  ! Compute the stress in the crust/mantle
  do ispec = 1, NSPEC_CRUST_MANTLE  
    if (USE_DEVILLE_PRODUCTS_VAL) then
      call compute_element_stress_cm_Dev(ispec,NGLOB_CRUST_MANTLE,NSPEC_CRUST_MANTLE, &
                                        displ_crust_mantle,ibool_crust_mantle, &
                                        hprime_xx,hprime_xxT, &
                                        deriv_mapping_crust_mantle, &
                                        sigma_xx,sigma_yy,sigma_zz, &
                                        sigma_xy,sigma_xz,sigma_yz, &
                                        sigma_yx,sigma_zx,sigma_zy)

    else
      call compute_element_stress_cm_noDev(ispec,NGLOB_CRUST_MANTLE,NSPEC_CRUST_MANTLE, &
                                          displ_crust_mantle, &
                                          hprime_xx,hprime_yy,hprime_zz,ibool_crust_mantle, &
                                          xix_crust_mantle,xiy_crust_mantle,xiz_crust_mantle, &
                                          etax_crust_mantle,etay_crust_mantle,etaz_crust_mantle, &
                                          gammax_crust_mantle,gammay_crust_mantle,gammaz_crust_mantle, &
                                          sigma_xx,sigma_yy,sigma_zz, &
                                          sigma_xy,sigma_xz,sigma_yz, &
                                          sigma_yx,sigma_zx,sigma_zy)
    endif
  
    ! compute the norm of stress, and track the maximum norm across all timesteps
    do k=1, NGLLZ
      do j=1, NGLLY
        do i=1, NGLLX
          nstress_cm=scale_bar*sqrt(sigma_xx(i,j,k)**2 & 
			  	  + sigma_yy(i,j,k)**2 &
				  + sigma_zz(i,j,k)**2 &
				  + sigma_xy(i,j,k)**2 &
				  + sigma_xz(i,j,k)**2 &
				  + sigma_yz(i,j,k)**2 & 
				  + sigma_xy(i,j,k)**2 &
				  + sigma_xz(i,j,k)**2 &
				  + sigma_yz(i,j,k)**2)
		  
          maxnormstress_cm(i,j,k,ispec)=max(nstress_cm,maxnormstress_cm(i,j,k,ispec))
        enddo
      enddo
    enddo
  enddo
  
  end subroutine max_norms_stress_cm  

!
!--------------------------------------------------------------------------------------------
!

  subroutine write_bin_ndispvel(maxnormdisp_cm,maxnormvel_cm,maxnormdisp_ic, &
				maxnormvel_ic,maxnormdisp_oc,maxnormvel_oc)
  
  use constants_solver							
  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore
  use specfem_par_outercore
  use specfem_par_movie

  implicit none
   
  ! input arguments, with the maximum norms of displacement and velocity throughout the interior
  real(kind=CUSTOM_REAL), dimension(NGLOB_CRUST_MANTLE), intent(in) :: maxnormdisp_cm,maxnormvel_cm
  real(kind=CUSTOM_REAL), dimension(NGLOB_OUTER_CORE), intent(in) :: maxnormdisp_oc,maxnormvel_oc 
  real(kind=CUSTOM_REAL), dimension(NGLOB_INNER_CORE), intent(in) :: maxnormdisp_ic,maxnormvel_ic
  
  ! bin file names
  character(len=MAX_STRING_LEN) :: dispfile_cm,dispfile_ic,dispfile_oc,velfile_cm,velfile_ic,velfile_oc
  
  integer :: ispec,iglob,i,j,k,ier
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: tmp_data


  ! Crust/mantle
  ! Use xcombine_vol_data_vtk to convert the resulting .bin files to VTK
  allocate(tmp_data(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE),stat=ier)
  if (ier /= 0 ) call exit_MPI(myrank,'Error allocating temporary array tmp_data')
  do ispec = 1, NSPEC_CRUST_MANTLE
    do k = 1, NGLLZ
      do j = 1, NGLLY
        do i = 1, NGLLX
          iglob = ibool_crust_mantle(i,j,k,ispec)
          ! norm
          tmp_data(i,j,k,ispec) = maxnormdisp_cm(iglob) 
        enddo
      enddo
    enddo
  enddo
  write(dispfile_cm, '(a,i6.6,a)') 'DATABASES_MPI/proc',myrank,'_reg1_maxnormdispl.bin'
  open(unit=IOUT,file=trim(dispfile_cm),status='unknown',form='unformatted',iostat=ier)
  if (ier /= 0 ) call exit_MPI(myrank,'Error opening file '//trim(dispfile_cm))
  write(IOUT) tmp_data
  close(IOUT)
  deallocate(tmp_data)
  !
  ! Velocity
  allocate(tmp_data(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE),stat=ier)
  if (ier /= 0 ) call exit_MPI(myrank,'Error allocating temporary array tmp_data')
  do ispec = 1, NSPEC_CRUST_MANTLE
    do k = 1, NGLLZ
      do j = 1, NGLLY
        do i = 1, NGLLX
          iglob = ibool_crust_mantle(i,j,k,ispec)
          ! norm
          tmp_data(i,j,k,ispec) = maxnormvel_cm(iglob) 
        enddo
      enddo
    enddo
  enddo
  write(velfile_cm, '(a,i6.6,a)') 'DATABASES_MPI/proc',myrank,'_reg1_maxnormveloc.bin'
  open(unit=IOUT,file=trim(velfile_cm),status='unknown',form='unformatted',iostat=ier)
  if (ier /= 0 ) call exit_MPI(myrank,'Error opening file '//trim(velfile_cm))
  write(IOUT) tmp_data
  close(IOUT)
  deallocate(tmp_data)
  
  
  ! Outer Core
  ! Potential and first derivative of potential, not displacement or velocity
  ! Displacement
  allocate(tmp_data(NGLLX,NGLLY,NGLLZ,NSPEC_OUTER_CORE),stat=ier)
  if (ier /= 0 ) call exit_MPI(myrank,'Error allocating temporary array tmp_data')
  do ispec = 1, NSPEC_OUTER_CORE
    do k = 1, NGLLZ
      do j = 1, NGLLY
        do i = 1, NGLLX
          iglob = ibool_outer_core(i,j,k,ispec)
          ! norm
          tmp_data(i,j,k,ispec) = maxnormdisp_oc(iglob) 
        enddo
      enddo
    enddo
  enddo
  write(dispfile_oc, '(a,i6.6,a)') 'DATABASES_MPI/proc',myrank,'_reg2_maxnormdispl.bin'
  open(unit=IOUT,file=trim(dispfile_oc),status='unknown',form='unformatted',iostat=ier)
  if (ier /= 0 ) call exit_MPI(myrank,'Error opening file '//trim(dispfile_oc))
  write(IOUT) tmp_data
  close(IOUT)
  deallocate(tmp_data)
  !
  ! Velocity
  allocate(tmp_data(NGLLX,NGLLY,NGLLZ,NSPEC_OUTER_CORE),stat=ier)
  if (ier /= 0 ) call exit_MPI(myrank,'Error allocating temporary array tmp_data')
  do ispec = 1, NSPEC_OUTER_CORE
    do k = 1, NGLLZ
      do j = 1, NGLLY
        do i = 1, NGLLX
          iglob = ibool_outer_core(i,j,k,ispec)
          ! norm
          tmp_data(i,j,k,ispec) = maxnormvel_oc(iglob) 
        enddo
      enddo
    enddo
  enddo
  write(velfile_oc, '(a,i6.6,a)') 'DATABASES_MPI/proc',myrank,'_reg2_maxnormveloc.bin'
  open(unit=IOUT,file=trim(velfile_oc),status='unknown',form='unformatted',iostat=ier)
  if (ier /= 0 ) call exit_MPI(myrank,'Error opening file '//trim(velfile_oc))
  write(IOUT) tmp_data
  close(IOUT)
  deallocate(tmp_data)
  
  
  ! Inner Core
  ! Displacement
  allocate(tmp_data(NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE),stat=ier)
  if (ier /= 0 ) call exit_MPI(myrank,'Error allocating temporary array tmp_data')
  do ispec = 1, NSPEC_INNER_CORE
    do k = 1, NGLLZ
      do j = 1, NGLLY
        do i = 1, NGLLX
          iglob = ibool_inner_core(i,j,k,ispec)
          ! norm
          tmp_data(i,j,k,ispec) = maxnormdisp_ic(iglob) 
        enddo
      enddo
    enddo
  enddo
  write(dispfile_ic, '(a,i6.6,a)') 'DATABASES_MPI/proc',myrank,'_reg3_maxnormdispl.bin'
  open(unit=IOUT,file=trim(dispfile_ic),status='unknown',form='unformatted',iostat=ier)
  if (ier /= 0 ) call exit_MPI(myrank,'Error opening file '//trim(dispfile_ic))
  write(IOUT) tmp_data
  close(IOUT)
  deallocate(tmp_data)
  !
  ! Velocity
  allocate(tmp_data(NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE),stat=ier)
  if (ier /= 0 ) call exit_MPI(myrank,'Error allocating temporary array tmp_data')
  do ispec = 1, NSPEC_INNER_CORE
    do k = 1, NGLLZ
      do j = 1, NGLLY
        do i = 1, NGLLX
          iglob = ibool_inner_core(i,j,k,ispec)
          ! norm
          tmp_data(i,j,k,ispec) = maxnormvel_ic(iglob) 
        enddo
      enddo
    enddo
  enddo
  write(velfile_ic, '(a,i6.6,a)') 'DATABASES_MPI/proc',myrank,'_reg3_maxnormveloc.bin'
  open(unit=IOUT,file=trim(velfile_ic),status='unknown',form='unformatted',iostat=ier)
  if (ier /= 0 ) call exit_MPI(myrank,'Error opening file '//trim(velfile_ic))
  write(IOUT) tmp_data
  close(IOUT)
  deallocate(tmp_data)	
								
  end subroutine write_bin_ndispvel		

!
!-------------------------------------------------------------------------------------------------
!

  subroutine write_bin_strain_cm(maxnormstrain_cm)
  
  use constants_solver
  use specfem_par
  use specfem_par_crustmantle
  
  implicit none
  
  ! Peak global norm of strain in the crust/mantle
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE),intent(inout) :: maxnormstrain_cm
  
  ! variables
  integer :: ier
  character(len=MAX_STRING_LEN) :: outputstrain
  
  ! Use xcombine_vol_data_vtk to convert the resulting .bin file to VTK
  write(outputstrain, '(a,i6.6,a)') 'DATABASES_MPI/proc',myrank,'_reg1_maxnormstrain.bin'
  open(unit=IOUT,file=trim(outputstrain),status='unknown',form='unformatted',iostat=ier)
  if (ier /= 0 ) call exit_MPI(myrank,'Error opening file '//trim(outputstrain))
  write(IOUT) maxnormstrain_cm
  close(IOUT)
  
  end subroutine write_bin_strain_cm


!
!-------------------------------------------------------------------------------------------------
!

  subroutine write_bin_stress_cm(maxnormstress_cm)
  
  use constants_solver
  use specfem_par
  use specfem_par_crustmantle
  
  implicit none
  
  ! Peak global norm of strain in the crust/mantle
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE),intent(inout) :: maxnormstress_cm
  
  ! variables
  integer :: ier
  character(len=MAX_STRING_LEN) :: outputstress
  
  ! Use xcombine_vol_data_vtk to convert the resulting .bin file to VTK
  write(outputstress, '(a,i6.6,a)') 'DATABASES_MPI/proc',myrank,'_reg1_maxnormstress.bin'
  open(unit=IOUT,file=trim(outputstress),status='unknown',form='unformatted',iostat=ier)
  if (ier /= 0 ) call exit_MPI(myrank,'Error opening file '//trim(outputstress))
  write(IOUT) maxnormstress_cm
  close(IOUT)
  
  
  end subroutine write_bin_stress_cm

