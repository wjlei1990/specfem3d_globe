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
  ! Temporary array to write the bin files
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
  
  ! Peak norm of strain in the crust/mantle
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE),intent(inout) :: maxnormstrain_cm
  
  ! Local variables
  ! Peak global/assembled norm of strain in the crust/mantle
  real(kind=CUSTOM_REAL), dimension(NGLOB_CRUST_MANTLE) :: maxnormstrain_cm_glob
  integer :: ispec,iglob,i,j,k,ier
  ! Temporary array to write the bin files
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: tmp_data
  ! bin file name
  character(len=MAX_STRING_LEN) :: outputstrain
  ! bit-mask
  logical, dimension(NGLOB_CRUST_MANTLE) :: mask_ibool
  
  ! Initialize global strain
  maxnormstrain_cm_glob(:)=0  
  ! Initialize the bit-mask
  mask_ibool(:) = .false.
  
  ! Local to global mapping---ignore
  do ispec = 1, NSPEC_CRUST_MANTLE
    do k = 1, NGLLZ
      do j = 1, NGLLY
        do i = 1, NGLLX
          iglob = ibool_crust_mantle(i,j,k,ispec)
          if (.not. mask_ibool(iglob)) then
             mask_ibool(iglob) = .true.
             maxnormstrain_cm_glob(iglob)=maxnormstrain_cm(i,j,k,ispec)
          endif
        enddo
      enddo
    enddo
  enddo
  
  ! Use xcombine_vol_data_vtk to convert the resulting .bin file to VTK
  allocate(tmp_data(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE),stat=ier)
  if (ier /= 0 ) call exit_MPI(myrank,'Error allocating temporary array tmp_data')
  do ispec = 1, NSPEC_CRUST_MANTLE
    do k = 1, NGLLZ
      do j = 1, NGLLY
        do i = 1, NGLLX
          iglob = ibool_crust_mantle(i,j,k,ispec)
          ! local to global mapping in combine_vol_data
          tmp_data(i,j,k,ispec) = maxnormstrain_cm(i,j,k,ispec)
        enddo
      enddo
    enddo
  enddo
  write(outputstrain, '(a,i6.6,a)') 'DATABASES_MPI/proc',myrank,'_reg1_maxnormstrain.bin'
  open(unit=IOUT,file=trim(outputstrain),status='unknown',form='unformatted',iostat=ier)
  if (ier /= 0 ) call exit_MPI(myrank,'Error opening file '//trim(outputstrain))
  write(IOUT) tmp_data
  close(IOUT)
  deallocate(tmp_data)

  
  end subroutine write_bin_strain_cm


!
!-------------------------------------------------------------------------------------------------
!

  subroutine write_bin_stress_cm(maxnormstress_cm)
  
  use constants_solver
  use specfem_par
  use specfem_par_crustmantle
  
  implicit none
  
  ! Peak norm of stress in the crust/mantle
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE),intent(inout) :: maxnormstress_cm
  
  ! Local variables
  ! Peak global/assembled norm of stress in the crust/mantle
  real(kind=CUSTOM_REAL), dimension(NGLOB_CRUST_MANTLE) :: maxnormstress_cm_glob  
  integer :: ispec,iglob,i,j,k,ier
  ! Temporary array to write the bin files
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: tmp_data
  ! bin file name
  character(len=MAX_STRING_LEN) :: outputstress
  ! bit-mask
  logical, dimension(NGLOB_CRUST_MANTLE) :: mask_ibool
  
  ! Initialize global stress
  maxnormstress_cm_glob(:)=0
  ! Initialize the bit-mask
  mask_ibool(:) = .false.
  
  ! Local to global mapping---ignore
  do ispec = 1, NSPEC_CRUST_MANTLE
    do k = 1, NGLLZ
      do j = 1, NGLLY
        do i = 1, NGLLX
          iglob = ibool_crust_mantle(i,j,k,ispec)
          if (.not. mask_ibool(iglob)) then
             mask_ibool(iglob) = .true.
             maxnormstress_cm_glob(iglob)=maxnormstress_cm(i,j,k,ispec)
          endif
        enddo
      enddo
    enddo
  enddo
  
  ! Use xcombine_vol_data_vtk to convert the resulting .bin file to VTK
  allocate(tmp_data(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE),stat=ier)
  if (ier /= 0 ) call exit_MPI(myrank,'Error allocating temporary array tmp_data')
  do ispec = 1, NSPEC_CRUST_MANTLE
    do k = 1, NGLLZ
      do j = 1, NGLLY
        do i = 1, NGLLX
          iglob = ibool_crust_mantle(i,j,k,ispec)
          ! local to global mapping in combine_vol_data
          tmp_data(i,j,k,ispec) = maxnormstress_cm(i,j,k,ispec)
        enddo
      enddo
    enddo
  enddo
  write(outputstress, '(a,i6.6,a)') 'DATABASES_MPI/proc',myrank,'_reg1_maxnormstress.bin'
  open(unit=IOUT,file=trim(outputstress),status='unknown',form='unformatted',iostat=ier)
  if (ier /= 0 ) call exit_MPI(myrank,'Error opening file '//trim(outputstress))
  write(IOUT) tmp_data
  close(IOUT)
  deallocate(tmp_data)
  
  
  end subroutine write_bin_stress_cm
