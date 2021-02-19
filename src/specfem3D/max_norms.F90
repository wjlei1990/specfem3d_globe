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
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ),intent(inout) :: sigma_xx,sigma_yy,sigma_zz, &
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
	sigma_xzl = c15*duxdxl + c56*duxdyl_plus_duydxl + c25*duydyl + &
                           c55*duzdxl_plus_duxdzl + c45*duzdyl_plus_duydzl + c35*duzdzl
        sigma_yzl = c14*duxdxl + c46*duxdyl_plus_duydxl + c24*duydyl + &
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
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ),intent(inout) :: sigma_xx,sigma_yy,sigma_zz, &
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
  
  
