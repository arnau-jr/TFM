      module md_utils
            use tools
            use ff_utils
            implicit none
            real*8 :: M_solvent

            contains

            function init_normalmodes_micro(Natoms,Base,T) result(vel_init)
                  implicit none
                  integer :: Natoms
                  real*8 :: Base(3*Natoms,3*Natoms),T
                  real*8 :: vel_init(3,Natoms)
                  real*8 :: unpacked_vel_init(3*Natoms)
                  integer :: i

                  unpacked_vel_init = 0.d0
                  do i=7,3*Natoms
                        ! unpacked_vel_init(i) = 2.d0*T
                        unpacked_vel_init(i) = sqrt(2.d0*T)
                  enddo

                  unpacked_vel_init = matmul(Base,unpacked_vel_init)
                  vel_init = pack_coords(unpacked_vel_init)

                  vel_init(1,:) = vel_init(1,:)/sqrt(M)
                  vel_init(2,:) = vel_init(2,:)/sqrt(M)
                  vel_init(3,:) = vel_init(3,:)/sqrt(M)
      end function init_normalmodes_micro

            function comp_kinetic_energy(vel,flag) result(K)
                  implicit none
                  real*8 :: vel(:,:),K
                  integer :: flag,i
                  if(flag==0) then
                        K = sum(0.5d0*M*sum(vel**2,1))
                  elseif(flag==1) then
                        K = sum(0.5d0*M_solvent*sum(vel**2,1))
                  endif
      end function comp_kinetic_energy

            function comp_normal_energy(Natoms,xyz_cm,xyz_eckart,vel,Base,freqs) result(NM_energies)
                  implicit none
                  integer :: Natoms
                  real*8  :: xyz_cm(3,Natoms),xyz_eckart(3,Natoms)
                  real*8  :: Base(3*Natoms,3*Natoms),freqs(3*Natoms),NM_energies(3*Natoms)
                  real*8  :: vel(3,Natoms),vel_nm(3*Natoms)
                  real*8  :: xyz_cm_nm(3*Natoms),xyz_eckart_nm(3*Natoms) 
                  integer :: i
                  
                  !Kinetic part
                  vel_nm = cart_to_normal(Natoms,vel,Base)

                  NM_energies = 0.5d0*vel_nm**2

                  !Potential part
                  xyz_cm_nm = cart_to_normal(Natoms,xyz_cm,Base)
                  xyz_eckart_nm = cart_to_normal(Natoms,xyz_eckart,Base)

                  NM_energies = NM_energies + 0.5d0*freqs*(xyz_cm_nm-xyz_eckart_nm)**2
      end function comp_normal_energy

      end module md_utils