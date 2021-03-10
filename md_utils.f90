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

            function comp_normal_energy(Natoms,vel,Base) result(NM_energies)
                  implicit none
                  integer :: Natoms
                  real*8 :: Base(3*Natoms,3*Natoms),NM_energies(3*Natoms)
                  real*8 :: vel(3,Natoms),mw_vel(3,Natoms),mw_vel_unpacked(3*Natoms)
                  real*8 :: nm_vel_unpacked(3*Natoms)
                  
                  mw_vel(1,:) = vel(1,:)*sqrt(M)
                  mw_vel(2,:) = vel(2,:)*sqrt(M)
                  mw_vel(3,:) = vel(3,:)*sqrt(M)

                  mw_vel_unpacked = unpack_coords(mw_vel)

                  nm_vel_unpacked = matmul(transpose(Base),mw_vel_unpacked)

                  NM_energies = 0.5d0*nm_vel_unpacked**2
      end function comp_normal_energy

      end module md_utils