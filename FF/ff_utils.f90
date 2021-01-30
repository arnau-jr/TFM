      module ff_utils
            use tools
            implicit none

            contains

            real function comp_bond_energy(Nbonds,bond_vals,req,k) result(E)
                  implicit none
                  integer :: Nbonds
                  real :: req(Nbonds),bond_vals(Nbonds),k(Nbonds)
                  integer :: i

                  do i=1,Nbonds
                        E = E + k(i)*(bond_vals(i)-req(i))**2
                  enddo
            end function comp_bond_energy


      end module ff_utils