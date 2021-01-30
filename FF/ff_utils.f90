      module ff_utils
            use tools
            implicit none

            contains

            real function comp_bond_energy(Nbonds,bond_vals,req,k) result(E)
                  implicit none
                  integer :: Nbonds
                  real :: bond_vals(Nbonds),req(Nbonds),k(Nbonds)
                  integer :: i
                  E = 0.
                  do i=1,Nbonds
                        E = E + k(i)*(bond_vals(i)-req(i))**2
                  enddo
      end function comp_bond_energy

            real function comp_angle_energy(Nangles,angle_vals,aeq,k) result(E)
                  implicit none
                  integer :: Nangles
                  real :: angle_vals(Nangles),aeq(Nangles),k(Nangles)
                  integer :: i
                  E = 0.
                  do i=1,Nangles
                        E = E + k(i)*(angle_vals(i)-aeq(i))**2
                  enddo
      end function comp_angle_energy

            real function comp_torsion_energy(Ntorsions,torsion_vals,An,n,delta) result(E)
                  implicit none
                  integer :: Ntorsions
                  real :: torsion_vals(Ntorsions),An(Ntorsions),n(Ntorsions),delta(Ntorsions)
                  integer :: i
                  E = 0.
                  do i=1,Ntorsions
                        E = E + An(i)*cos(n(i)*torsion_vals(i) - delta(i))
                  enddo
      end function comp_torsion_energy

            real function comp_energy(Nbonds,Nangles,Ntorsions,bond_vals,&
            angle_vals,torsion_vals,req,kb,aeq,ka,An,n,delta) result (E)
                  implicit none
                  integer :: Nbonds,Nangles,Ntorsions
                  real :: bond_vals(Nbonds),angle_vals(Nangles),torsion_vals(Ntorsions)
                  real :: req(Nbonds),kb(Nbonds)
                  real :: aeq(Nangles),ka(Nangles)
                  real :: An(Nangles),n(Nangles),delta(Nangles)
                  E = 0.
                  E = E + comp_bond_energy(Nbonds,bond_vals,req,kb)
                  E = E + comp_angle_energy(Nangles,angle_vals,aeq,ka)
                  E = E + comp_torsion_energy(Ntorsions,torsion_vals,An,n,delta)
      end function comp_energy


      end module ff_utils