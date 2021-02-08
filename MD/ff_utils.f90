      module ff_utils
            use tools
            implicit none
            real*8,allocatable :: req(:),kb(:)
            real*8,allocatable :: aeq(:),ka(:)
            real*8,allocatable :: An(:),n(:),delta(:) 

            contains

            subroutine get_param(port,Nbonds,Nangles,Ntorsions)
                  implicit none
                  integer :: port, Nbonds,Nangles,Ntorsions
                  character :: dummy*90,dummy2*90
                  integer :: i,dummy3

                  allocate(req(Nbonds),kb(Nbonds))
                  allocate(aeq(Nangles),ka(Nangles))
                  allocate(An(Ntorsions),n(Ntorsions),delta(Ntorsions))
            
                  do i=1,Nbonds
                        read(port,*)dummy,dummy2,kb(i),req(i)
                  enddo
                  read(port,*)
                  do i=1,Nangles
                        read(port,*)dummy,dummy2,ka(i),aeq(i)
                  enddo
                  read(port,*)
                  do i=1,Ntorsions
                        read(port,*)dummy,dummy2,dummy3,An(i),delta(i),n(i)
                  enddo
      end subroutine get_param
            real*8 function comp_bond_energy(r,req,k) result(E)
                  implicit none
                  real*8 :: r,req,k
                  E = k*(r-req)**2
      end function comp_bond_energy

            real*8 function comp_bonds_energy(Nbonds,bond_vals) result(E)
                  implicit none
                  integer :: Nbonds
                  real*8 :: bond_vals(Nbonds)
                  integer :: i
                  E = 0.
                  do i=1,Nbonds
                        E = E + kb(i)*(bond_vals(i)-req(i))**2
                  enddo
      end function comp_bonds_energy

                  real*8 function comp_angle_energy(a,aeq,k) result(E)
                  implicit none
                  real*8 :: a,aeq,k
                  E = k*((pi/180.d0)*(a-aeq))**2
      end function comp_angle_energy

            real*8 function comp_angles_energy(Nangles,angle_vals) result(E)
                  implicit none
                  integer :: Nangles
                  real*8 :: angle_vals(Nangles)
                  integer :: i
                  E = 0.
                  do i=1,Nangles
                        E = E + ka(i)*((pi/180.d0)*(angle_vals(i)-aeq(i)))**2
                  enddo
      end function comp_angles_energy

            real*8 function comp_torsion_energy(T,An,n,delta) result(E)
                  implicit none
                  real*8 :: T,An,n,delta
                  E = An*(1. + cos((pi/180.d0)*(n*T - delta)))
      end function comp_torsion_energy

            real*8 function comp_torsions_energy(Ntorsions,torsion_vals) result(E)
                  implicit none
                  integer :: Ntorsions
                  real*8 :: torsion_vals(Ntorsions)
                  integer :: i
                  E = 0.
                  do i=1,Ntorsions
                        E = E + An(i)*(1. + cos((pi/180.d0)*(n(i)*torsion_vals(i) - delta(i))))
                  enddo
      end function comp_torsions_energy

            real*8 function comp_energy(Nbonds,Nangles,Ntorsions,bond_vals,&
            angle_vals,torsion_vals) result (E)
                  implicit none
                  integer :: Nbonds,Nangles,Ntorsions
                  real*8 :: bond_vals(Nbonds),angle_vals(Nangles),torsion_vals(Ntorsions)
                  E = 0.
                  E = E + comp_bonds_energy(Nbonds,bond_vals)
                  E = E + comp_angles_energy(Nangles,angle_vals)
                  E = E + comp_torsions_energy(Ntorsions,torsion_vals)
      end function comp_energy

            function build_gradient(Natoms,xyz,Nbonds,Nangles,Ntorsions,&
            bond_pairs,angle_pairs,torsion_pairs) result(G)
                  implicit none
                  integer :: Natoms,Nbonds,Nangles,Ntorsions
                  real*8 :: xyz(3,Natoms),xyzmod(3,Natoms)
                  integer :: bond_pairs(Nbonds,2),angle_pairs(Nangles,3),torsion_pairs(Ntorsions,4)
                  real*8 :: bond_vals(Nbonds),angle_vals(Nangles),torsion_vals(Ntorsions)
                  real*8 :: G(3,Natoms)
                  real*8,parameter :: hi = 0.01
                  real*8 :: Vp,Vm
                  integer :: a,b,p,q,i,j
                  integer :: bond,angle,torsion,k,l,m


                  G = 0.d0
                  do a=1,Natoms
                        do p=1,3
                              xyzmod = xyz
                              xyzmod(p,a) = xyz(p,a) + hi

                              bond_vals = recomp_bonds(Natoms,Nbonds,xyzmod,bond_pairs)
                              angle_vals = recomp_angles(Natoms,Nangles,xyzmod,angle_pairs)
                              torsion_vals = recomp_torsions(Natoms,Ntorsions,xyzmod,torsion_pairs)

                              Vp = comp_energy(Nbonds,Nangles,Ntorsions,bond_vals,&
                              angle_vals,torsion_vals)

                              xyzmod = xyz
                              xyzmod(p,a) = xyz(p,a) - hi

                              bond_vals = recomp_bonds(Natoms,Nbonds,xyzmod,bond_pairs)
                              angle_vals = recomp_angles(Natoms,Nangles,xyzmod,angle_pairs)
                              torsion_vals = recomp_torsions(Natoms,Ntorsions,xyzmod,torsion_pairs)

                              Vm = comp_energy(Nbonds,Nangles,Ntorsions,bond_vals,&
                              angle_vals,torsion_vals)

                              G(p,a) = (Vp-Vm)/(2.d0*hi)
                        enddo

                  enddo
      end function build_gradient

            function build_hessian(Natoms,xyz,Nbonds,Nangles,Ntorsions,&
            bond_pairs,angle_pairs,torsion_pairs) result(H)
                  implicit none
                  integer :: Natoms,Nbonds,Nangles,Ntorsions
                  real*8 :: xyz(3,Natoms),xyzmod(3,Natoms)
                  integer :: bond_pairs(Nbonds,2),angle_pairs(Nangles,3),torsion_pairs(Ntorsions,4)
                  real*8 :: bond_vals(Nbonds),angle_vals(Nangles),torsion_vals(Ntorsions)
                  real*8 :: H(3*Natoms,3*Natoms)
                  real*8,parameter :: hi = 0.01,hj = 0.01
                  real*8 :: disppi(3),dispmi(3),disppj(3),dispmj(3)
                  real*8 :: V,Vpp,Vmm,Vpm,Vmp
                  integer :: a,b,p,q,i,j
                  integer :: bond,angle,torsion,k,l,m

                  H = 0.d0
                  do a=1,Natoms
                        do b=1,Natoms
                              do p=0,2
                              do q=0,2
                              i = 3*a - p
                              j = 3*b - q
                              
                              if(i/=j) then
                              xyzmod = xyz
                              xyzmod(3-p,a) = xyz(3-p,a) + hi
                              xyzmod(3-q,b) = xyz(3-q,b) + hj

                              bond_vals = recomp_bonds(Natoms,Nbonds,xyzmod,bond_pairs)
                              angle_vals = recomp_angles(Natoms,Nangles,xyzmod,angle_pairs)
                              torsion_vals = recomp_torsions(Natoms,Ntorsions,xyzmod,torsion_pairs)

                              Vpp = comp_energy(Nbonds,Nangles,Ntorsions,bond_vals,&
                              angle_vals,torsion_vals)

                              xyzmod = xyz
                              xyzmod(3-p,a) = xyz(3-p,a) - hi
                              xyzmod(3-q,b) = xyz(3-q,b) - hj

                              bond_vals = recomp_bonds(Natoms,Nbonds,xyzmod,bond_pairs)
                              angle_vals = recomp_angles(Natoms,Nangles,xyzmod,angle_pairs)
                              torsion_vals = recomp_torsions(Natoms,Ntorsions,xyzmod,torsion_pairs)

                              Vmm = comp_energy(Nbonds,Nangles,Ntorsions,bond_vals,&
                              angle_vals,torsion_vals)

                              xyzmod = xyz
                              xyzmod(3-p,a) = xyz(3-p,a) + hi
                              xyzmod(3-q,b) = xyz(3-q,b) - hj

                              bond_vals = recomp_bonds(Natoms,Nbonds,xyzmod,bond_pairs)
                              angle_vals = recomp_angles(Natoms,Nangles,xyzmod,angle_pairs)
                              torsion_vals = recomp_torsions(Natoms,Ntorsions,xyzmod,torsion_pairs)

                              Vpm = comp_energy(Nbonds,Nangles,Ntorsions,bond_vals,&
                              angle_vals,torsion_vals)

                              xyzmod = xyz
                              xyzmod(3-p,a) = xyz(3-p,a) - hi
                              xyzmod(3-q,b) = xyz(3-q,b) + hj

                              bond_vals = recomp_bonds(Natoms,Nbonds,xyzmod,bond_pairs)
                              angle_vals = recomp_angles(Natoms,Nangles,xyzmod,angle_pairs)
                              torsion_vals = recomp_torsions(Natoms,Ntorsions,xyzmod,torsion_pairs)

                              Vmp = comp_energy(Nbonds,Nangles,Ntorsions,bond_vals,&
                              angle_vals,torsion_vals)

                              H(i,j) = (Vpp - Vpm - Vmp + Vmm)/(4.d0*hi*hj)
                              else

                              xyzmod = xyz
                              xyzmod(3-p,a) = xyz(3-p,a) + hi

                              bond_vals = recomp_bonds(Natoms,Nbonds,xyzmod,bond_pairs)
                              angle_vals = recomp_angles(Natoms,Nangles,xyzmod,angle_pairs)
                              torsion_vals = recomp_torsions(Natoms,Ntorsions,xyzmod,torsion_pairs)

                              Vpp = comp_energy(Nbonds,Nangles,Ntorsions,bond_vals,&
                              angle_vals,torsion_vals)

                              xyzmod = xyz
                              xyzmod(3-p,a) = xyz(3-p,a) - hi

                              bond_vals = recomp_bonds(Natoms,Nbonds,xyzmod,bond_pairs)
                              angle_vals = recomp_angles(Natoms,Nangles,xyzmod,angle_pairs)
                              torsion_vals = recomp_torsions(Natoms,Ntorsions,xyzmod,torsion_pairs)

                              Vmm = comp_energy(Nbonds,Nangles,Ntorsions,bond_vals,&
                              angle_vals,torsion_vals)

                              xyzmod = xyz

                              bond_vals = recomp_bonds(Natoms,Nbonds,xyzmod,bond_pairs)
                              angle_vals = recomp_angles(Natoms,Nangles,xyzmod,angle_pairs)
                              torsion_vals = recomp_torsions(Natoms,Ntorsions,xyzmod,torsion_pairs)

                              V = comp_energy(Nbonds,Nangles,Ntorsions,bond_vals,&
                              angle_vals,torsion_vals)

                              H(i,j) = (Vpp - 2.d0*V + Vmm)/(hi*hj)
                              
                              endif
                              enddo
                              enddo
                        enddo
                  enddo
      end function build_hessian

      end module ff_utils