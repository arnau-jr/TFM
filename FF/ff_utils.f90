      module ff_utils
            use tools
            implicit none

            contains

            real function comp_bond_energy(r,req,k) result(E)
                  implicit none
                  real :: r,req,k
                  E = k*(r-req)**2
      end function comp_bond_energy

            real function comp_bonds_energy(Nbonds,bond_vals,req,k) result(E)
                  implicit none
                  integer :: Nbonds
                  real :: bond_vals(Nbonds),req(Nbonds),k(Nbonds)
                  integer :: i
                  E = 0.
                  do i=1,Nbonds
                        E = E + k(i)*(bond_vals(i)-req(i))**2
                  enddo
      end function comp_bonds_energy

                  real function comp_angle_energy(a,aeq,k) result(E)
                  implicit none
                  real :: a,aeq,k
                  E = k*(a-aeq)**2
      end function comp_angle_energy

            real function comp_angles_energy(Nangles,angle_vals,aeq,k) result(E)
                  implicit none
                  integer :: Nangles
                  real :: angle_vals(Nangles),aeq(Nangles),k(Nangles)
                  integer :: i
                  E = 0.
                  do i=1,Nangles
                        E = E + k(i)*(angle_vals(i)-aeq(i))**2
                  enddo
      end function comp_angles_energy

            real function comp_torsion_energy(T,An,n,delta) result(E)
                  implicit none
                  real :: T,An,n,delta
                  E = An*(1. + cos(n*T - delta))
      end function comp_torsion_energy

            real function comp_torsions_energy(Ntorsions,torsion_vals,An,n,delta) result(E)
                  implicit none
                  integer :: Ntorsions
                  real :: torsion_vals(Ntorsions),An(Ntorsions),n(Ntorsions),delta(Ntorsions)
                  integer :: i
                  E = 0.
                  do i=1,Ntorsions
                        E = E + An(i)*(1. + cos(n(i)*torsion_vals(i) - delta(i)))
                  enddo
      end function comp_torsions_energy

            real function comp_energy(Nbonds,Nangles,Ntorsions,bond_vals,&
            angle_vals,torsion_vals,req,kb,aeq,ka,An,n,delta) result (E)
                  implicit none
                  integer :: Nbonds,Nangles,Ntorsions
                  real :: bond_vals(Nbonds),angle_vals(Nangles),torsion_vals(Ntorsions)
                  real :: req(Nbonds),kb(Nbonds)
                  real :: aeq(Nangles),ka(Nangles)
                  real :: An(Nangles),n(Nangles),delta(Nangles)
                  E = 0.
                  E = E + comp_bonds_energy(Nbonds,bond_vals,req,kb)
                  E = E + comp_angles_energy(Nangles,angle_vals,aeq,ka)
                  E = E + comp_torsions_energy(Ntorsions,torsion_vals,An,n,delta)
      end function comp_energy


            function build_hessian(Natoms,xyz,Nbonds,Nangles,Ntorsions,&
            bond_pairs,angle_pairs,torsion_pairs,req,kb,aeq,ka,An,n,delta) result(H)
                  implicit none
                  integer :: Natoms,Nbonds,Nangles,Ntorsions
                  real :: xyz(3,Natoms)
                  integer :: bond_pairs(Nbonds,2),angle_pairs(Nangles,3),torsion_pairs(Ntorsions,4)
                  real :: req(Nbonds),kb(Nbonds)
                  real :: aeq(Nangles),ka(Nangles)
                  real :: An(Nangles),n(Nangles),delta(Nangles)
                  real*8 :: H(3*Natoms,3*Natoms)
                  real,parameter :: hi = 0.01,hj = 0.01
                  real :: disppi(3),dispmi(3),disppj(3),dispmj(3)
                  real :: Vpp,Vmm,Vpm,Vmp
                  integer :: a,b,p,q,i,j
                  integer :: bond,angle,torsion,k,l,m
                  real :: r12(3)

                  H = 0.d0
                  !TODO: TREAT DIAGONAL PROPERLY
                  do a=1,Natoms
                        do b=1,Natoms
                              do p=0,2
                              do q=0,2
                              i = 3*a - p
                              j = 3*b - q

                              disppi = xyz(:,a)
                              dispmi = xyz(:,a)
                              disppj = xyz(:,b)
                              dispmj = xyz(:,b)

                              disppi(p+1) = xyz(p+1,a) + hi
                              dispmi(p+1) = xyz(p+1,a) - hi
                              disppj(q+1) = xyz(q+1,b) + hj
                              dispmj(q+1) = xyz(q+1,b) - hj
                        !Check all bonds
                        do bond=1,Nbonds
                              !If both atoms are involved in the bond
                              if((bond_pairs(bond,1)==a .and. bond_pairs(bond,2)==b)&
                              .or. (bond_pairs(bond,2)==a .and. bond_pairs(bond,1)==b)) then
                                    Vpp = comp_bond_energy(get_dist(disppi,disppj),req(bond),kb(bond))
                                    Vmm = comp_bond_energy(get_dist(dispmi,dispmj),req(bond),kb(bond))
                                    Vpm = comp_bond_energy(get_dist(disppi,dispmj),req(bond),kb(bond))
                                    Vmp = comp_bond_energy(get_dist(dispmi,disppj),req(bond),kb(bond))
                              !If only a is involved
                              elseif(any(bond_pairs(bond,:)==a)) then
                                    if(bond_pairs(bond,1)==a) k=bond_pairs(bond,2) 
                                    if(bond_pairs(bond,2)==a) k=bond_pairs(bond,1)
                                    Vpp = comp_bond_energy(get_dist(disppi,xyz(:,k)),req(bond),kb(bond))
                                    Vmm = comp_bond_energy(get_dist(dispmi,xyz(:,k)),req(bond),kb(bond))
                                    Vpm = comp_bond_energy(get_dist(disppi,xyz(:,k)),req(bond),kb(bond))
                                    Vmp = comp_bond_energy(get_dist(dispmi,xyz(:,k)),req(bond),kb(bond))
                              !If only b is involved
                              elseif(any(bond_pairs(bond,:)==b)) then
                                    if(bond_pairs(bond,1)==b) k=bond_pairs(bond,2)
                                    if(bond_pairs(bond,2)==b) k=bond_pairs(bond,1)
                                    Vpp = comp_bond_energy(get_dist(disppj,xyz(:,k)),req(bond),kb(bond))
                                    Vmm = comp_bond_energy(get_dist(dispmj,xyz(:,k)),req(bond),kb(bond))
                                    Vpm = comp_bond_energy(get_dist(dispmj,xyz(:,k)),req(bond),kb(bond))
                                    Vmp = comp_bond_energy(get_dist(disppj,xyz(:,k)),req(bond),kb(bond))
                              else
                                    Vpp = 0.
                                    Vpm = 0.
                                    Vmp = 0.
                                    Vmm = 0.
                              endif
                              H(i,j) = H(i,j) + (Vpp - Vpm - Vmp + Vmm)/(4.d0*hi*hj)
                        enddo

                        !For all angles
                        do angle=1,Nangles
                              !If both a and b are involved
                              if(any(angle_pairs(angle,:)==a) .or. any(angle_pairs(angle,:)==b)) then
                                    !If a is an extreme
                                    if(angle_pairs(angle,1)==a .or. angle_pairs(angle,3)==a) then
                                          !If b is also an extreme
                                          if(angle_pairs(angle,1)==b .or. angle_pairs(angle,3)==b) then
                                                k = angle_pairs(angle,2)
                                                Vpp = comp_angle_energy(get_angle(disppi,xyz(:,k),disppj),aeq(angle),ka(angle))
                                                Vmm = comp_angle_energy(get_angle(dispmi,xyz(:,k),dispmj),aeq(angle),ka(angle))
                                                Vpm = comp_angle_energy(get_angle(disppi,xyz(:,k),dispmj),aeq(angle),ka(angle))
                                                Vmp = comp_angle_energy(get_angle(dispmi,xyz(:,k),disppj),aeq(angle),ka(angle))
                                          else !If b is the middle atom
                                                if(angle_pairs(angle,1)==a) k=angle_pairs(angle,3)
                                                if(angle_pairs(angle,3)==a) k=angle_pairs(angle,1)
                                                Vpp = comp_angle_energy(get_angle(disppi,disppj,xyz(:,k)),aeq(angle),ka(angle))
                                                Vmm = comp_angle_energy(get_angle(dispmi,dispmj,xyz(:,k)),aeq(angle),ka(angle))
                                                Vpm = comp_angle_energy(get_angle(disppi,dispmj,xyz(:,k)),aeq(angle),ka(angle))
                                                Vmp = comp_angle_energy(get_angle(dispmi,disppj,xyz(:,k)),aeq(angle),ka(angle))
                                          endif
                                    else !If a is the middle atom (b is then an extreme)
                                          if(angle_pairs(angle,1)==b) k=angle_pairs(angle,3)
                                          if(angle_pairs(angle,3)==b) k=angle_pairs(angle,1)
                                          Vpp = comp_angle_energy(get_angle(disppj,disppi,xyz(:,k)),aeq(angle),ka(angle))
                                          Vmm = comp_angle_energy(get_angle(dispmj,dispmi,xyz(:,k)),aeq(angle),ka(angle))
                                          Vpm = comp_angle_energy(get_angle(disppj,disppi,xyz(:,k)),aeq(angle),ka(angle))
                                          Vmp = comp_angle_energy(get_angle(dispmj,dispmi,xyz(:,k)),aeq(angle),ka(angle))
                                    endif

                              !If a is an extreme
                              elseif(any(angle_pairs(angle,:)==a).and.angle_pairs(angle,2)/=a) then
                                    if(angle_pairs(angle,1)==a) then 
                                          k=angle_pairs(angle,2)
                                          l=angle_pairs(angle,3)
                                    elseif(angle_pairs(angle,3)==a) then 
                                          k=angle_pairs(angle,2)
                                          l=angle_pairs(angle,1)
                                    endif
                                    Vpp = comp_angle_energy(get_angle(disppi,xyz(:,k),xyz(:,l)),aeq(angle),ka(angle))
                                    Vmm = comp_angle_energy(get_angle(dispmi,xyz(:,k),xyz(:,l)),aeq(angle),ka(angle))
                                    Vpm = comp_angle_energy(get_angle(disppi,xyz(:,k),xyz(:,l)),aeq(angle),ka(angle))
                                    Vmp = comp_angle_energy(get_angle(dispmi,xyz(:,k),xyz(:,l)),aeq(angle),ka(angle))
                              !If a is the middle atom
                              elseif(angle_pairs(angle,2)==a) then
                                    k=angle_pairs(angle,1)
                                    l=angle_pairs(angle,2)
                                    Vpp = comp_angle_energy(get_angle(xyz(:,k),disppi,xyz(:,l)),aeq(angle),ka(angle))
                                    Vmm = comp_angle_energy(get_angle(xyz(:,k),dispmi,xyz(:,l)),aeq(angle),ka(angle))
                                    Vpm = comp_angle_energy(get_angle(xyz(:,k),disppi,xyz(:,l)),aeq(angle),ka(angle))
                                    Vmp = comp_angle_energy(get_angle(xyz(:,k),dispmi,xyz(:,l)),aeq(angle),ka(angle))
                              !Same for b
                              elseif(any(angle_pairs(angle,:)==b).and.angle_pairs(angle,2)/=b) then
                                    if(angle_pairs(angle,1)==b) then
                                          k=angle_pairs(angle,2)
                                          l=angle_pairs(angle,3)
                                    elseif(angle_pairs(angle,3)==b) then
                                          k=angle_pairs(angle,2)
                                          l=angle_pairs(angle,1)
                                    endif
                                    Vpp = comp_angle_energy(get_angle(disppj,xyz(:,k),xyz(:,l)),aeq(angle),ka(angle))
                                    Vmm = comp_angle_energy(get_angle(dispmj,xyz(:,k),xyz(:,l)),aeq(angle),ka(angle))
                                    Vpm = comp_angle_energy(get_angle(dispmj,xyz(:,k),xyz(:,l)),aeq(angle),ka(angle))
                                    Vmp = comp_angle_energy(get_angle(disppj,xyz(:,k),xyz(:,l)),aeq(angle),ka(angle))
                              elseif(angle_pairs(angle,2)==b) then
                                    k=angle_pairs(angle,1)
                                    l=angle_pairs(angle,2)
                                    Vpp = comp_angle_energy(get_angle(xyz(:,k),disppj,xyz(:,l)),aeq(angle),ka(angle))
                                    Vmm = comp_angle_energy(get_angle(xyz(:,k),dispmj,xyz(:,l)),aeq(angle),ka(angle))
                                    Vpm = comp_angle_energy(get_angle(xyz(:,k),dispmj,xyz(:,l)),aeq(angle),ka(angle))
                                    Vmp = comp_angle_energy(get_angle(xyz(:,k),disppj,xyz(:,l)),aeq(angle),ka(angle))

                              else
                                    Vpp = 0.
                                    Vpm = 0.
                                    Vmp = 0.
                                    Vmm = 0.
                              endif
                              H(i,j) = H(i,j) + (Vpp - Vpm - Vmp + Vmm)/(4.d0*hi*hj)
                        enddo

                        !For all torsions
                        do torsion=1,Ntorsions

                              !If both a and b are extremes and a is the first one
                              if(torsion_pairs(torsion,1)==a .and. torsion_pairs(torsion,4)==b) then
                                    k = torsion_pairs(torsion,2)
                                    l = torsion_pairs(torsion,3)

                                    Vpp = comp_torsion_energy(get_torsion(disppi,xyz(:,k)&
                                    ,xyz(:,l),disppj),An(torsion),n(torsion),delta(torsion))
                                    Vmm = comp_torsion_energy(get_torsion(dispmi,xyz(:,k)&
                                    ,xyz(:,l),dispmj),An(torsion),n(torsion),delta(torsion))
                                    Vpm = comp_torsion_energy(get_torsion(disppi,xyz(:,k)&
                                    ,xyz(:,l),dispmj),An(torsion),n(torsion),delta(torsion))
                                    Vmp = comp_torsion_energy(get_torsion(dispmi,xyz(:,k)&
                                    ,xyz(:,l),disppj),An(torsion),n(torsion),delta(torsion))
                              !If both a and b are extremes and a is the last one
                              elseif(torsion_pairs(torsion,4)==a .and. torsion_pairs(torsion,1)==b) then
                                    k = torsion_pairs(torsion,2)
                                    l = torsion_pairs(torsion,3)

                                    Vpp = comp_torsion_energy(get_torsion(disppj,xyz(:,k)&
                                    ,xyz(:,l),disppi),An(torsion),n(torsion),delta(torsion))
                                    Vmm = comp_torsion_energy(get_torsion(dispmj,xyz(:,k)&
                                    ,xyz(:,l),dispmi),An(torsion),n(torsion),delta(torsion))
                                    Vpm = comp_torsion_energy(get_torsion(dispmj,xyz(:,k)&
                                    ,xyz(:,l),disppi),An(torsion),n(torsion),delta(torsion))
                                    Vmp = comp_torsion_energy(get_torsion(disppj,xyz(:,k)&
                                    ,xyz(:,l),dispmi),An(torsion),n(torsion),delta(torsion))
                              !If a is an extreme and b is a middle atom next to a
                              elseif((torsion_pairs(torsion,1)==a .and. torsion_pairs(torsion,2)==b) .or.&
                                    (torsion_pairs(torsion,4)==a .and. torsion_pairs(torsion,3)==b)) then

                                    if(torsion_pairs(torsion,1)==a) then
                                          k=torsion_pairs(torsion,3)
                                          l=torsion_pairs(torsion,4)
                                    elseif(torsion_pairs(torsion,4)==a) then
                                          k=torsion_pairs(torsion,2)
                                          l=torsion_pairs(torsion,1)
                                    endif
                                    Vpp = comp_torsion_energy(get_torsion(disppi,disppj,xyz(:,k)&
                                    ,xyz(:,l)),An(torsion),n(torsion),delta(torsion))
                                    Vmm = comp_torsion_energy(get_torsion(dispmi,dispmj,xyz(:,k)&
                                    ,xyz(:,l)),An(torsion),n(torsion),delta(torsion))
                                    Vpm = comp_torsion_energy(get_torsion(disppi,dispmj,xyz(:,k)&
                                    ,xyz(:,l)),An(torsion),n(torsion),delta(torsion))
                                    Vmp = comp_torsion_energy(get_torsion(dispmi,disppj,xyz(:,k)&
                                    ,xyz(:,l)),An(torsion),n(torsion),delta(torsion))
                              !If a is an extreme and b is a middle atom not next to a
                              elseif((torsion_pairs(torsion,1)==a .and. torsion_pairs(torsion,3)==b) .or.&
                                    (torsion_pairs(torsion,4)==a .and. torsion_pairs(torsion,2)==b)) then

                                    if(torsion_pairs(torsion,1)==a) then
                                          k=torsion_pairs(torsion,2)
                                          l=torsion_pairs(torsion,4)
                                    elseif(torsion_pairs(torsion,4)==a) then
                                          k=torsion_pairs(torsion,3)
                                          l=torsion_pairs(torsion,1)
                                    endif
                                    Vpp = comp_torsion_energy(get_torsion(disppi,xyz(:,k),disppj&
                                    ,xyz(:,l)),An(torsion),n(torsion),delta(torsion))
                                    Vmm = comp_torsion_energy(get_torsion(dispmi,xyz(:,k),dispmj&
                                    ,xyz(:,l)),An(torsion),n(torsion),delta(torsion))
                                    Vpm = comp_torsion_energy(get_torsion(disppi,xyz(:,k),dispmj&
                                    ,xyz(:,l)),An(torsion),n(torsion),delta(torsion))
                                    Vmp = comp_torsion_energy(get_torsion(dispmi,xyz(:,k),disppj&
                                    ,xyz(:,l)),An(torsion),n(torsion),delta(torsion))
                              !If b is an extreme and a is a middle atom next to b
                              elseif((torsion_pairs(torsion,1)==b .and. torsion_pairs(torsion,2)==a) .or.&
                                    (torsion_pairs(torsion,4)==b .and. torsion_pairs(torsion,3)==a)) then

                                    if(torsion_pairs(torsion,1)==b) then
                                          k=torsion_pairs(torsion,3)
                                          l=torsion_pairs(torsion,4)
                                    elseif(torsion_pairs(torsion,4)==b) then
                                          k=torsion_pairs(torsion,2)
                                          l=torsion_pairs(torsion,1)
                                    endif
                                    Vpp = comp_torsion_energy(get_torsion(disppj,disppi,xyz(:,k)&
                                    ,xyz(:,l)),An(torsion),n(torsion),delta(torsion))
                                    Vmm = comp_torsion_energy(get_torsion(dispmj,dispmi,xyz(:,k)&
                                    ,xyz(:,l)),An(torsion),n(torsion),delta(torsion))
                                    Vpm = comp_torsion_energy(get_torsion(dispmj,disppi,xyz(:,k)&
                                    ,xyz(:,l)),An(torsion),n(torsion),delta(torsion))
                                    Vmp = comp_torsion_energy(get_torsion(disppj,dispmi,xyz(:,k)&
                                    ,xyz(:,l)),An(torsion),n(torsion),delta(torsion))
                              !If b is an extreme and a is a middle atom not next to b
                              elseif((torsion_pairs(torsion,1)==b .and. torsion_pairs(torsion,3)==a) .or.&
                                    (torsion_pairs(torsion,4)==b .and. torsion_pairs(torsion,2)==a)) then

                                    if(torsion_pairs(torsion,1)==b) then
                                          k=torsion_pairs(torsion,2)
                                          l=torsion_pairs(torsion,4)
                                    elseif(torsion_pairs(torsion,4)==b) then
                                          k=torsion_pairs(torsion,3)
                                          l=torsion_pairs(torsion,1)
                                    endif
                                    Vpp = comp_torsion_energy(get_torsion(disppj,xyz(:,k),disppi&
                                    ,xyz(:,l)),An(torsion),n(torsion),delta(torsion))
                                    Vmm = comp_torsion_energy(get_torsion(dispmj,xyz(:,k),dispmi&
                                    ,xyz(:,l)),An(torsion),n(torsion),delta(torsion))
                                    Vpm = comp_torsion_energy(get_torsion(dispmj,xyz(:,k),disppi&
                                    ,xyz(:,l)),An(torsion),n(torsion),delta(torsion))
                                    Vmp = comp_torsion_energy(get_torsion(disppj,xyz(:,k),dispmi&
                                    ,xyz(:,l)),An(torsion),n(torsion),delta(torsion))
                                          
                              !If a is an extreme
                              elseif(any(torsion_pairs(torsion,:)==a).and.torsion_pairs(torsion,2)/=a &
                              .and.torsion_pairs(torsion,3)/=a) then
                                    if(torsion_pairs(torsion,1)==a) then
                                          k=torsion_pairs(torsion,2)
                                          l=torsion_pairs(torsion,3)
                                          m=torsion_pairs(torsion,4)
                                    elseif(torsion_pairs(torsion,4)==a) then
                                          k=torsion_pairs(torsion,2)
                                          l=torsion_pairs(torsion,3)
                                          m=torsion_pairs(torsion,1)
                                    endif

                                    Vpp = comp_torsion_energy(get_torsion(disppi,xyz(:,k)&
                                    ,xyz(:,l),xyz(:,m)),An(torsion),n(torsion),delta(torsion))
                                    Vmm = comp_torsion_energy(get_torsion(dispmi,xyz(:,k)&
                                    ,xyz(:,l),xyz(:,m)),An(torsion),n(torsion),delta(torsion))
                                    Vpm = comp_torsion_energy(get_torsion(disppi,xyz(:,k)&
                                    ,xyz(:,l),xyz(:,m)),An(torsion),n(torsion),delta(torsion))
                                    Vmp = comp_torsion_energy(get_torsion(dispmi,xyz(:,k)&
                                    ,xyz(:,l),xyz(:,m)),An(torsion),n(torsion),delta(torsion))

                              !If a is the first middle atom
                              elseif(torsion_pairs(torsion,2)==a) then
                                    k=torsion_pairs(torsion,1)
                                    l=torsion_pairs(torsion,3)
                                    m=torsion_pairs(torsion,4)
                                    Vpp = comp_torsion_energy(get_torsion(xyz(:,k),disppi&
                                    ,xyz(:,l),xyz(:,m)),An(torsion),n(torsion),delta(torsion))
                                    Vmm = comp_torsion_energy(get_torsion(xyz(:,k),dispmi&
                                    ,xyz(:,l),xyz(:,m)),An(torsion),n(torsion),delta(torsion))
                                    Vpm = comp_torsion_energy(get_torsion(xyz(:,k),disppi&
                                    ,xyz(:,l),xyz(:,m)),An(torsion),n(torsion),delta(torsion))
                                    Vmp = comp_torsion_energy(get_torsion(xyz(:,k),dispmi&
                                    ,xyz(:,l),xyz(:,m)),An(torsion),n(torsion),delta(torsion))
                              !If a is the second middle atom
                              elseif(torsion_pairs(torsion,3)==a) then
                                    k=torsion_pairs(torsion,1)
                                    l=torsion_pairs(torsion,2)
                                    m=torsion_pairs(torsion,4)
                                    Vpp = comp_torsion_energy(get_torsion(xyz(:,k),xyz(:,l)&
                                    ,disppi,xyz(:,m)),An(torsion),n(torsion),delta(torsion))
                                    Vmm = comp_torsion_energy(get_torsion(xyz(:,k),xyz(:,l)&
                                    ,dispmi,xyz(:,m)),An(torsion),n(torsion),delta(torsion))
                                    Vpm = comp_torsion_energy(get_torsion(xyz(:,k),xyz(:,l)&
                                    ,disppi,xyz(:,m)),An(torsion),n(torsion),delta(torsion))
                                    Vmp = comp_torsion_energy(get_torsion(xyz(:,k),xyz(:,l)&
                                    ,dispmi,xyz(:,m)),An(torsion),n(torsion),delta(torsion))

                              !If b is an extreme
                              elseif(any(torsion_pairs(torsion,:)==b).and.torsion_pairs(torsion,2)/=b &
                              .and.torsion_pairs(torsion,3)/=b) then
                                    if(torsion_pairs(torsion,1)==b) then
                                          k=torsion_pairs(torsion,2)
                                          l=torsion_pairs(torsion,3)
                                          m=torsion_pairs(torsion,4)
                                    elseif(torsion_pairs(torsion,4)==b) then
                                          k=torsion_pairs(torsion,2)
                                          l=torsion_pairs(torsion,3)
                                          m=torsion_pairs(torsion,1)
                                    endif

                                    Vpp = comp_torsion_energy(get_torsion(disppj,xyz(:,k)&
                                    ,xyz(:,l),xyz(:,m)),An(torsion),n(torsion),delta(torsion))
                                    Vmm = comp_torsion_energy(get_torsion(dispmj,xyz(:,k)&
                                    ,xyz(:,l),xyz(:,m)),An(torsion),n(torsion),delta(torsion))
                                    Vpm = comp_torsion_energy(get_torsion(dispmj,xyz(:,k)&
                                    ,xyz(:,l),xyz(:,m)),An(torsion),n(torsion),delta(torsion))
                                    Vmp = comp_torsion_energy(get_torsion(disppj,xyz(:,k)&
                                    ,xyz(:,l),xyz(:,m)),An(torsion),n(torsion),delta(torsion))
                              !If b is the first middle atom
                              elseif(torsion_pairs(torsion,2)==b) then
                                    k=torsion_pairs(torsion,1)
                                    l=torsion_pairs(torsion,3)
                                    m=torsion_pairs(torsion,4)
                                    Vpp = comp_torsion_energy(get_torsion(xyz(:,k),disppj&
                                    ,xyz(:,l),xyz(:,m)),An(torsion),n(torsion),delta(torsion))
                                    Vmm = comp_torsion_energy(get_torsion(xyz(:,k),dispmj&
                                    ,xyz(:,l),xyz(:,m)),An(torsion),n(torsion),delta(torsion))
                                    Vpm = comp_torsion_energy(get_torsion(xyz(:,k),dispmj&
                                    ,xyz(:,l),xyz(:,m)),An(torsion),n(torsion),delta(torsion))
                                    Vmp = comp_torsion_energy(get_torsion(xyz(:,k),disppj&
                                    ,xyz(:,l),xyz(:,m)),An(torsion),n(torsion),delta(torsion))
                              !If b is the second middle atom
                              elseif(torsion_pairs(torsion,3)==b) then
                                    k=torsion_pairs(torsion,1)
                                    l=torsion_pairs(torsion,2)
                                    m=torsion_pairs(torsion,4)
                                    Vpp = comp_torsion_energy(get_torsion(xyz(:,k),xyz(:,l)&
                                    ,disppj,xyz(:,m)),An(torsion),n(torsion),delta(torsion))
                                    Vmm = comp_torsion_energy(get_torsion(xyz(:,k),xyz(:,l)&
                                    ,dispmj,xyz(:,m)),An(torsion),n(torsion),delta(torsion))
                                    Vpm = comp_torsion_energy(get_torsion(xyz(:,k),xyz(:,l)&
                                    ,dispmj,xyz(:,m)),An(torsion),n(torsion),delta(torsion))
                                    Vmp = comp_torsion_energy(get_torsion(xyz(:,k),xyz(:,l)&
                                    ,disppj,xyz(:,m)),An(torsion),n(torsion),delta(torsion))
                              else
                                    Vpp = 0.
                                    Vpm = 0.
                                    Vmp = 0.
                                    Vmm = 0.
                              endif
                              H(i,j) = H(i,j) + (Vpp - Vpm - Vmp + Vmm)/(4.d0*hi*hj)
                        enddo




                              enddo
                              enddo
                        enddo
                  enddo
      end function build_hessian

      end module ff_utils