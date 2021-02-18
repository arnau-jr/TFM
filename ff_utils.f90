      module ff_utils
            use tools
            implicit none
            integer :: Nimpropers
            integer,allocatable :: improper_pairs(:,:)
            real*8,allocatable :: improper_vals(:)
            real*8,allocatable :: req(:),kb(:) !Harmonic stretching
            real*8,allocatable :: De(:),beta(:) !Morse stretching 
            character :: pot_type*1 !Decides if harmonic or morse is used,
            !this decision is all or nothing: all must be harmonics or all must be morses.
            real*8,allocatable :: aeq(:),ka(:) !Harmonic angles
            real*8,allocatable :: An(:),n(:),delta(:) !Cosine torsions
            real*8,allocatable :: impAn(:),impn(:),impdelta(:) !Cosine impropers
            real*8,parameter :: kcal_to_kj = 4.184d0
            real*8,parameter :: h_cm_dps = 8065.6d0*4.135667696d-2
            real*8,parameter :: hbar_cm_dps = 8065.6d0*6.582119569d-3

            contains

            subroutine get_param(port,Nbonds,Nangles,Ntorsions)
                  implicit none
                  integer :: port, Nbonds,Nangles,Ntorsions
                  character :: dummy*90,dummy2*90,units*90
                  integer :: i,a,b,c,d,dummy3

                  Nimpropers = 0
                  allocate(aeq(Nangles),ka(Nangles))
                  allocate(An(Ntorsions),n(Ntorsions),delta(Ntorsions))
            
                  read(port,*)pot_type
                  if(pot_type=="H") allocate(req(Nbonds),kb(Nbonds))
                  if(pot_type=="M") allocate(De(Nbonds),beta(Nbonds),req(Nbonds))
                  do i=1,Nbonds
                        if(pot_type=="H") then
                              read(port,*)dummy,dummy2,kb(i),req(i)
                              kb(i) = kb(i)*kcal_to_kj
                        elseif(pot_type=="M") then
                              read(port,*)dummy,dummy2,De(i),beta(i),req(i)
                              !De should be in kJ/mol
                        endif
                  enddo
                  read(port,*,end=10)
                  read(port,*)units
                  do i=1,Nangles
                        if(units=="KC") then
                              read(port,*)dummy,dummy2,ka(i),aeq(i)
                              ka(i) = ka(i)*kcal_to_kj
                        elseif(units=="KJ") then
                              read(port,*)dummy,dummy2,ka(i),aeq(i)
                        endif
                  enddo
                  read(port,*,end=10)
                  read(port,*)units
                  do i=1,Ntorsions
                        if(units=="KC") then
                              read(port,*)dummy,dummy2,dummy3,An(i),delta(i),n(i)
                              An(i) = An(i)*kcal_to_kj
                        elseif(units=="KJ") then
                              read(port,*)dummy,dummy2,dummy3,An(i),delta(i),n(i)
                        endif
                  enddo
                  read(port,*,end=10)Nimpropers
                  read(port,*)units
                  allocate(improper_pairs(Nimpropers,4),improper_vals(Nimpropers))
                  allocate(impAn(Nimpropers),impn(Nimpropers),impdelta(Nimpropers))
                  do i=1,Nimpropers
                        if(units=="KC") then
                              read(port,*)dummy,a,b,c,d,impAn(i),impdelta(i),impn(i)
                              An(i) = An(i)*kcal_to_kj
                        elseif(units=="KJ") then
                              read(port,*)dummy,a,b,c,d,impAn(i),impdelta(i),impn(i)
                        endif
                        improper_pairs(i,:) = (/a,b,c,d/)
                        improper_vals(i) = get_improper(xyz(:,a),xyz(:,b),&
                        xyz(:,c),xyz(:,d))
                  enddo
                  10 continue
      end subroutine get_param


            real*8 function comp_bonds_energy(Nbonds,bond_vals) result(E)
                  implicit none
                  integer :: Nbonds
                  real*8 :: bond_vals(Nbonds)
                  integer :: i
                  E = 0.
                  if(pot_type=="H") then
                        do i=1,Nbonds
                              E = E + kb(i)*(bond_vals(i)-req(i))**2
                        enddo
                  elseif(pot_type=="M") then
                        do i=1,Nbonds
                              E = E + De(i)*((1.d0-exp(-beta(i)*(bond_vals(i)-req(i))))**2-1.d0)
                        enddo
                  endif
      end function comp_bonds_energy


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

            real*8 function comp_impropers_energy() result(E)
                  implicit none 
                  integer :: i
                  E = 0.
                  do i=1,Nimpropers
                        ! I = get_improper(xyz(:,improper_pairs(k,1)),xyz(:,improper_pairs(k,2)),&
                        ! xyz(:,improper_pairs(k,3)),xyz(:,improper_pairs(k,4)))
                        E = E + impAn(i)*(1. + cos((pi/180.d0)*(impn(i)*improper_vals(i) - impdelta(i))))
                  enddo
      end function comp_impropers_energy

            real*8 function comp_energy(Nbonds,Nangles,Ntorsions,bond_vals,&
            angle_vals,torsion_vals) result (E)
                  implicit none
                  integer :: Nbonds,Nangles,Ntorsions
                  real*8 :: bond_vals(Nbonds),angle_vals(Nangles),torsion_vals(Ntorsions)
                  E = 0.
                  E = E + comp_bonds_energy(Nbonds,bond_vals)
                  E = E + comp_angles_energy(Nangles,angle_vals)
                  E = E + comp_torsions_energy(Ntorsions,torsion_vals)
                  E = E + comp_impropers_energy()
      end function comp_energy

            function build_gradient(Natoms,xyz,Nbonds,Nangles,Ntorsions,&
            bond_pairs,angle_pairs,torsion_pairs) result(G)
                  implicit none
                  integer :: Natoms,Nbonds,Nangles,Ntorsions
                  real*8 :: xyz(3,Natoms),xyzmod(3,Natoms)
                  integer :: bond_pairs(Nbonds,2),angle_pairs(Nangles,3),torsion_pairs(Ntorsions,4)
                  real*8 :: bond_vals(Nbonds),angle_vals(Nangles),torsion_vals(Ntorsions)
                  real*8 :: G(3,Natoms)
                  real*8,parameter :: hi = 5.d-6
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
                              improper_vals = recomp_impropers(Natoms,Nimpropers,xyzmod,improper_pairs)

                              Vp = comp_energy(Nbonds,Nangles,Ntorsions,bond_vals,&
                              angle_vals,torsion_vals)

                              xyzmod = xyz
                              xyzmod(p,a) = xyz(p,a) - hi

                              bond_vals = recomp_bonds(Natoms,Nbonds,xyzmod,bond_pairs)
                              angle_vals = recomp_angles(Natoms,Nangles,xyzmod,angle_pairs)
                              torsion_vals = recomp_torsions(Natoms,Ntorsions,xyzmod,torsion_pairs)
                              improper_vals = recomp_impropers(Natoms,Nimpropers,xyzmod,improper_pairs)

                              Vm = comp_energy(Nbonds,Nangles,Ntorsions,bond_vals,&
                              angle_vals,torsion_vals)

                              G(p,a) = (Vp-Vm)/(2.d0*hi)
                        enddo

                  enddo
                  bond_vals = recomp_bonds(Natoms,Nbonds,xyz,bond_pairs)
                  angle_vals = recomp_angles(Natoms,Nangles,xyz,angle_pairs)
                  torsion_vals = recomp_torsions(Natoms,Ntorsions,xyz,torsion_pairs)
                  improper_vals = recomp_impropers(Natoms,Nimpropers,xyzmod,improper_pairs)
      end function build_gradient

            function build_hessian(Natoms,xyz,Nbonds,Nangles,Ntorsions,&
            bond_pairs,angle_pairs,torsion_pairs) result(H)
                  implicit none
                  integer :: Natoms,Nbonds,Nangles,Ntorsions
                  real*8 :: xyz(3,Natoms),xyzmod(3,Natoms)
                  integer :: bond_pairs(Nbonds,2),angle_pairs(Nangles,3),torsion_pairs(Ntorsions,4)
                  real*8 :: bond_vals(Nbonds),angle_vals(Nangles),torsion_vals(Ntorsions)
                  real*8 :: H(3*Natoms,3*Natoms)
                  real*8,parameter :: hi = 5.d-6,hj = 5.d-6
                  ! real*8,parameter :: hi = 1.d-5,hj = 1.d-5
                  real*8 :: V,Vpp,Vmm,Vpm,Vmp
                  integer :: a,b,p,q,i,j

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
                              improper_vals = recomp_impropers(Natoms,Nimpropers,xyzmod,improper_pairs)

                              Vpp = comp_energy(Nbonds,Nangles,Ntorsions,bond_vals,&
                              angle_vals,torsion_vals)

                              xyzmod = xyz
                              xyzmod(3-p,a) = xyz(3-p,a) - hi
                              xyzmod(3-q,b) = xyz(3-q,b) - hj

                              bond_vals = recomp_bonds(Natoms,Nbonds,xyzmod,bond_pairs)
                              angle_vals = recomp_angles(Natoms,Nangles,xyzmod,angle_pairs)
                              torsion_vals = recomp_torsions(Natoms,Ntorsions,xyzmod,torsion_pairs)
                              improper_vals = recomp_impropers(Natoms,Nimpropers,xyzmod,improper_pairs)

                              Vmm = comp_energy(Nbonds,Nangles,Ntorsions,bond_vals,&
                              angle_vals,torsion_vals)

                              xyzmod = xyz
                              xyzmod(3-p,a) = xyz(3-p,a) + hi
                              xyzmod(3-q,b) = xyz(3-q,b) - hj

                              bond_vals = recomp_bonds(Natoms,Nbonds,xyzmod,bond_pairs)
                              angle_vals = recomp_angles(Natoms,Nangles,xyzmod,angle_pairs)
                              torsion_vals = recomp_torsions(Natoms,Ntorsions,xyzmod,torsion_pairs)
                              improper_vals = recomp_impropers(Natoms,Nimpropers,xyzmod,improper_pairs)

                              Vpm = comp_energy(Nbonds,Nangles,Ntorsions,bond_vals,&
                              angle_vals,torsion_vals)

                              xyzmod = xyz
                              xyzmod(3-p,a) = xyz(3-p,a) - hi
                              xyzmod(3-q,b) = xyz(3-q,b) + hj

                              bond_vals = recomp_bonds(Natoms,Nbonds,xyzmod,bond_pairs)
                              angle_vals = recomp_angles(Natoms,Nangles,xyzmod,angle_pairs)
                              torsion_vals = recomp_torsions(Natoms,Ntorsions,xyzmod,torsion_pairs)
                              improper_vals = recomp_impropers(Natoms,Nimpropers,xyzmod,improper_pairs)

                              Vmp = comp_energy(Nbonds,Nangles,Ntorsions,bond_vals,&
                              angle_vals,torsion_vals)

                              H(i,j) = (Vpp - Vpm - Vmp + Vmm)/(4.d0*hi*hj)
                              else

                              xyzmod = xyz
                              xyzmod(3-p,a) = xyz(3-p,a) + hi

                              bond_vals = recomp_bonds(Natoms,Nbonds,xyzmod,bond_pairs)
                              angle_vals = recomp_angles(Natoms,Nangles,xyzmod,angle_pairs)
                              torsion_vals = recomp_torsions(Natoms,Ntorsions,xyzmod,torsion_pairs)
                              improper_vals = recomp_impropers(Natoms,Nimpropers,xyzmod,improper_pairs)

                              Vpp = comp_energy(Nbonds,Nangles,Ntorsions,bond_vals,&
                              angle_vals,torsion_vals)

                              xyzmod = xyz
                              xyzmod(3-p,a) = xyz(3-p,a) - hi

                              bond_vals = recomp_bonds(Natoms,Nbonds,xyzmod,bond_pairs)
                              angle_vals = recomp_angles(Natoms,Nangles,xyzmod,angle_pairs)
                              torsion_vals = recomp_torsions(Natoms,Ntorsions,xyzmod,torsion_pairs)
                              improper_vals = recomp_impropers(Natoms,Nimpropers,xyzmod,improper_pairs)

                              Vmm = comp_energy(Nbonds,Nangles,Ntorsions,bond_vals,&
                              angle_vals,torsion_vals)

                              xyzmod = xyz

                              bond_vals = recomp_bonds(Natoms,Nbonds,xyzmod,bond_pairs)
                              angle_vals = recomp_angles(Natoms,Nangles,xyzmod,angle_pairs)
                              torsion_vals = recomp_torsions(Natoms,Ntorsions,xyzmod,torsion_pairs)
                              improper_vals = recomp_impropers(Natoms,Nimpropers,xyzmod,improper_pairs)

                              V = comp_energy(Nbonds,Nangles,Ntorsions,bond_vals,&
                              angle_vals,torsion_vals)

                              H(i,j) = (Vpp - 2.d0*V + Vmm)/(hi*hj)
                              
                              endif
                              enddo
                              enddo
                        enddo
                  enddo
                  bond_vals = recomp_bonds(Natoms,Nbonds,xyz,bond_pairs)
                  angle_vals = recomp_angles(Natoms,Nangles,xyz,angle_pairs)
                  torsion_vals = recomp_torsions(Natoms,Ntorsions,xyz,torsion_pairs)
                  improper_vals = recomp_impropers(Natoms,Nimpropers,xyzmod,improper_pairs)
      end function build_hessian

      end module ff_utils