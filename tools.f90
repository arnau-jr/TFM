      module tools
            implicit none
            real*8,parameter :: pi=4.*atan(1.) 
            real*8,parameter :: cov(10) = (/0.37,0.,0.,0.,0.,0.77,0.75,0.73,0.,0./)
            real*8,parameter :: mass(10) = (/1.008,0.,0.,0.,0.,12.01,14.01,16.00,0.,0./)

            real*8,allocatable :: xyz(:,:)

            contains

            subroutine parse_atomic_symbol(N,S,Z,M)
                  implicit none
                  integer :: N
                  character :: S(N)*2
                  real*8 :: M(N)
                  integer :: Z(N),i

                  do i=1,N
                        if(S(i)=="H") then
                              Z(i) = 1
                        elseif(S(i)=="C") then
                              Z(i) = 6
                        elseif(S(i)=="N") then
                              Z(i) = 7
                        elseif(S(i)=="O") then
                              Z(i) = 8
                        else
                              write(*,*)"Parsing: atom not recognized"
                        endif
                        M(i) = mass(Z(i))
                  enddo
            end

            function unit_cross(u1,u2) result(u3)
                  implicit none
                  real*8 :: u1(3),u2(3)
                  real*8 :: u3(3)
      
                  u3(1) = u1(2)*u2(3) - u1(3)*u2(2)
                  u3(2) = - u1(1)*u2(3) + u1(3)*u2(1)
                  u3(3) = u1(1)*u2(2) - u1(2)*u2(1)
      
                  u3 = u3/sqrt(sum(u3**2))
            end function unit_cross

            subroutine get_local_coords(c1,c2,c3,basis)
                  implicit none
                  real*8 :: c1(3),c2(3),c3(3)
                  real*8 :: basis(3,3) !First coordinates, second vectors
                  real*8 :: u1(3),u2(3)
                  real*8 :: n(3),nn(3)
      
                  u1 = c2-c1
                  u2 = c3-c2
      
                  u1 = u1/sqrt(sum(u1**2))
                  u2 = u2/sqrt(sum(u2**2))
      
                  n = unit_cross(u2,u1)
                  nn = unit_cross(u1,n)
                  
                  basis(:,3) = u1
                  basis(:,2) = nn
                  n = unit_cross(basis(:,2),basis(:,3))
                  basis(:,1) = n
            end

            function get_dist(c1,c2) result (d)
                  implicit none
                  real*8 :: c1(3),c2(3)
                  real*8 :: d

                  d = sqrt(sum((c1-c2)**2))
            end function get_dist

            function get_angle(c1,c2,c3) result(A)
                  implicit none
                  real*8 :: c1(3),c2(3),c3(3)
                  real*8 :: u1(3),u2(3),proj
                  real*8 :: A

                  u1 = c2-c1
                  u2 = c3-c2
      
                  u1 = u1/sqrt(sum(u1**2))
                  u2 = u2/sqrt(sum(u2**2))

                  proj = sum(u1*u2)
                  if(proj>=1.d0) then
                        A = 180.d0
                  elseif(proj<=-1.d0) then
                        A = 0.d0
                  else
                        A = 180.d0 - (180.d0/pi)*acos(proj)
                  endif
            end function get_angle

            function get_torsion(c1,c2,c3,c4) result(T)
                  implicit none
                  real*8 :: c1(3),c2(3),c3(3),c4(3)
                  real*8 :: u12(3),u23(3),u32(3),u43(3)
                  real*8 :: u1232(3),u2343(3)
                  real*8 :: proj,proj2
                  real*8 :: T

                  u12 = c2-c1
                  u23 = c3-c2
                  u32 = c2-c3
                  u43 = c3-c4
      
                  u12 = u12/sqrt(sum(u12**2))
                  u23 = u23/sqrt(sum(u23**2))
                  u32 = u32/sqrt(sum(u32**2))
                  u43 = u43/sqrt(sum(u43**2))

                  u1232 = unit_cross(u12,u32)
                  u2343 = unit_cross(u23,u43)
                  
                  proj = sum(u1232*u2343)
                  proj2 = sum(u1232*u43)

                  ! if(proj>=1.d0) then
                  !       T = 180.d0
                  ! elseif(proj<=-1.d0) then
                  !       T = 180.d0 - 180.d0*sign(1.d0,proj2)
                  ! else
                  !       T = 180.d0 - (180.d0/pi)*sign(1.d0,proj2)*acos(proj)
                  ! endif
                  if(proj>=1.d0) then
                        T = 0.d0*sign(1.d0,-proj2)
                  elseif(proj<=-1.d0) then
                        T = 180.d0*sign(1.d0,-proj2)
                  else
                        T = (180.d0/pi)*sign(1.d0,-proj2)*acos(proj)
                  endif
            end function get_torsion

            function get_improper(c4,c1,c3,c2) result(T)
                  implicit none
                  real*8 :: c1(3),c2(3),c3(3),c4(3)
                  real*8 :: u12(3),u23(3),u32(3),u43(3)
                  real*8 :: u1232(3),u2343(3)
                  real*8 :: proj,proj2
                  real*8 :: T

                  u12 = c2-c1
                  u23 = c3-c2
                  u32 = c2-c3
                  u43 = c3-c4
      
                  u12 = u12/sqrt(sum(u12**2))
                  u23 = u23/sqrt(sum(u23**2))
                  u32 = u32/sqrt(sum(u32**2))
                  u43 = u43/sqrt(sum(u43**2))

                  u1232 = unit_cross(u12,u32)
                  u2343 = unit_cross(u23,u43)
                  
                  proj = sum(u1232*u2343)
                  proj2 = sum(u1232*u43)

                  ! if(proj>=1.d0) then
                  !       T = 180.d0
                  ! elseif(proj<=-1.d0) then
                  !       T = 180.d0 - 180.d0*sign(1.d0,proj2)
                  ! else
                  !       T = 180.d0 - (180.d0/pi)*sign(1.d0,proj2)*acos(proj)
                  ! endif
                  if(proj>=1.d0) then
                        T = 0.d0*sign(1.d0,-proj2)
                  elseif(proj<=-1.d0) then
                        T = 180.d0*sign(1.d0,-proj2)
                  else
                        T = (180.d0/pi)*sign(1.d0,-proj2)*acos(proj)
                  endif
            end function get_improper

            function get_improper_bac(c4,c1,c3,c2) result(I)
                  implicit none
                  real*8 :: c1(3),c2(3),c3(3),c4(3)
                  real*8 :: u24(3),u34(3),u14(3)
                  real*8 :: u2434(3)
                  real*8 :: proj
                  real*8 :: I

                  u24 = c4-c2
                  u34 = c4-c3
                  u14 = c4-c1
      
                  u24 = u24/sqrt(sum(u24**2))
                  u34 = u34/sqrt(sum(u34**2))
                  u14 = u14/sqrt(sum(u14**2))

                  u2434 = unit_cross(u24,u34)
                  
                  proj = sum(u2434*u14)

                  if(proj>=1.d0) then
                        I = 90.d0
                  elseif(proj<=-1.d0) then
                        I = -90.d0
                  else
                        I = (180.d0/pi)*asin(proj)
                  endif
            end function get_improper_bac


            subroutine get_xyz(port,N,S,xyz)
                  implicit none
                  integer :: port,N
                  character,allocatable :: S(:)*2
                  real*8,allocatable :: xyz(:,:)
                  integer :: i

                  read(port,*)N
                  read(port,*)

                  allocate(S(N),xyz(3,N))

                  do i = 1,N
                        read(port,*)S(i),xyz(:,i)
                  enddo
            end subroutine get_xyz

            function get_dist_matrix(N,xyz) result(dmat)
                  implicit none
                  integer :: N
                  real*8 :: xyz(3,N)
                  real*8,allocatable :: dmat(:,:)
                  integer :: i,j

                  allocate(dmat(N,N))
                  dmat = 0.
                  do i = 1,N
                        do j = i+1,N
                              dmat(i,j) = sqrt(sum((xyz(:,i)-xyz(:,j))**2))
                              dmat(j,i) = dmat(i,j)
                        enddo
                  enddo
            end function get_dist_matrix

            function get_bond_graph(N,Z,dmat) result(bond_graph)
                  implicit none
                  integer :: N
                  integer :: Z(N)
                  real*8 :: dmat(N,N)
                  logical,allocatable :: bond_graph(:,:)
                  real*8,parameter :: thresh = 1.2
                  integer :: counter,i,j

                  
                  allocate(bond_graph(N,N))
                  bond_graph = .false.
                  do i = 1,N
                        do j = i+1,N
                              if(dmat(i,j) < thresh*(cov(Z(i))+cov(Z(j)))) then
                                    bond_graph(i,j) = .true.
                                    bond_graph(j,i) = .true.
                              endif
                        enddo
                  enddo
            end function get_bond_graph

            subroutine get_bonds(N,Z,dmat,bond_pairsbis,bond_valsbis)
                  implicit none
                  integer :: N
                  integer :: Z(N)
                  real*8 :: dmat(N,N)
                  integer,allocatable :: bond_pairs(:,:),bond_pairsbis(:,:)
                  real*8,allocatable :: bond_vals(:),bond_valsbis(:)
                  real*8,parameter :: thresh = 1.2
                  integer :: Nbonds,i,j
                  
                  allocate(bond_pairs(N*N,2),bond_vals(N*N))
                  Nbonds = 0
                  do i = 1,N
                        do j = i+1,N
                              if(dmat(i,j) < thresh*(cov(Z(i))+cov(Z(j)))) then
                                    Nbonds = Nbonds + 1
                                    bond_pairs(Nbonds,:) = (/i,j/)
                                    bond_vals(Nbonds) = dmat(i,j)
                              endif
                        enddo
                  enddo
                  if(allocated(bond_pairsbis).eqv..false.) then
                        allocate(bond_pairsbis(Nbonds,2),bond_valsbis(Nbonds))
                  endif
                  bond_pairsbis = bond_pairs(:Nbonds,:)
                  bond_valsbis = bond_vals(:Nbonds)
            end subroutine get_bonds

            subroutine get_angles(N,xyz,bond_graph,angle_pairsbis,angle_valsbis)
                  implicit none
                  integer :: N
                  logical :: bond_graph(N,N)
                  real*8 :: xyz(3,N)
                  integer,allocatable :: angle_pairs(:,:),angle_pairsbis(:,:)
                  real*8,allocatable :: angle_vals(:),angle_valsbis(:)
                  integer :: i,j,k
                  integer :: a,b,c
                  integer :: Nangles

                  allocate(angle_pairs(N*N,3),angle_vals(N*N))
                  Nangles = 0
                  do i = 1,N
                        a = i
                        do j = 1,N
                              if(bond_graph(a,j).eqv..true.) then
                                    b = j
                                    do k = a+1,N
                                          if(bond_graph(b,k).eqv..true.) then
                                                c = k
                                                Nangles = Nangles + 1
                                                angle_pairs(Nangles,:) = (/a,b,c/)
                                                angle_vals(Nangles) = get_angle(xyz(:,a),xyz(:,b),xyz(:,c))
                                          endif
                                    enddo
                              endif
                        enddo
                  enddo
                  if(allocated(angle_pairsbis).eqv..false.) then
                        allocate(angle_pairsbis(Nangles,3),angle_valsbis(Nangles))
                  endif
                  angle_pairsbis = angle_pairs(:Nangles,:)
                  angle_valsbis = angle_vals(:Nangles)
            end subroutine get_angles

            subroutine get_torsions(N,xyz,bond_graph,torsion_pairsbis,torsion_valsbis)
                  implicit none
                  integer :: N
                  logical :: bond_graph(N,N)
                  real*8 :: xyz(3,N)
                  integer,allocatable :: torsion_pairs(:,:),torsion_pairsbis(:,:)
                  real*8,allocatable :: torsion_vals(:),torsion_valsbis(:)
                  integer :: i,j,k,l
                  integer :: a,b,c,d
                  integer :: Ntorsions

                  allocate(torsion_pairs(N*N,4),torsion_vals(N*N))
                  Ntorsions = 0
                  do i = 1,N
                        a = i
                        do j = 1,N
                              if(bond_graph(a,j).eqv..true.) then
                                    ! if(j>a) then
                                    !       cycle
                                    ! endif
                                    b = j
                                    do k = 1,N
                                          if(bond_graph(b,k).eqv..true.) then
                                                ! if(k>=a) then
                                                if(k==a) then
                                                      cycle
                                                endif
                                                c = k
                                                do l=1,N
                                                      if(bond_graph(c,l).eqv..true.) then
                                                            ! if(l >= a .or. l == b) then
                                                            !       cycle
                                                            ! endif
                                                            if(l >= a .or. l == b) then
                                                                  cycle
                                                            endif
                                                            d = l
                                                            Ntorsions = Ntorsions + 1
                                                            torsion_pairs(Ntorsions,:) = (/a,b,c,d/)
                                                            torsion_vals(Ntorsions) &
                                                      = get_torsion(xyz(:,a),xyz(:,b),xyz(:,c),xyz(:,d))
                                                      endif
                                                enddo
                                          endif
                                    enddo
                              endif
                        enddo
                  enddo
                  if(allocated(torsion_pairsbis).eqv..false.) then
                        allocate(torsion_pairsbis(Ntorsions,4),torsion_valsbis(Ntorsions))
                  endif
                  torsion_pairsbis = torsion_pairs(:Ntorsions,:)
                  torsion_valsbis = torsion_vals(:Ntorsions)
            end subroutine get_torsions

            function recomp_bonds(Natoms,Nbonds,xyz,bond_pairs) result(bond_vals)
                  implicit none
                  integer :: Natoms,Nbonds
                  integer :: bond_pairs(Nbonds,2)
                  real*8 :: xyz(3,Natoms),bond_vals(Nbonds)
                  integer :: bond

                  do bond=1,Nbonds
                        bond_vals(bond) = get_dist(xyz(:,bond_pairs(bond,1)),&
                        xyz(:,bond_pairs(bond,2)))
                  enddo
      end function recomp_bonds

             function recomp_angles(Natoms,Nangles,xyz,angle_pairs) result(angle_vals)
                  implicit none
                  integer :: Natoms,Nangles
                  integer :: angle_pairs(Nangles,3)
                  real*8 :: xyz(3,Natoms),angle_vals(Nangles)
                  integer :: angle

                  do angle=1,Nangles
                        angle_vals(angle) = get_angle(xyz(:,angle_pairs(angle,1)),&
                        xyz(:,angle_pairs(angle,2)),xyz(:,angle_pairs(angle,3)))
                  enddo
      end function recomp_angles

                  function recomp_torsions(Natoms,Ntorsions,xyz,torsion_pairs) result(torsion_vals)
                  implicit none
                  integer :: Natoms,Ntorsions
                  integer :: torsion_pairs(Ntorsions,4)
                  real*8 :: xyz(3,Natoms),torsion_vals(Ntorsions)
                  integer :: torsion

                  do torsion=1,Ntorsions
                        torsion_vals(torsion) = get_torsion(xyz(:,torsion_pairs(torsion,1)),&
                        xyz(:,torsion_pairs(torsion,2)),&
                        xyz(:,torsion_pairs(torsion,3)),&
                        xyz(:,torsion_pairs(torsion,4)))
                  enddo
      end function recomp_torsions

                  function recomp_impropers(Natoms,Nimpropers,xyz,improper_pairs) result(improper_vals)
                  implicit none
                  integer :: Natoms,Nimpropers
                  integer :: improper_pairs(Nimpropers,4)
                  real*8 :: xyz(3,Natoms),improper_vals(Nimpropers)
                  integer :: improper

                  do improper=1,Nimpropers
                        improper_vals(improper) = get_improper(xyz(:,improper_pairs(improper,1)),&
                        xyz(:,improper_pairs(improper,2)),&
                        xyz(:,improper_pairs(improper,3)),&
                        xyz(:,improper_pairs(improper,4)))
                  enddo
      end function recomp_impropers

            subroutine save_zmat(port,N,Nangle,Ntorsion,S,dmat,bond_graph,angle_pairs,angle_vals,torsion_pairs,torsion_vals)
                  implicit none
                  integer :: port,N,Nangle,Ntorsion
                  character :: S(N)*2
                  integer :: angle_pairs(Nangle,3),torsion_pairs(Ntorsion,4)
                  real*8 :: dmat(N,N),angle_vals(Nangle),torsion_vals(Ntorsion)
                  logical :: bond_graph(N,N)

                  integer :: Batom,Aatom,Tatom
                  real*8 :: B,A,T
                  integer :: i,j,k

                  character :: fmtlabel*90

                  write(fmtlabel,"(A)")"(A,1X,I2,1X,F4.2,1X,I2,1X,F6.1,1X,I2,1X,F6.1)"

                  write(port,"(A)")"#"
                  write(port,*)
                  write(port,*)"Comment"
                  write(port,*)
                  write(port,"(A)")"0,1"

                  write(port,"(A)")S(1)


                  Batom = 1
                  B = dmat(2,Batom)
                  write(port,"(A,1X,I2,1X,F4.2)")S(2),Batom,B


                  do j=1,3
                        if(bond_graph(3,j).eqv..true.) then
                              Batom = j
                              exit
                        endif
                  enddo
                  B = dmat(3,Batom)
                  do i=1,Nangle
                        if(angle_pairs(i,1)==3 .and. angle_pairs(i,3)<3) then
                              Aatom = angle_pairs(i,3)
                              A = angle_vals(i)
                        elseif(angle_pairs(i,3)==3 .and. angle_pairs(i,1)<3) then
                              Aatom = angle_pairs(i,1)
                              A = angle_vals(i)
                        endif
                  enddo
                  write(port,"(A,1X,I2,1X,F4.2,1X,I2,1X,F6.1)")S(3),Batom,B,Aatom,A

                  do i=4,N
                        do j=1,i
                              if(bond_graph(i,j).eqv..true.) then
                                    Batom = j
                                    exit
                              endif
                        enddo
                        B = dmat(i,Batom)

                        do j=1,Nangle
                              if(angle_pairs(j,1)==i .and. angle_pairs(j,3)<i &
                              .and. angle_pairs(j,2)<i) then
                                    Aatom = angle_pairs(j,3)
                                    A = angle_vals(j)
                                    exit
                              elseif(angle_pairs(j,3)==i .and. angle_pairs(j,1)<i &
                                    .and. angle_pairs(j,2)<i) then
                                    Aatom = angle_pairs(j,1)
                                    A = angle_vals(j)
                                    exit
                              endif
                        enddo

                        do j=1,Ntorsion
                              if(torsion_pairs(j,1)==i .and. torsion_pairs(j,4)<i) then
                                    Tatom = torsion_pairs(j,4)
                                    T = torsion_vals(j)
                                    exit
                              elseif(torsion_pairs(j,4)==i .and. torsion_pairs(j,1)<i) then
                                    Tatom = torsion_pairs(j,1)
                                    T = torsion_vals(j)
                              endif
                        enddo

                        write(port,fmtlabel)S(i),Batom,B,Aatom,A,Tatom,T
                  enddo
            end subroutine save_zmat

      end module tools