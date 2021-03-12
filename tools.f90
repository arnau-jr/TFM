      module tools
            implicit none
            real*8,parameter :: pi=4.*atan(1.) 
            real*8,parameter :: cov(10) = (/0.37,0.,0.,0.,0.,0.77,0.75,0.73,0.,0./)
            ! real*8,parameter :: mass(10) = (/1.008,0.,0.,0.,0.,12.01,14.01,16.00,0.,0./)
            real*8,parameter :: mass(10) = (/1.00782522,0.,0.,0.,0.,12.01,14.01,15.99491502,0.,0./)


            character,allocatable :: S(:)*2
            integer,allocatable :: Z(:)
            real*8,allocatable ::  M(:),xyz(:,:),xyz_mw(:,:)

            contains

            subroutine parse_atomic_symbol(N)
                  implicit none
                  integer :: N
                  integer :: i
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

            function get_norm(u) result(a)
                  implicit none
                  real*8 :: u(:),a
                  a = sqrt(sum(u**2))
      end function get_norm

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
                  ! proj2 = -1.d0

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


            subroutine get_xyz(port,N)
                  implicit none
                  integer :: port,N
                  integer :: i

                  read(port,*)N
                  read(port,*)

                  allocate(S(N),M(N),Z(N),xyz(3,N),xyz_mw(3,N))

                  do i = 1,N
                        read(port,*)S(i),xyz(:,i)
                  enddo
                  call parse_atomic_symbol(N)
                  do i=1,N
                        xyz_mw(:,i) = sqrt(M(i))*xyz(:,i)
                  enddo
            end subroutine get_xyz

            subroutine write_conf(D,N,r,port)
                  implicit none
                  integer :: N,D,port,i
                  real*8 :: r(D,N)
      
                  write(port,"(I5)")N
                  write(port,*)""
                  do i=1,N
                        write(port,"(A,2X,E20.10,2X,E20.10,2X,E20.10)")S(i),r(:,i)
                  enddo
            end

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

            function get_bond_graph(N,dmat) result(bond_graph)
                  implicit none
                  integer :: N
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

            subroutine get_bonds(N,dmat,bond_pairsbis,bond_valsbis)
                  implicit none
                  integer :: N
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

            function unpack_coords(packed) result (unpacked)
                  implicit none
                  real*8 :: packed(:,:)
                  real*8 :: unpacked(3*size(packed,2))
                  integer :: i
                  do i=1,size(packed,2)
                        unpacked(3*i-2) = packed(1,i)
                        unpacked(3*i-1) = packed(2,i)
                        unpacked(3*i  ) = packed(3,i)
                  enddo
      end function unpack_coords

            function pack_coords(unpacked) result(packed)
                  implicit none
                  real*8 :: unpacked(:)
                  real*8 :: packed(3,int(size(unpacked)/3))
                  integer :: i
                  do i=1,int(size(unpacked)/3)
                        packed(1,i) = unpacked(3*i-2)
                        packed(2,i) = unpacked(3*i-1)
                        packed(3,i) = unpacked(3*i  )
                  enddo
      end function pack_coords

            subroutine get_cm_coords(Natoms,xyz,cm_pos,xyz_cm)
                  implicit none
                  integer,intent(in)  :: Natoms
                  real*8 ,intent(in)  :: xyz(3,Natoms)
                  real*8 ,intent(out) :: cm_pos(3),xyz_cm(3,Natoms)
                  real*8              :: total_mass
                  integer             :: i,j
                  total_mass = 0.d0
                  cm_pos = 0.d0
                  xyz_cm = 0.d0
                  do i=1,Natoms
                        total_mass = total_mass + M(i)
                        cm_pos(:) = cm_pos(:) + M(i)*xyz(:,i)
                  end do
                  cm_pos = cm_pos/total_mass
                  do i=1,Natoms
                        xyz_cm(:,i) = xyz(:,i) - cm_pos(:)
                  end do
      end subroutine get_cm_coords

            function cart_to_normal(Natoms,xyz,Base) result(xyz_nm)
                  implicit none
                  integer :: Natoms
                  real*8  :: xyz(3,Natoms),Base(3*Natoms,3*Natoms),xyz_nm(3*Natoms)
                  real*8  :: xyz_mw(3,Natoms),xyz_mw_unpacked(3*Natoms)

                  xyz_mw(1,:) = xyz(1,:)*sqrt(M)
                  xyz_mw(2,:) = xyz(2,:)*sqrt(M)
                  xyz_mw(3,:) = xyz(3,:)*sqrt(M)

                  xyz_mw_unpacked = unpack_coords(xyz_mw)

                  xyz_nm = matmul(transpose(Base),xyz_mw_unpacked)
      end function cart_to_normal


            function build_eckart_matrix(Natoms,xyz_eq,xyz_cm) result(EM)
                  implicit none
                  integer :: Natoms
                  real*8  :: xyz_eq(3,Natoms),xyz_cm(3,Natoms)
                  real*8  :: xpa,xma,ypa,yma,zpa,zma
                  real*8  :: EM(4,4)
                  integer :: a
                  EM = 0.d0
                  do a=1,Natoms
                        xpa = xyz_eq(1,a) + xyz_cm(1,a)
                        xma = xyz_eq(1,a) - xyz_cm(1,a)

                        ypa = xyz_eq(2,a) + xyz_cm(2,a)
                        yma = xyz_eq(2,a) - xyz_cm(2,a)

                        zpa = xyz_eq(3,a) + xyz_cm(3,a)
                        zma = xyz_eq(3,a) - xyz_cm(3,a)

                        !Diagonal
                        EM(1,1) = EM(1,1) + M(a)*(xpa**2 + yma**2 + zma**2)
                        EM(2,2) = EM(2,2) + M(a)*(xma**2 + ypa**2 + zpa**2)
                        EM(3,3) = EM(3,3) + M(a)*(xpa**2 + yma**2 + zpa**2)
                        EM(4,4) = EM(4,4) + M(a)*(xpa**2 + ypa**2 + zma**2)

                        !Rest of 1st row
                        EM(1,2) = EM(1,2) + M(a)*(ypa*zma - yma*zpa)
                        EM(1,3) = EM(1,3) + M(a)*(xma*zpa - xpa*zma)
                        EM(1,4) = EM(1,4) + M(a)*(xpa*yma - xma*ypa)

                        !Rest of 2nd row
                        EM(2,3) = EM(2,3) + M(a)*(xma*yma - xpa*ypa)
                        EM(2,4) = EM(2,4) + M(a)*(xma*zma - xpa*zpa)

                        !Rest of 3rd row
                        EM(3,4) = EM(3,4) + M(a)*(yma*zma - ypa*zpa)
                  end do

                  !Symmetrize
                  EM(2,1) = EM(1,2)
                  EM(3,1) = EM(1,3)
                  EM(4,1) = EM(1,4)

                  EM(3,2) = EM(2,3)
                  EM(4,2) = EM(2,4)

                  EM(4,3) = EM(3,4)
      end function build_eckart_matrix

            function build_direction_cosine_matrix(V) result (U)
                  implicit none
                  real*8 :: V(4)
                  real*8 :: U(3,3)
                  real*8 :: q0,q1,q2,q3
                  q0 = V(1)
                  q1 = V(2)
                  q2 = V(3)
                  q3 = V(4)

                  !1st row
                  U(1,1) = q0**2 + q1**2 - q2**2 - q3**2
                  U(1,2) = 2.d0*(q1*q2 + q0*q3)
                  U(1,3) = 2.d0*(q1*q3 - q0*q2)
                  
                  !2nd row
                  U(2,1) = 2.d0*(q1*q2 - q0*q3)
                  U(2,2) = q0**2 - q1**2 + q2**2 - q3**2
                  U(2,3) = 2.d0*(q2*q3 + q0*q1)

                  !3rd row
                  U(3,1) = 2.d0*(q1*q3 + q0*q2)
                  U(3,2) = 2.d0*(q2*q3 - q0*q1)
                  U(3,3) = q0**2 - q1**2 - q2**2 + q3**2
      end function build_direction_cosine_matrix

            subroutine get_eckart_state(Natoms,xyz,xyz_eq,xyz_cm,xyz_eckart)
                  implicit none
                  integer,intent(in) :: Natoms
                  real*8,intent(in)  :: xyz(3,Natoms),xyz_eq(3,Natoms)
                  real*8,intent(out) :: xyz_cm(3,Natoms),xyz_eckart(3,Natoms)
                  real*8             :: cm_pos(3),EM(4,4),U(3,3)
                  real*8             :: work(100),d(4)
                  integer            :: i,ierror

                  call get_cm_coords(Natoms,xyz,cm_pos,xyz_cm)

                  EM = build_eckart_matrix(Natoms,xyz_eq,xyz_cm)

                  call dsyev("V","U",4,EM,4,d,work,100,ierror)

                  U = build_direction_cosine_matrix(EM(:,1))

                  do i=1,Natoms
                        xyz_eckart(:,i) = matmul(U,xyz_cm(:,i))
                  end do
      end subroutine get_eckart_state

            function build_normal_base(Natoms,xyz,M) result (Base)
                  implicit none
                  integer :: Natoms
                  real*8 :: xyz(3,Natoms),M(Natoms),Base(3*Natoms,3*Natoms)
                  real*8 :: a(3*Natoms),b(3*Natoms),c(3*Natoms)
                  integer :: i,j
                  
                  Base = 0.d0
                  do i=1,Natoms
                        !TX,TY,TZ
                        Base(3*i-2,1) = sqrt(M(i))
                        Base(3*i-1,2) = sqrt(M(i))
                        Base(3*i  ,3) = sqrt(M(i))

                        !RX
                        Base(3*i-1,4) = -xyz(3,i)*sqrt(M(i))
                        Base(3*i  ,4) =  xyz(2,i)*sqrt(M(i))
                        !RY
                        Base(3*i-2,5) = -xyz(3,i)*sqrt(M(i))
                        Base(3*i  ,5) =  xyz(1,i)*sqrt(M(i))
                        !RZ
                        Base(3*i-2,6) = -xyz(2,i)*sqrt(M(i))
                        Base(3*i-1,6) =  xyz(1,i)*sqrt(M(i))
                  enddo

                  do i=1,6
                        Base(:,i) = Base(:,i)/get_norm(base(:,i))
                  enddo

                  
                  do i=4,3*Natoms
                        if(i<7) then
                              a = Base(:,i)
                        else
                              a = 0.d0
                              a(i) = 1.d0
                              ! a = 1.d0
                              ! a = a/get_norm(a)
                        endif

                        ! do j=1,i-1
                        !       c = (sum(a*Base(:,j)))/sum(Base(:,j)**2)*Base(:,j)
                        ! enddo
                        ! b = a - c

                        b = a
                        do j=1,i-1
                              c = (sum(b*Base(:,j)))/sum(Base(:,j)**2)*Base(:,j)
                              b = b - c
                        enddo
                        

                        b = b/get_norm(b)
                        do j=1,i-1
                              print*,sum(b*Base(:,j))
                        enddo
                        print*,""
                        Base(:,i) = b
                  enddo
      end function build_normal_base

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