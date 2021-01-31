      module tools
            implicit none
            real,parameter :: pi=4.*atan(1.) 
            real,parameter :: cov(10) = (/0.37,0.,0.,0.,0.,0.77,0.75,0.73,0.,0./)

            contains

            subroutine parse_atomic_symbol_to_num(N,S,Z)
                  implicit none
                  integer :: N
                  character :: S(N)*2
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
                  enddo
            end

            function unit_cross(u1,u2) result(u3)
                  implicit none
                  real :: u1(3),u2(3)
                  real :: u3(3)
      
                  u3(1) = u1(2)*u2(3) - u1(3)*u2(2)
                  u3(2) = - u1(1)*u2(3) + u1(3)*u2(1)
                  u3(3) = u1(1)*u2(2) - u1(2)*u2(1)
      
                  u3 = u3/sqrt(sum(u3**2))
            end function unit_cross

            subroutine get_local_coords(c1,c2,c3,basis)
                  implicit none
                  real :: c1(3),c2(3),c3(3)
                  real :: basis(3,3) !First coordinates, second vectors
                  real :: u1(3),u2(3)
                  real :: n(3),nn(3)
      
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

            function get_angle(c1,c2,c3) result(A)
                  implicit none
                  real :: c1(3),c2(3),c3(3)
                  real :: u1(3),u2(3),proj
                  real :: A

                  u1 = c2-c1
                  u2 = c3-c2
      
                  u1 = u1/sqrt(sum(u1**2))
                  u2 = u2/sqrt(sum(u2**2))

                  proj = sum(u1*u2)
                  A = 180. - (180./pi)*acos(proj)
            end function get_angle

            function get_torsion(c1,c2,c3,c4) result(T)
                  implicit none
                  real :: c1(3),c2(3),c3(3),c4(3)
                  real :: u1(3),u2(3),u3(3),u4(3)
                  real :: u12(3),u34(3)
                  real :: proj,proj2
                  real :: T

                  u1 = c2-c1
                  u2 = c3-c2
                  u3 = c2-c3
                  u4 = c4-c3
      
                  u1 = u1/sqrt(sum(u1**2))
                  u2 = u2/sqrt(sum(u2**2))
                  u3 = u3/sqrt(sum(u3**2))
                  u4 = u4/sqrt(sum(u4**2))

                  u12 = unit_cross(u1,u2)
                  u34 = unit_cross(u3,u4)
                  
                  proj = sum(u12*u34)
                  proj2 = sum(u12*u4)

                  T = 180. - (180./pi)*sign(1.,proj2)*acos(proj)
            end function get_torsion

            subroutine get_xyz(port,N,S,xyz)
                  implicit none
                  integer :: port,N
                  character,allocatable :: S(:)*2
                  real,allocatable :: xyz(:,:)
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
                  real :: xyz(3,N)
                  real,allocatable :: dmat(:,:)
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

            function get_bond_graph(N,S,dmat) result(bond_graph)
                  implicit none
                  integer :: N
                  character :: S(N)*2
                  real :: dmat(N,N)
                  logical,allocatable :: bond_graph(:,:)
                  real,parameter :: thresh = 1.2
                  integer :: Z(N)
                  integer :: counter,i,j

                  call parse_atomic_symbol_to_num(N,S,Z)
                  
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

            subroutine get_bonds(N,S,dmat,bond_pairsbis,bond_valsbis)
                  implicit none
                  integer :: N
                  character :: S(N)*2
                  real :: dmat(N,N)
                  integer,allocatable :: bond_pairs(:,:),bond_pairsbis(:,:)
                  real,allocatable :: bond_vals(:),bond_valsbis(:)
                  real,parameter :: thresh = 1.2
                  integer :: Z(N)
                  integer :: Nbonds,i,j

                  call parse_atomic_symbol_to_num(N,S,Z)
                  
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

                  allocate(bond_pairsbis(Nbonds,2),bond_valsbis(Nbonds))
                  bond_pairsbis = bond_pairs(:Nbonds,:)
                  bond_valsbis = bond_vals(:Nbonds)
            end subroutine get_bonds

            subroutine get_angles(N,xyz,bond_graph,angle_pairsbis,angle_valsbis)
                  implicit none
                  integer :: N
                  logical :: bond_graph(N,N)
                  real :: xyz(3,N)
                  integer,allocatable :: angle_pairs(:,:),angle_pairsbis(:,:)
                  real,allocatable :: angle_vals(:),angle_valsbis(:)
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
                  allocate(angle_pairsbis(Nangles,3),angle_valsbis(Nangles))
                  angle_pairsbis = angle_pairs(:Nangles,:)
                  angle_valsbis = angle_vals(:Nangles)
            end subroutine get_angles

            subroutine get_torsions(N,xyz,bond_graph,torsion_pairsbis,torsion_valsbis)
                  implicit none
                  integer :: N
                  logical :: bond_graph(N,N)
                  real :: xyz(3,N)
                  integer,allocatable :: torsion_pairs(:,:),torsion_pairsbis(:,:)
                  real,allocatable :: torsion_vals(:),torsion_valsbis(:)
                  integer :: i,j,k,l
                  integer :: a,b,c,d
                  integer :: Ntorsions

                  allocate(torsion_pairs(N*N,4),torsion_vals(N*N))
                  Ntorsions = 0
                  do i = 1,N
                        a = i
                        do j = 1,N
                              if(bond_graph(a,j).eqv..true.) then
                                    if(j>a) then
                                          cycle
                                    endif
                                    b = j
                                    do k = 1,N
                                          if(bond_graph(b,k).eqv..true.) then
                                                if(k==a) then
                                                      cycle
                                                endif
                                                c = k
                                                do l=1,N
                                                      if(bond_graph(c,l).eqv..true.) then
                                                            if(l == a .or. l == b) then
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
                  allocate(torsion_pairsbis(Ntorsions,4),torsion_valsbis(Ntorsions))
                  torsion_pairsbis = torsion_pairs(:Ntorsions,:)
                  torsion_valsbis = torsion_vals(:Ntorsions)
            end subroutine get_torsions

            subroutine save_zmat(port,N,Nangle,Ntorsion,S,dmat,bond_graph,angle_pairs,angle_vals,torsion_pairs,torsion_vals)
                  implicit none
                  integer :: port,N,Nangle,Ntorsion
                  character :: S(N)*2
                  integer :: angle_pairs(Nangle,3),torsion_pairs(Ntorsion,4)
                  real :: dmat(N,N),angle_vals(Nangle),torsion_vals(Ntorsion)
                  logical :: bond_graph(N,N)

                  integer :: Batom,Aatom,Tatom
                  real :: B,A,T
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