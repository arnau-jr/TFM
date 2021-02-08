      use tools
      use ff_utils
      implicit none
      integer :: Natoms,Nbonds,Nangles,Ntorsions
      character,allocatable :: S(:)*2
      integer,allocatable :: Z(:)
      real*8,allocatable :: M(:),xyz(:,:),dmat(:,:)
      logical,allocatable :: bond_graph(:,:)
      integer,allocatable :: bond_pairs(:,:),angle_pairs(:,:),torsion_pairs(:,:)
      real*8,allocatable :: bond_vals(:),angle_vals(:),torsion_vals(:)
      real*8,allocatable :: G(:,:)

      real*8,allocatable :: v(:,:)
      real*8,parameter :: dt=0.01d0
      ! integer,parameter :: Nt = 10000

      integer :: i

      character :: input_filename*90,output_filename*90

      call get_command_argument(1,input_filename)

      if(index(input_filename,".xyz")==0) then
            print*, "Error : input file does not have xyz extension"
            stop
      endif


      open(1,file=input_filename)
      open(2,file="param.dat")
      open(3,file="trajectory.xyz")


      call get_xyz(1,Natoms,S,xyz)
      allocate(Z(Natoms),M(Natoms),v(3,Natoms))
      call parse_atomic_symbol(Natoms,S,Z,M)


      dmat =  get_dist_matrix(Natoms,xyz)
      bond_graph = get_bond_graph(Natoms,Z,dmat)

      call get_bonds(Natoms,Z,dmat,bond_pairs,bond_vals)
      Nbonds = size(bond_vals)


      call get_angles(Natoms,xyz,bond_graph,angle_pairs,angle_vals)
      Nangles = size(angle_vals)


      call get_torsions(Natoms,xyz,bond_graph,torsion_pairs,torsion_vals)
      Ntorsions = size(torsion_vals)
      

      call get_param(2,Nbonds,Nangles,Ntorsions)
      
      allocate(G(3,Natoms))
      G = -build_gradient(Natoms,xyz,Nbonds,Nangles,Ntorsions,&
      bond_pairs,angle_pairs,torsion_pairs)

      print*,sum(G*G)

      call write_conf(3,Natoms,S,xyz,3)
      v = 0.d0

      i = 0
      do while(sum(G*G)>2.5d-5) !2.5d-5
            i = i + 1
            call verletvel_step(Natoms,dt,M,xyz,v,G,Nbonds,Nangles,Ntorsions,&
            bond_pairs,angle_pairs,torsion_pairs)

            if(sum(G*G)>5.d-5) then
                  v = 0.9*v
            else
                  v = 0.8*v
            endif
            if(mod(i,100)==0) then
                  print*,sum(G*G)

                  call write_conf(3,Natoms,S,xyz,3)
            endif
      enddo

      close(3)
      open(3,file="DPNA_opti.xyz")
      call write_conf(3,Natoms,S,xyz,3)

      end




      subroutine verletvel_step(Natoms,dt,M,r,v,f,Nbonds,Nangles,Ntorsions,&
            bond_pairs,angle_pairs,torsion_pairs)
            use tools
            use ff_utils
            implicit none
            integer :: Natoms,Nbonds,Nangles,Ntorsions
            integer :: bond_pairs(Nbonds,2),angle_pairs(Nangles,3),torsion_pairs(Ntorsions,4)
            real*8 :: dt,M(Natoms),r(3,Natoms),v(3,Natoms),f(3,Natoms)
            real*8 :: fold(3,Natoms)
            integer :: i

            fold = f
            do i=1,Natoms
                  r(:,i) = r(:,i) + v(:,i)*dt + 0.5d0*f(:,i)*dt**2/M(i)
            enddo

            f = -build_gradient(Natoms,r,Nbonds,Nangles,Ntorsions,&
            bond_pairs,angle_pairs,torsion_pairs)

            do i=1,Natoms
                  v(:,i) = v(:,i) + 0.5d0*(fold(:,i)+f(:,i))*dt/M(i)
            enddo
      end

      subroutine write_conf(D,N,S,r,port)
            implicit none
            integer :: N,D,port,i
            character :: S(N)*2
            real*8 :: r(D,N)

            write(port,"(I5)")N
            write(port,*)""
            do i=1,N
                  write(port,"(A,2X,E14.7,2X,E14.7,2X,E14.7)")S(i),r(:,i)
            enddo
      end