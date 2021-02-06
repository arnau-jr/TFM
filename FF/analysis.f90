      use tools
      use ff_utils
      implicit none

      integer :: Natoms,Nbonds,Nangles,Ntorsions
      character,allocatable :: S(:)*2
      real*8,allocatable :: xyz(:,:),dmat(:,:)
      logical,allocatable :: bond_graph(:,:)
      integer,allocatable :: bond_pairs(:,:),angle_pairs(:,:),torsion_pairs(:,:)
      real*8,allocatable :: bond_vals(:),angle_vals(:),torsion_vals(:)
      ! real*8,allocatable :: req(:),kb(:)
      ! real*8,allocatable :: aeq(:),ka(:)
      ! real*8,allocatable :: An(:),n(:),delta(:)
      real*8,allocatable :: H(:,:),G(:,:)
      integer :: i
      
      character :: input_filename*90,output_filename*90

      call get_command_argument(1,input_filename)

      if(index(input_filename,".xyz")==0) then
            print*, "Error : input file does not have xyz extension"
            stop
      endif


      open(1,file=input_filename)
      open(2,file="param.dat")


      call get_xyz(1,Natoms,S,xyz)

      dmat =  get_dist_matrix(Natoms,xyz)
      bond_graph = get_bond_graph(Natoms,S,dmat)

      call get_bonds(Natoms,S,dmat,bond_pairs,bond_vals)
      Nbonds = size(bond_vals)


      call get_angles(Natoms,xyz,bond_graph,angle_pairs,angle_vals)
      Nangles = size(angle_vals)


      call get_torsions(Natoms,xyz,bond_graph,torsion_pairs,torsion_vals)
      Ntorsions = size(torsion_vals)
      

      call get_param(2,Nbonds,Nangles,Ntorsions)

      print*,comp_energy(Nbonds,Nangles,Ntorsions,bond_vals,&
            angle_vals,torsion_vals)
      
      allocate(H(3*Natoms,3*Natoms))
      H = build_hessian(Natoms,xyz,Nbonds,Nangles,Ntorsions,&
            bond_pairs,angle_pairs,torsion_pairs)

      allocate(G(3,Natoms))
      G = build_gradient(Natoms,xyz,Nbonds,Nangles,Ntorsions,&
      bond_pairs,angle_pairs,torsion_pairs)
      
      do i=1,3*Natoms
            print"(20(F7.2,2X))",H(i,:20)
      enddo

      end