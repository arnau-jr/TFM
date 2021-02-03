      use tools
      use ff_utils
      implicit none

      integer :: Natoms,Nbonds,Nangles,Ntorsions
      character,allocatable :: S(:)*2
      real,allocatable :: xyz(:,:),dmat(:,:)
      logical,allocatable :: bond_graph(:,:)
      integer,allocatable :: bond_pairs(:,:),angle_pairs(:,:),torsion_pairs(:,:)
      real,allocatable :: bond_vals(:),angle_vals(:),torsion_vals(:)
      real,allocatable :: req(:),kb(:)
      real,allocatable :: aeq(:),ka(:)
      real,allocatable :: An(:),n(:),delta(:)
      real*8,allocatable :: H(:,:)
      integer :: i
      
      character :: input_filename*90,output_filename*90

      call get_command_argument(1,input_filename)

      if(index(input_filename,".xyz")==0) then
            print*, "Error : input file does not have xyz extension"
            stop
      endif
      if(command_argument_count() == 2) then
            call get_command_argument(2,output_filename)
      else
            write(output_filename,"(A,A)")input_filename(:index(input_filename,".gzmat")),"xyz"  
      endif


      open(1,file=input_filename)
      open(2,file=output_filename)


      call get_xyz(1,Natoms,S,xyz)

      dmat =  get_dist_matrix(Natoms,xyz)
      bond_graph = get_bond_graph(Natoms,S,dmat)

      call get_bonds(Natoms,S,dmat,bond_pairs,bond_vals)
      Nbonds = size(bond_vals)


      call get_angles(Natoms,xyz,bond_graph,angle_pairs,angle_vals)
      Nangles = size(angle_vals)


      call get_torsions(Natoms,xyz,bond_graph,torsion_pairs,torsion_vals)
      Ntorsions = size(torsion_vals)

      allocate(req(Nbonds),kb(Nbonds))
      allocate(aeq(Nangles),ka(Nangles))
      allocate(An(Ntorsions),n(Ntorsions),delta(Ntorsions))

      req = 1.4
      kb = 1.
      aeq = 0.
      ka = 1.
      An = 10.
      n = 1.
      delta = 0.
      print*,comp_energy(Nbonds,Nangles,Ntorsions,bond_vals,&
            angle_vals,torsion_vals,req,kb,aeq,ka,An,n,delta)
      
      allocate(H(3*Natoms,3*Natoms))
      H = build_hessian(Natoms,xyz,Nbonds,Nangles,Ntorsions,&
            bond_pairs,angle_pairs,torsion_pairs,req,kb,aeq,ka,An,n,delta)

      do i=1,3*Natoms
            print*,H(i,i)
      enddo

      end