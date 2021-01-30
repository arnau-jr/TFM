      use tools
      implicit none

      integer :: N,Nangle,Ntorsion
      character,allocatable :: S(:)*2
      real,allocatable :: xyz(:,:),dmat(:,:)
      integer,allocatable :: bond_graph(:,:)
      real,allocatable :: angles(:,:),torsions(:,:)
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


      call get_xyz(1,N,S,xyz)

      dmat =  get_dist_matrix(N,xyz)
      bond_graph = get_bond_graph(N,dmat)

      angles = get_angles(N,xyz,bond_graph)
      Nangle = size(angles,1)

      torsions = get_torsions(N,xyz,bond_graph)
      Ntorsion = size(torsions,1)

      call save_zmat(2,N,Nangle,Ntorsion,S,dmat,bond_graph,angles,torsions)

      end

