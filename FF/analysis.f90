      use tools
      use ff_utils
      implicit none

      integer :: N,Nbond,Nangle,Ntorsion
      character,allocatable :: S(:)*2
      real,allocatable :: xyz(:,:),dmat(:,:)
      logical,allocatable :: bond_graph(:,:)
      integer,allocatable :: bond_pairs(:,:),angle_pairs(:,:),torsion_pairs(:,:)
      real,allocatable :: bond_vals(:),angle_vals(:),torsion_vals(:)
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
      bond_graph = get_bond_graph(N,S,dmat)

      call get_bonds(N,S,dmat,bond_pairs,bond_vals)
      Nbond = size(bond_vals)
      print*,Nbond
      do i=1,Nbond
            print*,bond_pairs(i,:),bond_vals(i)
      enddo
      print*,""


      call get_angles(N,xyz,bond_graph,angle_pairs,angle_vals)
      Nangle = size(angle_vals)
      print*,Nangle
      do i=1,Nangle
            print*,angle_pairs(i,:),angle_vals(i)
      enddo
      print*,""


      call get_torsions(N,xyz,bond_graph,torsion_pairs,torsion_vals)
      Ntorsion = size(torsion_vals)

      print*,Ntorsion
      do i=1,Ntorsion
            print*,torsion_pairs(i,:),torsion_vals(i)
      enddo

      end