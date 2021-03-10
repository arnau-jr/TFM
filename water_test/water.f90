      use tools
      use ff_utils
      implicit none

      integer :: Natoms,Nbonds,Nangles,Ntorsions
      real*8,allocatable :: dmat(:,:)
      logical,allocatable :: bond_graph(:,:)
      integer,allocatable :: bond_pairs(:,:),angle_pairs(:,:),torsion_pairs(:,:)
      real*8,allocatable :: bond_vals(:),angle_vals(:),torsion_vals(:)
      real*8,allocatable :: H(:,:),Hm(:,:),Hmcopy(:,:),G(:,:)
      real*8,allocatable :: d(:),v(:,:)
      real*8 :: work(100)

      real*8,allocatable :: norm_osc(:),cart_osc(:),xyz_osc(:,:),NBase(:,:)
      integer :: Nosc
      real*8 :: dosc,omega

      integer :: i,j,k,a,b,p,q
      
      character :: input_filename*90,param_filename*90

      call get_command_argument(1,input_filename)
      call get_command_argument(2,param_filename)

      if(index(input_filename,".xyz")==0) then
            print*, "Error : input file does not have xyz extension"
            stop
      endif


      open(1,file=input_filename)
      open(2,file=param_filename)


      call get_xyz(1,Natoms)

      dmat =  get_dist_matrix(Natoms,xyz)
      bond_graph = get_bond_graph(Natoms,dmat)

      call get_bonds(Natoms,dmat,bond_pairs,bond_vals)
      Nbonds = size(bond_vals)


      call get_angles(Natoms,xyz,bond_graph,angle_pairs,angle_vals)
      Nangles = size(angle_vals)


      call get_torsions(Natoms,xyz,bond_graph,torsion_pairs,torsion_vals)
      Ntorsions = size(torsion_vals)
      

      call get_param(2,Nbonds,Nangles,Ntorsions)

      close(1)
      close(2)
      
      allocate(H(3*Natoms,3*Natoms),Hm(3*Natoms,3*Natoms))
      H = build_hessian(Natoms,xyz,Nbonds,Nangles,Ntorsions,&
            bond_pairs,angle_pairs,torsion_pairs)
      
      do a=1,Natoms
            do b=1,Natoms
                  do p=0,2
                        do q=0,2
                              i = 3*a - p
                              j = 3*b - q

                              Hm(i,j) = H(i,j)/sqrt(M(a)*M(b))
                        enddo
                  enddo
            enddo
      enddo

      do i=1,3*Natoms
            ! print"(20(F7.2,2X))",Hm(i,:20)
            ! print"(9(F14.7,2X))",H(i,:)
      enddo
      ! print*,""


      allocate(G(3,Natoms))
      G = build_gradient(Natoms,xyz,Nbonds,Nangles,Ntorsions,&
      bond_pairs,angle_pairs,torsion_pairs)
      print"(A,E14.7,A)","Total gradient: ",sum(G*G)," (kJ/mol/A)^2"

      print*,""
      print*,"E = ",comp_energy(Nbonds,Nangles,Ntorsions,bond_vals,&
      angle_vals,torsion_vals),"kJ/mol"
      print*,"Bonds: ",comp_bonds_energy(Nbonds,bond_vals)
      print*,"Angles: ",comp_angles_energy(Nangles,angle_vals)
      print*,"Torsions: ",comp_torsions_energy(Ntorsions,torsion_vals)
      print*,"Impropers: ",comp_impropers_energy()


      allocate(d(3*Natoms),v(3*Natoms,3*Natoms))
      Hmcopy = Hm
      call dsyev("V","U",3*Natoms,Hm,3*Natoms,d,work,100,i)

      print*,"Eigenvalues"
      do i=1,6
            print"(F20.15)",d(i)
      enddo
      print*,""
      do i=7,3*Natoms
            print"(F20.15)",d(i)
      enddo

      print*,"Energy"
      do i=7,3*Natoms
            print"(F20.12,2X,I4)",sqrt(d(i))*hbar_cm_dps,nint(sqrt(d(i))*hbar_cm_dps)
      enddo
      print*,""

      ! j = 1
      ! do j = 1,3*Natoms
      ! ! do j = j,j
      ! print*,"Eigenvector",j
      ! do i=1,Natoms
      !       print"(2(F20.15,2X))",Hm(3*i-2,j)
      !       print"(2(F20.15,2X))",Hm(3*i-1,j)
      !       print"(2(F20.15,2X))",Hm(3*i,j)
      !       print*,""
      ! enddo
      ! print*,""
      ! enddo

      allocate(norm_osc(3*Natoms),cart_osc(3*Natoms),xyz_osc(3,Natoms),NBase(3*Natoms,3*Natoms))

      norm_osc = 0.d0
      k = 9
      omega = sqrt(d(k))
      dosc = 0.01/omega
      Nosc = nint(((2.d0*pi)/omega)/dosc)

      print*,"NM frequency:",omega
      open(1,file="normal_oscillation.xyz")
      do i=1,Nosc

            norm_osc(k) = 0.5*sin(omega*dosc*(i-1))

            cart_osc = matmul(Hm,norm_osc)
            
            do j=1,Natoms
                  xyz_osc(1,j) = xyz(1,j) + cart_osc(3*j-2)/sqrt(M(j))
                  xyz_osc(2,j) = xyz(2,j) + cart_osc(3*j-1)/sqrt(M(j))
                  xyz_osc(3,j) = xyz(3,j) + cart_osc(3*j)/sqrt(M(j))
            enddo
            call write_conf(3,Natoms,xyz_osc,1)
      enddo
      print*,""

      NBase = build_normal_base(Natoms,xyz,M)

      do i=1,3*Natoms
            print"(9(F14.7,X),A,F14.7)",Nbase(:,i),"|",sum(Nbase(:,i)**2)
      enddo

      end