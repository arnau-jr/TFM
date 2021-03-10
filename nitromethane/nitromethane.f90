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
      integer :: i,j,k,a,b,p,q,nrot
      
      character :: input_filename*90,param_filename*90,output_filename*90

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

      ! print*,Nbonds
      ! do i=1,Nbonds
      !       print*,bond_pairs(i,:),bond_vals(i)
      ! enddo
      ! print*,""
      ! print*,Nangles
      ! do i=1,Nangles
      !       print*,angle_pairs(i,:),angle_vals(i)
      ! enddo
      ! print*,""
      ! print*,Ntorsions
      ! do i=1,Ntorsions
      !       print*,torsion_pairs(i,:),torsion_vals(i)
      ! enddo
      

      call get_param(2,Nbonds,Nangles,Ntorsions)

      print*,Nimpropers
      do i=1,Nimpropers
            print*,improper_pairs(i,:),improper_vals(i)
      enddo
      
      allocate(H(3*Natoms,3*Natoms),Hm(3*Natoms,3*Natoms),Hmcopy(3*Natoms,3*Natoms))
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
      print*,""


      allocate(G(3,Natoms))
      G = build_gradient(Natoms,xyz,Nbonds,Nangles,Ntorsions,&
      bond_pairs,angle_pairs,torsion_pairs)

      print"(A,E14.7,A)","Total gradient: ",sum(G*G)," (kJ/mol/A)^2"
      do i=1,Natoms
            print"(3(F14.7,2X))",G(:,i)
      enddo

      print*,""
      print*,"E = ",comp_energy(Nbonds,Nangles,Ntorsions,bond_vals,&
      angle_vals,torsion_vals),"kJ/mol"
      print*,"Bonds: ",comp_bonds_energy(Nbonds,bond_vals),-sum(De)
      print*,"Angles: ",comp_angles_energy(Nangles,angle_vals)
      print*,"Torsions: ",comp_torsions_energy(Ntorsions,torsion_vals)
      print*,"Impropers: ",comp_impropers_energy()

      allocate(d(3*Natoms),v(3*Natoms,3*Natoms))
      ! call jacobi(Hm,1,3*Natoms,d,v,nrot)
      ! call sort_ev(d,v,3*Natoms)


      Hmcopy = Hm
      call dsyev("V","U",3*Natoms,Hm,3*Natoms,d,work,100,i)

      ! print*,"Jacobi finished, took",nrot,"rotations"
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

      ! j = 8
      ! print*,"Eigenvector",j
      ! do i=1,Natoms
      !       print"(2(F20.15,2X))",Hm(3*i-2,j)
      !       print"(2(F20.15,2X))",Hm(3*i-1,j),sqrt(M(i))
      !       print"(2(F20.15,2X))",Hm(3*i,j)
      !       print*,""
      ! enddo

      ! print*,sum(matmul(Hmcopy,Hm(:,1))/Hm(:,1))/(3.d0*Natoms)


      close(1)
      do k=7,3*Natoms
            print*,"Testing normal mode",k
            if(k<10) then
                  write(output_filename,"(A,I1,A)")"normal_oscillations/normalmode_",k,".xyz"
            elseif(k<100) then
                  write(output_filename,"(A,I2,A)")"normal_oscillations/normalmode_",k,".xyz"
            endif
            open(1,file=output_filename)

            call test_normal_mode(Natoms,xyz,Hm,k,sqrt(d(k)),1)
            close(1)
      enddo
      


      end