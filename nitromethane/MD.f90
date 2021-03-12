      use tools
      use ff_utils
      use md_utils
      implicit none
      integer :: Natoms,Nbonds,Nangles,Ntorsions
      real*8,allocatable :: dmat(:,:)
      logical,allocatable :: bond_graph(:,:)
      integer,allocatable :: bond_pairs(:,:),angle_pairs(:,:),torsion_pairs(:,:)
      real*8,allocatable :: bond_vals(:),angle_vals(:),torsion_vals(:)
      real*8,allocatable :: G(:,:)
      real*8,allocatable :: H(:,:),Hm(:,:),Hmcopy(:,:)
      real*8,allocatable :: d(:)
      real*8 :: work(100)
      real*8 :: E,K,T
      real*8,allocatable :: NM_energies(:)
      character :: format_label*90

      real*8,allocatable :: v(:,:)
      real*8,parameter :: dt=1d-3 !0.1 fs
      integer,parameter :: Nt = 10000 !With 0.1 fs dt it should be 1 ps

      real*8 :: cm_pos(3)
      real*8,allocatable :: xyz_cm(:,:),xyz_eq(:,:),xyz_eckart(:,:)

      integer :: i,j,a,b,p,q

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
      
      ! v = 0.d0

      ! print*,E,sum(v*v),sum(G*G)

      ! i = 0
      ! do while(sum(G*G)>1.d-14) !1.d-10
      !       i = i + 1
      !       call verletvel_step(Natoms,dt,M,xyz,v,G,Nbonds,Nangles,Ntorsions,&
      !       bond_pairs,angle_pairs,torsion_pairs)

      !       v = 0.99d0*v
            
      !       bond_vals = recomp_bonds(Natoms,Nbonds,xyz,bond_pairs)
      !       angle_vals = recomp_angles(Natoms,Nangles,xyz,angle_pairs)
      !       torsion_vals = recomp_torsions(Natoms,Ntorsions,xyz,torsion_pairs)
      !       improper_vals = recomp_impropers(Natoms,Nimpropers,xyz,improper_pairs)

      !       E = comp_energy(Nbonds,Nangles,Ntorsions,bond_vals,&
      !       angle_vals,torsion_vals)


      !       if(mod(i,100)==0) then
      !             print*,E,sum(v*v),sum(G*G)

      !             call write_conf(3,Natoms,S,xyz,3)
      !       endif
      ! enddo

      ! bond_vals = recomp_bonds(Natoms,Nbonds,xyz,bond_pairs)
      ! angle_vals = recomp_angles(Natoms,Nangles,xyz,angle_pairs)
      ! torsion_vals = recomp_torsions(Natoms,Ntorsions,xyz,torsion_pairs)
      ! improper_vals = recomp_impropers(Natoms,Nimpropers,xyz,improper_pairs)

      ! print*,"Final energy and gradient:",comp_energy(Nbonds,Nangles,Ntorsions,bond_vals,&
      ! angle_vals,torsion_vals),sum(G*G)


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

      allocate(d(3*Natoms))


      Hmcopy = Hm
      call dsyev("V","U",3*Natoms,Hm,3*Natoms,d,work,100,i)

      open(3,file="results/trajectory.xyz")
      open(4,file="results/thermodynamics.dat")
      open(5,file="results/nm_ener.dat")
      open(10,file="results/eckart_states.xyz")
      write(format_label,"(A,I2,A)")"(I5,2X,",3*Natoms,"(E14.7,2X))"


      T = 300.d0*R_kJ_mol_K !300K ~ 2.5 kJ/mol
      v = init_normalmodes_micro(Natoms,Hm,T)

      allocate(G(3,Natoms),NM_energies(3*Natoms))
      allocate(xyz_cm(3,Natoms),xyz_eq(3,Natoms),xyz_eckart(3,Natoms))

      G = -build_gradient(Natoms,xyz,Nbonds,Nangles,Ntorsions,&
      bond_pairs,angle_pairs,torsion_pairs)

      E = comp_energy(Nbonds,Nangles,Ntorsions,bond_vals,&
      angle_vals,torsion_vals)      

      K = comp_kinetic_energy(v,0)

      call write_conf(3,Natoms,xyz,3)

      write(4,"(I5,2X,4(E14.7,2X))")0,K,E,K+E,sum(NM_energies)


      call get_cm_coords(Natoms,xyz,cm_pos,xyz_eq)

      call get_eckart_state(Natoms,xyz,xyz_eq,xyz_cm,xyz_eckart)
      call write_conf(3,Natoms,xyz_eckart,10)

      NM_energies = comp_normal_energy(Natoms,xyz_cm,xyz_eckart,v,Hm,d)
      write(5,format_label)0,NM_energies

      do i=1,Nt
            call verletvel_step(Natoms,dt,xyz,v,G,Nbonds,Nangles,Ntorsions,&
            bond_pairs,angle_pairs,torsion_pairs)
            
            bond_vals = recomp_bonds(Natoms,Nbonds,xyz,bond_pairs)
            angle_vals = recomp_angles(Natoms,Nangles,xyz,angle_pairs)
            torsion_vals = recomp_torsions(Natoms,Ntorsions,xyz,torsion_pairs)
            improper_vals = recomp_impropers(Natoms,Nimpropers,xyz,improper_pairs)

            E = comp_energy(Nbonds,Nangles,Ntorsions,bond_vals,&
            angle_vals,torsion_vals)

            K = comp_kinetic_energy(v,0)
            NM_energies = comp_normal_energy(Natoms,xyz_cm,xyz_eckart,v,Hm,d)

            call write_conf(3,Natoms,xyz,3)

            call get_eckart_state(Natoms,xyz,xyz_eq,xyz_cm,xyz_eckart)
            call write_conf(3,Natoms,xyz_eckart,10)

            write(4,"(I5,2X,4(E14.7,2X))")i,K,E,K+E,sum(NM_energies)
            write(5,format_label)i,NM_energies
      enddo

      end




      subroutine verletvel_step(Natoms,dt,r,v,f,Nbonds,Nangles,Ntorsions,&
            bond_pairs,angle_pairs,torsion_pairs)
            use tools
            use ff_utils
            implicit none
            integer :: Natoms,Nbonds,Nangles,Ntorsions
            integer :: bond_pairs(Nbonds,2),angle_pairs(Nangles,3),torsion_pairs(Ntorsions,4)
            real*8 :: dt,r(3,Natoms),v(3,Natoms),f(3,Natoms)
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