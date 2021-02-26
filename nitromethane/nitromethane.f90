      use tools
      use ff_utils
      implicit none

      integer :: Natoms,Nbonds,Nangles,Ntorsions
      character,allocatable :: S(:)*2
      integer,allocatable :: Z(:)
      real*8,allocatable :: M(:),dmat(:,:)
      logical,allocatable :: bond_graph(:,:)
      integer,allocatable :: bond_pairs(:,:),angle_pairs(:,:),torsion_pairs(:,:)
      real*8,allocatable :: bond_vals(:),angle_vals(:),torsion_vals(:)
      real*8,allocatable :: H(:,:),Hm(:,:),G(:,:)
      real*8,allocatable :: d(:),v(:,:)
      real*8 :: work(100)
      real*8 :: mu,freq2
      integer :: i,j,a,b,p,q,nrot
      
      character :: input_filename*90,param_filename*90

      call get_command_argument(1,input_filename)
      call get_command_argument(2,param_filename)

      if(index(input_filename,".xyz")==0) then
            print*, "Error : input file does not have xyz extension"
            stop
      endif


      open(1,file=input_filename)
      open(2,file=param_filename)


      call get_xyz(1,Natoms,S,xyz)
      allocate(Z(Natoms),M(Natoms))
      call parse_atomic_symbol(Natoms,S,Z,M)

      dmat =  get_dist_matrix(Natoms,xyz)
      bond_graph = get_bond_graph(Natoms,Z,dmat)

      call get_bonds(Natoms,Z,dmat,bond_pairs,bond_vals)
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
      print*,Ntorsions
      do i=1,Ntorsions
            print*,torsion_pairs(i,:),torsion_vals(i)
      enddo
      

      call get_param(2,Nbonds,Nangles,Ntorsions)

      
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

      call dsyev("N","U",3*Natoms,Hm,3*Natoms,d,work,100,i)

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

      print*,""

      end

            subroutine jacobi(a,n,np,d,v,nrot)
!     Finds the eigenvalues and eigenvectors of the symmetric matrix a using the Jacobi diagonalization method.
!           Input
!           -----
!                 a : real*8,dimension(n,n)
!                       Matrix to diagonalize, the upper triangle will be set to 0 during the diagonalization process but the
!                       original diagonal and lower triangle will not be changed.
!                 n : integer
!                       Box parameter, if 1 all elements are eliminated, if different than 1 then the the elements of
!                       of the first np-n box are not eliminated.                              
!                 np : integer
!                       Dimension of the matrix to diagonalize.
!           Output
!           ------
!                 d : real*8,dimension(n)
!                       Array with the (not sorted) eigenvalues of the input matrix.
!                 v : real*8,dimension(n,n)
!                       Array with the eignevectors (not normalized) arranged by columns 
!                       with the same ordering as the eigenvalue array.
!                 nrot : integer
!                       Number of Jacobi rotations used to diagonalize the matrix.
            implicit none
            integer :: n,np
            real*8 :: a(np,np),d(np),v(np,np)
            integer :: nrot
            real*8 :: maxelement,theta,t,c,s,tau,sum
            real*8 :: ap,aq,vp,vq,diff
            integer,parameter :: maxiter = 1e6
            real*8,parameter :: eps = 1.d-12
            integer :: i,j,maxi,maxj,p,q

            d = 0.d0 !Initialize d to all zeros
            
            !Initialize v to identity
            do i=1,np
                  do j=1,np
                        if(i.eq.j) then
                              v(i,j) = 1.d0
                        else
                              v(i,j) = 0.d0
                        endif
                  enddo
            enddo

            !Convergence checks
            sum = 0.d0
            maxelement = 0.d0

            !Find greatest element in the upper triangle, element a(p,q)
            !Because it is the upper triangle q => p always
            if(n == 1) then
                  do i=1,np
                        do j=i+1,np
                              sum = sum + a(i,j)**2
                              if(maxelement < abs(a(i,j))) then
                                    maxelement = abs(a(i,j))
                                    maxi = i
                                    maxj = j
                              endif
                        enddo
                  enddo
            else
                  do i=1,np-n
                        do j=n+1,np
                              sum = sum + a(i,j)**2
                              if(maxelement < abs(a(i,j))) then
                                    maxelement = abs(a(i,j))
                                    maxi = i
                                    maxj = j
                              endif
                        enddo
                  enddo
            endif
            p = maxi
            q = maxj

            nrot = 0
            do while(sum > eps .and. maxelement > eps) !Main loop
                  !Compute all the relevant quantities
                  diff = a(q,q) - a(p,p)
                  !Machine epsilon check for theta, within 100 times(because the square)
                  if(abs(diff)+100.d0*a(p,q) .eq. abs(diff)) then
                        t = a(p,q)/diff
                  else
                        theta = (a(q,q) - a(p,p))/(2.d0*a(p,q))
                        t = sign(1.d0,theta)*(1.d0/(abs(theta) + sqrt(theta**2 + 1.d0)))
                  endif
                  c = 1.d0/sqrt((1.d0 + t**2))
                  s = t*c
                  tau = s/(1.d0 + c)

                  !Update app aqq and apq
                  a(q,q) = a(q,q) + t*a(p,q)
                  a(p,p) = a(p,p) - t*a(p,q)
                  a(p,q) = 0.d0

                  do i=1,np !Loop over all rows/columns
                        if(i<p) then !For all elements below the row/column p
                              !Store old elements of the column p and q
                              ap = a(i,p) 
                              aq = a(i,q)
                              !Modify all elements of the column p and q
                              a(i,p) = ap - s*(aq + ap*tau) 
                              a(i,q) = aq + s*(ap - aq*tau)
                        elseif(p<i .and. i<q) then !For all elements above p but below q (remember q => p always)
                              !Store old elements for the row p and column q
                              ap = a(p,i)
                              aq = a(i,q)
                              !Modify all elements of row p and column q
                              a(p,i) = ap - s*(aq + ap*tau)
                              a(i,q) = aq + s*(ap - aq*tau)
                        elseif(i>q) then !For the remaining elements (above q)
                              !Store old rows p and q
                              ap = a(p,i)
                              aq = a(q,i)
                              !Modify row p and q
                              a(p,i) = ap - s*(aq + ap*tau)
                              a(q,i) = aq + s*(ap - aq*tau)
                        endif
                        !Note that in the previous loop we have not modified elements pp,qq or pq, which are already modified

                        !Update eigenvector matrix
                        vp = v(i,p)
                        vq = v(i,q)
                        v(i,p) = vp - s*(vq + vp*tau)
                        v(i,q) = vq + s*(vp - vq*tau)
                        !Add eigenvalue to d
                        d(i) = a(i,i)
                  enddo
                  
                  
                  sum = 0.d0
                  maxelement = 0.d0
                  !Refind greatest value and convergence
                  if(n == 1) then
                        do i=1,np
                              do j=i+1,np
                                    sum = sum + a(i,j)**2
                                    if(maxelement < abs(a(i,j))) then
                                          maxelement = abs(a(i,j))
                                          maxi = i
                                          maxj = j
                                    endif
                              enddo
                        enddo
                  else
                        do i=1,np-n
                              do j=n+1,np
                                    sum = sum + a(i,j)**2
                                    if(maxelement < abs(a(i,j))) then
                                          maxelement = abs(a(i,j))
                                          maxi = i
                                          maxj = j
                                    endif
                              enddo
                        enddo
                  endif
                  
                  p = maxi
                  q = maxj
                  nrot = nrot + 1
                  if(nrot>maxiter) then
                        write(*,*)"Jacobi: maximum iterations exceeded"
                        exit
                  endif
            enddo
      end

            
      subroutine sort_ev(d,v,n)
            implicit none
            integer :: n
            real*8 :: d(n),v(n,n)
            real*8 :: dcheck
            real*8 :: vcheck(n)
            integer :: i,j,k
            k = 1
            do i=1,n-1
                  k = i
                  dcheck = d(i)
                  vcheck = v(:,i)
                  do j=i+1,n
                        if(d(j) < dcheck) then
                              k = j
                              dcheck = d(j)
                              vcheck = v(:,j)
                        endif
                  enddo

                  if(k /= i) then
                        d(k) = d(i)
                        d(i) = dcheck
                        v(:,k) = v(:,i)
                        v(:,i) = vcheck

                        ! vcheck = v(:,i)
                        ! v(:,i) = v(:,k)
                        ! v(:,k) = vcheck
                  endif
            enddo
      end