      use tools
      use ff_utils
      implicit none


      end




      subroutine verletvel_step(D,N,dt,r,v,L,rc,force,f,U,P)
            implicit none
            integer :: D,N
            real*8 :: dt,r(N,D),v(N,D),L,rc,f(N,D),U,P
            real*8 :: fold(N,D)

            fold = f
            r = r + v*dt + 0.5d0*f*dt**2
            G = build_gradient(Natoms,xyz,Nbonds,Nangles,Ntorsions,&
            bond_pairs,angle_pairs,torsion_pairs,req,kb,aeq,ka,An,n,delta)
            v = v + 0.5d0*(fold+f)*dt
      end