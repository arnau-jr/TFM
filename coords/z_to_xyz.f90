      use tools
      implicit none

      integer,parameter :: maxatom = 100
      character :: S(maxatom)*2
      real :: xyz(3,maxatom)

      integer :: Batom,Aatom,Datom
      real :: B,A,D
      real :: x,y,z
      real :: localv(3),basis(3,3)
      integer :: i,n
      
      character :: input_filename*90,output_filename*90

      call get_command_argument(1,input_filename)

      if(index(input_filename,".gzmat")==0) then
            print*, "Error : input file does not have a Gaussian Z-matrix extension"
            stop
      endif
      if(command_argument_count() == 2) then
            call get_command_argument(2,output_filename)
      else
            write(output_filename,"(A,A)")input_filename(:index(input_filename,".gzmat")),"xyz"  
      endif


      open(1,file=input_filename)
      open(2,file=output_filename)

      read(1,*)
      read(1,*)
      read(1,*)
      read(1,*)
      read(1,*)

      read(1,*)S(1)
      x = 0.
      y = 0.
      z = 0.
      xyz(:,1) = (/x,y,z/)

      read(1,*)S(2),Batom,B
      x = 0.
      y = 0.
      z = B
      xyz(:,2) = (/x,y,z/)

      read(1,*)S(3),Batom,B,Aatom,A
      x = xyz(1,Batom) + 0.
      y = xyz(2,Batom) + B*sin((pi/180)*A)
      z = xyz(3,Batom) - B*cos((pi/180)*A)
      xyz(:,3) = (/x,y,z/)

      n = 3
      do while(.true.)
            read(1,*,END=10)S(n+1),Batom,B,Aatom,A,Datom,D
            ! print*,n+1,S(n+1),Batom,Aatom,Datom
            call get_local_coords(xyz(:,Batom),xyz(:,Aatom),xyz(:,Datom),basis)

            x = B*sin((pi/180)*A)*sin((pi/180)*D)
            y = B*sin((pi/180)*A)*cos((pi/180)*D)
            z = B*cos((pi/180)*A)
            localv = (/x,y,z/)
            xyz(:,n+1) = xyz(:,Batom) + matmul(basis,localv)
            n = n + 1
      enddo
      10 continue

      write(2,"(I2)")n
      write(2,*)""
      
      do i=1,n
            write(2,"(A,7X,F8.5,7X,F8.5,7X,F8.5)")S(i),xyz(:,i)
      enddo
      end
