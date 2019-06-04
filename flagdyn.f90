!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!% 2017, Teresa Ruiz, truiz@seas.harvard.edu
!% Dynamics of rosette formation
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   

module structures
       !flagellate class contains:
       ! sec: 0 for non-secreting cells, 1 for secreting cells
       ! div: 0 for non-dividing cells, j for dividing cells, being j the id of the particle dividing from
       ! divN:   ! steps into division, 0 for non-dividing cells, k for dividing cells
       ! divT:  ~cell cycle, time before division
       ! rSEC:   distance between secreted ECM and cell head, rSEC starts being 0.5*(sig(1)-sig_ecm)  and increases each timestep dr=(sig_ecm+2**(1.0/6)-1)/gSEC
       ! n:    normal axis of plane form by two cells. division plane if dividing, previous final division plane if not dividing
       ! X(3,3)   position of beads forming cell , dim (3,3)  3 beads x 3D
       !F(3,3)    forces on  beads forming cell , dim (3,3)  3 beads x 3D


    type flagellate
       integer:: sec=0,div=0,divN
       real(8):: divT=0,rSEC=0
       real(8):: n(3),X(3,3),F(3,3)
     
    end type flagellate

end module structures




module functions

contains

real(8) function gasdev(dseed)
!returns number from a Normal distribution with 0 mean an standard deviation = 1
!dseed is the seed

!     Returns a Normal(0,1) random variate (idum < 0 for first call).
!     This program was copied from Press, Flannery, Teukolsky, and
!     Vetterling, Numerical Recipies, Cambridge University Press, 1988.

      implicit none
     
      real(8)::v1,v2,r,fac,gset,dseed,ran1b
      integer::iset

      iset=0
      if (iset.eq.0) then
  1      call ggub(dseed,ran1b)
        v1 = 2.0d0*ran1b - 1.0d0
        call ggub(dseed,ran1b)
        v2 = 2.0d0*ran1b - 1.0d0
        r = v1**2 + v2**2
        if (r.ge.1) goto 1
        fac = dsqrt(-2.0d0*dlog(r)/r)
        gset = v1*fac
        gasdev= v2*fac
        iset = 1 
      else
        gasdev = gset
      
        iset = 0 
      end if
 
      end function


function cross(a, b)
!cross product function used to calculate plane of division
  real(8), DIMENSION(3) :: cross
  real(8), DIMENSION(3), INTENT(IN) :: a, b

  cross(1) = a(2) * b(3) - a(3) * b(2)
  cross(2) = a(3) * b(1) - a(1) * b(3)
  cross(3) = a(1) * b(2) - a(2) * b(1)

end function cross


real(8) function CDF(x,mu,sigma)  
!cumulative density function used to calculate probability of division and secretion
  real(8),intent(in)::x,mu,sigma

  CDF=1d0/(1+exp(-(x-mu)/sigma))

end function CDF


end module functions


module parameters
use structures
use functions
implicit none
save
real(8),parameter::e0=1,l0=1,m0=1,kT=2,t0=1   !untis of energy, length, thermal energy and time. This can be set to real values to run simulations with real units instead of system units.  
integer::ncell,n_ecm,nsteps    !number of cells, ECM particles and steps
real(8):: tau_div,sigma_div,kFENE,kLEN,kSEC,kDIV,dDIV,gSEC,gDIV,dSEC,dc,dc_ecm,fc,fc_ecm,dseed
real(8):: eps_adh,eps_ecm,sig_ecm,eta,dt,sig(3)
real(8),allocatable::X_ecm(:,:),F_ecm(:,:)  !position vector and force vector of ECM particles, arrays have dimensions (n_ecm,3)
real(8),allocatable:: CC(:)    !cell clock, it's reset after division, CC has dimensions (ncell)   
integer,allocatable::DIV(:,:)    !matrix will contain pairs of particles dividing, each line will contain 1 pair [i,j,steps into Division]
integer,allocatable::SEC(:,:)    !matrix will contain particle secretion, each line will contain 1 pair [icell,steps into secretion]
type(flagellate),allocatable::cell(:)   !flagellate structure initialization

contains

function normal_axis(icell,jcell)  
!returns normal axis of plane formed by two flagellates
  implicit none
  real(8), DIMENSION(3) :: normal_axis,body_axis_i,body_axis_j
  integer,intent(in):: icell,jcell
 
 
  body_axis_i = cell(icell)%X(1,:) - cell(icell)%X(3,:)
  body_axis_j = cell(jcell)%X(1,:) - cell(jcell)%X(3,:)
    
  normal_axis = cross(body_axis_i, body_axis_j)   !cross product

end function normal_axis


end module parameters



!---------------------MAIN PROGRAM-----------------------


program flagdyn
   use parameters
   use structures
   use functions
   implicit none

   integer::i,j,k,l,ksave,nsave
   real(8)::START, FINISH,pi,randR(3),sig1,sig2,sig3,force
   real(8)::epot_WCA,epot_FENE,epot_LEN,epot_SLJ,epot_ECM,epot_MIT,epot_SEC,epot


   call cpu_time(start)  !this is just to monitor performance 

   !read parameters from external file
   open(11,file='param.dat')
   read(11,*)eps_adh,eps_ecm,sig1,sig2,sig3,sig_ecm,tau_div,eta
   print*,eps_adh,eps_ecm,sig1,sig2,sig3,sig_ecm,tau_div,eta
   close(11)


   !------------initiate  variables

   !system starts with 1 cell and no ECM particles
   ncell = 1
   n_ecm = 0

   sig = (/sig1,sig2,sig3/)
  
   dseed=32513
   pi=4.0*atan(1.0)
          
   !constants controlling division time of cells  cell(i)%divT=sigma_div*R(0,1)+tau_div
   tau_div = tau_div*t0
   sigma_div = tau_div/10d0
    
   !I kept fixed the following numbers but you actually can play with them
   nsteps=5000000     !maximum number of steps
   dt=0.001*t0        !time step
   kFENE = 500*e0     !FENE constant to keep beads connected
   kLEN = 200*e0      !spring constant to keep cell elongated
   kDIV =50*e0        !spring constant to keep cells together during division 
   !dDIV=0.05          !
   gSEC=5000          ! maximum number of steps to complete ECM secretion
   gDIV=7000          ! maximum number of steps to complet cell division
   !dSEC=(1-(sig1-sig_ecm)/(sig1+sig_ecm))/gSEC     !  

    !friction coeff (fc) and constrant for browninan term (dc) for cell beads and ECM beads
   fc = 6*pi*eta*sig2              !simplification, same for all of them, not considering flagellum, mass differences...
   fc_ecm= 6*pi*eta*sig_ecm
   dc = sqrt(2*kT*dt/fc)
   dc_ecm = sqrt(2*kT*dt/fc_ecm)

   !variables for writing files
   ksave=10    
   nsave=0

   !initialize structures 
   call init()

   !dynamics loop    
   do k=1,nsteps
       
       CC=CC+dt   !increase cell clock for all cells at every time step

       call forces(epot_WCA,epot_FENE,epot_LEN,epot_SLJ,epot_ECM,epot_MIT,epot_SEC) !calculate forces
       
   
       do i=1,ncell         !update positions
          do l=1,3    !loop over the three particles comprising the flagellate

             randR= (/gasdev(dseed), gasdev(dseed), gasdev(dseed)/)
             cell(i)%X(l,:) = cell(i)%X(l,:) + cell(i)%F(l,:)*dt/fc  + dc*randR  
          end do
       end do
      

       do i=1,n_ecm  !loop over ECM particles

          randR= (/gasdev(dseed), gasdev(dseed), gasdev(dseed)/)
          X_ecm(i,:) = X_ecm(i,:) + F_ecm(i,:)*dt/fc_ecm  + dc_ecm*randR  

       end do

  
       call division(k)  !check for division

       call secretion(k)  !check for secretion
      
       
!-----------------------------------------------------------------------------------------------------
!save configurations
       if (k == ksave) then
          !print*,k,epot_WCA/ncell,epot_FENE/ncell,epot_LEN/ncell,epot_SLJ,epot_ECM,epot_MIT,epot_SEC
          write(600,'(i7,5f12.5)')k,epot_WCA/ncell,epot_FENE/ncell,epot_LEN/ncell,epot_SLJ/(max(n_ecm,1)),epot_ECM/(max(n_ecm,1))
          nsave=nsave+1
          write(55,*) nsave,ncell
          do i=1,ncell
              do j=1,3
                 write(55,'(3f12.6)')cell(i)%X(j,1),cell(i)%X(j,2),cell(i)%X(j,3)
              end do
          end do
          write(56,*)nsave,n_ecm
          do i=1,n_ecm
              write(56,'(3f12.6)')X_ecm(i,1),X_ecm(i,2),X_ecm(i,3)
          end do
          ksave=ksave+2000
       end if

  
       if (ncell ==30) exit   !ends dynamics when colony has 30 cells
!--------------------------------------------------------------------------------------


   end do

 
end program flagdyn

 subroutine init()
    use parameters
    use structures
    use functions
    implicit none
    real(8)::x2,x3

    !initializes 1st cell
   
    allocate(cell(1),CC(ncell),X_ecm(0,3), DIV(0,3),SEC(0,2),F_ecm(0,3))   

    x2=0.5*(sig(1)+sig(2))
    x3=x2+ 0.5*(sig(2)+sig(3))
 

    cell(1)%n= (/ 0d0, 0d0, 1d0 /)     !select xy as division plane for first particle
    cell(1)%div = 0        ! 0 for non-dividing cells
    cell(1)%divT=sigma_div*gasdev(dseed)+tau_div !normrnd(tau_div,sigma_div);
    cell(1)%divN=0      ! 0 for non-dividing cells
    cell(1)%sec = 0        !0 for non-secreting cells
   
    cell(1)%X = transpose(reshape ((/ 0d0, 0d0, 0d0, x2, 0d0, 0d0, x3, 0d0, 0d0 /), (/ 3, 3 /)))    !set initial positions for 1st cell 
    cell(1)%F=0

    CC=0  !set cell clock to 0
        
  end subroutine init


  subroutine forces(epot_WCA,epot_FENE,epot_LEN,epot_SLJ,epot_ECM,epot_MIT,epot_SEC)
    use parameters
    use structures
    use functions
    implicit none
    real(8),intent(out)::epot_WCA,epot_FENE,epot_LEN,epot_SLJ,epot_ECM,epot_MIT,epot_SEC
    integer::i,j,icell,jcell,l1,l2,l,k,iECM
    real(8)::r0,soft,softi,softj,sigma,eps
    real(8)::delta(3),cut(3),ff(3),rr(3)
  
    
    do i=1,ncell
       cell(i)%F=0
    end do

    if (n_ecm > 0) then
       F_ecm=0
    end if
   
    epot_WCA = 0
    eps=e0
  
! steric forces   between cells
    if (ncell>1) then
  
       do icell = 1, ncell-1

          do jcell =icell +1, ncell
             if (cell(icell)%div .ne. jcell) then

                if (cell(icell)%divN==0) then
                   softi=1
                else
                   softi=1d0*cell(icell)%divN/gDIV
                   !print*,'softi',softi
                end if
    
                if (cell(jcell)%divN==0) then
                   softj=1
                else
                   softj=1d0*cell(jcell)%divN/gDIV
      
               end if
      
                do l1 = 1,3
                   do l2 = 1,3

                      rr= cell(icell)%X(l1,:) - cell(jcell)%X(l2,:)
                      sigma = 0.5 * (sig(l1) + sig(l2))
      
                      call force_WCA(sigma,rr,eps,epot_WCA,ff,softi*softj)  !try softer potential

                      cell(icell)%F(l1,:) = cell(icell)%F(l1,:) + ff
                      cell(jcell)%F(l2,:) = cell(jcell)%F(l2,:) - ff

                   end do
                end do
             end if
          end do
       end do
    end if



!intramolecular + adhesive forces
    
    epot_FENE=0
    epot_LEN=0
    epot_SLJ=0
    eps=e0
    r0 = 2.0*(0.5*sig(1)+sig(2)+0.5*sig(3)) !two times real distance to ensure straightness of molecule
    delta=(/ (sig(1)+sig_ecm)/2.0-1 ,(1.5*sig(2)+sig_ecm)/2.0-1, (1.5*sig(3)+sig_ecm)/2.0-1 /)  ! shifted lenndar jones parameter for adhesion forces between cells and ECM 
    cut= (/2.5, 2**(1.0/6) , 2**(1.0/6)/)  ! 2.5 for the head to add attractive region in SLJ


    do icell= 1, ncell

       do i=1,3

          if (i.ne.3) then
             j=i+1
             rr = cell(icell)%X(i,:) - cell(icell)%X(j,:)

             call force_FENE(rr,sig(i),sig(j),eps,epot_FENE,ff,icell)  !spring force between connected springs , FENE contains WCA
      
             cell(icell)%F(i,:) = cell(icell)%F(i,:) + ff
             cell(icell)%F(j,:) = cell(icell)%F(j,:) - ff

          end if

          if (i==1) then
             j = 3

             rr = cell(icell)%X(i,:) - cell(icell)%X(j,:)
             
             call force_SPRING(rr,r0,kLEN,epot_LEN,ff)  ! force to keep elongation
             
             cell(icell)%F(i,:) = cell(icell)%F(i,:) + ff
             cell(icell)%F(j,:) = cell(icell)%F(j,:) - ff
             
          end if
          
          !adhesion to ECM, only for heads, repulsion for others

          if (cell(icell)%divN==0) then
             softi=1d0
          else
             softi=1d0*cell(icell)%divN/gDIV
          end if

          
          do iECM= 1,n_ecm
           
             if (SEC(iECM,1) .ne. 0) then  !iECM is being secreted
                soft = 1d0*SEC(iECM,2)/gSEC
             else
                soft =1d0
             end if

             rr = cell(icell)%X(i,:) - X_ecm(iECM,:)
             
             if ((SEC(iECM,1) .ne. icell) .or. (i.ne.1)) then      !if jecm~=jsec || l~=1
                    
                call  force_SLJ(rr,eps_adh,delta(i),cut(i),epot_SLJ,ff,soft*softi)   !interaction between ECM and icell, ecm is not being secreted by icell
               
             else
                call  force_SLJ(rr,eps_adh,cell(icell)%rSEC-2**(1.0/6),cut(2),epot_SLJ,ff,soft*softi)  ! interaction between ECM and icell, ecm is being secreted by icell, this force will be close to 0, just added so when particle is finally secreted there is a smooth transition between forces before y after 
               
             end if
                        
             cell(icell)%F(i,:) = cell(icell)%F(i,:) + ff
             F_ecm(iECM,:) = F_ecm(iECM,:) - ff
    
          end do
       end do
       
    end do

    !ECM interaction
   
    epot_ECM=0
 
    if (n_ecm >1 ) then

       do i=1,n_ecm -1
    
          if (SEC(i,1) .ne. 0 ) then           !particle being secreted
             softi=1d0*SEC(i,2)/gSEC
          else
             softi=1
          end if

          do j = i+1, n_ecm

             if (SEC(j,1) .ne. 0 ) then          !particle being  secreted
                softj=1d0*SEC(j,2)/gSEC
             else
                softj=1
             end if

             rr= X_ecm(i,:) - X_ecm(j,:)

             call force_LJ42(rr,eps_ecm,sig_ecm,epot_ECM,ff,softi*softj)    !!!SAME eps_adh for ecm-ecm adhesion and ecm-cell adhesion ????

             F_ecm(i,:) = F_ecm(i,:) + ff
             F_ecm(j,:) = F_ecm(j,:) - ff
          end do
       end do
    end if

    !mitosis
    epot_MIT=0
    eps=e0
    
    if (size(DIV,1) > 0 )  then

       do k=1, size(DIV,1)
          icell = DIV(k,1)
          jcell = DIV(k,2)
          do l1=1,3
             do l2=1,3
                if (l1==l2) then
                   sigma = 0.5*sig(l1)+ 0.5*sig(l1)*DIV(k,3)/gDIV
               
                else
                   sigma = 0.5*(sig(l1)+sig(l2))
                end if

                rr = cell(icell)%X(l1,:) - cell(jcell)%X(l2,:)
                 
                  
                call force_WCA(sigma,rr,eps,epot_MIT,ff,1d0)                
                    
                cell(icell)%F(l1,:) = cell(icell)%F(l1,:) + ff
                cell(jcell)%F(l2,:) = cell(jcell)%F(l2,:) - ff


             end do
          end do
       end do
    end if

    
    !secretion 
    
    epot_SEC=0


       do j=1, n_ecm
          
          if (SEC(j,1)>0) then
          
             icell = SEC(j,1)
             r0 = cell(icell)%rSEC 
             rr = cell(icell)%X(1,:) - X_ecm(j,:)
   
             call force_SPRING(rr,r0,kSEC,epot_SEC,ff)
          
             cell(icell)%F(1,:) = cell(icell)%F(1,:) + ff
             F_ecm(j,:) = F_ecm(j,:) - ff
          end if

       end do
      
   
 end subroutine forces


subroutine force_WCA(sigma,rr,eps,epot,ff,soft)
  use parameters
  implicit none
  real(8),intent(in):: sigma,rr(3),eps,soft
  real(8),intent(inout)::epot
  real(8),intent(out)::ff(3)
  real(8)::cut,r2,r22,r1

  ff=0
  cut=sigma*2**(1.0/6)

  r2= dot_product(rr,rr)

  if (r2.lt.cut*cut) then

     r1 = (sigma**2)/r2
     r22 = r1**3

     epot= epot+ soft*4*eps*r22*(r22-1d0)

     ff = soft*rr*48*eps*r22*(r22-0.5d0)/r2
      
  end if
 
end subroutine FORCE_WCA



subroutine force_FENE(rr,s1,s2,eps,epot,ff,icell)
  use parameters
  implicit none
  real(8),intent(in):: rr(3),s1,s2,eps
  real(8),intent(inout)::epot
  real(8),intent(out)::ff(3)
  integer,intent(in)::icell
  real(8)::r1,rFENE,sigma,epot2,r0
  integer::i,j
  r1=norm2(rr)
  rFENE = 0.7*(s1+s2)  !distance between centers + 20% of this distance
  sigma = 0.5 * (s1 + s2)
  
  call force_WCA(sigma,rr,eps,epot,ff,1d0)
  
  epot = epot -(0.5*kFENE*rFENE**2)*log(1-(r1/rFENE)**2)
  ff= ff-kFENE*rr/(1-(r1/rFENE)**2) 
  if (r1>rFENE) then
     print*,'problem fene r1>rFENE',r1,rFENE,icell
         
         
         write(55,*) 'last',ncell
         do i=1,ncell
            do j=1,3
               write(55,'(3f12.6)')cell(i)%X(j,1),cell(i)%X(j,2),cell(i)%X(j,3)
            end do
         end do
         write(56,*)'last',n_ecm
         do i=1,n_ecm
            write(56,'(3f12.6)')X_ecm(i,1),X_ecm(i,2),X_ecm(i,3)
         end do
     stop

  end if
  !r0=s1+s2
  !epot = epot + 0.5*kFENE*(r1-r0)**2
  !ff = ff -kFENE*(r1-r0)*rr/r1
  
end subroutine force_FENE





subroutine force_SPRING(rr,r0,k,epot,ff)
  use parameters
  implicit none
  real(8),intent(in):: rr(3),r0,k
  real(8),intent(inout)::epot
  real(8),intent(out)::ff(3)
  real(8)::r1

  r1=norm2(rr)
 
  epot = epot + 0.5*k*(r1-r0)**2
  ff = -k*(r1-r0)*rr/r1
    
end subroutine force_SPRING


subroutine force_SLJ (rr,eps,delta,cut,epot,ff,soft)                                   
 use parameters
 implicit none
 real(8),intent(in):: rr(3),eps,delta,cut,soft
 real(8),intent(inout)::epot
 real(8),intent(out)::ff(3)
 real(8)::r1,cutslj

 r1=norm2(rr)
   
 cutslj=cut+delta
 ff=0
 if (r1 <= cutslj) then
    
    epot= epot+ soft*4*eps*((1.0/(r1-delta))**(12)-(1.0/(r1-delta))**(6))
    ff  = soft*24*eps*(2*(1.0/(r1-delta))**(13)-(1.0/(r1-delta))**(7))*rr/r1 
 end if
end subroutine force_SLJ



subroutine force_LJ42(rr,eps,sigma,epot,ff,soft)
use parameters
implicit none
real(8),intent(in):: rr(3),eps,sigma,soft
real(8),intent(inout)::epot
real(8),intent(out)::ff(3)
real(8):: r2,r1,cut
    
  r2= dot_product(rr,rr)
  cut=4
  ff=0

  if (r2.lt.cut*cut) then

     r1 = (sigma**2)/r2

     epot=epot+soft*4*eps*r1*(r1-1d0)

     ff = soft*rr*16*eps*r1*(r1-0.5d0)/r2
     
  end if
end subroutine force_LJ42


subroutine division(step)
  use parameters
  use structures
  use functions
  implicit none
  integer,intent(in)::step
  integer::i,j,nrem,icell,jcell
  integer,allocatable::DIVnew(:,:)
  real(8)::normal(3)
  real(8)::rand,acceptance_rate

  if (size(DIV,1)>0) then
     nrem=0
     i=1
     do
        if (size(DIV,1)<i) exit

        DIV(i,3)= DIV(i,3)+1
        icell=DIV(i,1)
        jcell=DIV(i,2)

        if (DIV(i,3) >= gDIV) then   !mitosis is complete between icell and jcell
               
           normal=normal_axis(icell,jcell)
           cell(icell)%n = normal
           cell(jcell)%n = normal

           cell(icell)%div=0
           cell(jcell)%div=0
           CC(icell) = 0
           CC(jcell) = 0
           DIV(i,:)=0
           nrem = nrem+1     
           
           cell(icell)%divT=sigma_div*gasdev(dseed)+tau_div    ! next time for division
           cell(jcell)%divT=sigma_div*gasdev(dseed)+tau_div
                
                   
        end if

        cell(icell)%divN = DIV(i,3)
        cell(jcell)%divN = DIV(i,3)

        i=i+1
     end do


     allocate(DIVnew(size(DIV,1)-nrem,3))

     j=0
     do i=1, size(DIV,1)
        if (DIV(i,1).ne.0) then
           j=j+1
           DIVnew(j,:)=DIV(i,:)
        end if
     end do

     if (j.ne.size(DIVnew,1)) then
        print*,'there is a problem with sizes'
     end if

     deallocate(DIV)
     allocate(DIV(size(DIVnew,1),3))

     DIV=DIVnew
  end if

  ! is there a new cell dividing?

  if ( mod(step,10)==0 ) then ! %mod(step,100) == 0
        
     do icell=1,ncell

        if (cell(icell)%div==0 .and. cell(icell)%sec==0) then
           acceptance_rate = 1-CDF(CC(icell),tau_div,sigma_div)
        
           call ggub(dseed,rand)
           if (acceptance_rate <= rand) then
            
              print*,icell,'starts division'

              call divide(icell)
           end if
        end if
     end do
      
  end if
    

 
    
end subroutine division
     



subroutine divide(icell)
  use parameters
  use structures
  use functions
  implicit none
  integer,intent(in)::icell
  integer::i,jcell,newcell,DIVnewsize,j
  real(8)::body_axis(3),division_axis(3),div_matrix(3,3)
  type(flagellate),allocatable::cellnew(:)
  integer,allocatable:: DIVnew(:,:)
  real(8),allocatable::CCnew(:)
  real(8)::rr(3),r1

  
  newcell = ncell +1
  DIVnewsize= size(DIV,1)+1
  allocate(cellnew(newcell),DIVnew(DIVnewsize,3),CCnew(newcell))

  cellnew(1:ncell)=cell
  cellnew(newcell)%X=cell(icell)%X
  
 !division axis perpendicular to body axis 
  body_axis =  cell(icell)%X(1,:)-cell(icell)%X(3,:)

  
  division_axis= cross(cell(icell)%n, body_axis)
  division_axis = division_axis / norm2(division_axis)

  do i=1,3
     div_matrix(i,:)=0.25*sig(i)*division_axis   
  
 end do

  cellnew(icell)%X = cellnew(icell)%X + div_matrix
  cellnew(newcell)%X = cellnew(newcell)%X - div_matrix

 
   
  DIVnew(1:DIVnewsize-1,:)=DIV
  DIVnew(DIVnewsize,:)=(/icell,newcell,1/)   !connection between new dividing cells
  
  CCnew(1:ncell)=CC
  CCnew(newcell)=CC(icell)
  deallocate(DIV,cell,CC)
  allocate(cell(newcell),DIV(DIVnewsize,3),CC(newcell))

  cell=cellnew
  DIV=DIVnew
  CC=CCnew

  cell(icell)%div = newcell
  cell(icell)%divT= 0
  cell(icell)%divN= 1


  cell(newcell)%n = cell(icell)%n
  cell(newcell)%div =icell
  cell(newcell)%sec=0
  cell(newcell)%divT=0
  cell(newcell)%divN=1


  ncell = ncell +1


end subroutine divide


subroutine secretion(step)
  use parameters
  use structures
  use functions

  implicit none
  integer,intent(in)::step
  integer::i,j,icell,jcell,nrem
  real(8),allocatable::SECnew(:,:)
  real(8)::rand,acceptance_rate
 

  do i=1,n_ecm

     if (SEC(i,1).ne. 0) then
        icell=SEC(i,1)

 !       SEC(i,3)= SEC(i,3)+(sig_ecm+2**(1.0/6)-1)/gSEC
        if (SEC(i,2)<=gSEC) then
           SEC(i,2)=SEC(i,2)+1
           cell(icell)%rSEC = cell(icell)%rSEC + (sig_ecm+2**(1.0/6)-1)/gSEC
        else               
           cell(icell)%sec = 0   !secretion is complete
           cell(icell)%rSEC=0
           SEC(i,:)= 0
        end if
     end if
  end do


  !secretion?
  if (mod(step,10) == 0) then
     call ggub(dseed,rand)
     icell = ceiling(ncell*rand)
            
     if (icell==0) then
        icell = icell +1
     end if
     
     if (cell(icell)%sec == 0  .and. cell(icell)%div ==0) then 
        acceptance_rate = CDF(CC(icell),tau_div,sigma_div) 
        
        call ggub(dseed,rand)
        if (acceptance_rate <= rand) then
           call secrete(icell)
        end if
     end if
        
  end if
end subroutine secretion



subroutine secrete(icell)
  use parameters
  use structures
  use functions
  implicit none
  integer,intent(in)::icell
  integer::newecm
  real(8),allocatable::X_ecmNew(:,:),SECnew(:,:)
  real(8)::rand,body_axis(3)


  newecm = n_ecm +1

print*,'secrete particle',newecm   
  allocate(X_ecmNew(newecm,3),SECnew(newecm,2))
  
  X_ecmNew(1:n_ecm,:)=X_ecm
  X_ecmNew(newecm,:) = cell(icell)%X(1,:) 

  body_axis = cell(icell)%X(1,:) - cell(icell)%X(3,:)
  body_axis = body_axis / norm2(body_axis)

  X_ecmNew(newecm,:) = X_ecmNew(newecm,:) + 0.5*(sig(1)-sig_ecm)*body_axis

  SECnew(1:n_ecm,:)=SEC
  SECnew(newecm,:)=(/icell,0/)  !connections between secretion particle and cell, and number of steps that has being secreted

  deallocate(SEC,X_ecm,F_ecm)

  allocate(SEC(newecm,3),X_ecm(newecm,3),F_ecm(newecm,3))

  SEC=SECnew
  X_ecm=X_ecmNew
  
  cell(icell)%sec = 1
  cell(icell)%rSEC = 0.5*(sig(1)-sig_ecm)
  n_ecm = n_ecm +1

end subroutine secrete


!
!  Generacion numeros aleatorios
!
      subroutine ggub(dseed,r)

      real(8):: z,d2p31m,d2pn31,dseed,r
      d2p31m=2147483!647
      d2pn31=2147483!648

      z = dseed
      z = mod(16807*z,d2p31m)
      r = z / d2pn31
      dseed = z

      end subroutine
      
