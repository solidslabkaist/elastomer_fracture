************************************************************************
! 
! Contributing authors: Jaehee Lee (KAIST), Jeongun Lee (KAIST), 
!                       Seounghee Yun (KAIST), Sanha Kim (KAIST),    
!                       Shawn A. Chester (NJIT), Hansohl Cho (KAIST)
! 
! Size-dependent fracture in elastomers: experiments and continuum modeling, 2024
! 
! For further information, contact: hansohl at kaist dot ac dot kr
! 
! User element for damage evolution, and large elastic deformation
! in elastomeric materials in 2D plane-stress settings.
! 
! Solution variables (or nodal variables) are the displacements and the
! damage (phase-) field.
! 
! This subroutine is for the following element types
!  > two-dimensional 4 node isoparametric (plane-stress) element as shown below
!       with 4pt (full) gauss integration.
! 
!     
!              A eta (=xi_2)
!  4-node      |
!   quad       |Face 3
!        4-----------3
!        |     |     |
!        |     |     |
!  Face 4|     ------|---> xi (=xi_1)
!        |           | Face2
!        |           |
!        1-----------2
!          Face 1
! 
! 
***********************************************************************
! 
! User element statement in the input file:
! 
!  2D Plane-stress elements
!  *User Element,Nodes=4,Type=U1,Iproperties=2,Properties=7,Coordinates=2,Variables=16,Unsymm
!  1,2,11
! 
!     State Variables
!     --------------------------------------------------------------
!     Global SDV's (used for visualization)
!       1) Bond-stretch (blambda)
!       2) History function (histmax)
!       3) True stress (2,2) component (S22)
!       4) Deformation gradient (3,3) component (F33)
! 
!     Local SDV's (used for the solution procedure)
!       j = 0
!       do k = 1,nIntPt
!          svars(1+j) = blambda_tau
!          svars(2+j) = Histmax_tau         
!          svars(3+j) = S22_tau
!          svars(4+j) = F_tau(3,3)
!          j = j + nlSdv
!       end loop over k
!
!     In the input file, set 'User output variables'= number of global SDV's
!
!     In the input file, set 'ngSdv'= number of global SDV's
!
!     In the input file, set 'nlSdv'= number of local SDV's
!
!     In the input file, set 'varibles'=(nlSdv*nIntPt)
!
!
!     Material Properties Vector
!     --------------------------------------------------------------
!     Gshear0 = props(1)       ! Shear modulus
!     ENK     = props(2)       ! number of Kuhn segments in a chain
!     Eb_bar  = props(3)       ! Bond stiffness
!     Kbulk   = props(4)       ! Bulk modulus
!     erff    = props(5)       ! chain scission energy per unit
!     lc      = props(6)       ! characteristic length scale of the graident theory
!     zeta    = props(7)       ! kinetic modulus for the evolution of the damage
!     nlSdv  = jprops(1)       ! Number of local sdv's per integ pt
!     ngSdv  = jprops(2)       ! Number of global sdv's per integ pt
!
!***********************************************************************

      module global

      ! This module is used to transfer SDV's from the UEL
      !  to the UVARM so that SDV's can be visualized on a
      !  dummy mesh
      !
      !  globalSdv(X,Y,Z)
      !   X - element pointer
      !   Y - integration point pointer
      !   Z - SDV pointer
      !
      !  numElem
      !   Total number of elements in the real mesh, the dummy
      !   mesh needs to have the same number of elements, and 
      !   the dummy mesh needs to have the same number of integ
      !   points.  You must set that parameter value here.
      !
      !  ElemOffset
      !   Offset between element numbers on the real mesh and
      !    dummy mesh.  That is set in the input file, and 
      !    that value must be set here the same.

      integer numElem,ElemOffset,err

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Set the number of UEL elements used here
       parameter(numElem=10982)    ! Single-edge-notched
      ! parameter(numElem=16615)   ! Randomly perforated
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Set the offset here for UVARM plotting, must match input file!
       parameter(ElemOffset=100000)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      real*8, allocatable :: globalSdv(:,:,:)

      end module global

***********************************************************************

      SUBROUTINE UVARM(UVAR,DIRECT,T,TIME,DTIME,CMNAME,ORNAME,
     1 NUVARM,NOEL,NPT,LAYER,KSPT,KSTEP,KINC,NDI,NSHR,COORD,
     2 JMAC,JMATYP,MATLAYO,LACCFLA)

      ! This subroutine is used to transfer SDV's from the UEL
      !  onto the dummy mesh for viewing.  Note that an offset of
      !  ElemOffset is used between the real mesh and the dummy mesh.
      !  If your model has more than ElemOffset UEL elements, then
      !  this will need to be modified.
     
      use global
     
      include 'ABA_PARAM.INC'

      CHARACTER*80 CMNAME,ORNAME
      CHARACTER*3 FLGRAY(15)
      DIMENSION UVAR(NUVARM),DIRECT(3,3),T(3,3),TIME(2)
      DIMENSION ARRAY(15),JARRAY(15),JMAC(*),JMATYP(*),COORD(*)

C     The dimensions of the variables FLGRAY, ARRAY and JARRAY
C     must be set equal to or greater than 15.

      uvar(1) = globalSdv(noel-ElemOffset,npt,1)
      uvar(2) = globalSdv(noel-ElemOffset,npt,2)
      uvar(3) = globalSdv(noel-ElemOffset,npt,3)
      uvar(4) = globalSdv(noel-ElemOffset,npt,4)

      return
      end subroutine uvarm

****************************************************************************

      SUBROUTINE UEL(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     +     PROPS,NPROPS,coords,MCRD,NNODE,Uall,DUall,Vel,Accn,JTYPE,
     +     TIME,DTIME,KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,
     +     PREDEF,NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,
     +     NJPROP,PERIOD)

      use global
*
      IMPLICIT NONE
*
*     VARIABLES DEFINED IN UEL, PASSED BACK TO ABAQUS
*
      REAL(8) :: RHS,AMATRX,SVARS,ENERGY
*
*     VARIABLES PASSED INTO UEL 
*
      REAL(8) :: PROPS,coords,Uall,DUall,Vel,Accn,TIME,
     1  DTIME,PARAMS,ADLMAG,PREDEF,DDLMAG,PNEWDT,PERIOD
      INTEGER :: NDOFEL,NRHS,NSVARS,NPROPS,MCRD,NNODE,JTYPE,KSTEP,KINC,
     1  JELEM,NDLOAD,JDLTYP,NPREDF,LFLAGS,MLVARX,MDLOAD,JPROPS,NJPROP
*
      DIMENSION RHS(MLVARX,*),AMATRX(NDOFEL,NDOFEL),PROPS(*),
     1 SVARS(*),ENERGY(8),coords(MCRD,NNODE),Uall(NDOFEL),
     2 DUall(MLVARX,*),Vel(NDOFEL),Accn(NDOFEL),TIME(2),PARAMS(*),
     3 JDLTYP(MDLOAD,*),ADLMAG(MDLOAD,*),DDLMAG(MDLOAD,*),
     4 PREDEF(2,NPREDF,NNODE),LFLAGS(*),JPROPS(*)

      integer lenJobName,lenOutDir,nDim,nInt,nIntS
      character*256 jobName,outDir,fileName



      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      parameter(nInt=4)  ! number of volume integration pionts
      parameter(nIntS=1) ! number of surface integration points
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




      !----------------------------------------------------------------
      ! 
      ! Perform initial checks
      !
      !
      ! Open the debug/error message file
      !
      call getJobName(jobName,lenJobName)
      call getOutDir(outDir,lenOutDir)
      fileName = outDir(1:lenOutDir)//'\aaMSGS_'//
     +     jobName(1:lenJobName)//'.dat'
      open(unit=80,file=fileName,status='unknown')



      ! Make sure Abaqus knows you are doing a large
      !  deformation problem, I think this only matters
      !  when it comes to output in viewer
      !
      if(lflags(2).eq.0) then
         !
         ! lflags(2)=0 -> small disp.
         ! lflags(2)=1 -> large disp.
         !
         write(*,*) 'Abaqus thinks you are doing'
         write(*,*) 'a small displacement analysis'
         write(*,*) 'go in and set nlgeom=yes'
         write(80,*) 'Abaqus thinks you are doing'
         write(80,*) 'a small displacement analysis'
         write(80,*) 'go in and set nlgeom=yes'
         call xit
      endif


      ! Check to see if you are doing a general
      !  step or a linear perturbation step
      !
      if(lflags(4).eq.1) then
         !
         ! lflags(4)=0 -> general step
         ! lflags(4)=1 -> linear perturbation step
         !
         write(*,*) 'Abaqus thinks you are doing'
         write(*,*) 'a linear perturbation step'
         write(80,*) 'Abaqus thinks you are doing'
         write(80,*) 'a linear perturbation step'
         call xit         
      endif


      ! Do nothing if a ``dummy'' step
      !
      if(dtime.eq.0.0) return
      !
      ! Done with initial checks
      !
      !----------------------------------------------------------------



      !----------------------------------------------------------------
      ! 
      ! Call the paricular element to perform the analysis
      !
      !
      ! This is a plane-stress analysis
      !
      nDim = 2
      call UPS4(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     +     PROPS,NPROPS,coords,MCRD,NNODE,Uall,DUall,Vel,Accn,JTYPE,
     +     TIME,DTIME,KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,
     +     PREDEF,NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,
     +     NJPROP,PERIOD,
     +     nDim,nInt,nIntS)

      !
      ! Done with this element, RHS and AMATRX already returned
      !  as output from the specific element routine called
      !
      !----------------------------------------------------------------


      return
      end subroutine uel



           
************************************************************************
      subroutine UPS4(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     +     PROPS,NPROPS,coords,MCRD,NNODE,Uall,DUall,Vel,Accn,JTYPE,
     +     TIME,DTIME,KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,
     +     PREDEF,NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,
     +     NJPROP,PERIOD,
     +     nDim,nInt,nIntS)

      use global
*
      IMPLICIT NONE
*
*     VARIABLES DEFINED IN UEL, PASSED BACK TO ABAQUS
*
      REAL(8) :: RHS,AMATRX,SVARS,ENERGY
*
*     VARIABLES PASSED INTO UEL 
*
      REAL(8) :: PROPS,coords,Uall,DUall,Vel,Accn,TIME,
     1  DTIME,PARAMS,ADLMAG,PREDEF,DDLMAG,PNEWDT,PERIOD
      INTEGER :: NDOFEL,NRHS,NSVARS,NPROPS,MCRD,NNODE,JTYPE,KSTEP,KINC,
     1  JELEM,NDLOAD,JDLTYP,NPREDF,LFLAGS,MLVARX,MDLOAD,JPROPS,NJPROP
*
      DIMENSION RHS(MLVARX,*),AMATRX(NDOFEL,NDOFEL),PROPS(*),
     1 SVARS(*),ENERGY(8),coords(MCRD,NNODE),Uall(NDOFEL),
     2 DUall(MLVARX,*),Vel(NDOFEL),Accn(NDOFEL),TIME(2),PARAMS(*),
     3 JDLTYP(MDLOAD,*),ADLMAG(MDLOAD,*),DDLMAG(MDLOAD,*),
     4 PREDEF(2,NPREDF,NNODE),LFLAGS(*),JPROPS(*)

      real*8 u(nNode,2),du(nNode,ndofel),thetaNew(nNode)
      real*8 thetaOld(nNode),dtheta(nNode),dmNew(nNode)
      real*8 dmOld(nNode),dDM(nNode),uNew(nNode,ndofel)
      real*8 uOld(nNode,ndofel),u_t(nNode,ndofel),v(nNode,2)
      real*8 coordsC(mcrd,nNode)

      integer i,j,k,l,m,n,nIntPt,nDim,intpt,pOrder,a1,b1,a11,b11,face
      integer nInt,ii,jj,a12,b12,pe,nSdv,stat,q,nIntV,nIntPtV,p,ngSdv
      integer nlSdv,kk,lenJobName,lenOutDir,faceFlag,nIntS

      real*8 Iden(3,3),Le,theta0,blambda0,Ru(2*nNode,1),Rd(nNode,1)
      real*8 Kuu(2*nNode,2*nNode),Kdd(nNode,nNode),sh0(nNode),detMapJ0
      real*8 dshxi(nNode,2),dsh0(nNode,2),dshC0(nNode,2),detMapJ0C,Vmol
      real*8 Fc_tau(3,3),Fc_t(3,3),detFc_tau,detFc_t,detFc,w(nInt),ds
      real*8 sh(nNode),detMapJ,blambda_t,dsh(nNode,2),detMapJC,blamLmt
      real*8 dshC(nNode,2),dm_tau,dm_t,dDMdX(2,1),dDMdt,F_tau(3,3)
      real*8 F_t(3,3),detF_tau,xi(nInt,2),detF,TR_tau(3,3),T_tau(3,3)
      real*8 SpTanMod(3,3,3,3),blambda_tau,body(3),umeror
      real*8 Smat(4,1),Bmat(4,2*nNode),BodyForceRes(2*nNode,1),Qmat(4,4)
      real*8 Gmat(4,2*nNode),G0mat(4,2*nNode),Amat(4,4),wS(nIntS)
      real*8 xLocal(nIntS),yLocal(nIntS),Kdu(nNode,2*nNode),detF_t
      real*8 Kud(2*nNode,nNode),Nvec(1,nNode),ResFac,AmatUD(4,1),TanFac
      real*8 SpUDMod(3,3),SpDUMod(3,3),AmatDU(1,4)
      real*8 Histmax_t,Histmax_tau,dTRdF(3,3,3,3)
      real*8 EN,ENK,Eb,Eb_bar,lc,zeta,Gshear0,erff,Kbulk
      real*8 C_tau(3,3),C2_tau(2,2),S220,S22_t,S22_tau
      real*8 Amat13(1,4),Amat31(1,4),delF(3,3),Finv(3,3),F_tau2(3,3),F33
      real*8 avglambda_t,avglambda,root,Hist1,Hist2,Hist3 
      real*8 delta,per(3,3),Fper(3,3),blambda_per,Histmax_per,Tper(3,3)
      real*8 TRper(3,3),crit,sol
     
            
      real*8 zero,one,two,half,Pi,three,third,T33,tol,del

      integer maxit,nn
   
      parameter(maxit=1000)            
      parameter(zero=0.d0,one=1.d0,two=2.d0,half=0.5d0,Pi=3.141592653d0,
     +     three=3.d0,third=1.d0/3.d0)


      ! Get element parameters
      !
      nlSdv = jprops(1) ! number of local sdv's per integ point
      ngSdv = jprops(2) ! number of global sdv's per integ point


      ! Allocate memory for the globalSdv's used for viewing
      !  results on the dummy mesh
      !
      if(.not.allocated(globalSdv)) then

         !
         ! allocate memory for the globalSdv's
         !
         ! numElem needs to be set in the MODULE
         ! nInt needs to be set in the UEL
         !
         stat=0
c         allocate(globalSdv(numElem,nInt,ngSdv))
c         deallocate(globalSdv)
         allocate(globalSdv(numElem,nInt,ngSdv),stat=err)
         if(stat.ne.0) then
            write(*,*) '//////////////////////////////////////////////'
            write(*,*) 'error when allocating globalSdv'
            write(*,*) '//////////////////////////////////////////////'
            write(*,*) '   stat=',stat
            write(*,*) '  ngSdv=',ngSdv
            write(*,*) '   nInt=',nInt
            write(*,*) 'numElem=',numElem
            write(*,*) '  nNode=',nNode
            write(*,*) 'lbound(globalSdv)',lbound(globalSdv)
            write(*,*) 'ubound(globalSdv)',ubound(globalSdv)
            write(*,*) '//////////////////////////////////////////////'
            write(80,*) '//////////////////////////////////////////////'
            write(80,*) 'error when allocating globalSdv'
            write(80,*) '//////////////////////////////////////////////'
            write(80,*) '   stat=',stat
            write(80,*) '  ngSdv=',ngSdv
            write(80,*) '   nInt=',nInt
            write(80,*) 'numElem=',numElem
            write(80,*) '  nNode=',nNode
            write(80,*) 'lbound(globalSdv)=',lbound(globalSdv)
            write(80,*) 'ubound(globalSdv)=',ubound(globalSdv)
            write(80,*) '//////////////////////////////////////////////'
            call xit
         endif
         write(*,*) '-------------------------------------------------'
         write(*,*) '----------- globalSDV ALLOCATED -----------------'
         write(*,*) '-------------------------------------------------'
         write(*,*) '---------- YOU PUT NUMBER OF ELEMENTS -----------'
         write(*,*) '---------- numElem=',numElem
         write(*,*) '---------- UPS4 ELEMENTS ------------------------'
         write(*,*) '-------------------------------------------------'
         write(*,*) '---------- YOU PUT NUMBER OF POINTS -------------'
         write(*,*) '---------- nInt =',nInt
         write(*,*) '---------- nIntS=',nIntS
         write(*,*) '-------------------------------------------------'
         write(*,*) '---------- YOU PUT NUMBER OF SDVs ---------------'
         write(*,*) '---------- ngSdv=',ngSdv
         write(*,*) '---------- nlSdv=',nlSdv
         write(*,*) '-------------------------------------------------'
      endif


      ! Identity tensor
      !
      call onem(Iden)


      ! Obtain initial conditions
      !
      theta0   = 298d0
      blambda0 = one
      
      erff    = props(5)       ! chain scission energy per unit
      lc      = props(6)       ! characteristic length scale of the graident theory
      zeta    = props(7)       ! kinetic modulus for the evolution of the damage
      

      ! Initialize the residual and tangent matrices to zero.
      !
      Ru  = zero
      Rd = zero
      Kuu = zero
      Kdd = zero
      Kud = zero
      Kdu = zero
      Energy = zero
      
c      pnewdt = one
      



      ! Body forces
      ! 
      body(1:3) = zero


      ! Obtain nodal displacements and damage
      !
      k = 0
      do i=1,nNode
         do j=1,nDim
            k = k + 1
            u(i,j) = Uall(k)
            du(i,j) = DUall(k,1)
            uOld(i,j) = u(i,j) - du(i,j)
         enddo
         k = k + 1
         dmNew(i) = Uall(k)
         dDM(i) = DUall(k,1)
         dmOld(i) = dmNew(i) - dDM(i)       
      enddo
      




      ! Obtain current nodal coordinates
      !
      do i=1,nNode
         do j=1,nDim
            coordsC(j,i) = coords(j,i) + u(i,j)
         enddo
      enddo


      ! Impose any time-stepping changes on the increments of
      !  damage or displacement if you wish
      !
      ! damage increment
      !
      do i=1,nNode
         if(dabs(dDM(i)).gt.1d-1) then
            pnewdt = 0.5
            write(*,*) '---------- damage cut-back ----------'
            write(*,*) 'time = ', time
            return
         endif
      enddo
      !
      ! displacement increment, based on element diagonal
      !
      Le = dsqrt(((coordsC(1,1)-coordsC(1,3))**two) + 
     +     ((coordsC(2,1)-coordsC(2,3))**two))
      !
      do i=1,nNode
         do j=1,nDim
            if(dabs(du(i,j)).gt.50.0*Le) then
               pnewdt = 0.5
               write(*,*) '---------- disp. cut-back ----------'  
               write(*,*) 'time = ', time                            
               return
            endif
         enddo
      enddo








      !----------------------------------------------------------------
      ! Begin the loop over body integration points
      !
      ! Obtain integration point local coordinates and weights
      !
      if(nNode.eq.4) then
         !
         ! gauss integration for a rectangular element
         !
         if(nInt.eq.4) then
            call xint2D4pt(xi,w,nIntPt) ! 4-pt integration, nInt=4 above
         elseif(nInt.eq.1) then
            call xint2D1pt(xi,w,nIntPt) ! 1-pt integration, nInt=1 above
         else
            write(*,*) 'Invalid number of int points, nInt=',nInt
            write(80,*) 'Invalid number of int points, nInt=',nInt
            call xit
         endif
      else
         write(*,*) 'Incorrect number of nodes: nNode.ne.4'
         write(80,*) 'Incorrect number of nodes: nNode.ne.4'
         call xit
      endif


      ! Loop over integration points
      !
      jj = 0 ! jj is used for tracking the state variables
      do intpt=1,nIntPt


         ! Obtain shape functions and their local gradients
         !
         if(nNode.eq.4) then
            call calcShape2DLinear(nIntPt,xi,intpt,sh,dshxi)
         else
            write(*,*) 'Incorrect number of nodes: nNode.ne.4'
            write(80,*) 'Incorrect number of nodes: nNode.ne.4'
            call xit
         endif
         

         ! Map shape functions from local to global reference coordinate system
         !
         if(mcrd.eq.2) then
            call mapShape2Da(nNode,dshxi,coords,dsh,detMapJ,stat)
            if(stat.eq.0) then
               pnewdt = 0.5
               write(*,*) 'detF.lt.zero in mapShape2Da reference'
               write(*,*) 'time = ', time
               return
            endif
         elseif(mcrd.eq.3) then
            call mapShape2D(nNode,dshxi,coords,dsh,detMapJ,stat)
            if(stat.eq.0) then
               pnewdt = 0.5
               write(*,*) 'detF.lt.zero in mapShape2D reference'
               write(*,*) 'time = ', time               
               return
            endif
         else
            ! big problem
            write(*,*) 'Unexpected error, mcrd=',mcrd
            write(80,*) 'Unexpected error, mcrd=',mcrd
            call xit
         endif




         ! Obtain the damage and its derivative's at 
         !  this intPt at the begining and end of the incrment
         !
         dm_tau = zero
         dm_t = zero
         dDMdt = zero
         dDMdX = zero
         do k=1,nNode
            dm_tau = dm_tau + dmNew(k)*sh(k)
            dm_t   = dm_t + dmOld(k)*sh(k)
            do i=1,nDim
               dDMdX(i,1) = dDMdX(i,1) + dmNew(k)*dsh(k,i)
            enddo
         enddo
         

         dDMdt = (dm_tau - dm_t)/dtime
         



         ! Obtain the deformation gradient at this integration point. 
         !
         F_tau = Iden
         F_t = Iden
         do i=1,nDim
            do j=1,nDim
               do k=1,nNode
                  F_tau(i,j) = F_tau(i,j) + dsh(k,j)*u(k,i)
                  F_t(i,j) = F_t(i,j) + dsh(k,j)*uOld(k,i)
               enddo
            enddo
         enddo




         ! Obtain state variables from previous increment
         !
         if((kinc.le.1).and.(kstep.eq.1)) then
            !
            ! this is the first increment, of the first step
            !  give initial conditions
            !
            blambda_t  = blambda0   ! Bond-stretch
            Histmax_t  = zero       ! History function
            F_t(3,3) = one          ! deformation gradient (3,3) component at time t
            F_tau(3,3) = one        ! deformation gradient (3,3) component at time tau
            !
         else
            !
            ! this is not the first increment, read old values
            !
            blambda_t  = svars(1+jj)   ! Bond-stretch
            Histmax_t  = svars(2+jj)   ! History function                  
            F_t(3,3) = svars(4+jj)     ! deformation gradient (3,3) component at time t
            F_tau(3,3) = svars(4+jj)   ! deformation gradient (3,3) component at time tau            
            !
         endif
         

         call mdet(F_t,detF_t) 
         call matInv3D(F_tau,Finv,detF,stat)


         ! Cut-back algorithm if det(F) < 0
         ! 
         if (stat.eq.0) then
            pnewdt = half
            write(*,*) 'detF.lt.zero before plane stress algorithm'
            write(*,*) 'time = ', time
            write(*,*) 'element = ', jelem
            write(*,*) 'F_tau = ', F_tau
            write(*,*) 'detF = ', detF
            return
         end if
         
         

         !----------------------------------------------------------------
         ! 
         !  Reference for the plane-stress algorithm:
         !  Klinkel and Govindjee, 2002.
         !  Using finite strain 3D-material models in beam and shell 
         !  elements. Engineering Computations, 19, 254-271. 
         ! 
         ! 
         ! Calculate F(3,3) at time tau satisfying |TR(3,3)| < tolerance
         ! 
         
         call planestress(sol,time,jelem,props,nprops,
     +     dtime,pnewdt,F_t,F_tau,dm_t,dm_tau,
     +     blambda_t,detF_t,theta0,Histmax_t)
     
     
     
         ! Cut-back algorithm if planestress algorithm exceeds maximum iteration
         !     
         if (pnewdt.eq.half) then
            return
         end if
     
         
         ! Update F(3,3) obtained from the planestress algorithm         
         !          
         F_tau(3,3) = sol


         
     
         ! Calculate stress, tangent modulus, and internal variables
         ! (bond-stretch and history function) at each integration point
         ! 
       
         call integAB(props,nprops,dtime,time,jelem,
     +     F_t,F_tau,dm_t,dm_tau,blambda_t,detF_t,theta0,
     +     TR_tau,T_tau,dTRdF,
     +     blambda_tau,root,
     +     Histmax_t,Histmax_tau)     
              
     

         ! Save the state variables at this integ point
         !  at the end of the increment
         !
         S22_tau = T_tau(2,2)
         svars(1+jj) = blambda_tau   ! bond stretch
         svars(2+jj) = Histmax_tau   ! History function         
         svars(3+jj) = S22_tau       ! True stress
         svars(4+jj) = F_tau(3,3)    ! deformation gradient (3,3) component       

       
         jj = jj + nlSdv ! setup for the next intPt


         ! Save the state variables at this integ point in the
         !  global array used for plotting field output
         !
         globalSdv(jelem,intPt,1) = blambda_tau   ! bond stretch
         globalSdv(jelem,intPt,2) = Histmax_tau   ! History field function
         globalSdv(jelem,intPt,3) = S22_tau       ! True stress
         globalSdv(jelem,intPt,4) = F_tau(3,3)    ! deformation gradient (3,3) component       

 




         ! Compute/update the displacement residual vector
         !
         Smat(1,1) = TR_tau(1,1)
         Smat(2,1) = TR_tau(2,1)
         Smat(3,1) = TR_tau(1,2)
         Smat(4,1) = TR_tau(2,2)
         !
         Gmat = zero
         do kk=1,nNode
            Gmat(1,1+nDim*(kk-1)) = dsh(kk,1)
            Gmat(2,2+nDim*(kk-1)) = dsh(kk,1)
            Gmat(3,1+nDim*(kk-1)) = dsh(kk,2)
            Gmat(4,2+nDim*(kk-1)) = dsh(kk,2)
         enddo      
         !
         BodyForceRes = zero
         do kk=1,nNode
            BodyForceRes(1+nDim*(kk-1),1) = sh(kk)*body(1)
            BodyForceRes(2+nDim*(kk-1),1) = sh(kk)*body(2)
         enddo
         !
         Ru = Ru + detmapJ*w(intpt)*
     +        (
     +        -matmul(transpose(Gmat),Smat)
     +        + BodyForceRes
     +        )




         ! Compute/update the damage residual vector
         !
         do kk=1,nNode
            Nvec(1,kk) = sh(kk)
         enddo


         ResFac = zeta*dDMdt+erff*dm_tau-two*(one-dm_tau)*Histmax_tau
         !
         Rd = Rd - detmapJ*w(intpt)*
     +        (
     +        transpose(Nvec)*ResFac + (erff*lc**two)*matmul(dsh,dDMdX)
     +        )
 
         
         
         

      
      
      

         ! Compute/update the displacement tangent matrix
         !
         Amat = zero
         Amat(1,1) = dTRdF(1,1,1,1)
         Amat(1,2) = dTRdF(1,1,2,1)
         Amat(1,3) = dTRdF(1,1,1,2)
         Amat(1,4) = dTRdF(1,1,2,2)
         Amat(2,1) = dTRdF(2,1,1,1)
         Amat(2,2) = dTRdF(2,1,2,1)
         Amat(2,3) = dTRdF(2,1,1,2)
         Amat(2,4) = dTRdF(2,1,2,2)
         Amat(3,1) = dTRdF(1,2,1,1)
         Amat(3,2) = dTRdF(1,2,2,1)
         Amat(3,3) = dTRdF(1,2,1,2)
         Amat(3,4) = dTRdF(1,2,2,2)
         Amat(4,1) = dTRdF(2,2,1,1)
         Amat(4,2) = dTRdF(2,2,2,1)
         Amat(4,3) = dTRdF(2,2,1,2)
         Amat(4,4) = dTRdF(2,2,2,2)
         
         
         Amat13(1,1) = dTRdF(1,1,3,3)
         Amat13(1,2) = dTRdF(2,1,3,3)
         Amat13(1,3) = dTRdF(1,2,3,3)
         Amat13(1,4) = dTRdF(2,2,3,3)
         
         
         Amat31(1,1) = dTRdF(3,3,1,1)
         Amat31(1,2) = dTRdF(3,3,2,1)
         Amat31(1,3) = dTRdF(3,3,1,2)
         Amat31(1,4) = dTRdF(3,3,2,2)
         

         
         ! Static condensation for displacement tangent modulus
         ! 
         Amat = Amat - matmul(transpose(Amat13),Amat31)/dTRdF(3,3,3,3)
         
         


         Kuu = Kuu + detMapJ*w(intpt)*
     +        (
     +        matmul(matmul(transpose(Gmat),Amat),Gmat)
     +        )




         ! Compute/update the damage tangent matrix
         !
         TanFac = zeta/dtime + erff + two*Histmax_tau
         
         Kdd = Kdd + detmapJ*w(intPt)*
     +        (
     +        TanFac*matmul(transpose(Nvec),Nvec)
     +        + (erff*lc**two)*matmul(dsh,transpose(dsh))
     +        )



         ! Compute/update the damage - displacement tangent matrix
         ! Here, we neglect the coupled tangents
         !    
         Kdu = zero     


         ! Compute/update the displacement - damage tangent matrix
         ! Here, we neglect the coupled tangents
         ! 
         Kud = zero

          
      enddo
      !
      ! End the loop over body integration points
      !----------------------------------------------------------------


            

      !----------------------------------------------------------------
      ! Return Abaqus the RHS vector and the Stiffness matrix.
      !
      call AssembleElement(nDim,nNode,nDofEl,
     +     Ru,Rd,Kuu,Kud,Kdu,Kdd,
     +     rhs,amatrx)
      ! 
      ! End return of RHS and AMATRX
      !----------------------------------------------------------------


      return
      end subroutine UPS4
      
      
     
           
 


************************************************************************

      ! (Nearly incompressible) Arruda-Boyce model with bond-stretching
      ! Talamini and Anand JMPS 2018
      ! 
      subroutine integAB(props,nprops,dtime,time,jelem,
     +     F_t,F_tau,dm_t,dm_tau,blambda_t,detF_t,theta,
     +     TR_tau,T_tau,dTRdF,
     +     blambda_tau,root,
     +     Histmax_t,Histmax_tau)

      ! This subroutine computes everything required for the time integration
      ! of the problem.
      ! 
      ! Inputs:
      !  1) material parameters, props(nprops)
      !  2) time increment, dtime
      !  3) deformation gradient, F_t(3,3) and F_tau(3,3)
      !  4) damage variable, d_t and d_tau
      !  5) bond-stretch at time t, blambda_t
      !  6) history function at time t, Histmax_t      
      ! 
      ! Outputs:
      !  1) first Piola stress, TR_tau(3,3) and Cauchy stress, T_tau(3,3)
      !  2) referential tangent modulus, dTRdF(3,3,3,3)
      !  3) bond-stretch at time tau, blambda_tau
      !  4) history function at time tau, Histmax_tau

      implicit none

      integer i,j,k,l,m,n,nprops,nargs,stat,jelem
      parameter(nargs=6)

      real*8 Iden(3,3),props(nprops),F_t(3,3),F_tau(3,3),dm_t,dm_tau
      real*8 theta,TR_tau(3,3),T_tau(3,3),dTRdF(3,3,3,3),Gshear,Kbulk
      real*8 spTanMod(3,3,3,3),detF,detF_t,FinvT(3,3),Finv(3,3)
      real*8 F_bar(3,3),B_bar(3,3),trB_bar,args(nargs),avglambda_t
      real*8 Gshear0,ENK,Eb_bar,erff,lc,zeta,kk
      real*8 avglambda,root,root_old,blambda_t,blambda_tau,PA,dPAdr
      real*8 dLb,dGdL,dGdLb,dLbdF,dLbddm,Histmax_t,Histmax_tau,dum1
      real*8 Ft_bar(3,3),B_t(3,3),trB_t,delta
      real*8 per(3,3),Fper(3,3),Finvper(3,3),FinvTper(3,3),detFper
      real*8 TRper(3,3),Fper_bar(3,3),Bper(3,3),trBper,avglambda_per
      real*8 Gshear_per,blambda_per,rootper,PAper,time,dtime
      
      real*8 zero,one,two,three,four,third,half
      parameter(zero=0.d0,one=1.d0,two=2.d0,three=3.d0,third=1.d0/3.d0,
     +     half=1.d0/2.d0,four=4.d0)


      ! Identity tensor
      !
      call onem(Iden)


      ! Obtain material properties
      !
      Gshear0 = props(1)       ! Shear modulus
      ENK     = props(2)       ! number of Kuhn segments in a chain
      Eb_bar  = props(3)       ! Bond stiffness
      Kbulk   = props(4)       ! Bulk modulus
      erff    = props(5)       ! chain scission energy per unit
      lc      = props(6)       ! characteristic length scale of the graident theory
      zeta    = props(7)       ! kinetic modulus for the evolution of the damage   
      
      kk     = 1d-4            ! degradation function is modified as (1-d)^2 + kk


      
      ! Compute the inverse of F, its determinant, and its transpose
      !
      call matInv3D(F_tau,Finv,detF,stat)
      FinvT = transpose(Finv)


      F_bar = (detF**(-third))*F_tau
      B_bar = matmul(F_bar,transpose(F_bar))
      trB_bar = B_bar(1,1) + B_bar(2,2) + B_bar(3,3)

      avglambda = dsqrt(trB_bar/three)      
      
      
      Ft_bar = (detF_t**(-third))*F_t
      B_t = matmul(Ft_bar,transpose(Ft_bar))
      trB_t = B_t(1,1) + B_t(2,2) + B_t(3,3)

      avglambda_t = dsqrt(trB_t/three)    



      
      ! Compute the bond-stretch 
      ! For convenience, define the root of implicit equation
      ! 
      root_old = avglambda_t/blambda_t/dsqrt(ENK)


      args(1)  = dm_tau
      args(2)  = avglambda
      args(3)  = ENK
      args(4)  = Gshear0
      args(5)  = Eb_bar
      args(6)  = kk
      call solveblambda(root,args,nargs,root_old,time,jelem)



      blambda_tau = avglambda/root/dsqrt(ENK)
      



      ! Update History function
      ! 
      dum1=half*(Eb_bar*(blambda_t-one)**two 
     +    + KBulk*(detF_t-one)**two)
     +    - erff/two

c      dum1 = half*(Eb_bar*(blambda_t-one)**two 
c     +    + KBulk*(detF_t-one)**two)


     
     
      if (dum1.gt.Histmax_t) then
          ! 
          ! if current history function is maximum, update the history function
          ! 
          Histmax_tau = dum1
      else 
          Histmax_tau = Histmax_t  
      end if




      ! Pade approaximation of inverse Langevin function
      ! 
      PA = root*(three-root**two)/(one-root**two)
      
      ! Shear modulus
      ! 
      Gshear = Gshear0/three*dsqrt(ENK)/(avglambda*blambda_tau)*PA

      ! Cauchy stress
      ! 
      T_tau = Gshear*(B_bar - avglambda**two*Iden)/detF
     +        + ((one-dm_tau)**two+kk)*Kbulk*(detF-one)*Iden
        
      ! first Piola stress
      !       
      TR_tau = Gshear*(detF**(-two*third)*F_tau-avglambda**two*FinvT)
     +         + ((one-dm_tau)**two+kk)*Kbulk*(detF**two-detF)*FinvT
     




      ! Derivatives for tangent modulus
      ! 
      dPAdr = (three+root**four)/(one-root**two)**two     
       
      dGdL = -Gshear/avglambda*(one - root/PA*dPAdr)    
     
      dGdLb = -Gshear/blambda_tau*(one + root/PA*dPAdr)
      
      dLb = ((one-dm_tau)**two+kk)*Eb_bar/Gshear0/ENK
     +       *(two*blambda_tau-one) 
     +       + (PA+root*dPAdr)*root/blambda_tau
     
     
      dLbdF=(PA+root*dPAdr)/dsqrt(ENK)/blambda_tau/dLb
     
     
      ! for use in coupled tangents
      ! 
c      dLbddm=two*(one-dm_tau)*Eb_bar/Gshear0/ENK*(blambda_tau-one)
c     +       *blambda_tau/dLb
     
     

      ! Compute material tangent modulus, dTRdF
      !
      dTRdF = zero
      do i=1,3
         do j = 1,3
            do k = 1,3
               do l = 1,3
                  dTRdF(i,j,k,l) = dTRdF(i,j,k,l)
c                 dTRdF|Lb,d
     +              + Gshear*(detF**(-two*third)*(Iden(i,k)*Iden(j,l)
     +              - two/three*F_tau(i,j)*Finv(l,k))
     +              + (trB_bar/three)*Finv(l,i)*Finv(j,k)
     +              - two/three*(detF**(-two/three)*F_tau(k,l)
     +              - (trB_bar/three)*Finv(l,k))*Finv(j,i))    
c     
     +              + dGdL*(detF**(-two/three)/three/avglambda
     +              *F_tau(k,l) - avglambda/three*Finv(l,k))
     +              *(detF**(-two/three)*F_tau(i,j)
     +              - (trB_bar/three)*Finv(j,i))
c     
     +              + ((one-dm_tau)**two+kk)*Kbulk*(
     +                (two*detF-one)*detF*Finv(l,k)*Finv(j,i)
     +              - (detF-one)*detF*Finv(l,i)*Finv(j,k))
c     
c                   dTRdLb|F,d * dLbdF|d  
     +              + dGdLb*(detF**(-two/three)*F_tau(i,j) 
     +              - (trB_bar/three)*Finv(j,i))*dLbdF
     +              *(detF**(-two/three)/three/avglambda*F_tau(k,l)
     +              - avglambda/three*Finv(l,k))
               enddo
           enddo
         enddo
      enddo



      

      return
      end subroutine integAB
      
      
      
      


****************************************************************************

      subroutine solveblambda(root,args,nargs,rootOld,time,jelem)

      ! This subroutine will numerically solve for the bond-stretch
      !  using rtsafe (in numerical recipes)

      implicit none

      integer maxit,j,nargs,jelem

      real*8 xacc,f,df,fl,fh,xl,xh,x1,x2,swap,root,dxold,one
      real*8 dx,args(nargs),zero,rootOld,temp,rootMax,rootMin
      real*8 time

      parameter(maxit=100)
      parameter(xacc=1.d-6,zero=0.d0,one=1.d0)

      rootMax = 1d0
      rootMin = 0.0001d0

      x1 = rootMin
      x2 = rootMax
      call blambdaFunc(x1,FL,DF,args,nargs)
      call blambdaFunc(x2,FH,DF,args,nargs)

      if(fl*fh.ge.zero) then
         root = rootOld
         write(*,*) 'time = ', time
         write(*,*) 'jelem = ', jelem         
         write(*,*) 'FYI, root not bracketed on blambda'
         write(*,*) 'fl=',fl
         write(*,*) 'fh=',fh
         write(*,*) 'rootOld=',rootOld
         write(80,*) 'FYI, the root is not bracketed on blambda'
         write(80,*) 'fl=',fl
         write(80,*) 'fh=',fh
         write(80,*) 'rootOld=',rootOld

         write(*,*) 'dm =',args(1)
         write(*,*) 'lambda =',args(2)
         write(*,*) 'ENK=',args(3)
         write(*,*) 'Gshear=',args(4)
         write(*,*) 'Eb_bar =',args(5)
         write(*,*) 'kk =',args(6)


         call xit
         return
      endif

C
C		ORIENT THE SEARCH SO THAT F(XL) < 0.
C
      IF( FL .LT. 0.D0 ) THEN
         XL = X1
         XH = X2
      ELSE
         XH = X1
         XL = X2
         SWAP = FL
         FL = FH
         FH = SWAP
      END IF
C
C		INITIALIZE THE GUESS FOR THE ROOT, THE ''STEP SIZE
C		BEFORE LAST'', AND THE LAST STEP
C
      if(rootOld.lt.rootMin) rootOld = rootMin
      if(rootOld.gt.rootMax) rootOld = rootMax
      ROOT = rootOld !0.5D0 *( X1 + X2)
      DXOLD = DABS(X2 - X1)
      DX = DXOLD
      
      call blambdaFunc(root,F,DF,args,nargs)

C
C			LOOP OVER ALLOWED ITERATIONS
C
      DO 10 J = 1,MAXIT
C
C			BISECT IF NEWTON OUT OF RANGE, OR NOT DECREASING
C			FAST ENOUGH.
C
         IF( ((ROOT-XH)*DF - F)*((ROOT - XL)*DF -F) .GE. 0.D0
     +        .OR. DABS(2.D0*F) .GT. DABS(DXOLD*DF) ) THEN

            DXOLD = DX
            DX = 0.5D0*(XH-XL)
            ROOT = XL + DX
            IF( XL .EQ. ROOT ) THEN
C
C			CHANGE IN ROOT IS NEGLIGIBLE
C
               RETURN
            END IF

         ELSE
C
C			NEWTON STEP IS ACCEPTABLE. TAKE IT.
C
            DXOLD = DX
            DX = F/DF
            TEMP = ROOT
            ROOT = ROOT - DX
            IF( TEMP .EQ. ROOT) THEN
C
C			 CHANGE IN ROOT IS NEGLIGIBLE
C
               RETURN
            END IF

         END IF
C
C		CONVERVEGENCE CRITERION
C
         IF( DABS(DX) .LT. XACC) RETURN

C
C			THE ONE NEW FUNCTION EVALUATION PER ITERATION
C
         call blambdaFunc(root,F,DF,args,nargs)

C
C		MAINTAIN THE BRACKET ON THE ROOT
C
         IF( F .LT. 0.D0) THEN
            XL = ROOT
            FL = F
         ELSE
            XH = ROOT
            FH = F
         END IF

 10   CONTINUE

      WRITE(*,'(/1X,A)') 'solveblambda EXCEEDING MAXIMUM ITERATIONS'
      WRITE(80,'(/1X,A)') 'solveblambda EXCEEDING MAXIMUM ITERATIONS'

      return
      end subroutine solveblambda
      



****************************************************************************

      subroutine blambdaFunc(root,f,df,args,nargs)

      ! This subroutine serves as the function we would like to solve for
      !  the bond-stretch by finding the root satisfying `f=0'
      
      implicit none

      integer nargs,NeoHookean,Langevin,material
      parameter(NeoHookean=1,Langevin=2)

      real*8 args(nargs),f,df,dm,lambda,ENK,Eb_Bar,Gshear,kk
      real*8 blambda,PA,dPAdr,root

      real*8 zero,one,two,three,four
      parameter(zero=0.d0,one=1.d0,two=2.d0,three=3.d0,four=4.d0)

      
      ! Obtain relevant quantities
      !
      dm      = args(1)
      lambda  = args(2)
      ENK     = args(3)
      Gshear  = args(4)
      Eb_bar  = args(5)
      kk      = args(6)
    


      blambda = lambda/root/dsqrt(ENK) 
      PA =  root*(three-root**two)/(one-root**two)
      
      dPAdr = (three+root**four)/(one-root**two)**two



      ! Compute the residual
      ! 
      f = ((one-dm)**two+kk)*Eb_bar*(blambda-one)*blambda
     +   - Gshear*ENK*root*PA


      ! Compute the tangent
      !
      df = -((one-dm)**two+kk)*Eb_bar*(two*blambda-one)*blambda/root
     +      - Gshear*ENK*(PA + root*dPAdr)


      return
      end subroutine blambdaFunc
      
      
      
      
      
      
      

      
      
****************************************************************************      
      subroutine planestress(root,time,jelem,props,nprops,dtime,pnewdt,
     +   F_t,F_tau,dm_t,dm_tau,blambda_t,detF_t,theta,Histmax_t)

      ! This subroutine will numerically solve for (3,3) component of 
      ! deformation gradient, F33, using rtsafe (in numerical recipes)
      
      implicit none

      integer maxit,j,jelem,nprops

      real*8 tol,f,df,fl,fh,xl,xh,x1,x2,swap,root,dxold,one
      real*8 dx,zero,half,rootOld,temp,rootMax,rootMin
      real*8 time,pnewdt,F_root(3,3)
      
      real*8 F_min(3,3),TR_min(3,3),T_min(3,3),dTRdF_min(3,3,3,3)
      real*8 F_max(3,3),TR_max(3,3),T_max(3,3),dTRdF_max(3,3,3,3)
      real*8 F_tau(3,3),TR_tau(3,3),T_tau(3,3),dTRdF(3,3,3,3)
      real*8 F_t(3,3),dm_t,dm_tau,blambda_t,detF_t,theta
      real*8 props(nprops),dtime,Histmax_t,dum1,dum2,dum3

      
      
      parameter(maxit=1000)
      parameter(tol=1.d-8,zero=0.d0,one=1.d0,half=0.5d0)
      
      ! Tolerance changes depending on the damage variable
      ! 
c      if (dm_tau.lt. 0.3d0) then
c         tol=1.d-8
c      elseif((dm_tau.ge. 0.3d0).and.(dm_tau.lt. 0.8d0)) then
c         tol=1.d-3
c      else
c         tol=1.d0
c      end if
      
c      pnewdt = one

      rootMax = 1d0
      rootMin = 0.0001d0

      x1 = rootMin
      x2 = rootMax
      
      F_min = F_tau
      F_max = F_tau
      
      F_min(3,3) = rootMin
      F_max(3,3) = rootMax
      

       
         call integAB(props,nprops,dtime,time,jelem,
     +        F_t,F_min,dm_t,dm_tau,blambda_t,detF_t,theta,
     +        TR_min,T_min,dTRdF_min,
     +        dum1,dum2,Histmax_t,dum3)     
     
            
         call integAB(props,nprops,dtime,time,jelem,
     +        F_t,F_max,dm_t,dm_tau,blambda_t,detF_t,theta,
     +        TR_max,T_max,dTRdF_max,
     +        dum1,dum2,Histmax_t,dum3)          
           
      
      FL = TR_min(3,3)
      FH = TR_max(3,3)      
     
     
     
C
C		ORIENT THE SEARCH SO THAT F(XL) < 0.
C

      IF (dabs(FH).LT.tol) then
         root = rootMax
         RETURN
      end if


      IF( FL .LT. 0.D0 ) THEN
         XL = X1
         XH = X2
      ELSE
         XH = X1
         XL = X2
         SWAP = FL
         FL = FH
         FH = SWAP
      END IF
C
C		INITIALIZE THE GUESS FOR THE ROOT, THE ''STEP SIZE
C		BEFORE LAST'', AND THE LAST STEP
C
      
      rootOld = F_tau(3,3) 
      
      if(rootOld.lt.rootMin) rootOld = rootMin
      if(rootOld.gt.rootMax) rootOld = rootMax
      
      
      ROOT = F_tau(3,3) 
      DXOLD = DABS(X2 - X1)
      DX = DXOLD
      



      F_root = F_tau  
           
         call integAB(props,nprops,dtime,time,jelem,
     +           F_t,F_root,dm_t,dm_tau,blambda_t,detF_t,theta,
     +           TR_tau,T_tau,dTRdF,
     +           dum1,dum2,Histmax_t,dum3)    
      
      F = TR_tau(3,3)
      DF = dTRdF(3,3,3,3)         
      
      

C
C			LOOP OVER ALLOWED ITERATIONS
C
      DO 10 J = 1,MAXIT
C
C			BISECT IF NEWTON OUT OF RANGE, OR NOT DECREASING
C			FAST ENOUGH.
C
c         IF( ((ROOT-XH)*DF - F)*((ROOT - XL)*DF -F) .GT. 0.D0
c     +        .OR. DABS(2.D0*F) .GT. DABS(DXOLD*DF) ) THEN
         IF(DABS(2.D0*F) .GT. DABS(DXOLD*DF)) THEN         


            
            DXOLD = DX
            DX = 0.5D0*(XH-XL)
            ROOT = XL + DX
            
            F_root = F_tau
            F_root(3,3) = root
            
                        
            call integAB(props,nprops,dtime,time,jelem,
     +           F_t,F_root,dm_t,dm_tau,blambda_t,detF_t,theta,
     +           TR_tau,T_tau,dTRdF,
     +           dum1,dum2,Histmax_t,dum3)
        
                 
            IF(dabs(F).lt.tol) THEN
               RETURN
            END IF

         ELSE
C
C			NEWTON STEP IS ACCEPTABLE. TAKE IT.
C
           
            DXOLD = DX
            DX = F/DF
            TEMP = ROOT
            ROOT = ROOT - DX

            IF(dabs(F).lt.tol) THEN            
               RETURN
            END IF

         END IF




         F_root = F_tau
         F_root(3,3) = root

       
         call integAB(props,nprops,dtime,time,jelem,
     +           F_t,F_root,dm_t,dm_tau,blambda_t,detF_t,theta,
     +           TR_tau,T_tau,dTRdF,
     +           dum1,dum2,Histmax_t,dum3)     
      
      
         F = TR_tau(3,3)
         DF = dTRdF(3,3,3,3)     

C
C		MAINTAIN THE BRACKET ON THE ROOT
C
         IF( F .LT. 0.D0) THEN
            XL = ROOT
            FL = F
         ELSE
            XH = ROOT
            FH = F
         END IF

 10   CONTINUE
 
 
      ! Cut-back algorithm when exceeding maximum iteration
      !  
      pnewdt = half      

  
      write(*,*) 'time = ', time 
      write(*,*) 'jelem = ', jelem
      write(*,*) 'F(3,3) = ', root      
      write(*,*) 'TR(3,3) = ', TR_tau(3,3)
      write(*,*) 'Planestress EXCEEDING MAXIMUM ITERATIONS'

      return
      end subroutine planestress
      
      
      
      
      
      
      


************************************************************************
************************************************************************
************************************************************************
************************************************************************

      subroutine AssembleElement(nDim,nNode,ndofel,
     +     Ru,Rd,Kuu,Kud,Kdu,Kdd,
     +     rhs,amatrx)
      
      !
      ! Subroutine to assemble the local elements residual and tangent
      !

      implicit none

      integer i,j,k,l,m,n,A11,A12,B11,B12,nDim,nNode,nDofEl,nDofN

      real*8 Ru(nDim*nNode,1),Rd(nNode,1),Kuu(nDim*nNode,nDim*nNode)
      real*8 Kdd(nNode,nNode),Kud(nDim*nNode,nNode),rhs(ndofel,1)
      real*8 Kdu(nNode,nDim*nNode),amatrx(ndofel,ndofel)


      ! Total number of degrees of freedom per node
      !
      nDofN = nDofEl/nNode


      ! init
      !
      rhs(:,1) = 0.d0
      amatrx = 0.d0

      if(nDim.eq.2) then
         !
         ! Assemble the element level residual
         !
         do i=1,nNode
            A11 = nDofN*(i-1)+1
            A12 = nDim*(i-1)+1
            !
            ! displacement
            !
            rhs(A11,1) = Ru(A12,1)
            rhs(A11+1,1) = Ru(A12+1,1)
            !
            ! damage
            !
            rhs(A11+2,1) = Rd(i,1)
         enddo
         !
         ! Assemble the element level tangent matrix
         !
         do i=1,nNode
            do j=1,nNode
               A11 = nDofN*(i-1)+1
               A12 = nDim*(i-1)+1
               B11 = nDofN*(j-1)+1
               B12 = nDim*(j-1)+1
               !
               ! displacement
               !
               amatrx(A11,B11) = Kuu(A12,B12)
               amatrx(A11,B11+1) = Kuu(A12,B12+1)
               amatrx(A11+1,B11) = Kuu(A12+1,B12)
               amatrx(A11+1,B11+1) = Kuu(A12+1,B12+1)
               !
               ! damage
               !
               amatrx(A11+2,B11+2) = Kdd(i,j)
               !
               ! displacement - damage
               !
               amatrx(A11,B11+2) = Kud(A12,j)
               amatrx(A11+1,B11+2) = Kud(A12+1,j)
               !
               ! damage - displacement
               !
               amatrx(A11+2,B11) = Kdu(i,B12)
               amatrx(A11+2,B11+1) = Kdu(i,B12+1)
               !
            enddo
         enddo
         !
      elseif(nDim.eq.3) then
         !
         ! Assemble the element level residual
         !
         do i=1,nNode
            A11 = nDofN*(i-1)+1
            A12 = nDim*(i-1)+1
            !
            ! displacement
            !
            rhs(A11,1)   = Ru(A12,1)
            rhs(A11+1,1) = Ru(A12+1,1)
            rhs(A11+2,1) = Ru(A12+2,1)
            !
            ! damage
            !
            rhs(A11+3,1) = Rd(i,1)
            !
         enddo
         !
         ! Assembly the element level tangent matrix
         !
         do i=1,nNode
            do j=1,nNode
               A11 = nDofN*(i-1)+1
               A12 = nDim*(i-1)+1
               B11 = nDofN*(j-1)+1
               B12 = nDim*(j-1)+1
               !
               ! displacement
               !
               amatrx(A11,B11)     = Kuu(A12,B12)
               amatrx(A11,B11+1)   = Kuu(A12,B12+1)
               amatrx(A11,B11+2)   = Kuu(A12,B12+2)
               amatrx(A11+1,B11)   = Kuu(A12+1,B12)
               amatrx(A11+1,B11+1) = Kuu(A12+1,B12+1)
               amatrx(A11+1,B11+2) = Kuu(A12+1,B12+2)
               amatrx(A11+2,B11)   = Kuu(A12+2,B12)
               amatrx(A11+2,B11+1) = Kuu(A12+2,B12+1)
               amatrx(A11+2,B11+2) = Kuu(A12+2,B12+2)
               !
               ! damage
               !
               amatrx(A11+3,B11+3) = Kdd(i,j)
               !
               ! displacement - damage
               !
               amatrx(A11,B11+3) = Kud(A12,j)
               amatrx(A11+1,B11+3) = Kud(A12+1,j)
               amatrx(A11+2,B11+3) = Kud(A12+2,j)
               !
               ! damage - displacement
               !
               amatrx(A11+3,B11) = Kdu(i,B12)
               amatrx(A11+3,B11+1) = Kdu(i,B12+1)
               amatrx(A11+3,B11+2) = Kdu(i,B12+2)
               !
            enddo
         enddo
         !
      else
         write(*,*) 'How did you get nDim=',nDim
         call xit
      endif

      return
      end subroutine AssembleElement

!****************************************************************************
!     Element subroutines
!****************************************************************************

      subroutine xint2D1pt(xi,w,nIntPt)
      !
      ! This subroutine will get the integration point locations
      !  and corresponding gauss quadrature weights for 2D elements
      !  using 1 gauss point for integration
      !
      !  xi(nIntPt,2): xi,eta coordinates for the integration pts
      !  w(nIntPt):    corresponding integration weights
      !
      implicit none
      !
      integer nIntPt,nDim
      !
      real*8 xi(1,2), w(1)


      ! Initialize
      !
      w = 0.d0
      xi = 0.d0


      ! Number of Gauss points
      !
      nIntPt = 1


      ! Gauss weights
      !
      w = 4.d0
      

      ! Gauss pt location in master element
      !
      xi(1,1) = 0.d0
      xi(1,2) = 0.d0


      return
      end subroutine xint2D1pt
      
!************************************************************************

      subroutine xint2D4pt(xi,w,nIntPt)
      !
      ! This subroutine will get the integration point locations
      !  and corresponding gauss quadrature weights for 2D elements
      !  using 4 gauss points for integration
      !
      !  xi(nIntPt,2): xi,eta coordinates for the integration pts
      !  w(nIntPt):    corresponding integration weights
      !
      implicit none
      !
      integer nIntPt,nDim
      !
      real*8 xi(4,2), w(4)


      ! Initialize
      !
      w = 0.d0
      xi = 0.d0


      ! Number of Gauss points
      !
      nIntPt = 4


      ! Gauss weights
      !
      w(1) = 1.d0
      w(2) = 1.d0
      w(3) = 1.d0
      w(4) = 1.d0
      

      ! Gauss pt locations in master element
      !
      xi(1,1) = -dsqrt(1.d0/3.d0)
      xi(1,2) = -dsqrt(1.d0/3.d0)
      xi(2,1) = dsqrt(1.d0/3.d0)
      xi(2,2) = -dsqrt(1.d0/3.d0)
      xi(3,1) = -dsqrt(1.d0/3.d0)
      xi(3,2) = dsqrt(1.d0/3.d0)
      xi(4,1) = dsqrt(1.d0/3.d0)
      xi(4,2) = dsqrt(1.d0/3.d0)


      return
      end subroutine xint2D4pt

************************************************************************

      subroutine xint3D1pt(xi,w,nIntPt)
      
      ! This subroutine will get the integration point locations
      !  and corresponding gauss quadrature weights for 3D elements
      !  using a 2 gauss points for integration
      !
      !  xi(nIntPt,3): xi,eta,zeta coordinates for the integration pts
      !  w(nIntPt):    corresponding integration weights
      
      implicit none

      integer nIntPt,nDim

      real*8 xi(1,3),w(1)


      ! Init
      !
      w = 0.d0
      xi = 0.d0


      ! Number of Gauss points
      !
      nIntPt = 1


      ! Gauss weights
      !
      w(1) = 8.d0
      

      ! Gauss pt locations in master element
      !
      xi(1,1) = 0.d0
      xi(1,2) = 0.d0
      xi(1,3) = 0.d0

      return
      end subroutine xint3D1pt
     
************************************************************************

      subroutine xint3D8pt(xi,w,nIntPt)
      
      ! This subroutine will get the integration point locations
      !  and corresponding gauss quadrature weights for 3D elements
      !  using 8 gauss points for integration
      !
      !  xi(nIntPt,3): xi,eta,zeta coordinates for the integration pts
      !  w(nIntPt):    corresponding integration weights
      
      implicit none

      integer nIntPt,nDim

      real*8 xi(8,3),w(8)


      ! Init
      !
      w = 0.d0
      xi = 0.d0


      ! Number of Gauss points
      !
      nIntPt = 8


      ! Gauss weights
      !
      w(1) = 1.d0
      w(2) = 1.d0
      w(3) = 1.d0
      w(4) = 1.d0
      w(5) = 1.d0
      w(6) = 1.d0
      w(7) = 1.d0
      w(8) = 1.d0
      

      ! Gauss pt locations in master element
      !
      xi(1,1) = -dsqrt(1.d0/3.d0)
      xi(1,2) = -dsqrt(1.d0/3.d0)
      xi(1,3) = -dsqrt(1.d0/3.d0)
      xi(2,1) = dsqrt(1.d0/3.d0)
      xi(2,2) = -dsqrt(1.d0/3.d0)
      xi(2,3) = -dsqrt(1.d0/3.d0)
      xi(3,1) = -dsqrt(1.d0/3.d0)
      xi(3,2) = dsqrt(1.d0/3.d0)
      xi(3,3) = -dsqrt(1.d0/3.d0)
      xi(4,1) = dsqrt(1.d0/3.d0)
      xi(4,2) = dsqrt(1.d0/3.d0)
      xi(4,3) = -dsqrt(1.d0/3.d0)
      xi(5,1) = -dsqrt(1.d0/3.d0)
      xi(5,2) = -dsqrt(1.d0/3.d0)
      xi(5,3) = dsqrt(1.d0/3.d0)
      xi(6,1) = dsqrt(1.d0/3.d0)
      xi(6,2) = -dsqrt(1.d0/3.d0)
      xi(6,3) = dsqrt(1.d0/3.d0)
      xi(7,1) = -dsqrt(1.d0/3.d0)
      xi(7,2) = dsqrt(1.d0/3.d0)
      xi(7,3) = dsqrt(1.d0/3.d0)
      xi(8,1) = dsqrt(1.d0/3.d0)
      xi(8,2) = dsqrt(1.d0/3.d0)
      xi(8,3) = dsqrt(1.d0/3.d0)


      return
      end subroutine xint3D8pt

************************************************************************

      subroutine xintSurf2D1pt(face,xLocal,yLocal,w)

      ! This subroutine will get the integration point locations
      !  and corresponding gauss quadrature weights for 2D elements
      !  using 2 gauss points for surface integration
      !
      !  xLocal(nIntPt): x coordinates for the integration pts
      !  yLocal(nIntPt): y coordinates for the integration pts
      !  w(nIntPt):    corresponding integration weights

      implicit none

      integer face

      real*8 xLocal(1),yLocal(1),w(1),zero,one,two
      parameter(zero=0.d0,one=1.d0,two=2.d0)


      ! Gauss weights
      !
      w(1) = two
      

      ! Gauss pt locations in master element
      !
      if(face.eq.1) then
         xLocal(1) = zero
         yLocal(1) = -one
      elseif(face.eq.2) then
         xLocal(1) = one
         yLocal(1) = zero
      elseif(face.eq.3) then
         xLocal(1) = zero
         yLocal(1) = one
      elseif(face.eq.4) then
         xLocal(1) = -one
         yLocal(1) = zero
      else
         write(*,*) 'face.ne.1,2,3,4'
         write(80,*) 'face.ne.1,2,3,4'
         call xit
      endif

      end subroutine xintSurf2D1pt

************************************************************************

      subroutine xintSurf2D2pt(face,xLocal,yLocal,w)

      ! This subroutine will get the integration point locations
      !  and corresponding gauss quadrature weights for 2D elements
      !  using 2 gauss points for surface integration
      !
      !  xLocal(nIntPt): x coordinates for the integration pts
      !  yLocal(nIntPt): y coordinates for the integration pts
      !  w(nIntPt):    corresponding integration weights

      implicit none

      integer face

      real*8 xLocal(2),yLocal(2),w(2),one,three
      parameter(one=1.d0,three=3.d0)


      ! Gauss weights
      !
      w(1) = one
      w(2) = one
      

      ! Gauss pt locations in master element
      !
      if(face.eq.1) then
         xLocal(1) = -dsqrt(one/three)
         yLocal(1) = -one
         xLocal(2) = dsqrt(one/three)
         yLocal(2) = -one
      elseif(face.eq.2) then
         xLocal(1) = one
         yLocal(1) = -dsqrt(one/three)
         xLocal(2) = one
         yLocal(2) = dsqrt(one/three)
      elseif(face.eq.3) then
         xLocal(1) = -dsqrt(one/three)
         yLocal(1) = one
         xLocal(2) = dsqrt(one/three)
         yLocal(2) = one
      elseif(face.eq.4) then
         xLocal(1) = -one
         yLocal(1) = dsqrt(one/three)
         xLocal(2) = -one
         yLocal(2) = -dsqrt(one/three)
      else
         write(*,*) 'face.ne.1,2,3,4'
         write(80,*) 'face.ne.1,2,3,4'
         call xit
      endif

      end subroutine xintSurf2D2pt

************************************************************************

      subroutine xintSurf2D3pt(face,xLocal,yLocal,w)

      ! This subroutine will get the integration point locations
      !  and corresponding gauss quadrature weights for 2D elements
      !  using 2 gauss points for surface integration
      !
      !  xLocal(nIntPt): x coordinates for the integration pts
      !  yLocal(nIntPt): y coordinates for the integration pts
      !  w(nIntPt):    corresponding integration weights

      implicit none

      integer face

      real*8 xLocal(3),yLocal(3),w(3),zero,one,two,three,five,eight,nine
      parameter(zero=0.d0,one=1.d0,two=2.d0,three=3.d0,five=5.d0,
     +     eight=8.d0,nine=9.d0)


      ! Gauss weights
      !
      w(1) = five/nine
      w(2) = eight/nine
      w(3) = five/nine
      

      ! Gauss pt locations in master element
      !
      if(face.eq.1) then
         xLocal(1) = -dsqrt(three/five)
         yLocal(1) = -one
         xLocal(2) = zero
         yLocal(2) = -one
         xLocal(2) = dsqrt(three/five)
         yLocal(2) = -one
      elseif(face.eq.2) then
         xLocal(1) = one
         yLocal(1) = -dsqrt(three/five)
         xLocal(2) = one
         yLocal(2) = zero
         xLocal(3) = one
         yLocal(3) = dsqrt(three/five)
      elseif(face.eq.3) then
         xLocal(1) = -dsqrt(three/five)
         yLocal(1) = one
         xLocal(2) = zero
         yLocal(2) = one
         xLocal(3) = dsqrt(three/five)
         yLocal(3) = one
      elseif(face.eq.4) then
         xLocal(1) = -one
         yLocal(1) = dsqrt(three/five)
         xLocal(2) = -one
         yLocal(2) = zero
         xLocal(3) = -one
         yLocal(3) = -dsqrt(three/five)
      else
         write(*,*) 'face.ne.1,2,3,4'
         write(80,*) 'face.ne.1,2,3,4'
         call xit
      endif

      end subroutine xintSurf2D3pt

************************************************************************

      subroutine xintSurf3D1pt(face,xLocal,yLocal,zLocal,w)

      ! This subroutine will get the integration point locations
      !  and corresponding gauss quadrature weights for 3D elements
      !  using 1 gauss point for surface integration
      !
      !  xLocal(nIntPt): x coordinates for the integration pts
      !  yLocal(nIntPt): y coordinates for the integration pts
      !  zLocal(nIntPt): z coordinates for the integration pts
      !  w(nIntPt):    corresponding integration weights

      implicit none

      integer face

      real*8 xLocal(1),yLocal(1),zLocal(1),w(1),zero,one,four
      parameter(zero=0.d0,one=1.d0,four=4.d0)


      ! Gauss weights
      !
      w(1) = four
      

      ! Gauss pt locations in master element
      !
      if(face.eq.1) then
         xLocal(1) = zero
         yLocal(1) = zero
         zLocal(1) = -one
      elseif(face.eq.2) then
         xLocal(1) = zero
         yLocal(1) = zero
         zLocal(1) = one
      elseif(face.eq.3) then
         xLocal(1) = zero
         yLocal(1) = -one
         zLocal(1) = zero
      elseif(face.eq.4) then
         xLocal(1) = one
         yLocal(1) = zero
         zLocal(1) = zero
      elseif(face.eq.5) then
         xLocal(1) = zero
         yLocal(1) = one
         zLocal(1) = zero
      elseif(face.eq.6) then
         xLocal(1) = -one
         yLocal(1) = zero
         zLocal(1) = zero
      else
         write(*,*) 'face.ne.1,2,3,4,5,6'
         write(80,*) 'face.ne.1,2,3,4,5,6'
         call xit
      endif

      end subroutine xintSurf3D1pt

************************************************************************

      subroutine xintSurf3D4pt(face,xLocal,yLocal,zLocal,w)

      ! This subroutine will get the integration point locations
      !  and corresponding gauss quadrature weights for 3D elements
      !  using 4 gauss points for surface integration
      !
      !  xLocal(nIntPt): x coordinates for the integration pts
      !  yLocal(nIntPt): y coordinates for the integration pts
      !  yLocal(nIntPt): z coordinates for the integration pts
      !  w(nIntPt):    corresponding integration weights

      implicit none

      integer face

      real*8 xLocal(4),yLocal(4),zLocal(4),w(4),one,three
      parameter(one=1.d0,three=3.d0)


      ! Gauss weights
      !
      w(1) = one
      w(2) = one
      w(3) = one
      w(4) = one
      

      ! Gauss pt locations in master element
      !
      if(face.eq.1) then
         xLocal(1) = -dsqrt(one/three)
         yLocal(1) = -dsqrt(one/three)
         zLocal(1) = -one
         xLocal(2) = dsqrt(one/three)
         yLocal(2) = -dsqrt(one/three)
         zLocal(2) = -one
         xLocal(3) = dsqrt(one/three)
         yLocal(3) = dsqrt(one/three)
         zLocal(3) = -one
         xLocal(4) = -dsqrt(one/three)
         yLocal(4) = dsqrt(one/three)
         zLocal(4) = -one
      elseif(face.eq.2) then
         xLocal(1) = -dsqrt(one/three)
         yLocal(1) = -dsqrt(one/three)
         zLocal(1) = one
         xLocal(2) = dsqrt(one/three)
         yLocal(2) = -dsqrt(one/three)
         zLocal(2) = one
         xLocal(3) = dsqrt(one/three)
         yLocal(3) = dsqrt(one/three)
         zLocal(3) = one
         xLocal(4) = -dsqrt(one/three)
         yLocal(4) = dsqrt(one/three)
         zLocal(4) = one
      elseif(face.eq.3) then
         xLocal(1) = -dsqrt(one/three)
         yLocal(1) = -one
         zLocal(1) = -dsqrt(one/three)
         xLocal(2) = dsqrt(one/three)
         yLocal(2) = -one
         zLocal(2) = -dsqrt(one/three)
         xLocal(3) = dsqrt(one/three)
         yLocal(3) = -one
         zLocal(3) = dsqrt(one/three)
         xLocal(4) = -dsqrt(one/three)
         yLocal(4) = -one
         zLocal(4) = dsqrt(one/three)
      elseif(face.eq.4) then
         xLocal(1) = one
         yLocal(1) = -dsqrt(one/three)
         zLocal(1) = -dsqrt(one/three)
         xLocal(2) = one
         yLocal(2) = dsqrt(one/three)
         zLocal(2) = -dsqrt(one/three)
         xLocal(3) = one
         yLocal(3) = dsqrt(one/three)
         zLocal(3) = dsqrt(one/three)
         xLocal(4) = one
         yLocal(4) = -dsqrt(one/three)
         zLocal(4) = dsqrt(one/three)
      elseif(face.eq.5) then
         xLocal(1) = -dsqrt(one/three)
         yLocal(1) = one
         zLocal(1) = -dsqrt(one/three)
         xLocal(2) = dsqrt(one/three)
         yLocal(2) = one
         zLocal(2) = -dsqrt(one/three)
         xLocal(3) = dsqrt(one/three)
         yLocal(3) = one
         zLocal(3) = dsqrt(one/three)
         xLocal(4) = -dsqrt(one/three)
         yLocal(4) = one
         zLocal(4) = dsqrt(one/three)
      elseif(face.eq.6) then
         xLocal(1) = -one
         yLocal(1) = -dsqrt(one/three)
         zLocal(1) = -dsqrt(one/three)
         xLocal(2) = -one
         yLocal(2) = dsqrt(one/three)
         zLocal(2) = -dsqrt(one/three)
         xLocal(3) = -one
         yLocal(3) = dsqrt(one/three)
         zLocal(3) = dsqrt(one/three)
         xLocal(4) = -one
         yLocal(4) = -dsqrt(one/three)
         zLocal(4) = dsqrt(one/three)
      else
         write(*,*) 'face.ne.1,2,3,4,5,6'
         write(80,*) 'face.ne.1,2,3,4,5,6'
         call xit
      endif

      end subroutine xintSurf3D4pt
     
!************************************************************************

      subroutine calcShape2DLinear(nIntPt,xi_int,intpt,sh,dshxi)
      !
      ! Calculate the shape functions and their derivatives at the
      ! given integration point in the master element


      ! Calculate the shape functions and their derivatives at the
      ! given integration point in the master element
      !
      !                          eta
      !   4-----------3          |
      !   |           |          |
      !   |           |          |
      !   |           |          |
      !   |           |          |
      !   |           |          O--------- xi
      !   1-----------2        origin at center
      !
      !
      ! sh(i) = shape function of node i at the intpt.
      ! dshxi(i,j) = derivative wrt j direction of shape fn of node i
      !
      implicit none
      !
      integer intpt,nDim,nIntPt
      !
      real*8 xi_int(nIntPt,2),sh(4),dshxi(4,2),xi,eta
      !
      real*8 zero,one,fourth
      parameter(zero=0.d0,one=1.d0,fourth=1.d0/4.d0)
      

      ! Location in the master element
      !
      xi = xi_int(intpt,1)
      eta = xi_int(intpt,2)
      
      
      ! The shape functions
      !
      sh(1) = fourth*(one - xi)*(one - eta)
      sh(2) = fourth*(one + xi)*(one - eta)
      sh(3) = fourth*(one + xi)*(one + eta)
      sh(4) = fourth*(one - xi)*(one + eta)
      
      
      ! The first derivatives
      !
      dshxi(1,1) = -fourth*(one - eta)
      dshxi(1,2) = -fourth*(one - xi)
      dshxi(2,1) = fourth*(one - eta)
      dshxi(2,2) = -fourth*(one + xi)
      dshxi(3,1) = fourth*(one + eta)
      dshxi(3,2) = fourth*(one + xi)
      dshxi(4,1) = -fourth*(one + eta)
      dshxi(4,2) = fourth*(one - xi)
      

      return
      end subroutine calcShape2DLinear

************************************************************************

      subroutine calcShape3DLinear(nIntPt,xi_int,intpt,sh,dshxi)
      !
      !
      ! Calculate the shape functions and their derivatives at the
      ! given integration point in the master element
      !
      ! This subroutine uses a 8-node linear 3D element as shown
      !
      !      8-----------7
      !     /|          /|       zeta
      !    / |         / |       
      !   5-----------6  |       |     eta
      !   |  |        |  |       |   /
      !   |  |        |  |       |  /
      !   |  4--------|--3       | /
      !   | /         | /        |/
      !   |/          |/         O--------- xi
      !   1-----------2        origin at cube center
      !
      !
      ! sh(i) = shape function of node i at the intpt.
      ! dshxi(i,j) = derivative wrt j direction of shape fn of node i
      ! d2shxi(i,j,k) = derivatives wrt j and k of shape fn of node i

      implicit none

      integer intpt,nDim,nIntPt,i,j

      real*8 xi_int(nIntPt,3),sh(8),dshxi(8,3)
      real*8 d2shxi(8,3,3),xi,eta,zeta

      real*8 zero,one,two,half,fourth,eighth
      parameter(zero=0.d0,one=1.d0,two=2.d0,half=0.5d0,fourth=0.25d0,
     +     eighth=1.d0/8.d0)
      

      ! Location in the master element
      !
      xi = xi_int(intpt,1)
      eta = xi_int(intpt,2)
      zeta = xi_int(intpt,3)
      !
      ! The shape functions
      !
      sh(1) = eighth*(one - xi)*(one - eta)*(one - zeta)
      sh(2) = eighth*(one + xi)*(one - eta)*(one - zeta)
      sh(3) = eighth*(one + xi)*(one + eta)*(one - zeta)
      sh(4) = eighth*(one - xi)*(one + eta)*(one - zeta)
      sh(5) = eighth*(one - xi)*(one - eta)*(one + zeta)
      sh(6) = eighth*(one + xi)*(one - eta)*(one + zeta)
      sh(7) = eighth*(one + xi)*(one + eta)*(one + zeta)
      sh(8) = eighth*(one - xi)*(one + eta)*(one + zeta)
      !
      ! The first derivatives
      !
      dshxi(1,1) = -eighth*(one - eta)*(one - zeta)
      dshxi(1,2) = -eighth*(one - xi)*(one - zeta)
      dshxi(1,3) = -eighth*(one - xi)*(one - eta)
      dshxi(2,1) = eighth*(one - eta)*(one - zeta)
      dshxi(2,2) = -eighth*(one + xi)*(one - zeta)
      dshxi(2,3) = -eighth*(one + xi)*(one - eta)
      dshxi(3,1) = eighth*(one + eta)*(one - zeta)
      dshxi(3,2) = eighth*(one + xi)*(one - zeta)
      dshxi(3,3) = -eighth*(one + xi)*(one + eta)
      dshxi(4,1) = -eighth*(one + eta)*(one - zeta)
      dshxi(4,2) = eighth*(one - xi)*(one - zeta)
      dshxi(4,3) = -eighth*(one - xi)*(one + eta)
      dshxi(5,1) = -eighth*(one - eta)*(one + zeta)
      dshxi(5,2) = -eighth*(one - xi)*(one + zeta)
      dshxi(5,3) = eighth*(one - xi)*(one - eta)
      dshxi(6,1) = eighth*(one - eta)*(one + zeta)
      dshxi(6,2) = -eighth*(one + xi)*(one + zeta)
      dshxi(6,3) = eighth*(one + xi)*(one - eta)
      dshxi(7,1) = eighth*(one + eta)*(one + zeta)
      dshxi(7,2) = eighth*(one + xi)*(one + zeta)
      dshxi(7,3) = eighth*(one + xi)*(one + eta)
      dshxi(8,1) = -eighth*(one + eta)*(one + zeta)
      dshxi(8,2) = eighth*(one - xi)*(one + zeta)
      dshxi(8,3) = eighth*(one - xi)*(one + eta)
      !
      ! The second derivatives
      !
      d2shxi = zero
      d2shxi(1,1,2) = eighth*(one - zeta)
      d2shxi(1,2,1) = d2shxi(1,1,2)
      d2shxi(1,1,3) = eighth*(one - eta)
      d2shxi(1,3,1) = d2shxi(1,1,3)
      d2shxi(1,2,3) = eighth*(one - xi)
      d2shxi(1,3,2) = d2shxi(1,2,3)
      d2shxi(2,1,2) = -eighth*(one - zeta)
      d2shxi(2,2,1) = d2shxi(2,1,2)
      d2shxi(2,1,3) = -eighth*(one - eta)
      d2shxi(2,3,1) = d2shxi(2,1,3)
      d2shxi(2,2,3) = eighth*(one + xi)
      d2shxi(2,3,2) = d2shxi(2,2,3)
      d2shxi(3,1,2) = eighth*(one - zeta)
      d2shxi(3,2,1) = d2shxi(2,1,2)
      d2shxi(3,1,3) = -eighth*(one + eta)
      d2shxi(3,3,1) = d2shxi(2,1,3)
      d2shxi(3,2,3) = -eighth*(one + xi)
      d2shxi(3,3,2) = d2shxi(2,2,3)
      d2shxi(4,1,2) = -eighth*(one - zeta)
      d2shxi(4,2,1) = d2shxi(2,1,2)
      d2shxi(4,1,3) = eighth*(one + eta)
      d2shxi(4,3,1) = d2shxi(2,1,3)
      d2shxi(4,2,3) = -eighth*(one - xi)
      d2shxi(4,3,2) = d2shxi(2,2,3)
      d2shxi(5,1,2) = eighth*(one + zeta)
      d2shxi(5,2,1) = d2shxi(2,1,2)
      d2shxi(5,1,3) = -eighth*(one - eta)
      d2shxi(5,3,1) = d2shxi(2,1,3)
      d2shxi(5,2,3) = -eighth*(one - xi)
      d2shxi(5,3,2) = d2shxi(2,2,3)
      d2shxi(6,1,2) = eighth*(one + zeta)
      d2shxi(6,2,1) = d2shxi(2,1,2)
      d2shxi(6,1,3) = eighth*(one - eta)
      d2shxi(6,3,1) = d2shxi(2,1,3)
      d2shxi(6,2,3) = -eighth*(one + xi)
      d2shxi(6,3,2) = d2shxi(2,2,3)
      d2shxi(7,1,2) = eighth*(one + zeta)
      d2shxi(7,2,1) = d2shxi(2,1,2)
      d2shxi(7,1,3) = eighth*(one + eta)
      d2shxi(7,3,1) = d2shxi(2,1,3)
      d2shxi(7,2,3) = eighth*(one + xi)
      d2shxi(7,3,2) = d2shxi(2,2,3)
      d2shxi(8,1,2) = -eighth*(one + zeta)
      d2shxi(8,2,1) = d2shxi(2,1,2)
      d2shxi(8,1,3) = -eighth*(one + eta)
      d2shxi(8,3,1) = d2shxi(2,1,3)
      d2shxi(8,2,3) = eighth*(one - xi)
      d2shxi(8,3,2) = d2shxi(2,2,3)
      
      return
      end subroutine calcShape3DLinear

!************************************************************************


      subroutine computeSurf(xLocal,yLocal,face,coords,sh,ds)

      ! This subroutine computes the shape functions, derivatives
      !  of shape functions, and the length ds, so that one can
      !  do the numerical integration on the boundary for fluxes 
      !  on the 4-node quadrilateral elements

      implicit none

      integer face

      real*8 xLocal,yLocal,ds,dshxi(4,2),sh(4),dXdXi,dXdEta,dYdXi
      real*8 dYdEta,one,coords(2,4),fourth,shape,normal(2,1)
      parameter(one=1.d0,fourth=1.d0/4.d0)

      sh(1) = fourth*(one - xLocal)*(one - yLocal)
      sh(2) = fourth*(one + xLocal)*(one - yLocal)
      sh(3) = fourth*(one + xLocal)*(one + yLocal)
      sh(4) = fourth*(one - xLocal)*(one + yLocal)
      
      dshxi(1,1) = -fourth*(one - yLocal)
      dshxi(1,2) = -fourth*(one - xLocal)
      dshxi(2,1) = fourth*(one - yLocal)
      dshxi(2,2) = -fourth*(one + xLocal)
      dshxi(3,1) = fourth*(one + yLocal)
      dshxi(3,2) = fourth*(one + xLocal)
      dshxi(4,1) = -fourth*(one + yLocal)
      dshxi(4,2) = fourth*(one - xLocal)

      dXdXi = dshxi(1,1)*coords(1,1)+dshxi(2,1)*coords(1,2)
     +     + dshxi(3,1)*coords(1,3)+dshxi(4,1)*coords(1,4)
      dXdEta = dshxi(1,2)*coords(1,1)+dshxi(2,2)*coords(1,2)
     +     + dshxi(3,2)*coords(1,3)+dshxi(4,2)*coords(1,4)
      dYdXi = dshxi(1,1)*coords(2,1)+dshxi(2,1)*coords(2,2)
     +     + dshxi(3,1)*coords(2,3)+dshxi(4,1)*coords(2,4)
      dYdEta = dshxi(1,2)*coords(2,1)+dshxi(2,2)*coords(2,2)
     +     + dshxi(3,2)*coords(2,3)+dshxi(4,2)*coords(2,4)


      ! Jacobian of the mapping
      !
      if((face.eq.2).or.(face.eq.4)) then
         ds = dsqrt(dXdEta*dXdEta + dYdEta*dYdEta)
      elseif((face.eq.1).or.(face.eq.3)) then
         ds = dsqrt(dXdXi*dXdXi + dYdXi*dYdXi)
      else
         write(*,*) 'never should get here'
         call xit
      endif


      ! Surface normal, outward pointing in this case. Useful for
      !  ``follower'' type loads. The normal is referential or spatial
      !  depending on which coords were supplied to this subroutine
      !  (NOT fully tested)
      !
      if((face.eq.2).or.(face.eq.4)) then
         normal(1,1) = dYdEta/dsqrt(dXdEta*dXdEta + dYdEta*dYdEta)
         normal(2,1) = -dXdEta/dsqrt(dXdEta*dXdEta + dYdEta*dYdEta)
         if(face.eq.4) normal = -normal
      elseif((face.eq.1).or.(face.eq.3)) then
         normal(1,1) = dYdXi/dsqrt(dXdXi*dXdXi + dYdXi*dYdXi)
         normal(2,1) = -dXdXi/dsqrt(dXdXi*dXdXi + dYdXi*dYdXi)
         if(face.eq.3) normal = -normal
      else
         write(*,*) 'never should get here'
         call xit
      endif

      return
      end subroutine computeSurf

************************************************************************

      subroutine computeSurf3D(xLocal,yLocal,zLocal,face,coords,sh,dA)

      ! This subroutine computes the shape functions, derivatives
      !  of shape functions, and the area dA, so that one can
      !  do the numerical integration on the boundary for fluxes 
      !  on the 8-node brick elements

      implicit none

      integer face,stat,i,j,k

      real*8 xLocal,yLocal,zLocal,dA,dshxi(8,3),sh(8),zero,dsh(8,3),one
      real*8 coords(3,8),two,eighth,mapJ(3,3),mag,normal(3,1)

      real*8 dXdXi,dXdEta,dXdZeta,dYdXi,dYdEta,dYdZeta,dZdXi,dZdEta
      real*8 dZdZeta

      parameter(one=1.d0,two=2.d0,eighth=1.d0/8.d0,zero=0.d0)

      ! The shape functions
      !
      sh(1) = eighth*(one - xLocal)*(one - yLocal)*(one - zLocal)
      sh(2) = eighth*(one + xLocal)*(one - yLocal)*(one - zLocal)
      sh(3) = eighth*(one + xLocal)*(one + yLocal)*(one - zLocal)
      sh(4) = eighth*(one - xLocal)*(one + yLocal)*(one - zLocal)
      sh(5) = eighth*(one - xLocal)*(one - yLocal)*(one + zLocal)
      sh(6) = eighth*(one + xLocal)*(one - yLocal)*(one + zLocal)
      sh(7) = eighth*(one + xLocal)*(one + yLocal)*(one + zLocal)
      sh(8) = eighth*(one - xLocal)*(one + yLocal)*(one + zLocal)


      ! Shape function derivatives
      !
      dshxi(1,1) = -eighth*(one - yLocal)*(one - zLocal)
      dshxi(1,2) = -eighth*(one - xLocal)*(one - zLocal)
      dshxi(1,3) = -eighth*(one - xLocal)*(one - yLocal)
      dshxi(2,1) = eighth*(one - yLocal)*(one - zLocal)
      dshxi(2,2) = -eighth*(one + xLocal)*(one - zLocal)
      dshxi(2,3) = -eighth*(one + xLocal)*(one - yLocal)
      dshxi(3,1) = eighth*(one + yLocal)*(one - zLocal)
      dshxi(3,2) = eighth*(one + xLocal)*(one - zLocal)
      dshxi(3,3) = -eighth*(one + xLocal)*(one + yLocal)
      dshxi(4,1) = -eighth*(one + yLocal)*(one - zLocal)
      dshxi(4,2) = eighth*(one - xLocal)*(one - zLocal)
      dshxi(4,3) = -eighth*(one - xLocal)*(one + yLocal)
      dshxi(5,1) = -eighth*(one - yLocal)*(one + zLocal)
      dshxi(5,2) = -eighth*(one - xLocal)*(one + zLocal)
      dshxi(5,3) = eighth*(one - xLocal)*(one - yLocal)
      dshxi(6,1) = eighth*(one - yLocal)*(one + zLocal)
      dshxi(6,2) = -eighth*(one + xLocal)*(one + zLocal)
      dshxi(6,3) = eighth*(one + xLocal)*(one - yLocal)
      dshxi(7,1) = eighth*(one + yLocal)*(one + zLocal)
      dshxi(7,2) = eighth*(one + xLocal)*(one + zLocal)
      dshxi(7,3) = eighth*(one + xLocal)*(one + yLocal)
      dshxi(8,1) = -eighth*(one + yLocal)*(one + zLocal)
      dshxi(8,2) = eighth*(one - xLocal)*(one + zLocal)
      dshxi(8,3) = eighth*(one - xLocal)*(one + yLocal)


      dXdXi = zero
      dXdEta = zero
      dXdZeta = zero
      dYdXi = zero
      dYdEta = zero
      dYdZeta = zero
      dZdXi = zero
      dZdEta = zero
      dZdZeta = zero
      do k=1,8
         dXdXi = dXdXi + dshxi(k,1)*coords(1,k)
         dXdEta = dXdEta + dshxi(k,2)*coords(1,k)
         dXdZeta = dXdZeta + dshxi(k,3)*coords(1,k)
         dYdXi = dYdXi + dshxi(k,1)*coords(2,k)
         dYdEta = dYdEta + dshxi(k,2)*coords(2,k)
         dYdZeta = dYdZeta + dshxi(k,3)*coords(2,k)
         dZdXi = dZdXi + dshxi(k,1)*coords(3,k)
         dZdEta = dZdEta + dshxi(k,2)*coords(3,k)
         dZdZeta = dZdZeta + dshxi(k,3)*coords(3,k)
      enddo


      ! Jacobian of the mapping
      !
      if((face.eq.1).or.(face.eq.2)) then
         ! zeta = constant on this face
         dA = dsqrt(
     +          (dYdXi*dZdEta - dYdEta*dZdXi)**two
     +        + (dXdXi*dZdEta - dXdEta*dZdXi)**two
     +        + (dXdXi*dYdEta - dXdEta*dYdXi)**two
     +        )
      elseif((face.eq.3).or.(face.eq.5)) then
         ! eta = constant on this face
         dA = dsqrt(
     +          (dYdXi*dZdZeta - dYdZeta*dZdXi)**two
     +        + (dXdXi*dZdZeta - dXdZeta*dZdXi)**two
     +        + (dXdXi*dYdZeta - dXdZeta*dYdXi)**two
     +        )
      elseif((face.eq.4).or.(face.eq.6)) then
         ! xi = constant on this face
         dA = dsqrt(
     +          (dYdEta*dZdZeta - dYdZeta*dZdEta)**two
     +        + (dXdEta*dZdZeta - dXdZeta*dZdEta)**two
     +        + (dXdEta*dYdZeta - dXdZeta*dYdEta)**two
     +        )
         else
            write(*,*) 'never should get here'
            call xit
      endif


      ! Surface normal, outward pointing in this case. Useful for
      !  ``follower'' type loads. The normal is referential or spatial
      !  depending on which coords were supplied to this subroutine
      !  (NOT fully tested)
      !
      if((face.eq.1).or.(face.eq.2)) then
         ! zeta = constant on this face
         normal(1,1) = dYdXi*dZdEta - dYdEta*dZdXi
         normal(2,1) = dXdXi*dZdEta - dXdEta*dZdXi
         normal(3,1) = dXdXi*dYdEta - dXdEta*dYdXi
         if(face.eq.1) normal = -normal
      elseif((face.eq.3).or.(face.eq.5)) then
         ! eta = constant on this face
         normal(1,1) = dYdXi*dZdZeta - dYdZeta*dZdXi
         normal(2,1) = dXdXi*dZdZeta - dXdZeta*dZdXi
         normal(3,1) = dXdXi*dYdZeta - dXdZeta*dYdXi
         if(face.eq.5) normal = -normal
      elseif((face.eq.4).or.(face.eq.6)) then
         ! xi = constant on this face
         normal(1,1) = dYdEta*dZdZeta - dYdZeta*dZdEta
         normal(2,1) = dXdEta*dZdZeta - dXdZeta*dZdEta
         normal(3,1) = dXdEta*dYdZeta - dXdZeta*dYdEta
         if(face.eq.6) normal = -normal
      else
         write(*,*) 'never should get here'
         call xit
      endif
      mag = dsqrt(normal(1,1)**two+normal(2,1)**two+normal(3,1)**two)
      normal(1,1) = normal(1,1)/mag
      normal(2,1) = normal(2,1)/mag
      normal(3,1) = normal(3,1)/mag

      end subroutine computeSurf3D

************************************************************************

      subroutine mapShape2D(nNode,dshxi,coords,dsh,detMapJ,stat)
      !
      ! Map derivatives of shape fns from xi-eta-zeta domain
      !  to x-y-z domain.
      !
      implicit none
      !
      integer i,j,k,nNode,ieror,stat
      !
      real*8 dshxi(nNode,2),dsh(nNode,2),coords(3,nNode),mapJ(2,2),
     +  mapJ_inv(2,2),detmapJ
      !
      real*8 zero,one,two,half,fourth,eighth
      parameter(zero=0.d0,one=1.d0,two=2.d0,half=0.5d0,fourth=0.25d0,
     +     eighth=1.d0/8.d0)


      ! Calculate the mapping Jacobian matrix:
      !
      mapJ = zero
      do i=1,2
        do j=1,2
          do k=1,nNode
              mapJ(i,j) = mapJ(i,j) + dshxi(k,i)*coords(j,k)
          end do
        end do
      end do


      ! Calculate the inverse and the determinant of Jacobian
      !
      call matInv2D(mapJ,mapJ_inv,detMapJ,stat)
      if(stat.eq.0) then
         write(*,*) 'Problem: detF.lt.zero in mapShape2D'
c         call xit
      endif


      ! Calculate first derivatives wrt x, y, z
      !
      dsh = transpose(matmul(mapJ_inv,transpose(dshxi)))
      

      return
      end subroutine mapShape2D

!*************************************************************************

      subroutine mapShape2Da(nNode,dshxi,coords,dsh,detMapJ,stat)
      !
      ! Map derivatives of shape fns from xi-eta-zeta domain
      !  to x-y-z domain.
      !
      ! This subroutine is exactly the same as the regular mapShape2D
      !  with the exception that coords(2,nNode) here and coords(3,nNode)
      !  in the regular.  I have noticed that a "heat transfer" and 
      !  "static" step uses MCRD=2, but for "coupled-temperature-displacement"
      !  you will get MCRD=3, even for a plane analysis.
      !
      implicit none
      !
      integer i,j,k,nNode,ieror,stat
      !
      real*8 dshxi(nNode,2),dsh(nNode,2),coords(2,nNode),mapJ(2,2),
     +  mapJ_inv(2,2),detmapJ
      !
      real*8 zero,one,two,half,fourth,eighth
      parameter(zero=0.d0,one=1.d0,two=2.d0,half=0.5d0,fourth=0.25d0,
     +     eighth=1.d0/8.d0)


      ! Calculate the mapping Jacobian matrix:
      !
      mapJ = zero
      do i=1,2
        do j=1,2
          do k=1,nNode
              mapJ(i,j) = mapJ(i,j) + dshxi(k,i)*coords(j,k)
          end do
        end do
      end do


      ! Calculate the inverse and the determinant of Jacobian
      !
      call matInv2D(mapJ,mapJ_inv,detMapJ,stat)
      if(stat.eq.0) then
         write(*,*) 'Problem: detF.lt.zero in mapShape2Da'
c         call xit
      endif


      ! Calculate first derivatives wrt x, y, z
      !
      dsh = transpose(matmul(mapJ_inv,transpose(dshxi)))
      

      return
      end subroutine mapShape2Da

************************************************************************

      subroutine mapShape3D(nNode,dshxi,coords,dsh,detMapJ,stat)
      !
      ! Map derivatives of shape fns from xi-eta-zeta domain
      !  to x-y-z domain.  This subroutine works for both 8-node
      !  linear and 20-node quadratic 3D elements.
      !
      implicit none

      integer i,j,k,nNode,ieror,stat

      real*8 dshxi(nNode,3),dsh(nNode,3),coords(3,nNode)
      real*8 mapJ(3,3),mapJ_inv(3,3),detmapJ

      real*8 zero,one,two,half,fourth,eighth
      parameter(zero=0.d0,one=1.d0,two=2.d0,half=0.5d0,fourth=0.25d0,
     +     eighth=1.d0/8.d0)
      

      ! Calculate the mapping Jacobian matrix:
      !
      mapJ = zero
      do i=1,3
        do j=1,3
          do k=1,nNode
              mapJ(i,j) = mapJ(i,j) + dshxi(k,i)*coords(j,k)
          end do
        end do
      end do


      ! Calculate the inverse and the determinant of Jacobian
      !
      call matInv3D(mapJ,mapJ_inv,detMapJ,stat)
      if(stat.eq.0) then
         write(*,*) 'Problem: detF.lt.zero in mapShape3D'
c         call xit
      endif


      ! Calculate first derivatives wrt x, y, z
      !
      dsh = transpose(matmul(mapJ_inv,transpose(dshxi)))


      ! The second derivatives may be calculated.
      !

      return
      end subroutine mapShape3D

!****************************************************************************
!     Utility subroutines
!****************************************************************************

      subroutine matInv3D(A,A_inv,det_A,istat)
      !
      ! Returns A_inv, the inverse and det_A, the determinant
      ! Note that the det is of the original matrix, not the
      ! inverse
      !
      implicit none
      !
      integer istat
      !
      real*8 A(3,3),A_inv(3,3),det_A,det_A_inv


      istat = 1
      
      det_A = A(1,1)*(A(2,2)*A(3,3) - A(3,2)*A(2,3)) -
     +        A(2,1)*(A(1,2)*A(3,3) - A(3,2)*A(1,3)) +
     +        A(3,1)*(A(1,2)*A(2,3) - A(2,2)*A(1,3))
      
      if (det_A .le. 0.d0) then
        write(*,*) 'WARNING: subroutine matInv3D:'
        write(*,*) 'WARNING: det of mat=',det_A
        istat = 0
        return
      end if
          
      det_A_inv = 1.d0/det_A
        
      A_inv(1,1) = det_A_inv*(A(2,2)*A(3,3)-A(3,2)*A(2,3))
      A_inv(1,2) = det_A_inv*(A(3,2)*A(1,3)-A(1,2)*A(3,3))
      A_inv(1,3) = det_A_inv*(A(1,2)*A(2,3)-A(2,2)*A(1,3))
      A_inv(2,1) = det_A_inv*(A(3,1)*A(2,3)-A(2,1)*A(3,3))
      A_inv(2,2) = det_A_inv*(A(1,1)*A(3,3)-A(3,1)*A(1,3))
      A_inv(2,3) = det_A_inv*(A(2,1)*A(1,3)-A(1,1)*A(2,3))
      A_inv(3,1) = det_A_inv*(A(2,1)*A(3,2)-A(3,1)*A(2,2))
      A_inv(3,2) = det_A_inv*(A(3,1)*A(1,2)-A(1,1)*A(3,2))
      A_inv(3,3) = det_A_inv*(A(1,1)*A(2,2)-A(2,1)*A(1,2))
      

      return
      end subroutine matInv3D

!****************************************************************************

      subroutine matInv2D(A,A_inv,det_A,istat)
      !
      ! Returns A_inv, the inverse, and det_A, the determinant
      ! Note that the det is of the original matrix, not the
      ! inverse
      !
      implicit none
      !
      integer istat
      !
      real*8 A(2,2),A_inv(2,2),det_A,det_A_inv

      
      istat = 1
      
      det_A = A(1,1)*A(2,2) - A(1,2)*A(2,1)
        
      if (det_A .le. 0.d0) then
        write(*,*) 'WARNING: subroutine matInv2D:'
        write(*,*) 'WARNING: det of mat=',det_A
        istat = 0
        return
      end if
            
      det_A_inv = 1.d0/det_A
          
      A_inv(1,1) =  det_A_inv*A(2,2)
      A_inv(1,2) = -det_A_inv*A(1,2)
      A_inv(2,1) = -det_A_inv*A(2,1)
      A_inv(2,2) =  det_A_inv*A(1,1)


      return
      end subroutine matInv2D

!****************************************************************************

      subroutine mdet(A,det)
      !
      ! This subroutine calculates the determinant
      ! of a 3 by 3 matrix [A]
      !
      implicit none
      !
      real*8  A(3,3),det


      det = A(1,1)*A(2,2)*A(3,3) 
     +	  + A(1,2)*A(2,3)*A(3,1)
     +	  + A(1,3)*A(2,1)*A(3,2)
     +	  - A(3,1)*A(2,2)*A(1,3)
     +	  - A(3,2)*A(2,3)*A(1,1)
     +	  - A(3,3)*A(2,1)*A(1,2)


      return
      end subroutine mdet
	
!****************************************************************************

      subroutine onem(A)
      !
      ! This subroutine stores the identity matrix in the
      ! 3 by 3 matrix [A]
      !
      implicit none
      !
      integer i,j
      !
      real*8 A(3,3)


      do i=1,3
         do J=1,3
	    if (i .eq. j) then
              A(i,j) = 1.0
            else
              A(i,j) = 0.0
            end if
         end do
      end do


      return
      end subroutine onem

****************************************************************************
