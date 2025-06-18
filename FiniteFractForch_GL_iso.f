!cccRoutine to solve the fractional porepressure diffusion equation to use in a couple thermomechanic problem
!ccc    Gunda Sachin 26/7/23
!ccc    Fractional Forchheimers correction model. Refer Fract Forchy
! Permeability in current configuration is scalar*unitmat
C!======================================================================
      module global
      !  globalSdv(X,Y,Z)
      !   X - element pointer
      !   Y - integration point pointer
      !   Z - SDV pointer
      real*8, allocatable :: Rqdot(:,:,:,:)
      real*8, allocatable :: globalSdv(:,:,:)
      end module global
C=======================================================================
      SUBROUTINE UMATHT(U,DUDT,DUDG,FLUX,DFDT,DFDG,
     1 STATEV,TEMP,DTEMP,DTEMDX,TIME,DTIME,PREDEF,DPRED,
     2 CMNAME,NTGRD,NSTATV,PROPS,NPROPS,COORDS,PNEWDT,
     3 NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
C
      use global
      INCLUDE 'ABA_PARAM.INC'
C
      CHARACTER*80 CMNAME
      DIMENSION DUDG(NTGRD),FLUX(NTGRD),DFDT(NTGRD),
     1 DFDG(NTGRD,NTGRD),STATEV(NSTATV),DTEMDX(NTGRD),
     2 TIME(2),PREDEF(1),DPRED(1),PROPS(NPROPS),COORDS(3)
      real*8 bina,eps, k0, m0, k0R  !
      dimension bina(kinc+1),eps(kinc,ntgrd) 
      dimension dfgrd(3,3), C(3,3),Cinv(3,3), UnitMat(3,3), Finv(3,3)
c
      Dimension Flux1(ntgrd,1), Flux2(ntgrd,1), dfdq(ntgrd,ntgrd)
      Dimension Feq(ntgrd,1), gradp(ntgrd,1), dfdqinv(ntgrd,ntgrd)
      Dimension errorn(ntgrd,1), dfdg_coeff(ntgrd,ntgrd)
      Dimension dfdg_coeffinv(ntgrd,ntgrd), RF(3,3), gradpi(ntgrd,1)
      Dimension	per(3,3), perinv(3,3), Qdot(ntgrd,1), hist(ntgrd,1)
      Dimension RFQ(ntgrd,1), Finv_t(3,3), Ft(3,3), Finvi(3,3)
      DIMENSION Forchh_u(3,3),Forchh(ntgrd,ntgrd), dfgrdi(3,3)
      DIMENSION flux1i(3,1), q_coeff(3,3), q_coeffiv(3,3), q_rhs(3,1)
      dimension dqdnormdgradp(3,1),dffdgradp(3,1), dRFdgradp(3,1), 
      dimension dAdgradpQ(3,1), dbdgradp(3,1)
C
      PARAMETER(ZERO=0.D0, ONE=1.D0, TWO=2.D0, THREE=3.D0, 
     1            FOUR=4.D0, FIVE = 5.d0, SIX=6.D0)
C
      alpha   = PROPS(1)
      k0      = props(2)
      m0      = props(3)
      k0R     = props(4)
      phi_s   = props(5)
      c0      = props(6)
      c1      = props()
      c2      = props()
      f_t     = props()
      R_o     = props()
      rho_f   = props()   ! Fluid density
      vis     = props()
      tc      = props()
      numelem = props()
      nintp   = props()
c
      Unitmat = ZERO
      do i = 1,3
        UnitMat(i,i) = ONE
      end do

!      write(*,*)"UMATHT started"
CCC   ANALOGY with the thermal problem
Ccc   U is  dP/dT= U+(du/dp)dp
CCCC  TEMP=P
CCC   DTEMP=DP
CCC   DUDT=DUDP
      if(.not.allocated(globalSdv)) then
         allocate(globalSdv(numElem,nIntp,21))
      endif
	  
      if(.not.allocated(Rqdot)) then
         allocate(Rqdot(10000,numelem,nintp,NTGRD))
      endif
	  
C
      DUDT = 1.0d0
      DU = DUDT*DTEMP 
      U = U+DU
C
      DFGRD(1,1) = globalSdv(noel,npt,3)
      DFGRD(2,2) = globalSdv(noel,npt,4)
      DFGRD(3,3) = globalSdv(noel,npt,5)
      DFGRD(2,3) = globalSdv(noel,npt,6)
      DFGRD(1,3) = globalSdv(noel,npt,7)
      DFGRD(1,2) = globalSdv(noel,npt,8)
      DFGRD(2,1) = globalSdv(noel,npt,9)
      DFGRD(3,1) = globalSdv(noel,npt,10)
      DFGRD(3,2) = globalSdv(noel,npt,11)
      NSHR       = globalSdv(noel,npt,1) 
      DFGRDi(1,1) = globalSdv(noel,npt,12)
      DFGRDi(2,2) = globalSdv(noel,npt,13)
      DFGRDi(3,3) = globalSdv(noel,npt,14)
      DFGRDi(2,3) = globalSdv(noel,npt,15)
      DFGRDi(1,3) = globalSdv(noel,npt,16)
      DFGRDi(1,2) = globalSdv(noel,npt,17)
      DFGRDi(2,1) = globalSdv(noel,npt,18)
      DFGRDi(3,1) = globalSdv(noel,npt,19)
      DFGRDi(3,2) = globalSdv(noel,npt,20)
      deti        = globalSdv(noel,npt,21)       
C
C JACOBIAN AND DISTORTION TENSOR
C
      DET=DFGRD(1, 1)*DFGRD(2, 2)*DFGRD(3, 3)-DFGRD(1, 2)*
     1      DFGRD(2, 1)*DFGRD(3, 3)
     
C Dimension Check - 3D or 2D elements 			

      IF(NSHR.EQ.3) THEN
        DET = DET+DFGRD(1, 2)*DFGRD(2, 3)*DFGRD(3, 1)+DFGRD(1, 3)*
     1         DFGRD(3, 2)*DFGRD(2, 1)-DFGRD(1, 3)*DFGRD(3,1)*
     2         DFGRD(2, 2)-DFGRD(2, 3)*DFGRD(3, 2)*DFGRD(1, 1)
      END IF

      call matinv3(DFGRD,Finv)
      Finv_t = transpose(Finv)
      Ft = transpose(dfgrd)
      
      C = matmul(Ft,DFGRD)
      call matinv3(C, Cinv)

c     Permeability
      
      per_sc = (k0R*((det-phi_S)/(one-phi_s))**k0*exp(one/two*m0*
     1          (det**two-one)))
      per = per_sc*det*Cinv
!      call matinv3(per,perinv) 
      
c     Forchheimer's constants
      phi_f = (det - phi_s)/det
      A_eq = c0*rho_f*vis**c2*(per_sc/det)**(c2+one)*(phi_f/det)**(c1) 
      !      A_eq = (one-f_t)/f_t**two*phi_f**(-five)/R_o*rho_f*sqrt(per_sc/
!     1       (vis*phi_f))
!      Forchh = Unitmat*A_eq

      do i = 1, ntgrd
        FLUX1(i,1) = flux(i)
        FLUX2(i,1) = flux(i)
        gradp(i,1) = dtemdx(i)
      end do
      
      Flux2 = det*matmul(Finv,flux2)
      call matinv3(dfgrdi,Finvi)
      flux1i = deti*matmul(Finvi,flux2)
      gradpi = matmul(Ft,gradp)
!      Forchh_u = matmul(matmul(Ft,Forchh),Finv_t)
        
!      do while (any(abs(errorN)>tol))
       qd = -matmul(per,gradpi)/vis
       qdnorm = one/det*sqrt(sum(matmul(transpose(matmul(dfgrd,qd)),
     1          matmul(dfgrd,qd))))
       fric_fac = two/(one + sqrt(one + four*A_eq*qdnorm))
       RF = phi_f*vis*(one+qdnorm*A_eq*firc_fac)/per_sc
      
c!    History calculations
      eps(1:kinc,1:ntgrd) = RQdot(1:kinc,noel,npt,1:ntgrd)
      
      bina(1)=one
      hist = zero
      do k = 1,kinc-1  ! initial flux is assumed to be zero
        do j = 1,ntgrd
          hist(j,1) = hist(j,1) + bina(k)*eps(kinc-k,j)	
        end do
        bina(k+1) = (k+1-alpha)*bina(k)/(k+1)
      end do
      hist = alpha*tc**alpha*dtime**(one-alpha)*det*hist
      
      q_coeff = Rf/det*dfgrd + alpha*tc**alpha*dtime**(-alpha)/det*
     1          dfgrd
      q_rhs = phi_f*vis/per_sc*matmul(dfgrd,qd)/det - hist - alpha*
     1        tc**alpha*dtime**(-alpha)*RF*matmul(dfgrd,flux1i)/det
      call matinv3(q_coeff, q_coeffinv)
      flux2 = matmul(q_coeffinv,q_rhs)
      
c!    Calculate Qdot for current time step
      Qdot = (Flux2 - flux1i)/dtime

c!    Calculate SDV RQdot
      RF = phi_f*vis*(one+qdnorm*A_eq*fric_fac)/per_sc
      Rqdot(kinc,noel,npt,:) = zero
      do i = 1,ntgrd
        do j = 1,ntgrd
          Rqdot(kinc,noel,npt,i) = Rqdot(kinc,noel,npt,i) + one/det*
     1                             dfgrd(i,j)*qdot(j,1)*RF
        end do
      end do
      
C!   dq/dgrap calculation
      dfdg_coeff = Rf*(one + alpha*tc**alpha*dtime**(-alpha))*dfgrd/
     1             det      
      
      dqdnormdgradp = -one/(vis)*per_sc*matmul(dfgrd,qd)
      dffdgradp = one/(vis*det)*per_sc*matmul(dfgrd,qd)*four*A_eq/
     1            ((on + sqrt(one+four*a_eq*qdnorm))**two*sqrt(one +
     2             four*A_eq*qdnorm))
      dRFdgradp = phi_f*vis/per_sc*(dffdgradp*a_eq*qdnorm+ firc_fac*
     1            A_eq*dqdnormdgradp)
      dAdgradpQ = one/det*(one + alpha*tc**alpha*dtime**(-alpha))*
     1            matmul(matmul(dfgrd,flux2),transpose(dRFdgradp)) 
      dbdgradp = -alpha*tc**alpha*dtime**(-alpha)*matmul(matmul(
     1           dfgrd,flux1i),transpose(dRFdgradp)) + phi_f*vis/
     2           (det*per_sc)*matmul(dfgrd,dqdnormdgradp)
      
      call matinv3(dfdg_coeff, dfdg_coeffinv)
            
      DFDG = one/det*matmul(dfgrd,matmul(dfdg_coeffinv,dbdgradp-
     1      dAdgradpQ))      
        
      STATEV(1) = sqrt(FLUX2(1,1)**2+FLUX2(2,1)**2+FLUX2(3,1)**2) 
      
      FLUX1 = 1/det*matmul(DFGRD,FLUX2)

      do k = 1,ntgrd
        flux(k) = Flux1(k,1)
      end do

      STATEV(3:5) = FLUX(1:3)
      STATEV(2) = sqrt(FLUX(1)**2+FLUX(2)**2+FLUX(3)**2)
      
      if(noel.gt.numelem) then
        write(*,*) "Error: Given numelem",numelem
        write(*,*)  "noel", noel 
      end if   
      

      RETURN
      END

c ======================================================================
C ABAQUS Subroutine for Neo-Hookean Hyperelasticity
C-----------------------------------------------------------------------
      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,JSTEP,KINC)
C
      use global
      INCLUDE 'ABA_PARAM.INC'
C
      CHARACTER*80 CMNAME
      DIMENSION STRESS(NTENS),STATEV(NSTATV),
     1 DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3 PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3),
     4 JSTEP(4)
C-----------------------------------------------------------------------
C
        DIMENSION B(6), UnitMat(3,3)
C
        PARAMETER(ZERO=0.D0, ONE=1.D0, TWO=2.D0, THREE=3.D0, 
     1            FOUR=4.D0, SIX=6.D0)
     
        Real*8 lambda, mu
C
C ----------------------------------------------------------------
!      write(*,*)"umat started"
C
C INPUT
C
      mu      = PROPS(1)
      lambda  = PROPS(2)
      phi_S   = props(3)
      numelem = props(4)
      nintp   = props(5) 
      
      if(.not.allocated(globalSdv)) then
         allocate(globalSdv(numElem,nIntp,21))
      endif
      
      globalSdv(noel,npt,3) = DFGRD1(1,1)
      globalSdv(noel,npt,4) = DFGRD1(2,2)
      globalSdv(noel,npt,5) = DFGRD1(3,3)
      globalSdv(noel,npt,6) = DFGRD1(2,3)
      globalSdv(noel,npt,7) = DFGRD1(1,3)
      globalSdv(noel,npt,8) = DFGRD1(1,2)
      globalSdv(noel,npt,9) = DFGRD1(2,1)
      globalSdv(noel,npt,10) = DFGRD1(3,1)
      globalSdv(noel,npt,11) = DFGRD1(3,2)
      globalSdv(noel,npt,1) = NSHR
      globalSdv(noel,npt,12) = DFGRD0(1,1)
      globalSdv(noel,npt,13) = DFGRD0(2,2)
      globalSdv(noel,npt,14) = DFGRD0(3,3)
      globalSdv(noel,npt,15) = DFGRD0(2,3)
      globalSdv(noel,npt,16) = DFGRD0(1,3)
      globalSdv(noel,npt,17) = DFGRD0(1,2)
      globalSdv(noel,npt,18) = DFGRD0(2,1)
      globalSdv(noel,npt,19) = DFGRD0(3,1)
      globalSdv(noel,npt,20) = DFGRD0(3,2)
      
      Unitmat = ZERO
      do i = 1,3
        UnitMat(i,i) = ONE
      end do
C
C JACOBIAN 
C
        DET=DFGRD1(1, 1)*DFGRD1(2, 2)*DFGRD1(3, 3)-DFGRD1(1, 2)*
     1      DFGRD1(2, 1)*DFGRD1(3, 3)
C Dimension Check - 3D or 2D elements 			
        IF(NSHR.EQ.3) THEN
          DET=DET+DFGRD1(1, 2)*DFGRD1(2, 3)*DFGRD1(3, 1)+DFGRD1(1, 3)*
     1        DFGRD1(3, 2)*DFGRD1(2, 1)-DFGRD1(1, 3)*DFGRD1(3,1)*
     2        DFGRD1(2, 2)-DFGRD1(2, 3)*DFGRD1(3, 2)*DFGRD1(1, 1)
        END IF
        
        DET0=DFGRD0(1, 1)*DFGRD0(2, 2)*DFGRD0(3, 3)-DFGRD0(1, 2)*
     1      DFGRD0(2, 1)*DFGRD0(3, 3)
C Dimension Check - 3D or 2D elements 			
        IF(NSHR.EQ.3) THEN
          DET0=DET0+DFGRD0(1, 2)*DFGRD0(2, 3)*DFGRD0(3, 1)+DFGRD0(1, 3)
     1        *DFGRD0(3, 2)*DFGRD0(2, 1)-DFGRD0(1, 3)*DFGRD0(3,1)*
     2        DFGRD0(2, 2)-DFGRD0(2, 3)*DFGRD0(3, 2)*DFGRD0(1, 1)
        END IF
        
      globalSdv(noel,npt,21) = Det0
C 			
C CALCULATE DEVIATORIC LEFT CAUCHY-GREEN DEFORMATION TENSOR
C
        B(1)=DFGRD1(1, 1)**2+DFGRD1(1, 2)**2+DFGRD1(1, 3)**2
        B(2)=DFGRD1(2, 1)**2+DFGRD1(2, 2)**2+DFGRD1(2, 3)**2
        B(3)=DFGRD1(3, 3)**2+DFGRD1(3, 1)**2+DFGRD1(3, 2)**2
        B(4)=DFGRD1(1, 1)*DFGRD1(2, 1)+DFGRD1(1, 2)*DFGRD1(2, 2)+
     1          DFGRD1(1, 3)*DFGRD1(2, 3)
C Dimension Check - 3D or 2D elements ?
 	
        IF(NSHR.EQ.3) THEN 
          B(5)=DFGRD1(1, 1)*DFGRD1(3, 1)+DFGRD1(1, 2)*DFGRD1(3, 2)+
     1            DFGRD1(1, 3)*DFGRD1(3, 3)
          B(6)=DFGRD1(2, 1)*DFGRD1(3, 1)+DFGRD1(2, 2)*DFGRD1(3, 2)+
     1            DFGRD1(2, 3)*DFGRD1(3, 3)
        END IF			

C
C CALCULATE THE STRESS
C
        TRB=(B(1)+B(2)+B(3))/THREE
        
c        write(*,*)'stress',stress
c        write(*,*)'stress',DTEMP
        DO K1=1,NDI
          STRESS(K1) = phi_s*(mu*B(K1)-mu+lambda*log(det))/det
     1                 - (TEMP+DTEMP)      !EG*(BBAR(K1)-TRBBAR)+PR !- TEMP!!!!Need to check
        END DO
c        write(*,*)'stress2',stress
        DO K1=NDI+1,NDI+NSHR
          STRESS(K1) = phi_s*(mu*B(K1))/det !EG*BBAR(K1)
        END DO	
c        globalSdv(noel,npt,2) = det
        
C 			
C CALCULATE THE STIFFNESS
C
        EG   = phi_S*mu/(det)
        EG23 = phi_s*mu*TWO/(THREE*det)
        EK   = phi_S*lambda/(det)-(TEMP + DTEMP)
        DDSDDE(1, 1)= EG23*(B(1)+TRB)+EK
        DDSDDE(2, 2)= EG23*(B(2)+TRB)+EK
        DDSDDE(3, 3)= EG23*(B(3)+TRB)+EK
        DDSDDE(1, 2)=-EG23*(B(1)+B(2)-TRB)+EK
        DDSDDE(1, 3)=-EG23*(B(1)+B(3)-TRB)+EK
        DDSDDE(2, 3)=-EG23*(B(2)+B(3)-TRB)+EK
        DDSDDE(1, 4)= EG23*B(4)/TWO
        DDSDDE(2, 4)= EG23*B(4)/TWO
        DDSDDE(3, 4)=-EG23*B(4)
        DDSDDE(4, 4)= EG*(B(1)+B(2))/TWO
C  Dimension Check - 3D or 2D elements ? 	
        IF(NSHR.EQ.3) THEN
          DDSDDE(1, 5)= EG23*B(5)/TWO
          DDSDDE(2, 5)=-EG23*B(5)
          DDSDDE(3, 5)= EG23*B(5)/TWO
          DDSDDE(1, 6)=-EG23*B(6)
          DDSDDE(2, 6)= EG23*B(6)/TWO
          DDSDDE(3, 6)= EG23*B(6)/TWO
          DDSDDE(5, 5)= EG*(B(1)+B(3))/TWO
          DDSDDE(6, 6)= EG*(B(2)+B(3))/TWO
          DDSDDE(4,5)= EG*B(6)/TWO
          DDSDDE(4,6)= EG*B(5)/TWO
          DDSDDE(5,6)= EG*B(4)/TWO
        END IF		
        DO K1=1, NTENS
          DO K2=1, K1-1
            DDSDDE(K1, K2)=DDSDDE(K2, K1)
          END DO
        END DO    
        
        DDSDDT = ZERO
        DDSDDT(1) = -ONE
        DDSDDT(2) = -one
        if(NSHR.eq.3) THEN
          DDSDDT(3) = -one
        END IF
        
        RPL = -(DET - det0)/dtime
       
        DRPLDE = ZERO
        DRPLDE(1) = -det
        DRPLDE(2) = -det
        if(NSHR.eq.3) THEN
          DRPLDE(3) = -det
        END IF 
       
        statev(6) = det
        
      RETURN
      END      
c ======================================================================
      Subroutine matinv2(A,B)
        !! Performs a direct calculation of the inverse of a 2×2 matrix.
        Real*8 A(2,2)   !! Matrix
        Real*8 B(2,2)   !! Inverse matrix
        Real*8 detinv

        ! Calculate the inverse determinant of the matrix
        detinv = 1/(A(1,1)*A(2,2) - A(1,2)*A(2,1))

        ! Calculate the inverse of the matrix
        B(1,1) = +detinv * A(2,2)
        B(2,1) = -detinv * A(2,1)
        B(1,2) = -detinv * A(1,2)
        B(2,2) = +detinv * A(1,1)
      end Subroutine

c ======================================================================
      Subroutine matinv3(A,B)
        !! Performs a direct calculation of the inverse of a 3×3 matrix.
        Real*8 A(3,3)   !! Matrix
        Real*8 B(3,3)   !! Inverse matrix
        Real*8 detinv

        ! Calculate the inverse determinant of the matrix
        detinv = 1/(A(1,1)*A(2,2)*A(3,3) - A(1,1)*A(2,3)*A(3,2)
     1         - A(1,2)*A(2,1)*A(3,3) + A(1,2)*A(2,3)*A(3,1)
     2         + A(1,3)*A(2,1)*A(3,2) - A(1,3)*A(2,2)*A(3,1))

        ! Calculate the inverse of the matrix
        B(1,1) = +detinv * (A(2,2)*A(3,3) - A(2,3)*A(3,2))
        B(2,1) = -detinv * (A(2,1)*A(3,3) - A(2,3)*A(3,1))
        B(3,1) = +detinv * (A(2,1)*A(3,2) - A(2,2)*A(3,1))
        B(1,2) = -detinv * (A(1,2)*A(3,3) - A(1,3)*A(3,2))
        B(2,2) = +detinv * (A(1,1)*A(3,3) - A(1,3)*A(3,1))
        B(3,2) = -detinv * (A(1,1)*A(3,2) - A(1,2)*A(3,1))
        B(1,3) = +detinv * (A(1,2)*A(2,3) - A(1,3)*A(2,2))
        B(2,3) = -detinv * (A(1,1)*A(2,3) - A(1,3)*A(2,1))
        B(3,3) = +detinv * (A(1,1)*A(2,2) - A(1,2)*A(2,1))
      end subroutine
