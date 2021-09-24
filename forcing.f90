!=======================================================================

#include "config.h"

  MODULE mod_forcing

    USE mod_param,  ONLY :  gammaDot,vis,pi!,UinfS,waveNS,waveNS2
    USE mod_common, ONLY :  vol,radMin

!-----------------------------------------------------------------------

    IMPLICIT NONE

    PRIVATE  lubCoeff,binarySearch,crossProd
    PUBLIC   nInteractionDyn        &
            ,yInteractionDyn

!-----------------------------------------------------------------------

  CONTAINS

!=======================================================================

  SUBROUTINE nInteractionDyn(rad,vel,ome                             &
                            ,velInf,omeInf,eetInf,velInfLap          &
                            ,for,tor,sig)

    USE mod_param,  ONLY :  waveNS2

!-----------------------------------------------------------------------

    IMPLICIT NONE

    REAL,                 INTENT(IN)   ::  rad
    REAL, DIMENSION(3),   INTENT(IN)   ::  vel,ome                   &
                                          ,velInf,omeInf,velInfLap
    REAL, DIMENSION(3,3), INTENT(IN)   ::  eetInf
    REAL, DIMENSION(3),   INTENT(OUT)  ::  for,tor
    REAL, DIMENSION(3,3), INTENT(OUT)  ::  sig

    REAL  ::  alpha1,alpha2

!-----------------------------------------------------------------------

    alpha1 = -6*pi*vis*rad
    alpha2 = -8*pi*vis*rad*rad*rad

    for    = alpha1*(vel-(velInf+rad*rad*velInfLap/6.))

    tor    = alpha2*(ome-omeInf)

    !Avoiding to compute explicitly the laplacian of eetInf
    sig      = (1.0 - 0.1*rad*rad*waveNS2)*eetInf
    sig(1,2) =  sig(1,2) + 0.1*rad*rad*waveNS2*0.5*gammaDot
    sig(2,1) =  sig(1,2)
    sig      = -2.5*alpha2/(3.*vol) * sig

!-----------------------------------------------------------------------

    RETURN

!-----------------------------------------------------------------------

  END SUBROUTINE nInteractionDyn

!=======================================================================

  SUBROUTINE yInteractionDyn(n,l                     &
                            ,rad1,rad2               &
                            ,deltaX,pos1             &
                            ,vel1,ome1               &
                            ,vel2,ome2,dist          &
                            ,velInf1,omeInf1,eetInf1 &
                            ,velInf2,omeInf2,eetInf2 &
                            ,forIJ,forJI             &
                            ,torIJ,torJI             &
                            ,sigL,sigCN,sigCT        &
                            ,sigR,sigA)

    USE mod_param,  ONLY : dtime
#ifdef _LUB_
    USE mod_param,  ONLY : lubOuterDist
#endif
#ifdef _CONTACT_
    USE mod_param,  ONLY : kN,gammaN
#endif
#ifdef _ESR_
    USE mod_param,  ONLY : fER,kappaDb,lambdaDb 
#endif
#ifdef _VDW_
    USE mod_param,  ONLY : hamA,eps
#endif

#if defined _FRICTION_ || defined _ROLLING_
    USE mod_common, ONLY : nAroundMax,idxCtcArr       &
                          ,ctcArr,ctcArrNew
#ifdef _FRICTION_
    USE mod_param,  ONLY : kT,muC
    USE mod_common, ONLY : csiX,csiY,csiZ             &
                          ,csiXNew,csiYNew,csiZNew
#endif
#ifdef _ROLLING_
    USE mod_param,  ONLY : kR,muR
    USE mod_common, ONLY : csiRX,csiRY,csiRZ          &
                          ,csiRXNew,csiRYNew,csiRZNew
#endif
#endif

!-----------------------------------------------------------------------

    IMPLICIT NONE

    INTEGER,              INTENT(IN)   ::  n,l
    REAL,                 INTENT(IN)   ::  rad1,rad2,dist
    REAL, DIMENSION(3),   INTENT(IN)   ::  deltaX,pos1             &
                                          ,vel1,ome1               &
                                          ,vel2,ome2               &
                                          ,velInf1,omeInf1         &
                                          ,velInf2,omeInf2
    REAL, DIMENSION(3,3), INTENT(IN)   ::  eetInf1,eetInf2
    REAL, DIMENSION(3),   INTENT(OUT)  ::  forIJ,forJI,torIJ,torJI
    REAL, DIMENSION(3,3), INTENT(OUT)  ::  sigL,sigCN,sigCT        &
                                          ,sigR,sigA


    INTEGER               ::  i,j,idx 
    REAL                  ::  invDist,h12,radM
    REAL, DIMENSION(3)    ::  pos2,nor1,aux,aux2

#ifdef _LUB_
    REAL                 ::  x11a,x12a,x22a       &
                            ,y11a,y12a,y22a       &
                            ,y11b,y12b,y21b,y22b  &
                            ,y11c,y12c,y22c       &
                            ,x11g,x12g,x21g,x22g  &
                            ,y11g,y12g,y21g,y22g  &
                            ,y11h,y12h,y21h,y22h  &
                            ,pE1,pE2,tE1,tE2
                            !,uInf1,uInf2          &
                            !,omeIn1,omeIn2        &
                            !,eIn121,eIn122
    REAL, DIMENSION(3)   ::  forL,forL2
#endif

#ifdef _CONTACT_
    REAL, DIMENSION(3)   ::  forC
#if defined _FRICTION_ || defined _ROLLING_
    REAL                 ::  csi1N,csi1Max      &
                            ,rad1c,rad2c
    REAL, DIMENSION(3)   ::  csi1
#ifdef _FRICTION_
    REAL, DIMENSION(3)   ::  csiOld
#endif
#ifdef _ROLLING_
    REAL                 ::  radMc
    REAL, DIMENSION(3)   ::  csiROld
#endif
#endif
#endif

#ifdef _ESR_
    REAL, DIMENSION(3)   ::  forE
#endif
#ifdef _VDW_
    REAL, DIMENSION(3)   ::  forA
#endif

!----------------------------------------------------------------------

    forIJ = 0.
    forJI = 0.
    torIJ = 0.
    torJI = 0.
    sigL  = 0.
    sigCN = 0.
    sigCT = 0.
    sigR = 0.
    sigA = 0.

    invDist = 1/dist

    nor1 = deltaX*invDist

    pos2 = pos1 + deltaX

    h12 = dist - (rad1+rad2)

    radM = 2*rad1*rad2/(rad1+rad2)

!----------------------------------------------------------------------

#ifdef _LUB_
    !Lubrication forcing
    IF (h12.LE.lubOuterDist*radMin .AND. h12.GT.0.) THEN
      CALL lubCoeff(rad1,rad2,h12          &
                   ,x11a,x12a,x22a         &
                   ,y11a,y12a,y22a         &
                   ,y11b,y12b,y21b,y22b    &
                   ,y11c,y12c,y22c         &
                   ,x11g,x12g,x21g,x22g    &
                   ,y11g,y12g,y21g,y22g    &
                   ,y11h,y12h,y21h,y22h)

      !uInf1  =  gammaDot*pos1(2) + UinfS*SIN(waveNS*pos1(2))
      !uInf2  =  gammaDot*pos2(2) + UinfS*SIN(waveNS*pos2(2))
      !omeIn1 = -0.5*(gammaDot+UinfS*waveNS*COS(waveNS*pos1(2)))
      !omeIn2 = -0.5*(gammaDot+UinfS*waveNS*COS(waveNS*pos2(2)))
      !eIn121 =  0.5*(gammaDot+UinfS*waveNS*COS(waveNS*pos1(2)))
      !eIn122 =  0.5*(gammaDot+UinfS*waveNS*COS(waveNS*pos2(2)))

      aux    = y11b*(ome1-omeInf1) + y21b*(ome2-omeInf2)
      forL   = crossProd(aux,nor1)
      aux    = y12b*(ome1-omeInf1) + y22b*(ome2-omeInf2)
      forL2  = crossProd(aux,nor1)

      aux    = vel1-velInf1
      aux2   = DOT_PRODUCT(nor1,aux)*nor1
      forL   = forL  - (x11a*aux2 + y11a*(aux - aux2))
      forL2  = forL2 - (x12a*aux2 + y12a*(aux - aux2))
      torIJ  = torIJ + crossProd(y11b*aux,nor1)
      torJI  = torJI + crossProd(y21b*aux,nor1)

      aux    = vel2-velInf2
      aux2   = DOT_PRODUCT(nor1,aux)*nor1
      forL   = forL  - (x12a*aux2 + y12a*(aux - aux2))
      forL2  = forL2 - (x22a*aux2 + y22a*(aux - aux2))
      torIJ  = torIJ + crossProd(y12b*aux,nor1)
      torJI  = torJI + crossProd(y22b*aux,nor1)

      pE1 = 0.
      pE2 = 0.
      tE1 = 0.
      tE2 = 0.
      DO i=1,3
        pE1 = pE1 + DOT_PRODUCT(eetInf1(:,i),nor1)*nor1(i)
        pE2 = pE2 + DOT_PRODUCT(eetInf2(:,i),nor1)*nor1(i)
        tE1 = tE1 + eetInf1(i,i)
        tE2 = tE2 + eetInf2(i,i)
      ENDDO

      aux(1) = 2*DOT_PRODUCT(eetInf1(:,1),nor1) !the tensor is symmetric
      aux(2) = 2*DOT_PRODUCT(eetInf1(:,2),nor1) !the tensor is symmetric
      aux(3) = 2*DOT_PRODUCT(eetInf1(:,3),nor1) !the tensor is symmetric
      forL   = forL  + (x11g*(pE1-tE1/3.)-2*y11g*pE1)*nor1   &
             + y11g*aux
      forL2  = forL2 + (x12g*(pE1-tE1/3.)-2*y12g*pE1)*nor1   &
             + y12g*aux
      torIJ  = torIJ + 2*y11h*crossProd(nor1,aux)
      torJI  = torJI + 2*y12h*crossProd(nor1,aux)

      aux(1) = 2*DOT_PRODUCT(eetInf2(:,1),nor1) !the tensor is symmetric
      aux(2) = 2*DOT_PRODUCT(eetInf2(:,2),nor1) !the tensor is symmetric
      aux(3) = 2*DOT_PRODUCT(eetInf2(:,3),nor1) !the tensor is symmetric
      forL   = forL  + (x21g*(pE2-tE2/3.)-2*y21g*pE2)*nor1   &
             + y21g*aux
      forL2  = forL2 + (x22g*(pE2-tE2/3.)-2*y22g*pE2)*nor1   &
             + y22g*aux
      torIJ  = torIJ + 2*y21h*crossProd(nor1,aux)
      torJI  = torJI + 2*y22h*crossProd(nor1,aux)

      forIJ  = forIJ + forL
      forJI  = forJI + forL2
      forL   = forL  - forL2
      DO j = 1,3
        DO i = 1,3
          sigL(i,j) = (forL(i)*deltaX(j) +        &
                       forL(j)*deltaX(i))*0.5/vol
        ENDDO
      ENDDO

      aux    = ome1-omeInf1
      aux    = aux - DOT_PRODUCT(nor1,aux)*nor1
      torIJ  = torIJ - y11c*aux
      torJI  = torJI - y12c*aux

      aux    = ome2-omeInf2
      aux    = aux - DOT_PRODUCT(nor1,aux)*nor1
      torIJ  = torIJ - y12c*aux
      torJI  = torJI - y22c*aux

    ENDIF
#endif

!----------------------------------------------------------------------

#ifdef _CONTACT_
    !Contact forcing
    IF (h12.LE.0.) THEN !On Ge's, the sign of h12 is the opposite
      forC = kN*h12*nor1

      aux = vel2 - vel1 !Sign useful for DSYMV/MATMUL
      forC = forC + gammaN*DOT_PRODUCT(nor1,aux)*nor1
      DO j = 1,3
        DO i = 1,3
          sigCN(i,j) = 2*forC(i)*deltaX(j)/vol
        ENDDO
      ENDDO

#if defined _FRICTION_ || defined _ROLLING_
      rad1c = rad1 + h12*0.5 !h12 is negative
      rad2c = rad2 + h12*0.5

      csi1Max = SQRT(forC(1)*forC(1) + forC(2)*forC(2) + forC(3)*forC(3))

      CALL binarySearch(nAroundMax,ctcArr(:,n),l,idx)
      idxCtcArr(n) = idxCtcArr(n) + 1

      IF (idx .EQ. -1) THEN
#ifdef _FRICTION_
        csiOld(1) = 0.
        csiOld(2) = 0.
        csiOld(3) = 0.
#endif
#ifdef _ROLLING_
        csiROld(1) = 0.
        csiROld(2) = 0.
        csiROld(3) = 0.
#endif
      ELSE
#ifdef _FRICTION_
        csiOld(1) = csiX(idx,n)
        csiOld(2) = csiY(idx,n)
        csiOld(3) = csiZ(idx,n)
        csiOld    = csiOld - DOT_PRODUCT(csiOld,nor1)*nor1
#endif
#ifdef _ROLLING_
        csiROld(1) = csiRX(idx,n)
        csiROld(2) = csiRY(idx,n)
        csiROld(3) = csiRZ(idx,n)
        csiROld    = csiROld - DOT_PRODUCT(csiROld,nor1)*nor1
#endif
      ENDIF

#ifdef _FRICTION_
      aux  = rad1c*ome1 + rad2c*ome2
      csi1 = crossProd(aux,nor1)
      aux  = vel2 - vel1 - csi1 !Sign useful for DSYMV
      csi1 = aux - DOT_PRODUCT(aux,nor1)*nor1

      csi1(1) = csiOld(1) + dtime*csi1(1)
      csi1(2) = csiOld(2) + dtime*csi1(2)
      csi1(3) = csiOld(3) + dtime*csi1(3)

      csi1N   = SQRT(csi1(1)*csi1(1) + csi1(2)*csi1(2) + csi1(3)*csi1(3))

      IF (csi1N.GE.csi1Max*muC/kT) csi1 = csi1/csi1N*csi1Max*muC/kT
      csiXNew(idxCtcArr(n),n) = csi1(1)
      csiYNew(idxCtcArr(n),n) = csi1(2)
      csiZNew(idxCtcArr(n),n) = csi1(3)

      forC  = forC  + kT*csi1
      DO j = 1,3
        DO i = 1,3
          sigCT(i,j) = 2*kT*csi1(i)*deltaX(j)/vol
        ENDDO
      ENDDO
      torIJ = torIJ + crossProd(rad1c*kT*nor1,csi1)
      torJI = torJI + crossProd(rad2c*kT*nor1,csi1)
#endif
#ifdef _ROLLING_
      radMc = 2*rad1c*rad2c/(rad1c+rad2c)
      aux   = ome1 - ome2
      csi1  = crossProd(aux,nor1)
      aux   = radMc*csi1 !Sign useful for DSYMV
      csi1  = aux - DOT_PRODUCT(aux,nor1)*nor1

      csi1(1) = csiROld(1) + dtime*csi1(1)
      csi1(2) = csiROld(2) + dtime*csi1(2)
      csi1(3) = csiROld(3) + dtime*csi1(3)

      csi1N   = SQRT(csi1(1)*csi1(1) + csi1(2)*csi1(2) + csi1(3)*csi1(3))

      IF (csi1N.GE.csi1Max*muR/kR) csi1 = csi1/csi1N*csi1Max*muR/kR
      csiRXNew(idxCtcArr(n),n) = csi1(1)
      csiRYNew(idxCtcArr(n),n) = csi1(2)
      csiRZNew(idxCtcArr(n),n) = csi1(3)

      torIJ = torIJ - crossProd(radMc*kR*nor1,csi1)
      torJI = torJI + crossProd(radMc*kR*nor1,csi1)
#endif
      ctcArrNew(idxCtcArr(n),n) = l
#endif

      forIJ = forIJ + forC
      forJI = forJI - forC

    ENDIF

#endif

!----------------------------------------------------------------------

#ifdef _ESR_
    !Electrostatic repulsion
    IF (h12.LE.20*lambdaDb*radMin) THEN
      IF (h12.GE.0) THEN
        forE = -fER*radM/radMin*EXP(-kappaDb/radMin*h12)*nor1
      ELSE
        forE = -fER*radM/radMin*nor1
      ENDIF
      forIJ = forIJ + forE
      forJI = forJI - forE
      DO j = 1,3
        DO i = 1,3
          sigR(i,j) = 2*forE(i)*deltaX(j)/vol
        ENDDO
      ENDDO
    ENDIF
#endif

!----------------------------------------------------------------------

#ifdef _VDW_
    !Van der Waals attraction
    IF (dist.LE.4*radMin) THEN
      IF (h12.GE.0.) THEN
        forA = hamA*radM/(12.*(h12**2+(eps*radM)**2))*nor1
      ELSE
        forA = hamA*radM/(12.*(eps*radM)**2)*nor1
      ENDIF
      forIJ = forIJ + forA
      forJI = forJI - forA
      DO j = 1,3
        DO i = 1,3
          sigA(i,j) = 2*forA(i)*deltaX(j)/vol
        ENDDO
      ENDDO
    ENDIF
#endif

!----------------------------------------------------------------------

    RETURN

!----------------------------------------------------------------------

  END SUBROUTINE yInteractionDyn

!=======================================================================

  SUBROUTINE lubCoeff(rad1,rad2,h12          &
                     ,x11a,x12a,x22a         &
                     ,y11a,y12a,y22a         &
                     ,y11b,y12b,y21b,y22b    &
                     ,y11c,y12c,y22c         &
                     ,x11g,x12g,x21g,x22g    &
                     ,y11g,y12g,y21g,y22g    &
                     ,y11h,y12h,y21h,y22h)

!----------------------------------------------------------------------

    IMPLICIT NONE

    REAL, INTENT(IN)   ::  rad1,rad2,h12
    REAL, INTENT(OUT)  ::  x11a,x12a,x22a         &
                          ,y11a,y12a,y22a         &
                          ,y11b,y12b,y21b,y22b    &
                          ,y11c,y12c,y22c         &
                          ,x11g,x12g,x21g,x22g    &
                          ,y11g,y12g,y21g,y22g    &
                          ,y11h,y12h,y21h,y22h


    REAL  ::  alpha,lambda,lambda2           &
             ,lambdaM1,lambdaM12             &
             ,inv1Plambda,inv1Plambda2       &
             ,inv1Plambda3,inv1Plambda5      &
             ,inv1PlambdaM1,inv1PlambdaM12   &
             ,inv1PlambdaM13,inv1PlambdaM15  &
             ,delta,invDelta                 &
             ,rad12,rad22,rad13,rad23        &
             ,r1r2,r1r22,r1r23

!----------------------------------------------------------------------

    alpha          =  6*pi*vis
    lambda         =  rad2/rad1
    lambda2        =  lambda*lambda
    lambdaM1       =  1/lambda
    lambdaM12      =  lambdaM1*lambdaM1
    inv1Plambda    =  1/(1+lambda)
    inv1Plambda2   =  inv1Plambda*inv1Plambda
    inv1plambda3   =  inv1plambda2*inv1Plambda
    inv1plambda5   =  inv1plambda3*inv1Plambda2
    inv1PlambdaM1  =  1/(1+lambdaM1)
    inv1PlambdaM12 =  inv1PlambdaM1*inv1PlambdaM1
    inv1PlambdaM13 =  inv1PlambdaM12*inv1PlambdaM1
    inv1PlambdaM15 =  inv1PlambdaM13*inv1PlambdaM12
    delta          =  2*(h12+1e-3*radMin)/(rad1+rad2)
    invDelta       =  1/delta

!----------------------------------------------------------------------

    rad12 = rad1*rad1
    rad22 = rad2*rad2
    rad13 = rad12*rad1
    rad23 = rad22*rad2
    r1r2  = rad1+rad2
    r1r22 = r1r2*r1r2
    r1r23 = r1r22*r1r2

!----------------------------------------------------------------------

    x11a =  2*rad1*lambda2*inv1Plambda3*invDelta*alpha
    x12a = -x11a/rad1*r1r2*inv1Plambda
    x22a =  2*rad2*lambdaM12*inv1PlambdaM13*invDelta*alpha

    y11a =  4*rad1*lambda*(2+lambda+2*lambda2)*inv1Plambda3/15.*LOG(invDelta)*alpha
    y12a = -y11a/rad1*r1r2*inv1Plambda
    y22a =  4*rad2*lambdaM1*(2+lambdaM1+2*lambdaM12)*inv1PlambdaM13/15.*LOG(invDelta)*alpha

    y11b = -2*rad12*lambda*(4+lambda)*inv1Plambda2/15.*LOG(invDelta)*alpha
    y12b = -y11b/rad12*r1r22*inv1Plambda2
    y22b =  2*rad22*lambdaM1*(4+lambdaM1)*inv1PlambdaM12/15.*LOG(invDelta)*alpha
    y21b = -y22b/rad22*r1r22*inv1PlambdaM12

    y11c =  8*rad13*lambda*inv1Plambda/15.*LOG(invDelta)*alpha
    y12c =  2*r1r23*lambda2*inv1Plambda2*inv1Plambda2/15.*LOG(invDelta)*alpha
    y22c =  lambda2*y11c

    x11g =  2*rad12*lambda2*inv1Plambda3*invDelta*alpha
    x12g = -x11g/rad12*r1r22*inv1Plambda2
    x22g = -2*rad22*lambdaM12*inv1PlambdaM13*invDelta*alpha
    x21g = -x22g/rad22*r1r22*inv1PlambdaM12

    y11g =  rad12*lambda*(4-lambda+7*lambda2)*inv1Plambda3/15.*LOG(invDelta)*alpha
    y12g = -y11g/rad12*r1r22*inv1Plambda2
    y22g = -rad22*lambdaM1*(4-lambdaM1+7*lambdaM12)*inv1PlambdaM13/15.*LOG(invDelta)*alpha
    y21g = -y22g/rad22*r1r22*inv1PlambdaM12

    y11h =  2*rad13*lambda*(2-lambda)*inv1Plambda2/15.*LOG(invDelta)*alpha
    y12h =  r1r23*lambda2*(1+7*lambda)*inv1Plambda5/15.*LOG(invDelta)*alpha
    y22h =  2*rad23*lambdaM1*(2-lambdaM1)*inv1PlambdaM12/15.*LOG(invDelta)*alpha
    y21h =  r1r23*lambdaM12*(1+7*lambdaM1)*inv1PlambdaM15/15.*LOG(invDelta)*alpha

!----------------------------------------------------------------------

    RETURN

!----------------------------------------------------------------------

  END SUBROUTINE lubCoeff

!=======================================================================

  SUBROUTINE binarySearch(n,arr,Tar,iT)

    !Binary search for an array sorted in a descending order

!----------------------------------------------------------------------

    IMPLICIT NONE

    INTEGER,               INTENT(IN)  :: n,Tar
    INTEGER, DIMENSION(:), INTENT(IN)  :: arr
    INTEGER,               INTENT(OUT) :: iT

    INTEGER :: l,r,idx

!----------------------------------------------------------------------

    l = 1
    r = n
    iT = -1

    DO WHILE (l .LE. r)
        idx = (l + r)/2
        IF (arr(idx) .GT. Tar) THEN
            r = idx - 1
        ELSEIF (arr(idx) .EQ. Tar) THEN
            iT = idx
            l = r + 1
        ELSE
            l = idx + 1
        ENDIF
    ENDDO

!----------------------------------------------------------------------
            
    RETURN

!----------------------------------------------------------------------

  END SUBROUTINE binarySearch

!=======================================================================

  FUNCTION crossProd(v1,v2) RESULT(v3)

    !Cross product v3 = v1 x v2

!----------------------------------------------------------------------

    IMPLICIT NONE

    REAL, DIMENSION(:)  ::  v1,v2
    REAL, DIMENSION(3)  ::  v3

!----------------------------------------------------------------------

    v3(1) = v1(2)*v2(3)-v1(3)*v2(2)
    v3(2) = v1(3)*v2(1)-v1(1)*v2(3)
    v3(3) = v1(1)*v2(2)-v1(2)*v2(1)

!----------------------------------------------------------------------

    RETURN

!----------------------------------------------------------------------

  END FUNCTION crossProd

!=======================================================================

  END MODULE mod_forcing

!=======================================================================
