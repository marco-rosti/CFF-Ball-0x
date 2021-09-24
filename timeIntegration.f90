!=======================================================================

#include "config.h"

  MODULE mod_timeInt

    USE mod_param,   ONLY :  nPar,dtime,tmax        &
                            ,lx,ly,lz               &
                            ,nxBin,nyBin,nzBin      &
                            ,xBin,yBin,zBin         &
                            ,rSearch                &
                            ,gammaDot,vis           &
                            ,itSave,itField         &
                            ,itRunTime,itPrint      &
                            ,itStat,tStat           &
                            ,nthreads

     USE mod_common, ONLY :  posArr,velArr,accArr   &
                            ,omeArr,omedArr         &
                            ,forArr,torArr          &
                            ,rArr,masArr,ineArr     &
                            ,sigmaT,sigmaL          &
                            ,sigmaCN,sigmaCT        &
                            ,sigmaR,sigmaA          &
                            ,ntStat                 &
                            ,sigmaT_,sigmaL_        &
                            ,sigmaCN_,sigmaCT_      &
                            ,sigmaR_,sigmaA_

#if defined _FRICTION_ || defined _ROLLING_
    USE mod_common, ONLY :   nAroundMax,idxCtcArr   &
                            ,ctcArr,ctcArrNew
#ifdef _FRICTION_
    USE mod_common, ONLY :   csiX,csiY,csiZ         &
                            ,csiXNew,csiYNew,csiZNew
#endif
#ifdef _ROLLING_
    USE mod_common, ONLY :   csiRX,csiRY,csiRZ      &
                            ,csiRXNew,csiRYNew,csiRZNew
#endif
#endif

    USE mod_forcing, ONLY :  nInteractionDyn        &
                            ,yInteractionDyn

    USE mod_stats,  ONLY :   updateStats

    IMPLICIT NONE

    PRIVATE  writeRestart,writeField                &
            ,compInfVal                             &
#ifndef _FLOW_SHEAR_
            ,compInfVelLap                          &
#endif
            ,heapsort,siftdown
    PUBLIC   timeIntegration

  CONTAINS

!=======================================================================

  SUBROUTINE timeIntegration

    USE omp_lib
    USE mod_common, ONLY :  it,time

    IMPLICIT NONE

    INTEGER  ::  ntmax              &
                ,i,j,k,n            &
                ,ii,jj,kk           &
                ,l,lS,lE            &
                ,nBins              &
                ,iB,jB,kB,iB2       &
                ,idx,idx2           &
                ,idx2isN,idx2Loc    &
                ,alpha,beta,ithread &
                ,neighBinMax        &
                ,sumAuxPS,offsetPS  &
                ,iaux

    INTEGER, DIMENSION(:), ALLOCATABLE  ::  iBin,nArr     &
                                           ,cnt,aux       &
                                           ,sumPS,iBS,iBE
    INTEGER, DIMENSION(:,:), ALLOCATABLE  ::  cnt2

    REAL                  ::  dist,auxPos1
    REAL, DIMENSION(3)    ::  auxDX,auxV2,forIJ,torIJ,forJI,torJI  &
                             ,velInf1,omeInf1,velInfLap1           &
                             ,velInf2,omeInf2
    REAL, DIMENSION(3,3)  ::  sigL,sigCN,sigCT,sigR,sigA           &
                             ,eetInf1,eetInf2

    REAL, DIMENSION(:,:),   ALLOCATABLE  ::  posArrOld             &
                                            ,velArrOld,accArrOld   &
                                            ,omeArrOld,omedArrOld
    REAL, DIMENSION(:),     ALLOCATABLE  ::  invMass,invIner  

    REAL, DIMENSION(:,:,:), ALLOCATABLE  ::  forAux,torAux

    REAL  ::  dtime2,start,finish,cpuTimeCumul,cpuTimeIter,tWall

#if defined _FRICTION_ || defined _ROLLING_
    INTEGER, DIMENSION(:,:), ALLOCATABLE  ::  ctcArrIdx
    INTEGER, DIMENSION(:),   ALLOCATABLE  ::  idxTest
#endif

    LOGICAL  ::  flagRes = .TRUE.

    ALLOCATE(posArrOld(3,nPar),velArrOld(3,nPar),accArrOld(3,nPar)  &
            ,omeArrOld(3,nPar),omedArrOld(3,nPar)                   &
            ,invMass(nPar),invIner(nPar)                            &
            ,forAux(3,nthreads,nPar),torAux(3,nthreads,nPar))
#if defined _FRICTION_ || defined _ROLLING_
    ALLOCATE(ctcArrIdx(nAroundMax,nPar),idxTest(nAroundMax))
    DO i = 1,nAroundMax
      idxTest(i) = i
    ENDDO
#endif

    invMass = 1./masArr
    invIner = 1./ineArr

    dtime2  = dtime*dtime
    ntmax   = NINT(tmax/dtime)

    tWall = 1.90

    nBins = nxBin*nyBin*nzBin
    ALLOCATE(iBin(nPar),nArr(nPar),cnt(nBins),aux(nPar)     &
            ,sumPS(nthreads+1),iBS(nthreads),iBE(nthreads),cnt2(nthreads,nBins))

    idx2isN = 5
    neighBinMax = 3

    jj = nBins/nthreads
    ii = MOD(nBins,nthreads)

    j=0
    DO i=1,nthreads
      iBS(i) = (i-1)*jj+1+j
      IF (i.LE.ii) THEN
          j=j+1
      ENDIF
    ENDDO
    DO i=1,nthreads-1
      iBE(i) = iBS(i+1)-1
    ENDDO
    iBE(nthreads)=nBins

    WRITE(*,9999)

    OPEN(1,FILE='checkConvergence.dat',STATUS='unknown' &
          ,ACTION='write',POSITION='append')
 
    CALL OMP_SET_DYNAMIC(.FALSE.)
    CALL OMP_SET_NUM_THREADS(nthreads)

!$OMP PARALLEL PRIVATE(sigmaT,start,finish,ithread)        &
!$OMP& PRIVATE(sumAuxPS,offsetPS,cpuTimeIter,cpuTimeCumul) &
!$OMP& PRIVATE(velInf1,omeInf1,eetInf1,velInfLap1)         &
!$OMP& PRIVATE(velInf2,omeInf2,eetInf2)                    &
!$OMP& PRIVATE(iaux) FIRSTPRIVATE(it,time)

    ithread = OMP_GET_THREAD_NUM()

    start = OMP_GET_WTIME()

    cpuTimeCumul = 0.

    sumPS(ithread+1) = 0
!$OMP SINGLE
    sumPS(nthreads+1) = 0
!$OMP END SINGLE NOWAIT
 
!$OMP DO SCHEDULE(static)
    DO i=1,nBins
      cnt(i) = 0
      cnt2(:,i) = 0
    ENDDO
!$OMP END DO

    velInf1 = 0.
    omeInf1 = 0.
    eetInf1 = 0.
    velInf2 = 0.
    omeInf2 = 0.
    eetInf2 = 0.

    velInfLap1 = 0.
    
!$OMP BARRIER
    DO WHILE (it.LT.ntmax .AND. time.LE.tmax)

      it = it + 1
      time = time + dtime

      sumAuxPS = 0
      offsetPS = 0

!$OMP DO SCHEDULE(static)
      DO i=1,nPar
        forAux(:,:,i) = 0.
        torAux(:,:,i) = 0.

        posArrOld(:,i)  = posArr(:,i)
        velArrOld(:,i)  = velArr(:,i)
        accArrOld(:,i)  = accArr(:,i)
        omeArrOld(:,i)  = omeArr(:,i)
        omedArrOld(:,i) = omedArr(:,i)
       
        posArr(:,i) = posArrOld(:,i) + dtime*velArrOld(:,i) + 0.5*dtime2*accArrOld(:,i)
        velArr(:,i) = velArrOld(:,i) + 0.5*dtime*accArrOld(:,i)
        omeArr(:,i) = omeArrOld(:,i) + 0.5*dtime*omedArrOld(:,i)
       
        posArrOld(:,i)  = posArr(:,i)

        IF (posArr(2,i).LT.0.) THEN 
          posArr(1,i) = posArr(1,i) + gammaDot*ly*time
          posArr(1,i) = posArr(1,i) - INT(posArr(1,i)/lx)*lx
          velArr(1,i) = velArr(1,i) + gammaDot*ly
          posArr(2,i) = posArr(2,i) + ly
        ENDIF
        IF (posArr(2,i).GE.ly) THEN 
          posArr(1,i) = posArr(1,i) - gammaDot*ly*time
          posArr(1,i) = (1-INT(posArr(1,i)/lx))*lx + posArr(1,i)
          velArr(1,i) = velArr(1,i) - gammaDot*ly
          posArr(2,i) = posArr(2,i) - ly
        ENDIF
        IF (posArr(1,i).LT.0.) THEN 
          posArr(1,i) = posArr(1,i) + lx
        ENDIF
        IF (posArr(1,i).GE.lx) THEN 
          posArr(1,i) = posArr(1,i) - lx
        ENDIF
        IF (posArr(3,i).LT.0.) THEN 
          posArr(3,i) = posArr(3,i) + lz
        ENDIF
        IF (posArr(3,i).GE.lz) THEN 
          posArr(3,i) = posArr(3,i) - lz
        ENDIF

        !Friction and/or Rolling init
#if defined _FRICTION_ || defined _ROLLING_
        idxCtcArr(i) = 0
#endif
      ENDDO
!$OMP END DO NOWAIT

!-----------------------------------------------------------------------

!$OMP DO SCHEDULE(static) PRIVATE(i,j,k,idx)
      DO n=1,nPar

        forArr(:,n) = 0.
        torArr(:,n) = 0.
 
        i = INT(posArr(1,n)/xBin)+1
        IF (i.GT.nxBin) i = i-nxBin
        j = INT(posArr(2,n)/yBin)+1
        IF (j.GT.nyBin) j = j-nyBin
        k = INT(posArr(3,n)/zBin)+1
        IF (k.GT.nzBin) k = k-nzBin
 
        idx = i+(j-1)*nxBin+(k-1)*nxBin*nyBin
        aux(n) = idx
        cnt2(ithread+1,idx) = cnt2(ithread+1,idx)+1
 
      ENDDO
!$OMP ENDDO
!$OMP DO SCHEDULE(static) PRIVATE(j)
      DO n=1,nBins
        DO j=1,nthreads
          cnt(n) = cnt(n)+cnt2(j,n)
        ENDDO
      ENDDO
!$OMP ENDDO

!$OMP DO SCHEDULE(STATIC)
      DO i=1,nBins
        sumAuxPS = sumAuxPS + cnt(i)
        cnt(i) = sumAuxPS
        sumPS(ithread+2) = sumAuxPS
      ENDDO
!$OMP END DO
      DO i=1,ithread+1
        offsetPS = offsetPS+sumPS(i)
      ENDDO
!$OMP DO SCHEDULE(STATIC)
      DO i = 1,nthreads
        cnt(iBS(i):iBE(i))=cnt(iBS(i):iBE(i))+offsetPS
      ENDDO
!$OMP END DO
          
      DO iaux=1,nPar
        IF (aux(iaux).GE.iBS(ithread+1) .AND.  &
            aux(iaux).LE.iBE(ithread+1)) THEN
          iBin(cnt(aux(iaux))) = aux(iaux)
          nArr(cnt(aux(iaux))) = iaux
          cnt(aux(iaux)) = cnt(aux(iaux))-1
        ENDIF
      ENDDO

      sigmaT      = 0.
      sigmaT(1,2) = 2*vis*0.5*gammaDot
      sigmaT(2,1) = 2*vis*0.5*gammaDot
!$OMP SINGLE
      !Stress tensor init --> Sigma = 2*vis*E^{\infty}
      sigmaL      = 0.
      sigmaCN     = 0.
      sigmaCT     = 0.
      sigmaR      = 0.
      sigmaA      = 0.
!$OMP END SINGLE
! Implicit barrier at the end of single

! !$OMP BARRIER

!$OMP DO SCHEDULE(static)                      &
!$OMP& PRIVATE(i,j,k,ii,jj,kk,iB,iB2,jB,kB)    &
!$OMP& PRIVATE(l,lS,lE,idx,idx2,idx2Loc)       &
!$OMP& PRIVATE(alpha,beta)                     &
!$OMP& PRIVATE(auxPos1,auxV2,auxDX,dist)       &
!$OMP& PRIVATE(forIJ,torIJ,forJI,torJI)        &
!$OMP& PRIVATE(sigL,sigCN,sigCT,sigR,sigA)     & 
!$OMP& REDUCTION(+:sigmaL,sigmaCN,sigmaCT)     &
!$OMP& REDUCTION(+:sigmaR,sigmaA)

      DO n=1,nPar
        idx = iBin(n)
        alpha = 0
        beta  = 0
        iB = idx - idx/nxBin*nxBin !INTEGERS DIVISION
        IF (iB.EQ.0) THEN
          iB = nxBin
          alpha = 1
        ENDIF
    
        jB = idx/nxBin + 1 - idx/(nxBin*nyBin)*nyBin - alpha !INTEGERS DIVISION
        IF (jB.EQ.0) THEN
          jB = nyBin
          beta = 1
        ENDIF
    
        kB = idx/(nxBin*nyBin) + 1 - beta !INTEGERS DIVISION

        CALL compInfVal(posArr(:,nArr(n)),velInf1,omeInf1,eetInf1)
#ifndef _FLOW_SHEAR_
        CALL compInfVelLap(posArr(:,nArr(n)),velInfLap1)
#endif
    
        CALL nInteractionDyn(rArr(nArr(n))                          &
                            ,velArr(:,nArr(n)),omeArr(:,nArr(n))    &
                            ,velInf1,omeInf1,eetInf1,velInfLap1     &
                            ,forIJ,torIJ,sigL)
        !NOTE: it should be actually forI and torI as, here, there 
        !is no interaction. Chosen forIJ and torIJ for memory reasons
        forArr(:,nArr(n)) = forArr(:,nArr(n)) + forIJ
        torArr(:,nArr(n)) = torArr(:,nArr(n)) + torIJ
        sigmaL = sigmaL + sigL

        DO k = kB,kB+1
          kk = k
          IF (k.GT.nzBin) kk = kk - nzBin
    
          DO j = jB-1,jB+1
            jj  = j
            iB2 = iB
            IF (j.LT.1) THEN
              jj  = jj + nyBin
              auxPos1 = posArr(1,nArr(n)) + gammaDot*ly*time
              auxPos1 = auxPos1 - INT(auxPos1/lx)*lx
              iB2 = INT(auxPos1/xBin)+1 !calcolo a che iBin sono (mettiamo iBin=10, con nBx=4)
            ENDIF
            IF (j.GT.nyBin) THEN
              jj  = jj - nyBin
              auxPos1 = posArr(1,nArr(n)) - gammaDot*ly*time
              auxPos1 = (1-INT(auxPos1/lx))*lx + auxPos1
              iB2 = INT(auxPos1/xBin)+1
            ENDIF
            iB2 = MOD(iB2,nxBin)        !mi riporto al iBin corretto (MOD(iBin,nBx)=2)
            IF (iB2.EQ.0) iB2 = nxBin
    
            DO i = iB2-1,iB2+1
              ii = i
              IF (i.LT.1)     ii = ii + nxBin
              IF (i.GT.nxBin) ii = ii - nxBin
    
              idx2 = ii + (jj-1)*nxBin + (kk-1)*nxBin*nyBin
              idx2Loc = i-(iB2-2) + ((j-(jB-2)-1)  & 
                                  +  (k-(kB-1)-1)*neighBinMax)*neighBinMax
              IF (idx2Loc.GE.idx2isN) THEN
    
              lS = cnt(idx2)+1
              IF (idx2.LT.nBins) THEN
                lE = cnt(idx2+1)
              ELSE
                lE = nPar
              ENDIF
             
              DO l = lS,lE
                IF ( (nArr(l).EQ.nArr(n)) .OR.                        &
                     (idx2Loc.EQ.idx2isN .AND. nArr(l).LT.nArr(n)) ) CYCLE
                auxPos1  = posArr(1,nArr(l))
                auxV2    = velArr(:,nArr(l))
                auxDX(2) = posArr(2,nArr(l))-posArr(2,nArr(n))
                IF (auxDX(2).GE. 2*yBin) THEN
                  auxDX(2) = auxDX(2)-ly
                  auxV2(1) = velArr(1,nArr(l))-gammaDot*ly
                  auxPos1 = auxPos1 - gammaDot*ly*time
                  auxPos1 = auxPos1 + (1-INT(auxPos1/lx))*lx
                  IF (auxPos1.GE.lx) auxPos1 = auxPos1 - lx
                ENDIF
                IF (auxDX(2).LE.-2*yBin) THEN
                  auxDX(2) = auxDX(2)+ly
                  auxV2(1) = velArr(1,nArr(l))+gammaDot*ly
                  auxPos1 = auxPos1 + gammaDot*ly*time
                  auxPos1 = auxPos1 - INT(auxPos1/lx)*lx
                ENDIF
                auxDX(1) =  auxPos1-posArr(1,nArr(n))
                IF (auxDX(1).GE. 2*xBin) auxDX(1) = auxDX(1)-lx
                IF (auxDX(1).LE.-2*xBin) auxDX(1) = auxDX(1)+lx
                auxDX(3) =  posArr(3,nArr(l))-posArr(3,nArr(n))
                IF (auxDX(3).GE. 2*zBin) auxDX(3) = auxDX(3)-lz
                IF (auxDX(3).LE.-2*zBin) auxDX(3) = auxDX(3)+lz
                dist = SQRT(auxDX(1)**2 + auxDX(2)**2 + auxDX(3)**2)
                IF (dist .LE. rSearch) THEN
                  CALL compInfVal(posArr(:,nArr(n))+auxDX,velInf2,omeInf2,eetInf2)
                  CALL yInteractionDyn(nArr(n),nArr(l)                      &
                                      ,rArr(nArr(n)),rArr(nArr(l)),auxDX    &
                                      ,posArr(:,nArr(n)),velArr(:,nArr(n))  &
                                      ,omeArr(:,nArr(n))                    &
                                      ,auxV2,omeArr(:,nArr(l)),dist         &
                                      ,velInf1,omeInf1,eetInf1              &
                                      ,velInf2,omeInf2,eetInf2              &
                                      ,forIJ,forJI,torIJ,torJI              &
                                      ,sigL,sigCN,sigCT,sigR,sigA)
                  forArr(:,nArr(n)) = forArr(:,nArr(n)) + forIJ
                  torArr(:,nArr(n)) = torArr(:,nArr(n)) + torIJ
                  forAux(:,ithread+1,nArr(l)) = forAux(:,ithread+1,nArr(l)) + forJI
                  torAux(:,ithread+1,nArr(l)) = torAux(:,ithread+1,nArr(l)) + torJI
                  sigmaL  = sigmaL  + sigL
                  sigmaCN = sigmaCN + sigCN
                  sigmaCT = sigmaCT + sigCT
                  sigmaR  = sigmaR  + sigR
                  sigmaA  = sigmaA  + sigA
                ENDIF
              ENDDO
              ENDIF
    
            ENDDO
          ENDDO
        ENDDO
    
      ENDDO
!$OMP END DO

!-----------------------------------------------------------------------

      sigmaT = sigmaT  + sigmaL + sigmaCN    &
             + sigmaCT + sigmaR + sigmaA
 
!-----------------------------------------------------------------------

!$OMP DO SCHEDULE(static) PRIVATE(j)
      DO i=1,nPar
        DO j = 1,nthreads
          forArr(:,i)  = forArr(:,i) + forAux(:,j,i)
          torArr(:,i)  = torArr(:,i) + torAux(:,j,i)
        ENDDO

        accArr(:,i)  = invMass(i)*forArr(:,i)
        omedArr(:,i) = invIner(i)*torArr(:,i)

        velArr(:,i) = velArrOld(:,i) + 0.5*dtime*(accArr(:,i)+accArrOld(:,i))
        omeArr(:,i) = omeArrOld(:,i) + 0.5*dtime*(omedArr(:,i)+omedArrOld(:,i))

        IF (posArrOld(2,i).LT.0.) THEN 
          velArr(1,i) = velArr(1,i) + gammaDot*ly
        ELSEIF (posArrOld(2,i).GE.ly) THEN 
          velArr(1,i) = velArr(1,i) - gammaDot*ly
        ENDIF

#if defined _FRICTION_ || defined _ROLLING_
        ctcArr(:,i)    = ctcArrNew(:,i)
        ctcArrNew(:,i) = 0

        CALL heapsort(ctcArr(:,i),ctcArrIdx(:,i),idxTest)
#ifdef _FRICTION_
        DO j = 1,nAroundMax
          csiX(j,i) = csiXNew(ctcArrIdx(j,i),i)
          csiY(j,i) = csiYNew(ctcArrIdx(j,i),i)
          csiZ(j,i) = csiZNew(ctcArrIdx(j,i),i)
        ENDDO
#endif
#ifdef _ROLLING_
        DO j = 1,nAroundMax
          csiRX(j,i) = csiRXNew(ctcArrIdx(j,i),i)
          csiRY(j,i) = csiRYNew(ctcArrIdx(j,i),i)
          csiRZ(j,i) = csiRZNew(ctcArrIdx(j,i),i)
        ENDDO
#endif

#ifdef _FRICTION_
        csiXNew(:,i) = 0.
        csiYNew(:,i) = 0.
        csiZNew(:,i) = 0.
#endif
#ifdef _ROLLING_
        csiRXNew(:,i) = 0.
        csiRYNew(:,i) = 0.
        csiRZNew(:,i) = 0.
#endif
#endif
      ENDDO
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(static)
      DO i=1,nBins
        cnt(i) = 0
        cnt2(:,i) = 0
      ENDDO
!$OMP END DO

 
!-----------------------------------------------------------------------

      finish = OMP_GET_WTIME()

      cpuTimeIter  = finish-start
      cpuTimeCumul = cpuTimeCumul + cpuTimeIter/3600.

      start = finish

!$OMP SINGLE
      IF (MOD(it,itRunTime).EQ.0)                                      &
        WRITE(1,9997)  time,sigmaT(1,2),sigmaL(1,2),sigmaCN(1,2)       &
                           ,sigmaCT(1,2),sigmaR(1,2),sigmaA(1,2)

      IF (MOD(it,itStat).EQ.0 .AND. time.GE.tStat) CALL updateStats

      IF (MOD(it,itSave).EQ.0 .OR. time.GE.tmax) THEN
        CALL writeRestart(time,it)
        WRITE(*,*) "Restart file written"
      ELSEIF (cpuTimeCumul.GE.tWall .AND. flagRes) THEN
        CALL writeRestart(time,it)
        flagRes = .FALSE.
        WRITE(*,*) "Restart file written"
      ENDIF

      IF (MOD(it,itField).EQ.0) THEN
        CALL writeField(time,it)
        WRITE(*,*) "Field file written"
      ENDIF

      IF (MOD(it,itPrint).EQ.0)   &
        WRITE(*,9998)  time,sigmaT(1,2)/(vis*gammaDot)                          &
                      ,-(sigmaT(1,1)+sigmaT(2,2)+sigmaT(3,3))/(3*vis*gammaDot)  &
                      ,cpuTimeCumul,cpuTimeIter,it
!$OMP END SINGLE NOWAIT

    ENDDO
!$OMP END PARALLEL

    CLOSE(1)

    OPEN(1,FILE='stressTensor.dat')
    WRITE(1,*)  sigmaT_
    WRITE(1,*)  sigmaL_
    WRITE(1,*)  sigmaCN_
    WRITE(1,*)  sigmaCT_
    WRITE(1,*)  sigmaR_
    WRITE(1,*)  sigmaA_
    CLOSE(1)
 
    DEALLOCATE(aux)

    DEALLOCATE(iBin,nArr,cnt,sumPS,iBS,iBE,cnt2)

    DEALLOCATE(posArrOld,velArrOld,accArrOld  &
              ,omeArrOld,omedArrOld           &
              ,invMass,invIner                &
              ,forAux,torAux)
#if defined _FRICTION_ || defined _ROLLING_
     DEALLOCATE(ctcArrIdx,idxTest)
#endif

    RETURN

 9999 FORMAT (/'-- Time Unit --  ---- eta_r ----  ---- eta_n ----' &
           '  ||-- CPUt [h]  --  - CPUt/it [s] -  ---- it ----')
 9998 FORMAT (7X, f8.5, 8X, f9.5, 8X, f9.5, 2X, '||'               &
             ,6X, f9.5, 8X, f9.5, 5X, i9)
 9997 FORMAT (2X, f12.7, 6(4X, f25.12))

  END SUBROUTINE timeIntegration

!======================================================================

  SUBROUTINE compInfVal(pos,velInf,omeInf,eetInf)

    USE mod_param,   ONLY :  UinfS,waveNS

!----------------------------------------------------------------------

    IMPLICIT NONE

    REAL, DIMENSION(:),   INTENT(IN)    :: pos
    REAL, DIMENSION(:),   INTENT(INOUT) :: velInf,omeInf
    REAL, DIMENSION(:,:), INTENT(INOUT) :: eetInf

#if defined(_FLOW_SHEAR_)
    velInf(1)   =  gammaDot*pos(2)
    omeInf(3)   = -0.5*gammaDot
    eetInf(1,2) =  0.5*gammaDot
    eetInf(2,1) =  eetInf(1,2)
#elif defined(_FLOW_SHEAR_WAVEXY_)
    velInf(1)   =  gammaDot*pos(2) + UinfS*SIN(waveNS*pos(2))
    omeInf(3)   = -0.5*(gammaDot + UinfS*waveNS*COS(waveNS*pos(2)))
    eetInf(1,2) =  0.5*(gammaDot + UinfS*waveNS*COS(waveNS*pos(2)))
    eetInf(2,1) =  eetInf(1,2)
#elif defined(_FLOW_SHEAR_WAVEYX_)
    velInf(1)   =  gammaDot*pos(2) 
    velInf(2)   =  UinfS*SIN(waveNS*pos(1))
    omeInf(3)   = -0.5*(gammaDot - UinfS*waveNS*COS(waveNS*pos(1)))
    eetInf(1,2) =  0.5*(gammaDot + UinfS*waveNS*COS(waveNS*pos(1)))
    eetInf(2,1) =  eetInf(1,2)
#elif defined(_FLOW_SHEAR_WAVEXZ_)
    velInf(1)   =  gammaDot*pos(2) + UinfS*SIN(waveNS*pos(3))
    omeInf(2)   =  0.5*UinfS*waveNS*COS(waveNS*pos(3))
    omeInf(3)   = -0.5*gammaDot
    eetInf(1,2) =  0.5*gammaDot 
    eetInf(1,3) =  0.5*UinfS*waveNS*COS(waveNS*pos(3))
    eetInf(2,1) =  eetInf(1,2)
    eetInf(3,1) =  eetInf(1,3)
#elif defined(_FLOW_SHEAR_WAVEZX_)
    velInf(1)   =  gammaDot*pos(2) 
    velInf(3)   =  UinfS*SIN(waveNS*pos(1))
    omeInf(2)   = -0.5*UinfS*waveNS*COS(waveNS*pos(1))
    omeInf(3)   = -0.5*gammaDot
    eetInf(1,2) =  0.5*gammaDot 
    eetInf(1,3) =  0.5*UinfS*waveNS*COS(waveNS*pos(1))
    eetInf(2,1) =  eetInf(1,2)
    eetInf(3,1) =  eetInf(1,3)
#elif defined(_FLOW_SHEAR_WAVEYZ_)
    velInf(1)   =  gammaDot*pos(2)
    velInf(2)   =  UinfS*SIN(waveNS*pos(3))
    omeInf(1)   = -0.5*UinfS*waveNS*COS(waveNS*pos(3))
    omeInf(3)   = -0.5*gammaDot
    eetInf(1,2) =  0.5*gammaDot
    eetInf(2,3) =  0.5*UinfS*waveNS*COS(waveNS*pos(3))
    eetInf(2,1) =  eetInf(1,2)
    eetInf(3,2) =  eetInf(2,3)
#elif defined(_FLOW_SHEAR_WAVEZY_)
    velInf(1)   =  gammaDot*pos(2)
    velInf(3)   =  UinfS*SIN(waveNS*pos(2))
    omeInf(1)   =  0.5*UinfS*waveNS*COS(waveNS*pos(2))
    omeInf(3)   = -0.5*gammaDot
    eetInf(1,2) =  0.5*gammaDot
    eetInf(2,3) =  0.5*UinfS*waveNS*COS(waveNS*pos(2))
    eetInf(2,1) =  eetInf(1,2)
    eetInf(3,2) =  eetInf(2,3)
#endif

!----------------------------------------------------------------------

    RETURN

!----------------------------------------------------------------------

  END SUBROUTINE compInfVal

!=======================================================================
#ifndef _FLOW_SHEAR_
  SUBROUTINE compInfVelLap(pos,velInfLap)

    USE mod_param,   ONLY :  UinfS,waveNS,waveNS2

!----------------------------------------------------------------------

    IMPLICIT NONE

    REAL, DIMENSION(:), INTENT(IN)    :: pos
    REAL, DIMENSION(:), INTENT(INOUT) :: velInfLap

#if defined(_FLOW_SHEAR_WAVEXY_)
    velInfLap(1) = -UinfS*waveNS2*SIN(waveNS*pos(2))
#elif defined(_FLOW_SHEAR_WAVEYX_)
    velInfLap(2) = -UinfS*waveNS2*SIN(waveNS*pos(1))
#elif defined(_FLOW_SHEAR_WAVEXZ_)
    velInfLap(1) = -UinfS*waveNS2*SIN(waveNS*pos(3))
#elif defined(_FLOW_SHEAR_WAVEZX_)
    velInfLap(3) = -UinfS*waveNS2*SIN(waveNS*pos(1))
#elif defined(_FLOW_SHEAR_WAVEYZ_)
    velInfLap(2) = -UinfS*waveNS2*SIN(waveNS*pos(3))
#elif defined(_FLOW_SHEAR_WAVEZY_)
    velInfLap(3) = -UinfS*waveNS2*SIN(waveNS*pos(2))
#endif

!----------------------------------------------------------------------

    RETURN

!----------------------------------------------------------------------

  END SUBROUTINE compInfVelLap
#endif
!=======================================================================

  SUBROUTINE writeRestart(time,it)

    USE mod_param,   ONLY :  itSave,nSave

    IMPLICIT NONE

    CHARACTER(LEN=3)  ::  iRestart
    INTEGER           ::  it
    REAL              ::  time

    WRITE(iRestart,'(i3.3)') MOD((it-1)/itSave,nSave)

    OPEN(999,FILE='simRestart-'//iRestart//'.bin',FORM='unformatted')
    WRITE(999) time,it
    WRITE(999) posArr
    WRITE(999) velArr
    WRITE(999) accArr
    WRITE(999) omeArr
    WRITE(999) omedArr
#if defined _FRICTION_ || defined _ROLLING_
    WRITE(999) ctcArr
#ifdef _FRICTION_
    WRITE(999) csiX,csiY,csiZ
#endif
#ifdef _ROLLING_
    WRITE(999) csiRX,csiRY,csiRZ
#endif
#endif
    WRITE(999) ntStat
    WRITE(999) sigmaT_
    WRITE(999) sigmaL_
    WRITE(999) sigmaCN_
    WRITE(999) sigmaCT_
    WRITE(999) sigmaR_
    WRITE(999) sigmaA_
    CLOSE(999)

    OPEN(999,FILE='simRestart.lst')
    WRITE(999,*) iRestart
    CLOSE(999)

    RETURN

  END SUBROUTINE writeRestart

!=======================================================================

  SUBROUTINE writeField(time,it)

    IMPLICIT NONE

    CHARACTER(LEN=9)  ::  fIt
    INTEGER           ::  it
    REAL              ::  time

    WRITE(fIt,'(i9.9)') it
    OPEN(998,FILE='simField-'//fIt//'.bin',FORM='unformatted')
    WRITE(998) time
    WRITE(998) posArr
    WRITE(998) velArr
    WRITE(998) accArr
    WRITE(998) omeArr
    WRITE(998) omedArr
    CLOSE(998)

    RETURN

  END SUBROUTINE writeField

!=======================================================================

  SUBROUTINE heapsort(arr,idx,idxTest)

    IMPLICIT NONE
 
    INTEGER, INTENT(INOUT)  ::  arr(0:),idx(0:)
    INTEGER, INTENT(IN)     ::  idxTest(0:)

    INTEGER  ::  start,N,bottom,temp,tempI
 
    N   = SIZE(arr)
    idx = idxTest
    DO start = (N-2)/2,0,-1
      CALL siftdown(arr,start,N,idx)
    ENDDO
 
    DO bottom = N-1,1,-1
      temp        = arr(0)
      tempI       = idx(0)
      arr(0)      = arr(bottom)
      idx(0)      = idx(bottom)
      arr(bottom) = temp
      idx(bottom) = tempI
      CALL siftdown(arr,0,bottom,idx)
    ENDDO
 
  END SUBROUTINE heapsort

!=======================================================================
 
  SUBROUTINE siftdown(arr,start,bottom,idx)

    IMPLICIT NONE
 
    INTEGER, INTENT(INOUT)  ::  arr(0:),idx(0:)
    INTEGER, INTENT(IN)     ::  start,bottom

    INTEGER  ::  child,root,temp,tempI
 
    root = start
    DO WHILE(root*2 + 1 < bottom)
      child = root*2 + 1
 
      IF (child + 1 < bottom) THEN
        IF (arr(child) < arr(child+1)) child = child + 1
      ENDIF
 
      IF (arr(root) < arr(child)) THEN
        temp       = arr(child)
        tempI      = idx(child)
        arr(child) = arr(root)
        idx(child) = idx(root)
        arr(root)  = temp
        idx(root)  = tempI
        root       = child
      ELSE
        RETURN
      ENDIF  
    ENDDO      
 
  END SUBROUTINE siftdown

!=======================================================================

  END MODULE mod_timeInt

!=======================================================================
