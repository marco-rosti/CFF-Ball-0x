!=======================================================================

#include "config.h"

  MODULE mod_init
 
    USE mod_param,   ONLY :  nPar,modop,lx,ly,lz
    USE mod_common,  ONLY :  posArr,rArr,radMin,vol

    IMPLICIT NONE

    PUBLIC   initParticlesArr

  CONTAINS

!=======================================================================

  SUBROUTINE initParticlesArr
    
    !Initialize the numerical variables

    USE mod_param,  ONLY :  rhoPar,gammaDot,pi  &
                           ,UinfS,waveNS,tstat  &
                           ,radFilename,iniFilename
    USE mod_common, ONLY :  velArr,accArr       &
                           ,omeArr,omedArr      &
                           ,forArr,torArr       &
                           ,masArr,ineArr       &
                           ,nAroundMax          &
                           ,ntStat              &
                           ,sigmaT_,sigmaL_     &
                           ,sigmaCN_,sigmaCT_   &
                           ,sigmaR_,sigmaA_     &
                           ,it,time
#if defined _FRICTION_ || defined _ROLLING_
    USE mod_common, ONLY :  idxCtcArr,ctcArr,ctcArrNew
#endif
#ifdef _FRICTION_
    USE mod_common, ONLY :  csiX,csiY,csiZ      &
                           ,csiXNew,csiYNew,csiZNew
#endif
#ifdef _ROLLING_
    USE mod_common, ONLY :  csiRX,csiRY,csiRZ   &
                           ,csiRXNew,csiRYNew,csiRZNew
#endif
    USE mod_stats,  ONLY :  initStats

    IMPLICIT NONE

    LOGICAL           ::  restart_exists
    CHARACTER(LEN=3)  ::  iRestart
    INTEGER           ::  i
    REAL              ::  radMax
    REAL, DIMENSION(:), ALLOCATABLE  ::  volArr

    ALLOCATE(posArr(3,nPar),rArr(nPar)        &
            ,velArr(3,nPar),accArr(3,nPar)    &
            ,omeArr(3,nPar),omedArr(3,nPar)   &
            ,forArr(3,nPar),torArr(3,nPar)    &
            ,masArr(nPar),ineArr(nPar)        &
            ,volArr(nPar))

    vol = lx*ly*lz
    OPEN(1,FILE=radFilename)
    DO i=1,nPar
      READ(1,*) rArr(i)
      volArr(i) = 4*pi*rArr(i)**3/3.
      masArr(i) = rhoPar*volArr(i)
      ineArr(i) = 2*masArr(i)*rArr(i)**2*0.2
    ENDDO

    DEALLOCATE(volArr)

    radMin = MINVAL(rArr)
    radMax = MAXVAL(rArr)

    nAroundMax = INT(4*pi*(radMax+radMin)**2 /          &
                 (2.*2.*SQRT(3.)*radMin**2)) + 1

#if defined _FRICTION_ || defined _ROLLING_
    ALLOCATE(idxCtcArr(nPar),ctcArr(nAroundMax,nPar)    &
            ,ctcArrNew(nAroundMax,nPar))
#endif
#ifdef _FRICTION_
    ALLOCATE(csiX(nAroundMax,nPar),csiY(nAroundMax,nPar)   &
            ,csiZ(nAroundMax,nPar)                         &
            ,csiXNew(nAroundMax,nPar)                      &
            ,csiYNew(nAroundMax,nPar)                      &
            ,csiZNew(nAroundMax,nPar))
#endif

#ifdef _ROLLING_
    ALLOCATE(csiRX(nAroundMax,nPar),csiRY(nAroundMax,nPar) &
            ,csiRZ(nAroundMax,nPar)                        &
            ,csiRXNew(nAroundMax,nPar)                     &
            ,csiRYNew(nAroundMax,nPar)                     &
            ,csiRZNew(nAroundMax,nPar))
#endif

    !Modus operandi 'modop':
    !modop = 0, new simulation;
    !modop = 1, new simulation; initial condition given by 
    !           the restart file
    !modop = 2, continue the simulation; new statistics
    !modop = 3, conitnue the simulation and statistics
    IF (modop.GT.0) THEN

      OPEN(1,FILE='simRestart.lst')
      READ(1,*) i
      CLOSE(1)
      WRITE(iRestart,'(i3.3)') i

      INQUIRE(FILE="simRestart-"//iRestart//".bin"         &
             ,EXIST=restart_exists)
      IF (.NOT.restart_exists) THEN
        WRITE(*,*) 'Error: Restart file missing'
        STOP
      ENDIF

      OPEN(1,FILE='simRestart-'//iRestart//'.bin',FORM='unformatted')
      READ(1) time,it
      READ(1) posArr
      READ(1) velArr
      READ(1) accArr
      READ(1) omeArr
      READ(1) omedArr
#if defined _FRICTION_ || defined _ROLLING_
      READ(1) ctcArr
#ifdef _FRICTION_
      READ(1) csiX,csiY,csiZ
#endif
#ifdef _ROLLING_
      READ(1) csiRX,csiRY,csiRZ
#endif
#endif
      IF (modop.EQ.3 .AND. time.GE.tstat) THEN
        READ(1) ntStat
        READ(1) sigmaT_
        READ(1) sigmaL_
        READ(1) sigmaCN_
        READ(1) sigmaCT_
        READ(1) sigmaR_
        READ(1) sigmaA_
      ELSE
        CALL initStats
        IF (modop.EQ.1) THEN
          time = 0.0
          it   = 0
        ENDIF
      ENDIF
      CLOSE(1)
      
    ELSE

      it        = 0
      time      = 0.0
      velArr    = 0.0
      accArr    = 0.0
      omeArr    = 0.0
      omedArr   = 0.0
#if defined _FRICTION_ || defined _ROLLING_
      ctcArr    = 0
#ifdef _FRICTION_
      csiX      = 0.0
      csiY      = 0.0
      csiZ      = 0.0
#endif
#ifdef _ROLLING_
      csiRX     = 0.0
      csiRY     = 0.0
      csiRZ     = 0.0
#endif
#endif

      CALL initStats

      OPEN(1,FILE=iniFilename)
      DO i=1,nPar
        READ(1,*) posArr(1,i),posArr(2,i),posArr(3,i)
      ENDDO
      CLOSE(1)

      velArr(1,:) = gammaDot*posArr(2,:) +  &
                    UinfS*SIN(waveNS*posArr(2,:))
      omeArr(3,:) = -0.5*( gammaDot +       &
                    UinfS*waveNS*COS(waveNS*posArr(2,:)) )
 
    ENDIF

#if defined _FRICTION_ || defined _ROLLING_
    ctcArrNew = 0
#ifdef _FRICTION_
    csiXNew   = 0.0
    csiYNew   = 0.0
    csiZNew   = 0.0
#endif
#ifdef _ROLLING_
    csiRXNew  = 0.0
    csiRYNew  = 0.0
    csiRZNew  = 0.0
#endif
#endif

    RETURN

  END SUBROUTINE initParticlesArr

!=======================================================================

  END MODULE mod_init

!=======================================================================
