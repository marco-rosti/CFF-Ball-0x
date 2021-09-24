!=======================================================================

#include "config.h"

  MODULE mod_param

!----------------------------------------------------------------------
    
    IMPLICIT NONE

    !Signal handling: terminate process(es)
    INTEGER, PARAMETER  ::  KILLSIGNUM=15

    !Area of the unit circle
    REAL,    PARAMETER  ::  pi=ACOS(-1.)

    CHARACTER(LEN=11) ::  radFilename & !Filename of the radii distrib.
                         ,iniFilename   !Filename of the initial condition

    INTEGER  ::  nPar                 & !Number of particles
                ,modop                & !Modus operandi (see init.f90)
                ,nxBin,nyBin,nzBin    & !nBins
                ,nWaveNS              & !Sin shear flow - number of waves
                ,nSave                & !Max number of restart files saved
                ,itSave,itField       & !Saving frequencies
                ,itRuntime,itPrint    &
                ,itStat               &
                ,nThreads               !Number of threads

    REAL     ::  dtime,tmax           & !Time parameters
                ,lx,ly,lz             & !Domain parameters
                ,rhoPar               & !Particles density
                ,xBin,yBin,zBin       & !Bins size
                ,rSearch              & !Searching neigh. size
                ,gammaDot,vis         & !Uniform shear flow parameters
                ,UinfS,waveNS         & !Sin shear flow - ampl,waveNum
                ,waveNS2              & !Sin shear flow - square waveNum
                ,lubOuterDist         & !Lubrication outer influence
                ,kN,gammaN            & !Contact parameters
                ,kT,muC               &
                ,kR,muR               &
                ,fER,kappaDb          & !Electrostatic parameters
                ,lambdaDb             &
                ,hamA,eps             & !Van der Waals parameters
                ,tStat

    PUBLIC assignParam

!----------------------------------------------------------------------

  CONTAINS

!=======================================================================

  SUBROUTINE assignParam

!----------------------------------------------------------------------
    
    IMPLICIT NONE

    REAL  ::  kTnum,kTden          & !Contact parameters
             ,kRnum,kRden

!----------------------------------------------------------------------

    !Read the parameters from file    
    OPEN(1,FILE='p.in',FORM='formatted')
    READ(1,*) radFilename,iniFilename
    READ(1,*) nPar
    READ(1,*) modop
    READ(1,*) dtime,tmax
    READ(1,*) lx,ly,lz
    READ(1,*) rhoPar
    READ(1,*) nxBin,nyBin,nzBin
    READ(1,*) gammaDot,vis
    READ(1,*) UinfS,nWaveNS
    READ(1,*) lubOuterDist
    READ(1,*) kN,gammaN
    READ(1,*) kTnum,kTden,muC
    READ(1,*) kRnum,kRden,muR
    READ(1,*) fER,lambdaDb
    READ(1,*) hamA,eps
    READ(1,*) nSave
    READ(1,*) itSave
    READ(1,*) itField
    READ(1,*) itRuntime
    READ(1,*) itPrint
    READ(1,*) itStat
    READ(1,*) tStat
    READ(1,*) nThreads
    CLOSE(1)

    xBin = lx/nxBin
    yBin = ly/nyBin
    zBin = lz/nzBin
    rSearch = MINVAL((/xBin,yBin,zBin/))
#if defined (_FLOW_SHEAR_WAVEYX_) || defined (_FLOW_SHEAR_WAVEZX_)
    waveNS  = 2*pi*nWaveNS/lx
#elif defined (_FLOW_SHEAR_WAVEXY_) || defined (_FLOW_SHEAR_WAVEZY_)
    waveNS  = 2*pi*nWaveNS/ly
#elif defined (_FLOW_SHEAR_WAVEYZ_) || defined (_FLOW_SHEAR_WAVEXZ_)
    waveNS  = 2*pi*nWaveNS/lz
#endif
    waveNS2 = waveNS*waveNS
    kT = kTnum/kTden*kN
    kR = kRnum/kRden*kN
    kappaDb = 1./lambdaDb

!----------------------------------------------------------------------

    RETURN

!----------------------------------------------------------------------

  END SUBROUTINE assignParam

!=======================================================================

  END MODULE mod_param

!=======================================================================
