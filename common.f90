!=======================================================================

#include "config.h"

  MODULE mod_common
    
    !Define common variables

    IMPLICIT NONE

    !Particles position, velocities and accelarations
    REAL, DIMENSION(:,:),ALLOCATABLE  ::  posArr                 &
                                         ,velArr,accArr          &
                                         ,omeArr,omedArr         &
                                         ,forArr,torArr!(3,nPar)

    !Particles radius, mass and inertia moment
    REAL, DIMENSION(:), ALLOCATABLE    ::  rArr,masArr,ineArr

    !Particles minimum radius
    REAL  ::  radMin

    !Max number of the smallest particles around the largest one
    INTEGER  ::  nAroundMax

    !Domain volume
    REAL  ::  vol

    !Time
    INTEGER  ::  it !time iteration
    REAL     ::  time

    !Stress tensor
    REAL, DIMENSION(3,3)     ::  sigmaT,sigmaL,sigmaCN  &
                                ,sigmaCT,sigmaR,sigmaA
    
    !Tangential contact
#if defined _FRICTION_ || defined _ROLLING_
    INTEGER, DIMENSION(:)  , ALLOCATABLE  ::  idxCtcArr
    INTEGER, DIMENSION(:,:), ALLOCATABLE  ::  ctcArr,ctcArrNew
#endif
#ifdef _FRICTION_
    REAL, DIMENSION(:,:), ALLOCATABLE  ::  csiX,csiY,csiZ      &!Friction
                                          ,csiXNew,csiYNew,csiZNew
#endif
#ifdef _ROLLING_
    REAL, DIMENSION(:,:), ALLOCATABLE  ::  csiRX,csiRY,csiRZ   &!Rolling
                                          ,csiRXNew,csiRYNew,csiRZNew
#endif

    !REAL, DIMENSION(3,25,nPar)  ::  csiMat
    !see "https://mathematica.stackexchange.com/questions/ ...
    !101350/packing-of-smaller-spheres-onto-a-bigger-one"

    !Statistics
    REAL, DIMENSION(3,3)        ::  sigmaT_,sigmaL_,sigmaCN_ &
                                   ,sigmaCT_,sigmaR_,sigmaA_
    INTEGER                     ::  ntStat

  END MODULE mod_common

!=======================================================================
