!=======================================================================

#include "config.h"

  MODULE mod_dest
 
    IMPLICIT NONE

    PUBLIC   destroyArr

  CONTAINS

!=======================================================================

  SUBROUTINE destroyArr
    
    !Initialize the numerical variables

    USE mod_common,  ONLY :  posArr,rArr         &
                            ,velArr,accArr       &
                            ,omeArr,omedArr      &
                            ,forArr,torArr       &
                            ,masArr,ineArr
#if defined _FRICTION_ || defined _ROLLING_
    USE mod_common,  ONLY :  idxCtcArr,ctcArr,ctcArrNew
#endif
#ifdef _FRICTION_
    USE mod_common,  ONLY :  csiX,csiY,csiZ      &
                            ,csiXNew,csiYNew,csiZNew
#endif
#ifdef _ROLLING_
    USE mod_common,  ONLY :  csiRX,csiRY,csiRZ   &
                            ,csiRXNew,csiRYNew,csiRZNew
#endif

    IMPLICIT NONE

    DEALLOCATE(posArr,rArr,velArr,accArr      &
              ,omeArr,omedArr,forArr,torArr   &
              ,masArr,ineArr)

#if defined _FRICTION_ || defined _ROLLING_
    DEALLOCATE(idxCtcArr,ctcArr,ctcArrNew)
#endif
#ifdef _FRICTION_
    DEALLOCATE(csiX,csiY,csiZ,csiXNew,csiYNew,csiZNew)
#endif

#ifdef _ROLLING_
    DEALLOCATE(csiRX,csiRY,csiRZ,csiRXNew,csiRYNew,csiRZNew)
#endif

    RETURN

  END SUBROUTINE destroyArr

!=======================================================================

  END MODULE mod_dest

!=======================================================================
