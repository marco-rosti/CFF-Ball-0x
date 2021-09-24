!=======================================================================

  PROGRAM testParticles

    USE mod_param,   ONLY : assignParam,KILLSIGNUM
    USE mod_init,    ONLY : initParticlesArr
    USE mod_timeInt, ONLY : timeIntegration
    USE mod_dest,    ONLY : destroyArr

    IMPLICIT NONE

    CALL signal(KILLSIGNUM,destroyArr)
    CALL assignParam
    CALL initParticlesArr
    CALL timeIntegration
    CALL destroyArr

  END PROGRAM

!=======================================================================
