!=======================================================================

  MODULE mod_stats
    
    !Define common variables
    
    USE mod_common, ONLY :  sigmaT,sigmaT_         &
                           ,sigmaL,sigmaL_         &
                           ,sigmaCN,sigmaCN_       &
                           ,sigmaCT,sigmaCT_       &
                           ,sigmaR,sigmaR_         &
                           ,sigmaA,sigmaA_         &
                           ,ntStat

    IMPLICIT NONE

    PRIVATE
    PUBLIC   initStats,updateStats

    CONTAINS

!=======================================================================

    SUBROUTINE initStats

      IMPLICIT NONE

      ntStat = 0

      sigmaL_    = 0.
      sigmaCN_   = 0.
      sigmaCT_   = 0.
      sigmaR_    = 0.
      sigmaA_    = 0.
      sigmaT_    = 0.

      RETURN

    END SUBROUTINE initStats

!=======================================================================

    SUBROUTINE updateStats

      IMPLICIT NONE

      REAL  ::  qnsamp

      ntStat = ntStat + 1
      qnsamp = 1./ntStat

      sigmaL_    = sigmaL_    + (sigmaL-sigmaL_)*qnsamp
      sigmaCN_   = sigmaCN_   + (sigmaCN-sigmaCN_)*qnsamp
      sigmaCT_   = sigmaCT_   + (sigmaCT-sigmaCT_)*qnsamp
      sigmaR_    = sigmaR_    + (sigmaR-sigmaR_)*qnsamp
      sigmaA_    = sigmaA_    + (sigmaA-sigmaA_)*qnsamp
      sigmaT_    = sigmaT_    + (sigmaT-sigmaT_)*qnsamp

      RETURN

    END SUBROUTINE updateStats

!=======================================================================

  END MODULE mod_stats

!=======================================================================
