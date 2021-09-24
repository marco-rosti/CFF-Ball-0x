!=======================================================================

! Choose the pairwise intra-particle forces

! Lubrication    > _LUB_
! Lubrication NO > _LUB_NO_
#define _LUB_

! Frictonless contact    > _CONTACT_
! Frictonless contact NO > _CONTACT_NO_
#define _CONTACT_

! Frictonal contact    > _FRICTION_
! Frictonal contact NO > _FRICTION_NO_
#define _FRICTION_NO_

! Rolling contact    > _ROLLING_
! Rolling contact NO > _ROLLING_NO_
#define _ROLLING_NO_

! Electrostatic repulsive    > _ESR_
! Electrostatic repulsive NO > _ESR_NO_
#define _ESR_

! Van der Walls attractive    > _VDW_
! Van der Walls attractive NO > _VDW_NO_
#define _VDW_NO_

!-----------------------------------------------------------------------

! Choose the flow type

! Plain shear flow > _FLOW_SHEAR_
! Plain shear flow + wave xy (ux = ux(y)) > _FLOW_SHEAR_WAVEXY_
! Plain shear flow + wave yx (uy = uy(x)) > _FLOW_SHEAR_WAVEYX_
! Plain shear flow + wave xz (ux = ux(z)) > _FLOW_SHEAR_WAVEXZ_
! Plain shear flow + wave zx (uz = uz(x)) > _FLOW_SHEAR_WAVEZX_
! Plain shear flow + wave yz (uy = uy(z)) > _FLOW_SHEAR_WAVEYZ_
! Plain shear flow + wave zy (uz = uz(y)) > _FLOW_SHEAR_WAVEZY_
#define _FLOW_SHEAR_

!=======================================================================
