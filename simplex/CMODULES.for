       MODULE contact

        IMPLICIT NONE

        REAL*8 :: normDirFlag, f_ext
        INTEGER  ::  totNumElem, nPtsBM
        INTEGER, DIMENSION(:), ALLOCATABLE  :: BM_num, mark, cnt_xh
        INTEGER, DIMENSION(:,:), ALLOCATABLE :: ElemNeig
        REAL*8, DIMENSION(:), ALLOCATABLE :: BM_x, BM_y, BM_z
        REAL*8, DIMENSION(:), ALLOCATABLE :: BM_x0, BM_y0, BM_z0
        REAL*8, DIMENSION(:), ALLOCATABLE :: BM_x1, BM_y1, BM_z1
        REAL*8, DIMENSION(:), ALLOCATABLE :: BM_x2, BM_y2, BM_z2
        REAL*8, DIMENSION(:), ALLOCATABLE :: BM_x3, BM_y3, BM_z3
        REAL*8, DIMENSION(:), ALLOCATABLE :: ElemNorm_x, ElemNorm_y, 
     &                                       ElemNorm_z
        REAL*8, DIMENSION(:), ALLOCATABLE :: ElemCentx, ElemCenty,
     &                                       ElemCentz

        REAL, DIMENSION(:), ALLOCATABLE :: xNormBM, yNormBM,
     &                                     zNormBM

        REAL*8  ::  pointOutsideX,pointOutsideY,
     &              pointOutsideZ

        !REAL*8, DIMENSION(:), ALLOCATABLE :: cnt_d
        !INTEGER, DIMENSION(:), ALLOCATABLE :: ElemNum

        REAL, PARAMETER :: cnt_h  = 0.08
!        REAL, PARAMETER :: cnt_alpha  = 0.7
        REAL, PARAMETER :: cnt_c  = 0.24
        REAL, PARAMETER :: cnt_k  = 0.40

        INTEGER, PARAMETER :: cnt_n  = 2709

        INTEGER, DIMENSION(:,:), ALLOCATABLE  :: ElemC

        REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: pks2, strn
        REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: cauchy

        INTEGER  ::  flag12, flag13, flag23

        INTEGER, DIMENSION(:), ALLOCATABLE  :: cnt_flg12,cnt_flg21
        INTEGER, DIMENSION(:), ALLOCATABLE  :: cnt_flg13,cnt_flg31
        INTEGER, DIMENSION(:), ALLOCATABLE  :: cnt_flg23,cnt_flg32

        REAL*8, DIMENSION(:), ALLOCATABLE  :: BI_normX
        REAL*8, DIMENSION(:), ALLOCATABLE  :: BI_normY
        REAL*8, DIMENSION(:), ALLOCATABLE  :: BI_normZ

        INTEGER  :: flg

       END MODULE
