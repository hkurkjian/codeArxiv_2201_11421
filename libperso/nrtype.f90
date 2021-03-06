MODULE nrtype
        INTEGER, PARAMETER :: I4B = SELECTED_INT_KIND(9)
        INTEGER, PARAMETER :: I2B = SELECTED_INT_KIND(4)
        INTEGER, PARAMETER :: I1B = SELECTED_INT_KIND(2)
        INTEGER, PARAMETER :: SP = KIND(1.0)
        INTEGER, PARAMETER :: DP = KIND(1.0D0)
        INTEGER, PARAMETER :: QP = SELECTED_REAL_KIND(33,4931)
        INTEGER, PARAMETER :: SPC = KIND((1.0,1.0))
        INTEGER, PARAMETER :: DPC = KIND((1.0D0,1.0D0))
        INTEGER, PARAMETER :: QPC = KIND((1.0_qp,1.0_qp))
        INTEGER, PARAMETER :: LGT = KIND(.true.)
        REAL(QP), PARAMETER :: PI=3.141592653589793238462643383279502884197_qp
        REAL(QP), PARAMETER :: PIs2=1.57079632679489661923132169163975144209858_qp
        REAL(SP), PARAMETER :: PIO2=1.57079632679489661923132169163975144209858_sp
        REAL(SP), PARAMETER :: TWOPI=6.283185307179586476925286766559005768394_sp
        REAL(SP), PARAMETER :: SQRT2=1.41421356237309504880168872420969807856967_sp
        REAL(SP), PARAMETER :: EULER=0.5772156649015328606065120900824024310422_sp
        REAL(DP), PARAMETER :: PI_D=3.141592653589793238462643383279502884197_dp
        REAL(DP), PARAMETER :: PIO2_D=1.57079632679489661923132169163975144209858_dp
        REAL(DP), PARAMETER :: TWOPI_D=6.283185307179586476925286766559005768394_dp
        COMPLEX(DPC), PARAMETER :: ii =dcmplx(0.0_dp,1.0_dp)
        COMPLEX(QPC), PARAMETER :: iiq=cmplx(0.0_qp,1.0_qp,QPC)
        COMPLEX(QPC), PARAMETER :: pis2c=cmplx(PI/2.0_qp,0.0_qp,QPC)
        REAL(QP), PARAMETER :: bidon=-1.0e300_qp
        INTEGER,  PARAMETER :: intbidon=-10000000
        CHARACTER(len=250)  ::     chainebidon
        TYPE sprs2_sp
                INTEGER(I4B) :: n,len
                REAL(SP), DIMENSION(:), POINTER :: val
                INTEGER(I4B), DIMENSION(:), POINTER :: irow
                INTEGER(I4B), DIMENSION(:), POINTER :: jcol
        END TYPE sprs2_sp
        TYPE sprs2_dp
                INTEGER(I4B) :: n,len
                REAL(DP), DIMENSION(:), POINTER :: val
                INTEGER(I4B), DIMENSION(:), POINTER :: irow
                INTEGER(I4B), DIMENSION(:), POINTER :: jcol
        END TYPE sprs2_dp
END MODULE nrtype
