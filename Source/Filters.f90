!-------------------------------------------------------------------------------------------------------------------------------
! This module contains all the filters
MODULE Filters
!...............................................................................................................................

    IMPLICIT NONE

CONTAINS

    !-------------------------------------------------------------------------------------------------------------------------------
    ! Discrete time Low-Pass Filter
    REAL FUNCTION LPFilter(InputSignal, DT, CornerFreq, iStatus, reset, inst, FilterVar)
    !...............................................................................................................................
        USE DRC_Types, ONLY : FilterVariables

            ! Inputs

        REAL(4), INTENT(IN)         :: InputSignal
        REAL(4), INTENT(IN)         :: DT                       ! time step [s]
        REAL(4), INTENT(IN)         :: CornerFreq               ! corner frequency [rad/s]
        INTEGER(4), INTENT(IN)      :: iStatus                  ! A status flag set by the simulation as follows: 0 if this is the first call, 1 for all subsequent time steps, -1 if this is the final call at the end of the simulation.
        INTEGER(4), INTENT(INOUT)   :: inst                     ! Instance number. Every instance of this function needs to have an unique instance number to ensure instances don't influence each other.
        LOGICAL(4), INTENT(IN)      :: reset                    ! Reset the filter to the input signal
        TYPE(FilterVariables), INTENT(INOUT) :: FilterVar

            ! Local

        ! REAL(4), DIMENSION(99), SAVE    :: InputSignalLast      ! Input signal the last time this filter was called. Supports 99 separate instances.
        ! REAL(4), DIMENSION(99), SAVE    :: OutputSignalLast ! Output signal the last time this filter was called. Supports 99 separate instances.

            ! Initialization

        IF ((iStatus == 0) .OR. reset) THEN
            FilterVar%LPFOutputSignalLast(inst) = InputSignal
            FilterVar%LPFInputSignalLast(inst) = InputSignal
        ENDIF

            ! Body

        LPFilter = (DT*CornerFreq*InputSignal + DT*CornerFreq*FilterVar%LPFInputSignalLast(inst) - (DT*CornerFreq-2.0)*FilterVar%LPFOutputSignalLast(inst))/(DT*CornerFreq+2.0)

            ! Save signals for next time step

        FilterVar%LPFInputSignalLast(inst)  = InputSignal
        FilterVar%LPFOutputSignalLast(inst) = LPFilter
        inst = inst + 1

    END FUNCTION LPFilter
    !-------------------------------------------------------------------------------------------------------------------------------
    ! Discrete time second order Low-Pass Filter
    REAL FUNCTION SecLPFilter(InputSignal, DT, CornerFreq, Damp, iStatus, reset, inst, FilterVar)
    !...............................................................................................................................

        USE DRC_Types, ONLY : FilterVariables
        IMPLICIT NONE

            ! Inputs

        REAL(4), INTENT(IN)         :: InputSignal
        REAL(4), INTENT(IN)         :: DT                       ! time step [s]
        REAL(4), INTENT(IN)         :: CornerFreq               ! corner frequency [rad/s]
        REAL(4), INTENT(IN)         :: Damp                     ! Dampening constant
        INTEGER(4), INTENT(IN)      :: iStatus                  ! A status flag set by the simulation as follows: 0 if this is the first call, 1 for all subsequent time steps, -1 if this is the final call at the end of the simulation.
        INTEGER(4), INTENT(INOUT)   :: inst                     ! Instance number. Every instance of this function needs to have an unique instance number to ensure instances don't influence each other.
        LOGICAL(4), INTENT(IN)      :: reset                    ! Reset the filter to the input signal
        TYPE(FilterVariables), INTENT(INOUT) :: FilterVar

            ! Local

        ! REAL(4), DIMENSION(99), SAVE :: InputSignalLast1    ! Input signal the last time this filter was called. Supports 99 separate instances.
        ! REAL(4), DIMENSION(99), SAVE :: InputSignalLast2    ! Input signal the next to last time this filter was called. Supports 99 separate instances.
        ! REAL(4), DIMENSION(99), SAVE :: OutputSFilterVar%SecLPFignalLast1   ! Output signal the last time this filter was called. Supports 99 separate instances.
        ! REAL(4), DIMENSION(99), SAVE :: OutputSignalLast2   ! Output signal the next to last time this filter was called. Supports 99 separate instances.

            ! Initialization

        IF ((iStatus == 0) .OR. reset )  THEN
            FilterVar%SecLPFOutputSignalLast1(inst)  = InputSignal
            FilterVar%SecLPFOutputSignalLast2(inst)  = InputSignal
            FilterVar%SecLPFInputSignalLast1(inst)   = InputSignal
            FilterVar%SecLPFInputSignalLast2(inst)   = InputSignal
        ENDIF

            ! Body

        SecLPFilter = 1/(4+4*DT*Damp*CornerFreq+DT**2*CornerFreq**2) * ( (8-2*DT**2*CornerFreq**2)*FilterVar%SecLPFOutputSignalLast1(inst) &
                        + (-4+4*DT*Damp*CornerFreq-DT**2*CornerFreq**2)*FilterVar%SecLPFOutputSignalLast2(inst) + (DT**2*CornerFreq**2)*InputSignal &
                            + (2*DT**2*CornerFreq**2)*FilterVar%SecLPFInputSignalLast1(inst) + (DT**2*CornerFreq**2)*FilterVar%SecLPFInputSignalLast2(inst) )

            ! Save signals for next time step

        FilterVar%SecLPFInputSignalLast2(inst)   = FilterVar%SecLPFInputSignalLast1 (inst)
        FilterVar%SecLPFInputSignalLast1(inst)   = InputSignal
        FilterVar%SecLPFOutputSignalLast2(inst)  = FilterVar%SecLPFOutputSignalLast1 (inst)
        FilterVar%SecLPFOutputSignalLast1(inst)  = SecLPFilter
        inst = inst + 1

    END FUNCTION SecLPFilter
    !-------------------------------------------------------------------------------------------------------------------------------
    ! Discrete time High-Pass Filter
    REAL FUNCTION HPFilter( InputSignal, DT, CornerFreq, iStatus, reset, inst, FilterVar)
    !...............................................................................................................................

        USE DRC_Types, ONLY : FilterVariables    

        ! Inputs

        REAL(4), INTENT(IN)     :: InputSignal
        REAL(4), INTENT(IN)     :: DT                       ! time step [s]
        REAL(4), INTENT(IN)     :: CornerFreq               ! corner frequency [rad/s]
        INTEGER, INTENT(IN)     :: iStatus                  ! A status flag set by the simulation as follows: 0 if this is the first call, 1 for all subsequent time steps, -1 if this is the final call at the end of the simulation.
        INTEGER, INTENT(INOUT)  :: inst                     ! Instance number. Every instance of this function needs to have an unique instance number to ensure instances don't influence each other.
        LOGICAL(4), INTENT(IN)  :: reset                    ! Reset the filter to the input signal
        TYPE(FilterVariables), INTENT(INOUT) :: FilterVar

            ! Local

        REAL(4)                         :: K                        ! Constant gain
        ! REAL(4), DIMENSION(99), SAVE    :: InputSignalLast      ! Input signal the last time this filter was called. Supports 99 separate instances.
        ! REAL(4), DIMENSION(99), SAVE    :: OutputSignalLast ! Output signal the last time this filter was called. Supports 99 separate instances.

            ! Initialization

        IF ((iStatus == 0) .OR. reset)  THEN
            FilterVar%HPFOutputSignalLast(inst) = InputSignal
            FilterVar%HPFInputSignalLast(inst) = InputSignal
        ENDIF

        K = 2.0 / DT

            ! Body

        HPFilter = K/(CornerFreq + K)*InputSignal - K/(CornerFreq + K)*FilterVar%HPFInputSignalLast(inst) - (CornerFreq - K)/(CornerFreq + K)*FilterVar%HPFOutputSignalLast(inst)

            ! Save signals for next time step

        FilterVar%HPFInputSignalLast(inst)   = InputSignal
        FilterVar%HPFOutputSignalLast(inst)  = HPFilter
        inst = inst + 1

    END FUNCTION HPFilter
    !-------------------------------------------------------------------------------------------------------------------------------
    ! Discrete time inverted Notch Filter with descending slopes, G = CornerFreq*s/(Damp*s^2+CornerFreq*s+Damp*CornerFreq^2)
    REAL FUNCTION NotchFilterSlopes(InputSignal, DT, CornerFreq, Damp, iStatus, reset, inst, FilterVar)
    !...............................................................................................................................

        USE DRC_Types, ONLY : FilterVariables
        
        ! Inputs

        REAL(4), INTENT(IN)     :: InputSignal
        REAL(4), INTENT(IN)     :: DT                       ! time step [s]
        REAL(4), INTENT(IN)     :: CornerFreq               ! corner frequency [rad/s]
        REAL(4), INTENT(IN)     :: Damp                     ! Dampening constant
        INTEGER, INTENT(IN)     :: iStatus                  ! A status flag set by the simulation as follows: 0 if this is the first call, 1 for all subsequent time steps, -1 if this is the final call at the end of the simulation.
        INTEGER, INTENT(INOUT)  :: inst                     ! Instance number. Every instance of this function needs to have an unique instance number to ensure instances don't influence each other.
        LOGICAL(4), INTENT(IN)  :: reset                    ! Reset the filter to the input signal
        TYPE(FilterVariables), INTENT(INOUT) :: FilterVar

            ! Local

        ! REAL(4), DIMENSION(99), SAVE :: InputSignalLast1    ! Input signal the last time this filter was called. Supports 99 separate instances.
        ! REAL(4), DIMENSION(99), SAVE :: InputSignalLast2    ! Input signal the next to last time this filter was called. Supports 99 separate instances.
        ! REAL(4), DIMENSION(99), SAVE :: OutputSignalLast1   ! Output signal the last time this filter was called. Supports 99 separate instances.
        ! REAL(4), DIMENSION(99), SAVE :: OutputSignalLast2   ! Output signal the next to last time this filter was called. Supports 99 separate instances.

            ! Initialization

        IF ((iStatus == 0) .OR. reset) THEN
            FilterVar%NFSlopesOutputSignalLast1(inst)  = InputSignal
            FilterVar%NFSlopesOutputSignalLast2(inst)  = InputSignal
            FilterVar%NFSlopesInputSignalLast1(inst)   = InputSignal
            FilterVar%NFSlopesInputSignalLast2(inst)   = InputSignal
        ENDIF

            ! Body

        NotchFilterSlopes = 1.0/(4.0+2.0*DT*Damp*CornerFreq+DT**2.0*CornerFreq**2.0) * ( (8.0-2.0*DT**2.0*CornerFreq**2.0)*FilterVar%NFSlopesOutputSignalLast1(inst) &
                        + (-4.0+2.0*DT*Damp*CornerFreq-DT**2.0*CornerFreq**2.0)*FilterVar%NFSlopesOutputSignalLast2(inst) + &
                            (2.0*DT*Damp*CornerFreq)*InputSignal + (-2.0*DT*Damp*CornerFreq)*FilterVar%NFSlopesInputSignalLast2(inst) )

            ! Save signals for next time step

        FilterVar%NFSlopesInputSignalLast2(inst)   = FilterVar%NFSlopesInputSignalLast1(inst)
        FilterVar%NFSlopesInputSignalLast1(inst)   = InputSignal          !Save input signal for next time step
        FilterVar%NFSlopesOutputSignalLast2(inst)  = FilterVar%NFSlopesOutputSignalLast1(inst)      !Save input signal for next time step
        FilterVar%NFSlopesOutputSignalLast1(inst)  = NotchFilterSlopes
        inst = inst + 1

    END FUNCTION NotchFilterSlopes
    !-------------------------------------------------------------------------------------------------------------------------------
    ! Discrete time Notch Filter, G = (s^2 + 2*omega*betaNum*s + omega^2)/(s^2 + 2*omega*betaDen*s + omega^2)
    REAL FUNCTION NotchFilter(InputSignal, DT, omega, betaNum, betaDen, iStatus, reset, inst, FilterVar)
    !...............................................................................................................................

        USE DRC_Types, ONLY : FilterVariables
        
        ! Inputs

        REAL(4), INTENT(IN)     :: InputSignal
        REAL(4), INTENT(IN)     :: DT                       ! time step [s]
        REAL(4), INTENT(IN)     :: omega                    ! corner frequency [rad/s]
        REAL(4), INTENT(IN)     :: betaNum                  ! Dampening constant in numerator of filter transfer function
        REAL(4), INTENT(IN)     :: betaDen                  ! Dampening constant in denominator of filter transfer function
        INTEGER, INTENT(IN)     :: iStatus                  ! A status flag set by the simulation as follows: 0 if this is the first call, 1 for all subsequent time steps, -1 if this is the final call at the end of the simulation.
        INTEGER, INTENT(INOUT)  :: inst                     ! Instance number. Every instance of this function needs to have an unique instance number to ensure instances don't influence each other.
        LOGICAL(4), INTENT(IN)  :: reset                    ! Reset the filter to the input signal
        TYPE(FilterVariables), INTENT(INOUT) :: FilterVar

            ! Local
        REAL(4)                         :: K, P1, P2, P3, P4, P5    ! Constant gain
        ! REAL(4), DIMENSION(99), SAVE    :: InputSignalLast1         ! Input signal the last time this filter was called. Supports 99 separate instances.
        ! REAL(4), DIMENSION(99), SAVE    :: InputSignalLast2         ! Input signal the next to last time this filter was called. Supports 99 separate instances.
        ! REAL(4), DIMENSION(99), SAVE    :: OutputSignalLast1        ! Output signal the last time this filter was called. Supports 99 separate instances.
        ! REAL(4), DIMENSION(99), SAVE    :: OutputSignalLast2        ! Output signal the next to last time this filter was called. Supports 99 separate instances.

            ! Initialization

        IF ((iStatus == 0) .OR. reset) THEN
            FilterVar%NFOutputSignalLast1(inst)  = InputSignal
            FilterVar%NFOutputSignalLast2(inst)  = InputSignal
            FilterVar%NFInputSignalLast1(inst)   = InputSignal
            FilterVar%NFInputSignalLast2(inst)   = InputSignal
        ENDIF

        K = 2/DT
        P1 = (K**2 + 2*omega*BetaNum*K + omega**2)/(K**2 + 2*omega*BetaDen*K + omega**2)
        P2 = (2*omega**2 - 2*K**2)  / (K**2 + 2*omega*BetaDen*K + omega**2);
        P3 = (K**2 - 2*omega*BetaNum*K + omega**2) / (K**2 + 2*omega*BetaDen*K + omega**2)
        P4 = (2*omega**2 - 2*K**2)  / (K**2 + 2*omega*BetaDen*K + omega**2)
        P5 = (K**2 - 2*omega*BetaDen*K + omega**2)/ (K**2 + 2*omega*BetaDen*K + omega**2)
        
        ! Body
        NotchFilter = P1*InputSignal + P2*FilterVar%NFInputSignalLast1(inst) + P3*FilterVar%NFInputSignalLast2(inst) - P4*FilterVar%NFOutputSignalLast1(inst) - P5*FilterVar%NFOutputSignalLast2(inst)

        ! Save signals for next time step
        FilterVar%NFInputSignalLast2(inst)   = FilterVar%NFInputSignalLast1(inst)
        FilterVar%NFInputSignalLast1(inst)   = InputSignal          !Save input signal for next time step
        FilterVar%NFOutputSignalLast2(inst)  = FilterVar%NFOutputSignalLast1(inst)      !Save input signal for next time step
        FilterVar%NFOutputSignalLast1(inst)  = NotchFilter
        inst = inst + 1

    END FUNCTION NotchFilter
    !-------------------------------------------------------------------------------------------------------------------------------
    ! Prefilter measured wind turbine signals to separate the filtering from the actual control actions
    SUBROUTINE PreFilterMeasuredSignals(CntrPar, LocalVar, objInst, FilterVar)
    !...............................................................................................................................

        USE DRC_Types, ONLY : ControlParameters, LocalVariables, ObjectInstances, FilterVariables
        
        TYPE(ControlParameters), INTENT(INOUT)  :: CntrPar
        TYPE(LocalVariables), INTENT(INOUT)     :: LocalVar
        TYPE(ObjectInstances), INTENT(INOUT)    :: objInst
        TYPE(FilterVariables), INTENT(INOUT)    :: FilterVar

        ! Filter the HSS (generator) speed measurement:
        ! Apply Low-Pass Filter (choice between first- and second-order low-pass filter)
        IF (CntrPar%F_LPFType == 1) THEN
            LocalVar%GenSpeedF = LPFilter(LocalVar%GenSpeed, LocalVar%DT, CntrPar%F_LPFCornerFreq, LocalVar%iStatus, .FALSE., objInst%instLPF, FilterVar)
        ELSEIF (CntrPar%F_LPFType == 2) THEN   
            LocalVar%GenSpeedF = SecLPFilter(LocalVar%GenSpeed, LocalVar%DT, CntrPar%F_LPFCornerFreq, CntrPar%F_LPFDamping, LocalVar%iStatus, .FALSE., objInst%instSecLPF, FilterVar) ! Second-order low-pass filter on generator speed
        END IF
        
        IF (CntrPar%F_NotchType == 1) THEN
            LocalVar%GenSpeedF = NotchFilter(LocalVar%GenSpeedF, LocalVar%DT, CntrPar%F_NotchCornerFreq, CntrPar%F_NotchBetaNumDen(1), CntrPar%F_NotchBetaNumDen(2), LocalVar%iStatus, .FALSE., objInst%instNotch, FilterVar) ! Notch filter on generator speed
        END IF
        
        ! Filtering the tower fore-aft acceleration signal 
        LocalVar%FA_AccHPF = HPFilter(LocalVar%FA_Acc, LocalVar%DT, CntrPar%FA_HPFCornerFreq, LocalVar%iStatus, .FALSE., objInst%instHPF, FilterVar)
        
        END SUBROUTINE PreFilterMeasuredSignals
    END MODULE Filters
