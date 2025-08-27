MODULE ReadSetParameters

    USE, INTRINSIC :: ISO_C_Binding
    IMPLICIT NONE

CONTAINS
    !..............................................................................................................................
    ! Read all constant control parameters from DISCON.IN parameter file
    !..............................................................................................................................
    SUBROUTINE ReadControlParameterFileSub(avrSWAP, accINFILE, CntrPar)
        USE DRC_Types, ONLY : ControlParameters

        INTEGER(4), PARAMETER :: UnControllerParameters = 89
        REAL(C_FLOAT), INTENT(IN)            :: avrSWAP(*)
        CHARACTER(C_CHAR), INTENT(IN) :: accINFILE(NINT(avrSWAP(50)))        ! The name of the parameter input file
        TYPE(ControlParameters), INTENT(INOUT) :: CntrPar

        CHARACTER(LEN=256) :: filename
        INTEGER :: i

        ! Initialize the filename variable
        filename = ' '
        
        ! Convert the character array to a scalar string
        DO i = 1, MIN(LEN(filename), SIZE(accINFILE))
            filename(i:i) = accINFILE(i)
        END DO
        
        OPEN(unit=UnControllerParameters, file=TRIM(filename), status='old', action='read')
        
        !----------------------- HEADER ------------------------
        READ(UnControllerParameters, *)
        READ(UnControllerParameters, *)
        READ(UnControllerParameters, *)

        !----------------------- DEBUG --------------------------
        READ(UnControllerParameters, *) 
        READ(UnControllerParameters, *) CntrPar%LoggingLevel
        READ(UnControllerParameters, *) 

        !----------------- CONTROLLER FLAGS ---------------------
        READ(UnControllerParameters, *)
        READ(UnControllerParameters, *) CntrPar%F_LPFType
        READ(UnControllerParameters, *) CntrPar%F_NotchType
        READ(UnControllerParameters, *) CntrPar%IPC_ControlMode
        READ(UnControllerParameters, *) CntrPar%VS_ControlMode
        READ(UnControllerParameters, *) CntrPar%Y_ControlMode        
        READ(UnControllerParameters, *)

        !----------------- FILTER CONSTANTS ---------------------
        READ(UnControllerParameters, *)        
        READ(UnControllerParameters, *) CntrPar%F_LPFCornerFreq
        READ(UnControllerParameters, *) CntrPar%F_LPFDamping
        READ(UnControllerParameters, *) CntrPar%F_NotchCornerFreq
        ALLOCATE(CntrPar%F_NotchBetaNumDen(2))
        READ(UnControllerParameters,*) CntrPar%F_NotchBetaNumDen
        READ(UnControllerParameters, *)

        !----------- BLADE PITCH CONTROLLER CONSTANTS -----------
        READ(UnControllerParameters, *)
        READ(UnControllerParameters, *) CntrPar%PC_GS_n
        ALLOCATE(CntrPar%PC_GS_angles(CntrPar%PC_GS_n))
        READ(UnControllerParameters,*) CntrPar%PC_GS_angles
        ALLOCATE(CntrPar%PC_GS_KP(CntrPar%PC_GS_n))
        READ(UnControllerParameters,*) CntrPar%PC_GS_KP
        ALLOCATE(CntrPar%PC_GS_KI(CntrPar%PC_GS_n))
        READ(UnControllerParameters,*) CntrPar%PC_GS_KI
        ALLOCATE(CntrPar%PC_GS_KD(CntrPar%PC_GS_n))
        READ(UnControllerParameters,*) CntrPar%PC_GS_KD
        ALLOCATE(CntrPar%PC_GS_TF(CntrPar%PC_GS_n))
        READ(UnControllerParameters,*) CntrPar%PC_GS_TF
        READ(UnControllerParameters, *) CntrPar%PC_MaxPit
        READ(UnControllerParameters, *) CntrPar%PC_MinPit
        READ(UnControllerParameters, *) CntrPar%PC_MaxRat
        READ(UnControllerParameters, *) CntrPar%PC_MinRat
        READ(UnControllerParameters, *) CntrPar%PC_RefSpd
        READ(UnControllerParameters, *) CntrPar%PC_FinePit
        READ(UnControllerParameters, *) CntrPar%PC_Switch
        !-- sine pitch excitiation --
        READ(UnControllerParameters, *) CntrPar%Z_EnableSine
        READ(UnControllerParameters, *) CntrPar%Z_PitchAmplitude
        READ(UnControllerParameters, *) CntrPar%Z_PitchFrequency
        READ(UnControllerParameters, *)

        !------------------- IPC CONSTANTS -----------------------
        READ(UnControllerParameters, *) 
        READ(UnControllerParameters, *) CntrPar%IPC_IntSat
        ALLOCATE(CntrPar%IPC_KI(2))
        READ(UnControllerParameters,*) CntrPar%IPC_KI
        ALLOCATE(CntrPar%IPC_aziOffset(2))
        READ(UnControllerParameters,*) CntrPar%IPC_aziOffset
        READ(UnControllerParameters, *) CntrPar%IPC_CornerFreqAct
        READ(UnControllerParameters, *)

        !------------ VS TORQUE CONTROL CONSTANTS ----------------
        READ(UnControllerParameters, *)
        READ(UnControllerParameters, *) CntrPar%VS_GenEff
        READ(UnControllerParameters, *) CntrPar%VS_ArSatTq
        READ(UnControllerParameters, *) CntrPar%VS_MaxRat
        READ(UnControllerParameters, *) CntrPar%VS_MaxTq
        READ(UnControllerParameters, *) CntrPar%VS_MinTq
        READ(UnControllerParameters, *) CntrPar%VS_MinOMSpd
        READ(UnControllerParameters, *) CntrPar%VS_Rgn2K
        READ(UnControllerParameters, *) CntrPar%VS_RtPwr
        READ(UnControllerParameters, *) CntrPar%VS_RtTq
        READ(UnControllerParameters, *) CntrPar%VS_RefSpd
        READ(UnControllerParameters, *) CntrPar%VS_n
        ALLOCATE(CntrPar%VS_KP(CntrPar%VS_n))
        READ(UnControllerParameters,*) CntrPar%VS_KP
        ALLOCATE(CntrPar%VS_KI(CntrPar%VS_n))
        READ(UnControllerParameters,*) CntrPar%VS_KI
        READ(UnControllerParameters, *)

        !------------ WIND SPEED ESTIMATOR CONTANTS --------------
        READ(UnControllerParameters, *)
        READ(UnControllerParameters, *) CntrPar%WE_BladeRadius
        READ(UnControllerParameters, *) CntrPar%WE_CP_n
        ALLOCATE(CntrPar%WE_CP(CntrPar%WE_CP_n))
        READ(UnControllerParameters, *) CntrPar%WE_CP
        READ(UnControllerParameters, *) CntrPar%WE_Gamma
        READ(UnControllerParameters, *) CntrPar%WE_GearboxRatio
        READ(UnControllerParameters, *) CntrPar%WE_Jtot
        READ(UnControllerParameters, *) CntrPar%WE_RhoAir
        READ(UnControllerParameters, *)

        !-------------- YAW CONTROLLER CONSTANTS -----------------
        READ(UnControllerParameters, *)
        READ(UnControllerParameters, *) CntrPar%Y_ErrThresh
        READ(UnControllerParameters, *) CntrPar%Y_IPC_IntSat
        READ(UnControllerParameters, *) CntrPar%Y_IPC_n
        ALLOCATE(CntrPar%Y_IPC_KP(CntrPar%Y_IPC_n))
        READ(UnControllerParameters,*) CntrPar%Y_IPC_KP
        ALLOCATE(CntrPar%Y_IPC_KI(CntrPar%Y_IPC_n))
        READ(UnControllerParameters,*) CntrPar%Y_IPC_KI
        READ(UnControllerParameters, *) CntrPar%Y_IPC_omegaLP
        READ(UnControllerParameters, *) CntrPar%Y_IPC_zetaLP
        READ(UnControllerParameters, *) CntrPar%Y_MErrSet
        READ(UnControllerParameters, *) CntrPar%Y_omegaLPFast
        READ(UnControllerParameters, *) CntrPar%Y_omegaLPSlow
        READ(UnControllerParameters, *) CntrPar%Y_Rate
        READ(UnControllerParameters, *)

        !------------ FORE-AFT TOWER DAMPER CONSTANTS ------------
        READ(UnControllerParameters, *)      
        READ(UnControllerParameters, *) CntrPar%FA_KI  
        READ(UnControllerParameters, *) CntrPar%FA_HPFCornerFreq
        READ(UnControllerParameters, *) CntrPar%FA_IntSat

        !------------ CONTROLLER RESTART -------------------------
        IF (avrSWAP(1) == 2) THEN
            READ(UnControllerParameters, *)      
            READ(UnControllerParameters, *) CntrPar%RestartFile
        END IF
        
        ! END OF INPUT FILE    
        
        !------------------- CALCULATED CONSTANTS -----------------------
        CntrPar%PC_RtTq99 = CntrPar%VS_RtTq*0.99
        CntrPar%VS_MinOMTq = CntrPar%VS_Rgn2K*CntrPar%VS_MinOMSpd**2
        CntrPar%VS_MaxOMTq = CntrPar%VS_Rgn2K*CntrPar%VS_RefSpd**2
        CntrPar%VS_Rgn3Pitch = CntrPar%PC_FinePit + CntrPar%PC_Switch
        
        CLOSE(UnControllerParameters)
    END SUBROUTINE ReadControlParameterFileSub
    
    SUBROUTINE ComputeVariablesSetpoints(CntrPar, LocalVar)
        USE DRC_Types, ONLY : ControlParameters, LocalVariables
        
        ! Local variables
        TYPE(ControlParameters), INTENT(INOUT)  :: CntrPar
        TYPE(LocalVariables), INTENT(INOUT)     :: LocalVar
        
        ! Calculate yaw misalignment error
        LocalVar%Y_MErr = LocalVar%Y_M + CntrPar%Y_MErrSet ! Yaw-alignment error
        
        ! Compute the current speed AND power error
        LocalVar%PC_SpdErr = CntrPar%PC_RefSpd - LocalVar%GenSpeedF            ! Speed error
        LocalVar%PC_PwrErr = CntrPar%VS_RtPwr - LocalVar%VS_GenPwr             ! Power error
        
        ! XXX
        LocalVar%VS_SpdErrAr = CntrPar%VS_RefSpd - LocalVar%GenSpeedF       ! Current speed error - Above-rated PI-control
        LocalVar%VS_SpdErrBr = CntrPar%VS_MinOMSpd - LocalVar%GenSpeedF     ! Current speed error - Below-rated PI-control
    END SUBROUTINE ComputeVariablesSetpoints
    
    SUBROUTINE ReadAvrSWAP(avrSWAP, LocalVar)
        USE DRC_Types, ONLY : LocalVariables
    
        REAL(C_FLOAT), INTENT(INOUT) :: avrSWAP(*)   ! The swap array, used to pass data to, and receive data from, the DLL controller.
        TYPE(LocalVariables), INTENT(INOUT) :: LocalVar
        
        ! Load variables from calling program (See Appendix A of Bladed User's Guide):
        LocalVar%iStatus = NINT(avrSWAP(1))
        LocalVar%Time = avrSWAP(2)
        LocalVar%DT = avrSWAP(3)
        LocalVar%BlPitch(1) = avrSWAP(4)
        LocalVar%VS_MechGenPwr = avrSWAP(14)
        LocalVar%VS_GenPwr = avrSWAP(15)
        LocalVar%GenSpeed = avrSWAP(20)
        LocalVar%RotSpeed = avrSWAP(21)
        LocalVar%GenTqMeas = avrSWAP(23)
        LocalVar%Y_M = avrSWAP(24)
        LocalVar%HorWindV = avrSWAP(27)
        LocalVar%rootMOOP(1) = avrSWAP(30)
        LocalVar%rootMOOP(2) = avrSWAP(31)
        LocalVar%rootMOOP(3) = avrSWAP(32)
        LocalVar%BlPitch(2) = avrSWAP(33)
        LocalVar%BlPitch(3) = avrSWAP(34)
        LocalVar%FA_Acc = avrSWAP(53)
        LocalVar%Azimuth = avrSWAP(60)
        LocalVar%NumBl = NINT(avrSWAP(61))
    END SUBROUTINE ReadAvrSWAP
    
    SUBROUTINE Assert(LocalVar, CntrPar, avrSWAP, aviFAIL, ErrMsg, size_avcMSG)
        USE, INTRINSIC :: ISO_C_Binding
        USE DRC_Types, ONLY : LocalVariables, ControlParameters
    
        IMPLICIT NONE
    
        ! Inputs
        TYPE(ControlParameters), INTENT(IN)     :: CntrPar
        TYPE(LocalVariables), INTENT(IN)        :: LocalVar
        INTEGER(4), INTENT(IN)                  :: size_avcMSG
        REAL(C_FLOAT), INTENT(IN)               :: avrSWAP(*)          ! The swap array, used to pass data to, and receive data from, the DLL controller.
        
        ! Outputs
        INTEGER(C_INT), INTENT(OUT)             :: aviFAIL             ! A flag used to indicate the success of this DLL call set as follows: 0 if the DLL call was successful, >0 if the DLL call was successful but cMessage should be issued as a warning messsage, <0 if the DLL call was unsuccessful or for any other reason the simulation is to be stopped at this point with cMessage as the error message.
        CHARACTER(size_avcMSG-1), INTENT(OUT)   :: ErrMsg              ! a Fortran version of the C string argument (not considered an array here) [subtract 1 for the C null-character]
        
        ! Local
        
        !..............................................................................................................................
        ! Check validity of input parameters:
        !..............................................................................................................................
        
        IF ((CntrPar%F_LPFType > 2.0) .OR. (CntrPar%F_LPFType < 1.0)) THEN
            aviFAIL = -1
            ErrMsg  = 'F_LPFType must be 1 or 2.'
        ENDIF
        
        IF ((CntrPar%F_LPFDamping > 1.0) .OR. (CntrPar%F_LPFDamping < 0.0)) THEN
            aviFAIL = -1
            ErrMsg  = 'Filter damping coefficient must be between [0, 1]'
        ENDIF
        
        IF (CntrPar%IPC_CornerFreqAct < 0.0) THEN
            aviFAIL = -1
            ErrMsg  = 'Corner frequency of IPC actuator model must be positive, or set to 0 to disable.'
        ENDIF
        
        IF (CntrPar%F_LPFCornerFreq <= 0.0) THEN
            aviFAIL = -1
            ErrMsg  = 'CornerFreq must be greater than zero.'
        ENDIF
        
        IF ((CntrPar%IPC_ControlMode > 0) .AND. (CntrPar%Y_ControlMode > 1)) THEN
            aviFAIL = -1
            ErrMsg  = 'IPC control for load reductions and yaw-by-IPC cannot be activated simultaneously'
        ENDIF
        
        IF (LocalVar%DT <= 0.0) THEN
            aviFAIL = -1
            ErrMsg  = 'DT must be greater than zero.'
        ENDIF
        
        IF (CntrPar%VS_MaxRat <= 0.0) THEN
            aviFAIL =  -1
            ErrMsg  = 'VS_MaxRat must be greater than zero.'
        ENDIF
        
        IF (CntrPar%VS_RtTq < 0.0) THEN
            aviFAIL = -1
            ErrMsg  = 'VS_RtTq must not be negative.'
        ENDIF
        
        IF (CntrPar%VS_Rgn2K < 0.0) THEN
            aviFAIL = -1
            ErrMsg  = 'VS_Rgn2K must not be negative.'
        ENDIF
        
        IF (CntrPar%VS_MaxTq < CntrPar%VS_RtTq) THEN
            aviFAIL = -1
            ErrMsg  = 'VS_RtTq must not be greater than VS_MaxTq.'
        ENDIF
        
        IF (CntrPar%VS_KP(1) > 0.0) THEN
            aviFAIL = -1
            ErrMsg  = 'VS_KP must be greater than zero.'
        ENDIF
        
        IF (CntrPar%VS_KI(1) > 0.0) THEN
            aviFAIL = -1
            ErrMsg  = 'VS_KI must be greater than zero.'
        ENDIF
        
        IF (CntrPar%PC_RefSpd <= 0.0) THEN
            aviFAIL = -1
            ErrMsg  = 'PC_RefSpd must be greater than zero.'
        ENDIF
        
        IF (CntrPar%PC_MaxRat <= 0.0) THEN
            aviFAIL = -1
            ErrMsg  = 'PC_MaxRat must be greater than zero.'
        ENDIF
        
        IF (CntrPar%PC_MinPit >= CntrPar%PC_MaxPit)  THEN
            aviFAIL = -1
            ErrMsg  = 'PC_MinPit must be less than PC_MaxPit.'
        ENDIF
        
        IF (CntrPar%IPC_KI(1) < 0.0)  THEN
            aviFAIL = -1
            ErrMsg  = 'IPC_KI(1) must be zero or greater than zero.'
        ENDIF
        
        IF (CntrPar%IPC_KI(2) < 0.0)  THEN
            aviFAIL = -1
            ErrMsg  = 'IPC_KI(2) must be zero or greater than zero.'
        ENDIF
        
        IF (CntrPar%Y_IPC_omegaLP <= 0.0)  THEN
            aviFAIL = -1
            ErrMsg  = 'Y_IPC_omegaLP must be greater than zero.'
        ENDIF
        
        IF (CntrPar%Y_IPC_zetaLP <= 0.0)  THEN
            aviFAIL = -1
            ErrMsg  = 'Y_IPC_zetaLP must be greater than zero.'
        ENDIF
        
        IF (CntrPar%Y_ErrThresh <= 0.0)  THEN
            aviFAIL = -1
            ErrMsg  = 'Y_ErrThresh must be greater than zero.'
        ENDIF
        
        IF (CntrPar%Y_Rate <= 0.0)  THEN
            aviFAIL = -1
            ErrMsg  = 'CntrPar%Y_Rate must be greater than zero.'
        ENDIF
        
        IF (CntrPar%Y_omegaLPFast <= 0.0)  THEN
            aviFAIL = -1
            ErrMsg  = 'Y_omegaLPFast must be greater than zero.'
        ENDIF
        
        IF (CntrPar%Y_omegaLPSlow <= 0.0)  THEN
            aviFAIL = -1
            ErrMsg  = 'Y_omegaLPSlow must be greater than zero.'
        ENDIF
        
        ! Abort if the user has not requested a pitch angle actuator (See Appendix A
        ! of Bladed User's Guide):
        IF (NINT(avrSWAP(10)) /= 0)  THEN ! .TRUE. if a pitch angle actuator hasn't been requested
            aviFAIL = -1
            ErrMsg  = 'Pitch angle actuator not requested.'
        ENDIF
        
        IF (NINT(avrSWAP(28)) == 0 .AND. ((CntrPar%IPC_ControlMode > 0) .OR. (CntrPar%Y_ControlMode > 1))) THEN
            aviFAIL = -1
            ErrMsg  = 'IPC enabled, but Ptch_Cntrl in ServoDyn has a value of 0. Set it to 1.'
        ENDIF
    END SUBROUTINE Assert
    
    SUBROUTINE SetParameters(avrSWAP, accINFILE, aviFAIL, ErrMsg, size_avcMSG, CntrPar, LocalVar, objInst, FilterVar, CntrVar)
        USE DRC_Types, ONLY : ControlParameters, LocalVariables, ObjectInstances, FilterVariables, ControllerVariables
        
        INTEGER(4), INTENT(IN) :: size_avcMSG
        TYPE(ControlParameters), INTENT(INOUT) :: CntrPar
        TYPE(LocalVariables), INTENT(INOUT) :: LocalVar
        TYPE(ObjectInstances), INTENT(INOUT) :: objInst
        TYPE(FilterVariables), INTENT(INOUT) :: FilterVar
        TYPE(ControllerVariables), INTENT(INOUT) :: CntrVar
        
        REAL(C_FLOAT), INTENT(INOUT) :: avrSWAP(*)          ! The swap array, used to pass data to, and receive data from, the DLL controller.
        CHARACTER(C_CHAR), INTENT(IN) :: accINFILE(NINT(avrSWAP(50)))        ! The name of the parameter input file
        INTEGER(C_INT), INTENT(OUT) :: aviFAIL              ! A flag used to indicate the success of this DLL call set as follows: 0 if the DLL call was successful, >0 if the DLL call was successful but cMessage should be issued as a warning messsage, <0 if the DLL call was unsuccessful or for any other reason the simulation is to be stopped at this point with cMessage as the error message.
        CHARACTER(size_avcMSG-1), INTENT(OUT) :: ErrMsg     ! a Fortran version of the C string argument (not considered an array here) [subtract 1 for the C null-character]
        INTEGER(4) :: K                                     ! Index used for looping through blades.
        INTEGER(4) :: fid = 98                              ! File ID for the controller input state file
        CHARACTER(LEN=256) :: fileName
        CHARACTER(LEN=256) :: directoryPath
        CHARACTER(LEN=256) :: FinalName
        CHARACTER(LEN=5)   :: file_id
        INTEGER(4) :: lastSlash                            ! Last occurrence of the directory separator in the filename
        INTEGER(4) :: i                                    ! Loop index
        
        ! Set aviFAIL to 0 in each iteration:
        aviFAIL = 0
        
        ! Initialize all filter instance counters at 1
        objInst%instLPF = 1
        objInst%instSecLPF = 1
        objInst%instHPF = 1
        objInst%instNotchSlopes = 1
        objInst%instNotch = 1
        objInst%instPI = 1
        
        ! Set unused outputs to zero (See Appendix A of Bladed User's Guide):
        avrSWAP(35) = 1.0 ! Generator contactor status: 1=main (high speed) variable-speed generator
        avrSWAP(36) = 0.0 ! Shaft brake status: 0=off
        avrSWAP(41) = 0.0 ! Demanded yaw actuator torque
        avrSWAP(46) = 0.0 ! Demanded pitch rate (Collective pitch)
        avrSWAP(55) = 0.0 ! Pitch override: 0=yes
        avrSWAP(56) = 0.0 ! Torque override: 0=yes
        avrSWAP(65) = 0.0 ! Number of variables returned for logging
        avrSWAP(72) = 0.0 ! Generator start-up resistance
        avrSWAP(79) = 0.0 ! Request for loads: 0=none
        avrSWAP(80) = 0.0 ! Variable slip current status
        avrSWAP(81) = 0.0 ! Variable slip current demand
        
        ! Read any External Controller Parameters specified in the User Interface
        !   and initialize variables:
        IF (LocalVar%iStatus == 0) THEN ! .TRUE. if we're on the first call to the DLL
            
            ! Inform users that we are using this user-defined routine:
            aviFAIL = 1
            ErrMsg = '                                                          '//NEW_LINE('A')// &
                     'Running the Delft Research Controller (DRC)               '//NEW_LINE('A')// &
                     'A wind turbine controller for use in the scientific field '//NEW_LINE('A')// &
                     'Written by S.P. Mulders, Jan-Willem van Wingerden         '//NEW_LINE('A')// &
                     'Delft University of Technology, The Netherlands           '//NEW_LINE('A')// &
                     'Visit our GitHub-page to contribute to this project:      '//NEW_LINE('A')// &
                     'https://github.com/TUDelft-DataDrivenControl              '
            
            CALL ReadControlParameterFileSub(avrSWAP, accINFILE, CntrPar)
            
            ! Initialize testValue (debugging variable)
            LocalVar%TestType = 0
        
            ! Initialize the SAVED variables:

            ! DO K = 1,LocalVar%NumBl
            LocalVar%PitCom = LocalVar%BlPitch ! This will ensure that the variable speed controller picks the correct control region and the pitch controller picks the correct gain on the first call
            ! END DO
            
            LocalVar%Y_AccErr = 0.0  ! This will ensure that the accumulated yaw error starts at zero
            LocalVar%Y_YawEndT = -1.0 ! This will ensure that the initial yaw end time is lower than the actual time to prevent initial yawing
            
            ! Wind speed estimator initialization, we always assume an initial wind speed of 10 m/s
            LocalVar%WE_Vw = 10
            LocalVar%WE_VwI = LocalVar%WE_Vw - CntrPar%WE_Gamma*LocalVar%RotSpeed
            
            ! Check validity of input parameters:
            CALL Assert(LocalVar, CntrPar, avrSWAP, aviFAIL, ErrMsg, size_avcMSG)

        ELSEIF (LocalVar%iStatus == 2) THEN
            
            CALL ReadControlParameterFileSub(avrSWAP, accINFILE, CntrPar)

            ! Initialize filename and directoryPath with spaces
            fileName = ' '
            directoryPath = ' '
            
            ! Convert the character array to a scalar string
            DO i = 1, MIN(LEN(fileName), SIZE(accINFILE))
                fileName(i:i) = accINFILE(i)
            END DO
            
            ! Find the last occurrence of the directory separator
            lastSlash = INDEX(fileName, '/', BACK=.TRUE.)
            IF (lastSlash == 0) THEN
                lastSlash = INDEX(fileName, '\', BACK=.TRUE.)
            END IF
            
            ! Extract the directory path
            IF (lastSlash > 0) THEN
                directoryPath = fileName(1:lastSlash)
            ELSE
                directoryPath = './' ! Default to current directory if no separator is found
            END IF

            write(file_id, '(I3.3)') NINT(avrSWAP(120))
            write(FinalName,'(6A)')  TRIM(directoryPath), '/', TRIM(CntrPar%RestartFile), '_turb_', TRIM(file_id), '.in'

            ! Read controller state input file
            open(unit=fid, file=trim(FinalName), status='old', form='unformatted', access='stream')
                read(fid) LocalVar, objInst, FilterVar, CntrVar
            close(fid)      
            
            ! Check validity of input parameters:
            CALL Assert(LocalVar, CntrPar, avrSWAP, aviFAIL, ErrMsg, size_avcMSG)
        ENDIF
    END SUBROUTINE SetParameters
END MODULE ReadSetParameters