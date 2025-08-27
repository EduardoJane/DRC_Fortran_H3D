!=======================================================================
! SUBROUTINE DISCON(avrSWAP, from_SC, to_SC, aviFAIL, accINFILE, avcOUTNAME, avcMSG) BIND (C, NAME='DISCON')
SUBROUTINE DISCON(avrSWAP, aviFAIL, accINFILE, avcOUTNAME, avcMSG) BIND (C, NAME='DISCON')
! DO NOT REMOVE or MODIFY LINES starting with "!DEC$" or "!GCC$"
! !DEC$ specifies attributes for IVF and !GCC$ specifies attributes for gfortran
!DEC$ ATTRIBUTES DLLEXPORT :: DISCON

USE, INTRINSIC  :: ISO_C_Binding
USE             :: DRC_Types
USE             :: ReadSetParameters
USE             :: Controllers
USE             :: Constants
USE             :: Filters

IMPLICIT NONE
#ifndef IMPLICIT_DLLEXPORT
!GCC$ ATTRIBUTES DLLEXPORT :: DISCON
#endif

!------------------------------------------------------------------------------------------------------------------------------
! Variable declaration and initialization
!------------------------------------------------------------------------------------------------------------------------------

! Passed Variables:
!REAL(C_FLOAT), INTENT(IN)      :: from_SC(*)       ! DATA from the super controller
!REAL(C_FLOAT), INTENT(INOUT)   :: to_SC(*)         ! DATA to the super controller

REAL(C_FLOAT), INTENT(INOUT)            :: avrSWAP(*)                       ! The swap array, used to pass data to, and receive data from, the DLL controller.
INTEGER(C_INT), INTENT(INOUT)           :: aviFAIL                          ! A flag used to indicate the success of this DLL call set as follows: 0 if the DLL call was successful, >0 if the DLL call was successful but cMessage should be issued as a warning messsage, <0 if the DLL call was unsuccessful or for any other reason the simulation is to be stopped at this point with cMessage as the error message.
CHARACTER(KIND=C_CHAR), INTENT(IN)      :: accINFILE(NINT(avrSWAP(50)))     ! The name of the parameter input file
CHARACTER(KIND=C_CHAR), INTENT(IN)      :: avcOUTNAME(NINT(avrSWAP(51)))    ! OUTNAME (Simulation RootName)
CHARACTER(KIND=C_CHAR), INTENT(INOUT)   :: avcMSG(NINT(avrSWAP(49)))        ! MESSAGE (Message from DLL to simulation code [ErrMsg])  The message which will be displayed by the calling program if aviFAIL <> 0.
CHARACTER(SIZE(avcOUTNAME)-1)           :: RootName                         ! a Fortran version of the input C string (not considered an array here)    [subtract 1 for the C null-character]
CHARACTER(SIZE(avcMSG)-1)               :: ErrMsg                           ! a Fortran version of the C string argument (not considered an array here) [subtract 1 for the C null-character]
CHARACTER(LEN=256)                      :: fileName                         ! Filename for the controller input parameters file
CHARACTER(LEN=256)                      :: directoryPath                    ! Directory path for the controller input parameters file
CHARACTER(LEN=256)                      :: FinalName                        ! Final filename for the controller output state file
CHARACTER(LEN=512)                      :: file_id                          ! File ID used in controller output state file name
INTEGER(4)                              :: fid = 99                         ! File ID for the controller output state file
INTEGER(4)                              :: lastSlash                        ! Last occurrence of the directory separator in the filename
INTEGER(4)                              :: i                                ! Loop index

TYPE(ControlParameters), SAVE         :: CntrPar
TYPE(LocalVariables), SAVE            :: LocalVar
TYPE(ObjectInstances), SAVE           :: objInst
TYPE(FilterVariables), SAVE           :: FilterVar
TYPE(ControllerVariables), SAVE       :: CntrVar

!------------------------------------------------------------------------------------------------------------------------------
! Main control calculations
!------------------------------------------------------------------------------------------------------------------------------
! Read avrSWAP array into derived types/variables
CALL ReadAvrSWAP(avrSWAP, LocalVar)

! Write controller state (derived types) to output file for future restart
IF (LocalVar%iStatus == 3) THEN
    ! Initialize filename and directoryPath with spaces
    fileName = ' '
    directoryPath = ' '
    
    ! Convert the character array to a scalar string
    DO i = 1, MIN(LEN(filename), SIZE(accINFILE))
        filename(i:i) = accINFILE(i)
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
    
    if (avrSWAP(121) .lt. 0.0) then
        write(FinalName,'(4A)')  TRIM(directoryPath), 'DISCON_state_turb_', TRIM(file_id), '.in'
    else
        write(FinalName,'(2A,I10.10,3A)')  TRIM(directoryPath), 'DISCON_state_', NINT(avrSWAP(121)), '_turb_', TRIM(file_id), '.in'
    end if

    open(unit=fid, file=trim(FinalName), status='replace', form='unformatted', access='stream')
        write(fid) LocalVar, objInst, FilterVar, CntrVar
    close(fid)
    
    RETURN
ELSE
    CALL SetParameters(avrSWAP, accINFILE, aviFAIL, ErrMsg, SIZE(avcMSG), CntrPar, LocalVar, objInst, FilterVar, CntrVar)
    CALL PreFilterMeasuredSignals(CntrPar, LocalVar, objInst, FilterVar)

    IF ((LocalVar%iStatus >= 0) .AND. (aviFAIL >= 0))  THEN  ! Only compute control calculations if no error has occurred and we are not on the last time step
        CALL ComputeVariablesSetpoints(CntrPar, LocalVar)
        
        CALL StateMachine(CntrPar, LocalVar)
        CALL WindSpeedEstimator(LocalVar, CntrPar)
        
        CALL VariableSpeedControl(avrSWAP, CntrPar, LocalVar, objInst, CntrVar)
        CALL PitchControl(avrSWAP, CntrPar, LocalVar, objInst, FilterVar, CntrVar)
        CALL YawRateControl(avrSWAP, CntrPar, LocalVar, objInst, FilterVar)

        CALL Debug(LocalVar, CntrPar, avrSWAP, RootName, SIZE(avcOUTNAME))
    END IF

    avcMSG = TRANSFER(TRIM(ErrMsg)//C_NULL_CHAR, avcMSG, SIZE(avcMSG))
    RETURN
END IF
END SUBROUTINE DISCON
