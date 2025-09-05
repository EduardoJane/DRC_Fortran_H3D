#================================================================================#
# This makefile created by B. Jonkman on 2-Apr-2013,                             #
# adapted from Crunch (M. Buhl on 25-Jan-2013).                                  #
# (c) 2013 National Renewable Energy Laboratory                                  #
#                                                                                #
# This makefile has been tested on Windows 7 with gfortran.                      #
# This makefile works with mingw32-make.exe.                                     #
#                                                                                #
# It was designed to be used with discon.f90 for the Bladed DLL Interface        #
#================================================================================#

   # 32-bit or 64-bit?
# BITS = 32
BITS = 64
   # Location of source files for the DLL.
   # You may need to change these for your DLL.

DLL_DIR     = .

SOURCE_FILE = Constants.f90 DRC_Types.f90 Functions.f90 Filters.f90 ReadSetParameters.f90 Controllers.f90 DISCON.f90

########################################
## User defined compilation flags 
########################################
GNU_HPC_FLAGS= -cpp -fallow-argument-mismatch -std=legacy -ffree-line-length-0 -O3 -ftree-vectorize \
		   -ftree-vectorizer-verbose=0 -fbackslash -D_has_Quad -march=znver3 -funroll-loops
GNU_RELEASE_FLAGS= -cpp -ffree-line-length-0 -O3 -ftree-vectorize  -ftree-vectorizer-verbose=0 \
		   -fbackslash -D_has_Quad
GNU_DEBUG_FLAGS=  -ffree-line-length-0 -cpp -O0 -g -fcheck=all -Wno-unused-variable \
		  -fbacktrace -fbounds-check -finit-real=snan -Wall \
		  -ffpe-trap=invalid,zero,overflow -fcheck-array-temporaries \
		  -fbackslash -fcheck=pointer -D_has_Quad
INTEL_HPC_FLAGS= -cpp -Ofast -msse4.2 -axAVX,CORE-AVX2 -assume bscc -D_has_Quad  -fp-model=fast=2 -ipo
INTEL_RELEASE_FLAGS= -cpp -O3 -assume bscc -D_has_Quad -fp-model precise
INTEL_DEBUG_FLAGS= -cpp -O0 -g -warn all -fpscomp logicals -traceback \
		   -check bounds -check uninit -ftrapuv -debug all \
		   -ftrace=full -ftz -fp-model precise \
		   -fp-speculation=off -assume bscc -D_has_Quad

##########################################
## Select compiler 
##########################################
ifeq ($(COMPILER),gfortran)
  FC:=gfortran
else ifeq ($(COMPILER),gnu)
  FC:=gfortran
else ifeq ($(COMPILER),GNU)
  FC:=gfortran
else ifeq ($(COMPILER),ifort)
  FC:=ifort
else ifeq ($(COMPILER),intel)
  FC:=ifort
else ifeq ($(COMPILER),INTEL)
  FC:=ifort
else
   ifneq ($(MAKECMDGOALS),clean)
      $(info COMPILER value unspecified or not recognised. Using default gfortran)
   endif
   FC:=gfortran
endif

##########################################
## Select mode
##########################################
ifeq ($(MODE),RELEASE)
     ifeq ($(FC),gfortran)
         FFLAGS=$(GNU_RELEASE_FLAGS)
     else ifeq ($(FC),ifort)
         FFLAGS=$(INTEL_RELEASE_FLAGS)
     endif
else ifeq ($(MODE),HPC)
     ifeq ($(FC),gfortran)
         FFLAGS=$(GNU_HPC_FLAGS)
     else ifeq ($(FC),ifort)
         FFLAGS=$(INTEL_HPC_FLAGS)
     endif	 
else ifeq ($(MODE),DEBUG)
     ifeq ($(FC),gfortran)
         FFLAGS=$(GNU_DEBUG_FLAGS)
     else ifeq ($(FC),ifort)
         FFLAGS=$(INTEL_DEBUG_FLAGS)
     endif
else
   ifneq ($(MAKECMDGOALS),clean)
      $(info MODE argument unspecified or not recognised. Using RELEASE compiler flags)
   endif
   ifeq ($(FC),gfortran)
       FFLAGS=$(GNU_RELEASE_FLAGS)
   else ifeq ($(FC),ifort)
       FFLAGS=$(INTEL_RELEASE_FLAGS)
   endif
endif

LDFLAGS = -shared -O2 -m$(BITS) -fbacktrace

# other useful gfortran options:
#  -fdefault-real-8 -fcheck=bounds -std=f2003 -O0 -v
# -Wl,--out-implib,DISCON.lib


   # Destination and RootName for executable

OUTPUT_NAME = DISCON
DEST_DIR    = ../DISCON

   #==========================================================#
   # You should not need to change anything beyond this point #
   #==========================================================#

   # System-specific settings.

ifeq ($(OS),Windows_NT)
      # Windows
   DEL_CMD   = del
   EXE_EXT   = _gwin$(BITS).dll
   INTER_DIR = Obj_win$(BITS)
   MD_CMD    = @mkdir
   OBJ_EXT   = .obj
   PATH_SEP  = \\
   SYS_FILE  = SysGnuWin
else

   FFLAGS := $(FFLAGS) -DIMPLICIT_DLLEXPORT

      # Linux
   DEL_CMD   = rm -f
   EXE_EXT   = _glin$(BITS).so
   INTER_DIR = Obj_lin$(BITS)
   MD_CMD    = @mkdir -p
   OBJ_EXT   = .o
   PATH_SEP  = /
   SYS_FILE  = SysGnuLinux

   UNAME := $(shell uname -s)

   ifneq ($(UNAME), Darwin)
      FFLAGS     := $(FFLAGS) -fPIC
      LDFLAGS    := $(LDFLAGS) -fPIC
   endif
endif


   # Source files (by module)

DLL_SOURCES =          \
	$(SOURCE_FILE)

vpath %.f90 $(DLL_DIR)
vpath %.mod $(INTER_DIR)
vpath %.obj $(INTER_DIR)

ALL_OBJS      = $(DLL_SOURCES:.f90=.obj)

   # Rule to do everything.

all:     default
default:  $(INTER_DIR) $(DEST_DIR)/$(OUTPUT_NAME)$(EXE_EXT)

   # General rules for compliling the files.

%.obj: %.f90
	$(FC)  $(FFLAGS) -c $< -o $(INTER_DIR)/$@ -J $(INTER_DIR) -B $(INTER_DIR)


   #  Dependency rules.
#$(SYS_FILE).obj:                   $(PREC).obj


   # Make sure the destination directory for the intermediate files exist.

$(INTER_DIR):
	$(MD_CMD) $(INTER_DIR)


   # For linking DLL.

$(DEST_DIR)/$(OUTPUT_NAME)$(EXE_EXT): $(ALL_OBJS) | $(INTER_DIR)
	$(FC) $(LDFLAGS) -I $(INTER_DIR) -o $(DEST_DIR)/$(OUTPUT_NAME)$(EXE_EXT) \
	$(foreach src, $(ALL_OBJS), $(addprefix $(INTER_DIR)/,$(src)))

   # Cleanup afterwards.

clean:
	$(DEL_CMD) $(INTER_DIR)$(PATH_SEP)*.mod $(INTER_DIR)$(PATH_SEP)*.obj

#gfortran -shared -O2 -m32 -fbacktrace -ffree-line-length-none -x f95-cpp-input -C ../CertTest/5MW_Baseline/ServoData/Source/DISCON.f90