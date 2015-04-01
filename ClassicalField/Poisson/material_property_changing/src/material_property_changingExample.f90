!> \file
!> \author Chris Bradley
!> \brief This is an example program to solve a constant source Poisson equation using OpenCMISS calls.
!>
!> \section LICENSE
!>
!> Version: MPL 1.1/GPL 2.0/LGPL 2.1
!>
!> The contents of this file are subject to the Mozilla Public License
!> Version 1.1 (the "License"); you may not use this file except in
!> compliance with the License. You may obtain a copy of the License at
!> http://www.mozilla.org/MPL/
!>
!> Software distributed under the License is distributed on an "AS IS"
!> basis, WITHOUT WARRANTY OF ANY KIND, either express or implied. See the
!> License for the specific language governing rights and limitations
!> under the License.
!>
!> The Original Code is OpenCMISS
!>
!> The Initial Developer of the Original Code is University of Auckland,
!> Auckland, New Zealand and University of Oxford, Oxford, United
!> Kingdom. Portions created by the University of Auckland and University
!> of Oxford are Copyright (C) 2007 by the University of Auckland and
!> the University of Oxford. All Rights Reserved.
!>
!> Contributor(s):
!>
!> Alternatively, the contents of this file may be used under the terms of
!> either the GNU General Public License Version 2 or later (the "GPL"), or
!> the GNU Lesser General Public License Version 2.1 or later (the "LGPL"),
!> in which case the provisions of the GPL or the LGPL are applicable instead
!> of those above. If you wish to allow use of your version of this file only
!> under the terms of either the GPL or the LGPL, and not to allow others to
!> use your version of this file under the terms of the MPL, indicate your
!> decision by deleting the provisions above and replace them with the notice
!> and other provisions required by the GPL or the LGPL. If you do not delete
!> the provisions above, a recipient may use your version of this file under
!> the terms of any one of the MPL, the GPL or the LGPL.
!>

!> \example ClassicalField/Poisson/42Master/src/42MasterExample.f90
!! Example program to solve a constant source Poisson equation using openCMISS calls.
!! \htmlinclude ClassicalField/Poisson/42Master/history.html
!<

!> Main program
PROGRAM LINEARPOISSONEXAMPLE

  USE OPENCMISS
  USE MPI

#ifdef WIN32
  USE IFQWIN
#endif

  IMPLICIT NONE

  !Test program parameters

  REAL(CMISSDP), PARAMETER :: HEIGHT=4.0_CMISSDP
  REAL(CMISSDP), PARAMETER :: WIDTH=6.0_CMISSDP
  REAL(CMISSDP), PARAMETER :: LENGTH=0.5_CMISSDP

  INTEGER(CMISSIntg), PARAMETER :: CoordinateSystemUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: RegionUserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: BasisUserNumber=3
  INTEGER(CMISSIntg), PARAMETER :: GeneratedMeshUserNumber=4
  INTEGER(CMISSIntg), PARAMETER :: MeshUserNumber=5
  INTEGER(CMISSIntg), PARAMETER :: DecompositionUserNumber=6
  INTEGER(CMISSIntg), PARAMETER :: GeometricFieldUserNumber=7
  INTEGER(CMISSIntg), PARAMETER :: DependentFieldUserNumber=8
  INTEGER(CMISSIntg), PARAMETER :: MaterialsFieldUserNumber=9
  INTEGER(CMISSIntg), PARAMETER :: SourceFieldUserNumber=10
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetUserNumber=11
  INTEGER(CMISSIntg), PARAMETER :: ProblemUserNumber=12
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetFieldUserNumber=13

  !Program variables

  INTEGER(CMISSIntg) :: NUMBER_DIMENSIONS,INTERPOLATION_TYPE,NUMBER_OF_GAUSS_XI
  INTEGER(CMISSIntg) :: NUMBER_GLOBAL_X_ELEMENTS,NUMBER_GLOBAL_Y_ELEMENTS,NUMBER_GLOBAL_Z_ELEMENTS
  INTEGER(CMISSIntg) :: component_idx,my_material_components
  INTEGER(CMISSIntg) :: i
  CHARACTER(LEN=255) :: Filename
 ! LOGICAL :: EXPORT_FIELD
  !CMISS variables

  TYPE(CMISSBasisType) :: Basis
  TYPE(CMISSBoundaryConditionsType) :: BoundaryConditions
  TYPE(CMISSCoordinateSystemType) :: CoordinateSystem,WorldCoordinateSystem
  TYPE(CMISSDecompositionType) :: Decomposition
  TYPE(CMISSEquationsType) :: Equations
  TYPE(CMISSEquationsSetType) :: EquationsSet
  TYPE(CMISSFieldType) :: GeometricField,DependentField,MaterialsField,SourceField
  TYPE(CMISSFieldsType) :: Fields
  TYPE(CMISSGeneratedMeshType) :: GeneratedMesh
  TYPE(CMISSMeshType) :: Mesh
  TYPE(CMISSProblemType) :: Problem
  TYPE(CMISSRegionType) :: Region,WorldRegion
  TYPE(CMISSSolverType) :: Solver,LinearSolver
  TYPE(CMISSSolverEquationsType) :: SolverEquations
  TYPE(CMISSFieldType) :: EquationsSetField
  TYPE(CMISSControlLoopType) :: ControlLoop
!==========================================================================================================
  !Generic CMISS variables

  INTEGER(CMISSIntg) :: EquationsSetIndex
  INTEGER(CMISSIntg) :: Err
  INTEGER(CMISSIntg) :: NumberOfComputationalNodes,ComputationalNodeNumber
  INTEGER(CMISSIntg), PARAMETER :: ControlLoopNode=0

#ifdef WIN32
  !Quickwin type
  LOGICAL :: QUICKWIN_STATUS=.FALSE.
  TYPE(WINDOWCONFIG) :: QUICKWIN_WINDOW_CONFIG

  !Initialise QuickWin
  QUICKWIN_WINDOW_CONFIG%TITLE="General Output" !Window title
  QUICKWIN_WINDOW_CONFIG%NUMTEXTROWS=-1 !Max possible number of rows
  QUICKWIN_WINDOW_CONFIG%MODE=QWIN$SCROLLDOWN
  !Set the window parameters
  QUICKWIN_STATUS=SETWINDOWCONFIG(QUICKWIN_WINDOW_CONFIG)
  !If attempt fails set with system estimated values
  IF(.NOT.QUICKWIN_STATUS) QUICKWIN_STATUS=SETWINDOWCONFIG(QUICKWIN_WINDOW_CONFIG)
#endif
!==========================================================================================================
                    ! Retinal modeling parameters
!==========================================================================================================

!======= Retinal layers conductivities 
  
  REAL(CMISSDP), PARAMETER :: Conductivity_Vitr = 1.25_CMISSDP
  REAL(CMISSDP), PARAMETER :: Conductivity_IPL = 0.0281_CMISSDP
  REAL(CMISSDP), PARAMETER :: Conductivity_OPL = 0.0151_CMISSDP
  REAL(CMISSDP), PARAMETER :: Conductivity_ONL = 0.0124_CMISSDP
  REAL(CMISSDP), PARAMETER :: Conductivity_RPE = 0.0483_CMISSDP
  REAL(CMISSDP), PARAMETER :: Conductivity_Electrode = 1.0_CMISSDP
  INTEGER(CMISSIntg),DIMENSION(8) :: injecting_electrode,returning_electrode ! 8 stands for the number of nodes

!======= the amount of injected current by the electrodes
 
  REAL(CMISSDP), PARAMETER :: injecting = 5.0_CMISSDP
  REAL(CMISSDP), PARAMETER :: returning = -5.0_CMISSDP

!======= the nodes of injecting and returning electrode
  
  injecting_electrode=(/1283,1284,1285,1272,1261,1262,1263,1274/)
  returning_electrode=(/1257,1258,1259,1268,1270,1279,1280,1281/)

!==========================================================================================================\
                 ! Geometry parameters
!==========================================================================================================

  NUMBER_DIMENSIONS = 3
  NUMBER_GLOBAL_X_ELEMENTS = 10
  NUMBER_GLOBAL_Y_ELEMENTS = 10
  NUMBER_GLOBAL_Z_ELEMENTS = 10
  INTERPOLATION_TYPE = 1

!==========================================================================================================

  !Intialise OpenCMISS
  CALL CMISSInitialise(WorldCoordinateSystem,WorldRegion,Err)

  !Trap all errors
  CALL CMISSErrorHandlingModeSet(CMISS_ERRORS_TRAP_ERROR,Err)
  CALL CMISSRandomSeedsSet(9999,Err)
  
  CALL CMISSDiagnosticsSetOn(CMISS_IN_DIAG_TYPE,[1,2,3,4,5],"Diagnostics",["DOMAIN_MAPPINGS_LOCAL_FROM_GLOBAL_CALCULATE"],Err)

  WRITE(Filename,'(A,"_",I0,"x",I0,"x",I0,"_",I0)') "Poisson",NUMBER_GLOBAL_X_ELEMENTS,NUMBER_GLOBAL_Y_ELEMENTS, &
    & NUMBER_GLOBAL_Z_ELEMENTS,INTERPOLATION_TYPE
  
  CALL CMISSOutputSetOn(Filename,Err)

  !CALL CMISSOutputSetOn("poisson",Err)
  !Output to a file
  !CALL CMISSOutputSetOn("Poisson",Err)

  !Get the computational nodes information
  CALL CMISSComputationalNumberOfNodesGet(NumberOfComputationalNodes,Err)
  CALL CMISSComputationalNodeNumberGet(ComputationalNodeNumber,Err)

  !Start the creation of a new RC coordinate system
  CALL CMISSCoordinateSystem_Initialise(CoordinateSystem,Err)
  CALL CMISSCoordinateSystem_CreateStart(CoordinateSystemUserNumber,CoordinateSystem,Err)
  CALL CMISSCoordinateSystem_DimensionSet(CoordinateSystem,NUMBER_DIMENSIONS,Err)
  CALL CMISSCoordinateSystem_CreateFinish(CoordinateSystem,Err)
!==========================================================================================================
  !Start the creation of the region
!==========================================================================================================
  CALL CMISSRegion_Initialise(Region,Err)
  CALL CMISSRegion_CreateStart(RegionUserNumber,WorldRegion,Region,Err)
  CALL CMISSRegion_LabelSet(Region,"PoissonRegion",Err)
  !Set the regions coordinate system to the 2D RC coordinate system that we have created
  CALL CMISSRegion_CoordinateSystemSet(Region,CoordinateSystem,Err)
  CALL CMISSRegion_CreateFinish(Region,Err)
!==========================================================================================================
  !Start the creation of a basis (default is trilinear lagrange)
!==========================================================================================================
  CALL CMISSBasis_Initialise(Basis,Err)
  CALL CMISSBasis_CreateStart(BasisUserNumber,Basis,Err)
  CALL CMISSBasis_NumberOfXiSet(Basis,NUMBER_DIMENSIONS,Err)
  CALL CMISSBasis_TypeSet(Basis,CMISS_BASIS_LAGRANGE_HERMITE_TP_TYPE,Err)
  NUMBER_OF_GAUSS_XI=2
  !Set the basis to be a tri-interpolation basis
  CALL CMISSBasis_NumberOfXiSet(Basis,3,Err)
  CALL CMISSBasis_InterpolationXiSet(Basis,[INTERPOLATION_TYPE,INTERPOLATION_TYPE,INTERPOLATION_TYPE],Err)
  IF(NUMBER_OF_GAUSS_XI>0) THEN
    CALL CMISSBasis_QuadratureNumberOfGaussXiSet(Basis,[NUMBER_OF_GAUSS_XI,NUMBER_OF_GAUSS_XI,NUMBER_OF_GAUSS_XI],Err)
  ENDIF
  !Finish the creation of the basis
  CALL CMISSBasis_CreateFinish(Basis,Err)
!==========================================================================================================
  !Start the creation of a generated mesh in the region
!==========================================================================================================
  CALL CMISSGeneratedMesh_Initialise(GeneratedMesh,Err)
  CALL CMISSGeneratedMesh_CreateStart(GeneratedMeshUserNumber,Region,GeneratedMesh,Err)
  !Set up a regular x*y*z mesh
  CALL CMISSGeneratedMesh_TypeSet(GeneratedMesh,CMISS_GENERATED_MESH_REGULAR_MESH_TYPE,Err)
  !Set the default basis
  CALL CMISSGeneratedMesh_BasisSet(GeneratedMesh,Basis,Err)   
  !Define the mesh on the region
  CALL CMISSGeneratedMesh_ExtentSet(GeneratedMesh,[WIDTH,HEIGHT,LENGTH],Err)
  CALL CMISSGeneratedMesh_NumberOfElementsSet(GeneratedMesh,[NUMBER_GLOBAL_X_ELEMENTS,NUMBER_GLOBAL_Y_ELEMENTS, &
    & NUMBER_GLOBAL_Z_ELEMENTS],Err)
  !Finish the creation of a generated mesh in the region
  CALL CMISSMesh_Initialise(Mesh,Err)
  CALL CMISSGeneratedMesh_CreateFinish(GeneratedMesh,MeshUserNumber,Mesh,Err)

!==========================================================================================================
  !Create a decomposition
!==========================================================================================================

  CALL CMISSDecomposition_Initialise(Decomposition,Err)
  CALL CMISSDecomposition_CreateStart(DecompositionUserNumber,Mesh,Decomposition,Err)
  !Set the decomposition to be a general decomposition with the specified number of domains
  CALL CMISSDecomposition_TypeSet(Decomposition,CMISS_DECOMPOSITION_CALCULATED_TYPE,Err)
  CALL CMISSDecomposition_NumberOfDomainsSet(Decomposition,NumberOfComputationalNodes,Err)
  CALL CMISSDecomposition_CreateFinish(Decomposition,Err)

  !Start to create a default (geometric) field on the region
  CALL CMISSField_Initialise(GeometricField,Err)
  CALL CMISSField_CreateStart(GeometricFieldUserNumber,Region,GeometricField,Err)
  CALL CMISSField_MeshDecompositionSet(GeometricField,Decomposition,Err)
  DO component_idx=1,NUMBER_DIMENSIONS
    CALL CMISSField_ComponentMeshComponentSet(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,component_idx,1,Err)
  ENDDO
  CALL CMISSField_CreateFinish(GeometricField,Err)
 
  !Update the geometric field parameters
  CALL CMISSGeneratedMesh_GeometricParametersCalculate(GeneratedMesh,GeometricField,Err)
!==========================================================================================================
!                        !Create the equations_set
!==========================================================================================================

  CALL CMISSEquationsSet_Initialise(EquationsSet,Err)
  CALL CMISSField_Initialise(EquationsSetField,Err)
  CALL CMISSEquationsSet_CreateStart(EquationsSetUserNumber,Region,GeometricField,CMISS_EQUATIONS_SET_CLASSICAL_FIELD_CLASS, &
    & CMISS_EQUATIONS_SET_POISSON_EQUATION_TYPE,CMISS_EQUATIONS_SET_TRANSIENT_SOURCE_POISSON_SUBTYPE,EquationsSetFieldUserNumber, &
    & EquationsSetField,EquationsSet,Err)
  CALL CMISSEquationsSet_CreateFinish(EquationsSet,Err)

!==========================================================================================================
                      !Create the equations set dependent field variables
!==========================================================================================================

  CALL CMISSField_Initialise(DependentField,Err)
  CALL CMISSEquationsSet_DependentCreateStart(EquationsSet,DependentFieldUserNumber,DependentField,Err)
  CALL CMISSField_VariableLabelSet(DependentField,CMISS_FIELD_U_VARIABLE_TYPE,'Phi',Err)
  CALL CMISSEquationsSet_DependentCreateFinish(EquationsSet,Err)

  !Initialise the field values
  CALL CMISSField_ComponentValuesInitialise(DependentField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,0.0_CMISSDP, &
    & Err)
!==========================================================================================================
!                      !Create the equations set material field variables 
!==========================================================================================================

  my_material_components=7; ! see the poisson routine- the first 3 components are the diagoanl which we need 

  CALL CMISSField_Initialise(MaterialsField,Err)
  CALL CMISSEquationsSet_MaterialsCreateStart(EquationsSet,MaterialsFieldUserNumber,MaterialsField,Err)
  CALL CMISSField_VariableLabelSet(MaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,'Materials',Err)
  Do component_idx=1,my_material_components
    CALL CMISSField_ComponentInterpolationSet(MaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,component_idx, &
      & CMISS_FIELD_ELEMENT_BASED_INTERPOLATION,Err)
  ENDDO
  CALL CMISSEquationsSet_MaterialsCreateFinish(EquationsSet,Err)

!====== MATERIAL PROPERTIES OF RETINAL LAYERS ======

  Do component_idx=1,my_material_components
    IF (component_idx<4) THEN
      Do i=1,1000
        IF(i<201) THEN
            CALL CMISSField_ParameterSetUpdateElement(MaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,i, &
              & component_idx,Conductivity_Vitr,Err)
        ENDIF 
        IF(i>200 .AND. i<401) THEN
            CALL CMISSField_ParameterSetUpdateElement(MaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,i, &
              & component_idx,Conductivity_IPL,Err)
        ENDIF
        IF(i>400 .AND. i<601) THEN
            CALL CMISSField_ParameterSetUpdateElement(MaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,i, &
              & component_idx,Conductivity_OPL,Err)
        ENDIF
        IF(i>600 .AND. i<801) THEN
            CALL CMISSField_ParameterSetUpdateElement(MaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,i, &
              & component_idx,Conductivity_ONL,Err)
        ENDIF
        IF (i>800) THEN
            CALL CMISSField_ParameterSetUpdateElement(MaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,i, &
              & component_idx,Conductivity_RPE,Err)
        ENDIF
      ENDDO
    ENDIF
    IF (component_idx>3) THEN
      Do i=1,1000
        CALL CMISSField_ParameterSetUpdateElement(MaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,i, &
          & component_idx,0.0_CMISSDP,Err)  !==== the x12, x21 conductivities are zere+ the last component which is the Source term  
      ENDDO
    ENDIF
  ENDDO

  CALL CMISSField_ParameterSetUpdateStart(MaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,Err)
  CALL CMISSField_ParameterSetUpdateFinish(MaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,Err)

!==========================================================================================================
!                                   !Create the source field
!==========================================================================================================
  
  CALL CMISSField_Initialise(SourceField,Err)
  CALL CMISSEquationsSet_SourceCreateStart(EquationsSet,SourceFieldUserNumber,SourceField,Err)
  CALL CMISSField_VariableLabelSet(SourceField,CMISS_FIELD_U_VARIABLE_TYPE,'Source',Err)
  CALL CMISSEquationsSet_SourceCreateFinish(EquationsSet,Err)

!==========================================================================================================
!                            !Create the equations set equations/Problem
!==========================================================================================================
 
  CALL CMISSEquations_Initialise(Equations,Err)
  CALL CMISSEquationsSet_EquationsCreateStart(EquationsSet,Equations,Err)
  CALL CMISSEquations_SparsityTypeSet(Equations,CMISS_EQUATIONS_SPARSE_MATRICES,Err)
  CALL CMISSEquations_OutputTypeSet(Equations,CMISS_EQUATIONS_NO_OUTPUT,Err)
  CALL CMISSEquationsSet_EquationsCreateFinish(EquationsSet,Err)
!==========================================================================================================
                                 !Start the creation of a problem/ control loop
!==========================================================================================================

  CALL CMISSProblem_Initialise(Problem,Err)
  CALL CMISSControlLoop_Initialise(ControlLoop,Err)
  CALL CMISSProblem_CreateStart(ProblemUserNumber,Problem,Err)
  CALL CMISSProblem_SpecificationSet(Problem,CMISS_PROBLEM_CLASSICAL_FIELD_CLASS,CMISS_PROBLEM_POISSON_EQUATION_TYPE, &
    & CMISS_PROBLEM_TRANSIENT_SOURCE_POISSON_SUBTYPE,Err)
  !Finish the creation of a problem.
  CALL CMISSProblem_CreateFinish(Problem,Err)

  CALL CMISSProblem_ControlLoopCreateStart(Problem,Err)
  !Get the control loop
  CALL CMISSProblem_ControlLoopGet(Problem,CMISS_CONTROL_LOOP_NODE,ControlLoop,Err)
  !Set the times
  CALL CMISSControlLoop_TimesSet(ControlLoop,0.0_CMISSDP,3.0_CMISSDP,1.0_CMISSDP,Err)
  !Finish creating the problem control loop
  CALL CMISSProblem_ControlLoopCreateFinish(Problem,Err)

!==========================================================================================================
  !Start the creation of the problem solvers
!==========================================================================================================
  CALL CMISSSolver_Initialise(Solver,Err)
  CALL CMISSSolver_Initialise(LinearSolver,Err)
  CALL CMISSProblem_SolversCreateStart(Problem,Err)
  CALL CMISSProblem_SolverGet(Problem,CMISS_CONTROL_LOOP_NODE,1,Solver,Err)
  CALL CMISSSolver_OutputTypeSet(Solver,CMISS_SOLVER_PROGRESS_OUTPUT,Err)
  CALL CMISSSolver_DynamicLinearSolverGet(Solver,LinearSolver,Err)
  CALL CMISSSolver_LinearIterativeMaximumIterationsSet(LinearSolver,300,Err)
  !Finish the creation of the problem solver
  CALL CMISSProblem_SolversCreateFinish(Problem,Err)

  !Start the creation of the problem solver equations
  CALL CMISSSolver_Initialise(Solver,Err)
  CALL CMISSSolverEquations_Initialise(SolverEquations,Err)
  CALL CMISSProblem_SolverEquationsCreateStart(Problem,Err)
  CALL CMISSProblem_SolverGet(Problem,CMISS_CONTROL_LOOP_NODE,1,Solver,Err)
  CALL CMISSSolver_SolverEquationsGet(Solver,SolverEquations,Err)
  CALL CMISSSolverEquations_SparsityTypeSet(SolverEquations,CMISS_SOLVER_SPARSE_MATRICES,Err)
  CALL CMISSSolverEquations_EquationsSetAdd(SolverEquations,EquationsSet,EquationsSetIndex,Err)
  CALL CMISSProblem_SolverEquationsCreateFinish(Problem,Err)

!==========================================================================================================
!                               !Set up the boundary conditions
!==========================================================================================================
 
  CALL CMISSBoundaryConditions_Initialise(BoundaryConditions,Err)
  CALL CMISSSolverEquations_BoundaryConditionsCreateStart(SolverEquations,BoundaryConditions,Err)

!=== Drichlet at the center of the guard:
  CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,1,1269,1, &
    & CMISS_BOUNDARY_CONDITION_FIXED,0.0_CMISSDP,Err)
!====
  CALL CMISSSolverEquations_BoundaryConditionsCreateFinish(SolverEquations,Err)
 
!==========================================================================================================
!                                 !Solve the problem
!==========================================================================================================

  CALL CMISSProblem_Solve(Problem,Err)

  !EXPORT_FIELD=.TRUE.
  !IF(EXPORT_FIELD) THEN
  CALL CMISSFields_Initialise(Fields,Err)
  CALL CMISSFields_Create(Region,Fields,Err)
  CALL CMISSFields_NodesExport(Fields,"Poisson","FORTRAN",Err)
  CALL CMISSFields_ElementsExport(Fields,"Poisson","FORTRAN",Err)
  CALL CMISSFields_Finalise(Fields,Err)
  !ENDIF

  !Finialise CMISS
  CALL CMISSFinalise(Err)

  WRITE(*,'(A)') "Program successfully completed."

  STOP

CONTAINS

  SUBROUTINE HANDLE_ERROR(ERROR_STRING)

    CHARACTER(LEN=*), INTENT(IN) :: ERROR_STRING

    WRITE(*,'(">>ERROR: ",A)') ERROR_STRING(1:LEN_TRIM(ERROR_STRING))
    STOP

  END SUBROUTINE HANDLE_ERROR

END PROGRAM LINEARPOISSONEXAMPLE
