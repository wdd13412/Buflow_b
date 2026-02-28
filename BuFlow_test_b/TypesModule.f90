module TypesModule
	implicit none
  ! 统一定义64位整数和实数常量
	integer, parameter :: int64 = 8       ! 8字节整数:real64可精确存储 15-17 位十进制有效数字,指数范围是 10⁻³⁰⁸~10³⁰⁸.因此完全可存储4.06575815E-20
	integer, parameter :: real64 = 8      ! 8字节实数（双精度）
	character(len=256), allocatable :: lines_defor(:)	
  
!!!!!!!!!!!!!!!!Buflow
!    real(kind=8), allocatable :: first5_refs(:,:)  ! 形状(5,3)!!!!!!!!!!!!
	real(kind=8), allocatable :: temp_cellPts(:,:) 
	integer(kind=8), allocatable :: currentBoundary(:)
	real(kind=8), allocatable :: P_eps(:), gradP(:,:), P_matrix(:,:), temp_grad(:,:,:), faceVals(:,:)
	real(kind=8), allocatable :: sj(:), sjCount(:), rj(:), rjsjF(:,:)
	real(kind=8), allocatable :: eps2(:), eps4(:)
	real(kind=8), allocatable :: timeFluxes(:), surfaceAreas(:)
	real(kind=8), allocatable :: timeFluxess(:), surfaceAreass(:)
	real(kind=8), allocatable :: faceRhoT(:,:), cellRhoT(:,:)
	real(kind=8), allocatable :: faceVel(:), positionn(:)
	real(kind=8), allocatable :: eps22(:), eps44(:), diffusionFlux(:), unitFA(:), fD(:), farOwnerfD(:), farNeighbourfD(:), epss(:,:)
	real(kind=8), allocatable :: fDeltas(:,:), fDGrads(:,:,:)
	real(kind=8), allocatable :: faceVals_dd(:,:), faceVals_d(:,:) 
	real(kind=8), allocatable :: fluxResiduals(:,:), dt_val(:)
	character(len=256), allocatable :: bLines(:)
	character(len=256), allocatable :: lines_neighbour(:)
	character(len=256), allocatable :: lines_owner(:)
	integer(kind=8), allocatable :: pts_Faces(:)
	character(len=1000), allocatable :: lines_Faces(:)
!	real(kind=8), allocatable :: first5_refs(:,:)
	type FaceType
		integer(kind=8), allocatable :: points(:)  ! 面的节点索引
	end type FaceType
	type FaceArray
		type(FaceType), allocatable :: faces(:)
	end type FaceArray
	type(FaceType), allocatable :: tmp_faces(:)
	real(kind=8), allocatable :: encodePrimitives3DD(:,:)
	integer(kind=8), allocatable :: boundaryNumFacess(:), boundaryStartFacess(:)
	real(kind=8), allocatable :: initialValues(:)
	real(kind=8), allocatable :: cellState_sln(:,:)
	real(kind=8), allocatable :: cellFluxes_sln(:,:)
	real(kind=8), allocatable :: cellPrimitives_sln(:,:)
	real(kind=8), allocatable :: fluxResiduals_sln(:,:)
	real(kind=8), allocatable :: faceFluxes_sln(:,:)		
	integer(kind=8), allocatable :: cells_mesh(:,:)
	real(kind=8), allocatable :: cVols_mesh(:)
	real(kind=8), allocatable :: cCenters_mesh(:,:)
	real(kind=8), allocatable :: cellSizes_mesh(:,:)		
	integer(kind=8), allocatable :: faces_mesh(:,:)
	real(kind=8), allocatable :: fAVecs_mesh(:,:)
	real(kind=8), allocatable :: fCenters_mesh(:,:)
	integer(kind=8), allocatable :: boundaryFaces_mesh(:,:)
	real(kind=8), allocatable :: dt_solve(:)
	real(kind=8), allocatable :: sol(:,:)
	real(kind=8), allocatable :: facePts(:,:)
   	real(kind=8), allocatable :: fCs(:,:), cell_fAVecs(:,:)
    real(kind=8), allocatable :: cellPts(:,:)
    integer(kind=8), allocatable :: cellFaceCount(:)
    real(kind=8), allocatable :: UunitVec(:)     
    real(kind=8), allocatable :: points_tempMesh(:,:)
!!!!!!!!!meshdeformation
    real(kind=8), allocatable :: points_meshdefor(:,:), mat_inv(:,:)
    real(kind=8), allocatable :: wing_coords(:,:), wing_update_coords(:,:), inoutput_coords(:,:), r_row(:), phi_row(:)
    real(kind=8), allocatable :: control_points(:,:)
    real(kind=8), allocatable :: phi(:,:), distances(:)
    integer(kind=8), allocatable :: owner_tempMesh(:), neighbour_tempMesh(:)
    integer(kind=8), allocatable :: facePoints_tempMesh(:,:)
    real(kind=8), allocatable :: point_update(:,:)
    real(kind=8), allocatable :: wing(:,:), inoutput(:,:)

end module TypesModule
