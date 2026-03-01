!1-10全来自src文件
!1.dataStructures.jl:数据结构   2.FvCFD.jl:主求解函数(subroutine solve)writeoutput_d
!3.boundaryConditions.jl:边界条件
!4.constitutiveRelations.jl:状态变量与通量推导公式   5.JST.jl:通量残差   
!6.mesh.jl:网格    7.numerics.jl:梯度,面值插值,面通量到单元格通量
!8.output.jl:输出文件    9.timeDiscretizations.jl:更新解    
!10.vectorFunctions.jl:向量计算函数
!11.NACA0012.jl:初始化 (来自Examples的文件)
!1.dataStructures.jl:数据结构 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module BuFlowModule	
!    use iso_fortran_env
    ! 2. 再禁用默认类型
    use TypesModule
    use meshdeformationn
    implicit none
!    include 'tapenade_includes.f90'  ! <<< 仅 Tapenade 用，原编译时注释   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!openfoam
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!		
    integer(kind=8), parameter :: N_PRIM_VARS = 5  
! 全局常量：cellFluxes第二维固定为15（3维*5守恒变量）
	integer(kind=8), parameter :: N_FLUX_VARS = 15  
	integer(kind=8), parameter :: max_face_num = 341510
	integer(kind=8), parameter :: DIM_1    = 1_8
	integer(kind=8), parameter :: DIM_2    = 2_8    ! 维度大小2（INTEGER(8)类型）
    integer(kind=8), parameter :: DIM_3    = 3_8    ! 维度大小3
    integer(kind=8), parameter :: DIM_4    = 4_8    ! 维度大小4
    integer(kind=8), parameter :: DIM_5    = 5_8    ! 维度大小5
    integer(kind=8), parameter :: DIM_8    = 8_8    ! 维度大小5
    integer(kind=8), parameter :: DIM_15   = 15_8   ! 维度大小15
    integer(kind=8), parameter :: one_1   = 1_8
	type Celll
		integer(kind=8), allocatable :: faceIndices(:)  ! 面索引列表（动态数组）
		!$AD NONDIFF
		integer(kind=8), allocatable :: pointIndices(:) ! 点索引列表（动态数组）
		!$AD NONDIFF
	end type Celll

	type Fluidd
		real(kind=8) :: Cp
		real(kind=8) :: R
		real(kind=8) :: gammaa
	end type Fluidd
	type SolverStatus
		real(kind=8) :: currentTime         ! 当前仿真时间
		integer(kind=8) :: nTimeSteps       ! 当前已执行的时间步数
		real(kind=8) :: nextOutputTime      ! 下一次写输出的时间点
		real(kind=8) :: endTime             ! 总仿真时间
	end type SolverStatus
	! 或者如果 FaceType 未定义，同时添加 FaceType 定义!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	type FaceType
!		integer(kind=8), allocatable :: points(:)  ! 面的节点索引
!	end type FaceType
!	type FaceArray
!		type(FaceType), allocatable :: faces(:)
!	end type FaceArray
	type BoundaryCondition
		integer(kind=8) :: type  ! 边界类型（用整数标识）
		real(kind=8) :: params(5)  ! 边界参数（如压力、速度等）
	end type BoundaryCondition
	! 定义边界类型常量
    integer(kind=8), parameter :: wallBoundary = 1
    integer(kind=8), parameter :: emptyBoundary = 2
    integer(kind=8), parameter :: InletBoundary = 3
    integer(kind=8), parameter :: OutletBoundary = 4
    type CelllArray
        type(Celll), allocatable :: cells(:)
        !$AD NONDIFF
    end type CelllArray
    ! 在模块开头的类型定义部分添加
	type RestrictResult
		logical :: needTruncate       ! 第1列：是否截断（逻辑值）
		real(kind=8) :: actualDt      ! 第2列：实际时间步长（实数）
	end type RestrictResult
    type MeshData
!		real(kind=8), allocatable :: points(:,:)
!		integer(kind=8), allocatable :: facePoints(:,:)! 面-点索引：(nFaces,4)，四边形面固定4个点
!		integer(kind=8), allocatable :: owner(:), neighbour(:)
		character(len=100), allocatable :: boundaryNames(:)
		integer(kind=8), allocatable :: boundaryNumFaces(:), boundaryStartFaces(:)
	end type MeshData
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
contains  
!2.FvCFD.jl:主求解函数(subroutine solve)!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
    !######################### 初始化 ###########################
    !# 创建一个均匀初始解（所有单元的原始变量均相同）
	subroutine initializeUniformSolution3D(P, T, Ux, Uy, Uz) 
		implicit none
!		type(Meshh), intent(in) :: mesh
		real(kind=8), intent(in) :: P, T, Ux, Uy, Uz
!		real(kind=8), allocatable, intent(inout) :: sol(:,:)!本质是cellPrimitives
		integer(kind=8) :: nCells, meshInfo(4) 
		integer(kind=8) :: nVars, c, nDims
!		real(kind=8), allocatable :: initialValues(:)
		nDims = 3
		call unstructuredMeshInfo(meshInfo)
		nCells = meshInfo(1)
		nVars = 2 + nDims  ! P, T, 和每个方向的速度
!		allocate(sol(nCells, nVars))
!		allocate(initialValues(nVars))

		! 修改1: 修正数组构造器语法
		initialValues(1) = P       ! 第1维：P
		initialValues(2) = T       ! 第2维：T
		initialValues(3) = Ux      ! 第3维：Ux
		initialValues(4) = Uy      ! 第4维：Uy
		initialValues(5) = Uz      ! 第5维：Uz
		
		do c = 1, nCells
			sol(c, :) = initialValues
		end do
		
	end subroutine initializeUniformSolution3D

	! ######################### CFL计算 ###########################
	! 根据当前速度、温度估算各单元的CFL数
	!julia程序里CFL!(CFL, mesh::Mesh, sln::SolutionState, fluid::Fluid, dt=1)的第一项fill!(CFL, 0.0)被强制为0
	subroutine CFL(nCells, fluid, dt, CFLL)
		implicit none
		integer(kind=8), intent(in) :: nCells
		real(kind=8), intent(out) :: CFLL(nCells) ! 函数返回值（原输出参数）
!		type(Meshh), intent(in) :: mesh
!		type(SolutionState), intent(in) :: sln
		type(Fluidd), intent(in) :: fluid
		real(kind=8), intent(in) :: dt
		integer(kind=8) :: nFaces, nBoundaries, nBdryFaces, meshInfo(4)
		integer(kind=8) :: f, ownerCell, neighbourCell
!		real(kind=8), allocatable :: faceRhoT(:,:), cellRhoT(:,:)
		real(kind=8) :: faceRho, faceT, flux, a, dt_val, mag_val, dott
!		real(kind=8), allocatable :: faceVel(:), positionn(:)
		call unstructuredMeshInfo(meshInfo)
		nFaces = meshInfo(2)
		nBoundaries = meshInfo(3)
		nBdryFaces = meshInfo(4)
		CFLL = 0.0_8
		
		! 插值得到面上的密度和温度
!		allocate(faceRhoT(nFaces, 2), source=0.0_8)		
!		allocate(cellRhoT(nCells, 2), source=0.0_8)
		! ② 手动拼接两列数据（第1列：密度；第2列：温度）
		cellRhoT(:, 1) = cellState_sln(:, 1)    ! 第1列：单元格密度（来自 cellState 第1列）
		cellRhoT(:, 2) = cellPrimitives_sln(:, 2)! 第2列：单元格温度（来自 cellPrimitives 第2
		! ③ 调用 linInterp_3D 插值到面
		call linInterp_3D(cellRhoT, faceRhoT)
!		allocate(faceVel(DIM_3), positionn(DIM_3))

		do f = 1, nFaces
			ownerCell = faces_mesh(f,1)
			neighbourCell = faces_mesh(f,2)
			
			faceRho = faceRhoT(f, 1)
			faceT = faceRhoT(f, 2)
			
			if(neighbourCell == -1) then  ! 边界面，用cell值替代
			    faceRho = cellState_sln(ownerCell, 1)
			    faceT = cellPrimitives_sln(ownerCell, 2)
			end if
			
			faceVel(:) = faceFluxes_sln(f, 1:3) / faceRho
			call dot_product_real(faceVel, fAVecs_mesh(f,:), dott, DIM_3)
			flux = abs(dott) * dt
			
			if(faceT <= 0.0_8) then
			    positionn(:) = fCenters_mesh(f,:)
!			    print*, "Warning: Negative temperature at face ", f, ": ", positionn
			end if
			
			! 计算声速
			a = sqrt(fluid%gammaa * fluid%R * faceT)
			call mag(fAVecs_mesh(f,:), mag_val)
			flux = flux + mag_val * a * dt
			
			! 将该face的流量贡献给相邻两个cell
			CFLL(ownerCell) = CFLL(ownerCell) + flux
			if(neighbourCell > -1) then
			    CFLL(neighbourCell) = CFLL(neighbourCell) + flux
			end if
		end do
!		deallocate(faceVel, positionn)
		! 虽然没有return,但在 Julia 中，函数的返回值默认是最后一个表达式的值
		CFLL = CFLL / (2.0_8 * cVols_mesh)  ! 除以两倍体积，符合CFL定义
!		deallocate(faceRhoT)
!		deallocate(cellRhoT)
	end subroutine CFL

	! ######################### 解结构初始化 ###########################
	! 构建SolutionState结构体，并初始化各变量矩阵
	subroutine populateSolution(cellPrimitives, nCells, nFaces, fluid, nDims)
		implicit none
		real(kind=8), intent(in) :: cellPrimitives(:,:)
		integer(kind=8), intent(in) :: nCells, nFaces
		type(Fluidd), intent(in) :: fluid
		integer(kind=8), intent(in), optional :: nDims
!		type(SolutionState), intent(inout) :: sln

		integer(kind=8) :: nConservedVars, nFluxes, nDims_val

		! 默认维数
		nDims_val = 3
		if (present(nDims)) nDims_val = nDims

		nConservedVars = 2 + nDims_val
		nFluxes = nConservedVars * nDims_val

		! 分配结构体数组
!		allocate(sln%cellPrimitives(nCells, 5))
!		allocate(sln%cellState(nCells, nConservedVars))
!		allocate(sln%cellFluxes(nCells, nFluxes))
!		allocate(sln%fluxResiduals(nCells, nConservedVars))
!		allocate(sln%faceFluxes(nFaces, nFluxes))

		! 赋初值
		cellPrimitives_sln = cellPrimitives
		call encodePrimitives3D(cellPrimitives, fluid, cellState_sln)
		cellFluxes_sln     = 0.0_8
		fluxResiduals_sln  = 0.0_8
		faceFluxes_sln     = 0.0_8

		! 用解码函数填充原始变量和通量
		call decodeSolution_3D(fluid)
	end subroutine populateSolution



	! ######################### 时间步长调整 ###########################
	! 判断是否需要限制当前时间步大小$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	subroutine restrictTimeStep(status_val, desiredDt, res)
		implicit none
		type(SolverStatus), intent(in) :: status_val
		real(kind=8), intent(in) :: desiredDt
		type(RestrictResult), intent(inout) :: res  ! 使用派生类型作为返回值
		real(kind=8) :: maxStep
		
		maxStep = min(status_val%endTime - status_val%currentTime, status_val%nextOutputTime - status_val%currentTime)
		
		if (desiredDt > maxStep) then
		    res%needTruncate = .true.    ! 第1列：需要截断
		    res%actualDt = maxStep       ! 第2列：实际时间步长
		else
		    res%needTruncate = .false.   ! 第1列：不需要截断
		    res%actualDt = desiredDt     ! 第2列：实际时间步长
		end if
	end subroutine restrictTimeStep
	! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	! targetCFL初始值为0.2
	! 用于局部时间步长方法 (Local Time Stepping)
	subroutine adjustTimeStep_LTS(targetCFL, dt, status_val, output) 
		implicit none
		real(kind=8), intent(in) :: targetCFL       ! 目标CFL数
		real(kind=8), intent(inout) :: dt(:)        ! 时间步长数组（修改第一个元素）
		type(SolverStatus), intent(in) :: status_val    ! 求解器状态
		real(kind=8), intent(inout) :: output(3)                   ! 返回数组：[writeOutputThisIteration(逻辑值转为实数), dt(1), CFL]		
		type(RestrictResult) :: res                            
		real(kind=8) :: CFLL                         ! 实际CFL数

		! 根据当前时间步数调整CFL（前10步逐渐增加到目标值）
		if (status_val%nTimeSteps < 10) then
		    CFLL = (status_val%nTimeSteps + 1) * targetCFL / 10
		else
		    CFLL = targetCFL
		end if

		! 调用restrictTimeStep限制时间步长
		call restrictTimeStep(status_val, CFLL, res)
		CFLL = res%actualDt  ! 直接将限制后的时间步长赋值给dt(1)
		dt(1) = CFLL
		print *, 'bdt', dt(1)
		! 修正：数组构造器使用逗号分隔，续行符正确
		output = [merge(1.0_real64, 0.0_real64, res%needTruncate), &  ! 第1元素：是否输出
		          dt(1), &        ! 这里事实上dt，但只能输出标量，那就只能输出dt(1),也即CFLL了
		          CFLL]
		! 新增：打印函数内部的output(3)，验证是否为0.5
	end subroutine adjustTimeStep_LTS

	! 更新时间步状态
	subroutine advanceStatus(status_val, dt, CFLL, timeIntegrationFn, silent)
		implicit none
		type(SolverStatus), intent(inout) :: status_val
		real(kind=8), intent(in) :: dt(:), CFLL
		character(len=*), intent(in) :: timeIntegrationFn
		logical, intent(in) :: silent

		! 更新时间步数量
		status_val%nTimeSteps = status_val%nTimeSteps + 1
		! 如果是LTSEuler时间推进方式，则累加模拟时间
		!if (timeIntegrationFn == 'LTSEuler')
		status_val%currentTime = status_val%currentTime + CFLL
		! 输出信息（可选）
		if (.not. silent) then
		    write(*,'(A,I5,A,F9.4,A,F9.4)') 'Timestep: ', status_val%nTimeSteps, &
		        ', simTime: ', status_val%currentTime, ', Max CFLL: ', CFLL
		end if
	end subroutine advanceStatus

    
    ! 默认空气属性初始化函数
    function fluiddd() result(fluid)
	    type(Fluidd) :: fluid
	    fluid%Cp = 1005.0_8
	    fluid%R = 287.05_8
	    fluid%gammaa = 1.4_8
	end function fluiddd
	! ######################### 主求解函数 ###########################
	! CFD主循环
	! 原行过长，拆分为多行
	subroutine solve(meshPath, cellPrimitives, boundaryConditions, point_update, cellPrimitivess)
		implicit none
		! 仅保留必要输入参数（移除已知参数的输入）
!		type(Meshh), intent(in) :: mesh
		character(len=*), intent(in) :: meshPath
		real*8, allocatable, intent(in) :: point_update(:,:)  ! 输入：变形网格点
		real(kind=8), intent(in) :: cellPrimitives(:,:)
		real(kind=8), intent(inout) :: cellPrimitivess(:,:)
		type(BoundaryCondition), intent(in) :: boundaryConditions(:) 
		! 函数返回值（重命名避免与函数名冲突，不改变功能）
		!real(kind=8), allocatable :: solve_result(:,:)		
		! 内部变量声明（完全保留原逻辑）
		integer(kind=8) :: nCells, nFaces, nBoundaries, nBdryFaces, meshInfo(4), tsResult(3), idx_row, idx_col, nFacePerBdry
!		type(SolutionState), allocatable, intent(inout) :: sln
!		real(kind=8), allocatable :: dt_solve(:)
		type(SolverStatus) :: status_val
		real(kind=8) :: CFLL, LTSResult(3)
		logical :: writeOutputThisIteration, silent, restart, createRestartFile, createVTKOutput
		character(len=32) :: timeIntegrationFn, fluxFunction
		real(kind=8) :: initDt, endTime, outputInterval, targetCFL
		type(Fluidd) :: fluid
		character(len=256) :: restartFile
		! 新增：用于输出的变量!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
		integer(kind=8) :: i, j, nPts, nPointsPerCell, l  ! 循环变量、总点数、单元顶点数
		integer(kind=8) :: pointIndices(8)  ! 存储单元的顶点索引（假设最多8个顶点）
		integer(kind=8) :: nConservedVars, nFluxes, nDims_val, print_num 
    	
		! 设置已知参数的默认值（直接固定，不再通过输入修改）
		timeIntegrationFn = 'LTSEuler'  ! 已知时间积分方法
		fluxFunction = 'unstructured_JSTFlux'  ! 已知通量函数
		initDt = 0.0000001                ! 已知初始时间步
		endTime = 5!36.5       ! 已知总仿真时间!!!!!
		outputInterval = 25         ! 已知输出间隔
		targetCFL = 0.5               ! 已知目标CFL数
		fluid = fluiddd()       ! 已知流体属性
		silent = .false.                 ! 已知静默模式
		restart = .false.               ! 已知不重启
		createRestartFile = .true.      ! 已知创建重启文件
		createVTKOutput = .true.        ! 已知创建VTK输出
		restartFile = 'FvCFDRestart.txt'! 已知重启文件路径
		
		! 以下代码完全保留原逻辑，不做任何修改
		if(.not. silent) then
		    print*, 'Initializing Simulation'
		end if

		call unstructuredMeshInfo(meshInfo)
		nCells = meshInfo(1)
		nFaces = meshInfo(2)
		nBoundaries = meshInfo(3)
		nBdryFaces = meshInfo(4)
!!!!!!!!!!!!!!!		

!!!!!!!!!!!!!!!
			
!!!!!!!!!!!!!!!!!!!!!!!!!!!!这些代码原来是populateSolution里面的!!!!!!!!!!!!!!!!!!!!!!!
		! 默认维数
		nDims_val = 3
!		if (present(nDims)) nDims_val = nDims

		nConservedVars = 2 + nDims_val
		nFluxes = nConservedVars * nDims_val
 !       allocate(sln)
		! 分配结构体数组
!		allocate(sln%cellPrimitives(nCells, nConservedVars), source=0.0d0)
!		allocate(sln%cellState(nCells, nConservedVars), source=0.0d0)
!		allocate(sln%cellFluxes(nCells, nFluxes), source=0.0d0)
!		allocate(sln%fluxResiduals(nCells, nConservedVars), source=0.0d0)
!		allocate(sln%faceFluxes(nFaces, nFluxes), source=0.0d0)
!!!!!!!!!!		
!		sln%cellPrimitives(1:nCells, 1:nConservedVars) = cellPrimitives(1:nCells, 1:nConservedVars)
!		!$AD OUTPUT sln%cellPrimitives(1:nCells, 1:5)  ! 输出可微（添加切片）
!		sln%cellState(1:nCells, 1:nConservedVars) = 0.0d0  ! 显式初始化+切片
!		sln%cellFluxes(1:nCells, 1:nFluxes) = 0.0d0        ! 显式初始化+切片
!		sln%fluxResiduals(1:nCells, 1:nConservedVars) = 0.0d0  ! 显式初始化+切片
!		sln%faceFluxes(1:nFaces, 1:nFluxes) = 0.0d0        ! 显式初始化+切片
!!!!!!!!!!!!!!
		! 赋初值
		cellPrimitives_sln = cellPrimitives
!		!$AD OUTPUT sln%cellPrimitives  ! 输出可微
		call encodePrimitives3D(cellPrimitives, fluid, cellState_sln)
		cellFluxes_sln    = 0.0d0
		fluxResiduals_sln  = 0.0d0
		faceFluxes_sln     = 0.0d0
		call decodeSolution_3D(fluid)
!		call populateSolution(cellPrimitives, nCells, nFaces, fluid, 3_int64, sln)
	
        !这里的cellstate没问题,cellprimite和cellfulex也基本没问题
!        allocate(dt_solve(nCells))
		if(timeIntegrationFn == 'LTSEuler') then		    
		    dt_solve = 0.0_8
		else
		    dt_solve = initDt
		end if
		
		status_val%currentTime = 0.0_8
		status_val%nTimeSteps = 0_8
		status_val%nextOutputTime = outputInterval
		status_val%endTime = endtime
		
		if(.not. silent) then
		    print*, 'Starting iterations'
		end if
		print *, "a2"
        !这里primitives没问题
		do while(status_val%currentTime < status_val%endTime)		    		     
	        call adjustTimeStep_LTS(targetCFL, dt_solve, status_val, LTSResult)
			writeOutputThisIteration = (abs(ltsResult(1) - 1.0_8) < 1d-10)  ! 转换逻辑值
			dt_solve(1) = ltsResult(2)
			CFLL = ltsResult(3)
			print *, 'adt', dt_solve(1)
			call LTSEuler(boundaryConditions, fluid, dt_solve)
	        !#dt和CFL都是向量，列数为1，行数为cells	    	    		    
            ! CFLL*ones(size(dt)将标量CFLL扩展为与dt同维度的数组
			! 原错误行：merge(dt, CFLL*ones(size(dt)), ...)
			! 修正为：使用repeat函数创建标量数组或直接初始化数组
			!原julia和现在的Fortran：dt都是向量，CFL是标量
			call advanceStatus(status_val, dt_solve, CFLL, timeIntegrationFn, silent)
			print *, "wh"
	    	if(writeOutputThisIteration) then
		        status_val%nextOutputTime = status_val%nextOutputTime + outputInterval
		    end if
		    print *, "whh"
		end do
		! 简化后的输出：仅保留单元索引和流场数据!!!!!!!!!!!!!!!!!!!!
		print *, "=== 前100个单元的流场数据 ==="
		print *, "总单元数: ", nCells
		print *, "流场数据维度: ", size(cellPrimitives_sln,1), "×", size(cellPrimitives_sln,2), "（行×列）"
		print *, "格式：单元索引 | P | T | Ux | Uy | Uz"

		! 输出前100个单元的流场数据（仅整数和实数，避免格式冲突）
		do i = 1, 15
			WRITE(*, '(I5, 5ES14.6)') i, &
			cellPrimitives_sln(i, 1),  &  ! P（科学计数法，14位宽，6位小数）
			cellPrimitives_sln(i, 2),  &  ! T
			cellPrimitives_sln(i, 3),  &  ! Ux
			cellPrimitives_sln(i, 4),  &  ! Uy
			cellPrimitives_sln(i, 5)     ! Uz  ! U分量（实数）
		end do!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!		cellPrimitives(1:nCells, 1:5) = sln%cellPrimitives(1:nCells, 1:5)
        do idx_row = 1, nCells
			do idx_col = 1, 5
				cellPrimitivess(idx_row, idx_col) = cellPrimitives_sln(idx_row, idx_col)
			end do
		end do
		print_num = 5  ! 最多输出5个单元，避免越界
		PRINT *, '=== 前5个单元的 cellprimitivesd 数据 ==='
		PRINT *, '格式：单元索引 | dP | dT | dUx | dUy | dUz'
!		DO i = 1, print_num
		    ! 输出每个单元的 cellprimitivesd 分量（假设为5列，与 cellprimitives 一致）
!		    WRITE(*, '(I5, 5F12.6)') i, &
!		        cellprimitivesd(i, 1),  &  ! dP
!		        cellprimitivesd(i, 2),  &  ! dT
!		        cellprimitivesd(i, 3),  &  ! dUx
!		        cellprimitivesd(i, 4),  &  ! dUy
!		        cellprimitivesd(i, 5)     ! dUz
!		END DO

!		call writeoutput_d(cellPrimitives_sln, cellprimitivesd, restartfile, meshpath, &  ! 新增cellprimitivesd
!&              createrestartfile, createvtkoutput, point_update)
        call writeOutput(cellPrimitives_sln, restartFile, meshPath, createRestartFile, createVTKOutput, point_update)	        
!		deallocate(dt_solve)
!!!!!!!!!!!		

!!!!!!!!!!!		
	end subroutine solve

	! end

	! 3.boundaryConditions.jl!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! 为一些特殊的边界（如入口，出口，壁面，0梯度边界等）设置特殊的条件即边界条件（例如壁面质量通量为0）

	! ######################### 边界条件实现 #########################
	! 这些边界条件函数的作用是：
	!  1. 在求解器中补充内部面通量计算之外的“边界面通量”；
	!  2. 在 JST 等数值格式中，内部面通量由通量函数计算，但边界面需要专门强制设置；
	!  3. 每个函数实现一种常见边界条件，如超音速入口、出口、壁面、对称等。

	! ------------ 亚音速入口边界条件（使用总压/总温）------------
	! 输入变量为总压、总温、入口速度方向单位向量（nx, ny, nz）
	! 算法来自 FUN3D 求解器：在边界计算静态压强和温度，再重建变量
	subroutine updateInletBoundary(boundaryNumber, inletConditions, fluid, currentBoundary)
		implicit none
!		type(Meshh), intent(in) :: mesh
!		type(SolutionState), intent(inout) :: sln
		integer(kind=8), intent(in) :: boundaryNumber
		real(kind=8), intent(in) :: inletConditions(DIM_5)  ! [总压Pt, 总温Tt, 速度方向nx, ny]
		type(Fluidd), intent(in) :: fluid

		integer(kind=8) :: face, ownerCell
		integer(kind=8), allocatable, intent(inout) :: currentBoundary(:)
		real(kind=8) :: Pt, Tt, nx, ny, nz, gammaa, R
		real(kind=8) :: adjustedVelocity, machNum, e
		real(kind=8) :: primitives(5), state(5), velUnitVector(3), boundaryFluxes(15)
		integer(kind=8) :: nCells, nFaces, nFacePerBdry, face_idx

		nFaces = size(faces_mesh, 1)
		nCells = size(cellState_sln, 1)
		gammaa = fluid%gammaa
		R = fluid%R

		Pt = inletConditions(1)
		Tt = inletConditions(2)
		nx = inletConditions(3)
		ny = inletConditions(4)
		nz = inletConditions(5)

		velUnitVector = [nx, ny, nz]

!		nFacePerBdry = size(mesh%boundaryFaces, 2)	
!		allocate(currentBoundary(nFacePerBdry))
		currentBoundary = boundaryFaces_mesh(boundaryNumber, :)
		nFacePerBdry = size(currentBoundary)
		do face = one_1, nFacePerBdry
		    face_idx = currentBoundary(face)
		    !exit 语句的作用是立即终止当前所在的循环
		    !nFacePerBdry是所有的边界面,currentBoundary
		    if (face_idx == 0_8) exit
		    ownerCell = faces_mesh(face_idx, 1)


		    call dot_product_real(velUnitVector, cellPrimitives_sln(ownerCell, 3:5), adjustedVelocity, 3_8)
		    primitives(2) = Tt - (gammaa - 1.0_8) / (2.0_8 * gammaa * R) * adjustedVelocity**2
		    machNum = abs(adjustedVelocity) / sqrt(gammaa * R * primitives(2))
		    primitives(1) = Pt * (1.0_8 + (gammaa - 1.0_8) / 2.0_8 * machNum**2)**(-gammaa / (gammaa - 1.0_8))
		    primitives(3:5) = adjustedVelocity * velUnitVector

		    call idealGasRho(primitives(2), primitives(1), R, state(1))
		    state(2:4) = primitives(3:5) * state(1)
		    call calPerfectEnergy(primitives(2), fluid, e)
		    state(5) = state(1) * (e + 0.5_8 * sum(primitives(3:5)**2))

		    boundaryFluxes = 0.0_8
		    call calculateFluxes3D(boundaryFluxes, primitives, state)

		    faceFluxes_sln(face_idx, :) = boundaryFluxes
		end do

!		deallocate(currentBoundary)
	end subroutine updateInletBoundary


	! ------------ 压力出口边界条件 ------------
	! 设定出口压力（压强），修正通量中包含压强项的部分
	subroutine updateOutletBoundary(boundaryNumber, outletPressure, fluid, currentBoundary)
		implicit none
!		type(Meshh), intent(in) :: mesh
!		type(SolutionState), intent(inout) :: sln
		integer(kind=8), intent(in) :: boundaryNumber
		integer(kind=8), allocatable, intent(inout) :: currentBoundary(:)  ! 新增
		real(kind=8), intent(in) :: outletPressure(:)
		type(Fluidd), intent(in) :: fluid

		integer(kind=8) :: nFluxes, face, ownerCell, flux, d
!		integer(kind=8), allocatable :: currentBoundary(:)
		real(kind=8) :: outletPressuree, origP
		integer(kind=8) :: nCells, nFaces, nFacePerBdry, face_idx

		nFaces = size(faces_mesh, 1)
		nCells = size(cellState_sln, 1)
		nFluxes = size(cellFluxes_sln, 2)
		nFacePerBdry = size(currentBoundary)
!		nFacePerBdry = size(mesh%boundaryFaces, 2)
!		allocate(currentBoundary(nFacePerBdry))
!		if (stat /= 0) stop "Error: 分配 currentBoundary 失败"
		currentBoundary = boundaryFaces_mesh(boundaryNumber, :)
		outletPressuree = outletPressure(1)
		do face = 1, nFacePerBdry
		    face_idx = currentBoundary(face)
		    if (face_idx == 0) exit

		    ownerCell = faces_mesh(face_idx, 1)
		    ! 复制通量（零梯度假设）
		    do flux = 1, nFluxes
		        faceFluxes_sln(face_idx, flux) = cellFluxes_sln(ownerCell, flux)
		    end do

		    origP = cellPrimitives_sln(ownerCell, 1)
		    do flux = 1, 3
		        faceFluxes_sln(face_idx, 4 * flux) = faceFluxes_sln(face_idx, 4 * flux) + (outletPressuree - origP)
		    end do

		    do d = 1, 3
		        faceFluxes_sln(face_idx, 12 + d) = faceFluxes_sln(face_idx, 12 + d) + &
		                                          cellPrimitives_sln(ownerCell, 2 + d) * (outletPressuree - origP)
		    end do
		end do

!		deallocate(currentBoundary)
	end subroutine updateOutletBoundary


	! ------------ 壁面边界条件（不可穿透）------------
	! 设置动量通量 = 压力，质量和能量通量 = 0
	subroutine updateWallBoundary(boundaryNumber, dummy1, dummy2, currentBoundary)
		implicit none
!		type(Meshh), intent(in) :: mesh
!		type(SolutionState), intent(inout) :: sln
		integer(kind=8), intent(in) :: boundaryNumber
		integer(kind=8), allocatable, intent(inout) :: currentBoundary(:)  
		real(kind=8), intent(in), optional :: dummy1(:)
		type(Fluidd), intent(in) :: dummy2

		integer(kind=8) :: face, stat = 0, i
!		integer(kind=8), allocatable :: currentBoundary(:)
		real(kind=8) :: faceP
		integer(kind=8) :: nCells, nFaces, nFacePerBdry, face_idx, ownerCell
		integer(kind=8) :: realFaceCount

		! 注意：假设mesh中存储了每个边界的真实面数（需在OpenFOAMMesh中保存）
		! 若未保存，需通过其他方式获取（如传入参数）
		
		! 变量初始化
		nFaces = size(faces_mesh, 1)
		nCells = size(cellState_sln, 1)
!		nFacePerBdry = size(boundaryFaces_mesh, 2)
		nFacePerBdry = size(currentBoundary) 
		! 分配currentBoundary并复制数据
!		allocate(currentBoundary(nFacePerBdry))
!		if (stat /= 0) stop "Error: 分配 currentBoundary 失败"
		currentBoundary = boundaryFaces_mesh(boundaryNumber, :)

		! 计算当前边界的真实面数
		realFaceCount = 0
		do i = 1, nFacePerBdry
		    if (currentBoundary(i) == 0) exit
		    realFaceCount = realFaceCount + 1
		end do

		! ====================== 遍历有效面并检查索引 ======================
		 ! 修正：使用realFaceCount作为循环上限，而非遍历到0
		do face = 1, realFaceCount
		    face_idx = currentBoundary(face)
		    ! 检查face_idx有效性（防止随机值越界）
		    if (face_idx < 1 .or. face_idx > nFaces) then
		        cycle  ! 跳过无效面，继续处理下一个
		    end if
		    ! 获取 ownerCell 并检查有效性
		    ownerCell = max(faces_mesh(face_idx, 1), faces_mesh(face_idx, 2))

		    ! ====================== 执行边界条件逻辑 ======================
		    faceP = cellPrimitives_sln(ownerCell, 1)  ! 压力

		    ! 检查 faceFluxes 赋值是否越界
		    faceFluxes_sln(face_idx, 4) = faceP    ! x动量
		    faceFluxes_sln(face_idx, 8) = faceP    ! y动量
		    faceFluxes_sln(face_idx, 12) = faceP   ! z动量
		    faceFluxes_sln(face_idx, 1:3) = 0.0_8   ! 质量通量
		    faceFluxes_sln(face_idx, 13:15) = 0.0_8 ! 能量通量

		end do

!		deallocate(currentBoundary)
	end subroutine updateWallBoundary


	! ------------ Empty 边界条件（2D 网格边界等）------------
	! 什么也不做，占位函数（如周期边界或计算域外）
	subroutine updateEmptyBoundary(boundaryNumber, dummy1, dummy2)
		implicit none
!		type(Meshh), intent(in) :: mesh
!		type(SolutionState), intent(inout) :: sln
		integer(kind=8), intent(in) :: boundaryNumber
		real(kind=8), intent(in), optional :: dummy1(:)
!		integer(kind=8), allocatable :: currentBoundary(:)
		type(Fluidd), intent(in) :: dummy2

		! 空边界什么都不改，直接返回
	end subroutine updateEmptyBoundary


	! 原调用处跟下面变量名一样
	! ------------ 应用所有边界条件的统一入口 ------------
	subroutine updateBoundaryConditions(boundaryConditions, nBoundaries, fluid)
		implicit none
		! 输入输出参数
!		type(Meshh), intent(in) :: mesh
!		type(SolutionState), intent(inout) :: sln
		type(BoundaryCondition), intent(in) :: boundaryConditions(:)
		integer(kind=8), intent(in) :: nBoundaries
		type(Fluidd), intent(in) :: fluid

		! 局部变量
		integer(kind=8) :: boundaryNumber, nFacePerBdry
		integer(kind=8) :: nCells, nFaces, nVars, dim1, dim2

		! 边界类型枚举定义
		integer(kind=8), parameter :: wallBoundary = 1
		integer(kind=8), parameter :: emptyBoundary = 2
		integer(kind=8), parameter :: InletBoundary = 3
		integer(kind=8), parameter :: OutletBoundary = 4

		! 获取网格信息
		nCells = size(cellState_sln, 1)
		nFaces = size(faces_mesh, 1)
		nVars = size(cellState_sln, 2)
		dim1 = 0; dim2 = 0
		if (allocated(faceFluxes_sln)) then
		    dim1 = size(faceFluxes_sln, 1)
		    dim2 = size(faceFluxes_sln, 2)
		end if
		
!        nFacePerBdry = size(mesh%boundaryFaces, 2)
!    	allocate(currentBoundary(nFacePerBdry)) 
		! 遍历边界
		do boundaryNumber = 1, nBoundaries	
			
!		    currentBoundary(1:nFacePerBdry) = mesh%boundaryFaces(boundaryNumber, 1:nFacePerBdry)  
		      
		    select case (boundaryConditions(boundaryNumber)%type)
		    case (wallBoundary)
		        call updateWallBoundary(boundaryNumber, boundaryConditions(boundaryNumber)%params, fluid, currentBoundary)
		    case (emptyBoundary)
		        call updateEmptyBoundary(boundaryNumber, boundaryConditions(boundaryNumber)%params, fluid)

		    case (InletBoundary)
		        call updateInletBoundary(boundaryNumber, boundaryConditions(boundaryNumber)%params, fluid, currentBoundary)

		    case (OutletBoundary)
		        call updateOutletBoundary(boundaryNumber, boundaryConditions(boundaryNumber)%params, fluid, currentBoundary)
		    case default	
		        stop
		    end select
		end do
!		deallocate(currentBoundary)
	end subroutine updateBoundaryConditions
	! 4.constitutiveRelations.jl!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! ######################### Internal/Private functions #######################
	! 内部函数：理想气体状态方程和能量计算

	! 根据理想气体状态方程计算密度 rho = P / (R*T)
	subroutine idealGasRho(T, P, R, rho_val)
		implicit none
		real(kind=8), intent(in) :: T, P, R
		real(kind=8), intent(out) :: rho_val		
		rho_val = P/(R*T)
	end subroutine idealGasRho

	! 根据理想气体状态方程计算压强 P = rho * R * T
	subroutine idealGasP(rho, T, R, P)
		implicit none
		real(kind=8), intent(in) :: rho, T
		real(kind=8), intent(in) :: R
		real(kind=8), intent(out) :: P
		! PV = mRT 推导得到 P = rho * R * T
		P = rho*R*T
	end subroutine idealGasP

	! 计算理想气体的比内能 e = T * (Cp - R)
	subroutine calPerfectEnergy(T, fluid, e_val)
		implicit none
		real(kind=8), intent(in) :: T
		type(Fluidd), intent(in) :: fluid
		real(kind=8), intent(out) :: e_val  ! 输出参数：原返回值
		e_val = T * (fluid%Cp - fluid%R)
	end subroutine calPerfectEnergy    

	! 反过来，根据比内能 e 计算温度 T = e / (Cp - R)
	subroutine calPerfectT(e, fluid, T)
		implicit none
		real(kind=8), intent(in) :: e
		type(Fluidd), intent(in) :: fluid
		real(kind=8), intent(out) :: T
		
		T = e / (fluid%Cp - fluid%R)
	end subroutine calPerfectT

	! ######################### Public functions (called by FvCFD.jl) #######################
	! 公共函数：供外部模块（如 FvCFD.jl）调用

	! =
	!     从守恒变量 cellState 解码出原始变量 primitives（压强、温度、速度）
	!     输入：
	!         cellState: [rho, rho*Ux, rho*Uy, rho*Uz, rho*eV2] （守恒变量）
	!         fluid: Fluid 类型，包含 Cp 和 R
	!     输出：
	!         primitives: [P, T, Ux, Uy, Uz]
	! =#
    subroutine decodePrimitives3D(primitives, cellState, fluid)
		implicit none
		real(kind=8), intent(inout) :: primitives(5)  ! 输入输出：原始变量（将被更新）
		real(kind=8), intent(in)    :: cellState(5)   ! 输入：守恒变量
		type(Fluidd), intent(in)    :: fluid

		real(kind=8) :: e, vel_mag, temp_T

		! 速度分量 Ux, Uy, Uz
		primitives(3) = cellState(2) / cellState(1)  ! Ux
		primitives(4) = cellState(3) / cellState(1)  ! Uy
		primitives(5) = cellState(4) / cellState(1)  ! Uz

		! 动能与比内能
		vel_mag = sqrt(primitives(3)**2 + primitives(4)**2 + primitives(5)**2)
		e = (cellState(5) / cellState(1)) - 0.5 * vel_mag**2

		! 温度 T
		call calPerfectT(e, fluid, primitives(2))
        ! 关键修改：用临时变量 temp_T 存储 primitives(2) 的值
    	temp_T = primitives(2)
		! 压强 P = rho * R * T
		call idealGasP(cellState(1), temp_T, fluid%R, primitives(1))
	end subroutine decodePrimitives3D

	! =
	!     从原始变量 primitives 计算守恒变量 cellState
	!     输入：
	!         cellPrimitives: [P, T, Ux, Uy, Uz]，每行代表一个 cell
	!         fluid: 流体属性，包含 Cp 和 R
	!     输出：
	!         cellState: [rho, rho*Ux, rho*Uy, rho*Uz, rho*eV2]
	! =#
	subroutine encodePrimitives3D(cellPrimitives, fluid, encodePrimitives3DD)
		implicit none
		real(kind=8), intent(in) :: cellPrimitives(:,:)
		type(Fluidd), intent(in) :: fluid
		real(kind=8), allocatable, intent(inout) :: encodePrimitives3DD(:,:)
		integer(kind=8) :: nCells, c
		real(kind=8) :: e, mag_vel
		
		nCells = size(cellPrimitives, 1)
		
		! 初始化守恒变量矩阵
!		allocate(encodePrimitives3D(nCells, 5))
		
		do c = 1, nCells
			! 根据 T 和 P 计算密度
			call idealGasRho(cellPrimitives(c,2), cellPrimitives(c,1), fluid%R, encodePrimitives3DD(c,1))
			! 计算动量项 rho * U
			encodePrimitives3DD(c,2:4) = cellPrimitives(c,3:5) * encodePrimitives3DD(c,1)
			! 计算内能密度 e = T*(Cp - R)
			call calPerfectEnergy(cellPrimitives(c,2), fluid, e)
			! 总能量 = rho * (e + 动能)
			call mag(cellPrimitives(c,3:5), mag_vel)  ! 计算速度模长
    		! 用临时变量 mag_vel 替代原 mag(...)
    		encodePrimitives3DD(c,5) = encodePrimitives3DD(c,1)*(e + (mag_vel**2)/2.0_8)
		end do
		
	end subroutine encodePrimitives3D

	! =
	!     根据 cell 的状态和原始变量，计算通量 fluxes（用于求解方程的显式格式）
	!     输入：
	!         fluxes: 输出通量向量（长度为15）
	!         prim: 原始变量 [P, T, Ux, Uy, Uz]
	!         state: 守恒变量 [rho, rho*Ux, rho*Uy, rho*Uz, rho*eV2]
	!     输出：
	!         无（直接修改 fluxes 向量）
	! 
	!     通量结构如下：
	!         1-3: 质量通量 (rho*U)
	!         4-12: 动量通量张量（3x3）
	!         13-15: 能量通量 (rho*eV2*U + P*U)
	! =#
	! 原始变量primitives = Vector{Float64}(undef, 5)    5：[P, T, Ux, Uy, Uz]
	! 守恒变量state = Vector{Float64}(undef, 5)   5:[ρ, ρUx, ρUy, ρUz, ρe]
	! 函数名以 ! 结尾，它直接修改传入的 fluxes 数组，而不是返回一个新的数组。且在调用这个函数时可不初始化fluex的值,且在调用这个函数时，可以直接使用输出的结果(无return)
	subroutine calculateFluxes3D(out_fluxes, prim, state)
		implicit none
		real(kind=8), intent(out) :: out_fluxes(15)  ! 输出：计算后的通量
		real(kind=8), intent(in) :: prim(5)          ! 原始变量
		real(kind=8), intent(in) :: state(5)         ! 守恒变量

		! 质量通量
		out_fluxes(1) = state(2)           ! rho*Ux
		out_fluxes(2) = state(3)           ! rho*Uy
		out_fluxes(3) = state(4)           ! rho*Uz

		! 动量通量
		out_fluxes(4) = state(2)*prim(3) + prim(1)  ! x动量x通量
		out_fluxes(7) = state(2)*prim(4)             ! y动量x通量
		out_fluxes(10) = state(2)*prim(5)            ! z动量x通量
		out_fluxes(5) = out_fluxes(7)                 ! x动量y通量（对称）
		out_fluxes(8) = state(3)*prim(4) + prim(1)   ! y动量y通量
		out_fluxes(11) = state(3)*prim(5)             ! z动量y通量
		out_fluxes(6) = out_fluxes(10)                ! x动量z通量（对称）
		out_fluxes(9) = out_fluxes(11)                 ! y动量z通量（对称）
		out_fluxes(12) = state(4)*prim(5) + prim(1)  ! z动量z通量

		! 能量通量
		out_fluxes(13) = prim(3)*state(5) + prim(1)*prim(3)
		out_fluxes(14) = prim(4)*state(5) + prim(1)*prim(4)
		out_fluxes(15) = prim(5)*state(5) + prim(1)*prim(5)
	end subroutine calculateFluxes3D

! 5.JST.jl!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! 4-55(函数unstructured_JSTEps)计算 JST 方法中的二阶（eps2）和四阶（eps4）人工扩散系数
	! 56-108(函数unstructured_JSTFlux)执行 JST 方法计算通量残差
	! 包装函数：兼容ALLOCATABLE要求+Tapenade维度推导
	
	subroutine unstructured_JSTEps(fluid, eps)
		implicit none
		! 输入参数
!		type(Meshh), intent(in) :: mesh
!		type(SolutionState), intent(in) :: sln
		type(Fluidd), intent(in) :: fluid		
		! 输出结果：二维数组，第1列=eps2，第2列=eps4
		real(kind=8), intent(inout) :: eps(:,:)		
		! 局部变量
		integer(kind=8) :: nCells, nFaces, nBoundaries, nBdryFaces, meshInfo(4)
		
		integer(kind=8) :: f, ownerCell, neighbourCell, c, nVars 
		real(kind=8) :: d(3), oP, nP, farOwnerP, farNeighbourP, epsilonn, dot_ownerCell, dot_neighbourCell
		real(kind=8) :: k2_val, k4_val, c4_val, vec_mag   ! JST系数默认值
		integer(kind=8), parameter :: DIM_SPACE = 3  ! 空间维度（x/y/z）
	    integer(kind=8), parameter :: P_var = 1    ! 压力变量数（仅1个）
	    integer(kind=8), parameter :: N_PRIM_VEL = 3 ! 速度分量（Ux/Uy/Uz）
		integer(kind=8), parameter :: N_EPS_COLS = 2 
		! 新增：静态中间数组（编译期确定维度，Tapenade可推导）
		! 1. 设置JST系数默认值（可选参数未输入时使用）
		k2_val = 0.5_8       ! 默认k2
		k4_val = 1.0_8/32.0_8 ! 默认k4
		c4_val = 1.0_8       ! 默认c4
		
		! 2. 获取网格信息（单元数、面数等）
		call unstructuredMeshInfo(meshInfo)
		nCells = meshInfo(1)
		nFaces = meshInfo(2)
		nBoundaries = meshInfo(3)
		nBdryFaces = meshInfo(4)

		! 3. 提取压力并计算压力梯度（用于激波检测）
!		allocate(P_eps(nCells))
!		allocate(temp_grad(nCells, 1, 3))  ! 临时存储梯度（单元数×3方向×1变量）
		print *, 'hhh'
!		allocate(gradP(nCells, 3))
!		allocate(sj(nCells), sjCount(nCells))
!		allocate(rj(nCells)) 
!		allocate(rjsjF(nFaces, 2), source=0.0_8)  ! 第1列=rj，第2列=sj
		! 7. 计算JST格式的耗散系数eps2和eps4
!		allocate(eps2(nFaces), source=0.0_8)
!		allocate(eps4(nFaces), source=0.0_8)
!!		allocate(eps(nFaces, 2), source=0.0d0)
!		P_eps = cellPrimitives_sln  ! 压力（原始变量第1列）
!		P_matrix(:, 1) = P_eps
		
		! 转换为二维矩阵（适配greenGaussGrad函数输入要求）
!		allocate(P_matrix(nCells, P_var))
		P_matrix(:, 1) = cellPrimitives_sln(:,1)
       
!		call greenGaussGrad(mesh, P_matrix, .false., temp_grad)  ! 计算压力梯度
!!!!!!!!!!!!!        
!        allocate(faceVals(nFaces, 5))
        faceVals = 0.0_8
	    call greenGaussGrad(P_matrix, .false., temp_grad, faceVals) 
!!!!!!!!!!!!!  
		gradP = temp_grad(:, :, 1) 	     
!		gradP(:, 1) = temp_grad(:, 1, 1)  ! 明确提取x方向梯度（第2维第1个元素）
!		gradP(:, 2) = temp_grad(:, 1, 2)  ! 明确提取y方向梯度（第2维第2个元素）
!		gradP(:, 3) = temp_grad(:, 1, 3) 
		! 4. 计算激波传感器（shock sensor）sj
		sj = 0.0_8         ! 激波传感器值
		sjCount = 0.0_8    ! 计数：每个单元被多少面访问
		
		epsilonn = 1.0d-10 ! 防止除零的小量
		print *, 'ha'
		do f = 1, nFaces - nBdryFaces  ! 只处理内部面（排除边界面）
		    ownerCell = faces_mesh(f, 1)
		    neighbourCell = faces_mesh(f, 2)
		    
		    ! 单元质心向量差（neighbour - owner）
		    d = cCenters_mesh(neighbourCell, :) - cCenters_mesh(ownerCell, :)
			    
		    ! 提取压力值
		    oP = P_matrix(ownerCell, 1)
		    nP = P_matrix(neighbourCell, 1)
		    call dot_product_real(d, gradP(ownerCell, :), dot_ownerCell, 3_8)
		    call dot_product_real(d, gradP(neighbourCell, :), dot_neighbourCell, 3_8)
		    ! 构造远端虚拟单元压力（基于梯度外推）
		    farOwnerP = nP - 2.0_8 * dot_ownerCell
		    farNeighbourP = oP + 2.0_8 * dot_neighbourCell
		    
		    ! 计算激波传感器值（累加）
		    sj(ownerCell) = sj(ownerCell) + (abs(nP - 2.0_8*oP + farOwnerP) &
		        / max(abs(nP - oP) + abs(oP - farOwnerP), epsilonn))**2
		    sjCount(ownerCell) = sjCount(ownerCell) + 1
		    
		    sj(neighbourCell) = sj(neighbourCell) + (abs(oP - 2.0_8*nP + farNeighbourP) &
		        / max(abs(farNeighbourP - nP) + abs(nP - oP), epsilonn))**2
		    sjCount(neighbourCell) = sjCount(neighbourCell) + 1
		end do
		print *, 'h'
		! 5. 计算谱半径（rj）和平均激波传感器（sj）
		 ! 谱半径 = 流速 + 声速
		do c = 1, nCells
		    ! 谱半径：|速度| + 声速
		    ! 原错误行拆分为两步：1. 调用mag获取模长；2. 计算rj(c)
			call mag(cellPrimitives_sln(c,3:5), vec_mag)  ! 第一步：调用mag，结果存到vec_mag
			rj(c) = vec_mag + real(sqrt(fluid%gammaa * fluid%R * cellPrimitives_sln(c,2)), kind=8)
		    
		    ! 平均sj（按面数平均）
		    if(sjCount(c) > 0) then
		        sj(c) = sj(c) / sjCount(c)
		    end if
		end do
		
		! 6. 插值到面（取两侧单元最大值）		
		call maxInterp(sj, rj, rjsjF)  ! 调用插值函数				
		
		do f = 1, nFaces - nBdryFaces  ! 内部面计算耗散
		    eps2(f) = k2_val * rjsjF(f,2) * rjsjF(f,1)  ! eps2 = k2 * sj * rj
		    eps4(f) = max(0.0_8, k4_val * rjsjF(f,1) - c4_val * eps2(f))  ! eps4 = max(0, k4*rj - c4*eps2)
		end do
		
		! 8. 输出结果：合并eps2和eps4为二维数组
!		allocate(eps(nFaces, 2), source=0.0d0)
		eps(:,1) = eps2  ! 第1列=eps2
		eps(:,2) = eps4  ! 第2列=eps4
		print *, 'h'
		! 9. 释放临时内存
!		deallocate(P, P_matrix, temp_grad, gradP)
!		deallocate(sj, sjCount, rj, rjsjF, eps2, eps4)
!!!!!!!!!!		
!		deallocate(faceVals)
!!!!!!!!!!			
        ! 释放前检查数组是否已分配，避免释放未分配内存导致double free
		! 先释放临时数组，再释放核心数组，避免依赖
		print *, 'hhh'
	end subroutine unstructured_JSTEps 

	! = 
	!     输入：sln.cellState, sln.cellPrimitives, sln.cellFluxes 应该是最新的
	!     输出：更新 sln.faceFluxes 和 sln.cellResiduals
	!     返回：更新后的 sln.cellResiduals
	! 
	!     实现 JST 方法：中央差分 + 人工粘性（JST扩散），每个面作为一维问题处理
	!     参考：https://aero-comlab.stanford.edu/Papers/jst_2015_updated_07_03_2015.pdf （特别是第5-6页）
	! =#
	! sln的初始值在第206行(主求解函数function solve)定义的
	! JST扩散系数计算的包装函数
	subroutine unstructured_JSTFlux(boundaryConditions, fluid, unstructured_JSTFluxx)
		implicit none
!		type(Meshh), intent(in) :: mesh
!		type(SolutionState), intent(inout) :: sln
		type(BoundaryCondition), intent(in) :: boundaryConditions(:)
		type(Fluidd), intent(in) :: fluid
		real(kind=8) :: d(3), grad_v(3), dot_grad_v
		real(kind=8), intent(inout) :: unstructured_JSTFluxx(:,:)
!		real(kind=8), allocatable :: eps22(:), eps44(:), diffusionFlux(:), unitFA(:), fD(:), farOwnerfD(:), farNeighbourfD(:), epss(:,:)
!		real(kind=8), allocatable :: fDeltas(:,:), fDGrads(:,:,:)
		integer(kind=8) :: f, ownerCell, neighbourCell, v, i1, i2, nFaces, nBoundaries, nBdryFaces, meshInfo(4), i, max_col
		INTEGER(kind=8) :: print_eps, print_face, nCells, nVars! 控制输出面数量	
		integer(kind=8), parameter :: DIM_SPACE = 3  ! 空间维度（x/y/z）
   		integer(kind=8), parameter :: N_EPS_COLS = 2
		call unstructuredMeshInfo(meshInfo)
		nCells = meshInfo(1)
		nFaces = meshInfo(2)
		nBoundaries = meshInfo(3)
		nBdryFaces = meshInfo(4)
		nVars = size(cellState_sln, 2) 
		
   ! 维度2：列数
		! #### 1. 应用边界条件 ####		
		call updateBoundaryConditions(boundaryConditions, nBoundaries, fluid)
		! #### 2. 计算中央差分通量 ####
		call linInterp_3D(cellFluxes_sln, faceFluxes_sln)
		! 输出前10个cellFluxes
!		write(*,*) "First 10 cellFluxes:"
!		do i = 1, 10
!		    write(*,'(I5,*(E12.5))') i, cellFluxes_sln(i,1:size(cellFluxes_sln,2))
!		end do
		
		! 输出前10个faceFluxes - 修正这里使用size获取行数
!		write(*,*) "First 10 faceFluxes:"
!		do f = 1, 10
!		    write(*,'(I5,*(E12.5))') f, faceFluxes_sln(f,:)
!		end do
		! #### 3. 计算 JST 人工扩散通量 ####
!		allocate(fDeltas(nFaces, nVars), source=0.0_8)		
		call faceDeltas(fDeltas)                      ! 单元变量在面的差值
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!		allocate(fDGrads(nCells, nVars, DIM_SPACE))
		print *, 'AAA'
		call greenGaussGradd(fDeltas, .false., fDGrads)    ! 差值的梯度
		print *, 'BBB'
!		allocate(eps22(nFaces), eps44(nFaces))
!		allocate(eps(nFaces, 2), source=0.0_8)
		call unstructured_JSTEps(fluid, epss)  ! 二阶和四阶人工扩散系数
		print *, 'CCC'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
		eps22 = epss(:,1)  ! 二阶耗散系数
		eps44 = epss(:,2)  ! 四阶耗散系数
!		deallocate(eps)  ! 及时释放临时数组eps（关键修正）
!		allocate(diffusionFlux(nVars), unitFA(3), fD(nVars), farOwnerfD(nVars), farNeighbourfD(nVars))
		do f = 1, nFaces-nBdryFaces
			ownerCell = faces_mesh(f, 1)
			neighbourCell = faces_mesh(f, 2)
			! 这个中心点cCenters是根据面中心点和体中心点加权得到
			d = cCenters_mesh(neighbourCell, :) - cCenters_mesh(ownerCell, :)

			! 获取该面的状态差值
			fD = fDeltas(f,:)
			! 虚拟单元的状态差值（逐个变量计算梯度点积）
			do v = 1, nVars  ! 遍历每个变量
				! 提取第v个变量在ownerCell处的梯度向量（1维数组：x,y,z方向梯度）				
				grad_v = real(fDGrads(ownerCell, v, :), kind=8)   ! 维度：(单元格, 变量v, 空间维度) → 1维
	            call dot_product_real(d, grad_v, dot_grad_v, 3_8)
				! 计算点积（d和grad_v都是1维向量，符合dot_product要求）
				farOwnerfD(v) = fD(v) - dot_grad_v
				farNeighbourfD(v) = fD(v) + dot_grad_v
			end do

			! 人工扩散项 = eps2 * 一阶差值 - eps4 * 二阶差值
!			diffusionFlux = eps2(f) * fD - eps44(f) * (farNeighbourfD - 2*fD + farOwnerfD)
            diffusionFlux = eps22(f) * fD 
			! 将扩散通量添加到面通量上，方向用单位面积向量调整
			! mesh.fAVecs是面法向量(通过面的子三角形的叉乘求法向量，然后相加作为面的法向量) ,normalize将法向量转化为单位向量
			call normalize(fAVecs_mesh(f,:), unitFA)
			
			do v = 1, nVars
			    i1 = (v-1)*3 + 1
			    i2 = i1 + 2
			    ! !!!!!!!!!线性加权计算出来的面通量就代表这两包这件通过这个面传递的通量残差!!!!!!!!!
			    faceFluxes_sln(f,i1:i2) = faceFluxes_sln(f,i1:i2) - (diffusionFlux(v) * unitFA)
			end do
		end do
		! 输出前10个faceFluxes - 修正这里使用size获取行数
		write(*,*) "First 10 faceFluxes:"
		do f = 1, 10
		    write(*,'(I5,*(E12.5))') f, faceFluxes_sln(f,:)
		end do
		! #### 4. 将面通量积分到单元内部，得到单元残差 ####
!		allocate(unstructured_JSTFluxx(nCells, nVars), source=0.0d0)
		call integrateFluxes_unstructured3D(unstructured_JSTFluxx)  
!		deallocate(fDeltas, fDGrads, eps22, eps44, diffusionFlux, unitFA, fD, farOwnerfD, farNeighbourfD)
	end subroutine unstructured_JSTFlux

	! 6.mesh.jl!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! 引用自 Moukalled 等人的《FVM》书籍——适用于 OpenFOAM、Matlab 等平台的计算方法
	! 1-98行(有限体积法准备工作): 计算单元/面的几何中心、面积向量、体积等关键几何属性
	! 98-310: 读取和处理OpenFOAM网格数据，包括点 (points)、面 (faces)、所有权 (owner)、邻接关系 (neighbour) 以及边界信息 (boundary)
	! 310-500: 将 OpenFOAM 网格文件夹中的面和点信息转换为单元（cell）结构，并确保其点的顺序满足 .vtk 文件格式的要求
	! 500-612: 主函数!!!!!!!，读取OpenFOAM 网格文件，并输出为FvCFD.jl可后续使用的Mesh 类型

	! ######################### 网格/单元几何相关函数 ###########################

	! 计算三角形（或多边体）的几何中心!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	function triangleCentroid(points)
		implicit none
		real(kind=8), intent(in) :: points(:,:)
		real(kind=8) :: triangleCentroid(3)
		integer(kind=8) :: nPts, pt
		
		triangleCentroid = [0.0_8, 0.0_8, 0.0_8]
		nPts = size(points, 1)
		
		do pt = 1, nPts
			triangleCentroid = triangleCentroid + points(pt,:)   ! 累加每个点坐标
		end do
		
		triangleCentroid = triangleCentroid / nPts   ! 求平均，即几何中心
	end function triangleCentroid

	! 为几何中心函数设置一个别名
	function geometricCenter(points)
		implicit none
		real(kind=8), intent(in) :: points(:,:)
		real(kind=8) :: geometricCenter(3)
		geometricCenter = triangleCentroid(points)
	end function geometricCenter

	! 计算三角形面积向量!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	function triangleArea(points)
		implicit none
		real(kind=8), intent(in) :: points(:,:)
		real(kind=8) :: triangleArea(3)
		real(kind=8) :: side1(3), side2(3)
		
		side1 = points(2,:) - points(1,:)
		side2 = points(3,:) - points(1,:)
		triangleArea = cross(side1, side2) / 2.0_8   ! 使用向量叉乘除以2
	end function triangleArea

	! 计算多边形面元的面积向量与几何中心!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine faceAreaCentroid(points, fAVec, centroid)
	!subroutine faceAreaCentroid(points, fAVec, centroid)
		implicit none
		real(kind=8), intent(in) :: points(:,:)
		real(kind=8), intent(out) :: fAVec(3), centroid(3)
		real(kind=8) :: gC(3), subTriPts(3,3), triCentroid(3), subFAVec(3)
		real(kind=8) :: mag_subFAVec, mag_fAVec  ! 新增：两个临时变量
		integer(kind=8) :: nPts, i
!		logical :: print_output  ! 标记是否需要输出
!		real(kind=8), intent(in) :: first5_refs(5,3)   ! 前5个面的参考顶点坐标!!!!!!!!!!!!!!
		integer(kind=8) :: current_face_id!, f ! 当前面的索引（1-5）!!!!!!!!!!!!!!!!!!!
!		logical :: is_in_first5!!!!!!!!!
		! 判断是否为前5个面（通过第一个顶点坐标匹配）!!!!!!!!
	
!		is_in_first5 = .false.
!		current_face_id = 0
!		do f = 1, 5
		    ! 允许微小数值误差（1e-6）
!		    if (abs(points(1,1)-first5_refs(f,1)) < 1d-6 .and. &
!		        abs(points(1,2)-first5_refs(f,2)) < 1d-6 .and. &
!		        abs(points(1,3)-first5_refs(f,3)) < 1d-6) then
!		        is_in_first5 = .true.
!		        current_face_id = f
!		        exit
!		    end if
!		end do!!!!!!!!!!!!!!!!
		gC = geometricCenter(points)       ! 得到整个面的几何中心
		nPts = size(points, 1) 
	      
		fAVec = [0.0_8, 0.0_8, 0.0_8]      ! 初始化面面积向量
		centroid = [0.0_8, 0.0_8, 0.0_8]   ! 初始化面中心
	      
	 	! 仅当处理第1个面时才输出
!		print_output = (f == 1)!!!!!!!!!!!!!!!
		! 将多边形面分割成多个由几何中心和两个相邻顶点组成的三角形!!!!!!!!!!!!!!!!!!!!
		do i = 1, nPts
			if (i < nPts) then
			    subTriPts(1,:) = gC
			    subTriPts(2,:) = points(i,:)
			    subTriPts(3,:) = points(i+1,:)
			else
			    subTriPts(1,:) = gC
			    subTriPts(2,:) = points(i,:)
			    subTriPts(3,:) = points(1,:)  ! 闭合面
			end if
            !subTriPts没问题
			triCentroid = triangleCentroid(subTriPts)    ! 子三角形中心
			subFAVec = triangleArea(subTriPts)           ! 子三角形面积向量
			! 前5个面输出subFAVec
!		    if (is_in_first5) then
!	        print *, "=== 前5个面的子三角形 ", i, " ==="
!	        print *, "  subFAVec (x,y,z): (", subFAVec(1), ",", subFAVec(2), ",", subFAVec(3), ")"
!	        print *, "  面积大小: ", mag_subFAVec  ! 面积大小（标量）
!		    end if					
			fAVec = fAVec + subFAVec			
			call mag(subFAVec, mag_subFAVec)  ! 计算 subFAVec 模长	
			centroid = centroid + triCentroid * mag_subFAVec
		end do

		call mag(fAVec, mag_fAVec)  ! 计算 fAVec 模长
		centroid = centroid / mag_fAVec  ! 用临时变   ! 总质心为加权平均
!		if (is_in_first5) then
!        print *, "centroid=", centroid 
 !       end if
	end subroutine faceAreaCentroid

	! 原调用处:cVols[c], cCenters[c] = cellVolCentroid(pts, cell_fAVecs, fCs),fCs为一个单元中各个面的几何中心,cell_fAVecs为这个包对应面的法向量
	! 计算单元体积和中心位置!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine cellVolCentroid(points, fAVecs, faceCentroids, vol, centroid)
		implicit none
		! 假设输入数组维度：
		! points: (nPoints, 3)  单元所有顶点坐标
		! fAVecs: (nFaces, 3)   每个面的面积向量（每行一个面）
		! faceCentroids: (nFaces, 3) 每个面的中心（每行一个面）
		real(kind=8), intent(in) :: points(:,:), fAVecs(:,:), faceCentroids(:,:)
		real(kind=8), intent(out) :: vol, centroid(3)
		real(kind=8) :: gC(3), cellCenterVec(3), subPyrVol, subPyrCentroid(3), dott
		integer(kind=8) :: nFaces, f
		
		! 计算单元几何中心（确保与Julia的geometricCenter逻辑一致）
		gC = geometricCenter(points)
		
		! 关键修正：正确获取面数（若fAVecs是(nFaces,3)，则用size(fAVecs,1)）
		! 若fAVecs是(3,nFaces)，则改为 nFaces = size(fAVecs, 2)
		nFaces = size(fAVecs, 1)  

		vol = 0.0_8
		centroid = [0.0_8, 0.0_8, 0.0_8]

		do f = 1, nFaces
		    ! 计算面中心到单元几何中心的向量（修正数组访问）
		    cellCenterVec = faceCentroids(f,:) - gC  ! 若faceCentroids是(3,nFaces)，则用faceCentroids(:,f)
		    call dot_product_real(fAVecs(f,:), cellCenterVec, dott, 3_8)
		    ! 计算子金字塔体积（点积+绝对值）
		    subPyrVol = abs( dott ) / 3.0_8  
		    ! 改用dot_product确保点积计算正确（与Julia的sum(.+)一致）
		    
		    ! 计算子金字塔中心
		    subPyrCentroid = 0.75_8 * faceCentroids(f,:) + 0.25_8 * gC
		    
		    ! 累加体积和中心
		    vol = vol + subPyrVol
		    centroid = centroid + subPyrCentroid * subPyrVol
		end do

		centroid = centroid / vol
	end subroutine cellVolCentroid

	! 计算每个面从其所属单元中心指向面中心的向量!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	function cellCentroidToFaceVec(faceCentroids, cellCentroids)
		implicit none
		real(kind=8), intent(in) :: faceCentroids(:,:), cellCentroids(3)
		real(kind=8) :: cellCentroidToFaceVec(size(faceCentroids,1), 3)
		integer(kind=8) :: nFaces, f
		
		nFaces = size(faceCentroids, 1)

		do f = 1, nFaces
			cellCentroidToFaceVec(f,:) = faceCentroids(f,:) - cellCentroids  ! 注意：cellCentroids 应是单个中心点？
		end do
	end function cellCentroidToFaceVec

	! ######################### 工具函数 ###########################

	! 返回网格的基本信息，包括单元数、面数、边界数量、边界面数!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine unstructuredMeshInfo(info)
		implicit none
!		type(Meshh), intent(in) :: mesh  ! 输入网格对象（与Julia的Mesh对应）
		integer(kind=8), intent(out) :: info(4)       ! 返回数组：[nCells, nFaces, nBoundaries, nBdryFaces]
		integer(kind=8) :: nCells, nFaces, nBoundaries, nBdryFaces
		integer(kind=8) :: bdry, countt, i         ! 循环变量（边界索引）
		
		! 计算单元总数（网格单元的行数）
		nCells = size(cells_mesh, 1)
		! 计算面总数（网格面的行数）
		nFaces = size(faces_mesh, 1)
		! 计算边界总数（边界面的行数）
		nBoundaries = size(boundaryFaces_mesh, 1)		
		! 统计边界面总数（累加每个边界的面数量）
		! 正确计算边界面数（忽略填充的0）
		nBdryFaces = 0
		do bdry = 1, nBoundaries
		    ! 统计当前边界的非零元素个数
		    i = 0
		    do while (i < size(boundaryFaces_mesh, 2) .and. boundaryFaces_mesh(bdry, i+1) /= 0)
		        i = i + 1
		    end do
		    nBdryFaces = nBdryFaces + i
		end do
		! 赋值返回结果（保持与Julia相同的顺序）
		info = [nCells, nFaces, nBoundaries, nBdryFaces]
	end subroutine unstructuredMeshInfo
	
    ! 核心函数：读取并处理OpenFOAM网格
    subroutine OpenFOAMMesh(polyMeshPath, point_update, tempMesh)
        character(len=*), intent(in) :: polyMeshPath
        real*8, allocatable, intent(in) :: point_update(:,:)  ! 输入：变形网格点
!        type(Meshh), intent(inout) :: mesh     	    
        type(MeshData), intent(in) :: tempMesh
        integer(kind=8) :: nCells, nFaces, nBoundaries, nPoints, new_size, count_cellPtss
        integer(kind=8) :: f = 0, c, b, i, j, k, startF, endF, faceCount, max_faces_per_cell, j_max, cell_faces(6), ff, fff
        real(kind=8) :: maxCoords(3), minCoords(3), cellPtss(8, 3)
		logical :: pt_exists
		integer(kind=8) :: pt_id
		print *, "b3"
!		allocate(first5_refs(5,3))
!		first5_refs = 0.0_8 !!!!!!!!!!
        print*, 'Reading mesh: ', polyMeshPath  
        print *, "b3"
        nPoints = size(points_tempMesh, 1)
        nFaces = size(owner_tempMesh)
        nCells = maxval(owner_tempMesh)
        nBoundaries = size(tempMesh%boundaryNames)     
!        allocate(mesh%faces(nFaces, 2))
!        allocate(fAVecs_mesh(nFaces, 3), mesh%fCenters(nFaces, 3))
        
!        allocate(cellFaceCount(nCells))

!        allocate(mesh%cells(nCells, max_faces_per_cell))
        print *, "b3"
        cells_mesh = 0
        print *, "b3"
!        allocate(mesh%cVols(nCells), mesh%cCenters(nCells, 3), mesh%cellSizes(nCells, 3))
        ! 计算每个面的面积向量和几何中心
!        allocate(facePts(size(tempMesh%faces%faces(1)%points), 3))
        do f = 1, nFaces
!            allocate(facePts(size(tempMesh%faces%faces(f)%points), DIM_3))
            do i = 1, 4
                pt_id = facePoints_tempMesh(f, i)
                facePts(i,:) = points_tempMesh(pt_id,:)
            end do
!            print *, "b4"
            ! 保存前5个面的第一个顶点坐标作为参考!!!!!!
!			if (f <= 5) then
!				first5_refs(f,:) = facePts(1,:)
!			end if!!!!!!!!!!!!!!!!
!			print *, "b4"
			! 调用面计算子程序，传递前5个面的参考特征
!			call faceAreaCentroid(facePts, first5_refs, mesh%fAVecs(f,:), mesh%fCenters(f,:))!!!!!!!
            call faceAreaCentroid(facePts, fAVecs_mesh(f,:), fCenters_mesh(f,:))
        end do
!        deallocate(facePts)
!        deallocate(first5_refs)!!!!!!!!!!!!!
        print *, "b4"
        cellFaceCount = 0
        do f = 1, nFaces
            c = owner_tempMesh(f)
            cellFaceCount(c) = cellFaceCount(c) + 1
            cells_mesh(c, cellFaceCount(c)) = f
            faces_mesh(f, 1) = c
            faces_mesh(f, 2) = -1
            
            if (allocated(neighbour_tempMesh) .and. f <= size(neighbour_tempMesh)) then
                c = neighbour_tempMesh(f)
                cellFaceCount(c) = cellFaceCount(c) + 1
                cells_mesh(c, cellFaceCount(c)) = f
                faces_mesh(f, 2) = c
            end if
        end do
!        deallocate(cellFaceCount)
 !       allocate(fCs(6, 3), cell_fAVecs(6, 3))
!        allocate(cellPts(0, 3))
!       allocate(cellPtss(1, 3))
        
        print *, "b4"
        !遍历网格中的每一个单元（cell），收集该单元所有面的所有顶点坐标，去除重复的顶点，最终得到每个单元的唯一顶点集合。
        do c = 1, nCells
            cellPtss = 0.0_8
   			count_cellPtss = 0
            faceCount = 0
            ! 步骤1：统计当前单元的有效面数量（faceCount）
            do i = 1, size(cells_mesh, 2)
                if (cells_mesh(c, i) == 0) exit
                faceCount = faceCount + 1
            end do  
            ! 步骤2：提取当前单元的所有有效面ID（存入cell_faces）                    
            cell_faces = 0  ! 初始化
			cell_faces(1:faceCount) = cells_mesh(c, 1:faceCount)
			! 步骤3：遍历当前单元的每一个面
            do i = 1, faceCount
            	
!                ff = cell_faces(i)
!                print *, "ff", ff
                
                fCs(i,:) = fCenters_mesh(cell_faces(i),:)
                cell_fAVecs(i,:) = fAVecs_mesh(cell_faces(i),:)
                ! 步骤4：遍历当前面的每一个顶点（j=1~4，对应面的4个顶点）
                do j = 1, 4
                    pt_id = facePoints_tempMesh(cell_faces(i), j)
                    pt_exists = .false.   ! 初始化为没找到
        
					! 步骤5：检查该顶点是否已存在于cellPtss中（去重核心逻辑）
					do k = 1, count_cellPtss  ! 只遍历有效行（1~count_cellPtss）
						! 浮点坐标比较：误差小于1e-12即认为是同一个点（和原函数一致）
						if (all(abs(cellPtss(k,:) - points_tempMesh(pt_id,:)) < 1d-12)) then
						    pt_exists = .true.  ! 找到该点，标记为true
						    exit  ! 找到后直接退出循环，无需继续比较
						end if
					end do
					! 步骤6：如果顶点不存在，就添加到cellPtss中）
					if (.not. pt_exists) then
						count_cellPtss = count_cellPtss + 1  ! 有效行数+1
						cellPtss(count_cellPtss, :) = points_tempMesh(pt_id,:) 
					end if
                end do
            end do

			! 仅对前5个单元输出输入变量的前5行
			if(c <= 5) then  ! 新增：限制仅前5个单元执行以下输出
				
				! 1. 输出输入变量cellPtss（当前单元的前5行顶点坐标）
				print *, "调用cellVolCentroid前 - 单元", c, "的输入cellPtss（前5行）："
				do k = 1, 5  ! 每个单元内最多输出5行
					print "(A, I0, A, F0.6, A, F0.6, A, F0.6)", "  第", k, "行: (", &
						  cellPtss(k,1), ", ", cellPtss(k,2), ", ", cellPtss(k,3)
				end do

				! 2. 输出输入变量cell_fAVecs（当前单元的前5行面面积向量）
				print *, "调用cellVolCentroid前 - 单元", c, "的输入cell_fAVecs（前5行）："
				do k = 1, 5
					print "(A, I0, A, F0.6, A, F0.6, A, F0.6)", "  第", k, "行: (", &
						  cell_fAVecs(k,1), ", ", cell_fAVecs(k,2), ", ", cell_fAVecs(k,3)
				end do

				! 3. 输出输入变量fCs（当前单元的前5行面中心坐标）
				print *, "调用cellVolCentroid前 - 单元", c, "的输入fCs（前5行）："
				do k = 1, 5
					print "(A, I0, A, F0.6, A, F0.6, A, F0.6)", "  第", k, "行: (", &
						  fCs(k,1), ", ", fCs(k,2), ", ", fCs(k,3)
				end do
				
				print *, "----------------------------------------"  ! 可选：分隔不同单元的输入信息

			end if  ! 新增：结束前5个单元的判
			! 对每个单元，利用其所有面的信息计算体积和中心
            call cellVolCentroid(cellPtss, cell_fAVecs, fCs, cVols_mesh(c), cCenters_mesh(c,:))
            ! 4. 输出输出变量mesh%cVols(c)（仅前5个单元）
			if(c <= 5) then  ! 新增：仅前5个单元输出
				print *, "调用cellVolCentroid后 - 输出mesh%cVols(c)（单元", c, "体积）："
				print "(A, F0.6)", "  单元体积: ", cVols_mesh(c)
			end if  ! 新增：结束前5个单元的判断

			! 5. 输出输出变量mesh%cCenters(c,:)（仅前5个单元）
			if(c <= 5) then  ! 新增：仅前5个单元输出
				print *, "调用cellVolCentroid后 - 输出mesh%cCenters(c,:)（单元", c, "中心）："
				print "(A, F0.6, A, F0.6, A, F0.6)", "  单元中心: (", &
					  cCenters_mesh(c,1), ", ", cCenters_mesh(c,2), ", ", cCenters_mesh(c,3)
				print *, "------------------------"  ! 可选：添加分隔线，区分不同单元输出
			end if  ! 新增：结束前5个单元的判断
        end do
!        deallocate(fCs, cell_fAVecs, cellPts, cellPtss)
        ! 生成 boundaryFaces 数组
		! 恢复原始分配方式
		! 生成 boundaryFaces 数组时，先全部初始化为0
!		allocate(mesh%boundaryFaces(nBoundaries, maxval(tempMesh%boundaryNumFaces)))
		boundaryFaces_mesh = 0  ! 关键：将所有元素初始化为0

		do b = 1, nBoundaries
			startF = tempMesh%boundaryStartFaces(b)
			endF = startF + tempMesh%boundaryNumFaces(b) - 1
			do f = startF, endF
				boundaryFaces_mesh(b, f - startF + 1) = f  ! 填充有效面索引（无效位置保持0）
			end do
		end do
        ! 计算每个单元在 x、y、z 方向上的尺寸（用边界包围盒法）
        
        do c = 1, nCells
            maxCoords = -huge(1.0_real64)
            minCoords = huge(1.0_real64)
            cell_faces(1:size(cells_mesh, 2)) = cells_mesh(c, 1:size(cells_mesh, 2))
            do i = 1, size(cells_mesh, 2)
!                fff = cell_faces(i)
                if (cell_faces(i) == 0) exit
                do j = 1, 4
                    pt_id = facePoints_tempMesh(cell_faces(i), j)
                    do k = 1, 3
                        maxCoords(k) = max(maxCoords(k), points_tempMesh(pt_id, k))
                        minCoords(k) = min(minCoords(k), points_tempMesh(pt_id, k))
                    end do
                end do
            end do
            cellSizes_mesh(c,:) = maxCoords - minCoords
        end do
        
!        call deallocateMeshData(tempMesh)
    end subroutine OpenFOAMMesh
    subroutine extend_array(old_arr, new_row, new_arr)
		implicit none
		real*8, allocatable, intent(inout) :: old_arr(:,:)  ! 原数组
		real*8, intent(in) :: new_row(3)                     ! 新增行
		real*8, allocatable, intent(out) :: new_arr(:,:)     ! 扩展后数组
		integer(kind=8) :: new_size

		! 计算新尺寸
		if (allocated(old_arr)) then
		    new_size = size(old_arr, 1) + 1
		    allocate(new_arr(new_size, 3))
		    new_arr(1:new_size-1, :) = old_arr  ! 复制旧数据
		else
		    new_size = 1
		    allocate(new_arr(new_size, 3))
		end if
		new_arr(new_size, :) = new_row  ! 添加新行
	end subroutine extend_array
    ! 读取OpenFOAM网格主函数
    subroutine readOpenFOAMMesh(polyMeshPath, point_update, tempMesh)
        character(len=*), intent(in) :: polyMeshPath
        real*8, allocatable, intent(in) :: point_update(:,:)  ! 输入：变形网格点
        type(MeshData), intent(inout) :: tempMesh
        character(len=256) :: pointsFilePath, facesFilePath, ownerFilePath, neighbourFilePath, boundaryFilePath

        pointsFilePath = trim(polyMeshPath) // "/points"
        facesFilePath = trim(polyMeshPath) // "/faces"
        ownerFilePath = trim(polyMeshPath) // "/owner"
        neighbourFilePath = trim(polyMeshPath) // "/neighbour"
        boundaryFilePath = trim(polyMeshPath) // "/boundary"
        
        points_tempMesh = point_update
        facePoints_tempMesh = readOFFacesFile(facesFilePath)
        owner_tempMesh = readOFOwnerFile(ownerFilePath)
        neighbour_tempMesh = readOFNeighbourFile(neighbourFilePath)
        call readOFBoundaryFile(boundaryFilePath, tempMesh%boundaryNames, tempMesh%boundaryNumFaces, tempMesh%boundaryStartFaces)

    end subroutine readOpenFOAMMesh
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! 主函数：从OpenFOAM网格路径读取并处理单元点索引
    subroutine OpenFOAMMesh_findCellPts(polyMeshPath, pointLocations, cells, point_update)
		character(len=*), intent(in) :: polyMeshPath
		real*8, allocatable, intent(in) :: point_update(:,:)
		real(kind=8), allocatable, intent(out) :: pointLocations(:,:)
		type(CelllArray), intent(out) :: cells
		!$AD NONDIFF  ! 明确该变量为非微分
		
		! 主程序变量（包含pointIndicesByFace，避免全局变量）
		type(MeshData) :: tempMesh  ! 网格数据
		integer(kind=8) :: nCells, nFaces, nBoundaries
		integer(kind=8) :: i, j, fi, nFacesCell, quadFaceCount
		! 核心：点索引数据（对应Julia的pointIndicesByFace，在主程序内定义）
		integer(kind=8), allocatable :: pointIndicesByFace(:,:)  ! (面索引, 点序号)
		integer(kind=8) :: ierr, len_common
		integer(kind=8), allocatable :: commonPoints(:)		! 1. 读取网格数据
!		allocate(points_tempMesh(size(point_update,1), size(point_update,2)))
!		allocate(owner_tempMesh(0), owner_tempMesh(0))
!        allocate(tempMesh%boundaryNames(0), tempMesh%boundaryNumFaces(0), meshDataa%boundaryStartFaces(0))
		call readOpenFOAMMesh(polyMeshPath, point_update, tempMesh)
		
		! 2. 提取pointIndicesByFace（关键：在主程序内定义，供内置函数使用）
		! 假设从meshDataa中提取点索引（根据实际数据结构调整）
		nFaces = size(facePoints_tempMesh, 1)
		allocate(pointIndicesByFace(nFaces, 4))  ! 假设最多4个点（四边形）
		do i = 1, nFaces
		    if (size(facePoints_tempMesh(i,:)) /= 4) then
				print *, "Error: 面", i, "顶点数不为4（六面体要求），实际为", size(facePoints_tempMesh(i,:))
				stop
			end if
		    pointIndicesByFace(i, 1:4) = int(facePoints_tempMesh(i, 1:4))
		end do
		
		! 3. 初始化单元数组
		nCells = maxval(owner_tempMesh)
		allocate(cells%cells(nCells))
		do i = 1, nCells
		    allocate(cells%cells(i)%faceIndices(0))
		    allocate(cells%cells(i)%pointIndices(0))
		end do
		
		! 4. 添加面索引（调用外部的addCellFaceIndices）
		call addCellFaceIndices(owner_tempMesh, cells, nCells)
		call addCellFaceIndices(neighbour_tempMesh, cells, nCells)
		
		! 5. 处理单元点索引（调用内置函数）
		do i = 1, nCells
		    nFacesCell = size(cells%cells(i)%faceIndices)
		    select case(nFacesCell)
		        case(4)
		            call populatePointIndices_Tet(cells%cells(i))
		        case(5)
		            quadFaceCount = 0
		            do j = 1, nFacesCell
		                fi = cells%cells(i)%faceIndices(j)
		                if (size(pointIndicesByFace(fi,:)) == 4) quadFaceCount = quadFaceCount + 1
		            end do
		            if (quadFaceCount == 1) then
		                call populatePointIndices_Pyramid(cells%cells(i))
		            else if (quadFaceCount == 3) then
		                call populatePointIndices_Wedge(cells%cells(i))
		            else
		                print *, "Error: 无法识别的单元类型（单元", i, "）"
		                stop
		            end if
		        case(6)
		        ! 在case(6)中，单元1处理前添加调试
					if (i == 1) then
						print *, "调试：单元1的6个面索引：", cells%cells(i)%faceIndices
						do j = 1, 6
							fi = cells%cells(i)%faceIndices(j)
							print *, "面", fi, "的顶点：", pointIndicesByFace(fi, :)
						end do
						! 检查面1与其他面的公共点（基准面f1应与4个侧面相邻）
						print *, "面1（f1）与其他面的公共点数量："
						do j = 2, 6
							fi = cells%cells(i)%faceIndices(j)
							call get_intersection(pointIndicesByFace(1, :), pointIndicesByFace(fi, :), commonPoints, len_common)
							print *, "面1与面", fi, "：", len_common, "个公共点"
							deallocate(commonPoints)
						end do
					end if
		        ! 新增：检查单元所有面是否为4个顶点（六面体必需）
					do j = 1, size(cells%cells(i)%faceIndices)
						fi = cells%cells(i)%faceIndices(j)
						if (size(pointIndicesByFace(fi, :)) /= 4) then
							print *, "Error: 六面体单元", i, "的面", fi, "顶点数不为4（实际为", size(pointIndicesByFace(fi, :)), "）"
							stop
						end if
					end do
					! 新增：检查面数是否确实为6（防止数据错误）
					if (size(cells%cells(i)%faceIndices) /= 6) then
						print *, "Error: 单元", i, "声明为六面体但实际面数为", size(cells%cells(i)%faceIndices)
						stop
					end if
					call populatePointIndices_Hex(cells%cells(i), ierr)
					if (ierr /= 0) then
						! 增强错误信息：输出单元1的面索引，方便排查
						print *, "Error: 六面体单元", i, "顶点排序失败，错误代码", ierr
						print *, "单元", i, "的面索引列表：", cells%cells(i)%faceIndices
						stop
					end if
		        case default
		            print *, "Error: 不支持的单元面数（单元", i, "，面数", nFacesCell, "）"
		            stop
		    end select
		end do
		
		! 6. 输出点坐标
		if (allocated(pointLocations)) deallocate(pointLocations)
		allocate(pointLocations(size(points_tempMesh,1), 3))
		pointLocations = points_tempMesh


	! 内置所有点索引处理函数（通过contains关联，共享主程序变量）
	contains
		! 四面体单元点索引组织
		subroutine populatePointIndices_Tet(cell)
		    type(Celll), intent(inout) :: cell
		    ! 直接调用内置的addAllNewPoints（无需传递参数，共享pointIndicesByFace）
		    !call addAllNewPoints(cell, cell%faceIndices(1))
		    !call addAllNewPoints(cell, cell%faceIndices(2))
		    print *, "Aa"
		end subroutine populatePointIndices_Tet

		! 金字塔单元点索引组织
		subroutine populatePointIndices_Pyramid(cell)
		    type(Celll), intent(inout) :: cell
		    integer(kind=8) :: i, quadFaceIdx, otherFaceIndex
		    quadFaceIdx = -1
		    do i = 1, size(cell%faceIndices)
		        ! 使用主程序的pointIndicesByFace判断面类型
		        if (size(pointIndicesByFace(cell%faceIndices(i),:)) == 4) then
		            quadFaceIdx = cell%faceIndices(i)
		            exit
		        end if
		    end do
		    if (quadFaceIdx /= -1) then
		        !call addAllNewPoints(cell, quadFaceIdx)
		        ! 选择另一个面（与Julia逻辑一致）
		        if (quadFaceIdx /= cell%faceIndices(1)) then
		            otherFaceIndex = cell%faceIndices(1)
		        else
		            otherFaceIndex = cell%faceIndices(2)
		        end if
		        !call addAllNewPoints(cell, otherFaceIndex)
		    end if
		    print *, "Bb"
		end subroutine populatePointIndices_Pyramid		
		! 棱柱（楔形）单元：两个三角形底面，其余三个是连接面		
		subroutine populatePointIndices_Wedge(cell)
			type(Celll), intent(inout) :: cell
			integer(kind=8) :: i, j, t1, t2, nUnused, k
			integer(kind=8), allocatable :: unusedFaces(:)
			
			nUnused = size(cell%faceIndices)
			allocate(unusedFaces(nUnused))
			unusedFaces = cell%faceIndices
			
			t1 = -1
			do j = 1, 2
				do i = 1, nUnused
				    if (count(pointIndicesByFace(unusedFaces(i),:) > 0) == 3) then
				        if (t1 == -1) then
				            t1 = unusedFaces(i)
				            ! 修正：使用 k 作为循环变量
				            do k = i, nUnused-1
				                unusedFaces(k) = unusedFaces(k+1)
				            end do
				            nUnused = nUnused - 1
				            exit
				        else
				            t2 = unusedFaces(i)
				            ! 修正：使用 k 作为循环变量
				            do k = i, nUnused-1
				                unusedFaces(k) = unusedFaces(k+1)
				            end do
				            nUnused = nUnused - 1
				            exit
				        end if
				    end if
				end do
			end do
			
			!if (t1 /= -1) call addAllNewPoints(cell, t1)
			!if (t2 /= -1) call addAllNewPoints(cell, t2)
			
			do i = 1, nUnused
				!call addAllNewPoints(cell, unusedFaces(i))
			end do
			
			deallocate(unusedFaces)
		end subroutine populatePointIndices_Wedge
	    ! 六面体（立方体）单元：按照 VTK 的顺序组织 8 个顶点
        ! 原定义错误：嵌套contains且参数多余，修改为内部子程序（直接访问主程序变量）
		subroutine populatePointIndices_Hex(cell, ierr)
			type(Celll), intent(inout) :: cell
			integer(kind=8), intent(out) :: ierr
			integer(kind=8) :: f1, f3, lastFace, noResult
			integer(kind=8), allocatable :: f1Points(:), fiPoints(:)
			integer(kind=8), allocatable :: unusedFaces(:)
			integer(kind=8) :: i, fi, lastFaceIdx
			logical :: found

			ierr = 0
			found = .false.

			! 1. 复制所有面索引到unusedFaces
			allocate(unusedFaces(size(cell%faceIndices)))
			unusedFaces = cell%faceIndices
			!!print *, "初始面索引: ", unusedFaces

			! 2. 取第一个面作为f1（修复核心错误）
			f1 = unusedFaces(1)
			call remove_element(unusedFaces, 1_int64)  ! 移除第一个元素
			!!print *, "基准面f1 = 面", f1, "，剩余面: ", unusedFaces

			! 3. 获取f1的顶点
			allocate(f1Points(4))
			f1Points = pointIndicesByFace(f1, :)
			!!print *, "f1的顶点: ", f1Points

			! 4. 找对面（与f1无公共点的面，应为面2）
			do i = 1, size(unusedFaces)
				fi = unusedFaces(i)
				allocate(fiPoints(4), source=0_int64)
				fiPoints = pointIndicesByFace(fi, :)
				if (disjoint(f1Points, fiPoints)) then
				    !!print *, "找到对面: 面", fi
				    call remove_element(unusedFaces, i)
				    found = .true.
				    deallocate(fiPoints)
				    exit
				end if
				deallocate(fiPoints)
			end do
			if (.not. found) then
				ierr = 1
				!!print *, "错误代码1: 未找到对面"
				return
			end if

			! 检查剩余面是否为4个（关键验证）
			if (size(unusedFaces) /= 4) then
				print *, "错误：剩余面数应为4，实际为", size(unusedFaces)
				ierr = 1
				return
			end if
			!!print *, "删除对面后剩余面（4个侧面）: ", unusedFaces

			! 5. 初始化顶点数组
			if (allocated(cell%pointIndices)) deallocate(cell%pointIndices)
			allocate(cell%pointIndices(8))
			cell%pointIndices(1:4) = f1Points(1:4)
			cell%pointIndices(5:8) = [0,0,0,0]
			!!print *, "初始顶点索引: ", cell%pointIndices

			! 6. 取第一个剩余面作为f3
			f3 = unusedFaces(1)
			call remove_element(unusedFaces, 1_int64)
			!!print *, "处理面f3 = 面", f3, "，剩余面: ", unusedFaces

			! 7. 第一次调用addEdges
			lastFace = addEdges(cell%pointIndices, pointIndicesByFace(f3, :), f1Points, unusedFaces)
			!!print *, "第一次addEdges后顶点: ", cell%pointIndices
			!!print *, "addEdges返回的lastFace: ", lastFace

			! 8. 检查lastFace有效性
			if (lastFace < 1 .or. lastFace > size(unusedFaces)) then
				ierr = 3
				print *, "错误代码3: lastFace无效"
				return
			end if

			! 9. 处理最后一个面
			lastFaceIdx = unusedFaces(lastFace)
			call remove_element(unusedFaces, lastFace)
			!!print *, "处理最后一个面: 面", lastFaceIdx

			! 10. 第二次调用addEdges
			noResult = addEdges(cell%pointIndices, pointIndicesByFace(lastFaceIdx, :), f1Points, unusedFaces)
			!!print *, "第二次addEdges后顶点: ", cell%pointIndices

			! 11. 最终检查
			if (noResult /= -1) then
				ierr = 4
				return
			end if

			deallocate(f1Points, unusedFaces)
		end subroutine populatePointIndices_Hex

		function addEdges(cellPoints, facePoints, endFacePoints, unusedFaces) result(oppositeFaceIndex)
			integer(kind=8), intent(inout) :: cellPoints(:)
			integer(kind=8), intent(in) :: facePoints(:), endFacePoints(:)
			integer(kind=8), allocatable, intent(in) :: unusedFaces(:)
			integer(kind=8) :: oppositeFaceIndex
			integer(kind=8) :: i, fi, p1, p2, pos, len_common
			integer(kind=8), allocatable :: fiPoints(:), commonPoints(:)
			logical :: p1_in_end
			integer(kind=8) :: processedCount, expectedCount
			logical :: has_opposite  ! 新增：标记是否存在对面

			oppositeFaceIndex = -1
			processedCount = 0
			has_opposite = .false.  ! 初始化为无对面

			do i = 1, size(unusedFaces)
				fi = unusedFaces(i)
				if (allocated(fiPoints)) deallocate(fiPoints)
				allocate(fiPoints(size(pointIndicesByFace(fi, :))))
				fiPoints = pointIndicesByFace(fi, :)
				call get_intersection(facePoints, fiPoints, commonPoints, len_common)
				!print *, "addEdges: 面", fi, "与当前面公共点数量：", len_common

				if (len_common == 2) then
				    processedCount = processedCount + 1
				    ! 顶点填充逻辑（不变）
				    p1_in_end = any(endFacePoints == commonPoints(1))
				    if (p1_in_end) then
				        p1 = commonPoints(1)
				        p2 = commonPoints(2)
				    else
				        p1 = commonPoints(2)
				        p2 = commonPoints(1)
				    end if
				    pos = findloc(endFacePoints, p1, dim=1)
				    if (pos > 0) cellPoints(pos + 4) = p2

				else if (len_common == 0) then
				    oppositeFaceIndex = i
				    has_opposite = .true.  ! 标记存在对面
				end if
				deallocate(commonPoints)
			end do
			if (allocated(fiPoints)) deallocate(fiPoints)

			! 关键修改：根据是否有对面计算预期值
			if (has_opposite) then
				expectedCount = size(unusedFaces) - 1  ! 有对面：总面数-1
			else
				expectedCount = size(unusedFaces)      ! 无对面：总面数
			end if

			! 检查是否匹配
			if (processedCount /= expectedCount) then
				!print *, "Error: addEdges处理面数异常，预期", expectedCount, "实际", processedCount
				oppositeFaceIndex = -3
			end if
		end function addEdges
			
		! 辅助函数：计算两个数组的交集
		subroutine get_intersection(a, b, intersection, len)
			integer(kind=8), intent(in) :: a(:), b(:)
			integer(kind=8), allocatable, intent(out) :: intersection(:)
			integer(kind=8), intent(out) :: len
			integer(kind=8) :: i, j, count
			integer(kind=8), allocatable :: temp(:)  ! 临时数组

			count = 0
			allocate(intersection(min(size(a), size(b))))

			do i = 1, size(a)
				do j = 1, size(b)
				    if (a(i) == b(j)) then
				        count = count + 1
				        intersection(count) = a(i)
				        exit
				    end if
				end do
			end do

			len = count
			if (count < size(intersection)) then
				! 修复：先复制到临时可分配数组，再移动
				allocate(temp(count))
				temp = intersection(1:count)
				call move_alloc(temp, intersection)  ! 现在第一个参数是可分配数组
!				!$AD call move_alloc_int(temp, intersection)  ! Tapenade 使用（intersection 是整数一维数组）
			end if
		end subroutine get_intersection

		! 辅助函数：检查元素是否在数组中
		logical function is_in(val, arr)
			integer(kind=8), intent(in) :: val, arr(:)
			integer(kind=8) :: i
			is_in = .false.
			do i = 1, size(arr)
				if (arr(i) == val) then
				    is_in = .true.
				    return
				end if
			end do
		end function is_in

		! 辅助函数：检查两个点集是否不相交
		! 辅助函数：检查两个点集是否不相交（同级）
		logical function disjoint(a, b)
			integer(kind=8), intent(in) :: a(:), b(:)
			integer(kind=8) :: i, j
			disjoint = .true.
			do i = 1, size(a)
			    do j = 1, size(b)
			        if (a(i) == b(j)) then
			            disjoint = .false.
			            return
			        end if
			    end do
			end do
		end function disjoint
		! 辅助子程序：从数组中移除指定索引的元素（同级）
		subroutine remove_element(arr, idx)
		    integer(kind=8), allocatable, intent(inout) :: arr(:)
		    integer(kind=8), intent(in) :: idx
		    integer(kind=8), allocatable :: temp(:)
		    integer(kind=8) :: i, n

		    n = size(arr)
		    if (idx < 1 .or. idx > n) return

		    allocate(temp(n-1))
		    do i = 1, idx-1
		        temp(i) = arr(i)
		    end do
		    do i = idx+1, n
		        temp(i-1) = arr(i)
		    end do
		    !TAPENADE NOCHECKPOINT  
		    call move_alloc(temp, arr)
!		    !$AD call move_alloc_int(temp, arr)  ! Tapenade 使用（arr 是整数一维数组）
		end subroutine remove_element		
		subroutine swap(a, b)
			integer(kind=8), intent(inout) :: a, b
			integer(kind=8) :: temp
			!$AD NONDIFF  ! 替换原!$AD noalias a, b
			temp = a; a = b; b = temp
		end subroutine swap
		subroutine addCellFaceIndices(adjacentCells, cells, nCells)
			integer(kind=8), intent(in) :: adjacentCells(:)      ! 面所属单元索引
			integer(kind=8), intent(in) :: nCells                ! 总单元数
			type(CelllArray), intent(inout) :: cells    ! 单元数组
			!$AD NONDIFF  ! 明确该变量为非微分
			! 局部变量：统计用
			integer(kind=8) :: f, cellIdx, len, i, num_faces, j
			integer(kind=8) :: count_6_faces, tmp                    ! 6个面的单元数量
			integer(kind=8), allocatable :: other_counts(:), face_nums(:)  ! 其他面数统计
			integer(kind=8), allocatable :: temp(:)
			logical :: is_duplicate, found
			
			! 初始化统计变量
			count_6_faces = 0
			allocate(other_counts(0), face_nums(0))
			
			! 1. 为单元添加面索引（严格对应Julia逻辑）
			do f = 1, size(adjacentCells)
				cellIdx = adjacentCells(f)
				
				! 检查单元索引有效性（防止越界）
				if (cellIdx < 1 .or. cellIdx > nCells) cycle
				
				! 直接添加面索引（不检查重复，与Julia一致）
				len = size(cells%cells(cellIdx)%faceIndices)
				allocate(temp(len + 1))
				if (len > 0) temp(1:len) = cells%cells(cellIdx)%faceIndices
				temp(len + 1) = f
				! 在move_alloc前添加指令，跳过checkpointing
				!TAPENADE NOCHECKPOINT  
				call move_alloc(temp, cells%cells(cellIdx)%faceIndices)
!				!$AD call move_alloc_int(temp, cells%cells(cellIdx)%faceIndices)  ! Tapenade 使用（整数一维数组）
			end do
			
			! 2. 统计面数（与Julia逻辑完全一致）
			do cellIdx = 1, nCells
				num_faces = size(cells%cells(cellIdx)%faceIndices)
				
				! 统计6个面的单元
				if (num_faces == 6) then
				    count_6_faces = count_6_faces + 1
				else
				    ! 统计其他面数（动态数组模拟字典）
				    found = .false.
				    do i = 1, size(face_nums)
				        if (face_nums(i) == num_faces) then
				            other_counts(i) = other_counts(i) + 1
				            found = .true.
				            exit
				        end if
				    end do
				    
				    ! 新增面数类型
				    if (.not. found) then
				        allocate(temp(size(face_nums) + 1))
				        if (size(face_nums) > 0) temp(1:size(face_nums)) = face_nums
				        temp(size(face_nums) + 1) = num_faces
				        !TAPENADE NOCHECKPOINT
				        call move_alloc(temp, face_nums)
!				        !$AD call move_alloc_int(temp, face_nums)  ! Tapenade 使用（整数一维数组）
				        
				        allocate(temp(size(other_counts) + 1))
				        if (size(other_counts) > 0) temp(1:size(other_counts)) = other_counts
				        temp(size(other_counts) + 1) = 1
				        !TAPENADE NOCHECKPOINT
				        call move_alloc(temp, other_counts)
!				        !$AD call move_alloc_int(temp, face_nums)  ! Tapenade 使用（整数一维数组）
				    end if
				end if
			end do
			
			! 3. 输出统计结果（与Julia格式完全一致）
			print *, "包含6个面的单元数量: ", count_6_faces
			if (size(face_nums) > 0) then
				print *, "其他面数的单元统计:"
				! 排序（按面数升序）
				do i = 1, size(face_nums) - 1
				    do j = i + 1, size(face_nums)
				        if (face_nums(j) < face_nums(i)) then
				        	tmp = face_nums(i)
							face_nums(i) = face_nums(j)
							face_nums(j) = tmp

							! 交换 other_counts
							tmp = other_counts(i)
							other_counts(i) = other_counts(j)
							other_counts(j) = tmp
						end if
				    end do
				end do
				! 打印排序后的结果
				do i = 1, size(face_nums)
				    write(*, '(A, I0, A, I0, A)') "  ", face_nums(i), " 个面的单元: ", other_counts(i), " 个"
				end do
			else
				print *, "所有单元都包含6个面"
			end if
			
			! 释放临时数组
			deallocate(other_counts, face_nums)
		end subroutine addCellFaceIndices		
	end	subroutine OpenFOAMMesh_findCellPts
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! 读取点文件
    subroutine readOFPointsFile(filePath) 
        character(len=*), intent(in) :: filePath
!        real(kind=8), allocatable :: points(:,:)
        character(len=256), allocatable :: lines(:)
        integer(kind=8) :: startLine, pCount, i, nLines, iostat, file_unit
        character(len=256) :: line, bracketsRemoved
        real(kind=8) :: coords(3)
        
        file_unit = get_free_unit()
        open(unit=file_unit, file=filePath, status='old', action='read', iostat=iostat)
        
        nLines = 0
        do
            read(file_unit, '(a)', iostat=iostat)
            if (iostat /= 0) exit
            nLines = nLines + 1
        end do
        rewind(file_unit)
        
        allocate(lines(nLines))
        do i = 1, nLines
            read(file_unit, '(a)') lines(i)
        end do
        close(file_unit)
        
        call OFFile_FindNItems(lines, startLine, pCount)

        
!        allocate(points(pCount, 3))
        do i = 1, pCount
            line = trim(lines(startLine + i - 1))
            bracketsRemoved = line(2:len_trim(line)-1)
            read(bracketsRemoved, *) coords
            points_tempMesh(i,:) = coords
        end do
        deallocate(lines)
    end subroutine readOFPointsFile

    ! 读取面文件
    function readOFFacesFile(filePath) result(facePoints)
        implicit none
        integer(kind=8), parameter :: int64 = 8  ! 显式定义int64为8字节整数（等价于kind=8）        
        character(len=*), intent(in) :: filePath
        integer(kind=8) :: facePoints(74151,4)
!        type(FaceType), allocatable :: tmp_faces(:)
!        character(len=1000), allocatable :: lines_Faces(:)
        integer(kind=8) :: startLine, fCount, i, nLines, iostat, file_unit, bracketL, bracketR
        character(len=1000) :: line
!        integer(kind=8), allocatable :: pts(:)
        integer(kind=8) :: pt_count
!        integer(kind=8), allocatable :: pts_Faces(:)  
        PRINT *, "A"
        file_unit = get_free_unit()
        open(unit=file_unit, file=filePath, status='old', action='read', iostat=iostat)


        nLines = 0
        do
            read(file_unit, '(a)', iostat=iostat)
            if (iostat /= 0) exit
            nLines = nLines + 1
        end do
        rewind(file_unit)
        
!        allocate(lines_Faces(nLines))
        do i = 1, nLines
            read(file_unit, '(a)') lines_Faces(i)
        end do
        close(file_unit)

        call OFFile_FindNItems(lines_Faces, startLine, fCount)


!        allocate(tmp_faces(fCount))

        ! 循环外一次性分配
!		allocate(pts_Faces(max_pt_count), source=0_int64)
!		allocate(facePoints(fCount, 4)) 
		! 解析面数据
		PRINT *, "A"
		do i = 1, fCount
			line = trim(lines_Faces(startLine + i - 1))
			bracketL = index(line, '(')
			bracketR = index(line, ')')
			read(line(1:bracketL-1), *) pt_count
			! 读取括号内的点索引
			read(line(bracketL+1:bracketR-1), *) pts_Faces	
			! 0基转1基
			pts_Faces = pts_Faces + 1

			! 给面结构分配空间并赋值
!			allocate(tmp_faces(i)%points(pt_count))

			facePoints(i, :) = pts_Faces
		end do
		PRINT *, "A"
!		deallocate(pts_Faces)
!        deallocate(lines_Faces, tmp_faces)
    end function readOFFacesFile

    ! 读取owner文件
    function readOFOwnerFile(filePath) result(owner)
        implicit none
        character(len=*), intent(in) :: filePath
        integer(kind=8), allocatable :: owner(:)
!        character(len=256), allocatable :: lines(:)
        integer(kind=8) :: startLine, oCount, i, nLines, iostat, file_unit
        integer(kind=8) :: val
        
        file_unit = get_free_unit()
        open(unit=file_unit, file=filePath, status='old', action='read', iostat=iostat)
        if (iostat /= 0) then
            print *, "Error: Can't open owner file: ", trim(filePath)
            allocate(owner(0))
            return
        end if
        
        nLines = 0
        do
            read(file_unit, '(a)', iostat=iostat)
            if (iostat /= 0) exit
            nLines = nLines + 1
        end do
        rewind(file_unit)
        
!        allocate(lines(nLines))
        do i = 1, nLines
            read(file_unit, '(a)') lines_owner(i)
        end do
        close(file_unit)
        
        call OFFile_FindNItems(lines_owner, startLine, oCount)
        if (oCount <= 0) then
            print *, "Error: Invalid owner count in ", trim(filePath)
            deallocate(lines_owner)
            allocate(owner(0))
            return
        end if
        
        allocate(owner(oCount))
        do i = 1, oCount
            read(lines_owner(startLine + i - 1), *) val
            owner(i) = int(val) + 1
        end do
!        deallocate(lines)
    end function readOFOwnerFile

    ! 读取neighbour文件
    function readOFNeighbourFile(filePath) result(neighbour)
		implicit none
		character(len=*), intent(in) :: filePath  ! 输入文件路径
		integer(kind=8), allocatable :: neighbour(:)      ! 输出邻居单元索引数组
!		character(len=256), allocatable :: lines_neighbour(:)  ! 存储文件所有行
		integer(kind=8) :: startLine, nCount  ! 数据起始行、数据总数
		integer(kind=8) :: i, nLines, iostat, file_unit  ! 循环变量、总行数、I/O状态、文件单元
		integer(kind=8) :: val  ! 临时存储读取的大整数（兼容OpenFOAM格式）

		! 1. 初始化：不预先分配数组（与owner处理一致）
		if (allocated(neighbour)) deallocate(neighbour)

		! 2. 打开文件（与owner逻辑完全一致）
		file_unit = get_free_unit()
		open(unit=file_unit, file=filePath, status='old', action='read', iostat=iostat)
		if (iostat /= 0) then
		    print *, "Error: Can't open neighbour file: ", trim(filePath)  ! 错误信息与owner格式统一
		    allocate(neighbour(0))  ! 错误时返回空数组
		    return
		end if

		! 3. 读取文件总行数（与owner逻辑一致）
		nLines = 0
		do
		    read(file_unit, '(a)', iostat=iostat)  ! 逐行读取计数
		    if (iostat /= 0) exit
		    nLines = nLines + 1
		end do
		rewind(file_unit)  ! 重置文件指针

		! 4. 存储所有行内容（与owner逻辑一致）
!		allocate(lines(nLines))  ! 分配存储所有行的数组
		do i = 1, nLines
		    read(file_unit, '(a)') lines_neighbour(i)  ! 逐行读取内容
		end do
		close(file_unit)  ! 关闭文件（后续使用lines处理）

		! 5. 查找数据起始行和数量（核心：与owner调用同一函数）
		call OFFile_FindNItems(lines_neighbour, startLine, nCount)

		! 6. 检查数据数量有效性（与owner逻辑一致）
		if (nCount <= 0) then
		    print *, "Error: Invalid neighbour count in ", trim(filePath)  ! 错误信息格式统一
		    allocate(neighbour(0))  ! 无效数量时返回空数组
		    return
		end if

		! 7. 分配邻居数组并读取数据（与owner完全一致）
		allocate(neighbour(nCount))  ! 确认数量有效后再分配
		do i = 1, nCount
		    ! 严格按照owner的读取方式：startLine + i - 1索引
		    read(lines_neighbour(startLine + i - 1), *) val  ! 读取原始值（OpenFOAM格式）
		    neighbour(i) = int(val) + 1  ! 转换为Fortran索引（+1逻辑与owner一致）
		end do

		! 8. 释放临时资源（与owner一致）
!		deallocate(lines_neighbour)

	end function readOFNeighbourFile
    ! 读取边界文件
    subroutine readOFBoundaryFile(filePath, boundaryNames, boundaryNumFacess, boundaryStartFacess)
        implicit none
        character(len=*), intent(in) :: filePath
        character(len=100), allocatable, intent(inout) :: boundaryNames(:)
        integer(kind=8), allocatable, intent(inout) :: boundaryNumFacess(:), boundaryStartFacess(:)
        integer(kind=8) :: i, nLines, iostat, startLine, bCount, file_unit
        integer(kind=8) :: bNameLine, bNFacesLine, bStartFaceLine, pos
!        character(len=256), allocatable :: bLines(:)
        character(len=256) :: lineStr, dummy
        
        file_unit = get_free_unit()
        open(unit=file_unit, file=filePath, status='old', action='read', iostat=iostat)
!        if (iostat /= 0) then
!            allocate(boundaryNames(0), boundaryNumFacess(0), boundaryStartFacess(0))
!            return
!        end if
        print *, "A0"
        nLines = 0
        do
            read(file_unit, '(a)', iostat=iostat)
            if (iostat /= 0) exit
            nLines = nLines + 1
        end do
        rewind(file_unit)
        
!        allocate(bLines(nLines))
!        do i = 1, nLines
!            read(file_unit, '(a)') bLines(i)
!        end do
        close(file_unit)
        
        call OFFile_FindNItems(bLines, startLine, bCount)
!        if (bCount <= 0) then
!            deallocate(bLines)
!            allocate(boundaryNames(0), boundaryNumFacess(0), boundaryStartFacess(0))
!            return
!        end if
        
        allocate(boundaryNumFacess(bCount), boundaryStartFacess(bCount))
		allocate(boundaryNames(bCount))
		print *, "A0"
        do i = 1, bCount
            bNameLine = findInLines("{", bLines, startLine) - 1
            boundaryNames(i) = trim(adjustl(bLines(bNameLine)))
            
            bNFacesLine = findInLines("nFaces", bLines, startLine)
            pos = index(bLines(bNFacesLine), "nFaces") + 6
            print *, "A1"
            read(bLines(bNFacesLine)(pos:), *) boundaryNumFacess(i)
            print *, "A2"
            bStartFaceLine = findInLines("startFace", bLines, startLine)
            pos = index(bLines(bStartFaceLine), "startFace") + 9
            read(bLines(bStartFaceLine), *) dummy, boundaryStartFacess(i)
			boundaryStartFacess(i) = boundaryStartFacess(i) + 1
			!boundaryEndFaces(i) = boundaryStartFacess(i) + boundaryNumFacess(i) - 1
            
            startLine = findInLines("}", bLines, startLine) + 1
        end do
!        deallocate(bLines)
    end subroutine readOFBoundaryFile

    ! 辅助函数：查找数据数量和起始行
    subroutine OFFile_FindNItems(fileLines, startLine, itemCount)
        implicit none
        character(len=*), intent(in) :: fileLines(:)
        integer(kind=8), intent(out) :: startLine, itemCount
        integer(kind=8) :: i, pos, iostat
        character(len=256) :: line
        
        itemCount = 0
        startLine = 0
        do i = 1, size(fileLines)
            line = trim(fileLines(i))
            pos = index(line, "nPoints") + index(line, "nFaces") + index(line, "size")
            if (pos == 0) cycle
            
            read(line(pos:), *, iostat=iostat) itemCount
            if (itemCount > 0) then
                startLine = i + 2
                return
            end if
        end do
        
        do i = 1, size(fileLines)
            if (isNumber(trim(fileLines(i)))) then
                read(fileLines(i), *, iostat=iostat) itemCount
                if (itemCount > 0) then
                    startLine = i + 2
                    return
                end if
            end if
        end do
        
        do i = 1, size(fileLines)
            if (index(fileLines(i), "(") > 0) then
                startLine = i + 1
                exit
            end if
        end do
        itemCount = 0
        do i = startLine, size(fileLines)
            if (index(fileLines(i), ")") > 0) exit
            itemCount = itemCount + 1
        end do
    end subroutine OFFile_FindNItems

    ! 辅助函数：判断是否为数字
    logical function isNumber(str)
        character(len=*), intent(in) :: str
        real(kind=8) :: num
        integer(kind=8) :: iostat
        isNumber = .false.
        read(str, *, iostat=iostat) num
        if (iostat == 0) isNumber = .true.
    end function isNumber

    ! 辅助函数：查找包含子串的行
    integer(kind=8) function findInLines(substr, lines, startLine)
        character(len=*), intent(in) :: substr, lines(:)
        integer(kind=8), intent(in) :: startLine
        integer(kind=8) :: i
        findInLines = 0
        do i = startLine, size(lines)
            if (index(lines(i), substr) > 0) then
                findInLines = i
                return
            end if
        end do
    end function findInLines

    ! 辅助函数：获取空闲单元号
    integer(kind=8) function get_free_unit()
        implicit none
        integer(kind=8) :: i, iostat
        logical :: opened
        do i = 10, 999
            inquire(unit=i, opened=opened, iostat=iostat)
            if (.not. opened .and. iostat == 0) then
                get_free_unit = i
                return
            end if
        end do
    end function get_free_unit

    ! 其他必要函数（简化实现）
    logical function isInArray(pts, pt)
        real(kind=8), intent(in) :: pts(:,:), pt(3)
        integer(kind=8) :: i
        isInArray = .false.
        do i = 1, size(pts,1)
            if (all(abs(pts(i,:)-pt) < 1d-12)) then
                isInArray = .true.
                return
            end if
        end do
    end function isInArray

    subroutine deallocateMeshData(md)
        type(MeshData), intent(inout) :: md
        integer(kind=8) :: i
!        if (allocated(md%points)) deallocate(md%points)
!        if (allocated(md%owner)) deallocate(md%owner)
!        if (allocated(md%neighbour)) deallocate(md%neighbour)
        if (allocated(md%boundaryNames)) deallocate(md%boundaryNames)
        if (allocated(md%boundaryNumFaces)) deallocate(md%boundaryNumFaces)
        if (allocated(md%boundaryStartFaces)) deallocate(md%boundaryStartFaces)
!        if (allocated(md%facePoints)) deallocate(md%facePoints)
    end subroutine deallocateMeshData

!7.numerics.jl!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!######################### Gradient Computation #######################
!1.梯度计算：（1）leastSqGrad: 最小二乘法梯度计算（未完成）（2）greenGaussGrad: Green-Gauss梯度计算方法（基于面通量积分）    2.面值插值：（1）linInterp_3D: 线性插值到面中心（2）maxInterp: 取相邻单元最大值作为面值（3）faceDeltas: 计算面两侧变量变化量（用于激波捕捉）    3.时间步进相关：（1）decodeSolution_3D: 从守恒变量计算原始变量和单元通量（2）integrateFluxes_unstructured3D: 通过面通量积分计算残差
	subroutine greenGaussGrad(matrix, valuesAtFaces, grad, faceVals_d)
		implicit none
!		type(Meshh), intent(in) :: mesh  ! 网格信息
		real(kind=8), intent(in) :: matrix(:,:)  ! 单元变量矩阵（行：单元，列：变量）
		logical, intent(in), optional :: valuesAtFaces  ! 是否直接提供面值
		real(kind=8), allocatable, intent(inout) :: grad(:,:,:)  ! 梯度矩阵（单元, 变量, 空间维度）
		integer(kind=8) :: nCells, nFaces, nBoundaries, nBdryFaces, meshInfo(4)
		integer(kind=8) :: f, v, d, ownerCell, neighbourCell, nVars
		real(kind=8), intent(inout) :: faceVals_d(:,:)  ! 面插值结果

		! 解包网格信息：单元格数、面数、边界数、边界面数
		call unstructuredMeshInfo(meshInfo)
		nCells = meshInfo(1)
		nFaces = meshInfo(2)
		nBoundaries = meshInfo(3)
		nBdryFaces = meshInfo(4)
		nVars = size(matrix, 2)  ! 变量数（如速度、压力等）
        
!        allocate(faceVals_d(nFaces, 5))  ! 维度：(面数, 变量数)
!    	faceVals_d = 0.0_8  ! 初始化避免垃圾值
		! 若面值未直接提供，则通过线性插值计算面值
		if (present(valuesAtFaces)) then
		    if (valuesAtFaces) then
		        faceVals_d = matrix  ! 直接使用输入矩阵作为面值
		    else
		        call linInterp_3D(matrix, faceVals_d)  ! 调用线性插值函数计算面值
		    end if
		else
		    call linInterp_3D(matrix, faceVals_d)  ! 默认插值计算面值
		end if

		! 初始化梯度矩阵：[单元格, 变量, 空间维度(x/y/z)]
!		allocate(grad(nCells, nVars, 3))
		grad = 0.0_8

		! 遍历所有面，积分通量以计算梯度
		do f = 1, nFaces
		    ownerCell = faces_mesh(f,1)       ! 面所属的主单元格（owner cell）
		    neighbourCell = faces_mesh(f,2)    ! 面相邻的邻接单元格（neighbour cell，边界面为-1）

		    do v = 1, nVars
		        do d = 1, 3  ! 空间维度（x/y/z方向）
		            ! 主单元格：累加面通量（面面积矢量 × 面值）
		            grad(ownerCell, v, d) = grad(ownerCell, v, d) + fAVecs_mesh(f,d) * faceVals_d(f,v)

		            ! 邻接单元格（非边界面时）：减去面通量（通量方向相反）
		            if (neighbourCell > -1) then
		                grad(neighbourCell, v, d) = grad(neighbourCell, v, d) - fAVecs_mesh(f,d) * faceVals_d(f,v)
		            end if
		        end do
		    end do
		end do

		! 将积分结果除以单元格体积，得到梯度（单位体积的通量）
		do d = 1, 3
		    do v = 1, nVars
		        do f = 1, nCells  ! f循环实际为单元索引c
		            grad(f, v, d) = grad(f, v, d) / cVols_mesh(f)
		        end do
		    end do
		end do

!		deallocate(faceVals_d)
	end subroutine greenGaussGrad
	subroutine greenGaussGradd(matrix, valuesAtFaces, grad)
		implicit none
!		type(Meshh), intent(in) :: mesh  ! 网格信息
		real(kind=8), intent(in) :: matrix(:,:)  ! 单元变量矩阵（行：单元，列：变量）
		logical, intent(in), optional :: valuesAtFaces  ! 是否直接提供面值
		real(kind=8), allocatable, intent(inout) :: grad(:,:,:)  ! 梯度矩阵（单元, 变量, 空间维度）
		integer(kind=8) :: nCells, nFaces, nBoundaries, nBdryFaces, meshInfo(4)
		integer(kind=8) :: f, v, d, ownerCell, neighbourCell, nVars
!		real(kind=8), allocatable :: faceVals(:,:)  ! 面插值结果

		! 解包网格信息：单元格数、面数、边界数、边界面数
		call unstructuredMeshInfo(meshInfo)
		nCells = meshInfo(1)
		nFaces = meshInfo(2)
		nBoundaries = meshInfo(3)
		nBdryFaces = meshInfo(4)
		nVars = size(matrix, 2)  ! 变量数（如速度、压力等）
        
!        allocate(faceVals_dd(nFaces, 5))  ! 维度：(面数, 变量数)
    	faceVals_dd = 0.0_8  ! 初始化避免垃圾值
		! 若面值未直接提供，则通过线性插值计算面值
		if (present(valuesAtFaces)) then
		    if (valuesAtFaces) then
		        faceVals_dd = matrix  ! 直接使用输入矩阵作为面值
		    else
		        call linInterp_3D(matrix, faceVals_dd)  ! 调用线性插值函数计算面值
		    end if
		else
		    call linInterp_3D(matrix, faceVals_dd)  ! 默认插值计算面值
		end if

		! 初始化梯度矩阵：[单元格, 变量, 空间维度(x/y/z)]
!		allocate(grad(nCells, nVars, 3))
		grad = 0.0_8

		! 遍历所有面，积分通量以计算梯度
		do f = 1, nFaces
		    ownerCell = faces_mesh(f,1)       ! 面所属的主单元格（owner cell）
		    neighbourCell = faces_mesh(f,2)    ! 面相邻的邻接单元格（neighbour cell，边界面为-1）

		    do v = 1, nVars
		        do d = 1, 3  ! 空间维度（x/y/z方向）
		            ! 主单元格：累加面通量（面面积矢量 × 面值）
		            grad(ownerCell, v, d) = grad(ownerCell, v, d) + fAVecs_mesh(f,d) * faceVals_dd(f,v)

		            ! 邻接单元格（非边界面时）：减去面通量（通量方向相反）
		            if (neighbourCell > -1) then
		                grad(neighbourCell, v, d) = grad(neighbourCell, v, d) - fAVecs_mesh(f,d) * faceVals_dd(f,v)
		            end if
		        end do
		    end do
		end do

		! 将积分结果除以单元格体积，得到梯度（单位体积的通量）
		do d = 1, 3
		    do v = 1, nVars
		        do f = 1, nCells  ! f循环实际为单元索引c
		            grad(f, v, d) = grad(f, v, d) / cVols_mesh(f)
		        end do
		    end do
		end do

!		deallocate(faceVals_dd)
	end subroutine greenGaussGradd
	!####################### 面值插值 ######################
	!#=
	!   对所有内部面进行插值（边界面需单独处理）
	!   输入矩阵格式：
	!    行：单元格，列：物理变量（如u, v, w, p）
	!   输出矩阵格式：
	!   行：面，列：物理变量（插值后的面值）
	!=#
	subroutine linInterp_3D(matrix, faceVals)
    	implicit none
!		type(Meshh), intent(in) :: mesh  ! 网格信息
		real(kind=8), intent(in) :: matrix(:,:)  ! 单元变量矩阵（行：单元，列：变量）
		real(kind=8), intent(inout) :: faceVals(:,:)  ! 预分配的面值矩阵（输入输出）
		integer(kind=8) :: nCells, nFaces, nBoundaries, nBdryFaces, meshInfo(4)
		integer(kind=8) :: nVars, f, v, i, c1, c2, row, nVars_matrix, nVars_face
		real(kind=8) :: c1Dist, c2Dist, totalDist

		! 解包网格信息
		call unstructuredMeshInfo(meshInfo)
		nCells = meshInfo(1)
		nFaces = meshInfo(2)
		nBoundaries = meshInfo(3)
		nBdryFaces = meshInfo(4)
		nVars = size(matrix, 2)  ! 变量数
		nVars_matrix = size(matrix, 2)  ! 实际列数
    	nVars_face = size(faceVals, 2)
        ! LTSeuler给的输入cellprimitvits(matrix)没问题
		! 验证faceVals维度是否正确（避免越界）
		if (size(faceVals, 1) < nFaces .or. size(faceVals, 2) < nVars) then
		    print *, "Error: faceVals dimensions mismatch!"
		    return
		end if
		! 遍历内部面（排除边界面nBdryFaces）
		do f = 1, nFaces - nBdryFaces 
		    c1 = faces_mesh(f,1)  ! 主单元格
		    c2 = faces_mesh(f,2)  ! 邻接单元格

		    ! 计算单元格中心到面中心的距离平方（欧氏距离）
		    c1Dist = 0.0_8
		    c2Dist = 0.0_8
		    do i = 1, 3  ! x/y/z方向
		        c1Dist = c1Dist + (cCenters_mesh(c1,i) - fCenters_mesh(f,i))**2
		        c2Dist = c2Dist + (cCenters_mesh(c2,i) - fCenters_mesh(f,i))**2
		    end do
		    totalDist = c1Dist + c2Dist  ! 总距离
		    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        
		    IF (f <= 3) THEN
		   	  PRINT*, '===== Face/Cell No.', f, ' ====='
!			  PRINT*, 'c1distd   : ', c1distd
			  PRINT*, 'c1dist    : ', c1dist
!			  PRINT*, 'c2distd   : ', c2distd
			  PRINT*, 'c2dist    : ', c2dist
!			  PRINT*, 'totaldistd: ', totaldistd
			  PRINT*, 'totaldist : ', totaldist
			  PRINT*, '---------------------------------'
			END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		    ! 线性插值：根据距离权重计算面值
		    do v = 1, nVars
		        faceVals(f, v) = matrix(c1, v) * (c2Dist / totalDist) + matrix(c2, v) * (c1Dist / totalDist)
		    end do
		end do
	end subroutine linInterp_3D
	! 与linInterp_3D类似，但取相邻单元格的最大值作为面值
	subroutine maxInterp(sjj, rjj, faceVals)
!		type(Meshh), intent(in) :: mesh
		real(kind=8), intent(in) :: sjj(:), rjj(:)
		real(kind=8), allocatable, intent(inout) :: faceVals(:,:)
		
		integer(kind=8) :: nCells, nFaces, nBoundaries, nBdryFaces, meshInfo(4)
		integer(kind=8) :: f, c1, c2
		integer(kind=8) :: nVars = 2  ! 固定为2个变量（sjj和rj）
		
		! 获取网格信息
		call unstructuredMeshInfo(meshInfo)
		nCells = meshInfo(1)
		nFaces = meshInfo(2)
		nBoundaries = meshInfo(3)
		nBdryFaces = meshInfo(4)
		
		! 分配内存给输出数组
!		allocate(faceVals(nFaces, nVars))
		faceVals = 0.0_8  ! 初始化为0
		
		! 遍历内部面
		do f = 1, nFaces - nBdryFaces
		    c1 = faces_mesh(f, 1)  ! 主单元格
		    c2 = faces_mesh(f, 2)  ! 邻居单元格
		    
		    ! 对每个变量取最大值
		    faceVals(f, 1) = max(sjj(c1), sjj(c2))  ! 第一列存储sj的插值结果
		    faceVals(f, 2) = max(rjj(c1), rjj(c2))  ! 第二列存储rjj的插值结果
		end do		
	end subroutine maxInterp
	! 计算面两侧单元格的变量差值（邻接值 - 主值）
	subroutine faceDeltas(deltas)
		implicit none
!		type(Meshh), intent(in) :: mesh  ! 网格信息
!		type(SolutionState), intent(in) :: sln  ! 流场解状态
		real(kind=8), allocatable, intent(inout) :: deltas(:,:)  ! 面差值矩阵（行：面，列：变量）
		integer(kind=8) :: nCells, nFaces, nBoundaries, nBdryFaces, meshInfo(4)
		integer(kind=8) :: nVars, f, v, ownerCell, neighbourCell

		! 解包网格信息
		call unstructuredMeshInfo(meshInfo)
		nCells = meshInfo(1)
		nFaces = meshInfo(2)
		nBoundaries = meshInfo(3)
		nBdryFaces = meshInfo(4)
		nVars = size(cellState_sln, 2)  ! 单元格状态变量数

		! 初始化面差值矩阵
!		allocate(deltas(nFaces, nVars))
		deltas = 0.0_8

		! 遍历内部面（除去边界面）
		do f = 1, nFaces - nBdryFaces
		    ownerCell = faces_mesh(f,1)    ! 主单元格
		    neighbourCell = faces_mesh(f,2)  ! 邻接单元格

		    do v = 1, nVars
		        ! 计算邻接单元格与主单元格的变量差值
		        deltas(f, v) = cellState_sln(neighbourCell, v) - cellState_sln(ownerCell, v)
		    end do
		end do
	end subroutine faceDeltas
	!######################### TimeStepping #################################
	!######################### 时间步进相关 #########################
	!# 从单元格状态向量解码为原始变量和通量
	subroutine decodeSolution_3D(fluid) 
		implicit none
!		type(SolutionState), intent(inout) :: sln
		type(Fluidd), intent(in) :: fluid

		integer(kind=8) :: nCells, i

		nCells = size(cellState_sln, 1)
		! 解码原始变量并计算通量
		do i = 1, nCells		   
		    call decodePrimitives3D(cellPrimitives_sln(i,:), cellState_sln(i,:), fluid)
		      ! 新增：输出前50行sln%cellState和sln%cellPrimitives（在calculateFluxes3D前）
		    ! 输出前50行cellState和cellPrimitives（增大宽度避免溢出）
!        if (i <= 50) then
!            if (i == 1) then
!                print *, new_line('a')//"=== 前50行sln%cellState和sln%cellPrimitives ==="
!                print *, "格式：单元索引 | cellState(5列) | cellPrimitives(5列)"
!            end if
            ! 关键修改：将F12.6改为F18.6（宽度从12增至18）
!            write(*, '(I5, 10F18.6)') i, &
!                sln%cellState(i,:), &       ! cellState(5列)
!                sln%cellPrimitives(i,:)     ! cellPrimitives(5列)
!       end if

        	call calculateFluxes3D(cellFluxes_sln(i,:), cellPrimitives_sln(i,:), cellState_sln(i,:))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   			  
		end do
	end subroutine decodeSolution_3D

	!#=
	!   通过面通量计算单元格的通量残差（用于时间步进）
	!   残差 = -∇·F（F为通量矢量）
	!=#
	subroutine integrateFluxes_unstructured3D(fluxResiduals)
		implicit none
!		type(Meshh), intent(in) :: mesh
!		type(SolutionState), intent(inout) :: sln
		real(kind=8), intent(inout) :: fluxResiduals(:,:)
		integer(kind=8) :: nFaces, nBoundaries, nBdryFaces, meshInfo(4), nCells, nVars
		integer(kind=8) :: f, v, ownerCell, neighbourCell, i1, i2, i, maxOutput, j
		real(kind=8) :: flow
		INTEGER :: y1d, min1d  ! 局部变量：meshd%favecs的行数、实际输出数量
		INTEGER :: y2d, min2d  ! 局部变量：meshd%cvols的行数、实际输出数量
		
		! 1. 正确获取网格和变量信息（不依赖未分配的数组）
		call unstructuredMeshInfo(meshInfo)
		nCells = meshInfo(1)         ! 单元总数（有效，非0）
		nFaces = meshInfo(2)
		nBoundaries = meshInfo(3)
		nBdryFaces = meshInfo(4)
		nVars = size(cellState_sln, 2)  
		! 输出前5个mesh.fAVecs和mesh.cVols
		write(*,*) "First 5 mesh.fAVecs:"
		do f = 1, 5
		    write(*,'(A,I3,A,3(E12.6,1X))') "Face ", f, ": ", fAVecs_mesh(f,1), fAVecs_mesh(f,2), fAVecs_mesh(f,3)
		end do
		
		write(*,*) "First 5 mesh.cVols:"
		do f = 1, 5
		    write(*,'(A,I3,A,E12.6)') "Cell ", f, ": ", cVols_mesh(f)
		end do
		! ======================================================================
		! 输出前5个 meshd%favecs		
!		WRITE(*, *) NEW_LINE('a')//'First 5 meshd%fAVecs (differential):'
!		y1d = SIZE(meshd%favecs, 1)  ! 获取meshd%favecs的总行数（面数）
!		min1d = 5          ! 最多输出5个面，避免越界
!		DO f=1,min1d
!		  WRITE(*, '(A,I3,A,3(E12.6,1X))') 'Diff Face ', f, ': ', &
!			meshd%favecs(f, 1), meshd%favecs(f, 2), meshd%favecs(f, 3)
!		END DO

		! 输出前5个 meshd%cvols		
!		WRITE(*, *) 'First 5 meshd%cVols (differential):'
!		y2d = SIZE(meshd%cvols, 1)  ! 获取meshd%cvols的总行数（单元数）
!		min2d = 5         ! 最多输出5个单元，避免越界
!		DO f=1,min2d
!		  WRITE(*, '(A,I3,A,E12.6)') 'Diff Cell ', f, ': ', meshd%cvols(f)
!		END DO
		! ======================================================================
		! 2. 确保 sln%fluxResiduals 已分配（关键修正）
!		if (.not. allocated(sln%fluxResiduals) .or. &
!		    size(sln%fluxResiduals,1)/=nCells .or. &
!		    size(sln%fluxResiduals,2)/=nVars) then
		    ! 若未分配或维度不匹配，重新分配
!		    if (allocated(sln%fluxResiduals)) deallocate(sln%fluxResiduals)
!		    allocate(sln%fluxResiduals(nCells, nVars))
!		end if
		
		! 3. 重置残差矩阵
		fluxResiduals_sln = 0.0_8
		
		! 4. 通量积分（核心逻辑）
		do f = 1, nFaces
		    ownerCell = faces_mesh(f,1)
		    neighbourCell = faces_mesh(f,2)
		    
		    do v = 1, nVars
		        i1 = (v-1)*3 + 1
		        i2 = i1 + 2
		        call dot_product_real(faceFluxes_sln(f,i1:i2), fAVecs_mesh(f,:), flow, 3_8)  ! 通量点积
		        ! 添加输出 - 只输出前10个面的flow值
				if (f <= 10) then
					write(*,'(A,I3,A,I3,A,E12.5)') "Face ", f, " Var ", v, " Flow = ", flow
				end if       
		        fluxResiduals_sln(ownerCell, v) = fluxResiduals_sln(ownerCell, v) - flow
		   		        
		        ! 邻接单元格残差（非边界）
		        if (neighbourCell > 0 .and. neighbourCell <= nCells) then  ! 避免边界的-1
		            fluxResiduals_sln(neighbourCell, v) = fluxResiduals_sln(neighbourCell, v) + flow
		        end if
		    end do
		end do
		maxOutput = 15  ! 确保不超过数组大小
		write(*,*) "3D:Flux residuals for first 15 cells:"
		do i = 1, maxOutput
		    write(*, '(A,I3,A,*(E12.6,1X))') "Cell ", i, ": ", (fluxResiduals_sln(i,j), j=1,nVars)
		end do
		! 5. 除以单元格体积（单位体积残差）
		do v = 1, nVars
		    do f = 1, nCells  ! f 实际是单元索引（建议重命名为 c 更清晰）
		        fluxResiduals_sln(f, v) = fluxResiduals_sln(f, v) / cVols_mesh(f)
		    end do
		end do
		
		! 6. 分配并返回 fluxResiduals
		!allocate(fluxResiduals(nCells, nVars))
		fluxResiduals = fluxResiduals_sln
	end subroutine integrateFluxes_unstructured3D
! 8.output.jl
	! 将单元原始变量写入重启文件!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine writeRestartFile(cellPrimitives, path)
		real(kind=8), intent(in) :: cellPrimitives(:,:)
		character(len=*), intent(in), optional :: path
		character(len=256) :: filePath
		
		if (present(path)) then
			filePath = path
		else
			filePath = "FvCFDRestart.txt"
		end if
		
		open(unit=10, file=filePath, status='replace', action='write')
		write(10, *) cellPrimitives
		close(10)
	end subroutine writeRestartFile
	
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    ! 2. 修正调用writeVTKFile时的参数传递（确保points为二维数组）
	subroutine outputVTK(meshPath, cellPrimitives, fileName, point_update)
		! （参数声明不变）
		character(len=*), intent(in) :: meshPath
		real(kind=8), intent(in) :: cellPrimitives(:,:)
		character(len=*), intent(in), optional :: fileName
		real*8, allocatable, intent(in) :: point_update(:,:)
				
		! （局部变量声明不变）
		real(kind=8), allocatable :: pointLocations(:,:)
		type(CelllArray) :: cells
		!$AD NONDIFF  ! 明确该变量为非微分
		character(len=256) :: vtkFileName
		integer(kind=8) :: nCells, nPoints, i
		
		! 获取点坐标和单元数据（确保pointLocations是二维数组）
		call OpenFOAMMesh_findCellPts(meshPath, pointLocations, cells, point_update)
		nCells = size(cells%cells)
		nPoints = size(pointLocations, 1)  ! 原维度：nPoints×3
		
		! 转置点坐标（转为3×nPoints，与Julia一致）
		call transposePoints(pointLocations)  ! 转置后pointLocations为3×nPoints
		
		! 修正调用方式：传递转置后的pointLocations（二维数组）
		if (present(fileName)) then
		    vtkFileName = trim(fileName) // ".vtk"
		else
		    vtkFileName = "solution.vtk"
		end if
		call writeVTKFile(vtkFileName, pointLocations, cells, cellPrimitives, &
		                 nPoints, nCells, [1_int64,3_int64,5_int64,10_int64,14_int64,13_int64,-1_int64,12_int64])  ! vtkCellTypes数组

		! （其余代码不变）
	end subroutine outputVTK
   	subroutine writeOutput(cellPrimitives, restartFile, meshPath, createRestartFile, createVTKOutput, point_update)
		! 输入参数（与Julia函数完全对应）
		real(kind=8), intent(in) :: cellPrimitives(:,:)  ! 单元原始数据
		character(len=*), intent(in) :: restartFile      ! 重启文件名
		character(len=*), intent(in) :: meshPath         ! 网格路径
		logical, intent(in) :: createRestartFile         ! 是否创建重启文件
		logical, intent(in) :: createVTKOutput           ! 是否创建VTK输出
		real*8, allocatable, intent(in) :: point_update(:,:) 
		! 局部变量
		integer(kind=8) :: unit, iostat
		character(len=256) :: solnName
		logical :: exists
		
		! 1. 处理重启文件（对应Julia的if createRestartFile）
		if (createRestartFile) then
		    print *, 'Writing Restart File: ', trim(restartFile)
		    call writeRestartFile(cellPrimitives, restartFile)  ! 调用重启文件写入子程序
		end if
		
		! 2. 处理VTK输出（核心修改：删除系统命令、文件列表、序号查找，改用固定文件名）
		if (createVTKOutput) then
		    solnName = 'solution_001.vtk'  ! 固定文件名，无需序号查找
		    print *, 'Writing VTK File: ', trim(solnName)
		    call outputVTK(meshPath, cellPrimitives, solnName, point_update)  ! 原调用逻辑不变
		end if

	end subroutine writeOutput
	subroutine writeVTKFile(fileName, points, cells, primitives, nPts, nCells, vtkTypes)
		character(len=*), intent(in) :: fileName
		real(kind=8), intent(in) :: points(:,:)  ! 3×nPts（x,y,z）
		type(CelllArray), intent(in) :: cells
		!$AD NONDIFF  ! 明确该变量为非微分
		real(kind=8), intent(in) :: primitives(:,:)  ! 流场数据：[nCells×5] = [P, T, Ux, Uy, Uz]
		integer(kind=8), intent(in) :: nPts, nCells
		integer(kind=8), intent(in) :: vtkTypes(8)
		integer(kind=8) :: fileUnit, i, j, n, cellType, ptIdx, connOffset
		real(kind=8) :: P, T, Ux, Uy, Uz
		! Fortran 95兼容：手动定义无穷大阈值（替代isinf）
		real(kind=8), parameter :: INF_THRESHOLD = 1.0d300  ! 超过此值视为无穷大

		! 2. 预检查：流场数据维度匹配
		if (size(primitives, 1) /= nCells .or. size(primitives, 2) < 5) then
		    print *, "Error: 流场数据维度不匹配（单元数=", nCells, ", 数据行数=", size(primitives,1), "）"
		    return
		end if

		! 3. 写入VTK文件
		fileUnit = get_free_unit()  ! 显式初始化
		open(unit=fileUnit, file=fileName, status='replace', action='write', form='formatted')

		! 3.1 VTK文件头
		write(fileUnit, '(A)') '# vtk DataFile Version 3.0'
		write(fileUnit, '(A)') 'Generated by BuFlowModule'
		write(fileUnit, '(A)') 'ASCII'
		write(fileUnit, '(A)') 'DATASET UNSTRUCTURED_GRID'

		! 3.2 点坐标
		write(fileUnit, '(A,I0,A)') 'POINTS ', nPts, ' double'
		do i = 1, nPts
		    write(fileUnit, '(3F16.8)') points(1,i), points(2,i), points(3,i)
		end do

		! 3.3 单元连接性
		connOffset = 0
		do i = 1, nCells
		    connOffset = connOffset + size(cells%cells(i)%pointIndices) + 1
		end do
		write(fileUnit, '(A,I0,A,I0)') 'CELLS ', nCells, ' ', connOffset
		do i = 1, nCells
		    n = size(cells%cells(i)%pointIndices)
		    write(fileUnit, '(I4)', advance='no') n
		    do j = 1, n
		        ptIdx = cells%cells(i)%pointIndices(j) - 1  ! 0-based索引
		        write(fileUnit, '(I6)', advance='no') ptIdx
		    end do
		    write(fileUnit, *)
		end do

		! 3.4 单元类型
		write(fileUnit, '(A,I0)') 'CELL_TYPES ', nCells
		do i = 1, nCells
		    n = size(cells%cells(i)%pointIndices)
		    cellType = merge(vtkTypes(n), 1_int64, n>=1 .and. n<=8)
		    write(fileUnit, '(I4)') cellType
		end do

		! 3.5 流场数据（兼容Fortran 95，无isinf）
		write(fileUnit, '(A,I0)') 'CELL_DATA ', nCells

		! 3.5.1 压力场P
		write(fileUnit, '(A)') 'SCALARS P double 1'
		write(fileUnit, '(A)') 'LOOKUP_TABLE default'
		do i = 1, nCells
		    P = primitives(i,1)
		    ! Fortran 95兼容：处理nan和无穷大（nan的特征是不等于自身）
!		    if (P /= P .or. abs(P) > INF_THRESHOLD) P = 0.0d0  ! P /= P 表示nan
		    write(fileUnit, '(F16.8)') P
		end do

		! 3.5.2 温度场T
		write(fileUnit, '(A)') 'SCALARS T double 1'
		write(fileUnit, '(A)') 'LOOKUP_TABLE default'
		do i = 1, nCells
		    T = primitives(i,2)
!		    if (T /= T .or. abs(T) > INF_THRESHOLD) T = 0.0d0
		    write(fileUnit, '(F16.8)') T
		end do

		! 3.5.3 速度场U
		write(fileUnit, '(A)') 'VECTORS U double'
		do i = 1, nCells
		    Ux = primitives(i,3)
		    Uy = primitives(i,4)
		    Uz = primitives(i,5)
		    ! 处理无效值
!		    if (Ux /= Ux .or. abs(Ux) > INF_THRESHOLD) Ux = 0.0d0
!		    if (Uy /= Uy .or. abs(Uy) > INF_THRESHOLD) Uy = 0.0d0
!		    if (Uz /= Uz .or. abs(Uz) > INF_THRESHOLD) Uz = 0.0d0
		    write(fileUnit, '(3F16.8)') Ux, Uy, Uz
		end do

		close(fileUnit)
	end subroutine writeVTKFile
	subroutine writeOutput_d(cellPrimitives, cellPrimitivesd, restartFile, meshPath, createRestartFile, createVTKOutput, point_update)
		! 输入参数（与Julia函数完全对应）
		real(kind=8), intent(in) :: cellPrimitives(:,:), cellPrimitivesd(:,:)  ! 单元原始数据
		character(len=*), intent(in) :: restartFile      ! 重启文件名
		character(len=*), intent(in) :: meshPath         ! 网格路径
		logical, intent(in) :: createRestartFile         ! 是否创建重启文件
		logical, intent(in) :: createVTKOutput           ! 是否创建VTK输出
		real*8, allocatable, intent(in) :: point_update(:,:) 
		! 局部变量
		integer(kind=8) :: unit, iostat
		character(len=256) :: solnName
		logical :: exists
		
		! 1. 处理重启文件（对应Julia的if createRestartFile）
		if (createRestartFile) then
		    print *, 'Writing Restart File: ', trim(restartFile)
		    call writeRestartFile(cellPrimitives, restartFile)  ! 调用重启文件写入子程序
		end if
		
		! 2. 处理VTK输出（核心修改：删除系统命令、文件列表、序号查找，改用固定文件名）
		if (createVTKOutput) then
		    solnName = 'solution_001.vtk'  ! 固定文件名，无需序号查找
		    print *, 'Writing VTK File: ', trim(solnName)
		    call outputVTK_d(meshPath, cellPrimitives, cellPrimitivesd, solnName, point_update)  ! 原调用逻辑不变
		end if

	end subroutine writeOutput_d
	! 2. 修正调用writeVTKFile时的参数传递（确保points为二维数组）
	subroutine outputVTK_d(meshPath, cellPrimitives, cellPrimitivesd, fileName, point_update)
		! （参数声明不变）
		character(len=*), intent(in) :: meshPath
		real(kind=8), intent(in) :: cellPrimitives(:,:), cellPrimitivesd(:,:)
		character(len=*), intent(in), optional :: fileName
		real*8, allocatable, intent(in) :: point_update(:,:)
				
		! （局部变量声明不变）
		real(kind=8), allocatable :: pointLocations(:,:)
		type(CelllArray) :: cells
		!$AD NONDIFF  ! 明确该变量为非微分
		character(len=256) :: vtkFileName
		integer(kind=8) :: nCells, nPoints, i
		
		! 获取点坐标和单元数据（确保pointLocations是二维数组）
		call OpenFOAMMesh_findCellPts(meshPath, pointLocations, cells, point_update)
		nCells = size(cells%cells)
		nPoints = size(pointLocations, 1)  ! 原维度：nPoints×3
		
		! 转置点坐标（转为3×nPoints，与Julia一致）
		call transposePoints(pointLocations)  ! 转置后pointLocations为3×nPoints
		
		! 修正调用方式：传递转置后的pointLocations（二维数组）
		if (present(fileName)) then
		    vtkFileName = trim(fileName) // ".vtk"
		else
		    vtkFileName = "solution.vtk"
		end if
		call writeVTKFile_d(vtkFileName, pointLocations, cells, cellPrimitives, cellPrimitivesd, &  ! 新增cellPrimitivesd
                     nPoints, nCells, [1_int64,3_int64,5_int64,10_int64,14_int64,13_int64,-1_int64,12_int64])  ! vtkCellTypes数组

		! （其余代码不变）
	end subroutine outputVTK_d
	subroutine writeVTKFile_d(fileName, points, cells, primitives, cellPrimitivesd, nPts, nCells, vtkTypes)
		character(len=*), intent(in) :: fileName
		real(kind=8), intent(in) :: points(:,:)  ! 3×nPts（x,y,z）
		type(CelllArray), intent(in) :: cells
		!$AD NONDIFF  ! 明确该变量为非微分
		real(kind=8), intent(in) :: primitives(:,:), cellPrimitivesd(:,:)  ! 流场数据：[nCells×5] = [P, T, Ux, Uy, Uz]
		integer(kind=8), intent(in) :: nPts, nCells
		integer(kind=8), intent(in) :: vtkTypes(8)
		integer(kind=8) :: fileUnit, i, j, n, cellType, ptIdx, connOffset
		real(kind=8) :: P, T, Ux, Uy, Uz, P_t, T_t, Ux_t, Uy_t, Uz_t
		! Fortran 95兼容：手动定义无穷大阈值（替代isinf）
		real(kind=8), parameter :: INF_THRESHOLD = 1.0d300  ! 超过此值视为无穷大

		! 2. 预检查：流场数据维度匹配
		if (size(primitives, 1) /= nCells .or. size(primitives, 2) < 5) then
		    print *, "Error: 流场数据维度不匹配（单元数=", nCells, ", 数据行数=", size(primitives,1), "）"
		    return
		end if

		! 3. 写入VTK文件
		fileUnit = get_free_unit()  ! 显式初始化
		open(unit=fileUnit, file=fileName, status='replace', action='write', form='formatted')

		! 3.1 VTK文件头
		write(fileUnit, '(A)') '# vtk DataFile Version 3.0'
		write(fileUnit, '(A)') 'Generated by BuFlowModule'
		write(fileUnit, '(A)') 'ASCII'
		write(fileUnit, '(A)') 'DATASET UNSTRUCTURED_GRID'

		! 3.2 点坐标
		write(fileUnit, '(A,I0,A)') 'POINTS ', nPts, ' double'
		do i = 1, nPts
		    write(fileUnit, '(3F16.8)') points(1,i), points(2,i), points(3,i)
		end do

		! 3.3 单元连接性
		connOffset = 0
		do i = 1, nCells
		    connOffset = connOffset + size(cells%cells(i)%pointIndices) + 1
		end do
		write(fileUnit, '(A,I0,A,I0)') 'CELLS ', nCells, ' ', connOffset
		do i = 1, nCells
		    n = size(cells%cells(i)%pointIndices)
		    write(fileUnit, '(I4)', advance='no') n
		    do j = 1, n
		        ptIdx = cells%cells(i)%pointIndices(j) - 1  ! 0-based索引
		        write(fileUnit, '(I6)', advance='no') ptIdx
		    end do
		    write(fileUnit, *)
		end do

		! 3.4 单元类型
		write(fileUnit, '(A,I0)') 'CELL_TYPES ', nCells
		do i = 1, nCells
		    n = size(cells%cells(i)%pointIndices)
		    cellType = merge(vtkTypes(n), 1_int64, n>=1 .and. n<=8)
		    write(fileUnit, '(I4)') cellType
		end do

		! 3.5 流场数据（兼容Fortran 95，无isinf）
		write(fileUnit, '(A,I0)') 'CELL_DATA ', nCells

		! 3.5.1 压力场P
		write(fileUnit, '(A)') 'SCALARS P double 1'
		write(fileUnit, '(A)') 'LOOKUP_TABLE default'
		do i = 1, nCells
		    P = primitives(i,1)
		    ! Fortran 95兼容：处理nan和无穷大（nan的特征是不等于自身）
!		    if (P /= P .or. abs(P) > INF_THRESHOLD) P = 0.0d0  ! P /= P 表示nan
		    write(fileUnit, '(F16.8)') P
		end do

		! 3.5.2 温度场T
		write(fileUnit, '(A)') 'SCALARS T double 1'
		write(fileUnit, '(A)') 'LOOKUP_TABLE default'
		do i = 1, nCells
		    T = primitives(i,2)
!		    if (T /= T .or. abs(T) > INF_THRESHOLD) T = 0.0d0
		    write(fileUnit, '(F16.8)') T
		end do

		! 3.5.3 速度场U
		write(fileUnit, '(A)') 'VECTORS U double'
		do i = 1, nCells
		    Ux = primitives(i,3)
		    Uy = primitives(i,4)
		    Uz = primitives(i,5)
		    ! 处理无效值
!		    if (Ux /= Ux .or. abs(Ux) > INF_THRESHOLD) Ux = 0.0d0
!		    if (Uy /= Uy .or. abs(Uy) > INF_THRESHOLD) Uy = 0.0d0
!		    if (Uz /= Uz .or. abs(Uz) > INF_THRESHOLD) Uz = 0.0d0
		    write(fileUnit, '(3F16.8)') Ux, Uy, Uz
		end do
		! 4.1 P_tagent（修改格式为科学计数法）
		WRITE(fileunit, '(A)') 'SCALARS P_tagent double 1'
		WRITE(fileunit, '(A)') 'LOOKUP_TABLE default'
		DO i=1,ncells
		    p_t = cellprimitivesd(i, 1)
		  ! 限制超大值（避免格式溢出）
!		    IF (ABS(p_t) > 1d20) p_t = SIGN(1d20, p_t)  ! 超过1e20的用1e20替代
		    WRITE(fileunit, '(ES16.8)') p_t  ! ES格式：科学计数法，适配大数值
		END DO

		! 4.2 T_tagent（同上述修改）
		WRITE(fileunit, '(A)') 'SCALARS T_tagent double 1'
		WRITE(fileunit, '(A)') 'LOOKUP_TABLE default'
		DO i=1,ncells
		    t_t = cellprimitivesd(i, 2)
!		    IF (ABS(t_t) > 1d20) t_t = SIGN(1d20, t_t)
		    WRITE(fileunit, '(ES16.8)') t_t
		END DO

		! 4.3 U_tagent 向量（同上述修改）
		WRITE(fileunit, '(A)') 'VECTORS U_tagent double'
		DO i=1,ncells
		    ux_t = cellprimitivesd(i, 3)
		    uy_t = cellprimitivesd(i, 4)
		    uz_t = cellprimitivesd(i, 5)
!		    IF (ABS(ux_t) > 1d20) ux_t = SIGN(1d20, ux_t)
!		    IF (ABS(uy_t) > 1d20) uy_t = SIGN(1d20, uy_t)
!		    IF (ABS(uz_t) > 1d20) uz_t = SIGN(1d20, uz_t)
		    WRITE(fileunit, '(3ES16.8)') ux_t, uy_t, uz_t  ! 向量也用科学计数法
		END DO

		close(fileUnit)
	end subroutine writeVTKFile_d
	! 3. 确认transposePoints子程序的正确性（保持不变）
	subroutine transposePoints(points)
		real(kind=8), allocatable, intent(inout) :: points(:,:)
		real(kind=8), allocatable :: temp(:,:)
		integer(kind=8) :: i, j, nPoints
		
		nPoints = size(points, 1)  ! 原维度：nPoints×3
		allocate(temp(3, nPoints))  ! 转置为3×nPoints
		
		do i = 1, nPoints
		    do j = 1, 3
		        temp(j, i) = points(i, j)  ! 原(i,j) → 新(j,i)
		    end do
		end do
		
		deallocate(points)
		call move_alloc(temp, points)  ! 转置后points为3×nPoints
!		!$AD call move_alloc_int(temp, face_nums)  ! Tapenade 使用（整数一维数组）
	end subroutine transposePoints

!9.timeDiscretizations.jl:更新解 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! 局部时间步长的欧拉方法.局部时间步长（Local Time Stepping, LTS）：为每个单元格独立计算时间步长
	subroutine LTSEuler(boundaryConditions, fluid, dt)
		implicit none
!		type(Meshh), intent(in) :: mesh
!		type(SolutionState), intent(inout) :: sln
		type(BoundaryCondition), intent(in) :: boundaryConditions(:)
		type(Fluidd), intent(in) :: fluid
		real(kind=8), intent(inout) :: dt(:)
		real(kind=8) :: targetCFL
!		real(kind=8), allocatable :: fluxResiduals(:,:), dt_val(:)
		integer(kind=8) :: maxOutput, nCells, nVars
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        integer(kind=8) :: i, j, nPts, nPointsPerCell, l  ! 循环变量、总点数、单元顶点数
		integer(kind=8) :: pointIndices(8)  ! 存储单元的顶点索引（假设最多8个顶点
        INTEGER(kind=8) :: print_dtd  ! 控制输出单元数量
        INTEGER(kind=8) :: print_slnd  ! 控制输出单元数量
        INTEGER(kind=8) :: print_resd  ! 控制输出单元数量

		targetCFL = dt(1)
		print *, 'targetCFL ', targetCFL
		nCells = size(cellState_sln, 1)   ! 获取单元数
		nVars = size(cellState_sln, 2)    ! 获取变量数
!		allocate(fluxResiduals(nCells, nVars), source=0.0_8)  ! 显式分配并初始化
!		allocate(dt(nCells), source=0.0d0)
!        allocate(dt_val(nCells), source=0.0_8) 
		! 计算残差
		call unstructured_JSTFlux(boundaryConditions, fluid, fluxResiduals)                            ! 显式维度参数2
!		nCells = size(sln%cellState, 1)
!		nVars = size(fluxResiduals, 2)
        maxOutput = 15  ! 确保不超过数组大小    
!		write(*,*) "Flux residuals for first 5 cells:"
!		do i = 1, maxOutput
!		  write(*, '(A,I3,A,*(E12.6,1X))') "Cell ", i, ": ", (fluxResiduals(i,j), j=1,nVars)
!		end do
        print *, 'www'
		! 计算时间步长
		call CFL(nCells, fluid, 1.0_8, dt)
		print *, 'www'		
		! 输出前5个dt值
!		write(*,*) "First 5 dt values:"
!		do i = 1, 5
!		    write(*,'(A,I3,A,E12.6)') "dt(", i, ") = ", dt(i)
!		end do
		dt = targetCFL / dt		
		! 输出前5个dt值
!		write(*,*) "First 5 dt values:"
!		do i = 1, 5
!		    write(*,'(A,I3,A,E12.6)') "dt(", i, ") = ", dt(i)
!		end do
        print *, "AAz"
		call smoothTimeStep(dt, 0.1_8, dt_val)
        print *, "AAz"
		! 输出前5个dt值
!		write(*,*) "First 5 dt values:"
!		do i = 1, 5
!		    write(*,'(A,I3,A,E12.6)') "dt(", i, ") = ", dt(i)
!		end do
        call smoothTimeStepp(dt_val, 0.1_8, dt)	
		! 用残差更新守恒变量（按广播匹配dt维度）
		do j = 1, nVars          ! 遍历每个变量
		    do i = 1, nCells     ! 遍历每个单元
		        cellState_sln(i,j) = cellState_sln(i,j) + fluxResiduals(i,j) * dt(i)
		    end do
		end do
		! 输出前6个cellState
		write(*,*) "First 6 cellState values:"
		do i = 1, 6
		    write(*,'(A,I3,A,*(E12.6,1X))') "Cell ", i, ": ", (cellState_sln(i,j), j=1,size(cellState_sln,2))
		end do
        print *, 'www'
		call decodeSolution_3D(fluid)
        ! 释放临时数组（可选，避免内存泄漏）
!    	deallocate(fluxResiduals, dt_val)
    	
        print *, 'www'
	end subroutine LTSEuler

	! 时间步长平滑函数（用于局部时间步长方法）.局部时间步长在空间上可能剧烈变化，导致相邻单元格步长差异过大，引发数值振荡。
	subroutine smoothTimeStep(dt, diffusionCoefficient, updated_dt)
		real(kind=8), intent(in) :: dt(:)
!		type(Meshh), intent(in) :: mesh
		real(kind=8), intent(in) :: diffusionCoefficient
!		real(kind=8), intent(out) :: updated_dt(:)
		real(kind=8), intent(inout) :: updated_dt(:)  ! 返回值，与输入的 dt 类型一致		
		integer(kind=8) :: nFaces, nBoundaries, nBdryFaces, f, ownerCell, neighbourCell, i, meshInfo(4), nCells
		real(kind=8) :: coeff, timeFlux, mag_val	
!		real(kind=8), allocatable :: timeFluxes(:), surfaceAreas(:)	
		
		coeff = diffusionCoefficient		
		! 获取网格信息
		call unstructuredMeshInfo(meshInfo)
		nCells = meshInfo(1)
		nFaces = meshInfo(2)
		nBoundaries = meshInfo(3)
		nBdryFaces = meshInfo(4)
!        allocate(updated_dt(nCells))
		! 初始化时间通量和表面积数组
!		allocate(timeFluxes(nCells))
!		allocate(surfaceAreas(nCells))
		timeFluxes = 0.0_8
		surfaceAreas = 0.0_8
		
		! 遍历所有内部面
		do f = 1, nFaces-nBdryFaces
		    ownerCell = faces_mesh(f,1)      ! 主单元
		    neighbourCell = faces_mesh(f,2)  ! 邻居单元
		    
		    ! 计算时间步长差异通量
		    call mag(fAVecs_mesh(f,:), mag_val)
			timeFlux = (dt(ownerCell) - dt(neighbourCell)) * mag_val
		    
		    ! 累加表面积（用于归一化）
		    surfaceAreas(ownerCell) = surfaceAreas(ownerCell) + mag_val
			surfaceAreas(neighbourCell) = surfaceAreas(neighbourCell) + mag_val
		    
		    ! 更新通量（守恒形式）
		    timeFluxes(ownerCell) = timeFluxes(ownerCell) - timeFlux
		    timeFluxes(neighbourCell) = timeFluxes(neighbourCell) + timeFlux
		end do

		! 应用扩散系数并归一化
		timeFluxes = timeFluxes * (coeff / surfaceAreas)

		! 确保只应用平滑（不增加时间步长）
		do i = 1, nCells
		    timeFluxes(i) = min(0.0_8, timeFluxes(i))
		end do

		! 更新时间步长
		updated_dt = dt + timeFluxes
		! 在子程序末尾添加释放
!		deallocate(timeFluxes, surfaceAreas)
	end subroutine smoothTimeStep
	subroutine smoothTimeStepp(dt, diffusionCoefficient, updated_dt)
		real(kind=8), intent(in) :: dt(:)
!		type(Meshh), intent(in) :: mesh
		real(kind=8), intent(in) :: diffusionCoefficient
!		real(kind=8), intent(out) :: updated_dt(:)
		real(kind=8), intent(inout) :: updated_dt(:)  ! 返回值，与输入的 dt 类型一致
		
		integer(kind=8) :: nFaces, nBoundaries, nBdryFaces, f, ownerCell, neighbourCell, i, meshInfo(4), nCells
		real(kind=8) :: coeff, timeFlux, mag_val
!		real(kind=8), allocatable :: timeFluxess(:), surfaceAreass(:)

		coeff = diffusionCoefficient
	! 获取网格信息
		call unstructuredMeshInfo(meshInfo)
		nCells = meshInfo(1)
		nFaces = meshInfo(2)
		nBoundaries = meshInfo(3)
		nBdryFaces = meshInfo(4)
!        allocate(updated_dt(nCells))
		! 初始化时间通量和表面积数组
!		allocate(timeFluxess(nCells))
!		allocate(surfaceAreass(nCells))
		timeFluxess = 0.0_8
		surfaceAreass = 0.0_8
		
		! 遍历所有内部面
		do f = 1, nFaces-nBdryFaces
		    ownerCell = faces_mesh(f,1)      ! 主单元
		    neighbourCell = faces_mesh(f,2)  ! 邻居单元
		    
		    ! 计算时间步长差异通量
		    call mag(fAVecs_mesh(f,:), mag_val)
			timeFlux = (dt(ownerCell) - dt(neighbourCell)) * mag_val
		    
		    ! 累加表面积（用于归一化）
		    surfaceAreass(ownerCell) = surfaceAreass(ownerCell) + mag_val
			surfaceAreass(neighbourCell) = surfaceAreass(neighbourCell) + mag_val
		    
		    ! 更新通量（守恒形式）
		    timeFluxess(ownerCell) = timeFluxess(ownerCell) - timeFlux
		    timeFluxess(neighbourCell) = timeFluxess(neighbourCell) + timeFlux
		end do

		! 应用扩散系数并归一化
		timeFluxess = timeFluxess * (coeff / surfaceAreass)

		! 确保只应用平滑（不增加时间步长）
		do i = 1, nCells
		    timeFluxess(i) = min(0.0_8, timeFluxess(i))
		end do

		! 更新时间步长
		updated_dt = dt + timeFluxess
		! 在子程序末尾添加释放
!		deallocate(timeFluxess, surfaceAreass)
	end subroutine smoothTimeStepp
!###########################################################################################
!10.vectorFunctions.jl#####################################################################
!######################### 向量函数 ########################
!#函数：点积，叉积，模，转化为单位向量
	! 情况1：向量 × 向量（均为1维数组）
	function dot_vec_vec(arg1, arg2) result(out)
		real(kind=8), intent(in) :: arg1(:)  ! 1维向量
		real(kind=8), intent(in) :: arg2(:)  ! 1维向量
		real(kind=8), allocatable :: out(:)
		integer(kind=8) :: i, n

		n = size(arg1)
		if (n /= size(arg2)) then
		    print *, "Error: Vectors must have same length"
		    stop
		end if
		allocate(out(1))
		out(1) = 0.0_8
		do i = 1, n
		    out(1) = out(1) + arg1(i) * arg2(i)  ! 合法的1维数组索引
		end do
	end function dot_vec_vec

	! 情况2：向量 × 二维矩阵（含向量矩阵）
	function dot_vec_mat(arg1, arg2) result(out)
		real(kind=8), intent(in) :: arg1(:)    ! 1维向量
		real(kind=8), intent(in) :: arg2(:,:)  ! 二维矩阵（含向量矩阵）
		real(kind=8), allocatable :: out(:), temp(:)
		integer(kind=8) :: i, n_rows, arg1_size, arg2_dim1, arg2_dim2

		arg1_size = size(arg1)
		arg2_dim1 = size(arg2, 1)  ! 矩阵行数
		arg2_dim2 = size(arg2, 2)  ! 矩阵列数

		! 子情况A：向量 × 向量矩阵（N×M矩阵，每行是M维向量）
		if (arg2_dim2 > 1) then
		    if (arg1_size /= arg2_dim2) then
		        print *, "Error: Vector length must match sub-vector length"
		        stop
		    end if
		    allocate(out(arg2_dim1))
		    do i = 1, arg2_dim1
		        temp = dot_vec_vec(arg1, arg2(i,:))   ! 先存储函数返回的数组
				out(i) = temp(1)   ! 调用向量×向量函数
		    end do

		! 子情况B：向量 × 二维矩阵（N×1矩阵，每行是标量）
		else if (arg2_dim2 == 1) then
		    if (arg1_size /= arg2_dim1) then
		        print *, "Error: Vector length must match matrix rows"
		        stop
		    end if
		    allocate(out(1))
		    out(1) = 0.0_8
		    do i = 1, arg2_dim1
		        out(1) = out(1) + arg1(i) * arg2(i,1)  ! 合法的二维数组索引
		    end do
		end if
	end function dot_vec_mat
	! 3D向量叉积（仅适用于3元素向量）
	function cross(v1, v2) result(cross_vec)
		implicit none
		real(kind=8), intent(in) :: v1(3), v2(3)  ! v1=(x1,y1,z1), v2=(x2,y2,z2)
		real(kind=8) :: cross_vec(3)
		
		! 正确的叉乘公式（严格匹配Julia）
		cross_vec(1) = v1(2)*v2(3) - v1(3)*v2(2)        ! y1z2 - z1y2
		cross_vec(2) = v1(3)*v2(1) - v1(1)*v2(3)        ! z1x2 - x1z2（修正符号）
		cross_vec(3) = v1(1)*v2(2) - v1(2)*v2(1)        ! x1y2 - y1x2
	end function cross

	! 向量的模（2-范数）
	subroutine mag(vec, sqrSum)
		real(kind=8), intent(in) :: vec(:)
		real(kind=8), intent(out) :: sqrSum  ! 输出参数：原返回值
		integer(kind=8) :: i
		sqrSum = 0.0_8
		do i = 1, size(vec)
		    sqrSum = sqrSum + vec(i)**2  ! 累加元素平方
		end do
		sqrSum = sqrt(sqrSum)  ! 开平方得模长
	end subroutine mag
	! 向量归一化（转化为单位向量）
	subroutine normalize(vec, unit_vec)
		real(kind=8), intent(in) :: vec(:)
		real(kind=8), intent(out) :: unit_vec(:)
		real(kind=8) :: vec_mag
		call mag(vec, vec_mag)  ! 调用mag函数
		! 避免除以零
		if (vec_mag < 1.0d-12) then
		    print *, "Error in normalize: Zero magnitude vector"
		    stop
		end if		
		unit_vec = vec / vec_mag  ! 每个元素除以模长
	end subroutine normalize
	subroutine dot_product_real(A, B, result, len_vec)
		implicit none
		! 输入：A/B为1维向量；len_vec为向量长度（可选，适配Tapenade）
		! 输出：点积结果（标量）
		real*8, intent(in)  :: A(:), B(:)                ! 1维向量（现有调用的输入）
		real*8, intent(out) :: result                    ! 标量结果
		integer(8), intent(in), optional :: len_vec         ! 新增：向量长度（可选，Tapenade用）
		integer :: i, vec_len                            ! 局部：实际使用的向量长度

		! 1. 确定向量长度：优先用显式传递的len_vec，否则用size(A)（兼容现有调用）
		if (present(len_vec)) then
		    vec_len = len_vec  ! Tapenade微分时，显式传递长度（如3）
		else
		    vec_len = size(A)  ! 现有调用：默认用A的长度（不影响原有逻辑）
		end if

		! 2. 计算点积（直接用显式长度vec_len循环，避免Tapenade无法识别size(B)）
		result = 0.0_8
		do i = 1, vec_len
		    result = result + A(i) * B(i)
		end do

	end subroutine dot_product_real
	subroutine dot_product_reall(A, B, result, len_vec)
		implicit none
		! 输入：A/B为1维向量；len_vec为向量长度（可选，适配Tapenade）
		! 输出：点积结果（标量）
		real*8, intent(in)  :: A(:), B(:)                ! 1维向量（现有调用的输入）
		real*8, intent(out) :: result                    ! 标量结果
		integer(8), intent(in), optional :: len_vec         ! 新增：向量长度（可选，Tapenade用）
		integer :: i, vec_len                            ! 局部：实际使用的向量长度

		! 1. 确定向量长度：优先用显式传递的len_vec，否则用size(A)（兼容现有调用）
		if (present(len_vec)) then
		    vec_len = len_vec  ! Tapenade微分时，显式传递长度（如3）
		else
		    vec_len = size(A)  ! 现有调用：默认用A的长度（不影响原有逻辑）
		end if

		! 2. 计算点积（直接用显式长度vec_len循环，避免Tapenade无法识别size(B)）
		result = 0.0_8
		do i = 1, vec_len
		    result = result + A(i) * B(i)
		end do

	end subroutine dot_product_reall
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 11.NACA0012.jl!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! All file paths are relative to the repository's main directory, 'include' this script from there
	!-----------------------------------------------------------------------------------#
	! 变量赋值
	!-----------------------------------------------------------------------------------#
	subroutine compute_CFD(cellPrimitivess)
		type(Fluidd) :: fluid
!		type(Meshh), allocatable :: mesh
!		real*8, allocatable, intent(in) :: point_update(:,:) ! 输入：变形网格点
		real(kind=8) :: cellPrimitives(18513, 5)
!		type(SolutionState), allocatable, intent(inout) :: sln
		real(kind=8), allocatable, intent(inout) :: cellPrimitivess(:,:)
!		real(kind=8), intent(inout) :: cellPrimitivess_1
!		real(kind=8), allocatable :: UunitVec(:)  
		integer(kind=8) :: nCells, meshInfo(4), f, nFaces, nBoundaries, max_faces_per_cell
		real(kind=8) :: P, T, U(3), a, machNum, Pt, Tt, mag_U
		character(256) :: meshPath
		type(BoundaryCondition) :: boundaryConditions(4)
		type(MeshData) :: tempMesh
		integer(kind=8), parameter :: max_bc_params = 5
		! 核心新增：定义临时数组存储边界条件参数（避免嵌套allocate）
   		real(kind=8), allocatable :: bc3_params(:), bc4_params(:)
        ! 初始化流体属性
		fluid%Cp = 1005.0d0
		fluid%R = 287.05d0
		fluid%gammaa = 1.4d0					
		! 初始化边界条件
		P = 100000.0d0
		T = 300.0d0
		!U = [ 277.7091d0, 6.059633d0, 0d0 ]
		!U = [33.30d0, 0.7d0, 0.0d0] !mach=0.096
		U = [33.30d0, 0.0d0, 0.0d0] !mach=0.096		
		call normalize(U, UunitVec)
		a = sqrt(fluid%gammaa * fluid%R * T)
		call mag(U, mag_U)  ! 子程序调用：计算 U 的模长
		machNum = mag_U / a  ! 用临时变量代入表达式
		print *, "machNum = ", machNum
		Pt = P*(1.0_8 + ((fluid%gammaa-1.0_8)/2.0_8)*machNum**2)**(fluid%gammaa/(fluid%gammaa-1.0_8))
		Tt = T*(1.0_8 + ((fluid%gammaa-1.0_8)/2.0_8)*machNum**2)
				
		! 修正边界条件参数大小（原代码中params(4)分配与5个元素赋值不匹配）
!		allocate(boundaryConditions(4))
		boundaryConditions(1)%type = wallBoundary
!		allocate(boundaryConditions(1)%params(0))
		boundaryConditions(2)%type = emptyBoundary
!		allocate(boundaryConditions(2)%params(0))
        boundaryConditions(1)%params(1:5) = 0.0_8
		boundaryConditions(2)%params(1:5) = 0.0_8
		boundaryConditions(3)%type = InletBoundary
!		allocate(boundaryConditions(3)%params(5))  ! 改为5个参数（Pt, Tt, Ux, Uy, Uz）
!        allocate(bc3_params(5))
!        bc3_params = [Pt, Tt, UunitVec(1), UunitVec(2), UunitVec(3)]
!		boundaryConditions(3)%params = [Pt, Tt, UunitVec(1), UunitVec(2), UunitVec(3)] 
!        allocate(boundaryConditions(3)%params(size(bc3_params)))
!        boundaryConditions(3)%params = bc3_params
!        deallocate(bc3_params)
        boundaryConditions(3)%params(1:5) = [Pt, Tt, UunitVec(1), UunitVec(2), UunitVec(3)]
        
		boundaryConditions(4)%type = OutletBoundary
!		allocate(boundaryConditions(4)%params(1))
!		boundaryConditions(4)%params(1) = P		
!		allocate(bc4_params(1))
!	    bc4_params(1) = P
!  		allocate(boundaryConditions(4)%params(size(bc4_params)))
!        boundaryConditions(4)%params = bc4_params
!        deallocate(bc4_params) ! 及时释放临时数组    
            
		boundaryConditions(4)%params(1) = P
		boundaryConditions(4)%params(2:5) = 0.0_8
		
		! 读取网格
		meshPath = "mesh/OFairfoilMesh"
!!!!!!!!!!!	   
!		allocate(tempMesh%owner(0), neighbour_tempMesh(0))		    
        call readOpenFOAMMesh(meshPath, point_update, tempMesh)    
        nFaces = size(owner_tempMesh)
        nCells = maxval(owner_tempMesh)
		nBoundaries = size(tempMesh%boundaryNames)    

!        allocate(first5_refs(5,3))		  
		!call OpenFOAMMesh(meshPath, point_update, mesh)
		print *, "abc"	
		call OpenFOAMMesh(meshPath, point_update, tempMesh)
        print *, "abc"	              
!        deallocate(first5_refs)
        call deallocateMeshData(tempMesh)
!!!!!!!!!!!!!!!!!!!!!!		                 
		write(*,*) "test_First 5 mesh.cVols:"
		do f = 1, 5
		    write(*,'(A,I3,A,E12.6)') "Cell ", f, ": ", cVols_mesh(f)
		end do
		! 修改后的代码
		call unstructuredMeshInfo(meshInfo) 
		nCells = meshInfo(1)  
!		allocate(cellPrimitives(nCells, 5)) 
	    ! 初始化均匀解
	    print *, "a1"
		call initializeUniformSolution3D(P, T, U(1), U(2), U(3))
		cellPrimitives = sol
		print *, "a1"
		! 调用求解器（局部时间步长）
		call solve(meshPath, cellPrimitives, boundaryConditions, point_update, cellPrimitivess)
!		cellPrimitivess_1 = cellPrimitivess(1, 1)
!		if (allocated(boundaryConditions)) deallocate(boundaryConditions)
	end subroutine compute_CFD	
	subroutine compute_steady_residual(data_4D137, cellPrimitives, fluxResidualsOut)
		implicit none
		real(kind=8), intent(in) :: data_4D137(1, 4)
		real(kind=8), intent(in) :: cellPrimitives(:,:)
		real(kind=8), intent(inout) :: fluxResidualsOut(:,:)
		type(Fluidd) :: fluid
		type(BoundaryCondition) :: boundaryConditions(4)
		type(MeshData) :: tempMesh
		integer(kind=8) :: nCells, meshInfo(4)
		real(kind=8) :: P, T, U(3), a, machNum, Pt, Tt, mag_U
		character(256) :: meshPath
		call airfoil_deformation_HH(data_4D137)
		fluid%Cp = 1005.0d0
		fluid%R = 287.05d0
		fluid%gammaa = 1.4d0
		P = 100000.0d0
		T = 300.0d0
		U = [33.30d0, 0.0d0, 0.0d0]
		call normalize(U, UunitVec)
		a = sqrt(fluid%gammaa * fluid%R * T)
		call mag(U, mag_U)
		machNum = mag_U / a
		Pt = P*(1.0_8 + ((fluid%gammaa-1.0_8)/2.0_8)*machNum**2)**(fluid%gammaa/(fluid%gammaa-1.0_8))
		Tt = T*(1.0_8 + ((fluid%gammaa-1.0_8)/2.0_8)*machNum**2)
		boundaryConditions(1)%type = wallBoundary
		boundaryConditions(2)%type = emptyBoundary
		boundaryConditions(1)%params(1:5) = 0.0_8
		boundaryConditions(2)%params(1:5) = 0.0_8
		boundaryConditions(3)%type = InletBoundary
		boundaryConditions(3)%params(1:5) = [Pt, Tt, UunitVec(1), UunitVec(2), UunitVec(3)]
		boundaryConditions(4)%type = OutletBoundary
		boundaryConditions(4)%params(1) = P
		boundaryConditions(4)%params(2:5) = 0.0_8
		meshPath = "mesh/OFairfoilMesh"
		call readOpenFOAMMesh(meshPath, point_update, tempMesh)
		call OpenFOAMMesh(meshPath, point_update, tempMesh)
		call deallocateMeshData(tempMesh)
		call unstructuredMeshInfo(meshInfo)
		nCells = meshInfo(1)
		cellPrimitives_sln = cellPrimitives
		call encodePrimitives3D(cellPrimitives, fluid, cellState_sln)
		cellFluxes_sln = 0.0_8
		fluxResiduals_sln = 0.0_8
		faceFluxes_sln = 0.0_8
		call decodeSolution_3D(fluid)
!		allocate(fluxResidualss(nCells, 5), source=0.0_8)
		call unstructured_JSTFlux(boundaryConditions, fluid, fluxResidualsOut)
	end subroutine compute_steady_residual
end module BuFlowModule
