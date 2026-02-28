module main_module
    use BuFlowModule  
    use meshdeformationn  
    
    implicit none 
!!!!!!!!!!      
!!!!!!!!!!	 

contains  ! 模块内包含子程序，自动生成显式接口

    subroutine main(input_param, cellPrimitivesout, cellPrimitivesout_1)   
        use TypesModule
    	!$AD &input_param          ! 标记输入
        !$AD &cellPrimitives       ! 标记输出
        ! 输入参数：仅原data_4D137的第一个维度（标量）
		real(kind=8), intent(in) :: input_param  
		real(kind=8), allocatable, intent(out) :: cellPrimitivesout(:,:)  	
		real(kind=8), intent(out) :: cellPrimitivesout_1	      
        ! 局部变量：构造1×4变形参数
        real(kind=8) :: data_4D137(1,4)  
        ! 传递变形网格点的变量
        integer(kind=8) :: nCells, meshInfo(4), iostat, nLines_points, file_unit_points, nFaces, nFacePerBdry, nBoundaries&
        , nLines_boundary, file_unit_boundary, nLines_owner, file_unit_owner, nLines_neighbour, file_unit_neighbour&
        , nLines_Faces, file_unit_Faces, startLine_Faces, fCount, i, bracketL, bracketR, pt_count
        integer(kind=8) :: bCount, startLine_boundary, nPoints 
!        real(kind=8), allocatable :: point_update(:,:)  
        character(256) :: meshPath, pointsFilePath, boundaryFilePath, ownerFilePath, neighbourFilePath, facesFilePath, line
        type(MeshData) :: tempMesh

!!!!!!!!!!!!
       
!!!!!!!!!!!!
!!!!!!!!!!!!!!meshdeformation
		file_unit_points = get_free_unitt()
		meshPath = "mesh/OFairfoilMesh"
		pointsFilePath = trim(meshPath)//"/points"
        open(unit=file_unit_points, file=pointsFilePath, status='old', action='read', iostat=iostat)        
        nLines_points = 0
        do
            read(file_unit_points, '(a)', iostat=iostat)
            if (iostat /= 0) exit
            nLines_points = nLines_points + 1
        end do
        rewind(file_unit_points)
		allocate(lines_defor(nLines_points))
		
		boundaryFilePath = trim(meshPath) // "/boundary"
		file_unit_boundary = get_free_unit()
        open(unit=file_unit_boundary, file=boundaryFilePath, status='old', action='read', iostat=iostat)     
        nLines_boundary = 0
        do
            read(file_unit_boundary, '(a)', iostat=iostat)
            if (iostat /= 0) exit
            nLines_boundary = nLines_boundary + 1
        end do
        rewind(file_unit_boundary)       
        allocate(bLines(nLines_boundary))
        do i = 1, nLines_boundary
            read(file_unit_boundary, '(a)') bLines(i)
        end do
        close(file_unit_boundary)       
        call OFFile_FindNItemsss(bLines, startLine_boundary, bCount)
        allocate(boundaryNumFacess(bCount), boundaryStartFacess(bCount))
        
        ownerFilePath = trim(meshPath) // "/owner"
        file_unit_owner = get_free_unit()
        open(unit=file_unit_owner, file=ownerFilePath, status='old', action='read', iostat=iostat)
        if (iostat /= 0) then
            print *, "Error: Can't open owner file: ", trim(meshPath)
            return
        end if        
        nLines_owner = 0
        do
            read(file_unit_owner, '(a)', iostat=iostat)
            if (iostat /= 0) exit
            nLines_owner = nLines_owner + 1
        end do
        rewind(file_unit_owner)
        allocate(lines_owner(nLines_owner))
        
        ! 2. 打开文件（与owner逻辑完全一致）
        neighbourFilePath = trim(meshPath) // "/neighbour"
		file_unit_neighbour = get_free_unit()
		open(unit=file_unit_neighbour, file=neighbourFilePath, status='old', action='read', iostat=iostat)
		if (iostat /= 0) then
		    print *, "Error: Can't open neighbour file: ", trim(meshPath)  ! 错误信息与owner格式统
		    return
		end if
		nLines_neighbour = 0
		do
		    read(file_unit_neighbour, '(a)', iostat=iostat)  ! 逐行读取计数
		    if (iostat /= 0) exit
		    nLines_neighbour = nLines_neighbour + 1
		end do
		rewind(file_unit_neighbour)  ! 重置文件指针
		allocate(lines_neighbour(nLines_neighbour))  ! 分配存储所有行的数组

		facesFilePath = trim(meshPath) // "/faces"
		file_unit_Faces = get_free_unit()
        open(unit=file_unit_Faces, file=facesFilePath , status='old', action='read', iostat=iostat)  
        if (iostat /= 0) then
            print *, "Error: Can't open faces file: ", trim(meshPath)
!            allocate(faces%faces(0))
            return
        end if    
        nLines_Faces = 0
        do
            read(file_unit_Faces, '(a)', iostat=iostat)
            if (iostat /= 0) exit
            nLines_Faces = nLines_Faces + 1
        end do
        rewind(file_unit_Faces)
        allocate(lines_Faces(74151))
        do i = 1, nLines_Faces
            read(file_unit_Faces, '(a)') lines_Faces(i)
        end do
        close(file_unit_Faces)                
        call OFFile_FindNItemsss(lines_Faces, startLine_Faces, fCount)
        allocate(tmp_faces(fCount))
!       do i = 1, fCount
!      		allocate(tmp_faces(i)%points(4))
!       end do
		print *, "ABC"		
		do i = 1, 2
            line = trim(lines_Faces(startLine_Faces + i - 1))
            bracketL = index(line, '(')
            bracketR = index(line, ')')
            read(line(1:bracketL-1), *) pt_count            
        end do
        allocate(pts_Faces(pt_count), source=0_int64)
        
!		allocate(wing_update(99, 2), source=0.0_real64)
!!!!!!!!!!!!!!
!!!!!!!!!!!!!!Buflow
		allocate(owner_tempMesh(0), neighbour_tempMesh(0))	
		allocate(facePoints_tempMesh(0,0))
		PRINT *, "A"
		allocate(points_tempMesh(37243, 3))
        call readOpenFOAMMeshh(meshPath, tempMesh)   
        print *, "ABC"
        nPoints = size(points_tempMesh,1)
        nFaces = size(owner_tempMesh)
        nCells = maxval(owner_tempMesh)
		nBoundaries = size(tempMesh%boundaryNames) !四种边界
		nFacePerBdry = size(tempMesh%boundaryNumFaces)  !所有边界面的总数
    	allocate(currentBoundary(nFacePerBdry)) 
    	allocate(gradP(nCells, 3))
		allocate(sj(nCells), sjCount(nCells))
		allocate(rj(nCells)) 
		allocate(rjsjF(nFaces, 2), source=0.0_8) 
		allocate(eps2(nFaces), source=0.0_8)
		allocate(eps4(nFaces), source=0.0_8)
		allocate(P_matrix(nCells, 1))
		allocate(faceVals(nFaces, 5))
!		allocate(P_eps(nCells))
		allocate(temp_grad(nCells, 1, 3))!这个是正确写法，错误写法在程序运行完之后会报错。正确写法和错误写法结果一样
		allocate(timeFluxes(nCells))
		allocate(surfaceAreas(nCells))
		allocate(timeFluxess(nCells))
		allocate(surfaceAreass(nCells))
		allocate(faceRhoT(nFaces, 2), source=0.0_8)		
		allocate(cellRhoT(nCells, 2), source=0.0_8)
		allocate(faceVel(3), positionn(3))
		allocate(fDeltas(nFaces, 5), source=0.0_8)
		allocate(fDGrads(nCells, 5, 3))
		allocate(eps22(nFaces), eps44(nFaces))
		allocate(epss(nFaces, 2), source=0.0_8)
		allocate(diffusionFlux(5), unitFA(3), fD(5), farOwnerfD(5), farNeighbourfD(5))
        allocate(faceVals_dd(nFaces, 5)) 
        allocate(faceVals_d(nFaces, 5)) 
        allocate(fluxResiduals(nCells, 5), source=0.0_8) 
		allocate(dt_val(nCells), source=0.0_8) 
!		allocate(first5_refs(5,3))	
		allocate(encodePrimitives3DD(nCells, 5))
		allocate(initialValues(5))
		! 分配结构体数组
		allocate(cellPrimitives_sln(nCells, 5), source=0.0d0)
		allocate(cellState_sln(nCells, 5), source=0.0d0)
		allocate(cellFluxes_sln(nCells, 15), source=0.0d0)
		allocate(fluxResiduals_sln(nCells, 5), source=0.0d0)
		allocate(faceFluxes_sln(nFaces, 15), source=0.0d0)
        allocate(faces_mesh(nFaces, 2))
        allocate(fAVecs_mesh(nFaces, 3), fCenters_mesh(nFaces, 3))                
		allocate(cells_mesh(nCells, 6))
		allocate(cVols_mesh(nCells), cCenters_mesh(nCells, 3), cellSizes_mesh(nCells, 3))
		allocate(boundaryFaces_mesh(nBoundaries, maxval(tempMesh%boundaryNumFaces)))
		allocate(dt_solve(nCells))
		allocate(sol(nCells, 5))
        allocate(fCs(6,3), cell_fAVecs(6,3))
        allocate(cellPts(0,3))
!        allocate(cellPtss(1,3)) 
        allocate(facePts(size(facePoints_tempMesh(1, :)), 3))
        allocate(cellFaceCount(nCells))    
        allocate(UunitVec(3)) 
!        allocate(points_tempMesh(size(point_update,1), size(point_update,2))) !在readOFPointsFilef里面已经allocate
		allocate(points_meshdefor(37254,3), source=0.0_8) 
		allocate(inoutput(99, 2))
		allocate(wing(99, 2))
		allocate(wing_coords(99,2), inoutput_coords(99,2))
		allocate(control_points(198,2))
		allocate(r_row(198), phi_row(198))
		allocate(mat_inv(198, 198), source=0.0_8)
		allocate(distances(198))
		allocate(phi(198, 198))
		allocate(point_update(nPoints,3))
		
!!!!!!!!!!!!!!Buflow
 !       !$AD intent(out) :: cellPrimitives(nCells, 5)  ! 强制 Tapenade 识别维度（nCells 已定义）[-0.75d0, 0.485d0, -0.588d0, -0.912d0]
        ! 构造变形参数（仅第一个维度为输入，其余为0）
        data_4D137(1,1) = input_param  
!        data_4D137(1,2:4) = 0.0d0  ! 后三个维度固定为0
        data_4d137(1, 2) = 0.485d0
		data_4d137(1, 3) = -0.588d0
		data_4d137(1, 4) = -0.912d0  ! 后三个维度固定为0
        
        ! 读取网格
!		meshPath = "mesh/OFairfoilMesh"	
!		mesh = OpenFOAMMesh(meshPath, point_update)
		! 修改后的代码
!		meshInfo = unstructuredMeshInfo(mesh) 
!		nCells = meshInfo(1)  
!		allocate(cellPrimitives(nCells, 5)) 
        allocate(cellPrimitivesout(18513, 5))
        cellPrimitivesout = 0.0d0  ! 显式赋值，确保 Tapenade 识别为“活跃输出”
		print *, "ABC"
        ! 调用网格变形，获取变形后的网格点（point_update）
        call airfoil_deformation_HH(data_4D137)
        print *, "ABC"
        ! 调用CFD计算，传入变形网格点，输出结果到cellPrimitives
        call compute_CFD(cellPrimitivesout, cellPrimitivesout_1)  
!        call unstructured_JSTFlux_AD(boundaryConditions, fluid, unstructured_JSTFluxx)
        ! 释放point_update内存        
!!!!!!!!!!!!meshdeformation  
        deallocate(lines_defor)
        if (allocated(inoutput)) deallocate(inoutput)
		if (allocated(wing)) deallocate(wing)	
		if (allocated(point_update)) deallocate(point_update) 
        if (allocated(points_meshdefor)) deallocate(points_meshdefor)	
        if (allocated(mat_inv)) deallocate(mat_inv)
        if (allocated(wing_coords)) deallocate(wing_coords)
        if (allocated(inoutput_coords)) deallocate(inoutput_coords)
        if (allocated(r_row)) deallocate(r_row)
        if (allocated(phi_row)) deallocate(phi_row)
        if (allocated(control_points)) deallocate(control_points)
        if (allocated(phi)) deallocate(phi)
        if (allocated(distances)) deallocate(distances)
!		if (allocated(wing_update)) deallocate(wing_update)
!!!!!!!!!!!! meshdeformation
!!!!!!!!!!!! Buflow
        deallocate(currentBoundary)
!        deallocate(P_eps)
		deallocate(P_matrix, temp_grad, gradP)
		deallocate(sj, sjCount, rj, rjsjF, eps2, eps4)	
		deallocate(faceVals)
		deallocate(timeFluxes)
		deallocate(surfaceAreas)
		deallocate(timeFluxess)
		deallocate(surfaceAreass)
		deallocate(faceVel, positionn)
		deallocate(faceRhoT)
		deallocate(cellRhoT)
		deallocate(epss) 
		deallocate(fDeltas, fDGrads, eps22, eps44, diffusionFlux, unitFA, fD, farOwnerfD, farNeighbourfD)
		deallocate(faceVals_dd)
		deallocate(faceVals_d)
		deallocate(fluxResiduals, dt_val) 
		print *, "Aa"
		deallocate(bLines)
		print *, "c"
!		if (allocated(lines_neighbour)) deallocate(lines_neighbour)
!     	deallocate()
		print *, "c"
		deallocate(lines_owner)
		print *, "Aa"
		deallocate(lines_Faces)
		deallocate(pts_Faces)
		deallocate(tmp_faces)
		deallocate(dt_solve)
!		deallocate(first5_refs(5,3))
		deallocate(fCs, cell_fAVecs, cellPts, cellFaceCount, facePts)
		deallocate(UunitVec)
		
		print *, "Aa"
!!!!!!!!!!!! Buflow		
    end subroutine main
    subroutine readOpenFOAMMeshh(polyMeshPath, tempMesh)
        character(len=*), intent(in) :: polyMeshPath
        type(MeshData), intent(out) :: tempMesh
        character(len=256) :: pointsFilePath, facesFilePath, ownerFilePath, neighbourFilePath, boundaryFilePath
        
        pointsFilePath = trim(polyMeshPath) // "/points"
        facesFilePath = trim(polyMeshPath) // "/faces"
        ownerFilePath = trim(polyMeshPath) // "/owner"
        neighbourFilePath = trim(polyMeshPath) // "/neighbour"
        boundaryFilePath = trim(polyMeshPath) // "/boundary"
        
        call readOFPointsFile(pointsFilePath)
        facePoints_tempMesh = readOFFacesFile(facesFilePath)
        owner_tempMesh = readOFOwnerFile(ownerFilePath)
        neighbour_tempMesh = readOFNeighbourFile(neighbourFilePath)
        call readOFBoundaryFile(boundaryFilePath, tempMesh%boundaryNames, tempMesh%boundaryNumFaces, tempMesh%boundaryStartFaces)
        
    end subroutine  readOpenFOAMMeshh
    subroutine OFFile_FindNItemsss(fileLines, startLine, itemCount)
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
            if (isNumberr(trim(fileLines(i)))) then
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
    end subroutine OFFile_FindNItemsss
    logical function isNumberr(str)
        character(len=*), intent(in) :: str
        real(kind=8) :: num
        integer(kind=8) :: iostat
        isNumberr = .false.
        read(str, *, iostat=iostat) num
        if (iostat == 0) isNumberr = .true.
    end function isNumberr
    

end module main_module
