! 主程序：通过use模块获取main的显式接口
program run_main
    use main_module  ! 关键：使用包含main的模块，自动获取显式接口
    implicit none
    real(kind=8) :: data_4D137(1,4)  ! 输入参数
!    !$AD INPUT  ! 标记为微分的输入变量（自变量）
    real(kind=8), allocatable :: cellPrimitivesout(:,:) 
!    real(kind=8), allocatable :: fluxResidualsOut(:,:)

    data_4D137(1,1) = -0.75d0  
!   data_4D137(1,2:4) = 0.0d0  ! 后三个维度固定为0
    data_4d137(1, 2) = 0.485d0
    data_4d137(1, 3) = -0.588d0
    data_4d137(1, 4) = -0.912d0  ! 后三个维度固定为0 

    ! 调用main子程序（此时编译器已知晓其接口，支持可分配数组参数）
    call main(data_4D137, cellPrimitivesout)


    if (allocated(cellPrimitivesout)) then
        print *, "CFD计算完成，结果维度：", size(cellPrimitivesout, 1), "×", size(cellPrimitivesout, 2)
        deallocate(cellPrimitivesout)  ! 释放内存
    end if
	if (allocated(fluxResiduals)) deallocate(fluxResiduals)
end program run_main
