! 主程序：通过use模块获取main的显式接口
program run_main
    use main_module  ! 关键：使用包含main的模块，自动获取显式接口
    implicit none
    real(kind=8) :: input_param  ! 输入参数
!    !$AD INPUT  ! 标记为微分的输入变量（自变量）
    real(kind=8), allocatable :: cellPrimitivesout(:,:) 
    real(kind=8) :: cellPrimitivesout_1
!	!$AD OUTPUT  ! 标记为微分的输出变量（因变量） 
    ! 将AD指令放在变量声明之后，程序逻辑之前
    !$AD &input_param
    !$AD &cellPrimitivesout
    ! 设置输入参数（示例值）
    input_param = -0.75d0  
!    input_param = 0.0d0  
    ! 调用main子程序（此时编译器已知晓其接口，支持可分配数组参数）
    call main(input_param, cellPrimitivesout, cellPrimitivesout_1)


    if (allocated(cellPrimitivesout)) then
        print *, "CFD计算完成，结果维度：", size(cellPrimitivesout, 1), "×", size(cellPrimitivesout, 2)
        deallocate(cellPrimitivesout)  ! 释放内存
    end if

end program run_main
