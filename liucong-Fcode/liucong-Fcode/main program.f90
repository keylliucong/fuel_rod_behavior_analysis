    program main
    use variables
    use physical_properties
    use thermal_hydraulic
    use hcore
    use match
    use thermal_hydraulic_transient
    implicit none
    integer i,j,k,l

    !print*,"请输入计算类型，1--稳态；2--瞬态"
    !read *,model_identification
    !if (model_identification==2) then
    !    print*,"请选择已知参数类型，1--功率；2--包壳温度"
    !    read *,transient_mode
    !end if

    call model_select
    call open_file                          !!!!!!!!!!打开文件夹
    call read_file                          !!!!!!!!!!参数输入


    call fuel_physical_properties




    if (model_identification==1) then
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!功率计算!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!功率计算!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!功率计算!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        call power_calculation
        write(1002,*)"功率输出"
        write(1002,"(f9.3)")(p_line(i),i=1,n_axis)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!功率计算结束!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!功率计算结束!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!功率计算结束!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!冷却剂温度计算!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!冷却剂温度计算!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!冷却剂温度计算!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        write(1002,*)"冷却剂温度输出"
        do i=1,n_axis
            call coolant_physical_properties(i)
            call coolant_T_calculation(i)           !!!!!!冷却剂温度计算
            Temperature(:,i,n_radial+3)=coolant_T   !!!!!!冷却剂温度更新
            write(1002,"(f9.3)")coolant_T
        end do
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!冷却剂温度计算结束!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!冷却剂温度计算结束!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!冷却剂温度计算结束!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!轴向温度计算!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!轴向温度计算!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!轴向温度计算!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        do i=1,time                                                 !!!!!!!!!!!!!!!时间循环
            do j=1,n_axis                                           !!!!!!!!!!!!!!!轴向循环
                call axis_T_calculation(i,j)                        !!!!!!!!!!!!!!!轴向温度计算
                call axis_T_update(i,j)                             !!!!!!!!!!!!!!!温度更新
                call Bu_calculation(i,j)                            !!!!!!!!!!!!!!!燃耗计算
                call P_S_calculation(i,j)
                call fuel_t_calculation(i,j)                        !!!!!!!!!!!!!!!芯块温度计算
            end do
        enddo
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!轴向温度计算结束!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!轴向温度计算结束!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!轴向温度计算结束!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!数据输出!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!数据输出!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!数据输出!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        write(1002,*)"轴向各点温度"
        do i=1,time
            do j=1,n_axis
                write(1002,55)(Temperature(i,j,k),k=1,n_radial+3)
                write(1002,56)(d(i,j,k)/2,k=1,n_radial+2)
55              format(14f9.3)
56              format(13f9.6)
            end do
        end do
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!数据输出!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!数据输出!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!数据输出!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    else
        do i=1,t_number
            if (i==1) then
                AT=T_transient(1)
            else
                AT=T_transient(i)-T_transient(i-1)
            end if
            if (transient_mode==1) then
                p_line_average=p_average_transient(i)
                call power_calculation

                !!!!!!!!!!!计算初始条件时燃耗begin （后续修改直接调用稳态计算值
                do k=1,time
                    do l=1,n_axis
                        call Bu_calculation(k,l)
                    end do
                enddo
                !!!!!!!!!!!计算初始条件时燃耗begin （后续修改直接调用稳态计算值

                !!!!!!!!!瞬态温度计算部分begin

                do j=1,n_axis
                    call P_S_calculation(i,j)
                    call coolant_physical_properties(j)
                    heat_dv=100.
                    heat_old=0.
                    ! do while(abs(heat_dv)>0.1)          !迭代条件，储热收敛
                    call coolant_T_calculation(j)   !冷却剂温度计算
                    Temperature_transient(i,j,n_radial+3)=coolant_T_transient(j)        !赋值，存储冷却剂温度

                    call T_transient_calculation(i,j)   !温度计算，TDMA
                    !end do

                end do

                !!!!!!!!!瞬态温度计算部分end



                write(1004,60)Temperature_transient(i,11,1),Temperature_transient(i,11,n_radial)        !输出数据
60              format(2f9.3)


            else
                p_line_average=p_average_transient(i)
                call power_calculation
                !!!!!!!!!!!计算初始条件时燃耗begin （后续修改直接调用稳态计算值
                do k=1,time
                    do l=1,n_axis
                        call Bu_calculation(k,l)
                    end do
                enddo
                !!!!!!!!!!!计算初始条件时燃耗begin （后续修改直接调用稳态计算值

                !!!!!!!!!瞬态温度计算部分begin

                do j=1,n_axis
                    call P_S_calculation(i,j)
                    call coolant_physical_properties(j)
                    !heat_dv=100.
                    !heat_old=0.
                    !do while(abs(heat_dv)>0.1)          !迭代条件，储热收敛
                    call coolant_T_calculation(j)   !冷却剂温度计算
                    Temperature_transient(i,j,n_radial+3)=coolant_T_transient(j)        !赋值，存储冷却剂温度
                    call T_transient_calculation(i,j)   !温度计算，TDMA
                    !end do
                end do

                !!!!!!!!!瞬态温度计算部分end

                write(1006,88)Temperature_transient(i,11,1),Temperature_transient(i,11,n_radial)
88              format(2f9.3)

                write(1004,*)""
            end if
        end do


        do i=1,t_number
            write(1004,99) T_transient(i),Temperature_transient(i,11,1),Temperature_transient(i,11,n_radial)
99          format(3f9.3)

        end do


    end if


    end