    module variables
    implicit none
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!输入变量声明!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!输入变量声明!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!输入变量声明!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer::time                       !!!!时间参数，单位day
    integer::n_axis                     !!!!轴向划分网格数
    integer::n_radial                   !!!!径向划分网格数
    integer::coolant_kind               !!!!指定冷却剂种类，1--水；2--铅铋
    integer::module_identification      !!!!模式识别，1--稳态；2--瞬态
    integer::transient_mode             !!!!瞬时输入参数方式识别，1--功率；2--包壳温度
    real(8)::pi                         !!!!圆周率，常量
    real(8)::p_line_factor              !!!!功率因子
    real(8)::p_line_average             !!!!平均线功率密度，单位W/m
    real(8)::coolant_V_flow             !!!!冷却剂流量，单位m3/s
    real(8)::coolant_T_in               !!!!冷却剂进口温度，单位K
    real(8)::cladding_length            !!!!初始燃料元件高度，单位m
    real(8)::cladding_width             !!!!初始包壳厚度，单位m
    real(8)::gas_gap                    !!!!初始气隙间距，单位m
    real(8)::pitch                      !!!!栅距，单位m
    real(8)::pellet_diameter            !!!!芯块初始直径，单位m
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!输入变量声明!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!输入变量声明!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!输入变量声明!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!中间变量声明!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!中间变量声明!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!中间变量声明!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    parameter(pi=3.1415926)
    real,allocatable::p_line(:)                     !!!!线功率，单位W/m
    real,allocatable::Q_line(:)                     !!!!轴向各点热量，单位W
    real,allocatable::Temperature(:,:,:)            !!!!任意时刻，任意位置处温度，n_radial+1代表包壳内表面温度；n_radial+2代表包壳外表面温度；n_radial+3代表冷却剂温度，单位K
    real,allocatable::d(:,:,:)                      !!!!任意时刻，任意位置处直径，n_radial+1代表包壳内径；n_radial+2代表包壳外径，单位m
    real,allocatable::Bu(:,:,:)                     !!!!任意时刻，任意位置处，燃耗，单位MWD/kg
    real,allocatable::Bu_old(:,:,:)                 !!!!上一时刻，任意位置处，燃耗，单位MWD/kg
    real,allocatable::k_fuel(:,:,:)                 !!!!任意时刻，任意位置处，芯块热导率，单位W/m*K
    real,allocatable::P_SQUARE(:,:,:)               !!!!任意时刻，任意位置处，体积功率，单位W/m3
    real(8)::p_line_max                             !!!!最大线功率，单位W/m
    real(8)::x                                      !!!!功率函数截尾位置，单位m
    real(8)::sin_average                            !!!!余弦平均值
    real(8)::Fq                                     !!!!热工因子（迭代计算）
    real(8)::coolant_density                        !!!!冷却剂密度，单位kg/m3
    real(8)::coolant_T                              !!!!冷却剂温度，单位K
    real(8)::coolant_Cp                             !!!!冷却剂定压比热容，单位J/(kg*k)
    real(8)::fuel_percentage                        !!!!燃料百分比，单位%
    real(8)::De                                     !!!!当量直径，单位m
    real(8)::coolant_M_flow                         !!!!冷却剂质量流量，单位kg/s
    real(8)::coolant_S                              !!!!冷却剂流通面积，单位m2
    real(8)::coolant_speed                          !!!!冷却剂流速，单位m/s
    real(8)::Pr,Re,Pe,Nu                            !!!!普朗特数，雷诺数等无量纲数
    real(8)::cladding_T_surface                     !!!!包壳外表面温度，单位K
    real(8)::cladding_T,Tc,Tc1                      !!!!包壳迭代计算中间值，单位K
    real(8)::cladding_T_internal                    !!!!包壳内表面温度，单位K
    real(8)::pellet_T_surface                       !!!!芯块表面温度，单位K
    real(8)::Bu_average                             !!!!平均燃耗，单位MWD/kg
    real(8)::UO2_idealdensity                       !!!!二氧化铀理论密度，单位kg/m3
    real(8)::fule_percentage                        !!!!燃料理论密度，单位%
    real(8)::UO2_density                            !!!!二氧化铀密度，单位kg/m3
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!中间变量声明!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!中间变量声明!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!中间变量声明!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!瞬态变量声明!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!瞬态变量声明!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!瞬态变量声明!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer t_number                                !!!!瞬态计算点个数
    real(8)::AT                                     !!!!瞬态计算步长，单位s
    real(8)::UO2_Cp                                 !!!!二氧化铀定压比热容，单位J/(kg*K)
    real(8)::heat_dv,heat_old                       !!!!储热，储热循环变量
    real(8)::length_every                           !!!!轴向单元平均长度，单位m
    allocatable::Temperature_transient(:,:,:)       !!!!瞬态温度，单位K
    allocatable::T_transient(:)                     !!!!瞬态计算对应时刻，单位s
    allocatable::p_average_transient(:)             !!!!瞬态计算各时刻对应平均功率，单位W/m
    allocatable::coolant_T_transient(:)             !!!!瞬态计算各点冷却剂温度，单位K




    !定义数组为双精度变量
    double precision Temperature_transient
    double precision T_transient
    double precision p_average_transient
    double precision coolant_T_transient


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!瞬态变量声明!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!瞬态变量声明!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!瞬态变量声明!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    contains

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!打开文件夹!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!打开文件夹!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!打开文件夹!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine open_file
    if (module_identification==1) then
        open(1001,file='wentaishuru\canshushuru.txt')
        open(1002,file='wentaishuchu\canshushuchu.txt')
    else
        if (transient_mode==1)then
            open(1003,file='shuntaishuru\shuntaicanshushuru001.txt')
            open(1004,file='shuntaishuchu\canshushuchu001.txt')
        else
            open(1005,file='shuntaishuru\shuntaicanshushuru002.txt')
            open(1006,file='shuntaishuchu\canshushuchu002.txt')
        end if
    end if
    end subroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!打开文件夹!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!打开文件夹!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!打开文件夹!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!输入参数!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!输入参数!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!输入参数!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine read_file
    integer i,j
    if (module_identification==1) then
        read (1001,*)
        read (1001,*) time,n_axis,n_radial,coolant_kind,p_line_factor,p_line_average
        read (1001,*)
        read (1001,*) coolant_V_flow,coolant_T_in,cladding_length,cladding_width,gas_gap,pitch
        read (1001,*)
        read (1001,*) pellet_diameter
        !---------------------------------------------------------------------------------参数初始化
        call variable_return
    else
        if (transient_mode==1) then
            read (1003,*)
            read (1003,*) time,n_axis,n_radial,coolant_kind,p_line_factor,p_line_average
            read (1003,*)
            read (1003,*) coolant_V_flow,coolant_T_in,cladding_length,cladding_width,gas_gap,pitch
            read (1003,*)
            read (1003,*) pellet_diameter
            read (1003,*) t_number
            !-----------------------------------------------------------------------------参数初始化
            call variable_return
            do i=1,n_axis
                read (1003,*) (Temperature(1,i,j),j=1,n_radial+3)
                Temperature(1,:,:)=500
            end do
            do i=1,t_number
                read (1003,*) T_transient(i),p_average_transient(i)
            end do
        else

            read (1005,*)
            read (1005,*) time,n_axis,n_radial,coolant_kind,p_line_factor,p_line_average
            read (1005,*)
            read (1005,*) coolant_V_flow,coolant_T_in,cladding_length,cladding_width,gas_gap,pitch
            read (1005,*)
            read (1005,*) pellet_diameter,AT
            !-----------------------------------------------------------------------------参数初始化
            call variable_return
            do i=1,n_axis
                read (1005,*) Temperature(1,i,n_radial+3)
            end do
        end if
    end if
    end subroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!输入参数!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!输入参数!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!输入参数!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!参数初始化!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!参数初始化!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!参数初始化!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine variable_return
    integer i
    allocate(p_line(n_axis))
    allocate(Q_line(n_axis))
    allocate(Temperature(time,n_axis,n_radial+3))
    allocate(Bu(time,n_axis,n_radial))
    allocate(Bu_old(time,n_axis,n_radial))
    if (module_identification==2.AND.transient_mode==1) then
        allocate(Temperature_transient(t_number,n_axis,n_radial+3))
        allocate(T_transient(t_number))
        allocate(p_average_transient(t_number))
        allocate(d(t_number,n_axis,n_radial+2))
        allocate(coolant_T_transient(n_axis))
        allocate(P_SQUARE(t_number,n_axis,n_radial-1))
        allocate(k_fuel(t_number,n_axis,n_radial))
    else
        allocate(d(time,n_axis,n_radial+2))
        allocate(P_SQUARE(time,n_axis,n_radial-1))
        allocate(k_fuel(time,n_axis,n_radial))
    end if
    d(:,:,n_radial+2)=pellet_diameter+2*gas_gap+2*cladding_width    !!!!初始时刻，包壳外径，单位m
    d(:,:,n_radial+1)=pellet_diameter+2*gas_gap                     !!!!初始时刻，包壳内径，单位m
    coolant_S=pitch**2.-pi/4.*(pellet_diameter+2*gas_gap+&          !!!!流通面积，单位m2
    &2*cladding_width)**2.
    De=4.*(pitch**2.-pi/4.*(pellet_diameter+2*gas_gap+2*&           !!!!当量直径，单位m
    &cladding_width)**2.)/pi/(pellet_diameter+2*gas_gap&
        &+2*cladding_width)
    length_every=cladding_length/((n_axis-1)*1.0)
    Bu_old(:,:,:)=0
    do i=1,n_radial
        d(:,:,i)=(i-1)*pellet_diameter/(n_radial-1)
    enddo
    end subroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!参数初始化!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!参数初始化!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!参数初始化!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



    end module variables
















