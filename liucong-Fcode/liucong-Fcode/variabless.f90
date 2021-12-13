    module variables
    implicit none
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!�����������!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!�����������!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!�����������!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer::time                       !!!!ʱ���������λday
    integer::n_axis                     !!!!���򻮷�������
    integer::n_radial                   !!!!���򻮷�������
    integer::coolant_kind               !!!!ָ����ȴ�����࣬1--ˮ��2--Ǧ��
    integer::module_identification      !!!!ģʽʶ��1--��̬��2--˲̬
    integer::transient_mode             !!!!˲ʱ���������ʽʶ��1--���ʣ�2--�����¶�
    real(8)::pi                         !!!!Բ���ʣ�����
    real(8)::p_line_factor              !!!!��������
    real(8)::p_line_average             !!!!ƽ���߹����ܶȣ���λW/m
    real(8)::coolant_V_flow             !!!!��ȴ����������λm3/s
    real(8)::coolant_T_in               !!!!��ȴ�������¶ȣ���λK
    real(8)::cladding_length            !!!!��ʼȼ��Ԫ���߶ȣ���λm
    real(8)::cladding_width             !!!!��ʼ���Ǻ�ȣ���λm
    real(8)::gas_gap                    !!!!��ʼ��϶��࣬��λm
    real(8)::pitch                      !!!!դ�࣬��λm
    real(8)::pellet_diameter            !!!!о���ʼֱ������λm
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!�����������!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!�����������!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!�����������!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!�м��������!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!�м��������!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!�м��������!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    parameter(pi=3.1415926)
    real,allocatable::p_line(:)                     !!!!�߹��ʣ���λW/m
    real,allocatable::Q_line(:)                     !!!!���������������λW
    real,allocatable::Temperature(:,:,:)            !!!!����ʱ�̣�����λ�ô��¶ȣ�n_radial+1��������ڱ����¶ȣ�n_radial+2�������������¶ȣ�n_radial+3������ȴ���¶ȣ���λK
    real,allocatable::d(:,:,:)                      !!!!����ʱ�̣�����λ�ô�ֱ����n_radial+1��������ھ���n_radial+2��������⾶����λm
    real,allocatable::Bu(:,:,:)                     !!!!����ʱ�̣�����λ�ô���ȼ�ģ���λMWD/kg
    real,allocatable::Bu_old(:,:,:)                 !!!!��һʱ�̣�����λ�ô���ȼ�ģ���λMWD/kg
    real,allocatable::k_fuel(:,:,:)                 !!!!����ʱ�̣�����λ�ô���о���ȵ��ʣ���λW/m*K
    real,allocatable::P_SQUARE(:,:,:)               !!!!����ʱ�̣�����λ�ô���������ʣ���λW/m3
    real(8)::p_line_max                             !!!!����߹��ʣ���λW/m
    real(8)::x                                      !!!!���ʺ�����βλ�ã���λm
    real(8)::sin_average                            !!!!����ƽ��ֵ
    real(8)::Fq                                     !!!!�ȹ����ӣ��������㣩
    real(8)::coolant_density                        !!!!��ȴ���ܶȣ���λkg/m3
    real(8)::coolant_T                              !!!!��ȴ���¶ȣ���λK
    real(8)::coolant_Cp                             !!!!��ȴ����ѹ�����ݣ���λJ/(kg*k)
    real(8)::fuel_percentage                        !!!!ȼ�ϰٷֱȣ���λ%
    real(8)::De                                     !!!!����ֱ������λm
    real(8)::coolant_M_flow                         !!!!��ȴ��������������λkg/s
    real(8)::coolant_S                              !!!!��ȴ����ͨ�������λm2
    real(8)::coolant_speed                          !!!!��ȴ�����٣���λm/s
    real(8)::Pr,Re,Pe,Nu                            !!!!������������ŵ������������
    real(8)::cladding_T_surface                     !!!!����������¶ȣ���λK
    real(8)::cladding_T,Tc,Tc1                      !!!!���ǵ��������м�ֵ����λK
    real(8)::cladding_T_internal                    !!!!�����ڱ����¶ȣ���λK
    real(8)::pellet_T_surface                       !!!!о������¶ȣ���λK
    real(8)::Bu_average                             !!!!ƽ��ȼ�ģ���λMWD/kg
    real(8)::UO2_idealdensity                       !!!!�������������ܶȣ���λkg/m3
    real(8)::fule_percentage                        !!!!ȼ�������ܶȣ���λ%
    real(8)::UO2_density                            !!!!���������ܶȣ���λkg/m3
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!�м��������!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!�м��������!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!�м��������!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!˲̬��������!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!˲̬��������!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!˲̬��������!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer t_number                                !!!!˲̬��������
    real(8)::AT                                     !!!!˲̬���㲽������λs
    real(8)::UO2_Cp                                 !!!!�������˶�ѹ�����ݣ���λJ/(kg*K)
    real(8)::heat_dv,heat_old                       !!!!���ȣ�����ѭ������
    real(8)::length_every                           !!!!����Ԫƽ�����ȣ���λm
    allocatable::Temperature_transient(:,:,:)       !!!!˲̬�¶ȣ���λK
    allocatable::T_transient(:)                     !!!!˲̬�����Ӧʱ�̣���λs
    allocatable::p_average_transient(:)             !!!!˲̬�����ʱ�̶�Ӧƽ�����ʣ���λW/m
    allocatable::coolant_T_transient(:)             !!!!˲̬���������ȴ���¶ȣ���λK




    !��������Ϊ˫���ȱ���
    double precision Temperature_transient
    double precision T_transient
    double precision p_average_transient
    double precision coolant_T_transient


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!˲̬��������!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!˲̬��������!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!˲̬��������!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    contains

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!���ļ���!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!���ļ���!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!���ļ���!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!���ļ���!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!���ļ���!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!���ļ���!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!�������!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!�������!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!�������!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine read_file
    integer i,j
    if (module_identification==1) then
        read (1001,*)
        read (1001,*) time,n_axis,n_radial,coolant_kind,p_line_factor,p_line_average
        read (1001,*)
        read (1001,*) coolant_V_flow,coolant_T_in,cladding_length,cladding_width,gas_gap,pitch
        read (1001,*)
        read (1001,*) pellet_diameter
        !---------------------------------------------------------------------------------������ʼ��
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
            !-----------------------------------------------------------------------------������ʼ��
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
            !-----------------------------------------------------------------------------������ʼ��
            call variable_return
            do i=1,n_axis
                read (1005,*) Temperature(1,i,n_radial+3)
            end do
        end if
    end if
    end subroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!�������!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!�������!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!�������!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!������ʼ��!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!������ʼ��!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!������ʼ��!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
    d(:,:,n_radial+2)=pellet_diameter+2*gas_gap+2*cladding_width    !!!!��ʼʱ�̣������⾶����λm
    d(:,:,n_radial+1)=pellet_diameter+2*gas_gap                     !!!!��ʼʱ�̣������ھ�����λm
    coolant_S=pitch**2.-pi/4.*(pellet_diameter+2*gas_gap+&          !!!!��ͨ�������λm2
    &2*cladding_width)**2.
    De=4.*(pitch**2.-pi/4.*(pellet_diameter+2*gas_gap+2*&           !!!!����ֱ������λm
    &cladding_width)**2.)/pi/(pellet_diameter+2*gas_gap&
        &+2*cladding_width)
    length_every=cladding_length/((n_axis-1)*1.0)
    Bu_old(:,:,:)=0
    do i=1,n_radial
        d(:,:,i)=(i-1)*pellet_diameter/(n_radial-1)
    enddo
    end subroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!������ʼ��!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!������ʼ��!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!������ʼ��!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



    end module variables
















