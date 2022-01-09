    program main
    use variables
    use physical_properties
    use thermal_hydraulic
    use hcore
    use match
    use thermal_hydraulic_transient
    implicit none
    integer i,j,k,l

    !print*,"������������ͣ�1--��̬��2--˲̬"
    !read *,model_identification
    !if (model_identification==2) then
    !    print*,"��ѡ����֪�������ͣ�1--���ʣ�2--�����¶�"
    !    read *,transient_mode
    !end if

    call model_select
    call open_file                          !!!!!!!!!!���ļ���
    call read_file                          !!!!!!!!!!��������


    call fuel_physical_properties




    if (model_identification==1) then
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!���ʼ���!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!���ʼ���!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!���ʼ���!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        call power_calculation
        write(1002,*)"�������"
        write(1002,"(f9.3)")(p_line(i),i=1,n_axis)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!���ʼ������!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!���ʼ������!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!���ʼ������!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!��ȴ���¶ȼ���!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!��ȴ���¶ȼ���!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!��ȴ���¶ȼ���!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        write(1002,*)"��ȴ���¶����"
        do i=1,n_axis
            call coolant_physical_properties(i)
            call coolant_T_calculation(i)           !!!!!!��ȴ���¶ȼ���
            Temperature(:,i,n_radial+3)=coolant_T   !!!!!!��ȴ���¶ȸ���
            write(1002,"(f9.3)")coolant_T
        end do
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!��ȴ���¶ȼ������!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!��ȴ���¶ȼ������!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!��ȴ���¶ȼ������!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!�����¶ȼ���!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!�����¶ȼ���!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!�����¶ȼ���!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        do i=1,time                                                 !!!!!!!!!!!!!!!ʱ��ѭ��
            do j=1,n_axis                                           !!!!!!!!!!!!!!!����ѭ��
                call axis_T_calculation(i,j)                        !!!!!!!!!!!!!!!�����¶ȼ���
                call axis_T_update(i,j)                             !!!!!!!!!!!!!!!�¶ȸ���
                call Bu_calculation(i,j)                            !!!!!!!!!!!!!!!ȼ�ļ���
                call P_S_calculation(i,j)
                call fuel_t_calculation(i,j)                        !!!!!!!!!!!!!!!о���¶ȼ���
            end do
        enddo
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!�����¶ȼ������!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!�����¶ȼ������!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!�����¶ȼ������!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!�������!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!�������!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!�������!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        write(1002,*)"��������¶�"
        do i=1,time
            do j=1,n_axis
                write(1002,55)(Temperature(i,j,k),k=1,n_radial+3)
                write(1002,56)(d(i,j,k)/2,k=1,n_radial+2)
55              format(14f9.3)
56              format(13f9.6)
            end do
        end do
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!�������!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!�������!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!�������!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





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

                !!!!!!!!!!!�����ʼ����ʱȼ��begin �������޸�ֱ�ӵ�����̬����ֵ
                do k=1,time
                    do l=1,n_axis
                        call Bu_calculation(k,l)
                    end do
                enddo
                !!!!!!!!!!!�����ʼ����ʱȼ��begin �������޸�ֱ�ӵ�����̬����ֵ

                !!!!!!!!!˲̬�¶ȼ��㲿��begin

                do j=1,n_axis
                    call P_S_calculation(i,j)
                    call coolant_physical_properties(j)
                    heat_dv=100.
                    heat_old=0.
                    ! do while(abs(heat_dv)>0.1)          !������������������
                    call coolant_T_calculation(j)   !��ȴ���¶ȼ���
                    Temperature_transient(i,j,n_radial+3)=coolant_T_transient(j)        !��ֵ���洢��ȴ���¶�

                    call T_transient_calculation(i,j)   !�¶ȼ��㣬TDMA
                    !end do

                end do

                !!!!!!!!!˲̬�¶ȼ��㲿��end



                write(1004,60)Temperature_transient(i,11,1),Temperature_transient(i,11,n_radial)        !�������
60              format(2f9.3)


            else
                p_line_average=p_average_transient(i)
                call power_calculation
                !!!!!!!!!!!�����ʼ����ʱȼ��begin �������޸�ֱ�ӵ�����̬����ֵ
                do k=1,time
                    do l=1,n_axis
                        call Bu_calculation(k,l)
                    end do
                enddo
                !!!!!!!!!!!�����ʼ����ʱȼ��begin �������޸�ֱ�ӵ�����̬����ֵ

                !!!!!!!!!˲̬�¶ȼ��㲿��begin

                do j=1,n_axis
                    call P_S_calculation(i,j)
                    call coolant_physical_properties(j)
                    !heat_dv=100.
                    !heat_old=0.
                    !do while(abs(heat_dv)>0.1)          !������������������
                    call coolant_T_calculation(j)   !��ȴ���¶ȼ���
                    Temperature_transient(i,j,n_radial+3)=coolant_T_transient(j)        !��ֵ���洢��ȴ���¶�
                    call T_transient_calculation(i,j)   !�¶ȼ��㣬TDMA
                    !end do
                end do

                !!!!!!!!!˲̬�¶ȼ��㲿��end

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