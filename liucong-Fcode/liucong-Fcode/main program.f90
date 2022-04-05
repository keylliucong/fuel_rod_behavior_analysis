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

















    !��̬���㿪ʼ

    !�ȹ����㿪ʼ
    if (model_identification==1) then
        plastic=0.
        strain_z_lasttime=0.
        strain_z_lasttime111=0.
        press_inter=5.
        press_begin=5.

        d_length_spring=0.02    !m����ʼ��ѹ
        d_length_p=0.
        d_length_c=0.
        strain_creep=0.
        strain_creep111=0.
        strain_creep_nb=0.
        strain_creep_nb111=0.
        do i=1,n_axis
            do j=1,n_radial+N_CLAD
                if (j<=n_radial) then
                    yield(i,j)=5.d16        !о�鿴����������ֱ�ӿ���
                else
                    yield(i,j)=5.1d8        !������������
                end if
            end do
        end do

        X_FC=gas_gap

        !ʱ��ѭ����ʼ
        do i=1,time     !��λday

            day=5.*i
            time_total=day*24.
            time_increment=5.*24.


            do j=1,n_axis                                           !!!!!!!!!!!!!!!����ѭ��

                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!���ʼ���!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!���ʼ���!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!���ʼ���!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                call power_calculation
                ! write(1002,*)"�������"
                !write(1002,"(f9.3)")(p_line(l),l=1,n_axis)
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!���ʼ������!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!���ʼ������!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!���ʼ������!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!��ȴ���¶ȼ���!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!��ȴ���¶ȼ���!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!��ȴ���¶ȼ���!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                !write(1002,*)"��ȴ���¶����"



                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!��ȴ���¶ȼ������!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!��ȴ���¶ȼ������!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!��ȴ���¶ȼ������!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!






                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!�����¶ȼ���!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!�����¶ȼ���!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!�����¶ȼ���!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




                do k=1,n_radial+N_clad
                    yield_stress(k)=yield(j,k)
                    strain_plastic(k,1)=plastic(j,k,1)
                    strain_plastic(k,2)=plastic(j,k,2)
                    strain_plastic(k,3)=plastic(j,k,3)
                end do
                do k=1,N_CLAD
                    strain_creep(k,1)=strain_creep_nb(j,k,1)
                    strain_creep(k,2)=strain_creep_nb(j,k,2)
                    strain_creep(k,3)=strain_creep_nb(j,k,3)
                end do
                !�����������ֵ
                Temperature_gap(j)=200.
                Temperature_JS(j)=0.
h_gass=5678.
                do while(ABS(Temperature_gap(j)-Temperature_JS(j))/Temperature_gap(j)>0.1)
                  
                    Temperature_JS(j)=Temperature_gap(j)
                    do k=1,n_radial
                        d(i,j,k)=(k-1)*(d(i,j,12)-2.*X_FC(i,j))/(n_radial-1)

                    end do
                    call coolant_physical_properties(j)
                    call coolant_T_calculation(j)           !!!!!!��ȴ���¶ȼ���
                    Temperature(:,j,n_radial+3)=coolant_T   !!!!!!��ȴ���¶ȸ���
                    !write(1002,"(f9.3)")coolant_T

                    call axis_T_calculation(i,j)                        !!!!!!!!!!!!!!!�����¶ȼ���
                    call axis_T_update(i,j)                             !!!!!!!!!!!!!!!�¶ȸ���
                    call Bu_calculation(i,j)                            !!!!!!!!!!!!!!!ȼ�ļ���
                    call P_S_calculation(i,j)
                    call fuel_t_calculation(i,j)                        !!!!!!!!!!!!!!!о���¶ȼ���
                    do k=1,n_radial+2
                        if (i==1) then
                            T_pre(k)=Temperature(i,j,k)
                            T_now(k)=T_pre(k)
                        else
                            T_pre(k)=Temperature(i-1,j,k)
                            T_now(k)=Temperature(i,j,k)
                        end if
                     
                        BU_BEGIN(k)=4.*p_LINE(j)/(pi*(d(i,j,12)-2.*X_FC(I,j))**2.)*(day-5.)/(UO2_DENSITY*238./270.)/(10.**6.)
                        BU_END(k)=4.*p_LINE(j)/(pi*(d(i,j,12)-2.*X_FC(I,j))**2.)*day/(UO2_DENSITY*238./270.)/(10.**6.)

                    end do

                    !�ȹ��������

                    !��е���㿪ʼ


                    call  MECHMODEL_unrigid(time_total,time_increment,p_line(j),do_original,di_original,T_PRE,T_NOW,gas_gap,BU_BEGIN,BU_END,n_radial,PRESS_INTER,PRESS_BEGIN,UR_PRA,P_contact(j),&
                        &N_CLAD,d_length_spring,d_length_p(j),d_length_c(j),yield_stress,yield_stress111,strain_plastic,strain_plastic111&
                        &,aaa,strain_z_lasttime,strain_z_lasttime111,stress_equ,strain_fuel,strain_creep,strain_creep111&
                        &,stress_cladding_z111,strain_cladding_z111)

                    !if (2.*UR_PRA(12)<d(i,j,11))then
                    !    ur_pra(12)=d(i,j,11)
                    !    pause
                    !    end if

                    if (Bu_AVERAGE<=5.) then
                        B_A(j)=Bu_AVERAGE/5.
                    else
                        B_A(j)=1.
                    end if
                    if (p_line(j)<=20000.) then
                        X_FC(i,j)=UR_PRA(12)-UR_PRA(11)-GAS_GAP*(30.+10.*B_A(j)/5.)/100.
                    else if (p_line(j)<=40000.) then
                        X_FC(i,j)=UR_PRA(12)-UR_PRA(11)-GAS_GAP*(28.+(p_line(j)-20000.)/4000.+(12.+(p_line(j)-20000.)/4000.)*B_A(j)/5.)/100.
                    else
                        X_FC(i,j)=UR_PRA(12)-UR_PRA(11)-GAS_GAP*(32.+18.*B_A(j)/5.)/100.
                    end if
                    IF(X_FC(I,j).LE.(5.*10**(-9.))) THEN !�ڲ�ѹ�����󣬿��ܻ����PC�ڷ�������
                        X_FC(I,j)=0.
                    END IF
                    
call h_gas_calculation(i,j)

                    Temperature_gap(j)=q_line(j)/(pi*d(i,j,12))/h_gass(i,j)

                    d(i,j,12)=2.*UR_PRA(12)


                end do !end dowhile




                stress_cladding_z(j)=stress_cladding_z111
                strain_cladding_z(j)=strain_cladding_z111

                strain_creep=strain_creep111


                do k=1,N_CLAD
                    strain_creep_nb(j,k,1)=strain_creep(k,1)
                    strain_creep_nb(j,k,2)=strain_creep(k,2)
                    strain_creep_nb(j,k,3)=strain_creep(k,3)
                end do

                strain_z_lasttime=strain_z_lasttime111

                do k=1,N_radial+N_CLAD
                    yield(j,k)=yield_stress111(k)
                    plastic(j,k,1)=strain_plastic111(k,1)
                    plastic(j,k,2)=strain_plastic111(k,2)
                    plastic(j,k,3)=strain_plastic111(k,3)

                    stress_equ_wh(j,k)=stress_equ(k)
                    if (k<=n_radial) then
                        strain_fuel_wh(j,k,1)=strain_fuel(k,1)
                        strain_fuel_wh(j,k,2)=strain_fuel(k,2)
                        strain_fuel_wh(j,k,3)=strain_fuel(k,3)
                    end if

                end do





            end do

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!�����¶ȼ������!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!�����¶ȼ������!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!�����¶ȼ������!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



            d_length_spring=0.02
            do j=1,n_axis
                d_length_spring=d_length_spring+d_length_p(j)-d_length_c(j)
            end do









        enddo
        !ʱ��ѭ������

















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