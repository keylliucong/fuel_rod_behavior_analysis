    module thermal_hydraulic
    use variables
    use physical_properties
    use hcore
    use match
    implicit none
    contains



    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!���ʼ���!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!���ʼ���!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!���ʼ���!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine power_calculation
    integer i
    real(8)::a,b,add
add=0.
    p_line_max=p_line_average*p_line_factor
    p_line_max=16000.
    add=0.
    x=0.                                !!!!�ٶ����ʺ�����βλ��
    Fq=0.                               !!!!������ֵ
    do while(abs(Fq-p_line_factor)>0.005)
        sin_average=(cos(x)-cos(pi-x))/(pi-2*x)
        Fq=1./sin_average
        if (Fq-p_line_factor>0) then
            x=x+0.001
        else
            x=x-0.001
        end if
    end do
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !sin_average=(cos(0.)-cos(pi))/(pi)
    !    Fq=1./sin_average   !!!!!!!!!!!!!!!!!��֤˲̬���򣬹��ʷֲ���ȡ����
    !    x=0.
    !    p_line_max=p_line_average*Fq
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        do i=1,n_axis
        p_line(i)=p_line_max*sin(x+(i-1)*(pi-2*x)/(1.0*(n_axis-1)))
        p_line(i)=p_line_max*factor(i)
       add=add+p_line(i)
        end do
     do i=1,n_axis
         p_line(i)=p_line_max*p_line(i)*n_axis/add
         end do
    do i=1,n_axis
        a=0
        b=0
        a=(i-1)*cladding_length/((n_axis-1)*1.0)
        b=(i-2)*cladding_length/((n_axis-2)*1.0)
        Q_line(i)=p_line_max*cladding_length/(pi-2*x)*(cos(x+b*(pi-2*x)/cladding_length)-cos(x+a*(pi-2*x)/cladding_length))
    end do
    
    end subroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!���ʼ������!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!���ʼ������!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!���ʼ������!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!��ȴ���¶ȼ���!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!��ȴ���¶ȼ���!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!��ȴ���¶ȼ���!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine coolant_T_calculation(n_a)
    integer i
    integer n_a
    real(8)::a,b
    if (model_identification==1)then
        a=0
        a=(n_a-1)*cladding_length/((n_axis-1)*1.0)
        Q_line(n_a)=p_line_max*cladding_length/(pi-2*x)*(cos(x)-cos(x+a*(pi-2*x)/cladding_length))
        coolant_T=Q_line(n_a)/coolant_M_flow/coolant_Cp+coolant_T_in
    else
        if (n_a==1) then
            Q_line(n_a)=0.
            coolant_T_transient(n_a)=coolant_T_in
        else
            coolant_T_transient(n_a)=coolant_T_transient(n_a-1)+Q_line(n_a)/coolant_M_flow/coolant_Cp
            
        end if
    end if
    end subroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!��ȴ���¶ȼ������!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!��ȴ���¶ȼ������!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!��ȴ���¶ȼ������!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!�����¶ȼ���!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!�����¶ȼ���!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!�����¶ȼ���!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine axis_T_calculation(time,n_a)
    integer i
    integer time,n_a
    cladding_T_surface=Temperature(time,n_a,n_radial+3)+p_line(n_a)&
        &/(pi*d(time,n_a,n_radial+2)*h_coefficient(time,n_a))
    cladding_T=cladding_T_surface                                       !!!!!!��������ڱ����¶�
    cladding_T_internal=cladding_T_surface+P_line(n_a)*log(&            !!!!!�����ڱ����ʼ����ֵ
    &d(time,n_a,n_radial+2)/(d(time,n_a,n_radial+2)-2*cl&
        &adding_width))/(2*pi*Kc(cladding_T,cladding_T_surface))
    Tc1=cladding_T_internal-cladding_T                                  !!!!!��ʼ�²�
    Tc=0                                                                !!!!!!�����²�
    !���ǵ������㿪ʼ
    do while((Tc1-Tc)/Tc1>0.005)
        Tc=Tc1
        cladding_T_internal=cladding_T_surface+p_line(n_a)&
            &*log(d(time,n_a,n_radial+2)/(d(time,n_a,n_radial+2)&
            &-2*cladding_width))/(2*pi*Kc(cladding_T,cladding_T_surface))
        Tc1=abs(cladding_T-cladding_T_internal)
        cladding_T=cladding_T_internal
    enddo
    !���ǵ����������
    pellet_T_surface=cladding_T_internal+p_line(n_a)/pi*d(time,n_a,n_radial)/h_gass(time,n_a)
    
    
    end subroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!�����¶ȼ������!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!�����¶ȼ������!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!�����¶ȼ������!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!�¶ȸ���!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!�¶ȸ���!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!�¶ȸ���!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine axis_T_update(time,n_a)
    integer i
    integer time,n_a
    Temperature(time,n_a,n_radial+2)=cladding_T_surface
    Temperature(time,n_a,n_radial+1)=cladding_T_internal
    do i=1,n_radial
        Temperature(time,n_a,i)=pellet_T_surface
    end do
    end subroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!�¶ȸ���!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!�¶ȸ���!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!�¶ȸ���!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!�����¶ȼ���!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!�����¶ȼ���!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!�����¶ȼ���!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine fuel_t_calculation(time1,n_a)
    integer i
    integer time1,n_a
    dimension A(n_radial,n_radial),B(n_radial),X(n_radial),Tc(n_radial)
    double precision A,B,X,Tc
    real(8)::tt
    X(:)=0                                                  !!!!�����ʼ��
    tt=100.
    A(:,:)=0.
    do while(tt>0.01)
        do i=1,n_radial
            call fuel_coefficient_calculation(time1,n_a,i)

        end do
        !!!!!!!ϵ������ֵ��ʼ!!!!!!!!

        A(1,1)=1.
        A(1,2)=-1.
        B(1)=0.5*P_SQUARE(TIME1,n_a,1)*(d(time1,n_a,2)**2)/4./k_fuel(time1,n_a,1)
        do i=2,n_radial-1
            A(i,i-1)=(d(time1,n_a,i)+d(time1,n_a,i-1))/(d(time1,n_a,i)-d(time1,n_a,i-1))/2.*k_fuel(time1,n_a,i-1)
            A(i,i+1)=(d(time1,n_a,i)+d(time1,n_a,i+1))/(d(time1,n_a,i+1)-d(time1,n_a,i))/2.*k_fuel(time1,n_a,i)
            A(i,i)=-(A(i,i-1)+A(i,i+1))
            B(i)=-1.*P_SQUARE(TIME1,n_a,i)*d(time1,n_a,i)/2.*(d(time1,n_a,i+1)-d(time1,n_a,i-1))/4.
           
        end do
        A(n_radial,n_radial)=1.
        B(n_radial)=pellet_T_surface
        !!!!!!!ϵ������ֵ����!!!!!!!!
        call TDMA(n_radial,A,B,X)                            !!!!!TDMA���㣬����������ȡ�ü������о���¶�ֵ
        do i=1,n_radial
            Tc(i)=x(i)-Temperature(time1,n_a,i)
            Temperature(time1,n_a,i)=x(i)
        end do
        tt=maxval(abs(Tc))
    end do
    end subroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!�����¶ȼ������!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!�����¶ȼ������!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!�����¶ȼ������!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!













    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!�¶ȼ���˲̬!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!�¶ȼ���˲̬!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!�¶ȼ���˲̬!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine T_transient_calculation(time11,n_a)
    integer i,time11,n_a
    dimension A(n_radial+2,n_radial+2),B(n_radial+2),X(n_radial+2)
    dimension Tc(n_radial+2),transient(n_radial+2),Temperature_mid(t_number,n_axis,n_radial+3)
    double precision A,B,X,transient,Temperature_mid,Tc
    real(8)::tt,rw,re,d_rw,d_re,heat
    cladding_T_surface=Temperature_transient(time11,n_a,n_radial+3)+p_line(n_a)&
        &/(pi*d(time11,n_a,n_radial+2)*h_coefficient(time11,n_a))
    
    do i=1,n_radial+2
        Temperature_transient(time11,n_a,i)=cladding_T_surface

    end do








    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    X(:)=0                                                  !!!!�����ʼ��
    tt=100.
    A(:,:)=0.
    do while(tt>0.01)



        do i=1,n_radial
            call fuel_coefficient_calculation(time11,n_a,i)
        end do




        !!!!!!!ϵ������ֵ��ʼ!!!!!!!!
        do i=1,n_radial+3
            if (time11==1) then

                Temperature_mid(time11,n_a,i)=Temperature(time,n_a,i)
            else
                Temperature_mid(time11,n_a,i)=Temperature_transient(time11-1,n_a,i)
            end if
        end do
        transient(1)=UO2_density*UO2_Cp*(d(time11,n_a,2)/4.)**2./AT/2.
        A(1,1)=-(1.*k_fuel(time11,n_a,1)/2.+transient(1))
        A(1,2)=1.*k_fuel(time11,n_a,1)/2.
        B(1)=-(0.5*P_SQUARE(time11,n_a,1)*(d(time11,n_a,2)/4.)**2.+transient(1)*Temperature_mid(time11,n_a,1))
        do i=2,n_radial+1
            rw=(d(time11,n_a,i)+d(time11,n_a,i-1))/4.
            re=(d(time11,n_a,i)+d(time11,n_a,i+1))/4.
            d_rw=(d(time11,n_a,i)-d(time11,n_a,i-1))/2.
            d_re=(d(time11,n_a,i+1)-d(time11,n_a,i))/2.
            if (i<n_radial) then
                transient(i)=UO2_density*UO2_Cp/AT*(re-rw)*(re+rw)/2.
                A(i,i-1)=rw/d_rw*k_fuel(time11,n_a,i-1)
                A(i,i+1)=re/d_re*k_fuel(time11,n_a,i)
                A(i,i)=-(A(i,i-1)+A(i,i+1)+transient(i))
                B(i)=-1.*(P_SQUARE(time11,n_a,i)*(re-rw)*(re+rw)/2.+transient(i)*Temperature_mid(time11,n_a,i))
            else if (i==n_radial) then
                A(i,i-1)=k_fuel(time11,n_a,i)/d_re
                A(i,i+1)=h_gass(time11,n_a)
                A(i,i)=-(A(i,i-1)+A(i,i+1))
                B(i)=0
            else if (i==n_radial+1) then
                A(i,i-1)=h_gass(time11,n_a)*d(time11,n_a,n_radial)/d(time11,n_a,n_radial+1)
                A(i,i+1)=Kc(Temperature_transient(time11,n_a,i),cladding_T_surface)/d_rw
                A(i,i)=-(A(i,i-1)+A(i,i+1))
                B(i)=0
            end if

        end do
        
        
        A(n_radial+2,n_radial+2)=1.
        B(n_radial+2)=cladding_T_surface

        !!!!!!!ϵ������ֵ����!!!!!!!!
        call TDMA(n_radial+2,A,B,X)                            !!!!!TDMA���㣬����������ȡ�ü������о���¶�ֵ
        do i=1,n_radial+2
            Tc(i)=x(i)-Temperature_transient(time11,n_a,i)
            Temperature_transient(time11,n_a,i)=x(i)
        end do
        tt=maxval(abs(Tc))
    end do







    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!���㴢��
    Tc(:)=0.
    heat=0.
    do i=1,n_radial+3
        if (time11==1) then

            Temperature_mid(time11,n_a,i)=Temperature(time,n_a,i)
        else
            Temperature_mid(time11,n_a,i)=Temperature_transient(time11-1,n_a,i)
        end if
    end do
    do i=1,n_radial+2
        Tc(i)=Temperature_transient(time11,n_a,i)-Temperature_mid(time11,n_a,i)
        if (i==1) then
            heat=heat+Tc(i)*UO2_Cp*UO2_density*pi*((d(time11,n_a,i)/2.)**2)/2.
        else if(i<n_radial) then
            rw=(d(time11,n_a,i)+d(time11,n_a,i-1))/4.
            re=(d(time11,n_a,i)+d(time11,n_a,i+1))/4.
            heat=heat+Tc(i)*UO2_Cp*UO2_density*pi*(re**2-rw**2)
        else if(i==n_radial) then
            rw=(d(time11,n_a,i)+d(time11,n_a,i-1))/4.
            heat=heat+Tc(i)*UO2_Cp*UO2_density*pi*(d(time11,n_a,i)**2/4.-rw**2)
        else
            rw=d(time11,n_a,n_radial+1)/2.
            re=d(time11,n_a,n_radial+2)/2.
            heat=heat+Tc(i)*UO2_Cp*UO2_density*pi*(re**2-rw**2)/2.
        end if
    end do


    do i=2,n_axis
        Q_line(i)=Q_line(i)-heat*length_every        !��������
    end do
    heat_dv=heat-heat_old                       !������ֵ


    heat_old=heat




    end subroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!�¶ȼ���˲̬!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!�¶ȼ���˲̬!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!�¶ȼ���˲̬!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
    
    
    
    
    
    
    end module