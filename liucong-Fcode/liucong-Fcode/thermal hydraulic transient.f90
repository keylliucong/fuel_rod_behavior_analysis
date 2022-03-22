module thermal_hydraulic_transient
    use variables
    use physical_properties
    use hcore
    use match
    implicit none
    contains
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!��ȴ���¶ȼ���!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!��ȴ���¶ȼ���!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!��ȴ���¶ȼ���!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine coolant_T_calculation_transient(time11,n_a)
    integer i
    integer n_a,time11
    real(8)::a,b
        if (n_a==1) then
            Q_line(n_a)=0.
            coolant_T_transient(n_a)=coolant_T_in
        else
            coolant_T_transient(n_a)=coolant_T_transient(n_a-1)+Q_line(n_a)/coolant_M_flow/coolant_Cp
            

    end if
    end subroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!��ȴ���¶ȼ������!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!��ȴ���¶ȼ������!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!��ȴ���¶ȼ������!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    end module thermal_hydraulic_transient
    