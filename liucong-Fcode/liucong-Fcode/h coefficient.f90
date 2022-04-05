module hcore
    use variables 
    use physical_properties 
    implicit none
    contains
    
    
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!对流换热系数!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!对流换热系数!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!对流换热系数!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function h_coefficient(t,n_axis)                 !!!!t代表时间t，n_axis代表轴向划分网格数
    double precision h_coefficient
    integer t,n_axis
    if (model_identification==1) then
        if (coolant_kind==1) then
            Re=(11096.0-1.3236*Temperature(t,n_axis,n_radial+3))*coolant_speed*De/(4.94e-4*exp(757.1/&
                &Temperature(t,n_axis,n_radial+3)))
            Pr=(159.0-0.0272*Temperature(t,n_axis,n_radial+3)+7.12e-6*Temperature(t,n_axis,n_radial+3&
                &)**2)*(4.94e-4*exp(757.1/Temperature(t,n_axis,n_radial+3)))/(3.61+0.01517*Temperatur&
                &e(t,n_axis,n_radial+3))
            Pe=Re*Pr
            Nu=7.55*(pitch/d(t,n_axis,n_radial+2))-20.*((pitch/d(t,n_axis,n_radial+2))**-13)+3.67/90.&
                &*((pitch/d(t,n_axis,n_radial+2))**-2)*(Pe**(0.56+0.19*(pitch/d(t,n_axis,n_radial+2))))
            h_coefficient=Nu*(3.61+0.01517*Temperature(t,n_axis,n_radial+3))/De
        end if
    else
            Re=(11096.0-1.3236*Temperature_transient(t,n_axis,n_radial+3))*coolant_speed*De/(4.94e-4*exp(757.1/&
                &Temperature_transient(t,n_axis,n_radial+3)))
            Pr=(159.0-0.0272*Temperature_transient(t,n_axis,n_radial+3)+7.12e-6*Temperature_transient(t,n_axis,n_radial+3&
                &)**2)*(4.94e-4*exp(757.1/Temperature_transient(t,n_axis,n_radial+3)))/(3.61+0.01517*Temperatur&
                &e_transient(t,n_axis,n_radial+3))
            Pe=Re*Pr
            Nu=7.55*(pitch/d(t,n_axis,n_radial+2))-20.*((pitch/d(t,n_axis,n_radial+2))**-13)+3.67/90.&
                &*((pitch/d(t,n_axis,n_radial+2))**-2)*(Pe**(0.56+0.19*(pitch/d(t,n_axis,n_radial+2))))
            h_coefficient=Nu*(3.61+0.01517*Temperature_transient(t,n_axis,n_radial+3))/De
    end if
    end function
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!对流换热系数!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!对流换热系数!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!对流换热系数!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!包壳热导率!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!包壳热导率!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!包壳热导率!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function Kc(i,j)   
    real(8)::i,j,Kc
     Kc=(27.67 + 1.39d-2*((i+j)/2.) - 2.42d-5*((i+j)/2.)**2 + 8.29d-9*((i+j)/2.)**3)
    end function
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!包壳热导率!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!包壳热导率!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!包壳热导率!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    
    
   
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!气隙热导率Kg!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!气隙热导率Kg!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!气隙热导率Kg!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine h_gas_calculation(time1,n_a)
    integer::time1,n_a
    real(8)::C_HE,AHE,AXE,MHE,X_JUMP,RF,RC,X_DEFF,H_GAS
    real(8)::FA,EF,EC,RCI,RFS,FE,H_R
    real(8)::H_SOLID,PRESS_RATIO,FST_JS
    real(8)::S,P,B,k0,fd,fp,fm,fr,fuek
    real(8)::RMULT,CONDUCTIVITY_MEAN
    real(8)::T,TC
    
    
    
    S=1.5           !!!!!!!!!!形状因子
    P=0.05          !!!!!!!!!!芯块孔隙率
    
    
    
    
    
     C_HE=2.639*0.001*((Temperature(time1,n_a,11)+Temperature(time1,n_a,12))/2.)**0.7085
    !C_HE=0.158*0.01*((Temperature(time1,n_a,11)+Temperature(time1,n_a,12))/2.)**0.79                      ! Kampf和Karsten关系式 
    AHE=0.425-2.3*10.**(-4.)*((Temperature(time1,n_a,11)+Temperature(time1,n_a,12))/2.)                    ! accommodation coefficient for helium
    AXE=0.749-2.5*10.**(-4.)*((Temperature(time1,n_a,11)+Temperature(time1,n_a,12))/2.)                    ! accommodation coefficient for xenon
    MHE=4.                                                                           ! molecular weight of helium
    X_JUMP=0.024688*(C_HE*(((Temperature(time1,n_a,11)+Temperature(time1,n_a,12))/2.)**0.5)/(PRESS_BEGIN*(10.**6.))*1./AHE*MHE**0.5)    ! temperature jump distance (m)

    RF=2.0*10.**(-6.)                                                                 !fuel surface roughness
    RC=1.0*10.**(-6.)                                                                 !clad surface roughness
    !X_DEFF=RF+RC    
    X_DEFF=1.5*(RF+RC)                          
    IF(X_FC(time1,n_a)==0.) THEN
      !X_DEFF=EXP(-0.00125*PRESS_BEGIN*(10.**6.))*(RF+RC)                              !closed fuel-cladding gaps(FRAPCON)
      X_DEFF=1.98*EXP(-0.000088*(PRESS_BEGIN*(10.**6.))/(10.**6.)*145.)*(RF+RC)        !closed fuel-cladding gaps(GAPCON THERMAL 2:A COMPUTER PROGRAM FOR CALCULATING THE THERMAL BEHAVIOR OF AN OXIDE FUEL ROD;P20)
    END IF
    H_GAS=C_HE/(X_JUMP+X_DEFF+X_FC(time1,n_a))
    
    FA=1.0                                                                           ! configuration factor
    EF=0.7856+1.5263*10.**(-5.)*Temperature(time1,n_a,11)/2.                                       ! fuel emissivity
    IF(Temperature(time1,n_a,11)<1500.) THEN
      EC=0.325                                                                         ! cladding surface temperature has not exceeded 1500K
    ELSE 
      EC=0.325*EXP((Temperature(time1,n_a,11)-1500.)/300.)
    END IF
    RCI=d(time1,n_a,12)                                                                      ! 包壳内径
    RFS=RCI-2.*X_FC(time1,n_a)                                                               ! 芯块外径
    FE=1./(1./EF+(RFS/RCI)*(1./EC-1.))                                               ! emissivity factor 
    H_R=5.6697*10.**(-8.)*FE*FA*(Temperature(time1,n_a,12)**2.+(Temperature(time1,n_a,11))**2.)*2.*(Temperature(time1,n_a,11)+Temperature(time1,n_a,12))/2.  ! 由于导热，TGAS=(TCI+TFS)/2.


    ! SOLID CONTACT
    H_SOLID=0.
    IF(X_FC(time1,n_a)==0.)THEN
      PRESS_RATIO=(PRESS_BEGIN*(10.**6.))/680./(10.**6.)
      FST_JS=(Temperature(time1,n_a,11)+Temperature(time1,n_a,12))/2.
      
      
      
   k0=1./(0.0375+(2.165d-4)*FST_JS)+4.715d9*exp(-16361./FST_JS)/FST_JS**2
        !!!!!!!!!!!!!!!!k0-未考虑燃耗作用的二氧化铀导热系数
        FD=(1.09/(Bu_average/9.383)**3.265+0.0643/((Bu_average/9.383)**0.5)*(FST_JS**0.5)&
         &)*ATAN2(1.0,(1.09/((Bu_average/9.383)**3.265)+0.0643/((Bu_average/9.383)**0.5)*(FST_JS**0.5)))
        !!!!!!!!!!!!!!!!FD-溶解因子
        FP=1+(0.019*(Bu_average/9.383)/(3.-0.019*(Bu_average/9.383)))*(1./(1.+EXP(-(FST_JS-1200)/100.)))
        !!!!!!!!!!!!!!!!FP-沉淀因子
        FM=(1-P)/(1+(S-1)*P)!!!!!!!!!FM-孔隙因子
        FR=1.-0.2/(1.+EXP((FST_JS-900.)/80.))!!!!!!!!!!!!!FR辐照因子
        fuek=k0*FD*FM*FR               !!!!!!!!二氧化铀真实导热系数
        
         T=FST_JS-273.15D0
         TC=7.73D0+3.15D-2*T-2.87D-5*T*T+1.552D-8*T*T*T
         
        
        
    CONDUCTIVITY_MEAN=2.*fuek*TC/(TC+fuek)

    IF(PRESS_RATIO.GT.0.0087) THEN
      RMULT=2.9
    ELSE 
      RMULT=333.3*PRESS_RATIO
    END IF

    IF(PRESS_RATIO.GT.0.003) THEN
      H_SOLID=5.0*CONDUCTIVITY_MEAN*PRESS_RATIO*RMULT/((RF**2.+RC**2.)**0.5)/EXP(5.738)
    ELSE IF(PRESS_RATIO.LT.0.000009) THEN
      H_SOLID=5.0*CONDUCTIVITY_MEAN*PRESS_RATIO**0.5/((RF**2.+RC**2.)**0.5)/EXP(5.738)
    ELSE
      H_SOLID=0.0150*CONDUCTIVITY_MEAN/((RF**2.+RC**2.)**0.5)/EXP(5.738-0.528*log(1.*10.**(-6.)))
    END IF
    !H_SOLID=1.189*CONDUCTIVITY_MEAN*PRESS_RATIO/100./(((RF*100.)**2.+(RC*100.)**2.)**0.25)*10000. ! (Ross-Stoute Model)
    END IF
    ! SOLID CONTACT

    h_gass(time1,n_a)=H_R+H_GAS+H_SOLID
    
    end subroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!气隙热导率Kg!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!气隙热导率Kg!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!气隙热导率Kg!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!燃耗与体积功率计算Bu&P!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!燃耗与体积功率计算Bu&P!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!燃耗与体积功率计算Bu&P!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine Bu_calculation(time1,n_a)
    integer n_a,time1,n_r,day,i
    real(8)::sum_q
    sum_q=0.
    Bu_average=4.*p_LINE(n_a)/(pi*(pellet_diameter)**2.)*time1*5./(UO2_DENSITY*238./270.)/(10.**6.)   
  
    Bu(time1,n_a,:)=Bu_average
   
    
    if (time1==1) then
        Bu(time1,n_a,1)=Bu(time1,n_a,n_radial-2)*(1.-Bu(time1,n_a,n_radial-2)/300.)
        sum_q=sum_q+Bu(time1,n_a,1)*pi*(0.1**2.)
        do n_r=2,n_radial-2
            Bu(time1,n_a,n_r)=Bu(time1,n_a,n_r-1)+(Bu(time1,n_a,n_radial-2)-Bu(time1,n_a,1))/8.
            sum_q=sum_q+Bu(time1,n_a,n_r)*pi*((0.1*n_r)**2.-(0.1*(n_r-1))**2.)
        end do
        Bu(time1,n_a,n_radial-1)=(Bu(time1,n_a,n_radial-2)*pi*(1.0**2.)-sum_q)/(pi*(1.-0.9**2.))
     
    else
        do i=1,n_radial
            Bu_old(time1,n_a,i)=Bu(time1-1,n_a,i)
        end do
        Bu(time1,n_a,1)=Bu(time1,n_a,n_radial-2)*(1.-Bu(time1,n_a,n_radial-2)/300.)
        sum_q=sum_q+Bu(time1,n_a,1)*pi*(0.1**2.)
        do n_r=2,n_radial-2
            Bu(time1,n_a,n_r)=Bu(time1,n_a,n_r-1)+(Bu(time1,n_a,n_radial-2)-Bu(time1,n_a,1))/8.
            sum_q=sum_q+Bu(time1,n_a,n_r)*pi*((0.1*n_r)**2.-(0.1*(n_r-1))**2.)
        end do
        Bu(time1,n_a,n_radial-1)=(Bu(time1,n_a,n_radial-2)*pi*(1.0**2.)-sum_q)/(pi*(1.-0.9**2.))
     
        do i=1,n_radial
            Bu(time1,n_a,i)=Bu(time1,n_a,i)+Bu_old(time1,n_a,i)

        end do
        
    end if

    end subroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    
    subroutine P_S_calculation(time1,n_a)
    integer time1,n_a,n_r
    real(8)::P_SQUARE_AVERAGE,sum_p
    sum_p=0.
    p_square_average=4.*P_line(n_a)/(pi*(pellet_diameter**2))
    P_SQUARE(:,n_a,:)=p_square_average
    P_SQUARE(time1,n_a,1)=P_SQUARE(time1,n_a,n_radial-2)*(1-Bu_average/300.)
    !sum_p=sum_p+P_SQUARE(time1,n_a,1)*pi*(0.1**2.)
    do n_r=2,n_radial-2
        P_SQUARE(time1,n_a,n_r)=P_SQUARE(time1,n_a,n_r-1)+(P_SQUARE(time1,n_a,n_radial-2)-P_SQUARE(time1,n_a,1))/8.
        !sum_p=sum_p+P_SQUARE(time1,n_a,n_r)*pi*((0.1*n_r)**2.-(0.1*(n_r-1))**2.)
    end do
    P_SQUARE(time1,n_a,n_radial-1)=(P_SQUARE(time1,n_a,n_radial-2)*pi*(1.0**2.)-sum_p)/(pi*(1.-0.9**2.))
    end subroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!燃耗与体积功率计算Bu&P!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!燃耗与体积功率计算Bu&P!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!燃耗与体积功率计算Bu&P!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    
    
    
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!UO2芯块热导率计算!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!UO2芯块热导率计算!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!UO2芯块热导率计算!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine fuel_coefficient_calculation(time1,n_a,n_r)
    integer time1,n_a,n_r
    real(8)::S,P,B,k0,fd,fp,fm,fr
    S=1.5           !!!!!!!!!!形状因子
    P=0.05          !!!!!!!!!!芯块孔隙率
    if (model_identification==1) then
       
        k0=1./(0.0375+(2.165d-4)*Temperature(time1,n_a,n_r))+4.715d9*exp(-16361./Temperature(time1,n_a,n_r))/Temperature(time1,n_a,n_r)**2
        !!!!!!!!!!!!!!!!k0-未考虑燃耗作用的二氧化铀导热系数
        FD=(1.09/(Bu_average/9.383)**3.265+0.0643/((Bu_average/9.383)**0.5)*(Temperature(time1,n_a,n_r)**0.5)&
         &)*ATAN2(1.0,(1.09/((Bu_average/9.383)**3.265)+0.0643/((Bu_average/9.383)**0.5)*(Temperature(time1,n_a,n_r)**0.5)))
        !!!!!!!!!!!!!!!!FD-溶解因子
        FP=1+(0.019*(Bu_average/9.383)/(3.-0.019*(Bu_average/9.383)))*(1./(1.+EXP(-(Temperature(time1,n_a,n_r)-1200)/100.)))
        !!!!!!!!!!!!!!!!FP-沉淀因子
        FM=(1-P)/(1+(S-1)*P)!!!!!!!!!FM-孔隙因子
        FR=1.-0.2/(1.+EXP((Temperature(time1,n_a,n_r)-900.)/80.))!!!!!!!!!!!!!FR辐照因子
        k_fuel(time1,n_a,n_r)=k0*FD*FM*FR               !!!!!!!!二氧化铀真实导热系数
       
        
    else
        !k0=1./(0.0375+(2.165d-4)*Temperature_transient(time1,n_a,n_r))+4.715d9*exp(-16361./Temperature_transient(time1,n_a,n_r))/Temperature_transient(time1,n_a,n_r)**2
        !!!!!!!!!!!!!!!!!k0-未考虑燃耗作用的二氧化铀导热系数
        !FD=(1.09/(Bu(time,n_a,n_r)/9.383)**3.265+0.0643/((Bu(time,n_a,n_r)/9.383)**0.5)*(Temperature_transient(time1,n_a,n_r)**0.5)&
        ! &)*ATAN2(1.0,(1.09/((Bu(time,n_a,n_r)/9.383)**3.265)+0.0643/((Bu(time,n_a,n_r)/9.383)**0.5)*(Temperature_transient(time1,n_a,n_r)**0.5)))
        !!!!!!!!!!!!!!!!!FD-溶解因子
        !FP=1+(0.019*(Bu(time,n_a,n_r)/9.383)/(3.-0.019*(Bu(time,n_a,n_r)/9.383)))*(1./(1.+EXP(-(Temperature_transient(time1,n_a,n_r)-1200)/100.)))
        !!!!!!!!!!!!!!!!!FP-沉淀因子
        !FM=(1-P)/(1+(S-1)*P)!!!!!!!!!FM-孔隙因子
        !FR=1.-0.2/(1.+EXP((Temperature_transient(time1,n_a,n_r)-900.)/80.))!!!!!!!!!!!!!FR辐照因子
        !k_fuel(time1,n_a,n_r)=k0*FD*FM*FR               !!!!!!!!二氧化铀真实导热系数   
        k_fuel(:,:,:)=5678.
    end if
    end subroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!UO2芯块热导率计算!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!UO2芯块热导率计算!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!UO2芯块热导率计算!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    
end module hcore