SUBROUTINE MECHMODEL_unrigid(time_total,time_increment,Q_LINE,D0,DI,T_PRE,T_NOW,X_FC,BU_BEGIN,BU_END,NH,PRESS_INTER,PRESS_BEGIN,UR_PRA,P_contact,&
&N_CLAD,d_length_spring,d_length_p,d_length_c,yield_stress,yield_stress111,strain_plastic,strain_plastic111&
&,deng,strain_z_lasttime,strain_z_lasttime111,stress_equ,strain_fuel,strain_creep,strain_creep111&
&,stress_cladding_z,strain_cladding_z)


!&,yield_stress,yield_stress111,plastic_cladding_r,plastic_cladding_c,plastic_cladding_z,&
!&plastic_cladding_r111,plastic_cladding_c111,plastic_cladding_z111,&
!&strain_z_clad_s,strain_z_pellet_s,strain_z_clad_s111,strain_z_pellet_s111,p_contact,stress_equivalent_cladding) 
!d0,ȼ��Ԫ���⾶����ģ����δ���£���λm   (Ӧ�ô�����̬�ߴ�)
!di,�����ھ�����ģ����δ���£���λm      (Ӧ�ô�����̬�ߴ�)
!T_PRE,�¶����飨nh+2��,��ģ����δ���£���һȼ��ʱ�̵��¶�
!T_NOW,�¶����飨nh+2��,��ģ����δ���£���һȼ��ʱ�̵��¶�
!X_FC,��϶��ȣ���ģ����δ���£���λm   (Ӧ�ô�����̬�ߴ�)
!BU_BEGIN,ȼ�����飨nh+2������һʱ��ȼ�ģ���ģ�����и��£�����ȼ��Ϊƽ��ȼ�ģ���ģ���л����������·ֲ�
!BU_END,ȼ�����飨nh+2������ʱ��ȼ�ģ���ģ�����и��£�����ȼ��Ϊƽ��ȼ�ģ���ģ���л����������·ֲ�
!NH��о�黮�ֵľ���ڵ���Ŀ
!PRESS_INTER����һʱ��ȼ��Ԫ����ѹ��MPa��ģ�����޸���
!PRESS_BEGIN����һʱ��ȼ��Ԫ����ѹ��MPa��ģ�����޸���
!UR_PRA���뾶���飨nh+2����Ϊ���κ�İ뾶�����ڴ�����������ر�ע���ǰ뾶
!��111�Ķ�Ϊ��ģ����Ҫ���µı�����������111�Ķ�Ϊ��һʱ�̵ı�����ֵ���ǻ�׼��ģ���в�����



IMPLICIT none
integer,PARAMETER::NB=21 
integer::NH,i
real(8)::H
real(8)::D0,DI,X_FC
real(8)::T_PRE(nh+2),T_NOW(nh+2),BU_BEGIN(nh+2),BU_END(nh+2)
real(8)::PRESS_INTER,PRESS_BEGIN
real(8)::UR_PRA(nh+2)

integer::INCREMENT
real(8)::PIE,RFS,DH,TSINT,UO2_IDEALDENSITY,FUEL_percentage,UO2_DENSITY
real(8)::ER(NH+2),NUR(NH+2),EQSTR(NH+2),MO(NH+2) !Nh�㸳ֵо��ģ�δ�Ӵ�ʱ���ǽ���NH+1���Ӵ�������ʹ�ã�����NH��ֵ����
real(8)::l_densification_max,l_coordinate,BU_coordinate,BU_1,l_densification_begin
real(8)::DELT(NH+2)

real(8)::AJS
real(8)::BU_BEGINJS(NH),BU_ENDJS(NH),DELBU(NH)
real(8)::l_strain_pre(nh),l_strain(nh),l_Thermal_expansion_pre(nh),T_exchange(nh),l_Thermal_expansion(nh) !о��ڵ�
real(8)::l_densification_pre(nh),BU_exchange(nh),l_densification(nh)
real(8)::soldsw_pre(nh),soldsw(nh),swell_pre(nh),swell(nh)
real(8)::R(NH+2),DLTX   !ԭʼ�뾶
real(8)::p_contact,p_coolant=15.5    !mpa   
     



integer::N_CLAD   !���ǽڵ��� ������2��
real(8)::m_elastic(NH+N_CLAD),m_plastic(NH+N_CLAD),Poisson_ration(NH+N_CLAD)   !��е���Բ���
real(8)::d_length_spring,d_length_spring111
real(8)::d_length_p,d_length_c   !���㵯�ɵĳ��ȱ仯
real(8)::force_z_clad,force_z_pellet   !z�������ܵĺ���
real(8)::k_spring   !����ϵ����N/M
real(8)::r_mc(NH+N_CLAD)
real(8)::strain_intrinsic(NH+N_CLAD,3)  !����Ӧ����
real(8)::stress_equ(NH+N_CLAD),strain_equ(NH+N_CLAD)        !��ЧӦ��,��ЧӦ��
integer::factor     !�������Խ׶εĽڵ���




!о�鵥������ʱ������
real(8)::m_l(NH+N_CLAD-1,3,3),m_m(NH+N_CLAD-1,3)   !�������󣬹涨��i�ڵ㵽i+1�ڵ�Ĺ�ϵ
real(8)::m_A(NH+N_CLAD-1,3,3),m_B(NH+N_CLAD-1,3)   !ת�ƾ��󣬹涨��1�ڵ㵽i+1�ڵ�Ĺ�ϵ
real(8)::AA(3,3),BB(3,3),CC(3),DD(3),XX(3,3),YY(3)  !������þ��������ӳ���
real(8)::CCC(3,3),DDD(3)    !z�����Ϻ����ı߽�������ʾ
real(8)::AAA(3,3),BBB(3),XXX(3)    !��������һ���ڵ�������Ӧ�������ľ���
real(8)::stress_fuel(nh,3)  !о��Ӧ��
real(8)::strain_fuel(nh,3)  !о��Ӧ��

real(8)::yield_stress(NH+N_CLAD),yield_stress111(NH+N_CLAD)

real(8)::d_strain_plastic(NH+N_CLAD,3),d_strain_plastic_js(NH+N_CLAD,3)  
real(8)::d_strain_plastic_equ(NH+N_CLAD)
real(8)::stress_curve(NH+N_CLAD)
real(8)::s1,s2,s3
real(8)::factor2,factor3    !�����ж��Ƿ�������������
real(8)::stress_fuel_js(nh,3)
real(8)::strain_plastic(NH+N_CLAD,3),strain_plastic111(NH+N_CLAD,3)
real(8)::deng,force_total





real(8)::temperature_cladding
real(8)::strain_thermal_cladding_r,strain_thermal_cladding_c,strain_thermal_cladding_z
real(8)::stress_clad(N_CLAD,3),strain_clad(N_CLAD,3)
real(8)::stress_clad_js(N_CLAD,3)


real(8)::strain_z_lasttime(2)  !1Ϊо���ֵ��2Ϊ���ǵ�ֵ����ģ�鲻���£��������и��´��뱾ģ��
real(8)::strain_z_lasttime111(2)  !��ģ�����



real(8)::time_total,time_increment,Q_LINE,n_flux
real(8)::strain_creep(N_CLAD,3),strain_creep111(N_CLAD,3)

real(8)::creep_rate_th,creep_rate_irr
real(8)::creep_rate_add,creep_rate_total
real(8)::creep_saturated_primary,creep_total

real(8)::d_strain_creep(N_CLAD,3),d_strain_creep_js(N_CLAD,3)
integer::i2,i3,i4,i5,i6,i7
real(8)::factor44
real(8)::tem_cladding(N_CLAD)


real(8)::stress_cladding_z,strain_cladding_z


n_flux=1.*Q_LINE/1000./20.*9.2d17









!write(*,*)  D0,DI,x_fc,"������ݴ����Ƿ���ȷ"
!pause
!write(*,*) T_NOW
!pause

k_spring=200000.    !���ȡֵ�㣬����Ե���



do  i=1,NH+N_CLAD

  if (i<=NH) then
    m_elastic(i)=2.0*(10.**11.)
    m_plastic(i)=0.2*(10.**9.)          !о�����Ը�ֵ
    Poisson_ration(i)=0.25 
  else 
    m_elastic(i)=1.*(10.**11.)
    m_plastic(i)=0.5*(10.**9.)          !�������Ը�ֵ
    Poisson_ration(i)=0.316 
  end if

end do



INCREMENT=1
h=3.6576
p_coolant=15.5   !mpa


!о�鸳ֵ
PIE=3.1415926 
RFS=(DI-2.*X_FC)/2.             !о���⾶
DH=H/NB

TSINT=1300.+273.15              !�ս��¶�
UO2_IDEALDENSITY=10.96*10.**3   !��λ kg/m3
FUEL_percentage=0.95            !TD%
UO2_DENSITY=UO2_IDEALDENSITY*FUEL_percentage


 


     R(1)=0.                          !Initial radius,�ڳ����в��ܱ�
     R(NH+1)= DI/2.0                  !cladding radius 
    
     R(NH+2)= D0/2.0  
     DLTX=RFS/(NH-1) 
     DO I=2,NH
     R(I)=R(I-1)+DLTX 

     END DO
     !write(*,*) R
     !pause

     
     do i=1,NH+N_CLAD

       if (i<=NH+1) then

         r_mc(i)=R(i)

       else

         r_mc(i)=R(NH+1)+(i-NH-1)*(R(NH+2)-R(NH+1))/(N_CLAD-1)

       end if

          !�������帳ֵ������ֻ��һ��ʱ�����ֶ����Ҫ�޸�
     end do



     !write(*,*)  d_length_spring
     force_z_pellet=-1.*d_length_spring*k_spring
     force_z_clad=-force_z_pellet-(p_coolant-PRESS_INTER)*1.d6*PIE*R(NH+2)**2

     !write(*,*) (p_coolant-PRESS_INTER),PIE*R(NH+2)**2
     !write(*,*) force_z_pellet,force_z_clad
     !pause







     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!о��λ�Ƴ�������!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!о��λ�Ƴ�������!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!о��λ�Ƴ�������!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





!��ʵ�������ת��begin
l_densification_max=22.2*(100.-FUEL_percentage*100.)/(TSINT-1453.)
l_coordinate=-3.+l_densification_max   !�������ת��
BU_coordinate=0.005
BU_1=0.
do while(ABS(BU_coordinate-BU_1)/BU_coordinate>0.000002)
BU_1=BU_coordinate
BU_coordinate=BU_1-(l_coordinate+3.0-0.93*EXP(-BU_coordinate)-2.07*EXP(-35.*BU_coordinate))/(0.93*EXP(-BU_coordinate)+2.07*35.*EXP(-35.*BU_coordinate)) ! Newton's method
end do
l_densification_begin=-3.0+0.93*EXP(-BU_coordinate)+2.07*EXP(-35.*BU_coordinate)
!��ʵ�������ת��end

DO I=1,NH+2
DELT(I)=(T_NOW(I)-T_PRE(I))/INCREMENT
END DO 

! ȼ�ĵľ���ֲ� 
BU_END(1)=BU_END(9)*(1.-BU_END(9)/300.)
DO I=2,NH-3
BU_END(I)=BU_END(I-1)+(BU_END(9)-BU_END(1))/8.
END DO
AJS=0.
DO I=1,NH-2
AJS=AJS+BU_END(I)*pie*((0.1*I)**2-(0.1*(I-1))**2)
END DO
BU_END(10)=(BU_END(9)*PIE*(1.0**2.)-AJS)/(PIE*(1.0-0.9**2))
BU_END(NH)=BU_END(NH-1)

BU_BEGIN(1)=BU_BEGIN(9)*(1.-BU_BEGIN(9)/300.) 
DO I=2,NH-3
BU_BEGIN(I)=BU_BEGIN(I-1)+(BU_BEGIN(9)-BU_BEGIN(1))/8.
END DO
AJS=0.
DO I=1,NH-2
AJS=AJS+BU_BEGIN(I)*pie*((0.1*I)**2-(0.1*(I-1))**2)
END DO
BU_BEGIN(NH-1)=(BU_BEGIN(NH-2)*PIE*(1.0**2.)-AJS)/(PIE*(1.0-0.9**2))
BU_BEGIN(NH)=BU_BEGIN(NH-1)


DO I=1,NH
BU_BEGINJS(I)=BU_BEGIN(I)
BU_ENDJS(I)=BU_END(I)
END DO

DO I=2,NH-1
BU_BEGIN(I)=(BU_BEGINJS(I-1)+BU_BEGINJS(I))/2.
BU_END(I)=(BU_ENDJS(I-1)+BU_ENDJS(I))/2.
END DO
! ȼ�ĵľ���ֲ�

!write(*,*) BU_BEGIN
!write(*,*) BU_END

DO I=1,NH
DELBU(I)=(BU_END(I)-BU_BEGIN(I))/INCREMENT
END DO





     l_strain_pre=0.
     l_strain=0.


     DO I=1,NH
       l_Thermal_expansion_pre(I)=1.*10.**(-5.)*T_PRE(I)-3.0*10.**(-3.)     &
                                 +4.0*10.**(-2.)*exp(-6.9*10.**(-20.)/(1.38*10.**(-23.)*T_PRE(I)))
       T_exchange(I)=T_PRE(I)+DELT(I)
       l_Thermal_expansion(I)=1.*10.**(-5.)*T_exchange(I)-3.0*10.**(-3.)    &
	              +4.0*10.**(-2.)*exp(-6.9*10.**(-20.)/(1.38*10.**(-23.)*T_exchange(I)))

       IF(BU_BEGIN(I).LT.BU_coordinate) then
	      l_densification_pre(I)=0.
	   ELSE
	      l_densification_pre(I)=(-3.0+0.93*EXP(-BU_BEGIN(I))+2.07*EXP(-35.*BU_BEGIN(I))-l_densification_begin)*0.01
	   END IF	   	   
	   BU_exchange(I)=BU_BEGIN(I)+DELBU(I)
      
       IF(BU_exchange(I).LT.BU_coordinate) then
	      l_densification(I)=0.
	   ELSE
	      l_densification(I)=(-3.0+0.93*EXP(-BU_exchange(I))+2.07*EXP(-35.*BU_exchange(I))-l_densification_begin)*0.01
	   END IF
          

          soldsw_pre(i)=1.0*7.435*(10.**(-13.))*UO2_DENSITY*BU_BEGIN(I)*24.*3600.
          soldsw(I)=1.0*7.435*(10.**(-13.))*UO2_DENSITY*BU_exchange(I)*24.*3600.
          swell(I)=1.0/3.*soldsw(I)
          swell_pre(i)=1.0/3.*soldsw_pre(I)


          l_strain_pre(I)=l_Thermal_expansion_pre(I)+0.68*l_densification_pre(I)+swell_pre(i)
		  l_strain(I)=l_Thermal_expansion(I)+0.68*l_densification(I)+swell(I)
         
          !write(*,*) i,l_Thermal_expansion(I),l_densification(I),swell(I)
          !pause



     END DO
     !write(*,*) l_Thermal_expansion(1),l_densification(1),swell(1)
     !write(*,*) l_Thermal_expansion(11),l_densification(11),swell(11)
     !pause
     !!!!о����ε�����ɣ������ͬ�ԣ�

    

     d_strain_plastic=0.

  



     


      do i=1,NH

          strain_intrinsic(i,1)=l_strain(I)+strain_plastic(i,1)+d_strain_plastic(i,1)
          strain_intrinsic(i,2)=l_strain(I)+strain_plastic(i,2)+d_strain_plastic(i,2)
          strain_intrinsic(i,3)=l_strain(I)+strain_plastic(i,3)+d_strain_plastic(i,3)

      end do
     

 

     !!!!!!!!!!!!!!!!!!!!!�������ƾ������!!!!!!!!!!!!!!!!!!!!
     do i=1,NH-1

       m_l(i,1,1)=1.-(r_mc(i+1)-r_mc(i))/r_mc(i+1)
       m_l(i,1,2)=(r_mc(i+1)-r_mc(i))/r_mc(i+1)
       m_l(i,1,3)=0.
       m_m(i,1)=0.

       m_l(i,2,1)=(r_mc(i+1)-r_mc(i))/r_mc(i+1)+m_elastic(i+1)/(1.-Poisson_ration(i+1))&
       &*(Poisson_ration(i+1)/m_elastic(i+1)-Poisson_ration(i)/m_elastic(i))
       m_l(i,2,2)=1.-(r_mc(i+1)-r_mc(i))/r_mc(i+1)+m_elastic(i+1)/(1.-Poisson_ration(i+1)**2)&
       &*((Poisson_ration(i+1)*Poisson_ration(i+1)/m_elastic(i+1)-Poisson_ration(i+1)*Poisson_ration(i)/m_elastic(i))&
       &-(1./m_elastic(i+1)-1./m_elastic(i)))
       m_l(i,2,3)=m_elastic(i+1)/(1.-Poisson_ration(i+1)**2)*((Poisson_ration(i+1)/m_elastic(i+1)-Poisson_ration(i)/m_elastic(i))&
       &-(Poisson_ration(i+1)/m_elastic(i+1)-Poisson_ration(i+1)/m_elastic(i)))
       m_m(i,2)=m_elastic(i+1)/(1.-Poisson_ration(i+1)**2)*(((r_mc(i+1)-r_mc(i))/r_mc(i+1)*(strain_intrinsic(i+1,1)-strain_intrinsic(i+1,2)))&
       &-(strain_intrinsic(i+1,2)-strain_intrinsic(i,2))-Poisson_ration(i+1)*(strain_intrinsic(i+1,3)-strain_intrinsic(i,3)))

       m_l(i,3,1)=m_elastic(i+1)/(1.-Poisson_ration(i+1))*(Poisson_ration(i+1)/m_elastic(i+1)-Poisson_ration(i)/m_elastic(i))
       m_l(i,3,2)=m_elastic(i+1)/(1.-Poisson_ration(i+1)**2)*((Poisson_ration(i+1)/m_elastic(i+1)-Poisson_ration(i)/m_elastic(i))&
       &-(Poisson_ration(i+1)/m_elastic(i+1)-Poisson_ration(i+1)/m_elastic(i)))
       m_l(i,3,3)=1.+m_elastic(i+1)/(1.-Poisson_ration(i+1)**2)&
       &*((Poisson_ration(i+1)*Poisson_ration(i+1)/m_elastic(i+1)-Poisson_ration(i+1)*Poisson_ration(i)/m_elastic(i))&
       &-(1./m_elastic(i+1)-1./m_elastic(i)))
       m_m(i,3)=m_elastic(i+1)/(1.-Poisson_ration(i+1)**2)*(((r_mc(i+1)-r_mc(i))/r_mc(i+1)*&
       &Poisson_ration(i+1)*(strain_intrinsic(i+1,1)-strain_intrinsic(i+1,2)))&
       &-Poisson_ration(i+1)*(strain_intrinsic(i+1,2)-strain_intrinsic(i,2))-(strain_intrinsic(i+1,3)-strain_intrinsic(i,3)))

       !write(*,*) i
       !write(*,*) m_l(i,1,1),m_l(i,1,2),m_l(i,1,3)
       !write(*,*) m_l(i,2,1),m_l(i,2,2),m_l(i,2,3)
       !write(*,*) m_l(i,3,1),m_l(i,3,2),m_l(i,3,3)
       !write(*,*) m_m(i,1),m_m(i,2),m_m(i,3)
       !pause
     end do
     !!!!!!!!!!!!!!!!!!!!!�������ƾ������!!!!!!!!!!!!!!!!!!!!
   

     !!!!!!!!!!!!!!!!!!!!!ת�ƾ������!!!!!!!!!!!!!!!!!!!!
     do i=1,NH-1

       if (i==1) then

         m_A(i,1,1)=m_l(i,1,1)
         m_A(i,1,2)=m_l(i,1,2)
         m_A(i,1,3)=m_l(i,1,3)
         m_A(i,2,1)=m_l(i,2,1)
         m_A(i,2,2)=m_l(i,2,2)
         m_A(i,2,3)=m_l(i,2,3)
         m_A(i,3,1)=m_l(i,3,1)
         m_A(i,3,2)=m_l(i,3,2)
         m_A(i,3,3)=m_l(i,3,3)

         m_B(i,1)=m_m(i,1)
         m_B(i,2)=m_m(i,2)
         m_B(i,3)=m_m(i,3)

       else
         

         AA(1,1)=m_l(i,1,1)
         AA(1,2)=m_l(i,1,2)
         AA(1,3)=m_l(i,1,3)
         AA(2,1)=m_l(i,2,1)
         AA(2,2)=m_l(i,2,2)
         AA(2,3)=m_l(i,2,3)
         AA(3,1)=m_l(i,3,1)
         AA(3,2)=m_l(i,3,2)
         AA(3,3)=m_l(i,3,3)

         BB(1,1)=m_A(i-1,1,1)
         BB(1,2)=m_A(i-1,1,2)
         BB(1,3)=m_A(i-1,1,3)
         BB(2,1)=m_A(i-1,2,1)
         BB(2,2)=m_A(i-1,2,2)
         BB(2,3)=m_A(i-1,2,3)
         BB(3,1)=m_A(i-1,3,1)
         BB(3,2)=m_A(i-1,3,2)
         BB(3,3)=m_A(i-1,3,3)


         CC(1)=m_m(i,1)
         CC(2)=m_m(i,2)
         CC(3)=m_m(i,3)

         DD(1)=m_B(i-1,1)
         DD(2)=m_B(i-1,2)
         DD(3)=m_B(i-1,3)

         call matrix_multiply(AA,BB,XX)
         call series_multiply(AA,DD,YY)
         YY=YY+CC

         m_A(i,1,1)=XX(1,1)
         m_A(i,1,2)=XX(1,2)
         m_A(i,1,3)=XX(1,3)
         m_A(i,2,1)=XX(2,1)
         m_A(i,2,2)=XX(2,2)
         m_A(i,2,3)=XX(2,3)
         m_A(i,3,1)=XX(3,1)
         m_A(i,3,2)=XX(3,2)
         m_A(i,3,3)=XX(3,3)

         m_B(i,1)=YY(1)
         m_B(i,2)=YY(2)
         m_B(i,3)=YY(3)

       end if
      

       !write(*,*) i
       !write(*,*) m_l(i,1,1),m_l(i,1,2),m_l(i,1,3)
       !write(*,*) m_l(i,2,1),m_l(i,2,2),m_l(i,2,3)
       !write(*,*) m_l(i,3,1),m_l(i,3,2),m_l(i,3,3)
       !write(*,*) m_m(i,1),m_m(i,2),m_m(i,3)

       !write(*,*) m_A(i,1,1),m_A(i,1,2),m_l(i,1,3)
       !write(*,*) m_A(i,2,1),m_A(i,2,2),m_l(i,2,3)
       !write(*,*) m_A(i,3,1),m_A(i,3,2),m_l(i,3,3)
       !write(*,*) m_B(i,1),m_B(i,2),m_B(i,3)
       !pause
     


     end do
     !!!!!!!!!!!!!!!!!!!!!ת�ƾ������!!!!!!!!!!!!!!!!!!!!
     !PAUSE



     !!!!!!!!!!!!!!!!!!!!!z��������������!!!!!!!!!!!!!!!!!!!!


     CCC=0.
     DDD=0.
     do i=1,NH-1
       
       if (i==1) then
         
         ccc(1,1)=ccc(1,1)+0.5*pie*(r_mc(i+1)**2-r_mc(i)**2)*(1.+m_A(i,1,1))
         ccc(1,2)=ccc(1,2)+0.5*pie*(r_mc(i+1)**2-r_mc(i)**2)*(0.+m_A(i,1,2))
         ccc(1,3)=ccc(1,3)+0.5*pie*(r_mc(i+1)**2-r_mc(i)**2)*(0.+m_A(i,1,3))
         ccc(2,1)=ccc(2,1)+0.5*pie*(r_mc(i+1)**2-r_mc(i)**2)*(0.+m_A(i,2,1))
         ccc(2,2)=ccc(2,2)+0.5*pie*(r_mc(i+1)**2-r_mc(i)**2)*(1.+m_A(i,2,2))
         ccc(2,3)=ccc(2,3)+0.5*pie*(r_mc(i+1)**2-r_mc(i)**2)*(0.+m_A(i,2,3))
         ccc(3,1)=ccc(3,1)+0.5*pie*(r_mc(i+1)**2-r_mc(i)**2)*(0.+m_A(i,3,1))
         ccc(3,2)=ccc(3,2)+0.5*pie*(r_mc(i+1)**2-r_mc(i)**2)*(0.+m_A(i,3,2))
         ccc(3,3)=ccc(3,3)+0.5*pie*(r_mc(i+1)**2-r_mc(i)**2)*(1.+m_A(i,3,3))

         ddd(1)=ddd(1)+0.5*pie*(r_mc(i+1)**2-r_mc(i)**2)*m_B(i,1)
         ddd(2)=ddd(2)+0.5*pie*(r_mc(i+1)**2-r_mc(i)**2)*m_B(i,2)
         ddd(3)=ddd(3)+0.5*pie*(r_mc(i+1)**2-r_mc(i)**2)*m_B(i,3)

       else
         
         ccc(1,1)=ccc(1,1)+0.5*pie*(r_mc(i+1)**2-r_mc(i)**2)*(m_A(i-1,1,1)+m_A(i,1,1))
         ccc(1,2)=ccc(1,2)+0.5*pie*(r_mc(i+1)**2-r_mc(i)**2)*(m_A(i-1,1,2)+m_A(i,1,2))
         ccc(1,3)=ccc(1,3)+0.5*pie*(r_mc(i+1)**2-r_mc(i)**2)*(m_A(i-1,1,3)+m_A(i,1,3))
         ccc(2,1)=ccc(2,1)+0.5*pie*(r_mc(i+1)**2-r_mc(i)**2)*(m_A(i-1,2,1)+m_A(i,2,1))
         ccc(2,2)=ccc(2,2)+0.5*pie*(r_mc(i+1)**2-r_mc(i)**2)*(m_A(i-1,2,2)+m_A(i,2,2))
         ccc(2,3)=ccc(2,3)+0.5*pie*(r_mc(i+1)**2-r_mc(i)**2)*(m_A(i-1,2,3)+m_A(i,2,3))
         ccc(3,1)=ccc(3,1)+0.5*pie*(r_mc(i+1)**2-r_mc(i)**2)*(m_A(i-1,3,1)+m_A(i,3,1))
         ccc(3,2)=ccc(3,2)+0.5*pie*(r_mc(i+1)**2-r_mc(i)**2)*(m_A(i-1,3,2)+m_A(i,3,2))
         ccc(3,3)=ccc(3,3)+0.5*pie*(r_mc(i+1)**2-r_mc(i)**2)*(m_A(i-1,3,3)+m_A(i,3,3))

         ddd(1)=ddd(1)+0.5*pie*(r_mc(i+1)**2-r_mc(i)**2)*(m_B(i-1,1)+m_B(i,1))
         ddd(2)=ddd(2)+0.5*pie*(r_mc(i+1)**2-r_mc(i)**2)*(m_B(i-1,2)+m_B(i,2))
         ddd(3)=ddd(3)+0.5*pie*(r_mc(i+1)**2-r_mc(i)**2)*(m_B(i-1,3)+m_B(i,3))

       end if


     end do


     !!!!!!!!!!!!!!!!!!!!!z��������������!!!!!!!!!!!!!!!!!!!!(�ⲿ���Թ�)

       !write(*,*) CCC(1,1),CCC(1,2),CCC(1,3)
       !write(*,*) CCC(2,1),CCC(2,2),CCC(2,3)
       !write(*,*) CCC(3,1),CCC(3,2),CCC(3,3)
       !write(*,*) DDD(1),DDD(2),DDD(3)
       !pause





     !!!!!!!!!!!!!!!!!!����һ���ڵ��ϵ�Ӧ������!!!!!!!!!!!!!!!!!!!!!

     AAA(1,1)=1.
     AAA(1,2)=-1.
     AAA(1,3)=0.
     AAA(2,1)=m_A(NH-1,1,1)
     AAA(2,2)=m_A(NH-1,1,2)
     AAA(2,3)=m_A(NH-1,1,3)
     AAA(3,1)=ccc(3,1)
     AAA(3,2)=ccc(3,2)
     AAA(3,3)=ccc(3,3)

     BBB(1)=0.
     BBB(2)=-m_B(NH-1,1)-PRESS_INTER*1.D6
     BBB(3)=-ddd(3)+force_z_pellet

     call agaus2(AAA,bbb,XXX)

     stress_fuel(1,1)=XXX(1)
     stress_fuel(1,2)=XXX(2)
     stress_fuel(1,3)=XXX(3)

     stress_equ(1)=(0.5*((stress_fuel(1,1)-stress_fuel(1,2))**2+(stress_fuel(1,1)-stress_fuel(1,3))**2+(stress_fuel(1,2)-stress_fuel(1,3))**2))**0.5
 
     !write(*,*)  XXX
     !write(*,*)  stress_equ(1)/1.d6
     !pause
     !!!!!!!!!!!!!!!!!!����һ���ڵ��ϵ�Ӧ������!!!!!!!!!!!!!!!!!!!!!
     


     !!!!!!!!!!!!!!!!!!��������ڵ��ϵ�Ӧ������!!!!!!!!!!!!!!!!!!!!!
    
     do i=2,nh

       stress_fuel(i,1)=m_A(i-1,1,1)*stress_fuel(1,1)+m_A(i-1,1,2)*stress_fuel(1,2)+m_A(i-1,1,3)*stress_fuel(1,3)+m_b(i-1,1)
       stress_fuel(i,2)=m_A(i-1,2,1)*stress_fuel(1,1)+m_A(i-1,2,2)*stress_fuel(1,2)+m_A(i-1,2,3)*stress_fuel(1,3)+m_b(i-1,2)
       stress_fuel(i,3)=m_A(i-1,3,1)*stress_fuel(1,1)+m_A(i-1,3,2)*stress_fuel(1,2)+m_A(i-1,3,3)*stress_fuel(1,3)+m_b(i-1,3)

       stress_equ(i)=(0.5*((stress_fuel(i,1)-stress_fuel(i,2))**2+(stress_fuel(i,1)-stress_fuel(i,3))**2+(stress_fuel(i,2)-stress_fuel(i,3))**2))**0.5

       !write(*,*) i
       !write(*,*) stress_fuel(i,1)
       !write(*,*) stress_fuel(i,2)
       !write(*,*) stress_fuel(i,3)
       !write(*,*) stress_equ(i)/1.d6 
       !pause
     end do

     !!!!!!!!!!!!!!!!!!��������ڵ��ϵ�Ӧ������!!!!!!!!!!!!!!!!!!!!!
     !pause

     !force_total=0
     !do i=1,nh-1
     !force_total=force_total+0.5*pie*(r_mc(i+1)**2-r_mc(i)**2)*(stress_fuel(i,3)+stress_fuel(i+1,3))
     !end do
     !write(*,*)  force_z_pellet,force_total
     !pause




     !!!!!!!!!!!!!!!!!����Ӧ����Ӧ��!!!!!!!!!!!!!!!!!

     do i=1,nh

        strain_fuel(i,1)=1./m_elastic(i)*(stress_fuel(i,1)-Poisson_ration(i)*(stress_fuel(i,2)+stress_fuel(i,3)))+strain_intrinsic(i,1)
        strain_fuel(i,2)=1./m_elastic(i)*(stress_fuel(i,2)-Poisson_ration(i)*(stress_fuel(i,1)+stress_fuel(i,3)))+strain_intrinsic(i,2)
        strain_fuel(i,3)=1./m_elastic(i)*(stress_fuel(i,3)-Poisson_ration(i)*(stress_fuel(i,2)+stress_fuel(i,1)))+strain_intrinsic(i,3)

     end do


     !!!!!!!!!!!!!!!!!����Ӧ����Ӧ��!!!!!!!!!!!!!!!!!






     !!!!!!!!!!!!!!!!!!�жϽ����������!!!!!!!!!!!!!!!!

     factor=0
     do i=1,nh

       if (stress_equ(i)>yield_stress(i)) then

       !write(*,*) i,"��������"
       !pause

         !d_strain_plastic(i,1)=0.01
         !d_strain_plastic(i,2)=-0.005
         !d_strain_plastic(i,3)=-0.005

         factor=factor+1

       else
         
         !d_strain_plastic(i,1)=0.
         !d_strain_plastic(i,2)=0.
         !d_strain_plastic(i,3)=0.

         factor=factor+0

       end if


     end do

     !!!!!!!!!!!!!!!!!!�жϽ����������!!!!!!!!!!!!!!!!






if (factor==0) then




       !exit  !û���κ�һ��������Ա��Σ�������do while ѭ��

       do i=1,nh
       yield_stress111(i)=yield_stress(i)
       strain_plastic111(i,1)=strain_plastic(i,1)
       strain_plastic111(i,2)=strain_plastic(i,2)
       strain_plastic111(i,3)=strain_plastic(i,3)
       end do



else 
       

        !write(*,*)  "��������ѭ��"


        !pause

     !!!!!!!!!!!!!!!!!!��ֵ���Ա���!!!!!!!!!!!!!!!!

     !factor=0
     do i=1,nh

       if (stress_equ(i)>yield_stress(i)) then

       !write(*,*) i,"��������"
       !pause

         d_strain_plastic(i,1)=-0.0001
         d_strain_plastic(i,2)=0.
         d_strain_plastic(i,3)=+0.0001

         !factor=factor+1

       else
         
         d_strain_plastic(i,1)=0.
         d_strain_plastic(i,2)=0.
         d_strain_plastic(i,3)=0.

         !factor=factor+0

       end if

       !write(*,*) i
       !write(*,*) d_strain_plastic(i,1),d_strain_plastic(i,2),d_strain_plastic(i,3)

     end do

     !!!!!!!!!!!!!!!!!!��ֵ���Ա���!!!!!!!!!!!!!!!!







      
       

     
     factor2=1.
     factor3=1.





do while ((factor3>0.01))    !do while�����Ե���

     d_strain_plastic_js=d_strain_plastic


     do i=1,nh
     stress_fuel_js(i,1)=stress_fuel(i,1)
     stress_fuel_js(i,2)=stress_fuel(i,2)
     stress_fuel_js(i,3)=stress_fuel(i,3)
     !write(*,*) stress_fuel_js(i,1),stress_fuel_js(i,2),stress_fuel_js(i,3)
     end do

       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !!!!!!!!!!!!!���µ��Լ���!!!!!!!!!!!!!!
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             do i=1,NH

             strain_intrinsic(i,1)=l_strain(I)+strain_plastic(i,1)+d_strain_plastic(i,1)
             strain_intrinsic(i,2)=l_strain(I)+strain_plastic(i,2)+d_strain_plastic(i,2)
             strain_intrinsic(i,3)=l_strain(I)+strain_plastic(i,3)+d_strain_plastic(i,3)

             end do
     

 

          !!!!!!!!!!!!!!!!!!!!!�������ƾ������!!!!!!!!!!!!!!!!!!!!
          do i=1,NH-1

          m_l(i,1,1)=1.-(r_mc(i+1)-r_mc(i))/r_mc(i+1)
          m_l(i,1,2)=(r_mc(i+1)-r_mc(i))/r_mc(i+1)
          m_l(i,1,3)=0.
          m_m(i,1)=0.

          m_l(i,2,1)=(r_mc(i+1)-r_mc(i))/r_mc(i+1)+m_elastic(i+1)/(1.-Poisson_ration(i+1))&
          &*(Poisson_ration(i+1)/m_elastic(i+1)-Poisson_ration(i)/m_elastic(i))
          m_l(i,2,2)=1.-(r_mc(i+1)-r_mc(i))/r_mc(i+1)+m_elastic(i+1)/(1.-Poisson_ration(i+1)**2)&
          &*((Poisson_ration(i+1)*Poisson_ration(i+1)/m_elastic(i+1)-Poisson_ration(i+1)*Poisson_ration(i)/m_elastic(i))&
          &-(1./m_elastic(i+1)-1./m_elastic(i)))
          m_l(i,2,3)=m_elastic(i+1)/(1.-Poisson_ration(i+1)**2)*((Poisson_ration(i+1)/m_elastic(i+1)-Poisson_ration(i)/m_elastic(i))&
          &-(Poisson_ration(i+1)/m_elastic(i+1)-Poisson_ration(i+1)/m_elastic(i)))
          m_m(i,2)=m_elastic(i+1)/(1.-Poisson_ration(i+1)**2)*(((r_mc(i+1)-r_mc(i))/r_mc(i+1)*(strain_intrinsic(i+1,1)-strain_intrinsic(i+1,2)))&
          &-(strain_intrinsic(i+1,2)-strain_intrinsic(i,2))-Poisson_ration(i+1)*(strain_intrinsic(i+1,3)-strain_intrinsic(i,3)))

          m_l(i,3,1)=m_elastic(i+1)/(1.-Poisson_ration(i+1))*(Poisson_ration(i+1)/m_elastic(i+1)-Poisson_ration(i)/m_elastic(i))
          m_l(i,3,2)=m_elastic(i+1)/(1.-Poisson_ration(i+1)**2)*((Poisson_ration(i+1)/m_elastic(i+1)-Poisson_ration(i)/m_elastic(i))&
          &-(Poisson_ration(i+1)/m_elastic(i+1)-Poisson_ration(i+1)/m_elastic(i)))
          m_l(i,3,3)=1.+m_elastic(i+1)/(1.-Poisson_ration(i+1)**2)&
          &*((Poisson_ration(i+1)*Poisson_ration(i+1)/m_elastic(i+1)-Poisson_ration(i+1)*Poisson_ration(i)/m_elastic(i))&
          &-(1./m_elastic(i+1)-1./m_elastic(i)))
          m_m(i,3)=m_elastic(i+1)/(1.-Poisson_ration(i+1)**2)*(((r_mc(i+1)-r_mc(i))/r_mc(i+1)*&
          &Poisson_ration(i+1)*(strain_intrinsic(i+1,1)-strain_intrinsic(i+1,2)))&
          &-Poisson_ration(i+1)*(strain_intrinsic(i+1,2)-strain_intrinsic(i,2))-(strain_intrinsic(i+1,3)-strain_intrinsic(i,3)))

          !write(*,*) i
          !write(*,*) m_l(i,1,1),m_l(i,1,2),m_l(i,1,3)
          !write(*,*) m_l(i,2,1),m_l(i,2,2),m_l(i,2,3)
          !write(*,*) m_l(i,3,1),m_l(i,3,2),m_l(i,3,3)
          !write(*,*) m_m(i,1),m_m(i,2),m_m(i,3)
          !pause
        end do
        !!!!!!!!!!!!!!!!!!!!!�������ƾ������!!!!!!!!!!!!!!!!!!!!
   

        !!!!!!!!!!!!!!!!!!!!!ת�ƾ������!!!!!!!!!!!!!!!!!!!!
        do i=1,NH-1

        if (i==1) then

         m_A(i,1,1)=m_l(i,1,1)
         m_A(i,1,2)=m_l(i,1,2)
         m_A(i,1,3)=m_l(i,1,3)
         m_A(i,2,1)=m_l(i,2,1)
         m_A(i,2,2)=m_l(i,2,2)
         m_A(i,2,3)=m_l(i,2,3)
         m_A(i,3,1)=m_l(i,3,1)
         m_A(i,3,2)=m_l(i,3,2)
         m_A(i,3,3)=m_l(i,3,3)

         m_B(i,1)=m_m(i,1)
         m_B(i,2)=m_m(i,2)
         m_B(i,3)=m_m(i,3)

        else
         

         AA(1,1)=m_l(i,1,1)
         AA(1,2)=m_l(i,1,2)
         AA(1,3)=m_l(i,1,3)
         AA(2,1)=m_l(i,2,1)
         AA(2,2)=m_l(i,2,2)
         AA(2,3)=m_l(i,2,3)
         AA(3,1)=m_l(i,3,1)
         AA(3,2)=m_l(i,3,2)
         AA(3,3)=m_l(i,3,3)

         BB(1,1)=m_A(i-1,1,1)
         BB(1,2)=m_A(i-1,1,2)
         BB(1,3)=m_A(i-1,1,3)
         BB(2,1)=m_A(i-1,2,1)
         BB(2,2)=m_A(i-1,2,2)
         BB(2,3)=m_A(i-1,2,3)
         BB(3,1)=m_A(i-1,3,1)
         BB(3,2)=m_A(i-1,3,2)
         BB(3,3)=m_A(i-1,3,3)


         CC(1)=m_m(i,1)
         CC(2)=m_m(i,2)
         CC(3)=m_m(i,3)

         DD(1)=m_B(i-1,1)
         DD(2)=m_B(i-1,2)
         DD(3)=m_B(i-1,3)

         call matrix_multiply(AA,BB,XX)
         call series_multiply(AA,DD,YY)
         YY=YY+CC

         m_A(i,1,1)=XX(1,1)
         m_A(i,1,2)=XX(1,2)
         m_A(i,1,3)=XX(1,3)
         m_A(i,2,1)=XX(2,1)
         m_A(i,2,2)=XX(2,2)
         m_A(i,2,3)=XX(2,3)
         m_A(i,3,1)=XX(3,1)
         m_A(i,3,2)=XX(3,2)
         m_A(i,3,3)=XX(3,3)

         m_B(i,1)=YY(1)
         m_B(i,2)=YY(2)
         m_B(i,3)=YY(3)

        end if

        !write(*,*) i
        !write(*,*) m_l(i,1,1),m_l(i,1,2),m_l(i,1,3)
        !write(*,*) m_l(i,2,1),m_l(i,2,2),m_l(i,2,3)
        !write(*,*) m_l(i,3,1),m_l(i,3,2),m_l(i,3,3)
        !write(*,*) m_m(i,1),m_m(i,2),m_m(i,3)

        !write(*,*) m_A(i,1,1),m_A(i,1,2),m_l(i,1,3)
        !write(*,*) m_A(i,2,1),m_A(i,2,2),m_l(i,2,3)
        !write(*,*) m_A(i,3,1),m_A(i,3,2),m_l(i,3,3)
        !write(*,*) m_B(i,1),m_B(i,2),m_B(i,3)
        !pause
     


      end do
      !!!!!!!!!!!!!!!!!!!!!ת�ƾ������!!!!!!!!!!!!!!!!!!!!
      !PAUSE



      !!!!!!!!!!!!!!!!!!!!!z��������������!!!!!!!!!!!!!!!!!!!!


      CCC=0.
      DDD=0.
      do i=1,NH-1
       
        if (i==1) then
         
         ccc(1,1)=ccc(1,1)+0.5*pie*(r_mc(i+1)**2-r_mc(i)**2)*(1.+m_A(i,1,1))
         ccc(1,2)=ccc(1,2)+0.5*pie*(r_mc(i+1)**2-r_mc(i)**2)*(0.+m_A(i,1,2))
         ccc(1,3)=ccc(1,3)+0.5*pie*(r_mc(i+1)**2-r_mc(i)**2)*(0.+m_A(i,1,3))
         ccc(2,1)=ccc(2,1)+0.5*pie*(r_mc(i+1)**2-r_mc(i)**2)*(0.+m_A(i,2,1))
         ccc(2,2)=ccc(2,2)+0.5*pie*(r_mc(i+1)**2-r_mc(i)**2)*(1.+m_A(i,2,2))
         ccc(2,3)=ccc(2,3)+0.5*pie*(r_mc(i+1)**2-r_mc(i)**2)*(0.+m_A(i,2,3))
         ccc(3,1)=ccc(3,1)+0.5*pie*(r_mc(i+1)**2-r_mc(i)**2)*(0.+m_A(i,3,1))
         ccc(3,2)=ccc(3,2)+0.5*pie*(r_mc(i+1)**2-r_mc(i)**2)*(0.+m_A(i,3,2))
         ccc(3,3)=ccc(3,3)+0.5*pie*(r_mc(i+1)**2-r_mc(i)**2)*(1.+m_A(i,3,3))

         ddd(1)=ddd(1)+0.5*pie*(r_mc(i+1)**2-r_mc(i)**2)*m_B(i,1)
         ddd(2)=ddd(2)+0.5*pie*(r_mc(i+1)**2-r_mc(i)**2)*m_B(i,2)
         ddd(3)=ddd(3)+0.5*pie*(r_mc(i+1)**2-r_mc(i)**2)*m_B(i,3)

        else
         
         ccc(1,1)=ccc(1,1)+0.5*pie*(r_mc(i+1)**2-r_mc(i)**2)*(m_A(i-1,1,1)+m_A(i,1,1))
         ccc(1,2)=ccc(1,2)+0.5*pie*(r_mc(i+1)**2-r_mc(i)**2)*(m_A(i-1,1,2)+m_A(i,1,2))
         ccc(1,3)=ccc(1,3)+0.5*pie*(r_mc(i+1)**2-r_mc(i)**2)*(m_A(i-1,1,3)+m_A(i,1,3))
         ccc(2,1)=ccc(2,1)+0.5*pie*(r_mc(i+1)**2-r_mc(i)**2)*(m_A(i-1,2,1)+m_A(i,2,1))
         ccc(2,2)=ccc(2,2)+0.5*pie*(r_mc(i+1)**2-r_mc(i)**2)*(m_A(i-1,2,2)+m_A(i,2,2))
         ccc(2,3)=ccc(2,3)+0.5*pie*(r_mc(i+1)**2-r_mc(i)**2)*(m_A(i-1,2,3)+m_A(i,2,3))
         ccc(3,1)=ccc(3,1)+0.5*pie*(r_mc(i+1)**2-r_mc(i)**2)*(m_A(i-1,3,1)+m_A(i,3,1))
         ccc(3,2)=ccc(3,2)+0.5*pie*(r_mc(i+1)**2-r_mc(i)**2)*(m_A(i-1,3,2)+m_A(i,3,2))
         ccc(3,3)=ccc(3,3)+0.5*pie*(r_mc(i+1)**2-r_mc(i)**2)*(m_A(i-1,3,3)+m_A(i,3,3))

         ddd(1)=ddd(1)+0.5*pie*(r_mc(i+1)**2-r_mc(i)**2)*(m_B(i-1,1)+m_B(i,1))
         ddd(2)=ddd(2)+0.5*pie*(r_mc(i+1)**2-r_mc(i)**2)*(m_B(i-1,2)+m_B(i,2))
         ddd(3)=ddd(3)+0.5*pie*(r_mc(i+1)**2-r_mc(i)**2)*(m_B(i-1,3)+m_B(i,3))

        end if


      end do


      !!!!!!!!!!!!!!!!!!!!!z��������������!!!!!!!!!!!!!!!!!!!!(�ⲿ���Թ�)

       !write(*,*) CCC(1,1),CCC(1,2),CCC(1,3)
       !write(*,*) CCC(2,1),CCC(2,2),CCC(2,3)
       !write(*,*) CCC(3,1),CCC(3,2),CCC(3,3)
       !write(*,*) DDD(1),DDD(2),DDD(3)
       !pause





      !!!!!!!!!!!!!!!!!!����һ���ڵ��ϵ�Ӧ������!!!!!!!!!!!!!!!!!!!!!

      AAA(1,1)=1.
      AAA(1,2)=-1.
      AAA(1,3)=0.
      AAA(2,1)=m_A(NH-1,1,1)
      AAA(2,2)=m_A(NH-1,1,2)
      AAA(2,3)=m_A(NH-1,1,3)
      AAA(3,1)=ccc(3,1)
      AAA(3,2)=ccc(3,2)
      AAA(3,3)=ccc(3,3)

      BBB(1)=0.
      BBB(2)=-m_B(NH-1,1)-PRESS_INTER*1.D6
      BBB(3)=-ddd(3)+force_z_pellet

      call agaus2(AAA,bbb,XXX)

      stress_fuel(1,1)=XXX(1)
      stress_fuel(1,2)=XXX(2)
      stress_fuel(1,3)=XXX(3)

      stress_equ(1)=(0.5*((stress_fuel(1,1)-stress_fuel(1,2))**2+(stress_fuel(1,1)-stress_fuel(1,3))**2+(stress_fuel(1,2)-stress_fuel(1,3))**2))**0.5

      !write(*,*)  XXX
      !write(*,*)  stress_equ(1)/1.d6
      !pause
      !!!!!!!!!!!!!!!!!!����һ���ڵ��ϵ�Ӧ������!!!!!!!!!!!!!!!!!!!!!
     


      !!!!!!!!!!!!!!!!!!��������ڵ��ϵ�Ӧ������!!!!!!!!!!!!!!!!!!!!!
    
      do i=2,nh

       stress_fuel(i,1)=m_A(i-1,1,1)*stress_fuel(1,1)+m_A(i-1,1,2)*stress_fuel(1,2)+m_A(i-1,1,3)*stress_fuel(1,3)+m_b(i-1,1)
       stress_fuel(i,2)=m_A(i-1,2,1)*stress_fuel(1,1)+m_A(i-1,2,2)*stress_fuel(1,2)+m_A(i-1,2,3)*stress_fuel(1,3)+m_b(i-1,2)
       stress_fuel(i,3)=m_A(i-1,3,1)*stress_fuel(1,1)+m_A(i-1,3,2)*stress_fuel(1,2)+m_A(i-1,3,3)*stress_fuel(1,3)+m_b(i-1,3)

       stress_equ(i)=(0.5*((stress_fuel(i,1)-stress_fuel(i,2))**2+(stress_fuel(i,1)-stress_fuel(i,3))**2+(stress_fuel(i,2)-stress_fuel(i,3))**2))**0.5

       !write(*,*) i
       !write(*,*) stress_fuel(i,1)
       !write(*,*) stress_fuel(i,2)
       !write(*,*) stress_fuel(i,3)
       !write(*,*) stress_equ(i)/1.d6 
       !pause
      end do

      !!!!!!!!!!!!!!!!!!��������ڵ��ϵ�Ӧ������!!!!!!!!!!!!!!!!!!!!!
      !pause


      !!!!!!!!!!!!!!!!!����Ӧ����Ӧ��!!!!!!!!!!!!!!!!!

      do i=1,nh

        strain_fuel(i,1)=1./m_elastic(i)*(stress_fuel(i,1)-Poisson_ration(i)*(stress_fuel(i,2)+stress_fuel(i,3)))+strain_intrinsic(i,1)
        strain_fuel(i,2)=1./m_elastic(i)*(stress_fuel(i,2)-Poisson_ration(i)*(stress_fuel(i,1)+stress_fuel(i,3)))+strain_intrinsic(i,2)
        strain_fuel(i,3)=1./m_elastic(i)*(stress_fuel(i,3)-Poisson_ration(i)*(stress_fuel(i,2)+stress_fuel(i,1)))+strain_intrinsic(i,3)

      end do


      !!!!!!!!!!!!!!!!����Ӧ����Ӧ��!!!!!!!!!!!!!!!!!



       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !!!!!!!!!!!!!���µ��Լ���!!!!!!!!!!!!!!
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



     !!!!!!!!!!!!!!!!!!���¸�ֵ���Ա���!!!!!!!!!!!!!!!!

     !factor=0
     do i=1,nh

       if ((stress_equ(i)>yield_stress(i)).and.(abs(d_strain_plastic(i,1))+abs(d_strain_plastic(i,2))+abs(d_strain_plastic(i,3))==0.)) then


           d_strain_plastic(i,1)=-0.0001
           d_strain_plastic(i,2)=0.
           d_strain_plastic(i,3)=+0.0001

           !write(*,*) i,"���ĸ�ֵ"
           !write(*,*) d_strain_plastic(i,1),d_strain_plastic(i,2),d_strain_plastic(i,3)

       end if





     end do

     !!!!!!!!!!!!!!!!!!���¸�ֵ���Ա���!!!!!!!!!!!!!!!!







  



       factor2=0.
       factor3=0.
       !!!!!!!!!!!!!!!!!!���¼������Է���!!!!!!!!!!!!!!!!!!!
       do i=1,nh

         if (stress_equ(i)>yield_stress(i)) then
         !write(*,*)  i,yield_stress(i),stress_equ(i)/1.d6

           d_strain_plastic_equ(i)=2.**0.5/3.*((d_strain_plastic(i,1)-d_strain_plastic(i,2))**2+&
           &(d_strain_plastic(i,1)-d_strain_plastic(i,3))**2+(d_strain_plastic(i,2)-d_strain_plastic(i,2))**2)**0.5

           stress_curve(i)=yield_stress(i)+d_strain_plastic_equ(i)*m_plastic(i)

           s1=stress_fuel(i,1)-(stress_fuel(i,1)+stress_fuel(i,2)+stress_fuel(i,3))/3.  !pa
           s2=stress_fuel(i,2)-(stress_fuel(i,1)+stress_fuel(i,2)+stress_fuel(i,3))/3.  !pa
           s3=stress_fuel(i,3)-(stress_fuel(i,1)+stress_fuel(i,2)+stress_fuel(i,3))/3.  !pa


           d_strain_plastic(i,1)=3./2.*d_strain_plastic_equ(i)/stress_curve(i)*s1
           d_strain_plastic(i,2)=3./2.*d_strain_plastic_equ(i)/stress_curve(i)*s2
           d_strain_plastic(i,3)=3./2.*d_strain_plastic_equ(i)/stress_curve(i)*s3

           factor2=factor2+abs((d_strain_plastic(i,1)-d_strain_plastic_js(i,1))/d_strain_plastic(i,1))&
           &+abs((d_strain_plastic(i,2)-d_strain_plastic_js(i,2))/d_strain_plastic(i,2))&                    
           &+abs((d_strain_plastic(i,3)-d_strain_plastic_js(i,3))/d_strain_plastic(i,3))

           !write(*,*) d_strain_plastic(i,1),d_strain_plastic(i,2),d_strain_plastic(i,3)
           !write(*,*) d_strain_plastic_js(i,1),d_strain_plastic_js(i,2),d_strain_plastic_js(i,3)


           yield_stress111(i)=stress_equ(i)
           strain_plastic111(i,1)=strain_plastic(i,1)+d_strain_plastic(i,1)
           strain_plastic111(i,2)=strain_plastic(i,2)+d_strain_plastic(i,2)
           strain_plastic111(i,3)=strain_plastic(i,3)+d_strain_plastic(i,3)

           !write(*,*) i,stress_curve(i)/1.d6,stress_equ(i)/1.d6
           !write(*,*) d_strain_plastic(i,1),d_strain_plastic(i,2),d_strain_plastic(i,3)
           
 
         else
  

           d_strain_plastic(i,1)=0.5*d_strain_plastic(i,1)
           d_strain_plastic(i,2)=0.5*d_strain_plastic(i,2)      !���Ϊ�㣬����0.5����ν
           d_strain_plastic(i,3)=0.5*d_strain_plastic(i,3)      !�����Ϊ�㣬��˵����ʼ�ĳ�ֵ����ˣ���С

           yield_stress111(i)=yield_stress(i)
           strain_plastic111(i,1)=strain_plastic(i,1)
           strain_plastic111(i,2)=strain_plastic(i,2)
           strain_plastic111(i,3)=strain_plastic(i,3)

           factor2=factor2+0.

           !write(*,*) i,stress_curve(i)/1.d6,stress_equ(i)/1.d6
           !write(*,*) d_strain_plastic(i,1),d_strain_plastic(i,2),d_strain_plastic(i,3)
 
         end if

        

         factor3=factor3+abs((stress_fuel(i,1)-stress_fuel_js(i,1))/stress_fuel(i,1))&
           &+abs((stress_fuel(i,2)-stress_fuel_js(i,2))/stress_fuel(i,2))&
           &+abs((stress_fuel(i,3)-stress_fuel_js(i,3))/stress_fuel(i,3))
        
        !write(*,*) stress_fuel(i,1),stress_fuel(i,2),stress_fuel(i,3)
        !write(*,*) stress_fuel_js(i,1),stress_fuel_js(i,2),stress_fuel_js(i,3)
        !write(*,*) i,factor2,factor3
        !write(*,*) d_strain_plastic(i,1),d_strain_plastic(i,2),d_strain_plastic(i,3)
        !pause

       end do
       !!!!!!!!!!!!!!!!!!���¼������Է���!!!!!!!!!!!!!!!!!!!


        !write(*,*) factor2,factor3
        !write(*,*) 
        !pause

    end do !do while�����Ե���




end if    !�нڵ��������


deng=strain_fuel(11,2)




!write(*,*)  "һ������������"
!pause


     UR_PRA(1)=R(1)
     DO i=2,nh
     UR_PRA(i)=(1.+strain_fuel(i,2))*R(i)           !�ڶ�����������
     !UR_PRA(i)=UR_PRA(i-1)+DLTX*(1.+0.5(l_strain(I)+l_strain(I-1)))
     !UR_PRA(i)=UR_PRA(i-1)+DLTX+DLTX*(l_strain(I)+l_strain(I-1))/2.        !ԭ�����������о��λ��
     !write(*,*) i
     !write(*,*) UR_PRA(i),R(I),UR_PRA(i)-R(I)
     end do
     !write(*,*) UR_PRA(11),R(11),UR_PRA(11)-R(11)
     

     d_length_p=strain_fuel(11,3)*DH                          !��ʵ�κ�һ�������������Ӧ����һ����
     
     
     !write(*,*) 

     !pause

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!о��λ�Ƴ�������!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!о��λ�Ƴ�������!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!о��λ�Ƴ�������!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


     







!integer::N_CLAD   !���ǽڵ��� ������2��
!real(8)::m_elastic(NH+N_CLAD),m_plastic(NH+N_CLAD),Poisson_ration(NH+N_CLAD)   !��е���Բ���
!real(8)::r_mc(NH+N_CLAD)
!real(8)::strain_intrinsic(NH+N_CLAD,3)  !����Ӧ����
!real(8)::stress_equ(NH+N_CLAD),strain_equ(NH+N_CLAD)        !��ЧӦ��,��ЧӦ��
!real(8)::yield_stress(NH+N_CLAD),yield_stress111(NH+N_CLAD)
!real(8)::d_strain_plastic(NH+N_CLAD,3),d_strain_plastic_js(NH+N_CLAD,3)  
!real(8)::d_strain_plastic_equ(NH+N_CLAD)
!real(8)::stress_curve(NH+N_CLAD)
!real(8)::strain_plastic(NH+N_CLAD,3),strain_plastic111(NH+N_CLAD,3)




!о�鵥������ʱ������
!real(8)::m_l(NH+N_CLAD-1,3,3),m_m(NH+N_CLAD-1,3)   !�������󣬹涨��i�ڵ㵽i+1�ڵ�Ĺ�ϵ
!real(8)::m_A(NH+N_CLAD-1,3,3),m_B(NH+N_CLAD-1,3)   !ת�ƾ��󣬹涨��1�ڵ㵽i+1�ڵ�Ĺ�ϵ
!real(8)::AA(3,3),BB(3,3),CC(3),DD(3),XX(3,3),YY(3)  !������þ��������ӳ���
!real(8)::CCC(3,3),DDD(3)    !z�����Ϻ����ı߽�������ʾ
!real(8)::AAA(3,3),BBB(3),XXX(3)    !��������һ���ڵ�������Ӧ�������ľ���
!real(8)::stress_fuel(nh,3)  !о��Ӧ��
!real(8)::strain_fuel(nh,3)  !о��Ӧ��
!real(8)::stress_fuel_js(nh,3)



!real(8)::time_total,time_increment,Q_LINE,n_flux
!real(8)::strain_creep(N_CLAD,3),strain_creep111(N_CLAD,3)

!real(8)::creep_rate_th,creep_rate_irr
!real(8)::creep_rate_add,creep_rate_total
!real(8)::creep_saturated_primary,creep_total

!real(8)::d_strain_creep(N_CLAD,3),d_strain_creep_js(N_CLAD,3)
!integer:i2,i3,i4,i5,i6,i7
!real(8)::factor44









do i2=1,n_clad
d_strain_creep(i2,1)=1.d-4
d_strain_creep(i2,2)=5.d-5
d_strain_creep(i2,3)=5.d-5
d_strain_creep_js(i2,1)=0
d_strain_creep_js(i2,2)=0
d_strain_creep_js(i2,3)=0
end do

factor44=maxval(abs((d_strain_creep-d_strain_creep_js)/d_strain_creep))



do while (factor44>0.02)



!WRITE(*,*)  factor44
!PAUSE

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!����λ�Ƴ�������!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!����λ�Ƴ�������!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!����λ�Ƴ�������!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    
do i2=1,n_clad
d_strain_creep_js(i2,1)=d_strain_creep(i2,1)
d_strain_creep_js(i2,2)=d_strain_creep(i2,2)
d_strain_creep_js(i2,3)=d_strain_creep(i2,3)
end do


    

     !d_strain_plastic=0.     !��Ҫ���Σ���Ϊо����Ѿ�����ˣ�������������Ϊ���ˣ�Ҳ�����ٸ��㣬��Ϊǰ�渳ֵ����

  
  do i=1,N_CLAD
       
       tem_cladding(i)=T_NOW(nh+1)+(i-1)*(T_NOW(nh+2)-T_NOW(nh+1))/(N_CLAD-1)
     
       temperature_cladding=T_NOW(nh+1)+(i-1)*(T_NOW(nh+2)-T_NOW(nh+1))/(N_CLAD-1)
      
     if (temperature_cladding<1083.) then
       strain_thermal_cladding_r=1.*(4.95d-6*temperature_cladding-1.485d-3)
       strain_thermal_cladding_z=1.*(1.26d-5*temperature_cladding-3.78d-3)
       strain_thermal_cladding_c=strain_thermal_cladding_r
  
     else if (temperature_cladding<1244.) then
       strain_thermal_cladding_r=(2.77763+1.09822*cos(3.14159265*(temperature_cladding-1083.)/161.))/1000.
       strain_thermal_cladding_z=(2.77763+1.09822*cos(3.14159265*(temperature_cladding-1083.)/161.))/1000.
       strain_thermal_cladding_c=strain_thermal_cladding_r
     else if (temperature_cladding<2098.) then
       strain_thermal_cladding_r=9.7d-6*temperature_cladding-1.04d-2
       strain_thermal_cladding_z=9.76d-6*temperature_cladding-4.4d-3
       strain_thermal_cladding_c=strain_thermal_cladding_r
     else
       strain_thermal_cladding_r=9.7d-6*2098.-1.04d-2
       strain_thermal_cladding_z=9.76d-6*2098.-4.4d-3
       strain_thermal_cladding_c=strain_thermal_cladding_r
     end if

     !write(*,*) temperature_cladding
          !����Ӧ�������㣬ע���Ǵ����

          strain_intrinsic(NH+i,1)=strain_thermal_cladding_r+strain_plastic(NH+i,1)+d_strain_plastic(NH+i,1)+strain_creep(i,1)+d_strain_creep(i,1)
          strain_intrinsic(NH+i,2)=strain_thermal_cladding_c+strain_plastic(NH+i,2)+d_strain_plastic(NH+i,2)+strain_creep(i,2)+d_strain_creep(i,2)
          strain_intrinsic(NH+i,3)=strain_thermal_cladding_z+strain_plastic(NH+i,3)+d_strain_plastic(NH+i,3)+strain_creep(i,3)+d_strain_creep(i,3)
    
          !write(*,*) strain_intrinsic(NH+i,1),strain_intrinsic(NH+i,2),strain_intrinsic(NH+i,3)
          !write(*,*) strain_thermal_cladding_r,strain_thermal_cladding_c,strain_thermal_cladding_z
   end do
     

 !pause




     !!!!!!!!!!!!!!!!!!!!!�������ƾ������!!!!!!!!!!!!!!!!!!!!
     do i=NH+1,NH+N_CLAD-1           !о��Ͱ���δ�Ӵ�ʱ����NH��������Ҫ����Ϊо��Ͱ���֮���޹���

       m_l(i,1,1)=1.-(r_mc(i+1)-r_mc(i))/r_mc(i+1)
       m_l(i,1,2)=(r_mc(i+1)-r_mc(i))/r_mc(i+1)
       m_l(i,1,3)=0.
       m_m(i,1)=0.

       m_l(i,2,1)=(r_mc(i+1)-r_mc(i))/r_mc(i+1)+m_elastic(i+1)/(1.-Poisson_ration(i+1))&
       &*(Poisson_ration(i+1)/m_elastic(i+1)-Poisson_ration(i)/m_elastic(i))
       m_l(i,2,2)=1.-(r_mc(i+1)-r_mc(i))/r_mc(i+1)+m_elastic(i+1)/(1.-Poisson_ration(i+1)**2)&
       &*((Poisson_ration(i+1)*Poisson_ration(i+1)/m_elastic(i+1)-Poisson_ration(i+1)*Poisson_ration(i)/m_elastic(i))&
       &-(1./m_elastic(i+1)-1./m_elastic(i)))
       m_l(i,2,3)=m_elastic(i+1)/(1.-Poisson_ration(i+1)**2)*((Poisson_ration(i+1)/m_elastic(i+1)-Poisson_ration(i)/m_elastic(i))&
       &-(Poisson_ration(i+1)/m_elastic(i+1)-Poisson_ration(i+1)/m_elastic(i)))
       m_m(i,2)=m_elastic(i+1)/(1.-Poisson_ration(i+1)**2)*(((r_mc(i+1)-r_mc(i))/r_mc(i+1)*(strain_intrinsic(i+1,1)-strain_intrinsic(i+1,2)))&
       &-(strain_intrinsic(i+1,2)-strain_intrinsic(i,2))-Poisson_ration(i+1)*(strain_intrinsic(i+1,3)-strain_intrinsic(i,3)))

       m_l(i,3,1)=m_elastic(i+1)/(1.-Poisson_ration(i+1))*(Poisson_ration(i+1)/m_elastic(i+1)-Poisson_ration(i)/m_elastic(i))
       m_l(i,3,2)=m_elastic(i+1)/(1.-Poisson_ration(i+1)**2)*((Poisson_ration(i+1)/m_elastic(i+1)-Poisson_ration(i)/m_elastic(i))&
       &-(Poisson_ration(i+1)/m_elastic(i+1)-Poisson_ration(i+1)/m_elastic(i)))
       m_l(i,3,3)=1.+m_elastic(i+1)/(1.-Poisson_ration(i+1)**2)&
       &*((Poisson_ration(i+1)*Poisson_ration(i+1)/m_elastic(i+1)-Poisson_ration(i+1)*Poisson_ration(i)/m_elastic(i))&
       &-(1./m_elastic(i+1)-1./m_elastic(i)))
       m_m(i,3)=m_elastic(i+1)/(1.-Poisson_ration(i+1)**2)*(((r_mc(i+1)-r_mc(i))/r_mc(i+1)*&
       &Poisson_ration(i+1)*(strain_intrinsic(i+1,1)-strain_intrinsic(i+1,2)))&
       &-Poisson_ration(i+1)*(strain_intrinsic(i+1,2)-strain_intrinsic(i,2))-(strain_intrinsic(i+1,3)-strain_intrinsic(i,3)))

       !write(*,*) i
       !write(*,*) m_l(i,1,1),m_l(i,1,2),m_l(i,1,3)
       !write(*,*) m_l(i,2,1),m_l(i,2,2),m_l(i,2,3)
       !write(*,*) m_l(i,3,1),m_l(i,3,2),m_l(i,3,3)
       !write(*,*) m_m(i,1),m_m(i,2),m_m(i,3)
       !pause
     end do
     !!!!!!!!!!!!!!!!!!!!!�������ƾ������!!!!!!!!!!!!!!!!!!!!
   
   !pause



     !!!!!!!!!!!!!!!!!!!!!ת�ƾ������!!!!!!!!!!!!!!!!!!!!
     do i=NH+1,NH+N_CLAD-1

       if (i==NH+1) then   !����ط���Ҫ��

         m_A(i,1,1)=m_l(i,1,1)
         m_A(i,1,2)=m_l(i,1,2)
         m_A(i,1,3)=m_l(i,1,3)
         m_A(i,2,1)=m_l(i,2,1)
         m_A(i,2,2)=m_l(i,2,2)
         m_A(i,2,3)=m_l(i,2,3)
         m_A(i,3,1)=m_l(i,3,1)
         m_A(i,3,2)=m_l(i,3,2)
         m_A(i,3,3)=m_l(i,3,3)

         m_B(i,1)=m_m(i,1)
         m_B(i,2)=m_m(i,2)
         m_B(i,3)=m_m(i,3)

       else
         

         AA(1,1)=m_l(i,1,1)
         AA(1,2)=m_l(i,1,2)
         AA(1,3)=m_l(i,1,3)
         AA(2,1)=m_l(i,2,1)
         AA(2,2)=m_l(i,2,2)
         AA(2,3)=m_l(i,2,3)
         AA(3,1)=m_l(i,3,1)
         AA(3,2)=m_l(i,3,2)
         AA(3,3)=m_l(i,3,3)

         BB(1,1)=m_A(i-1,1,1)
         BB(1,2)=m_A(i-1,1,2)
         BB(1,3)=m_A(i-1,1,3)
         BB(2,1)=m_A(i-1,2,1)
         BB(2,2)=m_A(i-1,2,2)
         BB(2,3)=m_A(i-1,2,3)
         BB(3,1)=m_A(i-1,3,1)
         BB(3,2)=m_A(i-1,3,2)
         BB(3,3)=m_A(i-1,3,3)


         CC(1)=m_m(i,1)
         CC(2)=m_m(i,2)
         CC(3)=m_m(i,3)

         DD(1)=m_B(i-1,1)
         DD(2)=m_B(i-1,2)
         DD(3)=m_B(i-1,3)

         call matrix_multiply(AA,BB,XX)
         call series_multiply(AA,DD,YY)
         YY=YY+CC

         m_A(i,1,1)=XX(1,1)
         m_A(i,1,2)=XX(1,2)
         m_A(i,1,3)=XX(1,3)
         m_A(i,2,1)=XX(2,1)
         m_A(i,2,2)=XX(2,2)
         m_A(i,2,3)=XX(2,3)
         m_A(i,3,1)=XX(3,1)
         m_A(i,3,2)=XX(3,2)
         m_A(i,3,3)=XX(3,3)

         m_B(i,1)=YY(1)
         m_B(i,2)=YY(2)
         m_B(i,3)=YY(3)

       end if

       !write(*,*) i
       !write(*,*) m_l(i,1,1),m_l(i,1,2),m_l(i,1,3)
       !write(*,*) m_l(i,2,1),m_l(i,2,2),m_l(i,2,3)
       !write(*,*) m_l(i,3,1),m_l(i,3,2),m_l(i,3,3)
       !write(*,*) m_m(i,1),m_m(i,2),m_m(i,3)

       !write(*,*) m_A(i,1,1),m_A(i,1,2),m_l(i,1,3)
       !write(*,*) m_A(i,2,1),m_A(i,2,2),m_l(i,2,3)
       !write(*,*) m_A(i,3,1),m_A(i,3,2),m_l(i,3,3)
       !write(*,*) m_B(i,1),m_B(i,2),m_B(i,3)
       !pause
     


     end do
     !!!!!!!!!!!!!!!!!!!!!ת�ƾ������!!!!!!!!!!!!!!!!!!!!
     !PAUSE



     !!!!!!!!!!!!!!!!!!!!!z��������������!!!!!!!!!!!!!!!!!!!!


     CCC=0.
     DDD=0.
     do i=NH+1,NH+N_CLAD-1    !��Ҫ�޸�
       
       if (i==NH+1) then      !��Ҫ�޸�
         
         ccc(1,1)=ccc(1,1)+0.5*pie*(r_mc(i+1)**2-r_mc(i)**2)*(1.+m_A(i,1,1))
         ccc(1,2)=ccc(1,2)+0.5*pie*(r_mc(i+1)**2-r_mc(i)**2)*(0.+m_A(i,1,2))
         ccc(1,3)=ccc(1,3)+0.5*pie*(r_mc(i+1)**2-r_mc(i)**2)*(0.+m_A(i,1,3))
         ccc(2,1)=ccc(2,1)+0.5*pie*(r_mc(i+1)**2-r_mc(i)**2)*(0.+m_A(i,2,1))
         ccc(2,2)=ccc(2,2)+0.5*pie*(r_mc(i+1)**2-r_mc(i)**2)*(1.+m_A(i,2,2))
         ccc(2,3)=ccc(2,3)+0.5*pie*(r_mc(i+1)**2-r_mc(i)**2)*(0.+m_A(i,2,3))
         ccc(3,1)=ccc(3,1)+0.5*pie*(r_mc(i+1)**2-r_mc(i)**2)*(0.+m_A(i,3,1))
         ccc(3,2)=ccc(3,2)+0.5*pie*(r_mc(i+1)**2-r_mc(i)**2)*(0.+m_A(i,3,2))
         ccc(3,3)=ccc(3,3)+0.5*pie*(r_mc(i+1)**2-r_mc(i)**2)*(1.+m_A(i,3,3))

         ddd(1)=ddd(1)+0.5*pie*(r_mc(i+1)**2-r_mc(i)**2)*m_B(i,1)
         ddd(2)=ddd(2)+0.5*pie*(r_mc(i+1)**2-r_mc(i)**2)*m_B(i,2)
         ddd(3)=ddd(3)+0.5*pie*(r_mc(i+1)**2-r_mc(i)**2)*m_B(i,3)

       else
         
         ccc(1,1)=ccc(1,1)+0.5*pie*(r_mc(i+1)**2-r_mc(i)**2)*(m_A(i-1,1,1)+m_A(i,1,1))
         ccc(1,2)=ccc(1,2)+0.5*pie*(r_mc(i+1)**2-r_mc(i)**2)*(m_A(i-1,1,2)+m_A(i,1,2))
         ccc(1,3)=ccc(1,3)+0.5*pie*(r_mc(i+1)**2-r_mc(i)**2)*(m_A(i-1,1,3)+m_A(i,1,3))
         ccc(2,1)=ccc(2,1)+0.5*pie*(r_mc(i+1)**2-r_mc(i)**2)*(m_A(i-1,2,1)+m_A(i,2,1))
         ccc(2,2)=ccc(2,2)+0.5*pie*(r_mc(i+1)**2-r_mc(i)**2)*(m_A(i-1,2,2)+m_A(i,2,2))
         ccc(2,3)=ccc(2,3)+0.5*pie*(r_mc(i+1)**2-r_mc(i)**2)*(m_A(i-1,2,3)+m_A(i,2,3))
         ccc(3,1)=ccc(3,1)+0.5*pie*(r_mc(i+1)**2-r_mc(i)**2)*(m_A(i-1,3,1)+m_A(i,3,1))
         ccc(3,2)=ccc(3,2)+0.5*pie*(r_mc(i+1)**2-r_mc(i)**2)*(m_A(i-1,3,2)+m_A(i,3,2))
         ccc(3,3)=ccc(3,3)+0.5*pie*(r_mc(i+1)**2-r_mc(i)**2)*(m_A(i-1,3,3)+m_A(i,3,3))

         ddd(1)=ddd(1)+0.5*pie*(r_mc(i+1)**2-r_mc(i)**2)*(m_B(i-1,1)+m_B(i,1))
         ddd(2)=ddd(2)+0.5*pie*(r_mc(i+1)**2-r_mc(i)**2)*(m_B(i-1,2)+m_B(i,2))
         ddd(3)=ddd(3)+0.5*pie*(r_mc(i+1)**2-r_mc(i)**2)*(m_B(i-1,3)+m_B(i,3))

       end if


     end do


     !!!!!!!!!!!!!!!!!!!!!z��������������!!!!!!!!!!!!!!!!!!!!(�ⲿ���Թ�)

       !write(*,*) CCC(1,1),CCC(1,2),CCC(1,3)
       !write(*,*) CCC(2,1),CCC(2,2),CCC(2,3)
       !write(*,*) CCC(3,1),CCC(3,2),CCC(3,3)
       !write(*,*) DDD(1),DDD(2),DDD(3)
       !pause





     !!!!!!!!!!!!!!!!!!����һ���ڵ��ϵ�Ӧ������!!!!!!!!!!!!!!!!!!!!!

     AAA(1,1)=1.
     AAA(1,2)=0.
     AAA(1,3)=0.

     AAA(2,1)=m_A(NH+N_CLAD-1,1,1)
     AAA(2,2)=m_A(NH+N_CLAD-1,1,2)   !�޸�
     AAA(2,3)=m_A(NH+N_CLAD-1,1,3)

     AAA(3,1)=ccc(3,1)
     AAA(3,2)=ccc(3,2)
     AAA(3,3)=ccc(3,3)

     BBB(1)=-PRESS_INTER*1.D6   !�޸�
     BBB(2)=-m_B(NH+N_CLAD-1,1)-p_coolant*1.D6   !�޸�
     BBB(3)=-ddd(3)+force_z_clad   !�޸�

     call agaus2(AAA,bbb,XXX)

     stress_clad(1,1)=XXX(1)
     stress_clad(1,2)=XXX(2)   !�޸�
     stress_clad(1,3)=XXX(3)

     stress_equ(NH+1)=(0.5*((stress_clad(1,1)-stress_clad(1,2))**2+(stress_clad(1,1)-stress_clad(1,3))**2+(stress_clad(1,2)-stress_clad(1,3))**2))**0.5

     !write(*,*)  XXX
     !write(*,*)  stress_equ(NH+1)/1.d6
     !pause
     !!!!!!!!!!!!!!!!!!����һ���ڵ��ϵ�Ӧ������!!!!!!!!!!!!!!!!!!!!!
     


     !!!!!!!!!!!!!!!!!!��������ڵ��ϵ�Ӧ������!!!!!!!!!!!!!!!!!!!!!
    
     do i=2,N_CLAD    !�޸ģ��ӵڶ����ڵ㵽���Ľڵ�

       stress_clad(i,1)=m_A(NH+i-1,1,1)*stress_clad(1,1)+m_A(NH+i-1,1,2)*stress_clad(1,2)+m_A(NH+i-1,1,3)*stress_clad(1,3)+m_b(NH+i-1,1)
       stress_clad(i,2)=m_A(NH+i-1,2,1)*stress_clad(1,1)+m_A(NH+i-1,2,2)*stress_clad(1,2)+m_A(NH+i-1,2,3)*stress_clad(1,3)+m_b(NH+i-1,2)
       stress_clad(i,3)=m_A(NH+i-1,3,1)*stress_clad(1,1)+m_A(NH+i-1,3,2)*stress_clad(1,2)+m_A(NH+i-1,3,3)*stress_clad(1,3)+m_b(NH+i-1,3)

       stress_equ(NH+i)=(0.5*((stress_clad(i,1)-stress_clad(i,2))**2+(stress_clad(i,1)-stress_clad(i,3))**2+(stress_clad(i,2)-stress_clad(i,3))**2))**0.5

       !write(*,*) i
       !write(*,*) stress_clad(i,1)
       !write(*,*) stress_clad(i,2)
       !write(*,*) stress_clad(i,3)
       !write(*,*) stress_equ(NH+i)/1.d6 
       !pause
     end do

     !!!!!!!!!!!!!!!!!!��������ڵ��ϵ�Ӧ������!!!!!!!!!!!!!!!!!!!!!
     !pause



     !!!!!!!!!!!!!!!!!����Ӧ����Ӧ��!!!!!!!!!!!!!!!!!

     do i=1,N_CLAD

        strain_clad(i,1)=1./m_elastic(nh+i)*(stress_clad(i,1)-Poisson_ration(nh+i)*(stress_clad(i,2)+stress_clad(i,3)))+strain_intrinsic(nh+i,1)
        strain_clad(i,2)=1./m_elastic(nh+i)*(stress_clad(i,2)-Poisson_ration(nh+i)*(stress_clad(i,1)+stress_clad(i,3)))+strain_intrinsic(nh+i,2)
        strain_clad(i,3)=1./m_elastic(nh+i)*(stress_clad(i,3)-Poisson_ration(nh+i)*(stress_clad(i,2)+stress_clad(i,1)))+strain_intrinsic(nh+i,3)
        !if (i==1) then
        !write(*,*) stress_clad(i,1),stress_clad(i,2),stress_clad(i,3)
        !write(*,*) m_elastic(nh+i),Poisson_ration(nh+i)
        !write(*,*) strain_clad(i,1),strain_clad(i,2),strain_clad(i,3)
        !end if
        !pause
     end do


     !!!!!!!!!!!!!!!!!����Ӧ����Ӧ��!!!!!!!!!!!!!!!!!
     !pause





     !!!!!!!!!!!!!!!!!!�жϽ����������!!!!!!!!!!!!!!!!
 
     factor=0
     do i=nh+1,nh+N_CLAD
          
      !write(*,*)  stress_equ(i),yield_stress(i)
       if (stress_equ(i)>yield_stress(i)) then

       !write(*,*) i,"��������"
       !write(*,*)  i,stress_equ(i),yield_stress(i)
       !pause

         !d_strain_plastic(i,1)=0.01
         !d_strain_plastic(i,2)=-0.005
         !d_strain_plastic(i,3)=-0.005

         factor=factor+1

       else
         
         !d_strain_plastic(i,1)=0.
         !d_strain_plastic(i,2)=0.
         !d_strain_plastic(i,3)=0.

         factor=factor+0

       end if
       

     end do

     !!!!!!!!!!!!!!!!!!�жϽ����������!!!!!!!!!!!!!!!!

     !pause




if (factor==0) then




       !exit  !û���κ�һ��������Ա��Σ�������do while ѭ��

       do i=nh+1,nh+N_CLAD
              
       yield_stress111(i)=yield_stress(i)
       strain_plastic111(i,1)=strain_plastic(i,1)
       strain_plastic111(i,2)=strain_plastic(i,2)
       strain_plastic111(i,3)=strain_plastic(i,3)
       end do



else 
       

        !write(*,*)  "��������ѭ��"


        !pause

     !!!!!!!!!!!!!!!!!!��ֵ���Ա���!!!!!!!!!!!!!!!!

     !factor=0
     do i=nh+1,nh+N_CLAD
         
       if (stress_equ(i)>yield_stress(i)) then

       write(*,*) i,"��������"
   

         d_strain_plastic(i,1)=-0.0001
         d_strain_plastic(i,2)=0.
         d_strain_plastic(i,3)=+0.0001

         !factor=factor+1

       else
         
         d_strain_plastic(i,1)=0.
         d_strain_plastic(i,2)=0.
         d_strain_plastic(i,3)=0.

         !factor=factor+0

       end if

       !write(*,*) i
       !write(*,*) d_strain_plastic(i,1),d_strain_plastic(i,2),d_strain_plastic(i,3)

     end do

     !!!!!!!!!!!!!!!!!!��ֵ���Ա���!!!!!!!!!!!!!!!!







      
       

     
     factor2=1.
     factor3=1.





do while ((factor3>0.05))    !do while�����Ե���

     d_strain_plastic_js=d_strain_plastic


     do i=1,N_CLAD
     stress_clad_js(i,1)=stress_clad(i,1)
     stress_clad_js(i,2)=stress_clad(i,2)
     stress_clad_js(i,3)=stress_clad(i,3)
     !write(*,*) stress_clad_js(i,1),stress_clad_js(i,2),stress_clad_js(i,3)
     end do

       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !!!!!!!!!!!!!���µ��Լ���!!!!!!!!!!!!!!
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  do i=1,N_CLAD

       temperature_cladding=T_NOW(nh+1)+(i-1)*(T_NOW(nh+2)-T_NOW(nh+1))/(N_CLAD-1)
     if (temperature_cladding<1083.) then
       strain_thermal_cladding_r=1.*(4.95d-6*temperature_cladding-1.485d-3)
       strain_thermal_cladding_z=1.*(1.26d-5*temperature_cladding-3.78d-3)
       strain_thermal_cladding_c=strain_thermal_cladding_r
     else if (temperature_cladding<1244.) then
       strain_thermal_cladding_r=(2.77763+1.09822*cos(3.14159265*(temperature_cladding-1083.)/161.))/1000.
       strain_thermal_cladding_z=(2.77763+1.09822*cos(3.14159265*(temperature_cladding-1083.)/161.))/1000.
       strain_thermal_cladding_c=strain_thermal_cladding_r
     else if (temperature_cladding<2098.) then
       strain_thermal_cladding_r=9.7d-6*temperature_cladding-1.04d-2
       strain_thermal_cladding_z=9.76d-6*temperature_cladding-4.4d-3
       strain_thermal_cladding_c=strain_thermal_cladding_r
     else
       strain_thermal_cladding_r=9.7d-6*2098.-1.04d-2
       strain_thermal_cladding_z=9.76d-6*2098.-4.4d-3
       strain_thermal_cladding_c=strain_thermal_cladding_r
     end if

     !write(*,*) temperature_cladding
          !����Ӧ�������㣬ע���Ǵ����
          strain_intrinsic(NH+i,1)=strain_thermal_cladding_r+strain_plastic(NH+i,1)+d_strain_plastic(NH+i,1)+strain_creep(i,1)+d_strain_creep(i,1)
          strain_intrinsic(NH+i,2)=strain_thermal_cladding_c+strain_plastic(NH+i,2)+d_strain_plastic(NH+i,2)+strain_creep(i,2)+d_strain_creep(i,2)
          strain_intrinsic(NH+i,3)=strain_thermal_cladding_z+strain_plastic(NH+i,3)+d_strain_plastic(NH+i,3)+strain_creep(i,3)+d_strain_creep(i,3)
          !write(*,*) strain_intrinsic(NH+i,1),strain_intrinsic(NH+i,2),strain_intrinsic(NH+i,3)
          !write(*,*) strain_thermal_cladding_r,strain_thermal_cladding_c,strain_thermal_cladding_z
   end do
     

 

              !!!!!!!!!!!!!!!!!!!!!�������ƾ������!!!!!!!!!!!!!!!!!!!!
     do i=NH+1,NH+N_CLAD-1           !о��Ͱ���δ�Ӵ�ʱ����NH��������Ҫ����Ϊо��Ͱ���֮���޹���

       m_l(i,1,1)=1.-(r_mc(i+1)-r_mc(i))/r_mc(i+1)
       m_l(i,1,2)=(r_mc(i+1)-r_mc(i))/r_mc(i+1)
       m_l(i,1,3)=0.
       m_m(i,1)=0.

       m_l(i,2,1)=(r_mc(i+1)-r_mc(i))/r_mc(i+1)+m_elastic(i+1)/(1.-Poisson_ration(i+1))&
       &*(Poisson_ration(i+1)/m_elastic(i+1)-Poisson_ration(i)/m_elastic(i))
       m_l(i,2,2)=1.-(r_mc(i+1)-r_mc(i))/r_mc(i+1)+m_elastic(i+1)/(1.-Poisson_ration(i+1)**2)&
       &*((Poisson_ration(i+1)*Poisson_ration(i+1)/m_elastic(i+1)-Poisson_ration(i+1)*Poisson_ration(i)/m_elastic(i))&
       &-(1./m_elastic(i+1)-1./m_elastic(i)))
       m_l(i,2,3)=m_elastic(i+1)/(1.-Poisson_ration(i+1)**2)*((Poisson_ration(i+1)/m_elastic(i+1)-Poisson_ration(i)/m_elastic(i))&
       &-(Poisson_ration(i+1)/m_elastic(i+1)-Poisson_ration(i+1)/m_elastic(i)))
       m_m(i,2)=m_elastic(i+1)/(1.-Poisson_ration(i+1)**2)*(((r_mc(i+1)-r_mc(i))/r_mc(i+1)*(strain_intrinsic(i+1,1)-strain_intrinsic(i+1,2)))&
       &-(strain_intrinsic(i+1,2)-strain_intrinsic(i,2))-Poisson_ration(i+1)*(strain_intrinsic(i+1,3)-strain_intrinsic(i,3)))

       m_l(i,3,1)=m_elastic(i+1)/(1.-Poisson_ration(i+1))*(Poisson_ration(i+1)/m_elastic(i+1)-Poisson_ration(i)/m_elastic(i))
       m_l(i,3,2)=m_elastic(i+1)/(1.-Poisson_ration(i+1)**2)*((Poisson_ration(i+1)/m_elastic(i+1)-Poisson_ration(i)/m_elastic(i))&
       &-(Poisson_ration(i+1)/m_elastic(i+1)-Poisson_ration(i+1)/m_elastic(i)))
       m_l(i,3,3)=1.+m_elastic(i+1)/(1.-Poisson_ration(i+1)**2)&
       &*((Poisson_ration(i+1)*Poisson_ration(i+1)/m_elastic(i+1)-Poisson_ration(i+1)*Poisson_ration(i)/m_elastic(i))&
       &-(1./m_elastic(i+1)-1./m_elastic(i)))
       m_m(i,3)=m_elastic(i+1)/(1.-Poisson_ration(i+1)**2)*(((r_mc(i+1)-r_mc(i))/r_mc(i+1)*&
       &Poisson_ration(i+1)*(strain_intrinsic(i+1,1)-strain_intrinsic(i+1,2)))&
       &-Poisson_ration(i+1)*(strain_intrinsic(i+1,2)-strain_intrinsic(i,2))-(strain_intrinsic(i+1,3)-strain_intrinsic(i,3)))

       !write(*,*) i
       !write(*,*) m_l(i,1,1),m_l(i,1,2),m_l(i,1,3)
       !write(*,*) m_l(i,2,1),m_l(i,2,2),m_l(i,2,3)
       !write(*,*) m_l(i,3,1),m_l(i,3,2),m_l(i,3,3)
       !write(*,*) m_m(i,1),m_m(i,2),m_m(i,3)
       !pause
     end do
     !!!!!!!!!!!!!!!!!!!!!�������ƾ������!!!!!!!!!!!!!!!!!!!!
   
   !pause



     !!!!!!!!!!!!!!!!!!!!!ת�ƾ������!!!!!!!!!!!!!!!!!!!!
     do i=NH+1,NH+N_CLAD-1

       if (i==NH+1) then   !����ط���Ҫ��

         m_A(i,1,1)=m_l(i,1,1)
         m_A(i,1,2)=m_l(i,1,2)
         m_A(i,1,3)=m_l(i,1,3)
         m_A(i,2,1)=m_l(i,2,1)
         m_A(i,2,2)=m_l(i,2,2)
         m_A(i,2,3)=m_l(i,2,3)
         m_A(i,3,1)=m_l(i,3,1)
         m_A(i,3,2)=m_l(i,3,2)
         m_A(i,3,3)=m_l(i,3,3)

         m_B(i,1)=m_m(i,1)
         m_B(i,2)=m_m(i,2)
         m_B(i,3)=m_m(i,3)

       else
         

         AA(1,1)=m_l(i,1,1)
         AA(1,2)=m_l(i,1,2)
         AA(1,3)=m_l(i,1,3)
         AA(2,1)=m_l(i,2,1)
         AA(2,2)=m_l(i,2,2)
         AA(2,3)=m_l(i,2,3)
         AA(3,1)=m_l(i,3,1)
         AA(3,2)=m_l(i,3,2)
         AA(3,3)=m_l(i,3,3)

         BB(1,1)=m_A(i-1,1,1)
         BB(1,2)=m_A(i-1,1,2)
         BB(1,3)=m_A(i-1,1,3)
         BB(2,1)=m_A(i-1,2,1)
         BB(2,2)=m_A(i-1,2,2)
         BB(2,3)=m_A(i-1,2,3)
         BB(3,1)=m_A(i-1,3,1)
         BB(3,2)=m_A(i-1,3,2)
         BB(3,3)=m_A(i-1,3,3)


         CC(1)=m_m(i,1)
         CC(2)=m_m(i,2)
         CC(3)=m_m(i,3)

         DD(1)=m_B(i-1,1)
         DD(2)=m_B(i-1,2)
         DD(3)=m_B(i-1,3)

         call matrix_multiply(AA,BB,XX)
         call series_multiply(AA,DD,YY)
         YY=YY+CC

         m_A(i,1,1)=XX(1,1)
         m_A(i,1,2)=XX(1,2)
         m_A(i,1,3)=XX(1,3)
         m_A(i,2,1)=XX(2,1)
         m_A(i,2,2)=XX(2,2)
         m_A(i,2,3)=XX(2,3)
         m_A(i,3,1)=XX(3,1)
         m_A(i,3,2)=XX(3,2)
         m_A(i,3,3)=XX(3,3)

         m_B(i,1)=YY(1)
         m_B(i,2)=YY(2)
         m_B(i,3)=YY(3)

       end if

       !write(*,*) i
       !write(*,*) m_l(i,1,1),m_l(i,1,2),m_l(i,1,3)
       !write(*,*) m_l(i,2,1),m_l(i,2,2),m_l(i,2,3)
       !write(*,*) m_l(i,3,1),m_l(i,3,2),m_l(i,3,3)
       !write(*,*) m_m(i,1),m_m(i,2),m_m(i,3)

       !write(*,*) m_A(i,1,1),m_A(i,1,2),m_l(i,1,3)
       !write(*,*) m_A(i,2,1),m_A(i,2,2),m_l(i,2,3)
       !write(*,*) m_A(i,3,1),m_A(i,3,2),m_l(i,3,3)
       !write(*,*) m_B(i,1),m_B(i,2),m_B(i,3)
       !pause
     


     end do
     !!!!!!!!!!!!!!!!!!!!!ת�ƾ������!!!!!!!!!!!!!!!!!!!!
     !PAUSE



     !!!!!!!!!!!!!!!!!!!!!z��������������!!!!!!!!!!!!!!!!!!!!


     CCC=0.
     DDD=0.
     do i=NH+1,NH+N_CLAD-1    !��Ҫ�޸�
       
       if (i==NH+1) then      !��Ҫ�޸�
         
         ccc(1,1)=ccc(1,1)+0.5*pie*(r_mc(i+1)**2-r_mc(i)**2)*(1.+m_A(i,1,1))
         ccc(1,2)=ccc(1,2)+0.5*pie*(r_mc(i+1)**2-r_mc(i)**2)*(0.+m_A(i,1,2))
         ccc(1,3)=ccc(1,3)+0.5*pie*(r_mc(i+1)**2-r_mc(i)**2)*(0.+m_A(i,1,3))
         ccc(2,1)=ccc(2,1)+0.5*pie*(r_mc(i+1)**2-r_mc(i)**2)*(0.+m_A(i,2,1))
         ccc(2,2)=ccc(2,2)+0.5*pie*(r_mc(i+1)**2-r_mc(i)**2)*(1.+m_A(i,2,2))
         ccc(2,3)=ccc(2,3)+0.5*pie*(r_mc(i+1)**2-r_mc(i)**2)*(0.+m_A(i,2,3))
         ccc(3,1)=ccc(3,1)+0.5*pie*(r_mc(i+1)**2-r_mc(i)**2)*(0.+m_A(i,3,1))
         ccc(3,2)=ccc(3,2)+0.5*pie*(r_mc(i+1)**2-r_mc(i)**2)*(0.+m_A(i,3,2))
         ccc(3,3)=ccc(3,3)+0.5*pie*(r_mc(i+1)**2-r_mc(i)**2)*(1.+m_A(i,3,3))

         ddd(1)=ddd(1)+0.5*pie*(r_mc(i+1)**2-r_mc(i)**2)*m_B(i,1)
         ddd(2)=ddd(2)+0.5*pie*(r_mc(i+1)**2-r_mc(i)**2)*m_B(i,2)
         ddd(3)=ddd(3)+0.5*pie*(r_mc(i+1)**2-r_mc(i)**2)*m_B(i,3)

       else
         
         ccc(1,1)=ccc(1,1)+0.5*pie*(r_mc(i+1)**2-r_mc(i)**2)*(m_A(i-1,1,1)+m_A(i,1,1))
         ccc(1,2)=ccc(1,2)+0.5*pie*(r_mc(i+1)**2-r_mc(i)**2)*(m_A(i-1,1,2)+m_A(i,1,2))
         ccc(1,3)=ccc(1,3)+0.5*pie*(r_mc(i+1)**2-r_mc(i)**2)*(m_A(i-1,1,3)+m_A(i,1,3))
         ccc(2,1)=ccc(2,1)+0.5*pie*(r_mc(i+1)**2-r_mc(i)**2)*(m_A(i-1,2,1)+m_A(i,2,1))
         ccc(2,2)=ccc(2,2)+0.5*pie*(r_mc(i+1)**2-r_mc(i)**2)*(m_A(i-1,2,2)+m_A(i,2,2))
         ccc(2,3)=ccc(2,3)+0.5*pie*(r_mc(i+1)**2-r_mc(i)**2)*(m_A(i-1,2,3)+m_A(i,2,3))
         ccc(3,1)=ccc(3,1)+0.5*pie*(r_mc(i+1)**2-r_mc(i)**2)*(m_A(i-1,3,1)+m_A(i,3,1))
         ccc(3,2)=ccc(3,2)+0.5*pie*(r_mc(i+1)**2-r_mc(i)**2)*(m_A(i-1,3,2)+m_A(i,3,2))
         ccc(3,3)=ccc(3,3)+0.5*pie*(r_mc(i+1)**2-r_mc(i)**2)*(m_A(i-1,3,3)+m_A(i,3,3))

         ddd(1)=ddd(1)+0.5*pie*(r_mc(i+1)**2-r_mc(i)**2)*(m_B(i-1,1)+m_B(i,1))
         ddd(2)=ddd(2)+0.5*pie*(r_mc(i+1)**2-r_mc(i)**2)*(m_B(i-1,2)+m_B(i,2))
         ddd(3)=ddd(3)+0.5*pie*(r_mc(i+1)**2-r_mc(i)**2)*(m_B(i-1,3)+m_B(i,3))

       end if


     end do


     !!!!!!!!!!!!!!!!!!!!!z��������������!!!!!!!!!!!!!!!!!!!!(�ⲿ���Թ�)

       !write(*,*) CCC(1,1),CCC(1,2),CCC(1,3)
       !write(*,*) CCC(2,1),CCC(2,2),CCC(2,3)
       !write(*,*) CCC(3,1),CCC(3,2),CCC(3,3)
       !write(*,*) DDD(1),DDD(2),DDD(3)
       !pause





     !!!!!!!!!!!!!!!!!!����һ���ڵ��ϵ�Ӧ������!!!!!!!!!!!!!!!!!!!!!

     AAA(1,1)=1.
     AAA(1,2)=0.
     AAA(1,3)=0.

     AAA(2,1)=m_A(NH+N_CLAD-1,1,1)
     AAA(2,2)=m_A(NH+N_CLAD-1,1,2)   !�޸�
     AAA(2,3)=m_A(NH+N_CLAD-1,1,3)

     AAA(3,1)=ccc(3,1)
     AAA(3,2)=ccc(3,2)
     AAA(3,3)=ccc(3,3)

     BBB(1)=-PRESS_INTER*1.D6   !�޸�
     BBB(2)=-m_B(NH+N_CLAD-1,1)-p_coolant*1.D6   !�޸�
     BBB(3)=-ddd(3)+force_z_clad   !�޸�

     call agaus2(AAA,bbb,XXX)

     stress_clad(1,1)=XXX(1)
     stress_clad(1,2)=XXX(2)   !�޸�
     stress_clad(1,3)=XXX(3)

     stress_equ(NH+1)=(0.5*((stress_clad(1,1)-stress_clad(1,2))**2+(stress_clad(1,1)-stress_clad(1,3))**2+(stress_clad(1,2)-stress_clad(1,3))**2))**0.5

     !write(*,*)  XXX
     !write(*,*)  stress_equ(NH+1)/1.d6
     !pause
     !!!!!!!!!!!!!!!!!!����һ���ڵ��ϵ�Ӧ������!!!!!!!!!!!!!!!!!!!!!
     


     !!!!!!!!!!!!!!!!!!��������ڵ��ϵ�Ӧ������!!!!!!!!!!!!!!!!!!!!!
    
     do i=2,N_CLAD    !�޸ģ��ӵڶ����ڵ㵽���Ľڵ�

       stress_clad(i,1)=m_A(NH+i-1,1,1)*stress_clad(1,1)+m_A(NH+i-1,1,2)*stress_clad(1,2)+m_A(NH+i-1,1,3)*stress_clad(1,3)+m_b(NH+i-1,1)
       stress_clad(i,2)=m_A(NH+i-1,2,1)*stress_clad(1,1)+m_A(NH+i-1,2,2)*stress_clad(1,2)+m_A(NH+i-1,2,3)*stress_clad(1,3)+m_b(NH+i-1,2)
       stress_clad(i,3)=m_A(NH+i-1,3,1)*stress_clad(1,1)+m_A(NH+i-1,3,2)*stress_clad(1,2)+m_A(NH+i-1,3,3)*stress_clad(1,3)+m_b(NH+i-1,3)

       stress_equ(NH+i)=(0.5*((stress_clad(i,1)-stress_clad(i,2))**2+(stress_clad(i,1)-stress_clad(i,3))**2+(stress_clad(i,2)-stress_clad(i,3))**2))**0.5

       !write(*,*) i
       !write(*,*) stress_clad(i,1)
       !write(*,*) stress_clad(i,2)
       !write(*,*) stress_clad(i,3)
       !write(*,*) stress_equ(NH+i)/1.d6 
       !pause
     end do

     !!!!!!!!!!!!!!!!!!��������ڵ��ϵ�Ӧ������!!!!!!!!!!!!!!!!!!!!!
     !pause



     !!!!!!!!!!!!!!!!!����Ӧ����Ӧ��!!!!!!!!!!!!!!!!!

     do i=1,N_CLAD

        strain_clad(i,1)=1./m_elastic(nh+i)*(stress_clad(i,1)-Poisson_ration(nh+i)*(stress_clad(i,2)+stress_clad(i,3)))+strain_intrinsic(nh+i,1)
        strain_clad(i,2)=1./m_elastic(nh+i)*(stress_clad(i,2)-Poisson_ration(nh+i)*(stress_clad(i,1)+stress_clad(i,3)))+strain_intrinsic(nh+i,2)
        strain_clad(i,3)=1./m_elastic(nh+i)*(stress_clad(i,3)-Poisson_ration(nh+i)*(stress_clad(i,2)+stress_clad(i,1)))+strain_intrinsic(nh+i,3)
     
        !if (i==1) then
        !write(*,*) stress_clad(i,1),stress_clad(i,2),stress_clad(i,3)
        !write(*,*) m_elastic(nh+i),Poisson_ration(nh+i)
        !write(*,*) strain_clad(i,1),strain_clad(i,2),strain_clad(i,3)
        !end if
        !pause
     end do


     !!!!!!!!!!!!!!!!!����Ӧ����Ӧ��!!!!!!!!!!!!!!!!!



       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !!!!!!!!!!!!!���µ��Լ���!!!!!!!!!!!!!!
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



     !!!!!!!!!!!!!!!!!!���¸�ֵ���Ա���!!!!!!!!!!!!!!!!

     !factor=0
     do i=nh+1,nh+n_clad

       if ((stress_equ(i)>yield_stress(i)).and.(abs(d_strain_plastic(i,1))+abs(d_strain_plastic(i,2))+abs(d_strain_plastic(i,3))==0.)) then


           d_strain_plastic(i,1)=-0.0001
           d_strain_plastic(i,2)=0.
           d_strain_plastic(i,3)=+0.0001

           !write(*,*) i,"���ĸ�ֵ"
           !write(*,*) d_strain_plastic(i,1),d_strain_plastic(i,2),d_strain_plastic(i,3)

       end if





     end do

     !!!!!!!!!!!!!!!!!!���¸�ֵ���Ա���!!!!!!!!!!!!!!!!







  



       factor2=0.
       factor3=0.
       !!!!!!!!!!!!!!!!!!���¼������Է���!!!!!!!!!!!!!!!!!!!
       do i=nh+1,nh+n_clad

         if (stress_equ(i)>yield_stress(i)) then
         !write(*,*)  i,yield_stress(i),stress_equ(i)/1.d6

           d_strain_plastic_equ(i)=2.**0.5/3.*((d_strain_plastic(i,1)-d_strain_plastic(i,2))**2+&
           &(d_strain_plastic(i,1)-d_strain_plastic(i,3))**2+(d_strain_plastic(i,2)-d_strain_plastic(i,2))**2)**0.5

           stress_curve(i)=yield_stress(i)+d_strain_plastic_equ(i)*m_plastic(i)

           s1=stress_clad(i-nh,1)-(stress_clad(i-nh,1)+stress_clad(i-nh,2)+stress_clad(i-nh,3))/3.  !pa
           s2=stress_clad(i-nh,2)-(stress_clad(i-nh,1)+stress_clad(i-nh,2)+stress_clad(i-nh,3))/3.  !pa
           s3=stress_clad(i-nh,3)-(stress_clad(i-nh,1)+stress_clad(i-nh,2)+stress_clad(i-nh,3))/3.  !pa


           d_strain_plastic(i,1)=3./2.*d_strain_plastic_equ(i)/stress_curve(i)*s1
           d_strain_plastic(i,2)=3./2.*d_strain_plastic_equ(i)/stress_curve(i)*s2
           d_strain_plastic(i,3)=3./2.*d_strain_plastic_equ(i)/stress_curve(i)*s3

           factor2=factor2+abs((d_strain_plastic(i,1)-d_strain_plastic_js(i,1))/d_strain_plastic(i,1))&
           &+abs((d_strain_plastic(i,2)-d_strain_plastic_js(i,2))/d_strain_plastic(i,2))&                    
           &+abs((d_strain_plastic(i,3)-d_strain_plastic_js(i,3))/d_strain_plastic(i,3))

           !write(*,*) d_strain_plastic(i,1),d_strain_plastic(i,2),d_strain_plastic(i,3)
           !write(*,*) d_strain_plastic_js(i,1),d_strain_plastic_js(i,2),d_strain_plastic_js(i,3)


           yield_stress111(i)=stress_equ(i)
           strain_plastic111(i,1)=strain_plastic(i,1)+d_strain_plastic(i,1)
           strain_plastic111(i,2)=strain_plastic(i,2)+d_strain_plastic(i,2)
           strain_plastic111(i,3)=strain_plastic(i,3)+d_strain_plastic(i,3)

           !write(*,*) i,stress_curve(i)/1.d6,stress_equ(i)/1.d6
           !write(*,*) d_strain_plastic(i,1),d_strain_plastic(i,2),d_strain_plastic(i,3)
           
 
         else
  

           d_strain_plastic(i,1)=0.5*d_strain_plastic(i,1)
           d_strain_plastic(i,2)=0.5*d_strain_plastic(i,2)      !���Ϊ�㣬����0.5����ν
           d_strain_plastic(i,3)=0.5*d_strain_plastic(i,3)      !�����Ϊ�㣬��˵����ʼ�ĳ�ֵ����ˣ���С

           yield_stress111(i)=yield_stress(i)
           strain_plastic111(i,1)=strain_plastic(i,1)
           strain_plastic111(i,2)=strain_plastic(i,2)
           strain_plastic111(i,3)=strain_plastic(i,3)

           factor2=factor2+0.

           !write(*,*) i,stress_curve(i)/1.d6,stress_equ(i)/1.d6
           !write(*,*) d_strain_plastic(i,1),d_strain_plastic(i,2),d_strain_plastic(i,3)
 
         end if

        

         factor3=factor3+abs((stress_clad(i-nh,1)-stress_clad_js(i-nh,1))/stress_clad(i-nh,1))&
           &+abs((stress_clad(i-nh,2)-stress_clad_js(i-nh,2))/stress_clad(i-nh,2))&
           &+abs((stress_clad(i-nh,3)-stress_clad_js(i-nh,3))/stress_clad(i-nh,3))
        
        !write(*,*) stress_clad(i-nh,1),stress_clad(i-nh,2),stress_clad(i-nh,3)
        !write(*,*) stress_clad_js(i-nh,1),stress_clad_js(i-nh,2),stress_clad_js(i-nh,3)
        !write(*,*) i,factor2,factor3
        !write(*,*) d_strain_plastic(i,1),d_strain_plastic(i,2),d_strain_plastic(i,3)
        !pause

       end do
       !!!!!!!!!!!!!!!!!!���¼������Է���!!!!!!!!!!!!!!!!!!!


        !write(*,*) factor2,factor3
        !write(*,*) 
        !pause

    end do !do while�����Ե���




end if    !�нڵ��������






     !UR_PRA(1)=R(1)
     !!DO i=2,nh
     !UR_PRA(i)=(1.+strain_fuel(i,2))*R(i)           !�ڶ�����������
     !UR_PRA(i)=UR_PRA(i-1)+DLTX*(1.+0.5(l_strain(I)+l_strain(I-1)))
     !UR_PRA(i)=UR_PRA(i-1)+DLTX+DLTX*(l_strain(I)+l_strain(I-1))/2.        !ԭ�����������о��λ��
     !write(*,*) i
     !write(*,*) UR_PRA(i),R(I),UR_PRA(i)-R(I)
     !end do
     !write(*,*) UR_PRA(11),R(11),UR_PRA(11)-R(11

     !d_length_p=strain_fuel(11,3)*DH           



     UR_PRA(nh+1)=(1.+strain_clad(1,2))*R(nh+1) 
     UR_PRA(nh+2)=(1.+strain_clad(n_clad,2))*R(nh+n_clad) 
     d_length_c=strain_clad(1,3)*DH
    

 

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!����λ�Ƴ�������!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!����λ�Ƴ�������!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!����λ�Ƴ�������!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!real(8)::time_total,time_increment,Q_LINE,n_flux
!real(8)::strain_creep(N_CLAD,3),strain_creep111(N_CLAD,3)

!real(8)::creep_rate_th,creep_rate_irr
!real(8)::creep_rate_add,creep_rate_total
!real(8)::creep_saturated_primary,creep_total

!real(8)::d_strain_creep(N_CLAD,3),d_strain_creep_js(N_CLAD,3)
!integer:i2,i3,i4,i5,i6,i7
!real(8)::factor44




     !stress_clad(i,1/2/3)/1.d6
     !tem_cladding(i)
     !stress_equ(i+nh)/1.d6
     !m_elastic(nh+i)


     do i3=1,n_clad
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!�����μ���!!!!!!!!!!!!!!!!!     !����ѹ���߽Ӵ�ѹ��û�б������±�����һ�µ�
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
        creep_rate_th=5.47d8*m_elastic(nh+i3)/1.d6/tem_cladding(i3)*((sinh(320.*stress_equ(i3+nh)/m_elastic(nh+i3)))**3.5)&
        &*exp(-201./0.008314/tem_cladding(i3))

        if (tem_cladding(i3)<570.) then
        creep_rate_irr=1.87d-24*n_flux**0.85*stress_equ(i3+nh)*0.7944/1.d6
        else if (tem_cladding(i3)<625.) then
        creep_rate_irr=1.87d-24*n_flux**0.85*stress_equ(i3+nh)*(-3.18562+0.00699132*tem_cladding(i3))/1.d6
        else
        creep_rate_irr=1.87d-24*n_flux**0.85*stress_equ(i3+nh)*1.184/1.d6
        end if
        !write(*,*)  creep_rate_th,creep_rate_irr
        !pause

        creep_rate_add=creep_rate_th+creep_rate_irr
        creep_saturated_primary=0.0216*(creep_rate_add**0.109)*(2.-tanh(35500.*creep_rate_add))**(-2.05)
        creep_total=creep_saturated_primary*(1.-exp(-52.*(creep_rate_add*time_total)**0.5))+creep_rate_add*time_total
        creep_rate_total=52.*creep_saturated_primary*(creep_rate_add**0.5)/2./(time_total**0.5)*&
        &exp(-52.*(creep_rate_add*time_total)**0.5)+creep_rate_add
      
        !write(*,*)  creep_saturated_primary,creep_total,creep_rate_total
        !pause

        s1=stress_clad(i3,1)-(stress_clad(i3,1)+stress_clad(i3,2)+stress_clad(i3,3))/3.  !pa
        s2=stress_clad(i3,2)-(stress_clad(i3,1)+stress_clad(i3,2)+stress_clad(i3,3))/3.
        s3=stress_clad(i3,3)-(stress_clad(i3,1)+stress_clad(i3,2)+stress_clad(i3,3))/3.

        d_strain_creep(i3,1)=1.5*creep_rate_total*time_increment*s1/stress_equ(i3+nh)
        d_strain_creep(i3,2)=1.5*creep_rate_total*time_increment*s2/stress_equ(i3+nh)
        d_strain_creep(i3,3)=1.5*creep_rate_total*time_increment*s3/stress_equ(i3+nh)
    
        !write(*,*) d_creep_cladding_r,d_creep_cladding_c,d_creep_cladding_z
        !pause
            
        strain_creep111(i3,1)=strain_creep(i3,1)+d_strain_creep(i3,1)
        strain_creep111(i3,2)=strain_creep(i3,2)+d_strain_creep(i3,2)
        strain_creep111(i3,3)=strain_creep(i3,3)+d_strain_creep(i3,3)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!�����μ���!!!!!!!!!!!!!!!!!     !����ѹ���߽Ӵ�ѹ��û�б������±�����һ�µ�
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !write(*,*) i3
        !write(*,*) d_strain_creep(i3,1),d_strain_creep111(i3,2),d_strain_creep111(i3,3)
     end do

     !pause

     factor44=maxval(abs((d_strain_creep-d_strain_creep_js)/d_strain_creep))
     !write(*,*)  d_strain_creep(1,2),creep_rate_add,creep_rate_total
     !pause

end do    !!!!!!!!!!!!!!�����εĴ�ѭ��


     UR_PRA(nh+1)=(1.+strain_clad(1,2))*R(nh+1) 
     UR_PRA(nh+2)=(1.+strain_clad(n_clad,2))*R(nh+n_clad) 
     d_length_c=strain_clad(1,3)*DH

strain_z_lasttime111(1)=strain_fuel(nh,3)   !z��������һ���ڵ�
strain_z_lasttime111(2)=strain_clad(1,3)    !z��������һ���ڵ�

P_contact=0.










do i2=1,n_clad
d_strain_creep(i2,1)=1.d-4
d_strain_creep(i2,2)=5.d-5
d_strain_creep(i2,3)=5.d-5
d_strain_creep_js(i2,1)=0
d_strain_creep_js(i2,2)=0
d_strain_creep_js(i2,3)=0
end do

factor44=maxval(abs((d_strain_creep-d_strain_creep_js)/d_strain_creep))



do while (factor44>0.05)

do i2=1,n_clad
d_strain_creep_js(i2,1)=d_strain_creep(i2,1)
d_strain_creep_js(i2,2)=d_strain_creep(i2,2)
d_strain_creep_js(i2,3)=d_strain_creep(i2,3)
end do

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!�Ӵ������о�����ϼ���!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!�Ӵ������о�����ϼ���!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!�Ӵ������о�����ϼ���!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if (UR_PRA(nh+1)<=UR_PRA(nh)) then     
     
     !pause
     
     
     d_strain_plastic=0.
     

  
      do i=1,NH

          strain_intrinsic(i,1)=l_strain(I)+strain_plastic(i,1)+d_strain_plastic(i,1)
          strain_intrinsic(i,2)=l_strain(I)+strain_plastic(i,2)+d_strain_plastic(i,2)
          strain_intrinsic(i,3)=l_strain(I)+strain_plastic(i,3)+d_strain_plastic(i,3)

      end do
     

      do i=1,N_CLAD
             tem_cladding(i)=T_NOW(nh+1)+(i-1)*(T_NOW(nh+2)-T_NOW(nh+1))/(N_CLAD-1)
             temperature_cladding=T_NOW(nh+1)+(i-1)*(T_NOW(nh+2)-T_NOW(nh+1))/(N_CLAD-1)
             if (temperature_cladding<1083.) then
                strain_thermal_cladding_r=1.*(4.95d-6*temperature_cladding-1.485d-3)
                strain_thermal_cladding_z=1.*(1.26d-5*temperature_cladding-3.78d-3)
                strain_thermal_cladding_c=strain_thermal_cladding_r
             else if (temperature_cladding<1244.) then
                strain_thermal_cladding_r=(2.77763+1.09822*cos(3.14159265*(temperature_cladding-1083.)/161.))/1000.
                strain_thermal_cladding_z=(2.77763+1.09822*cos(3.14159265*(temperature_cladding-1083.)/161.))/1000.
                strain_thermal_cladding_c=strain_thermal_cladding_r
             else if (temperature_cladding<2098.) then
                strain_thermal_cladding_r=9.7d-6*temperature_cladding-1.04d-2
                strain_thermal_cladding_z=9.76d-6*temperature_cladding-4.4d-3
                strain_thermal_cladding_c=strain_thermal_cladding_r
             else
                strain_thermal_cladding_r=9.7d-6*2098.-1.04d-2
                strain_thermal_cladding_z=9.76d-6*2098.-4.4d-3
                strain_thermal_cladding_c=strain_thermal_cladding_r
             end if

          strain_intrinsic(NH+i,1)=strain_thermal_cladding_r+strain_plastic(NH+i,1)+d_strain_plastic(NH+i,1)+strain_creep(i,1)+d_strain_creep(i,1)
          strain_intrinsic(NH+i,2)=strain_thermal_cladding_c+strain_plastic(NH+i,2)+d_strain_plastic(NH+i,2)+strain_creep(i,2)+d_strain_creep(i,2)
          strain_intrinsic(NH+i,3)=strain_thermal_cladding_z+strain_plastic(NH+i,3)+d_strain_plastic(NH+i,3)+strain_creep(i,3)+d_strain_creep(i,3)
          !write(*,*) strain_intrinsic(NH+i,1),strain_intrinsic(NH+i,2),strain_intrinsic(NH+i,3)
          !write(*,*) strain_thermal_cladding_r,strain_thermal_cladding_c,strain_thermal_cladding_z
      end do





 

     !!!!!!!!!!!!!!!!!!!!!�������ƾ������!!!!!!!!!!!!!!!!!!!!
     do i=1,NH-1

       m_l(i,1,1)=1.-(r_mc(i+1)-r_mc(i))/r_mc(i+1)
       m_l(i,1,2)=(r_mc(i+1)-r_mc(i))/r_mc(i+1)
       m_l(i,1,3)=0.
       m_m(i,1)=0.

       m_l(i,2,1)=(r_mc(i+1)-r_mc(i))/r_mc(i+1)+m_elastic(i+1)/(1.-Poisson_ration(i+1))&
       &*(Poisson_ration(i+1)/m_elastic(i+1)-Poisson_ration(i)/m_elastic(i))
       m_l(i,2,2)=1.-(r_mc(i+1)-r_mc(i))/r_mc(i+1)+m_elastic(i+1)/(1.-Poisson_ration(i+1)**2)&
       &*((Poisson_ration(i+1)*Poisson_ration(i+1)/m_elastic(i+1)-Poisson_ration(i+1)*Poisson_ration(i)/m_elastic(i))&
       &-(1./m_elastic(i+1)-1./m_elastic(i)))
       m_l(i,2,3)=m_elastic(i+1)/(1.-Poisson_ration(i+1)**2)*((Poisson_ration(i+1)/m_elastic(i+1)-Poisson_ration(i)/m_elastic(i))&
       &-(Poisson_ration(i+1)/m_elastic(i+1)-Poisson_ration(i+1)/m_elastic(i)))
       m_m(i,2)=m_elastic(i+1)/(1.-Poisson_ration(i+1)**2)*(((r_mc(i+1)-r_mc(i))/r_mc(i+1)*(strain_intrinsic(i+1,1)-strain_intrinsic(i+1,2)))&
       &-(strain_intrinsic(i+1,2)-strain_intrinsic(i,2))-Poisson_ration(i+1)*(strain_intrinsic(i+1,3)-strain_intrinsic(i,3)))

       m_l(i,3,1)=m_elastic(i+1)/(1.-Poisson_ration(i+1))*(Poisson_ration(i+1)/m_elastic(i+1)-Poisson_ration(i)/m_elastic(i))
       m_l(i,3,2)=m_elastic(i+1)/(1.-Poisson_ration(i+1)**2)*((Poisson_ration(i+1)/m_elastic(i+1)-Poisson_ration(i)/m_elastic(i))&
       &-(Poisson_ration(i+1)/m_elastic(i+1)-Poisson_ration(i+1)/m_elastic(i)))
       m_l(i,3,3)=1.+m_elastic(i+1)/(1.-Poisson_ration(i+1)**2)&
       &*((Poisson_ration(i+1)*Poisson_ration(i+1)/m_elastic(i+1)-Poisson_ration(i+1)*Poisson_ration(i)/m_elastic(i))&
       &-(1./m_elastic(i+1)-1./m_elastic(i)))
       m_m(i,3)=m_elastic(i+1)/(1.-Poisson_ration(i+1)**2)*(((r_mc(i+1)-r_mc(i))/r_mc(i+1)*&
       &Poisson_ration(i+1)*(strain_intrinsic(i+1,1)-strain_intrinsic(i+1,2)))&
       &-Poisson_ration(i+1)*(strain_intrinsic(i+1,2)-strain_intrinsic(i,2))-(strain_intrinsic(i+1,3)-strain_intrinsic(i,3)))

       !write(*,*) i
       !write(*,*) m_l(i,1,1),m_l(i,1,2),m_l(i,1,3)
       !write(*,*) m_l(i,2,1),m_l(i,2,2),m_l(i,2,3)
       !write(*,*) m_l(i,3,1),m_l(i,3,2),m_l(i,3,3)
       !write(*,*) m_m(i,1),m_m(i,2),m_m(i,3)
       !pause
     end do







       m_l(NH,1,1)=1.
       m_l(NH,1,2)=0.
       m_l(NH,1,3)=0.
       m_m(NH,1)=0.

       m_l(NH,2,1)=m_elastic(nh+1)*(r_mc(nh+1)*(Poisson_ration(nh+1)**2)/m_elastic(nh+1)+r_mc(nh+1)*Poisson_ration(nh+1)/m_elastic(nh+1)&
                  &-r_mc(nh+1)*Poisson_ration(nh+1)*Poisson_ration(nh)/m_elastic(nh)-r_mc(nh)*Poisson_ration(nh)/m_elastic(nh))&
                  &/(r_mc(nh+1)*(1.-Poisson_ration(nh+1)**2))
       m_l(NH,2,2)=m_elastic(nh+1)*(r_mc(nh)/m_elastic(nh)-r_mc(nh+1)*Poisson_ration(nh+1)*Poisson_ration(nh)/m_elastic(nh))&
                  &/(r_mc(nh+1)*(1.-Poisson_ration(nh+1)**2))
       m_l(NH,2,3)=m_elastic(nh+1)*(r_mc(nh+1)*Poisson_ration(nh+1)/m_elastic(nh)-r_mc(nh)*Poisson_ration(nh)/m_elastic(nh))&
                  &/(r_mc(nh+1)*(1.-Poisson_ration(nh+1)**2))
       m_m(NH,2)=m_elastic(nh+1)*(r_mc(nh+1)*Poisson_ration(nh+1)*((strain_intrinsic(nh,3)-strain_intrinsic(nh+1,3))&
                 &+(strain_z_lasttime111(2)-strain_z_lasttime111(1)))+(r_mc(nh)*strain_intrinsic(nh,2)-r_mc(nh+1)*strain_intrinsic(nh+1,2))&
                 &-(r_mc(nh+1)-r_mc(nh)))/(r_mc(nh+1)*(1.-Poisson_ration(nh+1)**2))

       m_l(NH,3,1)=m_elastic(nh+1)*(r_mc(nh+1)*Poisson_ration(nh+1)/m_elastic(nh+1)+r_mc(nh+1)*(Poisson_ration(nh+1)**2)/m_elastic(nh+1)&
                  &-r_mc(nh+1)*Poisson_ration(nh)/m_elastic(nh)-r_mc(nh)*Poisson_ration(nh+1)*Poisson_ration(nh)/m_elastic(nh))&
                  &/(r_mc(nh+1)*(1.-Poisson_ration(nh+1)**2))
       m_l(NH,3,2)=m_elastic(nh+1)*(r_mc(nh)*Poisson_ration(nh+1)/m_elastic(nh)-r_mc(nh+1)*Poisson_ration(nh)/m_elastic(nh))&
                  &/(r_mc(nh+1)*(1.-Poisson_ration(nh+1)**2))
       m_l(NH,3,3)=m_elastic(nh+1)*(r_mc(nh+1)/m_elastic(nh)-r_mc(nh)*Poisson_ration(nh+1)*Poisson_ration(nh)/m_elastic(nh))&
                  &/(r_mc(nh+1)*(1.-Poisson_ration(nh+1)**2))
       m_m(NH,3)=m_elastic(nh+1)*(r_mc(nh+1)*((strain_intrinsic(nh,3)-strain_intrinsic(nh+1,3))&
                 &+(strain_z_lasttime111(2)-strain_z_lasttime111(1)))+Poisson_ration(nh+1)*((r_mc(nh)*strain_intrinsic(nh,2)&
                 &-r_mc(nh+1)*strain_intrinsic(nh+1,2))-(r_mc(nh+1)-r_mc(nh))))/(r_mc(nh+1)*(1.-Poisson_ration(nh+1)**2))








     do i=NH+1,NH+N_CLAD-1           !о��Ͱ���δ�Ӵ�ʱ����NH��������Ҫ����Ϊо��Ͱ���֮���޹���

       m_l(i,1,1)=1.-(r_mc(i+1)-r_mc(i))/r_mc(i+1)
       m_l(i,1,2)=(r_mc(i+1)-r_mc(i))/r_mc(i+1)
       m_l(i,1,3)=0.
       m_m(i,1)=0.

       m_l(i,2,1)=(r_mc(i+1)-r_mc(i))/r_mc(i+1)+m_elastic(i+1)/(1.-Poisson_ration(i+1))&
       &*(Poisson_ration(i+1)/m_elastic(i+1)-Poisson_ration(i)/m_elastic(i))
       m_l(i,2,2)=1.-(r_mc(i+1)-r_mc(i))/r_mc(i+1)+m_elastic(i+1)/(1.-Poisson_ration(i+1)**2)&
       &*((Poisson_ration(i+1)*Poisson_ration(i+1)/m_elastic(i+1)-Poisson_ration(i+1)*Poisson_ration(i)/m_elastic(i))&
       &-(1./m_elastic(i+1)-1./m_elastic(i)))
       m_l(i,2,3)=m_elastic(i+1)/(1.-Poisson_ration(i+1)**2)*((Poisson_ration(i+1)/m_elastic(i+1)-Poisson_ration(i)/m_elastic(i))&
       &-(Poisson_ration(i+1)/m_elastic(i+1)-Poisson_ration(i+1)/m_elastic(i)))
       m_m(i,2)=m_elastic(i+1)/(1.-Poisson_ration(i+1)**2)*(((r_mc(i+1)-r_mc(i))/r_mc(i+1)*(strain_intrinsic(i+1,1)-strain_intrinsic(i+1,2)))&
       &-(strain_intrinsic(i+1,2)-strain_intrinsic(i,2))-Poisson_ration(i+1)*(strain_intrinsic(i+1,3)-strain_intrinsic(i,3)))

       m_l(i,3,1)=m_elastic(i+1)/(1.-Poisson_ration(i+1))*(Poisson_ration(i+1)/m_elastic(i+1)-Poisson_ration(i)/m_elastic(i))
       m_l(i,3,2)=m_elastic(i+1)/(1.-Poisson_ration(i+1)**2)*((Poisson_ration(i+1)/m_elastic(i+1)-Poisson_ration(i)/m_elastic(i))&
       &-(Poisson_ration(i+1)/m_elastic(i+1)-Poisson_ration(i+1)/m_elastic(i)))
       m_l(i,3,3)=1.+m_elastic(i+1)/(1.-Poisson_ration(i+1)**2)&
       &*((Poisson_ration(i+1)*Poisson_ration(i+1)/m_elastic(i+1)-Poisson_ration(i+1)*Poisson_ration(i)/m_elastic(i))&
       &-(1./m_elastic(i+1)-1./m_elastic(i)))
       m_m(i,3)=m_elastic(i+1)/(1.-Poisson_ration(i+1)**2)*(((r_mc(i+1)-r_mc(i))/r_mc(i+1)*&
       &Poisson_ration(i+1)*(strain_intrinsic(i+1,1)-strain_intrinsic(i+1,2)))&
       &-Poisson_ration(i+1)*(strain_intrinsic(i+1,2)-strain_intrinsic(i,2))-(strain_intrinsic(i+1,3)-strain_intrinsic(i,3)))

       !write(*,*) i
       !write(*,*) m_l(i,1,1),m_l(i,1,2),m_l(i,1,3)
       !write(*,*) m_l(i,2,1),m_l(i,2,2),m_l(i,2,3)
       !write(*,*) m_l(i,3,1),m_l(i,3,2),m_l(i,3,3)
       !write(*,*) m_m(i,1),m_m(i,2),m_m(i,3)
       !pause
     end do
















     !!!!!!!!!!!!!!!!!!!!!�������ƾ������!!!!!!!!!!!!!!!!!!!!
   

     !!!!!!!!!!!!!!!!!!!!!ת�ƾ������!!!!!!!!!!!!!!!!!!!!
     do i=1,NH+N_CLAD-1

       if (i==1) then

         m_A(i,1,1)=m_l(i,1,1)
         m_A(i,1,2)=m_l(i,1,2)
         m_A(i,1,3)=m_l(i,1,3)
         m_A(i,2,1)=m_l(i,2,1)
         m_A(i,2,2)=m_l(i,2,2)
         m_A(i,2,3)=m_l(i,2,3)
         m_A(i,3,1)=m_l(i,3,1)
         m_A(i,3,2)=m_l(i,3,2)
         m_A(i,3,3)=m_l(i,3,3)

         m_B(i,1)=m_m(i,1)
         m_B(i,2)=m_m(i,2)
         m_B(i,3)=m_m(i,3)

       else
         

         AA(1,1)=m_l(i,1,1)
         AA(1,2)=m_l(i,1,2)
         AA(1,3)=m_l(i,1,3)
         AA(2,1)=m_l(i,2,1)
         AA(2,2)=m_l(i,2,2)
         AA(2,3)=m_l(i,2,3)
         AA(3,1)=m_l(i,3,1)
         AA(3,2)=m_l(i,3,2)
         AA(3,3)=m_l(i,3,3)

         BB(1,1)=m_A(i-1,1,1)
         BB(1,2)=m_A(i-1,1,2)
         BB(1,3)=m_A(i-1,1,3)
         BB(2,1)=m_A(i-1,2,1)
         BB(2,2)=m_A(i-1,2,2)
         BB(2,3)=m_A(i-1,2,3)
         BB(3,1)=m_A(i-1,3,1)
         BB(3,2)=m_A(i-1,3,2)
         BB(3,3)=m_A(i-1,3,3)


         CC(1)=m_m(i,1)
         CC(2)=m_m(i,2)
         CC(3)=m_m(i,3)

         DD(1)=m_B(i-1,1)
         DD(2)=m_B(i-1,2)
         DD(3)=m_B(i-1,3)

         call matrix_multiply(AA,BB,XX)
         call series_multiply(AA,DD,YY)
         YY=YY+CC

         m_A(i,1,1)=XX(1,1)
         m_A(i,1,2)=XX(1,2)
         m_A(i,1,3)=XX(1,3)
         m_A(i,2,1)=XX(2,1)
         m_A(i,2,2)=XX(2,2)
         m_A(i,2,3)=XX(2,3)
         m_A(i,3,1)=XX(3,1)
         m_A(i,3,2)=XX(3,2)
         m_A(i,3,3)=XX(3,3)

         m_B(i,1)=YY(1)
         m_B(i,2)=YY(2)
         m_B(i,3)=YY(3)

       end if

       !write(*,*) i
       !write(*,*) m_l(i,1,1),m_l(i,1,2),m_l(i,1,3)
       !write(*,*) m_l(i,2,1),m_l(i,2,2),m_l(i,2,3)
       !write(*,*) m_l(i,3,1),m_l(i,3,2),m_l(i,3,3)
       !write(*,*) m_m(i,1),m_m(i,2),m_m(i,3)

       !write(*,*) m_A(i,1,1),m_A(i,1,2),m_l(i,1,3)
       !write(*,*) m_A(i,2,1),m_A(i,2,2),m_l(i,2,3)
       !write(*,*) m_A(i,3,1),m_A(i,3,2),m_l(i,3,3)
       !write(*,*) m_B(i,1),m_B(i,2),m_B(i,3)
       !pause
     


     end do
     !!!!!!!!!!!!!!!!!!!!!ת�ƾ������!!!!!!!!!!!!!!!!!!!!
     !PAUSE



     !!!!!!!!!!!!!!!!!!!!!z��������������!!!!!!!!!!!!!!!!!!!!


     CCC=0.
     DDD=0.
     do i=1,NH+N_CLAD-1
       
       if (i==1) then
         
         ccc(1,1)=ccc(1,1)+0.5*pie*(r_mc(i+1)**2-r_mc(i)**2)*(1.+m_A(i,1,1))
         ccc(1,2)=ccc(1,2)+0.5*pie*(r_mc(i+1)**2-r_mc(i)**2)*(0.+m_A(i,1,2))
         ccc(1,3)=ccc(1,3)+0.5*pie*(r_mc(i+1)**2-r_mc(i)**2)*(0.+m_A(i,1,3))
         ccc(2,1)=ccc(2,1)+0.5*pie*(r_mc(i+1)**2-r_mc(i)**2)*(0.+m_A(i,2,1))
         ccc(2,2)=ccc(2,2)+0.5*pie*(r_mc(i+1)**2-r_mc(i)**2)*(1.+m_A(i,2,2))
         ccc(2,3)=ccc(2,3)+0.5*pie*(r_mc(i+1)**2-r_mc(i)**2)*(0.+m_A(i,2,3))
         ccc(3,1)=ccc(3,1)+0.5*pie*(r_mc(i+1)**2-r_mc(i)**2)*(0.+m_A(i,3,1))
         ccc(3,2)=ccc(3,2)+0.5*pie*(r_mc(i+1)**2-r_mc(i)**2)*(0.+m_A(i,3,2))
         ccc(3,3)=ccc(3,3)+0.5*pie*(r_mc(i+1)**2-r_mc(i)**2)*(1.+m_A(i,3,3))

         ddd(1)=ddd(1)+0.5*pie*(r_mc(i+1)**2-r_mc(i)**2)*m_B(i,1)
         ddd(2)=ddd(2)+0.5*pie*(r_mc(i+1)**2-r_mc(i)**2)*m_B(i,2)
         ddd(3)=ddd(3)+0.5*pie*(r_mc(i+1)**2-r_mc(i)**2)*m_B(i,3)

       else if (i==NH) then

         ccc(1,1)=ccc(1,1)+0.5*0.*(1.+m_A(i,1,1))
         ccc(1,2)=ccc(1,2)+0.5*0.*(0.+m_A(i,1,2))
         ccc(1,3)=ccc(1,3)+0.5*0.*(0.+m_A(i,1,3))
         ccc(2,1)=ccc(2,1)+0.5*0.*(0.+m_A(i,2,1))
         ccc(2,2)=ccc(2,2)+0.5*0.*(1.+m_A(i,2,2))
         ccc(2,3)=ccc(2,3)+0.5*0.*(0.+m_A(i,2,3))
         ccc(3,1)=ccc(3,1)+0.5*0.*(0.+m_A(i,3,1))
         ccc(3,2)=ccc(3,2)+0.5*0.*(0.+m_A(i,3,2))
         ccc(3,3)=ccc(3,3)+0.5*0.*(1.+m_A(i,3,3))

         ddd(1)=ddd(1)+0.5*0.*m_B(i,1)
         ddd(2)=ddd(2)+0.5*0.*m_B(i,2)
         ddd(3)=ddd(3)+0.5*0.*m_B(i,3)

       else
         
         ccc(1,1)=ccc(1,1)+0.5*pie*(r_mc(i+1)**2-r_mc(i)**2)*(m_A(i-1,1,1)+m_A(i,1,1))
         ccc(1,2)=ccc(1,2)+0.5*pie*(r_mc(i+1)**2-r_mc(i)**2)*(m_A(i-1,1,2)+m_A(i,1,2))
         ccc(1,3)=ccc(1,3)+0.5*pie*(r_mc(i+1)**2-r_mc(i)**2)*(m_A(i-1,1,3)+m_A(i,1,3))
         ccc(2,1)=ccc(2,1)+0.5*pie*(r_mc(i+1)**2-r_mc(i)**2)*(m_A(i-1,2,1)+m_A(i,2,1))
         ccc(2,2)=ccc(2,2)+0.5*pie*(r_mc(i+1)**2-r_mc(i)**2)*(m_A(i-1,2,2)+m_A(i,2,2))
         ccc(2,3)=ccc(2,3)+0.5*pie*(r_mc(i+1)**2-r_mc(i)**2)*(m_A(i-1,2,3)+m_A(i,2,3))
         ccc(3,1)=ccc(3,1)+0.5*pie*(r_mc(i+1)**2-r_mc(i)**2)*(m_A(i-1,3,1)+m_A(i,3,1))
         ccc(3,2)=ccc(3,2)+0.5*pie*(r_mc(i+1)**2-r_mc(i)**2)*(m_A(i-1,3,2)+m_A(i,3,2))
         ccc(3,3)=ccc(3,3)+0.5*pie*(r_mc(i+1)**2-r_mc(i)**2)*(m_A(i-1,3,3)+m_A(i,3,3))

         ddd(1)=ddd(1)+0.5*pie*(r_mc(i+1)**2-r_mc(i)**2)*(m_B(i-1,1)+m_B(i,1))
         ddd(2)=ddd(2)+0.5*pie*(r_mc(i+1)**2-r_mc(i)**2)*(m_B(i-1,2)+m_B(i,2))
         ddd(3)=ddd(3)+0.5*pie*(r_mc(i+1)**2-r_mc(i)**2)*(m_B(i-1,3)+m_B(i,3))

       end if


     end do


     !!!!!!!!!!!!!!!!!!!!!z��������������!!!!!!!!!!!!!!!!!!!!(�ⲿ���Թ�)

       !write(*,*) CCC(1,1),CCC(1,2),CCC(1,3)
       !write(*,*) CCC(2,1),CCC(2,2),CCC(2,3)
       !write(*,*) CCC(3,1),CCC(3,2),CCC(3,3)
       !write(*,*) DDD(1),DDD(2),DDD(3)
       !pause





     !!!!!!!!!!!!!!!!!!����һ���ڵ��ϵ�Ӧ������!!!!!!!!!!!!!!!!!!!!!

     AAA(1,1)=1.
     AAA(1,2)=-1.
     AAA(1,3)=0.
     AAA(2,1)=m_A(NH+N_CLAD-1,1,1)
     AAA(2,2)=m_A(NH+N_CLAD-1,1,2)
     AAA(2,3)=m_A(NH+N_CLAD-1,1,3)
     AAA(3,1)=ccc(3,1)
     AAA(3,2)=ccc(3,2)
     AAA(3,3)=ccc(3,3)

     BBB(1)=0.
     BBB(2)=-m_B(NH+N_CLAD-1,1)-p_coolant*1.D6
     BBB(3)=-ddd(3)-(p_coolant-PRESS_INTER)*1.d6*PIE*R(NH+2)**2

     call agaus2(AAA,bbb,XXX)

     stress_fuel(1,1)=XXX(1)
     stress_fuel(1,2)=XXX(2)
     stress_fuel(1,3)=XXX(3)

     stress_equ(1)=(0.5*((stress_fuel(1,1)-stress_fuel(1,2))**2+(stress_fuel(1,1)-stress_fuel(1,3))**2+(stress_fuel(1,2)-stress_fuel(1,3))**2))**0.5

     !write(*,*)  XXX
     !write(*,*)  stress_equ(1)/1.d6
     !pause
     !!!!!!!!!!!!!!!!!!����һ���ڵ��ϵ�Ӧ������!!!!!!!!!!!!!!!!!!!!!
     


     !!!!!!!!!!!!!!!!!!��������ڵ��ϵ�Ӧ������!!!!!!!!!!!!!!!!!!!!!
    
     do i=2,nh+N_CLAD 

       if (i<=nh) then
       stress_fuel(i,1)=m_A(i-1,1,1)*stress_fuel(1,1)+m_A(i-1,1,2)*stress_fuel(1,2)+m_A(i-1,1,3)*stress_fuel(1,3)+m_b(i-1,1)
       stress_fuel(i,2)=m_A(i-1,2,1)*stress_fuel(1,1)+m_A(i-1,2,2)*stress_fuel(1,2)+m_A(i-1,2,3)*stress_fuel(1,3)+m_b(i-1,2)
       stress_fuel(i,3)=m_A(i-1,3,1)*stress_fuel(1,1)+m_A(i-1,3,2)*stress_fuel(1,2)+m_A(i-1,3,3)*stress_fuel(1,3)+m_b(i-1,3)

       stress_equ(i)=(0.5*((stress_fuel(i,1)-stress_fuel(i,2))**2+(stress_fuel(i,1)-stress_fuel(i,3))**2+(stress_fuel(i,2)-stress_fuel(i,3))**2))**0.5
       
       !write(*,*) i
       !write(*,*) stress_fuel(i,1),stress_fuel(i,2),stress_fuel(i,3)
       !write(*,*) stress_equ(i)/1.d6 

       else

       stress_clad(i-nh,1)=m_A(i-1,1,1)*stress_fuel(1,1)+m_A(i-1,1,2)*stress_fuel(1,2)+m_A(i-1,1,3)*stress_fuel(1,3)+m_b(i-1,1)
       stress_clad(i-nh,2)=m_A(i-1,2,1)*stress_fuel(1,1)+m_A(i-1,2,2)*stress_fuel(1,2)+m_A(i-1,2,3)*stress_fuel(1,3)+m_b(i-1,2)
       stress_clad(i-nh,3)=m_A(i-1,3,1)*stress_fuel(1,1)+m_A(i-1,3,2)*stress_fuel(1,2)+m_A(i-1,3,3)*stress_fuel(1,3)+m_b(i-1,3)

       stress_equ(i)=(0.5*((stress_clad(i-nh,1)-stress_clad(i-nh,2))**2+(stress_clad(i-nh,1)-stress_clad(i-nh,3))**2&
       &+(stress_clad(i-nh,2)-stress_clad(i-nh,3))**2))**0.5

       !write(*,*) i
       !write(*,*) stress_clad(i-nh,1),stress_clad(i-nh,2),stress_clad(i-nh,3)
       !write(*,*) stress_equ(i)/1.d6 

       end if

       !pause
     end do

     !!!!!!!!!!!!!!!!!!��������ڵ��ϵ�Ӧ������!!!!!!!!!!!!!!!!!!!!!
     !pause

     !force_total=0
     !do i=1,nh+N_CLAD-1
     !if (i<=nh-1)  then
     !force_total=force_total+0.5*pie*(r_mc(i+1)**2-r_mc(i)**2)*(stress_fuel(i,3)+stress_fuel(i+1,3))
     !else if (i==nh) then
     !force_total=force_total+0.
     !else
     !force_total=force_total+0.5*pie*(r_mc(i+1)**2-r_mc(i)**2)*(stress_clad(i-nh,3)+stress_clad(i-nh+1,3))
     !end if
     !end do
     !write(*,*)  -(p_coolant-PRESS_INTER)*1.d6*PIE*R(NH+2)**2,force_total
     !pause




     !!!!!!!!!!!!!!!!!����Ӧ����Ӧ��!!!!!!!!!!!!!!!!!

     do i=1,nh

        strain_fuel(i,1)=1./m_elastic(i)*(stress_fuel(i,1)-Poisson_ration(i)*(stress_fuel(i,2)+stress_fuel(i,3)))+strain_intrinsic(i,1)
        strain_fuel(i,2)=1./m_elastic(i)*(stress_fuel(i,2)-Poisson_ration(i)*(stress_fuel(i,1)+stress_fuel(i,3)))+strain_intrinsic(i,2)
        strain_fuel(i,3)=1./m_elastic(i)*(stress_fuel(i,3)-Poisson_ration(i)*(stress_fuel(i,2)+stress_fuel(i,1)))+strain_intrinsic(i,3)

     end do


     do i=1,N_CLAD

        strain_clad(i,1)=1./m_elastic(nh+i)*(stress_clad(i,1)-Poisson_ration(nh+i)*(stress_clad(i,2)+stress_clad(i,3)))+strain_intrinsic(nh+i,1)
        strain_clad(i,2)=1./m_elastic(nh+i)*(stress_clad(i,2)-Poisson_ration(nh+i)*(stress_clad(i,1)+stress_clad(i,3)))+strain_intrinsic(nh+i,2)
        strain_clad(i,3)=1./m_elastic(nh+i)*(stress_clad(i,3)-Poisson_ration(nh+i)*(stress_clad(i,2)+stress_clad(i,1)))+strain_intrinsic(nh+i,3)
        !if (i==1) then
        !write(*,*) stress_clad(i,1),stress_clad(i,2),stress_clad(i,3)
        !write(*,*) m_elastic(nh+i),Poisson_ration(nh+i)
        !write(*,*) strain_clad(i,1),strain_clad(i,2),strain_clad(i,3)
        !end if
        !pause
     end do

     !write(*,*)  (1.+strain_fuel(nh,2))*r_mc(nh),(1.+strain_clad(1,2))*r_mc(nh+1)
     !write(*,*)  strain_fuel(nh,3)-strain_z_lasttime111(1),strain_clad(1,3)-strain_z_lasttime111(2)
     !write(*,*)  strain_z_lasttime111(2),strain_z_lasttime111(1),PRESS_INTER

     !!!!!!!!!!!!!!!!!����Ӧ����Ӧ��!!!!!!!!!!!!!!!!!(�������)

     !pause




     !!!!!!!!!!!!!!!!!!�жϽ����������!!!!!!!!!!!!!!!!

     factor=0
     do i=1,nh+N_CLAD

       if (stress_equ(i)>yield_stress(i)) then

       !write(*,*) i,"��������"
       !pause

         !d_strain_plastic(i,1)=0.01
         !d_strain_plastic(i,2)=-0.005
         !d_strain_plastic(i,3)=-0.005

         factor=factor+1

       else
         
         !d_strain_plastic(i,1)=0.
         !d_strain_plastic(i,2)=0.
         !d_strain_plastic(i,3)=0.

         factor=factor+0

       end if


     end do

     !!!!!!!!!!!!!!!!!!�жϽ����������!!!!!!!!!!!!!!!!






if (factor==0) then




       !exit  !û���κ�һ��������Ա��Σ�������do while ѭ��

       do i=1,nh+N_CLAD
       yield_stress111(i)=yield_stress(i)
       strain_plastic111(i,1)=strain_plastic(i,1)
       strain_plastic111(i,2)=strain_plastic(i,2)
       strain_plastic111(i,3)=strain_plastic(i,3)
       end do



else 
       

        !write(*,*)  "��������ѭ��"


        !pause

     !!!!!!!!!!!!!!!!!!��ֵ���Ա���!!!!!!!!!!!!!!!!

     !factor=0
     do i=1,nh+N_CLAD

       if (stress_equ(i)>yield_stress(i)) then

       !write(*,*) i,"��������"
       !pause

         d_strain_plastic(i,1)=-0.0001
         d_strain_plastic(i,2)=0.
         d_strain_plastic(i,3)=+0.0001

         !factor=factor+1

       else
         
         d_strain_plastic(i,1)=0.
         d_strain_plastic(i,2)=0.
         d_strain_plastic(i,3)=0.

         !factor=factor+0

       end if

       !write(*,*) i
       !write(*,*) d_strain_plastic(i,1),d_strain_plastic(i,2),d_strain_plastic(i,3)

     end do

     !!!!!!!!!!!!!!!!!!��ֵ���Ա���!!!!!!!!!!!!!!!!







      
       

     
     factor2=1.
     factor3=1.





do while ((factor3>0.05))    !do while�����Ե���

     d_strain_plastic_js=d_strain_plastic


     do i=1,nh
     stress_fuel_js(i,1)=stress_fuel(i,1)
     stress_fuel_js(i,2)=stress_fuel(i,2)
     stress_fuel_js(i,3)=stress_fuel(i,3)
     !write(*,*) stress_fuel_js(i,1),stress_fuel_js(i,2),stress_fuel_js(i,3)
     end do

     do i=1,N_CLAD
     stress_clad_js(i,1)=stress_clad(i,1)
     stress_clad_js(i,2)=stress_clad(i,2)
     stress_clad_js(i,3)=stress_clad(i,3)
     !write(*,*) stress_clad_js(i,1),stress_clad_js(i,2),stress_clad_js(i,3)
     end do


       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !!!!!!!!!!!!!���µ��Լ���!!!!!!!!!!!!!!
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      do i=1,NH

          strain_intrinsic(i,1)=l_strain(I)+strain_plastic(i,1)+d_strain_plastic(i,1)
          strain_intrinsic(i,2)=l_strain(I)+strain_plastic(i,2)+d_strain_plastic(i,2)
          strain_intrinsic(i,3)=l_strain(I)+strain_plastic(i,3)+d_strain_plastic(i,3)

      end do
     

      do i=1,N_CLAD

             temperature_cladding=T_NOW(nh+1)+(i-1)*(T_NOW(nh+2)-T_NOW(nh+1))/(N_CLAD-1)
             if (temperature_cladding<1083.) then
                strain_thermal_cladding_r=1.*(4.95d-6*temperature_cladding-1.485d-3)
                strain_thermal_cladding_z=1.*(1.26d-5*temperature_cladding-3.78d-3)
                strain_thermal_cladding_c=strain_thermal_cladding_r
             else if (temperature_cladding<1244.) then
                strain_thermal_cladding_r=(2.77763+1.09822*cos(3.14159265*(temperature_cladding-1083.)/161.))/1000.
                strain_thermal_cladding_z=(2.77763+1.09822*cos(3.14159265*(temperature_cladding-1083.)/161.))/1000.
                strain_thermal_cladding_c=strain_thermal_cladding_r
             else if (temperature_cladding<2098.) then
                strain_thermal_cladding_r=9.7d-6*temperature_cladding-1.04d-2
                strain_thermal_cladding_z=9.76d-6*temperature_cladding-4.4d-3
                strain_thermal_cladding_c=strain_thermal_cladding_r
             else
                strain_thermal_cladding_r=9.7d-6*2098.-1.04d-2
                strain_thermal_cladding_z=9.76d-6*2098.-4.4d-3
                strain_thermal_cladding_c=strain_thermal_cladding_r
             end if

          strain_intrinsic(NH+i,1)=strain_thermal_cladding_r+strain_plastic(NH+i,1)+d_strain_plastic(NH+i,1)+strain_creep(i,1)+d_strain_creep(i,1)
          strain_intrinsic(NH+i,2)=strain_thermal_cladding_c+strain_plastic(NH+i,2)+d_strain_plastic(NH+i,2)+strain_creep(i,2)+d_strain_creep(i,2)
          strain_intrinsic(NH+i,3)=strain_thermal_cladding_z+strain_plastic(NH+i,3)+d_strain_plastic(NH+i,3)+strain_creep(i,3)+d_strain_creep(i,3)
          !write(*,*) strain_intrinsic(NH+i,1),strain_intrinsic(NH+i,2),strain_intrinsic(NH+i,3)
          !write(*,*) strain_thermal_cladding_r,strain_thermal_cladding_c,strain_thermal_cladding_z
      end do





 

     !!!!!!!!!!!!!!!!!!!!!�������ƾ������!!!!!!!!!!!!!!!!!!!!
     do i=1,NH-1

       m_l(i,1,1)=1.-(r_mc(i+1)-r_mc(i))/r_mc(i+1)
       m_l(i,1,2)=(r_mc(i+1)-r_mc(i))/r_mc(i+1)
       m_l(i,1,3)=0.
       m_m(i,1)=0.

       m_l(i,2,1)=(r_mc(i+1)-r_mc(i))/r_mc(i+1)+m_elastic(i+1)/(1.-Poisson_ration(i+1))&
       &*(Poisson_ration(i+1)/m_elastic(i+1)-Poisson_ration(i)/m_elastic(i))
       m_l(i,2,2)=1.-(r_mc(i+1)-r_mc(i))/r_mc(i+1)+m_elastic(i+1)/(1.-Poisson_ration(i+1)**2)&
       &*((Poisson_ration(i+1)*Poisson_ration(i+1)/m_elastic(i+1)-Poisson_ration(i+1)*Poisson_ration(i)/m_elastic(i))&
       &-(1./m_elastic(i+1)-1./m_elastic(i)))
       m_l(i,2,3)=m_elastic(i+1)/(1.-Poisson_ration(i+1)**2)*((Poisson_ration(i+1)/m_elastic(i+1)-Poisson_ration(i)/m_elastic(i))&
       &-(Poisson_ration(i+1)/m_elastic(i+1)-Poisson_ration(i+1)/m_elastic(i)))
       m_m(i,2)=m_elastic(i+1)/(1.-Poisson_ration(i+1)**2)*(((r_mc(i+1)-r_mc(i))/r_mc(i+1)*(strain_intrinsic(i+1,1)-strain_intrinsic(i+1,2)))&
       &-(strain_intrinsic(i+1,2)-strain_intrinsic(i,2))-Poisson_ration(i+1)*(strain_intrinsic(i+1,3)-strain_intrinsic(i,3)))

       m_l(i,3,1)=m_elastic(i+1)/(1.-Poisson_ration(i+1))*(Poisson_ration(i+1)/m_elastic(i+1)-Poisson_ration(i)/m_elastic(i))
       m_l(i,3,2)=m_elastic(i+1)/(1.-Poisson_ration(i+1)**2)*((Poisson_ration(i+1)/m_elastic(i+1)-Poisson_ration(i)/m_elastic(i))&
       &-(Poisson_ration(i+1)/m_elastic(i+1)-Poisson_ration(i+1)/m_elastic(i)))
       m_l(i,3,3)=1.+m_elastic(i+1)/(1.-Poisson_ration(i+1)**2)&
       &*((Poisson_ration(i+1)*Poisson_ration(i+1)/m_elastic(i+1)-Poisson_ration(i+1)*Poisson_ration(i)/m_elastic(i))&
       &-(1./m_elastic(i+1)-1./m_elastic(i)))
       m_m(i,3)=m_elastic(i+1)/(1.-Poisson_ration(i+1)**2)*(((r_mc(i+1)-r_mc(i))/r_mc(i+1)*&
       &Poisson_ration(i+1)*(strain_intrinsic(i+1,1)-strain_intrinsic(i+1,2)))&
       &-Poisson_ration(i+1)*(strain_intrinsic(i+1,2)-strain_intrinsic(i,2))-(strain_intrinsic(i+1,3)-strain_intrinsic(i,3)))

       !write(*,*) i
       !write(*,*) m_l(i,1,1),m_l(i,1,2),m_l(i,1,3)
       !write(*,*) m_l(i,2,1),m_l(i,2,2),m_l(i,2,3)
       !write(*,*) m_l(i,3,1),m_l(i,3,2),m_l(i,3,3)
       !write(*,*) m_m(i,1),m_m(i,2),m_m(i,3)
       !pause
     end do







       m_l(NH,1,1)=1.
       m_l(NH,1,2)=0.
       m_l(NH,1,3)=0.
       m_m(NH,1)=0.

       m_l(NH,2,1)=m_elastic(nh+1)*(r_mc(nh+1)*(Poisson_ration(nh+1)**2)/m_elastic(nh+1)+r_mc(nh+1)*Poisson_ration(nh+1)/m_elastic(nh+1)&
                  &-r_mc(nh+1)*Poisson_ration(nh+1)*Poisson_ration(nh)/m_elastic(nh)-r_mc(nh)*Poisson_ration(nh)/m_elastic(nh))&
                  &/(r_mc(nh+1)*(1.-Poisson_ration(nh+1)**2))
       m_l(NH,2,2)=m_elastic(nh+1)*(r_mc(nh)/m_elastic(nh)-r_mc(nh+1)*Poisson_ration(nh+1)*Poisson_ration(nh)/m_elastic(nh))&
                  &/(r_mc(nh+1)*(1.-Poisson_ration(nh+1)**2))
       m_l(NH,2,3)=m_elastic(nh+1)*(r_mc(nh+1)*Poisson_ration(nh+1)/m_elastic(nh)-r_mc(nh)*Poisson_ration(nh)/m_elastic(nh))&
                  &/(r_mc(nh+1)*(1.-Poisson_ration(nh+1)**2))
       m_m(NH,2)=m_elastic(nh+1)*(r_mc(nh+1)*Poisson_ration(nh+1)*((strain_intrinsic(nh,3)-strain_intrinsic(nh+1,3))&
                 &+(strain_z_lasttime111(2)-strain_z_lasttime111(1)))+(r_mc(nh)*strain_intrinsic(nh,2)-r_mc(nh+1)*strain_intrinsic(nh+1,2))&
                 &-(r_mc(nh+1)-r_mc(nh)))/(r_mc(nh+1)*(1.-Poisson_ration(nh+1)**2))

       m_l(NH,3,1)=m_elastic(nh+1)*(r_mc(nh+1)*Poisson_ration(nh+1)/m_elastic(nh+1)+r_mc(nh+1)*(Poisson_ration(nh+1)**2)/m_elastic(nh+1)&
                  &-r_mc(nh+1)*Poisson_ration(nh)/m_elastic(nh)-r_mc(nh)*Poisson_ration(nh+1)*Poisson_ration(nh)/m_elastic(nh))&
                  &/(r_mc(nh+1)*(1.-Poisson_ration(nh+1)**2))
       m_l(NH,3,2)=m_elastic(nh+1)*(r_mc(nh)*Poisson_ration(nh+1)/m_elastic(nh)-r_mc(nh+1)*Poisson_ration(nh)/m_elastic(nh))&
                  &/(r_mc(nh+1)*(1.-Poisson_ration(nh+1)**2))
       m_l(NH,3,3)=m_elastic(nh+1)*(r_mc(nh+1)/m_elastic(nh)-r_mc(nh)*Poisson_ration(nh+1)*Poisson_ration(nh)/m_elastic(nh))&
                  &/(r_mc(nh+1)*(1.-Poisson_ration(nh+1)**2))
       m_m(NH,3)=m_elastic(nh+1)*(r_mc(nh+1)*((strain_intrinsic(nh,3)-strain_intrinsic(nh+1,3))&
                 &+(strain_z_lasttime111(2)-strain_z_lasttime111(1)))+Poisson_ration(nh+1)*((r_mc(nh)*strain_intrinsic(nh,2)&
                 &-r_mc(nh+1)*strain_intrinsic(nh+1,2))-(r_mc(nh+1)-r_mc(nh))))/(r_mc(nh+1)*(1.-Poisson_ration(nh+1)**2))








     do i=NH+1,NH+N_CLAD-1           !о��Ͱ���δ�Ӵ�ʱ����NH��������Ҫ����Ϊо��Ͱ���֮���޹���

       m_l(i,1,1)=1.-(r_mc(i+1)-r_mc(i))/r_mc(i+1)
       m_l(i,1,2)=(r_mc(i+1)-r_mc(i))/r_mc(i+1)
       m_l(i,1,3)=0.
       m_m(i,1)=0.

       m_l(i,2,1)=(r_mc(i+1)-r_mc(i))/r_mc(i+1)+m_elastic(i+1)/(1.-Poisson_ration(i+1))&
       &*(Poisson_ration(i+1)/m_elastic(i+1)-Poisson_ration(i)/m_elastic(i))
       m_l(i,2,2)=1.-(r_mc(i+1)-r_mc(i))/r_mc(i+1)+m_elastic(i+1)/(1.-Poisson_ration(i+1)**2)&
       &*((Poisson_ration(i+1)*Poisson_ration(i+1)/m_elastic(i+1)-Poisson_ration(i+1)*Poisson_ration(i)/m_elastic(i))&
       &-(1./m_elastic(i+1)-1./m_elastic(i)))
       m_l(i,2,3)=m_elastic(i+1)/(1.-Poisson_ration(i+1)**2)*((Poisson_ration(i+1)/m_elastic(i+1)-Poisson_ration(i)/m_elastic(i))&
       &-(Poisson_ration(i+1)/m_elastic(i+1)-Poisson_ration(i+1)/m_elastic(i)))
       m_m(i,2)=m_elastic(i+1)/(1.-Poisson_ration(i+1)**2)*(((r_mc(i+1)-r_mc(i))/r_mc(i+1)*(strain_intrinsic(i+1,1)-strain_intrinsic(i+1,2)))&
       &-(strain_intrinsic(i+1,2)-strain_intrinsic(i,2))-Poisson_ration(i+1)*(strain_intrinsic(i+1,3)-strain_intrinsic(i,3)))

       m_l(i,3,1)=m_elastic(i+1)/(1.-Poisson_ration(i+1))*(Poisson_ration(i+1)/m_elastic(i+1)-Poisson_ration(i)/m_elastic(i))
       m_l(i,3,2)=m_elastic(i+1)/(1.-Poisson_ration(i+1)**2)*((Poisson_ration(i+1)/m_elastic(i+1)-Poisson_ration(i)/m_elastic(i))&
       &-(Poisson_ration(i+1)/m_elastic(i+1)-Poisson_ration(i+1)/m_elastic(i)))
       m_l(i,3,3)=1.+m_elastic(i+1)/(1.-Poisson_ration(i+1)**2)&
       &*((Poisson_ration(i+1)*Poisson_ration(i+1)/m_elastic(i+1)-Poisson_ration(i+1)*Poisson_ration(i)/m_elastic(i))&
       &-(1./m_elastic(i+1)-1./m_elastic(i)))
       m_m(i,3)=m_elastic(i+1)/(1.-Poisson_ration(i+1)**2)*(((r_mc(i+1)-r_mc(i))/r_mc(i+1)*&
       &Poisson_ration(i+1)*(strain_intrinsic(i+1,1)-strain_intrinsic(i+1,2)))&
       &-Poisson_ration(i+1)*(strain_intrinsic(i+1,2)-strain_intrinsic(i,2))-(strain_intrinsic(i+1,3)-strain_intrinsic(i,3)))

       !write(*,*) i
       !write(*,*) m_l(i,1,1),m_l(i,1,2),m_l(i,1,3)
       !write(*,*) m_l(i,2,1),m_l(i,2,2),m_l(i,2,3)
       !write(*,*) m_l(i,3,1),m_l(i,3,2),m_l(i,3,3)
       !write(*,*) m_m(i,1),m_m(i,2),m_m(i,3)
       !pause
     end do
















     !!!!!!!!!!!!!!!!!!!!!�������ƾ������!!!!!!!!!!!!!!!!!!!!
   

     !!!!!!!!!!!!!!!!!!!!!ת�ƾ������!!!!!!!!!!!!!!!!!!!!
     do i=1,NH+N_CLAD-1

       if (i==1) then

         m_A(i,1,1)=m_l(i,1,1)
         m_A(i,1,2)=m_l(i,1,2)
         m_A(i,1,3)=m_l(i,1,3)
         m_A(i,2,1)=m_l(i,2,1)
         m_A(i,2,2)=m_l(i,2,2)
         m_A(i,2,3)=m_l(i,2,3)
         m_A(i,3,1)=m_l(i,3,1)
         m_A(i,3,2)=m_l(i,3,2)
         m_A(i,3,3)=m_l(i,3,3)

         m_B(i,1)=m_m(i,1)
         m_B(i,2)=m_m(i,2)
         m_B(i,3)=m_m(i,3)

       else
         

         AA(1,1)=m_l(i,1,1)
         AA(1,2)=m_l(i,1,2)
         AA(1,3)=m_l(i,1,3)
         AA(2,1)=m_l(i,2,1)
         AA(2,2)=m_l(i,2,2)
         AA(2,3)=m_l(i,2,3)
         AA(3,1)=m_l(i,3,1)
         AA(3,2)=m_l(i,3,2)
         AA(3,3)=m_l(i,3,3)

         BB(1,1)=m_A(i-1,1,1)
         BB(1,2)=m_A(i-1,1,2)
         BB(1,3)=m_A(i-1,1,3)
         BB(2,1)=m_A(i-1,2,1)
         BB(2,2)=m_A(i-1,2,2)
         BB(2,3)=m_A(i-1,2,3)
         BB(3,1)=m_A(i-1,3,1)
         BB(3,2)=m_A(i-1,3,2)
         BB(3,3)=m_A(i-1,3,3)


         CC(1)=m_m(i,1)
         CC(2)=m_m(i,2)
         CC(3)=m_m(i,3)

         DD(1)=m_B(i-1,1)
         DD(2)=m_B(i-1,2)
         DD(3)=m_B(i-1,3)

         call matrix_multiply(AA,BB,XX)
         call series_multiply(AA,DD,YY)
         YY=YY+CC

         m_A(i,1,1)=XX(1,1)
         m_A(i,1,2)=XX(1,2)
         m_A(i,1,3)=XX(1,3)
         m_A(i,2,1)=XX(2,1)
         m_A(i,2,2)=XX(2,2)
         m_A(i,2,3)=XX(2,3)
         m_A(i,3,1)=XX(3,1)
         m_A(i,3,2)=XX(3,2)
         m_A(i,3,3)=XX(3,3)

         m_B(i,1)=YY(1)
         m_B(i,2)=YY(2)
         m_B(i,3)=YY(3)

       end if

       !write(*,*) i
       !write(*,*) m_l(i,1,1),m_l(i,1,2),m_l(i,1,3)
       !write(*,*) m_l(i,2,1),m_l(i,2,2),m_l(i,2,3)
       !write(*,*) m_l(i,3,1),m_l(i,3,2),m_l(i,3,3)
       !write(*,*) m_m(i,1),m_m(i,2),m_m(i,3)

       !write(*,*) m_A(i,1,1),m_A(i,1,2),m_l(i,1,3)
       !write(*,*) m_A(i,2,1),m_A(i,2,2),m_l(i,2,3)
       !write(*,*) m_A(i,3,1),m_A(i,3,2),m_l(i,3,3)
       !write(*,*) m_B(i,1),m_B(i,2),m_B(i,3)
       !pause
     


     end do
     !!!!!!!!!!!!!!!!!!!!!ת�ƾ������!!!!!!!!!!!!!!!!!!!!
     !PAUSE



     !!!!!!!!!!!!!!!!!!!!!z��������������!!!!!!!!!!!!!!!!!!!!


     CCC=0.
     DDD=0.
     do i=1,NH+N_CLAD-1
       
       if (i==1) then
         
         ccc(1,1)=ccc(1,1)+0.5*pie*(r_mc(i+1)**2-r_mc(i)**2)*(1.+m_A(i,1,1))
         ccc(1,2)=ccc(1,2)+0.5*pie*(r_mc(i+1)**2-r_mc(i)**2)*(0.+m_A(i,1,2))
         ccc(1,3)=ccc(1,3)+0.5*pie*(r_mc(i+1)**2-r_mc(i)**2)*(0.+m_A(i,1,3))
         ccc(2,1)=ccc(2,1)+0.5*pie*(r_mc(i+1)**2-r_mc(i)**2)*(0.+m_A(i,2,1))
         ccc(2,2)=ccc(2,2)+0.5*pie*(r_mc(i+1)**2-r_mc(i)**2)*(1.+m_A(i,2,2))
         ccc(2,3)=ccc(2,3)+0.5*pie*(r_mc(i+1)**2-r_mc(i)**2)*(0.+m_A(i,2,3))
         ccc(3,1)=ccc(3,1)+0.5*pie*(r_mc(i+1)**2-r_mc(i)**2)*(0.+m_A(i,3,1))
         ccc(3,2)=ccc(3,2)+0.5*pie*(r_mc(i+1)**2-r_mc(i)**2)*(0.+m_A(i,3,2))
         ccc(3,3)=ccc(3,3)+0.5*pie*(r_mc(i+1)**2-r_mc(i)**2)*(1.+m_A(i,3,3))

         ddd(1)=ddd(1)+0.5*pie*(r_mc(i+1)**2-r_mc(i)**2)*m_B(i,1)
         ddd(2)=ddd(2)+0.5*pie*(r_mc(i+1)**2-r_mc(i)**2)*m_B(i,2)
         ddd(3)=ddd(3)+0.5*pie*(r_mc(i+1)**2-r_mc(i)**2)*m_B(i,3)

       else if (i==NH) then

         ccc(1,1)=ccc(1,1)+0.5*0.*(1.+m_A(i,1,1))
         ccc(1,2)=ccc(1,2)+0.5*0.*(0.+m_A(i,1,2))
         ccc(1,3)=ccc(1,3)+0.5*0.*(0.+m_A(i,1,3))
         ccc(2,1)=ccc(2,1)+0.5*0.*(0.+m_A(i,2,1))
         ccc(2,2)=ccc(2,2)+0.5*0.*(1.+m_A(i,2,2))
         ccc(2,3)=ccc(2,3)+0.5*0.*(0.+m_A(i,2,3))
         ccc(3,1)=ccc(3,1)+0.5*0.*(0.+m_A(i,3,1))
         ccc(3,2)=ccc(3,2)+0.5*0.*(0.+m_A(i,3,2))
         ccc(3,3)=ccc(3,3)+0.5*0.*(1.+m_A(i,3,3))

         ddd(1)=ddd(1)+0.5*0.*m_B(i,1)
         ddd(2)=ddd(2)+0.5*0.*m_B(i,2)
         ddd(3)=ddd(3)+0.5*0.*m_B(i,3)

       else
         
         ccc(1,1)=ccc(1,1)+0.5*pie*(r_mc(i+1)**2-r_mc(i)**2)*(m_A(i-1,1,1)+m_A(i,1,1))
         ccc(1,2)=ccc(1,2)+0.5*pie*(r_mc(i+1)**2-r_mc(i)**2)*(m_A(i-1,1,2)+m_A(i,1,2))
         ccc(1,3)=ccc(1,3)+0.5*pie*(r_mc(i+1)**2-r_mc(i)**2)*(m_A(i-1,1,3)+m_A(i,1,3))
         ccc(2,1)=ccc(2,1)+0.5*pie*(r_mc(i+1)**2-r_mc(i)**2)*(m_A(i-1,2,1)+m_A(i,2,1))
         ccc(2,2)=ccc(2,2)+0.5*pie*(r_mc(i+1)**2-r_mc(i)**2)*(m_A(i-1,2,2)+m_A(i,2,2))
         ccc(2,3)=ccc(2,3)+0.5*pie*(r_mc(i+1)**2-r_mc(i)**2)*(m_A(i-1,2,3)+m_A(i,2,3))
         ccc(3,1)=ccc(3,1)+0.5*pie*(r_mc(i+1)**2-r_mc(i)**2)*(m_A(i-1,3,1)+m_A(i,3,1))
         ccc(3,2)=ccc(3,2)+0.5*pie*(r_mc(i+1)**2-r_mc(i)**2)*(m_A(i-1,3,2)+m_A(i,3,2))
         ccc(3,3)=ccc(3,3)+0.5*pie*(r_mc(i+1)**2-r_mc(i)**2)*(m_A(i-1,3,3)+m_A(i,3,3))

         ddd(1)=ddd(1)+0.5*pie*(r_mc(i+1)**2-r_mc(i)**2)*(m_B(i-1,1)+m_B(i,1))
         ddd(2)=ddd(2)+0.5*pie*(r_mc(i+1)**2-r_mc(i)**2)*(m_B(i-1,2)+m_B(i,2))
         ddd(3)=ddd(3)+0.5*pie*(r_mc(i+1)**2-r_mc(i)**2)*(m_B(i-1,3)+m_B(i,3))

       end if


     end do


     !!!!!!!!!!!!!!!!!!!!!z��������������!!!!!!!!!!!!!!!!!!!!(�ⲿ���Թ�)

       !write(*,*) CCC(1,1),CCC(1,2),CCC(1,3)
       !write(*,*) CCC(2,1),CCC(2,2),CCC(2,3)
       !write(*,*) CCC(3,1),CCC(3,2),CCC(3,3)
       !write(*,*) DDD(1),DDD(2),DDD(3)
       !pause





     !!!!!!!!!!!!!!!!!!����һ���ڵ��ϵ�Ӧ������!!!!!!!!!!!!!!!!!!!!!

     AAA(1,1)=1.
     AAA(1,2)=-1.
     AAA(1,3)=0.
     AAA(2,1)=m_A(NH+N_CLAD-1,1,1)
     AAA(2,2)=m_A(NH+N_CLAD-1,1,2)
     AAA(2,3)=m_A(NH+N_CLAD-1,1,3)
     AAA(3,1)=ccc(3,1)
     AAA(3,2)=ccc(3,2)
     AAA(3,3)=ccc(3,3)

     BBB(1)=0.
     BBB(2)=-m_B(NH+N_CLAD-1,1)-p_coolant*1.D6
     BBB(3)=-ddd(3)-(p_coolant-PRESS_INTER)*1.d6*PIE*R(NH+2)**2

     call agaus2(AAA,bbb,XXX)

     stress_fuel(1,1)=XXX(1)
     stress_fuel(1,2)=XXX(2)
     stress_fuel(1,3)=XXX(3)

     stress_equ(1)=(0.5*((stress_fuel(1,1)-stress_fuel(1,2))**2+(stress_fuel(1,1)-stress_fuel(1,3))**2+(stress_fuel(1,2)-stress_fuel(1,3))**2))**0.5

     !write(*,*)  XXX
     !write(*,*)  stress_equ(1)/1.d6
     !pause
     !!!!!!!!!!!!!!!!!!����һ���ڵ��ϵ�Ӧ������!!!!!!!!!!!!!!!!!!!!!
     


     !!!!!!!!!!!!!!!!!!��������ڵ��ϵ�Ӧ������!!!!!!!!!!!!!!!!!!!!!
    
     do i=2,nh+N_CLAD 

       if (i<=nh) then
       stress_fuel(i,1)=m_A(i-1,1,1)*stress_fuel(1,1)+m_A(i-1,1,2)*stress_fuel(1,2)+m_A(i-1,1,3)*stress_fuel(1,3)+m_b(i-1,1)
       stress_fuel(i,2)=m_A(i-1,2,1)*stress_fuel(1,1)+m_A(i-1,2,2)*stress_fuel(1,2)+m_A(i-1,2,3)*stress_fuel(1,3)+m_b(i-1,2)
       stress_fuel(i,3)=m_A(i-1,3,1)*stress_fuel(1,1)+m_A(i-1,3,2)*stress_fuel(1,2)+m_A(i-1,3,3)*stress_fuel(1,3)+m_b(i-1,3)

       stress_equ(i)=(0.5*((stress_fuel(i,1)-stress_fuel(i,2))**2+(stress_fuel(i,1)-stress_fuel(i,3))**2+(stress_fuel(i,2)-stress_fuel(i,3))**2))**0.5
       
       !write(*,*) i
       !write(*,*) stress_fuel(i,1),stress_fuel(i,2),stress_fuel(i,3)
       !write(*,*) stress_equ(i)/1.d6 

       else

       stress_clad(i-nh,1)=m_A(i-1,1,1)*stress_fuel(1,1)+m_A(i-1,1,2)*stress_fuel(1,2)+m_A(i-1,1,3)*stress_fuel(1,3)+m_b(i-1,1)
       stress_clad(i-nh,2)=m_A(i-1,2,1)*stress_fuel(1,1)+m_A(i-1,2,2)*stress_fuel(1,2)+m_A(i-1,2,3)*stress_fuel(1,3)+m_b(i-1,2)
       stress_clad(i-nh,3)=m_A(i-1,3,1)*stress_fuel(1,1)+m_A(i-1,3,2)*stress_fuel(1,2)+m_A(i-1,3,3)*stress_fuel(1,3)+m_b(i-1,3)

       stress_equ(i)=(0.5*((stress_clad(i-nh,1)-stress_clad(i-nh,2))**2+(stress_clad(i-nh,1)-stress_clad(i-nh,3))**2&
       &+(stress_clad(i-nh,2)-stress_clad(i-nh,3))**2))**0.5

       !write(*,*) i
       !write(*,*) stress_clad(i-nh,1),stress_clad(i-nh,2),stress_clad(i-nh,3)
       !write(*,*) stress_equ(i)/1.d6 

       end if

       !pause
     end do

     !!!!!!!!!!!!!!!!!!��������ڵ��ϵ�Ӧ������!!!!!!!!!!!!!!!!!!!!!
     !pause

     !force_total=0
     !do i=1,nh+N_CLAD-1
     !if (i<=nh-1)  then
     !force_total=force_total+0.5*pie*(r_mc(i+1)**2-r_mc(i)**2)*(stress_fuel(i,3)+stress_fuel(i+1,3))
     !else if (i==nh) then
     !force_total=force_total+0.
     !else
     !force_total=force_total+0.5*pie*(r_mc(i+1)**2-r_mc(i)**2)*(stress_clad(i-nh,3)+stress_clad(i-nh+1,3))
     !end if
     !end do
     !write(*,*)  -(p_coolant-PRESS_INTER)*1.d6*PIE*R(NH+2)**2,force_total
     !pause




     !!!!!!!!!!!!!!!!!����Ӧ����Ӧ��!!!!!!!!!!!!!!!!!

     do i=1,nh

        strain_fuel(i,1)=1./m_elastic(i)*(stress_fuel(i,1)-Poisson_ration(i)*(stress_fuel(i,2)+stress_fuel(i,3)))+strain_intrinsic(i,1)
        strain_fuel(i,2)=1./m_elastic(i)*(stress_fuel(i,2)-Poisson_ration(i)*(stress_fuel(i,1)+stress_fuel(i,3)))+strain_intrinsic(i,2)
        strain_fuel(i,3)=1./m_elastic(i)*(stress_fuel(i,3)-Poisson_ration(i)*(stress_fuel(i,2)+stress_fuel(i,1)))+strain_intrinsic(i,3)

     end do


     do i=1,N_CLAD

        strain_clad(i,1)=1./m_elastic(nh+i)*(stress_clad(i,1)-Poisson_ration(nh+i)*(stress_clad(i,2)+stress_clad(i,3)))+strain_intrinsic(nh+i,1)
        strain_clad(i,2)=1./m_elastic(nh+i)*(stress_clad(i,2)-Poisson_ration(nh+i)*(stress_clad(i,1)+stress_clad(i,3)))+strain_intrinsic(nh+i,2)
        strain_clad(i,3)=1./m_elastic(nh+i)*(stress_clad(i,3)-Poisson_ration(nh+i)*(stress_clad(i,2)+stress_clad(i,1)))+strain_intrinsic(nh+i,3)
        !if (i==1) then
        !write(*,*) stress_clad(i,1),stress_clad(i,2),stress_clad(i,3)
        !write(*,*) m_elastic(nh+i),Poisson_ration(nh+i)
        !write(*,*) strain_clad(i,1),strain_clad(i,2),strain_clad(i,3)
        !end if
        !pause
     end do

     !write(*,*)  (1.+strain_fuel(nh,2))*r_mc(nh),(1.+strain_clad(1,2))*r_mc(nh+1)
     !write(*,*)  strain_fuel(nh,3)-strain_z_lasttime111(1),strain_clad(1,3)-strain_z_lasttime111(2)
     !write(*,*)  strain_z_lasttime111(2),strain_z_lasttime111(1),PRESS_INTER

     !!!!!!!!!!!!!!!!!����Ӧ����Ӧ��!!!!!!!!!!!!!!!!!(�������)



       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !!!!!!!!!!!!!���µ��Լ���!!!!!!!!!!!!!!
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!






     !!!!!!!!!!!!!!!!!!���¸�ֵ���Ա���!!!!!!!!!!!!!!!!

     !factor=0
     do i=1,nh+n_clad

       if ((stress_equ(i)>yield_stress(i)).and.(abs(d_strain_plastic(i,1))+abs(d_strain_plastic(i,2))+abs(d_strain_plastic(i,3))==0.)) then


           d_strain_plastic(i,1)=-0.0001
           d_strain_plastic(i,2)=0.
           d_strain_plastic(i,3)=+0.0001

           !write(*,*) i,"���ĸ�ֵ"
           !write(*,*) d_strain_plastic(i,1),d_strain_plastic(i,2),d_strain_plastic(i,3)

       end if





     end do

     !!!!!!!!!!!!!!!!!!���¸�ֵ���Ա���!!!!!!!!!!!!!!!!







  



       factor2=0.
       factor3=0.
       !!!!!!!!!!!!!!!!!!���¼������Է���!!!!!!!!!!!!!!!!!!!
       do i=1,nh+n_clad

       if (i<=nh) then

         if (stress_equ(i)>yield_stress(i)) then
         !write(*,*)  i,yield_stress(i),stress_equ(i)/1.d6

           d_strain_plastic_equ(i)=2.**0.5/3.*((d_strain_plastic(i,1)-d_strain_plastic(i,2))**2+&
           &(d_strain_plastic(i,1)-d_strain_plastic(i,3))**2+(d_strain_plastic(i,2)-d_strain_plastic(i,2))**2)**0.5

           stress_curve(i)=yield_stress(i)+d_strain_plastic_equ(i)*m_plastic(i)

           s1=stress_fuel(i,1)-(stress_fuel(i,1)+stress_fuel(i,2)+stress_fuel(i,3))/3.  !pa
           s2=stress_fuel(i,2)-(stress_fuel(i,1)+stress_fuel(i,2)+stress_fuel(i,3))/3.  !pa
           s3=stress_fuel(i,3)-(stress_fuel(i,1)+stress_fuel(i,2)+stress_fuel(i,3))/3.  !pa


           d_strain_plastic(i,1)=3./2.*d_strain_plastic_equ(i)/stress_curve(i)*s1
           d_strain_plastic(i,2)=3./2.*d_strain_plastic_equ(i)/stress_curve(i)*s2
           d_strain_plastic(i,3)=3./2.*d_strain_plastic_equ(i)/stress_curve(i)*s3

           factor2=factor2+abs((d_strain_plastic(i,1)-d_strain_plastic_js(i,1))/d_strain_plastic(i,1))&
           &+abs((d_strain_plastic(i,2)-d_strain_plastic_js(i,2))/d_strain_plastic(i,2))&                    
           &+abs((d_strain_plastic(i,3)-d_strain_plastic_js(i,3))/d_strain_plastic(i,3))

           !write(*,*) d_strain_plastic(i,1),d_strain_plastic(i,2),d_strain_plastic(i,3)
           !write(*,*) d_strain_plastic_js(i,1),d_strain_plastic_js(i,2),d_strain_plastic_js(i,3)


           yield_stress111(i)=stress_equ(i)
           strain_plastic111(i,1)=strain_plastic(i,1)+d_strain_plastic(i,1)
           strain_plastic111(i,2)=strain_plastic(i,2)+d_strain_plastic(i,2)
           strain_plastic111(i,3)=strain_plastic(i,3)+d_strain_plastic(i,3)

           !write(*,*) i,stress_curve(i)/1.d6,stress_equ(i)/1.d6
           !write(*,*) d_strain_plastic(i,1),d_strain_plastic(i,2),d_strain_plastic(i,3)
           
 
         else
  

           d_strain_plastic(i,1)=0.5*d_strain_plastic(i,1)
           d_strain_plastic(i,2)=0.5*d_strain_plastic(i,2)      !���Ϊ�㣬����0.5����ν
           d_strain_plastic(i,3)=0.5*d_strain_plastic(i,3)      !�����Ϊ�㣬��˵����ʼ�ĳ�ֵ����ˣ���С

           yield_stress111(i)=yield_stress(i)
           strain_plastic111(i,1)=strain_plastic(i,1)
           strain_plastic111(i,2)=strain_plastic(i,2)
           strain_plastic111(i,3)=strain_plastic(i,3)

           factor2=factor2+0.

           !write(*,*) i,stress_curve(i)/1.d6,stress_equ(i)/1.d6
           !write(*,*) d_strain_plastic(i,1),d_strain_plastic(i,2),d_strain_plastic(i,3)
 
         end if

        

         factor3=factor3+abs((stress_fuel(i,1)-stress_fuel_js(i,1))/stress_fuel(i,1))&
           &+abs((stress_fuel(i,2)-stress_fuel_js(i,2))/stress_fuel(i,2))&
           &+abs((stress_fuel(i,3)-stress_fuel_js(i,3))/stress_fuel(i,3))
        
        !write(*,*) stress_fuel(i,1),stress_fuel(i,2),stress_fuel(i,3)
        !write(*,*) stress_fuel_js(i,1),stress_fuel_js(i,2),stress_fuel_js(i,3)
        !write(*,*) i,factor2,factor3
        !write(*,*) d_strain_plastic(i,1),d_strain_plastic(i,2),d_strain_plastic(i,3)
        !pause




       else





         if (stress_equ(i)>yield_stress(i)) then
         !write(*,*)  i,yield_stress(i),stress_equ(i)/1.d6

           d_strain_plastic_equ(i)=2.**0.5/3.*((d_strain_plastic(i,1)-d_strain_plastic(i,2))**2+&
           &(d_strain_plastic(i,1)-d_strain_plastic(i,3))**2+(d_strain_plastic(i,2)-d_strain_plastic(i,2))**2)**0.5

           stress_curve(i)=yield_stress(i)+d_strain_plastic_equ(i)*m_plastic(i)

           s1=stress_clad(i-nh,1)-(stress_clad(i-nh,1)+stress_clad(i-nh,2)+stress_clad(i-nh,3))/3.  !pa
           s2=stress_clad(i-nh,2)-(stress_clad(i-nh,1)+stress_clad(i-nh,2)+stress_clad(i-nh,3))/3.  !pa
           s3=stress_clad(i-nh,3)-(stress_clad(i-nh,1)+stress_clad(i-nh,2)+stress_clad(i-nh,3))/3.  !pa


           d_strain_plastic(i,1)=3./2.*d_strain_plastic_equ(i)/stress_curve(i)*s1
           d_strain_plastic(i,2)=3./2.*d_strain_plastic_equ(i)/stress_curve(i)*s2
           d_strain_plastic(i,3)=3./2.*d_strain_plastic_equ(i)/stress_curve(i)*s3

           factor2=factor2+abs((d_strain_plastic(i,1)-d_strain_plastic_js(i,1))/d_strain_plastic(i,1))&
           &+abs((d_strain_plastic(i,2)-d_strain_plastic_js(i,2))/d_strain_plastic(i,2))&                    
           &+abs((d_strain_plastic(i,3)-d_strain_plastic_js(i,3))/d_strain_plastic(i,3))

           !write(*,*) d_strain_plastic(i,1),d_strain_plastic(i,2),d_strain_plastic(i,3)
           !write(*,*) d_strain_plastic_js(i,1),d_strain_plastic_js(i,2),d_strain_plastic_js(i,3)


           yield_stress111(i)=stress_equ(i)
           strain_plastic111(i,1)=strain_plastic(i,1)+d_strain_plastic(i,1)
           strain_plastic111(i,2)=strain_plastic(i,2)+d_strain_plastic(i,2)
           strain_plastic111(i,3)=strain_plastic(i,3)+d_strain_plastic(i,3)

           !write(*,*) i,stress_curve(i)/1.d6,stress_equ(i)/1.d6
           !write(*,*) d_strain_plastic(i,1),d_strain_plastic(i,2),d_strain_plastic(i,3)
           
 
         else
  

           d_strain_plastic(i,1)=0.5*d_strain_plastic(i,1)
           d_strain_plastic(i,2)=0.5*d_strain_plastic(i,2)      !���Ϊ�㣬����0.5����ν
           d_strain_plastic(i,3)=0.5*d_strain_plastic(i,3)      !�����Ϊ�㣬��˵����ʼ�ĳ�ֵ����ˣ���С

           yield_stress111(i)=yield_stress(i)
           strain_plastic111(i,1)=strain_plastic(i,1)
           strain_plastic111(i,2)=strain_plastic(i,2)
           strain_plastic111(i,3)=strain_plastic(i,3)

           factor2=factor2+0.

           !write(*,*) i,stress_curve(i)/1.d6,stress_equ(i)/1.d6
           !write(*,*) d_strain_plastic(i,1),d_strain_plastic(i,2),d_strain_plastic(i,3)
 
         end if   
         
               
           factor3=factor3+abs((stress_clad(i-nh,1)-stress_clad_js(i-nh,1))/stress_clad(i-nh,1))&
           &+abs((stress_clad(i-nh,2)-stress_clad_js(i-nh,2))/stress_clad(i-nh,2))&
           &+abs((stress_clad(i-nh,3)-stress_clad_js(i-nh,3))/stress_clad(i-nh,3))


       end if   !��о�黹�ǰ���




       end do
       !!!!!!!!!!!!!!!!!!���¼������Է���!!!!!!!!!!!!!!!!!!!


        !write(*,*) factor2,factor3
        !write(*,*) 
        !pause

    end do !do while�����Ե���




end if    !�нڵ��������


deng=strain_fuel(11,2)






     UR_PRA(1)=R(1)
     DO i=2,nh
     UR_PRA(i)=(1.+strain_fuel(i,2))*R(i)           !�ڶ�����������
     end do
     UR_PRA(nh+1)=(1.+strain_clad(1,2))*R(nh+1) 
     UR_PRA(nh+2)=(1.+strain_clad(n_clad,2))*R(nh+n_clad) 
   
     !write(*,*) UR_PRA(11),R(11),UR_PRA(11)-R(11)
     

     d_length_p=strain_fuel(11,3)*DH                          !��ʵ�κ�һ�������������Ӧ����һ����
     d_length_c=strain_clad(1,3)*DH
     !write(*,*) 

     !pause



     strain_z_lasttime111(1)=strain_fuel(nh,3)   !z��������һ���ڵ�
     strain_z_lasttime111(2)=strain_clad(1,3)    !z��������һ���ڵ�

     P_contact=stress_clad(1,1)/1.d6

 !UR_PRA=r

 !write(*,*)  UR_PRA
     !write(*,*)  P_contact,PRESS_INTER,UR_PRA(12)-UR_PRA(11)



     





 end if      !�ж��Ƿ����Ӵ��׶�






     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!�Ӵ������о�����ϼ���!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!�Ӵ������о�����ϼ���!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!�Ӵ������о�����ϼ���!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




!real(8)::time_total,time_increment,Q_LINE,n_flux
!real(8)::strain_creep(N_CLAD,3),strain_creep111(N_CLAD,3)

!real(8)::creep_rate_th,creep_rate_irr
!real(8)::creep_rate_add,creep_rate_total
!real(8)::creep_saturated_primary,creep_total

!real(8)::d_strain_creep(N_CLAD,3),d_strain_creep_js(N_CLAD,3)
!integer:i2,i3,i4,i5,i6,i7
!real(8)::factor44




     !stress_clad(i,1/2/3)/1.d6
     !tem_cladding(i)
     !stress_equ(i+nh)/1.d6
     !m_elastic(nh+i)






     do i3=1,n_clad
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!�����μ���!!!!!!!!!!!!!!!!!     !����ѹ���߽Ӵ�ѹ��û�б������±�����һ�µ�
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
        creep_rate_th=5.47d8*m_elastic(nh+i3)/1.d6/tem_cladding(i3)*((sinh(320.*stress_equ(i3+nh)/m_elastic(nh+i3)))**3.5)&
        &*exp(-201./0.008314/tem_cladding(i3))

        if (tem_cladding(i3)<570.) then
        creep_rate_irr=1.87d-24*n_flux**0.85*stress_equ(i3+nh)*0.7944/1.d6
        else if (tem_cladding(i3)<625.) then
        creep_rate_irr=1.87d-24*n_flux**0.85*stress_equ(i3+nh)*(-3.18562+0.00699132*tem_cladding(i3))/1.d6
        else
        creep_rate_irr=1.87d-24*n_flux**0.85*stress_equ(i3+nh)*1.184/1.d6
        end if
        !write(*,*)  creep_rate_th,creep_rate_irr
        !pause

        creep_rate_add=creep_rate_th+creep_rate_irr
        creep_saturated_primary=0.0216*(creep_rate_add**0.109)*(2.-tanh(35500.*creep_rate_add))**(-2.05)
        creep_total=creep_saturated_primary*(1.-exp(-52.*(creep_rate_add*time_total)**0.5))+creep_rate_add*time_total
        creep_rate_total=52.*creep_saturated_primary*(creep_rate_add**0.5)/2./(time_total**0.5)*&
        &exp(-52.*(creep_rate_add*time_total)**0.5)+creep_rate_add
        !write(*,*)  creep_saturated_primary,creep_total,creep_rate_total
        !pause

        s1=stress_clad(i3,1)-(stress_clad(i3,1)+stress_clad(i3,2)+stress_clad(i3,3))/3.  !pa
        s2=stress_clad(i3,2)-(stress_clad(i3,1)+stress_clad(i3,2)+stress_clad(i3,3))/3.
        s3=stress_clad(i3,3)-(stress_clad(i3,1)+stress_clad(i3,2)+stress_clad(i3,3))/3.

        d_strain_creep(i3,1)=1.5*creep_rate_total*time_increment*s1/stress_equ(i3+nh)
        d_strain_creep(i3,2)=1.5*creep_rate_total*time_increment*s2/stress_equ(i3+nh)
        d_strain_creep(i3,3)=1.5*creep_rate_total*time_increment*s3/stress_equ(i3+nh)
        !write(*,*) d_creep_cladding_r,d_creep_cladding_c,d_creep_cladding_z
        !pause
            
        strain_creep111(i3,1)=strain_creep(i3,1)+d_strain_creep(i3,1)
        strain_creep111(i3,2)=strain_creep(i3,2)+d_strain_creep(i3,2)
        strain_creep111(i3,3)=strain_creep(i3,3)+d_strain_creep(i3,3)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!�����μ���!!!!!!!!!!!!!!!!!     !����ѹ���߽Ӵ�ѹ��û�б������±�����һ�µ�
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !write(*,*) i3
        !write(*,*) d_strain_creep(i3,1),d_strain_creep111(i3,2),d_strain_creep111(i3,3)
     end do

     !pause

     factor44=maxval(abs((d_strain_creep-d_strain_creep_js)/d_strain_creep))
     !write(*,*)  d_strain_creep(1,2),creep_rate_add,creep_rate_total
     !pause

end do    !!!!!!!!!!!!!!�����εĴ�ѭ��










stress_cladding_z=stress_clad(1,3)
strain_cladding_z=strain_clad(1,3)








END SUBROUTINE MECHMODEL_unrigid
