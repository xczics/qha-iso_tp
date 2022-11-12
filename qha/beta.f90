subroutine GetLnbeta (Hfreq,Lfreq,Weights,Nq,Nfreq,An,T,NT,Lnbeta) bind(C, name="GetLnbeta")
    use iso_c_binding
    implicit none
    integer(c_int), intent(in),value :: Nq, Nfreq,An,NT
    real(c_double), intent(in) :: Lfreq(Nq,Nfreq), Hfreq(Nq,Nfreq), Weights(Nq)
    real(c_double), intent(in) :: T(NT)
    real(c_double), intent(inout) :: Lnbeta(NT)
    real(c_double),PARAMETER:: H=6.62606957D-34,K=1.3806488D-23,C=2.99792458D8*100
    real(c_double),PARAMETER:: CK=1.0D12
    !H是普朗克常数 K是玻尔兹曼常熟 C是光束，单位是cm/s CK是THZ单位换算成Hz单位的系数
    real(c_double)::MSH(Nq,NT),MSL(Nq,NT),NU(Nfreq),NDP(Nfreq)    
!NU 用来记录各Q点的每一个频率,U是频率，单位THz
!MS 用来记录各Q点经过计算出来的配分值,NDP 用来记录各Q点每一个频率的贡献
    integer(c_int)::L,I
    !write(*,*)Lfreq!,Hfreq,Weights,Nq,Nfreq,An,T,NT,Lnbeta
    DO L = 1,NT
        DO I = 1,Nq
            NU=Lfreq(I,:)*CK*H/K/T(L)
            NDP=NU*EXP(-NU/2.0)/(1-EXP(-NU))
            MSL(I,L)=PRODUCT(NDP)
        END DO
    END DO
    !write(*,*) MSL   
    DO L = 1,NT
        DO I = 1,Nq
            NU=Hfreq(I,:)*CK*H/K/T(L)
            NDP=NU*EXP(-NU/2.0)/(1-EXP(-NU))
            MSH(I,L)=PRODUCT(NDP)
        END DO
    END DO
    DO L= 1,NT
        Lnbeta(L)=1000*DOT_PRODUCT(Weights,LOG(MSH(:,L)/MSL(:,L)))/SUM(Weights)/real(An)
    ENDDO

end subroutine GetLnbeta   

