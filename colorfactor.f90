module colorfactor
    implicit none
    save
    double precision,dimension(3,3,10) :: colorfactor1F
    double precision,dimension(3,3) :: colorfactor3t5
    double precision,dimension(2,3) :: colorfactorqtqg
    double precision,dimension(2,3) :: colorfactorgtqq

end module colorfactor

subroutine ColorMatrix()
    use colorfactor

    implicit none

    colorfactor1F(1,1,:)=(/-2.0D0/3.0D0,-2.0D0/3.0D0,0.0D0,0.0D0,-2.0D0/3.0D0,0.0D0,0.0D0,0.0D0,0.0D0,4.0D0/3.0D0/)
    colorfactor1F(1,2,:)=(/1.0D0/6.0D0/sqrt(2.0D0),0.0D0,-1.0D0/3.0D0/sqrt(2.0D0),-1.0D0/3.0D0/sqrt(2.0D0),0.0D0,&
        & -1.0D0/3.0D0/sqrt(2.0D0),-1.0D0/3.0D0/sqrt(2.0D0),sqrt(2.0D0)/3.0D0,sqrt(2.0D0)/3.0D0,0.0D0/)
    colorfactor1F(1,3,:)=(/0.0D0,0.0D0,1.0D0/sqrt(6.0D0),1.0D0/sqrt(6.0D0),0.0D0,-1.0D0/sqrt(6.0D0),1.0D0/sqrt(6.0D0)&
        & ,0.0d0,0.0d0,0.0d0/)
    colorfactor1F(2,2,:)=(/-2.0D0/3.0D0,1.0D0/12.0D0,-7.0D0/12.0D0,1.0D0/6.0D0,1.0D0/12.0D0,-7.0D0/12.0D0,1.0D0/6.0D0&
        & ,-1.0d0/3.0D0,7.0d0/6.0D0,-1.0d0/6.0D0/)
    colorfactor1F(2,3,:)=(/0.0D0,-sqrt(3.0D0)/4.0D0,-1.0D0/4.0D0/sqrt(3.0D0),-1.0D0/sqrt(3.0D0),sqrt(3.0D0)/4.0D0,&
        & 1.0D0/4.0D0/sqrt(3.0D0),1.0D0/sqrt(3.0D0),0.0D0,0.0D0,0.0D0/)
    colorfactor1F(3,3,:)=(/1.0D0/3.0D0,-5.0D0/12.0D0,-5.0D0/12.0D0,5.0D0/6.0D0,-5.0D0/12.0D0,&
        & -5.0D0/12.0D0,5.0D0/6.0D0,-2.0D0/3.0D0,-1.0D0/6.0D0,-1.0D0/6.0D0/)
    colorfactor1F(2,1,:)=(/1.0D0/6.0D0/sqrt(2.0D0),0.0D0,-1.0D0/3.0D0/sqrt(2.0D0),-1.0D0/3.0D0/sqrt(2.0D0),0.0D0,&
        & -1.0D0/3.0D0/sqrt(2.0D0),-1.0D0/3.0D0/sqrt(2.0D0),sqrt(2.0D0)/3.0D0,sqrt(2.0D0)/3.0D0,0.0D0/)
    colorfactor1F(3,1,:)=(/0.0D0,0.0D0,1.0D0/sqrt(6.0D0),1.0D0/sqrt(6.0D0),0.0D0,-1.0D0/sqrt(6.0D0),1.0D0/sqrt(6.0D0)&
        & ,0.0d0,0.0d0,0.0d0/)
    colorfactor1F(3,2,:)=(/0.0D0,-sqrt(3.0D0)/4.0D0,-1.0D0/4.0D0/sqrt(3.0D0),-1.0D0/sqrt(3.0D0),sqrt(3.0D0)/4.0D0,&
        & 1.0D0/4.0D0/sqrt(3.0D0),1.0D0/sqrt(3.0D0),0.0D0,0.0D0,0.0D0/)

    colorfactor3t5(1,:)=(/0.0d0,0.0d0,0.0d0/)
    colorfactor3t5(2,:)=(/-sqrt(6.0d0)/6.0D0,-sqrt(6.0d0)/6.0d0,sqrt(6.0d0)/3.0d0/)
    colorfactor3t5(3,:)=(/sqrt(2.0d0)/2.0d0,-sqrt(2.0d0)/2.0d0,0.0d0/)
    ! colorfactor3t5(1,:)=(/1.0d0,1.0d0,1.0d0/)
    ! colorfactor3t5(2,:)=(/1.0d0,1.0d0,1.0d0/)
    ! colorfactor3t5(3,:)=(/1.0d0,1.0d0,1.0d0/)

    colorfactorqtqg(1,:)=(/2.0d0*sqrt(3.0d0)/3.0d0,-sqrt(3.0d0)/3.0d0,-sqrt(3.0d0)/3.0d0/)
    colorfactorqtqg(2,:)=(/0.0d0,1.0d0,-1.0d0/)

    colorfactorgtqq(1,:)=(/0.0D0,sqrt(2.0D0)/2.0D0,0.0D0/)
    colorfactorgtqq(2,:)=(/0.0D0,0.0D0,sqrt(2.0D0)/2.0D0/)


    return
end subroutine ColorMatrix