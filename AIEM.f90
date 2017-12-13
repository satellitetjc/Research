!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!c
!         Advanced Integral Equation Model                  c
!                                                           c
!  Program computes emissivity                              c
!  from 3d finitely conducting surface                      c
!  ( only for single scattering )                           c
!                                                           c
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!c
!                                                           c
!  Parameters:                                              c
!                                                           c
!     er: surface relative dielectric constant              c
!     kl: normalized surface correlation length             c
!     ks: normalized surface rms height                     c
!                                                           c
!     itype: select what type of surface correlation        c
!           =1  Gaussian correlation function               c
!           =2  exponential correlation function            c
!           =3  transformed exponential correlation         c
!                                                           c
!     theta: incident angle in deg                          c
!     (phi=0 : incident azimuth angle in deg)               c
!                                                           c
!                                                           c
!***********************************************************c
!        Approximations of Fresnel reflection coeff.        c
!        -------------------------------------------        c
!                                                           c
!    irc=1: Fresnel reflection coeff. approxmated by        c
!           R(incident_angle)                               c
!    irc=2: Fresnel reflection coeff. approxmated by        c
!           R(specular_angle)                               c
!                                                           c
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!c
!                                                           c
!                                   Sept. 18, 2001          c
!                                                           c
!                                             T. D. Wu      c
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!c 
      implicit real*8 (a-h,k,o-z)
      complex*16 er,ur,stem,rvi,rhi
      integer irc,npp,nss,itype
      dimension ww11(512),zz11(512),ww22(512),zz22(512)
      dimension wt(512),zt(512)
      dimension sigma0(4)
	  logical  back,shadowing
      common /surface/itype
      common /eu/er,ur
	  common /al/ effslop
      real*4 mv,fre, temperature,sig,cl
!-- for input parameter file name
        character * 255 fname
      npp=40
      nss=80
      ur=1.0
      pi=dacos(-1.0d0)
!
!-----------------------------------------------------------c
!        Approximations of Fresnel reflection coeff.        c
!-----------------------------------------------------------c
       irc=1     ! approxmated by R(theta)
!      irc=2     ! approxmated by R(theta_specular)
       itype=1
!      itype=1     ! Gaussian correlated surface 
!      itype=2     ! Exponential correlated surface
!      itype=3     ! Transformed exponential correlated surface
!------------------------------------------c
!    set up integration limits             c
!------------------------------------------c
       aa11=0.0
       bb11=3.14
       aa22=0.0
       bb22=1.57
!-----------------------------------------c
!  generates zeros and absissa            c
!-----------------------------------------c
      call quagen(zt,wt,npp)
      as11=(bb11-aa11)/2.0
      bs11=(bb11+aa11)/2.0
      do 11 i=1,npp
       ww11(i)=wt(i)*as11
       zz11(i)=zt(i)*as11+bs11
11     continue
      call quagen(zt,wt,nss)
      as22=(bb22-aa22)/2.0
      bs22=(bb22+aa22)/2.0
      do 12 i=1,nss
       ww22(i)=wt(i)*as22
       zz22(i)=zt(i)*as22+bs22
12     continue
!
!-----------------------------------------------------------c
!             open files                                    c
!-----------------------------------------------------------c   
!      open(10,file='iem.in',status='old',access='sequential')
!-----------------------------------------------------------c
!          input   Parameters     c
!-----------------------------------------------------------c
!先装VS2010再装Intel Visual Fortran2011，打开VS2010，新建Fortran项目console1，将AIEM2003里的后缀为f90文件（这些代码就是）添加到Console1-source Files在右侧“解决方案资源管理器”     
     open(11,file='AIEMOUT_89.0.txt')
      fre=89.0                ! 10.65 19.35   22.23  37   85.5卫星的波段，下面参数一致时，可以直接改.AMSR2的波段有6.93、7.3、10.65、18.7、23.8、36.5、89.0
	  theta=55                ! 52.8 卫星的入射角度
!	  sig=2
      sig=2.5                   ! 均方根高度
      cl=7.5                    ! 相关长度
      mv=0.4                    ! 土壤水分
      err=16                    ! 介电常数虚部&实部
	  eri=9                     ! 介电常数虚部&实部
	  temperature=20            ! 地表温度
   
                                  ! 下面do之后的2008、2010、2012、2013是指代码的标识，
                                  ! 下面continue语句前也有红色的标识，表示循环调到哪个代码。
                                  ! 通常情况下只改【步长】就有很多数据了

      do 2008 sig = 0.2,3,0.4     ! [三个数据分别代表（以下相同）：起始值，终值，步长]  sig:均方根高度
      do 2010 cl = 5,30,5         ! cl:相关长度
      do 2012 mv = 0.02,0.5,0.02  ! mv:土壤水分
      do 2013 temperature=1,40,3  ! 地表温度
!       do 2008 sig = 0.8,3.1,0.5
!       do 2010 cl = 8,29,6
!       do 2012  mv = 0.08,0.46,0.06
!   	   do 2013 temperature=5,46,8

!     do 2008 sig = 2.75,3.1,0.25

!       do 2008 sig = 0.2,3.1,0.3
!       do 2010 cl = 5,31,3
!       do 2012  mv = 0.02,0.45,0.05

!       do 2008 sig = 0.25,3.1,0.4
!       do 2010 cl = 5,30,3
!       do 2012  mv = 0.02,0.45,0.02
!	   do 2013 temperature=0,60,5
!     501   read(10,*,end=1998) freq,theta,cl,sig,err,eri
!
 

      call mois2dc(mv,fre,er,temperature)
      shadowing=.false.
      back=.false.
      k=pi*fre/15.0
      kl=k*cl
      ks=k*sig
	  err=dreal(er)
	  eri=dimag(er)
!      kl=ks*5.
!      er=dcmplx(err,eri)
!	
	  if(itype.eq.1) slope=1.414*sig/cl		
!	   if(itype.eq.1) slope=1.414*0.2	
      if(itype.eq.2) slope=2.0*sig/cl
      if(itype.eq.3) slope=1.732*sig/cl
      
      effslop=0.0+slope
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~c
! reflection coefficients based on the incident angle   c
! ==> R(theta)                                          c
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~c
      cs=dcosd(theta)
      si2=dsind(theta)**2.0
        stem=cdsqrt(er*ur-si2)
      rvi=(er*cs-stem)/(er*cs+stem)
      rhi=(ur*cs-stem)/(ur*cs+stem)
!
!---------------------------------------------------------------c
!              main program                                     c
!---------------------------------------------------------------c
      rh0=cdabs(rhi)**2
      rv0=cdabs(rvi)**2
!--------------------------------------------------------------c
!          coherent component                                  c
!--------------------------------------------------------------c
       if(.not.shadowing) then
       shfct=1.0
       else
        call shadow(back,dabs(theta),dabs(theta),shfct)
      endif

!       sh2=shfct**2
       sh2=shfct
	  cohh=(cdabs(rhi)**2.)*dexp(-(2.*dcosd(theta)*ks)**2.)*sh2
      covv=(cdabs(rvi)**2.)*dexp(-(2.*dcosd(theta)*ks)**2.)*sh2
!
!--------------------------------------------------------------c
!          non-coherent component                                  c
!--------------------------------------------------------------c
      reflvv=0.0
      reflhh=0.0
      
       do 30 j=1,npp
        phis=zz11(j)*180.0/pi
        if(phis.eq.180.0) phis=179.99
       do 20 is=1,nss
        thetas=zz22(is)*180.0/pi
        if(thetas.eq.90.0) thetas=89.999
        if(thetas.eq.theta) thetas=thetas+0.001
!----------------------------------------------------------------c
      call sigma(ks,kl,theta,thetas,phis,sigma0,irc)
!----------------------------------------------------------------c
       if(.not.shadowing) then
         shfct=1.0
       else
         call shadow(back,dabs(theta),dabs(thetas),shfct)
       endif
!       sh2=shfct**2
	   sh2=shfct

!
      if(sigma0(1).gt.0.0) vv=sigma0(1)*sh2
      if(sigma0(2).gt.0.0) hh=sigma0(2)*sh2
      if(sigma0(3).gt.0.0) hv=sigma0(3)*sh2
      if(sigma0(4).gt.0.0) vh=sigma0(4)*sh2
      reflhh=reflhh+2.0*ww11(j)*ww22(is)*(hh+vh)*dsind(thetas)&
           /(4.0*pi*dcosd(theta))
      reflvv=reflvv+2.0*ww11(j)*ww22(is)*(vv+hv)*dsind(thetas)&
           /(4.0*pi*dcosd(theta))

20    continue
30    continue
!
     
       emisvv=1.0-(reflvv+covv)
       emishh=1.0-(reflhh+cohh)

!
 !     open(unit=3,file='iem.out',status='unknown',&
 !             access='sequential',form='formatted')
!      write(6,117) emisvv(1),emisvv(2),emishh(1),emishh(2)
!     &     ,sig,cl,rv0,rh0
!       write(3,117) emisvv(1),emisvv(2),emishh(1),emishh(2) &
!          ,sig,cl,rv0,rh0
!       write(11,117) temperature,mv,err,eri,emisvv,emishh
       write(11,117) temperature,mv,err,eri,emisvv,emishh,sig, cl  ! emisvv:垂直极化发射率；emishh：水平极化发射率
                                                                    ! 计算亮温  亮温=（273k+temperature）*emisvv（emishh）
   	   write(*,117) temperature,mv,err,eri,emisvv,emishh,sig,cl

117       format(8(f10.4))
     
!
2013    continue
2012	continue
2010	continue
2008	continue   


!
1998  continue
!
      stop
      end
!
!
!
!*******************************************************************c
! subroutine sigma calculates the scattering coefficients           c
!*******************************************************************c
      subroutine sigma(ks,kl,theta,thetas,phis,sigma0,irc)
      implicit real*8 (a-h,k,o-z)
      integer irc,itype,ip
      complex*16 er,ur,rv,rh,rvh,stem,steml,fvv,fhh,fvh,fhv
      complex*16 rvi,rhi,rvhi,rvl,rhl,rvhl
      complex*16 Fahh,Favv,Favh,Fahv
        complex*16 Fbhh,Fbvv,Fbvh,Fbhv
        complex*16 expkc1,expkc2,expc1,expc2,expc3,expc4
        complex*16 expkcaup,expkcadn,expkcbup,expkcbdn
        complex*16 expcauau,expcadau,expcauad,expcadad
      complex*16 expcbuau,expcbdau,expcbuad,expcbdad
        complex*16 expcaubu,expcadbu,expcaubd,expcadbd
        complex*16 expcbubu,expcbdbu,expcbubd,expcbdbd
      complex*16 qq,qqs,qqt,qqts
      complex*16 kterm(4),kcterm(4),cterm(4)
        real*8 w(30000)
        dimension sigma0(4)
        dimension tempk_LOG(30000)
        common /cont1/ si,sis,cs,css,sfs,csfs,si2,sis2,cs2,css2
        common /cont2/ ks2
        common /cont3/ iterm
        common /spectra/ w
      common /eu/er,ur
      common /surface/itype
      common /rhhvv/rh,rv,rvh
!
!    
!
      torlant=1.0d-16
!
      cs=dcosd(theta)
      css=dcosd(thetas)
      si=dsind(theta)
      sis=dsind(thetas)
      csfs=dcosd(phis)
      sfs=dsind(phis)
        cs2=cs*cs
        css2=css*css
        si2=si*si
        sis2=sis*sis
!
      ks2=ks*ks
        kl2=kl*kl
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~c
! reflection coefficients based on the incident angle   c
! ==> R(theta)                                          c
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~c
      stem=cdsqrt(er*ur-si2)
      rvi=(er*cs-stem)/(er*cs+stem)
      rhi=(ur*cs-stem)/(ur*cs+stem)
        rvhi=(rvi-rhi)/2.0
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~c
! reflection coefficients based on the specular angle   c
! ==> R(theta_specular)                                 c
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~c
      csl=dsqrt(1.0d0+cs*css-si*sis*csfs)/dsqrt(2.0d0)
        sil=dsqrt(1.0d0-csl*csl)
        steml=cdsqrt(er*ur-sil*sil)
      rvl=(er*csl-steml)/(er*csl+steml)
      rhl=(ur*csl-steml)/(ur*csl+steml)
        rvhl=(rvl-rhl)/2.0
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~c
!     Reflection coefficients rv, rh, rvh               c
!           for kirchhoff field coefficients            c
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~c
        go to (101,102), irc
101          rh=rhi
             rv=rvi
       rvh=rvhi
        go to 100
102          rh=rhl
             rv=rvl
       rvh=rvhl
100   continue
!
!---------------------------------------------c
!     kirchhoff field coefficients            c
!---------------------------------------------c
      fvv=2.0d0*rv*(si*sis-(1.0d0+cs*css)*csfs)/(css+cs)
      fhh=-2.0d0*rh*(si*sis-(1.0d0+cs*css)*csfs)/(css+cs)
        fvh=-2.0d0*rvh*sfs
        fhv=2.0d0*rvh*sfs
!
        do 109,n=1,30000
      w(n)=0.0
109   continue
!
      if(ks.gt.5.)then
        tempk= ks2*(css+cs)**2.0
      do 111,n=1,1.5*tempk+100
        fn=float(n)
        ql=kl*dsqrt((-si+sis*csfs)**2.0+(0.+sis*sfs)**2.0)
!
        if((itype.gt.4).or.(itype.lt.1)) itype=1
!
        go to (911,912,913), itype
!------ Gaussian correlated surface -----------c  
911     w(n)=kl2*dexp(-ql*ql/(4.0*fn))/(2*fn)
      go to 910
!------ exponential correlated surface -----------c 
912     w(n)=(kl/fn)**2.0*(1.0+(ql/fn)**2.0)**(-1.5)
      go to 910
!------ 1.5-power correlated surface -----------c 
913     e = 1.5 * fn - 1.0d0
        y = 1.5 * fn
        gam = ( alogam(y) ) !gamma function (1.5n)
        if ( ql .eq. 0.0d0 ) then
        w(n) = kl*kl/(3.0*fn-2.0)
        go to 910
      endif
          if ( dmod( fn, 2.0d0) .eq. 0.0d0) then
           m = 1.5 * fn - 1.0            !fn is even;
            bk = log(BESSK( m, ql))           !integer order, check OK
             else
            m = 1.5 * fn - 1.0  - 0.5    !fn is odd
           bk = log(BesselK( m, ql))          !fractional order
          endif
          out = kl * kl * ( ql / 2.0d0) ** e
          w(n)=  out * dexp(bk - gam ) 
      go to 910
910   continue
!
111   continue
!------------------------------------------c
!    compute kirchhoff term                c
!------------------------------------------c
!      sum=0.0
!       temp=1.0
        
!       tempk= ks2*(css+cs)**2.0
      tempk_LOG(1)=dlog(tempk)

      n=2
          do  2 n=2,1.5*tempk+100
!-- revised by logarithm, then the code can run.
       fn=float(n)
       tempk_LOG(n)=tempk_LOG(1)+tempk_LOG(n-1)-dlog(fn)
2     continue
      mm=n-1
        
         sumk_LOG=0
         ep_LOG=-ks2*(css+cs)**2.0
        do 3 n=1,mm
          fn=float(n)
          sumk_LOG=sumk_LOG+dexp(tempk_LOG(n)+ep_LOG)*w(n)
!       print *, sumk_LOG, w(n), dexp(tempk_LOG(n)+ep_LOG)

3     continue  

!       do 3 n=1,iterm
        
!        fn=float(n)
!          temp=temp*(ks2*(cs+css)**2.0)/fn
!      temp=temp*ep*exp(-ep)/fn
!       sum=sum+temp*w(n)
!3      continue  
!      expk=dexp(-ks2*(css+cs)**2.0)*sum
      expk=sumk_LOG
!       print *, dexp(-ks2*(css+cs)**2.0),ks2*(css+cs)**2.0
      kterm(1)=0.5*expk*cdabs(fvv)**2
      kterm(2)=0.5*expk*cdabs(fhh)**2
      kterm(3)=0.5*expk*cdabs(fhv)**2
      kterm(4)=0.5*expk*cdabs(fvh)**2
!
!------------------------------------------c
!  end of kirchhoff term computation       c
!------------------------------------------c
!
      sigma0(1)=dreal(kterm(1))
      sigma0(2)=dreal(kterm(2))
        sigma0(3)=dreal(kterm(3))
        sigma0(4)=dreal(kterm(4))
!
!----------------------------------------------------------------c
      return

        else 
!---------------------------------------------c
! compute roughness spectrum w(n)             c
!---------------------------------------------c
      iterm=1
        tempold=0.0
        temp=(ks2*(cs+css)**2.0)
        do while(dabs(temp-tempold).gt.torlant)
           tempold=temp
           iterm=iterm+1
           fiterm=float(iterm)
                   temp=tempold*(ks2*(cs+css)**2.0)/fiterm
        end do
!
      
      do 110,n=1,iterm
        fn=float(n)
        ql=kl*dsqrt((-si+sis*csfs)**2.0+(0.+sis*sfs)**2.0)
!
        if((itype.gt.4).or.(itype.lt.1)) itype=1
!
        go to (901,902,903), itype
!------ Gaussian correlated surface -----------c  
901     w(n)=kl2*dexp(-ql*ql/(4.0*fn))/(2*fn)
      go to 900
!------ exponential correlated surface -----------c 
902     w(n)=(kl/fn)**2.0*(1.0+(ql/fn)**2.0)**(-1.5)
      go to 900
!------ 1.5-power correlated surface -----------c 
903     e = 1.5 * fn - 1.0d0
        y = 1.5 * fn
        gam = ( alogam(y) ) !gamma function (1.5n)
        if ( ql .eq. 0.0d0 ) then
        w(n) = kl*kl/(3.0*fn-2.0)
        go to 900
      endif
          if ( dmod( fn, 2.0d0) .eq. 0.0d0) then
           m = 1.5 * fn - 1.0            !fn is even;
            bk = log(BESSK( m, ql))           !integer order, check OK
             else
            m = 1.5 * fn - 1.0  - 0.5    !fn is odd
           bk = log(BesselK( m, ql))          !fractional order
          endif
          out = kl * kl * ( ql / 2.0d0) ** e
          w(n)=  out * dexp(bk - gam ) 
      go to 900
900   continue
!
110   continue
!------------------------------------------c
!    compute kirchhoff term                c
!------------------------------------------c
      sum=0.0
        temp=1.0
       do 30 n=1,iterm
         fn=float(n)
          temp=temp*(ks2*(cs+css)**2.0)/fn
        sum=sum+temp*w(n)
30      continue  
      expk=dexp(-ks2*(css+cs)**2.0)*sum
      kterm(1)=0.5*expk*cdabs(fvv)**2
      kterm(2)=0.5*expk*cdabs(fhh)**2
      kterm(3)=0.5*expk*cdabs(fhv)**2
      kterm(4)=0.5*expk*cdabs(fvh)**2
!
!------------------------------------------c
!  end of kirchhoff term computation       c
!------------------------------------------c
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~c
!     Reflection coefficients rv, rh, rvh               c
!           for complementary field coefficients        c
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~c
!--------------------------------------------------------------c
       ex1=dexp(-ks2*(cs2+css2+cs*css))
       ex2=dexp(-ks2*(cs2+css2))
       qq=cs
       qqt=cdsqrt(er-si2)
       qqs=css
       qqts=cdsqrt(er-sis2)
!--------------------------------------------------------------c
!------------------------------------------c
!   compute cross term                     c
!------------------------------------------c
      expkcaup=dreal(dconjg(fvv)*Favv(-si,0.0d0,qq,qq))*expkc1(qq)&
       +dreal(dconjg(fvv)*Favv(-sis*csfs,-sis*sfs,qqs,qqs))&
       *expkc2(qqs)
      expkcadn=dreal(dconjg(fvv)*Favv(-si,0.0d0,-qq,qq))*expkc1(-qq)&
       +dreal(dconjg(fvv)*Favv(-sis*csfs,-sis*sfs,-qqs,qqs))&
       *expkc2(-qqs)
      expkcbup=dreal(dconjg(fvv)*Fbvv(-si,0.0d0,qqt,qqt))*expkc1(qqt)&
       +dreal(dconjg(fvv)*Fbvv(-sis*csfs,-sis*sfs,qqts,qqts))&
       *expkc2(qqts)
      expkcbdn=dreal(dconjg(fvv)*Fbvv(-si,0.0d0,-qqt,qqt))*expkc1(-qqt)&
       +dreal(dconjg(fvv)*Fbvv(-sis*csfs,-sis*sfs,-qqts,qqts))&
       *expkc2(-qqts)
        kcterm(1)=0.25*ex1*(expkcaup+expkcadn+expkcbup+expkcbdn)
!
      expkcaup=dreal(dconjg(fhh)*Fahh(-si,0.0d0,qq,qq))*expkc1(qq)&
       +dreal(dconjg(fhh)*Fahh(-sis*csfs,-sis*sfs,qqs,qqs))&
       *expkc2(qqs)
      expkcadn=dreal(dconjg(fhh)*Fahh(-si,0.0d0,-qq,qq))*expkc1(-qq)&
       +dreal(dconjg(fhh)*Fahh(-sis*csfs,-sis*sfs,-qqs,qqs))&
       *expkc2(-qqs)
      expkcbup=dreal(dconjg(fhh)*Fbhh(-si,0.0d0,qqt,qqt))*expkc1(qqt)&
       +dreal(dconjg(fhh)*Fbhh(-sis*csfs,-sis*sfs,qqts,qqts))&
       *expkc2(qqts)
      expkcbdn=dreal(dconjg(fhh)*Fbhh(-si,0.0d0,-qqt,qqt))*expkc1(-qqt)&
       +dreal(dconjg(fhh)*Fbhh(-sis*csfs,-sis*sfs,-qqts,qqts))&
       *expkc2(-qqts)
        kcterm(2)=0.25*ex1*(expkcaup+expkcadn+expkcbup+expkcbdn)
!
      expkcaup=dreal(dconjg(fhv)*Fahv(-si,0.0d0,qq,qq))*expkc1(qq)&
       +dreal(dconjg(fhv)*Fahv(-sis*csfs,-sis*sfs,qqs,qqs))&
       *expkc2(qqs)
      expkcadn=dreal(dconjg(fhv)*Fahv(-si,0.0d0,-qq,qq))*expkc1(-qq)&
       +dreal(dconjg(fhv)*Fahv(-sis*csfs,-sis*sfs,-qqs,qqs))&
       *expkc2(-qqs)
      expkcbup=dreal(dconjg(fhv)*Fbhv(-si,0.0d0,qqt,qqt))*expkc1(qqt)&
       +dreal(dconjg(fhv)*Fbhv(-sis*csfs,-sis*sfs,qqts,qqts))&
       *expkc2(qqts)
      expkcbdn=dreal(dconjg(fhv)*Fbhv(-si,0.0d0,-qqt,qqt))*expkc1(-qqt)&
       +dreal(dconjg(fhv)*Fbhv(-sis*csfs,-sis*sfs,-qqts,qqts))&
       *expkc2(-qqts)
        kcterm(3)=0.25*ex1*(expkcaup+expkcadn+expkcbup+expkcbdn)
!
      expkcaup=dreal(dconjg(fvh)*Favh(-si,0.0d0,qq,qq))*expkc1(qq)&
       +dreal(dconjg(fvh)*Favh(-sis*csfs,-sis*sfs,qqs,qqs))&
       *expkc2(qqs)
      expkcadn=dreal(dconjg(fvh)*Favh(-si,0.0d0,-qq,qq))*expkc1(-qq)&
       +dreal(dconjg(fvh)*Favh(-sis*csfs,-sis*sfs,-qqs,qqs))&
       *expkc2(-qqs)
      expkcbup=dreal(dconjg(fvh)*Fbvh(-si,0.0d0,qqt,qqt))*expkc1(qqt)&
       +dreal(dconjg(fvh)*Fbvh(-sis*csfs,-sis*sfs,qqts,qqts))&
       *expkc2(qqts)
      expkcbdn=dreal(dconjg(fvh)*Fbvh(-si,0.0d0,-qqt,qqt))*expkc1(-qqt)&
       +dreal(dconjg(fvh)*Fbvh(-sis*csfs,-sis*sfs,-qqts,qqts))&
       *expkc2(-qqts)
        kcterm(4)=0.25*ex1*(expkcaup+expkcadn+expkcbup+expkcbdn)
!
!--------------------------------------c
! end of computation of cross terms    c
!--------------------------------------c
!
!--------------------------------------c
! evaluate  complementary term         c
!--------------------------------------c
      expcauau=Favv(-si,0.0d0,qq,qq)&
        *dconjg(Favv(-si,0.0d0,qq,qq))*expc1(qq,qq)&
        +Favv(-si,0.0d0,qq,qq)&
        *dconjg(Favv(-sis*csfs,-sis*sfs,qqs,qqs))*expc2(qq,qqs)&
        +Favv(-sis*csfs,-sis*sfs,qqs,qqs)&
        *dconjg(Favv(-si,0.0d0,qq,qq))*expc3(qqs,qq)&
        +Favv(-sis*csfs,-sis*sfs,qqs,qqs)&
        *dconjg(Favv(-sis*csfs,-sis*sfs,qqs,qqs))*expc4(qqs,qqs)
      expcadau=Favv(-si,0.0d0,-qq,qq)&
        *dconjg(Favv(-si,0.0d0,qq,qq))*expc1(-qq,qq)&
        +Favv(-si,0.0d0,-qq,qq)&
        *dconjg(Favv(-sis*csfs,-sis*sfs,qqs,qqs))*expc2(-qq,qqs)&
        +Favv(-sis*csfs,-sis*sfs,-qqs,qqs)&
        *dconjg(Favv(-si,0.0d0,qq,qq))*expc3(-qqs,qq)&
        +Favv(-sis*csfs,-sis*sfs,-qqs,qqs)&
        *dconjg(Favv(-sis*csfs,-sis*sfs,qqs,qqs))*expc4(-qqs,qqs)
      expcauad=Favv(-si,0.0d0,qq,qq)&
        *dconjg(Favv(-si,0.0d0,-qq,qq))*expc1(qq,-qq)&
        +Favv(-si,0.0d0,qq,qq)&
        *dconjg(Favv(-sis*csfs,-sis*sfs,-qqs,qqs))*expc2(qq,-qqs)&
        +Favv(-sis*csfs,-sis*sfs,qqs,qqs)&
        *dconjg(Favv(-si,0.0d0,-qq,qq))*expc3(qqs,-qq)&
        +Favv(-sis*csfs,-sis*sfs,qqs,qqs)&
        *dconjg(Favv(-sis*csfs,-sis*sfs,-qqs,qqs))*expc4(qqs,-qqs)
      expcadad=Favv(-si,0.0d0,-qq,qq)&
        *dconjg(Favv(-si,0.0d0,-qq,qq))*expc1(-qq,-qq)&
        +Favv(-si,0.0d0,-qq,qq)&
        *dconjg(Favv(-sis*csfs,-sis*sfs,-qqs,qqs))*expc2(-qq,-qqs)&
        +Favv(-sis*csfs,-sis*sfs,-qqs,qqs)&
        *dconjg(Favv(-si,0.0d0,-qq,qq))*expc3(-qqs,-qq)&
        +Favv(-sis*csfs,-sis*sfs,-qqs,qqs)&
        *dconjg(Favv(-sis*csfs,-sis*sfs,-qqs,qqs))*expc4(-qqs,-qqs)
      expcbuau=Fbvv(-si,0.0d0,qqt,qqt)&
        *dconjg(Favv(-si,0.0d0,qq,qq))*expc1(qqt,qq)&
        +Fbvv(-si,0.0d0,qqt,qqt)&
        *dconjg(Favv(-sis*csfs,-sis*sfs,qqs,qqs))*expc2(qqt,qqs)&
        +Fbvv(-sis*csfs,-sis*sfs,qqts,qqts)&
        *dconjg(Favv(-si,0.0d0,qq,qq))*expc3(qqts,qq)&
        +Fbvv(-sis*csfs,-sis*sfs,qqts,qqts)&
        *dconjg(Favv(-sis*csfs,-sis*sfs,qqs,qqs))*expc4(qqts,qqs)
      expcbdau=Fbvv(-si,0.0d0,-qqt,qqt)&
        *dconjg(Favv(-si,0.0d0,qq,qq))*expc1(-qqt,qq)&
        +Fbvv(-si,0.0d0,-qqt,qqt)&
        *dconjg(Favv(-sis*csfs,-sis*sfs,qqs,qqs))*expc2(-qqt,qqs)&
        +Fbvv(-sis*csfs,-sis*sfs,-qqts,qqts)&
        *dconjg(Favv(-si,0.0d0,qq,qq))*expc3(-qqts,qq)&
        +Fbvv(-sis*csfs,-sis*sfs,-qqts,qqts)&
        *dconjg(Favv(-sis*csfs,-sis*sfs,qqs,qqs))*expc4(-qqts,qqs)
      expcbuad=Fbvv(-si,0.0d0,qqt,qqt)&
        *dconjg(Favv(-si,0.0d0,-qq,qq))*expc1(qqt,-qq)&
        +Fbvv(-si,0.0d0,qqt,qqt)&
        *dconjg(Favv(-sis*csfs,-sis*sfs,-qqs,qqs))*expc2(qqt,-qqs)&
        +Fbvv(-sis*csfs,-sis*sfs,qqts,qqts)&
        *dconjg(Favv(-si,0.0d0,-qq,qq))*expc3(qqts,-qq)&
        +Fbvv(-sis*csfs,-sis*sfs,qqts,qqts)&
        *dconjg(Favv(-sis*csfs,-sis*sfs,-qqs,qqs))*expc4(qqts,-qqs)
      expcbdad=Fbvv(-si,0.0d0,-qqt,qqt)&
        *dconjg(Favv(-si,0.0d0,-qq,qq))*expc1(-qqt,-qq)&
        +Fbvv(-si,0.0d0,-qqt,qqt)&
        *dconjg(Favv(-sis*csfs,-sis*sfs,-qqs,qqs))*expc2(-qqt,-qqs)&
        +Fbvv(-sis*csfs,-sis*sfs,-qqts,qqts)&
        *dconjg(Favv(-si,0.0d0,-qq,qq))*expc3(-qqts,-qq)&
        +Fbvv(-sis*csfs,-sis*sfs,-qqts,qqts)&
        *dconjg(Favv(-sis*csfs,-sis*sfs,-qqs,qqs))*expc4(-qqts,-qqs)
      expcaubu=Favv(-si,0.0d0,qq,qq)&
        *dconjg(Fbvv(-si,0.0d0,qqt,qqt))*expc1(qq,qqt)&
        +Favv(-si,0.0d0,qq,qq)&
        *dconjg(Fbvv(-sis*csfs,-sis*sfs,qqts,qqts))*expc2(qq,qqts)&
        +Favv(-sis*csfs,-sis*sfs,qqs,qqs)&
        *dconjg(Fbvv(-si,0.0d0,qqt,qqt))*expc3(qqs,qqt)&
        +Favv(-sis*csfs,-sis*sfs,qqs,qqs)&
        *dconjg(Fbvv(-sis*csfs,-sis*sfs,qqts,qqts))*expc4(qqs,qqts)
      expcadbu=Favv(-si,0.0d0,-qq,qq)&
        *dconjg(Fbvv(-si,0.0d0,qqt,qqt))*expc1(-qq,qqt)&
        +Favv(-si,0.0d0,-qq,qq)&
        *dconjg(Fbvv(-sis*csfs,-sis*sfs,qqts,qqts))*expc2(-qq,qqts)&
        +Favv(-sis*csfs,-sis*sfs,-qqs,qqs)&
        *dconjg(Fbvv(-si,0.0d0,qqt,qqt))*expc3(-qqs,qqt)&
        +Favv(-sis*csfs,-sis*sfs,-qqs,qqs)&
        *dconjg(Fbvv(-sis*csfs,-sis*sfs,qqts,qqts))*expc4(-qqs,qqts)
      expcaubd=Favv(-si,0.0d0,qq,qq)&
        *dconjg(Fbvv(-si,0.0d0,-qqt,qqt))*expc1(qq,-qqt)&
        +Favv(-si,0.0d0,qq,qq)&
        *dconjg(Fbvv(-sis*csfs,-sis*sfs,-qqts,qqts))*expc2(qq,-qqts)&
        +Favv(-sis*csfs,-sis*sfs,qqs,qqs)&
        *dconjg(Fbvv(-si,0.0d0,-qqt,qqt))*expc3(qqs,-qqt)&
        +Favv(-sis*csfs,-sis*sfs,qqs,qqs)&
        *dconjg(Fbvv(-sis*csfs,-sis*sfs,-qqts,qqts))*expc4(qqs,-qqts)
      expcadbd=Favv(-si,0.0d0,-qq,qq)&
        *dconjg(Fbvv(-si,0.0d0,-qqt,qqt))*expc1(-qq,-qqt)&
        +Favv(-si,0.0d0,-qq,qq)&
        *dconjg(Fbvv(-sis*csfs,-sis*sfs,-qqts,qqts))*expc2(-qq,-qqts)&
        +Favv(-sis*csfs,-sis*sfs,-qqs,qqs)&
        *dconjg(Fbvv(-si,0.0d0,-qqt,qqt))*expc3(-qqs,-qqt)&
        +Favv(-sis*csfs,-sis*sfs,-qqs,qqs)&
        *dconjg(Fbvv(-sis*csfs,-sis*sfs,-qqts,qqts))*expc4(-qqs,-qqts)
      expcbubu=Fbvv(-si,0.0d0,qqt,qqt)&
        *dconjg(Fbvv(-si,0.0d0,qqt,qqt))*expc1(qqt,qqt)&
        +Fbvv(-si,0.0d0,qqt,qqt)&
        *dconjg(Fbvv(-sis*csfs,-sis*sfs,qqts,qqts))*expc2(qqt,qqts)&
        +Fbvv(-sis*csfs,-sis*sfs,qqts,qqts)&
        *dconjg(Fbvv(-si,0.0d0,qqt,qqt))*expc3(qqts,qqt)&
        +Fbvv(-sis*csfs,-sis*sfs,qqts,qqts)&
        *dconjg(Fbvv(-sis*csfs,-sis*sfs,qqts,qqts))*expc4(qqts,qqts)
      expcbdbu=Fbvv(-si,0.0d0,-qqt,qqt)&
        *dconjg(Fbvv(-si,0.0d0,qqt,qqt))*expc1(-qqt,qqt)&
        +Fbvv(-si,0.0d0,-qqt,qqt)&
        *dconjg(Fbvv(-sis*csfs,-sis*sfs,qqts,qqts))*expc2(-qqt,qqts)&
        +Fbvv(-sis*csfs,-sis*sfs,-qqts,qqts)&
        *dconjg(Fbvv(-si,0.0d0,qqt,qqt))*expc3(-qqts,qqt)&
        +Fbvv(-sis*csfs,-sis*sfs,-qqts,qqts)&
        *dconjg(Fbvv(-sis*csfs,-sis*sfs,qqts,qqts))*expc4(-qqts,qqts)
      expcbubd=Fbvv(-si,0.0d0,qqt,qqt)&
        *dconjg(Fbvv(-si,0.0d0,-qqt,qqt))*expc1(qqt,-qqt)&
        +Fbvv(-si,0.0d0,qqt,qqt)&
        *dconjg(Fbvv(-sis*csfs,-sis*sfs,-qqts,qqts))*expc2(qqt,-qqts)&
        +Fbvv(-sis*csfs,-sis*sfs,qqts,qqts)&
        *dconjg(Fbvv(-si,0.0d0,-qqt,qqt))*expc3(qqts,-qqt)&
        +Fbvv(-sis*csfs,-sis*sfs,qqts,qqts)&
        *dconjg(Fbvv(-sis*csfs,-sis*sfs,-qqts,qqts))*expc4(qqts,-qqts)
      expcbdbd=Fbvv(-si,0.0d0,-qqt,qqt)&
        *dconjg(Fbvv(-si,0.0d0,-qqt,qqt))*expc1(-qqt,-qqt)&
        +Fbvv(-si,0.0d0,-qqt,qqt)&
        *dconjg(Fbvv(-sis*csfs,-sis*sfs,-qqts,qqts))*expc2(-qqt,-qqts)&
        +Fbvv(-sis*csfs,-sis*sfs,-qqts,qqts)&
        *dconjg(Fbvv(-si,0.0d0,-qqt,qqt))*expc3(-qqts,-qqt)&
        +Fbvv(-sis*csfs,-sis*sfs,-qqts,qqts)&
        *dconjg(Fbvv(-sis*csfs,-sis*sfs,-qqts,qqts))*expc4(-qqts,-qqts)
      cterm(1)=0.03125*ex2*(expcauau+expcadau+expcauad+expcadad&
                           +expcbuau+expcbdau+expcbuad+expcbdad&
                           +expcaubu+expcadbu+expcaubd+expcadbd&
                           +expcbubu+expcbdbu+expcbubd+expcbdbd)
!
      expcauau=Fahh(-si,0.0d0,qq,qq)&
        *dconjg(Fahh(-si,0.0d0,qq,qq))*expc1(qq,qq)&
        +Fahh(-si,0.0d0,qq,qq)&
        *dconjg(Fahh(-sis*csfs,-sis*sfs,qqs,qqs))*expc2(qq,qqs)&
        +Fahh(-sis*csfs,-sis*sfs,qqs,qqs)&
        *dconjg(Fahh(-si,0.0d0,qq,qq))*expc3(qqs,qq)&
        +Fahh(-sis*csfs,-sis*sfs,qqs,qqs)&
        *dconjg(Fahh(-sis*csfs,-sis*sfs,qqs,qqs))*expc4(qqs,qqs)
      expcadau=Fahh(-si,0.0d0,-qq,qq)&
        *dconjg(Fahh(-si,0.0d0,qq,qq))*expc1(-qq,qq)&
        +Fahh(-si,0.0d0,-qq,qq)&
        *dconjg(Fahh(-sis*csfs,-sis*sfs,qqs,qqs))*expc2(-qq,qqs)&
        +Fahh(-sis*csfs,-sis*sfs,-qqs,qqs)&
        *dconjg(Fahh(-si,0.0d0,qq,qq))*expc3(-qqs,qq)&
        +Fahh(-sis*csfs,-sis*sfs,-qqs,qqs)&
        *dconjg(Fahh(-sis*csfs,-sis*sfs,qqs,qqs))*expc4(-qqs,qqs)
      expcauad=Fahh(-si,0.0d0,qq,qq)&
        *dconjg(Fahh(-si,0.0d0,-qq,qq))*expc1(qq,-qq)&
        +Fahh(-si,0.0d0,qq,qq)&
        *dconjg(Fahh(-sis*csfs,-sis*sfs,-qqs,qqs))*expc2(qq,-qqs)&
        +Fahh(-sis*csfs,-sis*sfs,qqs,qqs)&
        *dconjg(Fahh(-si,0.0d0,-qq,qq))*expc3(qqs,-qq)&
        +Fahh(-sis*csfs,-sis*sfs,qqs,qqs)&
        *dconjg(Fahh(-sis*csfs,-sis*sfs,-qqs,qqs))*expc4(qqs,-qqs)
      expcadad=Fahh(-si,0.0d0,-qq,qq)&
        *dconjg(Fahh(-si,0.0d0,-qq,qq))*expc1(-qq,-qq)&
        +Fahh(-si,0.0d0,-qq,qq)&
        *dconjg(Fahh(-sis*csfs,-sis*sfs,-qqs,qqs))*expc2(-qq,-qqs)&
        +Fahh(-sis*csfs,-sis*sfs,-qqs,qqs)&
        *dconjg(Fahh(-si,0.0d0,-qq,qq))*expc3(-qqs,-qq)&
        +Fahh(-sis*csfs,-sis*sfs,-qqs,qqs)&
        *dconjg(Fahh(-sis*csfs,-sis*sfs,-qqs,qqs))*expc4(-qqs,-qqs)
      expcbuau=Fbhh(-si,0.0d0,qqt,qqt)&
        *dconjg(Fahh(-si,0.0d0,qq,qq))*expc1(qqt,qq)&
        +Fbhh(-si,0.0d0,qqt,qqt)&
        *dconjg(Fahh(-sis*csfs,-sis*sfs,qqs,qqs))*expc2(qqt,qqs)&
        +Fbhh(-sis*csfs,-sis*sfs,qqts,qqts)&
        *dconjg(Fahh(-si,0.0d0,qq,qq))*expc3(qqts,qq)&
        +Fbhh(-sis*csfs,-sis*sfs,qqts,qqts)&
        *dconjg(Fahh(-sis*csfs,-sis*sfs,qqs,qqs))*expc4(qqts,qqs)
      expcbdau=Fbhh(-si,0.0d0,-qqt,qqt)&
        *dconjg(Fahh(-si,0.0d0,qq,qq))*expc1(-qqt,qq)&
        +Fbhh(-si,0.0d0,-qqt,qqt)&
        *dconjg(Fahh(-sis*csfs,-sis*sfs,qqs,qqs))*expc2(-qqt,qqs)&
        +Fbhh(-sis*csfs,-sis*sfs,-qqts,qqts)&
        *dconjg(Fahh(-si,0.0d0,qq,qq))*expc3(-qqts,qq)&
        +Fbhh(-sis*csfs,-sis*sfs,-qqts,qqts)&
        *dconjg(Fahh(-sis*csfs,-sis*sfs,qqs,qqs))*expc4(-qqts,qqs)
      expcbuad=Fbhh(-si,0.0d0,qqt,qqt)&
        *dconjg(Fahh(-si,0.0d0,-qq,qq))*expc1(qqt,-qq)&
        +Fbhh(-si,0.0d0,qqt,qqt)&
        *dconjg(Fahh(-sis*csfs,-sis*sfs,-qqs,qqs))*expc2(qqt,-qqs)&
        +Fbhh(-sis*csfs,-sis*sfs,qqts,qqts)&
        *dconjg(Fahh(-si,0.0d0,-qq,qq))*expc3(qqts,-qq)&
        +Fbhh(-sis*csfs,-sis*sfs,qqts,qqts)&
        *dconjg(Fahh(-sis*csfs,-sis*sfs,-qqs,qqs))*expc4(qqts,-qqs)
      expcbdad=Fbhh(-si,0.0d0,-qqt,qqt)&
        *dconjg(Fahh(-si,0.0d0,-qq,qq))*expc1(-qqt,-qq)&
        +Fbhh(-si,0.0d0,-qqt,qqt)&
        *dconjg(Fahh(-sis*csfs,-sis*sfs,-qqs,qqs))*expc2(-qqt,-qqs)&
        +Fbhh(-sis*csfs,-sis*sfs,-qqts,qqts)&
        *dconjg(Fahh(-si,0.0d0,-qq,qq))*expc3(-qqts,-qq)&
        +Fbhh(-sis*csfs,-sis*sfs,-qqts,qqts)&
        *dconjg(Fahh(-sis*csfs,-sis*sfs,-qqs,qqs))*expc4(-qqts,-qqs)
      expcaubu=Fahh(-si,0.0d0,qq,qq)&
        *dconjg(Fbhh(-si,0.0d0,qqt,qqt))*expc1(qq,qqt)&
        +Fahh(-si,0.0d0,qq,qq)&
        *dconjg(Fbhh(-sis*csfs,-sis*sfs,qqts,qqts))*expc2(qq,qqts)&
        +Fahh(-sis*csfs,-sis*sfs,qqs,qqs)&
        *dconjg(Fbhh(-si,0.0d0,qqt,qqt))*expc3(qqs,qqt)&
        +Fahh(-sis*csfs,-sis*sfs,qqs,qqs)&
        *dconjg(Fbhh(-sis*csfs,-sis*sfs,qqts,qqts))*expc4(qqs,qqts)
      expcadbu=Fahh(-si,0.0d0,-qq,qq)&
        *dconjg(Fbhh(-si,0.0d0,qqt,qqt))*expc1(-qq,qqt)&
        +Fahh(-si,0.0d0,-qq,qq)&
        *dconjg(Fbhh(-sis*csfs,-sis*sfs,qqts,qqts))*expc2(-qq,qqts)&
        +Fahh(-sis*csfs,-sis*sfs,-qqs,qqs)&
        *dconjg(Fbhh(-si,0.0d0,qqt,qqt))*expc3(-qqs,qqt)&
        +Fahh(-sis*csfs,-sis*sfs,-qqs,qqs)&
        *dconjg(Fbhh(-sis*csfs,-sis*sfs,qqts,qqts))*expc4(-qqs,qqts)
      expcaubd=Fahh(-si,0.0d0,qq,qq)&
        *dconjg(Fbhh(-si,0.0d0,-qqt,qqt))*expc1(qq,-qqt)&
        +Fahh(-si,0.0d0,qq,qq)&
        *dconjg(Fbhh(-sis*csfs,-sis*sfs,-qqts,qqts))*expc2(qq,-qqts)&
        +Fahh(-sis*csfs,-sis*sfs,qqs,qqs)&
        *dconjg(Fbhh(-si,0.0d0,-qqt,qqt))*expc3(qqs,-qqt)&
        +Fahh(-sis*csfs,-sis*sfs,qqs,qqs)&
        *dconjg(Fbhh(-sis*csfs,-sis*sfs,-qqts,qqts))*expc4(qqs,-qqts)
      expcadbd=Fahh(-si,0.0d0,-qq,qq)&
        *dconjg(Fbhh(-si,0.0d0,-qqt,qqt))*expc1(-qq,-qqt)&
        +Fahh(-si,0.0d0,-qq,qq)&
        *dconjg(Fbhh(-sis*csfs,-sis*sfs,-qqts,qqts))*expc2(-qq,-qqts)&
        +Fahh(-sis*csfs,-sis*sfs,-qqs,qqs)&
        *dconjg(Fbhh(-si,0.0d0,-qqt,qqt))*expc3(-qqs,-qqt)&
        +Fahh(-sis*csfs,-sis*sfs,-qqs,qqs)&
        *dconjg(Fbhh(-sis*csfs,-sis*sfs,-qqts,qqts))*expc4(-qqs,-qqts)
      expcbubu=Fbhh(-si,0.0d0,qqt,qqt)&
        *dconjg(Fbhh(-si,0.0d0,qqt,qqt))*expc1(qqt,qqt)&
        +Fbhh(-si,0.0d0,qqt,qqt)&
        *dconjg(Fbhh(-sis*csfs,-sis*sfs,qqts,qqts))*expc2(qqt,qqts)&
        +Fbhh(-sis*csfs,-sis*sfs,qqts,qqts)&
        *dconjg(Fbhh(-si,0.0d0,qqt,qqt))*expc3(qqts,qqt)&
        +Fbhh(-sis*csfs,-sis*sfs,qqts,qqts)&
        *dconjg(Fbhh(-sis*csfs,-sis*sfs,qqts,qqts))*expc4(qqts,qqts)
      expcbdbu=Fbhh(-si,0.0d0,-qqt,qqt)&
        *dconjg(Fbhh(-si,0.0d0,qqt,qqt))*expc1(-qqt,qqt)&
        +Fbhh(-si,0.0d0,-qqt,qqt)&
        *dconjg(Fbhh(-sis*csfs,-sis*sfs,qqts,qqts))*expc2(-qqt,qqts)&
        +Fbhh(-sis*csfs,-sis*sfs,-qqts,qqts)&
        *dconjg(Fbhh(-si,0.0d0,qqt,qqt))*expc3(-qqts,qqt)&
        +Fbhh(-sis*csfs,-sis*sfs,-qqts,qqts)&
        *dconjg(Fbhh(-sis*csfs,-sis*sfs,qqts,qqts))*expc4(-qqts,qqts)
      expcbubd=Fbhh(-si,0.0d0,qqt,qqt)&
        *dconjg(Fbhh(-si,0.0d0,-qqt,qqt))*expc1(qqt,-qqt)&
        +Fbhh(-si,0.0d0,qqt,qqt)&
        *dconjg(Fbhh(-sis*csfs,-sis*sfs,-qqts,qqts))*expc2(qqt,-qqts)&
        +Fbhh(-sis*csfs,-sis*sfs,qqts,qqts)&
        *dconjg(Fbhh(-si,0.0d0,-qqt,qqt))*expc3(qqts,-qqt)&
        +Fbhh(-sis*csfs,-sis*sfs,qqts,qqts)&
        *dconjg(Fbhh(-sis*csfs,-sis*sfs,-qqts,qqts))*expc4(qqts,-qqts)
      expcbdbd=Fbhh(-si,0.0d0,-qqt,qqt)&
        *dconjg(Fbhh(-si,0.0d0,-qqt,qqt))*expc1(-qqt,-qqt)&
        +Fbhh(-si,0.0d0,-qqt,qqt)&
        *dconjg(Fbhh(-sis*csfs,-sis*sfs,-qqts,qqts))*expc2(-qqt,-qqts)&
        +Fbhh(-sis*csfs,-sis*sfs,-qqts,qqts)&
        *dconjg(Fbhh(-si,0.0d0,-qqt,qqt))*expc3(-qqts,-qqt)&
        +Fbhh(-sis*csfs,-sis*sfs,-qqts,qqts)&
        *dconjg(Fbhh(-sis*csfs,-sis*sfs,-qqts,qqts))*expc4(-qqts,-qqts)
      cterm(2)=0.03125*ex2*(expcauau+expcadau+expcauad+expcadad&
                           +expcbuau+expcbdau+expcbuad+expcbdad&
                           +expcaubu+expcadbu+expcaubd+expcadbd&
                           +expcbubu+expcbdbu+expcbubd+expcbdbd)
!
      expcauau=Fahv(-si,0.0d0,qq,qq)&
        *dconjg(Fahv(-si,0.0d0,qq,qq))*expc1(qq,qq)&
        +Fahv(-si,0.0d0,qq,qq)&
        *dconjg(Fahv(-sis*csfs,-sis*sfs,qqs,qqs))*expc2(qq,qqs)&
        +Fahv(-sis*csfs,-sis*sfs,qqs,qqs)&
        *dconjg(Fahv(-si,0.0d0,qq,qq))*expc3(qqs,qq)&
        +Fahv(-sis*csfs,-sis*sfs,qqs,qqs)&
        *dconjg(Fahv(-sis*csfs,-sis*sfs,qqs,qqs))*expc4(qqs,qqs)
      expcadau=Fahv(-si,0.0d0,-qq,qq)&
        *dconjg(Fahv(-si,0.0d0,qq,qq))*expc1(-qq,qq)&
        +Fahv(-si,0.0d0,-qq,qq)&
        *dconjg(Fahv(-sis*csfs,-sis*sfs,qqs,qqs))*expc2(-qq,qqs)&
        +Fahv(-sis*csfs,-sis*sfs,-qqs,qqs)&
        *dconjg(Fahv(-si,0.0d0,qq,qq))*expc3(-qqs,qq)&
        +Fahv(-sis*csfs,-sis*sfs,-qqs,qqs)&
        *dconjg(Fahv(-sis*csfs,-sis*sfs,qqs,qqs))*expc4(-qqs,qqs)
      expcauad=Fahv(-si,0.0d0,qq,qq)&
        *dconjg(Fahv(-si,0.0d0,-qq,qq))*expc1(qq,-qq)&
        +Fahv(-si,0.0d0,qq,qq)&
        *dconjg(Fahv(-sis*csfs,-sis*sfs,-qqs,qqs))*expc2(qq,-qqs)&
        +Fahv(-sis*csfs,-sis*sfs,qqs,qqs)&
        *dconjg(Fahv(-si,0.0d0,-qq,qq))*expc3(qqs,-qq)&
        +Fahv(-sis*csfs,-sis*sfs,qqs,qqs)&
        *dconjg(Fahv(-sis*csfs,-sis*sfs,-qqs,qqs))*expc4(qqs,-qqs)
      expcadad=Fahv(-si,0.0d0,-qq,qq)&
        *dconjg(Fahv(-si,0.0d0,-qq,qq))*expc1(-qq,-qq)&
        +Fahv(-si,0.0d0,-qq,qq)&
        *dconjg(Fahv(-sis*csfs,-sis*sfs,-qqs,qqs))*expc2(-qq,-qqs)&
        +Fahv(-sis*csfs,-sis*sfs,-qqs,qqs)&
        *dconjg(Fahv(-si,0.0d0,-qq,qq))*expc3(-qqs,-qq)&
        +Fahv(-sis*csfs,-sis*sfs,-qqs,qqs)&
        *dconjg(Fahv(-sis*csfs,-sis*sfs,-qqs,qqs))*expc4(-qqs,-qqs)
      expcbuau=Fbhv(-si,0.0d0,qqt,qqt)&
        *dconjg(Fahv(-si,0.0d0,qq,qq))*expc1(qqt,qq)&
        +Fbhv(-si,0.0d0,qqt,qqt)&
        *dconjg(Fahv(-sis*csfs,-sis*sfs,qqs,qqs))*expc2(qqt,qqs)&
        +Fbhv(-sis*csfs,-sis*sfs,qqts,qqts)&
        *dconjg(Fahv(-si,0.0d0,qq,qq))*expc3(qqts,qq)&
        +Fbhv(-sis*csfs,-sis*sfs,qqts,qqts)&
        *dconjg(Fahv(-sis*csfs,-sis*sfs,qqs,qqs))*expc4(qqts,qqs)
      expcbdau=Fbhv(-si,0.0d0,-qqt,qqt)&
        *dconjg(Fahv(-si,0.0d0,qq,qq))*expc1(-qqt,qq)&
        +Fbhv(-si,0.0d0,-qqt,qqt)&
        *dconjg(Fahv(-sis*csfs,-sis*sfs,qqs,qqs))*expc2(-qqt,qqs)&
        +Fbhv(-sis*csfs,-sis*sfs,-qqts,qqts)&
        *dconjg(Fahv(-si,0.0d0,qq,qq))*expc3(-qqts,qq)&
        +Fbhv(-sis*csfs,-sis*sfs,-qqts,qqts)&
        *dconjg(Fahv(-sis*csfs,-sis*sfs,qqs,qqs))*expc4(-qqts,qqs)
      expcbuad=Fbhv(-si,0.0d0,qqt,qqt)&
        *dconjg(Fahv(-si,0.0d0,-qq,qq))*expc1(qqt,-qq)&
        +Fbhv(-si,0.0d0,qqt,qqt)&
        *dconjg(Fahv(-sis*csfs,-sis*sfs,-qqs,qqs))*expc2(qqt,-qqs)&
        +Fbhv(-sis*csfs,-sis*sfs,qqts,qqts)&
        *dconjg(Fahv(-si,0.0d0,-qq,qq))*expc3(qqts,-qq)&
        +Fbhv(-sis*csfs,-sis*sfs,qqts,qqts)&
        *dconjg(Fahv(-sis*csfs,-sis*sfs,-qqs,qqs))*expc4(qqts,-qqs)
      expcbdad=Fbhv(-si,0.0d0,-qqt,qqt)&
        *dconjg(Fahv(-si,0.0d0,-qq,qq))*expc1(-qqt,-qq)&
        +Fbhv(-si,0.0d0,-qqt,qqt)&
        *dconjg(Fahv(-sis*csfs,-sis*sfs,-qqs,qqs))*expc2(-qqt,-qqs)&
        +Fbhv(-sis*csfs,-sis*sfs,-qqts,qqts)&
        *dconjg(Fahv(-si,0.0d0,-qq,qq))*expc3(-qqts,-qq)&
        +Fbhv(-sis*csfs,-sis*sfs,-qqts,qqts)&
        *dconjg(Fahv(-sis*csfs,-sis*sfs,-qqs,qqs))*expc4(-qqts,-qqs)
      expcaubu=Fahv(-si,0.0d0,qq,qq)&
        *dconjg(Fbhv(-si,0.0d0,qqt,qqt))*expc1(qq,qqt)&
        +Fahv(-si,0.0d0,qq,qq)&
        *dconjg(Fbhv(-sis*csfs,-sis*sfs,qqts,qqts))*expc2(qq,qqts)&
        +Fahv(-sis*csfs,-sis*sfs,qqs,qqs)&
        *dconjg(Fbhv(-si,0.0d0,qqt,qqt))*expc3(qqs,qqt)&
        +Fahv(-sis*csfs,-sis*sfs,qqs,qqs)&
        *dconjg(Fbhv(-sis*csfs,-sis*sfs,qqts,qqts))*expc4(qqs,qqts)
      expcadbu=Fahv(-si,0.0d0,-qq,qq)&
        *dconjg(Fbhv(-si,0.0d0,qqt,qqt))*expc1(-qq,qqt)&
        +Fahv(-si,0.0d0,-qq,qq)&
        *dconjg(Fbhv(-sis*csfs,-sis*sfs,qqts,qqts))*expc2(-qq,qqts)&
        +Fahv(-sis*csfs,-sis*sfs,-qqs,qqs)&
        *dconjg(Fbhv(-si,0.0d0,qqt,qqt))*expc3(-qqs,qqt)&
        +Fahv(-sis*csfs,-sis*sfs,-qqs,qqs)&
        *dconjg(Fbhv(-sis*csfs,-sis*sfs,qqts,qqts))*expc4(-qqs,qqts)
      expcaubd=Fahv(-si,0.0d0,qq,qq)&
        *dconjg(Fbhv(-si,0.0d0,-qqt,qqt))*expc1(qq,-qqt)&
        +Fahv(-si,0.0d0,qq,qq)&
        *dconjg(Fbhv(-sis*csfs,-sis*sfs,-qqts,qqts))*expc2(qq,-qqts)&
        +Fahv(-sis*csfs,-sis*sfs,qqs,qqs)&
        *dconjg(Fbhv(-si,0.0d0,-qqt,qqt))*expc3(qqs,-qqt)&
        +Fahv(-sis*csfs,-sis*sfs,qqs,qqs)&
        *dconjg(Fbhv(-sis*csfs,-sis*sfs,-qqts,qqts))*expc4(qqs,-qqts)
      expcadbd=Fahv(-si,0.0d0,-qq,qq)&
        *dconjg(Fbhv(-si,0.0d0,-qqt,qqt))*expc1(-qq,-qqt)&
        +Fahv(-si,0.0d0,-qq,qq)&
        *dconjg(Fbhv(-sis*csfs,-sis*sfs,-qqts,qqts))*expc2(-qq,-qqts)&
        +Fahv(-sis*csfs,-sis*sfs,-qqs,qqs)&
        *dconjg(Fbhv(-si,0.0d0,-qqt,qqt))*expc3(-qqs,-qqt)&
        +Fahv(-sis*csfs,-sis*sfs,-qqs,qqs)&
        *dconjg(Fbhv(-sis*csfs,-sis*sfs,-qqts,qqts))*expc4(-qqs,-qqts)
      expcbubu=Fbhv(-si,0.0d0,qqt,qqt)&
        *dconjg(Fbhv(-si,0.0d0,qqt,qqt))*expc1(qqt,qqt)&
        +Fbhv(-si,0.0d0,qqt,qqt)&
        *dconjg(Fbhv(-sis*csfs,-sis*sfs,qqts,qqts))*expc2(qqt,qqts)&
        +Fbhv(-sis*csfs,-sis*sfs,qqts,qqts)&
        *dconjg(Fbhv(-si,0.0d0,qqt,qqt))*expc3(qqts,qqt)&
        +Fbhv(-sis*csfs,-sis*sfs,qqts,qqts)&
        *dconjg(Fbhv(-sis*csfs,-sis*sfs,qqts,qqts))*expc4(qqts,qqts)
      expcbdbu=Fbhv(-si,0.0d0,-qqt,qqt)&
        *dconjg(Fbhv(-si,0.0d0,qqt,qqt))*expc1(-qqt,qqt)&
        +Fbhv(-si,0.0d0,-qqt,qqt)&
        *dconjg(Fbhv(-sis*csfs,-sis*sfs,qqts,qqts))*expc2(-qqt,qqts)&
        +Fbhv(-sis*csfs,-sis*sfs,-qqts,qqts)&
        *dconjg(Fbhv(-si,0.0d0,qqt,qqt))*expc3(-qqts,qqt)&
        +Fbhv(-sis*csfs,-sis*sfs,-qqts,qqts)&
        *dconjg(Fbhv(-sis*csfs,-sis*sfs,qqts,qqts))*expc4(-qqts,qqts)
      expcbubd=Fbhv(-si,0.0d0,qqt,qqt)&
        *dconjg(Fbhv(-si,0.0d0,-qqt,qqt))*expc1(qqt,-qqt)&
        +Fbhv(-si,0.0d0,qqt,qqt)&
        *dconjg(Fbhv(-sis*csfs,-sis*sfs,-qqts,qqts))*expc2(qqt,-qqts)&
        +Fbhv(-sis*csfs,-sis*sfs,qqts,qqts)&
        *dconjg(Fbhv(-si,0.0d0,-qqt,qqt))*expc3(qqts,-qqt)&
        +Fbhv(-sis*csfs,-sis*sfs,qqts,qqts)&
        *dconjg(Fbhv(-sis*csfs,-sis*sfs,-qqts,qqts))*expc4(qqts,-qqts)
      expcbdbd=Fbhv(-si,0.0d0,-qqt,qqt)&
        *dconjg(Fbhv(-si,0.0d0,-qqt,qqt))*expc1(-qqt,-qqt)&
        +Fbhv(-si,0.0d0,-qqt,qqt)&
        *dconjg(Fbhv(-sis*csfs,-sis*sfs,-qqts,qqts))*expc2(-qqt,-qqts)&
        +Fbhv(-sis*csfs,-sis*sfs,-qqts,qqts)&
        *dconjg(Fbhv(-si,0.0d0,-qqt,qqt))*expc3(-qqts,-qqt)&
        +Fbhv(-sis*csfs,-sis*sfs,-qqts,qqts)&
        *dconjg(Fbhv(-sis*csfs,-sis*sfs,-qqts,qqts))*expc4(-qqts,-qqts)
      cterm(3)=0.03125*ex2*(expcauau+expcadau+expcauad+expcadad&
                           +expcbuau+expcbdau+expcbuad+expcbdad&
                           +expcaubu+expcadbu+expcaubd+expcadbd&
                           +expcbubu+expcbdbu+expcbubd+expcbdbd)
!
      expcauau=Favh(-si,0.0d0,qq,qq)&
        *dconjg(Favh(-si,0.0d0,qq,qq))*expc1(qq,qq)&
        +Favh(-si,0.0d0,qq,qq)&
        *dconjg(Favh(-sis*csfs,-sis*sfs,qqs,qqs))*expc2(qq,qqs)&
        +Favh(-sis*csfs,-sis*sfs,qqs,qqs)&
        *dconjg(Favh(-si,0.0d0,qq,qq))*expc3(qqs,qq)&
        +Favh(-sis*csfs,-sis*sfs,qqs,qqs)&
        *dconjg(Favh(-sis*csfs,-sis*sfs,qqs,qqs))*expc4(qqs,qqs)
      expcadau=Favh(-si,0.0d0,-qq,qq)&
        *dconjg(Favh(-si,0.0d0,qq,qq))*expc1(-qq,qq)&
        +Favh(-si,0.0d0,-qq,qq)&
        *dconjg(Favh(-sis*csfs,-sis*sfs,qqs,qqs))*expc2(-qq,qqs)&
        +Favh(-sis*csfs,-sis*sfs,-qqs,qqs)&
        *dconjg(Favh(-si,0.0d0,qq,qq))*expc3(-qqs,qq)&
        +Favh(-sis*csfs,-sis*sfs,-qqs,qqs)&
        *dconjg(Favh(-sis*csfs,-sis*sfs,qqs,qqs))*expc4(-qqs,qqs)
      expcauad=Favh(-si,0.0d0,qq,qq)&
        *dconjg(Favh(-si,0.0d0,-qq,qq))*expc1(qq,-qq)&
        +Favh(-si,0.0d0,qq,qq)&
        *dconjg(Favh(-sis*csfs,-sis*sfs,-qqs,qqs))*expc2(qq,-qqs)&
        +Favh(-sis*csfs,-sis*sfs,qqs,qqs)&
        *dconjg(Favh(-si,0.0d0,-qq,qq))*expc3(qqs,-qq)&
        +Favh(-sis*csfs,-sis*sfs,qqs,qqs)&
        *dconjg(Favh(-sis*csfs,-sis*sfs,-qqs,qqs))*expc4(qqs,-qqs)
      expcadad=Favh(-si,0.0d0,-qq,qq)&
        *dconjg(Favh(-si,0.0d0,-qq,qq))*expc1(-qq,-qq)&
        +Favh(-si,0.0d0,-qq,qq)&
        *dconjg(Favh(-sis*csfs,-sis*sfs,-qqs,qqs))*expc2(-qq,-qqs)&
        +Favh(-sis*csfs,-sis*sfs,-qqs,qqs)&
        *dconjg(Favh(-si,0.0d0,-qq,qq))*expc3(-qqs,-qq)&
        +Favh(-sis*csfs,-sis*sfs,-qqs,qqs)&
        *dconjg(Favh(-sis*csfs,-sis*sfs,-qqs,qqs))*expc4(-qqs,-qqs)
      expcbuau=Fbvh(-si,0.0d0,qqt,qqt)&
        *dconjg(Favh(-si,0.0d0,qq,qq))*expc1(qqt,qq)&
        +Fbvh(-si,0.0d0,qqt,qqt)&
        *dconjg(Favh(-sis*csfs,-sis*sfs,qqs,qqs))*expc2(qqt,qqs)&
        +Fbvh(-sis*csfs,-sis*sfs,qqts,qqts)&
        *dconjg(Favh(-si,0.0d0,qq,qq))*expc3(qqts,qq)&
        +Fbvh(-sis*csfs,-sis*sfs,qqts,qqts)&
        *dconjg(Favh(-sis*csfs,-sis*sfs,qqs,qqs))*expc4(qqts,qqs)
      expcbdau=Fbvh(-si,0.0d0,-qqt,qqt)&
        *dconjg(Favh(-si,0.0d0,qq,qq))*expc1(-qqt,qq)&
        +Fbvh(-si,0.0d0,-qqt,qqt)&
        *dconjg(Favh(-sis*csfs,-sis*sfs,qqs,qqs))*expc2(-qqt,qqs)&
        +Fbvh(-sis*csfs,-sis*sfs,-qqts,qqts)&
        *dconjg(Favh(-si,0.0d0,qq,qq))*expc3(-qqts,qq)&
        +Fbvh(-sis*csfs,-sis*sfs,-qqts,qqts)&
        *dconjg(Favh(-sis*csfs,-sis*sfs,qqs,qqs))*expc4(-qqts,qqs)
      expcbuad=Fbvh(-si,0.0d0,qqt,qqt)&
        *dconjg(Favh(-si,0.0d0,-qq,qq))*expc1(qqt,-qq)&
        +Fbvh(-si,0.0d0,qqt,qqt)&
        *dconjg(Favh(-sis*csfs,-sis*sfs,-qqs,qqs))*expc2(qqt,-qqs)&
        +Fbvh(-sis*csfs,-sis*sfs,qqts,qqts)&
        *dconjg(Favh(-si,0.0d0,-qq,qq))*expc3(qqts,-qq)&
        +Fbvh(-sis*csfs,-sis*sfs,qqts,qqts)&
        *dconjg(Favh(-sis*csfs,-sis*sfs,-qqs,qqs))*expc4(qqts,-qqs)
      expcbdad=Fbvh(-si,0.0d0,-qqt,qqt)&
        *dconjg(Favh(-si,0.0d0,-qq,qq))*expc1(-qqt,-qq)&
        +Fbvh(-si,0.0d0,-qqt,qqt)&
        *dconjg(Favh(-sis*csfs,-sis*sfs,-qqs,qqs))*expc2(-qqt,-qqs)&
        +Fbvh(-sis*csfs,-sis*sfs,-qqts,qqts)&
        *dconjg(Favh(-si,0.0d0,-qq,qq))*expc3(-qqts,-qq)&
        +Fbvh(-sis*csfs,-sis*sfs,-qqts,qqts)&
        *dconjg(Favh(-sis*csfs,-sis*sfs,-qqs,qqs))*expc4(-qqts,-qqs)
      expcaubu=Favh(-si,0.0d0,qq,qq)&
        *dconjg(Fbvh(-si,0.0d0,qqt,qqt))*expc1(qq,qqt)&
        +Favh(-si,0.0d0,qq,qq)&
        *dconjg(Fbvh(-sis*csfs,-sis*sfs,qqts,qqts))*expc2(qq,qqts)&
        +Favh(-sis*csfs,-sis*sfs,qqs,qqs)&
        *dconjg(Fbvh(-si,0.0d0,qqt,qqt))*expc3(qqs,qqt)&
        +Favh(-sis*csfs,-sis*sfs,qqs,qqs)&
        *dconjg(Fbvh(-sis*csfs,-sis*sfs,qqts,qqts))*expc4(qqs,qqts)
      expcadbu=Favh(-si,0.0d0,-qq,qq)&
        *dconjg(Fbvh(-si,0.0d0,qqt,qqt))*expc1(-qq,qqt)&
        +Favh(-si,0.0d0,-qq,qq)&
        *dconjg(Fbvh(-sis*csfs,-sis*sfs,qqts,qqts))*expc2(-qq,qqts)&
        +Favh(-sis*csfs,-sis*sfs,-qqs,qqs)&
        *dconjg(Fbvh(-si,0.0d0,qqt,qqt))*expc3(-qqs,qqt)&
        +Favh(-sis*csfs,-sis*sfs,-qqs,qqs)&
        *dconjg(Fbvh(-sis*csfs,-sis*sfs,qqts,qqts))*expc4(-qqs,qqts)
      expcaubd=Favh(-si,0.0d0,qq,qq)&
        *dconjg(Fbvh(-si,0.0d0,-qqt,qqt))*expc1(qq,-qqt)&
        +Favh(-si,0.0d0,qq,qq)&
        *dconjg(Fbvh(-sis*csfs,-sis*sfs,-qqts,qqts))*expc2(qq,-qqts)&
        +Favh(-sis*csfs,-sis*sfs,qqs,qqs)&
        *dconjg(Fbvh(-si,0.0d0,-qqt,qqt))*expc3(qqs,-qqt)&
        +Favh(-sis*csfs,-sis*sfs,qqs,qqs)&
        *dconjg(Fbvh(-sis*csfs,-sis*sfs,-qqts,qqts))*expc4(qqs,-qqts)
      expcadbd=Favh(-si,0.0d0,-qq,qq)&
        *dconjg(Fbvh(-si,0.0d0,-qqt,qqt))*expc1(-qq,-qqt)&
        +Favh(-si,0.0d0,-qq,qq)&
        *dconjg(Fbvh(-sis*csfs,-sis*sfs,-qqts,qqts))*expc2(-qq,-qqts)&
        +Favh(-sis*csfs,-sis*sfs,-qqs,qqs)&
        *dconjg(Fbvh(-si,0.0d0,-qqt,qqt))*expc3(-qqs,-qqt)&
        +Favh(-sis*csfs,-sis*sfs,-qqs,qqs)&
        *dconjg(Fbvh(-sis*csfs,-sis*sfs,-qqts,qqts))*expc4(-qqs,-qqts)
      expcbubu=Fbvh(-si,0.0d0,qqt,qqt)&
        *dconjg(Fbvh(-si,0.0d0,qqt,qqt))*expc1(qqt,qqt)&
        +Fbvh(-si,0.0d0,qqt,qqt)&
        *dconjg(Fbvh(-sis*csfs,-sis*sfs,qqts,qqts))*expc2(qqt,qqts)&
        +Fbvh(-sis*csfs,-sis*sfs,qqts,qqts)&
        *dconjg(Fbvh(-si,0.0d0,qqt,qqt))*expc3(qqts,qqt)&
        +Fbvh(-sis*csfs,-sis*sfs,qqts,qqts)&
        *dconjg(Fbvh(-sis*csfs,-sis*sfs,qqts,qqts))*expc4(qqts,qqts)
      expcbdbu=Fbvh(-si,0.0d0,-qqt,qqt)&
        *dconjg(Fbvh(-si,0.0d0,qqt,qqt))*expc1(-qqt,qqt)&
        +Fbvh(-si,0.0d0,-qqt,qqt)&
        *dconjg(Fbvh(-sis*csfs,-sis*sfs,qqts,qqts))*expc2(-qqt,qqts)&
        +Fbvh(-sis*csfs,-sis*sfs,-qqts,qqts)&
        *dconjg(Fbvh(-si,0.0d0,qqt,qqt))*expc3(-qqts,qqt)&
        +Fbvh(-sis*csfs,-sis*sfs,-qqts,qqts)&
        *dconjg(Fbvh(-sis*csfs,-sis*sfs,qqts,qqts))*expc4(-qqts,qqts)
      expcbubd=Fbvh(-si,0.0d0,qqt,qqt)&
        *dconjg(Fbvh(-si,0.0d0,-qqt,qqt))*expc1(qqt,-qqt)&
        +Fbvh(-si,0.0d0,qqt,qqt)&
        *dconjg(Fbvh(-sis*csfs,-sis*sfs,-qqts,qqts))*expc2(qqt,-qqts)&
        +Fbvh(-sis*csfs,-sis*sfs,qqts,qqts)&
        *dconjg(Fbvh(-si,0.0d0,-qqt,qqt))*expc3(qqts,-qqt)&
        +Fbvh(-sis*csfs,-sis*sfs,qqts,qqts)&
        *dconjg(Fbvh(-sis*csfs,-sis*sfs,-qqts,qqts))*expc4(qqts,-qqts)
      expcbdbd=Fbvh(-si,0.0d0,-qqt,qqt)&
        *dconjg(Fbvh(-si,0.0d0,-qqt,qqt))*expc1(-qqt,-qqt)&
        +Fbvh(-si,0.0d0,-qqt,qqt)&
        *dconjg(Fbvh(-sis*csfs,-sis*sfs,-qqts,qqts))*expc2(-qqt,-qqts)&
        +Fbvh(-sis*csfs,-sis*sfs,-qqts,qqts)&
        *dconjg(Fbvh(-si,0.0d0,-qqt,qqt))*expc3(-qqts,-qqt)&
        +Fbvh(-sis*csfs,-sis*sfs,-qqts,qqts)&
        *dconjg(Fbvh(-sis*csfs,-sis*sfs,-qqts,qqts))*expc4(-qqts,-qqts)
      cterm(4)=0.03125*ex2*(expcauau+expcadau+expcauad+expcadad&
                           +expcbuau+expcbdau+expcbuad+expcbdad&
                           +expcaubu+expcadbu+expcaubd+expcadbd&
                           +expcbubu+expcbdbu+expcbubd+expcbdbd)
!
!----------------------------------------------c
! end of computation of complementary terms    c
!----------------------------------------------c
!
      do 11 ip=1,4
!      if(dreal(kcterm(ip)).lt.0.) kcterm(ip)=(0.,0)
!      if(dreal(cterm(ip)).lt.0.) cterm(ip)=(0.,0)
      sigma0(ip)=dreal(kterm(ip)+kcterm(ip)+cterm(ip))
 11  continue      
!      sigma0(1)=dreal(kterm(1)+kcterm(1)+cterm(1))
!      sigma0(2)=dreal(kterm(2)+kcterm(2)+cterm(2))
!        sigma0(3)=dreal(kterm(3)+kcterm(3)+cterm(3))
!        sigma0(4)=dreal(kterm(4)+kcterm(4)+cterm(4))
!
!----------------------------------------------------------------c
      return

        endif

      end
!
!------------------------------------------------------------------c
      function expkc1(q)
      implicit real*8 (a-h,k,o-z)
        complex*16 expkc1,q,temp,sum
        real*8 w(30000)
        common /cont1/ si,sis,cs,css,sfs,csfs,si2,sis2,cs2,css2
        common /cont2/ ks2
        common /cont3/ iterm
        common /spectra/ w
      sum=0.0
        temp=1.0
       do 80 n=1,iterm
        fn=float(n)
         temp=temp*(ks2*(css-q)*(css+cs))/fn
         sum=sum+temp*w(n)
80      continue  
      expkc1=cdexp(-ks2*(q*q-q*css+q*cs))*sum
      return
      end
!
!------------------------------------------------------------------c
      function expkc2(q)
      implicit real*8 (a-h,k,o-z)
        complex*16 expkc2,q,temp,sum
        real*8 w(30000)
        common /cont1/ si,sis,cs,css,sfs,csfs,si2,sis2,cs2,css2
        common /cont2/ ks2
        common /cont3/ iterm
        common /spectra/ w
      sum=0.0
        temp=1.0
       do 81 n=1,iterm
         fn=float(n)
         temp=temp*(ks2*(cs+q)*(css+cs))/fn
         sum=sum+temp*w(n)
81     continue
      expkc2=cdexp(-ks2*(q*q-q*css+q*cs))*sum
      return
        end
!
!------------------------------------------------------------------c
      function expc1(q,qp)
      implicit real*8 (a-h,k,o-z)
        complex*16 expc1,q,qp,temp,sum
        real*8 w(30000)
        common /cont1/ si,sis,cs,css,sfs,csfs,si2,sis2,cs2,css2
        common /cont2/ ks2
        common /cont3/ iterm
        common /spectra/ w
      sum=0.0
        temp=1.0
       do 82 n=1,iterm
         fn=float(n)
         temp=temp*(ks2*(css-q)*(css-qp))/fn
        sum=sum+temp*w(n)
82     continue
      expc1=cdexp(-ks2*(q*q+qp*qp-css*(q+qp)+cs*(q+qp)))*sum
      return
      end
!
!------------------------------------------------------------------c
      function expc2(q,qp)
      implicit real*8 (a-h,k,o-z)
        complex*16 expc2,q,qp,temp,sum
        real*8 w(30000)
        common /cont1/ si,sis,cs,css,sfs,csfs,si2,sis2,cs2,css2
        common /cont2/ ks2
        common /cont3/ iterm
        common /spectra/ w
      sum=0.0
        temp=1.0
       do 83 n=1,iterm
        fn=float(n)
         temp=temp*(ks2*(css-q)*(cs+qp))/fn
         sum=sum+temp*w(n)
83     continue
      expc2=cdexp(-ks2*(q*q+qp*qp-css*(q+qp)+cs*(q+qp)))*sum
      return
        end
!
!------------------------------------------------------------------c
      function expc3(q,qp)
      implicit real*8 (a-h,k,o-z)
        complex*16 expc3,q,qp,temp,sum
        real*8 w(30000)
        common /cont1/ si,sis,cs,css,sfs,csfs,si2,sis2,cs2,css2
        common /cont2/ ks2
        common /cont3/ iterm
        common /spectra/ w
      sum=0.0
        temp=1.0
       do 80 n=1,iterm
         fn=float(n)
         temp=temp*(ks2*(cs+q)*(css-qp))/fn
         sum=sum+temp*w(n)
80     continue
      expc3=cdexp(-ks2*(q*q+qp*qp-css*(q+qp)+cs*(q+qp)))*sum
      return
        end
!
!------------------------------------------------------------------c
      function expc4(q,qp)
      implicit real*8 (a-h,k,o-z)
        complex*16 expc4,q,qp,temp,sum
        real*8 w(30000)
        common /cont1/ si,sis,cs,css,sfs,csfs,si2,sis2,cs2,css2
        common /cont2/ ks2
        common /cont3/ iterm
        common /spectra/ w
      sum=0.0
        temp=1.0
       do 80 n=1,iterm
        fn=float(n)
         temp=temp*(ks2*(cs+q)*(cs+qp))/fn
        sum=sum+temp*w(n)
80     continue
      expc4=cdexp(-ks2*(q*q+qp*qp-css*(q+qp)+cs*(q+qp)))*sum
      return
        end
!
!-------------------------------------------------------------------c
      function Favv(u,v,q,qfix)
      implicit real*8 (a-h,k,o-z)
      complex*16 er,ur,Favv,rh,rv,rvh,q,qfix
      complex*16 c1,c2,c3,c4,c5,c6
        complex*16 zx,zy,zxp,zyp
        complex*16 rpv,rmv,av,bv
        common /cont1/ si,sis,cs,css,sfs,csfs,si2,sis2,cs2,css2
      common /eu/er,ur
      common /rhhvv/rh,rv,rvh
      kxu=si+u
        ksxu=sis*csfs+u
        kyv=v
        ksyv=sis*sfs+v
        if((dabs(dreal(css-q)).lt.0.000001).or.&
         (dabs(dreal(cs+q)).lt.0.000001)) then
        c1=0.0
        c2=0.0
        c3=0.0
        c4=0.0
        c5=0.0
        c6=0.0
        else
      zx=(-ksxu)/(css-q)
      zy=(-ksyv)/(css-q)
      zxp=(kxu)/(cs+q)
      zyp=(kyv)/(cs+q)
      c1=-csfs*(-1.0-zx*zxp)+sfs*zxp*zy
      c2=-csfs*(-cs*q-cs*u*zx-q*si*zxp-si*u*zx*zxp-cs*v*zyp-si*v*zx*zyp)&
         +sfs*(cs*u*zy+si*u*zxp*zy+q*si*zyp-cs*u*zyp+si*v*zy*zyp)
      c3=-csfs*(si*u-q*si*zx-cs*u*zxp+cs*q*zx*zxp)&
         +sfs*(-si*v+cs*v*zxp+q*si*zy-cs*q*zxp*zy)
      c4=-css*sfs*(-si*zyp+cs*zx*zyp)-csfs*css*(-cs-si*zxp-cs*zy*zyp)&
         +sis*(-cs*zx-si*zx*zxp-si*zy*zyp)
      c5=-css*sfs*(-v*zx+v*zxp)-csfs*css*(q+u*zxp+v*zy)&
         +sis*(q*zx+u*zx*zxp+v*zxp*zy)
        c6=-css*sfs*(-u*zyp+q*zx*zyp)-csfs*css*(v*zyp-q*zy*zyp)&
         +sis*(v*zx*zyp-u*zy*zyp)
        end if
!
        rpv=1.0+rv
        rmv=1.0-rv
      av=rpv/qfix
        bv=rmv/qfix
      Favv=bv*(-rpv*c1+rmv*c2+rpv*c3)+av*(rmv*c4+rpv*c5+rmv*c6)
      return
      end
!
!-------------------------------------------------------------------c
      function Fahh(u,v,q,qfix)
      implicit real*8 (a-h,k,o-z)
      complex*16 er,ur,Fahh,rh,rv,rvh,q,qfix
      complex*16 c1,c2,c3,c4,c5,c6
        complex*16 zx,zy,zxp,zyp
        complex*16 rph,rmh,ah,bh
        common /cont1/ si,sis,cs,css,sfs,csfs,si2,sis2,cs2,css2
      common /eu/er,ur
      common /rhhvv/rh,rv,rvh
      kxu=si+u
        ksxu=sis*csfs+u
        kyv=v
        ksyv=sis*sfs+v
        if((dabs(dreal(css-q)).lt.0.000001).or.&
         (dabs(dreal(cs+q)).lt.0.000001)) then
        c1=0.0
        c2=0.0
        c3=0.0
        c4=0.0
        c5=0.0
        c6=0.0
        else
      zx=(-ksxu)/(css-q)
      zy=(-ksyv)/(css-q)
      zxp=(kxu)/(cs+q)
      zyp=(kyv)/(cs+q)
      c1=-csfs*(-1.0-zx*zxp)+sfs*zxp*zy
      c2=-csfs*(-cs*q-cs*u*zx-q*si*zxp-si*u*zx*zxp-cs*v*zyp-si*v*zx*zyp)&
         +sfs*(cs*u*zy+si*u*zxp*zy+q*si*zyp-cs*u*zyp+si*v*zy*zyp)
      c3=-csfs*(si*u-q*si*zx-cs*u*zxp+cs*q*zx*zxp)&
         +sfs*(-si*v+cs*v*zxp+q*si*zy-cs*q*zxp*zy)
      c4=-css*sfs*(-si*zyp+cs*zx*zyp)-csfs*css*(-cs-si*zxp-cs*zy*zyp)&
         +sis*(-cs*zx-si*zx*zxp-si*zy*zyp)
      c5=-css*sfs*(-v*zx+v*zxp)-csfs*css*(q+u*zxp+v*zy)&
         +sis*(q*zx+u*zx*zxp+v*zxp*zy)
        c6=-css*sfs*(-u*zyp+q*zx*zyp)-csfs*css*(v*zyp-q*zy*zyp)&
         +sis*(v*zx*zyp-u*zy*zyp)
        end if
!
        rph=1.0+rh
        rmh=1.0-rh
      ah=rph/qfix
        bh=rmh/qfix
      Fahh=-bh*(-rph*c1+rmh*c2+rph*c3)-ah*(rmh*c4+rph*c5+rmh*c6)
      return
      end
!
!-------------------------------------------------------------------c
      function Fbvv(u,v,q,qfix)
      implicit real*8 (a-h,k,o-z)
      complex*16 er,ur,Fbvv,rh,rv,rvh,q,qfix
      complex*16 c1,c2,c3,c4,c5,c6
        complex*16 zx,zy,zxp,zyp
        complex*16 rpv,rmv,av,bv
        common /cont1/ si,sis,cs,css,sfs,csfs,si2,sis2,cs2,css2
      common /eu/er,ur
      common /rhhvv/rh,rv,rvh
      kxu=si+u
        ksxu=sis*csfs+u
        kyv=v
        ksyv=sis*sfs+v
      zx=(-ksxu)/(css-q)
      zy=(-ksyv)/(css-q)
      zxp=(kxu)/(cs+q)
      zyp=(kyv)/(cs+q)
      c1=-csfs*(-1.0-zx*zxp)+sfs*zxp*zy
      c2=-csfs*(-cs*q-cs*u*zx-q*si*zxp-si*u*zx*zxp-cs*v*zyp-si*v*zx*zyp)&
         +sfs*(cs*u*zy+si*u*zxp*zy+q*si*zyp-cs*u*zyp+si*v*zy*zyp)
      c3=-csfs*(si*u-q*si*zx-cs*u*zxp+cs*q*zx*zxp)&
         +sfs*(-si*v+cs*v*zxp+q*si*zy-cs*q*zxp*zy)
      c4=-css*sfs*(-si*zyp+cs*zx*zyp)-csfs*css*(-cs-si*zxp-cs*zy*zyp)&
         +sis*(-cs*zx-si*zx*zxp-si*zy*zyp)
      c5=-css*sfs*(-v*zx+v*zxp)-csfs*css*(q+u*zxp+v*zy)&
         +sis*(q*zx+u*zx*zxp+v*zxp*zy)
        c6=-css*sfs*(-u*zyp+q*zx*zyp)-csfs*css*(v*zyp-q*zy*zyp)&
         +sis*(v*zx*zyp-u*zy*zyp)
        rpv=1.0+rv
        rmv=1.0-rv
      av=rpv/qfix
        bv=rmv/qfix
      Fbvv=av*(rpv*c1-rmv*c2-rpv*c3/er)-bv*(rmv*c4*er+rpv*c5+rmv*c6)
      return
      end
!
!-------------------------------------------------------------------c
      function Fbhh(u,v,q,qfix)
      implicit real*8 (a-h,k,o-z)
      complex*16 er,ur,Fbhh,rh,rv,rvh,q,qfix
      complex*16 c1,c2,c3,c4,c5,c6
        complex*16 zx,zy,zxp,zyp
        complex*16 rph,rmh,ah,bh
        common /cont1/ si,sis,cs,css,sfs,csfs,si2,sis2,cs2,css2
      common /eu/er,ur
      common /rhhvv/rh,rv,rvh
      kxu=si+u
        ksxu=sis*csfs+u
        kyv=v
        ksyv=sis*sfs+v
      zx=(-ksxu)/(css-q)
      zy=(-ksyv)/(css-q)
      zxp=(kxu)/(cs+q)
      zyp=(kyv)/(cs+q)
      c1=-csfs*(-1.0-zx*zxp)+sfs*zxp*zy
      c2=-csfs*(-cs*q-cs*u*zx-q*si*zxp-si*u*zx*zxp-cs*v*zyp-si*v*zx*zyp)&
         +sfs*(cs*u*zy+si*u*zxp*zy+q*si*zyp-cs*u*zyp+si*v*zy*zyp)
      c3=-csfs*(si*u-q*si*zx-cs*u*zxp+cs*q*zx*zxp)&
         +sfs*(-si*v+cs*v*zxp+q*si*zy-cs*q*zxp*zy)
      c4=-css*sfs*(-si*zyp+cs*zx*zyp)-csfs*css*(-cs-si*zxp-cs*zy*zyp)&
         +sis*(-cs*zx-si*zx*zxp-si*zy*zyp)
      c5=-css*sfs*(-v*zx+v*zxp)-csfs*css*(q+u*zxp+v*zy)&
         +sis*(q*zx+u*zx*zxp+v*zxp*zy)
        c6=-css*sfs*(-u*zyp+q*zx*zyp)-csfs*css*(v*zyp-q*zy*zyp)&
         +sis*(v*zx*zyp-u*zy*zyp)
        rph=1.0+rh
        rmh=1.0-rh
      ah=rph/qfix
        bh=rmh/qfix
      Fbhh=ah*(-rph*c1*er+rmh*c2+rph*c3)+bh*(rmh*c4+rph*c5+rmh*c6/er)
      return
      end
!
!-------------------------------------------------------------------c
      function Fahv(u,v,q,qfix)
      implicit real*8 (a-h,k,o-z)
      complex*16 er,ur,Fahv,rh,rv,rvh,q,qfix
      complex*16 b1,b2,b3,b4,b5,b6
        complex*16 zx,zy,zxp,zyp
        complex*16 rp,rm,a,b
        common /cont1/ si,sis,cs,css,sfs,csfs,si2,sis2,cs2,css2
      common /eu/er,ur
      common /rhhvv/rh,rv,rvh
      kxu=si+u
        ksxu=sis*csfs+u
        kyv=v
        ksyv=sis*sfs+v
        if((dabs(dreal(css-q)).lt.0.000001).or.&
         (dabs(dreal(cs+q)).lt.0.000001)) then
                b1=0.0
                b2=0.0
                b3=0.0
                b4=0.0
                b5=0.0
                b6=0.0
      else
      zx=(-ksxu)/(css-q)
      zy=(-ksyv)/(css-q)
      zxp=(kxu)/(cs+q)
      zyp=(kyv)/(cs+q)
      b1=-css*sfs*(-1.0-zx*zxp)-sis*zy-csfs*css*zxp*zy
        b2=-css*sfs*(-cs*q-cs*u*zx-q*si*zxp-si*u*zx*zxp-cs*v*zyp &
          -si*v*zx*zyp)+ &
       sis*(-cs*q*zy-q*si*zxp*zy+q*si*zx*zyp-cs*u*zx*zyp-cs*v*zy*zyp)- &
       csfs*css*(cs*u*zy+si*u*zxp*zy+q*si*zyp-cs*u*zyp+si*v*zy*zyp)
      b3=-css*sfs*(si*u-q*si*zx-cs*u*zxp+cs*q*zx*zxp)- &
         csfs*css*(-si*v+cs*v*zxp+q*si*zy-cs*q*zxp*zy)+ &
         sis*(-si*v*zx+cs*v*zx*zxp+si*u*zy-cs*u*zxp*zy)
        b4=-csfs*(-si*zyp+cs*zx*zyp)+sfs*(-cs-si*zxp-cs*zy*zyp)
      b5=-csfs*(-v*zx+v*zxp)+sfs*(q+u*zxp+v*zy)
        b6=-csfs*(-u*zyp+q*zx*zyp)+sfs*(v*zyp-q*zy*zyp)
        end if
!
        rp=1.0+rvh
        rm=1.0-rvh
      a=rp/qfix
        b=rm/qfix
      Fahv=b*(rp*b1-rm*b2-rp*b3)+a*(rm*b4+rp*b5+rm*b6)
      return
      end
!
!-------------------------------------------------------------------c
      function Favh(u,v,q,qfix)
      implicit real*8 (a-h,k,o-z)
      complex*16 er,ur,Favh,rh,rv,rvh,q,qfix
      complex*16 b1,b2,b3,b4,b5,b6
        complex*16 zx,zy,zxp,zyp
        complex*16 rp,rm,a,b
        common /cont1/ si,sis,cs,css,sfs,csfs,si2,sis2,cs2,css2
      common /eu/er,ur
      common /rhhvv/rh,rv,rvh
      kxu=si+u
        ksxu=sis*csfs+u
        kyv=v
        ksyv=sis*sfs+v
        if((dabs(dreal(css-q)).lt.0.000001).or.&
         (dabs(dreal(cs+q)).lt.0.000001)) then
        b1=0.0
        b2=0.0
        b3=0.0
        b4=0.0
        b5=0.0
        b6=0.0
        else
      zx=(-ksxu)/(css-q)
      zy=(-ksyv)/(css-q)
      zxp=(kxu)/(cs+q)
      zyp=(kyv)/(cs+q)
      b1=-css*sfs*(-1.0-zx*zxp)-sis*zy-csfs*css*zxp*zy
        b2=-css*sfs*(-cs*q-cs*u*zx-q*si*zxp-si*u*zx*zxp-cs*v*zyp &
          -si*v*zx*zyp)+ &
       sis*(-cs*q*zy-q*si*zxp*zy+q*si*zx*zyp-cs*u*zx*zyp-cs*v*zy*zyp)- &
       csfs*css*(cs*u*zy+si*u*zxp*zy+q*si*zyp-cs*u*zyp+si*v*zy*zyp)
      b3=-css*sfs*(si*u-q*si*zx-cs*u*zxp+cs*q*zx*zxp)- &
         csfs*css*(-si*v+cs*v*zxp+q*si*zy-cs*q*zxp*zy)+ &
         sis*(-si*v*zx+cs*v*zx*zxp+si*u*zy-cs*u*zxp*zy)
        b4=-csfs*(-si*zyp+cs*zx*zyp)+sfs*(-cs-si*zxp-cs*zy*zyp)
      b5=-csfs*(-v*zx+v*zxp)+sfs*(q+u*zxp+v*zy)
        b6=-csfs*(-u*zyp+q*zx*zyp)+sfs*(v*zyp-q*zy*zyp)
      end if
!
        rp=1.0+rvh
        rm=1.0-rvh
      a=rp/qfix
        b=rm/qfix
      Favh=b*(rp*b4+rm*b5+rp*b6)-a*(-rm*b1+rp*b2+rm*b3)
      return
      end
!
!-------------------------------------------------------------------c
      function Fbhv(u,v,q,qfix)
      implicit real*8 (a-h,k,o-z)
      complex*16 er,ur,Fbhv,rh,rv,rvh,q,qfix
      complex*16 b1,b2,b3,b4,b5,b6
        complex*16 zx,zy,zxp,zyp
        complex*16 rp,rm,a,b
        common /cont1/ si,sis,cs,css,sfs,csfs,si2,sis2,cs2,css2
      common /eu/er,ur
      common /rhhvv/rh,rv,rvh
      kxu=si+u
        ksxu=sis*csfs+u
        kyv=v
        ksyv=sis*sfs+v
      zx=(-ksxu)/(css-q)
      zy=(-ksyv)/(css-q)
      zxp=(kxu)/(cs+q)
      zyp=(kyv)/(cs+q)
      b1=-css*sfs*(-1.0-zx*zxp)-sis*zy-csfs*css*zxp*zy
        b2=-css*sfs*(-cs*q-cs*u*zx-q*si*zxp-si*u*zx*zxp-cs*v*zyp &
          -si*v*zx*zyp)+ &
       sis*(-cs*q*zy-q*si*zxp*zy+q*si*zx*zyp-cs*u*zx*zyp-cs*v*zy*zyp)- &
       csfs*css*(cs*u*zy+si*u*zxp*zy+q*si*zyp-cs*u*zyp+si*v*zy*zyp)
      b3=-css*sfs*(si*u-q*si*zx-cs*u*zxp+cs*q*zx*zxp)- &
         csfs*css*(-si*v+cs*v*zxp+q*si*zy-cs*q*zxp*zy)+ &
         sis*(-si*v*zx+cs*v*zx*zxp+si*u*zy-cs*u*zxp*zy)
        b4=-csfs*(-si*zyp+cs*zx*zyp)+sfs*(-cs-si*zxp-cs*zy*zyp)
      b5=-csfs*(-v*zx+v*zxp)+sfs*(q+u*zxp+v*zy)
        b6=-csfs*(-u*zyp+q*zx*zyp)+sfs*(v*zyp-q*zy*zyp)
        rp=1.0+rvh
        rm=1.0-rvh
      a=rp/qfix
        b=rm/qfix
      Fbhv=a*(-rp*b1+rm*b2+rp*b3/er)-b*(rm*b4*er+rp*b5+rm*b6)
      return
      end
!
!-------------------------------------------------------------------c
      function Fbvh(u,v,q,qfix)
      implicit real*8 (a-h,k,o-z)
      complex*16 er,ur,Fbvh,rh,rv,rvh,q,qfix
      complex*16 b1,b2,b3,b4,b5,b6
        complex*16 zx,zy,zxp,zyp
        complex*16 rp,rm,a,b
        common /cont1/ si,sis,cs,css,sfs,csfs,si2,sis2,cs2,css2
      common /eu/er,ur
      common /rhhvv/rh,rv,rvh
      kxu=si+u
        ksxu=sis*csfs+u
        kyv=v
        ksyv=sis*sfs+v
      zx=(-ksxu)/(css-q)
      zy=(-ksyv)/(css-q)
      zxp=(kxu)/(cs+q)
      zyp=(kyv)/(cs+q)
      b1=-css*sfs*(-1.0-zx*zxp)-sis*zy-csfs*css*zxp*zy
        b2=-css*sfs*(-cs*q-cs*u*zx-q*si*zxp-si*u*zx*zxp-cs*v*zyp &
          -si*v*zx*zyp)+ &
       sis*(-cs*q*zy-q*si*zxp*zy+q*si*zx*zyp-cs*u*zx*zyp-cs*v*zy*zyp)- &
       csfs*css*(cs*u*zy+si*u*zxp*zy+q*si*zyp-cs*u*zyp+si*v*zy*zyp)
      b3=-css*sfs*(si*u-q*si*zx-cs*u*zxp+cs*q*zx*zxp)- &
         csfs*css*(-si*v+cs*v*zxp+q*si*zy-cs*q*zxp*zy)+ &
         sis*(-si*v*zx+cs*v*zx*zxp+si*u*zy-cs*u*zxp*zy)
        b4=-csfs*(-si*zyp+cs*zx*zyp)+sfs*(-cs-si*zxp-cs*zy*zyp)
      b5=-csfs*(-v*zx+v*zxp)+sfs*(q+u*zxp+v*zy)
        b6=-csfs*(-u*zyp+q*zx*zyp)+sfs*(v*zyp-q*zy*zyp)
        rp=1.0+rvh
        rm=1.0-rvh
      a=rp/qfix
        b=rm/qfix
      Fbvh=-a*(rp*b4+rm*b5+rp*b6/er)+b*(-rm*b1*er+rp*b2+rm*b3)
      return
      end
!
!
!-------------------------------------------------------------------c
!
      function bessj0(x)
      real*8 x,y,p1,p2,p3,p4,p5,q1,q2,q3,q4,q5,r1,r2,r3,r4,r5,r6,&
      s1,s2,s3,s4,s5,s6,bessj0
      data r1,r2,r3,r4,r5,r6/57568490574.d0,-13362590354.d0,&
      651619640.7d0,-11214424.18d0,77392.33017d0,-184.9052456d0/,&
      s1,s2,s3,s4,s5,s6/57568490411.d0,1029532985.d0,9494680.718d0,&
      59272.64853d0,267.8532712d0,1.d0/
      data p1,p2,p3,p4,p5/1.d0,-.1098628627d-2,.2734510407d-4,&
      -.2073370639d-5,.2093887211d-6/,q1,q2,q3,q4,q5/-.1562499995d-1,&
      .1430488765d-3,-.6911147651d-5,.7621095161d-6,-.934945152d-7/
      if(abs(x).lt.8.)then
        y=x**2
      bessj0=(r1+y*(r2+y*(r3+y*(r4+y*(r5+y*r6)))))&
        /(s1+y*(s2+y*(s3+y*(s4+y*(s5+y*s6)))))
        else
        ax=abs(x)
        z=8./ax
        y=z**2
        xx=ax-.785398164
      bessj0=sqrt(.636619772/ax)*(cos(xx)*(p1+y*(p2+y*(p3+y*(p4+y&
        *p5))))-z*sin(xx)*(q1+y*(q2+y*(q3+y*(q4+y*q5)))))
        endif
      return
      end

      subroutine shadowg(ti,ts,s,shfct)
      implicit real*8 (a-h,o-z)
!       if(ts.ge.ti) then
!         arg=ts
!         else
!           arg=ti
!        endif
        arg=ti
        if(arg.eq.0.0) then
          shfct=1.0
           return
        endif
       u=1.0/dtand(arg)
       pi=acos(-1.0d0)
       et=u/(sqrt(2.0)*s)
       if(et.ge.20) then
          shfct=1.0
            return
        endif
       f1=dsqrt(2.0/pi)*s*dexp(-et*et)/u
        f2=erfc(et)
       f=(f1-f2)/2.0
       shfct=1.0/(1.0+f)
!       shfct=(1.0-0.5*erfc(et))/(1.0+f)
       return
       end
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~c
! subroutine calculates shadowing function     c
! (Tsang et al 1985, pp.95                     c
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~c
      subroutine shadow(back,ti,ts,shfct)
      implicit real*8 (a-h,o-z)
      logical back
      common /al/ effslop
      s=effslop
        if(ti.eq.0.0.and.ts.eq.0.0) then
          shfct=1.0
           return
        endif
       if(ti.eq.0.0) then
            ui=1.0d+30
          else
            ui=1.0/tand(ti)
       endif
       if(ts.eq.0.0) then
            us=1.0d+30
           else
            us=1.0/tand(ts)
       endif
       pi=acos(-1.0d0)
       eti=ui/(sqrt(2.0)*s)
       ets=us/(sqrt(2.0)*s)
      
       if(.not.back) then
        f1i=sqrt(2.0/pi)*s/ui*exp(-eti*eti)
       f1s=sqrt(2.0/pi)*s/us*exp(-ets*ets)
       f2i=erfc(eti)
       f2s=erfc(ets)
       fi=(f1i-f2i)/2.0
       fs=(f1s-f2s)/2.0
        shfct=1.0/(1.0+fi+fs)
       else
       if(ts.ge.ti) then
        f1s=sqrt(2.0/pi)*s/us*exp(-ets*ets)
        f2s=erfc(ets)
        fs=(f1s-f2s)/2.0
          shfct=(1.0-0.5*f2s)/(1.0+fs) 
        else
           f1i=sqrt(2.0/pi)*s/ui*exp(-eti*eti)
           f2i=erfc(eti)
           fi=(f1i-f2i)/2.0
          shfct=(1.0-0.5*f2i)/(fi+1.0) 
        endif
        endif
       return
       end

      FUNCTION erfcc(x)
      REAL*8 erfcc,x
      REAL*8 t,z
      z=dabs(x)
      t=1.d0/(1.d0+0.5d0*z)
      erfcc=t*exp(-z*z-1.26551223d0+t*(1.00002368d0+t*(.37409196d0+t*&
      (.09678418d0+t*(-.18628806d0+t*(.27886807d0+t*(-1.13520398d0+t*&
      (1.48851587d0+t*(-.82215223d0+t*.17087277d0)))))))))
      if (x.lt.0.d0) erfcc=2.d0-erfcc
      return
      END

      function erfc(x)
      implicit real*8 (a-h,o-z)
      if(x.lt.0.d0)then
        erfc=1.d0+gammp(.5d0,x**2.d0)
      else
        erfc=gammq(.5d0,x**2.d0)
      endif
      return
      end

      function gammp(a,x)
      implicit real*8 (a-h,o-z)
      if(x.lt.0.d0.or.a.le.0.d0) pause
      if(x.lt.a+1.d0)then
        call gser(gamser,a,x,gln)
        gammp=gamser
      else
        call gcf(gammcf,a,x,gln)
        gammp=1.d0-gammcf
      endif
      return
      end

      function gammq(a,x)
      implicit real*8 (a-h,o-z)
      if(x.lt.0.d0.or.a.le.0.d0)pause
      if(x.lt.a+1.d0)then
        call gser(gamser,a,x,gln)
        gammq=1.-gamser
      else
        call gcf(gammcf,a,x,gln)
        gammq=gammcf
      endif
      return
      end

      subroutine gser(gamser,a,x,gln)
      implicit real*8 (a-h,o-z)
      parameter (itmax=100,eps=1.d-18)
      gln=gammln(a)
      if(x.le.0.d0)then
        if(x.lt.0.d0)pause
        gamser=0.d0
        return
      endif
      ap=a
      sum=1.d0/a
      del=sum
      do 11 n=1,itmax
        ap=ap+1.
        del=del*x/ap
        sum=sum+del
        if(dabs(del).lt.dabs(sum)*eps)go to 1
11    continue
      pause 'a too large, itmax too small'
1     gamser=sum*dexp(-x+a*dlog(x)-gln)
      return
      end

      subroutine gcf(gammcf,a,x,gln)
      implicit real*8 (a-h,o-z)
      parameter (itmax=100,eps=1.d-18)
      gln=gammln(a)
      gold=0.0d0
      a0=1.0d0
      a1=x
      b0=0.0d0
      b1=1.0d0
      fac=1.0d0
      do 11 n=1,itmax
        an=dfloat(n)
        ana=an-a
        a0=(a1+a0*ana)*fac
        b0=(b1+b0*ana)*fac
        anf=an*fac
        a1=x*a0+anf*a1
        b1=x*b0+anf*b1
        if(a1.ne.0.d0)then
          fac=1.0d0/a1
          g=b1*fac
          if(dabs((g-gold)/g).lt.eps)go to 1
          gold=g
        endif
11    continue
      pause 'a too large, itmax too small'
1     gammcf=dexp(-x+a*dlog(x)-gln)*g
      return
      end

      function gammln(xx)
      implicit real*8 (a-h,o-z)
      dimension cof(6)
      data cof,stp/76.18009173D0,-86.50532033D0,24.01409822D0,&
          -1.231739516D0,.120858003D-2,-.536382D-5,2.50662827465D0/
      data half,one,fpf/0.5D0,1.0D0,5.5D0/
      x=xx-one
      tmp=x+fpf
      tmp=(x+half)*dlog(tmp)-tmp
      ser=one
      do 11 j=1,6
        x=x+one
        ser=ser+cof(j)/x
11    continue
      gammln=tmp+dlog(stp*ser)
      return
      end
!
!
!
      function BesselK(n,x)
!***************************************************
!  Modified Bessel function of order n+0.5
!
!  Input parameters:
!    n : integer part of order
!    x : real parameter
!
!***************************************************

      implicit real*8 (a - h, k, o - z)
      if(x.eq.0.0d0) then
        print*,'BesselK: Singularity encountered !'
        return
      endif
      PI = 4.0d0*atan(1.0d0)
      cons = dsqrt(PI/(2.0d0*x))
      K0 = cons*exp(-x)
      if(n.eq.0) then
        BesselK = K0
        return
      endif
      K1 = K0*(1.0d0+1.0d0/x)
      if(n.eq.1) then
        BesselK = K1
        return
      endif
      K2 = K0*(1.0d0+3.0d0/x+3.0d0/x/x)
      if(n.eq.2) then
        BesselK = K2
        return
      endif
      fn0 = K1*cons
      fn1 = -K2*cons
      do 10 i=2,n-1
        fn2 = fn0-float(2*i+1)/x*fn1
        fn0 = fn1
        fn1 = fn2
10    continue
      BesselK = dabs(fn2/cons)
      return
      end
!
!
!
      function alogam( x )
!
! evaluates natural logarithm of GAMMA(x)
! for x > 0
!
      implicit real*8 ( a - h, o - z)
      pi = acos (-1.0d0)
      a1 = dlog( 2 * pi) / 2.0d0
      a2 = 1.0 / 1680.0d0
      a3 = 1.0 / 1260.0d0
      a4 = 1.0 / 360.0d0
      a5 = 1.0 / 12.0d0
      alogam = 0.0d0
      ifault = 1
      if (x .le. 0.0d0) return
      ifault = 0
      y = x
      f = 0.0d0
      if (y .ge. 7.0d0) go to 30
      f = y  
10    y = y + 1.0d0
      if (y .ge. 7.0d0) go to 20
      f = f * y
      go to 10
20    f = -dlog(f)
30    z = 1.0d0 / ( y * y)
      alogam = f + (y - 0.5d0) * dlog(y) - y + a1&
               + ((( -a2 * z + a3) * z - a4) * z&
               + a5) / y
      return
      end

      FUNCTION BESSK(N,X)
      implicit real*8 (a-h,o-z)
      IF (N.ge.2) go to 1
      if(n.eq.0) bessk = bessk0(x)
      if(n.eq.1) bessk = bessk1(x)
      return
1     continue
      TOX=2.0/X
      BKM=BESSK0(X)
      BK=BESSK1(X)
      DO 11 J=1,N-1
        BKP=BKM+J*TOX*BK
        BKM=BK
        BK=BKP
11    CONTINUE
      BESSK=BK
      RETURN
      END

      FUNCTION BESSK0(X)
      implicit real*8 (a-h,o-z)
      DATA P1,P2,P3,P4,P5,P6,P7/-0.57721566D0,0.42278420D0,0.23069756D0,&
          0.3488590D-1,0.262698D-2,0.10750D-3,0.74D-5/
      DATA Q1,Q2,Q3,Q4,Q5,Q6,Q7/1.25331414D0,-0.7832358D-1,0.2189568D-1,&
          -0.1062446D-1,0.587872D-2,-0.251540D-2,0.53208D-3/
      IF (X.LE.2.0d0) THEN
        Y=X*X/4.0d0
        BESSK0=(-dLOG(X/2.0d0)*BESSI0(X))+(P1+Y*(P2+Y*(P3+&
              Y*(P4+Y*(P5+Y*(P6+Y*P7))))))
      ELSE
        Y=(2.0/X)
        BESSK0=(dEXP(-X)/dSQRT(X))*(Q1+Y*(Q2+Y*(Q3+&
              Y*(Q4+Y*(Q5+Y*(Q6+Y*Q7))))))
      ENDIF
      RETURN
      END

      FUNCTION BESSK1(X)
      implicit real*8 (a-h,o-z)
      DATA P1,P2,P3,P4,P5,P6,P7/1.0D0,0.15443144D0,-0.67278579D0,&
          -0.18156897D0,-0.1919402D-1,-0.110404D-2,-0.4686D-4/
      DATA Q1,Q2,Q3,Q4,Q5,Q6,Q7/1.25331414D0,0.23498619D0,-0.3655620D-1,&
          0.1504268D-1,-0.780353D-2,0.325614D-2,-0.68245D-3/
      IF (X.LE.2.0d0) THEN
        Y=X*X/4.0d0
        BESSK1=(dLOG(X/2.0)*BESSI1(X))+(1.0/X)*(P1+Y*(P2+&
            Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7))))))
      ELSE
        Y=2.0/X
        BESSK1=(dEXP(-X)/dSQRT(X))*(Q1+Y*(Q2+Y*(Q3+&
            Y*(Q4+Y*(Q5+Y*(Q6+Y*Q7))))))
      ENDIF
      RETURN
      END

      FUNCTION BESSI0(X)
      implicit real*8 (a-h,o-z)
      DATA P1,P2,P3,P4,P5,P6,P7/1.0D0,3.5156229D0,3.0899424D0,1.2067492D-1,&
          0.2659732D0,0.360768D-1,0.45813D-2/
      DATA Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9/0.39894228D0,0.1328592D-1,&
          0.225319D-2,-0.157565D-2,0.916281D-2,-0.2057706D-1,&
          0.2635537D-1,-0.1647633D-1,0.392377D-2/
      IF (dABS(X).LT.3.75) THEN
        Y=(X/3.75)**2
        BESSI0=P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7)))))
      ELSE
        AX=ABS(X)
        Y=3.75/AX
        BESSI0=(dEXP(AX)/dSQRT(AX))*(Q1+Y*(Q2+Y*(Q3+Y*(Q4&
            +Y*(Q5+Y*(Q6+Y*(Q7+Y*(Q8+Y*Q9))))))))
      ENDIF
      RETURN
      END

      FUNCTION BESSI1(X)
      implicit real*8 (a-h,o-z)
      DATA P1,P2,P3,P4,P5,P6,P7/0.5D0,0.87890594D0,0.51498869D0,&
          0.15084934D0,0.2658733D-1,0.301532D-2,0.32411D-3/
      DATA Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9/0.39894228D0,-0.3988024D-1,&
          -0.362018D-2,0.163801D-2,-0.1031555D-1,0.2282967D-1,&
          -0.2895312D-1,0.1787654D-1,-0.420059D-2/
      IF (dABS(X).LT.3.75) THEN
        Y=(X/3.75)**2
        BESSI1=X*(P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7))))))
      ELSE
        AX=dABS(X)
        Y=3.75/AX
        BESSI1=(dEXP(AX)/dSQRT(AX))*(Q1+Y*(Q2+Y*(Q3+Y*(Q4+&
            Y*(Q5+Y*(Q6+Y*(Q7+Y*(Q8+Y*Q9))))))))
      ENDIF
      RETURN
      END
!
!
!
      subroutine quagen(z,wt,n)     
      integer n
      real*8 z(n),wt(n)
      real*8 p(513),c1(513),c2(513)
      real*8 pi,pi4,pdir,xnow,const,epsilon
      logical found
      data p(1),epsilon,nmax,ip /1.d0,1.d-10,10,512/

      if(n.le.0.or.n.gt.ip) then
        write(*,*) 'input order out of range'
       stop
      endif
!
!     calculate the coefficients for the legendre poly. recursive formula
!   
      do 10 i=2,n
      c1(i)=(2*i-1)/dfloat(i)
      c2(i)=(i-1)/dfloat(i)
 10   continue

!
!     initial constants
!
      n1=n+1
      pi=acos(-1.d0)
      pi4=pi/4
      const=1.d0/(n+0.5)
!
!     determine the number of roots(nr) needed to be calculated
!
      n2=n/2
      if(n2*2.eq.n) then
       ncal=n2
      else
        ncal=n2+1
      endif
! 
!     main loop begins here
!   
      do i=1,ncal
       k=n-i+1
       ncount=0
       xinc=1.d0
!
!     use newton's method and a good initial guess to locate the root
!
       xnow=cos((i*pi-pi4)*const)
       found=.false.
      do while (.not.found)
        ncount=ncount+1
         p(2)=xnow
!
!     the following loop calculate p_n(x) using recursive formula
!
       do 20 j=2,n
        p(j+1)=c1(j)*xnow*p(j)-c2(j)*p(j-1) 
 20    continue

!
!     the derivate of p_n(x) can be calculated from p_n(x) and p_n-1(x)
!
       pdir=n*(p(n)-xnow*p(n1))/(1.d0-xnow*xnow)
        if(abs(xinc).le.epsilon.or.ncount.gt.nmax) then
      found=.true.
        z(k)=xnow
        z(i)=-xnow
       wt(k)=2.d0/(1.d0-z(k)*z(k))/(pdir*pdir)
       wt(i)=wt(k)
      endif
      xinc=-p(n1)/pdir
      xnow=xnow+xinc
      enddo
      enddo
      return
      end

!  /////////mao

!********************************************************************
!	 subroutine CALCULATE SOIL DIELECTRIC CONSTANT	OF SOIL		  *
!********************************************************************


      subroutine mois2dc(mv,ff,er,temperature)
	implicit real*4(a-h,o-z)
	implicit byte(n)
	integer i,j
	real*4 idcw,idcs,idcw_2,idcs_2,mv,mv2,ff
	complex*8 dcsoil,dcsoil_2
	complex*16 er
!	ff=fre
	rho_s=2.66
	dcsr=4.7	
!-----------------------------------------------------------c
!          input   Parameters                               c
!-----------------------------------------------------------c
!------------------------------------------------------------------
!**** input soil parameters:soil bulk density rho_b;volum fraction
!**** of sand sv; volum fraction of clay cv; soil temprature st
!
!------------------------------------------------------------------
	rho_b = 1.3
!-----------------------------------------------
!---- sand and clay volum fraction (eg:sand=0.50,clay=0.30)
!-----------------------------------------------
!	cv = 0.22
!	sv = 0.20

    cv = 0.42
	sv = 0.20
	st=temperature
!-----------------------------------------------	
!	st = 24            //mao revision
	
!	do 1998 mv1=0.05, 0.70, 0.01
!	do 1999 mv2=0.02, 0.44, 0.02
!	wvfrac water volume fraction
	wvfrac_1=mv
	wvfrac_2=mv

	fv=rho_b/rho_s	
	ewin = 4.9
	alpha=0.65

!----------------------------------------------------------

!	 calculate texture dependent parameters

  

	betar = 1.2748 - 0.519 * sv - 0.152 * cv
	betai = 1.33797 - 0.603 * sv - 0.166 * cv

	if (ff>=1.4) then
		sigmae=-1.645 + 1.939 * fv * 2.66 - 2.25622 * sv + 1.594 * cv	
	else	
		sigmae=0.0467 + 0.2204 * fv * 2.66 - 0.4111 * sv + 0.6614 *cv	
	endif
	
!---------------------------------------------------------------------- 
!	 read input for volume fraction of water
!------------------------------------------------------------------------

!	 calculate temperature dependence of static dielectric constant of pure
!	 water (dimensionless)

!--------------------------------------------------------------------------

		ew0 = 88.045 - 0.4147 * st + 6.295 * 0.0001 * st * st&
     	 	  + 1.075 * 0.00001 * st * st * st

!--------------------------------------------------------------------------
!	 calculate temperature dependence of relaxation time of pure water
!	 where ff in frequncy in GHz
!--------------------------------------------------------------------------

		tauf = ff * (0.11109 - 0.003824 * st + 0.00006938 * st * st&
    		   - 5.096 * 0.0000001 * st * st * st)

!------------------------------------------------------------------------- 
!	 calculate real part of dielectric constant of water
!-------------------------------------------------------------------------

		rdcw = ewin + ( ew0 - ewin ) / ( 1.0 + tauf * tauf )

!--------------------------------------------------------------------------- 
!	 calculate imag part of dielectric constant of water in temporal 1
!--------------------------------------------------------------------------

		idcw = tauf * ( ew0 - ewin ) / ( 1.0 + tauf * tauf )&
     			+ sigmae * ( 1.0 - fv ) / (8.0 * atan(1.0) * 8.854 *&
     			0.001 * ff * wvfrac_1)
		tmp_1 = 1.0 + fv * ((dcsr**alpha) - 1.0 )&
     			+ (wvfrac_1**betar) * (rdcw**alpha) - wvfrac_1

		rdcs = tmp_1**(1.0/alpha)
		idcs = wvfrac_1**(betai/alpha) * idcw

		if(wvfrac_1>0) then
			dcsoil = Cmplx(rdcs,idcs)		
		else
			dcsoil=1
		end if
!--------------------------------------------------------------------------- 
!	 calculate imag part of dielectric constant of water in temporal 2
!--------------------------------------------------------------------------

		idcw_2 = tauf * ( ew0 - ewin ) / ( 1.0 + tauf * tauf )&
     			+ sigmae * ( 1.0 - fv ) / (8.0 * atan(1.0) * 8.854 *&
     			0.001 * ff * wvfrac_2)

		tmp_2 = 1.0 + fv * ((dcsr**alpha) - 1.0 )&
     			+ (wvfrac_2**betar) * (rdcw**alpha) - wvfrac_2
		rdcs_2 = tmp_2**(1.0/alpha)
		idcs_2 = wvfrac_2**(betai/alpha) * idcw_2

		if(wvfrac_2>0) then
			dcsoil_2 = Cmplx(rdcs_2,idcs_2)		
		else
			dcsoil_2=1
		end if

		reflect_1=abs((csqrt(dcsoil)-1)/(csqrt(dcsoil)+1))**2
		reflect_2=abs((csqrt(dcsoil_2)-1)/(csqrt(dcsoil_2)+1))**2
		Ratio_mv=wvfrac_1/wvfrac_2

!		Ratio_ref=reflect_1/reflect_2
		er = dcsoil
	RETURN
	end
