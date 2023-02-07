!
!   ParticleRandomizer.f90
!   BaseLineCode
!
!   Created by Roger Isied on 6/8/18.
!   Copyright 2018 RogerIsied. All rights reserved.
!

program ParticleRandomizer

USE IFPOSIX
implicit none

integer :: i, s, j, Np, Nt, err, bin(50), xres, yres, xbin, ybin, Kd, K
real :: rp(50), dp(50), t, i1, rSeedx(4), rSeedy(4), rSeedz(4), rSeedr(4), tf, Mr, ek, p, start, finish, phi
integer, allocatable :: old(:), binBorder(:)
real, dimension(2*50,3) :: X, fxLast, Xtest, Xtestnew
real :: dz, dt, rem, TOL, w0, wk, Zk, dtmax, Lambdak
character(len=8) :: fmt, tw, file ! format descriptor
character(*), PARAMETER :: fileplace = '/Users/ruisied/Desktop/PHD/Research/BaseLineCode/OutputPos/'
character :: str, status
call cpu_time(start)
fmt = '(I7.7)' ! an integer of width 5 with zeros at the left

!str = 'delete'//fileplace
!call system(str)

!call PXFRMDIR(fileplace,1,err)
!print*, err
!CALL PXFMKDIR(fileplace,50,0,err)
!pause

! Defining simulation parameters
Np = 50 ! number of particles

xres = 10 ! Bin spacing in x direction
yres = 5 ! Bin spacing in y direction
binBorder = (/1:yres, xres*yres-yres:xres*yres, yres+1:xres*yres-(yres+1):yres, 2*yres:xres*yres-(yres):yres/)
!print*, binBorder
!pause
t = 0 ! Initial time (s)
dt = 1.0E-6 ! Initial Time step (s)
tf = 0.1 ! Final Time i.e. time span of simulation (s)
Nt = tf/dt
TOL = 1E-3 ! fixed point iteration tolerance
Kd = 10 ! Max fixed point iteration count
dtmax = 1.0E-4
p = 2
phi = 0
! Create random seeds for positions and particle radii
!rSeedx = (/1, 4, 6, 9/)
!rSeedy = (/2, 5, 10, 15/)
!rSeedz = (/2, 6, 19, 21/)
!rSeedr = (/10, 11, 19, 421/)

call random_seed()

! Create random unifrom positions and gaussian distributed radii
call random_number(X(Np+1:,:))
allocate (old(4))

call random_seed(Get = old(1:4))

call slarnv(3,old,Np,rp)
!call random_number(rp)

! Assign initial velocities
X(1:Np:1,:) = 0

! Scale range and std to match position of that in (Ganeriwala, Zohdi)
dp = ((4.24661E-6)*rp + (26.5E-6))
rp = dp/2.0
!print*, rp
!pause
Mr = maxval(rp)
X(Np+1:,1) = ((1.0* 1.0E-3)-(1.0* Mr))*X(Np+1:,1) + (Mr)
X(Np+1:,2) = ((1.0* 5E-4)-(1.0* Mr))*X(Np+1:,2) + (Mr)
X(Np+1:,3) = ((1E-4))*X(Np+1:,3) + (2E-4)

    write(tw, fmt) s


    open(s, status="unknown", file = fileplace//'t_'//'0000000'//'.txt', action="write")
    write(s,*) 'x,', 'y,', 'z,', 'd'!, 't'

    do i = 1,Np
            xbin = floor(X(Np+i,1)*(1E4)) + 1
            ybin = floor(X(Np+i,2)*(1E4)) + 1
            bin(i) = yres*(xbin-1) + ybin
            write(s,*) X(Np+i,1), ',', X(Np+i,2), ',', X(Np+i,3), ',', dp(i)!, ',', t
    end do
    close(s)

s = 1
do while (t < tf)

    ! Applying fixed point iteration for trapezoidal method
    !uguess = X
    fxLast = fxt1(X,rp,Np,bin,binBorder)

!    do while (ek > 1E-12)
!        ustar = X + 0.5*dt*(fxLast + fxt1(uguess,rp,Np,bin,binBorder))
!        ek = NORM2(ustar-uguess)
!        uguess = ustar
!    end do
    Xtest = X
    Zk = 10
    do while (Zk > 1)
        do K = 0,Kd
            Xtestnew = X + dt*(phi*fxlast + (1.0-phi)*fxt1(Xtest,rp,Np,bin,binBorder))
            if (K == 0) then
                Xtest = Xtestnew
            else if (K == 1) then
                w0 = sum(abs((Xtestnew-X)))
                wk = w0
                Zk = wk/TOL
                Lambdak = ((TOL/w0)**(1.0/(p*Kd)))/((wk/w0)**(1.0/(p*K)))
                Xtest = Xtestnew
            else
                wk = (sum(abs(Xtestnew-Xtest)))/(sum(abs(Xtestnew-X)))
                Zk = wk/TOL
                Lambdak = ((TOL/w0)**(1.0/(p*Kd)))/((wk/w0)**(1.0/(p*K)))
                Xtest = Xtestnew
                !print*, ((TOL/w0)**(1/(p*Kd))), ((wk/w0)**(1/(p*k)))
            end if
            if (Zk <= 1 .AND. K < Kd) then
                t = t + dt
                dt = Lambdak*dt
                dt = MINVAL((/dt,dtmax/))
                EXIT
            else if (Zk > 1 .AND. K == Kd) then
                dt = Lambdak*dt
                !dt = MAXVAL((/dt,dtmin/))
            else if (isnan(Zk)) then
                dt = 0.2*dt
                Zk = 10
            end if
        end do
    end do

    X = Xtestnew
    do i = 1,Np
            if (X(i+Np,1) <= 0.0) then
                X(i+Np,1) = rp(i)
                X(i,1) = 0.0
            else if (X(i+Np,1) >= (1.0*1.0E-3)) then
                X(i+Np,1) = 1.0*1.0E-3 - rp(i)
                X(i,1) = 0.0
            end if
            if (X(i+Np,2) <= 0.0) then
                X(i+Np,2) = rp(i)
                X(i,2) = 0.0
            else if (X(i+Np,2) >= (1.0*5.0E-4)) then
                X(i+Np,2) = 1.0*   5.0E-4 - rp(i)
                X(i,2) = 0.0
            end if
            if (X(i+Np,3) <= 0.0) then
                X(i+Np,3) = rp(i)
                !X(i,3) = 0.0
            else if (X(i+Np,3) >= 3.1E-4) then
                X(i,3) = 0.0
            end if
    end do
    !X = X + 0.5*dt*(fxLast + fxt1(uguess,rp,Np,bin,binBorder))

    !t = t + dt
    ! March forward with fourth order explicit Runge-Kutta scheme
    !Y1 = fxt1(X,rp,Np)
    !Y2 = X + 0.5*dt*fxt1(Y1,rp,Np)
    !Y3 = X + 0.5*dt*fxt1(Y2,rp,Np)
    !Y4 = X + dt*fxt1(Y3,rp,Np)
        !print*, Y1(Np+1,:), Y2(Np+1,:), Y3(Np+1,:), Y4(Np+1,:)
    !pause
    !FY2 = (1/6)*dt*(fxt1(Y1,rp)+(2*fxt1(Y2,rp))+ (2*fxt1(Y3,rp)) + fxt1(Y4,rp))
    !X = X + (dt/6.0)*(fxt1(Y1,rp,Np)+(2*fxt1(Y2,rp,Np))+ (1.0* fxt1(Y3,rp,Np)) + fxt1(Y4,rp,Np))
    !X = X + (dt/6.0)*(Y1+2*Y2+2*Y3+Y4)

        !print*, X(Np+1,:)
        !pause
    !print*, X(1,3)
    !pause
    write(tw, fmt) s
    rem = modulo(s,50)
    if (rem == 0) then
        print*, dt
        open(s, status="unknown", file = fileplace//'t_'//trim(tw)//'.txt', action="write")
        write(s,*) 'x,', 'y,', 'z,', 'd'!, 't'

        do i = 1,Np
            xbin = floor(X(Np+i,1)*(1E4)) + 1
            ybin = floor(X(Np+i,2)*(1E4)) + 1
            bin(i) = yres*(xbin-1) + ybin
            write(s,*) X(Np+i,1), ',', X(Np+i,2), ',', X(Np+i,3), ',', dp(i)!, ',', t
        end do
        close(s)
    end if


    s = s+1
end do

call cpu_time(finish)
print*, finish
contains



function fxt1(X1,rp1,Np1,bin1,binBorder1) result(Xdot1)

implicit none



real :: E, nu, zeta, zeta2, muD, muF, h, Tm, Tv, Lm, Lv, eps, alpha, gamma, Power, vlaser
real, dimension(10) :: Mass(Np1) !Tdat(10), Cdat(10), kDat(10), rhoDat(10), Tint, Cint, kint, rhoInt
real :: Cextrap, kextrap, rhoextrap, Estar1, Rstar1, Rstar2, Pi = 3.1415927, Pdist, rho0, omega,T0, Tenv
real :: deltaij, deltaDotij, deltaiw, deltaDotiw, damp, dampW, gravity, mstar1
real, dimension(3) :: vecij, veciw, nij, xw, niw, tij, tiw, vij, vtj, vti, vtw
real, dimension(1,3) :: Fconij, Fconiw, Ffricij, Ffriciw, Fenv, Fgrav, Ftot1, Ftot2
integer :: i, j, l, mybin(Np1)
integer, intent(IN) :: Np1, bin1(Np1), binBorder1(:)
integer, allocatable :: binPart(:)
logical :: wCon
real :: rp1(Np1), X1(2*Np1,3)
real, target, dimension(2*Np1,3) :: Xdot1
!real, pointer, dimension(:,:) :: fxt1
! Defining Constant Material parameters (316L stainless steel)
E = (193.0E9) ! Young's modulus (N/m^2)
nu = 0.26 ! Poisson ratio (unitless)
zeta = 3000 !Damping parameter (unitless)
zeta2 = 3000
muD = 0.1 ! Dynamic Friction Coefficient  (unitless)
muF = 2.1E-5 ! Viscosity of Argon (N/m^2)-s
!h = 40E-12 ! Heat Transfer Coefficient (W/um^2-K)
!Tm = 1700 !Melting T (K)
!Tv = 3130 ! Boiling T (K)
!Lm = 2.99E5 ! Latent heat of melting J/kg
!Lv = 6.09E6 ! Latent heat of vaporization J/kg
!eps = 0.33 ! Material Emissivity (unitless)
!alpha = 0.33 ! Material absorptivity (unitless)
!gamma = 0.55 ! Power bed porosity (unitless)
!Power = 200 ! Laser Power (W)
!vlaser = 2.0 ! Laser Scan Speed (m/s)
!omega = 54E-6 ! Laser Spot Size (m)
!T0 = 363 ! Preheat Temperature (K)
!Tenv = 363 ! Environment Temperature (Tenv)
rho0 = 7919 !(kg/m^3) density of steel particles
gravity = 9.81 ! Gravitational constant (m/s^2)
! Defining Temperature Dependent Material parameters (316L stainless steel)

!Tdat = (/293, 366, 478, 589, 700, 811, 922, 1033, 1144, 1255/) !Measured Temp (K)
!Cdat = (/452, 485, 527, 548, 565, 573, 586, 615, 649, 690/) !Measured Specific Heat (J/kg-K)
!kDat = (/13.3, 14.3, 15.9, 17.5, 19.0, 19.8, 21.9, 23.2, 24.6, 26.2/)! Measured Thermal Conductivity (W/m-K)
!rhoDat = (/7952, 7919, 7877, 7831, 7786, 7739, 7692, 7640, 7587, 7537/) ! Density (kg/m^3)
!Cextrap = (690-645)/(1255-1144)
!kextrap = (((32.4E-6)-(26.2E-6))/(1255-1144))
!rhoextrap = (((7300E-18)-(7537E-18))/(1255-1144))*(10**-18)


! Initialize Xdot
!Xdot1=>X1

Xdot1(1:Np1,:) = 0
Xdot1(Np+1:,:) = X1(1:Np1,:)

!Define effective parameters not dependent on particle geometry and position
Estar1 = E/(1.0* (1-(nu**2))) !Effective modulus between two powder particles

Mass = (1.33)*Pi*(rp1**3.0)*rho0


do i = 1,Np1
mybin(:) = bin1(i)
binPart = pack([(l,l=1,size(bin1))],mybin == bin1)
!print*, binPart
!pause
        do l = 1,size(binPart)
        j = binPart(l)
            if (i == j) then
                cycle
            end if
            vecij(1) = X1(j+Np1,1)-X1(i+Np1,1)
            vecij(2) = X1(j+Np1,2)-X1(i+Np1,2)
            vecij(3) = X1(j+Np1,3)-X1(i+Np1,3)
            Pdist = NORM2(vecij)
            !print*, rp1(i)+rp1(j)
            !print*, Pdist
            !pause
            if (Pdist <= (rp1(i)+rp1(j))) then ! Determine Particle Contact

            !print*, rp1(i)+rp1(j)
            !print*, Pdist
            !pause
            ! Calculate particle to particle contact force
                Rstar1 = (rp1(i)*rp1(j))/(rp1(i)+rp1(j)) !effective radius between two particles
                deltaij = abs(Pdist - (rp1(i)+rp1(j)))
                nij = vecij/Pdist
                vij(1) = X1(j,1)-X1(i,1)
                vij(2) = X1(j,2)-X1(i,2)
                vij(3) = X1(j,3)-X1(i,3)

                !print*, vij
                deltaDotij = dot_product(vij,nij)
                mstar1 = (Mass(i)*Mass(j)) / (Mass(i)+Mass(j))
                damp = 1.0* zeta2*(deltaij**(0.25))*sqrt(1.0* Estar1*mstar1*sqrt(Rstar1))
                if (deltaDotij == 0) then
                    Fconij(1,:) = (/0.0,0.0,0.0/)
                else
                    FConij(1,:) = (-1.33)*sqrt(Rstar1)*Estar1*(deltaij**(1.5))*nij + damp*deltaDotij*nij
                    Fconij(1,:) = Fconij(1,:)*(1.0E-10)
                end if
               ! print*, deltaij**(1.5), rp(i), rp(j)
                !pause
                !Xdot1=>FConij
                !Xdot1(i,:) = Xdot1(i,:) + Fconij(1,:)

                ! Calculate Particle Friction Force
                vtj = X1(j,:) - dot_product(X1(j,:),nij)*nij
                vti = X1(i,:) - dot_product(X1(i,:),nij)*nij
                if ((NORM2(vtj-vti)) == 0) then
                 tij = (/0.0,0.0,0.0/)
                else
                 tij = (vtj-vti)/(NORM2(vtj-vti))
                end if
                Ffricij(1,:) = muD * NORM2(FConij)*tij
                !print*, Ffricij
                !pause
                !Xdot1=>Ffricij
                !Xdot1(i,:) = Xdot1(i,:) + Ffricij(1,:)
                !pause
             else
                FConij(1,:) = (/0.0,0.0,0.0/)
                Ffricij(1,:) = (/0.0,0.0,0.0/)
                !pause
             end if
            Ftot1(1,:) =  (FConij(1,:) + Ffricij(1,:))
            Xdot1(i,:) = Xdot1(i,:) + Ftot1(1,:)

        end do
end do


! Determine particle wall contact
do i = 1,Np1
wCon = .false.
xw = X1(i+Np1,:)
    if (any(bin(i) == binBorder1)) then
             if (X1(i+Np1,1) <= rp1(i)) then
                xw(1) = 0.0
                wCon = .true.
             else if (X1(i+Np1,1) >= (1.0*   1.0E-3)-(rp1(i))) then
                xw(1) = 1.0*   1.0E-3
                wCon = .true.
            end if
            if (X1(i+Np1,2) <= rp1(i)) then
                xw(2) = 0.0
                wCon = .true.
            else if (X1(i+Np1,2) >= (1.0*   5.0E-4)-rp1(i)) then
                xw(2) = 1.0* 5.0E-4
                wCon = .true.
             end if
      end if
        ! Calculate environmental force
             Fenv(1,:) = 6.0*Pi*muF*rp1(i)*X1(i,:)
            ! print*, Fenv(1,:)
           ! pause
             !Xdot1=>Fenv
             !Xdot1(i,:) = Xdot1(i,:) + Fenv(1,:)

             ! Calculate gravitational gorce
             Fgrav(1,:) = -Mass(i)*gravity*(/0,0,1/)
             !if (X1(i+Np1,3) <= rp1(i)) then
                !Fgrav(1,:) = (/0,0,0/)
             !end if
             !print*, Fgrav(1,:)
             !Xdot1=>Fgrav
             !Xdot1(i,:) = Xdot1(i,:) + Fgrav(1,:)
            !print*, 'a' , Xdot1(1,:)
            !print*, 1/(Mass(i))
            !print*, Xdot1(1,:)
            !pause
            if (X1(i+Np1,3) <= rp1(i)) then
                xw(3) = 0.0
                wCon = .true.
            else if (X1(i+Np1,3) > (3.0E-3) - rp1(i)) then
                xw(3) = 3.0E-3
                wCon = .true.
            end if
            if (wCon == .true.) then
                !print*, xw
                !pause
                ! Calculate wall contact force
                niw = (xw-X1(i+Np1,:)) * (1/NORM2(xw-X1(i+Np1,:)))
                niw(3) = niw(3)
                !print*, NORM2(xw-X1(i+Np1,:))
                !pause
                deltaiw = abs(NORM2(xw-X1(i+Np1,:)) - rp1(i))
                !print*, deltaiw
                deltaDotiw = (dot_product((-X1(i,:)),niw))
                !print*, deltaDotiw
                Rstar2 = rp1(i)
                dampW = 2.0*zeta*(deltaiw**(0.25))*sqrt(2*Estar1*Mass(i)*sqrt(Rstar2))
                !print*, deltaDotiw
                !print*, dampW

                !print*, deltaiw
                !pause
                FConiw(1,:) = ((-1.33)*sqrt(Rstar2)*Estar1*(deltaiw**(1.5))*niw + dampW*deltaDotiw*niw)
                Fconiw = Fconiw*1E-9
                !print*, ((-1.33)*sqrt(Rstar2)*Estar1*(deltaiw**(1.5))*niw), dampW*deltaDotiw*niw
                !pause
                !Fconiw = abs(Fconiw)
                !print*, Fconiw
                !Xdot1=>FConiw
                !Xdot1(i,:) = Xdot1(i,:) + FConiw(1,:)

                 ! Calculate Particle wall Friction Force
                vtw = (/0,0,0/)
                vti = X1(i,:) - dot_product(X1(i,:),niw)*niw
                if (NORM2(vtw-vti) == 0.0) then
                    tiw = (/0.0,0.0,0.0/)
                else
                    tiw = (-vti)/NORM2(-vti)
                end if
                Ffriciw(1,:) = muD * NORM2(FConiw)*tiw
                !Xdot1=>Ffricij
                !Xdot1(i,:) = Xdot1(i,:) + Ffriciw(1,:)
            else
                    FConiw(1,:) = (/0.0,0.0,0.0/)
                    Ffriciw(1,:) = (/0.0,0.0,0.0/)
            end if
            Ftot2(1,:) =  Fgrav(1,:) + Fenv(1,:) + ((Fconiw(1,:) + Ffriciw(1,:)))
            Xdot1(i,:) = Xdot1(i,:) + Ftot2(1,:)
            Xdot1(i,:) = Xdot1(i,:)/Mass(i)
end do
!fxt1 => Xdot1
!print*, Xdot1(1,:)
!pause
!print*, Xdot1(1,:)
end function fxt1



end program ParticleRandomizer










