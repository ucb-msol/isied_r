%%%% Matlab Prototype of Depsoition code for analysis %%%%

close all
clear
clc
format long
rmdir('NewOutputPos','s')
mkdir('./','NewOutputPos')
% fileName = fopen('./NewOutputPos/tDat.dat','w');

% Np = 10*randi(1) + 20; % number of particle1s

Np = 100; % number of particle1s

xres = 20; % Bin spacing in x direction
yres = 10; % Bin spacing in y direction
binBorder = [1:yres, xres*yres-yres:xres*yres, yres+1:xres*yres-(yres+1):yres, 2*yres:xres*yres-(yres):yres];
dt = 0.000025; % Initial Time step (s)
t = 0; % Initial time (s)
tf = 0.4; % Final Time i.e. time span of simulation (s)
Nt = tf/dt;
TOL = 5E-4; % fixed point iteration tolerance
Kd = 10; % Max fixed point iteration count
dtmax = 0.00025;
p = 2;
phi = 0.3;
a = 0;

scale = 1E4;
% 
Kpo = 5E8/scale; % Normal contact parameter
Kwo = 5E8/scale; % Normal wall contact parameter

ccd = 5E8/scale; % Contact damping parameter
Kf = 1E6/scale; % Friction contact parameter

muS = 1; % Coefficent of static friction

Knb = 5E8/scale; % Normal bond parameter
pb = 2; % binding law exponent

Krb = 5E8/scale; % Rotational bond parameter

% Kpo = 1E7; % Normal contact parameter
% Kwo = 1E9; % Normal wall contact parameter
% 
% ccd = 1E5; % Contact damping parameter
% Kf = 1E7; % Friction contact parameter
% 
% muS = 0.4; % Coefficent of static friction
% 
% Knb = 1E6; % Normal bond parameter
% pb = 2; % binding law exponent
% 
% Krb = 1E3; % Rotational bond parameter
% 
epsStar = 0.8; % deformation metric threshold

a1 = 0.5; b1 = 1; a2 = 0.01; b2 = 2; % nearfield parameters

ce = 1; % Interstitial damping coefficient

% Defining Constant Material parameters (316L stainless steel)
E = (193E9); % Young's modulus (N/m^2)
nu = 0.26; % Poisson ratio (unitless)
Estar = (E)/(2* (1-(nu^2))); % Effective modulus
zeta = 3000; %Damping parameter (unitless)
muD = 0.3; % Dynamic Friction Coefficient  (unitless)
muF = 2.1E-5; % Viscosity of Argon (N/m^2)-s
rho0 = 2000; %(kg/m^3) density of steel particles
gravity = 9.81; % Gravitational constant (m/s^2)

% dRange = 4.24661E-6; % Range of diameters (m)
% dMin = 26.5E-6; % Minimum powder diameter

dRange = 0.03; % Range of diameters (m)
dMin = 0.1; % Minimum powder diameter

rp = (dRange * rand(Np,1) + dMin)/2; % radii of all powder particles

rij = rp + rp';

Mass = (1.33)*pi*(rp.^3.0)*rho0; % mass of powder particles

Rstar = zeros(Np,Np);
Mstar = zeros(Np,Np);

for i = 1:Np
    Rstar(i,:) = (rp(i)*rp)./(rp(i)+rp);
    Mstar(i,:) = (Mass(i)*Mass)./(Mass(i)+Mass);
end

rMax = max(rp); % Max Radius of all powder


% LxRange = (1.0* 1.0E-3)-(rMax); % Powder bed range in X
% LyRange = (1.0* 5E-4)-(rMax); % Powder bed range in Y

LxRange = (0.5)-(rMax); % Powder bed range in X
LyRange = (0.5)-(rMax); % Powder bed range in Y
LxMin = rMax; LyMin = rMax; % Initial point of bed in x-y plane

% LzRange = (1E-4); % Range fo particle positions in z
% LzMin = 2E-4; % Minimum height of particles above print bed

LzRange = 0.5; % Range fo particle positions in z
LzMin = 2*rMax; % Minimum height of particles above print bed

Lx = LxMin+LxRange; % Length of bed in x
Ly = LyMin+LyRange; % Length of bed in y
Lz = LzMin + LzRange; % Length of bed in z

hx = Lx/xres; hy = Ly/yres; % Bin spacing in x and y

%% Initializing powder positions and velocities
  
% Xp = [0.9*LxRange*rand(Np,1) + 1.1*LxMin, 0.9*LyRange*rand(Np,1) + 1.1*LyMin, LzRange*rand(Np,1) + LzMin]; % Initial Powder positions

[Xp1,Xp2,Xp3] = meshgrid(linspace(1.2*LxMin,0.8*Lx,4)',linspace(1.2*LyMin,0.8*Ly,5)',linspace(1.2*LzMin,0.6*Lz,5)'); % Initial Powder positions

Xp = [Xp1(:),Xp2(:),Xp3(:)]; 

    
Xp = [zeros(Np,3); Xp]; % Appending velocities to positions
%% Main Time loop

% [ePts] = rebin(xres,Xp,hx,hy,Np); % Bin particles initially

maxV = 1;
ek = 1;
uguess = zeros(2*Np,2);
bigRe = [];
bigdt = [];

c = 0;

tic
while t < tf %maxV > 4E-4

ce = 1; % Interstitial damping coefficient
    
    [fxLast, mRe] = fxt2(Xp,rp,Np,Estar,zeta,muD,Mstar,Rstar,gravity,Lx,Ly,Lz,Mass,rho0,muF,rij,Kpo,Kwo,ccd,Kf,muS,Knb,pb,Krb,a1,b1,a2,b2,dt,epsStar);
    
    bigRe = [bigRe, mRe];
    
    Xtest = Xp;
    Zk = 10;
    while Zk > 1
        for K = 0:Kd
            Xtestnew = Xp + dt*(phi*fxLast + (1-phi)*fxt2(Xp,rp,Np,Estar,zeta,muD,Mstar,Rstar,gravity,Lx,Ly,Lz,Mass,rho0,muF,rij,Kpo,Kwo,ccd,Kf,muS,Knb,pb,Krb,a1,b1,a2,b2,dt,epsStar));
            
            if K == 0
                Xtest = Xtestnew;
            elseif K == 1
                w0 = sum(vecnorm(Xtestnew-Xp,2,2));
                wk = w0;
                Zk = wk/TOL;
                Lambdak = ((TOL/w0).^(1/(p*Kd)))./((wk/w0).^(1.0/(p*K)));
                Xtest = Xtestnew;
            else
                wk = (sum(vecnorm(Xtestnew-Xtest,2,2)))/(sum(vecnorm(Xtestnew-Xp,2,2)));
                Zk = wk/TOL;
                Lambdak = ((TOL/w0).^(1.0/(p*Kd)))/((wk/w0).^(1.0/(p*K)));
                Xtest = Xtestnew;
            end
            
            if (Zk <= 1 && K < Kd)
                t = t+dt;
                dt = Lambdak*dt;
                dt = min(dt,dtmax);
                break;
            elseif (Zk > 1 && K == Kd)
                dt = Lambdak*dt;
            elseif (isnan(Zk))
                dt = 0.1*dt;
                Zk = 10;
                if dt < 1e-10
                    break;
                end
    
                break;
            end
            if dt < 1e-10
                    break;
            end
        end
        if dt < 1e-10
                    break;
        end
    end

    if dt < 1e-10
        break;
    end
    
    Xp = Xtestnew; %Updating Xp
    
%     Xp(Np+find(Xp(Np+1:end,1)<= 0),1) = rp(Xp(Np+1:end,1)<= 0);
%     Xp(Xp(Np+1:end,1)<= 0,1) = 0;
%     
%     Xp(Np+find(Xp(Np+1:end,1)>= Lx),1) = Lx - rp(Xp(Np+1:end,1)>= Lx);
%     Xp(Xp(Np+1:end,1)>= Lx,1) = 0;
%     
%     Xp(Np+find(Xp(Np+1:end,2)<= 0),2) = rp(Xp(Np+1:end,2)<= 0);
%     Xp(Xp(Np+1:end,2)<= 0,2) = 0;
%     
%     Xp(Np+find(Xp(Np+1:end,2)>= Ly),2) = Ly - rp(Xp(Np+1:end,2)>= Ly);
%     Xp(Xp(Np+1:end,2)>= Ly,2) = 0;
%         
%     Xp(Np+find(Xp(Np+1:end,3)<= 0),3) = rp(Xp(Np+1:end,3)<= 0);

%       Xp(Xp(Np+1:end,3)<= 0,3) = 0.5*Xp(Xp(Np+1:end,3)<= 0,3);

%     
%     Xp(Xp(Np+1:end,3)>= Lz,3) = 0;
    
    if t > a
                
        if a == 0
            fileID = fopen('./NOV22/finalOut.plt','w');
            fprintf(fileID,'%12s\n','VARIABLES = "X" "Y" "Z" "U" "V" "W" "r"');

        else
            fileID = fopen('./NOV22/finalOut.plt','a');
        end
        
        dataLine = [Xp(Np+1:end,1), Xp(Np+1:end,2),Xp(Np+1:end,3),Xp(1:Np,1), Xp(1:Np,2),Xp(1:Np,3),rp]';


        fprintf(fileID,'%12s\n%12s\n','ZONE',...
                ['DATAPACKING = POINT, SOLUTIONTIME = ' num2str(t)]);

        

        fprintf(fileID,'%12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f \n',dataLine);
        fclose(fileID);
    
%         iter = iter + 1;

        a = a+2E-3;
    end

     
    bigdt = [bigdt dt];
    
    t
    c = c+1;
%     maxV = max(vecnorm(Xp(1:Np),2,2))
end
toc
% fclose(fileName);
% 
close all
% plot3(Xp(Np+1:end,1),Xp(Np+1:end,2),Xp(Np+1:end,3),'.','MarkerSize',rp)
scatter3(Xp(Np+1:end,1),Xp(Np+1:end,2),Xp(Np+1:end,3),1000*rp)
% for i = 1:size(bigX,2)
% plot3(bigX(:,i),bigY(:,i),bigZ(:,i),'*')
% title(['t  = ' num2str(bigT(i))]);
% view(0,0)
% xlim([LxMin,LxMin+LxRange]);
% ylim([LyMin,LyMin+LyRange]);
% zlim([0,LzMin+LzRange]);
% pause(0.2)
% end
