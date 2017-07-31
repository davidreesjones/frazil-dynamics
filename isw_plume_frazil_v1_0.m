% MIT License
% 
% Copyright (c) [2017] [David Rees Jones]
%     David Rees Jones [david.reesjones@earth.ox.ac.uk] 
% 	  University of Oxford 	  
% 	  Earth Sciences
% 	  South Parks Road
% 	  Oxford, OX1 3AN
% 	  United Kingdom
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE.

function [vdydt,Y,vx,vz,vU,vD,vT,vS,vCi,vCik,vdelT,vm,vfp,vfpk,vpp,vppk,vTf,vdrho] = isw_plume_frazil_v1_0(x0,xmax,y0,ri_mm,opt,tilde_nmax)
%isw_plume_frazil_v1_0 Crystal-laden ISW plume calculation
%   Version used in Cryosphere calculations
%   Includes conduction into ice shelf

% Input: [x0,xmax] distance from grounding line (in meters).
% Input: y0 = y(x=x0) initial conditions
% Input: ri_mm = discrete crystal sizes (in mm)
% Input: opt = growth law option (1 = SJ04, 2 = SO94)

% Components: y(1:4,:) = [D*U,D*U^2, D*U*T,D*U*S)
% Components: y(5:end,:) = [D*U*Ci(k))


%% Set up geometry and ambient statification
Lmax=600000; %m Grounding-line to ice-front distance
H=1400; %m Depth of grounding line
HL=285;
DH=HL-H;
S0=34.5; % ppt Linear salinity gradient
SH=34.71;
DS=SH-S0;
T0=-1.9; % C Linear temperature gradient
TH=-2.18;
DT=TH-T0;
T00=-2;  %C Reference temperature in equation of state
S00=34.5; % ppt Reference salinity in equation of state
Ts=-15; %C (core ice shelf temperature)

%% Model and material properties
E0=0.036; % Entrainment coefficient
theta=atan(-DH/Lmax); % Slope angle
g=9.81; 
Eth=E0*sin(theta);
gth=g*sin(theta);
UT=0; %m/s Tidal velocity
Cd=2.5e-3; % Drag coefficient   

% material properties
a=-0.0573; %C/ppt %Freezing point relations
b=0.0832; %C
c=-7.61e-4; %C/m (n.b. we define z positive downwards)
L=3.35e5; %J/kg
cw=3974; %J/kg/C
ci=2009; %J/kg/C (affects heat flux to ice shelf, set to zero to match JB95) 
nu=1.95e-6; %m^2/s
betaS=7.86e-4; %/ppt
betaT=3.87e-5; %/C
Pr=13.8; % Prandtl number
Sc=2432; % Schmidt number
rho_w=1028; %kg/m^3 1030 in JB95, 1028 in SJ04
rho_i=917; %kg/m^3 920 in JB95, 917 in SJ04
drhoi=(rho_w-rho_i)/(rho_w);
gp=g*drhoi; % reduced gravity
Nu=1;
KT=1.4e-7; %m^2/s thermal diffusivity
KS=8e-10; %m^2/s haline diffusivity

%% Set up frazil crystals
% properties
Nclass=length(ri_mm);
Vri=ri_mm*1e-3; % convert mm to m
if opt==1 
    ar=0.02; % aspect ratio (n.b. ar=0.02 in SJ)
    Vrie=Vri*(3*ar/2)^(1/3); % equivalent radius
    Vvi=2*pi*ar*Vri.^3; % volume
    Gc=(16*Vri.^3*gp*ar)/nu^2; 
end
if (opt==2)||(opt==3)
    ti0=0.05*1e-3;
    Vti=ti0+zeros(1,Nclass);
    Vrie=Vri.^(2/3).*(3*Vti/4).^(1/3); % equivalent radius
    Vvi=pi*Vri.^2.*Vti; % volume
    Gc=(8*Vri.^2.*Vti*gp)/nu^2; 
end
DVvi=diff(Vvi);
VRei=zeros(1,Nclass);
Vwi=zeros(1,Nclass);
for J=1:Nclass
    rt=roots([0.111,1.108,1.386-log10(Gc(J))]);
    if isreal(rt(1))==1
        VRei(J)=10^(max(rt));
        Vwi(J)=VRei(J)*nu/(2*Vri(J));
    else
        VRei(J)=0;
        Vwi(J)=0;
    end
end

%% set up crystal precipitation law
Vpp0=rho_i*cos(theta)*Vwi/rho_w;
Sh_crit=0.01; % Critical Shields Number - 0.05 (JB) or 0.1 (SJ)
VUC2=Sh_crit*gp*Vrie/Cd;

%% set up secondary nucleation law
eps_turb=7.4e-6; % W/kg
Ur=((eps_turb/(15*nu))*(2*Vrie).^2+Vwi.^2).^(0.5);
vUr=pi*Vrie(1)^3*Ur./Vrie;

%% MAIN Ode Routine
% Initial conditions are a user input, apart from constant seeding rate 
r_seed=0; %constant seeding rate 
opts = odeset('RelTol',1e-10,'AbsTol',1e-13,'NonNegative',5:4+Nclass); %ensures Ci>0
xmax=min(xmax,Lmax); % stops integrating beyond end of shelf
[X,Y]=ode15s(@ode,[x0 xmax],y0,opts); % Integrate to nucleation site
[vdydt, vz, vU, vD, vT, vS, vCi, vCik, vdelT, vm, vfp, vfpk, vpp, vppk, vTf, vdrho]=ode(transpose(X),transpose(Y)); % Compute auxilliary variables
vx=transpose(X);
save(strcat('Growth_',num2str(opt),'_tilde_nmax_',num2str(tilde_nmax,'%1.0e')));

    function    [dydt, z, U, D, T, S, Ci, Cik, delT, m, fp, fpk, pp, ppk, Tf, drho]=ode(x,y)
      % Main ode function
        LL=length(x);
        Cik=zeros(Nclass,LL);
        
        z=1400+DH*(x/Lmax); %z in metres (measured positively)        
        
        % CALCULATE MACROSCOPIC VARIABLES
        U=y(2,:)./(y(1,:));
        D=(y(1,:)).^2./y(2,:);
        T=y(3,:)./y(1,:);
        S=y(4,:)./y(1,:);
        for k=1:Nclass
            Cik(k,:)=y(4+k,:)./(y(1,:));
        end
        Ci=sum(Cik);
      
        % CALCULATE ENTRAINMENT
        e=Eth*U;
        
        % CALCULATE HEAT TRANSFER TO SHELF
        Ue=(U.^2+UT.^2).^0.5;
        gamT=Cd^0.5*Ue./(2.12*log(Cd^0.5*U.*D/nu)+12.5*Pr^(2/3)-9);
        gamS=Cd^0.5*Ue./(2.12*log(Cd^0.5*U.*D/nu)+12.5*Sc^(2/3)-9);
        gamR=gamS./gamT;     
        p2=-a+a*gamR*ci/cw;      
        p1=T-b-c*z+gamR.*(L/cw-(ci/cw)*(Ts-b-c*z+a*S));
        p0=-gamR.*S.*(L/cw-(ci/cw)*(Ts-b-c*z));        
        Sb1=(-p1+(p1.^2-4.*p0.*p2).^0.5)./(2.*p2);
        
        Tb1=a*Sb1+b+c*z;
        m1=gamS.*(S-Sb1)./(Sb1+eps);    
        delT1=T-Tb1;
        
        p2=-a;      
        p1=T-b-c*z+gamR.*(L/cw);
        p0=-gamR.*S.*(L/cw);        
        Sb2=(-p1+(p1.^2-4.*p0.*p2).^0.5)./(2.*p2);
        Tb2=a*Sb2+b+c*z;
        m2=gamS.*(S-Sb2)./(Sb2+eps);    
        delT2=T-Tb2;
        Tb=Tb1.*heavi(m1)+Tb2.*heavi(-m1);
        Sb=Sb1.*heavi(m1)+Sb2.*heavi(-m1);
        m=m1.*heavi(m1)+m2.*heavi(-m1);      
        delT=delT1.*heavi(m1)+delT2.*heavi(-m1);

        % CALCULATE BUOYANCY FORCING
        Ta=T0+(z+D*cos(theta))*DT/H;
        Sa=S0+(z+D*cos(theta))*DS/H;        
        Sab=(S0+DS*(z+D*cos(theta)/2)/H)/cos(theta);
        Tab=(T0+DT*(z+D*cos(theta)/2)/H)/cos(theta);
        drho=Ci*(1-rho_i/rho_w)+betaS*(Sab-S)-betaT*(Tab-T);

        % CALCULATE CRYSTAL GROWTH AND PRECIPITATION
        % freezing/melting
        Tf=a*S+b+c*(z+D*cos(theta)/2);
        tilde_n=min(sum(Cik./kron(Vvi',ones(1,LL))),tilde_nmax*ones(1,LL));
        dCdtkN=-(vUr'*tilde_n).*Cik;
        dCdtkN(1,:)=-sum(dCdtkN(2:end,:));
        dCdtkF=zeros(Nclass,LL);
        wpk=zeros(Nclass,LL);
        for k=1:Nclass
            if opt==1
            dCdtkF(k,:)=(cw*Nu*KT*(Tf-T)/L).*(1+heavi(T-Tf)/(2*ar)).*Cik(k,:)*(2/Vri(k)^2);
            end
            if opt==2
                dCdtkF(k,:)=(cw*Nu*KT*(Tf-T)/L).*(1).*Cik(k,:)*(2.0/(Vri(k)*Vti(k)));
            end
            if opt==3
                ar_h=Vti(k)/(2*Vri(k));
                fr_h=2/(0.9008-0.2634*log(ar_h));
                dCdtkF(k,:)=(cw*Nu*KT*(Tf-T)/L).*(1).*Cik(k,:)*(fr_h/(Vri(k)*Vti(k)));
            end
        end
        dCdtk=dCdtkF;
        k=1; wpk(k,:)=Vvi(k)*((dCdtk(k+1,:)/DVvi(k)).*heavi(T-Tf)+(dCdtk(k,:)/DVvi(k)-r_seed).*heavi(Tf-T));
        for k=2:Nclass-1           
            wpk(k,:)=Vvi(k)*((dCdtk(k+1,:)/DVvi(k)-dCdtk(k,:)/DVvi(k-1)).*heavi(T-Tf)+(dCdtk(k,:)/DVvi(k)-dCdtk(k-1,:)/DVvi(k-1)).*heavi(Tf-T));
        end
        k=Nclass; wpk(k,:)=Vvi(k)*((-dCdtk(k,:)/DVvi(k-1)).*heavi(T-Tf)+(-dCdtk(k-1,:)/DVvi(k-1)).*heavi(Tf-T));
        wp=sum(wpk);
        fp=D.*wp;
        for k=1:Nclass
            fpk(k,:)=D.*wpk(k,:)-D*dCdtkN(k);
        end 
             
        % precipiation
        for k=1:Nclass
            ppk(k,:)=-Vpp0(k).*Cik(k,:).*max(1-Ue.^2./VUC2(k),0);
        end
        
        pp=sum(ppk);

                
        % ODE 
        dydt(1,:)=e+m+fp;
        dydt(2,:)=D.*drho*gth-Cd*U.*Ue;               
        dydt(3,:)=e.*Ta+m.*Tb-gamT.*(T-Tb)-(L/cw-Tf).*fp;
        dydt(4,:)=e.*Sa;        
        dydt(5:4+Nclass,:)=(rho_w/rho_i).*(ppk-fpk);       
    end

    function hh=heavi(xx)
        % heaviside step function
        hh=(sign(xx)+1)/2;
    end
end

