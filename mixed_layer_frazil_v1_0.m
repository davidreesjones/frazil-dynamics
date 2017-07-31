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

function [vt,vn,vT,sol,flag_crit,tsp,rsp] = mixed_layer_frazil_v1_0(tmax,ri_mm,ti_mm,T_init,n_init,eps_turb,Qext,D,tilde_nmax,growth_opt)
%mixed_layer_frazil_v1_0 Mixed layer frazil calculation
%   Version used in Cryosphere calculations
%   Includes three different parameterizations of frazil growth rate (growth_opt: note
%   different numbering relative to paper, see comments below)

% material properties
L=3.35e5; %J/kg
cw=3974; %J/kg/C
nu=1.95e-6; %m^2/s
rho_w=1030; %kg/m^3 1030 in JB95, 1028 in SJ04
rho_i=920; %kg/m^3 920 in JB95, 917 in SJ04
Nu=1;
KT=1.4e-7; %m^2/s thermal diffusivity
Tf=0; %degC

% CRYSTAL SET UP
% properties
Nclass=length(ri_mm);
Vri=ri_mm*1e-3; % convert mm to m
Vti=ti_mm*1e-3; % convert mm to m
if numel(Vti)==1
    Vti=Vti(1)+zeros(1,Nclass);
end
Vvi=pi*Vri.^2.*Vti; % volume
Vwi=16*Vri; % rise velocity

% crystal precipitation
g=Vwi/D;

% secondary nucleation
Ur=(4*eps_turb*Vri.^2/(15*nu)+Vwi.^2).^(0.5);
a=pi.*(Vvi(1)./Vvi).*Ur.*Vri.^2; %alpha

%floculation
b=zeros(size(a)); % set equal to zero

%growth
switch growth_opt
    case 1 % Corresponds to law f_2 in manuscript
        G=Nu*(rho_w*cw*KT/(rho_i*L))*(2*pi*Vri(1:Nclass-1)./diff(Vvi));
        Qt=Nu*(rho_w*cw*KT)*(2*pi*Vri);
    case 2 % Corresponds to law f_3 in manuscript
        G=Nu*(rho_w*cw*KT/(rho_i*L))*(2*pi*Vti(1:Nclass-1)./diff(Vvi));
        Qt=Nu*(rho_w*cw*KT)*(2*pi*Vti);
    case 3 % Corresponds to law f_1 in manuscript
        G=Nu*(rho_w*cw*KT/(rho_i*L))*(2*pi*Vri(1:Nclass-1)./diff(Vvi))./(0.9008-0.2634*log(Vti(1:Nclass-1)./(2*Vri(1:Nclass-1))));
        Qt=Nu*(rho_w*cw*KT)*(2*pi*Vri(1:Nclass-1))./(0.9008-0.2634*log(Vti(1:Nclass-1)./(2*Vri(1:Nclass-1))));
    case 4
        G=Nu*(rho_w*cw*KT/(rho_i*L))*(2*pi*0.04*Vri(1:Nclass-1)./diff(Vvi));
        Qt=Nu*(rho_w*cw*KT)*(2*pi*Vti);
    
    otherwise
        disp('growth_opt entered incorrectly')
        return
end

% Reorientate
b=b';
g=g';
G=G';
a=a';
Qt=Qt';

% MAIN Ode Routine
% Initial conditions are a user input, apart from seeding
%F0=5e-9;
y0(1:Nclass)=n_init;
y0(Nclass+1)=T_init;
opts = odeset('RelTol',1e-10,'AbsTol',1e-13,'NonNegative',1:Nclass); %ensures ni>0
sol=ode15s(@ode,[0 tmax],y0,opts); % Integrate to nucleation site
vt=sol.x;
vn=sol.y(1:Nclass,:);
vT=sol.y(Nclass+1,:);

rbar=sum(vn.*kron(Vri',ones(1,length(vt))))./sum(vn);


if abs(vT(end)/(-Qext*tmax/(rho_w*cw))-1)<0.5
    flag_crit=0;
    tsp=nan;
    rsp=nan;
else
    flag_crit=1;
    [~,tspI]=min(vT);
    tsp=vt(tspI);
    rsp=rbar(tspI);
end

    function    [dydt, n, T]=ode(x,y)
        % Main ode function
        LL=length(x);
        n=y(1:Nclass,:);
        T=y(Nclass+1,:);
        
        % CALCULATE CRYSTAL GROWTH AND PRECIPITATION
        % secondary nucleation
        tilde_n=min(sum(y(1:Nclass,:)),tilde_nmax*ones(1,LL));
        an=zeros(Nclass,LL);
        for k=2:Nclass
            an(k,:)=-a(k)*y(k,:).*tilde_n;
        end
        
        an(1,:)=-sum(an.*kron(Vvi'./Vvi(1),ones(1,LL)));
        dydt=an;
        
        % freezing/melting
        DT=Tf-T;
        for k=1:Nclass-1
            dydt(k,:)=dydt(k,:)-(b(k)+G(k)*DT).*y(k,:);
        end
        for k=2:Nclass
            dydt(k,:)=dydt(k,:)+(Vvi(k-1)*b(k-1)/Vvi(k)+G(k-1).*DT).*y(k-1,:);
        end
        
        % gravitational removal
        for k=1:Nclass
            dydt(k,:)=dydt(k,:)-g(k)*y(k,:);
        end
        
        % heat balance
        Q=zeros(Nclass,LL);
        for k=1:Nclass-1
            Q(k,:)=Qt(k)*DT.*y(k,:);
        end
        dydt(Nclass+1,:)=(-Qext+sum(Q))/(rho_w*cw);
        
    end



end

