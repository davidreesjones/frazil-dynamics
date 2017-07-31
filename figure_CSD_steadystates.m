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

N=128; 
ti_mm=0.05; % Thickness in MM
ri0=0.2/1000; % Initial average radius
ri_mm=logspace(log10(5e-3),log10(20),N); Vri=ri_mm/1000;
n_init=Vri*log(Vri(end)/Vri(1))/(N-1)/(2*ri0);
n_init(Vri>2*ri0)=0; %NB logical indexing used instead of heaviside function
Qext=1200; 
eps_turb=5e-3;
nu=1.95e-6; 
T_init=0;
D=1;
tilde_nmax =4e6;
Vti=ti_mm*1e-3; Vvi=pi*Vri.^2.*Vti; Nclass=N;
tmax=15000; 


dr=Vri*log(Vri(end)/Vri(1))/(length(Vri)-1);
figure(7); clf; hold on;

n1c=4e6; 

[vt,vn,vT,sol,flag_crit,tsp,rsp] = mixed_layer_frazil_v1_0(tmax,ri_mm,ti_mm,T_init,n1c*n_init,eps_turb,Qext,D,tilde_nmax,3);
plot(Vri,deval(sol,15000,1:Nclass)./dr','DisplayName','$f_1$');

[vt,vn,vT,sol,flag_crit,tsp,rsp] = mixed_layer_frazil_v1_0(tmax,ri_mm,ti_mm,T_init,n1c*n_init,eps_turb,Qext,D,tilde_nmax,1);
plot(Vri,deval(sol,15000,1:Nclass)./dr','DisplayName','$f_2$');

[vt,vn,vT,sol,flag_crit,tsp,rsp] = mixed_layer_frazil_v1_0(tmax,ri_mm,ti_mm,T_init,n1c*n_init,eps_turb,Qext,D,tilde_nmax,2);
plot(Vri,deval(sol,15000,1:Nclass)./dr','DisplayName','$f_3$');

set(gca,'XScale','log','XLim',[5e-6,1e-3])
set(gca,'YScale','log','YLim',[1e6,1e14],'YTick',logspace(6,14,5))
l = legend('show');
l.Interpreter='latex';
l.Location='southwest';
box on
set(gca,'TickLabelInterpreter', 'latex')
xlabel('crystal size $R$ (m)','interpreter','latex')
ylabel('number density $n$ (m$^{-4}$)','interpreter','latex')

drawnow