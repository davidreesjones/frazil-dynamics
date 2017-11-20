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
tmax=3000; 

figure(6); clf; 
subplot(2,2,1); hold on; 
set(gca,'YScale','log','YLim',[1e-8,1e-2],'YTick',logspace(-8,-2,4)); box on;
ylabel('$C$','interpreter','latex')
title('$D=1$ m','interpreter','latex')
text(50,2.5e-8,'(\textit{a})','interpreter','latex')

subplot(2,2,3); hold on; 
xlabel('$t$ (seconds)','interpreter','latex')
ylabel('$T$ ($^\circ$ C)','interpreter','latex')
set(gca,'YLim',[-0.45,0],'YTick',-.4:0.1:0); box on;
text(50,-.42,'(\textit{b})','interpreter','latex')

subplot(2,2,2); hold on; 
set(gca,'YScale','log','YLim',[1e-8,1e-2],'YTick',logspace(-8,-2,4)); box on;
ylabel('$C$','interpreter','latex')
title('$D=10$ m','interpreter','latex')
text(50,2.5e-8,'(\textit{c})','interpreter','latex')

subplot(2,2,4); hold on; 
xlabel('$t$ (seconds)','interpreter','latex')
ylabel('$T$ ($^\circ$ C)','interpreter','latex')
set(gca,'YLim',[-0.45,0],'YTick',-.4:0.1:0); box on;
text(50,-.42,'(\textit{d})','interpreter','latex')

n1c=1e6; 
growth_opt=3;
[vt,vn,vT,sol,flag_crit,tsp,rsp] = mixed_layer_frazil_v1_0(tmax,ri_mm,ti_mm,T_init,n1c*n_init,eps_turb,Qext,D,tilde_nmax,growth_opt);
subplot(2,2,1); plot(vt,sum(vn.*kron(Vvi',ones(1,length(vt)))),'-','DisplayName','$f_1$');
subplot(2,2,3); plot(vt,vT,'-','DisplayName','$f_1$');

growth_opt=1;
[vt,vn,vT,sol,flag_crit,tsp,rsp] = mixed_layer_frazil_v1_0(tmax,ri_mm,ti_mm,T_init,n1c*n_init,eps_turb,Qext,D,tilde_nmax,growth_opt);
subplot(2,2,1); plot(vt,sum(vn.*kron(Vvi',ones(1,length(vt)))),'-.','DisplayName','$f_2$');
subplot(2,2,3); plot(vt,vT,'-.','DisplayName','$f_2$');

growth_opt=2;
[vt,vn,vT,sol,flag_crit,tsp,rsp] = mixed_layer_frazil_v1_0(tmax,ri_mm,ti_mm,T_init,n1c*n_init,eps_turb,Qext,D,tilde_nmax,growth_opt);
subplot(2,2,1); plot(vt,sum(vn.*kron(Vvi',ones(1,length(vt)))),'--','DisplayName','$f_3$');
subplot(2,2,3); plot(vt,vT,'--','DisplayName','$f_3$');

D=10;
n1c=1e6; 
growth_opt=3;
[vt,vn,vT,sol,flag_crit,tsp,rsp] = mixed_layer_frazil_v1_0(tmax,ri_mm,ti_mm,T_init,n1c*n_init,eps_turb,Qext,D,tilde_nmax,growth_opt);
subplot(2,2,2); plot(vt,sum(vn.*kron(Vvi',ones(1,length(vt)))),'-','DisplayName','$f_1$');
subplot(2,2,4); plot(vt,vT,'-','DisplayName','$f_1$');

growth_opt=1;
[vt,vn,vT,sol,flag_crit,tsp,rsp] = mixed_layer_frazil_v1_0(tmax,ri_mm,ti_mm,T_init,n1c*n_init,eps_turb,Qext,D,tilde_nmax,growth_opt);
subplot(2,2,2); plot(vt,sum(vn.*kron(Vvi',ones(1,length(vt)))),'-.','DisplayName','$f_2$');
subplot(2,2,4); plot(vt,vT,'-.','DisplayName','$f_2$');

growth_opt=2;
[vt,vn,vT,sol,flag_crit,tsp,rsp] = mixed_layer_frazil_v1_0(tmax,ri_mm,ti_mm,T_init,n1c*n_init,eps_turb,Qext,D,tilde_nmax,growth_opt);
subplot(2,2,2); plot(vt,sum(vn.*kron(Vvi',ones(1,length(vt)))),'--','DisplayName','$f_3$');
subplot(2,2,4); plot(vt,vT,'--','DisplayName','$f_3$');

for I=1:4
    subplot(2,2,I)
    l = legend('show');
    l.Interpreter='latex';
    l.Location='southeast';
end

set(gca,'TickLabelInterpreter', 'latex')
drawnow

