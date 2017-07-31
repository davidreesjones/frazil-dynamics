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


figure(3); clf; 
colormap([0.8 0.8 0.8; [32,101,171]/255])
subplot(1,3,1); hold on; box on;
set(gca,'XScale','log','YScale','log')
set(gca,'XLim',[1e-6,1e0],'YLim',[1e4 1e8],'XTick',logspace(-6,0,4))
xlabel('$\epsilon$ (m$^2$ s$^{-3}$)','interpreter','latex')
ylabel('$N$ (m$^{-3}$)','interpreter','latex')

subplot(1,3,2); hold on; box on;
set(gca,'XScale','log','YScale','log')
set(gca,'XLim',[1e-1,1e2],'YLim',[1e2 1e8],'XTick',logspace(-1,2,4))
xlabel('$D$ (m)','interpreter','latex')
ylabel('$N$ (m$^{-3}$)','interpreter','latex')

subplot(1,3,3); hold on; box on; 
set(gca,'XScale','log','YScale','log')
set(gca,'XLim',[1e2,1e5],'YLim',[1e4 1e7],'XTick',logspace(2,5,4))
xlabel('Q (W m$^{-3}$)','interpreter','latex')
ylabel('N (m$^{-3}$)','interpreter','latex')

load('regime_data_E.mat')
subplot(1,3,1); 
%contourf(sE,sN,sC,[0 1],'LineStyle','none');
plot(sE(sC==1),sN(sC==1),'.','Color',[32,101,171]/255);
plot(sE(sC==0),sN(sC==0),'.','Color',[0.8 0.8 0.8]);
plot(vE,((4*vE/(15*nu)+16^2).^-0.5)*1e8/2.7/1.8,'--r','LineWidth',1.5)
 
load('regime_data_D.mat')
subplot(1,3,2);
%contourf(sD,sN,sC,[0 1],'LineStyle','none');
plot(sD(sC==1),sN(sC==1),'.','Color',[32,101,171]/255);
plot(sD(sC==0),sN(sC==0),'.','Color',[0.8 0.8 0.8]);
plot(vD,vD.^-(7/3)*1e6,'-r','LineWidth',1.5)
plot(vD,vD.^-1*0.8e7,'--r','LineWidth',1.5)

load('regime_data_Q.mat')
subplot(1,3,3);
%contourf(sQ,sN,sC,[0 1],'LineStyle','none');
plot(sQ(sC==1),sN(sC==1),'.','Color',[32,101,171]/255);
plot(sQ(sC==0),sN(sC==0),'.','Color',[0.8 0.8 0.8]);
plot(vQ,vQ.^-0*1e7/5,'--r','LineWidth',1.5)
plot(vQ,vQ.^(-2/3)*0.7e8,'-r','LineWidth',1.5)

subplot(1,3,1);
text(10^(-6+0.03*(0--6)),10^(4+0.05*(8-4)),'(\textit{a})','interpreter','latex')
text(1e-4,1e7,'Frazil','color','k','interpreter','latex')
text(5e-6,1e5,'No Frazil','color','k','interpreter','latex')
set(gca, 'Layer', 'top','TickLabelInterpreter', 'latex')

subplot(1,3,2);
text(10^(-1+0.03*(2--1)),10^(2+0.05*(8-2)),'(\textit{b})','interpreter','latex')
text(3e0,1e7,'Frazil','color','k','interpreter','latex')
text(1.5e-1,1e3   ,'No Frazil','color','k','interpreter','latex')
set(gca, 'Layer', 'top','TickLabelInterpreter', 'latex')

subplot(1,3,3);
text(10^(2+0.03*(5-2)),10^(4+0.05*(7-4)),'(\textit{c})','interpreter','latex')
text(1e3,3e6,'Frazil','color','k','interpreter','latex')
text(2e2,5e4   ,'No Frazil','color','k','interpreter','latex')
set(gca, 'Layer', 'top','TickLabelInterpreter', 'latex')

