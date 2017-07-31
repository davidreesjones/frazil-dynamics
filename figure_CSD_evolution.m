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

dr=Vri*log(Vri(end)/Vri(1))/(length(Vri)-1);
figure(4); clf; hold on;
set(gca,'XScale','log','XLim',[5e-6,10e-3])
set(gca,'YScale','log','YLim',[1,1e12],'YTick',logspace(0,12,4))
h0=plot(Vri,deval(sol,0,1:Nclass)./dr','k','DisplayName','$t=0$s');
plot(Vri,deval(sol,100,1:Nclass)./dr','DisplayName','$t=100$s');
plot(Vri,deval(sol,200,1:Nclass)./dr','DisplayName','$t=200$s');
plot(Vri,deval(sol,300,1:Nclass)./dr','DisplayName','$t=300$s');
plot(Vri,deval(sol,400,1:Nclass)./dr','DisplayName','$t=400$s');
plot(Vri,deval(sol,500,1:Nclass)./dr','DisplayName','$t=500$s');
plot(Vri,deval(sol,600,1:Nclass)./dr','DisplayName','$t=600$s');
l = legend('show');
l.Interpreter='latex';
l.Location='southwest';
uistack(h0,'top')

box on
set(gca,'TickLabelInterpreter', 'latex')
xlabel('crystal size $R$ (m)','interpreter','latex')
ylabel('number density $n$ (m$^{-4}$)','interpreter','latex')
