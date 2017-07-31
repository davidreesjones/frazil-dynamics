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

h=logspace(-3,0);
ar=0.02;
f2=ones(size(h));
f1=(0.9008-0.2634*log(h)).^(-1);
f3=2*h;
f4=2*ar*ones(size(h));


figure(1); clf; 
subplot(2,1,1); hold on; 
set(gca,'YScale','log','YLim',[1e-8,1e-2],'YTick',logspace(-8,-2,4)); box on;
ylabel('$C$','interpreter','latex')
text(50,4e-8,'(\textit{a})','interpreter','latex')

figure(1); clf; hold on;
plot(h,f1,'-','DisplayName','$f_1$');
plot(h,f2,'-.','DisplayName','$f_2$');
plot(h,f3,'--','DisplayName','$f_3$');
set(gca,'ColorOrderIndex',3)
plot(0.02,2*0.02,'s','DisplayName','$f_3(h=0.02)$');

set(gca,'XScale','log','YScale','log','XDir','reverse','YLim',[1e-2,2])
l = legend('show');
l.Interpreter='latex';
l.Location='southwest';
box on
set(gca,'TickLabelInterpreter', 'latex')
xlabel('Aspect ratio $h=H/2R$','interpreter','latex')
ylabel('Growth rate $f$','interpreter','latex')

