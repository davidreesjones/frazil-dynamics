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

figure(9); clf;

for I=1:8
    subplot(4,2,I); hold on;
    set(gca,'XLim',[420,520]); box on;
end

filename_list=cell(1,6);
I=0;
for opt=1:2:3
    for tilde_nmax=[0,500,4e6]
        I=I+1;
        filename_list{I}=strcat('Growth_',num2str(opt),'_tilde_nmax_',num2str(tilde_nmax,'%1.0e'));
    end
end

lin_sty=cell(1,6);
for I=1:3; lin_sty{I}='-'; end
for I=4:6; lin_sty{I}='--'; end
displayname=cell(1,6);

for I=1:numel(filename_list)
    load(filename_list{I})
    switch opt
      case 1
        growth_law='Slow (SJ04)';
      case 2
        growth_law='Fast (SO94)';
      case 3
    growth_law='Fast (RJW15)';
end
    if tilde_nmax==0
        tilde_nmax_str='0';
    else
        tilde_nmax_exp=floor(log10(tilde_nmax));
        tilde_nmax_str=strcat(num2str(tilde_nmax/10^tilde_nmax_exp,0),'$\times 10^{',num2str(tilde_nmax_exp),'}$');
    end
    displayname{I}=strcat(growth_law,', $\tilde{n}_\mathrm{max}=\,$',tilde_nmax_str);
    
    subplot(4,2,1); plot(vx/1000,vU,'LineStyle',lin_sty{I},'DisplayName',displayname{I})
    subplot(4,2,3); plot(vx/1000,vD,'LineStyle',lin_sty{I},'DisplayName',displayname{I})
    subplot(4,2,4); plot(vx/1000,vpp*(-365.25*24*3600),'LineStyle',lin_sty{I},'DisplayName',displayname{I})
    subplot(4,2,5); plot(vx/1000,vCi,'LineStyle',lin_sty{I},'DisplayName',displayname{I})
    subplot(4,2,6); plot(vx/1000,vfp*(-365.25*24*3600),'LineStyle',lin_sty{I},'DisplayName',displayname{I})
    subplot(4,2,7); plot(vx/1000,vT-vTf,'LineStyle',lin_sty{I},'DisplayName',displayname{I})
    subplot(4,2,8); plot(vx/1000,vm*(-365.25*24*3600),'LineStyle',lin_sty{I},'DisplayName',displayname{I})
end

subplot(4,2,2); ax=gca; ax.Visible='off';
for I=1:numel(filename_list)
plot([0 1],[0 0],'LineStyle',lin_sty{I},'DisplayName',displayname{I})
end
set(ax,'XTick',[],'YTick',[]); 
set(ax,'TickLabelInterpreter', 'latex')
    l = legend('show');
    l.Interpreter='latex';
    l.Location='west';
    l.Box='off';

subplot(4,2,1); ax=gca;
title('Plume speed (m/s)','interpreter','latex')
set(ax,'YLim',[0,.12],'YTick',[0:0.04:.12])
text(ax.XLim(1)+0.03*(ax.XLim(2)-ax.XLim(1)),ax.YLim(1)+0.86*(ax.YLim(2)-ax.YLim(1)),'(\textit{a})','interpreter','latex')
set(ax,'TickLabelInterpreter', 'latex')

subplot(4,2,3); ax=gca;
title('Plume depth (m)','interpreter','latex')
text(ax.XLim(1)+0.03*(ax.XLim(2)-ax.XLim(1)),ax.YLim(1)+0.86*(ax.YLim(2)-ax.YLim(1)),'(\textit{b})','interpreter','latex')
set(ax,'TickLabelInterpreter', 'latex')

subplot(4,2,5); ax=gca;
title('Frazil concentration','interpreter','latex')
text(ax.XLim(1)+0.03*(ax.XLim(2)-ax.XLim(1)),ax.YLim(1)+0.10*(ax.YLim(2)-ax.YLim(1)),'(\textit{c})','interpreter','latex')
set(ax,'TickLabelInterpreter', 'latex')
set(ax,'YScale','log','YLim',[1e-9,1e-3],'YTick',[1e-9 1e-6,1e-3])

subplot(4,2,7); ax=gca;
xlabel('Distance from grounding line (km)','interpreter','latex')
title('Supercooling ($^\circ$C)','interpreter','latex')
set(ax,'YLim',[-.05 0],'YTick',[-.05 0])
text(ax.XLim(1)+0.03*(ax.XLim(2)-ax.XLim(1)),ax.YLim(1)+0.76*(ax.YLim(2)-ax.YLim(1)),'(\textit{d})','interpreter','latex')
set(ax,'TickLabelInterpreter', 'latex')

subplot(4,2,4); ax=gca;
title('Frazil precipitation (m/yr)','interpreter','latex')
text(ax.XLim(1)+0.03*(ax.XLim(2)-ax.XLim(1)),ax.YLim(1)+0.86*(ax.YLim(2)-ax.YLim(1)),'(\textit{e})','interpreter','latex')
set(ax,'TickLabelInterpreter', 'latex')

subplot(4,2,6); ax=gca;
title('Frazil growth (m/yr)','interpreter','latex')
text(ax.XLim(1)+0.03*(ax.XLim(2)-ax.XLim(1)),ax.YLim(1)+0.86*(ax.YLim(2)-ax.YLim(1)),'(\textit{f})','interpreter','latex')
set(ax,'TickLabelInterpreter', 'latex')

subplot(4,2,8); ax=gca;
xlabel('Distance from grounding line (km)','interpreter','latex')
title('Basal freezing (m/yr)','interpreter','latex')
text(ax.XLim(1)+0.03*(ax.XLim(2)-ax.XLim(1)),ax.YLim(1)+0.86*(ax.YLim(2)-ax.YLim(1)),'(\textit{g})','interpreter','latex')

set(ax,'TickLabelInterpreter', 'latex')
drawnow

