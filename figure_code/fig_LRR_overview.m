function [] = fig_LRR_overview()

dbstop if error; close all;
h = tw_helpers;

winsize = 201; %local slope window size

N=7;
f = linspace(1,2,N);
ph = linspace(-1,1,N);
t = 0:1e-3:4;

for i=1:N
    fm = f(i) + 0.05*(1+sin(2*pi*f(i)*(t+ph(i))));
    x = sin(2*pi*fm.*t);
    [vals,ix] = findpeaks(x);
    X{i} = x(1:ix(3));
end

[maps,dist,mapinfo] = dtwm(X{1},X{2});
LRRs = lrr_pairwise(X,'winsize',winsize,'useparallel',true);
LRRs = LRRs';

%marginal means and sds
geom_mean = @(x)exp( nansum(log(x)) ./ sum(~isnan(x)) ); %#ok<*NANSUM> 
geom_std = @(x)exp( sqrt(   nansum((log(x)-log(geom_mean(x))).^2) ./ sum(~isnan(x))  )     );

for i=1:N
    LRR_gm{i} = geom_mean(cell2mat(LRRs(1:N,i)));
    LRR_gsd{i} = geom_std(cell2mat(LRRs(1:N,i)));
end

med_len = median(cellfun('length',LRRs),'all');
max_len = max(cellfun('length',LRRs),[],'all');
for a=1:N
    for b=1:N
        x = LRRs{a,b};
        LRRsn{a,b} = interp1(1:length(x),x,linspace(1,length(x),max_len));
    end
end
LRR_ggm{i} = geom_mean(vertcat(LRRsn{:}));
LRR_ggsd{i} = geom_std(vertcat(LRRsn{:}));

%%
ax = stf([1 1 2 2; 1 1 2 2; 3 3 2 2; 3 3 nan nan],[0.05 0.05 0.05 0.05],[0.15 0.10]);

colors = hsv(numel(X));
axes(ax(1));
for i=1:N
    plot(2*i+ X{i},'color',colors(i,:),'linew',3); hold on;
end

axes(ax(2));
axs = stfig_subaxpos(ax(2),[N+2 N+1],[0 0 0 0 0 0]);
axs = reshape(axs,N+1,[])';
delete(ax(2));

maxlen = max(cellfun('length',LRRs),[],'all');
for a=1:N
    for b=1:N
        axes(axs(a,b));
        plot(LRRs{a,b},'color',colors(a,:),'linew',3); hold on;
        plot([1 maxlen],[1 1],'k--');
    end
end
for i=1:N
    axes(axs(end-1,i));
    plot(LRR_gm{i},'color',[.5 .5 .5],'linew',3); hold on;
    plot([1 maxlen],[1 1],'k:');
    axes(axs(end,i));
    plot(LRR_gsd{i},'color','k','linew',3); hold on;   
end
axes(axs(end-1,end));
plot(LRR_ggm{i},'color',[.5 .5 .5],'linew',4); hold on;
plot([1 maxlen],[1 1],'k:');    

axes(axs(end,end));
plot(LRR_ggsd{i},'color','k','linew',4);


axis(axs,'tight');
AXS = {axs(1:end-2,1:end-1), axs(end-1,1:end), axs(end,1:end-1),axs(end,end)};
for i=1:length(AXS)
    lims = getlims(AXS{i},'xy');
    set(AXS{i},'xlim',lims(1,:),'ylim',lims(2,:),'XTick',[],'YTick',[]);
    axrescale(AXS{i},0.025,0.035);
end

shiftposy(axs(end-1:end,:),-0.015)
shiftposy(axs(end,:),-0.015)
adjheight(axs(end-1:end,:),0.01)
shiftposx(axs(:,end),0.01);
delete(axs(1:end-2,end));

%
axes(ax(3));
imagesc(mapinfo.costMatrix); hold on;
colormap(viridis(1000)); 
plot(mapinfo.index2s,mapinfo.index1s,'w-','linew',3)
set(gca,'YDir','normal');
adjwidth(ax(3),-0.1)
shiftposx(ax(3),0.125);

%
axbak = stbgax;
PP = [objpos(axs(1,1),'topleft'); objpos(axs(end-2,1),'botleft')];
PP(:,1) = PP(:,1)-0.02;
plot(PP(:,1),PP(:,2),'k','linew',2);
text(PP(1,1),mean(PP(:,2)),'LRR timeseries array', ...
    'rotation',90,'hori','center','verti','bot','fontsize',h.fs(2));

PP = objpos(axs(end-1,1),'left')-[0.01 0];
text(PP(1),PP(2),{'geometric mean','(LRR-gm)'},'hori','right','fontsize',h.fs(3));

PP = objpos(axs(end,1),'left')-[0.01 0];
text(PP(1),PP(2),{'geometric stdev','(LRR-gsd)'},'hori','right','fontsize',h.fs(3));

pp = objpos(ax(3),'top');
text(pp(1),pp(2),'DTW distance matrix','hori','center','verti','bot','fontsize',h.fs(2));
pp = objpos(ax(3),'c');
text(pp(1),pp(2),'warping curve', ...
    'hori','center','verti','bot','fontsize',h.fs(2),'color','w','rotation',40);


end






