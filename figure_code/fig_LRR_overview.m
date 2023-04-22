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
axis tight;
set(gca,'Visible','off');
adjheight(gca,-0.1);
shiftposy(gca,0.1);

axes(ax(2));
axs = stfig_subaxpos(ax(2),[N+2 N+1],[0 0 0 0 0 0]);
axs = reshape(axs,N+1,[])';
delete(ax(2));

maxlen = max(cellfun('length',LRRs),[],'all');
for a=1:N
    for b=1:N
        axes(axs(a,b));
        plot(LRRs{a,b},'color',colors(b,:),'linew',3); hold on;
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

shiftposy(axs(end-1:end,:),-0.15)
shiftposy(axs(end,:),-0.015)
adjheight(axs(end-1:end,:),0.01)
shiftposx(axs(:,end),0.01);
shiftposy(axs,-0.05);
delete(axs(1:end-2,end));


%
axes(ax(3));
imagesc(mapinfo.costMatrix); hold on;
colormap(viridis(1000)); 
plot(mapinfo.index2s,mapinfo.index1s,'w-','linew',3)
set(gca,'YDir','normal','XTickLabel',[],'YTickLabel',[]);
adjwidth(ax(3),-0.1)
shiftposx(ax(3),0.05);


%
axbak = stbgax;

PP{1} = [objpos(axs(1,1),'topleft'); objpos(axs(end-2,1),'botleft')] - [0.01 0];
PP{2} =  [objpos(ax(1),'topleft'); objpos(ax(1),'botleft')];
PP{3} =  [objpos(ax(1),'botleft'); objpos(ax(1),'botright')] - [0 0.02];
PP{4} =  [objpos(axs(end-2,1),'bl'); objpos(axs(end-2,end-1),'br')];

labs = {'array of LRR timeseries'
    'set of input signals'
    ''
    ''};

oo = [0 0 1 1];
for i=1:length(PP)
    pp = PP{i};
    bh(i) = drawbrace(pp(2,:),pp(1,:),0.0075,'color','k','linew',1);
    switch(oo(i))
        case 0
            th(i) = text(min(bh(i).XData),mean(bh(i).YData),labs{i}, ...
                'hori','center','verti','bot','fontsize',h.fs(2),'rotation',90);
        case 1
            th(i) = text(mean(bh(i).XData),min(bh(i).YData), labs{i},...
                'hori','center','verti','top','fontsize',h.fs(2),'rotation',0);
    end
end

labs = {{'geometric mean','(LRR-gm)'}
    {'geometric st.dev','(LRR-gsd)'}
    'warping curve'
    'apply DTW to each pair of inputs'
    {'calculate LRR','from warping curve'}
    'calculate marginal geometric means and st.devs'
    {'grand','mean/stdev'}};

PX(1,:) = objpos(axs(end-1,1),'left')-[0.01 0];
PX(2,:) = objpos(axs(end,1),'left')-[0.01 0];
PX(3,:) = objpos(ax(3),'c');
PX(4,:) = mean([objpos(ax(1),'bot'); objpos(ax(3),'top')]);
PX(5,:) = mean([objpos(th(1),'left'); objpos(ax(3),'right')]);
PX(6,:) = [mean(bh(end).XData) min(bh(end).YData)-0.05];
PX(7,:) = objpos(axs(end-1,end),'top');

for i=1:size(PX,1)
    thx(i) = text(PX(i,1),PX(i,2),labs{i},'fontsize',h.fs(3));
end
set(thx(1:2),'hori','right');
set(thx(3),'hori','center','verti','bot','rotation',35,'color','w');
set(thx(4:6),'hori','center','verti','mid','backgroundcolor','w');
set(thx(7),'hori','center','verti','bot','fontsize',h.fs(end));

%arrows
PA{1} = [objpos(ax(1),'bot')-[0 0.05]; objpos(ax(3),'top')+[0 0.02]];
PA{2} = [objpos(ax(3),'r')+[0.01 0]; objpos(th(1),'left')+[-0.01 0]];
PA{3} = [objpos(axs(end-2,round(end/2)),'bot')-[0 0.02]; objpos(axs(end-1,round(end/2)),'top')+[0 0.01]];

for i=1:length(PA)
    pp = PA{i};
    arrow(pp(1,:),pp(2,:),'color','k','linew',2);
end
   
LP = [objpos(ax(1),'topleft') - [0.01 -0.025];
    objpos(thx(4),'left') - [0.015 0];
    objpos(thx(5),'top') - [0 -0.015];
    objpos(axs(1,1),'topleft') - [0.01 -0.025];
    objpos(axs(end-1,1),'topleft') - [0.01 -0.025]];

labs = {'i','ii','iii','iv','v'};
for i=1:size(LP,1)
    text(LP(i,1),LP(i,2),labs{i}, ...
        'fontweight','bold','hori','center','fontsize',h.fs(2));
end

thx(4).Position(3) = 1;
thx(5).Position(3) = 1;
thx(6).Position(3) = 1;


%%
set(gcf,'InvertHardcopy','off','color','w');
h.printfig(mfilename);

end






