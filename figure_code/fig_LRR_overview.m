function [] = fig_LRR_overview()

dbstop if error; close all;
h = tw_helpers;

winsize = 101; %local slope window size
winc = 250; %window center for zoom

S = test_signals('growth_decay_example');

D.X = S.X([1 end]);
D.rates = S.rates([1 end]);
D.rates = cellfun(@(c){c/S.dt},D.rates);

D = repmat(D,2,1);
D(2).X = D(1).X([2 1]);
D(2).rates = D(1).rates([2 1]);
%
for i=1:length(D)
    [D(i).map,D(i).distance,D(i).map_info] = dtwm(D(i).X{1},D(i).X{2});
    [D(i).Xa(1,:),D(i).Xa(2,:),D(i).map_ref] = align_signals(D(i).map,D(i).X{1},D(i).X{2},'aligntype','reference');
    [D(i).LRR] = lrr(D(i).map,'winsize',winsize);
    D(i).relrate = D(i).rates{1}./D(i).rates{2}(D(i).map_ref(2,:));
end

end
