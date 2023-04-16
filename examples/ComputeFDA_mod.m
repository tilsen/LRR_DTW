function [tAvg,nData,stiA,stiP,vA,vP,wf,lin_wf] = ComputeFDA_mod(s, varargin)
%COMPUTEFDA  - non-linear trajectory normalization
%
%	usage:  [tAvg, nData, stiA, stiP] = ComputeFDA(s, ...)
%
% Given a set of signals S registered to common onset and offset events this procedure
% uses Lucero's nonlinear timewarping technique to compute their normalized alignment,
% and returns separate STI measures for amplitude and phase variability.  The algorithm
% determines a strictly increasing and reasonably smooth transformation of time (warping
% function) such that the distance of each signal as a function of the transformed time 
% (normalized signal) from the average across normalized signals is minimized.  Smoothness
% is controlled by including a penalty for warping function roughness (LAMBDA).
%
% Expected structure of S {nSignals x 1}[nSampsN x nDims].  If SPLIT is non-empty separate
% analyses are performed:  all signals, first SPLIT, and last SPLIT
%
% The following optional 'NAME',<VALUE> pairs are supported (defaults in {braces}):
%	DIM     - dimension for comparison ([] computes tangential velocity):  {1}
%	NBREAK	- number of spline breakpoints:  {number of inflections in S(1)}
%	NSAMPS	- number of normalized samples (NSAMPS-1 must be a multiple of NBREAK):
%				default is ceil(maxTrajectoryLength/nBreak)*nBreak + 1
%	LAMBDA	- warping function roughness penalty:  {.1}
%	NORMAMP	- normalize amplitude (s - mean(s))/std(s):  'T' | {'F'} (or 1 | 0)
%	PLOT	- plot results:  {'T'} | 'F' (or 1 | 0)
%	SPLIT	- split analysis size ([] disables):  {[]}
%	NOISY	- show intermediate convergence results:  'T' | {'F'} | 'P'(arams)
%
% Returns
%	TAVG	- averaged normalized data [NSAMPS x 1]
%	NDATA	- nonlinearly normalized data [NSAMPS x nSigs]
%	STIA	- summed std. deviations of amplitude variability
%	STIP	- summed std. deviations of phase variability
%
% STI measures computed at 2% intervals.  If nSignals > SPLIT*2 each is returned as a three
% element vector specifying results from analyses for all signals, first SPLIT and last SPLIT
%
% see also COMPUTESTI
%
% cf. Lucero, J., Munhall, K., Gracco, V., & Ramsay, J. (1997), "On the registration
%     of time and the patterning of speech movements", JSLHR, 40, pp. 1111-1117.
%

% mkt 10/04 from J. Lucero's TrajAlign
% mkt 08/21 facelift

% st 04/23: edited to return warping functions and supress output

% 	defaults

dim = 1;
nBreak = [];
nSamps = [];
lambda = .01;
normAmp = 0;
noPlot = 1;
noisy = 0;
showParam = 0;
split = [];

%	parse args

if nargin < 1
	eval('help ComputeFDA');
	return;
end
if ~iscell(s)
	error('expected structure of S is a cell array {nSignals x 1}[nSampsN x 1]');
end
for ai = 1 : 2 : length(varargin)
	switch upper(varargin{ai})
		case 'DIM', dim = varargin{ai+1};
		case 'NBREAK', nBreak = varargin{ai+1};
		case 'NSAMPS', nSamps = varargin{ai+1};
		case 'LAMBDA', lambda = varargin{ai+1};
		case 'PLOT',   k = varargin{ai+1}(1); if ischar(k), noPlot = (upper(k)=='F'); else, noPlot = k; end
		case 'SPLIT',  split = varargin{ai+1};
		case 'NORMAMP',k = varargin{ai+1}(1); if ischar(k), normAmp = (upper(k)=='T'); else, normAmp = k; end
		case 'NOISY',  k = upper(varargin{ai+1}(1)); if ischar(k), noisy = (k=='T'); showParam = (k=='P'|noisy); else, noisy = k; end
		otherwise, error(sprintf('unrecognized COMPUTEFDA option (%s)', varargin{ai}));
	end
end

% 	S specified as cell array of variable length signals

k = zeros(1,length(s));
for n = 1 : length(s), if isempty(s{n}), k(n) = 1; end; end
if any(k), s(find(k)) = []; end 
nSigs = length(s);
tLen = zeros(nSigs,1);
for ti = 1 : length(s)
	if isempty(dim), s{ti} = ComputeVel(s{ti}); else, s{ti} = s{ti}(:,dim); end
	tLen(ti) = size(s{ti},1);
end

%	reformat data (into a NaN-padded array)

maxLen = max(tLen);				% longest trajectory
uData = repmat(NaN, maxLen, nSigs);
for ti = 1 : nSigs
	uData(1:tLen(ti),ti) = s{ti};
end

%	determine number of breakpoints (number of inflections in 1st signal, less slop)

if isempty(nBreak)
	nBreak = length(find(diff([0 ; diff(s{1})] > 0))) - 2;
	if nBreak > 0
		nBreak = nBreak - round(log(nBreak)*2);
	end
	if nBreak < 3, nBreak = 3; end
end

%	determine number of normalized samples

if isempty(nSamps)
	nSamps = ceil(maxLen/nBreak)*nBreak + 1;
elseif mod(nSamps-1,nBreak)
	error(sprintf('NSAMPS-1 (%d) must be a multiple of NBREAK (%d)', nSamps-1, nBreak));
end

%	linearly normalize data to common length

lData = zeros(nSamps, nSigs);
for ti = 1 : nSigs
	q = s{ti};
	if normAmp		% cvt to z-scores
		q = (q - mean(q)) ./ std(q);
	end
	lData(:,ti) = interp1(q, linspace(1,tLen(ti),nSamps)', 'linear');

    %st linear warping factors
    lin_wf(ti) = (nSamps-1)/tLen(ti);
end

%	normalize

tol = max(max(uData) - min(uData)) * .005;
if noisy, ns = 'ITER'; else, ns = 'OFF'; end
options = optimset('DISPLAY',ns,'TOLX',tol,'TOLFUN',tol,'DIFFMAXCHANGE',1,'LARGESCALE','OFF');

[nData,h,tAvg] = fdaregstr(lData, zeros(nSigs,nBreak), nBreak, lambda, options);

%st return warping function
wf = h;


% separate analyses for first split and last split signals
if ~isempty(split)
	[nData1,h1,tAvg1] = fdaregstr(lData(:,1:split), zeros(split,nBreak), nBreak, lambda, options);
	[nData2,h2,tAvg2] = fdaregstr(lData(:,end-split+1:end), zeros(split,nBreak), nBreak, lambda, options);
% middle set
	n = floor(size(lData,2)/2); r2 = floor(split/2);
	[nData3,h3,tAvg3] = fdaregstr(lData(:,n-r2:n+r2), zeros(split,nBreak), nBreak, lambda, options);
end

if showParam
	nBreak,nSamps,tol
end

%	compute variability

vA = nData - repmat(tAvg, 1, nSigs);					% amplitude
vP = h - (0:nSamps-1)' / (nSamps-1)*ones(1,nSigs);		% phase

if ~isempty(split)
	vA1 = nData1 - repmat(tAvg1, 1, split);
	vP1 = h1 - (0:nSamps-1)' / (nSamps-1)*ones(1,split);
	vA2 = nData2 - repmat(tAvg2, 1, split);
	vP2 = h2 - (0:nSamps-1)' / (nSamps-1)*ones(1,split);
end
	
%	compute STI

nBins = 50;					% 2 percent increments
idx = round(linspace(1,nSamps+1,nBins+1));
spb = max(diff(idx));
stiA = NaN * zeros(spb,nSigs,nBins);
stiP = NaN * zeros(spb,nSigs,nBins);
for k = 1 : nBins
	stiA(1:idx(k+1)-idx(k),:,k) = vA(idx(k):idx(k+1)-1,:);
	stiP(1:idx(k+1)-idx(k),:,k) = vP(idx(k):idx(k+1)-1,:);
end
if ~isempty(split)
	stiA1 = nanstd(reshape(stiA(:,1:split,:),[split*spb,nBins]));
	stiP1 = nanstd(reshape(stiP(:,1:split,:),[split*spb,nBins]));
	stiA2 = nanstd(reshape(stiA(:,end-split+1:end,:),[split*spb,nBins]));
	stiP2 = nanstd(reshape(stiP(:,end-split+1:end,:),[split*spb,nBins]));
% middle set
	n = floor(size(stiA,2)/2); r2 = floor(split/2);
	stiA3 = nanstd(reshape(stiA(:,n-r2:n+r2,:),[split*spb,nBins]));
	stiP3 = nanstd(reshape(stiP(:,n-r2:n+r2,:),[split*spb,nBins]));
end
stiA = nanstd(reshape(stiA,[nSigs*spb,nBins]));
stiP = nanstd(reshape(stiP,[nSigs*spb,nBins]));

if noPlot
	stiA = nansum(stiA);
	stiP = nansum(stiP);
	if ~isempty(split)
		stiA = [stiA , nansum(stiA1) , nansum(stiA2) , nansum(stiA3)];
		stiP = [stiP , nansum(stiP1) , nansum(stiP2) , nansum(stiP3)];
	end
	return;
end

%	plot results

fh = colordef('new', 'white');		% figure handle
figPos = get(0, 'ScreenSize');
dh = min([figPos(4) 1000]);
dw = min([figPos(3) 1000]);
figPos = [6 , figPos(4)-dh+36 , dw-22 , dh-105];
set(fh, 'position', figPos, 'name', inputname(1), 'visible','on');

cols = hsv(nSigs);
grn = [0 .6 0];

% [1,1] linearly normalized
p = [.05 .7 .44 .25];
axes('position', p, 'box', 'on');
xi = linspace(0,1,nSamps);
if ~isempty(split)
	h3 = line(xi,lData(:,split+1:end-split),'color',grn);
	h1 = line(xi,lData1,'color','b');
	h2 = line(xi,lData2,'color','r');
	legend([h1(1),h3(1),h2(1)], 'First grp', 'Mid grp', 'Last grp');
	h = [h1;h3;h2];
else
	h = plot(xi,lData);
	for n = 1 : nSigs
		set(h(n),'color',cols(n,:));
	end
end
for n = 1 : nSigs
	mh = uicontextmenu; 
	uimenu(mh, 'label', sprintf('R%d', n));
	set(h(n), 'uicontextmenu', mh);
end
	
xlabel('linear time');
if normAmp
	ylabel('z-score');
else
	ylabel('mm');					% assumption!
end
ylim = get(gca,'ylim');
v = axis;
text(v(1)+(v(2)-v(1))*.05, v(4)+(v(4)-v(3))*.05, ' linearly normalized trajectories');
%set(gca, 'buttonDownFcn','expand');

% [2,1] amplitude variability
p(2) = p(2) - .325;
axes('position', p, 'box', 'on');
if ~isempty(split)
	h1 = line(xi,vA1,'color','b');
	h2 = line(xi,vA2,'color','r');
	legend([h1(1),h2(1)], 'First grp', 'Last grp');
else
	h = plot(xi, vA);
	for n = 1 : nSigs
		mh = uicontextmenu; 
		uimenu(mh, 'label', sprintf('R%d', n));
		set(h(n),'color',cols(n,:), 'uicontextmenu',mh);
	end
end
xlabel('nonlinear time');
ylabel('\Delta from mean');
v = axis;
text(v(1)+(v(2)-v(1))*.05, v(4)+(v(4)-v(3))*.05, ' amplitude variability');
%set(gca, 'buttonDownFcn','expand');

% [3,1] amplitude STI
p(2) = p(2) - .325;
axes('position', p, 'box', 'on');
if ~isempty(split)
	a = stiA1';
	b = stiA2';
	h = plot([a , b],'-o');
	lh = line([1;1]*[1:nBins],[a,b]','color',[0 .6 0],'linestyle',':');
	k = find(b>a);
	if ~isempty(k), set(lh(k), 'color', 'm'); end
	set(h(1), 'color', 'b', 'markerFaceColor','b');
	set(h(2), 'color', 'r', 'markerFaceColor','r');
	legend(h, sprintf('First grp  (%.1f)',nansum(a)), sprintf('Last grp  (%.1f)',nansum(b)));
else
	plot(stiA,'b-o');
end
set(gca,'xlim',[0 nBins+1]);
xlabel('bins (2% intervals)');
ylabel('S.D.');
v = axis;
text(v(1)+(v(2)-v(1))*.05, v(4)+(v(4)-v(3))*.05, sprintf(' amplitude STI (Overall = %.1f)',nansum(stiA)));
%set(gca, 'buttonDownFcn','expand');

% [1,2] non-linearly normalized
p = [.55 .7 .44 .25];
axes('position', p, 'box', 'on');
if ~isempty(split)
	h3 = line(xi,tAvg,'color',grn);
	h1 = line(xi,tAvg1,'color','b');
	h2 = line(xi,tAvg2,'color','r');
	legend([h1,h3,h2], 'First grp', 'Mid grp', 'Last grp');
else
	h = plot(xi,nData); 
	for n = 1 : nSigs
		set(h(n),'color',cols(n,:));
	end
	if size(nData,2) > 3, line(xi,tAvg,'color','k','linewidth',1.5'); end
end
for n = 1 : nSigs
	mh = uicontextmenu; 
	uimenu(mh, 'label', sprintf('R%d', n));
	set(h(n), 'uicontextmenu', mh);
end
set(gca,'ylim',ylim);
xlabel('nonlinear time');
if normAmp
	ylabel('z-score');
else
	ylabel('mm');					% assumption!
end
v = axis;
text(v(1)+(v(2)-v(1))*.05, v(4)+(v(4)-v(3))*.05, ' non-linearly normalized trajectories');
%set(gca, 'buttonDownFcn','expand');

% [2,2] phase variability
p(2) = p(2) - .325;
axes('position', p, 'box', 'on');
if ~isempty(split)
%	line(xi,vP(:,6:end-5),'color',grn);
	h1 = line(xi,vP1,'color','b');
	h2 = line(xi,vP2,'color','r');
	legend([h1(1),h2(1)], 'First grp', 'Last grp');
else
	h = plot(xi, vP);
	for n = 1 : nSigs
		mh = uicontextmenu; 
		uimenu(mh, 'label', sprintf('R%d', n));
		set(h(n),'color',cols(n,:), 'uicontextmenu',mh);
	end
end
xlabel('nonlinear time');
ylabel('\Delta from mean');
v = axis;
text(v(1)+(v(2)-v(1))*.05, v(4)+(v(4)-v(3))*.05, ' phase variability');
%set(gca, 'buttonDownFcn','expand');

% [3,2] phase STI
p(2) = p(2) - .325;
axes('position', p, 'box', 'on');
if ~isempty(split)
	a = stiP1';
	b = stiP2';
	h = plot([a , b],'-o');
	lh = line([1;1]*[1:nBins],[a,b]','color',[0 .6 0],'linestyle',':');
	k = find(b>a);
	if ~isempty(k), set(lh(k), 'color', 'm'); end
	set(h(1), 'color', 'b', 'markerFaceColor','b');
	set(h(2), 'color', 'r', 'markerFaceColor','r');
	legend(h, sprintf('First grp  (%.1f)',nansum(a)), sprintf('Last grp  (%.1f)',nansum(b)));
else
	plot(stiP,'b-o');
end
set(gca,'xlim',[0 nBins+1]);
xlabel('bins (2% intervals)');
ylabel('S.D.');
v = axis;
text(v(1)+(v(2)-v(1))*.05, v(4)+(v(4)-v(3))*.05, sprintf(' phase STI  (Overall = %.1f)',nansum(stiP)));
%set(gca, 'buttonDownFcn','expand');

%	clean up

stiA = nansum(stiA);
stiP = nansum(stiP);
if ~isempty(split)
	stiA = [stiA , nansum(stiA1) , nansum(stiA2) , nansum(stiA3)];
	stiP = [stiP , nansum(stiP1) , nansum(stiP2) , nansum(stiP3)];
end

if nargout < 1, clear tAvg; end

%====================================================================================
% COMPUTEVEL  - compute central difference

function v = ComputeVel(s)

if size(s,1) == 1, s = s(:); end
ds = [diff(s([1 3],:)) ; s(3:end,:) - s(1:end-2,:) ; diff(s([end-2 end],:))] ./ 2;
if size(s,2) > 1, v = sqrt(sum(ds.^2,2)); end


%====================================================================================
% FIXNANGAPS  - interpolate through missing data

function fs = FixNaNGaps(s)

fs = NaN(size(s));
for ti = 1 : size(s,3)
	for ci = 1 : size(s,2)
		q = s(:,ci,ti);
		good = find(~isnan(q));
		if isempty(good), continue; end		% all missing -- give up
		if good(1) > 1
			q(1:good(1)-1) = q(good(1));		% first valid value
		end
		if good(end) < length(q)
			q(good(end)+1:end) = q(good(end));	% last valid value
		end
		if sum(isnan(q)) > 0					% interpolate across internal missing data points
			q(good(1):good(end)) = interp1(good,q(good),linspace(good(1),good(end),good(end)-good(1)+1)','pchip');
		end
		fs(:,ci,ti) = q;
	end
end

%====================================================================================
% FDAREGSTR  - FDA registration
%
% Jorge C. Lucero, May 22, 2000
% mkt 10/04 mods for Optimization Tbx changes

function [wdat,h,daver,c2,crit] = fdaregstr(fdat1,c,nbreak,lambda,options)

% Initialization

global fdat wdat daver h reg	% needed by sregstr

fdat = fdat1;

[points,nrec] = size(fdat);

nbreak = size(c,2);

wdat = fdat;
h = zeros(points,nrec);
daver = (sum(wdat'))'/nrec;	
crit = realmax;
error = realmax;
ipass = 0;

% Registration loop

noisy = (options.Display(1) == 'i');

while ipass < 2
	
  ipass = ipass + 1;
  
	if noisy
		fprintf('\n<< PASS %d >>\n', ipass);
	end

	for idat = 1:nrec
		if noisy
			%fprintf('\nprocessing signal #%d\n', idat);
		else
			%fprintf('.');
		end
		c(idat,:)=fminunc(@fdasregstr,c(idat,:),options,lambda,idat,points);
	    tcrit(idat)=reg;		
	end

	daver = (sum(wdat'))'/nrec;
	
	if ~noisy, fprintf('\n'); end

end

c2 = c;


%====================================================================================
% FDASREGSTR  - fminunc optimization function for FDAREGSTR
%
% Jorge C. Lucero, May 22, 2000

% Initialization

function reg = fdasregstr(c,lambda,idat,points)

global fdat wdat daver h reg

nbreak=length(c);
cc=[0 c];
delta=1/nbreak;
deltap=(points-1)/nbreak;
h(1,idat)=0;

% Roughness penalty integral
 
ww=diff(cc)/delta;
penw= sum(ww.^2)/delta;

% Warping function

for i=2:points
	
	tki=floor((i-2)/deltap);
	tk=tki*delta;
	
	x=ww(tki+1)*((i-1)/(points-1)-tk);

	if abs(x)<0.0001
		gx=1+x/2+x^2/6+x^3/2;  % This avoids losing precision when x is too small
	else
		gx=(exp(x)-1)/x;
	end

	h(i,idat)=h(tki*deltap+1,idat)+((i-1)/(points-1)-tk)*exp(cc(tki+1))*gx;
end

h(:,idat)=h(:,idat)/h(points,idat); % Normalizes the warping function

% Registered data file

wdat(:,idat)=spline(1:points,fdat(:,idat),h(:,idat)*(points-1)+1);

% Optimizing criterion

reg=norm(wdat(:,idat)-daver)^2+lambda*penw;	


