function MUA_blocks(rootdir_ad,pat, block_ms, mua_ms,pk_thr,minprom,maxchan,baswin,postwin,varargin)

dbstop if error
% MUA plots
if ~isempty(varargin)
    eventinx = varargin{1};
else
    eventinx = [];
end



%% Plots MUA for one event of one patient (peaks & peak avgs)
if  ~iscell(pat) && ~isempty(eventinx)
    ti = varargin{2};
    chans = varargin{3};
    %     maxchan = varargin{4};
    %%
    MUA_blocks_one_event(rootdir_ad,pat,eventinx,ti,chans, block_ms, mua_ms,pk_thr,minprom,maxchan,baswin,postwin)
    
elseif  ~iscell(pat) && isempty(eventinx)
    ti = varargin{2};
    chans = varargin{3};
    MUA_blocks_one_event(rootdir_ad,{pat},eventinx,ti,chans, block_ms, mua_ms,pk_thr,minprom,maxchan,baswin,postwin)
    
else
    
    %% MUA plots -  LOOP OVER PATIENTS
    
    MUA_blocks_allpatients(rootdir_ad,pat,block_ms, mua_ms,pk_thr,minprom,maxchan,baswin,postwin);
    
    
    
    % MUA_interstim
end
end


function MUA_blocks_allpatients(rootdir_ad,patients, block_ms, mua_ms,pk_thr,minprom,maxchan,baswin,postwin)


for pp = 1:length(patients)
    patname = patients{pp};
    fprintf('%s...',patname);
    load(fullfile(rootdir_ad,'Patients',patname,'Events.mat'));
    
    if isfield(Events,['Th2_STIMdist_cm']); thnr = 2; else; thnr = 1; end;
    
    for eventinx = 1:size(Events,2)
        fprintf('%d...',eventinx);
        
        isadd = Events(eventinx).AD;
        for ti = 1:thnr
            fprintf('Th %d...',ti)
            %             if isfield(Events,['Th' num2str(ti) '_STIMdist_cm'])
            %                 if  ~isadd && isempty(Events(eventinx).(['Th' num2str(ti) '_STIMdist_cm']))
            %                     fprintf('not close\n'); continue;
            %                 else
            %                 end;
            %             else
            %                 continue;
            %             end
            
            if isempty(Events(eventinx).Stim_start);
                fprintf('No Stim start, %s %d\m',patname, eventinx ); continue;
            end;
            switch ti; case 1; chans = 1:23; case 2; chans = 25:47; end;
            MUA_blocks_one_event(rootdir_ad,patname,eventinx,ti,chans, block_ms, mua_ms,pk_thr,minprom,maxchan,baswin,postwin)
        end
    end
    
    
end
end

%--------------------------------------------------------------------------
function MUA_blocks_one_event(rootdir_ad,patname,eventinx,ti,chans, block_ms, mua_ms,pk_thr,minprom,maxchan,baswin,postwin)

load(fullfile(rootdir_ad,'Patients',patname,'Events.mat')) % load event-info containing struct

% Layer data
load(fullfile(rootdir_ad,'channels_by_layers.mat'));
if ~ismember(patname,channels_by_layers.Properties.VariableNames)
    layers = {};
    fprintf('No layer data\n')
else
    layers = channels_by_layers.(patname)(:,ti);
end


isAD = Events(eventinx).AD; % select events with AD

% Directory to save results (MUA blocks struct & figure)

matdir =  fullfile(rootdir_ad,'Figures','Micro','MUA','MUA_blocks');
% figdir = fullfile(rootdir_ad,'Figures','Micro','MUA','MUA_blocks',[patname '_Ev' num2str(eventinx) '_AD' num2str(isAD)]);
figdir = fullfile(rootdir_ad,'Figures','Micro','MUA','MUA_blocks',['AD' num2str(isAD)]);
if ~isfolder(figdir); mkdir(figdir); end;



matnm = ['MUA_stimblocks_' patname '_Ev' num2str(eventinx) '_AD' num2str(isAD) '.mat'];
if exist(fullfile(matdir,matnm))~=2
    
    % Load event
    [~,~,~,~,ns,h,~] = load_synch_EEG_thumb(rootdir_ad,patname,eventinx,...
        'loadrec',{'thumb_mua'});
    if isempty(ns)
        return;
    end
    
    % Timestamps of stim. start & end
    stimstart = Events(eventinx).Stim_start;
    stimend = Events(eventinx).Stim_end;
    stimL = stimend - stimstart + 0.2; % (0.2 puffer);
    
    guidata(ns,h)
    srate = h.srate;
    min_peak_width = floor((1/50)*srate)-70; % half-width of stim artefact peaks
    
    
    %% Find stim. artefacts + cut blocks
    [block_MUA, basblock_MUA, postblock_MUA, ...
        blocks_lims, Bblocks_lims, Pblocks_lims] = find_block_around_stimpeaks(ns,h,...
        block_ms, mua_ms,pk_thr,min_peak_width,minprom,maxchan,baswin,postwin,stimstart,stimL,matdir,matnm,patname);
    
    
    close(ns);
else
    load(fullfile(matdir,matnm));
    block_MUA = MUA_stimblocks.stimblock_MUA;
    basblock_MUA = MUA_stimblocks.basblock_MUA;
    postblock_MUA = MUA_stimblocks.postblock_MUA;
    
    blocks_lims = MUA_stimblocks.params.blocks_lims;
    Bblocks_lims = MUA_stimblocks.params.Bblocks_lims;
    Pblocks_lims = MUA_stimblocks.params.Pblocks_lims;
    
    if ~isfield( MUA_stimblocks.params,'srate')
        MUA_stimblocks.params.srate = 20000;
        srate = 20000;
    else
        srate = MUA_stimblocks.params.srate;
    end
end
%% Plot concat blocks (baseline - stim - poststim)
% figdir_ccbl = fullfile(figdir,'Concat_blocks');
% if ~isfolder(figdir_ccbl); mkdir(figdir_ccbl); end;
% 
% concat_blocks_fig(basblock_MUA, block_MUA, postblock_MUA, ...
%     blocks_lims, Bblocks_lims, Pblocks_lims, patname, eventinx,ti, chans,figdir_ccbl,srate,layers)

%% Average across blocks and concat

figdir_cca = fullfile(figdir,'Concat_blockAVGs');
if ~isfolder(figdir_cca); mkdir(figdir_cca); end;

[tm,liminx,limlabels] = concat_blockAVGs_fig(basblock_MUA, block_MUA, postblock_MUA, ...
    blocks_lims, Bblocks_lims, Pblocks_lims, patname, eventinx,ti, chans,figdir_cca,layers,rootdir_ad);

%% Smooth over block averages in time


mvw = 10;
figdir_sm = fullfile(figdir,'Concat_blockAVGs_smoothed', ['smoothed_' num2str(mvw)]);
if ~isfolder(figdir_sm); mkdir(figdir_sm); end;

smoother_layerplots(patname,eventinx,ti,tm,mvw,liminx,limlabels, figdir_sm,layers)

% mvw = 2;
% figdir_sm = fullfile(figdir,'Concat_blockAVGs_smoothed', ['smoothed_' num2str(mvw)]);
% if ~isfolder(figdir_sm); mkdir(figdir_sm); end;
% 
% smoother_layerplots(patname,eventinx,ti,tm,mvw,liminx,limlabels, figdir_sm,layers)
%% Single stim-pulse related MUA (average across blocks)
% figdir_ba = fullfile(figdir,'BlockAVG_single_pulse');
% if ~isfolder(figdir_ba); mkdir(figdir_ba); end;
% avg_across_blocks(basblock_MUA, block_MUA, postblock_MUA,patname, eventinx,ti, chans,figdir_ba,mua_ms,srate, layers)


end

%--------------------------------------------------------------------------
function avg_across_blocks(basblock_MUA, block_MUA, postblock_MUA,patname, eventinx,ti, chans,figdir,mua_ms,srate,layers)

tim = mua_ms(1):1000/srate:mua_ms(2);
tL = length(tim);
t = [1 round(tL/2) tL]+tL;
timL = [mua_ms(1) mua_ms(1)+diff(mua_ms)/2 mua_ms(2)];
% timL = repmat( [mua_ms(1) mua_ms(1)+diff(mua_ms)/2 mua_ms(2)],[1 3] );
% xti = [t t+tL t+ 2*tL];

fig = figure;
set(fig,'Visible','off')
avgbasMUA = mean(basblock_MUA(:,chans,:),3);
avgMUA = mean(block_MUA(:,chans,:),3);
avgpostMUA = mean(postblock_MUA(:,chans,:),3);
cc = cat(1, avgbasMUA, avgMUA);
cc = cat(1, cc, avgpostMUA);

pcolor(1:size(cc,1),1:length(chans),zscore(cc)'); shading interp; set(gca,'YDir','reverse'); colorbar;
colormap(jet)
xticks(t); xticklabels(arrayfun(@num2str, timL ,'uniformoutput',0) ); xlabel('Time relative to stim. peak (ms)');
ylabel('Channels')
yL = ylim;
line([tL tL],yL,'Color','k','LineStyle','-','LineWidth',3)
line([tL*2 tL*2],yL,'Color','k','LineStyle','-','LineWidth',3)
title({'single-stim. pulse related MUA (averaged across interstim-peak blocks)', [patname ', Th' num2str(ti)]});

if ~isempty(layers)
    label_layers(layers)
end

fnm = fullfile(figdir,['blockAVG_Th' num2str(ti) '_' patname '_' num2str(eventinx)]);
saveas(fig,[fnm '.jpg'])
saveas(fig,[fnm '.fig']);
close(fig);

end


%--------------------------------------------------------------------------
function smoother_layerplots(patname,eventinx,ti,tm,mvw,liminx,limlabels,figdir,layers)


% Smooth time-channel map
tm2 = smoothdata(tm,2,'movmedian',mvw);


fig = figure;
% set(fig,'Visible','off')
pcolor(tm2); shading interp;  set(gca,'YDir','reverse'); colorbar;
colormap(jet); caxis([-1.5 1.5]);
yL = ylim;
line([liminx(2) liminx(2)],yL,'Color','k','LineStyle','-')
line([liminx(4) liminx(4)],yL,'Color','k','LineStyle','--');
xlim([0 size(tm2,2)]);
xticks(liminx); xticklabels(limlabels);
title({'averaged across interstim-peak blocks',['smoothed (movmedian, window: ' num2str(mvw) ')']});
ylabel('Channels'); xlabel('Time relative to stim. start')

if ~isempty(layers)
    label_layers(layers)
end

fnm = fullfile(figdir,['TChmap_movmed_smoothed_' num2str(mvw) '_Th' num2str(ti) '_' patname '_' num2str(eventinx)]);
saveas(fig,[fnm '.jpg'])
saveas(fig,[fnm '.fig']);
close(fig);

% Average across layers
if isempty(layers)
    return;
end


fig = figure;
set(fig,'Visible','off')

figtit = {'averaged across interstim-peak blocks',['smoothed (movmedian, window: ' num2str(mvw) ')']};
plot_layerAVGs(tm2,layers,liminx,limlabels,liminx([2 4]),figtit,[],'k',false,fig);

figdir_lays = [figdir '_layers']; if ~isfolder(figdir_lays); mkdir(figdir_lays); end;
fnm = fullfile(figdir_lays,['LayerPLOT_movmed_smoothed_' num2str(mvw) '_Th' num2str(ti) '_' patname '_' num2str(eventinx)]);
saveas(fig,[fnm '.jpg'])
saveas(fig,[fnm '.fig']);
close(fig);
% 
% laymap = cat(1,laymean{:});
% fig = figure;
% set(fig,'Visible','off')
% imagesc(laymap);
% colormap(jet); caxis([-1 2]); colorbar
% ylabel('Layers');
% line([liminx(2) liminx(2)],yL,'Color','k','LineStyle','-')
% line([liminx(3) liminx(3)],yL,'Color','k','LineStyle','--');
% xlim([0 size(tm2,2)]);
% xticks(liminx); xticklabels(limlabels); xlabel('Time relative to stim. start (s) ')
% 
% title({'averaged across interstim-peak blocks',['smoothed (movmedian, window: ' num2str(mvw) ')'],...
%     ['Ev' num2str(eventinx) ', Th' num2str(ti)]});
% 
% fnm = fullfile(figdir_lays,['LayerMAP_movmed_smoothed_' num2str(mvw) '_Th' num2str(ti) '_' patname '_' num2str(eventinx)]);
% saveas(fig,[fnm '.jpg'])
% saveas(fig,[fnm '.fig']);
% close(fig);
end







function [block_MUA, basblock_MUA, postblock_MUA, ...
    blocks_lims, Bblocks_lims, Pblocks_lims] = find_block_around_stimpeaks(ns,h,...
    block_ms, mua_ms,pk_thr,min_peak_width,minprom,maxchan,baswin,postwin,stimstart,stimL,figdir,matnm,patname)
srate = h.srate;
pk_thr_orig = pk_thr;
minprom_orig = minprom;
%% Stimulation blocks

% Find peaks
ok = 0;
while ~ok
    [ns, h ,pklocs,stimdat,~] = find_stimpeaks2(ns,h,min_peak_width,pk_thr,stimstart-0.1,stimL,maxchan,minprom);
    ok = input('Peaks ok? 1/0');
    if ~ok
        pk_thr = input('New peak threshold?');
        minprom = input('New min. peak prominence?');
        maxchan = input('New channel index?');
    end
    
    close(gcf);
end

pk_ts = (stimstart-0.1) + (pklocs/h.srate);
% Cut blocks
[blockdat, peak_indeces] = cut_blocks(h,pklocs,block_ms,stimdat);
pk_ts_inblock = pk_ts(peak_indeces);

% MUA
[block_MUA, block_MUA_norm] = filtMUAblocks(blockdat,srate,block_ms,mua_ms,[],[]);
blocks_lims  = pk_ts_inblock + mua_ms/1000;
%% Baseline blocks


% Find peaks
pk_thr = pk_thr_orig;
minprom = minprom_orig;
ok = 0;
while ~ok
    %
    [ns, h ,bas_pklocs,basdat,~] = find_stimpeaks2(ns,h,min_peak_width,pk_thr,stimstart-baswin,baswin,maxchan,minprom);
    
    if strcmp(patname, 'Pt11')
        ok = input('Peaks ok? 1/0/9 if no stim artefact at all');
    else
        ok = 9;
    end
    if ok==0
        pk_thr = input('New peak threshold?');
        minprom = input('New min. peak prominence?');
        maxchan = input('New channel index?');
    elseif ok==9
        df = 370;
        bas_pklocs = 1:df:size(basdat,1);bas_pklocs = bas_pklocs';
        ok = 1;
    end
    
    
    close(gcf);
end
Bpk_ts = (stimstart-baswin) + (bas_pklocs/h.srate);


% Cut blocks
[bas_blockdat, Bpeak_indeces] = cut_blocks(h,bas_pklocs,block_ms,basdat);
Bpk_ts_inblock = Bpk_ts(Bpeak_indeces);



% MUA
[basblock_MUA, ~] = filtMUAblocks(bas_blockdat,srate,block_ms,mua_ms,[],[]);
Bblocks_lims = Bpk_ts_inblock + mua_ms/1000;


% basMUAavg = mean(basblock_MUA,[1 3]);
% basMUAsd = std(basblock_MUA,[],[1 3]);
% avgrep = repmat(basMUAavg, [size(block_MUA,1), 1]);
% sdrep = repmat(basMUAsd, [size(block_MUA,1) 1]);

%% Poststim blocks

% Find peaks
pk_thr = pk_thr_orig;
minprom = minprom_orig;
ok = 0;
while ~ok
    
    [ns, h ,post_pklocs,postdat,artend_post] = find_stimpeaks2(ns,h,min_peak_width,pk_thr,stimstart+stimL,postwin,maxchan,minprom);
    
    if strcmp(patname, 'Pt11')
        ok = input('Peaks ok? 1/0/9 if no stim artefact at all');
    else
        ok = 9;
    end;
    
    if ok==0
        pk_thr = input('New peak threshold?');
        minprom = input('New min. peak prominence?');
        maxchan = input('New channel index?');
    elseif ok==9
        df = 370;
        post_pklocs = 1:df:size(postdat,1); post_pklocs = post_pklocs';
        artend_post=length(post_pklocs);
        ok = 1;
    end
    close(gcf);
end


% Cut blocks
% post_blockdat = cut_blocks(h,post_pklocs(1:end),block_ms,postdat);
[post_blockdat, Ppeak_indeces] = cut_blocks(h,post_pklocs(1:artend_post),block_ms,postdat);

if artend_post~=length(post_pklocs)
    mean_df = round(mean(diff(post_pklocs(1:artend_post))));
    newlocs = post_pklocs(Ppeak_indeces(end)) + [mean_df:mean_df:size(postdat,1)- post_pklocs(Ppeak_indeces(end))];
    
    [post_blockdat2,Ppeak_indeces2] = cut_blocks(h,newlocs,block_ms,postdat);
    post_blockdat = cat(3,post_blockdat,post_blockdat2);
else
    newlocs = [];
    Ppeak_indeces2 = [];
end

Ppk_ts_inblock = (stimstart+stimL) + ( [post_pklocs(1:Ppeak_indeces(end)); newlocs(Ppeak_indeces2)'] /h.srate );




% MUA
[postblock_MUA, postblock_MUA_norm] = filtMUAblocks(post_blockdat,srate,block_ms,mua_ms,[],[]);

Pblocks_lims = Ppk_ts_inblock + mua_ms/1000;


%%
MUA_stimblocks.basblock_MUA = basblock_MUA;
MUA_stimblocks.stimblock_MUA = block_MUA;
MUA_stimblocks.postblock_MUA = postblock_MUA;

MUA_stimblocks.params.baswin = baswin;
MUA_stimblocks.params.stimL = stimL;
MUA_stimblocks.params.postwin = postwin;

MUA_stimblocks.params.min_peak_width = min_peak_width;
MUA_stimblocks.params.pk_thr = pk_thr;
MUA_stimblocks.params.maxchan = maxchan;

MUA_stimblocks.params.block_ms = block_ms;
MUA_stimblocks.params.mua_ms = mua_ms;

MUA_stimblocks.params.pk_ts_inblock = pk_ts_inblock;
MUA_stimblocks.params.blocks_lims = blocks_lims;

MUA_stimblocks.params.Bpk_ts_inblock = Bpk_ts_inblock;
MUA_stimblocks.params.Bblocks_lims = Bblocks_lims;

MUA_stimblocks.params.Ppk_ts_inblock = Ppk_ts_inblock;
MUA_stimblocks.params.Pblocks_lims = Pblocks_lims;


if ~isfolder(figdir); mkdir(figdir); end;
save(fullfile(figdir,matnm),'MUA_stimblocks');
end


function concat_blocks_fig(basblock_MUA, block_MUA, postblock_MUA, ...
    blocks_lims, Bblocks_lims, Pblocks_lims, patname, eventinx,ti, chans,figdir,srate,layers)
%% Concat MUA blocks (baseline - stim - poststim)
cc = [];
for j = 1:size(basblock_MUA,3)
    cc = cat(1,cc,basblock_MUA(:,chans,j));
end
for j = 1:size(block_MUA,3)
    cc = cat(1,cc,block_MUA(:,chans,j));
end
for j = 1:size(postblock_MUA,3)
    cc = cat(1,cc,postblock_MUA(:,chans,j));
end


%%
blockL = size(block_MUA,1);
[basblock_time, bascc_time] = calc_blocktime(Bblocks_lims,blockL,srate);
[block_time, cc_time] = calc_blocktime(blocks_lims,blockL,srate);
[postblock_time, postcc_time] = calc_blocktime(Pblocks_lims,blockL,srate);
alltime = [bascc_time; cc_time; postcc_time];

%%

ds = 20;
ccd = downsample(cc,ds); % downsample for plotting
zcc = zscore(ccd);

alltime_ds = downsample(alltime,ds);
alltime_ds = round(alltime_ds*100)/100;
tL = length(alltime_ds);


fig = figure;
set(fig,'Visible','off')
pcolor(1:tL, 1:length(chans), zcc'); shading interp; set(gca,'YDir','reverse'); colorbar;

hold on;
basend = floor((size(basblock_MUA,3)*size(basblock_MUA,1))/ds) ;
blend = floor(basend+ size(block_MUA,3)*size(block_MUA,1)/ds);
pmid = floor(   blend + (size(postblock_MUA,3)*size(postblock_MUA,1)/ds)/2   );

yL = ylim;
line([basend basend],yL,'Color','k','LineWidth',3)
line([blend blend],yL,'Color','k','LineWidth',3)
xti = [1 basend-40 basend+40 blend-40 blend+40 pmid size(zcc,1)];
xticks(xti);  xticklabels(arrayfun(@num2str, alltime_ds( xti) , 'UniformOutput' ,0))


colormap(jet); caxis([-3 3])
title({[patname ', Ev' num2str(eventinx)],['Th' num2str(ti)]});
ylabel('Channels'); xlabel('Time (s)');

if ~isempty(layers)
    label_layers(layers)
end

set(fig,'Position',get(0,'Screensize'))
fnm = fullfile(figdir,['Wholeblocks_Th' num2str(ti) '_' patname '_' num2str(eventinx)]);
saveas(fig,[fnm '.jpg'])
saveas(fig,[fnm '.fig'])
close(fig);
end


function [tm,liminx,limlabels] = concat_blockAVGs_fig(basblock_MUA, block_MUA, postblock_MUA, ...
    blocks_lims, Bblocks_lims, Pblocks_lims,patname, eventinx,ti, chans,figdir,layers,rootdir_ad)

allchnr = size(basblock_MUA,2);
hem = load(fullfile(rootdir_ad,'Patients',patname,['h' num2str(allchnr) '_hem_' patname '.ldr']));

basbla = squeeze(mean(basblock_MUA(:,chans,:), 1) );
stimbla = squeeze( mean(block_MUA(:,chans,:), 1) );
postbla = squeeze( mean(postblock_MUA(:,chans,:), 1) );

basbl_lim = round( diff(Bblocks_lims([1 end]))*10 )/10;
bl_lim = round( diff(blocks_lims([1 end]))*10 )/10;
blmid_lim = bl_lim/2;
postbl_lim = round( diff(Pblocks_lims([1 end]))*10)/10;
pmid_lim = postbl_lim/2;

basend =size(basblock_MUA,3);
blend = basend+ size(block_MUA,3);
blmid = basend + round(size(block_MUA,3)/2);
pmid = blend+ round(size(postblock_MUA,3)/2);


tm0 = cat(2,basbla,stimbla);
tm1 = zscore(cat(2,tm0,postbla),[],2);

tm =  (tm1'*hem(chans,chans)')' ;

    
fig = figure;
set(fig,'Visible','off')
pcolor(tm); shading interp;  set(gca,'YDir','reverse'); colorbar;
yL = ylim;
line([basend basend],yL,'Color','k','LineWidth',3)
line([blend blend],yL,'Color','k','LineWidth',3,'LineStyle','--')


liminx = [1 basend  blmid blend pmid size(tm,2)];
limlabels = arrayfun(@num2str,[0 basbl_lim blmid_lim bl_lim pmid_lim postbl_lim] ,'UniformOutput',0);
xticks(liminx); xticklabels(limlabels);
xlabel('Time relative to stim. start (s) ')
caxis([-1.5 1.5])
colormap(jet)
title({[patname ', Ev' num2str(eventinx)],['Th' num2str(ti)]});
ylabel('Channels');

if ~isempty(layers)
    label_layers(layers)
end


set(fig,'Position', get(0,'Screensize') );
fnm = fullfile(figdir,['BlockAVGs_Th' num2str(ti) '_' patname '_' num2str(eventinx)]);
saveas(fig,[fnm '.jpg'])
saveas(fig,[fnm '.fig'])
close(fig);

end


function [block_time, cc_time] = calc_blocktime(blocks_lims,blockL,srate)

blnr = size(blocks_lims,1);
block_time = nan(blnr,blockL);
for k = 1:blnr
    block_time(k,:) = blocks_lims(k,1):1/srate:blocks_lims(k,2);
end
cc_time = reshape(block_time',[blockL*blnr 1]);
end