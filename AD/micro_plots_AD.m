function micro_plots_AD


pat = {'Pt3', 'Pt8', 'Pt11', 'Pt20'}; % if loop over patients | {'Pt11','Pt3','Pt8','Pt20'}
% % evxx = {[86, 86, 68,68, 90, 90, 81,81]; [153,153, 172];54;39};
% evxx = {[86, 86, 68,68, 90, 90, 81,81]; [153,153, 172];54;39};
% tix = {[1, 2,1, 2, 1, 2, 1, 2]; [1,2,2]; 1;2};
vargin = {};

pat = 'Pt20'; % if only one patient
eventinx =16; % if only 1 event of one patient | leave empty for all events of one patient
ti = 2; % only if 1 event of one patient
if ti==1; chans = 1:23; else; chans = 25:47; end; % only if 1 event of one patient (vector of thumb channel indeces)
maxchan = [];
vargin = {eventinx, ti, chans, maxchan};

win = [-.2 .5]; %[-.2 .5] |  [-.1 .1]
freqs = [0.5 80];
isfig = true;
indivpeaks = false;
pktype = 'Peaktimes'; % Small_peaktimes | Peaktimes

lf_freq = 30;
%% One event/ patient
% RAW
raw_plots(rootdir_ad,pat,win,freqs,isfig,indivpeaks,vargin{:})

% CSD
% Prepare matrices for hemming filter and csd (.ldr files) for the patient manually!
 csd_plots(rootdir_ad,pat,win,freqs,isfig,indivpeaks,pktype,vargin{:})

% MUA
muafig = true; timingfig = false;
mua_plots(rootdir_ad,pat,win,muafig,indivpeaks,pktype,timingfig,vargin{:})


% PAT averages
%%
raw_patavg(rootdir_ad,patname, eventss,ti, chans, win, freqs, pktype)
%% CSD
 
% eventss = % Pt3: [133 150 153 157 172], Pt11: [68 71 78 81 86];
csd_patavg(rootdir_ad,patname, eventss,ti, chans, win, freqs, pktype)
%%
mua_patavg(rootdir_ad,patname, eventss,ti, chans, win, pktype,lf_freq)

%%
raw_csd_mua_layers(rootdir_ad,patname, eventss,ti, chans, win, pktype,lf_freq)
%%
%% ERSP, Channet-Time map
smwins = [-.2 -.05; -.05 .1; .1 .25]; % [-.2 -.05; -.05 .1; .1 .25; .25 .45] |  [-.2 -.05; -.05 .1; .1 .25];
freqs2plot = [150 300] ; % [1 4; 4 7; 8 13; 13 30; 30 50; 80 150] | [13 30];;
savewin = [-.3 .3];

savewav = false;
bychannel = false;
timemap = true;
freqmap = false;
freqtimng = false;
chantimeAVG = false;


for k = 1:length(pats)
    patname = pats{k};
    fprintf('%s...',patname);
    load(fullfile(rootdir_ad,'Patients',patname,'Events.mat'));
    
    if isfield(Events,['Th2_STIMdist_cm']); thnr = 2; else; thnr = 1; end;
    ADevxx = find([Events.AD]);
    %%
    for eventinx = ADevxx
        fprintf('%d...',eventinx);
        
        for ti = 1:thnr
            
            switch ti; case 1; chans = 1:23; case 2; chans = 25:47; end;
            
            maxchan = [];
            vargin = {eventinx, ti, chans, maxchan};
            
            
            ersp_plots(rootdir_ad,patname,win,smwins,savewin,freqs,freqs2plot,indivpeaks,...
                savewav,bychannel, timemap, freqmap, freqtimng, chantimeAVG,vargin{:})
        end
    end
end

%% Compare latencies

MUA_LFP_timing

%%
pca_components
%% Post stim windows

poststim_raw
poststim_wav
poststim_CSD

%%
ADtype = 'SW';

for ip = 1:length(patients)
    patname = patients{ip};
    
    patdir = fullfile(rootdir_ad,'Patients',patname);
    cd(patdir)
    AD_pk_freqs(rootdir_ad,patname,ADtype)
end

%% ADmap
draw_AD_gridmaps
end



%--------------------------------------------------------------------------
function csd_patavg(rootdir_ad,patname, eventss,ti, chans, win, freqs, pktype)

figdir = fullfile(rootdir_ad,'Figures','Micro','CSD');

adnr=1;
patdir = fullfile(rootdir_ad,'Patients',patname);
cd(patdir)

load(fullfile(rootdir_ad,'channels_by_layers.mat'))
layers = channels_by_layers.(patname)(:,ti);

ldrnm_csd = dir(fullfile(patdir,'*csd.ldr'));
m_csd = load(fullfile(patdir,ldrnm_csd.name));
ldrnm_hem = dir(fullfile(patdir,['*hem_' patname '.ldr']));
m_hem = load(fullfile(patdir,ldrnm_hem.name));

evnr = length(eventss);
csd_allpats = cell(1,evnr);
for k = 1:evnr
    eventinx = eventss(k);
    [normepoch, sr] = load_filt_epoch(rootdir_ad,patname,eventinx,adnr,win,freqs,pktype);
    
    % CSD
    
    time = win(1):1/sr:win(2); time = time(1:size(normepoch,1));
    
    
    epoch2csd = mean(normepoch,3);
    
    
    
    thdat_hem = epoch2csd*m_hem';
    csdm = thdat_hem*m_csd';
    csd_pat = csdm(:,chans);
    csd_allpats{k} = csd_pat;
end

csd_all = cat(3,csd_allpats{:});
csd_avg = mean( csd_all ,3);

basinx = 1:abs(win(1))*sr;

fig = csd_simple(csd_avg,time,[-0.5 0.5],2);

csd_all_stat = permute(csd_all,[2 1 3]);
[exactp_ersp,maskersp,alphafdr] = boostat_eeglab_J(csd_all_stat,1:length(chans),0.05,200,false,'fdr',[],basinx);

hold on; contour(time,1:length(chans),maskersp,'Color','white','LineWidth',.7)
elab = arrayfun(@(x) [num2str(x) ' '],eventss,'UniformOutput',false);
title({[patname ', n = ' num2str(evnr)],['Evs: ' elab{:}]});
xlabel('Time (s) ')
label_layers(layers);

fnm = fullfile(figdir,[patname '_AVG_selected_events_TH' num2str(ti)]);
saveas(fig,[fnm '.jpg'])
saveas(fig,[fnm '.fig'])
close(fig);
end


%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
function raw_patavg(rootdir_ad,patname, eventss,ti, chans, win, freqs, pktype)

figdir = fullfile(rootdir_ad,'Figures','Micro','rawAD');

adnr=1;
patdir = fullfile(rootdir_ad,'Patients',patname);
cd(patdir)

load(fullfile(rootdir_ad,'channels_by_layers.mat'))
layers = channels_by_layers.(patname)(:,ti);


evnr = length(eventss);
raw_allpats = cell(1,evnr);
for k = 1:evnr
    eventinx = eventss(k);
    [normepoch, sr] = load_filt_epoch(rootdir_ad,patname,eventinx,adnr,win,freqs,pktype);
    
    
    time = win(1):1/sr:win(2); time = time(1:size(normepoch,1));
    
        
    raw_allpats{k} = mean(normepoch(:,chans,:),3);
end

raw_all = cat(3,raw_allpats{:});
raw_avg = mean( raw_all ,3);

basinx = 1:abs(win(1))*sr;
ch = arrayfun(@(x) x, 1:length(chans),'UniformOutput',0);

fig = figure;
tinx = [1 basinx(end) length(time)];
tlab = arrayfun(@(x) num2str(round(x*10)/10) ,time(tinx), 'UniformOutput',0);
plot_layerAVGs(raw_avg',ch,tinx,tlab,[],'',[],'k',false,fig);


fig = figure;
pcolor(time, 1:length(chans), raw_avg'); shading interp; set(gca,'YDir','reverse'); colorbar; 
caxis([-1 1])
raw_all_stat = permute(raw_all,[2 1 3]);
[exactp_ersp,maskersp,alphafdr] = boostat_eeglab_J(raw_all_stat,1:length(chans),0.05,200,false,'fdr',[],basinx);

hold on; contour(time,1:length(chans),maskersp,'Color','red','LineWidth',.7)
elab = arrayfun(@(x) [num2str(x) ' '],eventss,'UniformOutput',false);
title({[patname ', n = ' num2str(evnr)],['Evs: ' elab{:}]});
xlabel('Time (s) ')
label_layers(layers);
colormap(bone)
fnm = fullfile(figdir,[patname '_AVG_selected_events_TH' num2str(ti)]);
saveas(fig,[fnm '.jpg'])
saveas(fig,[fnm '.fig'])
close(fig);
end


%--------------------------------------------------------------------------
function mua_patavg(rootdir_ad,patname, eventss,ti, chans, win, pktype,lf_freq)

figdir = fullfile(rootdir_ad,'Figures','Micro','MUA','MUA_AD');

adnr=1;
patdir = fullfile(rootdir_ad,'Patients',patname);
cd(patdir)

load(fullfile(rootdir_ad,'channels_by_layers.mat'))
layers = channels_by_layers.(patname)(:,ti);

evnr = length(eventss);
dat_allpats = cell(1,evnr);
for k = 1:evnr
    eventinx = eventss(k);
    [normepoch, sr] = load_filt_MUAepoch(rootdir_ad,patname,eventinx,adnr,win,pktype,lf_freq);
    
    % CSD
    
    time = win(1):1/sr:win(2); time = time(1:size(normepoch,1));
    
    
    epochavg = mean(normepoch,3);
    
    dat_pat = epochavg(:,chans);
    dat_allpats{k} = dat_pat;
end

dat_all = cat(3,dat_allpats{:});
dat_avg = mean( dat_all ,3);

basinx = 1:abs(win(1))*sr;


fig = figure;
pcolor(time,1:length(chans),dat_avg'); shading interp; colorbar; colormap(parula)
set(gca,'ydir','reverse')
xlabel('Time (s)'); ylabel('Channels');
caxis([-0.5 0.5]);


dat_all_stat = permute(dat_all,[2 1 3]);
% [exactp_ersp,maskersp,alphafdr] = boostat_eeglab_J(dat_all_stat,1:length(chans),0.05,200,false,'fdr',[], basinx);
[exactp_ersp,maskersp,alphafdr] = boostat_eeglab_J(dat_all_stat,1:length(chans),0.05,200,false,'fdr');

hold on; contour(time,1:length(chans),maskersp,'Color','white','LineWidth',.7)

elab = arrayfun(@(x) [num2str(x) ' '],eventss,'UniformOutput',false);
title({[patname ', n = ' num2str(evnr)],['Evs: ' elab{:}]});
xlabel('Time (s) ')
label_layers(layers);
fnm = fullfile(figdir,[patname '_AVG_selected_events_TH' num2str(ti)]);
saveas(fig,[fnm '.jpg'])
saveas(fig,[fnm '.fig'])
close(fig);
end



