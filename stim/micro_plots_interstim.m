% micro_plots_interstim%% MUA blocks


pats = {'Pt11','Pt3','Pt8','Pt20'}; % if loop over patients | {'Pt11','Pt3','Pt8','Pt20','Pt5'}
% % evxx = {[86, 86, 68,68, 90, 90, 81,81]; [153,153, 172];54;39};
% evxx = {[86, 86, 68,68, 90, 90, 81,81]; [153,153, 172];54;39};
% tix = {[1, 2,1, 2, 1, 2, 1, 2]; [1,2,2]; 1;2};
vargin = {};

pats = 'Pt3'; % if only one patient
eventinx = 116; % if only 1 event of one patient | leave empty for all events of one patient
ti = 2; % only if 1 event of one patient
if ti==1; chans = 1:23; else; chans = 25:47; end; % only if 1 event of one patient (vector of thumb channel indeces)
maxchan = [];
vargin = {eventinx, ti, chans, maxchan};

intens_correct = true;

%% MUA prepost_stim
%mua_prepoststim_ratio

%% Cut, plot and save interstim MUA blocks

%%%
block_ms = [5 16]; % prv: [5 15]; % window around peak to cut
mua_ms = [6 15]; % [7 13]; % window around peak to plots
baswin = 1;
postwin = 2;
%%%

pk_thr = 5; % 2; SD
minprom = 1; % min. peak prominence
maxchan = 1; % channel to detect peaks

MUA_blocks(rootdir_ad,pats, block_ms, mua_ms,pk_thr,minprom,maxchan,baswin,postwin,vargin{:});


%% Event types

patsss = {'Pt3','Pt3', 'Pt8', 'Pt11', 'Pt20','Pt15','Pt19'};
thbsss = [1 2 1 2 2 1 1];

%% Cluster threshold for permutation statistics with cluster correction?
generate_random_MUA(rootdir_ad,patsss,block_ms,mua_ms)

%%

period = 'stim'; % 'stim' | 'post'
chlab = {}; % {} | 'Supra' | 'Infra'

for p = 1:length(patsss)
    patname = patsss{p};
    disp(patname)
    ti = thbsss(p);
    %%
    [evinxu, evnru] = calc_evtype_nrs(rootdir_ad,patname,ti, 'mua_up',period,chlab);
    [evinxd, evnrd] = calc_evtype_nrs(rootdir_ad,patname,ti, 'mua_down',period,chlab);
    
%     [evinxd, evnrd] = calc_evtype_nrs(rootdir_ad,patname,ti, 'all');
end

%% Event summary tables
% extract_eventnrs_table(rootdir_ad,patsss,thbsss);

extract_eventnrs_table(rootdir_ad,patsss,thbsss,'muachange',true,'period',period,'chlab',chlab)

%% Chi2 stat
compare_evtypenrs_chi2stat(rootdir_ad,patsss,thbsss)

%% Stim. vs. Post-stim MUA increase/ decrease events

% Overlap - map and pie chart
stim_postMUA_evtype_overlap(rootdir_ad,patsss,thbsss)

% Correlation between stim. and post-stim MUA change (rel. to baseline)
stim_postMUA_corr(rootdir_ad,patsss,thbsss)
%% AD vs noAD interstim MUA within one patient

evtype =  {'Within','Close','MUA decrease','noAD, out of SOZ'}; % 'close' | 'mua_up' | 'mua_down'| 'all' | 
% {'noAD, out of SOZ'/ 'noAD, in SOZ'/ 'AD',
%'Within'/ 'Close'/ 'Far',
% 'MUA increase'/ 'MUA decrease' / 'no MUA change'};


%%
for k = 1:length(patsss)
   patname = patsss{k};
    disp(patname)
    ti = thbsss(k);
    %%
    MUA_blocks_AVG_onepatient(rootdir_ad,patname,ti,'evtype',evtype,'maxintens',false,'bnorm','Indiv','itnum',1000)
end

%%
MUA_allpats(rootdir_ad,patsss,thbsss,'evtype',evtype,'bnorm','Indiv','patavg',false,'maxintens',true)

%% Count channels with local min/ max
for k = 1:length(patsss)
    patname = patsss{k};
    disp(patname)
    ti = thbsss(k);
    MUA_ch_latency(rootdir_ad,patname,ti,'evtype',evtype,'maxintens',false,'bnorm','Indiv')
end
%%
compar_tags =   {'noADzone'}; % {'ADzone AD', 'ADzone noAD'} | {'ADzone noAD', 'noADzone'} | {'ADzone AD', 'noADzone'}

% c_pats = {'Pt3', 'Pt3', 'Pt8', 'Pt11', 'Pt20';};
% c_thumbs = [1 2 1  2 2];

for k = 1:length(patsss)
    patname = patsss{k};
    ti = thbsss(k);
    %%
    compare_MUA_blocks(rootdir_ad,patname,ti,compar_tags,'evtype',evtype,'intens_correct',true,'maxintens',true,...
        'bnorm','Common','stattype','my_permutation','mcorr','none');
end

%%
compare_MUA_allpats(rootdir_ad,patsss,thbsss,compar_tags,'evtype',evtype,'intens_correct',false,'maxintens',false,...
    'bnorm','Common','stattype','my_permutation','mcorr','none','patavg',false);


%% MUA stim param correlations


for k = 1:length(patsss)
    patname = patsss{k};
    ti = thbsss(k);
    %     MUA_blocks_stimintens_correlation(rootdir_ad,patname,'intensity',intens_correct);
    %     MUA_blocks_stimintens_correlation(rootdir_ad,patname,'distance',intens_correct);
    
    %%
    %%
    MUA_stimint_corr_onepatient_ONESITE(rootdir_ad,patname,ti,'intensity','maxintens',false,'corrmap',false,'corrLine',true,'evtype',evtype);
end


%%

for k = 4:length(patsss)
    patname = patsss{k};
    ti = thbsss(k);
    %%
    MUA_stimint_corr_evtypes(rootdir_ad,patname,ti,'maxintens',true, 'bnorm','Indiv','evtype',evtype,'maxintens',false)

end

%% Line correlation  - pool all patient
MUA_stimint_corr_evtypes_allpatients(rootdir_ad,patsss,thbsss,'evtype',evtype,'bnorm','Indiv','patavg',false,'maxintens',false)




%% CORRELATION MAP COMPARISON
evtype_c = {};

corrmap_compare(rootdir_ad,patsss,thbsss,evtype_c)

%% Correlation on MUA avg - comparison
for k = 1:length(patsss)
    patname = patsss{k};
    ti = thbsss(k);
    
    %%
    MUA_stimint_corr_compare_onepatient(rootdir_ad,patname,ti,compar_tags,...
        'evtype',evtype,'intens_correct',true,'maxintens',false,'bnorm','Indiv','corrperiod','stim')
    
end


%%
MUA_stimint_corr_compare_allpatients(rootdir_ad,patsss,thbsss,compar_tags,...
        'evtype',evtype,'intens_correct',true,'maxintens',true,'bnorm','Indiv',...
        'stattype','my_permutation','mcorr','none')

%%
pats = {'Pt3','Pt5','Pt8','Pt11','Pt12','Pt15','Pt19','Pt20'};
best_thumb = [2, 1, 1, 2, 1, 1, 1, 1, 2];


%% MUA map
draw_MUAmap

%% MUA stim - th distance correlation

for p = 1:length(patsss)
    patname= patsss{p};
    ti= thbsss(p);
    mua_distance_corr(rootdir_ad,patname,ti,'evtype','mua_up')
    mua_distance_corr(rootdir_ad,patname,ti,'evtype','mua_down')
    
    mua_distance_corr(rootdir_ad,patname,ti,'evtype','mua_both')
    
end



%% Compare Th-stim distance of MUA increase/ decrease events

mua_updown_dist_compare(rootdir_ad, pats, thumbs)