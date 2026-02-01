function FM_stim_sua_analysis_MAIN
% 
patname = 'Pt19';

suadir = fullfile(rootdir_ad,'Patients',patname,'Thumb','SUA_analysis');

%% SPIKE DETECTION & CLUSTERING

% save_segment2spikedetection(cnt_path,cnt_filename,time_start_sec, time_end_sec,chan)
% [savedir, fnm2] = save_fnm(cnt_path,cnt_filename,time_start_sec, time_end_sec,chan);
tic
for c = 1:5
    switch c; case 1; choi = 1:10; case 2; choi = 11:20; case 3; choi = 21:30; case 4; choi = 31:40; case 5; choi = 41:47; end;
    save_BLOCKsegment2spikedetection(rootdir_ad,patname,choi);
end
toc
%

% WAVE_CLUST
tic
for chan = 1:47
    wave_clust_sortcluster(rootdir_ad,patname,chan);
end
toc

%% PREPROCESS

% Saves unit table ('SUA_data.mat') with name and quality measures of units
save_units(rootdir_ad,patname,choi,sr) 

% Add properties to unit table (saves old table as '...backup...')
add_unit_prop(rootdir_ad,prop) % prop = 'spike_nr';
ISIcountAD
isolation_distanceAD % ????

%% INTRASTIM EVENTS


%% 
event_AD_spiking(rootdir_ad,patname,eventinx,chan,clus,sr,issave)
evtype_spiking(rootdir_ad,patname,ti)
event_spiking
allunits_perevent
%%
% Auto-corr
cnt_spike_autocorr(cnt_path,cnt_filename,time_start_sec, time_end_sec,chan)


end
