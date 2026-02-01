function AD_50Hz_MAIN

rootdir_ad = '';

%% Preprocess data for each patient
preprocess_patient_data

% Set timestamps for stimulation start and end for each event individually
save_stim_limits
%% Micro

% ERSP & CSD & MUA
micro_plots_AD

micro_plots_interstim

% Other
[allevnr, evnr] = calc_ADpeak_nr(rootdir_ad,adpats);


% SUA
FM_stim_sua_analysis_MAIN


%% LOAD synchronized EVENTS INTO NSVIEW

%% Select from list
clearvars -except rootdir_ad;
 [ns_eeg,h_eeg,ns_thumb,h_thumb,ns_thumb_mua,h_thumb_mua,patname,eventinx,Events,evlab] = ...
     go_raw_allpatients3(rootdir_ad,'loadrec',{'EEG','thumb_low','thumb_mua'});

%% Preselected event
[ns_eeg,h_eeg,ns_thumb,h_thumb,ns_thumb_mua,h_thumb_mua,~] = load_synch_EEG_thumb(rootdir_ad,patname,eventinx,...
    'loadrec',{'EEG','thumb_low','thumb_mua'});

adnr = 1;
[ns_eeg,h_eeg,ns_thumb,h_thumb,ns_thumb_mua,h_thumb_mua,~] = load_EEG_thumb_peaks(rootdir_ad,patname,eventinx,adnr,...
    ns_eeg,h_eeg,ns_thumb,h_thumb,ns_thumb_mua,h_thumb_mua);




