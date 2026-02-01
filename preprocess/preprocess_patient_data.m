% PREPROCESS PATIENT DATA


% % look through data & events in nswiew before proceeding
% % complete data in patient_info excel sheet

% pattable = readtable(fullfile(rootdir_ad,'patient_info.xlsx'),'Sheet','Sheet1');
% patlist = pattable.Patient;
%   pinx = find(pattable.AD==1);
%
% pats = pattable.Patient(pinx);
% pats = patlist(~isnan(pattable.PatientCode));
% pinx = find(strcmp(patname,pats));

% adsaved = pattable.ADsSavedToEvents_mat(pinx);
% labeltype = pattable.EventLabels{pinx};
%% Patient information
% for ip = 1:length(pats)
%     patname = pats{ip};

patname = 'Pt11';
patdir = fullfile(rootdir_ad,'Patients',patname);
cd(patdir)


%% In case of older data where stimulation start is marked manually
% Align manually placed marking to first peak of stimulation artefact

% if strcmp(labeltype,'old')
    correct_old_eventlabels
% end


%% Make Events.mat file + find ADs
% Save Events.mat file containing information on stimulation events
% Correct faulty event labels:
% - in case of newer data where stimulation start is marked automatically,
% automatic event labels might be faulty -> correct if necessary
% (Events.Orig_Label contain original event labels, Events.Label contain corrected labels to match channel names)
% - in case of older data, labels should be correct, in this case only copy
% Orig_Label column to Label column

ev_eeg(patdir)

%%
% find_EEG_stimend(rootdir_ad,patname) % add stimulation end timestamps for EEG
%% Add AD field to Events.mat (if epoch contained AD)
% Based on AD information in patient_info excel sheet
% Each patient should have a separate sheet with AD containing epochs:
%       B column: stimulation channel pair; C column: stimulation
%       intensity; D column: channels presenting ADs

%     if adsaved==0
saveADs2Events(rootdir_ad,tabledir,patname)
%     end


%% Synchronize ECoG & Thumb

synchronize_eeg_thumb

%% Mark individual AD peaks on grid recs manually -> correct to actual AD peak -> synchronize to thumb
ADpeaks_mark_save


%% Read distance of Thumb from AD-channels


read_thumb_ADdist(rootdir_ad,patname); % mark distances from AD-presenting 
read_thumb_STIMdist(rootdir_ad,patname) % mark distances only from stimulating channels
% end


%% Save stim limits, based on thumb rec (only for close events)
save_stim_limits(rootdir_ad,patname)

%% Clean stimulation artifacts + create epoched eeglab dataset
% % EEG
% win_sec = [-7 17];
% stim_sec = [0 4];
% bas_sec = [-5 -1];
% isfig = 0;
%
%
% % EEG
%
% EEG2set(patdir,win_sec,stim_sec,bas_sec,isfig)



% newEEGsets(patdir,1,0)
