function [evinx, evnr] = calc_evtype_nrs(rootdir_ad,patname,ti, evtype, period,chlab)

patdir = fullfile(rootdir_ad,'Patients',patname);
load(fullfile(patdir,'Events.mat'));
load(fullfile(rootdir_ad, 'good_channels.mat'))


if ~isempty(chlab)
    load(fullfile(rootdir_ad,'channels_by_layers.mat'));
    
    if ismember(patname,channels_by_layers.Properties.VariableNames)
        if contains(lower(chlab),'sup')
            chans = [channels_by_layers.(patname){1:3}];
        elseif contains(lower(chlab),'inf')
            chans = [channels_by_layers.(patname){5:6}];
        end
    else
        chans = nolaypats_if(patname, chlab);
    end
end



if strcmp(period,'stim')
    if isempty(chlab)
        matnm = 'evtypes.mat';
        chans = [];
    else
        matnm = ['evtypes_' chlab '.mat'];
    end
elseif strcmp(period,'post')
    if isempty(chlab)
        matnm = 'evtypes_post.mat';
        chans = [];
    else
        matnm = ['evtypes_post_' chlab '.mat'];
    end
end

if exist(fullfile(rootdir_ad,matnm))==2
    load(fullfile(rootdir_ad,matnm))
else
    evtypes = struct;
end

minL = .8;
timewin = [0 .9];



% Only events with max. stim. intensity at a given stimulated electrode
% pair (for all MUA average)
maxintens = true;


evinx = select_evinx(rootdir_ad,patname,Events,maxintens,evtype,ti, minL, timewin, chans, period);
evnr = length(evinx);
evtypes.(patname).(evtype)(ti).maxintens = evinx;


% All events (for all AD - noAD MUA comparison)
maxintens = false;


evinx = select_evinx(rootdir_ad,patname,Events,maxintens,evtype,ti, minL, timewin, chans, period);
evnr = length(evinx);
evtypes.(patname).(evtype)(ti).all = evinx;



save(fullfile(rootdir_ad,matnm),'evtypes');
