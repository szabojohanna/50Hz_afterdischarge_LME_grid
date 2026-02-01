%  ADpeaks_mark_save

%% Mark peaks on EEG rec!
% It will be converted and saved for thumb lfp and thumb mua files
corr_win = 0.03;
markEEG_saveEEG_Thumb(rootdir_ad,patname,corr_win)

%% Correct faulty labels on thumb if necessary + save for EEG if missing
corr_win = 0.1;
correct_onTHUMB_saveEEG_Thumb(rootdir_ad,patname,corr_win)

%% Mark small peaks

% patdir = fullfile(rootdir_ad,'Patients',patname);
% load(fullfile(patdir,'Events.mat'))
% AD_eventinxx = find(  [Events.AD]==1 & ~cellfun(@isempty,{Events.Thumbfile})  );

patname = 'Pt11';
% AD_eventinxx = [68 71 78 81 86];
AD_eventinxx = [71 78 81 86];
smallpeak_markEEG_saveETh(rootdir_ad,patname,AD_eventinxx)


function smallpeak_markEEG_saveETh(rootdir_ad,patname,AD_eventinxx)

patdir = fullfile(rootdir_ad,'Patients',patname);
load(fullfile(patdir,'Events.mat'))
adnr = 1;



pkfold_eeg = fullfile(rootdir_ad,'Patients',patname,'EEG','Peaktimes');
newpkfold_eeg = fullfile(rootdir_ad,'Patients',patname,'EEG','Small_peaktimes');
if ~isdir(newpkfold_eeg); mkdir(newpkfold_eeg); end;

for eventinx = AD_eventinxx
    
    fprintf('Ev%d, %s-%s\n',eventinx, Events(eventinx).StimChans{:});
    [ns_eeg,h_eeg,ns_thumb,h_thumb,ns_thumb_mua,h_thumb_mua] = load_synch_EEG_thumb(rootdir_ad,patname,eventinx,...
        'loadrec',{'EEG','thumb_low','thumb_mua'});
    h_eeg1 = h_eeg;
    srate = h_eeg.srate;
    EegTS = Events(eventinx).EegTS;
    thTS = Events(eventinx).ThumbTS;
    
    pkfile_eeg = dir(fullfile(pkfold_eeg,['*Ev' num2str(eventinx) '_manualcorr_AD'  num2str(adnr) '*']));
    
    if ~isempty(pkfile_eeg)
        evs = load(fullfile(pkfile_eeg.folder, pkfile_eeg.name));
    else
        fprintf('No events for Ev%d\n',eventinx);
        h_eeg.event = []; close(ns_eeg); close(ns_thumb); close(ns_thumb_mua);
        continue;
    end
    pks = sort([evs(:,end)]/srate,'ascend');
    
    h_eeg = set_nswiew_filter(h_eeg,[1 45], 'bandpass',2);
    h_eeg = set_nswiew_filter(h_eeg,[45 55], 'stop',2);
    h_eeg = set_nswiew_filter(h_eeg,[95 105], 'stop',2);
    guidata(ns_eeg,h_eeg)
    
    chanad = input('Give channel to correct peaks (ex: Gr2)','s');
    chnames = h_eeg.chnames;
    clean_chnames = clean_channames(chnames);
    chaninx = find(strcmpi(clean_chnames,chanad));
    polarity =  input('Give polarity of peaks (1 if upward, -1 if downward)'); % -1 if peaks are negative
    
    prepeak = 0.02;
    pkt = cell(length(pks),1);
    for k = 1:length(pks)
        h_eeg = nswiew('inport',pks(k)-prepeak,0.5,h_eeg);
        guidata(ns_eeg,h_eeg)
        
        ok = 0;
        while ok==0
            dat = h_eeg.data(:,chaninx);
            fig =figure;
            if polarity==1
                findpeaks(dat,'MinPeakProminence',25)
                [~, locs] =  findpeaks(dat,'MinPeakProminence',25);
            else
                findpeaks(-dat,'MinPeakProminence',25)
                [~, locs] =  findpeaks(-dat,'MinPeakProminence',25);
            end
            title(['PK' num2str(k)]);
            ok = input('Peaks OK? 1/0/9 if no small peaks at all\n');
            close(fig);
            if ok==0
                chanad = input('Give channel to correct peaks (ex: Gr2)','s');
                chaninx = find(strcmpi(clean_chnames,chanad));
                polarity =  input('Give polarity of peaks (1 if upward, -1 if downward)'); % -1 if peaks are negative
                
            end
        end
        if ok==1
            pkt{k} = round((pks(k)-prepeak)*srate)+locs;
        elseif ok==9
            pkt{k} = [];
        end
        
    end
    pktimes_eeg = cat(1,pkt{:});
    marksC = arrayfun(@(x) repmat(x, [length(pkt{x}), 1] ) ,1:length(pkt),'UniformOutput',0);
    marks = cat(1,marksC{:});
    
    
    
    h_eeg1.event = cell(length(pktimes_eeg),2);
    h_eeg1.event(:,1) = arrayfun(@(x) x,pktimes_eeg,'UniformOutput',0);
    h_eeg1.event(:,2) = arrayfun(@(x) num2str(x),marks,'UniformOutput',0);
    
    guidata(ns_eeg,h_eeg1);
    
    fprintf('Correct manually (delete small peaks at the end etc.)!\n'); keyboard;
    
    h_eeg = guidata(ns_eeg);
    pktimes_eeg_M = [h_eeg.event{:,1}]';
    marks_M = cellfun(@str2num, h_eeg.event(:,2));
    
    newpkfile_eeg = [pkfile_eeg.name(1:end-4) '_smallpeaks'];
    ev2writer_bh_j(pktimes_eeg_M,marks_M',{[newpkfold_eeg filesep], newpkfile_eeg});
    fpa = fullfile(patdir,'Thumb','Small_peaktimes'); if ~isdir(fpa); mkdir(fpa); end;
    %%
    convert_EEGtimes2Th(h_eeg,h_thumb,h_thumb_mua,pktimes_eeg,EegTS,thTS,marks',fpa,newpkfile_eeg)
    
    h_eeg.event = []; close(ns_eeg); close(ns_thumb); close(ns_thumb_mua);
end
end




%--------------------------------------------------------------------------
function markEEG_saveEEG_Thumb(rootdir_ad,patname, corr_win)
patdir = fullfile(rootdir_ad,'Patients',patname);

load(fullfile(patdir,'Events.mat'))
AD_eventinxx = find(  [Events.AD]==1 & ~cellfun(@isempty,{Events.Thumbfile})  );


fpa = fullfile(patdir,'EEG','Peaktimes'); if ~isdir(fpa); mkdir(fpa); end;

for eventinx = AD_eventinxx
    
    % Event name
    eventname = Events(eventinx).Orig_Label;
    evnm = eventname(1:strfind(eventname,'mA')+1);
    evfieldnm = ['Ev' num2str(eventinx)];
    disp(eventname)
    if ~contains(upper(eventname),'GR')
        continue;
    end
    
    
    fnmd = dir(fullfile(fpa,['*_' evnm(1:end-4) '_' evfieldnm '_manualcorr_AD1.ev2']));
     moread = 1; 
        adnr = 1;
    if exist(fullfile(fpa,fnmd.name))~=2
        
        %% Load synchronized EEG + thumb
        [ns_eeg,h_eeg,ns_thumb,h_thumb,ns_thumb_mua,h_thumb_mua,~] = load_synch_EEG_thumb(rootdir_ad,patname,eventinx,...
            'loadrec',{'EEG','thumb_low','thumb_mua'});
        
        %         keyboard;
        okp = input('If peaks to mark, press 1, else any key');
        if okp~=1
            continue;
        end
        %% Mark on EEG, save for thumb
        
        
       
        h_thumb_mua = guidata(ns_thumb_mua);
        h_eeg = guidata(ns_eeg);
        
        % Start point
        EegTS = h_eeg.time(1);
        thTS = h_thumb_mua.time(1);
        
        % File names
        thfname2 = h_thumb_mua.fname(1:end-4);
        eegfname = h_eeg.fname(1:end-4);
        
        
        
        
        %%
        while moread ==1
            %% Are there events?
            
            
            h_thumb_mua = guidata(ns_thumb_mua);
            h_eeg = guidata(ns_eeg);
            
            % Save peaks
            
            fnm = [thfname2 '_' eegfname(1:end-4) '_' evnm(1:end-4) '_' evfieldnm '_manualcorr_AD' num2str(adnr)];
            fnm = fnm(~ismember(fnm,'.'));
            fil = dir([fnm '*']);
            
            
            %                 mes = msgbox('Mark event peaks labeled as 1');
            keyboard;
            
            %%  Get peaks on EEG
            h_eeg = guidata(ns_eeg);  % Load events, if exist
            
            fnm = [thfname2 '_' eegfname '_' evnm(1:end-4) '_' evfieldnm '_manualcorr_AD' num2str(adnr)];
            fnm = fnm(~ismember(fnm,'.'));
            
            if ~isempty(fil)
                fprintf('Peaktimes already saved! Loaded.')
            else
                %         time = h_eeg.time;
                pkev = find(ismember(h_eeg.event(:,2),'1'));
                
                
                evtimes = cell2mat(h_eeg.event(pkev,1))/h_eeg.srate;
                
                
                
                
                %% Correct marked peaks + save
                
                chanad = input('Give channel to correct peaks (ex: Gr2)','s');
                chnames = h_eeg.chnames;
                clean_chnames = clean_channames(chnames);
                chaninx = find(strcmpi(clean_chnames,chanad));
                
                
                polarity =  input('Give polarity of peaks (1 if upward, -1 if downward)'); % -1 if peaks are negative
                [pktimes_eeg,h_eeg] = correct_events(corr_win,chaninx,h_eeg,[],polarity,evtimes);
                
                guidata(ns_eeg,h_eeg)
                % Save EEG peaks
                fpa = fullfile(patdir,'EEG','Peaktimes'); if ~isdir(fpa); mkdir(fpa); end;
                ev2writer_bh(pktimes_eeg,1,{[fpa filesep], fnm});
                %%
            end
            
            
            
            
            %% Save THUMB peaks
            fpat = fullfile(patdir,'Thumb','Peaktimes'); if ~isdir(fpat); mkdir(fpat); end;
            
            convert_EEGtimes2Th(h_eeg,h_thumb,h_thumb_mua,pktimes_eeg,EegTS,thTS,1,fpat,fnm)
            
            %                 close(mes);
            
            
            answer = questdlg('More AD?','',{'No','Yes'});
            if strcmp(answer,'No'); moread = 0; end;
            adnr = adnr+1;
        end
        
        
    else
        eeg_pkts = load(fullfile(fpa,fnmd.name));
        pktimes_eeg = eeg_pkts(:,end);
        [ns_eeg,h_eeg,ns_thumb,h_thumb,ns_thumb_mua,h_thumb_mua,~] = load_synch_EEG_thumb(rootdir_ad,patname,eventinx,...
            'loadrec',{'EEG','thumb_low','thumb_mua'});
        keyboard;
        h_eeg = guidata(ns_eeg);
        h_thumb_mua = guidata(ns_thumb_mua);
        
        EegTS = h_eeg.time(1);
        thTS = h_thumb_mua.time(1);
        
        fpat = fullfile(patdir,'Thumb','Peaktimes'); if ~isdir(fpat); mkdir(fpat); end;
        
        eegfname = h_eeg.fname(1:end-4);
        thfname2 = h_thumb_mua.fname(1:end-4);
        fnm = [thfname2 '_' eegfname '_' evnm(1:end-4) '_' evfieldnm '_manualcorr_AD' num2str(adnr)];
        fnm = fnm(~ismember(fnm,'.'));
        convert_EEGtimes2Th(h_eeg,h_thumb,h_thumb_mua,pktimes_eeg,EegTS,thTS,1,fpat,fnm)
    end
    % Close nswiew
        h_thumb.event = {}; guidata(ns_thumb,h_thumb);close(ns_thumb)
        h_thumb_mua.event = {}; guidata(ns_thumb_mua,h_thumb_mua);close(ns_thumb_mua)
        h_eeg.event = {}; guidata(ns_eeg,h_eeg);close(ns_eeg)
end
end





%--------------------------------------------------------------------------
function correct_onTHUMB_saveEEG_Thumb(rootdir_ad,patname,corr_win,eventinx,ns_thumb,h_thumb)
%% Correct event peaks if necessary AGAIN according to thumb

pkfold = fullfile(rootdir_ad,'Patients',patname,'Thumb','Peaktimes');
pkfiles = dir(fullfile(pkfold,['*Ev' num2str(eventinx) '_manualcorr_AD' num2str(adnr) 'e*']));

fnm = pkfiles.name(1:end-4);
h_thumb.event = {};
h_thumb = load_ev2event(ns_thumb,h_thumb,pkfold,fnm);

% if ismember(patname,{'LuEn','GaJu','MoJo'}) && h_thumb.srate ==2000;
%     h_thumb.event(:,1) = cellfun(@(x) round(x/10),h_thumb.event(:,1),'UniformOutput',0); % MUA sr alapján van lementve... (20000 Hz)
% end
pknr = size(h_thumb.event(:,1),1);

guidata(ns_thumb,h_thumb);

%%
% corr_win = 0.1;
chaninx  = input('Give channel index for correction!');
[corr_evtimes,h_thumb] = correct_events(corr_win,chaninx,h_thumb,1:pknr,-1);
guidata(ns_thumb,h_thumb);

keyboard;
pktimes_thumb = cell2mat(h_thumb.event(:,1));
ev2writer_bh(pktimes_thumb,1,{[pkfold filesep ] , fnm});

%% SAVE EVENTS FOR MUA (if there is sep mua file)
fnm = [pkfiles.name(1:end-4)] ;
pktimes_thumb = cell2mat(h_thumb_mua.event(:,1));
ev2writer_bh(round(pktimes_thumb/10),1,{[pkfold filesep ] , fnm});



fnm = [pkfiles.name(1:end-5) 'u'] ;
pktimes_thumb = cell2mat(h_thumb.event(:,1));
ev2writer_bh(round(pktimes_thumb*10),1,{[pkfold filesep ] , fnm});

%% Missing EEG peaktimes ?!
pktimes_eeg = floor(((pktimes_thumb/h_thumb.srate)-thTS+EegTS)*h_eeg.srate);

% Save EEG peaks
fpa = fullfile(patdir,'EEG','Peaktimes'); if ~isdir(fpa); mkdir(fpa); end;
fnm = [thfname2 '_' eegfname(1:end-4) '_' evnm(1:end-4) '_' evfieldnm '_manualcorr_AD' num2str(adnr)'];
fnm = fnm(~ismember(fnm,'.'));

ev2writer_bh(pktimes_eeg,1,{[fpa filesep], fnm});
end


function convert_EEGtimes2Th(h_eeg,h_thumb,h_thumb_mua,pktimes_eeg,EegTS,thTS,k,fpa,fnm)

%% Save THUMB peaks
pktimes_thumb = floor(((pktimes_eeg/h_eeg.srate)-EegTS+thTS)*h_thumb_mua.srate);
% h_thumb_mua.event  ={};
% for i = 1:length(pktimes_thumb)
%     h_thumb_mua.event{i,1} = pktimes_thumb(i);
%     h_thumb_mua.event{i,2} = num2str(k);
% end

% guidata(ns_thumb_mua,h_thumb_mua);


if h_thumb_mua.srate==20000 && h_thumb.srate==2000
    pktimes_thumbU = pktimes_thumb;
    pktimes_thumbE = round(pktimes_thumb/10);
else
    
    pktimes_thumbU = pktimes_thumb;
    pktimes_thumbE = pktimes_thumb;
end


ev2writer_bh_j(pktimes_thumbU,k,{[fpa filesep ] , [fnm 'u']});
ev2writer_bh_j(pktimes_thumbE,k,{[fpa filesep ] , [fnm 'e']});
end
