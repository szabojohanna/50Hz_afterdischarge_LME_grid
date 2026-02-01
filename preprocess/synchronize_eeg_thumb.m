% synchronize_eeg_thumb

% look through data & events in nswiew before proceeding

  
% Detect matching triggers visually 
ns1 = nswiew; % Import EEG
ns2 = nswiew; % Import Thumb

disp('Import macro & micro data, mark triggers/ artefacts used for synch! If ready, press Continue!')
keyboard;


%% 
thumb_eeg_sync_ns(rootdir_ad,patname,ns1,ns2)


%% if Sync.mat has already been saved
ev_coreg_eeg2thumb(rootdir_ad,patname) % saves Events.mat (melyik esemény melyik thumb-/eegfileban van benne)

%% 
% find_trigger_pattern % in progress





function ev_coreg_eeg2thumb(rootdir_ad,patname)


% Matching timestamps of stimulation events using linear regression
% saves Events struct with fields:
%      Label = event label
%      TS = timestamp of event in thumb file
%      Thumbfile = name of thumb file which contains the event
%      EEGfile = name of EEG file which contains the event
%%

patdir = fullfile(rootdir_ad,'Patients',patname);


load(fullfile(patdir,'Sync.mat'))
load(fullfile(patdir,'Events.mat'));


ns1 = nswiew; h1 = guidata(ns1);
ns2 = nswiew; h2 = guidata(ns2);


eegfiles = fieldnames(Sync);
for ei = 1:length(eegfiles)
    
    eegfname = eegfiles{ei};
    
    
    if contains(patdir,'Pt3'); eeg_ext = '.vhdr'; else eeg_ext = '.TRC'; end;
        h11 = nswiew('inport_menu_Callback',h1, [], h1, { fullfile(patdir,'EEG\'), [eegfname eeg_ext]});

        guidata(ns1, h11);
        h1 =guidata(ns1);
        
        if exist(fullfile(patdir,'EEG',[eegfname '_nswev.mat']))==2
            newevents = load(fullfile(patdir,'EEG',[eegfname '_nswev.mat']));
            h1.event = newevents.events;
        else
            continue;
        end
    
    thumbfiles = fieldnames(Sync.(eegfname));
    for ti = 1:length(thumbfiles)
        
        thfname = thumbfiles{ti};
        thfname2 = thfname(3:end); sL = strfind(thfname2,'_');
        
         tn = dir(fullfile(patdir,'Thumb','*.cnt'));
         if isempty(tn); 
             tn = dir(fullfile(patdir,'Thumb','*.rhd')); ext = 'rhd'; 
             h2.rhddata = []; h2.rhdtime = [];
         else; ext = 'cnt'; 
         end;
         
        if  length(tn)>1; tnm = tn(1).name;else; tnm = tn.name; end;
        
        if any(strfind(tnm,'_')) && strcmp(ext,'cnt') ;
            if strfind(tnm,'#'); thfname2(sL(2)) = '#'; end;
            if strfind(tnm,'@');
                try thfname2(sL(3)) = '@';
                catch thfname2(sL(2)) = '@';
                end;
            end;
            
        elseif ~any(strfind(tnm,'_')) && strcmp(ext,'cnt') ;
            if strfind(tnm,'#'); thfname2(sL(1)) = '#'; end;
            if strfind(tnm,'@'); thfname2(sL(2)) = '@';end;
        end
    
    
        h22 = nswiew('inport_menu_Callback',h2, [], h2, { fullfile(patdir,'Thumb\'), [thfname2 '.' ext]});
        guidata(ns2, h22); h22 = [];
        h2 =guidata(ns2);
        
        
        ecs = Sync.(eegfname).(thfname).eeg_markerindex';
        ecd = Sync.(eegfname).(thfname).thumb_markerindex';
        
        eix = find(contains(h1.event(:,2),'0.0Hz')); % index of stimulation events in EEG event list
        
        e2c = cell2mat(h1.event(eix,1))/h1.srate; % timestamp of stimulation events in EEG recording
        
        % Linear regression
        X=[ones(size(ecs)),ecs];
        linear_coeffs=X\ecd;
        
        ecdcontrol=ecs*linear_coeffs(2)+linear_coeffs(1);
        
        d=ecd-ecdcontrol;
        
        disp(['Dispersion: ',num2str(std(d))]);
        
        evo=e2c*linear_coeffs(2)+linear_coeffs(1);
        
        th_end = h2.maxsec;
        
        thinx = find(evo>0&evo<th_end);
        
        if isempty(thinx)
            fprintf('No events in this thumb\n')
            continue
        end
        
        thevo = evo(thinx); % index/timestamp of stimulation events in thumb recording
        eegevo = e2c(thinx);
        evlabs = h1.event(eix,2);
        thevlabs = evlabs(thinx);
        
        evonr = length(thevo);
        
        % save event list in thumb file
        h2.event(1:evonr,1) = num2cell(thevo*h2.srate);
        h2.event(1:evonr,2) = thevlabs;
        guidata(ns2,h2);
        
        ev2writer2(h2,'all',fullfile(patdir,'Thumb',[thfname2 '_' eegfname]));
        
        
        if ~isempty([Events.Label])
            en = length(Events);
        else
            en = 0;
        end
        
        
        
        Events_EEGfileinx = find(ismember({Events.EEGfile},[eegfname eeg_ext])); % find events in current EEGfile stored in Events.mat
        Events_thumbevents = Events_EEGfileinx(dsearchn([Events(Events_EEGfileinx).EegTS]',eegevo)); % find events in currenct Thumb file in Events.mat
    
        for j = 1:length(Events_thumbevents)
            Eventsj = Events_thumbevents(j);
            Events(Eventsj).ThumbTS = thevo(j);
            Events(Eventsj).Thumbfile = thfname;
        end
        
    end
end


% Add name of stimulation channels + stimulation intensity
for i = 1:size(Events,2)
    evo = Events(i).Label;
    lstr = strfind(evo,'-');pstr = strfind(evo,'.');
    sp = find(isspace(evo(1:pstr)));
    stim_chnm = {evo(1:lstr-1), evo(lstr+1:sp-1)};
    stim_chnm = cellfun(@(x) x(~isspace(x)),stim_chnm,'UniformOutput',false);
    Events(i).StimChans = stim_chnm;
    
    mastr = strfind(evo,'mA');
    ma = str2double(evo(sp+1:mastr-1));
    Events(i).StimIntensity = ma;
    
end

close(ns1); close(ns2)
save(fullfile(patdir,'Events.mat'),'Events')

end