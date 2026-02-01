function save_BLOCKsegment2spikedetection(rootdir_ad,patname,choi,savedata);


patdir = fullfile(rootdir_ad,'Patients',patname);
load(fullfile(patdir,'Events.mat'))

[Th_TS,ev_order] = sort([Events.ThumbTS],'ascend'); % legtöbbször sorrendben, de biztonság kedvéért
Thfiles = {Events.Thumbfile};
Thfiles(cellfun(@isempty,Thfiles)) = {''};
U_thfiles = unique(Thfiles);

buff = 0.01; % sec - buffer before next event start;
%%
evix = 0;

inx = 1;

EVNR = length(Events);
if exist(fullfile(rootdir_ad,'SUA_Analysis',patname,'EventInfo.mat'))~=2
    EventInfo = table(...
        cell(EVNR, 1), ... % Intrastim_BaselineInx
        cell(EVNR, 1), ... % Intrastim_StimInx
        cell(EVNR, 1), ... % Intrastim_PostInx
        cell(EVNR, 1), ... % ADperiodTime
        cell(EVNR, 1), ... % ADperiodInx
        'VariableNames', {'Intrastim_BaselineInx', 'Intrastim_StimInx', 'Intrastim_PostInx', 'ADperiodTime', 'ADperiodInx'});
else
    load(fullfile(rootdir_ad,'SUA_Analysis',patname,'EventInfo.mat'))
end

EventData = cell(length(Events),1);
%%
for u = 1:length(U_thfiles) % sorrendben fileokon, hogy id?rendbe lehessen tenni eseményeket
    thf = U_thfiles{u};
    if isempty(thf); continue; end
    thf_evs = find(ismember(Thfiles,thf));
    thf_ts = [Events(thf_evs).Stim_start];
    
    ev_datapoins = cell(length(Events),4);
    for eventinx = thf_evs
        
        evix = evix+1;
        if isempty(Events(eventinx).Stim_start);
            fprintf('No Stim start, %s %d\m',patname, eventinx ); continue;
        end;
        
        
        isAD = Events(eventinx).AD;
        matdir =  fullfile(rootdir_ad,'Figures','Micro','MUA','MUA_blocks');
        matnm = ['MUA_stimblocks_' patname '_Ev' num2str(eventinx) '_AD' num2str(isAD) '.mat'];
        try
            load(fullfile(matdir,matnm))
        catch
            fprintf('No EV data start, %s %d\m',patname, eventinx ); continue;
        end
        
        
        % Concatenate intrastim. artefact-free blocks
        stimpeak_ts{1} = MUA_stimblocks.params.Bpk_ts_inblock;
        stimpeak_ts{2} = MUA_stimblocks.params.pk_ts_inblock;
        stimpeak_ts{3} = MUA_stimblocks.params.Ppk_ts_inblock;
        
        %         pk_ts = cat(1,stimpeak_ts{:});
        
        block_ms = MUA_stimblocks.params.block_ms;
        plot_ms = MUA_stimblocks.params.mua_ms;
        
%                 [ns_eeg,h_eeg,ns_lfp,h_lfp,ns,h,~] = load_synch_EEG_thumb(rootdir_ad,patname,eventinx,...
%                     'loadrec',{'EEG','thumb_low','thumb_mua'});
        
        [~,~,~,~,ns,h,~] = load_synch_EEG_thumb(rootdir_ad,patname,eventinx,...
            'loadrec',{'thumb_mua'});
        sr = h.srate;
        
        blockL = diff(block_ms)/1000;
        blockL_dp = blockL*sr;
        
        
        blockdat = [];bends = nan(3,1);
        [blinx, bl_ts] = deal(cell(3,1));
        for z = 1:3
            pk_ts = stimpeak_ts{z}; % peak of stim. artefact
            pknr = length(pk_ts); % nr. of stim pulses
            
        [blinx{z}, bl_ts{z}] = deal(nan(pknr,2)); jj = 1; st = 1;
            for k = 1:pknr
                block_start = pk_ts(k)+block_ms(1)/1000;
                block_end = block_start+blockL;
                h = nswiew('inport',block_start,blockL,h); % window around a single stim pulse
                guidata(ns,h)
                dat = h.data; %
                blockdat = cat(1,blockdat, dat);
                blinx{z}(jj,:) = [st st+size(dat,1)-1];
                st =st+ size(dat,1);
                bl_ts{z}(jj,:) = [block_start block_end];
                jj = jj+1;
            end
            bends(z) = size(blockdat,1);
        end
        
        if savedata
            ch_blockdat = blockdat(:,choi);
        end
        
        % Add post-stim period until next stim. event
        
        ev_th_start = Events(eventinx).Stim_start;
        ev_th_end = Events(eventinx).Stim_end;
        nextev_th_start = min(thf_ts(thf_ts>ev_th_start)); % legközelebbi következ? esemény kezdete, ugyanabban a fileban
        
        if ~isempty(nextev_th_start)
            ev_st = ev_th_end+buff;
        else
            nextev_th_start =  min(h.maxsec-0.1,ev_st+30); % if there is no next event in the file, get 30 sec (or until the file ends)
        end
        
        ev_L = nextev_th_start-buff-ev_st;
        
        h = nswiew('inport',ev_st,ev_L,h);
        guidata(ns,h)
        dat = h.data;
        if savedata
            ch_blockdat = cat(1,ch_blockdat, dat(:,choi));
            EventData{eventinx} = ch_blockdat;
        end
        
        
        
        
        %         EventNr{eventinx} = eventinx;
        %         Intrastim_BlockStart_Sec{eventinx} = pk_ts+block_ms(1)/1000; % time of block start (in sec)
        %         Intrastim_BlockInx{eventinx} = blinx; % block index corresponding to each datapoint
        %         Intrastim_EventInx{eventinx} = (inx:inx+size(blockdat,1)-1)'; % indices of event within array of events (=within the saved EventData)
        %         inx=inx+size(blockdat,1);
        
        
        EventInfo.Intrastim_BaselineInx{eventinx} = [inx; inx+bends(1)-1];
        EventInfo.Intrastim_StimInx{eventinx} = [inx+bends(1); inx+bends(2)-1];
        EventInfo.Intrastim_PostInx{eventinx} = [inx+bends(2); inx+bends(3)-1];
        inx=inx+bends(3);
        
        EventInfo.Intrastim_BasLim_TS{eventinx} = bl_ts{1};
        EventInfo.Intrastim_BasLim_Inx{eventinx} = blinx{1};
        
        EventInfo.Intrastim_StimLim_TS{eventinx} = bl_ts{2};
        EventInfo.Intrastim_StimLim_Inx{eventinx} = blinx{2};
        
        EventInfo.Intrastim_PostLim_TS{eventinx} = bl_ts{3};
        EventInfo.Intrastim_PostLim_Inx{eventinx} = blinx{3};

                
        EventInfo.ADperiodTime{eventinx} = [ev_st; ev_st+ev_L]; % time vector of potential AD period
        EventInfo.ADperiodInx{eventinx} = [inx; inx+size(dat,1)-1]; % indices of potential AD period
        inx=inx+size(dat,1);
        close(ns);
        
        save(fullfile(rootdir_ad,'SUA_Analysis',patname,'EventInfo.mat'),'EventInfo');
        
    end
end
% EventInfo = table(EventNr,Intrastim_BlockStart_Sec,Intrastim_EventInx,Intrastim_BlockInx,ADperiodTime,ADperiodInx);
% EventInfo = table(Intrastim_BaselineInx,Intrastim_StimInx, Intrastim_PostInx, ADperiodTime,ADperiodInx);
allevdata = cat(1,EventData{:});

% allpks = sort(cat(1,Intrastim_BlockStart_Sec{:}));
% allt = sort(cat(1,ADperiodTime{:}));
% time_start_sec = min(allpks);
% time_end_sec = max(allt);

if savedata
    for c = 1:length(choi)
        chan =choi(c);
        
        savedir = fullfile(rootdir_ad,'SUA_Analysis',patname,['Chan' num2str(chan)]);
        if ~isdir(savedir); mkdir(savedir); end;
        data = allevdata(:,c);
        
        evfnm = ['Chan' num2str(chan) '_AllEvents_Raw.mat'];
        %     [savedir, fnm2] = save_fnm(fullfile(patdir,'Thumb'),evfnm,time_start_sec, time_end_sec,chan);
        %     save(fullfile(savedir,[fnm2 '.mat']),'data','sr','EventInfo');
        
        save(fullfile(savedir,evfnm),'data','sr');
        
    end
end

save(fullfile(rootdir_ad,'SUA_Analysis',patname,'EventInfo.mat'),'EventInfo');
end