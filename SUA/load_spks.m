function   [cl_evix,clus_inx,clus_ts_s,ev_ix,EventInfo, basend,stimend,postend,...
    cl_evixB,ev_ixB, cl_evixS,ev_ixS, cl_evixP,ev_ixP] = load_spks(rootdir_ad,Uname,sr,eventinx,evperiod,varargin);
% Loads spikes within stim period/ ad period
% Input:
%   Uname - char. array of unit label
%   sr - double, sampling rate (ex 20000)
%   eventinx - double, index of event (based on Events.mat)
%   evperiod - char. array, period to load; 
%               'intrastim' | 'adperiod' | 'adevents'


if strcmp(evperiod,'intrastim')
    block_inx = varargin{1};
else
    win = varargin{1};
end



[cl_evix,clus_inx,clus_ts_s,ev_ix,EventInfo, basend,stimend,postend] = deal([]);


[basend,stimend,postend] = deal([]);
[patname, ~, ~] = unpack_Uname(Uname);
% savedir = fullfile(rootdir_ad,'SUA_Analysis',patname,['Chan' num2str(chan)]);
% fnm = ['Chan' num2str(chan) '_AllEvents_Raw'];
% % load(fullfile(savedir,[fnm '.mat']));
% load(fullfile(savedir,[fnm '_spikes.mat']))
% load(fullfile(savedir,['times_' fnm '.mat']))
% 
% sp_ts_s = index/1000; % timestamp of spikes in sec
% spk_inx = round(sp_ts_s*sr); % indices of spikes
% 
% clus_inx = spk_inx(cluster_class(:,1)==clnr);
% clus_ts_s = spk_inx(cluster_class(:,1)==clnr);

load(fullfile(rootdir_ad,'SUA_Analysis','SUA_data_corrected'));
Unames = unit_table.Unit_name;
uix = ismember(Unames, Uname);
clus_inx = unit_table.Cluster_inx{uix};

if ~isempty(eventinx)
    load(fullfile(rootdir_ad,'SUA_Analysis',patname,'EventInfo.mat'));
   
    
    if strcmp(evperiod,'intrastim')
        
        if isempty(block_inx) % all data blocks included
            ev_LB = EventInfo.Intrastim_BaselineInx{eventinx};
            ev_ixB = ev_LB(1):ev_LB(2);
            ev_LS = EventInfo.Intrastim_StimInx{eventinx};
            ev_ixS = ev_LS(1):ev_LS(2);
            ev_LP = EventInfo.Intrastim_PostInx{eventinx};
            ev_ixP = ev_LP(1):ev_LP(2);
            
            
        else % selected data blocks (matching)
            
            match_blinx = EventInfo.Intrastim_BasLim_Inx{eventinx}(block_inx{1},:);
            mb=arrayfun(@(x) match_blinx(x,1):match_blinx(x,2),1:size(match_blinx,1),'UniformOutput',0);
            ev_ixB = [mb{:}]+EventInfo.Intrastim_BaselineInx{eventinx}(1)-1;
            
            match_blinx = EventInfo.Intrastim_StimLim_Inx{eventinx}(block_inx{2},:);
            mb=arrayfun(@(x) match_blinx(x,1):match_blinx(x,2),1:size(match_blinx,1),'UniformOutput',0);
            ev_ixS = [mb{:}]+EventInfo.Intrastim_StimInx{eventinx}(1)-1;
            
            match_blinx = EventInfo.Intrastim_PostLim_Inx{eventinx}(block_inx{3},:);
            mb=arrayfun(@(x) match_blinx(x,1):match_blinx(x,2),1:size(match_blinx,1),'UniformOutput',0);
            ev_ixP = [mb{:}]+EventInfo.Intrastim_PostInx{eventinx}(1)-1;
            
        end
        
%        
%         ev_ix = [ev_ixB ev_ixS ev_ixP];
%         cl_evix= intersect(ev_ix,clus_inx);
        
        
        [cl_evixB] = intersect(ev_ixB,clus_inx);
        [cl_evixS] = intersect(ev_ixS,clus_inx);
        [cl_evixP] = intersect(ev_ixP,clus_inx);
        
              
        cl_evix = [cl_evixB cl_evixS cl_evixP];
        
        basend = length(ev_ixB); stimend = basend+length(ev_ixS); postend = stimend+length(ev_ixP);

       
        
    elseif strcmp(evperiod,'adperiod')
        ev_ix = [EventInfo.ADperiodInx{eventinx}];
        cl_evix= intersect(ev_ix,clus_inx);
        
    elseif strcmp(evperiod,'adevents')
        adevs = load_ADpeaks(rootdir_ad,patname,eventinx);
        if isempty(adevs); return; end;
        aT = EventInfo.ADperiodTime{eventinx};
        adtime = aT(1):1/sr:aT(2);
        aX = EventInfo.ADperiodInx{eventinx};
        ev_ix = aX(1):aX(2);
        ad_ts = adevs(:,end)/sr;
        adix = dsearchn(adtime',ad_ts);
        anr = length(adix);
        
        win_dp = win*sr;
        ad_st = adix+ev_ix(1)-1;
        ad_w = ad_st+win_dp;
        for k = 1:anr
            [~,cl_evix{k},~] =intersect(ad_w(k,1):ad_w(k,2),clus_inx);
        end
        basend = abs(win_dp(1));
    end
    
else
    cl_evix = clus_inx;
end
end


function adevs = load_ADpeaks(rootdir_ad,patname,eventinx)
adevs = [];
pkfold = fullfile(rootdir_ad,'Patients',patname,'Thumb','Peaktimes');
pkfiles = dir(fullfile(pkfold,['*Ev' num2str(eventinx) '_manualcorr_AD1*']));

if length(pkfiles)>1
    pkfiles_mua = dir(fullfile(pkfold,['*Ev' num2str(eventinx) '_manualcorr_AD1*u.ev2']));
else
    pkfiles_mua = pkfiles;
end
if isempty(pkfiles_mua)
    return
end
adevs = load(fullfile(pkfold,pkfiles_mua.name));
end