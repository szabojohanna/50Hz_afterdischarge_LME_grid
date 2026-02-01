function [normepoch, sr,h_thumb] = load_filt_epoch(rootdir_ad,patname,eventinx,adnr,win,freqs,pktype,varargin)

if isempty(varargin)|| contains(varargin,'thumb_low');
    [~,~,ns_thumb,h_thumb,~,~,~] = load_synch_EEG_thumb(rootdir_ad,patname,eventinx,...
        'loadrec',{'thumb_low'});
    if isempty(ns_thumb); return; end;
    
    [~,~,ns_thumb,h_thumb,~,~,~] = load_EEG_thumb_peaks(rootdir_ad,patname,eventinx,adnr,...
        [],[],ns_thumb,h_thumb,[],[],pktype);
    
elseif contains(varargin,'EEG')
    [ns_thumb,h_thumb,~,~,~,~,~] = load_synch_EEG_thumb(rootdir_ad,patname,eventinx,...
        'loadrec',{'EEG'});
    
    if isempty(ns_thumb); return; end;
    
    [ns_thumb,h_thumb,~,~,~,~,~] = load_EEG_thumb_peaks(rootdir_ad,patname,eventinx,adnr,...
        ns_thumb,h_thumb,[],[],[],[],pktype);
    
elseif contains(varargin,'thumb_mua')
    [~,~,~,~,~,ns_thumb,h_thumb] = load_synch_EEG_thumb(rootdir_ad,patname,eventinx,...
        'loadrec',{'thumb_mua'});
    if isempty(ns_thumb); return; end;
    
    [~,~,~,~,~,ns_thumb,h_thumb] = load_EEG_thumb_peaks(rootdir_ad,patname,eventinx,adnr,...
        [],[],[],[],ns_thumb,h_thumb,pktype);
end



[normepoch, sr] = deal([]);
% if isempty(h_thumb.event) || (size(h_thumb.data,2)==24 && ti==2); close(ns_thumb); return; end;
if isempty(h_thumb.event) ; close(ns_thumb); return; end;
%% Select and filter epochs
sr = h_thumb.srate;
winlen = diff(win);

evtimes = [h_thumb.event{:,1}];
[evtimes,sortix] = sort(evtimes,'ascend');
if contains(pktype,'Small');
    evlabs = h_thumb.event(sortix,2);
    evlabs2 = cellfun(@str2double,evlabs);
    [pklabs, firstpk_inx ] = unique(evlabs2);
    evtimes(firstpk_inx) = [];
end
pknr = length(evtimes);

beg = (evtimes(1)/sr)+win(1)-winlen;
h_thumb = nswiew('inport',beg,winlen,h_thumb);

chnr = size(h_thumb.data,2);
signlen = size(h_thumb.data,1);

epoch = nan(signlen,chnr,pknr+2);
epoch(:,:,1) = h_thumb.data;

begg = (evtimes(end)/sr)+win(2);
h_thumb = nswiew('inport',begg,winlen,h_thumb);
epoch(:,:,end) = h_thumb.data;% buffer for filtering


for i = 1:pknr
    
    sig = [];
    
    begin = (evtimes(i)/sr)+win(1);
    h_thumb = nswiew('inport',begin,winlen,h_thumb);
    sig = h_thumb.data;
    
    epoch(:,:,i+1)=sig;
end

filtepoch = nan(size(epoch));
for j = 1:chnr
    longdat = reshape(epoch(:,j,:),[1 signlen*(pknr+2)]);
    
    
    %     Notch + bandpass
    [b a]=butter(1,[45/(sr/2) 55/(sr/2)],'stop');
    filtdat=filtfilt(b,a,longdat);
    
    [b a]=butter(1,[95/(sr/2) 105/(sr/2)],'stop');
    filtdat=filtfilt(b,a,filtdat);
    
    [b a]=butter(2,[freqs(1)/(sr/2) freqs(2)/(sr/2)],'bandpass');
    filtdat=filtfilt(b,a,filtdat);
    
    filtepoch(:,j,:) = reshape(filtdat,[size(epoch,1) 1 pknr+2]);
end
truepoch = filtepoch(:,:,2:end-1);



normepoch = nan(size(truepoch));
for j = 1:chnr
    longtru = reshape(truepoch(:,j,:),[1 signlen*(pknr)]);
    normepoch(:,j,:) = reshape( (longtru - nanmean(longtru)) / std( longtru,'omitnan'),[signlen 1 pknr]);
end


if sr==20000
    normepoch = downsample(normepoch,10);
    sr = 2000;
end

h_thumb.event = []; guidata(ns_thumb,h_thumb); close(ns_thumb)
