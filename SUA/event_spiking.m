function event_spiking
%%
load(fullfile(rootdir_ad,'SUA_Analysis','SUA_data'));
Unames = unit_table.Unit_name;

patdir = fullfile(rootdir_ad,'Patients',patname);
load(fullfile(patdir,'Events.mat'))

evs = find([Events.AD]==1);
binsize = 0.01;
for h = 1:length(Unames)
    Uname = Unames{h};
    [patname, chan, clnr] = unpack_Uname(Uname);
    
    if chan<25; ti = 1; else; ti = 2; end;
    
    %     pr.evtype = evtype; pr.patname = patname; pr.period = 'stim';
    %     pr.rootdir_ad = rootdir_ad; pr.maxintens = false;
    %     evs= evs_from_tbl(pr,ti)';
    %     evs = sort(evs);
    
    for k = 1:length(evs)
        eventinx = evs(k);
        
        if isempty(Events(eventinx).AD_ch_order1); fprintf('No %d ev\n',eventinx); continue; end;
        gridchan = Events(eventinx).AD_ch_order1{1}{1};
        %%
        stimevent_spiking_oneUnit_oneEvent(rootdir_ad,Uname,sr,eventinx,true,binsize)
        
        %%
        %         AD_spiking_oneUnit_oneEvent(rootdir_ad,Uname,sr,eventinx,gridchan, true)
        
        
    end
end

%%

for h = 1:length(Unames)
    Uname = Unames{h};
    [patname, chan, clnr] = unpack_Uname(Uname);
    
    if chan<25; ti = 1; else; ti = 2; end;
    
    pr.evtype = evtype; pr.patname = patname; pr.period = 'stim';
    pr.rootdir_ad = rootdir_ad; pr.maxintens = false;
    evs= evs_from_tbl(pr,ti)';
    evs = sort(evs);
    
    for k = 1:length(evs)
        eventinx = evs(k);
        
        if isempty(Events(eventinx).AD_ch_order1); fprintf('No %d ev\n',eventinx); continue; end;
        gridchan = Events(eventinx).AD_ch_order1{1}{1};
        
        
        AD_spiking_oneUnit_evtype(rootdir_ad,Uname,sr,evss)
        
        
    end
end

%%
binsize = 0.001;
ccx(rootdir_ad,Unm1,Unm2,sr,eventinx,binsize)


%% Count spikes

Uname = 'Pt3_3_3'; % 'Pt3_3_3'
evss = find(~cellfun(@isempty,{Events.Stim_start}));
AD_spiking_oneUnit_evtype(rootdir_ad,Uname,sr,evss)

end

%--------------------------------------------------------------------------
function ccx(rootdir_ad,Unm1,Unm2,sr,eventinx,binsize)


%%
sp_train = cell(2,1);
for j = 1:2
    switch j; case 1; Uname = Unm1; case 2; Uname = Unm2; end;
    
    [patname, chan, clnr] = unpack_Uname(Uname);
    
    patdir = fullfile(rootdir_ad,'Patients',patname);
    load(fullfile(patdir,'Events.mat'))
    gridchan = Events(eventinx).AD_ch_order1{1}{1};
    
    %%
    [spike_train_MX, edges, psth_FR_s,centers, LFPavg,lfp_time, ...
        GRavg,grid_time] = AD_spiking_oneUnit_oneEvent(rootdir_ad,Uname,sr,eventinx,gridchan, false);
    sp_ad = reshape(spike_train_MX',[size(spike_train_MX,1)* size(spike_train_MX,2),1]);
    
    %%
    [sp_stim, edges, fr,centers] = stimevent_spiking_oneUnit_oneEvent(rootdir_ad,Uname,sr,eventinx,false,binsize);
    sp_train{j} = cat(1,sp_stim',sp_ad);
end
%%
figure;
% subplot(2,1,1);
% plot(centers, psth_FR_s1); hold on;
% plot(centers, psth_FR_s2); hold on;
[ccg, lags] = xcorr( sp_train{1}, sp_train{2}, 0.3/binsize, 'coeff');
figure;
bar(lags*binsize, ccg);
xlabel('Lag (sec)');
ylabel('Correlation');
title('Cross-Correlogram');
end


%--------------------------------------------------------------------------
function [spike_train, edges, fr,centers] = stimevent_spiking_oneUnit_oneEvent(rootdir_ad,Uname,sr,eventinx,isplot,binsize)


[patname, chan, clnr] = unpack_Uname(Uname);

figdir = fullfile(rootdir_ad,'SUA_Analysis','Figures','intrastim_spiking',patname);
if ~isdir(figdir); mkdir(figdir); end;

% binsize = 0.01;
[cl_evix,~,~,ev_ix,EventInfo,basend,stimend,postend] = load_spks(rootdir_ad,Uname,sr,eventinx,'intrastim');
[fr,edges,spike_train] = calc_ratemod(cl_evix,ev_ix,sr,binsize);
centers = edges(1:end-1) + diff(edges)/2;

if isplot
    fig = figure;
    
    plot(centers, fr,'Color','k');
    plot(centers,fr,'Color','k');
    yL = ylim; hold on;
    line(([basend basend])./sr,yL,'Color','k')
    line(([stimend stimend])./sr,yL,'Color','k')
    xlim([0 (postend)./sr])
    title({Uname,['Spike count, bin size:' num2str(binsize) 'dp']});
    ylabel('Rate (Hz)');
    % xlabel('Time (s)')
    
    fnm = ['Ev' num2str(eventinx) '_ch' num2str(chan) '_CL' num2str(clnr)];
    
    saveas(fig,fullfile(figdir,['ratemod_STIMperiod_' fnm '.jpg']));
    saveas(fig,fullfile(figdir,['ratemod_STIMperiod_' fnm '.fig']));
    close(fig)
end
end


function [spike_train_MX, edges, psth_FR_s,centers, LFPavg,lfp_time, GRavg,grid_time] = AD_spiking_oneUnit_oneEvent(rootdir_ad,Uname,sr,eventinx,gridchan, isplot)


[patname, chan, clnr] = unpack_Uname(Uname);

figdir = fullfile(rootdir_ad,'SUA_Analysis','Figures','AD_spiking',patname);
if ~isdir(figdir); mkdir(figdir); end;

binsize = 0.01;
win = [-0.2 0.5];

[cl_evix,~,~,ev_ix,EventInfo] = load_spks(rootdir_ad,Uname,sr,eventinx,'adevents',win); % cell array: 1x nr. of AD events
if isempty(cl_evix); fprintf('No %d ev\n',eventinx); return; end;

for j = 1:length(cl_evix)
    [frA{j},edgesA{j},spike_trainA{j}] = calc_ratemod(cl_evix{j},1:diff(win)*sr,sr,binsize);
end


%%

[PL0, spike_train_MX, edges, psth_FR_s,centers] = raster_psth(spike_trainA,edgesA,abs(win(1)),[],[],binsize,Uname,sr,isplot);

hold on;

if isplot
    yyaxis right
end
[PL1, LFPavg, lfp_time] = plot_lfpAD(rootdir_ad,patname,eventinx,win,chan,'thumb_low','b',isplot);
hold on

[PL2, GRavg, grid_time] = plot_lfpAD(rootdir_ad,patname,eventinx,win,gridchan,'EEG','r',isplot);

if isplot
    legend([PL0,PL1,PL2],{'SUA','LFP', 'Grid'})
    fig = gcf;
    fnm = fullfile(figdir,['AD_' num2str(eventinx) '_' Uname]);
    saveas(fig, [fnm '.fig']);
    saveas(fig, [fnm '.jpg']);
    close(fig);
end

end

function [PL, avgdat,t] = plot_lfpAD(rootdir_ad,patname,eventinx,win,adchan,rectype, color,isplot)

PL = [];
binsize = 0.01;
[normepoch, sr,h] = load_filt_epoch(rootdir_ad,patname,eventinx,1,win,[0.5 45],'Peaktimes','EEG');

t = linspace(win(1),win(2),size(normepoch,1));

if strcmp(rectype,'EEG')
    gr_chnames = clean_channames(h.chnames);
    adchan = find(ismember(gr_chnames,upper(adchan)));
end

avgdat = mean(normepoch(:,adchan),3);
if isplot
    PL = plot(t,avgdat,'Color',color,'LineStyle','-','LineWidth',2);
end


end
