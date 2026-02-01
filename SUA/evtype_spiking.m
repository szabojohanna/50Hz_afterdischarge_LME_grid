function evtype_spiking
load(fullfile(rootdir_ad,'SUA_Analysis','SUA_data_corrected'));
Unames = unit_table.Unit_name;


%%
Uname = 'Pt3_5_1'; % 'Pt3_3_3' ; 'Pt3_10_4'; 'Pt3_36_2'; 'Pt3_22_4';'Pt3_5_1'
evss =[ 133 138 153]; % [85 89 138 150 157 172]; [85 117 133 138 150 172]; [ 133 138 153];[153 138 133];[ 133 138 153]
AD_spiking_oneUnit_evtype(rootdir_ad,Uname,sr,evss)



%%
% evtype = {'Within','Close'};
for k = 1:2
    switch k;
        case 1; evtype = {'Within','Close'};
        case 2; evtype = {'Far'};
        case 3; evtype = {'Within','Close', 'MUA increase'};
        case 4; evtype = {'Within','Close', 'MUA decrease'};
        case 5; evtype = {'Far', 'MUA increase'};
        case 6; evtype = {'Far', 'MUA decrease'};
        case 7; evtype = {'Within','Close','noAD, out of SOZ'};
        case 8; evtype = {'Far','noAD, out of SOZ'};
    end
    %%
    for h = 1:length(Unames)
        Uname = Unames{h};
        %%
        isplot = true;
        stimevent_spiking_oneUnit_evtype(rootdir_ad,Uname,sr,evtype,isplot)
    end
  
end
evs2compare = {[85 89 138 150 157 172], s};

%%
% evtype = {'Within','Close','Far'};
% compar_tags =   {'AD','noAD, in SOZ'}; % {'noAD, out of SOZ'/ 'noAD, in SOZ'/ 'AD'

% Uname ='Pt3_36_2'; evtype = {[133 138 153],[132 137 92]}; compar_tags ={'AD','noAD, in SOZ'};
% Uname ='Pt3_36_2';evtype = {[132 137 92],[96   106    99]}; compar_tags ={'noAD, in SOZ','noAD, out of SOZ'};
% Uname ='Pt3_3_3'; evtype = {[85 89 138 150 157 172], [132 137 92 83 152 170]};  compar_tags ={'AD','noAD, in SOZ'};
% Uname ='Pt3_3_3'; evtype = {[85 89 138 150 157 172], [96 78  99 74 106 128]};  compar_tags ={'AD','noAD, out of SOZ'};
% Uname ='Pt3_5_1'; evtype = {[133 138 153], [132 137 92]};  compar_tags ={'AD','noAD, in SOZ'};
% Uname ='Pt3_5_1'; evtype = {[132 137 92],[96   106    99]};  compar_tags ={'AD','noAD, out of SOZ'};

for k = 1:11
    switch k;
        case 1; evtype = {'Within','Close'}; compar_tags = {'noAD, in SOZ','noAD, out of SOZ'};
        case 2; evtype = {'Far'}; compar_tags =   {'AD','noAD, in SOZ'};
        case 3; evtype = {'Within','Close', 'MUA increase'};compar_tags = {'AD','noAD, in SOZ'};
        case 4; evtype = {'Within','Close', 'MUA increase'}; compar_tags ={'noAD, in SOZ','noAD, out of SOZ'};
        case 5; evtype = {'Within','Close', 'MUA decrease'};compar_tags = {'AD','noAD, in SOZ'};
        case 6; evtype = {'Within','Close', 'MUA decrease'};compar_tags ={'noAD, in SOZ','noAD, out of SOZ'};
        case 7; evtype = {'Far', 'MUA increase'};compar_tags = {'AD','noAD, in SOZ'};
        case 8; evtype = {'Far', 'MUA increase'};compar_tags = {'noAD, in SOZ','noAD, out of SOZ'};
        case 9; evtype = {'Far', 'MUA decrease'};compar_tags = {'AD','noAD, in SOZ'};
        case 10; evtype = {'Far', 'MUA decrease'};compar_tags = {'noAD, in SOZ','noAD, out of SOZ'};
        case 11; evtype = {'Far'}; compar_tags =   {'noAD, in SOZ','noAD, out of SOZ'};
    end
for h = 1:length(Unames)
    Uname = Unames{h};
    %%
    stimevent_spiking_oneUnit_evtype_COMPARE(rootdir_ad,Uname,sr,evtype,compar_tags)
end
end
%



end




%--------------------------------------------------------------------------
function AD_spiking_oneUnit_evtype(rootdir_ad,Uname,sr,evtype)


[patname, chan, clnr] = unpack_Uname(Uname);

patdir = fullfile(rootdir_ad,'Patients',patname);
load(fullfile(patdir,'Events.mat'))


if chan<25; ti = 1; else; ti = 2; end;

if isvector(evtype);
    evs = sort(evtype);
else
    pr.evtype = evtype; pr.patname = patname; pr.period = 'stim';
    pr.rootdir_ad = rootdir_ad; pr.maxintens = false;
    evs= evs_from_tbl(pr,ti)';
    evs = sort(evs);
end
evnr = length(evs);
binsize = 0.01;

win = [-0.2 0.5];

%%
[spike_train_MX,edges, psth_FR_s,centers,LFPavg,lfp_time,GRavg,grid_time] = deal(cell(evnr,1));

for ei = 1:evnr
    eventinx = evs(ei);
    
    if isempty(Events(eventinx).AD_ch_order1); fprintf('No %d ev\n',eventinx); continue; end;
    gridchan = Events(eventinx).AD_ch_order1{1}{1};
    %%
    [spike_train_MX{ei},edges{ei}, psth_FR_s{ei},centers{ei}, LFPavg{ei},lfp_time{ei},GRavg{ei},grid_time{ei}] = AD_spiking_oneUnit_oneEvent(rootdir_ad,Uname,sr,eventinx,gridchan, false);
end

fnn = find(~cellfun(@isempty,edges));
edg = edges{fnn(1)};
cnt = centers{fnn(1)};
lfpt= lfp_time{fnn(1)};
gridt= grid_time{fnn(1)};


L = cellfun(@(x) size(x,1),spike_train_MX);
sp_tr = cat(1,spike_train_MX{fnn});
[rows, cols] = find(sp_tr);


figure;
subplot(2,1,1)
plot(edg(cols), rows, 'k.', 'MarkerSize', 10);
xL = xlim;
Lpos = 0;
for k = fnn'
    line(xL,[Lpos + L(k)+.1 Lpos + L(k)+.1]);
    Lpos = Lpos + L(k);
end
yL = ylim;
line([0 0],yL,'Color','k'); ylim([0 sum(L)])
setmyplot_balazs(gca); ylabel('AD events'); xlabel('Time (s)')

subplot(2,1,2)
for ei = fnn'
    pl= plot(cnt,psth_FR_s{ei},'Color',[.5 .5 .5]); hold on;
    set(pl,'Marker','none');
end
pl0 =plot(cnt,mean(cat(1,psth_FR_s{:}),1),'Color','k','LineWidth',3);
set(pl0,'Marker','none');
ylabel('Rate (Hz)');; xlabel('Time (s)')
yL = ylim;
line([0 0],yL,'Color','k')

yyaxis right
% subplot(3,1,3)
for ei = fnn'
    %     plot(lfpt,LFPavg{ei},'Color',[.5 .5 1]); hold on;
    pl = plot(gridt,GRavg{ei},'Color',[1 .5 .5],'LineStyle','-','LineJoin','chamfer');
    hold on;
    set(pl,'Marker','none');
end
% pl1 = plot(lfpt,mean(cat(2,LFPavg{:}),2),'Color','b','LineWidth',3);
pl2 = plot(gridt,mean(cat(2,GRavg{:}),2),'Color','r','LineWidth',3,'LineStyle','-');
set(pl2,'Marker','none');
ylabel('Norm. amplitude'); xlabel('Time (s)');
% legend([pl1, pl2],{'LFP','Grid'})
legend([pl0, pl2],{'SUA','Grid'},'Autoupdate','off')
yL = ylim;
line([0 0],yL,'Color','k')
setmyplot_balazs(gca);


Un = Uname;
Un(strfind(Uname,'_')) = '-';

nms = arrayfun(@(x) [num2str(x) ' '],evs,'UniformOutput',false);
nmss = [nms{:}];


title({[patname ', Unit: ' Un ', Evs: ' nmss], ['Ev nr = ' num2str(length(fnn)) ', AD nr = ' num2str(size(sp_tr, 1))]});

figdir = fullfile(rootdir_ad,'SUA_Analysis','Figures','AD_spiking',patname);
fnm = fullfile(figdir,[Uname '_' nmss]);
saveas(gcf,[fnm '.jpg'])
saveas(gcf,[fnm '.fig'])
saveas(gcf,[fnm '.pdf'])
close(gcf);
end

function stimevent_spiking_oneUnit_evtype_COMPARE(rootdir_ad,Uname,sr,evtype,compar_tags)


minL=.8; % minimum length of stim periods in sec
maxL_post = 1;% max length of post-stim period

binsize = 0.01;
bas_win = [-.3 0];

[patname, chan, clnr] = unpack_Uname(Uname);

if chan<25; ti = 1; else; ti = 2; end;
patdir = fullfile(rootdir_ad,'Patients',patname);
load(fullfile(patdir,'Events.mat'))
%%
if ischar(evtype{1});
    evs2compare = cell(2,1);
    for k = 1:2
        evtype2compare= {evtype{:},compar_tags{k}};
        evs2compare{k} = select_evtype(evtype2compare,rootdir_ad,patname,ti);
    end
else
    evs2compare = evtype;
end

%%
[evs2compare,ints] = intenscorrect_ADcompare(Events,evs2compare);

%%
allevs = cat(1,evs2compare{:})';
if isempty(allevs); fprintf('No events\n'); return; end;
[allevs,frE,edgesE,spike_trainE,basend,stimend,postend] = match_load_evs(rootdir_ad,Uname,allevs,minL,sr,binsize);

%%
[evs2comp_N, spike_train_MX, edges, psth_FR_s,centers,frE_gr] = deal(cell(2,1));
bas_ts = basend/sr;
stim_ts = stimend/sr;
post_ts = stim_ts+maxL_post;
for k = 1:2
    %%
    [evs2comp_N{k},evinx2comp] = intersect(allevs,evs2compare{k});
    [~,~,~,~, spike_train_MX{k}, edges{k}, psth_FR_s{k},centers{k}] = raster_psth(spike_trainE(evinx2comp),...
        edgesE(evinx2comp),bas_ts,stim_ts,post_ts,binsize,Uname,false,bas_win);
    frE_gr{k} = cat(1,frE{evinx2comp});
end


binx = 1:dsearchn(centers{1}',0);
sinx = binx(end)+1:dsearchn(centers{1}',stim_ts-bas_ts);
pinx = sinx(end)+1:dsearchn(centers{1}',post_ts-bas_ts);

% psth2 = cat(1,psth_FR_s{:});
% bavg = mean(psth2(binx),'all'); bstd = std(psth2(binx),[],'all');
% psth_norm = cellfun(@(x) (x-bavg)./bstd,psth_FR_s,'UniformOutput',0);

bavg = cellfun(@(x) repmat( mean(x(:,binx),'all') ,[size(x,1) size(x,2)]),frE_gr,'UniformOutput',0);
bstd = cellfun(@(x) repmat( std(x(:,binx),[],'all') ,[size(x,1) size(x,2)]),frE_gr,'UniformOutput',0);
frE_norm = arrayfun(@(x) (frE_gr{x}-bavg{x})./bstd{x},1:2,'UniformOutput',0);

FR_avg = cellfun(@(x) mean(x,1),frE_norm,'UniformOutput',0);
%%
fig = figure;
for k = 1:2
    
%     PL(k) = plot_psth(centers{k},psth_norm{k},bas_ts,stim_ts,post_ts);
    PL(k) = plot_psth(centers{k},FR_avg{k},bas_ts,stim_ts,post_ts);
    
    if strcmp(compar_tags{k},'AD');
        PL(k).Color = 'r';
    elseif strcmp(compar_tags{k},'noAD, in SOZ');
        PL(k).Color = 'k';
    elseif strcmp(compar_tags{k},'noAD, out of SOZ');
        PL(k).Color = 'b';
    end
    hold on;
end
ylabel('Norm. rate')
legend(PL,compar_tags,'Autoupdate','off')

%STAT

[pcond, ~, ~,~, ~,~] = std_stat({frE_norm{1}';frE_norm{2}'},'condstats','on','groupstats','off',...
    'mode','fieldtrip','fieldtripmethod','montecarlo','mcorrect','cluster',...
    'fieldtripalpha',0.05,'fieldtripclusterparam',{'clusterstatistic','maxsize'});


if ~isempty(pcond)
    sign_p = pcond{1};
    hold on; yL = ylim;
   
    draw_signifpatch(centers{k},sign_p,[0.8 0.8 0.8])

end

% bas_win = 1:binx;
% test_win1= binx:sinx;
% test_win2= sinx:pinx;

for k = 1:2
    stim_spk{k}=sum(spike_train_MX{k}(:,sinx),2)/length(sinx);
    post_spk{k}=sum(spike_train_MX{k}(:,pinx),2)/length(pinx);
end


[pval1, ~, ~] = ranksum(stim_spk{1},stim_spk{2});
[pval2, ~, ~] = ranksum(post_spk{1},post_spk{2});

yL = ylim;
if pval1<0.05; col = 'r'; else; col = 'k'; end;
text(0,yL(2),['p = ' num2str(pval1)],'Color',col);

if pval2<0.05; col = 'r'; else; col = 'k'; end;
text(stim_ts-bas_ts,yL(2),['p = ' num2str(pval2)],'Color',col);


% Save fig

if iscell(evtype{1})
    et = cellfun(@(x) [x ' '],evtype,'UniformOutput',false);
    et2= cellfun(@(x) [x '_'],evtype,'UniformOutput',false);
    stit='';
else
    et = cellfun(@(x) [x ' '],arrayfun(@num2str,[evtype{:}],'UniformOutput',false),'UniformOutput',false);
    stit0= cellfun(@(x) [x '_'],arrayfun(@num2str,[evtype{:}],'UniformOutput',false),'UniformOutput',false);
    stit= [stit0{:}];
    et2 = {'Selected'};
end
evtit = [et{:}];
evftit = [et2{:}];

title({evtit, [compar_tags{1} ': n = ' num2str(length(evs2comp_N{1}))],...
   [compar_tags{2} ': n = ' num2str(length(evs2comp_N{2}))] });


ct= cellfun(@(x) [x '_VS_'],compar_tags,'UniformOutput',false);
ctit = [ct{:}];


figdir = fullfile(rootdir_ad,'SUA_Analysis','Figures','intrastim_spiking',...
    [patname '_evtypes'],['Compare_' evftit '_' ctit]);
if ~isdir(figdir); mkdir(figdir); end;

fnm = fullfile(figdir,[Uname stit]);
saveas(fig,[fnm '.jpg'])
saveas(fig,[fnm '.fig'])
close(fig)

if exist(fullfile(figdir,'Stat.mat'))==2
    load(fullfile(figdir,'Stat.mat'));
else
    Stat = struct;
end
Stat.(Uname).Stim = pval1;
Stat.(Uname).Post = pval2;
save(fullfile(figdir,'Stat.mat'),'Stat');
end

%--------------------------------------------------------------------------
function [psth_FR_s,centers] = stimevent_spiking_oneUnit_evtype(rootdir_ad,Uname,sr,evtype,isplot)

minL=.8; % minimum length of stim periods in sec
maxL_post = 1;% max length of post-stim period
[patname, chan, clnr] = unpack_Uname(Uname);

if chan<25; ti = 1; else; ti = 2; end;

evs = select_evtype(evtype,rootdir_ad,patname,ti);


binsize = 0.01;

[evs,frE,edgesE,spike_trainE,basend,stimend,postend] = match_load_evs(rootdir_ad,Uname,evs,minL,sr,binsize);
evnr = length(evs);

%%
bas_win = [-.3 0];
[fig,pval,type,~, ~, ~, psth_FR_s,centers] = raster_psth(spike_trainE,edgesE,basend/sr,stimend/sr,stimend/sr+maxL_post,binsize,Uname,isplot, bas_win);


% STAT 2

if iscell(evtype)
    et2= cellfun(@(x) [x '_'],evtype,'UniformOutput',false);
    et = cellfun(@(x) [x ' '],evtype,'UniformOutput',false);

else
    et2 = arrayfun(@(x) [num2str(x) '_'] ,evtype, 'UniformOutput',false);
    et = arrayfun(@(x) [num2str(x) ' '] ,evtype, 'UniformOutput',false);
    
end
evftit = [et2{:}];
evtit = [et{:}];


figdir = fullfile(rootdir_ad,'SUA_Analysis','Figures','intrastim_spiking',...
    [patname '_evtypes'],evftit);
if ~isdir(figdir); mkdir(figdir); end;

if isplot
title([evtit ', n = ' num2str(evnr) ' events']);


fnm = fullfile(figdir,Uname);
saveas(fig,[fnm '.jpg'])
saveas(fig,[fnm '.fig'])
close(fig)
end
if exist(fullfile(figdir,'Stat.mat'))==2
    load(fullfile(figdir,'Stat.mat'));
else
    Stat = struct;
end
Stat.(Uname).pval= pval;
Stat.(Uname).type = type;
save(fullfile(figdir,'Stat.mat'),'Stat');
end

%--------------------------------------------------------------------------
function [fig,pval,type,PL, spike_train_MX, edges, psth_FR_s,centers] = raster_psth(spike_trainE,edgesE,bas_ts,stim_ts,post_ts,binsize,Uname,isplot, bas_win)


[fig,pval1,pval2,PL, spike_train_MX, edges, psth_FR_s,centers] = deal([]);
evnr = length(spike_trainE);
mX = 1; miL = length(spike_trainE{mX});

spike_train_MX = cat(1,spike_trainE{:});
edges = edgesE{1}-bas_ts;
centers = edges(1:end-1) + diff(edges)/2;


psth = sum(spike_train_MX,1);
psth_FR = psth./evnr;
psth_FR = psth_FR./binsize;
psth_FR_s = smoothdata(psth_FR, 'gaussian', 10);

% STAT
binx = dsearchn(centers',bas_win');
baswin_ix = binx(1):binx(2);

stimL = stim_ts-bas_ts;
testL = stimL/2;

pval = nan(4,1);
type = nan(4,1);
for t = 1:4
    switch t ;
        case 1; test_win = [0 testL];
        case 2; test_win = [testL 2*testL];
        case 3; test_win = [stimL stimL+testL];
        case 4; test_win = [stimL+testL stimL+2*testL];
    end

tinx = dsearchn(centers',test_win');
test_ix= tinx(1):tinx(2);


bas_spk=sum(spike_train_MX(:,baswin_ix),2)/length(baswin_ix);
stim_spk=sum(spike_train_MX(:,test_ix),2)/length(test_ix);

% Bsl vs. Stim.
[pval(t), ~, ~] = ranksum(bas_spk,stim_spk);

if pval(t)<0.05
    if mean(bas_spk)>mean(stim_spk)
        type(t) = -1;
    elseif mean(bas_spk)<mean(stim_spk)
        type(t) = 1;
    end
else
    type(t) = 0;
end
end
%% FIGURE
if isplot
    fig = figure;
    set(fig,'Visible','off');
    subplot(2,1,1)
    hold on;
    
    for k = 1:evnr
        if ~isnan(spike_trainE{k});
            spike_times = (edgesE{k}(spike_trainE{k}==1))-bas_ts;
            y = ones(size(spike_times)) * k;
            plot(spike_times, y, '.', 'Color', 'k', 'MarkerSize', 10); hold on;
        end
    end
    
    yL = [.5 evnr+.5]; ylim(yL)
    if ~isempty(bas_ts)
        %         line([bas_ts bas_ts],yL,'Color','k')
        line([0 0],yL,'Color','k')
    end
    if ~isempty(stim_ts)
        line([stim_ts stim_ts]-bas_ts,yL,'Color','k')
    end
    if ~isempty(post_ts)
        xlim([-bas_ts post_ts-bas_ts]);
        
    end
    % xlabel('Time (s)');
    ylabel('Events');
    title({Uname,['Spike count, bin size:' num2str(binsize) 'dp']});
    
    subplot(2,1,2)
    
    PL = plot_psth(centers,psth_FR_s,bas_ts,stim_ts,post_ts);
    
    yL = ylim;
    xp = 0;
    for t = 1:4
        if pval(t)<0.05; col = 'r'; else; col = 'k'; end;
        text(xp,yL(2)*(1-t/10)+.1,['p = ' num2str(pval(t))],'Color',col);
        xp = xp+ testL;
    end
    
    
end

end

function evs = select_evtype(evtype,rootdir_ad,patname,ti)

if ~iscell(evtype);
    evs = sort(evtype);
else
    pr.evtype = evtype; pr.patname = patname; pr.period = 'stim';
    pr.rootdir_ad = rootdir_ad; pr.maxintens = false;
    evs= evs_from_tbl(pr,ti)';
    evs = sort(evs);
end
end

function [evs,frE,edgesE,spike_trainE,basend,stimend,postend] = match_load_evs(rootdir_ad,Uname,evs,minL,sr,binsize);

[patname, chan, clnr] = unpack_Uname(Uname);

% Check if length of stim. is too short
if ~isempty(minL)
    patdir = fullfile(rootdir_ad,'Patients',patname);
    load(fullfile(patdir,'Events.mat'));
    evinx0 = evs;
    stimLs = round(([Events(evs).Stim_end] -[Events(evs).Stim_start])*10)/10;
    goodinx = stimLs>=minL;
    evs = evinx0(goodinx);
    fprintf('excluded: Ev%d \n', evinx0(~goodinx)');
end
evnr = length(evs);

[~,~,mixiB,~]=matchblocklims(rootdir_ad,patname,evs,'baseline', 'sua');
[~,~,mixiS,~]=matchblocklims(rootdir_ad,patname,evs,'stim', 'sua');
[~,~,mixiP,~]=matchblocklims(rootdir_ad,patname,evs,'post', 'sua');

frE = cell(evnr,1);
spike_trainE = cell(evnr,1);
edgesE = cell(evnr,1);
for k = 1:evnr
    eventinx = evs(k);
    
    [~,~,~,~,EventInfo,basend,stimend,postend, ...
        cl_evixB,ev_ixB, cl_evixS,ev_ixS, cl_evixP,ev_ixP] = load_spks(rootdir_ad,...
        Uname,sr,eventinx,'intrastim',{mixiB{k} mixiS{k} mixiP{k}});
    for j = 1:3
        switch j
            case 1; ev_ix = ev_ixB; cl_evix = cl_evixB;
            case 2; ev_ix = ev_ixS; cl_evix = cl_evixS;
            case 3; ev_ix = ev_ixP; cl_evix = cl_evixP;
        end
    [frE0{j},edgesE0{j},spike_trainE0{j}] = calc_ratemod(cl_evix,ev_ix,sr,binsize);

    end
    frE{k} = cat(2,frE0{:});
    
%     edgesE{k}= cat(2,edgesE0{:});
edgesE{k} = 0:binsize:binsize*length(frE{k});
    spike_trainE{k}= cat(2,spike_trainE0{:});
end


end


%--------------------------------------------------------------------------
function PL = plot_psth(centers,psth_FR_s,bas_ts,stim_ts,post_ts)

PL = plot(centers,psth_FR_s,'Color','k','LineWidth',3); hold on;
%     plot(psth_FR_s,'Color','k','LineWidth',3); hold on;
yL = ylim;
if ~isempty(bas_ts)
    %         line([bas_ts bas_ts],yL,'Color','k')
    line([0 0],yL,'Color','k')
end
if ~isempty(stim_ts)
    line([stim_ts stim_ts]-bas_ts,yL,'Color','k')
end

if ~isempty(post_ts)
    xlim([-bas_ts post_ts-bas_ts]);
end
ylabel('Rate (Hz)');
xlabel('Time (s)')
end