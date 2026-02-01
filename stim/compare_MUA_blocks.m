function compare_MUA_blocks(rootdir_ad,patname,ti,compar_tags,varargin)

prs = inputParser;
addRequired(prs,'rootdir_ad',@isfolder);
addRequired(prs,'patname',@ischar);
addRequired(prs,'ti',@isnumeric);
addRequired(prs,'compar_tags',@iscell);
addParameter(prs,'intens_correct',false,@islogical);
addParameter(prs,'maxintens',true,@islogical);
addParameter(prs,'cLim',[-1.5 1.5],@isvector);
addParameter(prs,'cLim_diff',[-2 2],@isvector);
addParameter(prs,'blockavg',true,@islogical);
addParameter(prs,'smoothwith','smoothdata',@ischar); % 'none' | 'smoothdata' | 'lowpass'
addParameter(prs,'smparam',10,@isnumeric); % NaN if smoothdata is 'none'; 5-20 if smoothdata is 'smoothdata'; ~20 if smoothdata is 'lowpass';
addParameter(prs,'stattype','fieldtrip',@ischar); % 'fieldtrip' | 'eeglab' | 'my_permutation' | 'bayesian';
addParameter(prs,'mcorr','cluster',@ischar); % 'none' | 'cluster';
addParameter(prs,'itnum',500,@isnumeric);
addParameter(prs,'minL',.8,@isnumeric);
addParameter(prs,'evtype','all',@(x) ischar(x)|iscell(x));
addParameter(prs,'timewin',[0 .9],@isvector); % relevant if evtype is 'mua_up'/'mua_down'; % [0 .45] | [0.23 .67] | [.45 .9] |  sec relative to PERIOD start
addParameter(prs,'upchans',[],@(x) isnumeric(x)|isvector(x));  % relevant if evtype is 'mua_up'/'mua_down'; if empty, all channels
addParameter(prs,'period','stim',@ischar); % relevant if evtype is 'mua_up'/'mua_down'; bas | stim | post
addParameter(prs,'bnorm','Indiv',@ischar); % Normalization method: 'Zscore' (zscore for each channel, whole window) | 'Common' (comm. baseline across each events) | 'Indiv' (indiv. baseline for each event)
parse(prs,rootdir_ad,patname,ti,compar_tags,varargin{:});
pr = prs.Results;


muablockdir = fullfile(rootdir_ad,'Figures','Micro','MUA','MUA_blocks');

load(fullfile(rootdir_ad,'Patients',patname,'Events.mat'));


load(fullfile(rootdir_ad,'channels_by_layers.mat'));
islaypat  = ismember(patname,channels_by_layers.Properties.VariableNames);


if ti==1; chans = 1:23; else; chans = 25:47; end;


% Load data
[all_blocks,z_groups_cc, z_groups_avg, close_evinx2compare, ts] = load_evgroups(ti, chans,muablockdir,true,pr);

if isempty(all_blocks); return; end;
% Channel-time map
[dfmap, pmap] = ch_t_map(all_blocks,z_groups_cc,z_groups_avg,ti,chans,close_evinx2compare,muablockdir,ts,pr);


%     avgavg_boxplot(z_groups_cc,dfmap,ts)
    
if ~any(cellfun(@isempty,z_groups_cc))
    % Avg across all channels
    all_chavg_plot(z_groups_cc,z_groups_avg,ti,chans,muablockdir,ts,pr)
    
    [~,z_groups_ccL, z_groups_avgL, ~, ts] = load_evgroups(ti, chans,muablockdir,false,pr);
    
    
    if islaypat
        
        
        % Avg across infra/supragran. layers
        infra_supra_chavg_plot(z_groups_ccL,z_groups_avgL,ti,muablockdir,ts,pr)
        
        % Avg. across indiv. layers
%         layer_by_layer_plot(z_groups_ccL,z_groups_avgL,ti,muablockdir,ts,pr)
    else
        % Avg across upper/ lower half of thumb rec.
        half_chavg_plot(z_groups_ccL,z_groups_avgL,ti,muablockdir,ts,pr);
    end
    
end
%%


end

function avgavg_boxplot(z_groups_cc,dfmap,ts)

basend = ts.BL;
stimend = ts.BL + ts.SL;
stimhalf = ts.BL + floor(ts.SL/2);
postend = ts.BL + ts.SL + ts.PL;
% 
% load(fullfile(rootdir_ad, 'good_channels.mat'))
% chs = good_channels.(patname){ti}; if ti==2; chs = chs-24; end;

figure(1);
subplot(2,1,1);
imagesc(dfmap); hold on; colormap  jet; caxis([-1.5 1.5]); contour(pmap,'Color','w');


chavg = mean( cat(4,z_groups_cc{:}) , [2 3 4]);

[mx,chmx] = findpeaks(chavg,'NPeaks',1,'SortStr','descend');



for chmx = 1:length(chavg)
per = cellfun(@(x) squeeze( mean(x(chmx,basend+1:stimhalf,:),[1 2]) ), z_groups_cc,'UniformOutput',false);
bp = cat(2, per{:});


[pval(chmx), ~,~] = ranksum(bp(:,1), bp(:,2));
if pval(chmx)<0.05; col = 'r'; else; col = 'k'; end;


end

figure(2); ph=plot(pval'); xlabel('Channels'); ylabel('P values');
hold on; xL = xlim; line(xL,[0.05 0.05]);
hold on; yyaxis right; ch=plot(chavg); ylabel('Channel avg');
legend([ph,ch],{'Diff. in first stim. half', 'Channel avg'}); 

disp('Give chs and tx params')
keyboard;

close(figure(2));

hold on; 
rectangle('Position',[tx(1), chs(1), diff(tx([1 end])),diff(chs([1 end]))],'LineWidth',2)
% chs = chmx-1:chmx+1;

% figure; plot(chavg,1:length(chavg)); hold on; scatter(mx,chmx)

% tx = basend+1:stimhalf;
per = cellfun(@(x) squeeze( mean(x(chs,tx,:),[1 2]) ), z_groups_cc,'UniformOutput',false);
bp = cat(2, per{:});


[pval, ~,~] = ranksum(bp(:,1), bp(:,2));
if pval<0.05; col = 'r'; else; col = 'k'; end;


subplot(2,1,2);
boxplot(bp)
yL = ylim;
text(0.5,yL(2)*.9,['p = ' num2str(pval)],'Color',col);
set_my_boxplot(gca); set





end

%
% function [cc, cc_time] = load_concat_blocks(patname,eventinx,isAD,chans,w2cc,muablockdir)
%
%
% load(fullfile(muablockdir,['MUA_stimblocks_' patname '_Ev' num2str(eventinx) '_AD' num2str(isAD) '.mat']));
% block_MUA = MUA_stimblocks.stimblock_MUA;
% basblock_MUA = MUA_stimblocks.basblock_MUA;
% postblock_MUA = MUA_stimblocks.postblock_MUA;
%
% Bblocks_lims = MUA_stimblocks.params.Bblocks_lims;
% blocks_lims = MUA_stimblocks.params.blocks_lims;
% Pblocks_lims = MUA_stimblocks.params.Pblocks_lims;
%
% block_dp = size(block_MUA,1);
% cc = []; cc_time = [];
% if contains('baseline',w2cc)
%     for j = 1:size(basblock_MUA,3)
%         cc = cat(1,cc,basblock_MUA(:,chans,j));
%     end
%      for k = 1:size(basblock_MUA,3)
%          tim = linspace(Bblocks_lims(k,1),Bblocks_lims(k,2),block_dp);
%         cc_time = cat(1,cc_time,tim');
%     end
% end
% if contains('stim',w2cc)
%     for j = 1:size(block_MUA,3)
%         cc = cat(1,cc,block_MUA(:,chans,j));
%     end
%      for k = 1:size(block_MUA,3)
%          tim = linspace(blocks_lims(k,1),blocks_lims(k,2),block_dp);
%         cc_time = cat(1,cc_time,tim');
%     end
% end
% if contains('post',w2cc)
%     for j = 1:size(postblock_MUA,3)
%         cc = cat(1,cc,postblock_MUA(:,chans,j));
%     end
%     for k = 1:size(postblock_MUA,3)
%         tim = linspace(Pblocks_lims(k,1),Pblocks_lims(k,2),block_dp);
%         cc_time = cat(1,cc_time,tim');
%     end
% end
%
% end
%
% function [allblocks,all_lims, matchblock,matchblock_lims, matchblock_cc,matchblock_cc_time ] = load_matching_blocks(patname,ad_events,isAD_vec,w2cc,muablockdir)
%
% evnr = length(ad_events);
% k = 1; [allblocks, all_lims] = deal(cell(1,evnr));
% for eventinx = ad_events
%
%     load(fullfile(muablockdir,['MUA_stimblocks_' patname '_Ev' num2str(eventinx) '_AD' num2str(isAD_vec(k)) '.mat']));
%
%     switch w2cc;
%         case 'baseline'
%             block = MUA_stimblocks.basblock_MUA;
%             lims =  MUA_stimblocks.params.Bblocks_lims;
%         case 'stim'
%             block = MUA_stimblocks.stimblock_MUA;
%             lims = MUA_stimblocks.params.blocks_lims;
%         case 'post'
%             block = MUA_stimblocks.postblock_MUA;
%             lims =  MUA_stimblocks.params.Pblocks_lims;
%     end
%     allblocks{k} = block;
%     all_lims{k} = lims; % limits of each block (s) - relative to start of original thumb file
%     k = k+1;
% end
% block_dp = size(block,1);
%
%
% [mL, mLix] = min(cellfun(@(x) size(x,3),allblocks)); % event with least numbert of blocks (shortest)
%
% all_lims_c = cellfun(@(x) x(:,1) - x(1,1),all_lims,'UniformOutput',0); % start of each limit relative to period start (s)
% ref_lims = all_lims_c{mLix}; % time vector of shortest event
%
% % Find matching blocks (calculate differences in normalized block-starts
% % -> block with the closest starting value is the corresponding block in time)
% % ("extra" blocks might be present even in the middle of the block sequence due to faulty detection of stim./ artefact peaks)
% [mix, mixi] = deal(cell(1,evnr));
% for j = 1:length(ref_lims) % loop over blocks
%     for m = 1:evnr % loop over events
%         [mix{m}(j), mixi{m}(j)] = min(abs(ref_lims(j) - all_lims_c{m}));
%     end
% end
%
% matchblock = arrayfun(@(x) allblocks{x}(:,:,mixi{x})  ,1:evnr,'UniformOutput',0);
% matchblock_lims = arrayfun(@(x) all_lims{x}(mixi{x},:)  ,1:evnr,'UniformOutput',0);
% matchblock_cc = cellfun(@(x) reshape( permute(x,[1 3 2]) ,[size(x,1)*size(x,3) size(x,2)]) , matchblock,'UniformOutput',0);
%
% matchblock_cc_time = cell(1,evnr);
% for r = 1:evnr
%     matchblock_cc_time{r} = [];
%     for q =1:mL
%         blL = linspace(matchblock_lims{r}(q,1) , matchblock_lims{r}(q,2), block_dp);
%         matchblock_cc_time{r} = cat(1,matchblock_cc_time{r}, blL');
%     end
% end
%
%
% end
%
% function [ad_blocks, ad_blocks_time, ad_smoth] = ds_smooth(patname,...
%     ad_events,isAD_vec,chans,w2cc,smoothwith,smparam,blockavg,muablockdir)
%
% %%
% [~,~, ad_blocks_3d,ad_blocks_3d_time, ad_blocks_concat,ad_blocks_concat_time ] = load_matching_blocks(patname,ad_events,isAD_vec,w2cc,muablockdir);
%
% if blockavg
%     ad_blocks = cellfun(@(x) permute( mean(x,1) , [3 2 1]), ad_blocks_3d,'UniformOutput' ,false);
%     ad_blocks_time = cellfun(@(x) x(:,1), ad_blocks_3d_time,'UniformOutput' ,false );
% else
%
%     ad_blocks = ad_blocks_concat;
%     ad_blocks_time = ad_blocks_concat_time;
% end
%
% %%
% ad_blocks = cellfun(@(x) x(:,chans) ,ad_blocks,'UniformOutput' ,false);
% if strcmp(smoothwith,'smoothdata')
%     ad_smoth = cellfun(@(x) smoothdata(x,1,'movmedian',smparam) ,ad_blocks, 'UniformOutput' ,false);
% elseif strcmp(smoothwith,'lowpass')
%     srate = 20000;
%     [b,a]=butter(2,smparam/(srate/2),'low');
%     ad_smoth = cellfun(@(x) filtfilt(b,a,x) ,ad_blocks, 'UniformOutput' ,false);
% elseif strcmp(smoothwith,'none')
%     ad_smoth = ad_blocks;
% end
%
%
%
% end



function [dfmap,pmap] =ch_t_map(all_blocks,z_groups_cc,z_groups_avg,ti,chans,close_evinx2compare,muablockdir,ts,pr)



for p = 1:2
   evlab = arrayfun(@(x) ['Ev' num2str(x) '  '],close_evinx2compare{p},'UniformOutput',0);
    tits{p} = [pr.compar_tags{p} ', stim sites: ' evlab{:}];
end

pr.basend = ts.BL;
pr.stimend = ts.BL + ts.SL;
pr.postend = ts.BL + ts.SL + ts.PL;
pr.ti = ti;
xti = [1 ts.BL round(ts.BL+ts.SL/3) round(ts.BL+ts.SL*2/3) pr.stimend round(pr.stimend+ts.PL/2) pr.postend];
evtime2 = round(ts.evtime*10)/10;
xtilabs = arrayfun(@num2str, evtime2(xti), 'UniformOutput', false);

suptit = {[pr.patname 'TH' num2str(ti)],['MUA blocks smoothed (' pr.smoothwith ', param:' num2str(pr.smparam) ')']};


[figdir, fnm, fnm2] = save_fnms(muablockdir,ti,pr);

[dfmap, pmap] = plot_stim_statmaps(z_groups_cc, tits,suptit,xti, xtilabs, fnm, pr);


%% Channel-time map figure


StatMAT.zscored_groups = z_groups_cc;
StatMAT.all_blocks = all_blocks;
StatMAT.DIFFmap = dfmap;
StatMAT.Pmap = pmap;
StatMAT.events2compare = close_evinx2compare;
StatMAT.basline_end = pr.basend;
StatMAT.stim_end = pr.stimend;
StatMAT.data_smooth = [pr.smoothwith ', param: ' pr.smparam];
StatMAT.Bblock_starts = ts.BB;
StatMAT.Sblock_starts = ts.SB;
StatMAT.Pblock_starts = ts.PB;

save(fnm2, 'StatMAT')

end

function all_chavg_plot(z_groups_cc,z_groups_avg,ti,chans,muablockdir,ts,pr)


if any(cellfun(@(x) size(x,3), z_groups_cc)<2); return; end;

basend = ts.BL;
stimend = ts.BL + ts.SL;
postend = ts.BL + ts.SL + ts.PL;


figtit = {[pr.patname 'TH' num2str(ti)],['MUA blocks smoothed (' pr.smoothwith ', param:' num2str(pr.smparam) ')']};
cols = {'r','k'};

%% AVG all channels
fig = figure;
t = 1:ts.BL+ts.SL+ts.PL;
chs = 1:size(z_groups_avg{1},1);

for p = 1:2
    avg = mean(z_groups_avg{p}(chs,:),1);
    if isempty(avg); continue; end;
    sd = std( mean(z_groups_cc{p}(chs,:,:),1) ,[],3) / size(z_groups_cc{p},3);
    sdup = avg + sd; sddo = avg-sd;
    pL(p) = plot(t, avg,'Color',cols{p}); hold on;
    patch('XData',[t fliplr(t)], 'YData',[sddo fliplr(sdup)],...
        'FaceAlpha',0.2,'FaceColor',cols{p},'EdgeColor','none');
end
yL = ylim;
line([basend basend],yL,'Color','k','LineStyle','--')
line([stimend stimend],yL,'Color','k','LineStyle','--'); ylim(yL);


if ~any(cellfun(@isempty, z_groups_cc))
    for j = 1:3
        switch j; case 1; tinx = 1:ts.BL; case 2; tinx = ts.BL+1:ts.BL+ts.SL; case 3; tinx = ts.BL+ts.SL+1:ts.BL+ts.SL+ts.PL; end;
        % tinx = 1:ts.BL+ts.SL+ts.PL;
        z_groups_bl = cellfun(@(x) mean(x(chs,tinx,:),1) , z_groups_cc,'UniformOutput',0);
        
        if strcmp(pr.stattype,'bayesian')
            [~,weak, moderate, strong] = my_stats(z_groups_bl,pr);
            
            draw_signifpatch(tinx,weak',[.2 .5 .5]); hold on;
            draw_signifpatch(tinx,moderate',[.5 .5 .5]); hold on;
            draw_signifpatch(tinx,strong',[.9 .5 .5]); hold on;
        else
            [~,sign_p] = my_stats(z_groups_bl,pr);
            draw_signifpatch(tinx,sign_p,[.7 .5 .5])
        end
        
    end
end
xti = [1 ts.BL round(ts.BL+ts.SL/3) round(ts.BL+ts.SL*2/3) stimend round(stimend+ts.PL/2) postend];
evtime2 = round(ts.evtime*10)/10;
xtilabs = arrayfun(@num2str, evtime2(xti), 'UniformOutput', false);

xticks(xti); xticklabels(xtilabs); xlim(xti([1 end]))
xlabel('Time rel. to period start (s)'); ylabel('Norm. MUA');
title(figtit)
legend(pL,pr.compar_tags)

[~, fnm, ~] = save_fnms(muablockdir,ti,pr);

saveas(fig,[fnm '_ChAVG.jpg']);
saveas(fig,[fnm '_ChAVG.fig']);
close(fig);
end

function infra_supra_chavg_plot(z_groups_cc,z_groups_avg,ti,muablockdir,ts,pr)


if any(cellfun(@(x) size(x,3), z_groups_cc)<2); return; end;

basend = ts.BL;
stimend = ts.BL + ts.SL;
postend = ts.BL + ts.SL + ts.PL;


load(fullfile(pr.rootdir_ad,'channels_by_layers.mat'));
layers = channels_by_layers.(pr.patname)(:,ti);



cols = {'r','k'};


for k = 1:3
    %% AVG channels
    fig = figure;
    t = 1:ts.BL+ts.SL+ts.PL;
    switch k;
        case 1
            chs = [layers{1:3}];
            Ltit = 'Supragran';
        case 2
            chs = [layers{4}];
            Ltit = 'Gran';
        case 3
            chs = [layers{5:6}];
            Ltit = 'Infragran';
    end
    
    
    figtit = {[pr.patname 'TH' num2str(ti) ', ' Ltit ' layers: ' num2str(chs(1)) '-' num2str(chs(end)) ' channels']...
        ,['MUA blocks smoothed (' pr.smoothwith ', param:' num2str(pr.smparam) ')']};
    
    for p = 1:2
        avg = mean(z_groups_avg{p}(chs,:),1);
        sd = std( mean(z_groups_cc{p}(chs,:,:),1) ,[],3) / size(z_groups_cc{p},3);
        sdup = avg + sd; sddo = avg-sd;
        pL(p) = plot(t, avg,'Color',cols{p}); hold on;
        patch('XData',[t fliplr(t)], 'YData',[sddo fliplr(sdup)],...
            'FaceAlpha',0.2,'FaceColor',cols{p},'EdgeColor','none');
    end
    yL = ylim;
    ylim(yL);
    line([basend basend],yL,'Color','k','LineStyle','--')
    line([stimend stimend],yL,'Color','k','LineStyle','--');
    for j = 1:3
        switch j; case 1; tinx = 1:ts.BL; case 2; tinx = ts.BL+1:ts.BL+ts.SL; case 3; tinx = ts.BL+ts.SL+1:ts.BL+ts.SL+ts.PL; end;
        % tinx = 1:ts.BL+ts.SL+ts.PL;
        z_groups_bl = cellfun(@(x) mean(x(chs,tinx,:),1) , z_groups_cc,'UniformOutput',0);
        
        if strcmp(pr.stattype,'bayesian')
            [~,weak, moderate, strong] = my_stats(z_groups_bl,pr);
            
            draw_signifpatch(tinx,weak',[.2 .5 .5]); hold on;
            draw_signifpatch(tinx,moderate',[.5 .5 .5]); hold on;
            draw_signifpatch(tinx,strong',[.9 .5 .5]); hold on;
        else
            [~,sign_p] = my_stats(z_groups_bl,pr);
            draw_signifpatch(tinx,sign_p,[.7 .5 .5])
        end
    end
    
    xti = [1 ts.BL round(ts.BL+ts.SL/3) round(ts.BL+ts.SL*2/3) stimend round(stimend+ts.PL/2) postend];
    evtime2 = round(ts.evtime*10)/10;
    xtilabs = arrayfun(@num2str, evtime2(xti), 'UniformOutput', false);
    
    xticks(xti); xticklabels(xtilabs); xlim(xti([1 end]))
    xlabel('Time rel. to period start (s)'); ylabel('Norm. MUA');
    title(figtit)
    legend(pL,pr.compar_tags)
    
    [~, fnm, ~] = save_fnms(muablockdir,ti,pr);
    
    saveas(fig,[fnm '_' Ltit '.jpg']);
    saveas(fig,[fnm '_' Ltit '.fig']);
    close(fig);
    
end


fig = figure;
n = 1;
for k = 1:3
    
    switch k;
        case 1
            chs = [layers{1:3}];
            Ltit = 'Supragran';
        case 2
            chs = [layers{4}];
            Ltit = 'Gran';
        case 3
            chs = [layers{5:6}];
            Ltit = 'Infragran';
    end
    for j = 1:3
        switch j;
            case 1; tinx = 1:ts.BL; Ptit = 'Baseline';
            case 2; tinx = ts.BL+1:ts.BL+ts.SL; Ptit = 'Stim.';
            case 3; tinx = ts.BL+ts.SL+1:ts.BL+ts.SL+ts.PL;  Ptit = 'Poststim.';
        end;
        
        z_bp = cellfun(@(x) squeeze(mean( mean(x(chs,tinx,:),1),2)), z_groups_cc,'UniformOutput',0);
        
        subplot(3,3,n)
        mL = max( cellfun(@length,z_bp) );
bp = cellfun(@(x) [x; nan(mL-length(x),1)] ,z_bp, 'UniformOutput',0);
        boxplot([bp{:}]); set_my_boxplot(gca);
        [pval, ~, ~] = ranksum(z_bp{1},z_bp{2});
        if pval<0.05; col = 'r'; else; col = 'k'; end;
        yL = ylim;
        text(0.5,yL(2)*.9,['p = ' num2str(pval)], 'Color',col)
        title([Ltit ', ' Ptit]);
        n = n+1;
    end
end

[~, fnm, ~] = save_fnms(muablockdir,ti,pr);

saveas(fig,[fnm '_SI_boxplots.jpg']);
saveas(fig,[fnm '_SI_boxplots.fig']);
close(fig);
end




function half_chavg_plot(z_groups_cc,z_groups_avg,ti,muablockdir,ts,pr)

if any(cellfun(@(x) size(x,3), z_groups_cc)<2); return; end;

basend = ts.BL;
stimend = ts.BL + ts.SL;
postend = ts.BL + ts.SL + ts.PL;


cols = {'r','k'};
yL = [-1 2];

for k = 1:2
    %% AVG channels
    fig = figure;
    t = 1:ts.BL+ts.SL+ts.PL;
    switch k;
        case 1
            chs = 1:8;
            Ltit = 'Upper half';
        case 2
            chs = 15:23;
            Ltit = 'Lower half';
    end
    
    
    figtit = {[pr.patname 'TH' num2str(ti) ', ' Ltit ' layers: ' num2str(chs(1)) '-' num2str(chs(end)) ' channels']...
        ,['MUA blocks smoothed (' pr.smoothwith ', param:' num2str(pr.smparam) ')']};
    
    for p = 1:2
        avg = mean(z_groups_avg{p}(chs,:),1);
        sd = std( mean(z_groups_cc{p}(chs,:,:),1) ,[],3) / size(z_groups_cc{p},3);
        sdup = avg + sd; sddo = avg-sd;
        pL(p) = plot(t, avg,'Color',cols{p}); hold on;
        patch('XData',[t fliplr(t)], 'YData',[sddo fliplr(sdup)],...
            'FaceAlpha',0.2,'FaceColor',cols{p},'EdgeColor','none');
    end
    ylim(yL);
    line([basend basend],yL,'Color','k','LineStyle','--')
    line([stimend stimend],yL,'Color','k','LineStyle','--');
    for j = 1:3
        switch j; case 1; tinx = 1:ts.BL; case 2; tinx = ts.BL+1:ts.BL+ts.SL; case 3; tinx = ts.BL+ts.SL+1:ts.BL+ts.SL+ts.PL; end;
        % tinx = 1:ts.BL+ts.SL+ts.PL;
        z_groups_bl = cellfun(@(x) mean(x(chs,tinx,:),1) , z_groups_cc,'UniformOutput',0);
        
        if strcmp(pr.stattype,'bayesian')
            [~,weak, moderate, strong] = my_stats(z_groups_bl,pr);
            
            draw_signifpatch(tinx,weak',[.2 .5 .5]); hold on;
            draw_signifpatch(tinx,moderate',[.5 .5 .5]); hold on;
            draw_signifpatch(tinx,strong',[.9 .5 .5]); hold on;
        else
            [~,sign_p] = my_stats(z_groups_bl,pr);
            draw_signifpatch(tinx,sign_p,[.7 .5 .5])
        end
    end
    
    xti = [1 ts.BL round(ts.BL+ts.SL/3) round(ts.BL+ts.SL*2/3) stimend round(stimend+ts.PL/2) postend];
    evtime2 = round(ts.evtime*10)/10;
    xtilabs = arrayfun(@num2str, evtime2(xti), 'UniformOutput', false);
    
    xticks(xti); xticklabels(xtilabs); xlim(xti([1 end]))
    xlabel('Time rel. to period start (s)'); ylabel('Norm. MUA');
    title(figtit)
    legend(pL,pr.compar_tags)
    
    [~, fnm, ~] = save_fnms(muablockdir,ti,pr);
    
    saveas(fig,[fnm '_' Ltit '.jpg']);
    saveas(fig,[fnm '_' Ltit '.fig']);
    close(fig);
end

% 
% fig = figure;
% n = 1;
% for k = 1:2
%     
%     switch k;
%         case 1
%             chs = 1:8;
%             Ltit = 'Upper_half';
%         case 2
%             chs = 15:23;
%             Ltit = 'Lower_half';
%     end
%     for j = 1:3
%         switch j;
%             case 1; tinx = 1:ts.BL; Ptit = 'Baseline';
%             case 2; tinx = ts.BL+1:ts.BL+ts.SL; Ptit = 'Stim.';
%             case 3; tinx = ts.BL+ts.SL+1:ts.BL+ts.SL+ts.PL;  Ptit = 'Poststim.';
%         end;
%         
%         z_bp = cellfun(@(x) squeeze(mean( mean(x(chs,tinx,:),1),2)), z_groups_cc,'UniformOutput',0);
%         
%         subplot(3,3,n)
%         boxplot([z_bp{:}]); set_my_boxplot(gca);
%         [pval, ~, ~] = ranksum(z_bp{1},z_bp{2});
%         if pval<0.05; col = 'r'; else; col = 'k'; end;
%         yL = ylim;
%         text(0.5,yL(2)*.9,['p = ' num2str(pval)], 'Color',col)
%         title([Ltit ', ' Ptit]);
%         n = n+1;
%     end
% end
% 
% [~, fnm, ~] = save_fnms(muablockdir,ti,pr);
% 
% saveas(fig,[fnm '_SI_boxplots.jpg']);
% saveas(fig,[fnm '_SI_boxplots.fig']);
% close(fig);
end



function layer_by_layer_plot(z_groups_cc,z_groups_avg,ti,muablockdir,ts,pr)

figtit = {[pr.patname 'TH' num2str(ti)],['MUA blocks smoothed (' pr.smoothwith ', param:' num2str(pr.smparam) ')']};
cols = {'r','k'};

load(fullfile(pr.rootdir_ad,'channels_by_layers.mat'));
layers = channels_by_layers.(pr.patname)(:,ti);

basend = ts.BL;
stimend = ts.BL + ts.SL;
postend = ts.BL + ts.SL + ts.PL;


fig = figure;

xti = [1 ts.BL round(ts.BL+ts.SL/3) round(ts.BL+ts.SL*2/3) stimend round(stimend+ts.PL/2) postend];
evtime2 = round(ts.evtime*10)/10;
xtilabs = arrayfun(@num2str, evtime2(xti), 'UniformOutput', false);

for p = 1:2
    [~,laylabs] = plot_layerAVGs(z_groups_avg{p},layers,xti,xtilabs,[basend stimend],figtit,[],cols{p},true,fig);
end

pmap = [];
for j = 1:3
    switch j; case 1; tinx = 1:ts.BL; case 2; tinx = ts.BL+1:ts.BL+ts.SL; case 3; tinx = ts.BL+ts.SL+1:ts.BL+ts.SL+ts.PL; end;
    % tinx = 1:ts.BL+ts.SL+ts.PL;
    z_groups_bl = cellfun(@(x) x(:,tinx,:) , z_groups_cc,'UniformOutput',0);
    
    laymeans = cellfun(@(x) calc_layAVG(x,layers), z_groups_bl, 'UniformOutput',0);
    
    for jj = 1:length(layers)
        data1 = permute( laymeans{1}{jj} ,[3 1 2]);
        data2 = permute( laymeans{2}{jj} ,[3 1 2]);
        
        [permmaps,zmap_p,~] = my_permutationtest1(data1,data2,pr.itnum,0.05);
        
        if strcmp(pr.mcorr,'none')
            sign_z = zmap_p;
        elseif strcmp(pr.mcorr,'cluster')
            sign_z = my_cluster_correction(permmaps,zmap_p,pr.itnum,0.05);
        end
        sign_p = sign_z~=0;
        fig = gcf;
        fig.CurrentAxes = findobj(fig,'Type','Axes','Tag',laylabs{jj});
        draw_signifpatch(tinx,sign_p,[.7 .5 .5])
    end
end

[~, fnm, ~] = save_fnms(muablockdir,ti,pr);

saveas(fig,[fnm '_layers.jpg']);
saveas(fig,[fnm '_layers.fig']);
close(fig);

end


function [figdir, fnm, fnm2] = save_fnms(muablockdir,ti,pr);


if ~iscell(pr.evtype)
    evtylab = pr.evtype;
else
    c = cellfun(@(x) [x '_'],pr.evtype,'UniformOutput',0);
    evtylab = [c{:}];
end;


if pr.blockavg
    figdir = fullfile(muablockdir,[pr.compar_tags{1} '_vs_' pr.compar_tags{2} '_blockAVG'],evtylab); if ~isfolder(figdir); mkdir(figdir); end;
else
    figdir = fullfile(muablockdir,[pr.compar_tags{1} '_vs_' pr.compar_tags{2} '_concatblocks'],evtylab); if ~isfolder(figdir); mkdir(figdir); end;
end

if ~pr.intens_correct
    fnm = fullfile(figdir,[pr.patname '_TH' num2str(ti) '_' pr.smoothwith '_' num2str(pr.smparam) '_' evtylab 'evs']);
    fnm2 = fullfile(figdir,['StatMAT_' pr.patname '_TH' num2str(ti) '_' evtylab 'evs']);
else
    fnm = fullfile(figdir,[pr.patname '_TH' num2str(ti) '_' pr.smoothwith '_' num2str(pr.smparam) '_' evtylab 'evs_intenscorrect']);
    fnm2 = fullfile(figdir,['StatMAT_' pr.patname '_TH' num2str(ti) '_' evtylab 'evs_intenscorrect']);
end

if pr.maxintens
    fnm = [fnm  '_maxintens'];
    fnm2 = [fnm2  '_maxintens'];
end
fnm = [fnm '_STAT_' pr.stattype '_' pr.mcorr '_' pr.bnorm];
fnm2 = [fnm2 '_STAT_' pr.stattype '_' pr.mcorr '_' pr.bnorm '.mat'];
end