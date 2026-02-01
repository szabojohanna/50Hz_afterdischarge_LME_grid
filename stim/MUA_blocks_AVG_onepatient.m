function MUA_blocks_AVG_onepatient(rootdir_ad,patname,ti,varargin)

prs = inputParser;
addRequired(prs,'rootdir_ad',@isfolder);
addRequired(prs,'patname',@ischar);
addRequired(prs,'ti',@isnumeric);
addParameter(prs,'maxintens',true,@islogical);
addParameter(prs,'cLim',[-1.5 1.5],@isvector);
addParameter(prs,'blockavg',true,@islogical);
addParameter(prs,'smoothwith','smoothdata',@ischar); % 'none' | 'smoothdata' | 'lowpass'
addParameter(prs,'smparam',10,@isnumeric); % NaN if smoothdata is 'none'; 5-20 if smoothdata is 'smoothdata'; ~20 if smoothdata is 'lowpass';
addParameter(prs,'mcorr','fdr',@ischar); % 'none' | 'fdr';
addParameter(prs,'itnum',500,@isnumeric);
addParameter(prs,'minL',.8,@isnumeric);
addParameter(prs,'evtype','all',@(x) (ischar(x)|iscell(x))|isvector(x));
addParameter(prs,'timewin',[0 .9],@isvector); % relevant if evtype is 'mua_up'/'mua_down'; % [0 .45] | [0.23 .67] | [.45 .9] |  sec relative to PERIOD start
addParameter(prs,'upchans',[],@(x) isnumeric(x)|isvector(x));  % relevant if evtype is 'mua_up'/'mua_down'; if empty, all good channels
addParameter(prs,'period','stim',@ischar); % relevant if evtype is 'mua_up'/'mua_down'; bas | stim | post
addParameter(prs,'bnorm','Indiv',@ischar); % Normalization method: 'Zscore' (zscore for each channel, whole window) |  'Indiv' (indiv. baseline for each event)
parse(prs,rootdir_ad,patname,ti,varargin{:});
pr = prs.Results;




muablockdir = fullfile(rootdir_ad,'Figures','Micro','MUA','MUA_blocks');

load(fullfile(rootdir_ad,'Patients',patname,'Events.mat'));

load(fullfile(rootdir_ad,'channels_by_layers.mat'));
islaypat  = ismember(patname,channels_by_layers.Properties.VariableNames);

load(fullfile(rootdir_ad, 'good_channels.mat'))


%%

if ti==1; chans = 1:23; else; chans = 25:47; end;


[all_blocks,z_blocks_cc, z_blocks_avg, evinx2avg, ts] = load_evdat(ti, chans,muablockdir,true,pr);
if isempty(evinx2avg); return; end;

% avgavg_boxplot(all_blocks,z_blocks_cc,z_blocks_avg,ti,chans,evinx2avg,muablockdir,ts,pr)

ch_time_map(all_blocks,z_blocks_cc,z_blocks_avg,ti,chans,evinx2avg,muablockdir,ts,pr)


supra_infra_plot(patname,z_blocks_cc,evinx2avg,ti,muablockdir,ts,pr)
if islaypat;
    layer_by_layer_plot(z_blocks_cc,evinx2avg,ti,muablockdir,ts,pr)
end





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

function ch_time_map(all_blocks,z_blocks_cc,z_blocks_avg,ti,chans,evinx2avg,muablockdir,ts,pr)



%% Figure
basend = ts.BL;
stimend = ts.BL + ts.SL;
postend = ts.BL + ts.SL + ts.PL;


load(fullfile(pr.rootdir_ad,'channels_by_layers.mat'));
islaypat  = ismember(pr.patname,channels_by_layers.Properties.VariableNames);
if islaypat
    layers = channels_by_layers.(pr.patname)(:,ti);
end

fig = figure;

pcolor(z_blocks_avg); shading interp; set(gca,'YDir','reverse'); colorbar;

yL = ylim;
line([basend basend],yL,'Color','k','LineWidth',2,'LineStyle','-');
line([stimend stimend],yL,'Color','k','LineWidth',2,'LineStyle','--');
xti = [1 ts.BL round(ts.BL+ts.SL/3) round(ts.BL+ts.SL*2/3) stimend round(stimend+ts.PL/2) postend];
evtime2 = round(ts.evtime*10)/10;
xtiL = arrayfun(@num2str, evtime2(xti), 'UniformOutput', false);

xticks(xti)
xticklabels(xtiL);
xlabel('Time (s)');
ylabel('Channels')
caxis(pr.cLim)
evlab = arrayfun(@(x) ['Ev' num2str(x) '  '],evinx2avg,'UniformOutput',0);
if islaypat;  label_layers(layers); end;
colormap(jet)
% STAT
[exactp_ersp,maskersp,alpha2] = boostat_eeglab_J(z_blocks_cc,1:length(chans),0.05,pr.itnum,0,pr.mcorr,[],1:ts.BL);
zmap_p = norminv(exactp_ersp);
 bootstatFDR_clustercorr_AD(pr.rootdir_ad,pr.patname, pr.ti, zmap_p,maskersp,'maxsum')
% hold on; contour(exactp_ersp<0.05,'Color','m','LineWidth',2);

% if strcmp(pr.mcorr,'none'); mc = false; else; mc = true; end;
% data = permute(z_blocks_cc,[3 1 2]);
% [zmap_cl_p, ~, ~,~] = bas_permstat(data,1:ts.BL,pr.itnum,0.1,0.05,mc);
% hold on; contour(zmap_cl_p~=0,'Color','w','LineWidth',2);



%     set(fig,'Position',get(0,'Screensize'))

[figdir, fnm, fnm2,titlab] = save_fnms(muablockdir,ti,pr,evlab);

title({[pr.patname 'TH' num2str(ti) ', Ev nr = ' num2str(length(evinx2avg)) ],titlab,...
    ['MUA blocks smoothed (' pr.smoothwith ', param:' num2str(pr.smparam) ')']});

saveas(fig,[fnm '.jpg']);
saveas(fig,[fnm '.fig']);
saveas(fig,[fnm '.svg']);
close(fig);

StatMAT.z_blocks_cc = z_blocks_cc;
StatMAT.all_blocks = all_blocks;
StatMAT.exactp_ersp = exactp_ersp;
StatMAT.maskersp = maskersp;
StatMAT.clustform_alpha = 0.1;
StatMAT.clustthr_alpha = 0.05;
StatMAT.evinx2avg = evinx2avg;
StatMAT.basline_end = basend;
StatMAT.stim_end = stimend;
StatMAT.data_smooth = [pr.smoothwith ', param: ' num2str(pr.smparam)];
StatMAT.timevec = ts.evtime;

save(fnm2, 'StatMAT')
end

function [figdir, fnm, fnm2,titlab] = save_fnms(muablockdir,ti,pr,evlab)

if ischar(pr.evtype)
    evtylab = pr.evtype;
elseif iscell(pr.evtype)
    c = cellfun(@(x) [x '_'],pr.evtype,'UniformOutput',0);
    evtylab = [c{:}];
elseif isvector(pr.evtype)
    c = arrayfun(@num2str,pr.evtype,'UniformOutput',0);
    evtylab = [c{:}];
end;

if pr.blockavg
    figdir = fullfile(muablockdir,'allzones_blockAVG',evtylab); if ~isfolder(figdir); mkdir(figdir); end;
else
    figdir = fullfile(muablockdir,'allzones_concatblocks',evtylab); if ~isfolder(figdir); mkdir(figdir); end;
end


if ~pr.maxintens
    
    titlab =  ['All close stim sites: ' evlab{:}];
    
    fnm = fullfile(figdir,[pr.patname '_TH' num2str(ti) '_' pr.smoothwith '_' num2str(pr.smparam) '_' evtylab 'evs']);
    fnm2 = fullfile(figdir,['StatMAT_' pr.patname '_TH' num2str(ti) '_' evtylab 'evs']);
else
    
    titlab =  ['Close stim sites with max Int: ' evlab{:}];
    
    fnm = fullfile(figdir,[pr.patname '_TH' num2str(ti) '_' pr.smoothwith '_' num2str(pr.smparam) '_' evtylab 'evs_maxintens']);
    fnm2 = fullfile(figdir,['StatMAT_' pr.patname '_TH' num2str(ti) '_' evtylab 'evs_maxintens']);
end


fnm = [fnm '_STAT' pr.mcorr '_' pr.bnorm];
fnm2 = [fnm2 '_STAT' pr.mcorr '_' pr.bnorm '.mat'];



end

function layer_by_layer_plot(z_blocks_cc,evinx2avg,ti,muablockdir,ts,pr)

basend = ts.BL;
stimend = ts.BL + ts.SL;
postend = ts.BL + ts.SL + ts.PL;

load(fullfile(pr.rootdir_ad,'channels_by_layers.mat'));
layers = channels_by_layers.(pr.patname)(:,ti);


laynr = sum( ~cellfun(@isempty,layers) );
lay_cc = cell(1,laynr);
for n = 1:laynr
    chs = layers{n};
    lay_cc{n} =  z_blocks_cc(chs,:,:) ;
end
laymean = cellfun(@(x) mean(x,[1 3]) ,lay_cc, 'UniformOutput',0);
fig = figure;
nn = 1;
mxa = max(cat(1,laymean{:}),[],'all');

for n = laynr:-1:1
    yti(nn) = (nn-1)*2*mxa;
    plot(laymean{n}+yti(nn),'Color','k'); hold on;
    nn = nn+1;
end
if laynr==7
    laylabs = [ {'WM'}, arrayfun(@(x) ['Layer ' num2str(x)], 6:-1:1, 'UniformOutput',0)];
else
    laylabs = arrayfun(@(x) ['Layer ' num2str(x)], laynr:-1:1, 'UniformOutput',0);
end
yticks(yti); yticklabels(laylabs);
yL = ylim;
line([basend basend],yL,'Color','k','LineStyle','-')
line([stimend stimend],yL,'Color','k','LineStyle','--');
xlim([0 size(z_blocks_cc,2)]);
xti = [1 ts.BL round(ts.BL+ts.SL/3) round(ts.BL+ts.SL*2/3) stimend round(stimend+ts.PL/2) postend];
evtime2 = round(ts.evtime*10)/10;
xtiL = arrayfun(@num2str, evtime2(xti), 'UniformOutput', false);

xticks(xti)
xticklabels(xtiL);
xlabel('Time relative to stim. start (s) ')

evlab = arrayfun(@(x) ['Ev' num2str(x) '  '],evinx2avg,'UniformOutput',0);

[figdir, fnm, fnm2,titlab] = save_fnms(muablockdir,ti,pr,evlab);

title({[pr.patname 'TH' num2str(ti) ', Ev nr = ' num2str(length(evinx2avg))],titlab,...
    ['MUA blocks smoothed (' pr.smoothwith ', param:' num2str(pr.smparam) ')']});

saveas(fig,[fnm '_LAYERplot.jpg'])
saveas(fig,[fnm '_LAYERplot.fig']);
saveas(fig,[fnm '_LAYERplot.svg']);
close(fig);
end


function supra_infra_plot(patname,z_blocks_cc,evinx2avg,ti,muablockdir,ts,pr)

laylabs = {'Infragran','Gran','Supragran'};
basend = ts.BL;
stimend = ts.BL + ts.SL;
postend = ts.BL + ts.SL + ts.PL;

load(fullfile(pr.rootdir_ad,'channels_by_layers.mat'));
islaypat  = ismember(pr.patname,channels_by_layers.Properties.VariableNames);
if islaypat
    layers = channels_by_layers.(pr.patname)(:,ti);
else
    layers = {};
end


laynr = 3;
lay_cc = cell(1,laynr);
for k = 1:laynr
    if ~isempty(layers)
        switch k;
            case 1
                chs = [layers{1:3}];
            case 2
                chs = [layers{4}];
            case 3
                chs = [layers{5:6}];
        end
    else
        chs = nolaypats_if(patname,laylabs{k});
    end
    if ~isempty(chs)
        try
        lay_cc{k} =  z_blocks_cc(chs,:,:) ;
        catch
            disp()
        end
    end
end
nept = ~cellfun(@isempty,lay_cc);
laynr = sum(nept);
laymean = cellfun(@(x) mean(x,[1 3]) ,lay_cc(nept), 'UniformOutput',0);
laym = cellfun(@(x) mean(x,1) ,lay_cc(nept), 'UniformOutput',0);

laysd = cellfun(@(x) std( x,[],3) ./ size(x,3), laym, 'UniformOutput',0);
sdup = arrayfun(@(x) laymean{x}+laysd{x}, 1:laynr,'UniformOutput',0);
sddown = arrayfun(@(x) laymean{x}-laysd{x}, 1:laynr,'UniformOutput',0);

fig = figure;
nn = 1;
mxa = max(cat(1,laymean{:}),[],'all');
t = 1:length(laymean{1});
for n = laynr:-1:1
    yti(nn) = (nn-1)*2*mxa;
    plot(laymean{n}+yti(nn),'Color','k'); hold on;
    
    patch('XData',[t fliplr(t)], 'YData',[sddown{n}+yti(nn) fliplr(sdup{n}+yti(nn))],...
        'FaceAlpha',0.2,'FaceColor','k','EdgeColor','none');
    
%     [exactp_ersp,maskersp,alpha2] = boostat_eeglab_J(laym{n},1,0.05,pr.itnum,0,pr.mcorr,[],1:ts.BL);
data = permute(laym{n},[3 1 2]);
if strcmp(pr.mcorr,'none'); mc = false; else; mc = true; end;
[zmap_cl_p, ~, ~,~] = bas_permstat(data,1:ts.BL,pr.itnum,0.1,0.05,mc);


    hold on;
    signp = find(zmap_cl_p~=0);
    scatter(signp,zeros(1,length(signp))+yti(nn),[],'k','filled');
    nn = nn+1;
    
end

laylabs = laylabs(nept);
yticks(yti); yticklabels(laylabs);
yL = ylim;
line([basend basend],yL,'Color','k','LineStyle','-')
line([stimend stimend],yL,'Color','k','LineStyle','--');
xlim([0 size(z_blocks_cc,2)]);
xti = [1 ts.BL round(ts.BL+ts.SL/3) round(ts.BL+ts.SL*2/3) stimend round(stimend+ts.PL/2) postend];
evtime2 = round(ts.evtime*10)/10;
xtiL = arrayfun(@num2str, evtime2(xti), 'UniformOutput', false);

xticks(xti)
xticklabels(xtiL);

evlab = arrayfun(@(x) ['Ev' num2str(x) '  '],evinx2avg,'UniformOutput',0);

[figdir, fnm, fnm2,titlab] = save_fnms(muablockdir,ti,pr,evlab);

title({[pr.patname 'TH' num2str(ti) ', Ev nr = ' num2str(length(evinx2avg))],titlab,...
    ['MUA blocks smoothed (' pr.smoothwith ', param:' num2str(pr.smparam) ')']});

saveas(fig,[fnm '_SupraInfra.jpg'])
saveas(fig,[fnm '_SupraInfra.fig']);
close(fig);

% Compare domains
fig = figure;

cols = {'k','b','g'}; cols = cols(nept);
comb  = nchoosek(1:laynr,2);
cnr = size(comb,1);
laylabs2 = flip(laylabs);
for c = 1:cnr
    subplot(cnr,1,c)
    c1 = comb(c,1);
    c2 = comb(c,2);
    nn = 1; pL = [];
    for n = comb(c,:)
        
        pL(nn) = plot(laymean{n},'Color',cols{n}); hold on;
        
        patch('XData',[t fliplr(t)], 'YData',[sddown{n} fliplr(sdup{n})],...
            'FaceAlpha',0.2,'FaceColor',cols{n},'EdgeColor','none');
        nn = nn+1;
    end
    legend(pL,laylabs2(comb(c,:)),'AutoUpdate','off')
    data1 = permute( laym{c1} ,[3 1 2]);
    data2 = permute( laym{c2} ,[3 1 2]);
    
    [permmaps,zmap_p,~] = my_permutation_Ttest(data1,data2,pr.itnum,0.05);
    [signmap] = my_cluster_correction(permmaps,zmap_p,pr.itnum,0.05);
    sign_p = signmap~=0;
    yL = pr.cLim; ylim(yL)
    draw_signifpatch(t,sign_p,[.7 .5 .5])
    
    line([basend basend],yL,'Color','k','LineStyle','-')
    line([stimend stimend],yL,'Color','k','LineStyle','--');
    
    xlim([0 size(z_blocks_cc,2)]);
    xti = [1 ts.BL round(ts.BL+ts.SL/3) round(ts.BL+ts.SL*2/3) stimend round(stimend+ts.PL/2) postend];
    evtime2 = round(ts.evtime*10)/10;
    xtiL = arrayfun(@num2str, evtime2(xti), 'UniformOutput', false);
    
    xticks(xti)
    xticklabels(xtiL);
    ylabel('Norm. MUA'); xlabel('Time rel. to period start')
end
suptitle({[pr.patname 'TH' num2str(ti) ', Ev nr = ' num2str(length(evinx2avg))],titlab,...
    ['MUA blocks smoothed (' pr.smoothwith ', param:' num2str(pr.smparam) ')']});

saveas(fig,[fnm '_SupraInfra_comp.jpg'])
saveas(fig,[fnm '_SupraInfra_comp.fig']);
close(fig);

% fig = figure;
% for j = 1:3
%     switch j; case 1; tinx = 1:ts.BL; case 2; tinx = ts.BL+1:ts.BL+ts.SL; case 3; tinx = ts.BL+ts.SL+1:ts.BL+ts.SL+ts.PL; end;
%
%     bp = cellfun(@(x) squeeze(mean(x(:,tinx,:),2)),laym([1 3]), 'UniformOutput', false);
%     subplot(1,3,j)
%     boxplot([bp{:}])
%     [pval(j),~,~] = ranksum(bp{1},bp{2});
% %     [pval(j),anovatab,stats] = kruskalwallis([bp{:}]);
% %     multcompare(stats)
% end
end