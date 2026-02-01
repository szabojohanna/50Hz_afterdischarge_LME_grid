function MUA_stimint_corr_onepatient_ONESITE(rootdir_ad,patname,ti,what2corr,varargin)

prs = inputParser;
addRequired(prs,'rootdir_ad',@isfolder);
addRequired(prs,'patname',@ischar);
addRequired(prs,'ti',@isnumeric);
addRequired(prs,'what2corr',@ischar); %'int' | 'dist'
addParameter(prs,'maxintens',true,@islogical);
addParameter(prs,'corrmap',true,@islogical);
addParameter(prs,'corrLine',true,@islogical);
addParameter(prs,'evtype','all', @(x) ischar(x)|iscell(x) );
parse(prs,rootdir_ad,patname,ti,what2corr,varargin{:});
pr = prs.Results;

% what2corr: int | dist


%%
% [z_blocks_cc, ints, StatMAT] = load_dats(rootdir_ad,patname,ti,what2corr,pr);
[~,z_blocks_cc, ~, ~, ~,ints] = load_evdat(ti, chans,muablockdir,ishem,pr);

ChT_map_stimlocs(rootdir_ad,z_blocks_cc, ints,what2corr,patname,ti,StatMAT,pr.corrmap, pr.corrLine);


%%






end





function [z_blocks_cc, ints, StatMAT] = load_dats(rootdir_ad,patname,ti,what2corr,pr)


patdir = fullfile(rootdir_ad,'Patients',patname);
load(fullfile(patdir,'Events.mat'))

muablockdir = fullfile(rootdir_ad,'Figures','Micro','MUA','MUA_blocks');
admap_dir = fullfile(muablockdir,'allzones_blockAVG');

if pr.maxintens
    fnm = fullfile(admap_dir,['StatMAT_' patname '_TH' num2str(ti) '_' pr.evtype 'evs_maxintens_STATfdr.mat']);
else
    fnm = fullfile(admap_dir,['StatMAT_' patname '_TH' num2str(ti) '_' pr.evtype 'evs_STATfdr.mat']);
end


load(fnm);

all_blocks = StatMAT.all_blocks;
basend = StatMAT.basline_end;
stimend = StatMAT.stim_end;

aevents = StatMAT.evinx2avg;

[~, z_blocks_cc ]  = basnorm_blocks(all_blocks,basend);


if contains(what2corr,'int')
    ints = [Events(aevents).StimIntensity];
elseif contains(what2corr,'dist')
    ints = {Events(aevents).(['Th' num2str(ti) '_STIMdist_cm'])};
end
end




function ChT_map_stimlocs(rootdir_ad,z_blocks_cc, ints,what2corr,patname,ti,StatMAT,corrmap, corrLine)

patdir = fullfile(rootdir_ad,'Patients',patname);
load(fullfile(patdir,'Events.mat'))

muablockdir = fullfile(rootdir_ad,'Figures','Micro','MUA','MUA_blocks');

if exist(fullfile(rootdir_ad,'stimsite_types.mat'))==2
    load(fullfile(rootdir_ad,'stimsite_types.mat'))
else
    stimsite_types = struct;
    stimsite_types.(patname) = struct;
end



% Directory to save figure
figdir = fullfile(muablockdir,['MUA_' what2corr 'correlations'],'Site_corr');
if ~isdir(figdir); mkdir(figdir); end;

load(fullfile(rootdir_ad,'channels_by_layers.mat'));
islaypat  = ismember(patname,channels_by_layers.Properties.VariableNames);
if islaypat;
    layers = channels_by_layers.(patname)(:,ti);
else
    layers = [];
end;

basend = StatMAT.basline_end;
stimend = StatMAT.stim_end;

aevents = StatMAT.evinx2avg;


stim_pairs = arrayfun(@(x) [Events(x).StimChans{1} '_' Events(x).StimChans{2}], aevents, 'UniformOutput',0 );
[uniq_sp, ~, tag_occur] = unique(stim_pairs,'stable');
uniq_sp_nr = length(uniq_sp);

mups = {};mdos = {}; mnos = {};
for m = 1:uniq_sp_nr
    evx = tag_occur==m; % events of one stim. loc
    
    
    stimp_nm = uniq_sp{m};
    stimp_evs = aevents(evx);
    
    
    evlab = arrayfun(@(x) ['Ev' num2str(x) '  '],stimp_evs,'UniformOutput',0);
    
    if corrmap
        chnr = size(z_blocks_cc,1); tnr = size(z_blocks_cc,2);
        [coef, pval] = deal(nan(chnr,tnr));
        for j = 1:chnr % loop over channels
            for k = 1:tnr % loop over time points
                muav = squeeze(z_blocks_cc(j,k,evx));
                [coef(j,k), pval(j,k)] =  corr(muav,ints(evx)','type','Spearman');
                
            end
            
            
        end
        
        act_islands = bwconncomp(pval<0.05);
        
        pval_corr = pval;
        if numel(act_islands.PixelIdxList)>0
            act_clustsizes = cellfun(@length,act_islands.PixelIdxList);
            
            % Exclude clusters that are too small (<8 points)
            f = find(act_clustsizes<8);
            
            for h = f
                pval_corr(act_islands.PixelIdxList{h}) = 1;
            end
        else;
            pval_corr(pval_corr<0.05) = 1; % Exclude single-point significant p values
        end
        pval_corr(pval_corr>0.05) = 1;
        
        stp = Events(stimp_evs(1)).StimChans;
        if any(pval_corr<0.05,'all')&& length(stimp_evs)>1
            cM = mean(coef(pval_corr<0.05));
            if cM>0 % events with overall positive MUA-stim int. correlation
                mups = cat(1,mups,{stp,stimp_evs});
            else % events with overall negative MUA-stim int. correlation
                mdos = cat(1,mdos,{stp,stimp_evs});
            end
        else
            mnos = cat(1,mnos,{stp,stimp_evs});
        end
        
        
        fig = coefmap(coef, pval_corr, basend, stimend, stimp_evs,islaypat,layers);
        title({['Stim. site: ' stimp_nm],['events: ' evlab{:}]});
        
        fnm = fullfile(figdir,[patname '_TH' num2str(ti) '_' stimp_nm]);
        
        saveas(fig,[fnm '.jpg']);
        saveas(fig,[fnm '.fig']);
        close(fig);
        
        CorrMAPs.(stimp_nm).coef = coef;
        CorrMAPs.(stimp_nm).pval = pval;
        CorrMAPs.(stimp_nm).pval_corr = pval_corr;
        CorrMAPs.(stimp_nm).evs = stimp_evs;
        CorrMAPs.(stimp_nm).stim_pair =  Events(stimp_evs(1)).StimChans;
        CorrMAPs.(stimp_nm).basend =  basend;
        CorrMAPs.(stimp_nm).stimend =  stimend;
    end
    
    if corrLine
        
        % Line plot
        
        avg_muav = squeeze( mean(z_blocks_cc(:,basend+1:stimend,evx),[1 2]) );
        if length(avg_muav)<4; continue; end;
        
        
        [avg_coef, avg_pval] =  corr(ints(evx)',avg_muav,'type','Spearman');
        
        
        fig = figure;
        scatter(ints(evx)',avg_muav,[],'k','filled'); hold on;
        
        if avg_coef<.05; col = 'r'; else; col = 'k'; end;
        yL = ylim; xL = xlim; text(xL(1),yL(2)*.9,['p = ' num2str(avg_pval) ...
            ', R = ' num2str(avg_coef)],'Color',col);
        
        mdl_ad = fitlm(ints(evx)',avg_muav);
        P_h = plot(mdl_ad); hold on;
        
        set(findobj(P_h,'Type','Line'),'LineWidth',2,'Color','k')
        xlabel('Avg. norm. MUA'); ylabel('Stim. intensity (mA)');
        title({['Stim. site: ' stimp_nm],['events: ' evlab{:}]});

        
        fnm = fullfile(figdir,[patname '_TH' num2str(ti) '_' stimp_nm '_STIMperiodAVG']);
        
        saveas(fig,[fnm '.jpg']);
        saveas(fig,[fnm '.fig']);
        close(fig);
        
        
    end
    
end


if corrmap
    stimsite_types.(patname)(ti).mua_up = mups;
    stimsite_types.(patname)(ti).mua_down  = mdos;
    stimsite_types.(patname)(ti).no_change = mnos;
    save(fullfile(rootdir_ad,'stimsite_types.mat'),'stimsite_types');
    save(fullfile(figdir,[patname '_Th' num2str(ti) '_CorrMAPs.mat']), 'CorrMAPs' )
end
end






