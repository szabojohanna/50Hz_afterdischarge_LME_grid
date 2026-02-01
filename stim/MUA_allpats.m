function MUA_allpats(rootdir_ad,pats,thumbs,varargin)

prs = inputParser;
addRequired(prs,'rootdir_ad',@isfolder);
addRequired(prs,'pats',@iscell);
addRequired(prs,'thumbs',@isvector);
addParameter(prs,'maxintens',true,@islogical);
addParameter(prs,'cLim',[-1.5 1.5],@isvector);
addParameter(prs,'blockavg',true,@islogical);
addParameter(prs,'smoothwith','smoothdata',@ischar); % 'none' | 'smoothdata' | 'lowpass'
addParameter(prs,'smparam',10,@isnumeric); % NaN if smoothdata is 'none'; 5-20 if smoothdata is 'smoothdata'; ~20 if smoothdata is 'lowpass';
addParameter(prs,'mcorr','fdr',@ischar); % 'none' | 'fdr';
addParameter(prs,'itnum',1000,@isnumeric);
addParameter(prs,'minL',.8,@isnumeric);
addParameter(prs,'evtype','all',@(x) ischar(x)|iscell(x));
addParameter(prs,'timewin',[0 .9],@isvector); % relevant if evtype is 'mua_up'/'mua_down'; % [0 .45] | [0.23 .67] | [.45 .9] |  sec relative to PERIOD start
addParameter(prs,'upchans',[],@(x) isnumeric(x)|isvector(x));  % relevant if evtype is 'mua_up'/'mua_down'; if empty, all channels
addParameter(prs,'period','stim',@ischar); % relevant if evtype is 'mua_up'/'mua_down'; bas | stim | post
addParameter(prs,'bnorm','Indiv',@ischar); % Normalization method: 'Zscore' (zscore for each channel, whole window) |  'Indiv' (indiv. baseline for each event)
addParameter(prs,'patavg',true,@islogical); % Normalization method: 'Zscore' (zscore for each channel, whole window) |  'Indiv' (indiv. baseline for each event)
parse(prs,rootdir_ad,pats,thumbs,varargin{:});
pr = prs.Results;



[z_blM, mb,ms,mp,figdir] = load_allpat_dat(rootdir_ad,pats,thumbs,pr);


supra_infra_plot(rootdir_ad,z_blM, pats,thumbs,mb,ms,mp,pr,figdir)
end




function supra_infra_plot(rootdir_ad,z_blM, pats,thumbs,mb,ms,mp,pr,figdir)

[z_avgs, laylabs, laynr] = select_layers_by_patients(rootdir_ad,pats,thumbs,z_blM);


nempt = ~cellfun(@isempty,z_avgs);

zp_avgs = cell(size(z_avgs));
if pr.patavg
zp_avgs = cellfun(@(x) mean(x,3) ,z_avgs, 'UniformOutput',0); % average across events within patients
else
    zp_avgs(nempt) = cellfun(@(x) permute(x,[3 2 1]), z_avgs(nempt), 'UniformOutput',0);
end

% poola = cat(1,zp_avgs{:,[1 2]}); % pool for baseline corr
% bas_a = repmat( mean( poola(:,1:mb), 'all' ) ,[1 mb+ms+mp]);
% bas_sd = repmat( std( poola(:,1:mb),[], 'all' ) ,[1 mb+ms+mp]);
poola = arrayfun(@(x)  cat(1,zp_avgs{:,x}), 1:laynr, 'UniformOutput',0);
bas_a = cellfun(@(x)  repmat( mean( x(:,1:mb), 'all' ) ,[1 mb+ms+mp])  ,poola, 'UniformOutput',0);
bas_sd = cellfun(@(x)  repmat( std( x(:,1:mb),[], 'all' ) ,[1 mb+ms+mp])  ,poola, 'UniformOutput',0);
a_bascor = cell(length(pats), laynr);


for n = 1:laynr    
    a_bascor(nempt(:,n),n) = cellfun(@(x) (x-bas_a{n})./bas_sd{n} ,zp_avgs(nempt(:,n),n), 'UniformOutput',0);
end


fig = figure;
yL = [-3 3];
% yL = [-2 2];
for n = 1:laynr
    
    switch n;
        case 1; Ltit = 'Supragran'; spnr = 1;
        case 2; Ltit = 'Infragran'; spnr = 3;
        case 3; Ltit = 'Gran'; spnr = 2;
    end
    laymm = zp_avgs(nempt(:,n),n);
    xx = cat(1,laymm{:});
    laymean = mean( xx ,1);
    
    laysd= std( xx ,[],1)./ sqrt(size(xx,1));
    sdup = laymean+laysd;
    sddown = laymean-laysd;
    
    laypats = cat(1,laymm{:});
    t = 1:length(laymean);
    
    subplot(3,1,spnr);
    plot(laymean,'Color','k'); hold on;
    
    patnr_si = size(laypats,1);
    patnr_g = sum(nempt(:,3));
    
    patch('XData',[t fliplr(t)], 'YData',[sddown fliplr(sdup)],...
        'FaceAlpha',0.2,'FaceColor','k','EdgeColor','none');
    
        [exactp_ersp,maskersp,alpha2] = boostat_eeglab_J(permute(laypats,[3 2 1]),1,0.05,pr.itnum,0,pr.mcorr,[],1:mb);
        signp = (exactp_ersp<0.05)';
%     if strcmp(pr.mcorr,'none'); mc = false; else; mc = true; end;
%     data = permute(laypats,[3 1 2]);
% data = permute(laypats,[1 3 2]);
%     [zmap_cl_p, ~, ~,~] = bas_permstat(data,1:mb,pr.itnum,0.05,0.05,mc);
% signp = zmap_cl_p~=0;
    
    hold on;
    ylim(yL);
    draw_signifpatch(t,signp,[.5 .5 .5])
    title([Ltit ', infra/supra: n = ' num2str(patnr_si) ', gran: n = ' num2str(patnr_g) ]);
    line([mb mb],yL,'Color','k','LineStyle','-')
    line([mb+ms mb+ms],yL,'Color','k','LineStyle','--');
    xlim([0 size(t,2)]);
    
    
    if ~iscell(pr.evtype)
        evtylab = pr.evtype;
    else
        c = cellfun(@(x) [x '_'],pr.evtype,'UniformOutput',0);
        evtylab = [c{:}];
    end;
    
    
    laymm = a_bascor(nempt(:,n),n);
    xxx = cat(1,laymm{:});
    laymeana = mean( xxx ,1);
    laysda = std( xxx ,[],1)./ sqrt( size(xxx,1) );
    sdupa = laymeana+laysda;
    sddowna = laymeana-laysda;
    
    
    LM{n} = laymeana;
    LP{n} = laypats;
    LSd{n} = sddowna;
    LSu{n} = sdupa;
end

fnm = fullfile(figdir,['PATS_' pats{:} '_' pr.smoothwith '_' num2str(pr.smparam) '_' evtylab 'evs']);
if pr.maxintens
    fnm = [fnm '_maxintens'];
end
fnm = [fnm '_STAT' pr.mcorr '_SupraGranInfra'];
saveas(fig,[fnm '.jpg'])
saveas(fig,[fnm '.fig'])
saveas(fig,[fnm '.pdf'])
close(fig)



% Compare domains
yL = [-2 4];
fig = figure;

cols = {'k','g','b'};
comb  = nchoosek(1:laynr,2);
cnr = size(comb,1);
laylabs = {'Supragran.','Infagran','Gran.'};
for c = 1:cnr
    subplot(cnr,1,c)
    c1 = comb(c,1);
    c2 = comb(c,2);
    nn = 1; pL = [];
    for n = comb(c,:)
        
        pL(nn) = plot(LM{n},'Color',cols{n}); hold on;
        
        patch('XData',[t fliplr(t)], 'YData',[LSd{n} fliplr(LSu{n})],...
            'FaceAlpha',0.2,'FaceColor',cols{n},'EdgeColor','none');
        nn = nn+1;
    end
    legend(pL,laylabs(comb(c,:)),'AutoUpdate','off')
    data1 = permute( LP{c1} ,[1 3 2]);
    data2 = permute( LP{c2} ,[1 3 2]);
    
    [permmaps,zmap_p,~] = my_permutationtest1(data1,data2,pr.itnum,0.05);
    [signmap] = my_cluster_correction(permmaps,zmap_p,pr.itnum,0.05);
    sign_p = signmap~=0;
    %     yL = [-4.5 4.5];
    %     yL = ylim;
    %     ylim(yL)
    draw_signifpatch(t,sign_p,[.5 .5 .5])
    
    line([mb mb],yL,'Color','k','LineStyle','-')
    line([mb+ms mb+ms],yL,'Color','k','LineStyle','--');
    xlim([0 size(t,2)]);
    ylim(yL);
    
    ylabel('Norm. MUA'); xlabel('Time rel. to period start')
end
suptitle({['Patients: ' pats{nempt(:,1)} ],['infra/supra: n = ' num2str(patnr_si) ', gran: n = ' num2str(patnr_g) ],...
    'MUA avg. across all channels'});

fnm = fullfile(figdir,['PATS_' pats{:} '_' pr.smoothwith '_' num2str(pr.smparam) '_' evtylab 'evs']);
if pr.maxintens
    fnm = [fnm '_maxintens'];
end
fnm = [fnm '_STAT' pr.mcorr];
saveas(fig,[fnm '_SupraInfra_comp.jpg'])
saveas(fig,[fnm '_SupraInfra_comp.fig']);
saveas(fig,[fnm '_SupraInfra_comp.pdf']);
close(fig);




end


