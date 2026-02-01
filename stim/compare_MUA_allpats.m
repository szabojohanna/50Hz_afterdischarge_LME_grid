function compare_MUA_allpats(rootdir_ad,pats,thumbs,compar_tags,varargin)


prs = inputParser;
addRequired(prs,'rootdir_ad',@isfolder);
addRequired(prs,'pats',@iscell);
addRequired(prs,'thumbs',@isvector);
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
addParameter(prs,'itnum',1000,@isnumeric);
addParameter(prs,'minL',.8,@isnumeric);
addParameter(prs,'evtype','all',@(x) ischar(x)|iscell(x));
addParameter(prs,'timewin',[0 .9],@isvector); % relevant if evtype is 'mua_up'/'mua_down'; % [0 .45] | [0.23 .67] | [.45 .9] |  sec relative to PERIOD start
addParameter(prs,'upchans',[],@(x) isnumeric(x)|isvector(x));  % relevant if evtype is 'mua_up'/'mua_down'; if empty, all channels
addParameter(prs,'period','stim',@ischar); % relevant if evtype is 'mua_up'/'mua_down'; bas | stim | post
addParameter(prs,'bnorm','Indiv',@ischar); % Normalization method: 'Zscore' (zscore for each channel, whole window) | 'Common' (comm. baseline across each events) | 'Indiv' (indiv. baseline for each event)
addParameter(prs,'patavg',true,@islogical);
parse(prs,rootdir_ad,pats,thumbs,compar_tags,varargin{:});
pr = prs.Results;


muablockdir = fullfile(rootdir_ad,'Figures','Micro','MUA','MUA_blocks');

[a_grs, ispat,mb,ms,mp,pr,figdir,a_lays,laylabs] = load_groups_allpats(rootdir_ad,pats,thumbs,compar_tags,'patavg',pr.patavg,varargin{:});

allch_plot(a_grs, pats(ispat),mb,ms,mp,pr,figdir);

suprainfra_plot(a_grs, pats(ispat),mb,ms,mp,a_lays,laylabs,pr,figdir)

% laych_plot(z_grsM(nemp), pats(nemp),mb,ms,mp,layers(ispat),pr,figdir)

end


function allch_plot(a_grs, pats,mb,ms,mp,pr,figdir)


if ~iscell(pr.evtype)
    evtylab = pr.evtype;
else
    c = cellfun(@(x) [x '_'],pr.evtype,'UniformOutput',0);
    evtylab = [c{:}];
end;

% All ch avg
z_avgs =  cellfun(@(y) permute( nanmean(y,1) ,[3 2 1]),a_grs,'UniformOutput',0);


% Plot all ch avg
cols = {'r','k'};
fig = figure;
pL =plotstat(z_avgs, pats,mb,ms,mp,pr,cols)
legend(pL,pr.compar_tags)
fnm = fullfile(figdir,['PATS_' pats{:} '_' pr.smoothwith '_' num2str(pr.smparam) '_' evtylab 'evs']);
if pr.intens_correct
    fnm = [fnm '_intenscorrect'];
end
fnm = [fnm '_STAT' pr.mcorr '_ChAVG'];
if pr.patavg; fnm  = [fnm '_patAVG']; end;
saveas(fig,[fnm '.jpg'])
saveas(fig,[fnm '.fig'])
close(fig)

end

function suprainfra_plot(a_grs, pats,mb,ms,mp,a_lays,laylabs,pr,figdir)


if ~iscell(pr.evtype)
    evtylab = pr.evtype;
else
    c = cellfun(@(x) [x '_'],pr.evtype,'UniformOutput',0);
    evtylab = [c{:}];
end;


laynr = 3;
patnr = length(pats);


for L = 1:laynr
    
   Ltit = laylabs{L};
   
   grdat = cell(1,2);
   for g = 1:2
       gr = a_grs{g};
           gr_avgs = cell(size(gr,3),1);
       for j = 1:size(gr,3)
           gr_L = a_lays{g}{L,j};
           if ~isempty(gr_L)
               gr_avgs{j} = nanmean(  gr(gr_L,:,j) ,1);
           else
               gr_avgs{j} = [];
           end
       end
       grdat{g} = cat(1,gr_avgs{:});
       
   end
%    poola = cat(1, grdat{:});
   
   grdat_N = cell(1,2);
   for g = 1:2
%        bas_a = repmat( mean( poola(:,1:mb), 'all') ,[size(grdat{g},1) size(poola,2)]);
%        bas_sd = repmat( std( poola(:,1:mb),[], 'all') ,[size(grdat{g},1) size(poola,2)]);
       
       bas_a = repmat( nanmean( grdat{g}(:,1:mb), 'all') ,[size(grdat{g},1) size(grdat{g},2)]);
       bas_sd = repmat( nanstd( grdat{g}(:,1:mb),[], 'all') ,[size(grdat{g},1) size(grdat{g},2)]);
       grdat_N{g} = (grdat{g}-bas_a)./bas_sd;
   end
   
   
    cols = {'r','k'};
    fig = figure;
    pL = plotstat(grdat_N, pats,mb,ms,mp,pr,cols);
    legend(pL,pr.compar_tags)
    ylim([-6 6])
    fnm = fullfile(figdir,['PATS_' pats{:} '_' pr.smoothwith '_' num2str(pr.smparam) '_' evtylab 'evs']);
    if pr.intens_correct
        fnm = [fnm '_intenscorrect'];
    end
    fnm = [fnm '_STAT' pr.mcorr '_' Ltit];
    if pr.maxintens
        fnm = [fnm  '_maxintens'];
    end
    if pr.patavg; fnm  = [fnm '_patAVG']; end;
    saveas(fig,[fnm '.jpg'])
    saveas(fig,[fnm '.fig'])
    close(fig)
    
    LZ(L,1:2) = grdat;
end




% Compare domains
cols = {'k','g'};
fig = figure;
for g = 1:2
    Ldats = LZ(:,g);
%     poola = cat(1, Ldats{1:2});
    
    Ldats_N = cell(2,1);
    for Lc = 1:2
%         bas_a = repmat( mean( poola(:,1:mb), 'all') ,[size(Ldats{Lc},1) size(poola,2)]);
%         bas_sd = repmat( std( poola(:,1:mb),[], 'all') ,[size(Ldats{Lc},1) size(poola,2)]);
        
        
        bas_a = repmat( nanmean( Ldats{Lc}(:,1:mb), 'all') ,[size(Ldats{Lc},1) size(Ldats{Lc},2)]);
        bas_sd = repmat( nanstd( Ldats{Lc}(:,1:mb),[], 'all') ,[size(Ldats{Lc},1) size(Ldats{Lc},2)]);
        
        Ldats_N{Lc} = (Ldats{Lc}-bas_a)./bas_sd;
    end
    
    
    
    subplot(2,1,g)
    pls = plotstat(Ldats_N, pats,mb,ms,mp,pr,cols);
    ylim([-6 6])
    title(pr.compar_tags{g})
    
    legend(pls,laylabs(1:2),'AutoUpdate','off')
end

fnm = fullfile(figdir,['PATS_' pats{:} '_' pr.smoothwith '_' num2str(pr.smparam) '_' evtylab 'evs']);
if pr.intens_correct
    fnm = [fnm '_intenscorrect'];
end
if pr.maxintens
    fnm = [fnm  '_maxintens'];
end
fnm = [fnm '_STAT' pr.mcorr '_supra_infra_comp'];
saveas(fig,[fnm '.jpg'])
saveas(fig,[fnm '.fig'])
close(fig)

end



function laych_plot(z_grsM, pats,mb,ms,mp,layers,pr,figdir)


if ~iscell(pr.evtype)
    evtylab = pr.evtype;
else
    c = cellfun(@(x) [x '_'],pr.evtype,'UniformOutput',0);
    evtylab = [c{:}];
end;
laynr = 7;
% Supra/infragranular layers avg
z_avgs = cell(length(pats), 2);
for p = 1:length(pats)
    for k = 1:laynr
        if ~isempty(layers{p})
            chs = [layers{p}{k}];
        else
            continue;
        end
        
        z_avgs{p,k} =  cellfun(@(y) nanmean(y(chs,:,:),1),z_grsM{p},'UniformOutput',0);
    end
end
nempt = ~cellfun(@isempty, z_avgs);
z1c = cell(1,2);
cols = {'r','k'};
fig = figure;
for L = 1:laynr
    Ltit = ['Layer' num2str(L)];
    
    for k =1:2
        z1= cellfun(@(x) nanmean(x{k},3),z_avgs(nempt(:,L),L),'UniformOutput',0);
        z1c{k} = cat(1,z1{:});
    end
    datnr = size(z1c{1},1);
    subplot(laynr,1,L)
    plotstat(z1c, pats,mb,ms,mp,pr,cols);
    title([Ltit ', pat nr = ' num2str(datnr)]);
    
    
end

suptitle(['Patients: ' pats{:} ]);
fnm = fullfile(figdir,['PATS_' pats{:} '_' pr.smoothwith '_' num2str(pr.smparam) '_' evtylab 'evs']);
if pr.intens_correct
    fnm = [fnm '_intenscorrect'];
end
if pr.maxintens
    fnm = [fnm  '_maxintens'];
end
fnm = [fnm '_STAT' pr.mcorr '_layers'];
saveas(fig,[fnm '.jpg'])
saveas(fig,[fnm '.fig'])
close(fig)
end



function pL = plotstat(z1c, pats,mb,ms,mp,pr,cols)


t = 1:size(z1c{1},2);
for p = 1:2
    avg = nanmean(z1c{p},1);
    sd = nanstd( z1c{p},[],1) / sqrt( size(z1c{p},1) );
    
%     sd = std( z1c{p},[],1);
    sdup = avg + sd; sddo = avg-sd;
    pL(p) = plot(t, avg,'Color',cols{p}); hold on;
    patch('XData',[t fliplr(t)], 'YData',[sddo fliplr(sdup)],...
        'FaceAlpha',0.2,'FaceColor',cols{p},'EdgeColor','none');
    
datnr(p) = size(z1c{p},1);
end
title({['Patients: ' pats{:} ', pat nrs = ' num2str(datnr(1)) ' vs ' num2str(datnr(2)) ],'MUA avg. across all channels'});
% yL = ylim;
yL = [-2 2]; ylim(yL);
line([mb mb],yL,'Color','k','LineStyle','--')
line([mb+ms mb+ms],yL,'Color','k','LineStyle','--'); ylim(yL);
xlim(t([1 end]));

for j = 1:3
    switch j; case 1; tinx = 1:mb; case 2; tinx = mb+1:mb+ms; case 3; tinx = mb+ms+1: mb+ms+mp; end;
    % tinx = 1:mb+ms+mp;
    z_groups_bl = cellfun(@(x) x(:,tinx) , z1c,'UniformOutput',0);
    
    data1 = z_groups_bl{1};
    data2 = z_groups_bl{2};
    
    if sum(~isnan(data1(:,1)))>2||sum(~isnan(data2(:,1)))>2
        %     [permmaps,zmap_p,~] = my_permutationtest1(data1,data2,pr.itnum,0.05);
        [permmaps,zmap_p,~] = my_permutation_Ttest(data1,data2,pr.itnum,0.05);
        
        if strcmp(pr.mcorr,'none')
            sign_z = zmap_p;
        elseif strcmp(pr.mcorr,'cluster')
            sign_z = my_cluster_correction(permmaps,zmap_p,pr.itnum,0.05);
        end
        sign_p = sign_z~=0;
        draw_signifpatch(tinx,sign_p',[.7 .5 .5])
    end
    
end
set(gca,'TickDir','out','box','off');


end