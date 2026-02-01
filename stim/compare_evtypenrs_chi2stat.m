function compare_evtypenrs_chi2stat(rootdir_ad,patsss,thbsss,period,chlab)

THdist_labs = {'Close','Far'};
ADnoAD_labs = {'noAD, out of SOZ','noAD, in SOZ'}; % 'noAD, out of SOZ','noAD, in SOZ','AD'
MUAch_labs = {'MUA increase','MUA decrease'};

%%


nm = 'Event_summary_Zone_Thloc_MUA';
if strcmp(period,'stim')
    if isempty(chlab)
        exnm = nm ;
    else
        exnm = [nm '_' chlab ];
    end
elseif strcmp(period,'post')
    if isempty(chlab)
        exnm = [nm '_post'];
    else
        exnm = [nm '_post_' chlab];
    end
end

maxint= true;
%%
sheetnm = 'All pats';
barplot_chistat(rootdir_ad,exnm,sheetnm,THdist_labs,MUAch_labs,{'noAD'},maxint)

%%
barplot_chistat(rootdir_ad,exnm,sheetnm,{'AD','noAD, in SOZ','noAD, out of SOZ'},MUAch_labs,{'Close','Far'},maxint)

%%


for p = 1:length(patsss)
    patname = patsss{p};
    ti = thbsss(p);
    sheetnm = [patname ' Th' num2str(ti)];
    barplot_chistat(rootdir_ad,exnm,sheetnm,{'Close','Far'},MUAch_labs,{})
end

end



%---------------------------------------------------------------------------
function barplot_chistat(rootdir_ad,exnm,sheetnm,group2comp1,group2comp2,group2split,maxint)

exnm =  fullfile(rootdir_ad,[exnm '.xlsx']);
tab = readtable(exnm, 'Sheet',sheetnm);

%%
fieldnm1 = sel_fieldnm(group2comp1);
fieldnm2 = sel_fieldnm(group2comp2);
if ~isempty(group2split)
    fieldnm3 = sel_fieldnm(group2split);
else
    fieldnm3 = '';
end


% All patients, all zones max int. events: AD groups: MUA increase vs decrease
gr1_inx = cellfun(@(x) find(ismember(tab.(fieldnm1),x)), group2comp1,'UniformOutput',0);
gr2_inx = cellfun(@(x) find(ismember(tab.(fieldnm2),x)), group2comp2,'UniformOutput',0);


conttab = create_conttab(gr1_inx,gr2_inx,group2comp1,group2comp2,tab,maxint);

%

%%
[chi2_val, pval] = my_chi2_test(conttab);

% Bar plot
barplot(conttab,group2comp1,group2comp2,pval,chi2_val);
title([sheetnm ',all events'])
%%
if ~isempty(group2split)
    locgr_inx = cellfun(@(x) find(ismember(tab.(fieldnm3),x)), group2split,'UniformOutput',0);
    
    for k = 1:length(group2split)
        gr1_inx2 = cellfun(@(x) intersect(locgr_inx{k},x) , gr1_inx,'UniformOutput',0);
        gr2_inx2 = cellfun(@(x) intersect(locgr_inx{k},x) , gr2_inx,'UniformOutput',0);
        
        conttab = create_conttab(gr1_inx2,gr2_inx2,group2comp1,group2comp2,tab,maxint);
        
        [chi2_val, pval] = my_chi2_test(conttab);
        barplot(conttab,group2comp1,group2comp2,pval,chi2_val);
        title([sheetnm ', ' group2split{k} ' events'])
        
    end
end
end

%---------------------------------------------------------------------------
function conttab = create_conttab(gr1_inx,gr2_inx,ADnoAD_labs,MUAch_labs,tab,maxint)

ad_grnr = length(ADnoAD_labs);
m_grnr =  length(MUAch_labs);
conttab = nan(m_grnr,ad_grnr);
for a = 1:ad_grnr
    for m = 1:m_grnr
        r = intersect(gr1_inx{a},gr2_inx{m});
        if maxint
            evnr = sum(tab.('MaxInt_EventNr_')(r));
        else
            evnr = sum(tab.('EventNr_')(r));
        end
        conttab(m,a) = evnr;
    end
end

end

%---------------------------------------------------------------------------
function barplot(conttab,ADnoAD_labs,MUAch_labs,pval,chi2)

fig = figure;
% bh = bar(conttab','stacked');
bh = stacked_bar_perc(conttab',[1 0 0; 0 0 1],'k')
xticks([1 2]); xticklabels(ADnoAD_labs)

% set(bh, 'FaceColor', 'Flat');
% bh(1).CData = [1 0 0];
% bh(2).CData = [0 0 1];
legend(bh,MUAch_labs);
yL = ylim;
if pval<0.05; col = 'r'; else col = 'k'; end;
text(0,yL(2)*.9,['p = ' num2str(pval) ', chi2 = ' num2str(chi2)],'Color',col);
end

function fieldnm = sel_fieldnm(group2comp1)

if contains(group2comp1,'AD'); fieldnm = 'AD_SOZ';
elseif contains(group2comp1,'MUA'); fieldnm = 'MUAChange';
else fieldnm = 'THLoc'; end

end