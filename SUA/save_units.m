function save_units(rootdir_ad,patname,choi,sr)
% Saves table with name and quality measures of units

unit_names = {}; uq = {}; cr = {}; 

if exist(fullfile(rootdir_ad,'SUA_Analysis','SUA_data'))==2
    load(fullfile(rootdir_ad,'SUA_Analysis','SUA_data'))
    n = size(unit_table,1)+1;
else
    unit_table = table;
    n = 1;
end
for k = 1:length(choi)
    chan = choi(k);
    savedir = fullfile(rootdir_ad,'SUA_Analysis',patname,['Chan' num2str(chan)]);
    fnm = ['Chan' num2str(chan) '_AllEvents_Raw'];
    % load(fullfile(savedir,[fnm '_spikes.mat']));
    load(fullfile(savedir,['times_' fnm '.mat']))
    cluss = unique(cluster_class(cluster_class(:,1)~=0,1)); % unit clusters
    
    
    [metrics, outputWaves, clusterIDs,unitQuality, contaminationRate] = spik_waveform_metrics(rootdir_ad,patname,chan,1,sr);

    
    for j = cluss'
        unit_name = [patname '_' num2str(chan) '_' num2str(j)];
        unit_names{n} = unit_name;
        uq{n} = unitQuality(clusterIDs==j);
        cr{n} = contaminationRate(clusterIDs==j);
        n = n+1;
    end
    
end

unit_quality = [uq{:}]';
contamination_rate = [cr{:}]';

unit_table0 = table(unit_names',unit_quality',contamination_rate');
unit_table0.Properties.VariableNames = {'Unit_name'  'Unit_quality'  'Contamination_rate'};
unit_table = [unit_table; unit_table0];

save(fullfile(rootdir_ad,'SUA_Analysis','SUA_data'),'unit_table');

