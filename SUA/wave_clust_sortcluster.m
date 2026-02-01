function wave_clust_sortcluster(rootdir_ad,patname,chan)


% dirnm = dir(fullfile(suadir,[patname '_Ev*_Chan' num2str(chan)]));
% savedir = fullfile(suadir,dirnm(1).name);
% matfnm =  dir(fullfile(savedir,[patname '_Ev*_Chan' num2str(chan) '.mat']));
% fnm2 = [matfnm.name(1:end-4) ];
savedir = fullfile(rootdir_ad,'SUA_Analysis',patname,['Chan' num2str(chan)]);
fnm2 = ['Chan' num2str(chan) '_AllEvents_Raw'];

cd(savedir);
par = set_parameters;
par.sr = 20000; par.detection = 'both'; par.stdmin = 3.5; par.segments_length = 5; par.features = 'pca';
Get_spikes([fnm2 '.mat'],'par',par);
load([fnm2 '_spikes.mat'])
Do_clustering([fnm2 '_spikes.mat']);

% wave_clus([fnm2 '.mat'],'par',par)
% wave_clus
% keyboard;
% ok = input('Okay clusters? 1/0');

% if ok
%     save_clusters(rootdir_ad,savedir,patname,chan)
% end
end