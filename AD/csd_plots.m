function csd_plots(rootdir_ad,pat,win,freqs,isfig,indivpeaks,pktype,varargin)

% CSD PLOTS
% Prepare matrices for hemming filter and csd (.ldr files) for the patient manually!

if ~isempty(varargin)
    eventinx = varargin{1};
else
    eventinx = [];
end


if ~iscell(pat) && ~isempty(eventinx)
    %% Plot CSD - one event of one patient (peaks & peak avgs)
    
    ti = varargin{2};
    chans = varargin{3};
    %%
    csd_filtered_bypeak(rootdir_ad,pat,eventinx,ti,chans,win,freqs,isfig,indivpeaks,pktype)
    
elseif  ~iscell(pat) && isempty(eventinx)
    
    CSD_allpatients2(rootdir_ad,{pat},win,freqs,indivpeaks,pktype,isfig)
else
    %% Plot CSD -all patients
    CSD_allpatients2(rootdir_ad,pat,win,freqs,indivpeaks,pktype,isfig)
    %% Plot CSD - all patients
end


end

%--------------------------------------------------------------------------
function CSD_allpatients2(rootdir_ad,patients,win,freqs,indivpeaks,pktype,isfig)

for pp = 1:length(patients)
    patname = patients{pp};
    fprintf('%s...',patname);
    load(fullfile(rootdir_ad,'Patients',patname,'Events.mat'));
    ADevs = find([Events.AD]==1);
    
    if isfield(Events,'Th2_STIMdist_cm'); thnr = 2; else; thnr = 1; end;

        
    for evi = 1:length(ADevs)
        eventinx = ADevs(evi);
        fprintf('Ev%d...\n', eventinx);
        for ti = 1:thnr
            switch ti; case 1; chans = 1:23; case 2; chans = 25:47; end;
            csd_filtered_bypeak(rootdir_ad,patname,eventinx,ti,chans,win,freqs,isfig,indivpeaks,pktype)
        end
        
    end
    
end
end


function [csdavg] = csd_filtered_bypeak(rootdir_ad,patname,eventinx,ti,chans,win,freqs,isfig,indivpeaks,pktype)

csdavg = [];
load(fullfile(rootdir_ad,'channels_by_layers.mat'));
islaypat  = ismember(patname,channels_by_layers.Properties.VariableNames);

avgfigdir = fullfile(rootdir_ad,'Figures','Micro','CSD',['PkAVGs_W'  num2str(win) '_' pktype]);
indivfigdir = fullfile(rootdir_ad,'Figures','Micro','CSD',['Indiv_pks_W'  num2str(win) '_' pktype],patname);

if ~isdir(avgfigdir); mkdir(avgfigdir); end;
if ~isdir(indivfigdir); mkdir(indivfigdir); end;



adnr = 1;
set(0, 'DefaultFigureRenderer', 'painters');
%% Load data in nswiew

[normepoch, sr] = load_filt_epoch(rootdir_ad,patname,eventinx,adnr,win,freqs,pktype);
if isempty(normepoch); return; end;

% CSD

time = win(1):1/sr:win(2); time = time(1:size(normepoch,1));

patdir = fullfile(rootdir_ad,'Patients',patname);
cd(patdir)

ldrnm_csd = dir(fullfile(patdir,'*csd.ldr'));
m_csd = load(fullfile(patdir,ldrnm_csd.name));
ldrnm_hem = dir(fullfile(patdir,['*hem_' patname '.ldr']));
m_hem = load(fullfile(patdir,ldrnm_hem.name));




%% By peaks
pknr = size(normepoch,3);
pks = 1:pknr;
if strcmp(patname,'Pt8')&&eventinx==54;
    pks = 3:pknr;
elseif strcmp(patname,'Pt8')&&eventinx==103;
    pks = 2:pknr;
end
inpknr = floor(linspace(pks(1),pknr,min(pknr,5)));

if indivpeaks
    clim = [-2 2];
    for k = inpknr
        
        thdat_hem = normepoch(:,:,k)*m_hem';
        
        csdm = thdat_hem*m_csd';
        csdm = csdm(:,chans);
        
        if isfig
            fig1 = csd_simple(csdm,time,clim,ti);
            hold on;
            if islaypat
                layers = channels_by_layers.(patname)(:,ti); 
                label_layers(layers);
            end
            
            title({['Ev' num2str(eventinx) ' CDS'], ['TH' num2str(ti) ', PK ' num2str(k)]})
            pause(1)
            
            fignm1 = ['Ev' num2str(eventinx) '_TH' num2str(ti) 'CSD_PK' num2str(k)];
            
          
            saveas(fig1,fullfile(indivfigdir,[fignm1 '.jpg']));
            saveas(fig1,fullfile(indivfigdir,[fignm1 '.fig']));
            close(fig1)
        end
        
    end
end
%% Average

clim = [-1 1];
% clim = [-0.3 0.3];


% Average block of peaks, it there are many peaks

if pknr>30 && ~contains(pktype,'Small')
    pk_bl = 1:19:pknr;
    if pk_bl(end)~=pknr; pk_bl = [pk_bl pknr]; end;
else
    pk_bl = pks([1 end]);
end

basinx = 1:abs(win(1))*sr;
st = pks(1);
for k = 2:length(pk_bl)
    
    th_pks = normepoch(:,:,st:pk_bl(k));
    
    thdat_hem = nanmean(th_pks,3)*m_hem';
    csdavg = thdat_hem*m_csd';
    csdavg = csdavg(:,chans);
    
    if isfig
        fig1 = csd_simple(csdavg,time,clim,ti);
        hold on;
        if islaypat
            layers = channels_by_layers.(patname)(:,ti);
            label_layers(layers);
        end
        % CSD by peaks for stat
        th_pks2 = nan(size(csdavg,2),size(csdavg,1),size(th_pks,3));
        for j = 1:size(th_pks,3)
            thpk_hem = th_pks(:,:,j)*m_hem';
            
            csdpk = thpk_hem*m_csd';
            th_pks2(:,:,j) = csdpk(:,chans)';
        end
        
%         % Bootstrap stat without correction
%         [exactp_ersp,maskersp,alphafdr] = boostat_eeglab_J(th_pks2,1:23,0.05,200,false,'none',[],basinx);
%         
%         hold on; contour(time,1:length(chans),maskersp,'Color',[1 1 1],'LineWidth',0.1);
        
        
        % Bootstrap stat with FDR correction
        [exactp_ersp,maskersp,alphafdr] = boostat_eeglab_J(th_pks2,1:23,0.05,200,false,'fdr',[],basinx);
        
        hold on; contour(time,1:length(chans),maskersp,'Color','white','LineWidth',0.7)
        
        
        
        
        title({['Ev' num2str(eventinx) ' CDS'], ['TH' num2str(ti) ', AVG (PK' num2str(st) '-PK' num2str(pk_bl(k)) ')']})
        pause(1)
        
        fignm1 = [patname '_Ev' num2str(eventinx) '_TH' num2str(ti) 'CSD_AVG_' num2str(k-1)];
        
        saveas(fig1,fullfile(avgfigdir,[fignm1 '.jpg']));
        saveas(fig1,fullfile(avgfigdir,[fignm1 '.fig']));
        close(fig1)
    end
    st = pk_bl(k)+1;
end

end