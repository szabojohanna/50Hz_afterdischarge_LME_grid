
function mua_plots(rootdir_ad,pat,win,muafig,indivpeaks,pktype,timingfig,varargin)

% MUA plots
if ~isempty(varargin)
    eventinx = varargin{1};
else
    eventinx = [];
end



%% Plots MUA for one event of one patient (peaks & peak avgs)
if  ~iscell(pat) && ~isempty(eventinx)
    ti = varargin{2};
    chans = varargin{3};
%     maxchan = varargin{4};
    
    if muafig
        [th_avg, ~] = mua_bypeak(rootdir_ad,pat, eventinx,ti, chans,win,muafig,indivpeaks,pktype);
    end
    
    
elseif  ~iscell(pat) && isempty(eventinx)
    
    MUA_allpatients(rootdir_ad,{pat},win,indivpeaks,muafig,timingfig,pktype);
else
    
    %% MUA plots -  LOOP OVER PATIENTS
    
    MUA_allpatients(rootdir_ad,pat,win,indivpeaks,muafig,timingfig,pktype);
   
% 
% MUA_interstim
end
end



%--------------------------------------------------------------------------
function MUA_allpatients(rootdir_ad,patients,win,indivpeaks, muafig,timingfig,pktype)

for pp = 1:length(patients)
    patname = patients{pp};
    fprintf('%s...',patname);
    load(fullfile(rootdir_ad,'Patients',patname,'Events.mat'));
    ADevs = find([Events.AD]==1);
    for evi = 1:length(ADevs)
        eventinx = ADevs(evi);
        fprintf('%d...',eventinx);
        for ti = 1:2
            switch ti; case 1; chans = 1:23; case 2; chans = 25:47; end;
            
            if muafig
                mua_bypeak(rootdir_ad,patname, eventinx,ti, chans,win,true,indivpeaks,pktype);
            end
            if timingfig
                MUAmax_timing(rootdir_ad,patname,eventinx,ti, chans,win,pktype)
            end
            
        end
        
    end
    
end
end

%---------------------------------------------------------------------------
function [th_avg, th_pks] = mua_bypeak(rootdir_ad,patname, eventinx,ti, chans,win,isfig,indivpeaks,pktype)

load(fullfile(rootdir_ad,'channels_by_layers.mat'));
islaypat  = ismember(patname,channels_by_layers.Properties.VariableNames);

%% Load data in nswiew

[th_avg, th_pks] = deal([]);
adnr = 1;

%% Select and filter epochs

if ~contains(pktype,'Small')
    lf_freq = 30;
else
    lf_freq = 40;
end
    

[normepoch, sr] = load_filt_MUAepoch(rootdir_ad,patname,eventinx,adnr,win,pktype,lf_freq);
if isempty(normepoch); return; end;


avgfigdir = fullfile(rootdir_ad,'Figures','Micro','MUA','MUA_AD',['PkAVGs_W'  num2str(win)]);
indivfigdir = fullfile(rootdir_ad,'Figures','Micro','MUA','MUA_AD',['Indiv_pks_W'  num2str(win)],patname);

if ~isdir(avgfigdir); mkdir(avgfigdir); end;
if ~isdir(indivfigdir); mkdir(indivfigdir); end;


set(0, 'DefaultFigureRenderer', 'painters');


time = win(1):1/sr:win(2); time = time(1:size(normepoch,1));


% clim  = [-.4 .4];
clim = [-2 2];
pknr = size(normepoch,3);
pks = 1:pknr;
if strcmp(patname,'Pt8')&&eventinx==54;
    pks = 3:pknr;
elseif strcmp(patname,'Pt8')&&eventinx==103;
    pks = 2:pknr;
end

%% By peaks
if indivpeaks
    for k = 1:pknr
        
        th_avg = normepoch(:,chans,k);
        
        if isfig
            fig1 = figure; pcolor(time,1:length(chans),th_avg')
            shading interp
            set(gca,'YDir','reverse');
            colorbar; caxis(clim)
            if islaypat
                layers = channels_by_layers.(patname)(:,ti); 
                label_layers(layers);
            end
            
            
            
            title({['Ev' num2str(eventinx) ' MUA'], ['TH' num2str(ti) ', PK ' num2str(k)]})
            pause(1)
            
            fignm1 = ['Ev' num2str(eventinx) '_TH' num2str(ti) 'MUA_PK' num2str(k)];
            
            saveas(fig1,fullfile(indivfigdir,[fignm1 '.jpg']));
            saveas(fig1,fullfile(indivfigdir,[fignm1 '.fig']));
            close(fig1)
        end
        
    end
end
%% Average

% Average block of peaks, it there are many peaks
if pknr>30 && ~contains(pktype,'Small')
    pk_bl = 1:19:pknr;
    if pk_bl(end)~=pknr; pk_bl = [pk_bl pknr]; end;
else
    pk_bl = pks([1 end]);
end
clim = [-1 1];
basinx = 1:abs(win(1))*sr;

st = pks(1);
for k = 2:length(pk_bl)
    th_pks = normepoch(:,chans,st:pk_bl(k));
    th_avg = nanmean(th_pks,3);
    
    if isfig
        fig1 = figure; pcolor(time,1:length(chans),th_avg')
        shading interp
        set(gca,'YDir','reverse');
        colorbar; caxis(clim)
        if islaypat
            layers = channels_by_layers.(patname)(:,ti);
            label_layers(layers);
        end
        th_pks2 = permute(th_pks,[2 1 3]);
%         
%         % Bootstrap stat without correction
%         [exactp_ersp,maskersp,alphafdr] = boostat_eeglab_J(th_pks2,1:23,0.05,200,false,'none',[],basinx);
%         
%         hold on; contour(time,1:length(chans),maskersp,'Color',[0.5 0.5 0.5],'LineWidth',0.2);
        
        
        % Bootstrap stat with FDR correction
        [exactp_ersp,maskersp,alphafdr] = boostat_eeglab_J(th_pks2,1:23,0.05,200,false,'fdr',[],basinx);
        
        hold on; contour(time,1:length(chans),maskersp,'Color','white','LineWidth',0.7)
        
        
        
        
        title({['Ev' num2str(eventinx) ' MUA'], ['TH' num2str(ti) ', AVG (PK' num2str(st) '-PK' num2str(pk_bl(k)) ')']},'LineWidth',2)
%         pause(1)
        
        fignm1 = [patname '_Ev' num2str(eventinx) '_TH' num2str(ti) 'MUA_AVG_' num2str(k-1)];
        
        saveas(fig1,fullfile(avgfigdir,[fignm1 '.jpg']));
        saveas(fig1,fullfile(avgfigdir,[fignm1 '.fig']));
        close(fig1)
    end
    
    st = pk_bl(k) +1;
end


end
