function raw_csd_mua_layers(rootdir_ad,patname,eventinx,ti, chans,win,freqs,pktype,lf_freq)



[raw, csd, mua, Rtime, Mtime,Rsr] = load_dats(rootdir_ad,patname,ti,chans,eventinx,win,freqs,pktype,lf_freq);

% suprainfra_plot_same(rootdir_ad,patname,eventinx,ti,raw,csd,mua, Rtime, Mtime)

PL_win = 0.1;
%% First (positive) phase of spike
% r_layer = 3;
% sm_win = [-40 -20]; % ms
% 
% tit = 'P1';
% rpol = 1;
% cpol = 1;
% difflayers_plot(rootdir_ad,patname,eventinx,ti,raw,csd,mua, Rtime, Mtime,r_layer,sm_win,PL_win,Rsr,tit,rpol, cpol)


%% Initial negative spike
r_layer = 4;
sm_win = [-70 70];

tit = 'N1';
rpol = 0;
cpol = 1;
difflayers_plot(rootdir_ad,patname,eventinx,ti,raw,csd,mua, Rtime, Mtime,r_layer,sm_win,PL_win,Rsr,tit,rpol, cpol)


% Deep positive wave
% r_layer = 5;
% sm_win = [-30 120];
% 
% tit = 'P2';
% rpol = 1;
% cpol = 0;
% difflayers_plot(rootdir_ad,patname,eventinx,ti,raw,csd,mua, Rtime, Mtime,r_layer,sm_win,PL_win,Rsr,tit,rpol, cpol)



%% Superficial negative
r_layer = 2;
sm_win = [20 150];

tit = 'N2';
rpol = 0;
cpol = 1;
difflayers_plot(rootdir_ad,patname,eventinx,ti,raw,csd,mua, Rtime, Mtime,r_layer,sm_win,PL_win,Rsr,tit,rpol, cpol)


%% Deep negative
r_layer = 6;
sm_win = [50 220];

tit = 'N3';
rpol = 0;
cpol = 1;
difflayers_plot(rootdir_ad,patname,eventinx,ti,raw,csd,mua, Rtime, Mtime,r_layer,sm_win,PL_win,Rsr,tit,rpol, cpol)


%% Superficial positive
r_layer = 2;
sm_win = [60 260];

tit = 'P3';
rpol = 1;
cpol = 0;
difflayers_plot(rootdir_ad,patname,eventinx,ti,raw,csd,mua, Rtime, Mtime,r_layer,sm_win,PL_win,Rsr,tit,rpol, cpol)
% 
% 
% % Superficial negative
% r_layer = 2;
% sm_win = [200 500];
% 
% tit = 'N4';
% rpol = 0;
% cpol = 1;
% difflayers_plot(rootdir_ad,patname,eventinx,ti,raw,csd,mua, Rtime, Mtime,r_layer,sm_win,PL_win,Rsr,tit,rpol, cpol)



end

function difflayers_plot(rootdir_ad,patname,eventinx,ti,raw,csd,mua, Rtime, Mtime,r_layer,sm_win,PL_win,sr,tit,rpol, cpol)

load(fullfile(rootdir_ad,'channels_by_layers.mat'))
layers = channels_by_layers.(patname)(:,ti);

figdir = fullfile(rootdir_ad,'Figures','Micro','Raw_CSD_MUA');
if ~isdir(figdir); mkdir(figdir); end;

pknr = size(raw,3);

% Find LFP peak
Rchs = [layers{r_layer}]; % Channels of selected LFP layer
raw_L = mean( mean( raw(:,Rchs,:) ,2) ,3); % Average across channels within layer

raw_L_sd = std( mean( raw(:,Rchs,:) ,2) ,[], 3)./ pknr;
raw_L_sdup = raw_L_sd+raw_L;
raw_L_sddo = raw_L - raw_L_sd;

rlims = dsearchn(Rtime',sm_win'); % Time limits of window to find LFP peak
rinx = rlims(1):rlims(2); % Time indeces



% rpol = input('Positive (1) or negative (0)?');

fig = figure;

ok = 0;
while ~ok
    
    if rpol
        findpeaks(raw_L(rinx),'Npeaks',1,'SortStr','descend'); title(['small win: ' num2str(sm_win) ' s'])
        [pks, Rlocs] = findpeaks(raw_L(rinx),'Npeaks',1,'SortStr','descend');
    else
        findpeaks(-raw_L(rinx),'Npeaks',1,'SortStr','descend'); title(['small win: ' num2str(sm_win) ' s'])
        [pks, Rlocs] = findpeaks(-raw_L(rinx),'Npeaks',1,'SortStr','descend');
    end
    ok = input('OK? 0/1');
    if ok
        close(fig);
    else
        sm_win = input('new small win: ');
        rlims = dsearchn(Rtime',sm_win');
        rinx = rlims(1):rlims(2);
    end
end

pk = Rlocs+rinx(1);
PL_dp = PL_win*sr;
pk_inx = pk-PL_dp:min(pk+PL_dp,length(Rtime));
PL_Rtime = Rtime(pk_inx);

if rpol
    [~, R_mch] = max( mean( raw(Rlocs,Rchs,:),3 ) );
else
    [~, R_mch] = max( -mean( raw(Rlocs,Rchs,:),3 ) );
end
Rchinx = Rchs(1) + R_mch-1;


% Find CSD peak

% cpol = input('Sink (1) or source (0)?');

csd_avg = mean(csd,3);

fig = figure;

subplot(121)
if cpol
    dat = -csd_avg;
else
    dat = csd_avg;
end

findpeaks(dat(pk,:),'Npeaks',1,'SortStr','descend');
[~,csd_chan] = findpeaks(dat(pk,:),'Npeaks',1,'SortStr','descend');

ok = 0;
while ~ok
    
    
    subplot(122)
    
    [~, csdmax] = findpeaks(dat(rinx,csd_chan),'Npeaks',1,'SortStr','descend');
    
    imagesc( csd_avg(rinx,:)' ); colormap(flipud(jet)); hold on;
    scatter(Rlocs,csd_chan,[],'w','filled');
    scatter(csdmax,csd_chan,50,'w','+');
    
    
    
    ok = input('OK? 0/1');
    if ok
        close(fig);
    else
        csd_chan = input('peak channel manually: ');
    end
end

c_layer = find( cellfun(@(x) ismember(csd_chan,x), layers) );
Cchs = [layers{c_layer}];
csd_L = mean( mean( csd(:,Cchs,:) ,2) ,3);
csd_L_sd = std( mean( csd(:,Cchs,:) ,2) ,[], 3)./pknr;
csd_L_sdup = csd_L_sd + csd_L;
csd_L_sddo = csd_L - csd_L_sd;

if cpol
    [~, csdmax] = findpeaks(-csd_L(rinx),'Npeaks',1,'SortStr','descend');
else
    [~, csdmax] = findpeaks(csd_L(rinx),'Npeaks',1,'SortStr','descend');
end
csdM = csdmax + rinx(1);

mua_L = mean( mean( mua(:,Cchs,:) ,2) ,3);
mua_L_sd = std( mean( mua(:,Cchs,:) ,2) ,[], 3)./pknr;
mua_L_sdup = mua_L_sd + mua_L;
mua_L_sddo = mua_L - mua_L_sd;

mx = dsearchn(Mtime',PL_Rtime');
PL_Mtime = Mtime(mx);


xti = PL_Rtime([1 round(length(PL_Rtime)/2) end]);
% yLL = [-5 5];
% yLR = [-1.5 1.5];

% Max chan figure
fig = figure;
yyaxis left
pl(1) = plot(PL_Rtime,mean( raw(pk_inx,Rchinx,:) ,3),'Color','k','LineStyle','-'); hold on; % LFP
% ylim(yLL);
yLL = ylim;
line([PL_Rtime(PL_dp) PL_Rtime(PL_dp)],yLL,'Color','k','LineStyle','-') % Line for LFP peak


yyaxis right
pl(2) = plot(PL_Rtime,mean( csd(pk_inx,csd_chan,:), 3),'Color','b','LineStyle','-'); hold on; % CSD
% ylim(yLR);
yLR = ylim;
line([Rtime(csdM) Rtime(csdM)],yLR,'Color','b','LineStyle','--') % Line for CSD peak

set(gca, 'XColor','k', 'YColor','k')

yyaxis left
pl(3) = plot(PL_Mtime,mean( mua(mx,csd_chan,:), 3),'Color','r','LineStyle','-'); hold on; % MUA
ylim(yLL);

set(gca, 'XColor','k', 'YColor','k')

xlim(PL_Rtime([1 end]))

xlabel('Time (s)');
ylabel('Norm. unit')


legend(pl,{['LFP: channel ' num2str(Rchinx)],...
    ['CSD: channel ' num2str(csd_chan)],...
    ['MUA: channel ' num2str(csd_chan)]});

xticks(xti); xticklabels(arrayfun(@num2str, xti, 'UniformOutput',0))

title([tit ', ' patname ', Ev' num2str(eventinx)]);
set(gca,'TickDir','out','box','off');
fnm = fullfile(figdir,[patname '_' tit '_Ev' num2str(eventinx) '_R' num2str(r_layer) '_C' num2str(c_layer) '_M' num2str(c_layer) '_MAXCHAN'] );

saveas(fig,[fnm '.jpg'])
saveas(fig,[fnm '.fig'])
close(fig);




% Layers figure
% yLL = [-4 4];
fig = figure;
yyaxis left
pl(1) = plot(PL_Rtime,raw_L(pk_inx),'Color','k','LineStyle','-'); hold on; % LFP
patch('XData',[PL_Rtime fliplr(PL_Rtime)], 'YData',[raw_L_sdup(pk_inx)' fliplr(raw_L_sddo(pk_inx)')],...
    'FaceAlpha',0.2,'FaceColor','k','EdgeColor','none');
% ylim(yLL); 
yLL = ylim;
line([PL_Rtime(PL_dp) PL_Rtime(PL_dp)],yLL,'Color','k','LineStyle','-') % Line for LFP peak


yyaxis right
pl(2) = plot(PL_Rtime,csd_L(pk_inx),'Color','b','LineStyle','-'); hold on; % CSD

patch('XData',[PL_Rtime fliplr(PL_Rtime)], 'YData',[csd_L_sdup(pk_inx)' fliplr(csd_L_sddo(pk_inx)')],...
    'FaceAlpha',0.2,'FaceColor','b','EdgeColor','none');
yLR = ylim;
% ylim(yLR); 
line([Rtime(csdM) Rtime(csdM)],yLR,'Color','b','LineStyle','--') % Line for CSD peak


yyaxis left
pl(3) = plot(PL_Mtime,mua_L(mx),'Color','r','LineStyle','-'); hold on; % MUA

patch('XData',[PL_Mtime fliplr(PL_Mtime)], 'YData',[mua_L_sdup(mx)' fliplr(mua_L_sddo(mx)')],...
    'FaceAlpha',0.2,'FaceColor','r','EdgeColor','none');
ylim(yLL);


xlim(PL_Rtime([1 end]))

xlabel('Time (s)');
ylabel('Norm. unit')
xticks(xti); xticklabels(arrayfun(@num2str, xti, 'UniformOutput',0))


legend(pl,{['LFP: layer ' num2str(r_layer)],...
    ['CSD: layer ' num2str(c_layer)],...
    ['MUA: layer ' num2str(c_layer)]});


title([tit ', ' patname ', Ev' num2str(eventinx)]);

set(gca,'TickDir','out','box','off');
fnm = fullfile(figdir,[patname '_' tit '_Ev' num2str(eventinx) '_R' num2str(r_layer) '_C' num2str(c_layer) '_M' num2str(c_layer) '_LAYERS'] );

saveas(fig,[fnm '.jpg'])
saveas(fig,[fnm '.fig'])
saveas(fig,[fnm '.pdf'])
close(fig);

end


function suprainfra_plot_same(rootdir_ad,patname,eventinx,ti,raw,csd,mua, Rtime, Mtime)


load(fullfile(rootdir_ad,'channels_by_layers.mat'))
layers = channels_by_layers.(patname)(:,ti);



figdir = fullfile(rootdir_ad,'Figures','Micro','Raw_CSD_MUA');
if ~isdir(figdir); mkdir(figdir); end;


% Supra-infra plot
laynr = 3;
for k = 1:laynr
    
    if ~isempty(layers)
        switch k;
            case 1
                chs = [layers{1:3}];
                tit = 'Supra';
            case 2
                chs = [layers{5:6}];
                tit = 'Infra';
            case 3
                chs = [layers{4}];
                tit = 'Gran';
        end
    else
        if strcmp(patname,'Pt20')
            switch k;
                case 1
                    chs = 1:8;
                case 2
                    chs = 15:23;
                case 3
                    continue;
            end
        end
    end
    
    raw_L = mean( raw(:,chs) ,2);
    csd_L = mean( csd(:,chs) ,2);
    mua_L = mean( mua(:,chs) ,2);
    
    
    
    fig = figure;
    plot(Rtime,raw_L,'Color','k'); hold on;
    plot(Rtime,csd_L,'Color','g'); hold on;
    plot(Mtime,mua_L,'Color','r'); hold on;
    
    xlabel('Time (s)');
    ylabel('Norm. unit')
    legend({'LFP','CSD','MUA'})
    title([tit ', ' patname ', Ev' num2str(eventinx)]);
    fnm = fullfile(figdir,[patname '_Ev' num2str(eventinx) '_' tit] );
    
    saveas(fig,[fnm '.jpg'])
    saveas(fig,[fnm '.fig'])
    close(fig);
end


end



function [raw, csd, mua, Rtime, Mtime,Rsr] = load_dats(rootdir_ad,patname,ti,chans,eventinx,win,freqs,pktype,lf_freq);

adnr=1;
patdir = fullfile(rootdir_ad,'Patients',patname);
cd(patdir)

ldrnm_csd = dir(fullfile(patdir,'*csd.ldr'));
m_csd = load(fullfile(patdir,ldrnm_csd.name));
ldrnm_hem = dir(fullfile(patdir,['*hem_' patname '.ldr']));
m_hem = load(fullfile(patdir,ldrnm_hem.name));



[normepoch, Rsr] = load_filt_epoch(rootdir_ad,patname,eventinx,1,win,freqs,pktype);

% % Raw
% epochmean= mean(normepoch,3);
% raw = epochmean(:,chans,:);
raw = normepoch(:,chans,:);

% CSD
csdm = nan(size(normepoch));
for k = 1:size(normepoch,3)
    thdat_hem = normepoch(:,:,k)*m_hem';
    csdm(:,:,k) = thdat_hem*m_csd';
end
csd = csdm(:,chans,:);


Rtime = win(1):1/Rsr:win(2); Rtime = Rtime(1:size(normepoch,1))*1000;

% MUA
[normepoch_mua, Msr] = load_filt_MUAepoch(rootdir_ad,patname,eventinx,adnr,win,pktype,lf_freq);


Mtime = win(1):1/Msr:win(2); Mtime = Mtime(1:size(normepoch_mua,1))*1000;



mua = normepoch_mua(:,chans,:);
end