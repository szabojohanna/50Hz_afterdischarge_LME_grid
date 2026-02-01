function fig = csd_simple(csdm,win_ts,clim,ti)
% csdm: time x channel
% win_ts: row vector
nbchan = size(csdm,2);
fig = figure;
% set(fig,'Visible','off');
% yyaxis left

try
    pcolor(win_ts,1:nbchan,csdm'); shading interp; colorbar
catch
    keyboard;
end
set(gca,'ydir','reverse')
if isempty(clim)
    clim = prctile(csdm(:),[0 100]);
end
set(gca,'clim',clim)
colormap(flipud(jet))
if ~isempty(ti)
    ylabel(['Channels (Th' num2str(ti) ')'])
end