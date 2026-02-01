function [block_MUA, block_MUA_norm] = filtMUAblocks(blockdat,srate,block_ms,mua_ms,avgrep,sdrep)


tms = block_ms(1):1000/srate:block_ms(2); tms = tms(1:end-1); % window cut around stim. peak
tmsmu = mua_ms(1):1000/srate:mua_ms(2); % window around stim. peak to use for MUA calculation
muinx = dsearchn(tms',tmsmu');

% Filter for MUA
[b, a] = butter(2,500/(srate/2),'high'); % high-pass
lpF = 300;
[bb, aa] = butter(2,lpF/(srate/2),'low'); % low-pass for smoothing

% timeL = size(blockdat,1);
% halfL = round(timeL/2);

% blockdat_pad = cat(1, zeros(size(blockdat)), blockdat,   zeros(size(blockdat)) );
[block_MUA, block_MUA_norm] = deal(nan(length(muinx),size(blockdat,2)));
for k = 1:size(blockdat,3); % loop over blocks
    
    fdata = filtfilt(b,a,blockdat(:,:,k));
    rectdat= abs(fdata);
    %     eMUA = envelope(rectdat);
    eMUA_pad = filtfilt(bb,aa,envelope(rectdat));
    eMUA = eMUA_pad(muinx,:,:);
    %     eMUA = eMUA_pad(halfL+1:end-halfL,:);
    
    
    %
    
    if isempty(sdrep) &&  isempty(avgrep)
        
        basMUAavg = mean(eMUA,1);
        basMUAsd = std(eMUA,[],1);
        avgrep = repmat(basMUAavg, [size(eMUA,1) 1]);
        sdrep = repmat(basMUAsd, [size(eMUA,1) 1]);
    end
    
    eMUA_norm = (eMUA - avgrep)./sdrep;
    block_MUA(:,:,k) = eMUA;
    block_MUA_norm(:,:,k) = eMUA_norm;
    %
    % %     Check blocks
    %     maxch = 32;
    %     figure;
    %
    %     subplot(311);plot(tms,fdata(:,maxch)); title('High-pass filtered 500 Hz'); xlim(block_ms)
    %
    %     subplot(312); plot(tms,rectdat(:,maxch));xlim(block_ms)
    %     hold on; plot(tms,envelope(rectdat(:,maxch)));  title('Abs + Envelope')
    %     subplot(313);  plot(tmsmu, eMUA(:,maxch));  title(['Low-pass filtered ' lpF ' Hz']); xlim(block_ms)
    %     suptitle(['Channel: Ch' num2str(maxch)])
    %     chansi = 25:47;
    %     figure
    %     pcolor(tmsmu,chansi,zscore(eMUA(:,chansi))'); shading interp; set(gca,'YDir','reverse'); colorbar;
    %     colormap(pink); xlabel('Time rel. to stim. peak (ms)')
    %     caxis([-3 3])
    %     pause(1); close(gcf)
    %
end
end
