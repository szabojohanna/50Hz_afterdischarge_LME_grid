function [mean2bas, win2bas,stim_tm, allwin,win2basS] = calc_muamean(rootdir_ad,patname,evinx,chans,period,timewin)

[mean2bas, win2bas,stim_tm] = deal([]);
mua_blockdir = fullfile(rootdir_ad,'Figures','Micro','MUA','MUA_blocks');

 mbfnm = dir(fullfile(mua_blockdir,['MUA_stimblocks_' patname '_Ev' num2str(evinx) '_*.mat']));
    
    if isempty(mbfnm)
        fprintf('No MUA blocks for Ev%d\n', evinx); return;
    end
    
    load(fullfile(mua_blockdir,mbfnm.name))
    
    block_MUA = MUA_stimblocks.stimblock_MUA;
    basblock_MUA = MUA_stimblocks.basblock_MUA;
    postblock_MUA = MUA_stimblocks.postblock_MUA;
    
    basend =size(basblock_MUA,3);
    stimend = basend+ size(block_MUA,3);
    
    
    basbla = squeeze(mean(basblock_MUA(:,chans,:), 1) );
    
    stimbla = squeeze( mean(block_MUA(:,chans,:), 1) );
    postbla = squeeze( mean(postblock_MUA(:,chans,:), 1) );
    
    tm0 = cat(2,basbla,stimbla);
    tm = cat(2,tm0,postbla);
    
    if ~iscell(period)
        period = {period};
    end
    if ~iscell(timewin)
        timewin = {timewin};
    end
    
    for j = 1:1+length(period)
        switch j;
            case 1; per = 'bas'; tw = [];
            otherwise; per = period{j-1};
                tw = timewin{j-1};
        end;
        
        if strcmp(per,'bas')
            stim_lims1 = MUA_stimblocks.params.Bblocks_lims(:,1);
            periodst = 1;
            perioden = basend-20;
        elseif strcmp(per,'stim')
            stim_lims1 = MUA_stimblocks.params.blocks_lims(:,1);
            periodst = basend+1;
            perioden = stimend;
        elseif strcmp(per,'post')
            stim_lims1 = MUA_stimblocks.params.Pblocks_lims(:,1);
            periodst = stimend+1;
            perioden = periodst+size(postblock_MUA,3)-1;
        end
        
        stim_lims1_c = stim_lims1 - stim_lims1(1);
        
        if ~isempty(tw)&& length(tw)==2
            if min(tw)>max(stim_lims1_c)
                fprintf('Ev%d too short\n',evinx)
                return
            end
            win_lims = dsearchn(stim_lims1_c,tw')-1;
            stim_tm = tm(:, periodst+win_lims(1):periodst+win_lims(2) );
        elseif ~isempty(tw)&& length(tw)==1 % for slide
            stim_tm = tm(:,periodst+tw);
        else
            stim_tm = tm(:,periodst:perioden);
        end
        
        
        if strcmp(per,'bas')
            bas_tm = stim_tm;
        else strcmp(per,'stim')
            per_tmS{j-1} = stim_tm;
        end
        
        
    end
    basm = mean(bas_tm,2); bassd = std(bas_tm,[],2);
    bas2bas = (bas_tm - repmat(basm,[1 size(bas_tm,2)]) )./repmat(bassd,[1 size(bas_tm,2)]);
    
    for h = 1:length(period)
        win2basS{h} = (per_tmS{h} - repmat(basm,[1 size(per_tmS{h},2)]) )./repmat(bassd,[1 size(per_tmS{h},2)]);
    end
    win2bas = [win2basS{:}];
     mean2bas = mean(win2bas,'all');
    allwin= cat(2,bas2bas, win2bas);
end