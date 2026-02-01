function [basend_s, stimend_s, postend_s] = find_evlimits_suamua(rootdir_ad,patname,eventinx,sr,EventInfo)

matdir =  fullfile(rootdir_ad,'Figures','Micro','MUA','MUA_blocks');
matnm = dir(fullfile(matdir,['MUA_stimblocks_' patname '_Ev' num2str(eventinx) '_AD*.mat']));
load(fullfile(matdir,matnm.name))

blockL = size(MUA_stimblocks.stimblock_MUA,1);

basend = size(MUA_stimblocks.params.Bblocks_lims,1);
stimend = basend +  size(MUA_stimblocks.params.blocks_lims,1);
postend = stimend + size(MUA_stimblocks.params.Pblocks_lims,1);


nne = find(~cellfun(@isempty,EventInfo.EventNr));
evifx =nne((find([EventInfo.EventNr{:}]==eventinx)));
ev_ix = [EventInfo.Intrastim_EventInx{evifx} ];

BL = diff((MUA_stimblocks.params.block_ms))*sr/1000;

basend_s = basend*BL;
stimend_s = stimend*BL;
postend_s = postend*BL;