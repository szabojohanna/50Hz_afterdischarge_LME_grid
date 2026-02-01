function [mua_up, mua_down] = get_MUAchange_direction(rootdir_ad, patname, timewin, chans, period)


mua_blockdir = fullfile(rootdir_ad,'Figures','Micro','MUA','MUA_blocks');

patdir = fullfile(rootdir_ad, 'Patients', patname);

load(fullfile(patdir,'Events.mat'));

allevnr = length(Events);

evss = select_evinx(rootdir_ad,patname,Events,false,'all',[], [], [], [], []);

mua_up = []; mua_down = [];
for ei = 1:length(evss)
    evinx = evss(ei);
    
    [mean2bas, ~,~] = calc_muamean(rootdir_ad,patname,evinx,chans,period,timewin);
    
    if round(mean2bas*100)/100>=0.02
        mua_up = [mua_up evinx];
    elseif round(mean2bas*100)/100<=-0.02
        mua_down = [mua_down evinx];
    end
    
end
end