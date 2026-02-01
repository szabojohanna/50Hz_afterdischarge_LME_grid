function evinx = select_evinx(rootdir_ad,patname,Events,maxintens,evtype,ti, minL, timewin, upchans, period)


isgr = find( cellfun(@(x) contains(lower(x),'gr')|contains(lower(x),'cov'),{Events.Label}));
is_stimlim = find( ~cellfun(@isempty,{Events.Stim_start}) );

xx = intersect(isgr,is_stimlim);

is_thumbrec = find(cellfun(@(x) ~isempty(x),{Events.Thumbfile}));
goodx = intersect(xx,is_thumbrec);

if strcmp(patname,'Pt20')
    goodx(goodx==1) = []; % noisy event
end


if isempty(upchans)&& contains(evtype,'mua');
    load(fullfile(rootdir_ad, 'good_channels.mat'));
    upchans = good_channels.(patname){ti};
end;

if strcmp(evtype,'close')
    clos_events = [];clos_dists = {};
    for eventinx = 1:size(Events,2)
        
        if isfield(Events,['Th' num2str(ti) '_STIMdist_cm'])
            
            
            if ~isempty(Events(eventinx).(['Th' num2str(ti) '_STIMdist_cm']))
                clos_events = [clos_events eventinx];
                clos_dists{eventinx} = Events(eventinx).(['Th' num2str(ti) '_STIMdist_cm']);
                
            end;
        else
            continue;
        end
    end
    
elseif contains(evtype, 'mua_up')
    [clos_events, ~] = get_MUAchange_direction(rootdir_ad, patname, timewin, upchans, period);
    
elseif contains(evtype, 'mua_down')
    
    [~, clos_events] = get_MUAchange_direction(rootdir_ad, patname, timewin, upchans, period);
    
else strcmp(evtype, 'all')
    clos_events = 1:size(Events,2);
end




if maxintens
    %     [evinx, ~] = find_maxintens_evs(Events, evinx);
    [maxx, ~] = find_maxintens_evs(Events, goodx);
    
    evinx = intersect(maxx,clos_events);
else
    evinx = intersect(goodx,clos_events);
end


if ~isempty(minL)
    evinx0 = evinx;
    stimLs = [Events(evinx).Stim_end] -[Events(evinx).Stim_start];
    goodinx = stimLs>=minL;
    evinx = evinx0(goodinx);
    fprintf('excluded: Ev%d \n', evinx0(~goodinx)')
end
end