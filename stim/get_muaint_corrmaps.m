function [maps,stimch_pairs, zone_compartags, ...
    basend, stimend,postend,intdir] = get_muaint_corrmaps(rootdir_ad,patname,ti,evtype,zone_compartags)

patdir = fullfile(rootdir_ad,'Patients',patname);
load(fullfile(patdir,'Events.mat'))


load_muaint_corrmap_pat(rootdir_ad,patname,ti)

uniq_sp = fieldnames(CorrMAPs);
uniq_sp_nr = length(uniq_sp);

[ADzone_labels, ADzone_tags, uniq_stimch, stimch_inx] = ADstate_stimpairs(Events);



% AD zone/ no AD zone tags
if ~isempty(zone_compartags)
    ct_inx  = cellfun(@(x) find(ismember({'AD zone', 'no AD zone'},x )),  zone_compartags);
    zone_compartags = zone_compartags(ct_inx);
    if isempty(ct_inx);
        keyboard;
    end
    zone = []; stimp_nms = cell(uniq_sp_nr,1);
    for m = 1:uniq_sp_nr
        
        stimp_nm = uniq_sp{m};
        snp =stimp_nm; snp(strfind(stimp_nm,'_')) = '-';
        stimp_nms{m} = snp;
        
        shix = find( ismember( upper( uniq_stimch) , upper(snp) ));
        if any( ADzone_tags(  stimch_inx==shix ) ); % AD zone
            zone = [zone 1];
        else % no AD zone
            zone = [zone 2];
        end
    end
end

% Get maps
[maps, stimch_pairs] = deal( cell(1,2) );
[Bends, Sends] = deal(nan(uniq_sp_nr,1));
for m = 1:uniq_sp_nr
    
    sp_nm = uniq_sp{m};
    if length( CorrMAPs.(sp_nm).evs )>3
        
        
        [~, ~, tags] = get_thzone(rootdir_ad,patname,ti,CorrMAPs.(sp_nm).evs(1));
        
        if ismember(upper(tags),upper(evtype)) || isempty(evtype)
            maps{zone(m)} = cat( 3, maps{zone(m)}, CorrMAPs.(sp_nm).coef );
            if isempty(stimch_pairs{zone(m)})
                stimch_pairs{zone(m)}  = stimp_nms(m);
            else
                stimch_pairs{zone(m)}  = [stimch_pairs{zone(m)}, stimp_nms{m}];
            end
            
        end
        
        
    end
    
    
    Bends(m) = CorrMAPs.(uniq_sp{m}).basend;
    Sends(m) = CorrMAPs.(uniq_sp{m}).stimend;
    Pends(m) = size( CorrMAPs.(uniq_sp{m}).coef ,2);
    
end

if length(unique(Bends))==1&& length(unique(Sends))==1 && length(unique(Pends))==1
    basend = Bends(1); stimend = Sends(1); postend = Pends(1);
else
    keyboard;
end

end
