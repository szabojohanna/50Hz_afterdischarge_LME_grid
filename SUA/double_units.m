function double_units(rootdir_ad)
sr = 20000;
load(fullfile(rootdir_ad,'SUA_Analysis','SUA_data'));
Unames = unit_table.Unit_name;
unit_table0 = unit_table;
Unr = length(Unames);
clus_inxs = cell(Unr,1);
for u = 1:Unr
    Uname = Unames{u};
    [patname, chan, clnr] = unpack_Uname(Uname);

    savedir = fullfile(rootdir_ad,'SUA_Analysis',patname,['Chan' num2str(chan)]);
    fnm = ['Chan' num2str(chan) '_AllEvents_Raw'];
    % load(fullfile(savedir,[fnm '.mat']));
    load(fullfile(savedir,[fnm '_spikes.mat']))
    load(fullfile(savedir,['times_' fnm '.mat']))
    
    sp_ts_s = index/1000; % timestamp of spikes in sec
    spk_inx = round(sp_ts_s*sr); % indices of spikes
    
    clus_inxs{u} = spk_inx(cluster_class(:,1)==clnr);
end

uL = nan(Unr,Unr);
uP = nan(Unr,Unr);
uP2 = nan(Unr,Unr);
for u1 = 1:Unr
    for u2 = u1+1:Unr
%         ui1 = find(ismember(Unames,Uname1));
%         ui2 = find(ismember(Unames,Uname2));
        
        [udub]=intersect(clus_inxs{u1},clus_inxs{u2});
        uL(u1,u2) =length(udub);
        uP(u1,u2) =(length(udub)*100)/(length([clus_inxs{u1} ]));
        uP2(u1,u2) =(length(udub)*100)/(length([clus_inxs{u2}]));
    end
end
dubnr = sum(uP(:)>50|uP2(:)>50);
[ia,ib]=find(uP>50|uP2>50);

Unit_name = Unames;
clus_inxs2 = clus_inxs;
dnr = 0;
for h = 1:dubnr
    uu(1) = ia(h); uu(2) = ib(h);
    [~, chan1, ~] = unpack_Uname(Unames{uu(1)});
    [~, chan2, ~] = unpack_Uname(Unames{uu(2)});
    if abs(diff([chan1 chan2]))==1
        
        new_clusinx = unique(([clus_inxs{uu(1)} clus_inxs{uu(2)}]));
        [~,mi] = min([length(clus_inxs{uu(1)}) length(clus_inxs{uu(2)})]);
        [~,ma] = max([length(clus_inxs{uu(1)}) length(clus_inxs{uu(2)})]);
        Unit_name{uu(mi)} = [];
        clus_inxs2{uu(mi)} =  [];
        clus_inxs2{uu(ma)} = new_clusinx;
        dnr = dnr+1;
    end
end

emp = cellfun(@isempty,Unit_name);
Unit_name(emp)= [];
clus_inxs2(cellfun(@isempty,clus_inxs2))= [];

unit_table = table(Unit_name);
unit_table.Contamination_rate = cell(height(unit_table),1);
unit_table.Unit_quality = cell(height(unit_table),1);
unit_table.Cluster_inx = cell(height(unit_table),1);

unit_table.Unit_quality = unit_table0.Unit_quality(~emp);
unit_table.Contamination_rate = unit_table0.Contamination_rate(~emp);
unit_table.Cluster_inx = clus_inxs2;

save(fullfile(rootdir_ad,'SUA_Analysis','SUA_data_corrected'),'unit_table');


end
