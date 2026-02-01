function SUAstatSUM
%%
evtype = {'Far'}; %{'Within','Close'}; {'Far'}

SUAresp_laminar(rootdir_ad,patname,evtype)
end

function SUAresp_laminar(rootdir_ad,patname,evtype)

et2= cellfun(@(x) [x '_'],evtype,'UniformOutput',false);
evftit = [et2{:}];

figdir = fullfile(rootdir_ad,'SUA_Analysis','Figures','intrastim_spiking',...
    [patname '_evtypes'],evftit);

load(fullfile(figdir,'Stat.mat'));
fdnms = fieldnames(Stat);

% Stat_c = structfun(@(x) bonf_holm(x.pval, 0.05),Stat,'UniformOutput',0); % correct for multiple comparison
Stat_c = structfun(@(x) x.pval,Stat,'UniformOutput',0);

signU{1} = structfun(@(x) x(1)<=0.05&x(2)>0.05,Stat_c);
signU{2} = structfun(@(x) x(2)<=0.05&x(1)>0.05,Stat_c);
signU{3} = structfun(@(x) x(2)<=0.05&x(1)<=0.05,Stat_c);
signU{4} = signU{1}+signU{2}+signU{3};
U1 = sum(signU{1});
U2 =  sum(signU{2});
U3 =  sum(signU{3});
U4 =  sum(signU{4});
figure;
pie([U1 U2 U3])

us  = find(signU{4}); types = nan(U4,2);
for u = 1:U4
    types(u,:) = Stat.(fdnms{us(u)}).type(1:2);
end
act = find(types(:,1)==1|types(:,2)==1);
inh = find(types(:,1)==-1|types(:,2)==-1);

load(fullfile(rootdir_ad,'channels_by_layers.mat'))


% for k = 4
Us0 = find(signU{4});
for m = 1:2
    switch m
        case 1;    Us = Us0(act); tit = 'act';
        case 2;  Us = Us0(inh);tit = 'inh';
    end
    
    fig = figure;
    col = [.5 .5 .5];
    Unr = length(Us);
    chans = nan(Unr,1);
    for j = 1:Unr
        Uname = fdnms{Us(j)};
        [~, chans(j), ~] = unpack_Uname(Uname);
    end
    
    Ls = nan(Unr,1);
    for c = 1:Unr
        if chans(c)<25
            ti = 1;
        else
            ti = 2; chans(c) = chans(c)-24;
        end
        lays = channels_by_layers.(patname)(:,ti);
        
        Ls(c) = find(cellfun(@(y) any(ismember(chans(c),y)),lays));
        
    end
    edges = .5:6.5;
    hold on;
    h=histogram(Ls,edges,'FaceColor',col);
    view(90,90)
    fnm = fullfile(figdir,['sign_units_hist_' tit]);
    saveas(fig,[fnm '.jpg'])
    saveas(fig,[fnm '.fig'])
    saveas(fig,[fnm '.pdf']);
    close(fig);
end

end