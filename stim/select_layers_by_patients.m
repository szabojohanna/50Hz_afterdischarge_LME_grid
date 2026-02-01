function [z_avgs, laylabs, laynr] = select_layers_by_patients(rootdir_ad,pats,thumbs,z_blM)

laylabs = {'Supragran', 'Infragran','Gran'};
laynr = 3;


load(fullfile(rootdir_ad,'channels_by_layers.mat'));
patnr = length(pats);

for p = 1:patnr
    patname = pats{p};
    ti = thumbs(p);

    islaypat  = ismember(patname,channels_by_layers.Properties.VariableNames);
    if islaypat
        layers{p} = channels_by_layers.(patname)(:,ti);
    end
    
end


z_avgs = cell(length(pats), laynr);
for p = 1:length(pats)
    patname = pats{p};
    if isempty(z_blM{p}); continue; end;
    for k = 1:laynr
        if ~isempty(layers{p})
            switch k;
                case 1
                    chs = [layers{p}{1:3}];
                case 2
                    chs = [layers{p}{5:6}];
                case 3
                    chs = [layers{p}{4}];
            end
        else
            chs = nolaypats_if(patname,laylabs{k});
        end
        
        if ~isempty(chs)
            z_avgs{p,k} =   mean(z_blM{p}(chs,:,:),1); % average across channels within domain
        end
        chs = [];
    end
end


end
