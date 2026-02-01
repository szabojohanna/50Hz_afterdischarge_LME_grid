function label_layers(layers)

xL = xlim; 
Laynr = length(layers);
laylim = [];
for g = 1:Laynr
    if ~isempty(layers{g})
        laylim = [laylim (layers{g}(1))];
        Lg = g;
    end
end
laylim = [laylim (layers{Lg}(end))];

for a = 1:Lg+1
line([xL(1) xL(1)+diff(xL)*.15], [laylim(a); laylim(a)], 'Color','k','LineStyle','-')    
end

laylab = diff(laylim)/2 + laylim(1:end-1);

if Lg==7
    laylabtxt = [arrayfun(@num2str, 1:6,'UniformOutput',0), {'WM'}];
else
    laylabtxt = arrayfun(@num2str, 1:Lg,'UniformOutput',0);
end
text(repmat(xL(1)+diff(xL)*.07,[1 length(laylab)]),laylab,laylabtxt)