function struct2uitable(htable, s,flip)
if nargin<3
    fliph=false;
else
    if strcmp(flip,'flip')
        fliph=true;
    else 
        fliph=false;
    end
end
tab=struct2table(s);


if flip
    htable.Data=table2cell(tab)';
    htable.RowName=tab.Properties.VariableNames;
else
    htable.Data=table2cell(tab);
    htable.ColumnName=tab.Properties.VariableNames;
end


