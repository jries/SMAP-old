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
% tab.names=tab.Properties.VariableNames;

if fliph
    data=table2cell(tab)';
    data=[tab.Properties.VariableNames' data];
%     horzcat(tab.Properties.VariableNames, data)
    htable.Data=data;
    
    htable.RowName={};
else
    htable.Data=table2cell(tab);
    htable.ColumnName=tab.Properties.VariableNames;
end


