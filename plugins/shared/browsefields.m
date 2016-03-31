function str=browsefields(prop,field,recursions,recursive)
if nargin<4
    recursive=0;
end
if nargin<3
    recursions=inf;
end
if nargin<2||isempty(field)
    fn=fieldnames(prop);
    field='';
    psub=prop;
elseif strcmp(field,'..')
    fn=fieldnames(prop);
    psub=prop;
elseif isstruct(prop.(field))
    fn=fieldnames(prop.(field));
    fn={'..',fn{:}};
    psub=prop.(field);
else
    str=field;
    return;
end
pos=get(0,'PointerLocation');
pos(1)=max(1,pos(1)-50);
pos(2)=max(1,pos(2)-50);
pos(3:4)=[100 100];
    answ=mylistdlg('ListString',fn,'SelectionMode','single','FontSize',14,'InitialValue',min(2,length(fn)),'Position',pos);
    if isempty(answ)
        str='abortbutton';
    elseif recursions<=0
        str=[field '.' fn{answ}];
    else
%         field
        field2=fn{answ};
        if strcmp(field2,'..')
            str=browsefields(prop,field2);
        elseif strcmp(field,'..')
            str=browsefields(prop,field2);
        else
            str=[field '.' browsefields(psub,field2,recursions-1,1)];
        end
    end
    if  ~recursive&&~isempty(strfind(str,'abortbutton'))
        str=[];
    elseif str(1)=='.'
        str=str(2:end);
    end
end