function str=browsefields(prop,field,recursions,recursive,parsefield)
if nargin<5
    parsefield=false;
end
if nargin<4||isempty(recursive)
    recursive=0;
end
if nargin<3||isempty(recursions)
    recursions=inf;
end
if nargin<2||isempty(field)
    fn=fieldnames(prop);
    fnt=fn;
    field='';
    psub=prop;
elseif strcmp(field,'..')
    fn=fieldnames(prop);
    fnt=fn;
     recursions=recursions+1;
    psub=prop;
    field='';
elseif isstruct(prop.(field))
    fn=fieldnames(prop.(field));
    fnt=fn;
    if parsefield %&&recursions==0
        for k=1:length(fn)
            ph=prop.(field).(fn{k});
            if iscell(ph)
            fnt{k}=ph{end};
            
            elseif isfield(ph,'module')
                fnt{k}=ph.module{end};
            end
        end
    end
    
    fn={'..',fn{:}};
    fnt={'..',fnt{:}};
    psub=prop.(field);
else
    str=field;
    return;
end
pos=get(0,'PointerLocation');
pos(1)=max(1,pos(1)-50);
pos(2)=max(1,pos(2)-50);
pos(3:4)=[100 100];

answ=mylistdlg('ListString',fnt,'SelectionMode','single','FontSize',14,'InitialValue',min(2,length(fn)),'Position',pos,'ListSize',[240,300]);
if isempty(answ)
    str='abortbutton';
elseif recursions<=0 && ~strcmp(fn{answ},'..')
    %return final field
    str=[field '.' fn{answ}];    
else
%         field
    field2=fn{answ};
    if strcmp(field2,'..')
        str=browsefields(prop,field2,recursions,true,parsefield);
    else
        str=[field '.' browsefields(psub,field2,recursions-1,1,parsefield)];
    end
end
if  ~recursive&&~isempty(strfind(str,'abortbutton'))
    str=[];
elseif str(1)=='.'
    str=str(2:end);
end
end