classdef MathParser<interfaces.DialogProcessor
    methods
        function obj=MathParser(varargin)      
            obj@interfaces.DialogProcessor(varargin{:}) ;  
        end
        
        function out=run(obj,p)
            out=[];
            obj.setPar('undoModule','Math Parser');
            notify(obj.P,'backup4undo');
            locs = obj.locData.loc;
            evalstr=p.equation;
            fn=fieldnames(locs);
            for k = 1:length (fn)
                ind=strfind(evalstr,fn{k});
                for m=1:length(ind)
                    if ind(m)>1
                        lbefore=evalstr(ind(m)-1);
                    else
                        lbefore='';
                    end
                    if lbefore~='.' && ~isstrprop(lbefore,'alpha')
                    evalstr=[evalstr(1:ind(m)-1) 'locs.' evalstr(ind(m):end)];
                    end
%                     evalstr = strrep(evalstr,[lbefore fn{k}],[lbefore 'locs.'  fn{k}]);
                end
            end
             try
                newval=eval(evalstr);
                 obj.locData.setloc(p.resultfield,newval)
                 obj.locData.filter(p.resultfield)
                 obj.locData.regroup;
                 obj.setPar('locFields',fieldnames(obj.locData.loc))
             catch
                 disp('could not evaluate equation')
             end
             
        end
        function pard=guidef(obj)
            pard=guidef(obj);
        end
        function dobutton(obj,callobj,event,inp1)
            disp('dobutton')
           x= obj.locData.loc.xnm;
           mean(x)
        end
    end
end

function pard=guidef(obj)
pard.t1.object=struct('String','result field','Style','text');
pard.t1.position=[2,1];

pard.resultfield.object=struct('String','within','Style','edit');
pard.resultfield.position=[3,1];
pard.resultfield.Width=0.8;


pard.t2.object = struct('String','=','Style','text');
pard.t2.position=[3,1.8];
pard.t2.Width=0.2;

pard.t3.object=struct('String','Equation. Use fieldnames (e.g. xnm, phot) as variables','Style','text');
pard.t3.position=[2,2];
pard.t3.Width=3;

pard.equation.object = struct('String','(locprecnm<25 & PSFxnm>100) | numberInGroup>1','Style','edit');
pard.equation.position=[3,2];
pard.equation.Width=3;
pard.plugininfo.type='ProcessorPlugin';
end