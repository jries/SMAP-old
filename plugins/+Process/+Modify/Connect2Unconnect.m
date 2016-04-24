classdef Connect2Unconnect<interfaces.DialogProcessor
    methods
        function obj=Connect2Unconnect(varargin)     
            obj@interfaces.DialogProcessor(varargin{:}) ;  
        end
        
        function out=run(obj,p) 
            out=[];
           switch p.connectmode.selection
               case 'connect->unconnect'
                   disp('if problems, tell Jonas')
                   indunc=obj.locData.loc.channel==p.channel;
                   indcon=obj.locData.grouploc.channel==p.channel;
                   obj.locData.removeind(indunc);
                   loccopy=obj.locData.copy;
                   loccopy.removeind(~indcon,'grouploc');
                   loccopy.loc=loccopy.grouploc;
                   obj.locData.addLocData(loccopy);
                   obj.locData.regroup;

               otherwise 
                   disp('not implemented')
           end  
        end
        function pard=guidef(obj)
            pard=guidef;
        end
    end
end




function pard=guidef



pard.textb.object=struct('String','Channel','Style','text');
pard.textb.position=[1,1];
pard.channel.object=struct('String','1','Style','edit');
pard.channel.position=[1,2];



pard.connectmode.object=struct('String','connect->unconnect|unconnect->connect','Style','popupmenu','Value',1);
pard.connectmode.position=[2,1];
end