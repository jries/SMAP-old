function l=MakeEVCoordinates_side
    %   Set attributes in pp
    pp1 = [];
    pp1.depth = 40; % the inward depth of the endocytic site.
    pp1.dia = 30; % the outer radius.
    pp1.thickness = 10; % distance between inner and outer cylinder.
    pp1.numMol = 150; % the number of molecular in the region of interest
    pp1.view = 'side';
    
    pp1.shapeIn = 'columnNT';
    pp1.shapeOut = 'columnNT';

    pp1.comp = 50;
    
    l1 = pointsIntheSpace(pp1);
    l1.channel = repelem(1, length(l1.x))';
    l1 = struct2table(l1);
    %   Set attributes in pp
    pp2 = [];
    pp2.depth = 120; % the inward depth of the endocytic site.
    pp2.dia = 50; % the outer radius.
    pp2.thickness = 30; % distance between inner and outer cylinder.
    pp2.numMol = 300; % the number of molecular in the region of interest
    pp2.view = 'side';
    
    pp2.shapeIn = 'columnRT';
    pp2.shapeOut = 'parabola';
        
    pp2.comp = 50;

    l2 = pointsIntheSpace(pp2);
    l2.channel = repelem(2, length(l2.x))';
    l2 = struct2table(l2);
    l = vertcat(l1, l2);
	l=table2struct(l, 'ToScalar',true);
end