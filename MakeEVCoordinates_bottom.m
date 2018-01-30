function l=MakeEVCoordinates_bottom
    %   Set attributes in pp
    pp = [];
    pp.depth = 100; % the inward depth of the endocytic site.
    pp.dia = 30; % the outer radius.
    pp.thickness = 15; % distance between inner and outer cylinder.
    pp.numMol = 150; % the number of molecular in the region of interest
    pp.view = 'bottom';

    pp.comp = 50;
    
    l = pointsIntheSpace(pp)
end