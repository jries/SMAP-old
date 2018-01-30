function l=MakeEVCoordinates_bottom2
    %   Set attributes in pp
    pp = [];
    pp.depth = 120; % the inward depth of the endocytic site.
    pp.dia = 50; % the outer radius.
    pp.thickness = 30; % distance between inner and outer cylinder.
    pp.numMol = 300; % the number of molecular in the region of interest
    pp.view = 'bottom';
    
    pp.comp = 50;

    l = pointsIntheSpace(pp)
end