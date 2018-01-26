function l=MakeEVSideCoordinates
    %   Set attributes in pp
    pp = [];
    pp.size = 100; % size of the cooridnates, e.g., 60 means a 60*60 coordinates
    pp.nFinalDot = 100; % the final number of molecular of interest in the ring region.
    pp.inRadius = 10; % the iner radius or the ring.
    pp.outRadius = 20; % the outer radius or the ring.
    pp.depth = 100; % the inward depth of the endocytic site.

    l = MakeEVcoordinates('side', pp)
end