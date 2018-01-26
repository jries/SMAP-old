function l=MakeEVcoordinates(view, pp)
    %   Defining the ring region
    factor = (1/2*pi*(pp.outRadius)^2-1/2*pi*(pp.inRadius)^2)/((pp.size/2)^2); % area(ring)/area(coodinates)
    centerShift = (pp.size/2)-pp.outRadius
    nInputDot = round(pp.nFinalDot/factor)*2; % make sure there will be enough dots line within the ring
    rx = pp.size.*rand(nInputDot,1); % random sampled nInputDot float numbers from 0:1 as dot.x
    ry = pp.size.*rand(nInputDot,1); % random sampled nInputDot float numbers from 0:1 as dot.y
    rz = pp.depth.*rand(nInputDot,1)
    %	Filtering by the ring region
    idx = (rx-pp.outRadius).^2+(ry-pp.outRadius).^2<pp.outRadius^2 & (rx-pp.outRadius).^2+(ry-pp.outRadius).^2>pp.inRadius^2;
    rx = rx(idx);
    ry = ry(idx);
    rz = rz(idx);
    rx = rx(1:pp.nFinalDot,:);
    ry = ry(1:pp.nFinalDot,:);
    rz = rz(1:pp.nFinalDot,:);
    % figure(20)
    % scatter(rx,ry)
    rx = rx+centerShift;
    ry = ry+centerShift;
switch view    
    case 'bottom'
    l.x = rx';
    l.y = ry';
    
    case 'side'
    l.x = rx';
    l.y = rz';
    
    otherwise
end