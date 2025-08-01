function velRes = getVelocity(loc, goal)

locsz = max(size(loc));
runidx = 1:locsz;                                 % Index
runidxq = linspace(min(runidx), max(runidx), goal);    % Interpolation Vector
locDown = interp1(runidx, loc, runidxq, 'linear');       % Downsampled Vector

% Normalize loc
locDown = locDown - min(locDown);
locDown = locDown/max(locDown) *99 + 1;

vel = diff(locDown);
vel(vel>5) = 1; % Put this in, because got big upward spikes at transition of env
vel (vel< -0.1*max(vel) )= 0;
vel = vel/max(vel);

pupilsz = size(vel, 2);
runidx = 1:pupilsz;                                 % Index
runidxq = linspace(min(runidx), max(runidx), goal);    % Interpolation Vector
vel = interp1(runidx, vel, runidxq, 'linear');       % Downsampled Vector
velRes = smooth(vel,10)';
