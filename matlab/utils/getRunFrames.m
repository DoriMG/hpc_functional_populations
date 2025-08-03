function runFrames = getRunFrames(allLoc, track_length, frame_rate)
locDown = allLoc;
locDown = locDown - min(locDown);
locDown = locDown/max(locDown) *99 + 1;

vel = diff(locDown);
vel(vel>5) = 1; % Put this in, because got big upward spikes at transition of env
vel (vel< -0.1*max(vel) )= 0;
vel = vel/max(vel);
vel = [0;vel];
vel = vel*track_length/frame_rate;
runFrames = vel>1;