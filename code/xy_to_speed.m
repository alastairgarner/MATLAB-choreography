%% generate_plots.m

% Other m-files required: 
% Subfunctions: 
% MAT-files required: 
% See also:

% Author: Alastair Garner
% email: alastairgarner@outlook.com
% Website: https://github.com/alastairgarner/
% August 2019; Last revision: 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function distance_per_frame = xy_to_speed(x,y,et,frame)

    xDiff = [0 diff(x)];
    yDiff = [0 diff(y)];
    xyDistRaw = sqrt([xDiff.^2]+[yDiff.^2]);
    xSqDist = movmean(xDiff,frame).^2;
    ySqDist = movmean(yDiff,frame).^2;
    xyDist = sqrt(xSqDist + ySqDist);
    distance_per_frame = movmean(xyDist,21) ./ mean(diff(et));

end