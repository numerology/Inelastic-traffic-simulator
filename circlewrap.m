function [newTrace] = circlewrap(trace, rad)
% circlewrap: Use wrap-around method to limit arbitrary user trace within a
% circle.
% Params:
%   trace: U x simulationTime x 2, user position as time changes, trace(u, t, 1)
%   is the x-coord of user u at t, trace(u, t, 2) is the y-coord.
%   rad: the radius of the circle.
% Return:
%   newTrace: U x simulationTime x 2, trace after wrapping around.
newTrace = zeros(size(trace));
nUsers = size(trace, 1);
T = size(trace, 2);

for u = 1:nUsers
    for t = 1:T
        trueRad = sqrt(trace(u, t, 1)^2 + trace(u, t, 2)^2);
        if (trueRad > rad) % out-of-bound
            nRad = floor(trueRad / rad);
            leftOver = trueRad - nRad * rad;
            if (mod(nRad, 2) == 1)  % opposite side
                newTrace(u, t, 1) = -(rad - leftOver) * trace(u, t, 1) / trueRad;
                newTrace(u, t, 2) = -(rad - leftOver) * trace(u, t, 2) / trueRad;
            else % same side
                newTrace(u, t, 1) = leftOver * trace(u, t, 1) / trueRad;
                newTrace(u, t, 2) = leftOver * trace(u, t, 2) / trueRad;
            end
        else
            newTrace(u, t, :) = trace(u, t, :);
        end
    end
end

end

