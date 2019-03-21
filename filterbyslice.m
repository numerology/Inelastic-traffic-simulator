function [rateOfSlice] = filterbyslice(rates, opBelongs, v)
% filterbyslice function used to filter out the observation belongs to a
% specific slice v.
% rates: 1 x simulationTime cell array
% opBelongs: 1 x simulationTime cell array

simulationTime = length(rates);
rateOfSlice = cell(1, simulationTime);
for t = 1:simulationTime
    rateOfSlice{t} = rates{t}(opBelongs{t} == v);
end

end