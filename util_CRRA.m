function [util] = util_CRRA(consumption, gamma)
util = consumption^(1-gamma)/(1-gamma);
end

