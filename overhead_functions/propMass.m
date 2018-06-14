function [m] = propMass(t,T)

Isp = 60; % [s]
g0 = 9.81; % [m/s2]

m = zeros(1,length(t));
for ii = 1:length(t)-1
    mfr = T(ii+1)/(Isp*g0);
    m(ii) = mfr*(t(ii+1) - t(ii));
end
m = sum(m);

end