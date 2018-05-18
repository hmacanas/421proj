% Code by Michael Johnston, May 18, 2018
% Spacecraft center of mass calculation
% 

function COM = getCOM(masses, Rs)
	len = length(Rs);
	rTimesM = zeros(3,len);
	for i = 1:len
		rTimesM(:,i) = Rs(:,i).*masses(i);
	end
	COM = sum(rTimesM,2)/sum(masses);
end