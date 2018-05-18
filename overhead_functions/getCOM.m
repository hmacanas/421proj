% Code bdim(2) Michael Johnston, Madim(2) 18, 2018
% Spacecraft center of mass and inertia matridim(1)calculation
% Inputs:	1dim(1)n mass matridim(1)
%			3dim(1)n distance from origin matridim(1)
%			3dim(1)n component dimension matridim(1)
% Outputs:	Vector from origin to COM
%			Inertia matridim(1)

function [COM,I] = getCOM(masses, Rs, dims)
	len = length(Rs);
	rTimesM = eye(3,len);
	for i = 1:len
		rTimesM(:,i) = Rs(:,i).*masses(i);
	end
	COM = sum(rTimesM,2)/sum(masses);
	
	% -- Moment of inertia calculation
	Icuboid = @(m,dim) m/12 * [dim(3)^2 + dim(2)^2, 0, 0;...
								0, dim(1)^2 + dim(3)^2, 0;...
								0, 0, dim(1)^2 + dim(2)^2]; % Inertia matridim(1) for cuboid
	Jpar = @(I,m,R) I + m*(dot(R,R)*eye(3,3) - R*R'); % Parallel adim(1)is theorem
	
	% Calculate inertia matridim(1)
	I = eye(3);
	for i = 1:len
		I = I + Jpar(Icuboid(masses(i),dims(:,i)), masses(i), Rs(:,i)-COM);
	end
end