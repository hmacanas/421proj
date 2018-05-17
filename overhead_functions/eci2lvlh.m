function [C] = eci2lvlh(r,v)
    % returns the transformation matrix from lvlh to eci
    
    % defining lvlh basis
    z2 = -r/norm(r);
    y2 = -(cross(r,v))/norm((cross(r,v)));
    x2 = cross(y2,z2);
    
    % defining eci basis
    x1 = [1;0;0];
    y1 = [0;1;0];
    z1 = [0;0;1];
    
    C = [dot(x2,x1) dot(x2,y1) dot(x2,z1);...
        dot(y2,x1) dot(y2,y1) dot(y2,z1);...
        dot(z2,x1) dot(z2,y1) dot(z2,z1)];
    
    
end

