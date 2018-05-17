function [C] = eci2lvlh(r,v)
    % returns the transformation matrix from lvlh to eci
    
    % defining lvlh basis
    z1 = -r/norm(r);
    y1 = -(cross(r,v))/norm((cross(r,v)));
    x1 = cross(y1,z1);
    
    % defining eci basis
    x2 = [1;0;0];
    y2 = [0;1;0];
    z2 = [0;0;1];
    
    C = [dot(x2,x1) dot(x2,y1) dot(x2,z1);...
        dot(y2,x1) dot(y2,y1) dot(y2,z1);...
        dot(z2,x1) dot(z2,y1) dot(z2,z1)];
    
    
end

