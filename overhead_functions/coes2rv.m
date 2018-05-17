function [ r2,v2 ] = coes2rv(ecc,inc,RAAN,W,h,nu )
    mu = 398600;
    r=(h^2/mu)*(1/(1+ecc*cos(nu)));
    rx=r*cos(nu);
    ry=r*sin(nu);
    rperifocal=[rx;ry;0];

    vx=mu/h*(-sin(nu));
    vy=mu/h*(ecc+cos(nu));
    vperifocal=[vx;vy;0];
    
    r2=(cz(W)*cx(inc)*cz(RAAN))'*rperifocal;%km
    v2=(cz(W)*cx(inc)*cz(RAAN))'*vperifocal;%km/s


end

