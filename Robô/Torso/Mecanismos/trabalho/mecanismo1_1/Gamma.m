function b=Gamma(x,xp,t)
l1 = 5;
l2 = 5;
b=[0; 
    0; 
    l1*cos(x(3))*xp(6)^2 + l2*cos(x(6))*xp(6)^2; 
    l1*sin(x(3))*xp(6)^2 + l2*sin(x(6))*xp(6)^2; 
    -cos(t)*(l1 + l2 - sqrt(l1^2 + l2^2))/2; 
    -sin(t)*(l1 + l2 - sqrt(l1^2 + l2^2))/2];