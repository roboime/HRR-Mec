function b=nu(x,t)
l1 = 5;
l2 =5;
b = [zeros(4,1); 
    -sin(t)*(l1 + l2 - sqrt(l1^2 + l2^2))/2; 
    cos(t)*(l1 + l2 - sqrt(l1^2 + l2^2))/2];