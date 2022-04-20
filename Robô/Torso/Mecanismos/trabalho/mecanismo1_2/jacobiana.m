function J=jacobiana(x,t)
l1 = 5;
l2 = 5;
J=[1 0 0 0 0 0;
    0 1 0 0 0 0;
    1 0 -l1*sin(x(3)) -1 0 -l2*sin(x(6));
    0 0 l1*cos(x(3)) 0 -1 l2*cos(x(6)); 
    0 0 0 1 0 0;
    0 0 0 0 1 0];