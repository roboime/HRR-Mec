function b=nu(x,t)
l1 = 5;
l2 =5;
b = [zeros(18,1);
     -1;
     1; 
     -l1*sin(t)/10;
     l1*cos(t)/10];