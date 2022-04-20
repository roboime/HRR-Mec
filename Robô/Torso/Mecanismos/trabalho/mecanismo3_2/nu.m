function b=nu(x,t)
l1 = 5;
l2 =3;
b = [zeros(13,1);
     -sin(t)*(l1/4 + l2/4);
     cos(t)*(l1/4 + l2/4);
     1];