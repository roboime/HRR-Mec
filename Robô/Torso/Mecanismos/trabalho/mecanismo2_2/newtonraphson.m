function x=newtonraphson(x0,t)
    %definindo indicador e par�metro de converg�ncia
    tol=0.0001;
    x=x0;
    z=10000*x0;
    while abs(norm(z))>tol
        %constru��o da Jacobiana (A) e do vetor -F
        A=jacobiana(x,t);
        b=-Phi(x,t);
        z=eliminacaogauss(A,b);
        x=x+z;
    end