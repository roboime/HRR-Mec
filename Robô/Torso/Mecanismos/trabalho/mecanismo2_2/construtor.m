function [Phiq, nu, Gamma] = construtor (q, qp, t, Phi)

Phiq = jacobian (Phi,q);

nu = -jacobian (Phi,t);

Phiqqp = Phiq*qp;

Phiqqpq = jacobian (Phiqqp,q);

det(Phiq);

Gamma = - Phiqqpq*qp - jacobian(-nu,t);

