function Finf = fraction_released_core_shell(d,Ds,ks,R,Rc,P,coating_type)

lambdas = sqrt(ks*R^2/Ds);
Rch = Rc/R;
Rsh = 1 - Rch;

if isequal(coating_type,'fully-permeable')
    if d == 1 % slab
        Finf = 1 / (cosh(lambdas*Rsh));
    elseif d == 2 % cylinder
        Finf = 1 / (Rch*lambdas*(besselk(0,lambdas)*besseli(1,Rch*lambdas) + ...
            besseli(0,lambdas)*besselk(1,Rch*lambdas)));
    elseif d == 3 % sphere
        Finf = lambdas / (sinh(lambdas*Rsh) + lambdas*Rch*cosh(lambdas*Rsh));
    end
elseif isequal(coating_type,'semi-permeable')
    Psh = P*R/Ds;
    if d == 1 % slab
        Finf = Psh / (Psh*cosh(lambdas*Rsh)+lambdas*sinh(lambdas*Rsh));
    elseif d == 2 % cylinder
        Finf = Psh / (Rch*lambdas*((Psh*besselk(0,lambdas)-lambdas*besselk(1,lambdas))*besseli(1,Rch*lambdas) ...
            + (Psh*besseli(0,lambdas) + lambdas*besseli(1,lambdas))*besselk(1,Rch*lambdas)));
    elseif d == 3 % sphere
        Finf = Psh*lambdas / ((Psh-1+Rch*lambdas^2)*sinh(lambdas*Rsh) + lambdas*(Rsh+Psh*Rch)*cosh(lambdas*Rsh));
    end
end
