function Finf = fraction_released_monolithic(d,D,k,R,P,coating_type)

lambda = sqrt(k*R^2/D);

if isequal(coating_type,'fully-permeable')
    if d == 1 % slab
        Finf = sinh(lambda)/(lambda*cosh(lambda));
    elseif d == 2 % cylinder
        Finf = 2*besseli(1,lambda)/(lambda*besseli(0,lambda));
    elseif d == 3 % sphere
        Finf = 3*(lambda*cosh(lambda)-sinh(lambda))/(lambda^2*sinh(lambda));
    end
elseif isequal(coating_type,'semi-permeable')
    Ph = P*R/D;
    if d == 1 % slab
        Finf = Ph*sinh(lambda)/(lambda*(lambda*sinh(lambda)+Ph*cosh(lambda)));
    elseif d == 2 % cylinder
        Finf = 2*Ph*besseli(1,lambda)/(lambda*(lambda*besseli(1,lambda)+Ph*besseli(0,lambda)));
    elseif d == 3 % sphere
        Finf = 3*Ph*(lambda*cosh(lambda)-sinh(lambda))/(lambda^2*(lambda*cosh(lambda)+(Ph-1)*sinh(lambda)));
    end
end