function [t,c] = reaction_diffusion_model(d,Dc,Ds,kc,ks,Rc,R,P,c0,Nrc,Nrs,T,Nt,system_type,coating_type)
% Solves the governing reaction-diffusion model

hc = Rc/Nrc;
hs = (R-Rc)/Nrs;
Nr = Nrc+Nrs+1;
dt = T/Nt; % time step size
t = linspace(0,T,Nt+1)'; % discrete times
rc = linspace(0,Rc,Nrc+1);
rs = linspace(Rc,R,Nrs+1);
r = [rc,rs(2:end)]';
rw(1) = r(1); rw(2:Nr) = (r(1:Nr-1)+r(2:Nr))/2; % west boundaries 
re(1:Nr-1) = rw(2:Nr); re(Nr) = r(Nr); % east boundaries

% Finite volume space discretisation
A = zeros(Nr,Nr); b = zeros(Nr,1);

% reflecting inner boundary
A(1,1) = -(Dc*re(1)^(d-1)/hc) - kc*hc/2*r(1)^(d-1);
A(1,2) = Dc*re(1)^(d-1)/hc;
if d == 1
    A(1,1) = A(1,1)/(hc/2*r(1)^(d-1));
    A(1,2) = A(1,2)/(hc/2*r(1)^(d-1));
    b(1) = b(1)/(hc/2*r(1)^(d-1));
end
for i = 2:Nrc % interior nodes
    A(i,i-1) = Dc*rw(i)^(d-1)/hc;
    A(i,i) = -(Dc*re(i)^(d-1) + Dc*rw(i)^(d-1))/hc - kc*hc*r(i)^(d-1);
    A(i,i+1) = Dc*re(i)^(d-1)/hc;
    A(i,i-1) = A(i,i-1)/(hc*r(i)^(d-1));
    A(i,i) = A(i,i)/(hc*r(i)^(d-1));
    A(i,i+1) = A(i,i+1)/(hc*r(i)^(d-1));
end
for i = Nrc+1 % interface node
    A(i,i-1) = Dc*rw(i)^(d-1)/hc;
    A(i,i) = -(Ds*re(i)^(d-1)/hs + Dc*rw(i)^(d-1)/hc) - (kc*hc+ks*hs)/2*r(i)^(d-1);
    A(i,i+1) = Ds*re(i)^(d-1)/hs;
    A(i,i-1) = A(i,i-1)/((hc+hs)/2*r(i)^(d-1));
    A(i,i) = A(i,i)/((hc+hs)/2*r(i)^(d-1));
    A(i,i+1) = A(i,i+1)/((hc+hs)/2*r(i)^(d-1));
end
for i = Nrc+2:Nr-1 % interior nodes
    A(i,i-1) = Ds*rw(i)^(d-1)/hs;
    A(i,i) = -(Ds*re(i)^(d-1) + Ds*rw(i)^(d-1))/hs - ks*hs*r(i)^(d-1);
    A(i,i+1) = Ds*re(i)^(d-1)/hs;
    A(i,i-1) = A(i,i-1)/(hs*r(i)^(d-1));
    A(i,i) = A(i,i)/(hs*r(i)^(d-1));
    A(i,i+1) = A(i,i+1)/(hs*r(i)^(d-1));
end
if isequal(coating_type,'semi-permeable') % reflecting or semi-absorbing boundary
    A(Nr,Nr-1) = Ds*rw(Nr)^(d-1)/hs;
    A(Nr,Nr) = -(Ds*rw(Nr)^(d-1)/hs + re(Nr)^(d-1)*P) - ks*hs/2*r(Nr)^(d-1);
    A(Nr,Nr-1) = A(Nr,Nr-1)/(hs/2*r(Nr)^(d-1));
    A(Nr,Nr) = A(Nr,Nr)/(hs/2*r(Nr)^(d-1));
    b(Nr) = b(Nr)/(hs/2*r(Nr)^(d-1));
end

% Backward Euler time descretisation (Crank Nicolson gave odd behaviour at interface)
theta = 1;
At(1,:) = -A(1,:);
if d == 1
    At(1,:) = [1,zeros(1,Nr-1)]-theta*dt*A(1,:);
end
for i = 2:Nr-1 % interior nodes
    At(i,i-1) = -theta*dt*A(i,i-1);
    At(i,i) = 1-theta*dt*A(i,i);
    At(i,i+1) = -theta*dt*A(i,i+1);
end
if isequal(coating_type,'fully-permeable')
    At(Nr,:) = [zeros(1,Nr-1),1];
else
    At(Nr,:) = [zeros(1,Nr-1),1]-theta*dt*A(Nr,:);
end
At = sparse(At);

c = zeros(Nr,Nt+1);
if isequal(system_type,'monolithic')
    c(:,1) = c0*ones(Nr,1); % initial uniform concentration
elseif isequal(system_type,'core-shell')
    c(:,1) = c0*ones(Nr,1).*(r<Rc) + 0.5*c0*ones(Nr,1).*(r==Rc);
end

% Time stepping
for n = 1:Nt
    bt = c(:,n) + (1-theta)*dt*A*c(:,n) + dt*b;
    if (d ~= 1)
        bt(1) = b(1);
    end
    if isequal(coating_type,'fully-permeable')
        bt(Nr) = 0;
    end
    c(:,n+1) = At\bt;
end
