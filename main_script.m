close all
clear all
clc

path_name = cd;
%save_figs = true;
save_figs = false;

%% Plot formatting
font_size = 26;
line_width = 4;
line_width_axis = 1.0;
colors = [0, 47, 108; 200, 16, 46]/255;
font_interpreter = 'LaTeX';

%% Problem setup
% Select test case
Case = 'A'; system_type = 'monolithic'; coating_type = 'fully-permeable';
% Case = 'B'; system_type = 'monolithic'; coating_type = 'semi-permeable';
% Case = 'C'; system_type ='core-shell'; coating_type = 'fully-permeable';
% Case = 'D'; system_type = 'core-shel'; coating_type = 'semi-permeable';

R = 1e-4;
if isequal(system_type,'monolithic')
    D = 1e-12; % diffusivity
    k = 5e-5; % reaction rate
    P = 5e-8; % mass transfer coefficient
    c0 = 0.4; % initial concentration (doesn't affect Finf)
    T = 5e4; % end time
    Nt = 1e4; % number of time steps
elseif isequal(system_type,'core-shell')
    Rc = R/2; % core radius
    Dc = 1e-11; % core diffusivity (doesn't affect Finf)
    Ds = 1e-13; % shell diffusivity
    ks = 5e-5; % shell reaction rate
    P = 5e-8; % mass transfer coefficient
    c0 = 0.4; % initial concentration (doesn't affect Finf)
    T = 3e5; % end time
    Nt = 1e4; % number of time steps
end
Nrc = 1000; % number of nodes in core
Nrs = 1000; % number of nodes in shell
Nr = Nrc+Nrs+1; % total number of nodes

for d = 1:3 % slab, cylinder, sphere
    
    %% Reaction-diffusion model
    r = linspace(0,R,Nr)'; % node locations
    if isequal(system_type,'monolithic')
        Dc = D; Ds = D;
        kc = k; ks = k;
        Rc = R/2;
    elseif isequal(system_type,'core-shell')
        kc = 0;
        D = max(Dc,Ds);
    end
    [Dch,Dsh,kch,ksh,Rh,Th,Rch,Ph,c0h] = deal(Dc/D,Ds/D,R^2*kc/D,R^2*ks/D,R/R,D*T/R^2,Rc/R,R*P/D,c0/c0); % non-dimensionalise
    [th,ch] = reaction_diffusion_model(d,Dch,Dsh,kch,ksh,Rch,Rh,Ph,c0h,Nrc,Nrs,Th,Nt,system_type,coating_type);
    t = R^2*th/D; % dimensional time
    c = c0*ch; % dimensional concentration
    
    %% Released mass over time
    Rt = zeros(length(t),1);
    Rs = R-Rc; % shell radius
    hs = Rs/Nrs; % shell node spacing
    % q = -Ds*(c(Nr,:)-c(Nr-1,:))/h; q = q';
    % q = -Ds*(3/2*c(Nr,:)-2*c(Nr-1,:)+1/2*c(Nr-2,:))/h; q = q';
    % q = -Ds*(11/6*c(Nr,:)-3*c(Nr-1,:)+3/2*c(Nr-2,:)-1/3*c(Nr-3,:))/h; q = q';
    q = -Ds*(25/12*c(Nr,:)-4*c(Nr-1,:)+3*c(Nr-2,:)-4/3*c(Nr-3,:)+1/4*c(Nr-4,:))/hs; q = q';
    dt = T/Nt; % time step size
    t = linspace(0,T,Nt+1)'; % discrete times
    Rt(1) = 0;
    for i = 2:length(t)
        if isequal(system_type,'monolithic')
            Rt(i) = d/(R*c0)*trapz(t(1:i),q(1:i));
        elseif isequal(system_type,'core-shell')
            Rt(i) = d*R^(d-1)/(Rc^d*c0)*trapz(t(1:i),q(1:i));
        end
    end
    
    %% Caculate total fraction of drug released (Finf)
    if isequal(system_type,'monolithic')
        Finf = fraction_released_monolithic(d,D,k,R,P,coating_type);
    elseif isequal(system_type,'core-shell')
        Finf = fraction_released_core_shell(d,Ds,ks,R,Rc,P,coating_type);
    end
    
    %% Plots
    figure;
    set(gcf,'Position',[434   239   560*0.8   520*0.7]);
    set(gcf,'Color','w')
    set(gcf,'Renderer','Painters');
    hold on
    plot(t/T,Rt,'-','LineWidth',line_width,'Color',colors(1,:));
    plot([0,t(end)/T],[Finf,Finf],'--','LineWidth',line_width,'Color',colors(2,:))
    if isequal(font_interpreter,'LaTeX')
        xl = xlabel('Time [$t$]','Interpreter',font_interpreter);
    else
        xl = xlabel('t','Interpreter',font_interpreter);
    end
    set(xl,'Position',[1/2,-0.029])

    if isequal(font_interpreter,'LaTeX')
        yl = ylabel('Fraction released [$F(t)$]','Interpreter',font_interpreter);
    else
        yl = ylabel('F(t)','Interpreter',font_interpreter);
    end
    posyl = get(yl,'Position');
    xlim([0,1])
    ylim([0,1])
    leg.ItemTokenSize = [60,18];
    set(gca,'FontSize',font_size,'LineWidth',line_width_axis,'FontName','Arial',...
        'XTick',[0,1],'XTickLabel',{'$0$','$T$'},...
        'YTick',[0,1],'YTickLabel',{'0','1'},...
        'TickLabelInterpreter',font_interpreter,'Layer','Bottom','Color',1.0*ones(1,3))
    if d == 1
        text(0.4,0.36,['Slab [$d = ',num2str(d,'%i'),'$]'],'Interpreter','LaTeX','FontSize',font_size-2)
    elseif d == 2
        text(0.4,0.36,['Cylinder [$d = ',num2str(d,'%i'),'$]'],'Interpreter','LaTeX','FontSize',font_size-2)
    elseif d == 3
        text(0.4,0.36,['Sphere [$d = ',num2str(d,'%i'),'$]'],'Interpreter','LaTeX','FontSize',font_size-2)
    end
    text(0.4,0.24,['$F_{\infty} = ',num2str(Finf,'%1.4f'),'$'],'Interpreter','LaTeX',...
        'FontSize',font_size,'Color',colors(2,:))
    text(0.4,0.12,['$F(T) = ',num2str(Rt(end),'%1.4f'),'$'],'Interpreter','LaTeX',...
        'FontSize',font_size,'Color',colors(1,:))
    box on
    drawnow
    
    if save_figs
        pause(1);
        print(gcf,[path_name,'Case',Case,num2str(d)],'-depsc2')
    end
    
end
