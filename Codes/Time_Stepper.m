function [runtime,w_old] = Time_Stepper(pars,Sol,hands,folder)
q=pars.q; x=Sol.x; y=Sol.y; usol=Sol.usol; vsol=Sol.vsol; musol=Sol.musol;
fx=hands.f0; gx=hands.g0; k=pars.k;


dt=1/80;
% dt: time step
% steps: number of spatial points in each coordinate direction

%% Model Paramters and initial conditions
Dv=q(1);beta=q(2);
epsln1 = 1; epsln2 = Dv;

musol = real(musol);
% create nodes
steps=length(x);
h = abs(x(1,1)-x(1,2));
fprintf('dt=%f h=%f\n',dt,h)
nnodes = steps^2;
nodes = zeros(nnodes,2);
j = 1;
for p = 1 : steps
        for i = 1:steps
               nodes(j,:) = [x(i) y(p)];
            j = j+1;
        end
end
nb = 2*nnodes; % becuase we are solving a syste of 2 RDE

% discretize time interval
tmax=600;
t = 0:dt:tmax; tlen = length(t);

% initial condition for u
u_old = usol; 


% initial condition for v
v_old = vsol;

% Stacking nodes for evolution
w_old = zeros(nb,1);
w_old(1:2:nb-1) = u_old; w_old(2:2:nb) = v_old; 

%% Block matrix Assembly
C = zeros(2);
C(1,1) = (epsln1*dt)/h^2;
C(2,2) = (epsln2*dt)/h^2;
Q = blktridiag(2,-1,-1,steps);
Q(1,2) = -2; Q(steps,steps-1) = -2; I = speye(steps);
A1 = kron(kron(I,Q),C);A2 = kron(kron(Q,I),C); 
I = speye(nb); 
M1 = sparse(I+A1); M2 = sparse(I+(1/3)*A1); M3 = sparse(I+(1/4)*A1);
M11 = sparse(I+A2); M22 = sparse(I+(1/3)*A2); M33 = sparse(I+(1/4)*A2);

%% Time Evolution 
[L1,U1] =lu(M1);
[L2,U2] =lu(M2);
[L3,U3] =lu(M3);
 [L11,U11] =lu(M11); 
 [L22,U22] =lu(M22);
 [L33,U33] =lu(M33);

 Usoln = reshape(u_old,steps,steps); 
  Vsoln = reshape(v_old,steps,steps);
  Uplot =real(Usoln(:,:));
  Vplot = real(Vsoln(:,:));

if isempty("Veg.mat")==1
Veg=flipud(summer);
else
load("Veg.mat","Veg")
end
if isempty("Water.mat")==1
Water=flipud(winter);
else
load("Water.mat","Water")
end


scrsz = get(0,'ScreenSize');  
figure(1);close gcf;figure(2);close gcf; 

        
    figure('Position',[0*scrsz(3)/4 scrsz(4)/4 scrsz(3)/2 scrsz(3)/3]);
    surf(x,y,Uplot')
    title(['u at t=' num2str(0)]);
    xlabel('x')
    ylabel('y')
    colormap(Veg);
    colorbar;
    caxis([0 max(max(Uplot))]);
    set(gca,'LineWidth', 1);
    set(gca,'FontSize',10);
    set(gca,'FontWeight','bold');
    pbaspect(gca,[1 1 1]);view(0,90);shading interp;
    xlim([-x(end) x(end)]);ylim([-x(end) x(end)]);
    set(gca,'XTick',[-1 -0.8 -0.6 -0.4 -0.2 0 0.2 0.4 0.6 0.8 1]*x(end));
    set(gca,'YTick',[-1 -0.8 -0.6 -0.4 -0.2 0 0.2 0.4 0.6 0.8 1]*x(end));
    FigNameu = [folder,'\sliceu',num2str(0)];
    SolnNameu = [FigNameu,'.mat'];
    U=Uplot';
    print(FigNameu,'-depsc2');
    save(SolnNameu,'x','y','U');
    
    figure('Position',[2*scrsz(3)/4 scrsz(4)/4 scrsz(3)/2 scrsz(3)/3]);
    surf(x,y,Vplot' + (beta*Dv/(Dv-1)).*Uplot')
    title(['v at t=' num2str(0)]);
    xlabel('x')
    ylabel('y')
    colormap(Water);
    colorbar;
    caxis([0 max(max(Vplot + (beta*Dv/(Dv-1)).*Uplot))]);
    set(gca,'LineWidth', 1);
    set(gca,'FontSize',10);
    set(gca,'FontWeight','bold');
    pbaspect(gca,[1 1 1]);view(0,90);shading interp;
    xlim([-x(end) x(end)]);ylim([-x(end) x(end)]);
    set(gca,'XTick',[-1 -0.8 -0.6 -0.4 -0.2 0 0.2 0.4 0.6 0.8 1]*x(end));
    set(gca,'YTick',[-1 -0.8 -0.6 -0.4 -0.2 0 0.2 0.4 0.6 0.8 1]*x(end));
    V=Vplot' + (beta*Dv/(Dv-1)).*Uplot';
    FigNamev = [folder,'\slicev',num2str(0)];
    SolnNamev = [FigNamev,'.mat'];
    print(FigNamev,'-depsc2');
    save(SolnNamev,'x','y','V');
 pause;
tic
for i = 2:tlen
     F_old = F(w_old,musol,fx,gx);
     w_star = U11\(L11\(w_old + dt*F_old));
     w_star = U1\(L1\w_star);
     F_star = F(w_star,musol,fx,gx);
     a_1 = U2\(L2\w_old);
     b_1 = U3\(L3\w_old);
     c_1 = 9*a_1 - 8*b_1;
     a_2 = U2\(L2\F_old);
     b_2 = U3\(L3\F_old);
     c_2 = 9*a_2-8*b_2;
     d_1 = U22\(L22\(9*c_1+2*dt*c_2 + dt*F_star));
     d_2 = U33\(L33\(-8*c_1-(3/2)*dt*c_2-(dt/2)*F_star));
     w_old = d_1+d_2;              
     
    if mod(i*dt,1)==0 
        u_soln = w_old(1:2:nb-1); 
 v_soln = w_old(2:2:nb);
 Usoln = reshape(u_soln,steps,steps); 
  Vsoln = reshape(v_soln,steps,steps);
  Uplot = real(Usoln(:,:));
  Vplot = real(Vsoln(:,:));
        
     figure(1)
    surf(x,y,Uplot')
    title(['u at t=' num2str(i*dt)]);
    xlabel('x')
    ylabel('y')
    colormap(Veg);
    colorbar;
    caxis([0 max(max(Uplot))]);
    set(gca,'LineWidth', 1);
    set(gca,'FontSize',10);
    set(gca,'FontWeight','bold');
    pbaspect(gca,[1 1 1]);view(0,90);shading interp;
    xlim([-x(end) x(end)]);ylim([-x(end) x(end)]);
    set(gca,'XTick',[-1 -0.8 -0.6 -0.4 -0.2 0 0.2 0.4 0.6 0.8 1]*x(end));
    set(gca,'YTick',[-1 -0.8 -0.6 -0.4 -0.2 0 0.2 0.4 0.6 0.8 1]*x(end));
    G1(i*dt) = getframe(figure(1)) ;drawnow;
    if mod(i*dt,100)==0 
    FigNameu = [folder,'\sliceu',num2str(i*dt)];
    SolnNameu = [FigNameu,'.mat'];
    U=Uplot';
    print(FigNameu,'-depsc2');
    save(SolnNameu,'x','y','U');
    end

    figure(2)
    surf(x,y,Vplot' + (beta*Dv/(Dv-1)).*Uplot')
    title(['v at t=' num2str(i*dt)]);
    xlabel('x')
    ylabel('y')
    colormap(Water);
    colorbar;
    caxis([0 max(max(Vplot + (beta*Dv/(Dv-1)).*Uplot))]);
    set(gca,'LineWidth', 1);
    set(gca,'FontSize',10);
    set(gca,'FontWeight','bold');
    pbaspect(gca,[1 1 1]);view(0,90);shading interp;
    xlim([-x(end) x(end)]);ylim([-x(end) x(end)]);
    set(gca,'XTick',[-1 -0.8 -0.6 -0.4 -0.2 0 0.2 0.4 0.6 0.8 1]*x(end));
    set(gca,'YTick',[-1 -0.8 -0.6 -0.4 -0.2 0 0.2 0.4 0.6 0.8 1]*x(end));
    G2(i*dt) = getframe(figure(2)) ;drawnow;
    if mod(i*dt,100)==0 
    FigNamev = [folder,'\slicev',num2str(i*dt)];
    SolnNamev = [FigNamev,'.mat'];
    V=Vplot' + (beta*Dv/(Dv-1)).*Uplot';
    print(FigNamev,'-depsc2');
    save(SolnNamev,'x','y','V');
    end
    end
end
u_soln = w_old(1:2:nb-1); 
 v_soln = w_old(2:2:nb);
runtime = toc;

VideoName1 = [folder,'\Veg_Video'];
 % create the video writer with 1 fps
  writerObj1 = VideoWriter(VideoName1,'MPEG-4');
  writerObj1.FrameRate = 30;
  % set the seconds per image
% open the video writer
open(writerObj1);
% write the frames to the video
for n=1:length(G1)
    % convert the image to a frame
    frame1 = G1(n) ;    
    writeVideo(writerObj1, frame1);
end
% close the writer object
close(writerObj1);

VideoName2 = [folder,'\Water_Video'];
 % create the video writer with 1 fps
  writerObj2 = VideoWriter(VideoName2,'MPEG-4');
  writerObj2.FrameRate = 30;
  % set the seconds per image
% open the video writer
open(writerObj2);
% write the frames to the video
for n=1:length(G2)
    % convert the image to a frame
    frame2 = G2(n) ;    
    writeVideo(writerObj2, frame2);
end
% close the writer object
close(writerObj2);

Usoln = reshape(u_soln,steps,steps); 
  Vsoln = reshape(v_soln,steps,steps);
  Uplot = real(Usoln(:,:));
  Vplot = real(Vsoln(:,:));

    figure(1)
    surf(x,y,Uplot')
    title(['u at t=' num2str(tmax)]);
    xlabel('x')
    ylabel('y')
    colormap(Veg);
    colorbar;
    caxis([0 max(max(Uplot))]);
    set(gca,'LineWidth', 1);
    set(gca,'FontSize',10);
    set(gca,'FontWeight','bold');
    pbaspect(gca,[1 1 1]);view(0,90);shading interp;
    xlim([-x(end) x(end)]);ylim([-x(end) x(end)]);
    set(gca,'XTick',[-1 -0.8 -0.6 -0.4 -0.2 0 0.2 0.4 0.6 0.8 1]*x(end));
    set(gca,'YTick',[-1 -0.8 -0.6 -0.4 -0.2 0 0.2 0.4 0.6 0.8 1]*x(end));
    FigNameu = [folder,'\sliceu',num2str(tmax)];
    print(FigNameu,'-depsc2');    

    
    figure(2)
    surf(x,y,Vplot' + (beta*Dv/(Dv-1)).*Uplot')
    title(['v at t=' num2str(tmax)]);
    xlabel('x')
    ylabel('y')
    colormap(Water);
    colorbar;
    caxis([0 max(max(Vplot + (beta*Dv/(Dv-1)).*Uplot))]);
    set(gca,'LineWidth', 1);
    set(gca,'FontSize',10);
    set(gca,'FontWeight','bold');
    pbaspect(gca,[1 1 1]);view(0,90);shading interp;
    xlim([-x(end) x(end)]);ylim([-x(end) x(end)]);
    set(gca,'XTick',[-1 -0.8 -0.6 -0.4 -0.2 0 0.2 0.4 0.6 0.8 1]*x(end));
    set(gca,'YTick',[-1 -0.8 -0.6 -0.4 -0.2 0 0.2 0.4 0.6 0.8 1]*x(end));
    FigNamev = [folder,'\slicev',num2str(tmax)];
    print(FigNamev,'-depsc2');

function Fr = F(U,mu,f,g)
 Fr = zeros(nb,1);
 u = U(1:2:nb-1); v = U(2:2:nb);
 f1 = f(mu,u,v);
 f2 = g(mu,u,v);
 Fr(1:2:nb-1) = -f1; Fr(2:2:nb) = -f2;
end

end
