%% **************************************************************
% A simple script to solve the 1D equation for the
% upper overturning cell, as discussed in Jansen (submitted to JPO)
%
% Written by Malte Jansen (mfj@uchicago.edu) in May/June of 2017 
%
% Notice that in the main loop everything is non-dimensionalized
% using the depth of  the cell, H, (which is a parameter that
% is solved for), and the Coriolis parameter, f.
%
% The script includes a vertically constant and a vertically
% varying  diffusivity profile, but as of now that needs to be
% changed by   commenting out (or uncommenting) the definitions
% of kappa and d(kappa)/dz in the main loop.
% The solver converges less well for vertically varying
% diffusivity profiles. (I believe this is may be due to a
% positive feed-back for downward increasing diff. profiles.)
% It often helps to simply change the initial guess a bit.
% However, for some diffusivity profiles and parameters I haven't
% been able to obtain sensible solutions at all. And in some of
% these cases there may in fact not be a solution; but as of now
% I have no means to tell whether or not a solution exists.
% (This is a challenge for the Mathematicians out there...;)
% 
% As always with code, there may well be errors, and it's
% not recommended to use and trust this as a "black box"
%
% 
% If you find any errors, please get in touch!
% If you make improvements, I would also be happy if you get in touch.
%
% I'm also hoping to write (or even better: have somebody write ;)
% a Python version of this. If you are interested in that, please
% get in touch.
%
%% ***************************************************************

clear;

f=1.2e-4; % Coriolis parameter (in the north of basin)
bs=0.025; % surface buoyancy in the basin
% kappa is the value used for constant diffusivity:
kappa=6e-5; 
% the folowing three are parameters for 
% vertically varying diffusivity
% add or remove comments in main loop to use constant or varying diff.
kappa_back=1e-5;
kappa_s=3e-5;
kappa_4k=3e-4;
a=6.37e6;   % planetary radius (used in basin area calculation)
% Area over which diapycnal upwelling occurs:
A=2*pi*a^2*59/360*(sind(69)-sind(-48));
% Diffusive buoyancy loss to abyss (vector to compute multiple cases):
Bint=[3e3,1.2e4,3e3,1.2e4];  
% Depth of SO upwelling (vector to compute multiple cases):
Hmaxso=[2000,2000,1500,1500]; 
% Maximum value of SO upwelling (in SV):
psimaxso=4;

 
figure(1) % initialize figure
clf
zz = linspace(-1,0,100); % this is just for plotting

N=length(Bint); % number of cases to be computed

for i=1:N % this loop is over the differnt "cases" to be considered
    
% For constant kappa, use the next two lines
kap=@(z,H)kappa/H^2/f;% nondimensional kappa
dkapdz=@(z,H)0;
% For vertically varying kappa, use the next two lines
%kap=@(z,H) 1/H^2/f*(kappa_back+kappa_s*exp(z*H/100)+kappa_4k*exp(-z*H/1000-4));% nondim kappa
%dkapdz=@(z,H) 1/H/f*(      kappa_s/100*exp(z*H/100)-kappa_4k/1000*exp(-z*H/1000-4));%nondimensional d(kappa)/dz  

a=@(z,H) 1/(kap(z,H)*A)*H^2; % 1/(kappa A) nondimensionalized 

% nondim. surface buoyancy: (still needs to be div. by H, which 
% will be done in BC itself, so it doesn't need to be function)
b=-bs/f^2; 
% (Nondim.) stratification at bottom of cell:
bz=@(H)Bint(i)/(A*kap(-1,H))/f^3/H^2;

% Southern ocean upwelling:
psiso=@(z,H) psimaxso*1e6/f/H^3*sin(-pi*max(z*H,-Hmaxso(i))/Hmaxso(i)).^2;

% Initial conditions for solver:
% if the solver doesn't converge, it often helps to play a little with the initial guess
% the number of z points here specifies (among other things) the minimum number of mesh points used
solinit = bvpinit(linspace(-1,0,50),[1,0,-0.1,-bz(Hmaxso(i)+200)],Hmaxso(i)+200);

% This is the actual differential equation we are solving
% (see Jansen ????)
% (y(1)=psi,y(2)=psi_z,y(3)=psi_zz,y(4)=psi_zzz) 
ode = @(z,y,H)[y(2);y(3);y(4);...
               a(z,H)*(y(1)- psiso(z,H) - A/H^2*dkapdz(z,H))*y(4)];        

% Boundary conditions for differential equation:
% (Notice that the eq. is 4th order, but we are also solving for H
% hence 4 BCs are needed)
bc = @(ya,yb,H)[ ya(1); ya(2);ya(4)+bz(H); yb(1) ;yb(3)-b/H]; % fixed bottom b_z,psi, and psi_z

% This is where the equation actually gets solved:
sol = bvp4c(ode,bc,solinit);

% evaluate solution at points in z for plotting
y = deval(sol,zz);
H=sol.parameters % H is the depth of the overturning cell


% This is for plotting:
if i==1
    line(y(1,:)*H^3*f/1e6,zz*H,'linestyle','--','color',[Bint(i)./max(Bint),0.5*Bint(i)./max(Bint),1-Bint(i)./max(Bint)]); hold on;
    ax1=gca;
    xlabel('$\Psi_{N} \; [SV]$','interpreter','Latex','fontsize',20)
    ylabel('depth $[m]$','interpreter','Latex','fontsize',20)
    ax1_pos = ax1.Position; % position of first axes
    set(ax1,'xlim',[-5 20],'ylim',[-3000 0],'fontsize',15)
    ax2 = axes('Position',ax1_pos,'XAxisLocation','top',...
    'YAxisLocation','right','XTick',(-0.01:0.01:0.04),...
    'Color','none','xlim',[-0.01 0.04],'ylim',[-3000 0],'fontsize',15);
    xlabel(ax2,'$b \; [ms^{-2}]$','interpreter','Latex','fontsize',20)
else
    if i<3
    line(y(1,:)*H^3*f/1e6,zz*H,'parent',ax1,'linestyle','--','color',[Bint(i)./max(Bint),0.5*Bint(i)./max(Bint),1-Bint(i)./max(Bint)]); hold on;
    else
    line(y(1,:)*H^3*f/1e6,zz*H,'parent',ax1,'linestyle','-.','color',[Bint(i)./max(Bint),0.5*Bint(i)./max(Bint),1-Bint(i)./max(Bint)]); hold on;
    end
end
if i<3
    line(-y(3,:)*f^2*H,zz*H,'parent',ax2,'linestyle','-','color',[Bint(i)./max(Bint),0.5*Bint(i)./max(Bint),1-Bint(i)./max(Bint)]);
else
    line(-y(3,:)*f^2*H,zz*H,'parent',ax2,'linestyle',':','color',[Bint(i)./max(Bint),0.5*Bint(i)./max(Bint),1-Bint(i)./max(Bint)]);
end

end % of main loop

% Add legend to plot (if desired). Note that labelling currently not automated!
legend(ax1,{'$\mathcal{B}=3\times10^3$m$^4$s$^{-3}$, $H_{S}=2\,$km',...
            '$\mathcal{B}=1.2\times10^4$m$^4$s$^{-3}$, $H_{S}=2\,$km',...
            '$\mathcal{B}=3\times10^3$m$^4$s$^{-3}$, $H_{S}=1.5\,$km',...
            '$\mathcal{B}=1.2\times10^4$m$^4$s$^{-3}$, $H_{S}=1.5\,$km',...
            },'Location','Southeast','interpreter','Latex','fontsize',16);
legend(ax1,'boxoff')

%% This can be used to plot diffusivity profile:
% H=4000;
% figure(99)
% clf
% semilogx(kap(zz,H)*H^2*f,zz*H,'k')
% set(gca,'fontsize',15,'xlim',[1e-5,5e-4],'xtick',[1e-5,1e-4])
% xlabel('$\kappa(z)$ [m$^2$s$^{-1}$]','interpreter','Latex','fontsize',20)
% ylabel('depth [m]','interpreter','Latex','fontsize',20)

