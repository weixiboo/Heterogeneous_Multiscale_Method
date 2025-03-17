
function  [delx,dely,delt,m_delx,m_dely,X_dual,Y_dual,H_z_new] = hmm_first_order_2d(k_mac,rho_constant,rho_div,m_Nx)
%[X_dual_Y_main,Y_main_X_dual,E_x_new]
%[X_main_Y_dual,Y_dual_X_main,E_y_new]
%[X_dual,Y_dual,H_z_new]

mu0=1;       % Free space permeability
%eps0 defined in hmm_E_flux
c = sqrt(1);

ND=20; delx=1/ND;     % Avoid dispersion   
Nx= round(1/delx);
Ny = Nx;
delx = 1/Nx;
dely = delx;

Final_T = 0.5;
delt= (0.9/c)*(1/sqrt((1/delx^2)+(1/dely^2))); 
Nt = round((Final_T)/delt); 
delt = (Final_T)/Nt;

delx = delx/k_mac;
dely = delx;
delt = delt/k_mac;
Nt = Nt*k_mac;
Nx = Nx*k_mac;
Ny = Ny*k_mac;


%rho shrink with delx
%micro_box = 1.25e-4;
% micro_box = delx/16;

if rho_constant == 0
    micro_box = delx/rho_div;
else
    micro_box = rho_constant/rho_div;
end

m_Nx = m_Nx/2;
m_Ny = m_Nx;

m_delx = micro_box/(2*m_Nx);
m_dely = m_delx;


if(micro_box > delx)
    warning('microbox > delx or microbox time > delt')
end


%X's and T's
x = 0:delx:Nx*delx;
y = 0:dely:Ny*dely;
x_dual = x(2:end) - delx/2;
y_dual = y(2:end) - dely/2;

%[X,Y] = meshgrid(x,y);
[X_dual,Y_dual] = meshgrid(x_dual,y_dual);
[X_dual_Y_main,Y_main_X_dual] = meshgrid(x_dual,y);
[X_main_Y_dual,Y_dual_X_main] = meshgrid(x,y_dual);


%create dictionary to pass parameters more easily
macro_grids = containers.Map({'x','y','x_dual','y_dual','X_dual','Y_dual',...
    'X_dual_Y_main','Y_main_X_dual',...
    'X_main_Y_dual','Y_dual_X_main'},...
    {x,y,x_dual,y_dual,X_dual,Y_dual,...
    X_dual_Y_main,Y_main_X_dual,...
    X_main_Y_dual,Y_dual_X_main});

step_size = containers.Map({'delx','dely','delt',...
    'm_delx','m_dely',...
    'm_Nx','m_Ny'},...
    {delx,dely,delt,...
    m_delx,m_dely,...
    m_Nx,m_Ny});

micro_grids = micro_grid(macro_grids,step_size);

%Initial Data
H_z_old = zeros(size(X_dual));
E_x_old = zeros(size(Y_main_X_dual));
E_y_old = zeros(size(Y_dual_X_main));

%Material parameters
a = 0.0005;
mat_vol_frac = 0.5;
func = @(x) vol_frac(x,a) - mat_vol_frac;

if mat_vol_frac ~= 1
    r = fzero(func,[0,a/sqrt(2)]);
else
    r = a/sqrt(2);
end

[Ex_x_ind,Ex_y_ind,eps_inf_mat_Ex] = epsilon_inf_mat(x_dual,y,a,r,m_delx,m_dely,m_Nx,m_Ny,'Ex');
[Ey_x_ind,Ey_y_ind,eps_inf_mat_Ey] = epsilon_inf_mat(x,y_dual,a,r,m_delx,m_dely,m_Nx,m_Ny,'Ey');

material_parameters = containers.Map({'Ex_x_ind','Ex_y_ind','eps_inf_mat_Ex',...
    'Ey_x_ind','Ey_y_ind','eps_inf_mat_Ey','a'},...
    {Ex_x_ind,Ex_y_ind,eps_inf_mat_Ex,Ey_x_ind,Ey_y_ind,eps_inf_mat_Ey,a});

%Source
source = b((delt/2):delt:(Final_T-delt/2),4);
src_loc = round(size(X_dual,1)*4/10);


for i = 1:Nt
    if(mod(i,round(Nt/100)) == 0)
        disp(round(i/round(Nt/100)))
        
        % figure(1)
        % im = imagesc(x_dual,y_dual,H_z_old,[-5 5]);
        % im.AlphaData = 0.8;
        % hold on
        % im2 = imagesc(x_dual,y,eps_inf_Ex);
        % im2.AlphaData = 0.25;
        % hold off
        % pause(0.01)

        %max(H_z_old,[],'all')
    end

        % figure(1)
        % imagesc(x_dual,y_dual,H_z_old,[-5 5]);
        % pause(0.01)

    %E is at main grid in space, dual grid in time
    %H is at dual grid in space, main grid in time
    %Frist loop calculalte E at time 1/2 and H at time 1

    macro_data = containers.Map({'E_x','E_y','H_z'},{E_x_old,E_y_old,H_z_old});

    [Ex_t,Ey_t] = hmm_E_flux(macro_data,macro_grids,micro_grids,step_size,material_parameters); 

    E_x_new = E_x_old + delt.*Ex_t;
    E_y_new = E_y_old - delt.*Ey_t;
    H_z_new = H_z_old - delt./mu0.*(diff(E_y_new,1,2)./delx - diff(E_x_new,1,1)./dely);

    H_z_new(src_loc:src_loc+1,src_loc:src_loc+1) = ...
        H_z_new(src_loc:src_loc+1,src_loc:src_loc+1)...
        + (delt)./(delx*dely)./mu0.*source(i)/2;
    
    %Update variables
    E_x_old = E_x_new;
    E_y_old = E_y_new;
    H_z_old = H_z_new;

    

end



end

function [x_ind_func,y_ind_func,out] = epsilon_inf_mat(x,y,a,r,m_delx,m_dely,m_Nx,m_Ny,Ex_or_Ey)
    
    x_ind = unique(ceil(mod(x,a)/m_delx));
    y_ind = unique(ceil(mod(y,a)/m_dely));

    x_ind_func = @(x) find(x_ind == ceil(mod(x,a)./m_delx),1);
    y_ind_func = @(y) find(y_ind == ceil(mod(y,a)./m_delx),1);

    switch Ex_or_Ey
        case 'Ex'
            out = zeros(2*m_Ny - 1,2*m_Nx,length(x_ind),length(y_ind));
        
            for i = 1:length(x_ind)
                for j = 1:length(y_ind)
                    center_x = (((x_ind(i)-1)*m_delx)+m_delx/2);
                    x = (center_x + m_delx/2 - m_delx*m_Nx):m_delx:(center_x - m_delx/2 + m_delx*m_Nx);
        
                    center_y = (((y_ind(j)-1)*m_dely)+m_dely/2);
                    y = (center_y - m_dely*(m_Ny-1)):m_dely:(center_y + m_dely*(m_Ny-1));
        
                   out(:,:,i,j) = epsilon_inf(x,y,r,a,m_delx,m_dely);
                end
            end

        case 'Ey'
            out = zeros(2*m_Ny,2*m_Nx-1,length(x_ind),length(y_ind));
        
            for i = 1:length(x_ind)
                for j = 1:length(y_ind)
                    center_x = (((x_ind(i)-1)*m_delx)+m_delx/2);
                    x = (center_x - m_delx*(m_Nx-1)):m_delx:(center_x+ m_delx*(m_Nx-1));
        
                    center_y = (((y_ind(j)-1)*m_dely)+m_dely/2);
                    y = (center_y + m_dely/2 - m_dely*m_Ny):m_dely:(center_y - m_dely/2 + m_dely*m_Ny);
        
                   out(:,:,i,j) = epsilon_inf(x,y,r,a,m_delx,m_dely);
                end
            end
    end

end


function out = epsilon_inf(x,y,r,a,delx,dely)

    [X,Y] = meshgrid(x,y);
    
    out = zeros(size(X));

    for i = [-1,1]
        for j = [-1,1]

            D = sqrt((mod(X+i*delx/4,a)-a/2).^2 + (mod(Y+j*dely/4,a)-a/2).^2);

            %D = max(abs(mod(X+i*delx/4,a)-a/2),abs(mod(Y+j*dely/4,a)-a/2));

            out(D<=r) = out(D<=r) + 2.7;
            out(D>r) = out(D>r) + 1.03;

        end
    end

    out = out/4;


    % out(X<=(0.5)) = 2.7;
    % out(X>(0.5)) = 1.03;
    % out(X==0.5) = mean([2.7,1.03]);
    
end

function out = vol_frac(r,a)
    
    out = length(r);

    in_ind = r<=a/2;
    out_ind = r>a/2;

    out(in_ind) = (pi*r(in_ind).^2)./(a^2);

    b = sqrt(r(out_ind).^2-(a.^2./4));
    out(out_ind) = (pi*r(out_ind).^2 - (4*r(out_ind).^2).*atan(2*b./a) + 2.*a.*b)./(a^2) ;

end

function out = b(t,f_bw)

    an = [0.353222222,-0.488,0.145,-0.010222222];
    out = zeros(size(t));
    ind = t > 0 & t < 1.55/f_bw;

    for i = 1:length(an)
        %out(ind) = out(ind) + an(i).*cos(2*pi*(i-1)*(f_bw/1.55).*t(ind));
        out(ind) = out(ind) - (i-1).*an(i).*sin(2*pi*(i-1)*(f_bw/1.55).*t(ind));
    end

end

function micro_grid = micro_grid(macro_grids,step_size)

    m_Nx = step_size('m_Nx');
    m_Ny = step_size('m_Ny');
    m_Nt = 0;
    m_delx = step_size('m_delx');
    m_dely = step_size('m_dely');

    x = macro_grids('x');
    y = macro_grids('y');
    x_dual = macro_grids('x_dual');
    y_dual = macro_grids('y_dual');


    Ex_mx = zeros([length(x_dual), 1+ 2*m_Nx + 2*m_Nt]);
    Ex_my = zeros([length(y), 1+ 2*m_Ny + 2*m_Nt]);

    Ex_mx_dual = zeros([length(x_dual), 2*m_Nx + 2*m_Nt]);
    Ex_my_dual = zeros([length(y), 2*m_Ny + 2*m_Nt]);
    
    for i = 1:length(x_dual)
            Ex_mx(i,:) = (x_dual(i) - m_delx*m_Nx - m_delx*m_Nt):...
                m_delx:(x_dual(i) + m_delx*m_Nx + m_delx*m_Nt);
            Ex_mx_dual(i,:) = (x_dual(i) + m_delx/2 - m_delx*m_Nx - m_delx*m_Nt):...
                m_delx:(x_dual(i) - m_delx/2 + m_delx*m_Nx + m_delx*m_Nt);
    end

    for i = 1:length(y)
            Ex_my(i,:) = (y(i) - m_dely*m_Ny - m_dely*m_Nt):...
                m_dely:(y(i) + m_dely*m_Ny + m_dely*m_Nt);
            Ex_my_dual(i,:) = (y(i) + m_dely/2 - m_dely*m_Ny - m_dely*m_Nt):...
                m_dely:(y(i) - m_dely/2 + m_dely*m_Ny + m_dely*m_Nt);
    end

    
    Ey_mx = zeros([length(x), 1+ 2*m_Nx + 2*m_Nt]);
    Ey_my = zeros([length(y_dual), 1+ 2*m_Ny + 2*m_Nt]);

    Ey_mx_dual = zeros([length(x), 2*m_Nx + 2*m_Nt]);
    Ey_my_dual = zeros([length(y_dual), 2*m_Ny + 2*m_Nt]);
    
    for i = 1:length(x)
            Ey_mx(i,:) = (x(i) - m_delx*m_Nx - m_delx*m_Nt):...
                m_delx:(x(i) + m_delx*m_Nx + m_delx*m_Nt);
            Ey_mx_dual(i,:) = (x(i) + m_delx/2 - m_delx*m_Nx - m_delx*m_Nt):...
                m_delx:(x(i) - m_delx/2 + m_delx*m_Nx + m_delx*m_Nt);
    end

    for i = 1:length(y_dual)
            Ey_my(i,:) = (y_dual(i) - m_dely*m_Ny - m_dely*m_Nt):...
                m_dely:(y_dual(i) + m_dely*m_Ny + m_dely*m_Nt);
            Ey_my_dual(i,:) = (y_dual(i) + m_dely/2 - m_dely*m_Ny - m_dely*m_Nt):...
                m_dely:(y_dual(i) - m_dely/2 + m_dely*m_Ny + m_dely*m_Nt);
    end

    micro_grid = containers.Map({'Ex_mx','Ex_mx_dual','Ex_my','Ex_my_dual',...
        'Ey_mx','Ey_mx_dual','Ey_my','Ey_my_dual'},...
        {Ex_mx,Ex_mx_dual,Ex_my,Ex_my_dual,...
        Ey_mx,Ey_mx_dual,Ey_my,Ey_my_dual});
end
