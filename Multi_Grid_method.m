%多格点有限差分方法求解泊松方程

clear;
clc;
ni = 128;
L = 1;
dx = L/(ni-1);

x = (0:ni-1)'*dx;

A = 10;
k = 4;
b = A*sin(k*2*pi*x/L);

C1 = A/(k*2*pi/L);
C2 = -C1*L;
phi_true = (-A*sin(k*2*pi*x/L)*(1/(k*2*pi/L))^2+C1*x+C2);


%MGsolver
ni_c = ni/2;
dx_c = 2*dx;
solver_max = 100000;    %max iteration times


R_f = zeros(ni,1);
R_c = zeros(ni_c,1);
% eps_c = zeros(ni_c,1);
eps_f = zeros(ni,1);
phi = zeros(ni,1);

%start iteration
for it = 1:solver_max
    
   inner_its = 1;   %fine grid iteration time
   inner2_its = 50; %coarse grid iteration time
   w = 1.4;         %sor factor
%    phi = zeros(ni,1);
   % 1) perform one or more iterations on fine mesh
   phi(1,1) = phi(2,1);    %Neuman boundary condition
   for i = 2:ni-1
       
       g = 0.5*(phi(i-1,1)+phi(i+1,1)-dx*dx*b(i,1));
       phi(i,1) = phi(i,1)+w*(g-phi(i,1));
   end
    
   % 2) compute residue on the fine mesh, R = A*phi-b
   R_f(2:ni-1,1) = (phi(1:ni-2,1)-2*phi(2:ni-1,1)+phi(3:ni,1))/(dx*dx)-b(2:ni-1,1);
   R_f(1,1) = (phi(1,1)-phi(2,1))/dx;  %neumann boundary
   R_f(ni,1) = phi(ni,1)-0;            %Dirichlet boundary
   
   % 2b) check for termiantion
%    r_sum = 0;
   r_sum = sum(R_f.^2);
   norm = sqrt(r_sum)/ni;
   if norm < 1e-4
       fprintf('Converged after %d iterations at 2b)\n', it);
       plot(x,phi,'b');
%        legend('calculated');
       hold on
       plot(x,phi_true,'-r');
       legend('real','calculated');
       hold off
       break;
   end
   fprintf("norm = %f\n",norm);
   
%    % 3) restrict residue to the coarse mesh
%    R_c(0.5*((3:2:ni-2)+1),1) = 0.25*(R_f(2:2:ni-3,1)+R_f(3:2:ni-2,1)+R_f(4:2:ni-1,1));
%    R_c(1,1) = R_f(1,1);
%    
%    % 4) perform few iteration of the correction vector on the coarse mesh
%    for n = 1:inner2_its
%        eps_c = zeros(ni_c,1);
%        eps_c(1,1) = eps_c(2,1)+dx_c*R_c(1,1);
%        g = 0.5*(eps_c(1:ni_c-2,1)+eps_c(3:ni_c,1)-dx_c*dx_c*R_c(2:ni_c-1,1));
%        eps_c(2:ni_c-1,1) = eps_c(2:ni_c-1,1)+w*(g-eps_c(2:ni_c-1,1));
%    end
%    
%    % 5) interpolate eps to fine mesh
%    eps_f(3:2:ni-2,1) = eps_c(0.5*((3:2:ni-2)+1),1);
%    eps_f(2:2:ni-1,1) = 0.5*(eps_c(0.5*(2:2:ni-1),1)+eps_c(0.5*(2:2:ni-1)+1,1));
%    eps_f(1,1) = eps_c(1,1);
%    % 6) update solution on the fine mesh
%    phi(1:ni-1,1) = phi(1:ni-1,1)-eps_f(1:ni-1,1);
   
   
   
   
   
   
end