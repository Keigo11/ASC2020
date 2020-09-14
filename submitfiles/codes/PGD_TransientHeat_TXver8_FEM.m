function [X,T] = PGD_TransientHeat_TXver8_FEM(t,T_in,x,Max_terms,Max_fp_iter,epsilon,epsilon_tilde,Ft,Fx,k,q)
    %PGD solution of problem 4.1 with homogeneous Dirichlet boundary conditions
    %and homogeneous initial condition
    %For this example, we use simple second order FD with trapezoidal
    %integration rule along the space dimension and an implicit euler scheme
    %along the time dimension
    %Outputs:
    %       X,T :  computed PGD solution. The ith column of X (resp. T) contains the nodal values of X_i (resp. T_i)
    %Inputs:
    %       t,x: uniform 1D grids used along each dimension
    %       Max_terms: Maximum number of enrichments
    %       Max_fp_iter: Maximum number of fixed point iterations
    %       epsilon: termination criterion for the fixed point (Eq. (2.8))
    %       epsilon_tilde: termination criterion used for the enrichment process (Eq. (2.26))
    %       Ft,Fx: separated representation of the source term
    %       k: thermal diffusivity (const)
    %       q: heat flux (function of time)            
    %
    %Copyright (c) 2013, Francisco Chinesta (Ecole Centrale de Nantes), Roland Keunings (Universite catholique de Louvain), Adrien Leygue(CNRS)
    %Author: Adrien Leygue.
    %All rights reserved.
    %See License file for more details

    %thermal diffusivity
    % k=1;
    
    %Mesh definition for each dimension
    Nt = numel(t);
    Nx = numel(x);
    %reshape to ensure that x & t are column vectors
    t = t(:);
    x = x(:);
    %Mesh size
    dt = (t(2)-t(1));
    hx = (x(2)-x(1));
   
    %elemental stiffness matrices
    Kex = reshape([-1.0./hx,1.0./hx,1.0./hx,-1.0./hx],[2,2]);
%     Ket = reshape([1.0./dt,-1.0./dt,-1.0./dt,1.0./dt],[2,2]);
    
    %elemental mass matrices
    Mex = reshape([hx.*(1.0./3.0),hx.*(1.0./6.0),hx.*(1.0./6.0),hx.*(1.0./3.0)],[2,2]);
    Met = reshape([dt.*(1.0./3.0),dt.*(1.0./6.0),dt.*(1.0./6.0),dt.*(1.0./3.0)],[2,2]);
    
    %Mass matrix along each dimension
    Mx = sparse(Nx,Nx);
    Mt = sparse(Nt,Nt);
    %Stiffness matrix along each dimension
    Kx = sparse(Nx,Nx);
%     Kt = sparse(Nt,Nt);
    Ct = sparse(Nt,Nt);
    
    %Assembly of the FE matrices
    for el = 1:Nx-1
        Mx([el el+1],[el el+1]) = Mx([el el+1],[el el+1]) + Mex;
        Kx([el el+1],[el el+1]) = Kx([el el+1],[el el+1]) + Kex;
    end
    for el = 1:Nt-1
        Mt([el el+1],[el el+1]) = Mt([el el+1],[el el+1]) + Met;
        Ct(el+1,el+1) = 1;
        Ct(el,el+1) = -1;
    end
    
    Ct = Ct/2.0;
    Ct(1,1) = 0;
   
    %Boundary condition
    Tmax = max(T_in);
    Tmin = min(T_in);
        
    X = zeros(Nx,1);
    X(1,:) = 1;
    T = (T_in - Tmin)/(Tmax - Tmin);
    q = q/(Tmax - Tmin);%Normalization
    
    %main enrichment loop
    for term=2:Max_terms
        %initialization of the fixed point loop
        St = randn(Nt,1);
        Sx = randn(Nx,1);
        
        %Satisfaction of the homogeneous Dirichlet Boundary conditions for the
        %enrichments
%         if term==1
%             Sx(1) = 1;
%         else
%             Sx(1) = 0;
%         end
        
        Sx(1) = 0;
        St(1) = 0;
        
        %fixed point iterations
        for iter=1:Max_fp_iter
            %Store the old values of Sx & St for later comparison
            St_old = St;
            Sx_old = Sx;
            
            %Solve for Sx
            %construction of the boundary value problem along x
            %LHS coefficients
            alpha_x = k*St'*Mt*St;
            beta_x  = St'*Ct*St;
            
            %Source term coefficient
            ksi_x = Fx*St'*Mt*Ft;
            
            %Boundary condition
            mu_x = k*St'*Mt*q;
%             mu_x = St(end)*q;

            %Construction of the RHS
            RHS = Mx*ksi_x;
            
            %In case this is not the first enrichment, previous terms are added
            %to the RHS 
            if (term>1)
                %RHS coefficients
                delta_x_i = St'*Ct*T;
                gamma_x_i = k*St'*Mt*T;
                                
                RHS = RHS - (Mx*X)*delta_x_i' + (Kx*X)*gamma_x_i' ;
            end
            
            RHS(end) = RHS(end) + mu_x;
            
            %construction of the FD boundary value problem
            A = -alpha_x*Kx + beta_x*Mx;
            %solution with homogeneous boundary conditions
            Sx(2:end) = A(2:end,2:end)\RHS(2:end);
%             Sx(1:end) = A(1:end,1:end)\RHS(1:end);
            
            %Solve for St
            %construction of the initial value problem along t
            %LHS coefficients
            alpha_t = -k*Sx'*Kx*Sx;
            beta_t = Sx'*Mx*Sx;
            
            %Source term coefficient
            ksi_t = Ft*Sx'*Mx*Fx;
            
            %Boundary condition
            mu_t = k*Sx(end)*q;
            
            %Construction of the RHS
            RHS = Mt*(ksi_t+mu_t);
                      
            %In case this is not the first enrichment, previous terms are added
            %to the RHS
            if (term>1)
                %RHS coefficients
                delta_t_i = Sx'*Mx*X;
                gamma_t_i = -k*Sx'*Kx*X;
                
                
                RHS = RHS - (Ct*T)*delta_t_i' - (Mt*T)*gamma_t_i';
            end
            
%             RHS = RHS;
%             disp(mu_t)
            
            %construction of the initial value problem
            A = alpha_t*Mt + beta_t*Ct;
            %solution with homogeneous initial condition.
%             St(2:end) = A(2:end,2:end)\RHS(2:end);
            St(1:end) = A(1:end,1:end)\RHS(1:end);
                        
            %Norm of the difference between the 2 fixed point iterations  (Eq. (2.8))
            S_difference = sqrt(trapz(x,Sx.^2)*trapz(t,St.^2) + trapz(x,Sx_old.^2)*trapz(t,St_old.^2) - 2*trapz(x,Sx.*Sx_old)*trapz(t,St.*St_old));
            %fixed point exit test
            if(S_difference < epsilon), break; end
            
            
        end
        %New normalized enrichments are added to the  existing ones
        fact_x = sqrt(trapz(x,Sx.^2)/(x(end)-x(1))^2);
        fact_t = sqrt(trapz(St.*St)/(t(end)-t(1))^2);
        fact_xt = sqrt(fact_x*fact_t);
        X = [X fact_xt*Sx/fact_x];
        T = [T fact_xt*St/fact_t];
                
        %simplified stopping criterion (Eq. (2.26))
        E = sqrt(trapz(x,Sx.^2)*trapz(t,St.^2)) / sqrt(trapz(x,X(:,1).^2)*trapz(t,T(:,1).^2));
        
        disp(['ŒJ‚è•Ô‚µ”=', num2str(iter)]);
        
        if(E<epsilon_tilde), break; end
        
    end