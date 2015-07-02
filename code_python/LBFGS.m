function [m_final misfit_data misfit_model,model]=LBFGS(m_initial,niteration,m)
% to test the performance of different algorithms in finding the minimum of
% the n-dimensional Rosenbrock valley function
% INPUT:
% m_initial : starting point
% niteration : maximum iterations
% m : how many previous iterations used for L-BFGS hessian approximation
% OUTPUT:
% m_final : final solution
% misfit_data : data space misfit
% misfit_model : model space misfit

%clear all
%close all

defval('m_initial',[-1.9 2]')
defval('niteration',30)
defval('m',1)

misfit_data=zeros(niteration+1,1);
misfit_model=zeros(niteration+1,1);
% dimensionality
n=length(m_initial);
% true solution
m_target=ones(n,1);

% line search parameters
step=10^(-4);  % step size
Nstep=10^6;   % Max steps

% iteration stop thresholds
threshold_model=10^(-6);
threshold_misfit=10^(-8);

disp 'L-BFGS method:'
        % intialization
        dm=zeros(n,m);
        dg=zeros(n,m);
        alpha=zeros(m,1);
        beta=zeros(m,1);
        s=zeros(n,1);
        y=zeros(n,1);
        
        % starting point
        m_current=m_initial;
        % initial misfit
        misfit0=func(m_initial);
        % iterations
        for iter=1:niteration
            % current misfit
            misfit_current=func(m_current)/misfit0;
            % data and model misfit
            misfit_data(iter)=misfit_current;
            misfit_model(iter)=norm(m_current-m_target);
            % current grad
            g_current=grad(m_current);
            % search direction
            if(iter==1) % steepest descent
                p_current=-g_current;
            else % L-BFGS
                % remember to calculate model change and grad change
                % shift one column to save storage
                for i=1:m-1
                    dm(:,i)=dm(:,i+1);
                    dg(:,i)=dg(:,i+1);
                end
                % store the last column
                dm(:,m)=m_current-m_old;
                dg(:,m)=g_current-g_old;
                % two-loop recursion to approximate Hessian
                q=g_current;
                % index
                j=m;
                for i=iter-1:-1:max(iter-m,1)
                    s=dm(:,j);
                    y=dg(:,j);
                    alpha(i)=dot(s,q)/dot(y,s);
                    q=q-alpha(i)*y;
                    j=j-1;
                end
                % preconditioning
                z=q;
                % index
                for i=max(iter-m,1):iter-1
                    j=j+1;
                    s=dm(:,j);
                    y=dg(:,j);
                    beta(i)=dot(y,z)/dot(y,s);
                    z=z+(alpha(i)-beta(i))*s;
                end
                
                % search direction
                p_current=-z;
            end
            
            % store current model, grad, search direction before updating
            m_old=m_current;
            g_old=g_current;
            
            % normalize the search direction
            norm_factor=norm(p_current);
            p_current=p_current/norm_factor;
            
            % line search
            for istep=1:Nstep
                % test model
                m_update=m_current+step*p_current;
                % test misfit
                misfit_update=func(m_update)/misfit0;
                if(misfit_update<misfit_current)
                    % update model
                    m_current=m_update;
                    % update misfit
                    misfit_current=misfit_update;
                else
                    break
                end
            end
            
            % stop creteria
            deltm=norm(m_current-m_old)/norm(m_current);
            if(deltm<=threshold_model)
                fprintf('iteration=%d: stop since model update %f smaller than threshold %f \n',...
                    iter,deltm,threshold_model);
                fprintf('final model = %f\n',m_current);
                fprintf('final misfit = %e\n',misfit_current);
                break;
            elseif(misfit_current<=threshold_misfit)
                fprintf('iteration=%d: stop since misfit %f smaller than threshold %f \n',...
                    iter,misfit_current,threshold_misfit);
                fprintf('final model = %f\n',m_current);
                fprintf('final misfit = %e\n',misfit_current);
                break;
            elseif(iter==niteration)
                fprintf('iteration=%d: stop since iteration reach the maximum %d \n',...
                    iter,niteration);
                fprintf('final model = %f\n',m_current);
                fprintf('final misfit = %e\n',misfit_current);
            end
        end
        % final solutions
        m_final=m_current;
        misfit_data(iter+1)=misfit_current;
        misfit_model(iter+1)=norm(m_final-m_target);



