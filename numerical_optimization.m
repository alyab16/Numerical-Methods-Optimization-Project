function [rcnd, x0, re0, rr0, ref] = numerical_optimization(n)

    rcnd=0;
    x0=zeros(n,n);
    re0=0;
    rr0=0;
    ref=0;

    %%%%%%%%%%%%%%%%%%%%%%%%%% DEFINING A, b, x %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % The matrix A as defined in the assignment
    A=zeros(n,n);
    for i = 1:n
        for j = 1:n
            A(i,j) = j^i;
        end
    end
    
    % The vector b as defined in the assignment
    b=zeros(n,1);
    for i = 1:n
        for j = 1:n
            b(i,1) = b(i,1)+(-1)^(j+1)*A(i,j);
        end
    end
    
    %True Solution to the system above
    x= zeros(n,1);
    for i = 1:n
        x(i,1) = (-1)^(i+1);
    end


    % estimate the reciprocal condition of A
    rcnd = rcond(A);
    
    %%%%%%%%%%%%%%%%%%%%%%% THE INITIAL SOLUTION X_%%%%%%%%%%%%%%%%%%%%%%%%
    % PA=LU Factorization
    [L,U,P] = lu(A);
    
    %Forward Solve
    d = inv(L)*(P*b);    
    
    %Backward Solve
    x0 = inv(U)*d;
        
    % The relative error in the initial solution
    re0 = norm(x-x0)/norm(x);

    % The initial relative residual
    rr0 = norm(b-(A*x0))/ norm(b);

    % Checking for Subtractive Cancellation In The Initial Solution x0
    check = isequal(b-(A*x0), zeros(n, 1));
    
    %%%%%%%%%%%%%%%%%%%%%% ITERATIVE IMPROVEMENT %%%%%%%%%%%%%%%%%%%%%%%%%%
    if check==0
        for i = 0:4
            if i == 0
                xi = x0;
              x_i_minus_1 = x0;
            end
        
            % Compute r_i= b -(A*xi)
            r_i = b - (A*xi);
            z_i = inv(A)*r_i;
            checking_value = norm(xi - x_i_minus_1);
            x_i_minus_1 = xi;
            xi = xi + z_i;
            

            
        %%%%%%%%%% Checking if iterative improvement is useless %%%%%%%%%%
            if (i>0) 
                
                if not(norm(xi - x_i_minus_1) < checking_value)
                    break
                end
            end
        end
        % The relative error in the final solution
        ref = norm(x-xi)/norm(x);
    end
end