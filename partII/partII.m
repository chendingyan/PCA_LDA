function [] = partII()

    % generate the data

    rng(1); 
    r = sqrt(rand(100,1)); 
    t = 2*pi*rand(100,1);  
    data1 = [r.*cos(t), r.*sin(t)]; 

    r2 = sqrt(3*rand(100,1)+1); 
    t2 = 2*pi*rand(100,1);      
    data2 = [r2.*cos(t2), r2.*sin(t2)]; 

    % plot the data

    figure;
    plot(data1(:,1),data1(:,2),'r.','MarkerSize',15)
    hold on
    plot(data2(:,1),data2(:,2),'b.','MarkerSize',15)
    axis equal
    hold on

    % work on class 1
    [a1, R1] = calcRandCentre(data1);

    fprintf('Class 1 center point[%f, %f]\n',a1(1), a1(2))
    fprintf('Class 1 Radii is %f\n', R1)
    % work on class 2
    [a2, R2] = calcRandCentre(data2);
    fprintf('Class 2 center point[%f, %f]\n',a2(1), a2(2))
    fprintf('Class 2 Radii is %f\n', R2)
    % plot centre and radius for class 1
    plot(a1(1), a1(2), 'rx', 'MarkerSize', 15);
    viscircles(a1', R1, 'Color', 'r', 'LineWidth', 1);
    hold on

    % plot centre and radius for class 2
    plot(a2(1), a2(2), 'bx', 'MarkerSize', 15);
    viscircles(a2', R2, 'Color', 'b', 'LineWidth', 1);

end

function [a, R] = calcRandCentre(data)

    % to be completed

    % Formulate the result from previous question, 
    % note that on my report, it is max.. but here quadprog solve a min
    % problem, so I need change my function to -(...)
    % quadprog(H,f,A,b,Aeq,beq,lb,ub,x0)
    % H, A, and Aeq are matrices, and f, b, beq, lb, ub, and x are vectors.
    % First, calculate H, which is K in my illustration and 
    % it is x_i .* x_j
    n = length(data);
    H = zeros(n);
    for i = 1:n
        for j = 1:n
            H(i,j) = data(i,:)*transpose(data(j,:));
        end
    end
    % Then f, which is x_i .* x_i in max problem, should change to -f^T
    f = zeros(n,1);
    for i = 1:n
        f(i,1)=data(i,:)*transpose(data(i,:));
    end
    f = -transpose(f);
    % other simple parameter, we assume C = 0.4
    C = 0.4;
    A = zeros(n);
    b = zeros(n,1);
    Aeq = ones(1,n);
    beq = 1;
    lb = zeros(n,1);
    ub = C*ones(n,1);
    x = quadprog(H,f,A,b,Aeq,beq,lb,ub);
    %now we get \mu and we can calculate center point and radii
    a_t =  transpose(x)*data/sum(x);
    a = transpose(a_t);
    distance = data-a_t;
    R_sum = 0;
    num = 0;
    for i = 1:length(x)
        % 0 <= x <= C, here we change 0 to a small number as a threshold
        % and get better result
        if x(i)>=0.00000006 && x(i)<= C
            % calculate the sum of radii, each radii is i * i
            R_sum = R_sum + sqrt(distance(i,:)*transpose(distance(i,:)));
            num=num+1;
        end
    end
    R = R_sum/num;

end