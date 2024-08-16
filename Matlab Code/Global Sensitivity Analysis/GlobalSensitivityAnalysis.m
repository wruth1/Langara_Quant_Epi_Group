%% 

%% Initialize data

numExp = 100;

numInput = 18;
numOutput = 9;


% A and B are matrices of data
A = zeros(numExp, numInput);
B = A;

for k=1:numExp
    A(k,:) = createBioparametersUniform();
    B(k,:) = createBioparametersUniform();
end

%% load initial data
TFP0 = 6186950;

X0 =    2.435279313142126e+06;
E0 = 1.484646245171584e+04;

L0 = 3.096493365566524e+06;
% T0 = 1.070931823865891e+03;
T0 = 1081;
R0 =      6.392498588396340e+05;
TFP0 = X0 + E0+L0+T0+R0;

IC = [X0, E0, L0, T0, R0];

%% 
% save outputs
YA = zeros(numExp, numOutput);
YB = YA;
YA_other = cell(numExp,1);
YB_other = cell(numExp,1);

for i=1:numExp
    paramsi = A(i,:);


    
    
    BPi = paramsi(1:11);
    ICi = paramsi(12:16);



    TFP0 = sum(ICi);

    % x_magnitude, the amount we will rescale by
    x_m = [BPi(9) BPi(10) BPi(11) BPi(8), ICi];
    

    
    % update initial conditions

    % % find initial conditions R0, E0, L0, T0, R0. use steady state
    ysteady = findSteadyState2(BPi, ICi, ReportedImmigration(1));
    % rescale to correct total population.  
    ysteady = ysteady*TFP0/sum(ysteady);
    localIC = ysteady;

    
    allParamsi = [BPi localIC paramsi(17:18)];



    
    % set up fmincon input
    x0 = x_m;

    
    numx=9;
    Aineq = zeros(1,numx);
    Aineq(1,1:3) = x_m(1:3); % q1+q2+q3 <= 100%
    bineq(1) = 1 ;


    Aeq = zeros(2,numx);
    Aeq(1,end-4:end) = x_m(end-4:end); % total foreign-born population
    
    beq = [TFP0];

    %T0 = 1081
    Aeq(2,1:numx) = [0 0 0 0 0 0 0 x_m(8) 0];
        prev_ampi = paramsi(17);

    beq(2) = [1081.*prev_ampi];


 
    lb = zeros(1, numx); 
    ub = lb + Inf;
    
    

    f2=@(y)IncidenceError6_rescale(y,allParamsi, x_m, ReportedImmigration, ReportedIncidence);
    y0 = x0./x_m; %q1 q2 q3

    
    [ymin2,fval2,exitflag,output,lambda,grad,hessian] = fmincon(f2, y0 , Aineq , bineq, Aeq, beq, lb, ub) ; % about 10 seconds to run

    xmin2 = x_m.*ymin2;
    
    YA(i,:) = xmin2;
YA_other{i} = {fval2,exitflag,output,lambda,grad,hessian};
end

% YB
for i=1:numExp
    paramsi = B(i,:);


    
    
    BPi = paramsi(1:11);
    ICi = paramsi(12:16);



    TFP0 = sum(ICi);

    % x_magnitude, the amount we will rescale by
    x_m = [BPi(9) BPi(10) BPi(11) BPi(8), ICi];
    

    
    % update initial conditions

    % % find initial conditions R0, E0, L0, T0, R0. use steady state
    ysteady = findSteadyState2(BPi, ICi, ReportedImmigration(1));
    % rescale to correct total population.  
    ysteady = ysteady*TFP0/sum(ysteady);
    localIC = ysteady;

    
    allParamsi = [BPi localIC paramsi(17:18)];



    
    % set up fmincon input
    x0 = x_m;

    
    numx=9;
    Aineq = zeros(1,numx);
    Aineq(1,1:3) = x_m(1:3); % q1+q2+q3 <= 100%
    bineq(1) = 1 ;


    Aeq = zeros(2,numx);
    Aeq(1,end-4:end) = x_m(end-4:end); % total foreign-born population
    
    beq = [TFP0];

    %T0 = 1081
    Aeq(2,1:numx) = [0 0 0 0 0 0 0 x_m(8) 0];
        prev_ampi = paramsi(17);

    beq(2) = [1081.*prev_ampi];


 
    lb = zeros(1, numx); 
    ub = lb + Inf;
    
    

    f2=@(y)IncidenceError6_rescale(y,allParamsi, x_m, ReportedImmigration, ReportedIncidence);
    y0 = x0./x_m; %q1 q2 q3

    
    [ymin2,fval2,exitflag,output,lambda,grad,hessian] = fmincon(f2, y0 , Aineq , bineq, Aeq, beq, lb, ub) ; % about 10 seconds to run

    xmin2 = x_m.*ymin2;
    
    YB(i,:) = xmin2;
    YB_other{i} = {fval2,exitflag,output,lambda,grad,hessian};

end

%% now compute YC
YC = zeros(numExp, numOutput, numInput);
YC_other = cell(numExp,numInput);

for j=1:numInput
    % C is B, but with jth column replaced by A
    C = B;
    C(:,j) = A(:,j);

    for i=1:numExp
        paramsi = C(i,:);
    
    
        
        
        BPi = paramsi(1:11);
        ICi = paramsi(12:16);
    
    
    
        TFP0 = sum(ICi);
    
        % x_magnitude, the amount we will rescale by
        x_m = [BPi(9) BPi(10) BPi(11) BPi(8), ICi];
        
    
        
        % update initial conditions
    
        % % find initial conditions R0, E0, L0, T0, R0. use steady state
        ysteady = findSteadyState2(BPi, ICi, ReportedImmigration(1));
        % rescale to correct total population.  
        ysteady = ysteady*TFP0/sum(ysteady);
        localIC = ysteady;
    
        
        allParamsi = [BPi localIC paramsi(17:18)];
    
    
    
        
        % set up fmincon input
        x0 = x_m;
    
        
        numx=9;
        Aineq = zeros(1,numx);
        Aineq(1,1:3) = x_m(1:3); % q1+q2+q3 <= 100%
        bineq(1) = 1 ;
    
    
        Aeq = zeros(2,numx);
        Aeq(1,end-4:end) = x_m(end-4:end); % total foreign-born population
        
        beq = [TFP0];
    
        %T0 = 1081
        Aeq(2,1:numx) = [0 0 0 0 0 0 0 x_m(8) 0];
            prev_ampi = paramsi(17);
    
        beq(2) = [1081.*prev_ampi];
    
    
     
        lb = zeros(1, numx); 
        ub = lb + Inf;
        
        
    
        f2=@(y)IncidenceError6_rescale(y,allParamsi, x_m, ReportedImmigration, ReportedIncidence);
        y0 = x0./x_m; %q1 q2 q3
    
        
        [ymin2,fval2,exitflag,output,lambda,grad,hessian] = fmincon(f2, y0 , Aineq , bineq, Aeq, beq, lb, ub) ; % about 10 seconds to run
    
        xmin2 = x_m.*ymin2;
        
        YC(i,:,j) = xmin2;
        YC_other{i,j} = {fval2,exitflag,output,lambda,grad,hessian};
    
    end
end
%% 
save('SensitivityIndices.mat')

save SensitivityIndices2.mat A B C YA YB YC YA_other YB_other YC_other;
%% 

