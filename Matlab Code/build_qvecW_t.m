function qvec_t = build_qvecW_t(piW_t,LTBIW, INFW, qLqR_ratio)
%{
INPUT:
    pi_t in R^numt
    LTBIW in R^6, 6 WHO geographic regions. Houben's 
    INFW in R^6
    qLqR_ratio scalar
%}

numq = 4;
numw = 6;

numt=size(piW_t,1);

qvec_t = zeros(numt,numq); 

for t = 1:numt
    qvec_local = zeros(1,numq);

    

    for w=1:numw
        
        % solve Ax = b
        A = zeros(4);
        b = zeros(4,1);
    
        % E + L + R = 22.4% * pi_w(2014)
        A(1,:) = [0 1 1 1];
        b(1) = LTBIW(w)*piW_t(t,w);
    
        % E/(X+E+L+R) = 1.31%
        A(2,:) = [0 1 0 0 ] - [1 1 1 1]*INFW(w);
        b(2) = 0;
  
        % L/R = 46%
        A(3,:) = [0 0 1 -qLqR_ratio];
        b(3)=0;
    
        % sigmaR / ( pwE + vL + sigma R ) = 110/1120;
        % A(3,:) = [0 0 0 sigmalocal] - [0 plocal*wlocal vlocal sigmalocal]*NgRelapseFraction;
        % b(3) = 0;
    
        % X+E+L+R = pi_c
        A(4,:) = [1 1 1 1];
        b(4) = piW_t(t,w);
    
        XELR_local = A\b;
    
        XELRN_w(w,:) = XELR_local';
    end

    for k=1:numq
        qvec_local(k) = sum(XELRN_w(:,k))/sum(piW_t(t,:));
    end
    
    % qvec_local(2) = sum(XELRN_w(:,2))/sum(pi_t(t,:));
    % qvec_local(3) = sum(XELRN_w(:,3))/sum(pi_t(t,:));
    % qvec_local(4) = sum(XELRN_w(:,4))/sum(pi_t(t,:));

    qvec_t(t,:) = qvec_local;
end