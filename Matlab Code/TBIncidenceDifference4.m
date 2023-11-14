%% Parameters Setup

%Parameters from York Paper
beta = 1*10^-8;        % Transmission rate within foreign-born population in Canada. 
pi = 223840;           % Average number of annual new immigrants to Canada     
w = 0.4;               % (q1+q2) Proportion of all LTBI immigrants pass through latent stage in the first 2.5 years
p = 0.05;              % P(progresses directly to active TB stage from early latent without treatment)
v = 0.0002;            % The rate of slow progression to active TB due to reactivation 
a = 0.06;              % TB-caused death rate 
d = 0.8;               % Constant rate of recovery by nature or treatment 
%q1 = 0.03;             % Percentages of early latent (high risk) new immigrants to develop TB
q1 = 0.05;
q2 = 0.37;             % Percentages of late latent latent (low risk) new immigrants to develop TB q
tspan = [0:1:100];     % Time span

n = 0.039;             %Natural removal rate, dX = dE = dL = dT = dR = 0.039

%Actual Pi, annual new immigration to Canada every year 2010-2020
ACTUAL_PI = [259110 260036 263101 267924 240763 323192 272707 303325 313601 284157 226309];



%Actual Pi, annual new immigration to Canada every year 2001-2020
%Actual_Pi = [223840 199170 239083 244578 254374 238125 249622 245289 270581 259110 260036 263101 267924 240763 323192 272707 303325 313601 284157 226309];

ACTUAL_TB_RATE = [14.1 14.7 14.6 15 14.3 15 15.5 15 14.8 15.9 14.3]; %Actual TB rate from 2010 - 2020
INIT_ACTUAL_TB_RATE = ACTUAL_TB_RATE(1);


% 2010 Initial Conditions from TB Surveilance by Stat Can
INIT_X0 = 5874633;
INIT_E0 = 20000; %22533 ; %14969; %12969;
%INIT_L0 =  3016687;
INIT_T0 = 1054;
INIT_R0 = 0;
y0 = [INIT_X0, INIT_E0, -1, INIT_T0, INIT_R0];

INIT_L0 = (INIT_ACTUAL_TB_RATE*(y0(1) + y0(2) + y0(4) + y0(5)) - 100000*p*w*y0(2))/(100000*v-INIT_ACTUAL_TB_RATE);

y0 = [INIT_X0, INIT_E0, INIT_L0, INIT_T0, INIT_R0];

INIT_ACTUAL_TB_RATE = ACTUAL_TB_RATE(1);




L0_LIMIT = 16000000;

%% Initialize time grid


% assume time is years after 2010.  create vector to match ACTUAL_TB_RATE
sizet = length(ACTUAL_TB_RATE);
sizet = sizet-1;

TimeGrid = 1:sizet;



%% Calculate TB Incididence 

% Initialize the final solution Y = (X, E, L, T, R).  it will be 5 by
% size_t. each row represents a time




    Tr = [getTBIncidenceRate(y0, p, w, v)];

    y0local = y0;
    for pi = ACTUAL_PI

        % y = (x,e,l,t,r) is solution to system of odes.
        [t,y] = ode23(@(t,y) odefcn(t, y, beta, pi, p, w, v, a ,d, q1, q2, n), [0,1], y0local); % [0,1] means 1 year forward

        y0local = y(end,:);
        Tri = getTBIncidenceRate(y0local, p, w, v);
        Tr = [Tr Tri(1,1)];


        %y0local 
    end
%% Define time grid
   
%LineSpec_Colors = ["-r" "-g" "-b" "-c" "-m" "-y"];
%YEAR_START = 2010;
%YEAR_END = 2020;

%timegrid = YEAR_START:YEAR_END;

TimeGrid0 = [0, TimeGrid];

    plot(TimeGrid0,ACTUAL_TB_RATE(1:end),'r--','LineWidth',1); 
    % m = "Actual TB rate";
    hold on

    plot(TimeGrid0,Tr(1:end-1), 'k')
    
    % plot immigration rates
    PI_mean = ACTUAL_TB_RATE(2);
    PI_shifted = (ACTUAL_PI-mean(ACTUAL_PI))/std(ACTUAL_PI) + PI_mean;
   plot(TimeGrid0,PI_shifted(1:end), 'g*-' );
    
hold off

title('Canadian Foreign-born Population TB Incidence Rate')
legend('Reported TB Incidence','Estimated TB Incidence', 'Scaled Immigration Rate','Location','northwest' )
xlabel('Years After 2010')
ylabel('New Active TB Cases (per 100,000 population)')

%ylim([0 20])

saveas(gcf, 'TB Incidence.png')
%% 
        %[t,y] = ode23(@(t,y) odefcn(t, y, beta, pi, p, w, v, a ,d, q1, q2, n), tspan, y0);
        %y0 = y(2,:);
        %Tri = getTBIncidenceRate(y, p, w, v);
        %Tr = [Tr Tri(1,1)];

    %% 

%% 

% lgd.FontSize = 14;

%% 

% legend_a = plot(TimeGrid,min_diff_tb_rate,'-r','LineWidth',1); m = "E_0 = " + min_diff_init_E_0 + ", L_0 = " + min_diff_init_L_0;
% A = [A, legend_a]; M = [M, m];
% 
% legend_a = plot(TimeGrid,ACTUAL_TB_RATE,'-k','LineWidth',1); m = "Actual TB rate";
% A = [A, legend_a]; M = [M, m];
% lgd = legend(A, M); %Labelling the lines
% lgd.FontSize = 14;
% 
% hold off



%% Functions

%% Setup the ODE System
%P.4 of York Paper

% X <-> y(1)
% E <-> y(2)
% L <-> y(3)
% T <-> y(4)
% R <-> y(5)


function dydt = odefcn(t, y, beta, pi, p, w, v, a ,d, q1, q2, n)

dydt = zeros(4,1);
dydt(1) = (1-q1-q2)*pi - beta*y(1)*y(4) - n*y(1);
dydt(2) = q1*pi + beta*y(1)*y(4)-(n+w)*y(2);
dydt(3) = q2*pi + (1-p)*w*y(2) - (n+v)*y(3);
dydt(4) = p*w*y(2) + v*y(3) - (n+a+d)*y(4);
dydt(5) = d*y(4) - n*y(5);
end

%% calculate the TB Incidence Rate
function Tri = getTBIncidenceRate(y, p, w, v)
    Xi = y(:,1);
    Ei = y(:,2);
    Li = y(:,3);
    Ti = y(:,4);
    Ri = y(:,5);

    FPi = Xi + Ei + Li + Ti + Ri; %Foreign Population
    
    Tri = p*w*Ei + v*Li; %Page 10(704) of the York paper
    Tri = Tri * 100000./FPi;
end

%% To-dos
% Make pi into a smooth function
% vary q1 and q2
% E0 and L0
