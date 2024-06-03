function plotPopulationVsTime(t,Y )
%{
  Input: 
    - bioParameters vector in R7, hard coded so position [beta p w v a d n q1 q2]
    - params2 vector in R5, hard coded so position [X0 E0 L0 T0 R0]
    - ImmigrationRate, which gives annual immigration rate, vector in R^N

  Output:
    - XELTR, an N+1 (N is size of ImmigrationRate) by 5 matrix where each row is a time
    - TBIncidence, vector in R^N+1 with calculated incidence
%}



% Graph settings
subplot(2, 3, 1);
plot(t,Y(:,1)); grid on;
title('X vs t');
xlabel('time'); ylabel('susceptible');

subplot(2, 3, 2);
plot(t,Y(:,2)); grid on;
title('E vs t');
xlabel('time'); ylabel('early latent');

subplot(2, 3, 3);
plot(t,Y(:,3)); grid on;
title('L vs t');
xlabel('time'); ylabel('late latent');

subplot(2, 3, 4);
plot(t,Y(:,4)); grid on;
title('T vs t');
xlabel('time'); ylabel('actively infectious');

subplot(2, 3, 5);
plot(t,Y(:,5)); grid on;
title('R vs t');
xlabel('time'); ylabel('recovered');

subplot(2,3,6)
% plot(t(1:end-1),ReportedImmigration)
% title('\pi vs t')
% xlabel('time'); ylabel('immigration rate');
% grid on;

end