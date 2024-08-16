function err = IncidenceError6_rescale(y, allParams, initalScales, ImmigrationRate, ReportedIncidence)
% input: x is in qE qL qR sigma X0 E0 L0 T0 R0
% bioParameters as described, R22
% paramMagnitude:
    % (1-9): x9. reflects magnitude of estimated param
%% unpack parameters

x = y.*initalScales;

localparams = allParams;
localparams(9)=x(1); % q1
localparams(10)=x(2); % q2
localparams(11) = x(3); %q3
localparams(8) = x(4); % u

localparams(12:16) = x(5:end);

% localIC(1) = y(5); %X0
% localIC(2) = y(6); %E0
% localIC(3) = y(7); %L0
% localIC(4) = y(8); %T0
% localIC(5) = y(9); %R0




err = IncidenceError6(x, localparams, ImmigrationRate, ReportedIncidence);




end
