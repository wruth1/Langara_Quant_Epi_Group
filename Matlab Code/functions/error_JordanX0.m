function err = error_JordanX0(XELTR,  JordanX0Ratio)

    % JordanEoverTB_mean = 1/488;
    
    % power_relapse = allParams(18);
    

    

    % E2021TB_ratio = XELTR(end,2)/sum(XELTR(end,[2:5]));

    % multiplying by 100 changes from decimal to percentages. 
    % drives us to be within 0.9*NgRelapseFraction to 1.1*NgRelapseFraction
    err = (( XELTR(1,1)/sum(XELTR(1,1:5)) - JordanX0Ratio )/JordanX0Ratio ).^2;
    % scalar

   
   

end