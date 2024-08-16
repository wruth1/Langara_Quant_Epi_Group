function err = error_E2021(XELTR,  JordanEoverTB_mean)

    % JordanEoverTB_mean = 1/488;
    
    % power_relapse = allParams(18);
    

    

    E2021TB_ratio = XELTR(end,2)/sum(XELTR(end,[2:5]));

    % multiplying by 100 changes from decimal to percentages. 
    % drives us to be within 0.9*NgRelapseFraction to 1.1*NgRelapseFraction
    err = ((E2021TB_ratio-JordanEoverTB_mean)/JordanEoverTB_mean)^2;
    % scalar

   
   

end