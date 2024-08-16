function TI = WHOincidence(REGION)
    TI = -1;
    if strcmp(REGION,'AFR')
        TI=281;
    elseif strcmp(REGION,'AMR')
        TI=28;
    elseif strcmp(REGION,'EMR')
        TI=117;
    elseif strcmp(REGION,'EUR')
        TI=37;
    elseif strcmp(REGION,'SEAR')
        TI=211;
    elseif strcmp(REGION,'WPR')
        TI=85;
    end
end