function burden_t = compute_LTBI(XELTR)

    numt=size(XELTR,1);

    burden_t=zeros(numt,1);

    for j =1:numt
        burden_t(j) = sum(XELTR(j,[2,3,5]))./sum(XELTR(j,:));
    end
   
   

end