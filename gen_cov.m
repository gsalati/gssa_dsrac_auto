function y = gen_cov(k,nh)
    %%
%     y = zeros(2*nh+1,2);
    j = 1;
    for i = -nh:nh
        foo = abs(k-i);
        if foo <= nh
            y(j,1) = k-i;
            y(j,2) = i;
            j = j+1;
       end
   end
end

