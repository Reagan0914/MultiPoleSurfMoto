function Y=polybase1(r,N)
L=1;
for i=0:1:N;
    for k=0:i;
        for l=0:k;
            XX(1,L)=r(1,1)^(i-k)*r(1,2)^(k-l)*r(1,3)^l;
            L=L+1;
         end
    end
end
Y=XX;
            