%Function that calculates an array of boolean complementairy codes

%Required inputs:

%N,M: N x M coded aperture size
%Q: Number of snapshots

%Output:
%code: Array of boolean patterns



function code=my_bool(N,M,Q)
for i=1:N
    for j=1:M
        b=1:Q;
        k=randi([1,Q]);
        b2=b(b~=k);
        code2(i,j,k)=1;
        code2(i,j,b2)=0;
    end
end

for i=1:Q
    code{i}=code2(:,:,i);
end