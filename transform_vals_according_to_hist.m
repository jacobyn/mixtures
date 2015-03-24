function v2n=transform_vals_according_to_hist(v1,v2)
%v1=[randn(100,1)*2;10+randn(100,1)*1];
%v2=[randn(200,1)];
assert(length(v1)==length(v2));
[~,idx1]=sort(v1);
[val2,~]=sort(v2);
v2n=val2(idx1);

%figure(2);clf;nhist({v1,v2,v3});