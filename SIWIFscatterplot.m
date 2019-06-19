%[X,Y,Z] = meshgrid(0:.1:1);%defines the domain and the amount of dots
%[x]=meshgrid(0:.1:1,0:.1:1,0:.1:1)
%[y]=meshgrid(0:0.0791:0.791,0:0.0791:0.791,0:0.0791:0.791)
%[z]=meshgrid(-22:2:-2,-22:2:-2,-22:2:-2)
X=0:.1:1;%parameters 
Y=0:0.0791:0.791;
Z=-22:2:-2;
[x,y,z] = meshgrid(X,Y,Z);
[q,s,e]=size([x,y,z]); 
%a0=2*1.315*rand(q,1);
%this wiggles the parameters
a1=2*0.573*rand(q,1);
a2=2*3.281*rand(q,1);
%a3=2*10.614*rand(q,1);
a4=2*1.034*rand(q,1);
a5=2*1.297*rand(q,1);
a6=2*1.26*rand(q,1);
a8=2*rand(q,1);
a9=2*rand(q,1);
a10=2*rand(q,1);
syms n i d%integral parameter
i=1.315*rand();
P = i*n^0.573*(1-n)^3.281*(1-d*n^1.034);
b=vpa(solve(int(P,n,[0 1]) == 1));
h=double(b);
a0=i.*(1:q)';
a11=h.*(1:q)';
w=(((x./y-a0).^a9+a10)+a0.*x.^a1.*(1-x).^a2.*(1+a11.*x.^a4)).*a5.*(1-z./(a6.^2)).^(-2) ; 
%above is the function that we wish to graph
%mw=mean(w)
%mmw= mean(mw)
%mmmw=mean(mmw)
 scatter3(x(:),y(:),z(:),[],w(:)); %creates the scatter plot
xlabel('x')%labels for each axis
ylabel('\xi')
zlabel('t')   
cb = colorbar; %creates the colorbar
cb.Label.String = 'h (units)';%label for the color bar
%We are ignoring possible complex values
%Other interesting things to note is that there is always a slice that has 
%different values that occurs at y=0
%in addition there is always small pockets of values that have a lower
%h value than the average