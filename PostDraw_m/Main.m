clear
clc
load('25143itokawa_obj.mat')
[azimuth,elevation,rr] = cart2sph(vertices(:,1),vertices(:,2),vertices(:,3));
%% 
Len = length(vertices(:,1));
for i=1:Len; 
x(i,1)=cos(elevation(i,1))*cos(azimuth(i,1));
x(i,2)=cos(elevation(i,1))*sin(azimuth(i,1));
x(i,3)=sin(elevation(i,1));
i=i+1;
end
axis equal
%%
for N = 1:1000
    if (N*(N^2+6*N+11)/6+1)>=Len
        break
    end
end
disp('N = '); 
disp(N);
prompt='Input order £¨Not larger than N£©£º'; 
N=input(prompt);
MX = zeros(Len,N*(N^2+6*N+11)/6+1);
A = zeros(N*(N^2+6*N+11)/6+1,1);
for k = 1:1:Len
    r = x(k,:);
    Y = polybase1(r,N);
    MX(k,:) = Y;
end
C=cond(MX,2);
fprintf('The condition number of the system is£º%f\n',C);
%% 
g = lamda(MX);
j = 1;
    A(:,j) = lsmr(MX,rr(:,j),g, 10^-8, 10^-8, 1.0000e+12,10^10, 0, 0);
    
dlmwrite('AP.TXT', A, 'delimiter', '\t','precision', 6,'newline', 'pc');

%% 
NN = 100;
[x1,y1,z1]=sphere(NN);
X = x1*0; Y = y1 * 0; Z = z1 * 0; 
for i = 1:1:NN+1
    for j = 1:1:NN+1
        r = [x1(i,j),y1(i,j),z1(i,j)];
        vec = polybase1(r,N);      
        [fai(i,j),theta(i,j),rho(i,j)] = cart2sph(x1(i,j),y1(i,j),z1(i,j));
        X(i,j) = vec*A(:,1)*cos(theta(i,j))*cos(fai(i,j));       
        Y(i,j) = vec*A(:,1)*cos(theta(i,j))*sin(fai(i,j));
        Z(i,j) = vec*A(:,1)*sin(theta(i,j));
    end
end
%%
figure; 
set(gcf, 'Position', [10 10 1150 650]);
hold on
axis equal
view(3);
light('Position',[-10,-10,10]);
lighting gouraud;
surf(X,Y,Z)
colormap gray
axis equal


