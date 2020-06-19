% 
clc
clear
Main
hold on
u = 1;  
C = cell(u, 1);
 axis([-3 3 -3 3 -2 2]);
tu = 0.002;
for i = 1:u
    C{i, 1} = load(['GLPATH',int2str(i),'.BT']);
end
%% calculate the position of equal time  by interpolation 
D = cell(u, 1);
E = cell(u, 1);
for j = 1:u
    time = C{j, 1}(:, 1);
    pos = [C{j, 1}(:, 2),C{j, 1}(:, 3),C{j, 1}(:, 4)];
    D{j, 1}(:, 1) = [time(1):0.002:time(end)];
    pos1 = interp1(time,pos,D{j, 1}(:, 1));
    D{j, 1}(:, 2) = pos1(:, 1);
    D{j, 1}(:, 3) = pos1(:, 2);
    D{j, 1}(:, 4) = pos1(:, 3);
    E{j, 1}(:, 1) = D{j, 1}(:, 2);
    E{j, 1}(:, 2) = D{j, 1}(:, 3);
    E{j, 1}(:, 3) = D{j, 1}(:, 4);
end   
F = cell(u, 1);
for i = 1:1:u
    ts = 2.2;
    time = D{i, 1}(:, 1);
    pos = [D{i, 1}(:, 2),D{i, 1}(:, 3),D{i, 1}(:, 4)];
    posx = pos(:, 1);
    posy = pos(:, 2);
    posz = pos(:, 3);
    posxe = posx(end);
    posye = posy(end);
    posze = posz(end);
    pose = [posx(end),posy(end),posz(end)];
    tn = time(end);
    time1= [];
    pos1 = [];
    posc = [];
    if tn < ts
        tc = [tn:0.002:ts];
        tct = tc';
        time1 = [time; tct];
        A = size(tc);
        N = A(2);
        posc(:, 1) = linspace(posxe,posxe,N);
        posc(:, 2) = linspace(posye,posye,N);
        posc(:, 3) = linspace(posze,posze,N);
        posct = posc';
        pos1 = [pos; posc];
        F{i, 1}(:, 1) = time1;
        F{i, 1}(:, 1) = pos1(:, 1);
        F{i, 1}(:, 2) = pos1(:, 2);
        F{i, 1}(:, 3) = pos1(:, 3);
    else
    end
end
%%
x = [];
y = [];
z = [];
for i = 1:1:1102
    for n = 1:1:u
    x(n, i) = F{n, 1}(i, 1);
    y(n, i) = F{n, 1}(i, 2);
    z(n, i) = F{n, 1}(i, 3);
    end
end
for m = 1:1:1102
        Hp = plot3(x(:, m), y(:, m), z(:, m),'r.','Marker','o','MarkerFaceColor','r','MarkerSize',3);
        rotate(Hp, [0 0 1], 0.72*m);
        rotate(p, [0 0 1], 0.72);
        shading interp
    drawnow;
if m < 1102;
   delete(Hp);
else
end

end
        

