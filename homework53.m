%% Higher-order Taylor Methods, 2d
b = 2;
a = 1;
h = 0.25;
N = (b-a) / h;
y(1) = 2;
t = a:h:b;
n = 2; % order two taylor method
for i=1:N
   y(i+1) = y(i) + h*(1/t(i)^2)*(sin(2*t(i))-2*t(i)*y(i)) + (.5*h^2)*((2*t(i)*y(i)+2*t(i)*cos(2*t(i))-2*sin(2*t(i)))/(t(i)^3) -2*(1/t(i)^2)*(sin(2*t(i))-2*t(i)*y(i)) / t(i));
end
figure
plot(t, y, 'DisplayName', 'order 2 taylor');
title('Higher Order Taylor Approximation H=0.25');
xlabel('1 < t < 2');
ylabel('y');
hold on
%% Exact Solution, Just for checking;

for i=1:N+1
    yexact(i) = (4 + cos(2) - cos(2*t(i)))/(2*t(i).^2);
end
plot(t, yexact, 'DisplayName', 'exact solution');
%% 4d, same problem but now T^4
y4(1) = 2; 
for i=1:N
     y4(i+1) = y4(i) + h*(1/t(i)^2)*(sin(2*t(i))-2*t(i)*y4(i))...
     + (.5*h^2)*((2*t(i)*y4(i)+2*t(i)*cos(2*t(i))-2*sin(2*t(i)))/(t(i)^3)...
     -2*(1/t(i)^2)*(sin(2*t(i))-2*t(i)*y4(i)) / t(i))...
     + (h^3/6)*((2*t(i)*sin(2*t(i))-24*t(i)*y4(i)+12*sin(2*t(i))-12*t(i)*cos(2*t(i)))/t(i)^4)...
     + (h^4/24)*(120*t(i)*y4(i)+24*(t(i)^2)*sin(2*t(i))+60*t(i)*cos(2*t(i))-6*t(i)*sin(2*t(i))-72*sin(2*t(i)))/t(i)^5
end
plot(t, y4, 'DisplayName','order 4 taylor');
legend
hold off
%% 10 A,B,C
b = 2;
a = 1;
t = a:h:b;
h = .05;
N = (b-a)/h; 
t = a:h:b;
y1(1) = -1;
for i = 1:N
    y1(i+1) = y1(i) + h * ((1/t(i)^2)-y1(i)/t(i) - y1(i)^2)...
    +((h^2)/2)*((y1(i)/t(i)^2)-2/(t(i)^3) +(-1/t(i)-2*y1(i))*((1/t(i)^2)-y1(i)/t(i) - y1(i)^2))...
    +((h^3)/6)*((9/t(i)^4)-((3*y1(i)^2)/t(i)^2) + ((1/t(i)^2)-(y1(i)/t(i)) - y1(i)^2)*(6*y1(i)^2 + 6*y1(i)/t(i)))...
    +((h^4)/24)*(((-36/t(i)^5)+((6*y1(i)^2)/(t(i)^3)) + ((12*(t(i)^3)*(y1(i)^3)-18*y1(i))/t(i)^4)...
    + ((1/t(i)^2)-(y1(i)/t(i))-y1(i)^2)*((6/t(i)^3) - (6*y1(i)/t(i)^2) - ((36*y1(i)^2)/t(i)) - 24*y1(i)^3)))
end
figure
plot(t, y1, 'DisplayName','order 4');
hold on
y2(1) = -1; 
for i = 1:N
    y2(i+1) = y2(i) + h * ((1/t(i)^2)-y2(i)/t(i) - y2(i)^2)...
    +((h^2)/2)*((y2(i)/t(i)^2)-2/(t(i)^3) +(-1/t(i)-2*y2(i))*((1/t(i)^2)-y2(i)/t(i) - y2(i)^2));
   
end
plot(t, y2, 'DisplayName', 'order 2');
title('Higher Order Taylor Approximation 4d');
xlabel('1 < t < 2');
ylabel('y');
legend

order2 = 1:3;
order4 = 1:3;
error = 1:3;

order2(1) = interp1(t, y2, 1.052);
order2(2) = interp1(t, y2, 1.555);
order2(3) = interp1(t, y2, 1.978);

order4(1) = interp1(t, y1, 1.052);
order4(2) = interp1(t, y1, 1.555);
order4(3) = interp1(t, y1, 1.978);

error(1) = order2(1) - order4(1);
error(2) = order2(2) - order4(2);
error(3) = order2(3) - order4(3);

myTable = table(order2, order4, error)











