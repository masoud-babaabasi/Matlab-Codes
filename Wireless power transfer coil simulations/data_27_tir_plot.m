clear
filename = 'data_27_tir.xlsx';
%% receiver center of transmiter power VS tansmiter voltage
start = 2;
Vs = 2:0.5:10;
H = 12:12:48;
Pr = zeros(length(H),length(Vs));
Effi = zeros(length(H),length(Vs));
n_i = length(Vs);
for i=1:length(H)
    rangex = 'F'+string(start+ ( i - 1 ) * n_i )+':F'+string(start + i * n_i );
    Pr(i,:) = table2array(readtable(filename,'Range',rangex));
    rangex = 'I'+string(start+ ( i - 1 ) * n_i )+':I'+string(start + i * n_i );
    Effi(i,:) = table2array(readtable(filename,'Range',rangex));
end
%% plot 3d
figure('Name','single transmiter recieved power')
surf(Vs,H,Pr)
xlabel('Transmitter Vltage(V)')
ylabel('Receiver Height(mm)')
zlabel('Received Power(w)')
figure('Name','single transmiter efficiendy')
surf(Vs,H,Effi)
xlabel('Transmitter Vltage(V)')
ylabel('Receiver Height(mm)')
zlabel('Efficiency(%)')
%% plot 2d
figure('Name','single transmiter recieved power')
v = 10;
plot( H , Pr(:,Vs == v) ) 
xlabel('Receiver Height(mm)')
ylabel('Received Power(W)')
xticks(H);
yticks(0:0.5:5);
xlim([H(1) H(end)])
legend("Vs = 10V")
grid on
%%
h = 36;
figure('Name','single transmiter recieved power')
plot( Vs , Pr(H == h,:)*1000) 
xlabel('Transmitter Vltage(V)')
ylabel('Received Power(mW)')
xticks(Vs);
xlim([Vs(1) Vs(end)])
legend("H = 36mm")
grid on
%%
figure('Name','single transmiter efficiendy')
plot( H , Effi(:,Vs == v))
xlabel('Receiver Height(mm)')
ylabel('Efficiency(%)')
xticks(H);
xlim([H(1) H(end)])
grid on