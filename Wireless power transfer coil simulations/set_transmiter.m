function  state = set_transmiter(s,transmiter_num,transmiter_state,state)
% s = serialport(serial_com,115200);
TXbuf(1) = double('T');
TXbuf(2) = transmiter_num + transmiter_state * 2 ^ 4;
write(s , TXbuf , "uint8");
 state(transmiter_num) = transmiter_state;
% clear s
end
