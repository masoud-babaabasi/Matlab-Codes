function [OUT] = FLASH2BIT (IN, REF)
if  IN<-0.5*REF
    OUT=0;
end
if -0.5*REF<=IN && IN<0 
    OUT=1;
end
if 0<=IN && IN<0.5*REF
    OUT=2;
end
if 0.5*REF<=IN
    OUT=3;
end
end