function [OUT, bi] = TOWHALF(IN, REF, N, Y)
if  IN<-5/8*REF
    A=0;
    OUT=(4+Y)*IN+(3+Y)*REF;
elseif -5/8*REF<=IN && IN<-3/8*REF
    A=1;
    OUT=(4+Y)*IN+(2+Y);
elseif -3/8*REF<=IN && IN<-1/8*REF
    A=2;
    OUT=(4+Y)*IN+(1+Y)*REF;
elseif -1/8*REF<=IN && IN<1/8*REF
    A=3;
    OUT=(4+Y)*IN;
elseif 1/8*REF<=IN && IN<3/8*REF
    A=4;
    OUT=(4+Y)*IN-(1+Y)*REF;
elseif 3/8*REF<=IN && IN<5/8*REF
    A=5;
    OUT=(4+Y)*IN-(2+Y)*REF;
elseif 5/8*REF<=IN
    A=6;
    OUT=(4+Y)*IN-(3+Y)*REF;
end
bi=A*2^(2*N-1);
end