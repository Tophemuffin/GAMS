$onText
Work with quick start
$offText

set tt  /0*9/;

Parameter
tax;
tax=0.25;

Parameter  taxx(tt)  tax range;

taxx(tt)= tax*(ord(tt)-1);

Display
taxx;


Positive Variables  P, Qd , Qs;

Equations     Pdemand, Psupply, Equilibrium;
Pdemand..     6 - 0.3*Qd            =g= P;
Psupply..     P =g= ( 1 + 0.2*Qs);
Equilibrium.. Qs           =g= Qd;

Model problem / Pdemand.Qd, Psupply.Qs, Equilibrium.P /;
*starting guesses
P.l=1;
Qd.l=1; 
Qs.l=1;

solve problem using MCP;

*post solution calculations
Parameters
CS_base
PS_base
;

CS_base = (1/2)*(Qd.l)*(6-P.l);
PS_base = (1/2)*(Qd.l)*(P.l-1);

display
CS_base
PS_base
;

*Policy Model
*policy Variable


Positive Variables  P_p, Qd_p, Qs_p;

Equations     Pdemand_p, Psupply_p, Equilibrium_p;
Pdemand_p(tt)..     6 - 0.3*Qd_p(tt)            =g= P_p(tt);
Psupply_p(tt)..     P_p(tt)-taxx(tt) =g= ( 1 + 0.2*Qs_p(tt));
Equilibrium_p(tt).. Qs_p(tt)           =g= Qd_p(tt);

Model problem_p / Pdemand_p.Qd_p, Psupply_p.Qs_p,
Equilibrium_p.P_p /;
*starting guesses
P_p.l(tt)=1;
Qd_p.l(tt)=1; 
Qs_p.l(tt)=1;

solve problem_p using MCP;

*post solution calculations
Parameters
CS_p
PS_p
taxrev
dwl
;

CS_p(tt) = (1/2)*(Qd_p.l(tt))*(6-P_p.l(tt));
PS_p(tt) = (1/2)*(Qd_p.l(tt))*(P_p.l(tt)-taxx(tt)-1);
taxrev(tt) = taxx(tt)*Qd_p.l(tt);
dwl(tt) = CS_base + PS_base - (CS_p(tt)+PS_p(tt)+taxrev(tt));


display
CS_base
PS_base
CS_p
PS_p
taxrev
dwl
;










File TH1EXT / TH1EXT.txt /;
put TH1EXT;
TH1EXT.nd=6;
put "Results, step, tax, CS, PS, taxrev, DWL" /;
loop((tt),
    put ord(tt), taxx(tt), CS_p(tt), PS_p(tt), taxrev(tt), dwl(tt) /
);
putclose;