$onText
Work with quick start
$offText


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

