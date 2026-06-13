$PROBLEM 1cmt iv example model

$INPUT ID TIME AMT EVID DV CMT ET1 ET2 ITER OID

$DATA ../data/sim_nm.csv IGNORE=@

$SUBR ADVAN1 TRANS2

$PK
CL = THETA(1) * EXP(ET1)
V = THETA(2) * EXP(ET2)
S1 = V

$ERROR
IPRED = F
Y = IPRED * (1 + THETA(3) * EPS(1)) 

$THETA
8.66 ; CL
100 ; V
0.2 FIX

$OMEGA BLOCK(2) FIX ; just a dummy for now
0.25 
0.00050625 0.09

$SIGMA 
1 FIX

; $EST METHOD=1 MAXEVAL=0 POSTHOC
$SIM (12345) (54321) ONLYSIM SUBPROBLEMS=1 ; use newly created simulation dataset with prop.score matching, sim only once!

$TABLE ID TIME EVID DV ITER OID
  NOAPPEND NOHEADER NOPRINT
  FILE=simtab1_pm
