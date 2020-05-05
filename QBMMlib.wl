(* ::Package:: *)

BeginPackage["QBMMlib`"];

pow::usage = "";
add::usage = "";
wheeler::usage = "";
project::usage = "";
project1::usage = "";
project2::usage = "";
momidx::usage = "";
cqmom12::usage = "";
cqmom21::usage = "";
cqmom12m::usage = "";
cqmom21m::usage = "";
hyqmom::usage = "";
chyqmom::usage = "";
chyqmom4::usage = "";
chyqmom9::usage = "";
pointer::usage="";
quad::usage="";

Begin["`Private`"];

pow[x_, 0] := 1 /; x == 0 || x == 0.
pow[x_, y_: 0] := x^y

add[x_, y_] := x + y;

quad[w_,xi_,xis_,method_,nr_,nrd_,perm_:0]:=Module[{momq},
	If[method=="CQMOM",
		If[perm==12,momq[p_,q_]:=Sum[w[j,i]xi[[j]]^p xis[j][[i]]^q,{j,nr},{i,nrd}]];
		If[perm==21,momq[p_,q_]:=Sum[w[j,i]xi[[j]]^q xis[j][[i]]^p,{j,nrd},{i,nr}]];
	];
	If[method=="CHyQMOM",
		momq[p_,q_]:=Sum[w[[i]]xi[[i]]^p xis[[i]]^q,{i,Length[w]}];
	];
	Return[momq,Module];
];

wheeler[m_, n_] := 
  Module[{nn = n, mm = m, \[Sigma], a, b, Ja, w, xi, eval, evec, 
    esys}, 
    \[Sigma] = Table[0., {i, 2 nn}, {j, 2 nn}]; 
    Do[\[Sigma][[2, i + 1]] = mm[[i + 1]], {i, 0, 2 nn - 1}]; 
    a = Table[0., {i, nn}]; b = a; a[[1]] = mm[[2]]/mm[[1]]; 
    b[[1]] = 0; 
    Do[Do[\[Sigma][[i + 2, j + 1]] = \[Sigma][[i + 1, j + 2]] - 
        a[[i]] \[Sigma][[i + 1, j + 1]] - b[[i]] \[Sigma][[i, j + 1]];
       a[[i + 1]] = -(\[Sigma][[i + 1, i + 1]]/\[Sigma][[i + 1, 
          i]]) + \[Sigma][[i + 2, i + 2]]/\[Sigma][[i + 2, i + 1]]; 
      b[[i + 1]] = \[Sigma][[i + 2, i + 1]]/\[Sigma][[i + 1, i]];, {j,
        i, 2 nn - i - 1}];, {i, nn - 1}];
   Ja = DiagonalMatrix[a]; 
   Do[
     Ja[[i, i + 1]] = -Sqrt[Abs[b[[i + 1]]]]; 
     Ja[[i + 1, i]] = -Sqrt[Abs[b[[i + 1]]]];
   ,{i,nn-1}]; 
   esys = Eigensystem[Ja]; eval = esys[[1]]; evec = esys[[2]]; 
   w = Table[evec[[i, 1]]^2 mm[[1]],{i,nn}]; 
   Return[{eval, w}, Module];
];
   
project[xix_, xiy_, w_, ks_, ksp_, nr_, nrd_] := 
  Module[{moms, mom, momsp}, 
   mom[p_, q_] := Sum[w[[i]] xix[[i]]^p xiy[[i]]^q, {i, nr nrd}]; 
   moms = Table[mom[ks[[i, 1]], ks[[i, 2]]], {i, 1, Length[ks]}]; 
   momsp = Table[mom[ksp[[i, 1]], ksp[[i, 2]]], {i, 1, Length[ksp]}]; 
  {moms,momsp}
];
    
project1[xi_, xis_, wtot_, ks_, ksp_, nr_, nrd_] := 
  Module[{moms, mom, momsp}, 
   mom[p_, q_] := Sum[wtot[j, i] xi[[j]]^p xis[j][[i]]^q, {j,nr}, {i,nrd}]; 
   moms = Table[mom[ks[[i, 1]], ks[[i, 2]]], {i, 1, Length[ks]}]; 
   momsp = Table[mom[ksp[[i, 1]], ksp[[i, 2]]], {i,Length[ksp]}]; 
   {moms, momsp}
];
    
project2[xi_, xis_, wtot_, ks_, ksp_, nr_, nrd_] := 
  Module[{moms, mom, momsp}, 
   mom[p_, q_] := Sum[wtot[j, i] xi[[j]]^q xis[j][[i]]^p, {j,nrd}, {i,nr}]; 
   moms = Table[mom[ks[[i, 1]], ks[[i, 2]]], {i, 1, Length[ks]}]; 
   momsp = Table[mom[ksp[[i, 1]], ksp[[i, 2]]], {i, 1, Length[ksp]}]; 
   {moms, momsp}
];
    
momidx[nr_, nrd_, method_, nro_: 0] := 
  Module[{ks, k1, k2}, 
   If[method == "CHyQMOM", 
    If[nr == 2, 
     If[nro == 0, 
       ks = {{0, 0}, {1, 0}, {0, 1}, {2, 0}, {1, 1}, {0, 2}};, 
       ks = Flatten[
          Table[{{0, 0, i}, {1, 0, i}, {0, 1, i}, {2, 0, i}, {1, 1, 
             i}, {0, 2, i}}, {i, 0, nro}], 1];];]; 
    If[nr == 3, 
     ks = {{0, 0}, {1, 0}, {0, 1}, {2, 0}, {1, 1}, {0, 2}, {3, 0}, {0,
         3}, {4, 0}, {0, 4}}]; Return[ks, Module];]; 
   If[method == "CQMOM", 
    If[nro == 0, 
      k1 = Flatten[{Table[{q, p}, {q, 0, nr - 1}, {p, 0, 2 nrd - 1}], 
         Table[{q, p}, {q, nr, 2 nr - 1}, {p, 0, 0}]}, 2]; 
      k2 = Flatten[{Table[{p, q}, {q, 0, nrd - 1}, {p, 0, 2 nr - 1}], 
         Table[{p, q}, {q, nrd, 2 nrd - 1}, {p, 0, 0}]}, 2]; 
      ks = DeleteDuplicates[Join[k1, k2]]; Return[ks, Module];
      ,
      k1 = Flatten[
        Table[Flatten[{Table[{q, p, i}, {q, 0, nr - 1}, {p, 0, 
             2 nrd - 1}], 
           Table[{q, p, i}, {q, nr, 2 nr - 1}, {p, 0, 0}]}, 2], {i, 0,
           nro - 1}], 1]; 
      k2 = Flatten[
        Table[Flatten[{Table[{p, q, i}, {q, 0, nrd - 1}, {p, 0, 
             2 nr - 1}], 
           Table[{p, q, i}, {q, nrd, 2 nrd - 1}, {p, 0, 0}]}, 2], {i, 
          0, nro - 1}], 1]; 
      ks = DeleteDuplicates[Join[k1, k2]]; 
      Return[ks, Module];];]; 
      Print["Error!"];
 ];
      
pointer[moms_,ks_,w_:0,ro_:0] := 
  Module[{eqns, vars, linsolv, pm, mymom}, 
   If[Length[ks[[1]]] == 2, 
    eqns = Table[
      moms[[i]] == pm[ks[[i, 1]], ks[[i, 2]]], {i, Length[ks]}]; 
    vars = DeleteDuplicates[
      Flatten[Table[pm[ks[[i, 1]], ks[[i, 2]]], {i, Length[ks]}]]]; 
    linsolv = First[vars /. Solve[eqns, vars]]; 
    mymom[q_, p_] := linsolv[[First[First[Position[ks, {q, p}]]]]]; 
    Return[mymom, Module];]; 
   If[Length[ks[[1]]] == 3, 
    eqns = Table[
      moms[[i]] == 
       Sum[
        w[[l + 1]] ro[[l + 1]]^
         ks[[i, 3]] pm[ks[[i, 1]], ks[[i, 2]], l], {l, 0, 
         Length[ro] - 1}], {i, Length[ks]}]; 
    vars = DeleteDuplicates[
      Flatten[Table[
        pm[ks[[i, 1]], ks[[i, 2]], l], {l, 0, Length[ro] - 1}, {i, 
         Length[ks]}]]]; linsolv = First[vars /. Solve[eqns, vars]]; 
    mymom[q_, p_, r_] := 
     linsolv[[First[First[Position[ks, {q, p, r}]]]]]; 
    Return[mymom, Module];]; 
    Print["Error, k index not 2D or 3D!"];
    ];
    
chyqmom[momin_, k_, q_, wRo_: 0, Ros_: 0] := 
  Module[{moms = momin, ks = k},
  If[wRo==0,
   If[Mod[Length[moms], 6] == 0, 
    Return[chyqmom4[moms, ks], Module]]; 
   If[Mod[Length[moms], 10] == 0, 
    Return[chyqmom9[moms, ks, q], Module]];
    ,
   If[Mod[Length[moms], 6] == 0, 
    Return[chyqmom4[moms, ks, wRo, Ros], Module]]; 
   If[Mod[Length[moms], 10] == 0, 
    Return[chyqmom9[moms, ks, q, wRo, Ros], Module]];   
   ];
   Print["Error!"];
 ];
   
chyqmom4[momin_, kk_] := Module[{moms = momin, ks = kk, eqns, vars, linsolv, pm, mymom, n, 
    bu, bv, d20, d11, d02, c20, c11, c02, M1, rho, up, Vf, mu2avg, 
    mu2, M3, rh3, up3, vp21, vp22, rho21, rho22, u, v}, 
   mymom = pointer[moms, ks]; 
   n = Table[0, {i, 4}]; u = n; v = n; 
   bu = mymom[1, 0]/mymom[0, 0]; 
   bv = mymom[0, 1]/mymom[0, 0]; 
   d20 = mymom[2, 0]/mymom[0, 0]; 
   d11 = mymom[1, 1]/mymom[0, 0]; 
   d02 = mymom[0, 2]/mymom[0, 0]; 
   c20 = d20 - bu^2; c11 = d11 - bu bv;
   c02 = d02 - bv^2; 
   M1 = {1, 0, c20}; 
   {rho, up} = hyqmom[M1]; 
   Vf = c11 up/c20; 
   mu2avg = c02 - Total[rho Vf^2]; 
   mu2avg = Max[mu2avg, 0]; mu2 = mu2avg; 
   M3 = {1, 0, mu2}; 
   {rh3, up3} = hyqmom[M3]; 
   vp21 = up3[[1]]; 
   vp22 = up3[[2]]; 
   rho21 = rh3[[1]]; 
   rho22 = rh3[[2]]; 
   n[[1]] = rho[[1]] rho21; 
   n[[2]] = rho[[1]] rho22; 
   n[[3]] = rho[[2]] rho21; 
   n[[4]] = rho[[2]] rho22; 
   n = mymom[0, 0] n; 
   u[[1]] = up[[1]]; 
   u[[2]] = up[[1]]; 
   u[[3]] = up[[2]]; 
   u[[4]] = up[[2]]; 
   u = bu + u; 
   v[[1]] = Vf[[1]] + vp21; 
   v[[2]] = Vf[[1]] + vp22; 
   v[[3]] = Vf[[2]] + vp21; 
   v[[4]] = Vf[[2]] + vp22; 
   v = bv + v; 
   Return[{n, u, v}, Module];
   ];
   
chyqmom9[momin_, kk_, qmax_] := 
  Module[{moms = momin, ks = kk, vars, linsolv, mymom, n, csmall, 
    verysmall, bu, bv, d20, d11, d02, d30, d03, d40, d04, c20, c11, 
    c02, c30, c03, c40, c04, M1, rho, up, Vf, M2, rho2, up2, vp21, 
    vp22, vp23, rho21, rho22, rho23, mu2avg, mu2, mu3, mu4, u, v, 
    eqns, pm, q, eta, slope, det, qm, qp, up3, M3, rh3}, 
   mymom = pointer[moms, ks]; n = Table[0, {i, 9}]; u = n; v = n; 
   csmall = 10^(-10); verysmall = 10^(-14); 
   bu = mymom[1, 0]/mymom[0, 0]; bv = mymom[0, 1]/mymom[0, 0]; 
   d20 = mymom[2, 0]/mymom[0, 0]; d11 = mymom[1, 1]/mymom[0, 0]; 
   d02 = mymom[0, 2]/mymom[0, 0]; d30 = mymom[3, 0]/mymom[0, 0]; 
   d03 = mymom[0, 3]/mymom[0, 0]; d40 = mymom[4, 0]/mymom[0, 0]; 
   d04 = mymom[0, 4]/mymom[0, 0]; c20 = d20 - bu^2; c11 = d11 - bu bv;
    c02 = d02 - bv^2; c30 = d30 - 3 bu d20 + 2 bu^3; 
   c03 = d03 - 3 bv d02 + 2 bv^3; 
   c40 = d40 - 4 bu d30 + 6 bu^2 d20 - 3 bu^4; 
   c04 = d04 - 4 bv d03 + 6 bv^2 d02 - 3 bv^4; 
   M1 = {1, 0, c20, c30, c40}; {rho, up} = hyqmom[M1]; 
   If[c20 <= csmall, rho[[1]] = 0; rho[[2]] = 1; rho[[3]] = 0; 
    Vf = 0 up; 
    M2 = {1, 0, c02, c03, c04}; {rho2, up2} = hyqmom[M2, qmax]; 
    vp21 = up2[[1]]; vp22 = up2[[2]]; vp23 = up2[[3]]; 
    rho21 = rho2[[1]]; rho22 = rho2[[2]]; rho23 = rho2[[3]];, 
    Vf = c11 up/c20; mu2avg = c02 - Total[rho Vf^2]; 
    mu2avg = Max[mu2avg, 0]; mu2 = mu2avg; mu3 = 0 mu2; mu4 = mu2^2; 
    If[mu2 > csmall, q = (c03 - Total[rho Vf^3])/mu2^(3/2); 
     eta = (c04 - Total[rho Vf^4] - 6 Total[rho Vf^2] mu2)/mu2^2; 
     If[eta < q^2 + 1, 
      If[Abs[q] > verysmall, slope = (eta - 3)/q; det = 8 + slope^2; 
       qp = 0.5 (slope + Sqrt[det]); qm = 0.5 (slope - Sqrt[det]); 
       If[Sign[q] == 1, q = qp, q = qm];, q = 0]; eta = q^2 + 1;]; 
     mu3 = q mu2^(3/2); mu4 = eta mu2^2;]; 
    M3 = {1, 0, mu2, mu3, mu4}; {rh3, up3} = hyqmom[M3, qmax]; 
    vp21 = up3[[1]]; vp22 = up3[[2]]; vp23 = up3[[3]]; 
    rho21 = rh3[[1]]; rho22 = rh3[[2]]; rho23 = rh3[[3]];]; 
   n[[1]] = rho[[1]] rho21; n[[2]] = rho[[1]] rho22; 
   n[[3]] = rho[[1]] rho23; n[[4]] = rho[[2]] rho21; 
   n[[5]] = rho[[2]] rho22; n[[6]] = rho[[2]] rho23; 
   n[[7]] = rho[[3]] rho21; n[[8]] = rho[[3]] rho22; 
   n[[9]] = rho[[3]] rho23; n = mymom[0, 0] n; u[[1]] = up[[1]]; 
   u[[2]] = up[[1]]; u[[3]] = up[[1]]; u[[4]] = up[[2]]; 
   u[[5]] = up[[2]]; u[[6]] = up[[2]]; u[[7]] = up[[3]]; 
   u[[8]] = up[[3]]; u[[9]] = up[[3]]; u = bu + u; 
   v[[1]] = Vf[[1]] + vp21; v[[2]] = Vf[[1]] + vp22; 
   v[[3]] = Vf[[1]] + vp23; v[[4]] = Vf[[2]] + vp21; 
   v[[5]] = Vf[[2]] + vp22; v[[6]] = Vf[[2]] + vp23; 
   v[[7]] = Vf[[3]] + vp21; v[[8]] = Vf[[3]] + vp22; 
   v[[9]] = Vf[[3]] + vp23; v = bv + v; 
   Return[{n, u, v}, Module];
 ];
   
hyqmom[momin_, qmax_: 0] := 
  Module[{moms = momin}, 
   If[Length[moms] == 3, Return[hyqmom2[moms], Module]]; 
   If[Length[moms] == 5, Return[hyqmom3[moms, qmax], Module]]; 
   Print["Error!"];];
   
hyqmom2[momin_] := 
  Module[{moms = momin, n, u, bu, d2, c2}, n = Table[0, {i, 2}]; 
   u = n; bu = moms[[2]]/moms[[1]]; d2 = moms[[3]]/moms[[1]]; 
   c2 = d2 - bu^2; n[[1]] = moms[[1]]/2; n[[2]] = moms[[1]]/2; 
   u[[1]] = bu - Sqrt[c2]; u[[2]] = bu + Sqrt[c2]; 
   Return[{n, u}, Module];
];
   
hyqmom3[momin_, qmax_] := 
  Module[{moms = momin, etasmall, verysmall, realizable, realsmall, n,
     u, bu, d2, d3, d4, c2, c3, c4, q, eta, slope, det, qp, qm, 
    scales, srho, err, mo, scale, dem, prod, rho, ups, up}, 
   etasmall = 10^(-10); verysmall = 10^(-14); realsmall = 10^(-14); 
   n = Table[0, {i, 3}]; u = n; 
   If[moms[[1]] <= verysmall, n[[2]] = moms[[1]]; 
    Return[{n, u}, Module];]; bu = moms[[2]]/moms[[1]]; 
   d2 = moms[[3]]/moms[[1]]; d3 = moms[[4]]/moms[[1]]; 
   d4 = moms[[5]]/moms[[1]]; c2 = d2 - bu^2; 
   c3 = d3 - 3*bu d2 + 2 bu^3; c4 = d4 - 4*bu d3 + 6 bu^2 d2 - 3 bu^4;
    realizable = c2 c4 - c2^3 - c3^2; 
   If[c2 < 0, 
    If[c2 < -verysmall, Print["c2 negative in three node hyqmom"]; 
     Abort[];]; c2 = 0; c3 = 0; c4 = 0;, 
    If[realizable < 0, 
      If[c2 >= etasmall, q = c3/Sqrt[c2]/c2; eta = c4/c2/c2; 
        If[Abs[q] > verysmall, slope = (eta - 3)/q; det = 8 + slope^2;
          qp = 0.5 (slope + Sqrt[det]); qm = 0.5 (slope - Sqrt[det]); 
         If[Sign[q] == 1, q = qp, q = qm];, q = 0;]; eta = q^2 + 1; 
        c3 = q Sqrt[c2] c2; c4 = eta c2^2; 
        If[realizable < -10^-6, 
         Print["c4 too small in three node hyqmom"]; Abort[];];, 
        c3 = 0; c4 = c2^2;];];]; scale = Sqrt[c2]; 
   If[c2 >= etasmall, q = c3/Sqrt[c2]/c2; eta = c4/c2/c2;, q = 0; 
    eta = 1;]; 
   If[q^2 > qmax^2, slope = (eta - 3)/q; q = qmax Sign[q]; 
    eta = 3 + slope q; realizable = eta - 1 - q^2; 
    If[realizable < 0, eta = 1 + q^2];]; 
   ups[1] = (q - Sqrt[4 eta - 3 q^2])/2; ups[2] = 0; 
   ups[3] = (q + Sqrt[4 eta - 3 q^2])/2; dem = 1/Sqrt[4 eta - 3 q^2]; 
   prod = -ups[1] ups[3]; prod = Max[prod, 1 + realsmall]; 
   rho[1] = -dem/ups[1]; rho[2] = 1 - 1/prod; rho[3] = dem/ups[3]; 
   srho = Sum[rho[i], {i, 3}]; Do[rho[i] = rho[i]/srho, {i, 3}]; 
   scales = Sum[rho[i] ups[i]^2, {i, 3}]/Sum[rho[i], {i, 3}]; 
   Do[up[i] = ups[i] scale/Sqrt[scales], {i, 3}]; 
   If[Min[Table[rho[i], {i, 3}]] < 0, 
    Print["Negative weight in HyQMOM"]; Abort[];]; 
    n[[1]] = rho[1]; 
   n[[2]] = rho[2]; n[[3]] = rho[3]; n = moms[[1]] n; u[[1]] = up[1]; 
   u[[2]] = up[2]; u[[3]] = up[3]; u = bu + u; 
   mo = Table[Sum[n[[i]] pow[u[[i]], j], {i, 1, 3}], {j, 0, 4}]; 
   err = mo - moms; 
   If[err[[3]] > 10^-6 && moms[[1]] > verysmall, Print["err ", err];];
   Return[{n, u}, Module];
];
    
cqmom12[moms_, ks_, nr_, nrd_, nro_:0,wro_:0,ro_:0] := 
  Module[{},
   If[Length[ks[[1]]] == 2, 
    Return[cqmom12m[moms, ks, nr, nrd], Module]]; 
   If[Length[ks[[1]]] == 3, 
    Return[cqmom12p[moms,ks,nr,nrd,nro,wro,ro], Module]];
   ];

cqmom21[moms_, ks_, nr_, nrd_, nro_:0,wro_:0,ro_:0] := 
  Module[{}, 
   If[Length[ks[[1]]] == 2, 
    Return[cqmom21m[moms, ks, nr, nrd], Module]]; 
   If[Length[ks[[1]]] == 3, 
    Return[cqmom21p[moms, ks, nr, nrd, nro,wro,ro], Module]];
];
    
cqmom12m[moms_, ks_, nr_, nrd_]:=Module[{pm, mymom, eqns, vars, linsolv, mRs, xi, w, v, r1, ms, 
    condmoms, condmomvec, mout, xis, ws, wtot, wout, xix, xiy}, 
   mymom=pointer[moms, ks]; 
   mRs=Table[mymom[i,0], {i,0,2nr-1}]; 
   {xi,w}=wheeler[mRs, nr]; 
   v=Table[Table[xi[[i]]^j,{i,nr}],{j,0,nr-1}]; 
   r1=DiagonalMatrix[w]; 
   Do[ms[i]=Table[mymom[j, i], {j, 0, nr - 1}], {i, 0, 2 nrd - 1}]; 
   Do[condmoms[i]=Re[LinearSolve[v.r1, ms[i]]], {i, 0, 2 nrd - 1}]; 
   Do[
   condmomvec[j]=Table[condmoms[i][[j]], {i, 0, 2 nrd - 1}]; 
   If[Norm[Im[condmomvec[j]]] > 0, Print["Imag cond momvec"];];
   , {j, nr}]; 
   Do[{xis[i], ws[i]} = wheeler[condmomvec[i], nrd], {i, nr}];
   Do[wtot[j,i]=w[[j]] ws[j][[i]]; 
   If[
   Re[wtot[j, i]]<0 || Im[wtot[j, i]] != 0, 
   Print["Negative weights"]; Abort[];
   ];
   ,{j,nr},{i,nrd}]; 
   wout = Flatten[Table[wtot[i, j], {i, nr}, {j, nrd}]]; 
   xix = Flatten[Table[xi[[i]], {i, nr}, {j, nrd}]]; 
   xiy = Flatten[Table[xis[i][[j]], {i, nr}, {j, nrd}]]; 
   Return[{wtot, xi, xis}, Module];
   ];
   
cqmom12p[moms_, ks_, nr_, nrd_, nro_,wRos_,Ros_] :=Module[{mymom, mRs, xi, w, v, r1, ms, 
	condmoms,condmomvec,xis,ws,
	wtot,wout,xix,xiy}, 
	Print["Start poly"];
	mymom=pointer[moms,ks,wRos,Ros]; 
	Do[
	mRs=Table[mymom[i,0,l-1],{i,0,2nr-1}];
	{xi[l],w}=wheeler[mRs, nr]; 
    v=Table[Table[xi[l][[i]]^j, {i,nr}],{j,0,nr-1}]; 
    r1=DiagonalMatrix[w];
    Do[ms[m]=Table[mymom[i,m,l-1],{i,0,nr-1}],{m,0,2nrd-1}]; 
    Do[condmoms[i]=LinearSolve[v.r1,ms[i]],{i,0,2 nrd-1}]; 
    Do[condmomvec[j]=Table[condmoms[i][[j]],{i,0,2nrd-1}],{j,nr}]; 
    Do[{xis[l][i],ws[i]}=wheeler[condmomvec[i],nrd],{i,nr}]; 
    Do[wtot[l][j,i]=w[[j]] ws[j][[i]],{i,nrd},{j,nr}];
    ,{l,nro}]; 
    Return[{wtot, xi, xis}, Module];
];
      
cqmom21m[moms_, ks_, nr_, nrd_] := 
  Module[{pm, mymom, eqns, vars, linsolv, mRds, xi, w, v, r1, ms, 
    condmoms, condmomvec, mout, xis, ws, wtot, wout, xix, xiy}, 
   mymom = pointer[moms, ks]; 
   mRds = Table[mymom[0, i], {i, 0, 2 nrd - 1}]; 
   If[Norm[Im[mRds]] > 0, Print["Imag momvec 2"];]; 
   {xi, w} = wheeler[mRds, nrd]; 
   v = Table[Table[xi[[i]]^j, {i, 1, nrd}], {j, 0, nrd - 1}]; 
   r1 = DiagonalMatrix[w]; 
   Do[ms[i] = Table[mymom[i, j], {j, 0, nrd - 1}], {i, 0, 2 nr - 1}]; 
   Do[condmoms[i] = Re[LinearSolve[v.r1, ms[i]]], {i, 0, 2 nr - 1}]; 
   Do[condmomvec[j] = Table[condmoms[i][[j]], {i, 0, 2 nr - 1}]; 
   If[Norm[Im[condmomvec[j]]] > 0, Print["Imag cond momvec"];];, {j, nrd}]; 
   Do[{xis[i], ws[i]} = wheeler[condmomvec[i], nr], {i, nrd}];
   Do[wtot[j, i] = w[[j]] ws[j][[i]]; 
   If[Re[wtot[j, i]] < 0 || Im[wtot[j, i]] != 0, 
   Print["Negative weights"]; Abort[];];, {j, nrd}, {i, nr}]; 
   wout = Flatten[Table[wtot[j, i], {i, nr}, {j, nrd}]]; 
   xix = Flatten[Table[xis[j][[i]], {i, nr}, {j, nrd}]]; 
   xiy = Flatten[Table[xi[[j]], {i, nr}, {j, nrd}]]; 
   Return[{wtot, xi, xis}, Module];
   ];

End[];
EndPackage[];
