(* ::Package:: *)

BeginPackage["QBMMlib`"];

pow::usage="
    pow[x,y] returns x^y where x^0 = 1, even if x=0
";

getterms::usage="
    getterms[eqn,invars,v4,idx] computes the coefficients and exponents of the internal coordinates that correspond to the RHS operator of the moment transport equations.
    'eqn' is the governing ODE, 'invars' are the internal coordinates, 'v4' is the acceleration term, and 'idx' are the indices.
";

momidx::usage = "
    momidx[nr,nrd,method,numperm,nro] computes the moments required by the 'method' (CQMOM or CHyQMOM) for a given number of nodes in the first, second, and third coordinate directions nr, nrd, and nro. 
    Extra moments are included if the number of permuations, 'numperm', requires it.
    Constraints: CHyQMOM requires nr=nrd and nr=2 or 3. 
    CQMOM does not have these constraints.
";

quad::usage="
    quad[w,xi,xis,method,nr,nrd,ks,perm,nro,wRos,Ros] computes the moments M_{lmn} = int[x^l y^m z^n P dx dy dz] with indices 'ks = {l,m,n}' via the weights 'w' and quadrature points in the first 'xi' and second 'xis' coordinate directions.
    An optional third coordinate direction Ro can be included, which has 'nro' number of nodes with weights 'wRos' and abscissa 'Ros'.
    Thus, 'nro', 'wRos', and 'Ros' are optional arguments.
    'perm' corresponds to the permutation ordering (12: v1 | v2) or (21: v2 | v1), which switches which abscissa are considered the first direction 'xi' versus 'xis'.
    The 'method' input takes either 'CQMOM' or 'CHyQMOM', though the operations performed are essentially identical.
";

cqmom::usage = "
    cqmom[perm,mom,ks,nr,nrd,nro,wro,ro] performs moment inversion for the optimal weights and abscissa as computed by the CQMOM algorithm [Yuan & Fox, JCP 2011].
    The input moments are 'mom', which correspond to moment indices 'ks'. 
    'perm' is the permutation [=12 for v1 | v2 or = 21 for v2 | v1].
    'nr' and 'nrd' are the number of nodes in the first and second internal coordinate directions.
    A third, static distribution can optionally be included via its number of quadrature points 'nro' and their locations 'ro' and weights 'wro'.
";

chyqmom::usage = "
    chyqmom[mom,k,q,wRo,Ros] performs moment inversion for the optimal weights and abscissa as computed by the CHyQMOM algorithm [Fox et al., JCP 2018; Patel et al., JCPX 2019].
    The input moments are 'mom', which correspond to moment indices 'k'. 
    For three-node closure in any of the coordinate directions (e.g. nr=3 or nrd=3), an optional parameter 'q' is available to limit the amount of skewness that direction. 
    A third, static distribution can optionally be included via its number of quadrature points nro and their locations 'ro' and weights 'wro'.
";

wheeler::usage=""

Begin["`Private`"];

pow[x_, 0] := 1 /; x == 0 || x == 0.
pow[x_, y_: 0] := x^y

getterms[eqn_,invars_,v4_,idx_]:=Module[{v,integrand,dim,mrdd,list,exp,coefs,exps,vars,l,m,n},
    dim=Length[invars];
    If[dim!=2&&dim!=3,Print["Incorrect dimesionallity ",dim];Abort[]];
    If[dim==2,vars={v[1],v[2]};{l,m}=idx;];
    If[dim==3,vars={v[1],v[2],v[3]};{l,m,n}=idx;];

    Do[v[i]=invars[[i]],{i,dim}];
    mrdd=v4/.Solve[eqn,v4][[1]];

    If[dim==2,integrand=mrdd v[1]^l v[2]^(m-1)];
    If[dim==3,integrand=mrdd v[1]^l v[2]^(m-1) v[3]^n];

    list=PowerExpand[List@@Distribute[integrand]];
    If[dim==2,
        exps=Table[Exponent[list[[i]],j],{i,Length[list]},{j,invars}];
        (* Pad with a zero to make 3D *)
        Do[AppendTo[exps[[i]],0],{i,Length[list]}];
    ];
    If[dim==3,exps=Table[Exponent[list[[i]],j],{i,Length[list]},{j,invars}]];
    coefs=m Table[Coefficient[list[[i]],Product[v[j]^exps[[i,j]],{j,dim}]],{i,Length[list]}];

    (* Add term that isn't in the integrand *)
    If[dim==2,AppendTo[exps,{l-1,m+1,0}]];
    If[dim==3,AppendTo[exps,{l-1,m+1,n}]];
    AppendTo[coefs,l];
    Return[{coefs,exps},Module];
];

quad[w_,xi_,xis_,method_,nr_,nrd_,ks_,perm_,wRos_:0,Ros_:0]:=Module[{momq,momc,moms,momsp,nro},
    nro=Length[wRos];
	If[method=="CQMOM",
		If[nro==0,
			If[perm==12,momq[p_,q_]:=Sum[w[j,i] xi[[j]]^p xis[j][[i]]^q,{j,nr},{i,nrd}]];
			If[perm==21,momq[p_,q_]:=Sum[w[j,i] xi[[j]]^q xis[j][[i]]^p,{j,nrd},{i,nr}]];
			moms=Table[momq[ks[[i,1]],ks[[i,2]]], {i,Length[ks]}]; 
			(* momsp=Table[momq[ksp[[i,1]],ksp[[i,2]]], {i,Length[ksp]}]; *)
		];
		If[nro>0,
			If[perm==12,
				momq[p_,q_,l_]:=Sum[w[l][j,i] xi[l][[j]]^p xis[l][j][[i]]^q,{i,nrd},{j,nr}];
			];
			If[perm==21,
				momq[p_,q_,l_]:=Sum[w[l][j,i] xi[l][[j]]^q xis[l][j][[i]]^p,{i,nr},{j,nrd}];
			];
			moms=Table[momq[ks[[i,1]],ks[[i,2]],1+ks[[i,3]]],{i,Length[ks]}]; 
			(* momsp=Table[momq[ksp[[i,1]],ksp[[i,2]],ks[[i,3]]],{i,Length[ksp]}]; *)
		];
	];
	If[method=="CHyQMOM",
		If[nro==0,
			momq[p_,q_]:=Sum[w[[i]]xi[[i]]^p xis[[i]]^q,{i,Length[w]}];
			moms=Table[momq[ks[[i,1]],ks[[i,2]]], {i,Length[ks]}]; 
			momsp=Table[momq[ksp[[i,1]],ksp[[i,2]]],{i,Length[ksp]}];
		];
		If[nro>0,
            (* gets full moments  *)
			(* 
            momc[p_,q_,l_]:=Sum[w[l][[i]]xi[l][[i]]^p xis[l][[i]]^q,{i,Length[w[l]]}];
			momq[p_,q_,r_]:=Sum[wRos[[l]]Ros[[l]]^r momc[p,q,l],{l,nro}];
			moms=Table[momq[ks[[i,1]],ks[[i,2]],ks[[i,3]]],{i,Length[ks]}]; 
			momsp=Table[momq[ksp[[i,1]],ksp[[i,2]],ks[[i,3]]],{i,Length[ksp]}];
            *) 

            (* get conditioned on R_o,l moments, +1 included on last index to promote 0 -> 1 *)
			momq[p_,q_,l_]:=Sum[w[l][[i]]xi[l][[i]]^p xis[l][[i]]^q,{i,Length[w[l]]}];
			moms=Table[momq[ks[[i,1]],ks[[i,2]],1+ks[[i,3]]],{i,Length[ks]}]; 
            (* these are no longer computable/correct. need to compute them after the fact *)
			(* momsp=Table[momq[ks[[i,1]],ks[[i,2]],1+ks[[i,3]]],{i,Length[ksp]}];  *)
		];
	];
    (* moms and momsp are projected moments, momq is moment function (via quadrature) for any power pqr *)
	Return[{moms,momq},Module];
];

(* wheeler::usage="
    wheeler[mom] computes the one-dimensional optimal weights and abscissa given moments 'mom' using Wheeler's algorithm [Wheeler, Rocky Mt. J. Math. 1974].
"; *)
wheeler[mom_] := Module[{mm = mom, nn, sig, a, b, Ja, w, xi, eval, evec, esys}, 
    nn=Length[mm]/2;
    sig = Table[0., {i, 2 nn}, {j, 2 nn}]; 

    Do[sig[[2, i + 1]] = mm[[i + 1]], {i, 0, 2 nn - 1}]; 
    a = Table[0., {i, nn}]; b = a; a[[1]] = mm[[2]]/mm[[1]]; 
    b[[1]] = 0; 
    Do[
        Do[
            sig[[i + 2, j + 1]] = sig[[i + 1, j + 2]] - a[[i]] sig[[i + 1, j + 1]] - b[[i]] sig[[i, j + 1]];
            a[[i + 1]] = -(sig[[i + 1, i + 1]]/sig[[i + 1, i]]) + sig[[i + 2, i + 2]]/sig[[i + 2, i + 1]]; 
            b[[i + 1]] = sig[[i + 2, i + 1]]/sig[[i + 1, i]];
        ,{j,i, 2 nn - i - 1}];
    ,{i,nn - 1}];
    Ja = DiagonalMatrix[a]; 
    Do[
        Ja[[i, i + 1]] = -Sqrt[Abs[b[[i + 1]]]]; 
        Ja[[i + 1, i]] = -Sqrt[Abs[b[[i + 1]]]];
    ,{i,nn-1}]; 
    esys = Eigensystem[Ja]; 
    eval = esys[[1]]; 
    evec = esys[[2]]; 
    w = Table[evec[[i,1]]^2 mm[[1]],{i,nn}]; 
    Return[{eval, w}, Module];
];

momidx[nr_, nrd_, method_, numperm_: 2, nro_: 0] := Module[{ks, k1, k2,kstemp}, 
    If[method == "CHyQMOM", 
        If[nr == 2, 
            If[
                nro == 0, 
                ks = {{0, 0}, {1, 0}, {0, 1}, {2, 0}, {1, 1}, {0, 2}};
            ,
                ks = Flatten[Table[{{0, 0, i}, {1, 0, i}, {0, 1, i}, {2, 0, i}, {1, 1, i}, {0, 2, i}}, {i, 0, nro-1}], 1];
            ];
        ]; 
        If[nr == 3, 
            If[
                nro == 0, 
                ks = {{0, 0}, {1, 0}, {0, 1}, {2, 0}, {1, 1}, {0, 2}, {3, 0}, {0,3}, {4, 0}, {0, 4}}
            ,
                ks = Flatten[Table[{{0, 0, i}, {1, 0, i}, {0, 1, i}, {2, 0, i}, {1, 1, i}, {0, 2, i}, {3, 0, i}, {0,3, i}, {4, 0, i}, {0, 4, i}}, {i,0,nro-1}],1];
            ];
        ]; 
    ]; 
    If[method == "CQMOM", 
        If[nro == 0, 
            k1 = Flatten[{Table[{q,p}, {q, 0, nr - 1}, {p, 0, 2 nrd - 1}], Table[{q, p}, {q, nr, 2 nr - 1}, {p, 0, 0}]}, 2]; 
            k2 = Flatten[{Table[{p,q}, {q, 0, nrd - 1}, {p, 0, 2 nr - 1}], Table[{p, q}, {q, nrd, 2 nrd - 1}, {p, 0, 0}]}, 2]; 
            If[
                numperm==2,
                ks = DeleteDuplicates[Join[k1, k2]]; 
            ,
                ks = DeleteDuplicates[k1]; 
            ];
        ,
            Do[
                k1[i] = Flatten[{Table[{q,p,i}, {q, 0, nr - 1}, {p, 0, 2 nrd - 1}], Table[{q,p,i}, {q, nr, 2 nr - 1}, {p, 0, 0}]}, 2]; 
                k2[i] = Flatten[{Table[{p,q,i}, {q, 0, nrd - 1}, {p, 0, 2 nr - 1}], Table[{p,q,i}, {q, nrd, 2 nrd - 1}, {p, 0, 0}]}, 2]; 
                If[
                    numperm==2,
                    kstemp[i] = DeleteDuplicates[Join[k1[i],k2[i]]]; 
                ,
                    kstemp[i] = DeleteDuplicates[k1[i]]; 
                ];
            ,{i,0,nro-1}];
            ks=Flatten[Join[Table[kstemp[i],{i,0,nro-1}]],1];
        ];
    ]; 
    Return[ks, Module];
];
      
pointer[moms_,ks_,w_:0,ro_:0] := Module[{eqns, vars, linsolv, pm, mymom}, 
    If[
        Length[ks[[1]]] == 2, 
        eqns = Table[moms[[i]] == pm[ks[[i, 1]], ks[[i, 2]]], {i,Length[ks]}]; 
        vars = DeleteDuplicates[ Flatten[Table[pm[ks[[i, 1]], ks[[i, 2]]], {i,Length[ks]}]]];
        linsolv = First[vars/.Solve[eqns, vars]];
        mymom[q_, p_] := linsolv[[First[First[Position[ks,{q,p}]]]]]; 
        Return[mymom, Module];
    ];
    If[
        Length[ks[[1]]] == 3,
        (* converts vector of total moments m_lmn to conditioned moments m_lm @ each R_o,k  *)
        (* eqns = Table[moms[[i]] == Sum[ w[[l + 1]] ro[[l + 1]]^ks[[i, 3]] pm[ks[[i, 1]], ks[[i, 2]], l], {l,0,Length[ro]-1}], {i,Length[ks]}]; 
        vars = DeleteDuplicates[Flatten[Table[pm[ks[[i, 1]], ks[[i, 2]], l], {l,0,Length[ro]-1},{i,Length[ks]}]]]; 
        linsolv = First[vars/.Solve[eqns, vars]]; 
        mymom[q_, p_, r_] := linsolv[[First[First[Position[ks,{q,p,r}]]]]]; 
        *)

        (* just tells you where to find each moments {qp,k}}, for each R_o,k  *)
        eqns = Table[moms[[i]] == pm[ks[[i, 1]], ks[[i, 2]], ks[[i, 3]]], {i,Length[ks]}]; 
        vars = DeleteDuplicates[ Flatten[Table[pm[ks[[i, 1]], ks[[i, 2]], ks[[i, 3]]], {i,Length[ks]}]]];
        linsolv = First[vars/.Solve[eqns, vars]];
        mymom[q_, p_, r_] := linsolv[[First[First[Position[ks,{q,p,r}]]]]]; 
        Return[mymom, Module];
    ]; 
    Print["Error, k index not 2D or 3D!"];
];
    
chyqmom[momin_, k_, q_, wRo_: 0, Ros_: 0] := Module[{moms = momin, ks = k},
    If[
        Length[wRo]==0,
        If[Mod[Length[moms], 6] == 0,Return[chyqmom4[moms, ks], Module]]; 
        If[Mod[Length[moms], 10] == 0,Return[chyqmom9[moms, ks, q], Module]];
    ,
        If[Mod[Length[moms], 6] == 0, Return[chyqmom4p[moms, ks, wRo, Ros], Module]]; 
        If[Mod[Length[moms], 10] == 0, Return[chyqmom9p[moms, ks, q, wRo, Ros], Module]];   
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
    c20 = d20 - bu^2; 
    c11 = d11 - bu bv;
    c02 = d02 - bv^2; 
    M1 = {1, 0, c20}; 
    {rho, up} = hyqmom[M1]; 
    Vf = c11 up/c20; 
    mu2avg = c02 - Total[rho Vf^2]; 
    mu2avg = Max[mu2avg, 0]; 
    mu2 = mu2avg; 
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
    Return[{n,u,v}, Module];
];

chyqmom4p[momin_,kk_,wRo_,Ros_] := Module[{moms = momin, ks = kk, eqns, vars, linsolv, pm, mymom, n, 
    bu, bv, d20, d11, d02, c20, c11, c02, M1, rho, up, Vf, mu2avg, 
    mu2, M3, rh3, up3, vp21, vp22, rho21, rho22, u, v, nro, nn, uu, vv}, 

    nro = Length[wRo];
    mymom = pointer[moms, ks, wRo, Ros]; 
    (* operating on conditioned moments for each Ro slice *)
    Do[
        n = Table[0, {i, 4}]; u = n; v = n; 
        bu  = mymom[1, 0, l-1]/mymom[0, 0, l-1]; 
        bv  = mymom[0, 1, l-1]/mymom[0, 0, l-1]; 
        d20 = mymom[2, 0, l-1]/mymom[0, 0, l-1]; 
        d11 = mymom[1, 1, l-1]/mymom[0, 0, l-1]; 
        d02 = mymom[0, 2, l-1]/mymom[0, 0, l-1]; 
        c20 = d20 - bu^2; 
        c11 = d11 - bu bv;
        c02 = d02 - bv^2; 
        M1 = {1, 0, c20}; 
        {rho, up} = hyqmom[M1]; 
        Vf = c11 up/c20; 
        mu2avg = c02 - Total[rho Vf^2]; 
        mu2avg = Max[mu2avg, 0]; 
        mu2 = mu2avg; 
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
        n = mymom[0, 0, l-1] n; 

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

        nn[l] = n;
        uu[l] = u;
        vv[l] = v;
    ,{l,nro}];
    Return[{nn,uu,vv}, Module];
];
 

chyqmom9p[momin_, kk_, qmax_, wRo_, Ros_] := 
  Module[{moms = momin, ks = kk, vars, linsolv, mymom, n, csmall, 
    verysmall, bu, bv, d20, d11, d02, d30, d03, d40, d04, c20, c11, 
    c02, c30, c03, c40, c04, M1, rho, up, Vf, M2, rho2, up2, vp21, 
    vp22, vp23, rho21, rho22, rho23, mu2avg, mu2, mu3, mu4, u, v, 
    eqns, pm, q, eta, slope, det, qm, qp, up3, M3, rh3, nro, nn, uu, vv}, 

    nro = Length[wRo];
    mymom = pointer[moms, ks, wRo, Ros]; 

    Do[
        n = Table[0, {i, 9}]; u = n; v = n; 
        csmall = 10^(-10); 
        verysmall = 10^(-14); 
        bu =  mymom[1, 0,l-1]/mymom[0, 0,l-1]; 
        bv =  mymom[0, 1,l-1]/mymom[0, 0,l-1]; 
        d20 = mymom[2, 0,l-1]/mymom[0, 0,l-1]; 
        d11 = mymom[1, 1,l-1]/mymom[0, 0,l-1]; 
        d02 = mymom[0, 2,l-1]/mymom[0, 0,l-1]; 
        d30 = mymom[3, 0,l-1]/mymom[0, 0,l-1]; 
        d03 = mymom[0, 3,l-1]/mymom[0, 0,l-1]; 
        d40 = mymom[4, 0,l-1]/mymom[0, 0,l-1]; 
        d04 = mymom[0, 4,l-1]/mymom[0, 0,l-1]; 
        c20 = d20 - bu^2; 
        c11 = d11 - bu bv;
        c02 = d02 - bv^2; 
        c30 = d30 - 3 bu d20 + 2 bu^3; 
        c03 = d03 - 3 bv d02 + 2 bv^3; 
        c40 = d40 - 4 bu d30 + 6 bu^2 d20 - 3 bu^4; 
        c04 = d04 - 4 bv d03 + 6 bv^2 d02 - 3 bv^4; 
        M1 = {1, 0, c20, c30, c40}; 
        {rho, up} = hyqmom[M1]; 
        If[c20 <= csmall, 
            rho[[1]] = 0; 
            rho[[2]] = 1; 
            rho[[3]] = 0; 
            Vf = 0 up; 
            M2 = {1, 0, c02, c03, c04}; 
            {rho2, up2} = hyqmom[M2, qmax]; 
            vp21 = up2[[1]]; 
            vp22 = up2[[2]]; 
            vp23 = up2[[3]]; 
            rho21 = rho2[[1]]; 
            rho22 = rho2[[2]]; 
            rho23 = rho2[[3]];
        , 
            Vf = c11 up/c20; 
            mu2avg = c02 - Total[rho Vf^2]; 
            mu2avg = Max[mu2avg, 0]; 
            mu2 = mu2avg; 
            mu3 = 0 mu2; 
            mu4 = mu2^2; 
            If[mu2 > csmall, 
                q = (c03 - Total[rho Vf^3])/mu2^(3/2); 
                eta = (c04 - Total[rho Vf^4] - 6 Total[rho Vf^2] mu2)/mu2^2; 
                If[eta < q^2 + 1, 
                    If[Abs[q] > verysmall, 
                        slope = (eta - 3)/q; 
                        det = 8 + slope^2; 
                        qp = 0.5 (slope + Sqrt[det]); 
                        qm = 0.5 (slope - Sqrt[det]); 
                        If[Sign[q] == 1, q = qp, q = qm];
                    , 
                        q = 0
                    ]; 
                    eta = q^2 + 1;
                ]; 
                mu3 = q mu2^(3/2); 
                mu4 = eta mu2^2;
            ]; 
            M3 = {1, 0, mu2, mu3, mu4}; 
            {rh3, up3} = hyqmom[M3, qmax]; 
            vp21 = up3[[1]]; 
            vp22 = up3[[2]]; 
            vp23 = up3[[3]]; 
            rho21 = rh3[[1]]; 
            rho22 = rh3[[2]]; 
            rho23 = rh3[[3]];
        ]; 

        n[[1]] = rho[[1]] rho21; 
        n[[2]] = rho[[1]] rho22; 
        n[[3]] = rho[[1]] rho23; 
        n[[4]] = rho[[2]] rho21; 
        n[[5]] = rho[[2]] rho22; 
        n[[6]] = rho[[2]] rho23; 
        n[[7]] = rho[[3]] rho21; 
        n[[8]] = rho[[3]] rho22; 
        n[[9]] = rho[[3]] rho23; 
        n = mymom[0, 0, l-1] n; 

        u[[1]] = up[[1]]; 
        u[[2]] = up[[1]];
        u[[3]] = up[[1]]; 
        u[[4]] = up[[2]]; 
        u[[5]] = up[[2]];
        u[[6]] = up[[2]];
        u[[7]] = up[[3]]; 
        u[[8]] = up[[3]];
        u[[9]] = up[[3]]; 
        u = bu + u; 

        v[[1]] = Vf[[1]] + vp21; 
        v[[2]] = Vf[[1]] + vp22; 
        v[[3]] = Vf[[1]] + vp23; 
        v[[4]] = Vf[[2]] + vp21; 
        v[[5]] = Vf[[2]] + vp22;
        v[[6]] = Vf[[2]] + vp23; 
        v[[7]] = Vf[[3]] + vp21;
        v[[8]] = Vf[[3]] + vp22; 
        v[[9]] = Vf[[3]] + vp23; 
        v = bv + v; 

        nn[l] = n;
        uu[l] = u;
        vv[l] = v;
    ,{l,nro}];
    Return[{nn,uu,vv}, Module];
];

chyqmom9[momin_, kk_, qmax_] := 
  Module[{moms = momin, ks = kk, vars, linsolv, mymom, n, csmall, 
    verysmall, bu, bv, d20, d11, d02, d30, d03, d40, d04, c20, c11, 
    c02, c30, c03, c40, c04, M1, rho, up, Vf, M2, rho2, up2, vp21, 
    vp22, vp23, rho21, rho22, rho23, mu2avg, mu2, mu3, mu4, u, v, 
    eqns, pm, q, eta, slope, det, qm, qp, up3, M3, rh3}, 
    mymom = pointer[moms, ks]; 
    n = Table[0, {i, 9}]; u = n; v = n; 
    csmall = 10^(-10); 
    verysmall = 10^(-14); 
    bu = mymom[1, 0]/mymom[0, 0]; 
    bv = mymom[0, 1]/mymom[0, 0]; 
    d20 = mymom[2, 0]/mymom[0, 0]; 
    d11 = mymom[1, 1]/mymom[0, 0]; 
    d02 = mymom[0, 2]/mymom[0, 0]; 
    d30 = mymom[3, 0]/mymom[0, 0]; 
    d03 = mymom[0, 3]/mymom[0, 0]; 
    d40 = mymom[4, 0]/mymom[0, 0]; 
    d04 = mymom[0, 4]/mymom[0, 0]; 
    c20 = d20 - bu^2; 
    c11 = d11 - bu bv;
    c02 = d02 - bv^2; 
    c30 = d30 - 3 bu d20 + 2 bu^3; 
    c03 = d03 - 3 bv d02 + 2 bv^3; 
    c40 = d40 - 4 bu d30 + 6 bu^2 d20 - 3 bu^4; 
    c04 = d04 - 4 bv d03 + 6 bv^2 d02 - 3 bv^4; 
    M1 = {1, 0, c20, c30, c40}; 
    {rho, up} = hyqmom[M1]; 
    If[c20 <= csmall, 
        rho[[1]] = 0; 
        rho[[2]] = 1; 
        rho[[3]] = 0; 
        Vf = 0 up; 
        M2 = {1, 0, c02, c03, c04}; 
        {rho2, up2} = hyqmom[M2, qmax]; 
        vp21 = up2[[1]]; 
        vp22 = up2[[2]]; 
        vp23 = up2[[3]]; 
        rho21 = rho2[[1]]; 
        rho22 = rho2[[2]]; 
        rho23 = rho2[[3]];
    , 
        Vf = c11 up/c20; 
        mu2avg = c02 - Total[rho Vf^2]; 
        mu2avg = Max[mu2avg, 0]; 
        mu2 = mu2avg; 
        mu3 = 0 mu2; 
        mu4 = mu2^2; 
        If[mu2 > csmall, 
            q = (c03 - Total[rho Vf^3])/mu2^(3/2); 
            eta = (c04 - Total[rho Vf^4] - 6 Total[rho Vf^2] mu2)/mu2^2; 
            If[eta < q^2 + 1, 
                If[Abs[q] > verysmall, 
                    slope = (eta - 3)/q; 
                    det = 8 + slope^2; 
                    qp = 0.5 (slope + Sqrt[det]); 
                    qm = 0.5 (slope - Sqrt[det]); 
                    If[Sign[q] == 1, q = qp, q = qm];
                , 
                    q = 0
                ]; 
                eta = q^2 + 1;
            ]; 
            mu3 = q mu2^(3/2); 
            mu4 = eta mu2^2;
        ]; 
        M3 = {1, 0, mu2, mu3, mu4}; 
        {rh3, up3} = hyqmom[M3, qmax]; 
        vp21 = up3[[1]]; 
        vp22 = up3[[2]]; 
        vp23 = up3[[3]]; 
        rho21 = rh3[[1]]; 
        rho22 = rh3[[2]]; 
        rho23 = rh3[[3]];
    ]; 

    n[[1]] = rho[[1]] rho21; 
    n[[2]] = rho[[1]] rho22; 
    n[[3]] = rho[[1]] rho23; 
    n[[4]] = rho[[2]] rho21; 
    n[[5]] = rho[[2]] rho22; 
    n[[6]] = rho[[2]] rho23; 
    n[[7]] = rho[[3]] rho21; 
    n[[8]] = rho[[3]] rho22; 
    n[[9]] = rho[[3]] rho23; 
    n = mymom[0, 0] n; 

    u[[1]] = up[[1]]; 
    u[[2]] = up[[1]];
    u[[3]] = up[[1]]; 
    u[[4]] = up[[2]]; 
    u[[5]] = up[[2]];
    u[[6]] = up[[2]];
    u[[7]] = up[[3]]; 
    u[[8]] = up[[3]];
    u[[9]] = up[[3]]; 
    u = bu + u; 

    v[[1]] = Vf[[1]] + vp21; 
    v[[2]] = Vf[[1]] + vp22; 
    v[[3]] = Vf[[1]] + vp23; 
    v[[4]] = Vf[[2]] + vp21; 
    v[[5]] = Vf[[2]] + vp22;
    v[[6]] = Vf[[2]] + vp23; 
    v[[7]] = Vf[[3]] + vp21;
    v[[8]] = Vf[[3]] + vp22; 
    v[[9]] = Vf[[3]] + vp23; 
    v = bv + v; 
    Return[{n, u, v}, Module];
 ];
   
hyqmom[momin_, qmax_: 0] := Module[{moms = momin}, 
   If[Length[moms] == 3, Return[hyqmom2[moms], Module]]; 
   If[Length[moms] == 5, Return[hyqmom3[moms, qmax], Module]]; 
   Print["Error!"];
];
   
hyqmom2[momin_] := Module[{moms = momin, n, u, bu, d2, c2}, 
    n = Table[0, {i, 2}]; 
    u = n; 
    bu = moms[[2]]/moms[[1]]; 
    d2 = moms[[3]]/moms[[1]]; 
    c2 = d2 - bu^2; 
    n[[1]] = moms[[1]]/2; 
    n[[2]] = moms[[1]]/2; 
    u[[1]] = bu - Sqrt[c2]; 
    u[[2]] = bu + Sqrt[c2]; 
    Return[{n, u}, Module];
];
   
hyqmom3[momin_, qmax_] := Module[{moms = momin, etasmall, verysmall, 
    realizable, realsmall, n,
    u, bu, d2, d3, d4, c2, c3, c4, q, eta, slope, det, qp, qm, 
    scales, srho, err, mo, scale, dem, prod, rho, ups, up}, 

    etasmall = 10^(-10); 
    verysmall = 10^(-14); 
    realsmall = 10^(-14); 
    n = Table[0, {i, 3}]; 
    u = n; 
    If[moms[[1]] <= verysmall, 
        n[[2]] = moms[[1]]; 
        Return[{n, u}, Module];
    ]; 
    bu = moms[[2]]/moms[[1]]; 
    d2 = moms[[3]]/moms[[1]]; 
    d3 = moms[[4]]/moms[[1]]; 
    d4 = moms[[5]]/moms[[1]]; 
    c2 = d2 - bu^2; 
    c3 = d3 - 3*bu d2 + 2 bu^3; 
    c4 = d4 - 4*bu d3 + 6 bu^2 d2 - 3 bu^4;
    realizable = c2 c4 - c2^3 - c3^2; 
    If[c2 < 0, 
        If[c2 < -verysmall, Print["c2 negative in three node hyqmom"];Abort[];]; 
        c2 = 0; 
        c3 = 0; 
        c4 = 0;
    , 
        If[realizable < 0, 
              If[c2 >= etasmall, 
                    q = c3/Sqrt[c2]/c2; 
                    eta = c4/c2/c2; 
                    If[
                        Abs[q] > verysmall, 
                        slope = (eta - 3)/q; 
                        det = 8 + slope^2;
                        qp = 0.5 (slope + Sqrt[det]); 
                        qm = 0.5 (slope - Sqrt[det]); 
                        If[Sign[q] == 1, q = qp, q = qm];
                    , 
                        q = 0;
                    ]; 
                    eta = q^2 + 1; 
                    c3 = q Sqrt[c2] c2; c4 = eta c2^2; 
                    If[realizable < -10^-6,Print["c4 small in hyqmom3"];Abort[];];
                , 
                    c3 = 0; 
                    c4 = c2^2;
              ];
        ];
    ]; 
    scale = Sqrt[c2]; 
    If[c2 >= etasmall, 
        q = c3/Sqrt[c2]/c2; 
        eta = c4/c2/c2;
    , 
        q = 0; 
        eta = 1;
    ]; 
    If[
        q^2 > qmax^2, 
        slope = (eta - 3)/q; 
        q = qmax Sign[q]; 
        eta = 3 + slope q; 
        realizable = eta - 1 - q^2; 
        If[realizable < 0, eta = 1 + q^2];
    ]; 
    ups[1] = (q - Sqrt[4 eta - 3 q^2])/2; 
    ups[2] = 0; 
    ups[3] = (q + Sqrt[4 eta - 3 q^2])/2; 
    dem = 1/Sqrt[4 eta - 3 q^2]; 
    prod = -ups[1] ups[3]; 
    prod = Max[prod, 1 + realsmall]; 
    rho[1] = -dem/ups[1]; 
    rho[2] = 1 - 1/prod; 
    rho[3] = dem/ups[3]; 
    srho = Sum[rho[i], {i, 3}]; 
    Do[rho[i] = rho[i]/srho, {i, 3}]; 
    scales = Sum[rho[i] ups[i]^2, {i, 3}]/Sum[rho[i], {i, 3}]; 
    Do[up[i] = ups[i] scale/Sqrt[scales], {i, 3}]; 
    If[
        Min[Table[rho[i], {i, 3}]] < 0, 
        Print["Negative weight in HyQMOM"]; Abort[];
    ]; 
    n[[1]] = rho[1]; 
    n[[2]] = rho[2]; 
    n[[3]] = rho[3]; 
    n = moms[[1]] n; 

    u[[1]] = up[1]; 
    u[[2]] = up[2]; 
    u[[3]] = up[3]; 
    u = bu + u; 

    mo = Table[Sum[n[[i]] pow[u[[i]], j], {i, 1, 3}], {j, 0, 4}]; 
    err = mo - moms; 
    If[err[[3]] > 10^-6 && moms[[1]] > verysmall, Print["err ", err];];
    Return[{n, u}, Module];
];
    
cqmom[perm_, moms_, ks_, nr_, nrd_, nro_:0,wro_:0,ro_:0] := Module[{},
    If[perm==12,
        If[Length[ks[[1]]] == 2,Return[cqmom12m[moms, ks, nr, nrd], Module]]; 
        If[Length[ks[[1]]] == 3,Return[cqmom12p[moms,ks,nr,nrd,nro,wro,ro], Module]];
    ];
    If[perm==21,
        If[Length[ks[[1]]] == 2,Return[cqmom21m[moms, ks, nr, nrd], Module]]; 
        If[Length[ks[[1]]] == 3,Return[cqmom21p[moms, ks, nr, nrd, nro,wro,ro], Module]];
    ];
    Print["Error!"];
];

cqmom12m[moms_, ks_, nr_, nrd_]:=Module[{pm, mymom, eqns, vars, linsolv, 
    mRs, xi, w, v, r1, ms, 
    condmoms, condmomvec, mout, xis, ws, wtot}, 
    mymom=pointer[moms, ks]; 
    mRs=Table[mymom[i,0], {i,0,2nr-1}]; 
    {xi,w}=wheeler[mRs]; 
    v=Table[Table[xi[[i]]^j,{i,nr}],{j,0,nr-1}]; 
    r1=DiagonalMatrix[w]; 
    Do[ms[i]=Table[mymom[j, i], {j, 0, nr - 1}], {i, 0, 2 nrd - 1}]; 
    Do[condmoms[i]=Re[LinearSolve[v.r1, ms[i]]], {i, 0, 2 nrd - 1}]; 
    Do[
        condmomvec[j]=Table[condmoms[i][[j]], {i, 0, 2 nrd - 1}]; 
        If[Norm[Im[condmomvec[j]]] > 0, Print["Imag cond momvec"];Abort[];];
    , {j, nr}]; 
    Do[{xis[i], ws[i]} = wheeler[condmomvec[i]], {i, nr}];
    Do[
        wtot[j,i]=w[[j]] ws[j][[i]]; 
        If[
            Re[wtot[j, i]]<0 || Im[wtot[j, i]] != 0, 
            Print["Negative weights"]; Abort[];
        ];
    ,{j,nr},{i,nrd}]; 
    Return[{wtot, xi, xis}, Module];
];
       
cqmom12p[moms_, ks_, nr_, nrd_, nro_,wRos_,Ros_] :=Module[{mymom, mRs, xi, w, v, r1, ms, 
	condmoms,condmomvec,xis,ws,
	wtot}, 
	mymom=pointer[moms,ks,wRos,Ros]; 
	Do[
        mRs=Table[mymom[i,0,l-1],{i,0,2nr-1}];
        {xi[l],w}=wheeler[mRs]; 
        v=Table[Table[xi[l][[i]]^j,{i,nr}],{j,0,nr-1}]; 
        r1=DiagonalMatrix[w];
        Do[ms[m]=Table[mymom[i,m,l-1],{i,0,nr-1}],{m,0,2nrd-1}]; 
        Do[condmoms[i]=LinearSolve[v.r1,ms[i]],{i,0,2 nrd-1}]; 
        Do[
            condmomvec[j]=Table[condmoms[i][[j]],{i,0,2nrd-1}];
            If[Norm[Im[condmomvec[j]]] > 0, Print["Imag cond momvec"];Abort[];];
        ,{j,nr}]; 
        Do[{xis[l][i],ws[i]}=wheeler[condmomvec[i]],{i,nr}]; 
        Do[wtot[l][j,i]=w[[j]] ws[j][[i]],{i,nrd},{j,nr}];
    ,{l,nro}]; 
    Return[{wtot,xi,xis},Module];
];

cqmom21m[moms_, ks_, nr_, nrd_] := Module[{pm, mymom, eqns, vars, 
    linsolv, mRds, xi, w, v, r1, ms, 
    condmoms, condmomvec, mout, xis, ws, wtot}, 
    mymom = pointer[moms, ks]; 
    mRds = Table[mymom[0, i], {i, 0, 2 nrd - 1}]; 
    If[Norm[Im[mRds]] > 0, Print["Imag momvec 2"]]; 
    {xi, w} = wheeler[mRds]; 
    v = Table[Table[xi[[i]]^j, {i, 1, nrd}], {j, 0, nrd - 1}]; 
    r1 = DiagonalMatrix[w]; 
    Do[ms[i] = Table[mymom[i, j], {j, 0, nrd - 1}], {i, 0, 2 nr - 1}]; 
    Do[condmoms[i] = Re[LinearSolve[v.r1, ms[i]]], {i, 0, 2 nr - 1}]; 
    Do[
        condmomvec[j] = Table[condmoms[i][[j]], {i, 0, 2 nr - 1}]; 
        If[Norm[Im[condmomvec[j]]] > 0, Print["Imag cond momvec"];Abort[];];
    , {j, nrd}]; 
    Do[{xis[i], ws[i]} = wheeler[condmomvec[i]], {i, nrd}];
    Do[wtot[j, i] = w[[j]] ws[j][[i]]; 
        If[Re[wtot[j, i]] < 0 || Im[wtot[j, i]] != 0, 
        Print["Negative weights"]; Abort[];];
    ,{j,nrd},{i,nr}]; 
    Return[{wtot, xi, xis},Module];
];
   
cqmom21p[moms_, ks_, nr_, nrd_, nro_,wRos_,Ros_] :=Module[{mymom, mRds, xi, w, v, r1, ms, 
	condmoms,condmomvec,xis,ws,wtot}, 
	mymom=pointer[moms,ks,wRos,Ros]; 
	Do[
        mRds=Table[mymom[0,i,l-1],{i,0,2nrd-1}];
        {xi[l],w}=wheeler[mRds]; 
        v=Table[Table[xi[l][[i]]^j, {i,nrd}],{j,0,nrd-1}]; 
        r1=DiagonalMatrix[w];
        Do[ms[i]=Table[mymom[i,j,l-1],{j,0,nrd-1}],{i,0,2nr-1}]; 
        Do[condmoms[i]=LinearSolve[v.r1,ms[i]],{i,0,2 nr-1}]; 
        Do[
            condmomvec[j]=Table[condmoms[i][[j]],{i,0,2nr-1}];
            If[Norm[Im[condmomvec[j]]] > 0, Print["Imag cond momvec"];Abort[];];
        ,{j,nrd}]; 
        Do[{xis[l][i],ws[i]}=wheeler[condmomvec[i]],{i,nr}]; 
        Do[wtot[l][j,i]=w[[j]] ws[j][[i]],{j,nrd},{i,nr}];
    ,{l,nro}]; 
    Return[{wtot, xi, xis}, Module];
];

vander[x_,q_]:=Module[{n,w,b,s,t,xx,c},
    n=Length[q];
    If[n==1,
        w={q[[1]]};
        Return[w,Module];
    ];
    c=Table[0,{i,n}];
    w=c;
    Do[
        xx=-x[[i]];
        Do[
            c[[j]]=c[[j]]+xx c[[j+1]];
        ,{j,n+1-i,n-1}];
        c[[n]]=c[[n]]+xx;
    ,{i,2,n}];
    Do[
        xx=x[[i]];
        t=1;
        b=1;
        s=q[[n]];
        Do[
            b=c[[k]]+xx b;
            s=s+q[[k-1]]b;
            t=xx t+b;
        ,{k,n,2,-1}];
        w[[i]]=s/t;
    ,{i,1,n}];
    Return[w,Module];
];

adaptwheeler[m_,n_,mins_] := Module[{nn=n,mom=m,rmin=mins[[1;;n]],
    eabs=10^-10,sig,a,b,Ja,w,x,xi,eval,evec,esys,myn,ind,nu,
    nout,cutoff,dab,mab,bmin,z,mindab,maxmab},
    (* rmin is minimum ratio wmin/wmax *)
    (* eabs is minimum distance between distinct abscissas *)
    cutoff=0.;
    If [mom[[1]]<0,
        Print["Negative number density"];
        Abort[];
    ];
    If[mom[[1]]==0,
        w=0.;
        x=0.;
        Return[{{x},{w},1},Module];
    ];
    If[
        nn==1||mom[[1]]<rmin[[1]],
        w=mom[[1]];
        x=mom[[2]]/mom[[1]];
        Return[{{x},{w},1},Module];
    ];

    (* Compute modified moments equal to moments *)
    nu=mom;
    (* Construct recurrence matrix *)
    ind=nn;
    a=Table[0.,{i,ind}];
    b=Table[0.,{i,ind}];
    sig=Table[0.,{i,2ind+1},{j,2ind+1}];
    Do[
        sig[[2,i]]=nu[[i-1]];
    ,{i,2,2ind+1}];
    a[[1]]=nu[[2]]/nu[[1]];
    b[[1]]=0;
    Do[
        Do[sig[[k,l]]=sig[[k-1,l+1]]-a[[k-2]]sig[[k-1,l]]-b[[k-2]]sig[[k-2,l]],{l,k,2ind-k+3}];
        a[[k-1]]=sig[[k,k+1]]/sig[[k,k]]-sig[[k-1,k]]/sig[[k-1,k-1]];
        b[[k-1]]=sig[[k,k]]/sig[[k-1,k-1]];
    ,{k,3,ind+1}];
    (* Determine maximum n using diag elements of sig *)
    myn=nn;
    Do[
        If[
            sig[[k,k]]<=cutoff,
            myn=k-2;
            If[
                myn==1,
                w=mom[[1]];
                x=mom[[2]]/mom[[1]];
                Return[{{x},{w},1},Module];
            ];
        ];
    ,{k,ind+1,3,-1}];

    (* Compute quadrature using maximum n *)
    a=Table[0,{i,myn}];
    b=Table[0,{i,myn}];
    w=Table[0,{i,myn}];
    x=Table[0,{i,myn}];
    sig=Table[0,{i,2 myn+1},{j,2myn+1}];
    Do[sig[[2,i]]=nu[[i-1]],{i,2,2 myn+1}];
    a[[1]]=nu[[2]]/nu[[1]];
    b[[1]]=0.;
    Do[
        Do[sig[[k,l]]=sig[[k-1,l+1]]-a[[k-2]] sig[[k-1,l]]-b[[k-2]] sig[[k-2,l]],{l,k,2 myn-k+3}];
        a[[k-1]]=sig[[k,k+1]]/sig[[k,k]]-sig[[k-1,k]]/sig[[k-1,k-1]];
        b[[k-1]]=sig[[k,k]]/sig[[k-1,k-1]];
    ,{k,3,myn+1}];
    (* Check if moments are not realizable (should never happen) *)
    bmin=Min[b];
    If[bmin<0,Print["Moments in Wheeler moments are not realizable!"];Abort[];];
    (* Setup Jacobi matrix for n-point quadrature, adapt n using rmax and eabs *)
    Do[
        If[nl==1,
        w=mom[[1]];
        x=mom[[2]]/mom[[1]];
        Return[{{x},{w},1},Module];
        ];
        z=Table[0,{i,1,nl},{j,1,nl}];
        Do[
            z[[i,i]]=a[[i]];
            z[[i,i+1]]=Sqrt[b[[i+1]]];
            z[[i+1,i]]=z[[i,i+1]];
        ,{i,1,nl-1}];
        z[[nl,nl]]=a[[nl]];
        (* Compute weights and abscissas *)
        esys=Eigensystem[z];
        eval=esys[[1]];evec=esys[[2]];
        w=Table[0,{i,1,nl}];
        x=w;
        dab=w;
        mab=w;
        x=eval;
        w=Table[evec[[i,1]]^2 mom[[1]],{i,1,nl}];
        Do[
            dab[[i]]=Min[Abs[x[[i]]-x[[1;;i-1]]]];
            mab[[i]]=Max[Abs[x[[i]]-x[[1;;i-1]]]];
        ,{i,nl,2,-1}];
        mindab=Min[dab[[2;;nl]]];
        maxmab=Max[mab[[2;;nl]]];
        If[nl==2,maxmab=1];
        (* Check conditions that weights and abscissas must both satisfy *)
        If[
            Min[w]/Max[w]>rmin[[nl]]&&mindab/maxmab>eabs,
            Return[{x,w,Length[w]},Module];
        ];
    ,{nl,myn,1,-1}];
];

End[];
EndPackage[];

project[xix_, xiy_, w_, ks_, ksp_, nr_, nrd_] := Module[{moms, mom, momsp}, 
    mom[p_, q_] := Sum[w[[i]] xix[[i]]^p xiy[[i]]^q, {i, nr nrd}]; 
    moms = Table[mom[ks[[i, 1]], ks[[i, 2]]], {i, 1, Length[ks]}]; 
    momsp = Table[mom[ksp[[i, 1]], ksp[[i, 2]]], {i, 1, Length[ksp]}]; 
    Return[{moms,momsp},Module];
];
     
project1[xi_, xis_, wtot_, ks_, ksp_, nr_, nrd_,nro_:0,wRo_:0,Ros_:0] := 
  Module[{moms, mom, momsp},
	If[nro==0,
		mom[p_, q_] := Sum[wtot[j, i] xi[[j]]^p xis[j][[i]]^q, {j,nr}, {i,nrd}]; 
		moms = Table[mom[ks[[i, 1]], ks[[i, 2]]], {i, 1, Length[ks]}]; 
		momsp = Table[mom[ksp[[i, 1]], ksp[[i, 2]]], {i,Length[ksp]}];
	];
	If[nro>0,
		mom[p_, q_] := Sum[wtot[j, i] xi[[j]]^p xis[j][[i]]^q, {j,nr}, {i,nrd}]; 
		moms = Table[mom[ks[[i, 1]], ks[[i, 2]]], {i, 1, Length[ks]}]; 
		momsp = Table[mom[ksp[[i, 1]], ksp[[i, 2]]], {i,Length[ksp]}];	
	];
	Return[{moms, momsp},Module];
];
    
project2[xi_, xis_, wtot_, ks_, ksp_, nr_, nrd_] := 
  Module[{moms, mom, momsp}, 
   mom[p_, q_] := Sum[wtot[j, i] xi[[j]]^q xis[j][[i]]^p, {j,nrd}, {i,nr}]; 
   moms = Table[mom[ks[[i, 1]], ks[[i, 2]]], {i, 1, Length[ks]}]; 
   momsp = Table[mom[ksp[[i, 1]], ksp[[i, 2]]], {i, 1, Length[ksp]}]; 
   Return[{moms, momsp},Module];
];
