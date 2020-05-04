(* ::Package:: *)

BeginPackage["TimeSteppers`"];

RK23::usage="";

Begin["`Private`"];

err[fine_,coarse_]:=Norm[fine-coarse]/Norm[fine];

RK23[mom_,myrhs_,t_,dt_]:=Module[{moms=mom,mome,momstemp1,momstemp2,rhs,error},
	(* SSP-RK2 *)
	{moms,rhs}=myrhs[moms,t];
	momstemp1=moms+dt rhs;
	{momstemp1,rhs}=myrhs[momstemp1,t+dt];
	mome=(1/2) moms+(1/2)(momstemp1+dt rhs);
	(* SSP-RK3 *)
	momstemp2=(3/4)moms+(1/4)(momstemp1+dt rhs);
	{momstemp2,rhs}=myrhs[momstemp2,t+dt/2];
	moms=(1/3)moms+(2/3)(momstemp2+dt rhs);
	Return[{moms,mome,err[moms,mome]},Module];
];


End[];
EndPackage[];
