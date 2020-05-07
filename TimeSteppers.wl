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

RK34[mom_,myrhs_,t_,dt_]:=Module[{moms=mom,mome,momstemp1,momstemp2,rhs,error,c},
    (* coefficients for SSP-RK4 *)

    Print["Not finished"];
    Abort[];
    (* TODO: Add time dependency into myrhs calls *)
    c[1]=37/378;c[2]=0;c[3]=250/621;c[4]=125/594;c[5]=0;c[6]=512/1771;
    (* SSP-RK3 *)
	{moms,rhs}=myrhs[moms,t];
    momstemp1=moms+dt rhs;
    {momstemp1,rhs}=myrhs[momstemp1,t];
    momstemp2=(3/4)moms+(1/4)(momstemp1+dt rhs);
    {momtemp2,rhs}=myrhs[momstemp2,t];
    mome=(1/3)moms+(2/3)(momtemp2+dt rhs);
    {mome,rhs}=myrhs[mome];

    (* SSP-RK4 *)
    {moms,rhs}=myrhs[moms];
    momstemp1=c[1,0]moms+c[1,1]dt rhs;
    {momstemp1,rhs}=myrhs[momstemp1];
    momstemp2=c[2,0]moms+c[2,1] momstemp1+c[2,2] dt rhs;
    {momstemp2,rhs}=myrhs[momstemp2];
    momstemp3=c[3,0]moms+c[3,1]momstemp2+c[3,2]dt rhs;
    {momstemp3,rhs3}=myrhs[momstemp3];
    momstemp4=c[4,0]moms+c[4,1]momstemp3+c[4,2] dt rhs3;
    {momstemp4,rhs4}=myrhs[momstemp4];
    moms=c[5,0]momstemp2+c[5,1]momstemp3+c[5,2]dt rhs3+c[5,3]momstemp4+c[5,4]dt rhs4;
    {moms,rhs}=myrhs[moms];*)
	Return[{moms,mome,err[moms,mome]},Module];
];

End[];
EndPackage[];


(*
cs[1]=2825/27648;cs[2]=0;cs[3]=18575/48384;cs[4]=13525/55296;cs[5]=277/14336;cs[6]=1/4;

b0[2,1]=1/5;
b0[3,1]=3/40;b0[3,2]=9/40;
b0[4,1]=3/10;b0[4,2]=-(9/10);b0[4,3]=6/5;
b0[5,1]=-(11/54);b0[5,2]=5/2;b0[5,3]=-(70/27);b0[5,4]=35/27;
b0[6,1]=1631/55296;b0[6,2]=175/512;b0[6,3]=575/13824;b0[6,4]=44275/110592;b0[6,5]=253/4096;
*)

(* SSP RK(5,4) coefficients *)
(*c[1,0]=1;c[1,1]=0.391752226571890;
c[2,0]=0.444370493651235;c[2,1]=0.555629506348765;c[2,2]=0.368410593050371;
c[3,0]=0.620101851488403;c[3,1]=0.379898148511597;c[3,2]=0.251891774271694;
c[4,0]=0.178079954393132;c[4,1]=0.821920045606868;c[4,2]=0.544974750228521;
c[5,0]=0.517231671970585;c[5,1]=0.096059710526147;c[5,2]=0.063692468666290;c[5,3]=0.386708617503269;
c[5,4]=0.226007483236906;*)

(* Embedded SSPRK32 from Gottleib, Ketchesen, Chu [doesn't work at all] *)
(*{moms,rhs}=myrhs[moms];
momstemp1=moms+(dt/2) rhs;
{momstemp1,rhs}=myrhs[momstemp1];
momstemp2=momstemp1+(dt/2)rhs;
{momstemp2,rhs}=myrhs[momstemp2];
mome=(1/3) moms+(2/3)(momstemp2+dt rhs);
momstemp3=(2/3) moms+(1/3)(momstemp2+dt rhs);
{momtemp3,rhs}=myrhs[momstemp3];
moms=momstemp3+(dt/2)rhs;*)

(* Four-stage SSP-RK3 *)
(*moms=moms0;
{moms,rhs}=myrhs[moms];
momstemp1=(1/2)moms+(1/2)(moms+dt rhs);
{momtemp1,rhs}=myrhs[momstemp1];
momstemp2=(1/2)momstemp1+(1/2)(momstemp1+dt rhs);
{momtemp2,rhs}=myrhs[momstemp2];
momstemp3=(2/3)moms+(1/6)momstemp2+(1/6)(momstemp2+dt rhs);
{momstemp3,rhs}=myrhs[momstemp3];
moms=(1/2)momstemp3+(1/2)(momstemp3+dt rhs);*)

(* Fifth-order RK *)
(*Do[b[i,j]=dt b0[i,j],{i,1,6},{j,1,6}];
{mtemp[1],myk[1]}=myrhs[moms];
{mtemp[2],myk[2]}=myrhs[mtemp[1]+b[2,1] myk[1]];
{mtemp[3],myk[3]}=myrhs[mtemp[2]+b[3,1] myk[1]+b[3,2] myk[2]];
{mtemp[4],myk[4]}=myrhs[mtemp[3]+b[4,1] myk[1]+b[4,2] myk[2]+b[4,3] myk[3]];
{mtemp[5],myk[5]}=myrhs[mtemp[4]+b[5,1] myk[1]+b[5,2] myk[2]+b[5,3] myk[3]+b[5,4] myk[4]];
{mtemp[6],myk[6]}=myrhs[mtemp[5]+b[6,1] myk[1]+b[6,2] myk[2]+b[6,2] myk[3]+b[6,4] myk[4]+b[6,5] myk[5]];

moms=mtemp[6]+dt Sum[c[i] myk[i],{i,1,6}];


(*(* Euler step *)
momstemp1=moms+dt rhs;
{momstemp1,rhs}=myrhs[momstemp1];
mome=momstemp1;
(* SSP-RK2 *)
moms=(1/2) moms+(1/2)(momstemp1+dt rhs);
{moms,rhs}=myrhs[moms];*)
