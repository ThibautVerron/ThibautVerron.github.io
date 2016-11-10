# File "water_computations.mpl";

read "functions.mpl";


# System

read "water_def.mpl":
DetD:=Determinant(MatrixD):
gradD_yz := [diff(DetD, y1), diff(DetD, z1), diff(DetD, y2), diff(DetD, z2)]:

r0 := 3:

Sys := [DetD,op(gradD_yz)]:


# Incidence varieties

F2,varsY2 := IncidenceVariety(MatrixD,2):
F3,varsY3 := IncidenceVariety(MatrixD,3):


# Singular points and critical values
# Rank(M) <= 2

Jac := Jacobian(F2,[z1,z2,y1,y2,op(varsY2)]):
Sys_rk2 := [op(F2),Minors(2*k+(k-r0+1)^2,Jac),u*g2-1]:

GB_rk2 := fgb_gbasis_elim(Sys_rk2,0,[u,op(varsY2),z1,z2,y1,y2],[G2,g2]):
lprint(factor(GB_rk2));

(* Results:
[(G2-1)^2*(-2*g2+1+G2)^2*(-g2-1+2*G2)^2*(2*G2-g2)^2]
*)


# Intersection with the boundary

# Side 1
Sys_bnd1 := [op(F2),H[1]]:
GB_bnd1 := fgb_gbasis_elim(Sys_bnd1,0,[op(varsY2),z1,z2,y1,y2],[G2,g2]):
lprint(factor(GB_bnd1));

(* Results:
[g2*(G2-1)^2*(-2*g2-1+3*G2)*(3*G2^2-5*G2*g2+g2^2+2*G2-2*g2+1)]
*)

# Side 2

Sys_bnd2 := [op(F2),H[2]]:
GB_bnd2 := fgb_gbasis_elim(Sys_bnd2,0,[op(varsY2),z1,z2,y1,y2],[G2,g2]):
lprint(factor(GB_bnd2));
(* Results:
[(G2-1)^2*(2*G2^2-5*G2*g2+2*g2^2-2*G2+3*g2)*(2*G2^3-3*G2^2*g2-3*G2*g2^2+2*g2^3+2*G2^2+9*G2*g2-11*g2^2-4*G2+6*g2)*(-g2+2*G2)^2]
*)



# Rank(M) = 3

allGB_rk3 := []:
MM := [Minors(k-1,MatrixD)]:

Sys_rk3 := [DetD,op(gradD_yz),op(F3),u1*g2-1]:
for i from 1 to k^2 do
    Sys_rk3_i := [op(Sys_rk3),op(MM[1..i-1]),u2*MM[i]-1]:
    GB := fgb_gbasis_elim(Sys_rk3_i,0,[u1,u2,op(varsY3),z1,z2,y1,y2],[G2,g2],{"verb"=3}):
    allGB_rk3 := [op(allGB_rk3),[i,GB]]:
od:

for g in allGB_rk3 do
    if g[2] <> [1] then
        printf("i%a Basis:%a\n",g[1],factor(g[2])):
    fi:
od:
(* Results
i=1 Basis:[(-g2+2*G2)*(g2-2+G2)*(2*G2^2-5*G2*g2+2*g2^2+1)]
*)


