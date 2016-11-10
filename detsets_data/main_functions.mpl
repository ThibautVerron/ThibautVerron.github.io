# read "main_functions.mpl";

read "functions.mpl";

DiscriminatingPolynomial := proc(M,r0,H,vars,params,eqs := [], neqs := [], elim:="")
    description "Compute a polynomial P as in section 3";
local res:
    t0 := time():
    res1 := RankExactly(M,r0,vars,params,eqs,neqs):
    t1 := time():
    res2 := DeterminantCritVals(M,r0,vars,params,eqs,neqs,true):
    t2 := time():
    res3 := DeterminantBoundary(M,r0,H,vars,params,eqs,neqs,true):
    t3 := time():
    printf("RankExactly\t%as\n",t1-t0):
    printf("DeterminantCritVals\t%as\n",t2-t1):
    printf("DeterminantBoundary\t%as\n",t3-t2):
    return res1*res2*res3, res1, res2, res3;
end:

IncidenceVariety := proc(M,r)
description "Compute a system of generators for the incidence variety "
    "of rank r of M";
local MatrixU, MatrixY,k,i,j,Prod1,Sys1,Prod2,Sys2,Sys:
    k := RowDimension(M):
    MatrixU := Matrix(k-r,k):
    # print(MatrixU);
    for i from 1 to k-r do
        for j from 1 to k do
            MatrixU[i,j] := randval():
            # print(i,j,MatrixU);
        od:
    od:
    # MatrixU1 := Matrix(k-r,k-r,shape=identity):
    # MatrixU2 := Matrix(k-r,r,shape=zero):
    # MatrixU := Matrix([MatrixU1,MatrixU2]):

    print(MatrixU);

    MatrixY := Matrix(k,k-r):
    for i from 1 to k do
        for j from 1 to k-r do
            MatrixY[i,j] := Y[i,j]:
        od:
    od:
    Prod1 := M . MatrixY:
    # print(Dimensions(Prod1));
    Sys1 := [seq(seq(Prod1[i][j],j=1..k-r),i=1..k)]:
    Prod2 := MatrixU . MatrixY - Matrix(k-r,k-r,shape=identity):
    # print(Dimensions(Prod2));
    Sys2 := [seq(seq(Prod2[i][j],j=1..k-r),i=1..k-r)]:
    Sys := [op(Sys1),op(Sys2)]:
    return Sys, indets(MatrixY):
end:

RankExactly := proc(M,r0,vars,params,eqs := [], neqs := [], {elim := ""})
    description "Algorithm RankExactly (section 3.4)";
local res,k,n,FVr,dim,JVr,FV,F0,varsY,Mins,i,Sysi,G, codim:
    k := RowDimension(M):
    n := nops(vars):
    res := 1:

    FVr := [Minors(r0+1,M)]:
    codim := (k-r0)^2:
    JVr := Jacobian(FVr,[op(vars)]):
    FV := [op(FVr), Minors(codim,JVr),op(eqs),u*mul(f, f in neqs)-1]:
    F0, varsY := IncidenceVariety(M,r0):

    Mins := [Minors(r0,M)]:
    for i from 1 to nops(Mins) do
        #i := 1:
        Sysi := [op(FV),op(F0),op(Mins[1..i-1]),uu*Mins[i]-1]:
        G := Elimination(Sysi,[uu,u,op(varsY),op(vars)],params):
        if res1 <> FAIL then
            res := res*G[1]:
        fi:
    od:
    printf("RankExactly: %a\n", factor(res));
    return res:
end:

DeterminantCritVals := proc(M,r0,vars,params,eqs := [], neqs := [], skipRankExactly := false, {elim:=""})
    description "Algorithm DeterminantCritVals (section 3.5)";
local res,F0,varsY,N,J,F1,G:
    if not skipRankExactly then
        res := RankExactly(M,r0,vars,params):
    else
        res := 1:
    fi:
    #res := 1;
    F0, varsY := IncidenceVariety(M,r0-1):
    N := nops(F0):
    J := Jacobian(F0,[op(varsY),op(vars)]):
    F1 := [op(F0),Minors(N,J),op(eqs),u*mul(f, f in neqs)-1]:
    G := Elimination(F1,[u,op(varsY),op(vars)],params):
    res := res*G[1]:
    printf("CritVals: %a\n", factor(res));
    return res:
end:

DeterminantBoundary := proc(M,r0,H,vars,params,
                            eqs := [], neqs := [], skipRankExactly := false,
                            {elim := ""})
    description "Algorithm DeterminantBoundary (section 3.6)";
local res, F0, varsY, h, F1, G:
    if not skipRankExactly then
        res := RankExactly(M,r0,vars,params):
    else
        res := 1:
    fi:
    F0, varsY := IncidenceVariety(M,r0-1):
    for h in H do
        F1 := [op(F0),h,op(eqs),u*mul(f, f in neqs)-1]:
        G := Elimination(F1,[u,op(varsY),op(vars)],params,algo=""):
        res := res * G[1]:
    od:
    printf("Boundary: %a\n", factor(res));
    return res:
end:
