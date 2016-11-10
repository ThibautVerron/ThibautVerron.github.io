# File "functions.mpl";

with(LinearAlgebra):
with(FGb):
with(VectorCalculus):
with(combinat):

randomize():
randval := rand(100..10000):

Minors := proc(s,K)
description "it computes the minors of size s of a (d1 X d2) matrix K output: the list of (s X s) minors of K";
local d1,d2,MIN,j,l,listofminors:
    (d1,d2):=Dimension(K):
    if s > min(d1,d2) then
        printf("ERROR: INPUT s TOO LARGE\n");
        # print(s,Dimensions(K));
        return;
    fi;

    MIN:=Matrix(binomial(d1,s),binomial(d2,s)):
    for j from 1 to binomial(d1,s) do
        for l from 1 to binomial(d2,s) do
            MIN[j,l]:=Determinant(SubMatrix(K,(choose([seq(i,i=1..d1)],s))[j],
                                            (choose([seq(i,i=1..d2)],s))[l])):
        od:
    od:

    listofminors:=seq(seq(MIN[j,l],l=1..binomial(d2,s)),j=1..binomial(d1,s)):
    return(listofminors):
end:

dimension := proc(sys, vars := indets(sys))
description "Find the dimension of the variety defined by the system "
    "by cutting it with random hyperplanes";
    local dim, hyp, randgen,hh, v, gg,sysCut,j, firstdone;

    dim := nops(vars):

    hyp := []:
    randgen := rand(-500..500):
    for j in seq(i,i=1..dim) do
        hh := randgen():
        for v in vars do
            hh := hh + randgen()*v:
        od;
        hyp := [op(hyp),hh]:
    od:

    sysCut := [op(sys),op(hyp)]:

    gg := [1]:
    firstdone := false;
    while gg = [1] and dim >= 0 do
        if firstdone then
            dim := dim -1:
            sysCut := sysCut[1..-2]:
        else
            firstdone := true:
        fi:
        printf("Trying dimension %d\n",dim):
        # printf("%a %a\n",nops(sysCut),nops(vars)):
        gg := fgb_gbasis(sysCut,65521,vars,[],{"index"=2000000}):
    od:
    return dim:
end:

sqfr := proc(f)
    description "Square-free reduction of f";
    return mul(ff,ff in map(x -> x[1],factors(f)[2])):
end:

GB_interpolate := proc(sys,charac,vars1,vars2,varinter,deg)
description "Use interpolation on varinter up to degree deg to compute a "
    "Groebner basis of sys for the order eliminating vars1";
local points, bases, v, i, gb, res, t0, safety, tgtdeg, m, bases_i, pol, frct;
    v := -1:
    points := []:
    bases := []:
    safety := 1:
    tgtdeg := 2*deg + 2 + 2*safety:
    for i from 1 to tgtdeg do
        printf("GB %a/%a: ", i,tgtdeg):
        t0 := time():
        while v = -1 or v in points do
            v := randval():
        od:
        gb := fgb_gbasis_elim(eval(sys,varinter=v),charac,vars1,vars2,{#"verb"=3,
                                                                       "index"=20000000}):
        points := [op(points),v]:
        bases := [op(bases),gb]:
        printf("done (%a) [v=%a]\n", time()-t0,v):
    od:
    res := []:

    m := mul(varinter-v,v in points):

    for i from 1 to nops(bases[1]) do
        bases_i := map(x -> x[i]/lcoeff(x[i])
                       , bases):
        # lprint(points,bases_i):
        pol := CurveFitting[PolynomialInterpolation](points,bases_i,varinter):
        frct := ratrecon(pol,m,varinter,deg+safety,deg+safety):
        res := [op(res),numer(frct)]:
    od:
    return res:
end:

GB_find_deg := proc(sys,vars1,vars2,varinter,ntrials := 10)
description "Find the maximal degree in variable varinter of the polynomials in "
    "a Groebner basis of sys eliminating vars1";
local deg, gb, ss, i, v, t0;
    deg := -1:
    printf("Finding degree (%a tests): ", ntrials):
    i := 0:
    facts := []:
    while i < ntrials do
        ss := {}:
        for v in vars2 do
            if v <> varinter then
                ss := {op(ss), v=randval()}:
            fi:
        od:
        gb := fgb_gbasis_elim(eval(sys,ss),0,vars1,vars2):
        #gb := sort(gb):
        for j from 1 to nops(gb) do
            if i = 0 then
                facts := [op(facts),0]:
            fi:
            #lprint(nops(facts),j):
            facts[j] := gcd(facts[j],gb[j]);
        od:
        # if nops(gb) >= idx then
        deg := max(map(degree,gb,varinter),deg):
        printf("|"):
        i := i+1:
        # else:
        #     printf("x"):
        # fi:
    od:
    printf("\n"):
    facts_prod := sqfr(mul(p,p in facts));
    return deg, facts_prod:
end:


Elimination_codim1 := proc(F,vars,params,{algo:=default})
description "Compute a polynomial whose zeroes cover the projection of V(F) on the parameter space"
    ""
    "Compute a system of generators for the elimination ideal "
    "of F obtained by eliminated the variables vars"
    "Uses the algorithm `algo` to compute the basis. Admissible values "
    "for `algo` are:"
    "- gb_direct : compute an elimination Gröbner basis directly "
    "- gb_interp : compute an elimination Gröbner basis using evaluation/interpolation"
    "Any other value defaults to gb_interp."
    ;
local var,deg_interp,G:
    # print(algo);
    if algo = gb_direct then
        ### Direct
        G := fgb_gbasis_elim(F,0,vars,params ,{"verb"=3}):
        return G[1];
    else
        ### Interpolation with first param
        var := params[1]:
        deg_interp, facts := GB_find_deg(F,vars,params,var):
        G := GB_interpolate(F,0,vars,params,var,deg_interp):
        # !!! Peut-être incorrect
        return facts*G[1];
    fi:
end:

DiscriminatingPolynomial := proc(M,r0,H,vars,params,eqs := [], neqs := [], {algo:=default})
description "Compute a polynomial P as in section 3."
    ""
    "The argument `eqs` is a set of extra equations and inequations restricting "
    "the solutions: the solutions returned are projections in the params "
    "of points at which all polynomials of eqs vanish."
    ""
    "The polynomials in `neqs` are saturated in the computations."
    ""
    "See Elimination for the parameter `algo`":
local res,res1,res2,res3,t0,t1,t2,t3:
    # print(algo);
    t0 := time():
    res1 := RankExactly(M,r0,vars,params,eqs,neqs,':-algo'=algo):
    t1 := time():
    res2 := DeterminantCritVals(M,r0,vars,params,eqs,neqs,true,':-algo'=algo):
    t2 := time():
    res3 := DeterminantBoundary(M,r0,H,vars,params,eqs,neqs,true,':-algo'=algo):
    t3 := time():
    printf("RankExactly\t%as\n",t1-t0):
    printf("DeterminantCritVals\t%as\n",t2-t1):
    printf("DeterminantBoundary\t%as\n",t3-t2):
    return res1*res2*res3, res1, res2, res3;
end:

# DiscPoly_factors := proc(M,r0,H,vars,params,eqs := [], neqs := [], {algo:=default})
# description "Compute a polynomial P as in section 3."
#     ""
#     "The argument `eqs` is a set of extra equations and inequations restricting "
#     "the solutions: the solutions returned are projections in the params "
#     "of points at which all polynomials of eqs vanish."
#     ""
#     "The polynomials in `neqs` are saturated in the computations."
#     ""
#     "See Elimination for the parameter `algo`":
# local res,res1,res2,res3,t0,t1,t2,t3:
#     # print(algo);
#     t0 := time():
#     res1 := RankExactly_factors(M,r0,vars,params,eqs,neqs,':-algo'=algo):
#     t1 := time():
#     res2 := DeterminantCritVals_factors(M,r0,vars,params,eqs,neqs,true,':-algo'=algo):
#     t2 := time():
#     res3 := DeterminantBoundary_factors(M,r0,H,vars,params,eqs,neqs,true,':-algo'=algo):
#     t3 := time():
#     printf("RankExactly\t%as\n",t1-t0):
#     printf("DeterminantCritVals\t%as\n",t2-t1):
#     printf("DeterminantBoundary\t%as\n",t3-t2):
#     return res1*res2*res3, res1, res2, res3;
# end:


IncidenceVariety := proc(M,r)
description "Compute a system of generators for the incidence variety "
    "of rank r of M"
    ;
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

RankExactly := proc(M,r0,vars,params,eqs := [], neqs := [], {algo := default})
description "Algorithm RankExactly (section 3.4)"
    ""
    "See DiscriminatingPolynomial for a description of parameters `eqs` and "
    "`neqs`."
    "See Elimination for a description of the parameter `algo`";
local res,k,n,FVr,dim,JVr,FV,F0,varsY,Mins,i,Sysi,G, codim:
    # print(algo);
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
        g := Elimination_codim1(Sysi,[uu,u,op(varsY),op(vars)],params,':-algo'=algo):
        if g <> FAIL then
            res := res*g:
        fi:
    od:
    printf("RankExactly: %a\n", factor(res));
    return res:
end:

DeterminantCritVals := proc(M,r0,vars,params,eqs := [], neqs := [], skipRankExactly := false, {algo:=default})
    description "Algorithm DeterminantCritVals (section 3.5)"
    ""
    "See DiscriminatingPolynomial for a description of parameters `eqs` and "
    "`neqs`."
    "See Elimination for a description of the parameter `algo`":
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
    g := Elimination_codim1(F1,[u,op(varsY),op(vars)],params,':-algo'=algo):
    res := res*g:
    printf("CritVals: %a\n", factor(res));
    return res:
end:

DeterminantBoundary := proc(M,r0,H,vars,params,
                            eqs := [], neqs := [], skipRankExactly := false,
                            {algo := default})
    description "Algorithm DeterminantBoundary (section 3.6)"
    ""
    "See DiscriminatingPolynomial for a description of parameters `eqs` and "
    "`neqs`."
    "See Elimination for a description of the parameter `algo`";
local res, F0, varsY, h, F1, G:
    if not skipRankExactly then
        res := RankExactly(M,r0,vars,params):
    else
        res := 1:
    fi:
    F0, varsY := IncidenceVariety(M,r0-1):
    for h in H do
        F1 := [op(F0),h,op(eqs),u*mul(f, f in neqs)-1]:
        g := Elimination_codim1(F1,[u,op(varsY),op(vars)],params,':-algo'=algo):
        res := res * g:
    od:
    printf("Boundary: %a\n", factor(res));
    return res:
end:


