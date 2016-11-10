# File "water_functions.mpl":

# Function definitions

read "functions.mpl";


# Definitions for water

read "water_def.mpl";


# Computations

t0 := time():
res_general,res1,res2,res3 := DiscriminatingPolynomial(MatrixD,3,H,vars,params,[],[]# , algo = gb_direct
                                                       , algo=gb_interp
                                                      ):
t1 := time():


# Results

# RankExactly
ff1 := map(x -> x[1],factors(res1)[2]):
nops(ff1);
# 6
map(degree,ff1);
# [1,1,1,2,1,1]

# Critical values
ff2 := map(x -> x[1],factors(res2)[2]):
nops(ff2);
# 5
map(degree,ff2);
# [1, 1, 1, 1, 1]

# Boundary
ff3 := map(x -> x[1],factors(res3)[2]):
nops(ff3);
# 7
map(degree,ff3);
# [1, 1, 1, 1, 2, 2, 3]

allfacts := [op({op(ff1),op(ff2),op(ff3)})]:
nops(allfacts);
# 13
map(degree,allfacts);
# [1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 3]

filename := "water_results.mpl";
fid := fopen(filename,WRITE):
fprintf(fid,"## RankExactly:\n##---------\n"):
fprintf(fid,"res1 := [\n"):
for i from 1 to nops(ff1)-1 do
    f := ff1[i]:
    fprintf(fid,"    %a,\n", f):
od:
fprintf(fid,"    %a\n]:\n\n",ff1[-1]):
fprintf(fid,"## CritVals:\n##---------\n"):
fprintf(fid,"res2 := [\n"):
for i from 1 to nops(ff2)-1 do
    f := ff2[i]:
    fprintf(fid,"    %a,\n", f):
od:
fprintf(fid,"    %a\n]:\n\n",ff2[-1]):
fprintf(fid,"## Boundary\n##---------\n"):
fprintf(fid,"res3 := [\n"):
for i from 1 to nops(ff3)-1 do
    f := ff3[i]:
    fprintf(fid,"    %a,\n", f):
od:
fprintf(fid,"    %a\n]:\n\n",ff3[-1]):
fprintf(fid,"## All (without duplicates)\n##---------\n"):
fprintf(fid,"allfacts := [\n"):
for i from 1 to nops(allfacts)-1 do
    f := allfacts[i]:
    fprintf(fid,"    %a,\n", f):
od:
fprintf(fid,"    %a\n]:\n\n",allfacts[-1]):
fclose(fid):
