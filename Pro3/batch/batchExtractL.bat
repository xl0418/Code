@echo off
set ticks=1e+07
(for %%s in (HH) do (

(for /l %%j in (1,1,5) do (

(for /l %%i in (1,1,6) do (

(for /l %%r in (1,1,100) do (

jcpp file=sce%%s/results/%ticks%/spatialpara%ticks%%%s%%j%%i/%%spsi%%js_phi%%irep%%r.m output=sce%%s/results/%ticks%/spatialpara%ticks%%%s%%j%%i/%%spsi%%js_phi%%irep%%rLtable 
echo sce %%s and psi %%j and phi %%i and replicate %%r
))
))
))
))