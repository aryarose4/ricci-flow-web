#!/bin/sh
# serial final matrix (clean wall times); run from the worktree root
E=./bin/crumple7.exe; R=artifacts/anti-crumple/final
B8='0.8,0.8,0.8'; P='0.06,0;-0.06,0;0,-0.7'
X1='0.25,0.4,0.5'; P1='0.7,0;-0.35,0.6;-0.35,-0.6'
X6='0.2,0.35,0.5'; P6='0.2,0;-0.2,0;0,-0.7'
B3='0.6,0.7,0.5,0.8'; P3='0.7,0.7;-0.7,0.7;-0.7,-0.7;0.7,-0.7'
N4='0.3,0.3,0.3,0.3'; SQ='0.7,0;0,0.7;-0.7,0;0,-0.7'
run() { name=$1; mode=$2; shift 2; $E $R/${name}_$mode.json $mode "$@" | sed "s/^/$name /"; }
for m in base A AB BC ABC; do run c2a $m gluing coarse $X1 "$P1"; done
for m in base A AB BC ABC; do run c2b $m gluing coarse $X6 "$P6"; done
for m in base A AB BC ABC; do run c3a $m polytope coarse $B3 "$P3"; done
for m in base A AB BC ABC; do run c3b $m polytope coarse $B8 "$P"; done
for m in base A B AB; do run c4 $m gluing coarse $N4 "$SQ"; done
for m in base BC ABC A; do run c1 $m gluing medium $B8 "$P"; done
