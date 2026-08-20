#!/usr/bin/env bash
# First MC sweep: disk n=5..10 (open values), plus square/triangle n=7..10 as extra
# Valtr anchors for the fast numba tester. nice'd; keeps clear of production load.
set -u
cd "$(dirname "$0")"
PY=../.venv/bin/python
OUT=../results
mkdir -p "$OUT"
log(){ echo "$(date '+%F %T') $*" | tee -a "$OUT/sweep.log"; }
log "sweep start (host $(hostname), nproc $(nproc))"
# cross-check both testers on a modest sample first (must agree exactly)
for n in 5 6 7 8; do
  nice -n 10 $PY convex_position.py --body disk --n $n --samples 2e6 --seed $((500+n)) --out "$OUT/disk_n${n}_xcheck_2e6.json" > "$OUT/disk_n${n}_xcheck_2e6.log" 2>&1
  log "xcheck disk n=$n: $(python3 -c "import json;d=json.load(open('$OUT/disk_n${n}_xcheck_2e6.json'));print(d['p_hat'],d.get('testers_agree'))")"
done
for n in 5 6 7 8 9 10; do
  nice -n 10 $PY convex_position.py --body disk --n $n --samples 2e8 --no-both --seed $((100+n)) --out "$OUT/disk_n${n}_2e8.json" > "$OUT/disk_n${n}_2e8.log" 2>&1
  log "disk n=$n: $(python3 -c "import json;d=json.load(open('$OUT/disk_n${n}_2e8.json'));print(d['p_hat'],'+-',d['std_err'],'(',round(d['seconds']),'s)')")"
done
for body in square triangle; do for n in 7 8 9 10; do
  nice -n 10 $PY convex_position.py --body $body --n $n --samples 1e8 --no-both --seed $((200+n)) --out "$OUT/${body}_n${n}_1e8.json" > "$OUT/${body}_n${n}_1e8.log" 2>&1
  log "$body n=$n: $(python3 -c "import json;d=json.load(open('$OUT/${body}_n${n}_1e8.json'));print(d['p_hat'],'+-',d['std_err'])")"
done; done
log "sweep done"
