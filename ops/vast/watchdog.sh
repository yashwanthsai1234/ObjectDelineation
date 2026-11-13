#!/usr/bin/env bash
# Onstart watchdog (unpacked by remote.onstart_cmd, which limits its size). Destroys (never stops) via
# CONTAINER_API_KEY when the uptime since the first start (persisted in WATCHDOG_STATE, or the earlier
# WATCHDOG_CREATED) reaches WATCHDOG_CAP_S, or when the GPU idled WATCHDOG_IDLE_S and every output is on the Mac (EXIT,
# PULLED_OK after the idle start, no file newer outside WATCHDOG_SKIP, remote.PULL_EXCLUDE), or idled
# 2 x WATCHDOG_IDLE_S and no job ever started (no HEARTBEAT).
set -euo pipefail
CAP_S=${WATCHDOG_CAP_S:?}; IDLE_S=${WATCHDOG_IDLE_S:-1800}; TICK=${WATCHDOG_TICK_S:-60}
MAX_TICKS=${WATCHDOG_MAX_TICKS:-0}; PULLED=${WATCHDOG_PULLED_OK:?}; API=${VAST_API_URL:-https://console.vast.ai/api/v0}
RUN=$(dirname "$PULLED"); SKIP=${WATCHDOG_SKIP:-}; STATE=${WATCHDOG_STATE:-/root/watchdog.start}; ENV_FILE=${WATCHDOG_ENV_FILE:-/etc/environment}
# ssh sees only /etc/environment: non-secret run variables, once
for v in RUN_ID CONTAINER_ID WATCHDOG_CAP_S WATCHDOG_PULLED_OK; do
    if [ -n "${!v:-}" ] && ! grep -qxF "$v=${!v}" "$ENV_FILE" 2>/dev/null; then echo "$v=${!v}" >> "$ENV_FILE"; fi
done
destroy() {
    echo "watchdog $(date -u +%FT%TZ): destroying instance ${CONTAINER_ID}: $1"
    until curl -fsS --max-time 30 -X DELETE -H "Authorization: Bearer ${CONTAINER_API_KEY}" "$API/instances/${CONTAINER_ID}/"; do
        sleep 10
    done
    exit 0
}
[ -f "$STATE" ] || date +%s > "$STATE"
start=$(cat "$STATE"); created=${WATCHDOG_CREATED:-$start}
if (( created < start )); then start=$created; fi
prune=-false; for d in ${SKIP//,/ }; do prune="$prune -o -name ${d%/}"; done
idle_since=; n=0
while :; do
    now=$(date +%s)
    if (( now - start >= CAP_S )); then destroy "cap ${CAP_S} s reached"; fi
    util=$(nvidia-smi --query-gpu=utilization.gpu --format=csv,noheader,nounits 2>/dev/null | sort -n | tail -1) || util=
    [[ $util =~ ^[0-9]+$ ]] || util=100  # no reading counts as busy; the cap is the backstop
    if (( util < 5 )); then idle_since=${idle_since:-$now}; else idle_since=; fi
    if [ -n "$idle_since" ] && [ ! -f "$RUN/HEARTBEAT" ] && (( now - idle_since >= 2 * IDLE_S )); then
        destroy "GPU idle $(( now - idle_since )) s and no job ever started"
    fi
    # on the Mac: EXIT and PULLED_OK exist, no file is newer than PULLED_OK, and the pull came after the idle start
    if [ -n "$idle_since" ] && [ -f "$PULLED" ] && [ -f "$RUN/EXIT" ] \
        && [ -z "$(find "$RUN" \( $prune \) -prune -o -type f -newer "$PULLED" -print -quit)" ] \
        && (( $(date -r "$PULLED" +%s) >= idle_since && now - idle_since >= IDLE_S )); then
        destroy "GPU idle $(( now - idle_since )) s, job exited and pulled after the idle start"
    fi
    n=$((n + 1))
    if (( MAX_TICKS && n >= MAX_TICKS )); then exit 3; fi
    sleep "$TICK"
done
