#!/usr/bin/env bash
# Onstart watchdog (PLAN §4.5, G5; unpacked by policy.onstart_cmd, size-limited: details in M7.md). Destroys (never
# stops, G4) via CONTAINER_API_KEY (V9; DELETE endpoint UNVERIFIED until V0) when the uptime since the first start
# (persisted in WATCHDOG_STATE, or the earlier WATCHDOG_CREATED; BK4-17) reaches WATCHDOG_CAP_S, or when the GPU
# idled WATCHDOG_IDLE_S and every output is on the Mac (EXIT, PULLED_OK after the idle start, nothing newer; BK4-11).
set -euo pipefail
CAP_S=${WATCHDOG_CAP_S:?}; IDLE_S=${WATCHDOG_IDLE_S:-1800}; TICK=${WATCHDOG_TICK_S:-60}
MAX_TICKS=${WATCHDOG_MAX_TICKS:-0}; PULLED=${WATCHDOG_PULLED_OK:?}; API=${VAST_API_URL:-https://console.vast.ai/api/v0}
RUN=$(dirname "$PULLED"); STATE=${WATCHDOG_STATE:-/root/watchdog.start}; ENV_FILE=${WATCHDOG_ENV_FILE:-/etc/environment}
# ssh sees only /etc/environment (V9): non-secret run variables, once
for v in RUN_ID CONTAINER_ID WATCHDOG_CAP_S WATCHDOG_PULLED_OK; do
    if [ -n "${!v:-}" ] && ! grep -qxF "$v=${!v}" "$ENV_FILE" 2>/dev/null; then echo "$v=${!v}" >> "$ENV_FILE"; fi
done
destroy() {
    echo "watchdog $(date -u +%FT%TZ): destroying instance ${CONTAINER_ID}: $1"
    until curl -fsS -X DELETE -H "Authorization: Bearer ${CONTAINER_API_KEY}" "$API/instances/${CONTAINER_ID}/"; do
        sleep 10
    done
    exit 0
}
[ -f "$STATE" ] || date +%s > "$STATE"
start=$(cat "$STATE"); created=${WATCHDOG_CREATED:-$start}
if (( created < start )); then start=$created; fi
idle_since=; n=0
while :; do
    now=$(date +%s)
    if (( now - start >= CAP_S )); then destroy "cap ${CAP_S} s reached"; fi
    util=$(nvidia-smi --query-gpu=utilization.gpu --format=csv,noheader,nounits 2>/dev/null | sort -n | tail -1) || util=
    [[ $util =~ ^[0-9]+$ ]] || util=100  # no reading counts as busy; the cap is the backstop
    if (( util < 5 )); then idle_since=${idle_since:-$now}; else idle_since=; fi
    # on the Mac: EXIT and PULLED_OK exist, nothing is newer than PULLED_OK, and the pull came after the idle start
    if [ -n "$idle_since" ] && [ -f "$PULLED" ] && [ -f "$RUN/EXIT" ] \
        && [ -z "$(find "$RUN" -newer "$PULLED" -print -quit)" ] \
        && (( $(date -r "$PULLED" +%s) >= idle_since && now - idle_since >= IDLE_S )); then
        destroy "GPU idle $(( now - idle_since )) s, job exited and pulled after the idle start"
    fi
    n=$((n + 1))
    if (( MAX_TICKS && n >= MAX_TICKS )); then exit 3; fi
    sleep "$TICK"
done
