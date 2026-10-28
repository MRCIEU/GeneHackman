#!/usr/bin/env bash
# Start run_build.sh in a detached tmux session so it survives logging out (the full build takes days).
#
#   run_in_tmux.sh phase0          run_in_tmux.sh full EUR          run_in_tmux.sh full
#
# Window "build" runs the build; window "status" refreshes run_build.sh status every minute.
#   attach:  tmux attach -t ldpack        detach: Ctrl-b then d        switch window: Ctrl-b then n
# If the machine reboots or the run dies, start the same command again: finished work is skipped.
# Without tmux it falls back to nohup and prints where the log goes.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=config.sh
source "$HERE/config.sh"
SESSION="${TMUX_SESSION:-ldpack}"
[[ $# -gt 0 ]] || { sed -n '2,10p' "$0"; exit 1; }
mkdir -p "$WORK_DIR/logs"
args="$(printf '%q ' "$@")"

if ! command -v tmux >/dev/null 2>&1; then
  log_file="$WORK_DIR/logs/nohup_$(date +%Y%m%d_%H%M%S).log"
  nohup "$HERE/run_build.sh" "$@" > "$log_file" 2>&1 &
  echo "tmux not found: started with nohup (pid $!). Follow it with: tail -f $log_file"
  exit 0
fi
if tmux has-session -t "$SESSION" 2>/dev/null; then
  echo "A tmux session '$SESSION' already exists. Attach with: tmux attach -t $SESSION"
  echo "(or choose another name: TMUX_SESSION=ldpack2 $0 $args)"
  exit 1
fi
tmux new-session -d -s "$SESSION" -n build \
  "bash -c '\"$HERE/run_build.sh\" $args; code=\$?; echo; echo \"run_build.sh $args finished with exit code \$code\"; exec bash'"
tmux new-window -t "$SESSION" -n status \
  "bash -c 'while true; do clear; \"$HERE/run_build.sh\" status; sleep 60; done'"
tmux select-window -t "$SESSION:build"
echo "Started 'run_build.sh $args' in tmux session '$SESSION'."
echo "Attach: tmux attach -t $SESSION   (detach: Ctrl-b d; next window: Ctrl-b n)"
