#!/usr/bin/env bash
# Opens the EVE curation page with a double click: picks the results folder, starts
# server.py on a free port, opens the browser and stops the server when the small
# "server running" window is closed.
#
#   ./open_curation.sh                      # asks for the results folder
#   ./open_curation.sh Species_phlebotominae
#
# EVE_PYTHON chooses another interpreter (default: the active conda env, or eve_env).

AQUI="$(cd "$(dirname "$(readlink -f "$0")")" && pwd)"
CONFIG="$HOME/.config/eve_tool"
mkdir -p "$CONFIG"

aviso() {   # message in a window when there is a desktop, otherwise in the terminal
    if command -v zenity >/dev/null; then zenity --error --title="EVE curation" --width=420 --text="$1"
    else echo "$1" >&2; fi
}

# Python with the pipeline tools (samtools, diamond, bedtools, blast) on the PATH
# Order: EVE_PYTHON, the active conda env, an env called eve_env, python3 on the PATH
PY="${EVE_PYTHON:-}"
for c in "$PY" "${CONDA_PREFIX:+$CONDA_PREFIX/bin/python3}" \
         "$(conda info --base 2>/dev/null)/envs/eve_env/bin/python3" \
         "$HOME/miniconda3/envs/eve_env/bin/python3" "$HOME/miniforge3/envs/eve_env/bin/python3" \
         "$(command -v python3)"; do
    [ -n "$c" ] && [ -x "$c" ] && PY="$c" && break
done
export PATH="$(dirname "$PY"):$PATH"
export PYTHONPATH="$AQUI${PYTHONPATH:+:$PYTHONPATH}"   # the evetool package

# Results folder: argument, otherwise a folder dialog that starts at the last one used
RES="${1:-}"
if [ -z "$RES" ] && ! command -v zenity >/dev/null; then
    echo "Usage: $0 <results folder>   (zenity is not installed, so there is no folder dialog)" >&2
    exit 1
fi
if [ -z "$RES" ]; then
    ULTIMO="$(cat "$CONFIG/ultima_pasta" 2>/dev/null)"
    [ -d "$ULTIMO" ] || ULTIMO="$(eval "$("$PY" -m evetool.config --shell)"; echo "$EVE_RESULTS_DIR")"
    [ -d "$ULTIMO" ] || ULTIMO="$AQUI"
    RES="$(zenity --file-selection --directory --title="EVE curation - choose the results folder" \
           --filename="$ULTIMO/" 2>/dev/null)" || exit 0
fi
RES="$(readlink -f "$RES")"
# Accept the 12_dashboard_web folder itself as well
[ "$(basename "$RES")" = "12_dashboard_web" ] && RES="$(dirname "$RES")"
if [ ! -f "$RES/eves_classified.tsv" ]; then
    aviso "No eves_classified.tsv in:\n$RES\n\nChoose the results folder of a pipeline run."
    exit 1
fi
echo "$RES" > "$CONFIG/ultima_pasta"

# Free port, so two result folders can be open at the same time
PORTA="$("$PY" -c 'import socket; s=socket.socket(); s.bind(("127.0.0.1",0)); print(s.getsockname()[1])')"
LOG="$RES/13_curadoria/server.log"
mkdir -p "$RES/13_curadoria"

"$PY" -m evetool.curation.server --dir "$RES/12_dashboard_web" --outdir "$RES" --port "$PORTA" \
    --threads "$(nproc 2>/dev/null || echo 4)" > "$LOG" 2>&1 &
PID=$!
trap 'kill $PID 2>/dev/null' EXIT

# Wait for the server to answer (it writes the snapshot before listening)
for _ in $(seq 1 60); do
    if ! kill -0 $PID 2>/dev/null; then
        aviso "The server stopped while starting.\n\nLog: $LOG\n\n$(tail -n 8 "$LOG")"
        exit 1
    fi
    "$PY" -c "import urllib.request; urllib.request.urlopen('http://127.0.0.1:$PORTA/api/state', timeout=1)" \
        2>/dev/null && break
    sleep 0.5
done

URL="http://127.0.0.1:$PORTA/"
if command -v xdg-open >/dev/null; then xdg-open "$URL" >/dev/null 2>&1 &
elif command -v open >/dev/null; then open "$URL" &          # macOS
fi

if command -v zenity >/dev/null; then
    # Closing this window (or pressing the button) stops the server
    zenity --info --title="EVE curation - server running" --width=440 --ok-label="Stop the server" \
        --text="The page is open at\n<b>$URL</b>\n\nResults: $RES\n\nKeep this window open while you work.\nSave your curation on the page before stopping." \
        2>/dev/null
else
    echo "Page: $URL  (Ctrl+C stops the server)"
    wait $PID
fi
