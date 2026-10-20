#!/usr/bin/env bash
# Runs the EVE pipeline with results and log in their own folder, never next to the code.
#
#   ./run_pipeline.sh <run_name> --species-tsv species.tsv --genomes-dir /path/genomes [more eve-tool run options]
#   ./run_pipeline.sh llong_chr3 -i genome.fna --gff genome.gff3
#
# Results: <results_dir>/<run_name>/  (log: pipeline.log; folders set in eve_config.ini, see ./eve-tool config)
# The databases below are the defaults; pass --dmnd-db, --hmm-db... to override them.
# CD-Search and Pfam run by default when installed; NO_CDD=1 or NO_PFAM=1 turn them off.
# EVE_PYTHON chooses another interpreter (default: the active conda env, or eve_env).

set -u
AQUI="$(cd "$(dirname "$(readlink -f "$0")")" && pwd)"
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

# Database and results folders from eve_config.ini / EVE_* variables (see ./eve-tool config)
CAMINHOS="$("$PY" -m evetool.config --shell)" || exit 1
eval "$CAMINHOS"
RESULTADOS="$EVE_RESULTS_DIR"
DMND_DB="$EVE_VIRAL_DMND"
HMM_DB="$EVE_HMM_DB"
DMND_GERAL="$EVE_GERAL_DMND"
if [ ! -f "$DMND_DB" ] || [ ! -f "$HMM_DB" ]; then
    echo "Viral databases not installed (./eve-tool databases). Check with: ./eve-tool check"
    exit 1
fi

if [ $# -lt 2 ] || [[ "$1" == -* ]]; then
    sed -n '2,10p' "$0" | sed 's/^# \{0,1\}//'
    exit 1
fi
NOME="$1"; shift
SAIDA="$RESULTADOS/$NOME"
if [ -e "$SAIDA" ]; then
    echo "Already exists: $SAIDA"
    echo "Choose another run name (results are never overwritten by this script)."
    exit 1
fi

# Default databases, unless the user passed their own
EXTRA=()
[[ " $* " == *" --dmnd-db "* || " $* " == *" --fasta-dir "* ]] || EXTRA+=(--dmnd-db "$DMND_DB")
[[ " $* " == *" --hmm-db "* || " $* " == *" --hmm-dir "* ]] || EXTRA+=(--hmm-db "$HMM_DB")
[[ " $* " == *" --diamond-db-geral "* || ! -f "$DMND_GERAL" ]] || EXTRA+=(--diamond-db-geral "$DMND_GERAL")
# CD-Search and Pfam only when their databases are installed (prepare.py --cdd --pfam)
[[ " $* " == *" --cdsearch "* || -n "${NO_CDD:-}" || ! -e "$EVE_CDD_RPSDB.pal" && ! -e "$EVE_CDD_RPSDB.rps" ]] \
    || EXTRA+=(--cdsearch --cdsearch-mode local)
[[ " $* " == *" --pfam "* || -n "${NO_PFAM:-}" || ! -f "$EVE_PFAM_HMM" ]] || EXTRA+=(--pfam)
[[ " $* " == *" -t "* || " $* " == *" --threads "* ]] || EXTRA+=(-t "$(nproc 2>/dev/null || echo 4)")

mkdir -p "$SAIDA"
LOG="$SAIDA/pipeline.log"
{
    echo "# $(date '+%F %T')  code: $AQUI"
    echo "# $PY -m evetool.main -o $SAIDA ${EXTRA[*]} $*"
} > "$LOG"
echo "Results: $SAIDA"
echo "Log:     $LOG   (follow it with: tail -f \"$LOG\")"

INICIO=$(date +%s)
"$PY" -u -m evetool.main -o "$SAIDA" "${EXTRA[@]}" "$@" 2>&1 | tee -a "$LOG"
CODIGO=${PIPESTATUS[0]}
echo "# exit code $CODIGO after $(( ($(date +%s) - INICIO) / 60 )) min  ($(date '+%F %T'))" | tee -a "$LOG"
[ "$CODIGO" -eq 0 ] && echo "Open the curation with:  $AQUI/open_curation.sh \"$SAIDA\""
exit "$CODIGO"
