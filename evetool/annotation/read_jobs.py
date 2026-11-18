#!/usr/bin/env python3
"""
Background read-mapping jobs started from the curation page (small RNA with bowtie, long and
other RNA with minimap2). One job runs at a time, so two mappings never share the CPUs; the
others wait in line. The page polls `situacao()` for the step each job is in.
"""
import threading
import time
import traceback

_trava = threading.Lock()          # guards the job list
_vez = threading.Lock()            # one mapping at a time
_jobs = []                         # newest last
_contador = [0]


def iniciar(tipo: str, species: str, nome: str, funcao) -> dict:
    """
    Run `funcao(passo)` in a thread; `passo(text)` reports the current step. Returns the job
    record ({id, kind, species, name, state, step, error, started, finished, result}).
    """
    with _trava:
        if any(j["state"] in ("queued", "running") and j["kind"] == tipo and j["species"] == species
               and j["name"] == nome for j in _jobs):
            raise RuntimeError(f"{nome} is already being mapped")
        _contador[0] += 1
        job = {"id": _contador[0], "kind": tipo, "species": species, "name": nome, "state": "queued",
               "step": "waiting for the previous mapping", "error": None, "started": time.time(),
               "finished": None, "result": None}
        _jobs.append(job)
        del _jobs[:-30]

    def correr():
        with _vez:
            job["state"], job["step"] = "running", "starting"
            try:
                job["result"] = funcao(lambda texto: job.__setitem__("step", str(texto)))
                job["state"], job["step"] = "done", "done"
            except Exception as erro:
                traceback.print_exc()
                job["state"], job["error"], job["step"] = "error", str(erro), "failed"
            finally:
                job["finished"] = time.time()

    threading.Thread(target=correr, daemon=True, name=f"map_{tipo}_{nome}").start()
    return dict(job)


def situacao() -> list:
    """Every job of this server session, newest first, with its running time in seconds."""
    agora = time.time()
    with _trava:
        return [{**{k: v for k, v in j.items() if k != "result"},
                 "seconds": round((j["finished"] or agora) - j["started"])} for j in reversed(_jobs)]


def ativo() -> bool:
    with _trava:
        return any(j["state"] in ("queued", "running") for j in _jobs)
