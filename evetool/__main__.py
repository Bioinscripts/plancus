"""
Command line of the EVE tool.

  eve-tool run ...            the pipeline (same options as before; eve-tool run --help)
  eve-tool curate ...         the curation server (open_curation.sh starts it with a folder dialog)
  eve-tool databases ...      download and install the databases
  eve-tool build-databases    rebuild the viral databases from new source releases, same criteria
  eve-tool pack-databases     pack this machine's curated viral databases (for Zenodo)
  eve-tool check              tools, Python libraries and databases of this machine
  eve-tool config             the paths in use

Any module can also be run directly: python -m evetool.annotation.repeats -h
"""
import runpy
import sys

COMANDOS = {
    "run": "evetool.main",
    "curate": "evetool.curation.server",
    "databases": "evetool.databases.prepare",
    "build-databases": "evetool.databases.build",
    "pack-databases": "evetool.databases.pack",
    "check": "evetool.check",
    "config": "evetool.config",
}


def main():
    if len(sys.argv) < 2 or sys.argv[1] in ("-h", "--help") or sys.argv[1] not in COMANDOS:
        print(__doc__.strip())
        sys.exit(0 if len(sys.argv) < 2 or sys.argv[1] in ("-h", "--help") else 2)
    comando = sys.argv.pop(1)
    sys.argv[0] = f"eve-tool {comando}"
    runpy.run_module(COMANDOS[comando], run_name="__main__")


if __name__ == "__main__":
    main()
