#!/usr/bin/env python3
"""
postprocessing/cli.py

Top‐level CLI for postprocessing pipeline:
  - loaders
  - transformers
  - exporters
  - transporters (backup + transport)
"""
import click
import os
from dotenv import load_dotenv
from pathlib import Path
import logging
import sys
import shutil
from ..LogFormatter import ColorFormatter


# ---------------------------- #
# SETUP
# ---------------------------- #

handler = logging.StreamHandler(sys.stdout)
handler.setFormatter(ColorFormatter())
root = logging.getLogger()
root.setLevel(logging.DEBUG)
root.addHandler(handler)

# __file__ is the location of this script; adjust if your .env lives elsewhere
env_path = Path(__file__).parent / "ppcfg.env"
load_dotenv(dotenv_path=env_path)


from .exporters.refertazione_exporters import main as refertazione_exporters_main
from .exporters.variant_exporters      import main as variant_exporters_main
from .exporters.checklist_exporters    import main as checklist_exporters_main
from .transporters.backup              import run_backup
from .transporters.transport           import run_transport

# ---------------------------- #
# CLI
# ---------------------------- #
@click.group()
def cli():
    """Pipeline postprocessing commands."""
    pass


# ─── TRANSFORMER SUBGROUP ───────────────────────────────────────────────────────
# @cli.group()
# def transform():
#     """Apply transforms to your data."""
#     pass

# @transform.command("all")
# @click.option("--config", default=None,
#               help="Path to transform config")
# def run_all_transforms(config):
#     """Run all data transformers."""
#     # run_transform(config)
#     click.echo(f"Would run transformers with config={config!r}")


# ─── EXPORTER SUBGROUP ──────────────────────────────────────────────────────────
@cli.group()
def export():
    """Export to flat-files or intermediate formats."""
    pass

@export.command("refertazione")
@click.option("--in-dir", default=None,
              help="Where to write refertazione CSVs")
def export_ref(in_dir):
    if in_dir is None:
        in_dir = os.environ["STAGING_PATH"]

    """Export refertazione tables."""
    refertazione_exporters_main(Path(in_dir))

@export.command("variants")
@click.option("--in-dir", default=None,
              help="Where to write variant CSVs")
def export_var(in_dir):
    if in_dir is None:
        in_dir = os.environ["STAGING_PATH"]

    """Export variant tables."""
    variant_exporters_main(Path(in_dir))

@export.command("checklist")
@click.option("--in-dir", default=None,
              help="Where to write checklist CSVs")
def export_checklist(in_dir):
    if in_dir is None:
        in_dir = os.environ["STAGING_PATH"]

    """Export checklist tables."""
    checklist_exporters_main(Path(in_dir))


# ─── TRANSPORTER SUBGROUP ───────────────────────────────────────────────────────
@cli.group()
def transport():
    """Move files around (backup, staging)."""
    pass

@transport.command("backup")
@click.option("--destination", default=None,
              help="Where to write; choices: r, b, z")
@click.option("--in-dir", default=None,
              help="Where to read results from.")
@click.option("--in-fq", default=None,
              help="Where to take fq from.")
def backup(destination, in_dir, in_fq): # TODO: this will use a dedicated backup function implemented in backup.py; right now it just calls run_transport
    """Copy raw results into backup location."""
    run_backup(destination, in_dir, in_fq)

@transport.command("stage")
@click.option("--in-dir", default=None,
              help="Where to read variant CSVs")
def stage(in_dir):
    """Copy filtered results into staging area."""
    dest = Path(os.environ["STAGING_PATH"])
    run_transport(Path(in_dir), dest)


# ─── CLEAN SUBGROUP ───────────────────────────────────────────────────────
@cli.group()
def clean():
    """Clean up."""
    pass

@clean.command("stage")
def clean_stage():
    """Clean up staging area."""
    dest = Path(os.environ["STAGING_PATH"])
    if not dest.is_dir():
        raise NotADirectoryError(f"{dest!r} is not a directory")
    
    for item in dest.iterdir():
        if item.is_dir():
            shutil.rmtree(item)
        else:
            item.unlink()

@clean.command("backup")
def clean_backup():
    """Clean up backup area."""
    pass

if __name__ == "__main__":
    cli()
