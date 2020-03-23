import subprocess as sp
import os


def unlock(args):
    smk_dir = os.path.join(os.path.dirname(__file__), "smk")
    snakefile = os.path.join(smk_dir, "Snakefile")
    cmd = f"snakemake --snakefile {snakefile}"
    if not os.path.exists(snakefile):
        raise ValueError(f"snakefile {snakefile} doesn't exist.")
    cmd += f" --configfile {args.configfile} --unlock"
    if args.cleanup_metadata:
        cmd += f" --cleanup-metadata {args.cleanup_metadata}"
    if args.rerun_incomplete:
        cmd += " --rerun-incomplete"
    if args.extra_args:
        cmd += f" {args.extra_args}"
    print(cmd)
    sp.run(cmd.split())
    # print("unlock pipeline")

def parse_args(subparsers):
    parser = subparsers.add_parser(
        name="unlock",
        description="unlock snakemake",
        help="unlock snakemake"
    )
    parser.add_argument(
        'configfile', type=str, help="""Path to the workflow config file. Use 
        `pipeline init` to generate a template."""
    )
    parser.add_argument(
        '--cleanup-metadata', '--cm', type=str, help="""Cleanup the metadata of given
        files. That means that snakemake removes any tracked version info, and
        any marks that files are incomplete."""
    )
    parser.add_argument(
        '--rerun-incomplete', '--ri',  action="store_true", help="""Re-run all
        jobs the output of which is recognized as incomplete. Default: False"""
    )
    parser.add_argument(
        '--extra-args', type=str, help="""Extra args to be parsed to
        snakemake."""
    )
    parser.set_defaults(func=unlock)