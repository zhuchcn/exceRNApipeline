import subprocess as sp
import os
import argparse


def run(args):
    smk_dir = os.path.join(os.path.dirname(__file__), "smk")
    snakefile = os.path.join(smk_dir, "Snakefile")
    cmd = f"snakemake --snakefile {snakefile}"
    
    if args.configfile is None:
        raise argparse.ArgumentError("--configfile can not be None")
    cmd += f" --configfile {args.configfile}"
    if args.directory:
        cmd += f" -d {args.directory}"
    if args.report:
        cmd += f" --report {args.report}"
    cmd += f" -j {args.jobs}"
    cmd += f" --latency-wait {60 if args.slurm_config else 5}"
    if args.use_singularity:
        cmd += " --use-singularity"
    if args.singularity_args:
        cmd += f" --singularity-args \"{args.singularity_args}\""
    if args.printshellcmds:
        cmd += " -p"
    if args.dry_run:
        cmd += " -n"
    cmd += f" {args.snakemake_args}"

    if args.slurm_config:

        sm_params = ["job-name","partition", "mail-user", "mail-type", 
                "output", "error", "time", "nodes", "mem", "exclude"]
        sm_args = "--ntasks {threads} --mem-per-cpu {cluster.mem-per-cpu}"
        for param in sm_params:
            sm_args += f" --{param} {{cluster.{param}}}"

        user_mail = os.environ["USER_EMAIL"] or "None"
        cwd = args.directory or "."
        slurm_cmd = f"srun -N 1 -n 1 -t 14-0 --mail-user {user_mail}" +\
            f" --mail-type ALL -o {cwd}/slurmout/snakemake_%A.log" +\
            f" -e {cwd}/slurmout/snakemake_%A.err "
        if args.slurm_args:
            slurm_cmd += args.slurm_args
        cmd = slurm_cmd + " " + cmd + f" --cluster-config {args.slurm_config}"+\
            f" --cluster \"sbatch {sm_args}\" &"
    
    print(cmd)
    sp.call(cmd, shell=True)

def parse_args(subparsers):
    parser = subparsers.add_parser(
        name="run",
        description="run pipeline",
        help="run pipeline"
    )
    parser.set_defaults(func=run)
    parser.add_argument(
        '-n', '--dry-run', action="store_true",
        help="""Do not execute anything, and display what would be done. If you
        have a very large workflow, use –dry-run –quiet to just print a summary
        of the DAG of jobs. Default: False"""
    )
    parser.add_argument(
        '-p', '--printshellcmds', action="store_true", help="""Print out the
        shell commands that will be executed. Default: False"""
    )
    parser.add_argument(
        '--slurm-config', type=str, help="""Path to the slurm config file.
        Use `pipeline init` to generate a template. Pipeline will be ran with
        slurm with this flag."""
    )
    parser.add_argument(
        '--configfile', type=str, help="""Path to the workflow config file.
        Use `pipeline init` to generate a template."""
    )
    parser.add_argument(
        '-d', '--directory', type=str, help="""Specify working
        directory (relative paths in the snakefile will use this as their
        origin)."""
    )
    parser.add_argument(
        '--report', nargs="?", const="report.html", help="""Create an HTML
        report with results and statistics. If no filename is given,
        report.html is the default."""
    )
    parser.add_argument(
        '-j', '--jobs', type=int, default = 999, help="""Use at most N cores 
        in parallel. If N is omitted or ‘all’, the limit is set to the number
        of available cores. Default: 999"""
    )
    parser.add_argument(
        '-w', '--latency-wait', type=int, help="""Wait given
        seconds if an output file of a job is not present after the job
        finished. This helps if your filesystem suffers from latency. 
        Default: 5 for local run, and 60 on slurm"""
    )
    parser.add_argument(
        '--use-singularity', action="store_true", help="""Run pipeline within
        the singularity (docker) container. Default: False"""
    )
    parser.add_argument(
        '--singularity-args', type=str, default="",
        help="Pass additional args to singularity. Must be quoted. " +\
        "Example: \"--bind /scratch:/scratch\""
    )
    parser.add_argument(
        '--snakemake-args', type=str, default="",
        help="Pass additional args to snakemake. Must be quoted. " +\
        "Example: \"--verbose\""
    )
    parser.add_argument(
        '--slurm-args', type=str, default="",
        help="Pass additional args to slurm. Must be quoted. " +\
        "Example: \"--mail-user jdoe@gmail.com\""
    )
