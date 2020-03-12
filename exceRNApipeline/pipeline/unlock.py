import subprocess as sp


def unlock(args):
    cmd = "snakemake --unlock"
    print(cmd)
    sp.run(cmd.split())
    # print("unlock pipeline")

def parse_args(subparsers):
    parser = subparsers.add_parser(
        name="unlock",
        description="unlock snakemake",
        help="unlock snakemake"
    )
    parser.set_defaults(func=unlock)