import argparse
import getpass
from mcod_prep.utils.mcod_cluster_tools import submit_mcod


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Launcher for map_spec_drug_bug.py")
    parser.add_argument("source", nargs="+")
    write_parser = parser.add_mutually_exclusive_group(required=True)
    write_parser.add_argument(
        '--write', dest='write', action='store_true',
        help='Write the output file from mapping '
        '(only use once your cleaning script has been approved)'
    )
    write_parser.add_argument(
        '--no-write', dest='write', action='store_false',
        help='Do not write the output file from mapping (used for testing purposes)'
    )
    args = parser.parse_args()
    assert {'all', 'new'}.isdisjoint(set(args.source)),\
        "The 'all' and 'new' options don't make sense in the launcher yet"
    print(f"Submitting jobs for sources {args.source}")

    mem = {'US_MedMined_BD': "200G"}
    runtime = {'US_MedMined_BD': "04:00:00"}
    cores = {'US_MedMined_BD': 50}

    worker = "FILEPATH"
    for source in args.source:
        jobname = f"map_amr_{source}"
        submit_mcod(
            jobname, language='python', worker=worker,
            cores=cores.get(source, 1), memory=mem.get(source, "10G"),
            runtime=runtime.get(source, '00:30:00'),
            params=[source, '--write'], logging=True
        )
