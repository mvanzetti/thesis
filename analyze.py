#! /usr/bin/env python
# -*- coding: utf-8 -*-
import logging
import sys
from overlaps_db.cli import cli
from overlaps_db.analysis.encode_analyzer import EncodeAnalyzer

log = cli.start_logging(sys.argv[0], level=logging.DEBUG)


def main(args):
    storage_dir = "/Users/manuel/development/thesis/storage"

    if args.source == 'ENCODE':
        log.info("Initializing EncodeOverlapper")
        analyzer = EncodeAnalyzer(storage_path=storage_dir)

        if args.target == 'FANTOM':
            log.info("Overlapping with FANTOM")

            if args.analysis == 'overlap':
                analyzer.perform_overlap_analysis_with_fantom_permissive(assembly=args.assembly, method=args.method,
                                                                         overlap_intervals=args.overlap_intervals,
                                                                         samples_num=args.samples_num,
                                                                         biosample_type=args.biosample_type)

            if args.analysis == 'reldist':
                analyzer.perform_reldist_analsys_with_fantom_permissive(assembly=args.assembly, method=args.method,
                                                                        biosample_type=args.biosample_type)

        elif args.target == 'RepeatMasker':
            log.info("Overlapping with RepeatMasker")

            if args.analysis == 'overlap':
                analyzer.perform_overlap_analysis_with_repeatmasker(assembly=args.assembly, method=args.method,
                                                                    overlap_intervals=args.overlap_intervals,
                                                                    samples_num=args.samples_num,
                                                                    biosample_type=args.biosample_type,
                                                                    repeat_class=args.repeat_class)

            if args.analysis == 'reldist':
                analyzer.perform_reldist_analsys_with_repeatmasker(assembly=args.assembly, method=args.method,
                                                                   biosample_type=args.biosample_type,
                                                                   repeat_class=args.repeat_class)

        else:
            print("The selected target is not available")

        log.info("%s analysis completed", args.analysis)
    else:
        print("selected source not implemented yet")


if __name__ == '__main__':
    cli.preamble("Analyze")
    # cli.print_banner("/Users/manuel/development/thesis/overlaps_db/resources/banner.txt")

    parser = cli.get_default_parser()

    parser.add_argument("source", help='Select one of the available sources',
                        choices=['ENCODE', 'FANTOM', 'ENCODE_FANTOM'])
    parser.add_argument("target",
                        help='Select the target to consider for the overlap analysis',
                        choices=['FANTOM', 'RepeatMasker'])
    parser.add_argument("--analysis",
                        help='Select the kind of analysis you want to perform', action="store", dest="analysis",
                        choices=['overlap', 'reldist'])
    parser.add_argument('--assembly', action="store", dest="assembly",
                        help='The assembly to use to build the overlaps files. \
                        For the human genome assembly, type hg19.',
                        default='hg19')
    parser.add_argument('--method', action="store", dest="method",
                        help='Filter the source enhancer list by the provided method',
                        default='DNase_H3K27ac')
    parser.add_argument('--type', action="store", dest="biosample_type",
                        help='If specified, consider a specific biosample type')
    parser.add_argument('--repeat_class', action="store", dest="repeat_class",
                        help='If target is RepeatMasker, please specify a repeat class family (e.g.: "SINE/Alu")')
    parser.add_argument('--intervals', action="store", dest="overlap_intervals",
                        help='The number of intervals to split the min overlap requested for the overlapping tests',
                        type=int, default=10)
    parser.add_argument('--samples', action="store", dest="samples_num",
                        help='The number of samples to consider while building random null models', type=int,
                        default=20)

    received_args = parser.parse_args()
    log.info("Performing %s analysis for %s and %s", received_args.analysis, received_args.source, received_args.target)
    log.info("Optional parameters: Assembly: %s, Method: %s", received_args.assembly,
             received_args.method, )

    if received_args.analysis == 'overlap':
        print("Intervals: %d, Samples: %d", received_args.overlap_intervals, received_args.samples_num)

    if received_args.biosample_type:
        log.info("Consider only %s", received_args.biosample_type)

    if received_args.target == 'RepeatMasker':
        if not received_args.repeat_class:
            log.error("Please specify a valid repeat class")
            sys.exit()

        log.info("Consider repeat class %s", received_args.repeat_class)

    sys.exit(main(received_args))
