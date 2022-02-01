#!/usr/bin/env python3
"""
Assembles a trace file for the chromium tracing tool using multiple log files and prints the result.

The tool is available through chromium browsers (e.g. Google Chrome) using the url chrome://tracing or by using the standalone.

Format reference: https://docs.google.com/document/d/1CvAClvFfyA5R-PhYUmn5OOQtYMH4h6I0nSsKchNAySU/preview
"""

import argparse, datetime, json, sys


class StoreDictKeyPair(argparse.Action):
    """
    Helper class used as an action in argparse to store a dictionary of the format
    KEY=VALUE and sets the according namespace attribute
    """
    def __init__(self, option_strings, dest, nargs=None, **kwargs):
        self._nargs = nargs
        super().__init__(option_strings, dest, nargs=nargs, **kwargs)

    def __call__(self, parser, namespace, values, option_string=None):
        my_dict = {}
        for kv in values:
            k, v = kv.split("=")
            my_dict[k] = v
        setattr(namespace, self.dest, my_dict)



def check_and_parse_args():
    """
    Handles the parsing of the command line arguments
    """
    parser = argparse.ArgumentParser(description="Assembles a trace file for "
                                     "the chromium tracing tool using multiple "
                                     "log files and prints it.")
    parser.add_argument("logs", action=StoreDictKeyPair,
                        nargs="+", metavar="PARTICIPANT=LOGFILE")
    parser.add_argument("-p", "--pretty",  action="store_true",
                        help="Print the JSON in a pretty format.")
    parser.add_argument("-d", "--default",  default="default", metavar="CATEGORY",
                        help="The default category for unknown events.")
    parser.add_argument("-m", "--mapping", metavar="FILE",
                        help="The file containing mappings from event-names to categories.")
    parser.add_argument("-g", "--noglobal", action="store_true",
                        help="Ignore the global event.")
    parser.add_argument("-k", "--ranks", type = int, nargs="+", metavar="RANK",
                        help="Only output the given ranks.")
    parser.add_argument("-t", "--maxtime", type = int, default = -1,
                        help = "Maximum time stamp to convert, milliseconds after init of first rank.")
    parser.add_argument("--no-normalize", action = "store_true",
                        help = "Disable time normalization amoung participants")


    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    return parser.parse_args()


def eprint(*args, **kwargs):
    """ Debug function, prints to stderr. """
    print(*args, file=sys.stderr, **kwargs)
    
def normalize_times(*args):
    """ Normalize times to first t0 amoung all participants. """
    # Find the minimum Initialized time among all participants
    fmt_str =  r"%Y-%m-%dT%H:%M:%S.%f" # replaces the fromisoformatm, not available in python 3.6
    minT = min([datetime.datetime.strptime(d["Initialized"], fmt_str) for d in args])
    
    for d in args:
        init = datetime.datetime.strptime(d["Initialized"], fmt_str)
        delta = init - minT
        for ranks in d["Ranks"]:
            for sc in ranks["StateChanges"]:
                sc["Timestamp"] = int(sc["Timestamp"] + (delta.total_seconds() * 1000))
                
    return args


def build_process_name_entry(name, pid):
    """
    Builds a dictionary entry representing a metadata event to name a process id.
    """
    return { "name": "process_name",
             "ph": "M",
             "pid": pid,
             "tid": 0,
             "args": {"name": name}
    }


def build_thread_name_entry(name, pid, tid):
    """
    Builds a dictionary entry representing a metadata event to name a thread id
    local to a process.
    """
    return { "name": "thread_name",
             "ph": "M",
             "pid": pid,
             "tid": tid,
             "args": {"name": name}
    }


def main():
    args = check_and_parse_args()

    event_mapping = {}
    if args.mapping:
        event_mapping = json.loads(args.mapping)

    logs = args.logs.items()
    pids = range(len(logs))
    jsons = [json.load(open(i[1])) for i in logs]
    if not args.no_normalize:
        jsons = normalize_times(*jsons)
    
    # The output will be in the JSONArray format described in the specification
    traces = []
    for pid, participant, data in zip(pids, logs, jsons):
        # The pid identifies each participant and is used as process id
        traces.append(build_process_name_entry(participant[0], pid))
        
        for rank, rank_data in enumerate(data["Ranks"]):
            if (args.ranks) and (rank not in args.ranks):
                continue
            traces.append(build_thread_name_entry("Rank {:4d}".format(rank), pid, rank))
            
            for sc in rank_data["StateChanges"]:
                if args.noglobal and sc["Name"] == "_GLOBAL":
                    continue
                if args.maxtime > -1 and sc["Timestamp"] > args.maxtime:
                    continue

                # The current log format contains begin and end timestamps of
                # events which corresponds to the specified duration events    
                event = {
                    "name": sc["Name"],
                    "cat": event_mapping.get(sc["Name"], args.default),
                    "tid": rank,
                    "pid": pid,
                    "ts": sc["Timestamp"] * 1000, # convert from ms to µs
                    "ph" : "B" if sc["State"] == 1 else "E"
                }
                traces.append(event)

    if args.pretty:
        print(json.dumps(traces, indent=2))
    else:
        print(json.dumps(traces))


if __name__ == "__main__":
    main()