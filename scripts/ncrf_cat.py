#!/usr/bin/env python
"""
Concatenate several output files from Noise Cancelling Repeat Finder.
"""

from sys   import argv,stdin,stdout,stderr,exit
from os    import path as os_path
from errno import EPIPE
from gzip  import open as gzip_open


def usage(s=None):
    message = """
usage: ncrf_cat <file1> [<file2> ...] [--markend]
  <file1>    an output file from Noise Cancelling Repeat Finder
  <file2>    another output file from Noise Cancelling Repeat Finder
  --markend  assume end-of-file markers are absent in the input, and add an
             end-of-file marker to the output
             (by default we require inputs to have proper end-of-file markers)

Concatenate several output files from Noise Cancelling Repeat Finder.  This
is little more than copying the files and adding a blank line between the
files.

It can also be used to verify that the input files contain end-of-file markers
i.e. that they were not truncated when created."""

    if (s == None): 
        exit (message)
    else: 
        exit("%s\n%s" % (s ,message))


def _open(filename):
    if (filename.endswith(".gz")) or (filename.endswith(".gzip")):
        return gzip_open(filename, "rt")
    else:
        return open(filename, "rt")

def main():

    # parse the command line
    requireEof    = True
    markEndOfFile = False
    filenames     = []

    for arg in argv[1:]:
        if (arg in ["--noendmark", "--noeof", "--nomark"]):   # (unadvertised)
            requireEof = False
        elif (arg in ["--markend]", "--markeof"]):
            requireEof = False
            markEndOfFile = True
        else:
            filenames += [arg]

    if not filenames:
        usage("you have to give me at least one file")

    # copy the files;  note that we don't bother (or care) to verify that they
    # are really output from ncrf
    for (ix, filename) in enumerate(filenames):
        if (ix > 0):
            print()
        eofMarkerSeen = False

        with _open(filename) as f:
            for line in f:
                line = line.rstrip("\n")
                if (eofMarkerSeen) and (line != ""):
                    exit(f"{os_path.basename(argv[0])}: \"{filename}\" contains additional stuff after end marker (starting with \"{line[:10]}\")")

                if (line == "# ncrf end-of-file"):
                    eofMarkerSeen = True
                    markEndOfFile = True
                    continue
                if (not eofMarkerSeen):
                    try:
                        print(line)
                    except IOError as ex:
                        # "Broken pipe" can happen when downstream tools reject
                        # our output as their input
                        if (ex.errno == EPIPE):
                            exit(f"{os_path.basename(argv[0])}: [Errno {ex.errno}] Broken pipe")

        if (requireEof) and (not eofMarkerSeen):
            exit(f"{os_path.basename(argv[0])}: \"{filename}\" may have been truncated (end marker is absent)")
           
    if (markEndOfFile):
        print("# ncrf end-of-file")

if __name__ == "__main__": 
    main()
