#!/usr/bin/python2.7
# =======================================================================================
# Copyright Asthma Collaboratory (2018)
# coded by Christopher R. Gignoux and Zach Pincus
# modified by Kevin L. Keys
#
# This script creates a barebones job scheduling mechanism.
# The jobs are read from a filename, one job per line.
# At any one time, up to MAX_JOBS will run at once.
# =======================================================================================

# =======================================================================================
# load libraries
# =======================================================================================
import threading
import subprocess
from optparse import OptionParser
from sys import stdout

# =======================================================================================
# command line options
# =======================================================================================
# add command line options
parser = OptionParser()
parser.add_option('-f', '--file', dest='filename', type = 'string', default='joblist.sh', help='text file of commands to multithread')
parser.add_option('-j', '--jobs', dest='MAX_JOBS', type = 'int', default=8, help='number of parallel jobs to run')

# parse command line options
(options, args) = parser.parse_args()
filename = options.filename
MAX_JOBS = options.MAX_JOBS

# =======================================================================================
# subroutines
# =======================================================================================

# subroutine to make a list of parallel jobs
def generate_jobs():
    # your code here:
    # for each job, a list of strings corresponding to the executable and command-line args to run:
    # e.g. ['ls', '-l', 'foo/bar']
    # here 'tis for this job
    return [line.strip() for line in file(filename).readlines()]

def run_jobs():
    njobs = len(JOB_LIST)
    while True:
        LIST_LOCK.acquire()
        try:
            job = JOB_LIST.pop()
            print job
            print '\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b%s left' % (len(JOB_LIST) + 1)
            stdout.flush()
        except IndexError:
            # return if the job list is empty
            return
        finally:
            LIST_LOCK.release()
        # block until the job ends
        val = subprocess.call(job, shell=True)
    if val != 0:
        # real crude error reporting
        print "Job '%s' failed." % (' '.join(job))

# =======================================================================================
# executable code
# =======================================================================================

# generate list of jobs to execute
JOB_LIST = list(generate_jobs())

# start a threadlock mechanism to control jobs
LIST_LOCK = threading.Lock()

# loop through jobs
for i in range(MAX_JOBS):
    thread = threading.Thread(target = run_jobs)
    thread.start()
