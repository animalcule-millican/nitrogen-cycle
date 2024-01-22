#!/usr/bin/env python3
#!/usr/bin/env python

import sys
import glob
import htcondor
from htcondor import JobEventType
from os.path import join


def print_and_exit(s):
    print(s)
    exit()


jobID, UUID, clusterID = sys.argv[1].split("_")
jobDir = f"/home/glbrc.org/millican/repos/nitrogen-cycle/workflow/log/*/{jobID}_{UUID}"
logdir = f"/home/glbrc.org/millican/repos/nitrogen-cycle/workflow/log/*/{jobID}_{UUID}/*.log"
#jobDir = "/home/glbrc.org/millican/repos/Slime_Py/workflow/log/*/{}_{}".format(jobID, UUID)
join(jobDir, f"{jobID}-{UUID}.$(Cluster).log")
jobLog = glob.glob(logdir)[-1]


failed_states = [
    JobEventType.JOB_HELD,
    JobEventType.JOB_ABORTED,
    JobEventType.EXECUTABLE_ERROR,
]

try:
    jel = htcondor.JobEventLog(join(jobLog))
    for event in jel.events(stop_after=5):
        if event.type in failed_states:
            print_and_exit("failed")
        if event.type is JobEventType.JOB_TERMINATED:
            if event["ReturnValue"] == 0:
                print_and_exit("success")
            print_and_exit("failed")
except OSError as e:
    print_and_exit("failed: {}".format(e))

print_and_exit("running")
