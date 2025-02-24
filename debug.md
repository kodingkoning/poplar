# How to Address Poplar Errors

This document will expand as more types of errors are documented.

## Common Errors

### Input Files

At the beginning of execution, Poplar will confirm that all the referenced files in the database catalog do exist. However, if they do not contain valid data for their file type, that will cause a failure of the app using the file as input.

### Dependencies

Poplar will also confirm that it can find the commands required in the bash apps. It only confirms that `which` results in something other than `None`, so if there is an issue with the binary or a system incompatibility, Poplar will not catch that issue.

### Configuration

If the configuration is mismatched with the system, then Poplar will not be able to launch apps. For example, if you are not using a slurm system, then a config using slurm cannot succeed. If you are using a job submission system, the name of the queue to use must be accurate.

## Navigating the Ouput

### Standard Out of Managing Process

The main Python script will print message such as "Building BLAST Database". In most of the cases, this output means that the apps have been sent to Parsl, not that they have necessarily started execution. The runinfo and monitoring info can show what has started and completed execution.

### Run Info

Parsl saves the output from the apps within `runinfo/##`. If the error occurred within a the execution of a tool that was called in a bash app, then `runinfo/##/submit_scripts/parsl*.sh.out` and `runinfo/##/submit_scripts/parsl*.sh.err` are the most informative. These files contain stdout and sterr. If the fatal error immediately stopped the workflow from running, then the error will be at the end. However, other apps can continue to execute until there is a dependency issue between a failed app and a queued app, so the error may be earlier in the output.

If the error is in the Python part of a bash app or in a python app, then the log in `runinfo/##/parsl.log` will have information on the error. This file is more expansive but can be quite informative as well. Searching for "ERROR" can help quickly find failures, while "error" will often be referencing the redirection of standard error and not indicate any failure.

Some commands within Poplar will have additional logs, and these will be referenced in stdout and will be located in the temporary directory.

### Intermediate Output

Noting the presence or absense of output in the temporary directory can be useful for identifying where the error occurred. For example, if it appears that the app to construct the gene tree is failing, but there are no alignment files in the directory, then the failure must have actually occurred before the gene tree construction.

### Checkpointing

By default, Poplar will use checkpointing. If you rerun Poplar with the same command line arguments (even if the files have changed internally) it will use the checkpoints from previous runs. Depending on the failure type, Poplar may rerun failed apps. For example, if the submitted job ran out of time, Poplar should rerun that app. Other times, such as if the failure cascaded through apps before Poplar exited, Poplar will not rerun the apps with invalid output, and the checkpoints need to be deleted from `runinfo/##/checkpoint`.

## Monitoring

### Launching Monitoring

Running the Parsl monitoring requires including monitoring in the config and then starting `parsl-visualize` from the directory where `runinfo` is located. The depenencies for `parsl-visualize` are installed via `setup.sh` but are not required for running Poplar. It will then tell you where it is running, and that will often be `http://127.0.0.1:8080`.

Launching monitoring from a remote connection requires connecting to ssh with `ssh -L 50000:127.0.0.1.:8080 username@cluster`, running `parsl-visualize` on the remote machine, and then viewing in the browser at `http://127.0.0.1:50000` on the local machine.

If you delete the files in `runinfo`, then you must restart `parsl-visualize` for it to update with the new files.

### Viewing Monitoring

The front page of monitoring will list each workflow with information available in `runinfo`. Clicking on one of them will show the summary and the apps that have been queued and their statuses. The monitoring information is clearly labeled. The workflow DAG colored by task states can be particularly useful for understanding the progress of Poplar or the point of failure.
