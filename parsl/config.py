import parsl
from parsl.config import Config
from parsl.providers import SlurmProvider
from parsl.executors import HighThroughputExecutor
from parsl.launchers import SrunLauncher

""" This config assumes that it is used to launch parsl tasks from the login nodes
of Sandia's Kahuna. Each job submitted to the scheduler will request 2 nodes for 10 minutes.
"""
config = Config(
     executors=[
          HighThroughputExecutor(
               label="kahuna",
               worker_debug=False,
               cores_per_worker=16.0,  # each worker uses a full node
               provider=SlurmProvider(
                    partition='compute',  # partition
                    nodes_per_block=2,  # number of nodes
                    init_blocks=1,
                    max_blocks=1,
                    scheduler_options='',
                    cmd_timeout=60,
                    walltime='02:10:00',
                    launcher=SrunLauncher(),
                    worker_init='conda activate parsl_py38; export PATH=$PATH:/home/erkonin/daily_notes/BLAST/ncbi-blast-2.14.1+/bin:/home/erkonin/.conda/envs/pipeline_env/bin:/home/erkonin/daily_notes/gene_trees/raxml/from_annotation/astralpro/ASTER-Linux/bin:/home/erkonin/daily_notes/raxml-ng/build/bin:/home/erkonin/daily_notes/gene_trees',  # requires conda environment with parsl
               ),
          )
     ],
)
