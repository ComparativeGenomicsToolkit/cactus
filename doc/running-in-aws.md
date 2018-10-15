Running Cactus on AWS
===
Cactus supports running on AWS with an auto-scaling cluster using [Toil](https://toil.readthedocs.io/en/latest/). Check out the [Toil docs on running in the cloud](https://toil.readthedocs.io/en/latest/running/cloud/cloud.html) for the full story, but here's a short walkthrough of running on AWS.

## AWS setup
If you have a fresh AWS account or haven't used EC2 before, you'll need to go through some initial setup.
### Keypair
Make sure you have an AWS keypair ready. [This document](https://docs.aws.amazon.com/AWSEC2/latest/UserGuide/ec2-key-pairs.html) will tell you how to create an AWS keypair that will allow you to log into the instances you create.
### Access keys
You'll also need to have your AWS access credentials set up in `~/.aws/credentials` or the typical AWS environment variables (`AWS_ACCESS_KEY_ID` and `AWS_SECRET_ACCESS_KEY`). [Here is AWS's documentation on setting up your access key](https://docs.aws.amazon.com/IAM/latest/UserGuide/id_credentials_access-keys.html), and [this guide](https://docs.aws.amazon.com/cli/latest/userguide/cli-chap-getting-started.html) will help you set up `~/.aws/credentials`.
### Instance limits
By default, AWS will restrict you to running only a few small instances at a time. If you have a new AWS account, or you're not sure what your limits are, you will probably need to increase them. See [this guide](https://docs.aws.amazon.com/AWSEC2/latest/UserGuide/ec2-resource-limits.html) for information on how to check your existing limits and how to increase them if necessary. AWS support takes only one or two business days to process your request if you have the default "basic" support package.

[Look below](#estimating-the-maximum-number-of-worker-instances-youll-need) for tips on how many instances you may need for your specific alignment. Keep in mind that spot-market instance limits are separate from "on-demand" (non-preemptable) instance limits. (You'll probably want to request a slightly higher limit than you think you need, just in case you want to tweak the number of instances later on.)
## Installing Cactus on your local machine
Follow the steps in the README, making sure to install toil with its extra AWS support:
```
git clone https://github.com/comparativegenomicstoolkit/cactus.git
cd cactus
virtualenv venv
source venv/bin/activate
pip install --upgrade .
pip install --upgrade toil[aws]
```
## Estimating the maximum number of worker instances you'll need
The cluster will automatically scale up and down, but you'll want to set a maximum number of nodes so the scaler doesn't get overly aggressive and waste money, or go over your AWS limits. We typically use `c4.8xlarge` on the spot market for most jobs, and `r4.8xlarge` on-demand for database jobs. Here are some very rough estimates of what we typically use for the maximum of each type (round up):

- `N` mammal-size genomes (~2-4Gb): `(N / 2) * 20` c4.8xlarge on the spot market, `(N / 2)` r3.8xlarge on-demand
- `N` bird-size genomes (~1-2Gb): `(N / 2) * 10` c4.8xlarge on the spot market, `(N / 4)` r3.8xlarge on-demand
- `N` nematode-size genomes (~100-300Mb): `(N / 2)` c4.8xlarge on the spot market, `(N / 10)` r3.8xlarge on-demand
- For anything less than 100Mb, the computational requirements are so small that you may be better off running it on a single machine than using an autoscaling cluster.
## Launching the "leader" instance
Make sure you have your AWS keypair active in your `ssh-agent`, unless you used your existing SSH key. 
If `ssh-agent` isn't started by your operating system, you may have to run `eval $(ssh-agent)` to start it. You can activate the AWS keypair for a session by running:
```
ssh-add path/to/your_aws_ssh_keypair
```


Then launch the leader like so:

```
toil launch-cluster -z us-west-2b <clusterName> --keyPairName <yourKeyPairName> --leaderNodeType t2.medium
```
## Transfer over local input data (or use URLs in your seqfile)
You need to get your actual data to the leader somehow. You can use `http://` and `s3://` to specify FASTA locations in your input file, or you can rsync your data over like so:
```
toil rsync-cluster -z us-west-2b my-cactus-cluster -avP seqFile.txt input1.fa input2.fa :/
```
## Log into the leader and set up the cluster's Cactus installation
Log into the leader by running:
```
toil ssh-cluster -z us-west-2b <clusterName>
```

You should now get a prompt like:
```
root@ip-172-31-34-148:/#
```
indicating that you're on the leader.

Install Cactus in a virtual environment on the leader:
```
apt update
apt install -y git tmux
virtualenv --system-site-packages venv
source venv/bin/activate
git clone https://github.com/comparativegenomicstoolkit/cactus.git
cd cactus
pip install --upgrade .
cd /
```
## Run the alignment (on the leader)
The key parameters you'll care about changing (besides the usual Cactus parameters) are the autoscaling parameters `--nodeTypes`, `--minNodes`, `--maxNodes` and the jobStore location, `aws:<region>:<jobStoreName>`.

You must use the AWS jobstore (or other cloud jobstore, though others may incur data egress charges), not a directory jobstore, because there is no shared filesystem in the cluster. Set the region to whatever region you're running the leader in. The jobStoreName must be globally unique.

The `--nodeTypes` option lets you specify the list of instances you want as well as the price you're willing to pay for spot instances. For example, `c4.8xlarge:0.6` says that we want a c4.8xlarge instance on the spot market, and we're willing to pay up to 60 cents an hour for it. A value like `r3.8xlarge` indicates that we also want on-demand r3.8xlarge instances (for which we pay exactly the on-demand price, which is usually substantially higher). The `--maxNodes` option will let you specify a list containing the maximum number of nodes of each type to use at any given time (in the same order as `--nodeTypes`. The `--minNodes` option should probably be left at 0 unless you really know what you're doing and you're having substantial difficulties with the autoscaler.

```
cactus --nodeTypes c4.8xlarge:0.6,r3.8xlarge --minNodes 0,0 --maxNodes 20,2 --provisioner aws --batchSystem mesos --metrics aws:us-west-2:<jobstoreName> seqFile.txt output.hal
```

This will take a while. You'll want to run this inside something that will preserve your session, like `tmux` or `screen`, so the command doesn't terminate when you disconnect.

## Restarting after failure
If you cancel your run, or it fails for some reason, you can start where you left off by running:
```
cactus [all your usual options...] --restart
```
## Shut down your leader
When the alignment is done, your leader will still be active, costing you a little bit of money every day. Make sure you get rid of it when you're done:
```
toil destroy-cluster -z us-west-2b <yourClusterName>
```
