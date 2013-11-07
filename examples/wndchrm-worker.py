watch tubes:
	jobs
	job-deps
	jobs-ready

jobs gets reserved, deps added to deps-id

import time # for sleep
import ConfigParser
import io
import socket # for socket.gethostbyaddr
import os     # for popen

main:
	conf_path = /etc/wndchrm/wndchrm-queue.conf
	config = ConfigParser.SafeConfigParser()
	config.readfp(io.BytesIO('[conf]\n'+open(conf_path, 'r').read()))
	config.get("conf", "wndchrm_executable")
		wndchrm_executable
		num_workers
		worker_PID_dir
		worker_log
		beanstalkd_tube
		beanstalkd_host
		beanstalkd_port
		beanstalkd_wait_retry
		beanstalkd_retries
	CPU cores & RAM:
		sysctl -n hw.ncpu
		sysctl -n hw.memsize
	FQDN for host:
		socket.getfqdn()
	wndchrm path:
		os.popen('command -v wndchrm').read().strip()
	
	while (1):
		try:
			queue = dependent_queue()
			retries = 0
			queue.run()
		except beanstalkc.SocketError:
			retries += 1
			if (retries > max_retries):
				print_log ("Giving up after "+retries+" retries. Worker exiting.")
				sys.exit(1)
			if queue:
				print_log ("beanstalkd on "+beanstalkd_host+" port "+beanstalkd_port+" not responding. Retry in "+wait_retry+" seconds.")
				queue = None
			time.sleep(wait_retry)
		except Exception e:
			retries += 1
			if (retries > max_retries):
				print_log ("Giving up after "+retries+" retries. Worker exiting.")
				sys.exit(1)
			print_log (str(e))
			if not queue:
				time.sleep(wait_retry)



def which(program):
    import os
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None

def add_job_deps (job, job_dep_tube):
	ndeps = 0
	old_tube = beanstalk.using()
	
	beanstalk.use (job_dep_tube)
	# ...
	# beanstalk.put (dep_job)
	# ndeps += 1
	beanstalk.use (old_tube)
	return ndeps

def run_job (job):
	# ...
	job.delete()

job = beanstalk.reserve()
job_tube = job.stats()['tube']
if (job_tube == deps_tube):
	(job_id, job_dep_tube) = job.body.split("\t")
	if (beanstalk.stats_tube(job_dep_tube)['current-jobs-ready'] > 0):
		job.release()
		if (job_dep_tube not in beanstalk.watching()): beanstalk.watch(job_dep_tube)
		beanstalk.ignore (deps_tube)
	else:
		beanstalk.watch(deps_tube)
		beanstalk.ignore(job_dep_tube)
		if (not beanstalk.stats_tube(job_dep_tube)['current-jobs-reserved'] > 0):
			ready_job = beanstalk.peek (job_id)
			if (ready_job):
				beanstalk.use (jobs_ready_tube)
				beanstalk.put (ready_job.body)
elif (job_tube == jobs_tube):
	job_id = job.stats()['id']
	job_dep_tube = dep_tube+'-'+str(job_id)
	ndeps = add_job_deps (job, job_dep_tube)
	if (ndeps)
		beanstalk.use (dep_tube)
		beanstalk.put (job_id+"\t"+job_dep_tube)
		job.bury()
else:
	run_job (job)
	tube_stats = beanstalk.stats_tube(job_tube)
	if (not tube_stats['current-jobs-reserved'] + tube_stats['current-jobs-ready'] > 0):
		beanstalk.watch(deps_tube)
		beanstalk.ignore(job_tube)
		beanstalk.use (deps_tube)
	
