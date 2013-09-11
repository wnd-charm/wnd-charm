#include "WORMfile.h"
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <errno.h>

#define N_FILES 10
int main (int argc, char **argv) {

	int i=0, count=0, pid;
	char path[255];
	pid = getpid();
	
	while (count < N_FILES) {
		snprintf (path,sizeof(path),"%03d",i++);
		WORMfile wf (path);
		if (wf.status == WORMfile::WORM_WR) {
			count++;
			fprintf (wf.fp(),"%d\n",pid);
			printf ("%5d created file #%3d %s\n",pid, count, path);
			sleep(1);
			fprintf (wf.fp(),"%d\n",pid);
			wf.finish();
		} else if (wf.status == WORMfile::WORM_BUSY) {
			printf ("%5d file %s BUSY\n", pid, path);
		} else if (wf.status == WORMfile::WORM_RD) {
			printf ("%5d file %s READ\n", pid, path);
		} else {
			printf ("%5d file %s ERROR:%s\n", pid, path, strerror(errno));
		}
	}
}
