// WORMfile: Write-once-read-many file for concurrent processes, where only a single process can write.
// synopsis:
// 		WORMfile wf (path);
// 		if (wf.status == WORMfile::WORM_WR) {
//		// File is open for writing with an exclusive write-lock
//		// Write to file using wf.fp() or wf.fd()
//		// Call finish before destructing or the file will be unlinked.
// 			wf.finish();
// 		} else if (wf.status == WORMfile::WORM_BUSY) {
//		// File has active writelock by another process. The file is closed, wf.fp() = NULL and wf.fd() = -1
// 		} else if (wf.status == WORMfile::WORM_RD) {
//		// File is open for reading and has a readlock (no process can write)
// 		} else {
//		// There was an I/O error opening/creating/locking the file.
//		// The file is closed, wf.fp() = NULL and wf.fd() = -1.
//		// wf.status_errno has information about the I/O error (wf.status = WORM_RDLK_ERR || WORM_WRLK_ERR || WORM_IO_ERR || WORM_STALE)
// 		}
//
// The file has four defined states (in order):
//   non-existant: unless the optional read-only parameter is true, the file will be created by open_rw()
//                 in read/write mode, possibly with a write-lock
//   busy: file has read/write permissions and an active write lock
//   stale: file has read/write permissions but no active write lock (process terminated before calling finish())
//          stale files will be opened with a write-lock and truncated barring any locking or I/O errors
//          N.B.: A file may still be invalid/stale even if its not empty, but this is up to the caller to determine.
//   read-only: file has read-only permissions and can be read-locked
// The file may be opened read-only, optionally blocking until a readlock can be established.
// In this case, no write-lock or file creation will be attempted.
//	ex:
//	WORMfile wf (path, true);  // read-only
//	WORMfile wf (path, true, true);  // read-only, wait for a readlock
// wf.status should be checked for WORM_STALE, indicating an empty file without an active write-lock.

#include "WORMfile.h"
#define OPEN_RETRIES 36
#define MAX_WAIT_MULT 8192 // must be smaller than RAND_MAX, which is a minimum of 32767

WORMfile::WORMfile () {
	status = WORM_UNDEF;
	status_errno = 0;
	path = "";
	_fd = -1;
	_fp = NULL;
	_read_mode = def_read_mode;
}
WORMfile::WORMfile (const char *p_path, bool readonly, bool wait) {
	status = WORM_UNDEF;
	status_errno = 0;
	path = p_path;
	_fd = -1;
	_fp = NULL;
	_read_mode = def_read_mode;

	reopen (readonly, wait);
}

WORMfile::~WORMfile() {
	if (_fd > -1) {
		if (status == WORM_WR) {
			unlink (path.c_str());
		}
		if (_fp) fclose(_fp);
		close (_fd);
	}
	status = WORM_UNDEF;
	_fd = -1;
	_fp = NULL;
}

void WORMfile::reopen (bool readonly, bool wait) {

	if (_fd > 0) return;

	// These retries are to resolve the race between the
	// read-only permissions set with finish() and the time between failing a readlock and opening/creating with read/write
	// on glusterfs at least and under heavy load, there are a variety of errors reported
	// by open with read/write/create, not limited to permission errors
	// maximum of 13 retries seen so far with a heavy load on a shared filesystem.
	int k = 1;
	int retries = OPEN_RETRIES;
	status = WORM_UNDEF;
	while ( !(status == WORM_WR || status == WORM_RD || status == WORM_BUSY || status == WORM_ENOENT || status == WORM_STALE) && retries--) {
		if (readonly) open_r (wait);
		else open_rw();
		if (status == WORM_WR || status == WORM_RD || status == WORM_BUSY || status == WORM_ENOENT || status == WORM_STALE) break;

	// exponential backoff ala ethernet protocol: wait for k * 109us, where k is a random number between 0 and 2^fails - 1
	// ethernet protocol uses 51.2 us, we use 109 us.
		k = k << 1;
		if (k > MAX_WAIT_MULT) k = MAX_WAIT_MULT;
		usleep (109 * (rand() % (k - 1)) );
	}
}

int WORMfile::get_read_mode () {
	return (_read_mode);
}

int WORMfile::set_read_mode (int p_read_mode) {
	_read_mode = p_read_mode;
	return (_read_mode);
}

int WORMfile::fd () {
	return (_fd);
};

FILE *WORMfile::fp () {
	if (_fp == NULL && _fd > -1) {
		errno = 0;
		if (status == WORM_RD) _fp = fdopen (_fd, "r");
		else if (status == WORM_WR) _fp = fdopen (_fd, "w"); // in POSIX, no truncation in fdopen
		status_errno = errno;
	} else if (_fd < 0) {
		_fp = NULL;
	}
	return (_fp);
}

void WORMfile::finish (bool reopen) {

	if (_fp < 0) return;

	if (status == WORM_WR) {
		// Change to read-only permissions
		// r--r--r--: mode_t mask = S_IRUSR | S_IRGRP | S_IROTH;
		// These read-only permissions race with the time between failing a readlock and opening/creating with read/write
		// Also, we are changing to read-only while having the file open with write access.
		// No other processes can have the file open for reading or writing, and we are immediately closing the file, so it should be OK
		// This can fail if we are not the file owner, but have write permissions
		//   e.g., file created by different owner, then failed before calling finish().
		// Because we can succeed finishing the file, but fail to change permissions, we can't use the write permission as a check
		// for stale files.
		fchmod (_fd, _read_mode);
	}

	if (_fp) fclose(_fp);
	close (_fd);
	_fp = NULL;
	_fd = -1;
	status = WORM_FINISHED;
	// FIXME: Could do more detailed error checking, but in principle, the caller does not expect to process errors.
	errno = 0;
	// reopen if requested, without waiting.
	if (reopen) open_r ();
}

void WORMfile::open_rw () {
	struct flock fl_w, fl_r;
	struct stat fstat_buf;
	memset (&fl_w, 0, sizeof (struct flock));
	memset (&fl_r, 0, sizeof (struct flock));
	memset (&fstat_buf, 0, sizeof (struct stat));
	fl_w.l_type = F_WRLCK;
	fl_r.l_type = F_RDLCK;
	fl_w.l_whence = fl_r.l_whence = SEEK_SET;
	fl_w.l_start = fl_r.l_start = 0;
	fl_w.l_len = fl_r.l_len = 0;

// The first attempt to get a read-only lock is here in case permissions prevent writing, but
// otherwise things are A-OK.
// Busy files with active write locks are A-OK too, but they are closed.
	open_r();
	if (status == WORM_RD) return;
	if (status == WORM_BUSY) {
		close (_fd);
		_fd = -1;
		status_errno = errno = 0;
		return;
	}

// unable to get a read lock on a non-empty, existing file (otherwise return above)
// Generally, the outcome here is either busy or a write lock (things are not A-OK)
// n.b.: the file is closed here, so we may still end up with a readlock if we lost the race
// have to check either way.
// finish() may have been called on the file by another process, so we may no longer have permissions
// to open with writing.  In this case, the calling method will retry from the beginning with a readlock.
//
// Open the file rw with creation, and try to get a readlock.
// If the readlock succeeds, then return.  Otherwise request a writelock on the readlocked file.
// If the file is busy, return.
	status_errno = errno = 0;
	_fd = open (path.c_str(), O_RDWR | O_CREAT, S_IWUSR |  S_IRUSR);

	if (_fd < 0) {
		// if the file was determined stale by open_r, keep that as the status
		if (status != WORM_STALE) status = WORM_IO_ERR;
		status_errno = errno;
	} else {
		status_errno = errno = 0;

		// Try to get a readlock, return if busy or succeeds
		if (fcntl(_fd, F_SETLK, &fl_r) != -1) {
			int fstat_err = fstat(_fd, &fstat_buf);
			if (fstat_buf.st_size > 0) {
				status_errno = errno = 0;
				status = WORM_RD;
				return;
			} else if (!fstat_err) {
				status = WORM_STALE;
			}

		} else if (errno == EACCES || errno == EAGAIN) {
		// other process has an active readlock, this is not an error
			status_errno = errno = 0;
			status = WORM_BUSY;
			close (_fd);
			_fd = -1;
			status_errno = errno = 0;
			return;
		} else {
		// Some other I/O error establishing the readlock - will try for a writelock anyway
			if (status != WORM_STALE) status = WORM_RDLK_ERR;
		}

		// locked for reading and empty
		// try for a write lock
		// has to be open with write permissions to succeed.
		status_errno = errno = 0;
		if (fcntl(_fd, F_SETLK, &fl_w) != -1) {
			status = WORM_WR;
			ftruncate(_fd,0);
			status_errno = errno = 0;
			return;

		} else if (errno == EACCES || errno == EAGAIN) {
		// busy is not an error, but the file is closed
			status = WORM_BUSY;
			close (_fd);
			_fd = -1;
			status_errno = errno = 0;
			return;

		} else {
		// some I/O error establishing the writelock
			if (status != WORM_STALE) status = WORM_WRLK_ERR;
		}
		
		// if we didn't return yet, then we have an error			
		// close the file if we can't read or write it.
		// possible fixes are to unlink the file and/or fix its permissions, and/or truncate it
		// n.b.: busy files are closed without error
		status_errno = errno;
		close (_fd);
		_fd = -1;
		errno = status_errno; // in case close() changes it
	}
}

void WORMfile::open_r (bool wait) {
	struct flock fl_r;
	struct stat fstat_buf;
	fl_r.l_type = F_RDLCK;
	fl_r.l_whence = SEEK_SET;
	fl_r.l_start = 0;
	fl_r.l_len = 0;
	int set_lock = wait ? F_SETLKW : F_SETLK;

	status_errno = errno = 0;
	_fd = open (path.c_str(), O_RDONLY);

	if (_fd > -1) {
		status_errno = errno = 0;
		if (fcntl(_fd, set_lock, &fl_r) != -1) {
		// We have a readlock, but the file is still invalid if its 0-length
			int fstat_err = fstat(_fd, &fstat_buf);
			if (fstat_buf.st_size > 0) {
				status = WORM_RD;
			} else if (!fstat_err) {
				// A file is invalid (stale) if we can get a readlock on an empty file.
				// If we were called from open_rw, then this will result in a writelock, assuming
				// permissions are OK and no I/O errors.
				// However, if we were called with a read-only option, then we will return with a WORM_STALE status and a closed file.
				close (_fd);
				_fd = -1;
				status = WORM_STALE;
			} else {
				close (_fd);
				_fd = -1;
				status = WORM_IO_ERR;
			}
		} else if (errno == EACCES || errno == EAGAIN) {
			status_errno = errno = 0;
			status = WORM_BUSY;
		} else {
			status_errno = errno;
			close (_fd);
			_fd = -1;
			errno = status_errno; // in case close() changes it
			status = WORM_RDLK_ERR;
		}
	} else {
	// Various I/O errors opening the file
		if (errno == ENOENT) status = WORM_ENOENT;
		else status = WORM_IO_ERR;
		status_errno = errno;
	}
	
	return;
}
int WORMfile::def_read_mode = S_IRUSR | S_IRGRP | S_IROTH;
