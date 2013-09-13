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
//		// File is open for reading and has a readlock (no other process can write)
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

#ifndef __WORMfile_H__
#define __WORMfile_H__
#include <string>
#include <stdio.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <string.h>
#include <errno.h>
#include <stdlib.h>

class WORMfile {
public:

	// file status
	enum Status {
		WORM_UNDEF,     // no I/O operations have been attempted yet
		WORM_BUSY,      // file has active write-lock by another process
		WORM_RD,        // file has readlock
		WORM_WR,        // file has writelock
		WORM_FINISHED,  // finish() was called and the file is closed
		WORM_RDLK_ERR,  // error acquiring readlock other than lock contention
		WORM_WRLK_ERR,  // error acquiring writelock other than lock contention
		WORM_IO_ERR,    // error opening/creating the file
		WORM_STALE,     // The file exists, and can be read-locked, but is empty 
		WORM_ENOENT     // The file does not exist while attempting to open read-only
	};
	
	WORMfile::Status status;
	int status_errno;
	std::string path;
	int read_mode;
	static int def_read_mode;

	WORMfile();
	WORMfile (const char *p_path, bool readonly = false, bool wait = false);

	// calling the destructor when the file is opened for writing causes the file to be unlinked!
	// must call finish() before the destructor to keep the file.
	~WORMfile();
	
	void reopen (bool readonly = false, bool wait = false);

	int get_read_mode ();
	int set_read_mode (int p_read_mode);

	int fd ();
	FILE *fp ();
	
	// This finishes a write operation and prepares the file for reading
	// passing the optional boolean parameter as true, causes the file to be reopened with a read-lock
	// It is not possible to change the open flags of an open file descriptor from O_WRONLY to O_RDONLY,
	// so the file must be closed and reopened.  The reopening of the file may lead to errors,
	// so the status should be checked if reopen is true.
	// N.B.: This must be called before calling the destructor if the file was opened for writing!
	void finish (bool reopen = false);

private:
	int _fd;
	FILE *_fp;
	int _read_mode;

	// open read-only or create write-only using locking and file permissions
	void open_rw ();
	void open_r (bool wait = false);
};

#endif // __WORMfile_H__
