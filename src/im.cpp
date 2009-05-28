#include "im.h"

int main(int argc, char *argv[]) {
	set_new_handler(alloc_error);
	
	if(argc < 6) {
    print_usage(cout);
    exit(0);
  }
	
	string seqfile;                       // file with sequences
	string exprfile;                      // file with expression data
	if ((! GetArg2(argc, argv, "-s", seqfile)) 
				|| (! GetArg2(argc, argv, "-e", exprfile)) 
				|| (! GetArg2(argc, argv, "-o", outfile))) {
		print_usage(cout);
    exit(0);
	}
	
	// Decide mode of running
	archive = false;
	if(! GetArg2(argc, argv, "-worker", worker)) {
		worker = -1;
		archive = true;
	}
	
	// Read parameters
	if(! GetArg2(argc, argv, "-numcols", ncol)) ncol = 10;
  vector<string> seqs, nameset1;
  cerr << "Reading sequence data from '" << seqfile << "'... ";
	get_fasta_fast(seqfile.c_str(), seqs, nameset1);
	cerr << "done.\n";
	
	vector<string> nameset2;
	cerr << "Reading expression data from '" << exprfile << "'... ";
  expr = get_expr(exprfile.c_str(), &npoints, nameset2);
	cerr << "done.\n";
	
	if(nameset1.size() == nameset2.size()) {
		ngenes = nameset1.size();
	} else {
		cerr << "Inconsistent sizes for sequence and expression datasets\n";
		exit(0);
	}
	
	bool match = true;
	for (int i = 0; i < ngenes; i++) {
		match = match && (nameset1[i] == nameset2[i]);
		if (nameset1[i] != nameset2[i]) {
			cerr << "Row " << i << ": Sequence: '" << nameset1[i] << "'  " << "Expression: '" << nameset2[i] << "'\n";
		}
	}
	if(! match) {
		cerr << "Inconsistent gene names for sequence and expression datasets\n";
		exit(0);
	}
	
	cerr << "Successfully read input files -- dataset size is " << ngenes << " genes X " << npoints << " timepoints\n";
	
	cerr << "Setting up SEModel... ";
	SEModel se(seqs, expr, npoints, nameset1, ncol);
	se.modify_params(argc, argv);
	se.set_final_params();
	se.ace_initialize();
	cerr << "done.\n";
	
	if(archive) {
		cerr << "Running in archive mode...\n";
		string archinstr(outfile);
		archinstr.append(".adj.ace");
		ifstream archin(archinstr.c_str());
		if(archin) {
			cerr << "Refreshing from existing archive file " << archinstr << "... ";
			se.get_archive().read(archin);
			cerr << "done.\n";
		}
		while(true) {
			int found = read_motifs(se);
			if(found > 0) output(se);
			sleep(60);
		}
	} else {
		cerr << "Running as worker " << worker << "...\n";
		int nruns = se.total_positions()
								/ se.get_params().expect
								/ ncol
								/ se.get_params().undersample
								* se.get_params().oversample;
		string archinstr(outfile);
		archinstr.append(".adj.ace");
		for(int j = 1; j <= 1; j++) {
			struct flock fl;
			int fd;
			fl.l_type   = F_RDLCK;
			fl.l_whence = SEEK_SET;
			fl.l_start  = 0;
			fl.l_len    = 0;
			fl.l_pid    = getpid();
			fd = open("arch.lock", O_RDONLY);
			if(fd == -1) {
				if(errno != ENOENT)
					cerr << "\t\tUnable to read lock file, error was " << strerror(errno) << "\n";
			} else {
				while(fcntl(fd, F_SETLK, &fl) == -1) {
					cerr << "\t\tWaiting for lock release on archive file... \n";
					sleep(10);
				}
				ifstream archin(archinstr.c_str());
				if(archin) {
					cerr << "\t\tRefreshing archive from " << archinstr << "... ";
					se.get_archive().clear();
					se.get_archive().read(archin);
					archin.close();
					cerr << "done.\n";
				}
				fl.l_type = F_UNLCK;
				fcntl(fd, F_SETLK, &fl);
				close(fd);
				cerr << "\t\tArchive now has " << se.get_archive().nmots() << " motifs\n";
			}
			cerr << "\t\tSearch restart #" << j << "/" << nruns << "\n";
			se.search_for_motif(worker, j);
		}
	}
	return 0;
}

int read_motifs(SEModel& se) {
	DIR* workdir;
	struct dirent* dirp;
	string filename;
	string extension;
	int nfound = 0;
	int nmot = 0;
	unsigned int len = 0;
	unsigned long pos = 0;
	workdir = opendir(".");
	while((dirp = readdir(workdir))) {
		filename = string(dirp->d_name);
		len = filename.length();
		pos = filename.find_last_of('.');
		if(len > 4 && pos != string::npos)
			extension = filename.substr(pos, len - pos);
		else
			extension = "";
		if(extension.compare(".mot") == 0) {
			cerr << "Reading from file " << filename << endl;
			// check motif against archive, then move to "mots" directory
			string newname("mots/");
			newname.append(filename.c_str());
			if(se.consider_motif(filename.c_str())) {
				cerr << "Motif was added\n";
				nmot++;
			} else {
				cerr << "Motif was not added\n";
			}
			cerr << "Moving motif to " << newname << "\n\n";
			rename(filename.c_str(), newname.c_str());
			nfound++;
		}
	}
	closedir(workdir);
	cerr << "Read " << nfound << " motif(s), added " << nmot << " motif(s)\n";
	return nmot;
}

void output(SEModel& se) {
	string tmpstr(outfile);
	string outstr(outfile);
	tmpstr.append(".tmp.ace");
	outstr.append(".adj.ace");
	ofstream tmp(tmpstr.c_str(), ios::trunc);
	print_ace(tmp, se);
	tmp.close();
	struct flock fl;
	int fd;
	fl.l_type   = F_WRLCK;
	fl.l_whence = SEEK_SET;
	fl.l_start  = 0;
	fl.l_len    = 0;
	fl.l_pid    = getpid();
	fd = open("arch.lock", O_WRONLY | O_CREAT, 0644);
	while(fcntl(fd, F_SETLK, &fl) == -1) {
		cerr << "Waiting for lock release on archive file...\n";
		sleep(10);
	}
	rename(tmpstr.c_str(), outstr.c_str());
	fl.l_type = F_UNLCK;
	fcntl(fd, F_SETLK, &fl);
	close(fd);
}

void print_full_ace(ostream& out, SEModel& se) {
	out << "Parameter values:\n";
  se.output_params(out);
  out << "\nInput sequences:\n";
  for(unsigned int x = 0; x < se.names().size(); x++)
		out << "#" << x << '\t' << (se.names())[x] << endl;
  out << '\n';
  se.full_output(out);
}

void print_ace(ostream& out, SEModel& se) {
	out << "Parameter values:\n";
  se.output_params(out);
  out << "\nInput sequences:\n";
  for(unsigned int x = 0; x < se.names().size(); x++)
		out << "#" << x << '\t' << se.names()[x] << endl;
  out << endl;
	se.full_output(out);
}

void print_usage(ostream& fout) {
	fout<<"Usage: im -s seqfile -e exprfile -o outputfile (options)\n";
  fout<<" Seqfile must be in FASTA format.\n";
	fout<<" Exprfile must be in tab-delimited format.\n";
  fout<<"Options:\n";
	fout<<" -numcols    \tnumber of columns to align (10)\n";
  fout<<" -expect     \tnumber of sites expected in model (10)\n";
  fout<<" -minpass    \tminimum number of non-improved passes in phase 1 (200)\n";
  fout<<" -seed       \tset seed for random number generator (time)\n";
  fout<<" -undersample\tpossible sites / (expect * numcols * seedings) (1)\n"; 
  fout<<" -oversample\t1/undersample (1)\n";
}
