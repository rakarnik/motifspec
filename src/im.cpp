#include "im.h"

int main(int argc, char *argv[]) {
	set_new_handler(out_of_memory);

	string seqfile;                       // file with sequences
	string exprfile;                      // file with expression data
	if(argc < 6) {
    print_usage(cout);
    exit(0);
  }
	
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
			cerr << "Row " << i << ": Sequence: '" << nameset1[i] << "'  " << "Expression: '" << nameset2[i] << "'" << endl;
		}
	}
	if(! match) {
		cerr << "Inconsistent gene names for sequence and expression datasets\n";
		exit(0);
	}
	
	cerr << "Successfully read input files -- dataset size is " << ngenes << " genes X " << npoints << " timepoints" << endl;
	
	cerr << "Setting up SEModel... ";
	se.init(seqs, expr, npoints, nameset1, ncol);
	se.modify_params(argc, argv);
	se.set_final_params();
	se.ace_initialize();
	cerr << "done." << endl;
	
	if(archive) {
		cerr << "Running in archive mode..." << endl;
		signal (SIGTERM, final_output);
		while(true) {
			bool found = false;
			found = read_motifs();
			if(found) output();
			sleep(300);
		}
	} else {
		cerr << "Running as worker " << worker << "..." << endl;
		for(int g = 0; g < ngenes; g++)
			se.add_possible(g);
		int nruns = se.possible_positions()
								/ se.get_params().expect
								/ ncol
								/ se.get_params().undersample
								* se.get_params().oversample;
		string archinstr(outfile);
		archinstr.append(".adj.ace");
		string workoutstr(outfile);
		workoutstr.append(".tmp.ace");
		for(int j = 1; j <= nruns; j++) {
			cerr << "\t\tSearch restart #" << j << "/" << nruns << endl;
			se.search_for_motif(worker, j);
			if(j % 100 == 0 && access(archinstr.c_str(), F_OK) == 0) {
				se.get_archive()->clear();
				struct flock fl;
				int fd;
				fl.l_type   = F_RDLCK;
				fl.l_whence = SEEK_SET;
				fl.l_start  = 0;
				fl.l_len    = 0;
				fl.l_pid    = getpid();
				fd = open(archinstr.c_str(), O_RDONLY);
				fcntl(fd, F_SETLKW, &fl);
				ifstream archin(archinstr.c_str());
				se.get_archive()->read(archin);
				archin.close();
				fl.l_type = F_UNLCK;
				fcntl(fd, F_SETLK, F_UNLCK);
				close(fd);
				ofstream workout(workoutstr.c_str());
				print_full_ace(workout, se);
			}
		}
	}
	
	return 0;
}

bool read_motifs() {
	DIR* workdir;
	struct dirent* dirp;
	string filename;
	string extension;
	bool found = false;
	workdir = opendir(".");
	while(dirp = readdir(workdir)) {
		filename = string(dirp->d_name);
		if(filename.length() > 4 && filename.find_last_of('.') != string::npos)
			extension = filename.substr(filename.find_last_of('.'), 4);
		else
			extension = "";
		if(extension.compare(".mot") == 0) {
			// check motif against archive, then delete
			se.consider_motif(filename.c_str());
			remove(filename.c_str());
			found = true;
		}
	}
	closedir(workdir);
	return found;
}

void output() {
	string outstr(outfile);
	outstr.append(".adj.ace");
	struct flock fl;
	int fd;
	fl.l_type   = F_WRLCK;
	fl.l_whence = SEEK_SET;
	fl.l_start  = 0;
	fl.l_len    = 0;
	fl.l_pid    = getpid();
	fd = open(outstr.c_str(), O_WRONLY);
	fcntl(fd, F_SETLKW, &fl);
	ofstream out(outstr.c_str(), ios::trunc);
	print_ace(out, se);
	out.close();
	fl.l_type = F_UNLCK;
	fcntl(fd, F_SETLK, F_UNLCK);
	close(fd);
}

void final_output(int param) {
	read_motifs();
	output();
	exit(0);
}

void print_full_ace(ostream& out, SEModel& se) {
	out << "Parameter values:\n";
  se.output_params(out);
  out << "\nInput sequences:\n";
  for(int x = 0; x < se.names().size(); x++) out << "#" << x << '\t' << (se.names())[x] << endl;
  out << '\n';
  se.full_output(out);
}

void print_ace(ostream& out, SEModel& se) {
	out << "Parameter values:\n";
  se.output_params(out);
  out << "\nInput sequences:\n";
  for(int x = 0; x < se.names().size(); x++) out << "#" << x << '\t' << se.names()[x] << endl;
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
  fout<<" -gcback     \tbackground fractional GC content of input sequence (0.38)\n";
  fout<<" -minpass    \tminimum number of non-improved passes in phase 1 (200)\n";
  fout<<" -seed       \tset seed for random number generator (time)\n";
  fout<<" -undersample\tpossible sites / (expect * numcols * seedings) (1)\n"; 
  fout<<" -oversample\t1/undersample (1)\n";
}
