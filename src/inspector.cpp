#include "inspector.h"

int main(int argc, char *argv[]) {
	set_new_handler(alloc_error);

	if(argc < 6) {
		print_usage(cout);
		exit(0);
	}
	
	string seqfile;                       // file with sequences
	string exprfile;                      // file with expression data
	string subsetfile;                    // file with subset of sequence names to search
	string scorefile;                     // file with scores

	if(! GetArg2(argc, argv, "-s", seqfile)) {
		cerr << "Please specify sequence file\n\n";
		print_usage(cout);
		exit(0);
	}
	if(! GetArg2(argc, argv, "-o", outfile)) {
		cerr << "Please specify output file\n\n";
		print_usage(cout);
		exit(0);
	}
	
	int search_type = UNDEFINED;
	if(GetArg2(argc, argv, "-ex", exprfile)) {
		search_type = EXPRESSION;
	}
	if(GetArg2(argc, argv, "-su", subsetfile)) {
		search_type = SUBSET;
	}
	if(GetArg2(argc, argv, "-sc", scorefile)) {
		search_type = SCORE;
	}
	if(search_type == UNDEFINED) {
		cerr << "Please specify either an expression data file, a file with a subset of sequence names, or a file with sequence scores.\n\n";
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
	vector<string> seqs, nameset1;
	cerr << "Reading sequence data from '" << seqfile << "'... ";
	get_fasta_fast(seqfile.c_str(), seqs, nameset1);
	cerr << "done.\n";
	ngenes = nameset1.size();
	
	vector<vector <float> > expr;
	vector<string> nameset2;
	vector<string> subset;
	vector<float> scores;
	npoints = 0;
	if(search_type == EXPRESSION) {
		cerr << "Reading expression data from '" << exprfile << "'... ";
		get_expr(exprfile.c_str(), expr, nameset2);
		cerr << "done.\n";
		npoints = expr[0].size();
		nsubset = 0;
	} else if(search_type == SUBSET) {
		cerr << "Reading subset of sequences to search from '" << subsetfile << "'... ";
		get_list(subsetfile.c_str(), subset);
		cerr << "done.\n";
		npoints = 0;
		for(int i = 0; i < ngenes; i++) {
			vector<float> row(0);
			expr.push_back(row);
		}
		nsubset = subset.size();
		sort(subset.begin(), subset.end());
	} else if(search_type == SCORE) {
		cerr << "Reading sequence scores from '" << scorefile << "'... ";
		get_scores(scorefile.c_str(), scores, nameset2);
		cerr << "done.\n";
		npoints = 1;
		nsubset = 0;
	}

	if(search_type == EXPRESSION || search_type == SCORE) {
		if(nameset1.size() != nameset2.size()) {
			if(search_type == EXPRESSION) {
				cerr << "Inconsistent sizes for sequence and expression datasets\n";
			} else {
				cerr << "Inconsistent sizes for sequence and score datasets\n";
			}
			exit(0);
		}

		bool match = true;
		for (int i = 0; i < ngenes; i++) {
			match = match && (nameset1[i] == nameset2[i]);
			if (nameset1[i] != nameset2[i]) {
				cerr << "Row " << i << ": Sequence: '" << nameset1[i] << "'  " << "Data: '" << nameset2[i] << "'\n";
			}
		}
		if(! match) {
			if(search_type == EXPRESSION) {
				cerr << "Inconsistent sequence names for sequence and expression datasets\n";
			} else {
				cerr << "Inconsistent sequence names for sequence and score datasets\n";
			}
			exit(0);
		}
	}

	if(search_type == EXPRESSION) {
		cerr << "Successfully read input files -- dataset size is " << ngenes << " sequences X " << npoints << " timepoints\n";
	} else if(search_type == SUBSET) {
		cerr << "Successfully read input files -- dataset size is " << ngenes << " sequences, with " << nsubset << " to be searched\n";
	} else if(search_type == SCORE) {
		cerr << "Successfully read input files -- dataset size is " << ngenes << " scored sequences\n";
	}

	cerr << "Setting up MotifSearch... ";
	if(! GetArg2(argc, argv, "-numcols", ncol)) ncol = 10;
	if(! GetArg2(argc, argv, "-order", order)) order = 0;
	if(! GetArg2(argc, argv, "-simcut", simcut)) simcut = 0.9;
	MotifSearch* ms;
	if(search_type == EXPRESSION) {
		ms = new MotifSearchExpr(nameset1, seqs, ncol, order, simcut, expr, npoints);
	} else if(search_type == SCORE) {
		ms = new MotifSearchScore(nameset1, seqs, ncol, order, simcut, scores);
	} else {
		ms = new MotifSearchSubset(nameset1, seqs, ncol, order, simcut, subset);
	}
	ms->modify_params(argc, argv);
	ms->set_final_params();
	ms->ace_initialize();
	cerr << "done.\n";
	cerr << "Random seed: " << ms->get_params().seed << '\n';

	if(archive) {
		cerr << "Running in archive mode...\n";
		string archinstr(outfile);
		archinstr.append(".adj.ace");
		ifstream archin(archinstr.c_str());
		if(archin) {
			cerr << "Refreshing from existing archive file " << archinstr << "... ";
			ms->get_archive().read(archin);
			cerr << "done.\n";
		}
		while(true) {
			int found = read_motifs(ms);
			if(found > 0) output(ms);
			sleep(60);
		}
	} else {
		cerr << "Running as worker " << worker << "...\n";
		int nruns = ms->positions_in_search_space()/(ms->get_params().expect * ncol);
		nruns *= ms->get_params().oversample;
		nruns /= ms->get_params().undersample;
		cerr << "Restarts planned: " << nruns << '\n';
		string archinstr(outfile);
		archinstr.append(".adj.ace");
		string lockstr(outfile);
		lockstr.append(".lock");
		for(int j = 1; j <= nruns; j++) {
			if(j == 1 || j % 50 == 0 || search_type == SUBSET) {
				struct flock fl;
				int fd;
				fl.l_type   = F_RDLCK;
				fl.l_whence = SEEK_SET;
				fl.l_start  = 0;
				fl.l_len    = 0;
				fl.l_pid    = getpid();
				fd = open(lockstr.c_str(), O_RDONLY);
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
						cerr << "\t\tRefreshing archive from " << archinstr << "...";
						ms->get_archive().clear();
						ms->get_archive().read(archin);
						archin.close();
						cerr << "done.\n";
					}
					fl.l_type = F_UNLCK;
					fcntl(fd, F_SETLK, &fl);
					close(fd);
					cerr << "\t\tArchive now has " << ms->get_archive().nmots() << " motifs\n";
				}
			}
			cerr << "\t\tSearch restart #" << j << "/" << nruns << "\n";
			ms->search_for_motif(worker, j, outfile);
		}
	}
	delete ms;
	return 0;
}

int read_motifs(MotifSearch* ms) {
	DIR* workdir;
	struct dirent* dirp;
	string filename;
	string extension;
	string match(outfile);
	match.append(".");
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
		if(filename.find(match) == 0 && extension.compare(".mot") == 0) {
			cerr << "Reading from file " << filename << endl;
			// check motif against archive, then move to "mots" directory
			string newname("mots/");
			newname.append(filename.c_str());
			if(ms->consider_motif(filename.c_str())) {
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

void output(MotifSearch* ms) {
	string tmpstr(outfile);
	string outstr(outfile);
	string lockstr(outfile);
	tmpstr.append(".tmp.ace");
	outstr.append(".adj.ace");
	lockstr.append(".lock");
	ofstream tmp(tmpstr.c_str(), ios::trunc);
	print_ace(tmp, ms);
	tmp.close();
	struct flock fl;
	int fd;
	fl.l_type   = F_WRLCK;
	fl.l_whence = SEEK_SET;
	fl.l_start  = 0;
	fl.l_len    = 0;
	fl.l_pid    = getpid();
	fd = open(lockstr.c_str(), O_WRONLY | O_CREAT, 0644);
	while(fcntl(fd, F_SETLK, &fl) == -1) {
		cerr << "Waiting for lock release on archive file...\n";
		sleep(10);
	}
	rename(tmpstr.c_str(), outstr.c_str());
	cerr << "Archive output completed.\n";
	fl.l_type = F_UNLCK;
	fcntl(fd, F_SETLK, &fl);
	close(fd);
}

void print_full_ace(ostream& out, MotifSearch* ms) {
	out << "Parameter values:\n";
	ms->output_params(out);
	out << "\nInput sequences:\n";
	for(unsigned int x = 0; x < ms->names().size(); x++)
		out << "#" << x << '\t' << (ms->names())[x] << endl;
	out << '\n';
	ms->full_output(out);
}

void print_ace(ostream& out, MotifSearch* ms) {
	ms->full_output(out);
}

void print_usage(ostream& fout) {
	fout << "Usage: inspector -s seqfile [-ex exprfile | -su subsetfile | -sc scorefile] -o outputfile (options)\n";
	fout << " Seqfile must be in FASTA format.\n";
	fout << " Exprfile must be in tab-delimited format.\n";
	fout << " Subsetfile has one sequence name per line.\n";
	fout << " Scorefile is tab-delimited, with one sequence name and score per line.";
	fout << "Options:\n";
	fout << " -numcols    \tnumber of columns to align (10)\n";
	fout << " -order      \torder of the background model (3, can be 0 to 5)\n";
	fout << " -simcut     \tsimilarity cutoff for motifs (0.8)\n"; 
	fout << " -expect     \tnumber of sites expected in model (10)\n";
	fout << " -minpass    \tminimum number of non-improved passes in phase 1 (200)\n";
	fout << " -seed       \tset seed for random number generator (time)\n";
	fout << " -undersample\tpossible sites / (expect * numcols * seedings) (1)\n"; 
	fout << " -oversample\t1/undersample (1)\n";
}